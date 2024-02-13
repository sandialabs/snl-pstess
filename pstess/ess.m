function ess(i,k,bus,flag,varargin)
% Syntax: ess(i,k,bus,flag) OR
%         ess(i,k,bus,flag,h_sol) OR
%         ess(i,k,bus,flag,vnc_ess)
%
% Purpose: energy storage system,
%          with vectorized computation option
%          NOTE - energy storage buses must be declared as a
%                 non-conforming load bus
%
% Input:   i - energy storage system number (index)
%              if i = 0, vectorized computation
%          k - integer time
%          bus - solved loadflow bus data
%          vnc_ess - voltage at nc_load buses for flag == 4
%                 pass vnc_ess = nan for other cases
%          flag - 0 - initialization
%                 1 - network interface computation (unused)
%                 2 - generator dynamics computation
%                 3 - generator dynamics computation for linearization
%                 4 - iterative network interface computation (for nc_load)
%
% ess_con matrix format
% col   data                                        units
%  1    energy storage system/converter number      integer
%  2    bus number                                  integer
%  3    voltage transducer time constant            sec
%  4    use Pade approximation flag (0 = bypass)    binary
%  5    voltage magnitude time delay (Pade)         sec
%  6    power capacity                              MVA
%  7    energy capacity                             MWh
%  8    converter voltage for current limit         pu
%  9    complex power modulation priority mode      integer
%           1 = active; 2 = reactive; 3 = propor.   --
% 10    converter allowable power factor            pf
% 11    initial state of charge [0,1]               fraction of cap.
% 12    minimum state of charge [0,1]               fraction of cap.
% 13    maximum state of charge [0,1]               fraction of cap.
% 14    converter interface time constant           sec
% 15    active current ramp rate limit              pu/sec on ess base
% 16    reactive current ramp rate limit            pu/sec on ess base
% 17    LVPL breakpoint current (1.22)              pu on ess base
% 18    LVPL zero crossing voltage (0.5)            pu
% 19    LVPL breakpoint voltage (0.9)               pu
% 20    charge/discharge only indicator             integer
%           1 = charge only; 2 = discharge only     --
% 21    charge/discharge efficiency [0,1]           pct as fraction
% 22    scaling on limit ipmin (ke) [0,1]           pu
%
% note 1: SOC limits are deactivated when ess_con(,7) <= 0
% note 2: LVPL is deactivated when ess_con(,17) <= 0
%
% ess_pot matrix format
% col   data                           units
%  1    power capacity                 pu on syst. base
%  2    energy base conversion factor  pu per sec
%  3    current capacity               pu on syst. base
%
% ess states
% var   description
% ess1  transducer for local voltage magnitude
% ess2  Pade approx. for local voltage magnitude
% ess3  active current converter interface
% ess4  reactive current converter interface
% ess5  state of charge integrator

%-----------------------------------------------------------------------------%
% Version history
%
% Version: 3.1
% Author:  Ryan T. Elliott
% Date:    November 2020
% Note:    Reorganized the data structures to separate states associated
%          with the control and device models. Added reactive modulation.
%
% Version: 3.0
% Author:  D. Trudnowski
% Date:    Sep 2020
% Note:    Separated the version 2.0 ess.m into the version 3.0 ess.m and
%          ess_anglcntrl (ess_anglcntrl is now called lsc).
%
% Version: 2.0
% Author:  Ryan T. Elliott
% Date:    September 2019
% Note:    Modified compensation paths and implemented setpoint tracking
%
% Version: 1.0
% Author:  Ryan T. Elliott
% Date:    August 2019
% Note:    Initial version
%-----------------------------------------------------------------------------%

global g;  % declaring struct of global variables

lbnd = 1e-3;  % lower bound to prevent division by zero

if (~isempty(g.ess.ess_con) && i ~= 0)
    % non-vectorized computation
    error('ess: all ess calculations must be vectorized.');
elseif (~isempty(g.ess.ess_con) && i == 0)
    if (flag == 0)  % initialization
        % vectorized computation (initialization)

        % busnum -- ess bus number vector
        busnum = g.bus.bus_int(g.ess.ess_con(:,2));

        % ess_vmag(:,1) -- initial ess bus voltage magnitude (pu)
        g.ess.ess_vmag(:,1) = bus(busnum,2);

        % ess_sinj(:,1) -- initial complex power inj. (pu on syst. base)
        g.ess.ess_sinj(:,1) = -(bus(busnum,6) + 1j*bus(busnum,7));
        i_inj = g.ess.ess_sinj(:,1)./max(g.ess.ess_vmag(:,1),lbnd);

        % ess_pot(:,1) -- power capacity (pu on syst. base)
        % ess_pot(:,2) -- energy base conversion factor
        % ess_pot(:,3) -- current capacity (pu on syst. base)
        g.ess.ess_pot(:,1) = g.ess.ess_con(:,6)/g.sys.basmva;
        g.ess.ess_pot(:,2) = g.sys.basmva./(g.ess.ess_con(:,7)*3600);
        g.ess.ess_pot(:,3) = g.ess.ess_pot(:,1)./max(g.ess.ess_con(:,8),lbnd);

        % translating the vdl curves into the system base
        for ii = 1:g.ess.n_ess
            vdl = g.ess.vdl{ii};
            vdl.ip = vdl.ip.*g.ess.ess_pot(:,1);
            vdl.iq = vdl.iq.*g.ess.ess_pot(:,1);
            g.ess.vdl{ii} = vdl;
        end

        % ess_soc(:,1) -- initial ess state of charge (pct of ess capacity)
        g.ess.ess_soc(:,1) = g.ess.ess_con(:,11);

        g.ess.ess1(:,1) = g.ess.ess_vmag(:,1);        % voltage transducer
        g.ess.ess2(:,1) = g.ess.ess_vmag(:,1);        % Pade approximation
        g.ess.ess3(:,1) = real(i_inj);                % active curr. interface
        g.ess.ess4(:,1) = imag(i_inj);                % reactive curr. interface
        g.ess.ess5(:,1) = 0.0;                        % ess SOC integrator

        g.ess.dess1(:,1) = 0.0;
        g.ess.dess2(:,1) = 0.0;
        g.ess.dess3(:,1) = 0.0;
        g.ess.dess4(:,1) = 0.0;
        g.ess.dess5(:,1) = 0.0;

        mask = ~ismember(g.ess.ess_con(:,9),[1,2,3]);
        if any(mask)
            estr = '\ness: invalid modulation priority mode selected at ';
            estr = [estr, 'bus %0.0f.'];
            error(sprintf(estr,bus(busnum(mask),1)));
        end

        mask = ((g.ess.ess_con(:,10) < 0) | (g.ess.ess_con(:,10) > 1));
        if any(mask)
            estr = '\ness: impermissible converter allowable power factor at ';
            estr = [estr, 'bus %0.0f.'];
            error(sprintf(estr,bus(busnum(mask),1)));
        end

        mask = (g.ess.ess_con(:,11) < g.ess.ess_con(:,12));
        if any(mask)
            estr = '\ness: initial SOC below minimum at ';
            estr = [estr, 'bus %0.0f.'];
            error(sprintf(estr,bus(busnum(mask),1)));
        end

        mask = (g.ess.ess_con(:,11) > g.ess.ess_con(:,13));
        if any(mask)
            estr = '\ness: initial SOC exceeds maximum at ';
            estr = [estr, 'bus %0.0f.'];
            error(sprintf(estr,bus(busnum(mask),1)));
        end

        mask = (g.ess.ess_con(:,18) >= g.ess.ess_con(:,19));
        if any(mask)
            estr = '\ness: impermissible LVPL parameters at ';
            estr = [estr, 'bus %0.0f, zerox >= brkpt.'];
            error(sprintf(estr,bus(busnum(mask),1)));
        end

        mask = ((g.ess.ess_con(:,21) <= 0) | (g.ess.ess_con(:,21) > 1));
        if any(mask)
            estr = '\ness: impermissible charge/discharge efficiency at ';
            estr = [estr, 'bus %0.0f.'];
            error(sprintf(estr,bus(busnum(mask),1)));
        end

        mask = (abs(g.ess.ess_sinj(:,1)) > g.ess.ess_pot(:,1));
        if any(mask)
            estr = '\ness: initial apparent power exceeds maximum at ';
            estr = [estr, 'bus %0.0f.'];
            error(sprintf(estr,bus(busnum(mask),1)));
        end

        mask = (abs(i_inj) > g.ess.ess_pot(:,3));
        if any(mask)
            estr = '\ness: initial current exceeds maximum at ';
            estr = [estr, 'bus %0.0f.'];
            error(sprintf(estr,bus(busnum(mask),1)));
        end
    end

    if (flag == 1)  % network interface computation
        % iterative interface calculation required, done in flag == 4
    end

    if (flag == 2 || flag == 3)  % dynamics calculation
        % vectorized computation (dynamics)

        % busnum -- ess bus number vector
        busnum = g.bus.bus_int(g.ess.ess_con(:,2));

        %-----------------------------------------------------------------%
        % voltage magnitude transducer

        % ess_vmag -- raw ess bus voltage magnitude
        g.ess.ess_vmag(:,k) = abs(g.bus.bus_v(busnum,k));

        g.ess.dess1(:,k) = (g.ess.ess_vmag(:,k) - g.ess.ess1(:,k)) ...
                           ./max(g.ess.ess_con(:,3),lbnd);

        mask = (g.ess.ess_con(:,3) < lbnd);
        if any(mask)                                  % integrator bypass
            g.ess.dess1(mask,k) = 0.0;
            g.ess.ess1(mask,k) = g.ess.ess_vmag(mask,k);
        end

        %-----------------------------------------------------------------%
        % time delay via Pade approximation

        % input is the local voltage measurement ess1
        g.ess.dess2(:,k) = (g.ess.ess1(:,k) - g.ess.ess2(:,k)) ...
                           ./max(g.ess.ess_con(:,5)/2,lbnd);

        mask = (g.ess.ess_con(:,5)/2 < lbnd);
        if any(mask)                                  % integrator bypass
            g.ess.dess2(mask,k) = 0.0;
            g.ess.ess2(mask,k) = g.ess.ess1(mask,k);  % no delay
        end

        g.ess.ess_vmag_pade(:,k) = 2*g.ess.ess2(:,k) - g.ess.ess1(:,k);

        v_term = g.ess.ess1(:,k);                     % local voltage measurement
        mask = (g.ess.ess_con(:,4) ~= 0);             % check whether to use Pade
        if any(mask)
            v_term(mask) = g.ess.ess_vmag_pade(mask,k);
        end

        %-----------------------------------------------------------------%
        % determining the active power order

        if (g.reec.n_reec ~= 0)                       % convert to power order
            g.ess.ess_scmd(g.reec.ess_idx,k) = ...
                g.ess.ess_scmd(g.reec.ess_idx,k) ...
                .*max(v_term(g.reec.ess_idx),lbnd) ...
                - g.ess.ess_sinj(g.reec.ess_idx,1);
        end

        p_ord = real(g.ess.ess_scmd(:,k) + g.ess.ess_sinj(:,1) ...
                     + g.ess.ess_sig(:,k) + g.ess.ess_dsig(:,k));

        p_ord_lb = -g.ess.ess_pot(:,1);
        p_ord_ub = g.ess.ess_pot(:,1);

        if (g.reec.n_reec ~= 0)                       % don't bind for reec
            p_ord_lb(g.reec.ess_idx) = -100*g.ess.ess_pot(g.reec.ess_idx,1);
            p_ord_ub(g.reec.ess_idx) = 100*g.ess.ess_pot(g.reec.ess_idx,1);
        end

        mask = (g.ess.ess_con(:,20) == 1);            % charge only mode
        if any(mask)
            p_ord_ub(mask) = real(g.ess.ess_sinj(mask,1));
        end

        mask = (g.ess.ess_con(:,20) == 2);            % discharge only mode
        if any(mask)
            p_ord_lb(mask) = real(g.ess.ess_sinj(mask,1));
        end

        % determining the soc-dependent bounds
        g.ess.ess_soc(:,k) = g.ess.ess_soc(:,1) + g.ess.ess5(:,k);

        mask = (g.ess.ess_con(:,7) > 0) ...
               & (g.ess.ess_soc(:,k) + 1e-6 >= g.ess.ess_con(:,13));

        if any(mask)
            p_ord_lb(mask) = 0.0;
        end

        mask = (g.ess.ess_con(:,7) > 0) ...
               & (g.ess.ess_soc(:,k) - 1e-6 <= g.ess.ess_con(:,12));

        if any(mask)
            p_ord_ub(mask) = 0.0;
        end

        p_ord = max(p_ord,p_ord_lb);                  % active power limits
        p_ord = min(p_ord,p_ord_ub);

        %-----------------------------------------------------------------%
        % determining the reactive power order

        q_ord = imag(g.ess.ess_scmd(:,k) + g.ess.ess_sinj(:,1) ...
                     + g.ess.ess_sig(:,k) + g.ess.ess_dsig(:,k));

        q_ord_lb = -g.ess.ess_pot(:,1);
        q_ord_ub = g.ess.ess_pot(:,1);

        if (g.reec.n_reec ~= 0)                       % don't bind for reec
            q_ord_lb(g.reec.ess_idx) = -100*g.ess.ess_pot(g.reec.ess_idx,1);
            q_ord_ub(g.reec.ess_idx) = 100*g.ess.ess_pot(g.reec.ess_idx,1);
        end

        q_ord = max(q_ord,q_ord_lb);                  % reactive power limits
        q_ord = min(q_ord,q_ord_ub);

        %-----------------------------------------------------------------%
        % real and reactive current interface

        if (nargin ~= 5)
            estr = 'ess: incorrect number of arguments passed to function ';
            estr = [estr, '(flag == 2 || flag == 3).'];
            error(estr);
        else
            h_sol = varargin{1};
        end

        [ip_ord,iq_ord] = ess_sat(k,p_ord,q_ord,v_term,flag,h_sol);

        g.ess.ess_iord(:,k) = ip_ord + 1j*iq_ord;

        g.ess.dess3(:,k) = (ip_ord - g.ess.ess3(:,k)) ...
                           ./max(g.ess.ess_con(:,14),lbnd);

        g.ess.dess4(:,k) = (iq_ord - g.ess.ess4(:,k)) ...
                           ./max(g.ess.ess_con(:,14),lbnd);

        mask = (g.ess.ess_con(:,14) < lbnd);
        if any(mask)                                  % integrator bypass
            g.ess.dess3(mask,k) = 0.0;
            g.ess.ess3(mask,k) = ip_ord(mask);        % active current
            %
            g.ess.dess4(mask,k) = 0.0;
            g.ess.ess4(mask,k) = iq_ord(mask);        % reactive current
        end

        %-----------------------------------------------------------------%
        % state of charge tracking

        p_ch = -min(real(g.ess.ess_sinj(:,k)),0);
        p_dis = max(real(g.ess.ess_sinj(:,k)),0);

        tmp_ch = p_ch.*g.ess.ess_con(:,21) - p_dis./g.ess.ess_con(:,21);

        g.ess.dess5(:,k) = g.ess.ess_pot(:,2).*tmp_ch;
    end

    if (flag == 4)  % current calculation for nc_load
        % vectorized computation (iterative interface calculation)

        % busnum -- ess bus number vector
        busnum = g.bus.bus_int(g.ess.ess_con(:,2));

        if (nargin ~= 5)
            estr = 'ess: incorrect number of arguments passed to function ';
            estr = [estr, '(flag == 4).'];
            error(estr);
        else
            vnc_ess = varargin{1};
            mask = isnan(vnc_ess);
            if any(mask)
                estr = '\ness: undefined voltage during nc_load calc. at ';
                estr = [estr, 'bus %0.0f.'];
                error(sprintf(estr,bus(busnum(mask),1)));
            end
        end

        p_inj = max(abs(vnc_ess).*g.ess.ess3(:,k),-g.ess.ess_pot(:,1));
        p_inj = min(p_inj,g.ess.ess_pot(:,1));
        q_inj = abs(vnc_ess).*g.ess.ess4(:,k);

        % to limit reactive injection based on the mva rating use:
        % q_inj = max(abs(vnc_ess).*g.ess.ess4(:,k),-g.ess.ess_pot(:,1));
        % q_inj = min(q_inj,g.ess.ess_pot(:,1));

        [ip_inj,iq_inj] = ess_sat(k,p_inj,q_inj,abs(vnc_ess),flag);

        % complex power *into* system
        g.ess.ess_sinj(:,k) = abs(vnc_ess).*(ip_inj + 1j*iq_inj);

        % complex power modulation
        tmp_smod = g.ess.ess_sinj(:,k) - g.ess.ess_sinj(:,1);
        vnc_ess = max(abs(vnc_ess),lbnd).*exp(1j*angle(vnc_ess));
        g.ess.ess_cur(:,k) = conj(tmp_smod./vnc_ess);
    end
end

end  % function end

% eof
