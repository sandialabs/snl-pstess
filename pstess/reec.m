function reec(i,k,bus,flag,varargin)
% Syntax: reec(i,k,bus,flag) OR
%         reec(i,k,bus,flag,h_sol)
%
% Purpose: Renewable energy electrical control
%          with vectorized computation
%          NOTE - there must be a matching ess instance to which
%                 reec connects, e.g., a storage system or PV plant
%
% Input:   i - reec number (index)
%              if i = 0, vectorized computation
%          k - integer time
%          bus - solved loadflow bus data
%          flag - 0 - initialization
%                 1 - network interface computation
%                 2 - generator dynamics computation
%                 3 - generator dynamics computation for linearization
%
% Note:    The reec power base (MVA) is inherited from the ess to which
%          the controller is connected.
%
% reec_con matrix format
%
% col   data                                                units
% 1     MVA base (mvab), inherited from ess                 MVA
% 2     bus number                                          integer
% 3     voltage dip logic upper limit (vdiph)               pu
% 4     voltage dip logic lower limit (vdipl)               pu
% 5     transducer time constant (trv)                      sec
% 6     lower deadband in voltage error (dbd1)              pu
% 7     upper deadband in voltage error (dbd2)              pu
% 8     reactive current injection gain (kqv)               pu
% 9     maximum limit of reactive current injection (iqh1)  pu
% 10    minimum limit of reactive current injection (iql1)  pu
% 11    value at which iqcmd is held for "thldq" seconds    pu
%           following a voltage dip (iqfrz)
% 12    post voltage dip hold time for iqcmd (thldq)        sec
% 13    post voltage dip hold time for ipcmd (thldp)        sec
% 14    active power transducer time constant (tp)          sec
% 15    reactive power transducer time constant (tq)        sec
% 16    reactive power voltage maximum limit (qvmax)        pu
% 17    reactive power voltage minimum limit (qvmin)        pu
% 18    voltage control maximum limit (vmax)                pu
% 19    voltage control minimum limit (vmin)                pu
% 20    anti-windup flag for Q, V control (=1 anti-windup)  binary
% 21    proportional gain (kqp)                             pu
% 22    integral gain (kqi)                                 pu
% 23    proportional gain (kvp)                             pu
% 24    integral gain (kvi)                                 pu
% 25    time constant (tiq)                                 sec
% 26    up ramp rate on power reference (dpmax)             pu/sec
% 27    down ramp rate on power reference (dpmin)           pu/sec
% 28    maximum power reference (pmax)                      pu
% 29    minimum power reference (pmin)                      pu
% 30    time constant (tpord)                               pu
% 31    voltage compensation flag (vcmpflag)                integer
%           1 = current compensation; 0 = reactive droop
% 32    power factor flag (pfflag)                          integer
%           1 = power factor control; 0 = Q control
% 33    voltage control flag (vflag)                        integer
%           1 = Q control; 0 = Voltage control
% 34    reactive power control flag (qflag)                 integer
%           1 = Voltage/Q control;
%           0 = Constant power factor or Q control
% 35    power reference flag (pflag) [UNUSED]               integer
%           1 = Pref*speed; 0 = Pref
% 36    current compensation resistance on gen. base (rc)   pu
% 37    current compensation reactance on gen. base (xc)    pu
% 38    filter time constant for voltage measurement (tr1)  sec
% 39    reactive current compensation gain (kc)             pu
% 40    maximum converter blocking voltage (vblkh)          pu
% 41    minimum converter blocking voltage (vblkl)          pu
% 42    voltage blocking delay (tblkdelay)                  sec
% 43    time constant for active current command filter     sec
% 44    time constant for reactive current command filter   sec

% reec states
%   var     description
%   reec1   terminal voltage transducer
%   reec2   active power transducer
%   reec3   reactive power transducer
%   reec4   reactive PI loop integrator
%   reec5   voltage PI loop integrator
%   reec6   reactive current order filter (PI bypass)
%   reec7   active power order filter
%   reec8   voltage compensation filter (vcmpflag)
%   reec9   active current command filter, ipcmd
%   reec10  reactive current command filter, iqcmd

%-----------------------------------------------------------------------------%
% Version history
%
% Version: 1.0
% Author:  Ryan T. Elliott and Sam Ojetola
% Date:    November 2022
% Note:    Initial version, based on REEC_D
%-----------------------------------------------------------------------------%

global g;  % declaring struct of global variables

lbnd = 1e-3;  % lower bound to prevent division by zero

if (~isempty(g.reec.reec_con) && i ~= 0)
    % non-vectorized computation
    error('reec: all reec calculations must be vectorized.');
elseif (~isempty(g.reec.reec_con) && i == 0)
    if (flag == 0)  % initialization
        % vectorized computation (initialization)

        % busnum -- reec bus number vector
        busnum = g.bus.bus_int(g.reec.reec_con(:,2));

        % reec_pot(:,1) -- factor to convert from system to ess mva base
        g.reec.reec_pot(:,1) = g.sys.basmva./g.ess.ess_con(g.reec.ess_idx,6);

        % sgen -- injected complex power pu on generator base
        sgen = -(bus(busnum,6) + 1j*bus(busnum,7)).*g.reec.reec_pot(:,1);

        % icmd -- commanded current pu on generator base
        g.reec.icmd(:,1) = sgen./max(bus(busnum,2),lbnd);

        % iqmin, iqmax -- PI loop current limits
        g.reec.iqmax = (g.ess.ess_con(g.reec.ess_idx,6)/g.sys.basmva) ...
                       ./max(g.ess.ess_con(g.reec.ess_idx,8),lbnd);
        g.reec.iqmin = -g.reec.iqmax;

        % vcomp -- compensated voltage for cascaded PI loop
        vcomp = bus(busnum,2) + g.reec.reec_con(:,39).*imag(sgen);

        vcmpflag = (g.reec.reec_con(:,31) == 1);
        if any(vcmpflag)
            Vt = bus(busnum,2).*exp(1j*bus(busnum,3)*pi/180);
            It = g.reec.icmd(:,1);
            Zc = g.reec.reec_con(:,36) + 1j*g.reec.reec_con(:,37);

            vcomp(vcmpflag) = abs(Vt(vcmpflag) - It(vcmpflag).*Zc(vcmpflag));
        end

        % state variables
        g.reec.reec1(:,1) = bus(busnum,2);               % terminal voltage
        g.reec.reec2(:,1) = real(sgen);                  % active power
        g.reec.reec3(:,1) = imag(sgen);                  % reactive power
        g.reec.reec4(:,1) = vcomp;                       % reactive PI integrator
        g.reec.reec5(:,1) = imag(sgen)./bus(busnum,2);   % voltage PI integrator
        g.reec.reec6(:,1) = imag(sgen)./bus(busnum,2);   % lower PI bypass
        g.reec.reec7(:,1) = real(sgen);                  % active power order
        g.reec.reec8(:,1) = vcomp;                       % compensated voltage
        g.reec.reec9(:,1) = real(sgen)./bus(busnum,2);   % active current command
        g.reec.reec10(:,1) = imag(sgen)./bus(busnum,2);  % reactive current command

        % state derivatives
        g.reec.dreec1(:,1) = 0.0;
        g.reec.dreec2(:,1) = 0.0;
        g.reec.dreec3(:,1) = 0.0;
        g.reec.dreec4(:,1) = 0.0;
        g.reec.dreec5(:,1) = 0.0;
        g.reec.dreec6(:,1) = 0.0;
        g.reec.dreec7(:,1) = 0.0;
        g.reec.dreec8(:,1) = 0.0;
        g.reec.dreec9(:,1) = 0.0;
        g.reec.dreec10(:,1) = 0.0;

        % control references
        g.reec.pref = diag(real(sgen))*ones(size(g.reec.pref));
        g.reec.qref = diag(imag(sgen))*ones(size(g.reec.qref));
        g.reec.pfaref = atan2(imag(sgen),real(sgen));

        qlim = g.reec.qref(:,1);

        pfflag = (g.reec.reec_con(:,32) == 1);
        if any(pfflag)
            qlim(pfflag) = g.reec.reec2(pfflag,1).*tan(g.reec.pfaref(pfflag));
        end

        qlim = min(qlim,g.reec.reec_con(:,16));
        qlim = max(qlim,g.reec.reec_con(:,17));

        g.reec.vref0 = g.reec.reec1(:,1);
        g.reec.vref1 = g.reec.reec8(:,1) - qlim;

        %-----------------------------------------------------------------%
        % parameter checks

        mask = (g.reec.reec_con(:,3) < g.reec.reec_con(:,4));
        if any(mask)
            estr = '\nreec: invalid voltage dip logic limits selected at ';
            estr = [estr, 'bus %0.0f '];
            error(sprintf(estr,bus(busnum(mask),1)));
        end

        mask = (g.reec.reec_con(:,6) > 0);
        if any(mask)
            estr = '\nreec: lower deadband voltage is positive at ';
            estr = [estr, 'bus %0.0f '];
            error(sprintf(estr,bus(busnum(mask),1)));
        end

        mask = (g.reec.reec_con(:,7) < g.reec.reec_con(:,6));
        if any(mask)
            estr = '\nreec: invalid deadband voltages selected at ';
            estr = [estr, 'bus %0.0f '];
            error(sprintf(estr,bus(busnum(mask),1)));
        end

        mask = (g.reec.reec_con(:,9) < g.reec.reec_con(:,10));
        if any(mask)
            estr = '\nreec: invalid current injection limits selected at ';
            estr = [estr, 'bus %0.0f '];
            error(sprintf(estr,bus(busnum(mask),1)));
        end

        mask = (g.reec.reec_con(:,13) < 0);
        if any(mask)
            estr = '\nreec: invalid post voltage dip time for ipcmd selected at ';
            estr = [estr, 'bus %0.0f '];
            error(sprintf(estr,bus(busnum(mask),1)));
        end

        mask = (g.reec.reec_con(:,16) < g.reec.reec_con(:,17));
        if any(mask)
            estr = '\nreec: invalid reactive power voltage limits selected at ';
            estr = [estr, 'bus %0.0f '];
            error(sprintf(estr,bus(busnum(mask),1)));
        end

        mask = (g.reec.reec_con(:,18) < g.reec.reec_con(:,19));
        if any(mask)
            estr = '\nreec: invalid voltage control limits selected at ';
            estr = [estr, 'bus %0.0f '];
            error(sprintf(estr,bus(busnum(mask),1)));
        end

        mask = (g.reec.reec_con(:,27) > 0);
        if any(mask)
            estr = '\nreec: downward ramp rate on power reference is positive at ';
            estr = [estr, 'bus %0.0f '];
            error(sprintf(estr,bus(busnum(mask),1)));
        end

        mask = (g.reec.reec_con(:,26) < g.reec.reec_con(:,27));
        if any(mask)
            estr = '\nreec: invalid power reference ramp rate limits selected at ';
            estr = [estr, 'bus %0.0f '];
            error(sprintf(estr,bus(busnum(mask),1)));
        end

        mask = (g.reec.reec_con(:,28) < g.reec.reec_con(:,29));
        if any(mask)
            estr = '\nreec: invalid power reference limits selected at ';
            estr = [estr, 'bus %0.0f '];
            error(sprintf(estr,bus(busnum(mask),1)));
        end

        mask = (g.reec.reec_con(:,40) < g.reec.reec_con(:,41));
        if any(mask)
            estr = '\nreec: invalid converter blocking voltages selected at ';
            estr = [estr, 'bus %0.0f '];
            error(sprintf(estr,bus(busnum(mask),1)));
        end

        mask = (g.reec.reec_con(:,42) < 0);
        if any(mask)
            estr = '\nreec: invalid voltage blocking time delay selected at ';
            estr = [estr, 'bus %0.0f '];
            error(sprintf(estr,bus(busnum(mask),1)));
        end

        %-----------------------------------------------------------------%
        % flag checks

        mask = ~ismember(g.reec.reec_con(:,31),[0,1]);
        if any(mask)
            estr = '\nreec: invalid voltage compensation flag selected at ';
            estr = [estr, 'bus %0.0f '];
            error(sprintf(estr,bus(busnum(mask),1)));
        end

        mask = ~ismember(g.reec.reec_con(:,32),[0,1]);
        if any(mask)
            estr = '\nreec: invalid power factor flag selected at ';
            estr = [estr, 'bus %0.0f '];
            error(sprintf(estr,bus(busnum(mask),1)));
        end

        mask = ~ismember(g.reec.reec_con(:,33),[0,1]);
        if any(mask)
            estr = '\nreec: invalid voltage control flag selected at ';
            estr = [estr, 'bus %0.0f '];
            error(sprintf(estr,bus(busnum(mask),1)));
        end

        mask = ~ismember(g.reec.reec_con(:,34),[0,1]);
        if any(mask)
            estr = '\nreec: invalid reactive power control flag selected at ';
            estr = [estr, 'bus %0.0f '];
            error(sprintf(estr,bus(busnum(mask),1)));
        end

        mask = ~ismember(g.reec.reec_con(:,35),[0,1]);
        if any(mask)
            estr = '\nreec: invalid power reference flag selected at ';
            estr = [estr, 'bus %0.0f (CAUTION UNUSED) '];
            error(sprintf(estr,bus(busnum(mask),1)));
        end
    end

    if (flag == 1)  % network interface computation
        % interface calculation handled by ess
    end

    if (flag == 2 || flag == 3)  % dynamics calculation
        % vectorized computation (dynamics calculation)

        % busnum -- ess bus number vector
        busnum = g.bus.bus_int(g.reec.reec_con(:,2));

        %-----------------------------------------------------------------%
        % voltage magnitude transducer

        g.reec.dreec1(:,k) = (abs(g.bus.bus_v(busnum,k)) - g.reec.reec1(:,k)) ...
                             ./max(g.reec.reec_con(:,5),lbnd);

        mask = (g.reec.reec_con(:,5) < lbnd);
        if any(mask)                                  % integrator bypass
            g.reec.dreec1(mask,k) = 0.0;
            g.reec.reec1(mask,k) = abs(g.bus.bus_v(busnum(mask),k));
        end

        % vt_filt > vdiph or vt_filt < vdipl
        voltage_dip = (g.reec.reec1(:,k) > g.reec.reec_con(:,3)) ...
                      | (g.reec.reec1(:,k) < g.reec.reec_con(:,4));

        %-----------------------------------------------------------------%
        % active power transducer

        pgen = real(g.ess.ess_sinj(g.reec.ess_idx,k)) ...
               .*g.reec.reec_pot(:,1);

        g.reec.dreec2(:,k) = (pgen - g.reec.reec2(:,k)) ...
                             ./max(g.reec.reec_con(:,14),lbnd);

        mask = (g.reec.reec_con(:,14) < lbnd);
        if any(mask)                                  % integrator bypass
            g.reec.dreec2(mask,k) = 0.0;
            g.reec.reec2(mask,k) = pgen(mask);
        end

        %-----------------------------------------------------------------%
        % reactive power transducer

        qgen = imag(g.ess.ess_sinj(g.reec.ess_idx,k)) ...
               .*g.reec.reec_pot(:,1);

        g.reec.dreec3(:,k) = (qgen - g.reec.reec3(:,k)) ...
                             ./max(g.reec.reec_con(:,15),lbnd);

        mask = (g.reec.reec_con(:,15) < lbnd);
        if any(mask)                                  % integrator bypass
            g.reec.dreec3(mask,k) = 0.0;
            g.reec.reec3(mask,k) = qgen(mask);
        end

        %-----------------------------------------------------------------%
        % local voltage control loop

        verr = g.reec.vref0 - g.reec.reec1(:,k) + real(g.reec.reec_sig(:,k));

        mask = (verr >= g.reec.reec_con(:,6)) & (verr <= g.reec.reec_con(:,7));
        if any(mask)
            verr(mask) = 0.0;
        end

        mask = (verr > 0);
        if any(mask)
            verr(mask) = verr(mask) - g.reec.reec_con(mask,7);
        end

        mask = (verr < 0);
        if any(mask)
            verr(mask) = verr(mask) - g.reec.reec_con(mask,6);
        end

        iqv = g.reec.reec_con(:,8).*verr;
        iqinj = min(iqv,g.reec.reec_con(:,9));
        iqinj = max(iqinj,g.reec.reec_con(:,10));

        %-----------------------------------------------------------------%
        % voltage compensation path

        vcomp = abs(g.bus.bus_v(busnum,k)) + g.reec.reec_con(:,39).*qgen;

        vcmpflag = (g.reec.reec_con(:,31) == 1);
        if any(vcmpflag)
            Vt = g.bus.bus_v(busnum,k);
            It = (pgen + 1j*qgen)./max(abs(Vt),0.01);
            Zc = g.reec.reec_con(:,36) + 1j*g.reec.reec_con(:,37);

            vcomp(vcmpflag) = abs(Vt(vcmpflag) - It(vcmpflag).*Zc(vcmpflag));
        end

        g.reec.dreec8(:,k) = (vcomp - g.reec.reec8(:,k)) ...
                             ./max(g.reec.reec_con(:,38),lbnd);

        mask = (g.reec.reec_con(:,38) < lbnd);
        if any(mask)                                  % integrator bypass
            g.reec.dreec8(mask,k) = 0.0;
            g.reec.reec8(mask,k) = vcomp(mask);
        end

        %-----------------------------------------------------------------%
        % first stage PI loop

        qlim = g.reec.qref(:,k);

        pfflag = (g.reec.reec_con(:,32) == 1);
        if any(pfflag)
            qlim(pfflag) = g.reec.reec2(pfflag,k).*tan(g.reec.pfaref(pfflag));
        end

        qlim = min(qlim,g.reec.reec_con(:,16));
        qlim = max(qlim,g.reec.reec_con(:,17));

        qerr = qlim - g.reec.reec3(:,k) + imag(g.reec.reec_sig(:,k));

        g.reec.dreec4(:,k) = g.reec.reec_con(:,22).*qerr;

        y_piq = g.reec.reec_con(:,21).*qerr + g.reec.reec4(:,k);

        % anti-windup limits
        mask = (y_piq > g.reec.reec_con(:,18));
        if any(mask)
            y_piq(mask) = g.reec.reec_con(mask,18);

            anti_windup = mask & (g.reec.reec_con(:,20) == 1) ...
                          & (g.reec.dreec4(:,k) > 0);

            if any(anti_windup)
                g.reec.dreec4(anti_windup,k) = 0;
            end
        end

        mask = (y_piq < g.reec.reec_con(:,19));
        if any(mask)
            y_piq(mask) = g.reec.reec_con(mask,19);

            anti_windup = mask & (g.reec.reec_con(:,20) == 1) ...
                          & (g.reec.dreec4(:,k) < 0);

            if any(anti_windup)
                g.reec.dreec4(anti_windup,k) = 0;
            end
        end

        g.reec.dreec4(voltage_dip,k) = 0;

        %-----------------------------------------------------------------%
        % second stage PI loop

        vlim = qlim + g.reec.vref1;

        vflag = (g.reec.reec_con(:,33) == 1);
        if any(vflag)
            vlim(vflag) = y_piq(vflag);
        end

        vlim = min(vlim,g.reec.reec_con(:,18));
        vlim = max(vlim,g.reec.reec_con(:,19));

        g.reec.dreec5(:,k) = g.reec.reec_con(:,24).*(vlim - g.reec.reec8(:,k));

        y_piv = g.reec.reec_con(:,23).*(vlim - g.reec.reec8(:,k)) ...
                + g.reec.reec5(:,k);

        % anti-windup limits
        mask = (y_piv > g.reec.iqmax);
        if any(mask)
            y_piv(mask) = g.reec.iqmax(mask);

            anti_windup = mask & (g.reec.reec_con(:,20) == 1) ...
                          & (g.reec.dreec5(:,k) > 0);

            if any(anti_windup)
                g.reec.dreec5(anti_windup,k) = 0;
            end
        end

        mask = (y_piv < g.reec.iqmin);
        if any(mask)
            y_piv(mask) = g.reec.iqmin(mask);

            anti_windup = mask & (g.reec.reec_con(:,20) == 1) ...
                          & (g.reec.dreec5(:,k) < 0);

            if any(anti_windup)
                g.reec.dreec5(anti_windup,k) = 0;
            end
        end

        g.reec.dreec5(voltage_dip,k) = 0;

        %-----------------------------------------------------------------%
        % lower PI bypass loop

        g.reec.dreec6(:,k) = (qlim./max(g.reec.reec1(:,k),0.01) ...
                              - g.reec.reec6(:,k))./max(g.reec.reec_con(:,25),lbnd);

        mask = (g.reec.reec_con(:,25) < lbnd);
        if any(mask)                                  % integrator bypass
            g.reec.dreec6(mask,k) = 0.0;
            g.reec.reec6(mask,k) = qlim(mask)./max(g.reec.reec1(mask,k),0.01);
        end

        g.reec.dreec6(voltage_dip,k) = 0;

        %-----------------------------------------------------------------%
        % forming iqcmd

        iq2 = g.reec.reec6(:,k);

        qflag = (g.reec.reec_con(:,34) == 1);
        if any(qflag)
            iq2(qflag) = y_piv(qflag);
        end

        iqcmd = iqinj + iq2;

        g.reec.dreec10(:,k) = (iqcmd - g.reec.reec10(:,k)) ...
                              ./max(g.reec.reec_con(:,44),lbnd);

        mask = (g.reec.reec_con(:,44) < lbnd);
        if any(mask)                                  % integrator bypass
            g.reec.dreec10(mask,k) = 0.0;
            g.reec.reec10(mask,k) = iqcmd(mask);
        end

        %-----------------------------------------------------------------%
        % rate limit for real power reference

        if (nargin > 4)
            h_sol = varargin{1};
        else
             error('\nreec: incorrect number of arguments passed to function.');
        end

        mask = ((g.reec.pref(:,k) - g.reec.pref(:,max(k-1,1)))/h_sol ...
                > g.reec.reec_con(:,26));

        if any(mask)
            g.reec.pref(mask,k) = g.reec.pref(mask,max(k-1,1)) ...
                                  + h_sol*g.reec.reec_con(mask,26);
        end

        mask = ((g.reec.pref(:,k) - g.reec.pref(:,max(k-1,1)))/h_sol ...
                < g.reec.reec_con(:,27));

        if any(mask)
            g.reec.pref(mask,k) = g.reec.pref(mask,max(k-1,1)) ...
                                  + h_sol*g.reec.reec_con(mask,27);
        end

        %-----------------------------------------------------------------%
        % active power control loop

        g.reec.dreec7(:,k) = (g.reec.pref(:,k) - g.reec.reec7(:,k)) ...
                             ./max(g.reec.reec_con(:,30),lbnd);

        mask = (g.reec.reec_con(:,30) < lbnd);
        if any(mask)                                  % integrator bypass
            g.reec.dreec7(mask,k) = 0.0;
            g.reec.reec7(mask,k) = g.reec.pref(mask,k);
        end

        g.reec.dreec7(voltage_dip,k) = 0;

        g.reec.reec7(:,k) = min(g.reec.reec7(:,k),g.reec.reec_con(:,28));
        g.reec.reec7(:,k) = max(g.reec.reec7(:,k),g.reec.reec_con(:,29));

        ipcmd = (g.reec.reec7(:,k) + g.reec.paux(:,k)) ...
                ./max(g.reec.reec1(:,k),0.01);

        g.reec.dreec9(:,k) = (ipcmd - g.reec.reec9(:,k)) ...
                             ./max(g.reec.reec_con(:,43),lbnd);

        mask = (g.reec.reec_con(:,43) < lbnd);
        if any(mask)                                  % integrator bypass
            g.reec.dreec9(mask,k) = 0.0;
            g.reec.reec9(mask,k) = ipcmd(mask);
        end

        %-----------------------------------------------------------------%
        % passing the current command to the generator/converter (system base)

        g.reec.icmd(:,k) = ipcmd + 1j*iqcmd;

        g.ess.ess_scmd(g.reec.ess_idx,k) = g.reec.icmd(:,k)./g.reec.reec_pot(:,1);
    end
end

end  % function end

% eof
