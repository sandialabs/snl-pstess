function mac_ivm(i,k,bus,flag)
% Syntax: mac_ivm(i,k,bus,flag)
%
% Purpose: Internal Voltage Model (IVM) type generator
%          for representing grid-forming inverters
%
% Input: i - generator number
%          - 0 for vectorized computation
%        k - integer time
%        bus - solved loadflow bus data
%        flag - 0 - initialization
%               1 - network interface computation
%               2 - generator dynamics computation and state
%                   state matrix building
%
% ivm_con format
% col   data                                      units
%  1    ivm generator number (*not* mac number)   integer
%  2    bus number                                integer
%  3    ivm controller type                       integer
%           0 = user-defined; 1 = gfma            --
%  4    ivm generator mva base                    MVA
%  5    inverter coupling resistance (typ. zero)  pu
%  6    inverter coupling reactance               pu
%  7    ivm comm. angle time constant, Td         sec
%  8    ivm comm. voltage mag. time constant, Tv  sec
%  9    virtual inertia time constant             sec
%
% mac_con format (imported from ivm_con)
% col   data                                      units
%  1    generator number (*not* ivm number)       integer
%  2    bus number                                integer
%  3    generator mva base                        MVA
%  5    inverter coupling resistance (typ. zero)  pu
%  7    inverter coupling reactance               pu
%  9    ivm comm. angle time constant, Td         sec
% 10    ivm comm. voltage mag. time constant, Tv  sec
% 16    virtual inertia time constant             sec
% 19    bus number                                integer
%
% mac_pot matrix format
% col   data                                      units
%  1    factor to convert from syst. to gen base  pu
%
% mac_ivm states
% var       description
% eqprime   internal voltage phasor magnitude
% mac_ang   internal voltage phasor angle

%-----------------------------------------------------------------------------%
% Version history
%
% Version: 1.1
% Date:    October 2022
% Author:  Ryan Elliott
%
% Version: 1.0
% Date:    June 2019
% Author:  D. Trudnowski
%-----------------------------------------------------------------------------%

global g;  % declaring struct of global variables

lbnd = 1e-3;  % lower bound to prevent division by zero

if (g.mac.n_ivm ~= 0 && i ~= 0)
    % non-vectorized computation
    error('mac_ivm: all ivm calculations must be vectorized.');
elseif (g.mac.n_ivm ~= 0 && i == 0)
    if (flag == 0)   % initialization
        % busnum -- vector of bus numbers where ivm's are connected
        busnum = g.bus.bus_int(g.mac.mac_con(g.mac.mac_ivm_idx,2));

        % mac_pot(,1) -- factor to convert from system to gen MVA base
        % mac_pot(,2) -- base kV
        g.mac.mac_pot(g.mac.mac_ivm_idx,1) = g.sys.basmva ...
                                             ./g.mac.mac_con(g.mac.mac_ivm_idx,3);
        g.mac.mac_pot(g.mac.mac_ivm_idx,2) = 1.0*ones(g.mac.n_ivm,1);

        % eterm -- terminal bus voltage magnitude
        % theta -- terminal bus angle in radians
        g.mac.eterm(g.mac.mac_ivm_idx,1) = bus(busnum,2);
        g.bus.theta(busnum,1) = bus(busnum,3)*pi/180;

        % pelect -- electrical output power, active
        % qelect -- electrical output power, reactive
        g.mac.pelect(g.mac.mac_ivm_idx,1) = ...
            bus(busnum,4).*g.mac.mac_con(g.mac.mac_ivm_idx,22);
        g.mac.qelect(g.mac.mac_ivm_idx,1) = ...
            bus(busnum,5).*g.mac.mac_con(g.mac.mac_ivm_idx,23);

        % terminal voltage phasor in system reference frame
        Vt = g.mac.eterm(g.mac.mac_ivm_idx,1).*exp(1j*g.bus.theta(busnum,1));

        % complex current in system reference frame (generator base)
        It = conj((g.mac.pelect(g.mac.mac_ivm_idx,1) ...
                   + 1j*g.mac.qelect(g.mac.mac_ivm_idx,1))./Vt) ...
             .*g.mac.mac_pot(g.mac.mac_ivm_idx,1);

        % voltage behind reactance in system reference frame
        E = Vt + (g.mac.mac_con(g.mac.mac_ivm_idx,5) ...
                  + 1j*g.mac.mac_con(g.mac.mac_ivm_idx,7)).*It;

        % psi_re -- real-part of internal voltage used for network interface
        % psi_im -- imag-part of internal voltage used for network interface
        g.mac.psi_re(g.mac.mac_ivm_idx,1) = real(E);
        g.mac.psi_im(g.mac.mac_ivm_idx,1) = imag(E);

        % mac_ang -- internal voltage angle (rad)
        % mac_spd -- angular velocity of the voltage phasor in steady state (pu)
        g.mac.mac_ang(g.mac.mac_ivm_idx,1) = atan2(imag(E),real(E));
        g.mac.mac_spd(g.mac.mac_ivm_idx,1) = ones(g.mac.n_ivm,1);

        % system reference frame rotation to Park's frame
        rot = 1j*exp(-1j*g.mac.mac_ang(g.mac.mac_ivm_idx,1));

        % eprime = E.*rot;
        g.mac.eqprime(g.mac.mac_ivm_idx,1) = abs(E);  % imag(eprime) = abs(E); real(eprime) = 0;

        % vex -- commanded voltage phasor magnitude (pu)
        % fldcur -- commanded voltage phasor angle (rad)
        g.mac.vex(g.mac.mac_ivm_idx,1) = g.mac.eqprime(g.mac.mac_ivm_idx,1);
        g.mac.fldcur(g.mac.mac_ivm_idx,1) = g.mac.mac_ang(g.mac.mac_ivm_idx,1);

        curr = It.*rot;
        g.mac.curdg(g.mac.mac_ivm_idx,1) = real(curr);
        g.mac.curqg(g.mac.mac_ivm_idx,1) = imag(curr);

        g.mac.curd(g.mac.mac_ivm_idx,1) = ...
            real(curr)./g.mac.mac_pot(g.mac.mac_ivm_idx,1);
        g.mac.curq(g.mac.mac_ivm_idx,1) = ...
            imag(curr)./g.mac.mac_pot(g.mac.mac_ivm_idx,1);

        v = Vt.*rot;
        g.mac.ed(g.mac.mac_ivm_idx,1) = real(v);
        g.mac.eq(g.mac.mac_ivm_idx,1) = imag(v);

        % pmech -- mech power = elec power + losses (on generator base)
        % g.mac.pmech(g.mac.mac_ivm_idx,1) = ...
        %     g.mac.pelect(g.mac.mac_ivm_idx,1) ...
        %     .*g.mac.mac_pot(g.mac.mac_ivm_idx,1) ...
        %     + sign(g.mac.pelect(g.mac.mac_ivm_idx,1)) ...
        %       .*(abs(curr).^2).*g.mac.mac_con(g.mac.mac_ivm_idx,5);
    end

    if (flag == 1)  % network interface computation
        % g.mac.fldcur -- commanded voltage angle (rad)
        d_ord = g.mac.fldcur(g.mac.mac_ivm_idx,k);

        mask = (g.mac.mac_con(g.mac.mac_ivm_idx,9) < lbnd);
        if any(mask)
            g.mac.mac_ang(g.mac.mac_ivm_idx(mask),k) = d_ord(mask);
        end

        % g.mac.vex -- commanded voltage magnitude (pu)
        e_ord = g.mac.vex(g.mac.mac_ivm_idx,k);

        mask = (g.mac.mac_con(g.mac.mac_ivm_idx,10) < lbnd);
        if any(mask)
            g.mac.eqprime(g.mac.mac_ivm_idx(mask),k) = e_ord(mask);
        end

        % mac_ang -- angle with respect to machine reference (rad)
        g.mac.mac_ang(g.mac.mac_ivm_idx,k) = ...
            g.mac.mac_ang(g.mac.mac_ivm_idx,k) ...
            - g.sys.mach_ref(k)*ones(g.mac.n_ivm,1);

        % psi_re -- real part of internal voltage
        g.mac.psi_re(g.mac.mac_ivm_idx,k) = ...
            cos(g.mac.mac_ang(g.mac.mac_ivm_idx,k)) ...
            .*g.mac.eqprime(g.mac.mac_ivm_idx,k);

        % psi_im -- imag part of internal voltage
        g.mac.psi_im(g.mac.mac_ivm_idx,k) = ...
            sin(g.mac.mac_ang(g.mac.mac_ivm_idx,k)) ...
            .*g.mac.eqprime(g.mac.mac_ivm_idx,k);
    end

    if (flag == 2)  % generator dynamics calculation
        % g.mac.fldcur -- for ivm, this is the commanded voltage angle (rad)
        d_ord = g.mac.fldcur(g.mac.mac_ivm_idx,k);

        g.mac.dmac_ang(g.mac.mac_ivm_idx,k) = ...
            (d_ord - g.mac.mac_ang(g.mac.mac_ivm_idx,k)) ...
            ./max(g.mac.mac_con(g.mac.mac_ivm_idx,9),lbnd);

        mask = (g.mac.mac_con(g.mac.mac_ivm_idx,9) < lbnd);
        if any(mask)
            g.mac.dmac_ang(g.mac.mac_ivm_idx(mask),k) = 0.0;
            % mac_ang is updated in the network interface when bypassed
        end

        % g.mac.vex -- for ivm, this is the commanded voltage magnitude (pu)
        e_ord = g.mac.vex(g.mac.mac_ivm_idx,k);

        g.mac.deqprime(g.mac.mac_ivm_idx,k) = ...
            (e_ord - g.mac.eqprime(g.mac.mac_ivm_idx,k)) ...
            ./max(g.mac.mac_con(g.mac.mac_ivm_idx,10),lbnd);

        mask = (g.mac.mac_con(g.mac.mac_ivm_idx,10) < lbnd);
        if any(mask)
            g.mac.deqprime(g.mac.mac_ivm_idx(mask),k) = 0.0;
            % eqprime is updated in the network interface when bypassed
        end

        g.mac.curd(g.mac.mac_ivm_idx,k) = ...
            sin(g.mac.mac_ang(g.mac.mac_ivm_idx,k)) ...
            .*g.mac.cur_re(g.mac.mac_ivm_idx,k) ...
            - cos(g.mac.mac_ang(g.mac.mac_ivm_idx,k)) ...
              .*g.mac.cur_im(g.mac.mac_ivm_idx,k);

        g.mac.curq(g.mac.mac_ivm_idx,k) = ...
            cos(g.mac.mac_ang(g.mac.mac_ivm_idx,k)) ...
            .*g.mac.cur_re(g.mac.mac_ivm_idx,k) ...
            + sin(g.mac.mac_ang(g.mac.mac_ivm_idx,k)) ...
              .*g.mac.cur_im(g.mac.mac_ivm_idx,k);

        g.mac.curdg(g.mac.mac_ivm_idx,k) = ...
            g.mac.curd(g.mac.mac_ivm_idx,k).*g.mac.mac_pot(g.mac.mac_ivm_idx,1);

        g.mac.curqg(g.mac.mac_ivm_idx,k) = ...
            g.mac.curq(g.mac.mac_ivm_idx,k).*g.mac.mac_pot(g.mac.mac_ivm_idx,1);

        g.mac.ed(g.mac.mac_ivm_idx,k) = ...
            -g.mac.mac_con(g.mac.mac_ivm_idx,5) ...
             .*g.mac.curdg(g.mac.mac_ivm_idx,k) ...
            + g.mac.mac_con(g.mac.mac_ivm_idx,7) ...
              .*g.mac.curqg(g.mac.mac_ivm_idx,k);

        g.mac.eq(g.mac.mac_ivm_idx,k) = ...
            g.mac.eqprime(g.mac.mac_ivm_idx,k) ...
            - g.mac.mac_con(g.mac.mac_ivm_idx,5) ...
              .*g.mac.curqg(g.mac.mac_ivm_idx,k) ...
            - g.mac.mac_con(g.mac.mac_ivm_idx,7) ...
              .*g.mac.curdg(g.mac.mac_ivm_idx,k);

        g.mac.eterm(g.mac.mac_ivm_idx,k) = ...
            sqrt(g.mac.ed(g.mac.mac_ivm_idx,k).^2 ...
                 + g.mac.eq(g.mac.mac_ivm_idx,k).^2);

        g.mac.pelect(g.mac.mac_ivm_idx,k) = ...
            g.mac.eq(g.mac.mac_ivm_idx,k).*g.mac.curq(g.mac.mac_ivm_idx,k) ...
            + g.mac.ed(g.mac.mac_ivm_idx,k).*g.mac.curd(g.mac.mac_ivm_idx,k);

        g.mac.qelect(g.mac.mac_ivm_idx,k) = ...
            g.mac.eq(g.mac.mac_ivm_idx,k).*g.mac.curd(g.mac.mac_ivm_idx,k) ...
            - g.mac.ed(g.mac.mac_ivm_idx,k).*g.mac.curq(g.mac.mac_ivm_idx,k);
    end
end

end  % function end

% eof
