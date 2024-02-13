function mac_ib(i,k,bus,flag)
% Syntax: mac_ib(i,k,bus,flag)
%
% Purpose: generator infinite bus model
%     - forms a constant internal voltage model of a synchronous generator
%     - takes impedance from mac_con
%
% Note: infinite buses are defined in the vector ibus_con, which has zero
%     entries for generators that are not to be inf. buses and unity
%     entries for generators which are to be taken as infinite buses
%     uses the appropriate initialization model to determine the initial
%     internal voltage

%-----------------------------------------------------------------------------%

global g;  % declaring struct of global variables

% exit if no infinite buses
if (g.mac.n_ib ~= 0)
    if (flag == 0)
        % initialization
        if (i ~= 0)  % non-vector initialization
            % check that machine i is defined as an infinite bus
            ib_chk = find(g.mac.mac_ib_idx == i);
            if ~isempty(ib_chk)
                g.mac.mac_pot(i,1) = g.sys.basmva/g.mac.mac_con(i,3);
                % scaled MVA base
                % check for em generator
                if ~isempty(g.mac.mac_ib_em)
                    ib_em_chk = find(g.mac.mac_ib_em == i);
                    if ~isempty(ib_em_chk)
                        % initialize as classical machine model

                        % busnum -- bus number
                        busnum = g.bus.bus_int(g.mac.mac_con(i,2));

                        % eterm -- terminal bus voltage
                        % theta -- bus angle in rad
                        g.mac.eterm(i,1) = bus(busnum,2);
                        g.bus.theta(busnum,1) = bus(busnum,3)*pi/180;

                        % pelect -- electrical output power, active
                        % qelect -- electrical output power, reactive
                        g.mac.pelect(i,1) = bus(busnum,4)*g.mac.mac_con(i,22);
                        g.mac.qelect(i,1) = bus(busnum,5)*g.mac.mac_con(i,23);

                        % curr -- current magnitude on generator base
                        % phi -- power factor angle (rad)
                        curr = sqrt(g.mac.pelect(i,1)^2 + g.mac.qelect(i,1)^2) ...
                               /g.mac.eterm(i,1)*g.mac.mac_pot(i,1);

                        phi = atan2(g.mac.qelect(i,1),g.mac.pelect(i,1));

                        % v -- complex voltage in system reference frame
                        % curr -- complex current in system reference frame
                        v = g.mac.eterm(i,1)*exp(1j*g.bus.theta(busnum,1));
                        curr = curr*exp(1j*(g.bus.theta(busnum,1) - phi));

                        % voltage behind transient reactance
                        eprime = v + 1j*g.mac.mac_con(i,7)*curr;
                        ei = eprime;

                        % mac_ang -- machine angle (delta)
                        g.mac.mac_ang(i,1) = atan2(imag(ei),real(ei));

                        % system reference frame rotation
                        rot = 1j*exp(-1j*g.mac.mac_ang(i,1));
                        eprime = eprime*rot;

                        g.mac.edprime(i,1) = real(eprime);
                        g.mac.eqprime(i,1) = imag(eprime);

                        g.mac.psi_re(i,1) = ...
                            sin(g.mac.mac_ang(i,1))*g.mac.edprime(i,1) ...
                            + cos(g.mac.mac_ang(i,1))*g.mac.eqprime(i,1);

                        g.mac.psi_im(i,1) = ...
                            -cos(g.mac.mac_ang(i,1))*g.mac.edprime(i,1) ...
                            + sin(g.mac.mac_ang(i,1))*g.mac.eqprime(i,1);
                    end
                end

                % check for transient generator
                if ~isempty(g.mac.mac_ib_tra)
                    ib_tra_chk = find(g.mac.mac_ib_tra == i);
                    if ~isempty(ib_tra_chk) ~= 0
                        % initialize as transient generator

                        % busnum -- bus number
                        busnum = g.bus.bus_int(g.mac.mac_con(i,2));

                        % eterm -- terminal bus voltage
                        % theta -- bus angle in rad
                        g.mac.eterm(i,1) = bus(busnum,2);
                        g.bus.theta(busnum,1) = bus(busnum,3)*pi/180;

                        % pelect -- electrical output power, active
                        % qelect -- electrical output power, reactive
                        g.mac.pelect(i,1) = bus(busnum,4)*g.mac.mac_con(i,22);
                        g.mac.qelect(i,1) = bus(busnum,5)*g.mac.mac_con(i,23);

                        % curr -- current magnitude on generator base
                        % phi -- power factor angle (rad)
                        curr = sqrt(g.mac.pelect(i,1)^2 + g.mac.qelect(i,1)^2) ...
                               /g.mac.eterm(i,1)*g.mac.mac_pot(i,1);

                        phi = atan2(g.mac.qelect(i,1),g.mac.pelect(i,1));

                        % complex voltage in system reference frame
                        v = g.mac.eterm(i,1)*exp(1j*g.bus.theta(busnum,1));

                        % complex current in system reference frame
                        curr = curr*exp(1j*(g.bus.theta(busnum,1) - phi));

                        % voltage behind transient reactance
                        eprime = v + (g.mac.mac_con(i,5) ...
                                      + 1j*g.mac.mac_con(i,7))*curr;

                        ei = v + (g.mac.mac_con(i,5) ...
                                  + 1j*g.mac.mac_con(i,11))*curr;

                        % mac_ang -- machine angle (delta)
                        g.mac.mac_ang(i,1) = atan2(imag(ei),real(ei));

                        % system reference frame rotation
                        rot = 1j*exp(-1j*g.mac.mac_ang(i,1));
                        eprime = eprime*rot;

                        g.mac.edprime(i,1) = real(eprime);
                        g.mac.eqprime(i,1) = imag(eprime);

                        g.mac.psi_re(i,1) = ...
                            sin(g.mac.mac_ang(i,1))*g.mac.edprime(i,1) ...
                            + cos(g.mac.mac_ang(i,1))*g.mac.eqprime(i,1);

                        g.mac.psi_im(i,1) = ...
                            -cos(g.mac.mac_ang(i,1))*g.mac.edprime(i,1) ...
                            + sin(g.mac.mac_ang(i,1))*g.mac.eqprime(i,1);
                    end
                end

                % check for subtransient generator model
                if ~isempty(g.mac.mac_ib_sub)
                    ib_sub_chk = find(g.mac.mac_ib_sub == i);
                    if ~isempty(ib_sub_chk)
                        % initialize as subtransient model

                        % busnum -- bus number
                        busnum = g.bus.bus_int(g.mac.mac_con(i,2));

                        % eterm -- terminal bus voltage
                        % theta -- bus angle in rad
                        g.mac.eterm(i,1) = bus(busnum,2);
                        g.bus.theta(busnum,1) = bus(busnum,3)*pi/180;

                        % pelect -- electrical output power, active
                        % qelect -- electrical output power, reactive
                        g.mac.pelect(i,1) = bus(busnum,4)*g.mac.mac_con(i,22);
                        g.mac.qelect(i,1) = bus(busnum,5)*g.mac.mac_con(i,23);

                        % curr -- current magnitude on generator base
                        % phi -- power factor angle (rad)
                        curr = sqrt(g.mac.pelect(i,1)^2 + g.mac.qelect(i,1)^2) ...
                               /g.mac.eterm(i,1)*g.mac.mac_pot(i,1);

                        phi = atan2(g.mac.qelect(i,1),g.mac.pelect(i,1));

                        % complex voltage in system reference frame
                        v = g.mac.eterm(i,1)*exp(1j*g.bus.theta(busnum,1));

                        % complex current in system reference frame
                        curr = curr*exp(1j*(g.bus.theta(busnum,1) - phi));

                        ei = v + (g.mac.mac_con(i,5) + 1j*g.mac.mac_con(i,11))*curr;

                        % mac_ang -- machine angle (delta)
                        g.mac.mac_ang(i,1) = atan2(imag(ei),real(ei));

                        % system reference frame rotation
                        rot = 1j*exp(-1j*g.mac.mac_ang(i,1));
                        curr = curr*rot;

                        g.mac.curdg(i,1) = real(curr);
                        g.mac.curqg(i,1) = imag(curr);

                        v = v*rot;
                        g.mac.ed(i,1) = real(v);
                        g.mac.eq(i,1) = imag(v);

                        eqra = g.mac.eq(i,1) + g.mac.mac_con(i,5)*g.mac.curqg(i,1);

                        g.mac.psidpp = eqra + g.mac.mac_con(i,8)*g.mac.curdg(i,1);

                        g.mac.psikd(i,1) = ...
                            eqra + g.mac.mac_con(i,4)*g.mac.curdg(i,1);

                        g.mac.eqprime(i,1) = ...
                            eqra + g.mac.mac_con(i,7)*g.mac.curdg(i,1);

                        edra = -g.mac.ed(i,1) - g.mac.mac_con(i,5)*g.mac.curdg(i,1);

                        g.mac.psiqpp = edra + g.mac.mac_con(i,13)*g.mac.curqg(i,1);

                        g.mac.psikq(i,1) = ...
                            edra + g.mac.mac_con(i,4)*g.mac.curqg(i,1);

                        g.mac.edprime(i,1) = ...
                            edra + g.mac.mac_con(i,12)*g.mac.curqg(i,1);

                        g.mac.psi_re(i,1) = ...
                            sin(g.mac.mac_ang(i,1)).*(-g.mac.psiqpp) ...
                            + cos(g.mac.mac_ang(i,1)).*g.mac.psidpp;

                        g.mac.psi_im(i,1) = ...
                            -cos(g.mac.mac_ang(i,1)).*(-g.mac.psiqpp) ...
                            + sin(g.mac.mac_ang(i,1)).*g.mac.psidpp;
                    end
                end

                % check for ivm generator model
                if ~isempty(g.mac.mac_ib_ivm)
                    ib_ivm_chk = find(g.mac.mac_ib_ivm == i);
                    if ~isempty(ib_ivm_chk)
                        error('mac_ib: all ivm calculations must be vectorized.');
                    end
                end
            end
        else
            % vector computation

            % mac_pot(,1) -- scaled MVA base
            g.mac.mac_pot(g.mac.mac_ib_idx,1) = ...
                g.sys.basmva*ones(g.mac.n_ib,1)./g.mac.mac_con(g.mac.mac_ib_idx,3);

            % check for em models
            if (g.mac.n_ib_em ~= 0)
                % busnum -- bus number vector
                busnum = g.bus.bus_int(g.mac.mac_con(g.mac.mac_ib_em,2));

                % eterm -- terminal bus voltage
                % theta -- bus angle in rad
                g.mac.eterm(g.mac.mac_ib_em,1) = bus(busnum,2);
                g.bus.theta(busnum,1) = bus(busnum,3)*pi/180;

                % pelect -- electrical output power, active
                % qelect -- electrical output power, reactive
                g.mac.pelect(g.mac.mac_ib_em,1) = ...
                    bus(busnum,4).*g.mac.mac_con(g.mac.mac_ib_em,22);

                g.mac.qelect(g.mac.mac_ib_em,1) = ...
                    bus(busnum,5).*g.mac.mac_con(g.mac.mac_ib_em,23);

                % curr -- current magnitude on generator base
                % phi -- power factor angle (rad)
                curr = sqrt(g.mac.pelect(g.mac.mac_ib_em,1).^2 ...
                            + g.mac.qelect(g.mac.mac_ib_em,1).^2) ...
                       ./g.mac.eterm(g.mac.mac_ib_em,1) ...
                       .*g.mac.mac_pot(g.mac.mac_ib_em,1);

                phi = atan2(g.mac.qelect(g.mac.mac_ib_em,1), ...
                            g.mac.pelect(g.mac.mac_ib_em,1));

                % complex voltage in system reference frame
                v = g.mac.eterm(g.mac.mac_ib_em,1).*exp(1j*g.bus.theta(busnum,1));

                % complex current in system reference frame
                curr = curr.*exp(1j*(g.bus.theta(busnum,1) - phi));

                % voltage behind transient reactance
                eprime = v + 1j*g.mac.mac_con(g.mac.mac_ib_em,7).*curr;
                ei = eprime;

                % mac_ang -- machine angle (delta)
                g.mac.mac_ang(g.mac.mac_ib_em,1) = atan2(imag(ei),real(ei));

                % system reference frame rotation
                rot = 1j*exp(-1j*g.mac.mac_ang(g.mac.mac_ib_em,1));
                eprime = eprime.*rot;

                g.mac.edprime(g.mac.mac_ib_em,1) = real(eprime);
                g.mac.eqprime(g.mac.mac_ib_em,1) = imag(eprime);

                g.mac.psi_re(g.mac.mac_ib_em,1) = ...
                    sin(g.mac.mac_ang(g.mac.mac_ib_em,1)) ...
                    .*g.mac.edprime(g.mac.mac_ib_em,1) ...
                    + cos(g.mac.mac_ang(g.mac.mac_ib_em,1)) ...
                      .*g.mac.eqprime(g.mac.mac_ib_em,1);

                g.mac.psi_im(g.mac.mac_ib_em,1) = ...
                    -cos(g.mac.mac_ang(g.mac.mac_ib_em,1)) ...
                     .*g.mac.edprime(g.mac.mac_ib_em,1) ...
                    + sin(g.mac.mac_ang(g.mac.mac_ib_em,1)) ...
                      .*g.mac.eqprime(g.mac.mac_ib_em,1);
            end

            % check for transient models
            if (g.mac.n_ib_tra ~= 0)
                % busnum -- bus number vector
                busnum = g.bus.bus_int(g.mac.mac_con(g.mac.mac_ib_tra,2));

                % eterm -- terminal bus voltage
                % theta -- bus angle in rad
                g.mac.eterm(g.mac.mac_ib_tra,1) = bus(busnum,2);
                g.bus.theta(busnum,1) = bus(busnum,3)*pi/180;

                % pelect -- electrical output power, active
                % qelect -- electrical output power, reactive
                g.mac.pelect(g.mac.mac_ib_tra,1) = ...
                    bus(busnum,4).*g.mac.mac_con(g.mac.mac_ib_tra,22);

                g.mac.qelect(g.mac.mac_ib_tra,1) = ...
                    bus(busnum,5).*g.mac.mac_con(g.mac.mac_ib_tra,23);

                % curr -- current magnitude on generator base
                % phi -- power factor angle (rad)
                curr = sqrt(g.mac.pelect(g.mac.mac_ib_tra,1).^2 ...
                            + g.mac.qelect(g.mac.mac_ib_tra,1).^2) ...
                       ./g.mac.eterm(g.mac.mac_ib_tra,1) ...
                       .*g.mac.mac_pot(g.mac.mac_ib_tra,1);

                phi = atan2(g.mac.qelect(g.mac.mac_ib_tra,1), ...
                            g.mac.pelect(g.mac.mac_ib_tra,1));

                % complex voltage in system reference frame
                v = g.mac.eterm(g.mac.mac_ib_tra,1).*exp(1j*g.bus.theta(busnum,1));

                % complex current in system reference frame
                curr = curr.*exp(1j*(g.bus.theta(busnum,1) - phi));

                % voltage behind transient reactance
                eprime = v + (g.mac.mac_con(g.mac.mac_ib_tra,5) ...
                              + 1j*g.mac.mac_con(g.mac.mac_ib_tra,7)).*curr;

                ei = v + (g.mac.mac_con(g.mac.mac_ib_tra,5) ...
                          + 1j*g.mac.mac_con(g.mac.mac_ib_tra,11)).*curr;

                % mac_ang -- machine angle (delta)
                g.mac.mac_ang(g.mac.mac_ib_tra,1) = atan2(imag(ei),real(ei));

                % system reference frame rotation
                rot = 1j*exp(-1j*g.mac.mac_ang(g.mac.mac_ib_tra,1));
                eprime = eprime.*rot;

                g.mac.edprime(g.mac.mac_ib_tra,1) = real(eprime);
                g.mac.eqprime(g.mac.mac_ib_tra,1) = imag(eprime);

                g.mac.psi_re(g.mac.mac_ib_tra,1) = ...
                    sin(g.mac.mac_ang(g.mac.mac_ib_tra,1)) ...
                    .*g.mac.edprime(g.mac.mac_ib_tra,1) ...
                    + cos(g.mac.mac_ang(g.mac.mac_ib_tra,1)) ...
                      .*g.mac.eqprime(g.mac.mac_ib_tra,1);

                g.mac.psi_im(g.mac.mac_ib_tra,1) = ...
                    -cos(g.mac.mac_ang(g.mac.mac_ib_tra,1)) ...
                     .*g.mac.edprime(g.mac.mac_ib_tra,1) ...
                    + sin(g.mac.mac_ang(g.mac.mac_ib_tra,1)) ...
                      .*g.mac.eqprime(g.mac.mac_ib_tra,1);
            end

            % check for subtransient models
            if (g.mac.n_ib_sub ~= 0)
                % busnum -- bus number vector
                busnum = g.bus.bus_int(g.mac.mac_con(g.mac.mac_ib_sub,2));

                % eterm -- terminal bus voltage
                % theta -- bus angle in rad
                g.mac.eterm(g.mac.mac_ib_sub,1) = bus(busnum,2);
                g.bus.theta(busnum,1) = bus(busnum,3)*pi/180;

                % pelect -- electrical output power, active
                % qelect -- electrical output power, reactive
                g.mac.pelect(g.mac.mac_ib_sub,1) = ...
                    bus(busnum,4).*g.mac.mac_con(g.mac.mac_ib_sub,22);

                g.mac.qelect(g.mac.mac_ib_sub,1) = ...
                    bus(busnum,5).*g.mac.mac_con(g.mac.mac_ib_sub,23);

                % curr -- current magnitude on generator base
                % phi -- power factor angle (rad)
                curr = sqrt(g.mac.pelect(g.mac.mac_ib_sub,1).^2 ...
                            + g.mac.qelect(g.mac.mac_ib_sub,1).^2) ...
                       ./g.mac.eterm(g.mac.mac_ib_sub,1) ...
                       .*g.mac.mac_pot(g.mac.mac_ib_sub,1);

                phi = atan2(g.mac.qelect(g.mac.mac_ib_sub,1), ...
                            g.mac.pelect(g.mac.mac_ib_sub,1));

                % complex voltage in system reference frame
                v = g.mac.eterm(g.mac.mac_ib_sub,1).*exp(1j*g.bus.theta(busnum,1));

                % complex current in system reference frame
                curr = curr.*exp(1j*(g.bus.theta(busnum,1) - phi));

                ei = v + (g.mac.mac_con(g.mac.mac_ib_sub,5) ...
                          + 1j*g.mac.mac_con(g.mac.mac_ib_sub,11)).*curr;

                % mac_ang -- machine angle (delta)
                g.mac.mac_ang(g.mac.mac_ib_sub,1) = atan2(imag(ei),real(ei));

                % system reference frame rotation
                rot = 1j*exp(-1j*g.mac.mac_ang(g.mac.mac_ib_sub,1));
                curr = curr.*rot;

                g.mac.curdg(g.mac.mac_ib_sub,1) = real(curr);
                g.mac.curqg(g.mac.mac_ib_sub,1) = imag(curr);

                g.mac.curd(g.mac.mac_ib_sub,1) = ...
                    real(curr)./g.mac.mac_pot(g.mac.mac_ib_sub,1);

                g.mac.curq(g.mac.mac_ib_sub,1) = ...
                    imag(curr)./g.mac.mac_pot(g.mac.mac_ib_sub,1);

                v = v.*rot;
                g.mac.ed(g.mac.mac_ib_sub,1) = real(v);
                g.mac.eq(g.mac.mac_ib_sub,1) = imag(v);

                eqra = g.mac.eq(g.mac.mac_ib_sub,1) ...
                       + g.mac.mac_con(g.mac.mac_ib_sub,5) ...
                         .*g.mac.curqg(g.mac.mac_ib_sub,1);

                g.mac.psidpp = eqra + g.mac.mac_con(g.mac.mac_ib_sub,8) ...
                                      .*g.mac.curdg(g.mac.mac_ib_sub,1);

                g.mac.psikd(g.mac.mac_ib_sub,1) = ...
                    eqra + g.mac.mac_con(g.mac.mac_ib_sub,4) ...
                           .*g.mac.curdg(g.mac.mac_ib_sub,1);

                g.mac.eqprime(g.mac.mac_ib_sub,1) = ...
                    eqra + g.mac.mac_con(g.mac.mac_ib_sub,7) ...
                           .*g.mac.curdg(g.mac.mac_ib_sub,1);

                edra = -g.mac.ed(g.mac.mac_ib_sub,1) ...
                       - g.mac.mac_con(g.mac.mac_ib_sub,5) ...
                         .*g.mac.curdg(g.mac.mac_ib_sub,1);

                g.mac.psiqpp = edra + g.mac.mac_con(g.mac.mac_ib_sub,13) ...
                                      .*g.mac.curqg(g.mac.mac_ib_sub,1);

                g.mac.psikq(g.mac.mac_ib_sub,1) = ...
                    edra + g.mac.mac_con(g.mac.mac_ib_sub,4) ...
                           .*g.mac.curqg(g.mac.mac_ib_sub,1);

                g.mac.edprime(g.mac.mac_ib_sub,1) = ...
                    edra + g.mac.mac_con(g.mac.mac_ib_sub,12) ...
                           .*g.mac.curqg(g.mac.mac_ib_sub,1);

                g.mac.psi_re(g.mac.mac_ib_sub,1) = ...
                    sin(g.mac.mac_ang(g.mac.mac_ib_sub,1)).*(-g.mac.psiqpp) ...
                    + cos(g.mac.mac_ang(g.mac.mac_ib_sub,1)).*g.mac.psidpp;

                g.mac.psi_im(g.mac.mac_ib_sub,1) = ...
                    -cos(g.mac.mac_ang(g.mac.mac_ib_sub,1)).*(-g.mac.psiqpp) ...
                    + sin(g.mac.mac_ang(g.mac.mac_ib_sub,1)).*g.mac.psidpp;
            end

            % check for ivm generator models
            if (g.mac.n_ib_ivm ~= 0)
                % busnum -- bus number vector
                busnum = g.bus.bus_int(g.mac.mac_con(g.mac.mac_ib_ivm,2));

                % eterm -- terminal bus voltage magnitude
                % theta -- terminal bus angle in radians
                g.mac.eterm(g.mac.mac_ib_ivm,1) = bus(busnum,2);
                g.bus.theta(busnum,1) = bus(busnum,3)*pi/180;

                % pelect -- electrical output power, active
                % qelect -- electrical output power, reactive
                g.mac.pelect(g.mac.mac_ib_ivm,1) = ...
                    bus(busnum,4).*g.mac.mac_con(g.mac.mac_ib_ivm,22);

                g.mac.qelect(g.mac.mac_ib_ivm,1) = ...
                    bus(busnum,5).*g.mac.mac_con(g.mac.mac_ib_ivm,23);

                % terminal voltage phasor in system reference frame
                Vt = g.mac.eterm(g.mac.mac_ib_ivm,1) ...
                     .*exp(1j*g.bus.theta(busnum,1));

                % complex current in system reference frame (generator base)
                It = conj((g.mac.pelect(g.mac.mac_ib_ivm,1) ...
                           + 1j*g.mac.qelect(g.mac.mac_ib_ivm,1))./Vt) ...
                     .*g.mac.mac_pot(g.mac.mac_ib_ivm,1);

                % voltage behind reactance in system reference frame
                E = Vt + (g.mac.mac_con(g.mac.mac_ib_ivm,5) ...
                          + 1j*g.mac.mac_con(g.mac.mac_ib_ivm,7)).*It;

                % psi_re -- real-part of internal voltage used for network interface
                % psi_im -- imag-part of internal voltage used for network interface
                g.mac.psi_re(g.mac.mac_ib_ivm,1) = real(E);
                g.mac.psi_im(g.mac.mac_ib_ivm,1) = imag(E);

                % mac_ang -- internal voltage angle (rad)
                % mac_spd -- angular velocity of the voltage phasor in steady state (pu)
                g.mac.mac_ang(g.mac.mac_ib_ivm,1) = atan2(imag(E),real(E));
                g.mac.mac_spd(g.mac.mac_ib_ivm,1) = ones(g.mac.n_ivm,1);

                % system reference frame rotation to Park's frame
                rot = 1j*exp(-1j*g.mac.mac_ang(g.mac.mac_ib_ivm,1));

                % eprime = E.*rot;
                g.mac.eqprime(g.mac.mac_ib_ivm,1) = abs(E);  % imag(eprime) = abs(E); real(eprime) = 0;

                % vex -- commanded voltage phasor magnitude (pu)
                % fldcur -- commanded voltage phasor angle (rad)
                g.mac.vex(g.mac.mac_ib_ivm,:) = g.mac.eqprime(g.mac.mac_ib_ivm,1);
                g.mac.fldcur(g.mac.mac_ib_ivm,:) = g.mac.mac_ang(g.mac.mac_ib_ivm,1);

                curr = It.*rot;
                g.mac.curdg(g.mac.mac_ib_ivm,1) = real(curr);
                g.mac.curqg(g.mac.mac_ib_ivm,1) = imag(curr);

                g.mac.curd(g.mac.mac_ib_ivm,1) = ...
                    real(curr)./g.mac.mac_pot(g.mac.mac_ib_ivm,1);

                g.mac.curq(g.mac.mac_ib_ivm,1) = ...
                    imag(curr)./g.mac.mac_pot(g.mac.mac_ib_ivm,1);

                v = Vt.*rot;
                g.mac.ed(g.mac.mac_ib_ivm,1) = real(v);
                g.mac.eq(g.mac.mac_ib_ivm,1) = imag(v);
            end
        end

        % end initialization
    end

    if (flag == 1)
        % network interface

        % for all infinite buses set psi to the initial value
        if (i ~= 0)
            % test for infinite bus
            ib_chk = find(g.mac.mac_ib_idx == i);
            if ~isempty(ib_chk)
                g.mac.psi_re(i,k) = g.mac.psi_re(i,1);
                g.mac.psi_im(i,k) = g.mac.psi_im(i,1);
            end
        else
            g.mac.psi_re(g.mac.mac_ib_idx,k) = g.mac.psi_re(g.mac.mac_ib_idx,1);
            g.mac.psi_im(g.mac.mac_ib_idx,k) = g.mac.psi_im(g.mac.mac_ib_idx,1);
        end
    end

    if (flag == 2)
        % no dynamics
    end
end

end  % function end

% eof
