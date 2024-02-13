function mac_em(i,k,bus,flag)
% Syntax: mac_em(i,k,bus,flag)
%
% Purpose: Electromechanical generator model, with
%          vectorized computation option
%
% Note:    State variables are: mac_ang, mac_spd
%
% Input:   i - generator number
%            - 0 for vectorized computation
%          k - integer time
%          bus - solved loadflow bus data
%          flag - 0 - initialization
%                 1 - network interface computation
%                 2 - generator dynamics computation and state
%                     state matrix building

%-----------------------------------------------------------------------------%
% Version history
%
% Version: 2.0
% Date:    June 1996
% Author:  Graham Rogers
% Purpose: add facility to allow different machine models in vector run
%
% Version: 1.0
% Author:  Joe H. Chow
% Date:    January 1991
%-----------------------------------------------------------------------------%

global g;  % declaring struct of global variables

if (g.mac.n_em ~= 0)
    if (flag == 0)  % initialization
        if (i ~= 0)
            % non-vector calculation
            % check for em model
            em = find(g.mac.mac_em_idx == i)
            if (length(em) ~= 0)
                % busnum -- bus number
                busnum = g.bus.bus_int(g.mac.mac_con(i,2));

                % mac_pot(i,1) -- scaled MVA base
                % mac_pot(i,2) -- base kV
                g.mac.mac_pot(i,1) = g.sys.basmva/g.mac.mac_con(i,3);
                g.mac.mac_pot(i,2) = 1.0;

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
                curr = sqrt(g.mac.pelect(i,1)^2+g.mac.qelect(i,1)^2) ...
                       /g.mac.eterm(i,1)*g.mac.mac_pot(i,1);

                phi = atan2(g.mac.qelect(i,1),g.mac.pelect(i,1));

                % complex voltage in system reference frame
                v = g.mac.eterm(i,1)*exp(1j*g.bus.theta(busnum,1));

                % complex current in system reference frame
                curr = curr*exp(1j*(g.bus.theta(busnum,1) - phi));

                % voltage behind transient reactance
                eprime = v + 1j*g.mac.mac_con(i,7)*curr;
                ei = eprime;

                % mac_ang -- machine angle (delta)
                % mac_spd -- machine speed at steady state
                g.mac.mac_ang(i,1) = atan2(imag(ei),real(ei));
                g.mac.mac_spd(i,1) = 1.0;

                % system reference frame rotation
                rot = 1j*exp(-1j*g.mac.mac_ang(i,1));
                g.mac.psi_re(i,1) = real(eprime);
                g.mac.psi_im(i,1) = imag(eprime);

                eprime = eprime*rot;
                g.mac.edprime(i,1) = real(eprime);
                g.mac.eqprime(i,1) = imag(eprime);

                curr = curr*rot;  % current on Park's frame
                g.mac.curdg(i,1) = real(curr);
                g.mac.curqg(i,1) = imag(curr);
                g.mac.curd(i,1) = real(curr)/g.mac.mac_pot(i,1);
                g.mac.curq(i,1) = imag(curr)/g.mac.mac_pot(i,1);

                % convert to system base
                v = v*rot;
                g.mac.ed(i,1) = real(v);
                g.mac.eq(i,1) = imag(v);  % in Park's frame
                g.mac.vex(i,1) = g.mac.eqprime(i,1);

                % set input mechanical power equal to electrical output power
                % since losses are zero for em model. On generator base
                g.mac.pmech(i,1) = g.mac.pelect(i,1)*g.mac.mac_pot(i,1);
            end
        else
            % vectorized computation
            busnum = g.bus.bus_int(g.mac.mac_con(g.mac.mac_em_idx,2));  % bus vector

            % mac_pot(,1) -- scaled MVA base
            % mac_pot(,2) -- base kV
            g.mac.mac_pot(g.mac.mac_em_idx,1) = g.sys.basmva*ones(g.mac.n_em,1) ...
                                                ./g.mac.mac_con(g.mac.mac_em_idx,3);

            g.mac.mac_pot(g.mac.mac_em_idx,2) = 1.0*ones(g.mac.n_em,1);

            % eterm -- terminal bus voltage
            % theta -- bus angle in rad
            g.mac.eterm(g.mac.mac_em_idx,1) = bus(busnum,2);
            g.bus.theta(busnum,1) = bus(busnum,3)*pi/180;

            % pelect -- electrical output power, active
            % qelect -- electrical output power, reactive
            g.mac.pelect(g.mac.mac_em_idx,1) = bus(busnum,4) ...
                                               .*g.mac.mac_con(g.mac.mac_em_idx,22);

            g.mac.qelect(g.mac.mac_em_idx,1) = bus(busnum,5) ...
                                               .*g.mac.mac_con(g.mac.mac_em_idx,23);

            % current magnitude on generator base
            curr = sqrt(g.mac.pelect(g.mac.mac_em_idx,1).^2 ...
                        + g.mac.qelect(g.mac.mac_em_idx,1).^2) ...
                   ./g.mac.eterm(g.mac.mac_em_idx,1) ...
                   .*g.mac.mac_pot(g.mac.mac_em_idx,1);

            % power factor angle
            phi = atan2(g.mac.qelect(g.mac.mac_em_idx,1), ...
                        g.mac.pelect(g.mac.mac_em_idx,1));

            % complex voltage in system reference frame
            v = g.mac.eterm(g.mac.mac_em_idx,1).*exp(1j*g.bus.theta(busnum,1));

            % complex current in system reference frame
            curr = curr.*exp(1j*(g.bus.theta(busnum,1) - phi));

            eprime = v + 1j*g.mac.mac_con(g.mac.mac_em_idx,7).*curr;
            ei = eprime;

            % mac_ang -- machine angle (delta)
            % mac_spd -- machine speed at steady state
            g.mac.mac_ang(g.mac.mac_em_idx,1) = atan2(imag(ei),real(ei));
            g.mac.mac_spd(g.mac.mac_em_idx,1) = ones(g.mac.n_em,1);

            % system reference frame rotation
            rot = 1j*exp(-1j*g.mac.mac_ang(g.mac.mac_em_idx,1));

            g.mac.psi_re(g.mac.mac_em_idx,1) = real(eprime);
            g.mac.psi_im(g.mac.mac_em_idx,1) = imag(eprime);

            eprime = eprime.*rot;  % in Park's frame
            g.mac.edprime(g.mac.mac_em_idx,1) = real(eprime);
            g.mac.eqprime(g.mac.mac_em_idx,1) = imag(eprime);

            curr = curr.*rot;      % in Park's frame
            g.mac.curdg(g.mac.mac_em_idx,1) = real(curr);
            g.mac.curqg(g.mac.mac_em_idx,1) = imag(curr);

            g.mac.curd(g.mac.mac_em_idx,1) = real(curr) ...
                                             ./g.mac.mac_pot(g.mac.mac_em_idx,1);
            g.mac.curq(g.mac.mac_em_idx,1) = imag(curr) ...
                                             ./g.mac.mac_pot(g.mac.mac_em_idx,1);

            v = v.*rot;            % in Park's frame
            g.mac.ed(g.mac.mac_em_idx,1) = real(v);
            g.mac.eq(g.mac.mac_em_idx,1) = imag(v);
            g.mac.vex(g.mac.mac_em_idx,1) = g.mac.eqprime(g.mac.mac_em_idx,1);

            % set input mechanical power equal to electrical output power
            % since losses are zero in em model (on generator base)
            g.mac.pmech(g.mac.mac_em_idx,1) = g.mac.pelect(g.mac.mac_em_idx,1) ...
                                              .*g.mac.mac_pot(g.mac.mac_em_idx,1);
        end

        % end initialization
    end

    if (flag == 1)  % network interface computation
        if (i ~= 0)
            % check for em machine
            em = find(g.mac.mac_em_idx == i);
            if (length(em) ~= 0)
                % angle wrt machine reference
                g.mac.mac_ang(i,k) = g.mac.mac_ang(i,k) - g.sys.mach_ref(k);

                % real part of psi
                g.mac.psi_re(i,k) = sin(g.mac.mac_ang(i,k))*g.mac.edprime(i,k) ...
                                    + cos(g.mac.mac_ang(i,k))*g.mac.eqprime(i,k);

                % imag part of psi
                g.mac.psi_im(i,k) = -cos(g.mac.mac_ang(i,k))*g.mac.edprime(i,k) ...
                                    + sin(g.mac.mac_ang(i,k))*g.mac.eqprime(i,k);
            end
        else  % vectorized computation
            % angle wrt machine reference
            g.mac.mac_ang(g.mac.mac_em_idx,k) = ...
                g.mac.mac_ang(g.mac.mac_em_idx,k) ...
                - g.sys.mach_ref(k)*ones(g.mac.n_em,1);

            % real part of psi
            g.mac.psi_re(g.mac.mac_em_idx,k) = ...
                sin(g.mac.mac_ang(g.mac.mac_em_idx,k)) ...
                .*g.mac.edprime(g.mac.mac_em_idx,k) ...
                + cos(g.mac.mac_ang(g.mac.mac_em_idx,k)) ...
                  .*g.mac.eqprime(g.mac.mac_em_idx,k);

            % imag part of psi
            g.mac.psi_im(g.mac.mac_em_idx,k) = ...
                -cos(g.mac.mac_ang(g.mac.mac_em_idx,k)) ...
                 .*g.mac.edprime(g.mac.mac_em_idx,k) ...
                + sin(g.mac.mac_ang(g.mac.mac_em_idx,k)) ...
                  .*g.mac.eqprime(g.mac.mac_em_idx,k);
        end

        % end interface
    end

    if (flag == 2)  % generator dynamics calculation
        if (i ~= 0)
            % check for em machine
            em = find(g.mac.mac_em_idx == i);
            if (length(em) ~= 0)
                g.mac.curd(i,k) = sin(g.mac.mac_ang(i,k))*g.mac.cur_re(i,k) ...
                                  - cos(g.mac.mac_ang(i,k))*g.mac.cur_im(i,k);
                g.mac.curq(i,k) = cos(g.mac.mac_ang(i,k))*g.mac.cur_re(i,k) ...
                                  + sin(g.mac.mac_ang(i,k))*g.mac.cur_im(i,k);

                g.mac.curdg(i,k) = g.mac.curd(i,k)*g.mac.mac_pot(i,1);
                g.mac.curqg(i,k) = g.mac.curq(i,k)*g.mac.mac_pot(i,1);

                g.mac.dedprime(i,k) = 0;
                g.mac.deqprime(i,k) = 0;

                g.mac.ed(i,k) = g.mac.edprime(i,k) ...
                                + g.mac.mac_con(i,7)*g.mac.curqg(i,k);
                g.mac.eq(i,k) = g.mac.eqprime(i,k) ...
                                - g.mac.mac_con(i,7)*g.mac.curdg(i,k);

                g.mac.eterm(i,k) = sqrt(g.mac.ed(i,k)^2+g.mac.eq(i,k)^2);

                g.mac.pelect(i,k) = g.mac.eq(i,k)*g.mac.curq(i,k) ...
                                    + g.mac.ed(i,k)*g.mac.curd(i,k);
                g.mac.qelect(i,k) = g.mac.eq(i,k)*g.mac.curd(i,k) ...
                                    - g.mac.ed(i,k)*g.mac.curq(i,k);

                g.mac.dmac_ang(i,k) = g.sys.basrad*(g.mac.mac_spd(i,k) - 1);

                g.mac.dmac_spd(i,k) = ...
                    (g.mac.pmech(i,k) + g.mac.pm_sig(i,k) ...
                     - g.mac.pelect(i,k)*g.mac.mac_pot(i,1) ...
                     - g.mac.mac_con(i,17)*(g.mac.mac_spd(i,k) - 1)) ...
                    /(2*g.mac.mac_con(i,16));
            end
        else
            % vectorized computation
            g.mac.curd(g.mac.mac_em_idx,k) = ...
                sin(g.mac.mac_ang(g.mac.mac_em_idx,k)) ...
                .*g.mac.cur_re(g.mac.mac_em_idx,k) ...
                - cos(g.mac.mac_ang(g.mac.mac_em_idx,k)) ...
                  .*g.mac.cur_im(g.mac.mac_em_idx,k);

            g.mac.curq(g.mac.mac_em_idx,k) = ...
                cos(g.mac.mac_ang(g.mac.mac_em_idx,k)) ...
                .*g.mac.cur_re(g.mac.mac_em_idx,k) ...
                + sin(g.mac.mac_ang(g.mac.mac_em_idx,k)) ...
                  .*g.mac.cur_im(g.mac.mac_em_idx,k);

            g.mac.curdg(g.mac.mac_em_idx,k) = ...
                g.mac.curd(g.mac.mac_em_idx,k).*g.mac.mac_pot(g.mac.mac_em_idx,1);

            g.mac.curqg(g.mac.mac_em_idx,k) = ...
                g.mac.curq(g.mac.mac_em_idx,k).*g.mac.mac_pot(g.mac.mac_em_idx,1);

            g.mac.dedprime(g.mac.mac_em_idx,k) = zeros(g.mac.n_em,1);
            g.mac.deqprime(g.mac.mac_em_idx,k) = zeros(g.mac.n_em,1);

            g.mac.ed(g.mac.mac_em_idx,k) = ...
                g.mac.edprime(g.mac.mac_em_idx,k) ...
                + g.mac.mac_con(g.mac.mac_em_idx,7).*g.mac.curqg(g.mac.mac_em_idx,k);

            g.mac.eq(g.mac.mac_em_idx,k) = ...
                g.mac.eqprime(g.mac.mac_em_idx,k) ...
                - g.mac.mac_con(g.mac.mac_em_idx,7).*g.mac.curdg(g.mac.mac_em_idx,k);

            g.mac.eterm(g.mac.mac_em_idx,k) = ...
                sqrt(g.mac.ed(g.mac.mac_em_idx,k).^2 ...
                     + g.mac.eq(g.mac.mac_em_idx,k).^2);

            g.mac.pelect(g.mac.mac_em_idx,k) = ...
                g.mac.eq(g.mac.mac_em_idx,k).*g.mac.curq(g.mac.mac_em_idx,k) ...
                + g.mac.ed(g.mac.mac_em_idx,k).*g.mac.curd(g.mac.mac_em_idx,k);

            g.mac.qelect(g.mac.mac_em_idx,k) = ...
                g.mac.eq(g.mac.mac_em_idx,k).*g.mac.curd(g.mac.mac_em_idx,k) ...
                - g.mac.ed(g.mac.mac_em_idx,k).*g.mac.curq(g.mac.mac_em_idx,k);

            g.mac.dmac_ang(g.mac.mac_em_idx,k) = ...
                g.sys.basrad ...
                *(g.mac.mac_spd(g.mac.mac_em_idx,k) - ones(g.mac.n_em,1));

            g.mac.dmac_spd(g.mac.mac_em_idx,k) = ...
                (g.mac.pmech(g.mac.mac_em_idx,k) ...
                 + g.mac.pm_sig(g.mac.mac_em_idx,k) ...
                 - g.mac.pelect(g.mac.mac_em_idx,k) ...
                   .*g.mac.mac_pot(g.mac.mac_em_idx,1) ...
                 - g.mac.mac_con(g.mac.mac_em_idx,17) ...
                   .*(g.mac.mac_spd(g.mac.mac_em_idx,k) - ones(g.mac.n_em,1))) ...
                ./(2*g.mac.mac_con(g.mac.mac_em_idx,16));
        end

        % end rate calculation
    end
end

end  % function end

% eof
