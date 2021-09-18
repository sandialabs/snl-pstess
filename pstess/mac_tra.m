function mac_tra(i,k,bus,flag)
% Syntax: mac_tra(i,k,bus,flag)
%
% Purpose: voltage-behind-transient-reactance generator
%          model, with vectorized computation option
%
% Input:   i - generator number
%            - 0 for vectorized computation
%          k - integer time
%          bus - solved loadflow bus data
%          flag - 0 - initialization
%                 1 - network interface computation
%                 2 - generator dynamics computation
%
% See Also: mac_em, mac_sub

%-----------------------------------------------------------------------------%
% Version history
%
% Version   2.0
% Date:     September 1997
% Purpose:  Saturation modified for low Eqprime
%
% Modified: G.J. Rogers May 1995
%
% Version:  1.0
% Author:   Joe H. Chow
% Date:     March 1991
%-----------------------------------------------------------------------------%

global g;  % declaring struct of global variables

if (g.mac.n_tra ~= 0)
    if (flag == 0)  % initialization
        if (i ~= 0)
            % check that machine i is a transient machine
            tra = find(g.mac.mac_tra_idx == i);
            if (length(tra) ~= 0)
                % make x'd = x'q
                if (g.mac.mac_con(i,7) ~= g.mac.mac_con(i,12))
                    g.mac.mac_con(i,12) = g.mac.mac_con(i,7);
                    %
                    wstr = '\nmac_tra: xqp made equal to xdp at generator %0.0f.';
                    warning(sprintf(wstr,i));
                end

                % check Tqo'
                if (g.mac.mac_con(i,14) == 0)
                    g.mac.mac_con(i,14) = 999.0;
                end

                busnum = g.bus.bus_int(g.mac.mac_con(i,2));  % bus number

                % mac_pot(i,1) -- scaled MVA base
                % mac_pot(i,2) -- base kV
                g.mac.mac_pot(i,1) = g.sys.basmva/g.mac.mac_con(i,3);
                g.mac.mac_pot(i,2) = 1.0;

                % eterm -- terminal bus voltage
                % theta -- terminal bus angle in radians
                g.mac.eterm(i,1) = bus(busnum,2);
                g.bus.theta(busnum,1) = bus(busnum,3)*pi/180;

                % pelect -- terminal bus angle in radians
                % qelect -- electrical output power, active
                g.mac.pelect(i,1) = bus(busnum,4)*g.mac.mac_con(i,22);
                g.mac.qelect(i,1) = bus(busnum,5)*g.mac.mac_con(i,23);

                % curr -- current magnitude
                % phi  -- power factor angle
                curr = sqrt(g.mac.pelect(i,1)^2+g.mac.qelect(i,1)^2) ...
                       /g.mac.eterm(i,1)*g.mac.mac_pot(i,1);
                phi = atan2(g.mac.qelect(i,1),g.mac.pelect(i,1));

                % v    -- complex voltage on the system reference frame
                % curr -- complex current on the system reference frame
                v = g.mac.eterm(i,1)*exp(1j*g.bus.theta(busnum,1));
                curr = curr*exp(1j*(g.bus.theta(busnum,1) - phi));

                eprime = v + (g.mac.mac_con(i,5) + 1j*g.mac.mac_con(i,7))*curr;
                ei = v + (g.mac.mac_con(i,5) + 1j*g.mac.mac_con(i,11))*curr;

                % mac_ang -- machine angle (delta)
                % mac_spd -- machine speed at steady state (omega)
                g.mac.mac_ang(i,1) = atan2(imag(ei),real(ei));
                g.mac.mac_spd(i,1) = 1;

                % reference frame rotation
                rot = 1j*exp(-1j*g.mac.mac_ang(i,1));
                g.mac.psi_re(i,1) = real(eprime);
                g.mac.psi_im(i,1) = imag(eprime);
                eprime = eprime*rot;
                g.mac.edprime(i,1) = real(eprime);
                g.mac.eqprime(i,1) = imag(eprime);
                curr = curr*rot;
                mcurmag = abs(curr);

                % pmech -- mech power = elec power + losses (on generator base)
                g.mac.pmech(i,1) = g.mac.pelect(i,1)*g.mac.mac_pot(i,1) ...
                                   + g.mac.mac_con(i,5)*mcurmag*mcurmag;

                g.mac.curdg(i,1) = real(curr);
                g.mac.curqg(i,1) = imag(curr);
                g.mac.curd(i,1) = real(curr)/g.mac.mac_pot(i,1);
                g.mac.curq(i,1) = imag(curr)/g.mac.mac_pot(i,1);
                v = v*rot;
                g.mac.ed(i,1) = real(v);
                g.mac.eq(i,1) = imag(v);

                % compute saturation
                inv_sat = inv([0.64 0.8 1; 1 1 1; 1.44 1.2 1]);
                b = [0.8 1+g.mac.mac_con(i,20) 1.2*(1+g.mac.mac_con(i,21))];
                g.mac.mac_pot(i,3) = b*inv_sat(1,:)';
                g.mac.mac_pot(i,4) = b*inv_sat(2,:)';
                g.mac.mac_pot(i,5) = b*inv_sat(3,:)';

                E_Isat = ...
                    g.mac.mac_pot(i,3)*g.mac.eqprime(i,1)^2 ...
                    + g.mac.mac_pot(i,4)*g.mac.eqprime(i,1) + g.mac.mac_pot(i,5);

                if (g.mac.eqprime(i,1) < 0.8)
                    E_Isat = g.mac.eqprime(i,1);
                end

                g.mac.vex(i,1) = E_Isat + (g.mac.mac_con(i,6) ...
                                           - g.mac.mac_con(i,7))*g.mac.curdg(i,1);

                g.mac.fldcur(i,1) = g.mac.vex(i,1);
            end
        else
            % vectorized computation
            % make xd' = xq'
            uexp_idx = find(g.mac.mac_con(g.mac.mac_tra_idx,7) ...
                            ~= g.mac.mac_con(g.mac.mac_tra_idx,12));
            if ~isempty(uexp_idx)
                g.mac.mac_con(g.mac.mac_tra_idx(uexp_idx),12) = ...
                    g.mac.mac_con(g.mac.mac_tra_idx(uexp_idx),7);
                %
                wstr = '\nmac_tra: xqp made equal to xdp at generator %0.0f.';
                warning(sprintf(wstr,g.mac.mac_tra_idx(uexp_idx)));
            end

            % make Tqo' non-zero
            notqp_idx = find(g.mac.mac_con(g.mac.mac_tra_idx,14) == 0);
            if ~isempty(notqp_idx)
                g.mac.mac_con(g.mac.mac_tra_idx(notqp_idx),14) = ...
                    999.0*ones(length(notqp_idx),1);
            end

            % busnum -- bus number vector
            busnum = g.bus.bus_int(g.mac.mac_con(g.mac.mac_tra_idx,2));

            % mac_pot(,1) -- scaled MVA base
            % mac_pot(,2) -- base kV
            g.mac.mac_pot(g.mac.mac_tra_idx,1) = ...
                g.sys.basmva*ones(g.mac.n_tra,1)./g.mac.mac_con(g.mac.mac_tra_idx,3);
            g.mac.mac_pot(g.mac.mac_tra_idx,2) = ones(g.mac.n_tra,1);

            % eterm -- terminal bus voltage
            % theta -- terminal bus angle in radians
            g.mac.eterm(g.mac.mac_tra_idx,1) = bus(busnum,2);
            g.bus.theta(busnum,1) = bus(busnum,3)*pi/180;

            % pelect -- electrical output power, active
            % qelect -- electrical output power, reactive
            g.mac.pelect(g.mac.mac_tra_idx,1) = ...
                bus(busnum,4).*g.mac.mac_con(g.mac.mac_tra_idx,22);
            g.mac.qelect(g.mac.mac_tra_idx,1) = ...
                bus(busnum,5).*g.mac.mac_con(g.mac.mac_tra_idx,23);

            % curr -- current magnitude
            curr = sqrt(g.mac.pelect(g.mac.mac_tra_idx,1).^2 ...
                        + g.mac.qelect(g.mac.mac_tra_idx,1).^2) ...
                   ./g.mac.eterm(g.mac.mac_tra_idx,1) ...
                   .*g.mac.mac_pot(g.mac.mac_tra_idx,1);

            % phi  -- power factor angle
            phi = atan2(g.mac.qelect(g.mac.mac_tra_idx,1), ...
                        g.mac.pelect(g.mac.mac_tra_idx,1));

            % v    -- complex voltage on the system reference frame
            % curr -- complex current on the system reference frame
            v = g.mac.eterm(g.mac.mac_tra_idx,1).*exp(1j*g.bus.theta(busnum,1));
            curr = curr.*exp(1j*(g.bus.theta(busnum,1) - phi));

            eprime = v + (g.mac.mac_con(g.mac.mac_tra_idx,5) ...
                          + 1j*g.mac.mac_con(g.mac.mac_tra_idx,7)).*curr;
            ei = v + (g.mac.mac_con(g.mac.mac_tra_idx,5) ...
                      + 1j*g.mac.mac_con(g.mac.mac_tra_idx,11)).*curr;

            % mac_ang -- machine angle (delta)
            % mac_spd -- machine speed at steady state (omega)
            g.mac.mac_ang(g.mac.mac_tra_idx,1) = atan2(imag(ei),real(ei));
            g.mac.mac_spd(g.mac.mac_tra_idx,1) = ones(g.mac.n_tra,1);

            % reference frame rotation
            rot = 1j*exp(-1j*g.mac.mac_ang(g.mac.mac_tra_idx,1));
            g.mac.psi_re(g.mac.mac_tra_idx,1) = real(eprime);
            g.mac.psi_im(g.mac.mac_tra_idx,1) = imag(eprime);
            eprime = eprime.*rot;
            g.mac.edprime(g.mac.mac_tra_idx,1) = real(eprime);
            g.mac.eqprime(g.mac.mac_tra_idx,1) = imag(eprime);
            curr = curr.*rot;
            mcurmag = abs(curr);

            % pmech -- mech power = elec power + losses (on generator base)
            g.mac.pmech(g.mac.mac_tra_idx,1) = ...
                g.mac.pelect(g.mac.mac_tra_idx,1) ...
                .*g.mac.mac_pot(g.mac.mac_tra_idx,1) ...
                + g.mac.mac_con(g.mac.mac_tra_idx,5).*mcurmag.*mcurmag;

            g.mac.curdg(g.mac.mac_tra_idx,1) = real(curr);
            g.mac.curqg(g.mac.mac_tra_idx,1) = imag(curr);

            g.mac.curd(g.mac.mac_tra_idx,1) = real(curr) ...
                                              ./g.mac.mac_pot(g.mac.mac_tra_idx,1);
            g.mac.curq(g.mac.mac_tra_idx,1) = imag(curr) ...
                                              ./g.mac.mac_pot(g.mac.mac_tra_idx,1);

            v = v.*rot;
            g.mac.ed(g.mac.mac_tra_idx,1) = real(v);
            g.mac.eq(g.mac.mac_tra_idx,1) = imag(v);

            % compute saturation
            inv_sat = inv([0.64 0.8 1; 1 1 1; 1.44 1.2 1]);

            b = [0.8*ones(g.mac.n_tra,1), ...
                 ones(g.mac.n_tra,1) + g.mac.mac_con(g.mac.mac_tra_idx,20), ...
                 1.2*(ones(g.mac.n_tra,1) + g.mac.mac_con(g.mac.mac_tra_idx,21))];

            g.mac.mac_pot(g.mac.mac_tra_idx,3) = b*inv_sat(1,:)';
            g.mac.mac_pot(g.mac.mac_tra_idx,4) = b*inv_sat(2,:)';
            g.mac.mac_pot(g.mac.mac_tra_idx,5) = b*inv_sat(3,:)';

            E_Isat = g.mac.mac_pot(g.mac.mac_tra_idx,3) ...
                     .*g.mac.eqprime(g.mac.mac_tra_idx,1).^2 ...
                     + g.mac.mac_pot(g.mac.mac_tra_idx,4) ...
                       .*g.mac.eqprime(g.mac.mac_tra_idx,1) ...
                     + g.mac.mac_pot(g.mac.mac_tra_idx,5);

            nosat_idx = find(g.mac.eqprime(g.mac.mac_tra_idx,1) < 0.8);
            if ~isempty(nosat_idx)
                E_Isat(nosat_idx)=g.mac.eqprime(g.mac.mac_tra_idx(nosat_idx),1);
            end

            g.mac.vex(g.mac.mac_tra_idx,1) = ...
                E_Isat + (g.mac.mac_con(g.mac.mac_tra_idx,6) ...
                          - g.mac.mac_con(g.mac.mac_tra_idx,7)) ...
                         .*g.mac.curdg(g.mac.mac_tra_idx,1);

            g.mac.fldcur(g.mac.mac_tra_idx,1) = g.mac.vex(g.mac.mac_tra_idx,1);

        end
        % end initialization
    end

    if (flag == 1)  % network interface computation
        if (i ~= 0)
            % check i for transient machine
            tra = find(g.mac.mac_tra_idx == i,1);
            if ~isempty(tra)
                % mac_ang -- wrt machine reference
                g.mac.mac_ang(i,k) = g.mac.mac_ang(i,k) - g.sys.mach_ref(k);

                % psr_re -- real part of psi
                % psr_im -- imaginary part of psi
                g.mac.psi_re(i,k) = sin(g.mac.mac_ang(i,k))*g.mac.edprime(i,k) ...
                                    + cos(g.mac.mac_ang(i,k))*g.mac.eqprime(i,k);
                g.mac.psi_im(i,k) = -cos(g.mac.mac_ang(i,k))*g.mac.edprime(i,k) ...
                                    + sin(g.mac.mac_ang(i,k))*g.mac.eqprime(i,k);
            end
        else  % vectorized computation
            % mac_ang -- wrt machine reference
            g.mac.mac_ang(g.mac.mac_tra_idx,k) = ...
                g.mac.mac_ang(g.mac.mac_tra_idx,k) ...
                - g.sys.mach_ref(k)*ones(g.mac.n_tra,1);

            % psi_re -- real part of psi
            % psi_im -- imaginary part of psi
            g.mac.psi_re(g.mac.mac_tra_idx,k) = ...
                sin(g.mac.mac_ang(g.mac.mac_tra_idx,k)) ...
                .*g.mac.edprime(g.mac.mac_tra_idx,k) ...
                + cos(g.mac.mac_ang(g.mac.mac_tra_idx,k)) ...
                  .*g.mac.eqprime(g.mac.mac_tra_idx,k);

            g.mac.psi_im(g.mac.mac_tra_idx,k) = ...
                -cos(g.mac.mac_ang(g.mac.mac_tra_idx,k)) ...
                 .*g.mac.edprime(g.mac.mac_tra_idx,k) ...
                + sin(g.mac.mac_ang(g.mac.mac_tra_idx,k)) ...
                  .*g.mac.eqprime(g.mac.mac_tra_idx,k);
        end

        % end of interface calculation
    end

    if (flag == 2 || flag == 3)  % generator dynamics calculation
        if (i ~= 0)
            % check that i is a transient machine
            tra = find(g.mac.mac_tra_idx == i);
            if (length(tra) ~= 0)
                % curd -- d-axis current
                % curq -- q-axis current
                g.mac.curd(i,k) = sin(g.mac.mac_ang(i,k))*g.mac.cur_re(i,k) ...
                                  - cos(g.mac.mac_ang(i,k))*g.mac.cur_im(i,k);
                g.mac.curq(i,k) = cos(g.mac.mac_ang(i,k))*g.mac.cur_re(i,k) ...
                                  + sin(g.mac.mac_ang(i,k))*g.mac.cur_im(i,k);

                g.mac.curdg(i,k) = g.mac.curd(i,k)*g.mac.mac_pot(i,1);
                g.mac.curqg(i,k) = g.mac.curq(i,k)*g.mac.mac_pot(i,1);

                E_Isat = g.mac.mac_pot(i,3)*g.mac.eqprime(i,k)^2 ...
                         + g.mac.mac_pot(i,4)*g.mac.eqprime(i,k) ...
                         + g.mac.mac_pot(i,5);

                if (g.mac.eqprime(i,1) < 0.8)
                    E_Isat = g.mac.eqprime(i,1);
                end

                g.mac.dedprime(i,k) = ...
                    (-g.mac.edprime(i,k) ...
                     + (g.mac.mac_con(i,11) - g.mac.mac_con(i,12)) ...
                       *g.mac.curqg(i,k))/g.mac.mac_con(i,14);

                g.mac.fldcur(i,k) = ...
                    E_Isat + (g.mac.mac_con(i,6) - g.mac.mac_con(i,7)) ...
                             *g.mac.curdg(i,k);

                g.mac.deqprime(i,k) = (g.mac.vex(i,k) - g.mac.fldcur(i,k)) ...
                                      /g.mac.mac_con(i,9);

                g.mac.ed(i,k) = g.mac.edprime(i,k) ...
                                - g.mac.mac_con(i,5)*g.mac.curdg(i,k) ...
                                + g.mac.mac_con(i,7)*g.mac.curqg(i,k);

                g.mac.eq(i,k) = g.mac.eqprime(i,k) ...
                                - g.mac.mac_con(i,5)*g.mac.curqg(i,k) ...
                                - g.mac.mac_con(i,7)*g.mac.curdg(i,k);

                g.mac.eterm(i,k) = sqrt(g.mac.ed(i,k)^2 + g.mac.eq(i,k)^2);

                g.mac.pelect(i,k) = g.mac.eq(i,k)*g.mac.curq(i,k) ...
                                    + g.mac.ed(i,k)*g.mac.curd(i,k);
                g.mac.qelect(i,k) = g.mac.eq(i,k)*g.mac.curd(i,k) ...
                                    - g.mac.ed(i,k)*g.mac.curq(i,k);

                curmag = abs(g.mac.curdg(i,k) + 1j*g.mac.curqg(i,k));

                Te = g.mac.pelect(i,k)*g.mac.mac_pot(i,1) ...
                     + g.mac.mac_con(i,5)*curmag*curmag;

                g.mac.dmac_ang(i,k) = g.sys.basrad*(g.mac.mac_spd(i,k) - 1);

                % JHC add missing pm_sig term per DKF (2011)
                g.mac.dmac_spd(i,k) = ...
                    (g.mac.pmech(i,k) + g.mac.pm_sig(i,k) - Te ...
                     - g.mac.mac_con(i,17)*(g.mac.mac_spd(i,k) - 1)) ...
                    /(2*g.mac.mac_con(i,16));
            end
        else  % vectorized computation
            g.mac.curd(g.mac.mac_tra_idx,k) = ...
                sin(g.mac.mac_ang(g.mac.mac_tra_idx,k)) ...
                .*g.mac.cur_re(g.mac.mac_tra_idx,k) ...
                - cos(g.mac.mac_ang(g.mac.mac_tra_idx,k)) ...
                  .*g.mac.cur_im(g.mac.mac_tra_idx,k);

            g.mac.curq(g.mac.mac_tra_idx,k) = ...
                cos(g.mac.mac_ang(g.mac.mac_tra_idx,k)) ...
                .*g.mac.cur_re(g.mac.mac_tra_idx,k) ...
                + sin(g.mac.mac_ang(g.mac.mac_tra_idx,k)) ...
                  .*g.mac.cur_im(g.mac.mac_tra_idx,k);

            g.mac.curdg(g.mac.mac_tra_idx,k) = ...
                g.mac.curd(g.mac.mac_tra_idx,k).*g.mac.mac_pot(g.mac.mac_tra_idx,1);

            g.mac.curqg(g.mac.mac_tra_idx,k) = ...
                g.mac.curq(g.mac.mac_tra_idx,k).*g.mac.mac_pot(g.mac.mac_tra_idx,1);

            E_Isat = g.mac.mac_pot(g.mac.mac_tra_idx,3) ...
                     .*g.mac.eqprime(g.mac.mac_tra_idx,k).^2 ...
                     + g.mac.mac_pot(g.mac.mac_tra_idx,4) ...
                       .*g.mac.eqprime(g.mac.mac_tra_idx,k) ...
                     + g.mac.mac_pot(g.mac.mac_tra_idx,5);

            nosat_idx = find(g.mac.eqprime(g.mac.mac_tra_idx,1) < 0.8);
            if ~isempty(nosat_idx)
                E_Isat(nosat_idx) = g.mac.eqprime(g.mac.mac_tra_idx(nosat_idx),k);
            end

            g.mac.fldcur(g.mac.mac_tra_idx,k) = ...
                E_Isat + (g.mac.mac_con(g.mac.mac_tra_idx,6) ...
                          - g.mac.mac_con(g.mac.mac_tra_idx,7)) ...
                         .*g.mac.curdg(g.mac.mac_tra_idx,k);

            g.mac.deqprime(g.mac.mac_tra_idx,k) = ...
                (g.mac.vex(g.mac.mac_tra_idx,k) ...
                 - g.mac.fldcur(g.mac.mac_tra_idx,k)) ...
                ./g.mac.mac_con(g.mac.mac_tra_idx,9);

            g.mac.dedprime(g.mac.mac_tra_idx,k) = ...
                (-g.mac.edprime(g.mac.mac_tra_idx,k) ...
                 + (g.mac.mac_con(g.mac.mac_tra_idx,11) ...
                    - g.mac.mac_con(g.mac.mac_tra_idx,12)) ...
                   .*g.mac.curqg(g.mac.mac_tra_idx,k)) ...
                ./g.mac.mac_con(g.mac.mac_tra_idx,14);

            g.mac.ed(g.mac.mac_tra_idx,k) = ...
                g.mac.edprime(g.mac.mac_tra_idx,k) ...
                - g.mac.mac_con(g.mac.mac_tra_idx,5) ...
                  .*g.mac.curdg(g.mac.mac_tra_idx,k) ...
                + g.mac.mac_con(g.mac.mac_tra_idx,7) ...
                  .*g.mac.curqg(g.mac.mac_tra_idx,k);

            g.mac.eq(g.mac.mac_tra_idx,k) = ...
                g.mac.eqprime(g.mac.mac_tra_idx,k) ...
                - g.mac.mac_con(g.mac.mac_tra_idx,5) ...
                  .*g.mac.curqg(g.mac.mac_tra_idx,k) ...
                - g.mac.mac_con(g.mac.mac_tra_idx,7) ...
                  .*g.mac.curdg(g.mac.mac_tra_idx,k);

            g.mac.eterm(g.mac.mac_tra_idx,k) = ...
                sqrt(g.mac.ed(g.mac.mac_tra_idx,k).^2 ...
                     + g.mac.eq(g.mac.mac_tra_idx,k).^2);

            g.mac.pelect(g.mac.mac_tra_idx,k) = ...
                g.mac.eq(g.mac.mac_tra_idx,k).*g.mac.curq(g.mac.mac_tra_idx,k) ...
                + g.mac.ed(g.mac.mac_tra_idx,k).*g.mac.curd(g.mac.mac_tra_idx,k);

            g.mac.qelect(g.mac.mac_tra_idx,k) = ...
                g.mac.eq(g.mac.mac_tra_idx,k).*g.mac.curd(g.mac.mac_tra_idx,k) ...
                - g.mac.ed(g.mac.mac_tra_idx,k).*g.mac.curq(g.mac.mac_tra_idx,k);

            curmag = abs(g.mac.curdg(g.mac.mac_tra_idx,k) ...
                         + 1j*g.mac.curqg(g.mac.mac_tra_idx,k));

            Te = g.mac.pelect(g.mac.mac_tra_idx,k) ...
                 .*g.mac.mac_pot(g.mac.mac_tra_idx,1) ...
                 + g.mac.mac_con(g.mac.mac_tra_idx,5).*curmag.*curmag;

            g.mac.dmac_ang(g.mac.mac_tra_idx,k) = ...
                g.sys.basrad ...
                .*(g.mac.mac_spd(g.mac.mac_tra_idx,k) - ones(g.mac.n_tra,1));

            g.mac.dmac_spd(g.mac.mac_tra_idx,k) = ...
                (g.mac.pmech(g.mac.mac_tra_idx,k) ...
                 + g.mac.pm_sig(g.mac.mac_tra_idx,k) - Te ...
                 - g.mac.mac_con(g.mac.mac_tra_idx,17) ...
                   .*(g.mac.mac_spd(g.mac.mac_tra_idx,k) - ones(g.mac.n_tra,1))) ...
                ./(2*g.mac.mac_con(g.mac.mac_tra_idx,16));
        end

        % end calculation of rates of change
    end
end

end  % function end

% eof
