function [bus_new] = mac_ind(i,k,bus,flag)
% Syntax: [bus_new] = mac_ind(i,k,bus,flag)
%
% Purpose: Simple Induction Motor Model (single cage)
%
% Input:   i is the motor number, 0 for vectorized computation
%          k is the time step
%          bus is the solved load flow bus data
%          flag is 0 for initialization
%                  1 for network interface
%                  2 for dertermination of rates of change of states
%                  3 for formation of linearized state matrix
%
% Output:  bus_new is bus with the power and reactive power loads
%          modified to subtract the motor loads
%          modification is made only when the motors are initialized
%          i.e. flag = 0
%
% Note:    No leakage inductance saturation
%
% Data format ind_con
%          1 - motor number
%          2 - busnumber
%          3 - base MVA
%          4 - rs
%          5 - xs -stator leakage reactance
%          6 - Xm - magnetizing reactance
%          7 - rr
%          8 - xr - rotor leakage reactance
%          9 - H  - inertia constant motor + load in sec
%         15 - fraction of bus load power taken by motor
%
% Note:    If entry 15 is zero, it is assumed that the motor is to be
%          started on the specified bus

%-----------------------------------------------------------------------------%
% Version history
%
% Author Graham Rogers
% Date November 1995
%-----------------------------------------------------------------------------%

global g;  % declaring struct of global variables

bus_new = bus;

if ~isempty(g.ind.ind_con)
    if (flag == 0)
        % initialization
        if (i == 0)
            % vector computation

            g.ind.motbus = g.bus.bus_int(g.ind.ind_con(:,2));

            % ind_pot(,1) -- scaled mva base
            % ind_pot(,2) -- base kv
            g.ind.ind_pot(:,1) = g.sys.basmva./g.ind.ind_con(:,3);
            g.ind.ind_pot(:,2) = ones(g.ind.n_mot,1);

            % mot_vm -- motor terminal voltage mag
            % mot_ang -- motor terminal voltage angle
            mot_vm(:,1) = bus(g.ind.motbus,2);
            mot_ang(:,1) = bus(g.ind.motbus,3)*pi/180;

            % complex voltage
            v = mot_vm(:,1).*exp(1j*mot_ang(:,1));
            g.ind.vdmot(:,1) = real(v);
            g.ind.vqmot(:,1) = imag(v);

            % p_mot -- motor power demand
            g.ind.p_mot(:,1) = bus(g.ind.motbus,6).*g.ind.ind_con(:,15);

            % modify bus load power
            bus_new(g.ind.motbus,6) = bus(g.ind.motbus,6) - g.ind.p_mot(:,1);

            % index of motors to be initialized for running
            % assumes motor starting is power fraction zero
            run_ind = find(g.ind.ind_con(:,15) ~= 0);

            % ind_pot(,3) -- Xs
            % ind_pot(,4) -- Xr
            g.ind.ind_pot(:,3) = g.ind.ind_con(:,5) + g.ind.ind_con(:,6);
            g.ind.ind_pot(:,4) = g.ind.ind_con(:,8) + g.ind.ind_con(:,6);

            % ind_pot(,5) -- Xsp
            % ind_pot(,6) -- Xs-Xsp
            % ind_pot(,7) -- 1/Tr
            g.ind.ind_pot(:,5) = ...
                g.ind.ind_con(:,5) ...
                + g.ind.ind_con(:,6).*g.ind.ind_con(:,8)./g.ind.ind_pot(:,4);

            g.ind.ind_pot(:,6) = g.ind.ind_pot(:,3) - g.ind.ind_pot(:,5);

            g.ind.ind_pot(:,7) = g.sys.basrad*g.ind.ind_con(:,7)./g.ind.ind_pot(:,4);

            rs = g.ind.ind_con(:,4);
            xs = g.ind.ind_con(:,5);
            Xm = g.ind.ind_con(:,6);
            rr = g.ind.ind_con(:,7);
            xr = g.ind.ind_con(:,8);

            % find initial slip
            slip_old = zeros(g.ind.n_mot,1);
            slip_new = ones(g.ind.n_mot,1);

            % set defaults for motor starting
            imot = zeros(g.ind.n_mot,1);
            pem = zeros(g.ind.n_mot,1);
            qem = zeros(g.ind.n_mot,1);

            g.ind.vdp(:,1) = zeros(g.ind.n_mot,1);
            g.ind.vqp(:,1) = zeros(g.ind.n_mot,1);

            % default for motor starting
            g.ind.t_init = ones(g.ind.n_mot,1);

            % Newton-Raphson iteration to determine initial slip
            motrun = length(run_ind);  % number of running motors
            if (motrun ~= 0)  % check that some motors are running
                rsrun = g.ind.ind_con(run_ind,4);
                xsrun = g.ind.ind_con(run_ind,5);
                Xmrun = g.ind.ind_con(run_ind,6);
                rrrun = g.ind.ind_con(run_ind,7);
                xrrun = g.ind.ind_con(run_ind,8);

                iter = 0;
                itermax = 50;
                err = max(abs(slip_new - slip_old));
                while ((err >= 1e-8) && (iter <= itermax))
                    iter = iter + 1;
                    y = g.sys.basrad.*slip_old(run_ind)./g.ind.ind_pot(run_ind,7);

                    denom = ones(motrun,1) + y.*y;
                    zr = rsrun + y.*g.ind.ind_pot(run_ind,6)./denom;
                    zi = g.ind.ind_pot(run_ind,5) + g.ind.ind_pot(run_ind,6)./denom;

                    dzr = g.ind.ind_pot(run_ind,6) ...
                          .*(ones(motrun,1) - y.*y)./denom./denom;

                    dzi = -2*g.ind.ind_pot(run_ind,6).*y./denom./denom;

                    zmod2 = zr.*zr + zi.*zi;

                    dp = v(run_ind).*conj(v(run_ind)) ...
                         .*(dzr.*zmod2 - 2*zr.*(dzr.*zr+dzi.*zi));

                    dp = dp./zmod2./zmod2;

                    pem(run_ind) = v(run_ind).*conj(v(run_ind)).*zr./zmod2;

                    ynew = y - (pem(run_ind) - g.ind.p_mot(run_ind) ...
                                               .*g.ind.ind_pot(run_ind,1))./dp;

                    % define the mismatch and update the slip
                    slip_new(run_ind) = ynew.*g.ind.ind_pot(run_ind,7)/g.sys.basrad;
                    err = max(abs(slip_new - slip_old));

                    slip_old = slip_new;
                end

                if (iter > itermax)
                    error('mac_ind: motor slip calculation failed to converge.');
                end
            end

            g.ind.slip(:,1) = slip_new;
            y = g.sys.basrad*g.ind.slip(:,1)./g.ind.ind_pot(:,7);

            denom = ones(g.ind.n_mot,1) + y.*y;
            zr = rs + y.*g.ind.ind_pot(:,6)./denom;
            zi = g.ind.ind_pot(:,5) + g.ind.ind_pot(:,6)./denom;

            imot(run_ind) = v(run_ind)./(zr(run_ind) + 1j*zi(run_ind));
            smot(run_ind) = v(run_ind).*conj(imot(run_ind));
            pem(run_ind) = real(smot(run_ind));
            qem(run_ind) = imag(smot(run_ind));

            % complex initial rotor states
            vp = v - (rs + 1j*g.ind.ind_pot(:,5)).*imot;
            g.ind.vdp(run_ind,1) = real(vp(run_ind));
            g.ind.vqp(run_ind,1) = imag(vp(run_ind));

            % initial motor torque
            g.ind.t_init(run_ind) = real(vp(run_ind).*conj(imot(run_ind)));
            ind_ldto(0,1);
            lzero_ind = find(g.ind.tload(:,1) == 0);
            if (length(lzero_ind) ~= 0)
                g.ind.tload(lzero_ind,1) = ones(length(lzero_ind),1);
            end

            g.ind.t_init(run_ind) = g.ind.t_init(run_ind)./g.ind.tload(run_ind,1);
            g.ind.idmot(:,1) = real(imot)./g.ind.ind_pot(:,1);
            g.ind.iqmot(:,1) = imag(imot)./g.ind.ind_pot(:,1);

            % modify qload
            bus_new(g.ind.motbus,7) = bus(g.ind.motbus,7) - qem./g.ind.ind_pot(:,1);
        else
            % non-vectorized initialization

            g.ind.motbus = g.bus.bus_int(g.ind.ind_con(i,2));

            % ind_pot(,1) -- scaled mva base
            % ind_pot(,2) -- base kv
            g.ind.ind_pot(i,1) = g.sys.basmva/g.ind.ind_con(i,3);
            g.ind.ind_pot(i,2) = 1.0;

            % mot_vm -- motor terminal voltage mag
            % mot_ang -- motor terminal voltage angle
            mot_vm(i,1) = bus(g.ind.motbus,2);
            mot_ang(i,1) = bus(g.ind.motbus,3)*pi/180;

            % complex voltage
            v = mot_vm(i,1)*exp(1j*mot_ang(i,1));
            g.ind.vdmot(i,1) = real(v);
            g.ind.vqmot(i,1) = imag(v);

            % p_mot -- motor power demand
            g.ind.p_mot(i,1) = bus(g.ind.motbus,6)*g.ind.ind_con(i,15);

            % modify motor bus load
            bus_new(g.ind.motbus,6) = bus(g.ind.motbus,6) - g.ind.p_mot(i,1);

            % ind_pot(,3) -- Xs
            % ind_pot(,4) -- Xr
            g.ind.ind_pot(i,3) = g.ind.ind_con(i,5) + g.ind.ind_con(i,6);
            g.ind.ind_pot(i,4) = g.ind.ind_con(i,8) + g.ind.ind_con(i,6);

            % ind_pot(,5) -- Xsp
            % ind_pot(,6) -- Xs-Xsp
            % ind_pot(,7) -- 1/Tr
            g.ind.ind_pot(i,5) = ...
                g.ind.ind_con(i,5) ...
                + g.ind.ind_con(i,6)*g.ind.ind_con(i,8)/g.ind.ind_pot(i,4);

            g.ind.ind_pot(i,6) = g.ind.ind_pot(i,3) - g.ind.ind_pot(i,5);

            g.ind.ind_pot(i,7) = g.sys.basrad*g.ind.ind_con(i,7)/g.ind.ind_pot(i,4);

            rs = g.ind.ind_con(i,4);
            xs = g.ind.ind_con(i,5);
            Xm = g.ind.ind_con(i,6);
            rr = g.ind.ind_con(i,7);
            xr = g.ind.ind_con(i,8);

            % find initial slip
            slip_old = 0.005;
            slip_new = 1;

            g.ind.slip(i,1) = 1;
            g.ind.vdp(i,1) = 0;
            g.ind.vqp(i,1) = 0;
            g.ind.t_init(i) = 1;
            g.ind.idmot(i,1) = 0;
            g.ind.iqmot(i,1) = 0;

            err = abs(slip_old - slip_new);
            if (g.ind.ind_con(i,15) > 1e-6)  % motor running
                % Newton-Raphson iteration to find initial slip
                iter = 0;
                itermax = 50;
                while ((err >= 1e-8) && (iter <= itermax))
                    iter = iter + 1;
                    y = g.sys.basrad*slip_old/g.ind.ind_pot(i,7);

                    denom = 1 + y*y;
                    zr = rs + y.*g.ind.ind_pot(i,6)/denom;
                    zi = g.ind.ind_pot(i,5) + g.ind.ind_pot(i,6)/denom;

                    dzr = g.ind.ind_pot(i,6)*(1-y*y)/denom/denom;
                    dzi = -2*g.ind.ind_pot(i,6)*y/denom/denom;

                    modz2 = zr*zr + zi*zi;
                    dp = v*conj(v)*(dzr*modz2 - 2*zr*(dzr*zr+dzi*zi));
                    dp = dp/modz2/modz2;

                    pem = v*conj(v)*zr/modz2;
                    ynew = y - (pem - g.ind.p_mot*g.ind.ind_pot(i,1))/dp;

                    % define the mismatch and update the slip
                    slip_new = ynew/g.sys.basrad*g.ind.ind_pot(i,7);
                    err = abs(slip_old - slip_new);

                    slip_old = slip_new;
                end

                if (iter > itermax)
                    error('mac_ind: motor slip calculation failed to converge.');
                end

                g.ind.slip(i,1) = slip_new;
                y = g.sys.basrad*g.ind.slip(i,1)/g.ind.ind_pot(i,7);

                denom = 1 + y*y;
                zr = rs + y*g.ind.ind_pot(i,6)/denom;
                zi = g.ind.ind_pot(i,5) + g.ind.ind_pot(i,6)/denom;

                imot = v/(zr + 1j*zi);  % compex motor current
                smot = v*conj(imot);
                pem = real(smot);
                qem = imag(smot);

                % complex initial rotor states
                vp = v - (rs + 1j*g.ind.ind_pot(i,5))*imot;

                g.ind.vdp(i,1) = real(vp);
                g.ind.vqp(i,1) = imag(vp);

                g.ind.t_init(i) = real(vp*conj(imot));
                ind_ldto(i,1);  % find the motor load torque
                if (g.ind.tload(i,1) == 0)
                    error('mac_ind: undefined load speed characteristic.');
                end

                % set the load multiplier
                g.ind.t_init(i) = g.ind.t_init(i)/g.ind.tload(i,1);

                % motor currents on system base
                g.ind.idmot(i,1) = real(imot)/g.ind.ind_pot(i,1);
                g.ind.iqmot(i,1) = imag(imot)/g.ind.ind_pot(i,1);

                % modify qload
                bus_new(g.ind.motbus,7) = bus(g.ind.motbus,7) ...
                                          - qem/g.ind.ind_pot(i,1);
            end
        end
    end

    if (flag == 1)
        % network interface
        % no interface required for induction motor
    end

    if (flag == 2)
        % motor dynamics calculation
        if (i == 0)
            % vector calculation

            ind_ldto(0,k);

            % convert to machine base
            idm = g.ind.idmot(:,k).*g.ind.ind_pot(:,1);
            iqm = g.ind.iqmot(:,k).*g.ind.ind_pot(:,1);

            % Brereton, Lewis and Young motor model
            g.ind.dvdp(:,k) = -(iqm.*g.ind.ind_pot(:,6) ...
                                + g.ind.vdp(:,k)).*g.ind.ind_pot(:,7) ...
                              + g.ind.vqp(:,k).*g.ind.slip(:,k)*g.sys.basrad;

            g.ind.dvqp(:,k) = (idm.*g.ind.ind_pot(:,6) ...
                               - g.ind.vqp(:,k)).*g.ind.ind_pot(:,7) ...
                              - g.ind.vdp(:,k).*g.ind.slip(:,k)*g.sys.basrad;

            g.ind.dslip(:,k) = (g.ind.tload(:,k).*g.ind.t_init(:) ...
                                - g.ind.vdp(:,k).*idm - g.ind.vqp(:,k).*iqm)/2 ...
                               ./g.ind.ind_con(:,9);
        else
            ind_ldto(i,k);

            idm = g.ind.idmot(i,k)*g.ind.ind_pot(i,1);
            iqm = g.ind.iqmot(i,k)*g.ind.ind_pot(i,1);

            g.ind.dvdp(i,k) = -(iqm*g.ind.ind_pot(i,6) ...
                                + g.ind.vdp(i,k))*g.ind.ind_pot(i,7) ...
                              + g.ind.vqp(i,k)*g.ind.slip(i,k)*g.sys.basrad;

            g.ind.dvqp(i,k) = (idm*g.ind.ind_pot(i,6) ...
                               - g.ind.vqp(i,k))*g.ind.ind_pot(i,7) ...
                              - g.ind.vdp(i,k)*g.ind.slip(i,k)*g.sys.basrad;

            g.ind.dslip(i,k) = (g.ind.tload(i,k)*g.ind.t_init(i) ...
                                - g.ind.vdp(i,k)*idm - g.ind.vqp(i,k)*iqm)/2 ...
                               /g.ind.ind_con(i,9);
        end
    end

    if (flag == 3)
        % linearize
        % add code later
    end
end

end  % function end

% eof
