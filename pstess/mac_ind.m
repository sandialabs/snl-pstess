function bus_new = mac_ind(i,k,bus,flag)
% Syntax: [bus_new] = mac_ind(i,k,bus,flag)
%
% Purpose: Induction motor model
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
% Data format ind_con
%          1 - motor number
%          2 - busnumber
%          3 - base MVA
%          4 - rs
%          5 - xs, stator leakage reactance
%          6 - Xm, magnetizing reactance
%          7 - rr
%          8 - xr, rotor leakage reactance
%          9 - H, inertia constant motor + load in sec
%         10 - r2, double cage resistance
%         11 - x2, intercage reactance
%         12 - dbf, deep bar factor
%         13 - isat, current at which leakage inductance starts to saturate
%         15 - fraction of bus load power taken by motor
%
% Note:    If entry 15 is zero, it is assumed that the motor is to be
%          started on the specified bus

%-----------------------------------------------------------------------------%
% Version history
%
% Version: 2.0
% Purpose: added deep bar, double cage and leakage inductance saturation
% Date:    June 2002
% Author:  Graham Rogers
%
% Version: 1.0 (initial version)
% Author:  Graham Rogers
% Date:    November 1995
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

            rs = g.ind.ind_con(:,4);
            xs = g.ind.ind_con(:,5);
            Xm = g.ind.ind_con(:,6);
            rr = g.ind.ind_con(:,7);
            xr = g.ind.ind_con(:,8);
            rr2 = g.ind.ind_con(:,10);
            xr2 = g.ind.ind_con(:,11);
            dbf = g.ind.ind_con(:,12);
            isat = g.ind.ind_con(:,13);

            % dbc_idx -- motors with double cage
            % db_idx  -- motors with deep bars
            % sat_idx -- motors with leakage inductance saturation
            g.ind.dbc_idx = find(rr2 ~= 0);
            g.ind.db_idx = find(dbf ~= 0);
            g.ind.sat_idx = find(isat ~= 0);

            % ind_pot(,3) -- Xs
            % ind_pot(,4) -- Xr
            g.ind.ind_pot(:,3) = xs + Xm;
            g.ind.ind_pot(:,4) = xr + Xm;

            % ind_pot(,5) -- Xsp
            % ind_pot(,6) -- Xs-Xsp
            % ind_pot(,7) -- 1/Tr
            g.ind.ind_pot(:,5) = xs + Xm.*xr./g.ind.ind_pot(:,4);
            g.ind.ind_pot(:,6) = g.ind.ind_pot(:,3) - g.ind.ind_pot(:,5);
            g.ind.ind_pot(:,7) = g.sys.basrad*rr./g.ind.ind_pot(:,4);

            % run_ind -- index of running motors
            % motrun -- number of running motors
            run_ind = find(g.ind.ind_con(:,15) ~= 0);
            motrun = length(run_ind);

            % start_ind -- index of starting motors
            % motstart -- number of starting motors
            start_ind = find(g.ind.ind_con(:,15) == 0);
            motstart = length(start_ind);

            % assumes motor starting if power fraction zero
            % find initial slip
            slip_old = zeros(g.ind.n_mot,1);
            slip_new = ones(g.ind.n_mot,1);

            % reset ind_pot for double cage and deepbar rotor machines
            s = 0.01*ones(g.ind.n_mot,1);
            s(start_ind) = ones(motstart,1);
            if ~isempty(g.ind.dbc_idx)
                [rdc,xdc] = ...
                    dbcage(rr(g.ind.dbc_idx),xr(g.ind.dbc_idx), ...
                           rr2(g.ind.dbc_idx),xr2(g.ind.dbc_idx),s(g.ind.dbc_idx));

                g.ind.ind_pot(g.ind.dbc_idx,4) = Xm(g.ind.dbc_idx) + xdc;

                g.ind.ind_pot(g.ind.dbc_idx,5) = xs(g.ind.dbc_idx) ...
                                                 + Xm(g.ind.dbc_idx).*xdc ...
                                                   ./g.ind.ind_pot(g.ind.dbc_idx,4);

                g.ind.ind_pot(g.ind.dbc_idx,6) = g.ind.ind_pot(g.ind.dbc_idx,3) ...
                                                 - g.ind.ind_pot(g.ind.dbc_idx,5);

                g.ind.ind_pot(g.ind.dbc_idx,7) = g.sys.basrad*rdc ...
                                                 ./g.ind.ind_pot(g.ind.dbc_idx,4);
            end

            if ~isempty(g.ind.db_idx)
                [rdb,xdb] = ...
                    deepbar(rr(g.ind.db_idx),dbf(g.ind.db_idx),s(g.ind.db_idx));

                g.ind.ind_pot(g.ind.db_idx,4) = ...
                    Xm(g.ind.db_idx) + xr(g.ind.db_idx) + xdb;

                g.ind.ind_pot(g.ind.db_idx,5) = ...
                    xs(g.ind.db_idx) ...
                    + Xm(g.ind.db_idx).*(xr(g.ind.db_idx) + xdb) ...
                      ./g.ind.ind_pot(g.ind.db_idx,4);

                g.ind.ind_pot(g.ind.db_idx,6) = g.ind.ind_pot(g.ind.db_idx,3) ...
                                                - g.ind.ind_pot(g.ind.db_idx,5);

                g.ind.ind_pot(g.ind.db_idx,7) = g.sys.basrad*rdb ...
                                                ./g.ind.ind_pot(g.ind.db_idx,4);
            end

            % Set defaults for motor starting
            imot = zeros(g.ind.n_mot,1);
            pem = zeros(g.ind.n_mot,1);
            qem = zeros(g.ind.n_mot,1);

            g.ind.vdp(:,1) = zeros(g.ind.n_mot,1);
            g.ind.vqp(:,1) = zeros(g.ind.n_mot,1);

            vp = imot;

            % default for motor starting
            g.ind.t_init = ones(g.ind.n_mot,1);

            % Newton-Raphson iteration to determine initial slip
            if (motrun ~= 0)  % check that some motors are running
                iter = 0;
                itermax = 50;
                err = max(abs(slip_new - slip_old));
                while ((err >= 1e-8) && (iter <= itermax))
                    iter = iter + 1;
                    y = g.sys.basrad.*slip_old(run_ind)./g.ind.ind_pot(run_ind,7);

                    denom = ones(motrun,1) + y.*y;
                    zr = rs(run_ind) + y.*g.ind.ind_pot(run_ind,6)./denom;
                    zi = g.ind.ind_pot(run_ind,5) + g.ind.ind_pot(run_ind,6)./denom;

                    dzr = g.ind.ind_pot(run_ind,6) ...
                          .*(ones(motrun,1) - y.*y)./denom./denom;

                    dzi = -2*g.ind.ind_pot(run_ind,6).*y./denom./denom;

                    zmod2 = zr.*zr + zi.*zi;

                    dp = v(run_ind).*conj(v(run_ind)) ...
                         .*(dzr.*zmod2 - 2*zr.*(dzr.*zr + dzi.*zi));

                    dp = dp./zmod2./zmod2;

                    pem(run_ind) = v(run_ind).*conj(v(run_ind)).*zr./zmod2;

                    ynew = y - (pem(run_ind) - g.ind.p_mot(run_ind,1) ...
                                               .*g.ind.ind_pot(run_ind,1))./dp;

                    % define the mismatch and update the slip
                    slip_new(run_ind) = ynew.*g.ind.ind_pot(run_ind,7)/g.sys.basrad;
                    err = max(abs(slip_new - slip_old));

                    slip_old = slip_new;

                    % checking for deepbar induction machines
                    if ~isempty(g.ind.dbc_idx)
                        [rdc,xdc] = dbcage(rr(g.ind.dbc_idx),xr(g.ind.dbc_idx), ...
                                           rr2(g.ind.dbc_idx),xr2(g.ind.dbc_idx), ...
                                           slip_new(g.ind.dbc_idx));

                        g.ind.ind_pot(g.ind.dbc_idx,4) = Xm(g.ind.dbc_idx) + xdc;

                        g.ind.ind_pot(g.ind.dbc_idx,5) = ...
                            xs(g.ind.dbc_idx) ...
                            + Xm(g.ind.dbc_idx).*xdc./g.ind.ind_pot(g.ind.dbc_idx,4);

                        g.ind.ind_pot(g.ind.dbc_idx,6) = ...
                            g.ind.ind_pot(g.ind.dbc_idx,3) ...
                            - g.ind.ind_pot(g.ind.dbc_idx,5);

                        g.ind.ind_pot(g.ind.dbc_idx,7) = ...
                            g.sys.basrad*rdc./g.ind.ind_pot(g.ind.dbc_idx,4);
                    end

                    if ~isempty(g.ind.db_idx)
                        [rdb,xdb] = deepbar(rr(g.ind.db_idx),dbf(g.ind.db_idx), ...
                                            slip_new(g.ind.db_idx));

                        g.ind.ind_pot(g.ind.db_idx,4) = ...
                            Xm(g.ind.db_idx) + xr(g.ind.db_idx) + xdb;

                        g.ind.ind_pot(g.ind.db_idx,5) = ...
                            xs(g.ind.db_idx) ...
                            + Xm(g.ind.db_idx).*(xr(g.ind.db_idx) + xdb) ...
                              ./g.ind.ind_pot(g.ind.db_idx,4);

                        g.ind.ind_pot(g.ind.db_idx,6) = ...
                            g.ind.ind_pot(g.ind.db_idx,3) ...
                            - g.ind.ind_pot(g.ind.db_idx,5);

                        g.ind.ind_pot(g.ind.db_idx,7) = ...
                            g.sys.basrad*rdb./g.ind.ind_pot(g.ind.db_idx,4);
                    end
                end

                if (iter > itermax)
                    error('mac_ind: motor slip calculation failed to converge.');
                end
            end

            g.ind.slip(:,1) = slip_new;

            % induction motor load torque calculation
            ind_ldto(0,1);
            y = g.sys.basrad*g.ind.slip(:,1)./g.ind.ind_pot(:,7);

            denom = ones(g.ind.n_mot,1) + y.*y;
            zr = rs + y.*g.ind.ind_pot(:,6)./denom;
            zi = g.ind.ind_pot(:,5) + g.ind.ind_pot(:,6)./denom;
            if ~isempty(run_ind)
                imot(run_ind) = v(run_ind)./(zr(run_ind) + 1j*zi(run_ind));
                smot(run_ind) = v(run_ind).*conj(imot(run_ind));
                pem(run_ind) = real(smot(run_ind));
                qem(run_ind) = imag(smot(run_ind));

                % complex initial rotor states
                vp(run_ind) = v(run_ind) ...
                              - (rs(run_ind) ...
                                 + 1j*g.ind.ind_pot(run_ind,5)).*imot(run_ind);

                g.ind.vdp(run_ind,1) = real(vp(run_ind));
                g.ind.vqp(run_ind,1) = imag(vp(run_ind));
            end

            g.ind.idmot(:,1) = real(imot)./g.ind.ind_pot(:,1);
            g.ind.iqmot(:,1) = imag(imot)./g.ind.ind_pot(:,1);

            % modify qload
            bus_new(g.ind.motbus,7) = bus(g.ind.motbus,7)-qem./g.ind.ind_pot(:,1);
            tlm = g.ind.vdp(:,k).*real(imot) + g.ind.vqp(:,k).*imag(imot);
            trat = tlm./g.ind.tload(:,1);

            % modify load specification to correct the initial load
            g.ind.mld_con(run_ind,[3 5]) = diag(trat(run_ind)) ...
                                           *g.ind.mld_con(run_ind,[3 5]);
        else
            error('mac_ind: initialization must be vectorized.');
        end
    end

    if (flag == 1)
        v = g.bus.bus_v(g.ind.motbus,k);
        g.ind.vdmot(:,k) = real(v);
        g.ind.vqmot(:,k) = imag(v);
    end

    if (flag == 2)
        % motor dynamics calculation
        if (i == 0)
            % vector calculation

            ind_ldto(0,k);

            % convert to machine base
            idm = g.ind.idmot(:,k).*g.ind.ind_pot(:,1);
            iqm = g.ind.iqmot(:,k).*g.ind.ind_pot(:,1);

            rs = g.ind.ind_con(:,4);
            xs = g.ind.ind_con(:,5);
            Xm = g.ind.ind_con(:,6);
            rr = g.ind.ind_con(:,7);
            xr = g.ind.ind_con(:,8);
            rr2 = g.ind.ind_con(:,10);
            xr2 = g.ind.ind_con(:,11);
            dbf = g.ind.ind_con(:,12);

            imot = abs(idm + 1j*iqm);
            if ~isempty(g.ind.sat_idx)
                % saturation of leakage inductance
                ism = imot(g.ind.sat_idx);
                isat = g.ind.ind_con(g.ind.sat_idx,13);

                ir = 1j*Xm(g.ind.sat_idx) ...
                     .*(idm(g.ind.sat_idx) + 1j*iqm(g.ind.sat_idx)) ...
                     ./(rr(g.ind.sat_idx) + 1j*g.ind.ind_pot(g.ind.sat_idx,4));

                gs = dessat(ism,isat);
                gr = dessat(abs(ir),isat);

                xs(g.ind.sat_idx) = xs(g.ind.sat_idx).*(1 + gs)/2;
                xr(g.ind.sat_idx) = xr(g.ind.sat_idx).*(1 + gr)/2;

                g.ind.ind_pot(g.ind.sat_idx,3) = Xm(g.ind.sat_idx) ...
                                                 + xs(g.ind.sat_idx);

                g.ind.ind_pot(g.ind.sat_idx,4) = Xm(g.ind.sat_idx) ...
                                                 + xr(g.ind.sat_idx);

                g.ind.ind_pot(g.ind.sat_idx,5) = ...
                    xs(g.ind.sat_idx) ...
                    + Xm(g.ind.sat_idx).*xr(g.ind.sat_idx) ...
                      ./g.ind.ind_pot(g.ind.sat_idx,4);

                g.ind.ind_pot(g.ind.sat_idx,6) = g.ind.ind_pot(g.ind.sat_idx,3) ...
                                                 - g.ind.ind_pot(g.ind.sat_idx,5);

                g.ind.ind_pot(g.ind.sat_idx,7) = g.sys.basrad*rr(g.ind.sat_idx) ...
                                                 ./g.ind.ind_pot(g.ind.sat_idx,4);
            end

            if ~isempty(g.ind.dbc_idx)
                % reset double cage
                [rdc,xdc] = dbcage(rr(g.ind.dbc_idx),xr(g.ind.dbc_idx), ...
                                   rr2(g.ind.dbc_idx),xr2(g.ind.dbc_idx), ...
                                   g.ind.slip(g.ind.dbc_idx,k));

                g.ind.ind_pot(g.ind.dbc_idx,4) = Xm(g.ind.dbc_idx) + xdc;

                g.ind.ind_pot(g.ind.dbc_idx,5) = ...
                    xs(g.ind.dbc_idx) + Xm(g.ind.dbc_idx).*xdc ...
                                        ./g.ind.ind_pot(g.ind.dbc_idx,4);

                g.ind.ind_pot(g.ind.dbc_idx,6) = g.ind.ind_pot(g.ind.dbc_idx,3) ...
                                                 - g.ind.ind_pot(g.ind.dbc_idx,5);

                g.ind.ind_pot(g.ind.dbc_idx,7) = g.sys.basrad*rdc ...
                                                 ./g.ind.ind_pot(g.ind.dbc_idx,4);
            end

            if ~isempty(g.ind.db_idx)
                % reset deepbar
                [rdb,xdb] = deepbar(rr(g.ind.db_idx),dbf(g.ind.db_idx), ...
                                    g.ind.slip(g.ind.db_idx,k));

                g.ind.ind_pot(g.ind.db_idx,4) = Xm(g.ind.db_idx) ...
                                                + xr(g.ind.db_idx) + xdb;

                g.ind.ind_pot(g.ind.db_idx,5) = ...
                    xs(g.ind.db_idx) ...
                    + Xm(g.ind.db_idx).*(xr(g.ind.db_idx) + xdb) ...
                      ./g.ind.ind_pot(g.ind.db_idx,4);

                g.ind.ind_pot(g.ind.db_idx,6) = g.ind.ind_pot(g.ind.db_idx,3) ...
                                                - g.ind.ind_pot(g.ind.db_idx,5);

                g.ind.ind_pot(g.ind.db_idx,7) = g.sys.basrad*rdb ...
                                                ./g.ind.ind_pot(g.ind.db_idx,4);
            end

            % Brereton, Lewis and Young motor model
            g.ind.dvdp(:,k) = -(iqm.*g.ind.ind_pot(:,6) ...
                                + g.ind.vdp(:,k)).*g.ind.ind_pot(:,7) ...
                              + g.ind.vqp(:,k).*g.ind.slip(:,k)*g.sys.basrad;

            g.ind.dvqp(:,k) = (idm.*g.ind.ind_pot(:,6) ...
                               - g.ind.vqp(:,k)).*g.ind.ind_pot(:,7) ...
                              - g.ind.vdp(:,k).*g.ind.slip(:,k)*g.sys.basrad;

            g.ind.t_mot(:,k) = g.ind.vdp(:,k).*idm + g.ind.vqp(:,k).*iqm;

            g.ind.dslip(:,k) = (g.ind.tload(:,k) - g.ind.t_mot(:,k))/2 ...
                               ./g.ind.ind_con(:,9);
        else
            error('mac_ind: dynamics calculation must be vectorized.');
        end
    end

    if (flag == 3)
        % linearize
        % add code later
    end
end

end  % function end

% eof
