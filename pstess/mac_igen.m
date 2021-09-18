function [bus_new] = mac_igen(i,k,bus,flag)
% Syntax: [bus_new] = mac_igen(i,k,bus,flag)
%
% Purpose: Simple induction generator model (single cage)
%
% Note:    No leakage inductance saturation. Induction generator pick up
%          power from negative load.
%
% Input:   i is the induction generator number, 0 for vectorized computation
%          k is the time step
%          bus is the solved load flow bus data
%          flag is 0 for initialization
%                  1 for network interface
%                  2 for dertermination of rates of change of states
%                  3 for formation of linearized state matrix
%
% Output:  bus_new is a modified bus matrix with the induction generator
%          active and reactive powers subtracted from the original loads
%          at the generator bus
%
% Data format igen_con
%          1 - induction generator number
%          2 - busnumber
%          3 - base MVA
%          4 - rs
%          5 - xs -stator leakage reactance
%          6 - Xm - magnetizing reactance
%          7 - rr
%          8 - xr - rotor leakage reactance
%          9 - H  - inertia constant generator + turbine in sec
%          15 - fraction of bus load power taken by generator

%-----------------------------------------------------------------------------%
% Version history
%
% Purpose: Induction Generator Model
% Version: 1.0
% Author: Graham Rogers
% Date July 1997
%-----------------------------------------------------------------------------%

global g;  % declaring struct of global variables

bus_new = bus;
if ~isempty(g.igen.igen_con)
    if (flag == 0)
        % initialization
        if (i == 0)
            % vector computation

            g.igen.igbus = g.bus.bus_int(g.igen.igen_con(:,2));
            g.igen.igen_pot = zeros(g.igen.n_ig,7);

            % igen_pot(,1) -- scaled mva base
            % igen_pot(,2) -- base kv
            g.igen.igen_pot(:,1) = g.sys.basmva./g.igen.igen_con(:,3);
            g.igen.igen_pot(:,2) = ones(g.igen.n_ig,1);

            % ig_vm -- ind gen terminal voltage mag
            % ig_ang -- ind gen term voltage angle
            ig_vm(:,1) = bus(g.igen.igbus,2);
            ig_ang(:,1) = bus(g.igen.igbus,3)*pi/180;

            % complex voltage
            v = ig_vm(:,1).*exp(1j*ig_ang(:,1));

            g.igen.vdig(:,1) = real(v);
            g.igen.vqig(:,1) = imag(v);

            % pig -- ind generator power
            g.igen.pig(:,1) = bus(g.igen.igbus,6).*g.igen.igen_con(:,15);

            % modify bus load power
            bus_new(g.igen.igbus,6) = bus(g.igen.igbus,6) - g.igen.pig(:,1);

            % igen_pot(,3) -- Xs
            % igen_pot(,4) -- Xr
            g.igen.igen_pot(:,3) = g.igen.igen_con(:,5) + g.igen.igen_con(:,6);
            g.igen.igen_pot(:,4) = g.igen.igen_con(:,8) + g.igen.igen_con(:,6);

            % igen_pot(,5) -- Xsp
            % igen_pot(,6) -- Xs-Xsp
            % igen_pot(,7) -- 1/Tr
            g.igen.igen_pot(:,5) = ...
                g.igen.igen_con(:,5) ...
                + g.igen.igen_con(:,6).*g.igen.igen_con(:,8)./g.igen.igen_pot(:,4);

            g.igen.igen_pot(:,6) = g.igen.igen_pot(:,3) - g.igen.igen_pot(:,5);

            g.igen.igen_pot(:,7) = ...
                g.sys.basrad*g.igen.igen_con(:,7)./g.igen.igen_pot(:,4);

            rs = g.igen.igen_con(:,4);
            xs = g.igen.igen_con(:,5);
            Xm = g.igen.igen_con(:,6);
            rr = g.igen.igen_con(:,7);
            xr = g.igen.igen_con(:,8);

            % find initial slip
            slip_old = zeros(g.igen.n_ig,1);
            slip_new = ones(g.igen.n_ig,1);

            % Newton-Raphson iteration to determine initial slip
            iter = 0;
            itermax = 50;
            err = max(abs(slip_new - slip_old));
            while ((err >= 1e-8) && (iter <= itermax))
                iter = iter + 1;
                y = g.sys.basrad.*slip_old./g.igen.igen_pot(:,7);

                denom = ones(g.igen.n_ig,1) + y.*y;
                zr = rs + y.*g.igen.igen_pot(:,6)./denom;
                zi = g.igen.igen_pot(:,5) + g.igen.igen_pot(:,6)./denom;

                dzr = g.igen.igen_pot(:,6) ...
                      .*(ones(g.igen.n_ig,1) - y.*y)./denom./denom;

                dzi = -2*g.igen.igen_pot(:,6).*y./denom./denom;

                zmod2 = zr.*zr + zi.*zi;
                dp = v.*conj(v).*(dzr.*zmod2 - 2*zr.*(dzr.*zr + dzi.*zi));
                dp = dp./zmod2./zmod2;

                peig = v.*conj(v).*zr./zmod2;
                ynew = y - (peig - g.igen.pig(:,1).*g.igen.igen_pot(:,1))./dp;

                % define the mismatch and update the slip
                slip_new = ynew.*g.igen.igen_pot(:,7)/g.sys.basrad;
                err = max(abs(slip_new - slip_old));

                slip_old = slip_new;
            end

            if (iter > itermax)
                error('mac_igen: generator slip calculation failed to converge.');
            end

            g.igen.slig(:,1) = slip_new;
            y = g.sys.basrad*g.igen.slig(:,1)./g.igen.igen_pot(:,7);

            denom = ones(g.igen.n_ig,1) + y.*y;
            zr = rs + y.*g.igen.igen_pot(:,6)./denom;
            zi = g.igen.igen_pot(:,5) + g.igen.igen_pot(:,6)./denom;

            iig = v./(zr + 1j*zi);
            sig = v.*conj(iig);

            peig = real(sig);
            qeig = imag(sig);

            % complex initial rotor states
            vp = v - (rs + 1j*g.igen.igen_pot(:,5)).*iig;
            g.igen.vdpig = real(vp);
            g.igen.vqpig = imag(vp);

            % initial prime mover torque
            g.igen.tmig(:,1) = real(vp.*conj(iig));
            g.igen.idig(:,1) = real(iig)./g.igen.igen_pot(:,1);
            g.igen.iqig(:,1) = imag(iig)./g.igen.igen_pot(:,1);

            % modify qload
            g.igen.qig(:,1) = qeig./g.igen.igen_pot(:,1);
            bus_new(g.igen.igbus,7) = bus(g.igen.igbus,7) - g.igen.qig(:,1);
        else
            % generator by generator initialization
            error('mac_igen: initialization must be vectorized.');
        end
    end

    if (flag == 1)
        % network interface
        % no interface required for induction generators
    end

    if (flag == 2)
        % induction generator dynamics calculation
        if (i == 0)
            % vector calculation

            % convert to machine base
            idigm = g.igen.idig(:,k).*g.igen.igen_pot(:,1);
            iqigm = g.igen.iqig(:,k).*g.igen.igen_pot(:,1);

            % Brereton, Lewis and Young motor model
            g.igen.dvdpig(:,k) = -(iqigm.*g.igen.igen_pot(:,6) ...
                                   + g.igen.vdpig(:,k)).*g.igen.igen_pot(:,7) ...
                                 + g.igen.vqpig(:,k).*g.igen.slig(:,k)*g.sys.basrad;

            g.igen.dvqpig(:,k) = (idigm.*g.igen.igen_pot(:,6) ...
                                  - g.igen.vqpig(:,k)).*g.igen.igen_pot(:,7) ...
                                 - g.igen.vdpig(:,k).*g.igen.slig(:,k)*g.sys.basrad;

            g.igen.dslig(:,k) = (g.igen.tmig(:,k) - g.igen.vdpig(:,k).*idigm ...
                                 - g.igen.vqpig(:,k).*iqigm)/2./g.igen.igen_con(:,9);
        else
            error('mac_igen: dynamics calculation must be vectorized.');
        end
    end

    if (flag == 3)
        % linearize
        % add code later
    end
end

end  % function end

% eof
