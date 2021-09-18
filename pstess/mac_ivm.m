function mac_ivm(i,k,bus,flag)
% Syntax: mac_ivm(i,k,bus,flag)
%
% Purpose: Internal Voltage Model type generator
%
% Input: i - generator number
%          - 0 for vectorized computation
%        k - integer time
%        bus - solved loadflow bus data
%        flag - 0 - initialization
%               1 - network interface computation
%               2 - generator dynamics computation and state
%                   state matrix building

%-----------------------------------------------------------------------------%
% Version history
%
% Version: 1.0
% Date:    June 2019
% Author:  D. Trudnowski
%-----------------------------------------------------------------------------%

global g;  % declaring struct of global variables

if (g.mac.n_ivm ~= 0)
    if (flag == 0)   % initialization
        if (i ~= 0)
            % non-vectorized computation
            error('mac_ivm: initialization must be vectorized.');
        else
            % vectorized computation

            % busnum -- vector of bus numbers where ivm's are connected
            busnum = g.bus.bus_int(g.mac.mac_con(g.mac.mac_ivm_idx,2));

            % mac_pot(,1) -- factor to convert between MVA bases
            % mac_pot(,2) -- base kv?
            % eterm -- terminal bus voltage magnitude
            % theta -- terminal bus angle in radians
            g.mac.mac_pot(g.mac.mac_ivm_idx,1) = ...
                g.sys.basmva*ones(g.mac.n_ivm,1)./g.mac.mac_con(g.mac.mac_ivm_idx,3);
            g.mac.mac_pot(g.mac.mac_ivm_idx,2) = 1.0*ones(g.mac.n_ivm,1);
            g.mac.eterm(g.mac.mac_ivm_idx,1) = bus(busnum,2);
            g.bus.theta(busnum,1) = bus(busnum,3)*pi/180;

            % pelect -- electrical output power, active
            % qelect -- electrical output power, reactive
            g.mac.pelect(g.mac.mac_ivm_idx,1) = ...
                bus(busnum,4).*g.mac.mac_con(g.mac.mac_ivm_idx,22);
            g.mac.qelect(g.mac.mac_ivm_idx,1) = ...
                bus(busnum,5).*g.mac.mac_con(g.mac.mac_ivm_idx,23);

            % Terminal voltage
            Vt = g.mac.eterm(g.mac.mac_ivm_idx,1).*exp(1j*g.bus.theta(busnum,1));

            % It -- Current out of terminal
            % E -- Internal voltage
            It = conj((g.mac.pelect(g.mac.mac_ivm_idx,1) ...
                       + 1j*g.mac.qelect(g.mac.mac_ivm_idx,1))./Vt);
            E = Vt + (g.mac.mac_con(g.mac.mac_ivm_idx,5) ...
                      + 1j*g.mac.mac_con(g.mac.mac_ivm_idx,7))*It;

            % edprime -- internal voltage magnitude
            % mac_ang -- internal voltage angle (rad.)
            g.mac.edprime(g.mac.mac_ivm_idx,1) = abs(E);
            g.mac.mac_ang(g.mac.mac_ivm_idx,1) = atan2(imag(E),real(E));

            % psi_re -- real-part of internal voltage used for network interface
            % psi_im -- imag-part of internal voltage used for network interface
            g.mac.psi_re(g.mac.mac_ivm_idx,1) = real(E);
            g.mac.psi_im(g.mac.mac_ivm_idx,1) = imag(E);

            % g.mac.pmech(g.mac.mac_ivm_idx,1) = ...
            %     g.mac.pelect(g.mac.mac_ivm_idx,1) ...
            %     + sign(g.mac.pelect(g.mac.mac_ivm_idx,1)) ...
            %       .*(abs(It).^2).*g.mac.mac_con(g.mac.mac_ivm_idx,5);
        end
    end

    if (flag == 1)  % network interface computation
        if (i ~= 0)
            % non-vectorized computation
            error('mac_ivm: network interface calculation must be vectorized.');
        else
            % vectorized computation

            % mac_ang -- angle with respect to machine reference
            g.mac.mac_ang(g.mac.mac_ivm_idx,k) = ...
                g.mac.mac_ang(g.mac.mac_ivm_idx,k) ...
                - g.sys.mach_ref(k)*ones(g.mac.n_ivm,1);
            E = g.mac.edprime(g.mac.mac_ivm_idx,k) ...
                .*exp(1j*g.mac.mac_ang(g.mac.mac_ivm_idx,k));

            % psi_re -- real-part of internal voltage used for network interface
            % psi_im -- imag-part of internal voltage used for network interface
            g.mac.psi_re(g.mac.mac_ivm_idx,k) = real(E);
            g.mac.psi_im(g.mac.mac_ivm_idx,k) = imag(E);
        end
    end

    if (flag == 2)  % generator dynamics calculation
        if (i ~= 0)
            % non-vectorized computation
            error('mac_ivm: dynamics calculation must be vectorized.');
        else
            % vectorized computation
            g.mac.dmac_ang(g.mac.mac_ivm_idx,k) = ...
                (-g.mac.mac_ang(g.mac.mac_ivm_idx,k) ...
                 + g.mac.ivmmod_d_sig(:,k))./g.mac.mac_con(g.mac.mac_ivm_idx,9);
            g.mac.dedprime(g.mac.mac_ivm_idx,k) = ...
                (-g.mac.edprime(g.mac.mac_ivm_idx,k) ...
                 + g.mac.ivmmod_e_sig(:,k))./g.mac.mac_con(g.mac.mac_ivm_idx,10);

            % busnum -- vector of bus numbers where ivm's are connected
            busnum = g.bus.bus_int(g.mac.mac_con(g.mac.mac_ivm_idx,2));
            Vt = g.bus.bus_v(busnum,k);
            E = g.mac.edprime(g.mac.mac_ivm_idx,k) ...
                .*exp(1j*g.mac.mac_ang(g.mac.mac_ivm_idx,k));
            It = (E - Vt)./(g.mac.mac_con(g.mac.mac_ivm_idx,5) ...
                            + 1j*g.mac.mac_con(g.mac.mac_ivm_idx,7));

            S = Vt.*conj(It);
            g.mac.pelect(g.mac.mac_ivm_idx,k) = real(S);
            g.mac.qelect(g.mac.mac_ivm_idx,k) = imag(S);
            g.mac.eterm(g.mac.mac_ivm_idx,k) = abs(Vt);
        end
    end
end

end  % function end

% eof
