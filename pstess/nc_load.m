function V_nc = nc_load(bus,flag,Y22,Y21,psi,V_o,tol,k,kdc)
% Syntax: V_nc = nc_load(bus,flag,Y22,Y21,psi,V_o,tol,k,kdc)
%
% Purpose: non-conforming load model, for constant power
%          and constant current loads; the nonlinear
%          equations are solved using a Newton's algorithm
%
% Input:   bus  - solved loadflow bus data
%          flag - 0 - initialization
%                 1 - network interface computation
%                 2 - computation not needed
%          Y22  - reduced self Y matrix of non-conforming loads
%          Y21  - reduced mutual Y matrix of generator internal
%                 node to non-conforming loads
%          psi  - machine internal node voltage,
%                 not used in initialization
%          V_o  - initial non-conforming load bus voltage
%                 values (vector), not used in initialization
%          tol  - tolerance for Newton's algorithm convergence,
%                 not used in initialization
%          k    - integer time (only for svc/facts models),
%                 not used in initialization
%          kdc  - dc integer time
%
% Output:  V_nc - solved non-conforming load bus voltage
%                 values (vector)
%
% See Also: red_Ybus, svc
%
% Calls: dc_load, dc_cur
%
% Called By: svm_mgen, s_simu

%-----------------------------------------------------------------------------%
% Version history
%
% Version   2.4
% Author:   Ryan T. Elliott
% Date:     September 2019
% Modified: Added functionality related to energy storage
%
% Version:  2.3
% Date:     Feb 2015
% Author:   D. Trudnowski
% Purpose:  Add power modulation
%
% Version:  2.2
% Date:     August 1997
% Author:   Graham Rogers
% Purpose:  Add load modulation
%
% Version:  2.1
% Date:     March 1997
% Author:   Graham Rogers
% Purpose:  Add support for hvdc models
%
% Version:  2.0
% Date:     November 1996
% Author:   Graham Rogers
% Purpose:  Improve convergence
% Modified: Reduction to constant impedance load for low voltages
%           modification of Jacobian
%
% Version:  1.0
% Author:   Joe H. Chow
% Date:     April 1991
%-----------------------------------------------------------------------------%

global g;  % declaring struct of global variables

if ~isempty(g.ncl.load_con)
    if (flag == 0)  % initialization
        % set up constant power and current load components in load_pot
        % load_pot(:,1) -- constant power component
        % load_pot(:,2) -- constant current component
        % load_pot(:,3) -- constant impedance equiv. of load_pot(:,1)
        % load_pot(:,4) -- constant impedance equiv. of load_pot(:,2)

        j = g.bus.bus_int(g.ncl.load_con(:,1));

        % no need for special treatment of dc buses on initialization
        V_nc = bus(j,2).*exp(1j*bus(j,3)*pi/180);

        % complex power from constant power injections
        g.ncl.load_pot(:,1) = bus(j,6).*g.ncl.load_con(:,2) ...
                              + 1j*bus(j,7).*g.ncl.load_con(:,3);

        % complex power from constant current injections
        S_cc = bus(j,6).*g.ncl.load_con(:,4) + 1j*bus(j,7).*g.ncl.load_con(:,5);

        % complex current from constant current injections
        g.ncl.load_pot(:,2) = S_cc./abs(V_nc);

        % conj. of impedance for constant power injections
        g.ncl.load_pot(:,3) = g.ncl.load_pot(:,1)./V_nc./conj(V_nc);

        % conj. of impedance for constant current injections
        g.ncl.load_pot(:,4) = S_cc./V_nc./conj(V_nc);
    end

    if (flag == 1)  % network interface computation
        if (nargin < 7)
            tol = 1e-6;
        end

        if (nargin < 6)
            V_o = ones(g.ncl.n_load,1);                       % set flat start
        end

        V_nc = V_o;                                           % set initial voltages

        if (g.svc.n_svc ~= 0)
            j = g.svc.svc_idx;
            % Y22 is a file local variable
            Y22(j,j) = Y22(j,j) + 1j*diag(g.svc.B_cv(:,k));
        end

        if (g.tcsc.n_tcsc ~= 0)
            j = g.tcsc.tcscf_idx;
            jj = g.tcsc.tcsct_idx;
            Y22(j,j) = Y22(j,j) + 1j*diag(g.tcsc.B_tcsc(:,k));
            Y22(j,jj) = Y22(j,jj) - 1j*diag(g.tcsc.B_tcsc(:,k));
            Y22(jj,j) = Y22(jj,j) - 1j*diag(g.tcsc.B_tcsc(:,k));
            Y22(jj,jj) = Y22(jj,jj) + 1j*diag(g.tcsc.B_tcsc(:,k));
        end

        if (g.lmod.n_lmod ~= 0)
            j = g.lmod.lmod_idx;
            Y22(j,j) = Y22(j,j) + diag(g.lmod.lmod_st(:,k));
        end

        if (g.rlmod.n_rlmod ~= 0)
            j = g.rlmod.rlmod_idx;
            Y22(j,j) = Y22(j,j) + 1j*diag(g.rlmod.rlmod_st(:,k));
        end

        if ((g.trip.enable) && (g.trip.n_trip_ncl ~= 0))
            j = g.trip.bus_lshed_idx(g.trip.bus_lshed_type);
            Y22(j,j) = Y22(j,j) + diag(g.trip.load_trip_ncl(:,k));  % complex
        end

        % pwrmod (modifies load_pot to account for modulation)
        if (g.pwr.n_pwrmod ~= 0)
            for index = 1:g.pwr.n_pwrmod
                if ((g.ncl.load_con(g.pwr.pwrmod_idx(index),2) == 1) ...
                    && (g.ncl.load_con(g.pwr.pwrmod_idx(index),3) == 1))
                    % power modulation
                    g.ncl.load_pot(g.pwr.pwrmod_idx(index),1) = ...
                        -(g.pwr.pwrmod_p_st(index,k) ...
                          + 1j*g.pwr.pwrmod_q_st(index,k));

                    % power modulation with v<0.5 - maybe disable?
                    % g.ncl.load_pot(g.pwr.pwrmod_idx(index),3) = ...
                    %     -(g.pwr.pwrmod_p_st(index,k) ...
                    %       + 1j*g.pwr.pwrmod_q_st(index,k)) ...
                    %      ./(V_nc(index)*conj(V_nc(index)));
                elseif ((g.ncl.load_con(g.pwr.pwrmod_idx(index),4) == 1) ...
                        && (g.ncl.load_con(g.pwr.pwrmod_idx(index),5) == 1))
                    % current modulation
                    g.ncl.load_pot(g.pwr.pwrmod_idx(index),2) = ...
                        -(g.pwr.pwrmod_p_st(index,k) ...
                          + 1j*g.pwr.pwrmod_q_st(index,k));

                    % current modulation with v<0.5 - maybe disable?
                    % g.ncl.load_pot(g.pwr.pwrmod_idx(index),4) = ...
                    %     -(g.pwr.pwrmod_p_st(index,k) ...
                    %       + 1j*g.pwr.pwrmod_q_st(index,k))./abs(V_nc(index));
                end
            end
        end

        hv_idx = find(abs(V_nc) > 0.5);    % non-faulted buses
        lv_idx = find(abs(V_nc) <= 0.5);   % faulted buses

        curr_nc = zeros(length(V_nc),1);   % initialization
        curr_mis = zeros(g.ncl.n_load,1);
        curr_load = Y21*psi + Y22*V_nc;

        if ~isempty(hv_idx)
            % S_nc_cp -- complex power from constant P
            % S_nc_cc -- complex power from constant I
            S_nc_cp = g.ncl.load_pot(hv_idx,1);
            S_nc_cc = g.ncl.load_pot(hv_idx,2).*abs(V_nc(hv_idx));
            curr_nc(hv_idx) = -conj((S_nc_cp + S_nc_cc)./V_nc(hv_idx));
        end

        % converts to constant impedance if abs(V_nc) <= 0.5
        if ~isempty(lv_idx)
            % Y_nc_conj -- conjugate of admittance matrix from constant Z
            Y_nc_conj = diag(g.ncl.load_pot(lv_idx,3) + g.ncl.load_pot(lv_idx,4));
            curr_nc(lv_idx) = -conj(Y_nc_conj)*V_nc(lv_idx);
        end

        if ~isempty(g.dc.ldc_idx)
            % modify nc-current to take account of ac current (for hvdc)
            i_ac = dc_cur(V_nc(g.dc.ldc_idx),kdc);
            curr_nc(g.dc.ldc_idx) = curr_nc(g.dc.ldc_idx) - i_ac;
        end

        if ~isempty(g.ess.ess_idx)
            % accounting for ess
            ess(0,k,bus,4,V_nc(g.ess.ess_idx));  % function call
            curr_nc(g.ess.ess_idx) = curr_nc(g.ess.ess_idx) + g.ess.ess_cur(:,k);
        end

        curr_mis = curr_nc - curr_load;    % current mismatch

        iter = 0;
        itermax = 50;
        Y22_real = real(Y22);
        Y22_imag = imag(Y22);

        % Newton's algorithm
        while(norm(curr_mis,'inf') > tol)
            if ~isempty(hv_idx)
                % form Jacobian
                v_re = real(V_nc(hv_idx));
                v_im = imag(V_nc(hv_idx));
                v_mag = abs(V_nc(hv_idx));
                v_mag2 = v_mag.*v_mag;
                v1 = v_re.*v_re - v_im.*v_im;
                v2 = v_re.*v_im;
                vp11 = -v1.*real(g.ncl.load_pot(hv_idx,1)) ...
                       - 2*v2.*imag(g.ncl.load_pot(hv_idx,1));
                vp12 = v1.*imag(g.ncl.load_pot(hv_idx,1)) ...
                       - 2*v2.*real(g.ncl.load_pot(hv_idx,1));
                vi11 = v_im.*v_im.*real(g.ncl.load_pot(hv_idx,2)) ...
                       - v2.*imag(g.ncl.load_pot(hv_idx,2));
                vi12 = v_re.*v_re.*imag(g.ncl.load_pot(hv_idx,2)) ...
                       - v2.*real(g.ncl.load_pot(hv_idx,2));
                vi21 = -v_im.*v_im.*imag(g.ncl.load_pot(hv_idx,2)) ...
                       - v2.*real(g.ncl.load_pot(hv_idx,2));
                vi22 = v_re.*v_re.*real(g.ncl.load_pot(hv_idx,2)) ...
                       + v2.*imag(g.ncl.load_pot(hv_idx,2));
                vp11 = vp11./v_mag2./v_mag2;
                vp12 = vp12./v_mag2./v_mag2;
                vi11 = vi11./v_mag./v_mag2;
                vi12 = vi12./v_mag./v_mag2;
                vi21 = vi21./v_mag./v_mag2;
                vi22 = vi22./v_mag./v_mag2;
                v_s11(hv_idx) = vi11 + vp11;
                v_s12(hv_idx) = vp12 + vi12;
                v_s21(hv_idx) = vp12 + vi21;
                v_s22(hv_idx) = -vp11 + vi22;
            end

            if ~isempty(lv_idx)
                v_s11(lv_idx) = -real(g.ncl.load_pot(lv_idx,3) ...
                                      + g.ncl.load_pot(lv_idx,4));
                v_s12(lv_idx) = imag(g.ncl.load_pot(lv_idx,3) ...
                                     + g.ncl.load_pot(lv_idx,4));
                v_s22(lv_idx) = v_s11(lv_idx);
                v_s21(lv_idx) = -v_s12(lv_idx);
            end

            % modify Jacobian for dc
            if (g.dc.n_conv ~= 0)
                V = V_o(g.dc.ldc_idx);
                [y1,y2,y3,y4] = dc_load(V,kdc);
                v_s11(g.dc.ldc_idx) = v_s11(g.dc.ldc_idx) + y1';
                v_s12(g.dc.ldc_idx) = v_s12(g.dc.ldc_idx) + y2';
                v_s21(g.dc.ldc_idx) = v_s21(g.dc.ldc_idx) + y3';
                v_s22(g.dc.ldc_idx) = v_s22(g.dc.ldc_idx) + y4';
            end

            % update voltage
            iter = iter + 1;

            Jac_nc = [Y22_real + diag(v_s11), diag(v_s12)-Y22_imag;
                      Y22_imag + diag(v_s21), Y22_real + diag(v_s22)];

            b = [real(curr_mis); imag(curr_mis)];
            x = Jac_nc\b;  % solve for voltage increment
            V_nc = V_nc + x(1:g.ncl.n_load,1) ...
                   + 1j*x(g.ncl.n_load + 1:2*g.ncl.n_load,1);

            hv_idx = find(abs(V_nc) > 0.5);
            lv_idx = find(abs(V_nc) <= 0.5);

            % compute non-conforming load current
            curr_load = Y21*psi + Y22*V_nc;
            if ~isempty(hv_idx)
                % S_nc_cp -- complex power from constant P
                % S_nc_cc -- complex power from constant I
                S_nc_cp = g.ncl.load_pot(hv_idx,1);
                S_nc_cc = g.ncl.load_pot(hv_idx,2).*abs(V_nc(hv_idx));
                curr_nc(hv_idx) = -conj((S_nc_cp + S_nc_cc)./V_nc(hv_idx));
            end

            if ~isempty(lv_idx)
                % Y_nc_conj -- conjugate of admittance matrix from constant Z
                Y_nc_conj = diag(g.ncl.load_pot(lv_idx,3) ...
                                 + g.ncl.load_pot(lv_idx,4));
                curr_nc(lv_idx) = -conj(Y_nc_conj)*V_nc(lv_idx);
            end

            if ~isempty(g.dc.ldc_idx)
                % modify nc-current to take account of ac current (for hvdc)
                i_ac = dc_cur(V_nc(g.dc.ldc_idx),kdc);
                curr_nc(g.dc.ldc_idx) = curr_nc(g.dc.ldc_idx) - i_ac;
            end

            if ~isempty(g.ess.ess_idx)
                % accounting for ess
                ess(0,k,bus,4,V_nc(g.ess.ess_idx));  % function call
                curr_nc(g.ess.ess_idx) = curr_nc(g.ess.ess_idx) + g.ess.ess_cur(:,k);
            end

            curr_mis = curr_nc - curr_load;    % current mismatch

            if (iter > itermax)
                estr = '\nnc_load: Newton algorithm ';
                estr = [estr, 'has not converged in %g iterations.'];
                estr = [estr, '\n         Current mismatch = %g with tol = %g.'];
                error(sprintf(estr,[itermax,norm(curr_mis,'inf'),tol]));
            end
        end
    end

    if (flag == 2)  % no dynamics calculation needed
    end
end

end  % function end

% eof
