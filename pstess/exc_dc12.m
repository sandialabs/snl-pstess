function exc_dc12(i,k,flag)
% Syntax: exc_dc12(i,k,flag)
%
% Purpose: excitation system, model DC1 (exc_con(i,1)=1)
%          and model DC2 (exc_con(i,1)=2)
%          with vectorized computation option
%
% Input:   i - generator number
%                 0 - vectorized computation
%          k - integer time
%          flag - 0 - initialization
%                 1 - network interface computation
%                 2 - generator dynamics computation
%
% See Also: exc_st3, smp_exc

%-----------------------------------------------------------------------------%
% Version history
%
% Purpose:  Corrected error on satuation block
% Date:     12/17/2015
% Author:   Joe chow
%
% Purpose:  Corrected error with no_TE
% Date:     2/7/98
% Author:   Graham Rogers
%
% Version:  2.2
% Date:     30/6/98
% Purpose:  Corrected satuaration Asat modified
%
% Version:  2.1
% Date:     14/8/97
% Author:   Graham Rogers
% Purpose:  Reverted to exciter number for exc_sig
%           Added pss_out so that exc_sig can be used for other
%           control inputs if desired
%           Corrected sign of exc_sig in error signal
%
% Version:  2.0
% Date:     22/6/96
% Author:   Graham Rogers
% Purpose:  To enable vector calculations with different
%           exciter models on generators
%           This means that the vector option (i=0)
%           is the mode normally used
%           Indexes formed in exc_indx are required
% Modified: Changed way exc_sig is referenced, i.e.,
%           in terms of the generator number and not the
%           exciter number
%
% Version:  1.1
% Date:     22/5/95
% Author:   GJR
% Purpose:  Correct errors
% Modified: Saturation model and Lead/Lag initialization
%
% Version:  1.0
% Author:   Joe H. Chow
% Date:     March 1991
%-----------------------------------------------------------------------------%

global g;  % declaring struct of global variables

[nexc,~] = size(g.exc.exc_con);

if (flag == 0)   % initialization
    if (i ~= 0)  % scalar computation
        n = g.mac.mac_int(g.exc.exc_con(i,2));  % machine number
        g.exc.Efd(i,1) = g.mac.vex(n,1);

        if (g.exc.exc_con(i,1) ~= 1 || g.exc.exc_con(i,1) ~= 2)
            estr = '\nexc_dc12: the specified exciter type does not match the ';
            estr = [estr, 'dynamic model at exciter index %0.0f.'];
            error(sprintf(estr,i));
        end

        if (g.exc.exc_con(i,11) == 0)           % T_E = 0
            g.exc.V_R(i,1) = g.mac.vex(n,1);
            g.exc.exc_pot(i,1) = 0;
            g.exc.exc_pot(i,2) = 0;
        else
            if (g.exc.exc_con(i,14) == 0)       % Asat and Bsat specified
                g.exc.exc_pot(i,1) = g.exc.exc_con(i,13);
                g.exc.exc_pot(i,2) = g.exc.exc_con(i,15);
            else
                % exc_pot(,1) -- A
                % exc_pot(,2) -- B
                g.exc.exc_pot(i,2) = ...
                    log(g.exc.exc_con(i,15)*g.exc.exc_con(i,14) ...
                        /g.exc.exc_con(i,12)/g.exc.exc_con(i,13)) ...
                    /(g.exc.exc_con(i,14) - g.exc.exc_con(i,12));

                g.exc.exc_pot(i,1) = g.exc.exc_con(i,13) ...
                                     /exp(g.exc.exc_pot(i,2)*g.exc.exc_con(i,12));
            end

            if (g.exc.exc_con(i,8) == 0)        % Vrmax = 0 and compute Vrmax
                % assuming E2 is max Efd
                g.exc.exc_con(i,8) = ...
                    g.exc.exc_con(i,10)*g.exc.exc_con(i,14) ...
                    + sign(g.exc.exc_con(i,14))*g.exc.exc_pot(i,1) ...
                      *exp(g.exc.exc_pot(i,2)*abs(g.exc.exc_con(i,14)));

                if (g.exc.exc_con(i,1) == 2)    % bus fed voltage limit
                    g.exc.exc_con(i,8) = g.exc.exc_con(i,8)/g.mac.eterm(n,1);
                end

                g.exc.exc_con(i,9) = -g.exc.exc_con(i,8);
            end
        end

        % KE = 0, calculate KE to make V_R zero
        if (g.exc.exc_con(i,10) == 0)
            g.exc.exc_con(i,10) = -g.exc.exc_pot(i,1) ...
                                   *exp(g.exc.exc_pot(i,2)*abs(g.exc.Efd(i,1))) ...
                                   /abs(g.exc.Efd(i,1));
        end

        g.exc.V_R(i,1) = g.exc.exc_con(i,10)*g.exc.Efd(i,1) ...
                         + g.exc.Efd(i,1)*g.exc.exc_pot(i,1) ...
                           *exp(g.exc.exc_pot(i,2)*(g.exc.Efd(i,1)));

        if (g.exc.exc_con(i,1) == 1)
            mult = 1;
        else
            mult = g.mac.eterm(n,1);
        end

        if (g.exc.V_R(i,1) > g.exc.exc_con(i,8)*mult)
            estr = '\nexc_dc12: V_R exceeds upper limit at machine index %0.0f.';
            error(sprintf(estr,n));
        elseif (g.exc.V_R(i,1) < g.exc.exc_con(i,9)*mult)
            estr = '\nexc_dc12: V_R under lower limit at machine index %0.0f.';
            error(sprintf(estr,n));
        end

        if (g.exc.exc_con(i,4) == 0)  % KA = 0
            g.exc.exc_con(i,4) == 1;  % reset to 1
        end

        g.exc.V_A(i,1) = g.exc.V_R(i,1)/g.exc.exc_con(i,4);  % leadlag
        g.exc.V_As(i,1) = g.exc.V_A(i,1);                    % leadlag state variable
        if (g.exc.exc_con(i,6) ~= 0)
            g.exc.exc_pot(i,4) = g.exc.exc_con(i,7)/g.exc.exc_con(i,6);
        end

        % err         -- summing junction error
        % R_f         -- rate feedback state
        % V_FB        -- feedback from stabilizing transformer
        % V_TR        -- input filter state
        % exc_pot(,3) -- reference voltage
        % exc_pot(,5) -- stabilizer gain conversion factor

        err = g.exc.V_A(i,1);

        g.exc.R_f(i,1) = g.exc.Efd(i,1);
        g.exc.V_FB(i,1) = 0;
        g.exc.V_TR(i,1) = g.mac.eterm(n,1);

        g.exc.exc_pot(i,3) = g.mac.eterm(n,1) + err;
        g.exc.exc_pot(i,5) = g.exc.exc_con(i,16)/g.exc.exc_con(i,17);

    else  % vectorized computation
        if (g.exc.n_dc ~= 0)
            % dc machine numbers
            n = g.mac.mac_int(g.exc.exc_con(g.exc.dc_idx,2));
            if (g.exc.n_dc2 ~= 0)
                % dc type 2 exciters gen numbers
                n2 = g.mac.mac_int(g.exc.exc_con(g.exc.dc2_idx,2));
            end

            g.exc.Efd(g.exc.dc_idx,1) = g.mac.vex(n,1);
            TE = g.exc.dc_TE(g.exc.dc_TE_idx);

            expsat = find(g.exc.exc_con(g.exc.dc_idx,14) == 0);
            if ~isempty(expsat)
                g.exc.exc_pot(g.exc.dc_idx(expsat),1) = ...
                    g.exc.exc_con(g.exc.dc_idx(expsat),13);

                g.exc.exc_pot(g.exc.dc_idx(expsat),2) = ...
                    g.exc.exc_con(g.exc.dc_idx(expsat),15);
            end

            sesat = find(g.exc.exc_con(g.exc.dc_idx,14) ~= 0);
            if ~isempty(sesat)
                g.exc.exc_pot(g.exc.dc_idx(sesat),2) = ...
                    log(g.exc.exc_con(g.exc.dc_idx(sesat),15) ...
                        .*g.exc.exc_con(g.exc.dc_idx(sesat),14) ...
                        ./g.exc.exc_con(g.exc.dc_idx(sesat),12) ...
                        ./g.exc.exc_con(g.exc.dc_idx(sesat),13)) ...
                    ./(g.exc.exc_con(g.exc.dc_idx(sesat),14) ...
                       - g.exc.exc_con(g.exc.dc_idx(sesat),12));     % Bsat

                g.exc.exc_pot(g.exc.dc_idx(sesat),1) = ...
                    g.exc.exc_con(g.exc.dc_idx(sesat),13) ...
                    ./exp(g.exc.exc_pot(g.exc.dc_idx(sesat),2) ...
                          .*g.exc.exc_con(g.exc.dc_idx(sesat),12));  % Asat
            end

            if ~isempty(g.exc.dc_noTE_idx)
                no_TE = g.exc.dc_noTE_idx;
                g.exc.V_R(g.exc.dc_idx(no_TE),1) = g.mac.vex(n(no_TE),1);

                g.exc.exc_pot(g.exc.dc_idx(no_TE),1) = zeros(length(no_TE),1);
                g.exc.exc_pot(g.exc.dc_idx(no_TE),2) = ...
                    g.exc.exc_pot(g.exc.dc_idx(no_TE),1);
            end

            % Vrmax = 0 and compute Vrmax assuming E2 is max Efd
            no_Vrmax = find(g.exc.exc_con(g.exc.dc_idx,8) == 0);
            if ~isempty(no_Vrmax)
                g.exc.exc_con(g.exc.dc_idx(no_Vrmax),8) = ...
                    g.exc.exc_con(g.exc.dc_idx(no_Vrmax),10) ...
                    .*g.exc.exc_con(g.exc.dc_idx(no_Vrmax),14) ...
                    + g.exc.exc_pot(g.exc.dc_idx(noVrmax),1) ...
                      .*exp(g.exc.exc_pot(g.exc.dc_idx(no_Vrmax),2) ...
                            .*abs(g.exc.exc_con(g.exc.dc_idx(no_Vrmax),14))) ...
                      .*sign(g.exc.exc_con(g.exc.dc_idx(no_Vrmax),14));

                g.exc.exc_con(g.exc.dc_idx(no_Vrmax),9) = ...
                    -g.exc.exc_con(g.exc.dc_idx(no_Vrmax),8);

                bus_fed = find(g.exc.exc_con(g.exc.dc_idx(no_Vrmax),1) == 2);
                if ~isempty(bus_fed)
                    g.exc.exc_con(g.exc.dc_idx(no_Vrmax(bus_fed)),8) = ...
                        g.exc.exc_con(g.exc.dc_idx(no_Vrmax(bus_fed)),8) ...
                        ./g.mac.eterm(n(no_Vrmax(bus_fed)),1);

                    g.exc.exc_con(g.exc.dc_idx(no_Vrmax(bus_fed)),9) = ...
                        -g.exc.exc_con(g.exc.dc_idx(no_Vrmax(bus_fed)),8);
                end
            end

            % calculate KE to make initial V_R zero
            no_KE = find(g.exc.exc_con(g.exc.dc_idx,10) == 0);       % KE = 0
            if ~isempty(no_KE)
                asat = g.exc.exc_pot(g.exc.dc_idx(no_KE),1);
                bsat = g.exc.exc_pot(g.exc.dc_idx(no_KE),2);
                vf = abs(g.exc.Efd(g.exc.dc_idx(no_KE),k));
                g.exc.exc_con(g.exc.dc_idx(no_KE),10) = -asat.*exp(bsat.*vf)./vf;
            end

            if ~isempty(g.exc.dc_TE_idx)
                R_idx = g.exc.dc_idx(g.exc.dc_TE_idx);
                g.exc.V_R(R_idx,1) = ...
                    g.exc.exc_con(R_idx,10).*g.exc.Efd(R_idx,1) ...
                    + g.exc.Efd(R_idx,1).*g.exc.exc_pot(R_idx,1) ...
                      .*exp(g.exc.exc_pot(R_idx,2).*abs(g.exc.Efd(R_idx,1)));
            end

            % check limits
            mult = ones(g.exc.n_dc,1);
            if (g.exc.n_dc2 ~= 0)
                mult(find(g.exc.exc_con(g.exc.dc_idx,1) == 2)) = g.mac.eterm(n2,1);
            end

            over_lmt = find(g.exc.V_R(g.exc.dc_idx,1) ...
                            > g.exc.exc_con(g.exc.dc_idx,8).*mult);

            if ~isempty(over_lmt)
                n_error = g.mac.mac_int(g.exc.exc_con(g.exc.dc_idx(over_lmt),2));
                estr = '\nexc_dc12: V_R exceeds upper limit at machine index %0.0f.';
                error(sprintf(estr,n_error));
            end

            under_lmt = find(g.exc.V_R(g.exc.dc_idx,1) ...
                             < g.exc.exc_con(g.exc.dc_idx,9).*mult);

            if ~isempty(under_lmt)
                n_error = g.mac.mac_int(g.exc.exc_con(g.exc.dc_idx(under_lmt),2));
                estr = '\nexc_dc12: V_R under lower limit at machine index %0.0f.';
                error(sprintf(estr,n_error));
            end

            no_KA = find(g.exc.exc_con(g.exc.dc_idx,4) == 0);        % KA = 0
            if ~isempty(no_KA)
                g.exc.exc_con(g.exc.dc_idx(no_KA),4) = ones(length(no_KA),1);
            end

            % V_A  -- leadlag
            % V_As -- leadlag state variable
            g.exc.V_A(g.exc.dc_idx,1) = g.exc.V_R(g.exc.dc_idx,1) ...
                                        ./g.exc.exc_con(g.exc.dc_idx,4);

            g.exc.V_As(g.exc.dc_idx,1) = g.exc.V_A(g.exc.dc_idx,1);

            TB = g.exc.dc_TB_idx;
            if ~isempty(TB)
                g.exc.exc_pot(g.exc.dc_idx(TB),4) = ...
                    g.exc.exc_con(g.exc.dc_idx(TB),7) ...
                    ./g.exc.exc_con(g.exc.dc_idx(TB),6);
            end

            % err         -- summing junction error
            % R_f         -- rate feedback state
            % V_FB        -- feedback from stabilizing transformer
            % V_TR        -- input filter state
            % exc_pot(,3) -- reference voltage
            % exc_pot(,5) -- stabilizer gain conversion factor

            err = g.exc.V_A(g.exc.dc_idx,1);

            g.exc.R_f(g.exc.dc_idx,1) = g.exc.Efd(g.exc.dc_idx,1);
            g.exc.V_FB(g.exc.dc_idx,1) = zeros(g.exc.n_dc,1);
            g.exc.V_TR(g.exc.dc_idx,1) = g.mac.eterm(n,1);

            g.exc.exc_pot(g.exc.dc_idx,3) = g.mac.eterm(n,1) + err;
            g.exc.exc_pot(g.exc.dc_idx,5) = g.exc.exc_con(g.exc.dc_idx,16) ...
                                            ./g.exc.exc_con(g.exc.dc_idx,17);
        end
    end
end

if (flag == 1)   % network interface computation
    if (i ~= 0)  % scalar computation
        n = g.mac.mac_int(g.exc.exc_con(i,2));                 % machine number
        g.mac.vex(n,k) = g.exc.Efd(i,k);                       % set field voltage
    else         % vectorized computation
        if (g.exc.n_dc ~= 0)
            n = g.mac.mac_int(g.exc.exc_con(g.exc.dc_idx,2));  % machine number
            g.mac.vex(n,k) = g.exc.Efd(g.exc.dc_idx,k);        % set field voltage
        end
    end
end

if (flag == 2)   % exciter dynamics calculation
    if (i ~= 0)  % scalar computation
        % machine number
        n = g.mac.mac_int(g.exc.exc_con(i,2));
        if (g.exc.exc_con(i,3) == 0)  % transducer time constant = 0
            g.exc.dV_TR(i,k) = 0;
            g.exc.V_TR(i,k) = g.mac.eterm(n,k);
        else
            g.exc.dV_TR(i,k) = (-g.exc.V_TR(i,k) + g.mac.eterm(n,k)) ...
                               /g.exc.exc_con(i,3);
        end

        g.exc.V_FB(i,k) = g.exc.exc_pot(i,5)*(g.exc.Efd(i,k) - g.exc.R_f(i,k));

        err = g.exc.exc_sig(i,k) - g.exc.V_FB(i,k) ...
              + g.exc.exc_pot(i,3) - g.exc.V_TR(i,k);

        err = err + g.exc.pss_out(i,k);

        if (g.exc.exc_con(i,6) == 0)  % no leadlag
            g.exc.dV_As(i,k) = 0;
            g.exc.V_As(i,k) = err;
            g.exc.V_A(i,k) = err;
        else
            g.exc.dV_As(i,k) = (-g.exc.V_As(i,k) + err)/g.exc.exc_con(i,6);

            g.exc.V_A(i,k) = g.exc.exc_pot(i,4)*err ...
                             + (1 - g.exc.exc_pot(i,4))*g.exc.V_As(i,k);
        end

        if (g.exc.exc_con(i,1) == 1)
            mult = 1;                 % solid fed
        else
            mult = g.mac.eterm(n,k);  % bus fed
        end

        if (g.exc.exc_con(i,5) == 0)  % no TA
            g.exc.dV_R(i,k) = 0.0;
            g.exc.V_R(i,k) = g.exc.exc_con(i,4)*g.exc.V_A(i,k);
            g.exc.V_R(i,k) = max(g.exc.exc_con(i,9)*mult, ...
                                 min(g.exc.V_R(i,k),g.exc.exc_con(i,8)*mult));
        else
            g.exc.dV_R(i,k) = ...
                (-g.exc.V_R(i,k) + g.exc.exc_con(i,4)*g.exc.V_A(i,k)) ...
                /g.exc.exc_con(i,5);

            % anti-windup reset
            if (g.exc.V_R(i,k) > g.exc.exc_con(i,8)*mult)
                g.exc.V_R(i,k) = g.exc.exc_con(i,8)*mult;
                if g.exc.dV_R(i,k)>0.0
                    g.exc.dV_R(i,k) = 0.0;
                end
            end

            if (g.exc.V_R(i,k) < g.exc.exc_con(i,9)*mult)
                g.exc.V_R(i,k) = g.exc.exc_con(i,9)*mult;
                if (g.exc.dV_R(i,k) < 0)
                    g.exc.dV_R(i,k) = 0.0;
                end
            end
        end

        if (g.exc.exc_con(i,11) == 0)  % no exciter dynamics
            g.exc.dEfd(i,k) = 0.0;
            g.exc.Efd(i,k) = g.exc.V_R(i,k);
        else
            SE = sign(g.exc.Efd(i,k))*g.exc.exc_pot(i,1) ...
                 *exp(g.exc.exc_pot(i,2)*abs(g.exc.Efd(i,k)));

            g.exc.dEfd(i,k) = ...
                (g.exc.V_R(i,k) - (g.exc.exc_con(i,10) + SE)*g.exc.Efd(i,k)) ...
                /g.exc.exc_con(i,11);
        end

        g.exc.dR_f(i,k) = (-g.exc.R_f(i,k) + g.exc.Efd(i,k)) ...
                          /g.exc.exc_con(i,17);

    else  % vectorized computation
        if (g.exc.n_dc ~= 0)
            n = g.mac.mac_int(g.exc.exc_con(g.exc.dc_idx,2));  % machine number
            if (g.exc.n_dc2 ~= 0)
                n2 = g.mac.mac_int(g.exc.exc_con(g.exc.dc2_idx,2));
            end

            TR = g.exc.dc_TR_idx;
            no_TR = g.exc.dc_noTR_idx;
            if ~isempty(no_TR)
                g.exc.dV_TR(g.exc.dc_idx(no_TR),k) = zeros(length(no_TR),1);
                g.exc.V_TR(g.exc.dc_idx(no_TR),k) = g.mac.eterm(n(no_TR),k);
            end

            if ~isempty(TR)
                g.exc.dV_TR(g.exc.dc_idx(TR),k) = ...
                    (-g.exc.V_TR(g.exc.dc_idx(TR),k) + g.mac.eterm(n(TR),k)) ...
                    ./g.exc.exc_con(g.exc.dc_idx(TR),3);
            end

            g.exc.V_FB(g.exc.dc_idx,k) = ...
                g.exc.exc_pot(g.exc.dc_idx,5) ...
                .*(g.exc.Efd(g.exc.dc_idx,k) - g.exc.R_f(g.exc.dc_idx,k));

            err = g.exc.exc_sig(g.exc.dc_idx,k) - g.exc.V_FB(g.exc.dc_idx,k) ...
                  + g.exc.exc_pot(g.exc.dc_idx,3) - g.exc.V_TR(g.exc.dc_idx,k);

            err = err + g.exc.pss_out(g.exc.dc_idx,k);

            no_TB = g.exc.dc_noTB_idx;
            if ~isempty(no_TB)
                g.exc.dV_As(g.exc.dc_idx(no_TB),k) = zeros(length(no_TB),1);
                g.exc.V_As(g.exc.dc_idx(no_TB),k) = err(no_TB);
                g.exc.V_A(g.exc.dc_idx(no_TB),k) = err(no_TB);
            end

            TB = g.exc.dc_TB_idx;
            if ~isempty(TB)
                g.exc.dV_As(g.exc.dc_idx(TB),k) = ...
                    (-g.exc.V_As(g.exc.dc_idx(TB),k) + err(TB)) ...
                    ./g.exc.exc_con(g.exc.dc_idx(TB),6);

                g.exc.V_A(g.exc.dc_idx(TB),k) = ...
                    g.exc.exc_pot(g.exc.dc_idx(TB),4).*err(TB) ...
                    + (ones(length(TB),1) - g.exc.exc_pot(g.exc.dc_idx(TB),4)) ...
                      .*g.exc.V_As(g.exc.dc_idx(TB),k);
            end

            mult = ones(g.exc.n_dc,1);
            if (g.exc.n_dc2 ~= 0)
                mult(find(g.exc.exc_con(g.exc.dc_idx,1) == 2)) = g.mac.eterm(n2,k);
            end

            no_TA = g.exc.dc_noTA_idx;
            if ~isempty(no_TA)
                g.exc.dV_R(g.exc.dc_idx(no_TA),k) = zeros(length(no_TA),1);

                g.exc.V_R(g.exc.dc_idx(no_TA),k) = ...
                    g.exc.exc_con(g.exc.dc_idx(no_TA),4) ...
                    .*g.exc.V_A(g.exc.dc_idx(no_TA),k);

                g.exc.V_R(g.exc.dc_idx(no_TA),k) = ...
                    max(g.exc.exc_con(g.exc.dc_idx(no_TA),9).*mult(no_TA), ...
                        min(g.exc.V_R(g.exc.dc_idx(no_TA),k), ...
                            g.exc.exc_con(g.exc.dc_idx(no_TA),8).*mult(no_TA)));
            end

            TA = g.exc.dc_TA_idx;
            if ~isempty(TA)
                g.exc.dV_R(g.exc.dc_idx(TA),k) = ...
                    (-g.exc.V_R(g.exc.dc_idx(TA),k) ...
                     + g.exc.exc_con(g.exc.dc_idx(TA),4) ...
                       .*g.exc.V_A(g.exc.dc_idx(TA),k)) ...
                    ./g.exc.exc_con(g.exc.dc_idx(TA),5);

                % anti-windup reset
                maxlmt = find(g.exc.V_R(g.exc.dc_idx(TA),k) ...
                              > g.exc.exc_con(g.exc.dc_idx(TA),8).*mult(TA));

                if ~isempty(maxlmt)
                    g.exc.V_R(g.exc.dc_idx(TA(maxlmt)),k) = ...
                        g.exc.exc_con(g.exc.dc_idx(TA(maxlmt)),8).*mult(TA(maxlmt));

                    pos_rate = find(g.exc.dV_R(g.exc.dc_idx(TA(maxlmt)),k) > 0);

                    prl = length(pos_rate);
                    if (prl ~= 0)
                        g.exc.dV_R(g.exc.dc_idx(TA(maxlmt(pos_rate))),k) = ...
                            zeros(prl,1);
                    end
                end

                minlmt = find(g.exc.V_R(g.exc.dc_idx(TA),k) ...
                              < g.exc.exc_con(g.exc.dc_idx(TA),9).*mult(TA));

                if ~isempty(minlmt)
                    g.exc.V_R(g.exc.dc_idx(TA(minlmt)),k) = ...
                        g.exc.exc_con(g.exc.dc_idx(TA(minlmt)),9).*mult(TA(minlmt));

                    neg_rate = find(g.exc.dV_R(g.exc.dc_idx(TA(minlmt)),k) < 0);

                    nrl = length(neg_rate);
                    if (nrl ~= 0)
                        g.exc.dV_R(g.exc.dc_idx(TA(minlmt(neg_rate))),k) = ...
                            zeros(nrl,1);
                    end
                end
            end

            no_TE = g.exc.dc_noTE_idx;
            if ~isempty(no_TE)
                g.exc.dEfd(g.exc.dc_idx(no_TE),k) = zeros(length(no_TE),1);
                g.exc.Efd(g.exc.dc_idx(no_TE),k) = g.exc.V_R(g.exc.dc_idx(no_TE),k);
            end

            TE = g.exc.dc_TE_idx;
            if ~isempty(TE)
                SE = sign(g.exc.Efd(g.exc.dc_idx(TE),k)) ...
                     .*g.exc.exc_pot(g.exc.dc_idx(TE),1) ...
                     .*exp(g.exc.exc_pot(g.exc.dc_idx(TE),2) ...
                           .*abs(g.exc.Efd(g.exc.dc_idx(TE),k)));

                g.exc.dEfd(g.exc.dc_idx(TE),k) = ...
                    (g.exc.V_R(g.exc.dc_idx(TE),k) ...
                     - (g.exc.exc_con(g.exc.dc_idx(TE),10) + SE) ...
                       .*g.exc.Efd(g.exc.dc_idx(TE),k)) ...
                    ./g.exc.exc_con(g.exc.dc_idx(TE),11);
            end

            % rate feedback state derivative
            g.exc.dR_f(g.exc.dc_idx,k) = ...
                (-g.exc.R_f(g.exc.dc_idx,k) + g.exc.Efd(g.exc.dc_idx,k)) ...
                ./g.exc.exc_con(g.exc.dc_idx,17);
        end
    end
end

end  % function end

% eof
