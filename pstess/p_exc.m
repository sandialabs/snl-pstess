% script file
%
% perturb the exciter states
% part of svm_mgen

%-----------------------------------------------------------------------------%

k_exc = find(mac_exc == k);
if ~isempty(k_exc)
    if (g.exc.n_smp ~= 0)
        chk_smp = find(g.exc.smp_idx == k_exc);
        if ~isempty(chk_smp)
            disp('disturbing simple exciter')
            k_smp = find(g.exc.smp_idx == k_exc);

            not_TR = 1;
            if ~isempty(g.exc.smp_TR_idx)
                k_exc_idx = find(g.exc.smp_TR_idx == k_smp);
                if ~isempty(k_exc_idx)
                    not_TR = 0;
                end
            end

            if (not_TR == 0)
                j = j + 1;
                pert = p_ratio*abs(g.exc.V_TR(k_exc,1));
                pert = max(pert,p_ratio);
                g.exc.V_TR(k_exc,3) = g.exc.V_TR(k_exc,1) + pert;
                g.exc.V_TR(k_exc,4) = g.exc.V_TR(k_exc,1) - pert;

                run('p_file');
                st_name(k,j) = 7;
            end

            not_TB = 1;
            if ~isempty(g.exc.smp_TB_idx)
                k_exc_idx = find(g.exc.smp_TB_idx == k_smp);
                if ~isempty(k_exc_idx)
                    not_TB = 0;
                end
            end

            if (not_TB == 0)
                j = j + 1;
                pert = p_ratio*abs(g.exc.V_As(k_exc,1));
                pert = max(pert,p_ratio);
                g.exc.V_As(k_exc,3) = g.exc.V_As(k_exc,1) + pert;
                g.exc.V_As(k_exc,4) = g.exc.V_As(k_exc,1) - pert;

                run('p_file');
                st_name(k,j) = 8;
            end

            not_TA = 1;
            if ~isempty(g.exc.smp_TA_idx)
                k_exc_idx = find(g.exc.smp_TA_idx == k_smp);
                if ~isempty(k_exc_idx)
                    not_TA = 0;
                end
            end

            if (not_TA == 0)
                j = j + 1;
                pert = p_ratio*abs(g.exc.Efd(k_exc,1));
                pert = max(pert,p_ratio);
                g.exc.Efd(k_exc,3) = g.exc.Efd(k_exc,1) + pert;
                g.exc.Efd(k_exc,4) = g.exc.Efd(k_exc,1) - pert;

                run('p_file');
                st_name(k,j) = 10;
            end
        end
    end

    if (g.exc.n_smppi ~= 0)
        chk_smppi = find(g.exc.smppi_idx == k_exc);
        if ~isempty(chk_smppi)
            disp('disturbing simple PI exciter')
            k_smppi = find(g.exc.smppi_idx == k_exc);

            not_TR = 1;
            if ~isempty(g.exc.smppi_TR_idx)
                k_exc_idx = find(g.exc.smppi_TR_idx == k_smppi);
                if ~isempty(k_exc_idx)
                    not_TR = 0;
                end
            end

            if (not_TR == 0)
                j = j + 1;
                pert = p_ratio*abs(g.exc.V_TR(k_exc,1));
                pert = max(pert,p_ratio);
                g.exc.V_TR(k_exc,3) = g.exc.V_TR(k_exc,1) + pert;
                g.exc.V_TR(k_exc,4) = g.exc.V_TR(k_exc,1) - pert;

                run('p_file');
                st_name(k,j) = 7;
            end

            j = j + 1;
            pert = p_ratio*abs(g.exc.V_As(k_exc,1));
            pert = max(pert,p_ratio);
            g.exc.V_As(k_exc,3) = g.exc.V_As(k_exc,1) + pert;
            g.exc.V_As(k_exc,4) = g.exc.V_As(k_exc,1) - pert;

            run('p_file');
            st_name(k,j) = 8;

            j = j + 1;
            pert = p_ratio*abs(g.exc.Efd(k_exc,1));
            pert = max(pert,p_ratio);
            g.exc.Efd(k_exc,3) = g.exc.Efd(k_exc,1) + pert;
            g.exc.Efd(k_exc,4) = g.exc.Efd(k_exc,1) - pert;

            run('p_file');
            st_name(k,j) = 10;
        end
    end

    if (g.exc.n_dc ~= 0)
        chk_dc = find(g.exc.dc_idx == k_exc);
        if ~isempty(chk_dc)
            disp('disturbing dc exciter')
            k_dc = find(g.exc.dc_idx == k_exc);

            not_TR = 1;
            if ~isempty(g.exc.dc_TR_idx)
                k_exc_idx = find(g.exc.dc_TR_idx == k_dc);
                if ~isempty(k_exc_idx)
                    not_TR = 0;
                end
            end

            if (not_TR == 0)
                j = j + 1;
                pert = p_ratio*abs(g.exc.V_TR(k_exc,1));
                pert = max(pert,p_ratio);
                g.exc.V_TR(k_exc,3) = g.exc.V_TR(k_exc,1) + pert;
                g.exc.V_TR(k_exc,4) = g.exc.V_TR(k_exc,1) - pert;

                run('p_file');
                st_name(k,j) = 7;
            end

            not_TB = 1;
            if ~isempty(g.exc.dc_TB_idx)
                k_exc_idx = find(g.exc.dc_TB_idx == k_dc);
                if ~isempty(k_exc_idx)
                    not_TB = 0;
                end
            end

            if (not_TB == 0)
                j = j + 1;
                pert = p_ratio*abs(g.exc.V_As(k_exc,1));
                pert = max(pert,p_ratio);
                g.exc.V_As(k_exc,3) = g.exc.V_As(k_exc,1) + pert;
                g.exc.V_As(k_exc,4) = g.exc.V_As(k_exc,1) - pert;

                run('p_file');
                st_name(k,j) = 8;
            end

            not_TA = 1;
            if ~isempty(g.exc.dc_TA_idx)
                k_exc_idx = find(g.exc.dc_TA_idx == k_dc);
                if ~isempty(k_exc_idx)
                    not_TA = 0;
                end
            end

            if (not_TA == 0)
                j = j + 1;
                pert = p_ratio*abs(g.exc.V_R(k_exc,1));
                pert = max(pert,p_ratio);
                g.exc.V_R(k_exc,3) = g.exc.V_R(k_exc,1) + pert;
                g.exc.V_R(k_exc,4) = g.exc.V_R(k_exc,1) - pert;

                run('p_file');
                st_name(k,j) = 9;
            end

            not_TE = 1;
            if ~isempty(g.exc.dc_TE_idx)
                k_exc_idx = find(g.exc.dc_TE_idx == k_dc);
                if ~isempty(k_exc_idx)
                    not_TE = 0;
                end
            end

            if (not_TE == 0)
                j = j + 1;
                pert = p_ratio*abs(g.exc.Efd(k_exc,1));
                pert = max(pert,p_ratio);
                g.exc.Efd(k_exc,3) = g.exc.Efd(k_exc,1) + pert;
                g.exc.Efd(k_exc,4) = g.exc.Efd(k_exc,1) - pert;

                run('p_file');
                st_name(k,j) = 10;
            end

            not_TF = 1;
            if ~isempty(g.exc.dc_TF_idx)
                k_exc_idx = find(g.exc.dc_TF_idx == k_dc);
                if ~isempty(k_exc_idx)
                    not_TF = 0;
                end
            end

            if (not_TF == 0)
                j = j + 1;
                pert = p_ratio*abs(g.exc.R_f(k_exc,1));
                pert = max(pert,p_ratio);
                g.exc.R_f(k_exc,3) = g.exc.R_f(k_exc,1) + pert;
                g.exc.R_f(k_exc,4) = g.exc.R_f(k_exc,1) - pert;

                run('p_file');
                st_name(k,j) = 11;
            end
        end
    end

    if (g.exc.n_st3 ~= 0)
        chk_st3 = find(g.exc.st3_idx == k_exc);
        if ~isempty(chk_st3)
            disp('disturbing st3 exciter')
            k_st3 = find(g.exc.st3_idx == k_exc);

            not_TR = 1;
            if ~isempty(g.exc.st3_TR_idx)
                k_exc_idx = find(g.exc.st3_TR_idx == k_st3);
                if ~isempty(k_exc_idx)
                    not_TR = 0;
                end
            end

            if (not_TR == 0)
                j = j + 1;
                pert = p_ratio*abs(g.exc.V_TR(k_exc,1));
                pert = max(pert,p_ratio);
                g.exc.V_TR(k_exc,3) = g.exc.V_TR(k_exc,1) + pert;
                g.exc.V_TR(k_exc,4) = g.exc.V_TR(k_exc,1) - pert;

                run('p_file');
                st_name(k,j) = 7;
            end

            not_TB = 1;
            if ~isempty(g.exc.st3_TB_idx)
                k_exc_idx = find(g.exc.st3_TB_idx == k_st3);
                if ~isempty(k_exc_idx)
                    not_TB = 0;
                end
            end

            if (not_TB == 0)
                j = j + 1;
                pert = p_ratio*abs(g.exc.V_As(k_exc,1));
                pert = max(pert,p_ratio);
                g.exc.V_As(k_exc,3) = g.exc.V_As(k_exc,1) + pert;
                g.exc.V_As(k_exc,4) = g.exc.V_As(k_exc,1) - pert;

                run('p_file');
                st_name(k,j) = 8;
            end

            not_TA = 1;
            if ~isempty(g.exc.st3_TA_idx)
                k_exc_idx = find(g.exc.st3_TA_idx == k_st3);
                if ~isempty(k_exc_idx)
                    not_TA = 0;
                end
            end

            if (not_TA == 0)
                j = j + 1;
                pert = p_ratio*abs(g.exc.V_R(k_exc,1));
                pert = max(pert,p_ratio);
                g.exc.V_R(k_exc,3) = g.exc.V_R(k_exc,1) + pert;
                g.exc.V_R(k_exc,4) = g.exc.V_R(k_exc,1) - pert;

                run('p_file');
                st_name(k,j) = 9;
            end
        end
    end
end

% eof
