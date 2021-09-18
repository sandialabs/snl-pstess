% script file for exciter state perturbation
% 5:08 PM 15/08/97
% perturb the exciter variable
% part of svm_mgen

k_exc = find(mac_exc==k);
if ~isempty(k_exc)
    if n_smp ~= 0
        chk_smp = find(smp_idx==k_exc);
        if ~isempty(chk_smp)~=0
            disp('disturb simple exciter')
            k_smp = find(smp_idx==k_exc);
            not_TR = 1;
            if ~isempty(smp_TR_idx);
                k_exc_idx = find(smp_TR_idx == k_smp);
                if ~isempty(k_exc_idx);not_TR=0;end
            end
            if not_TR ==0
                j=j+1;
                pert = 0.0001*abs(V_TR(k_exc,1));
                pert = max(pert,0.0001);
                V_TR(k_exc,2)=V_TR(k_exc,1)+pert;
                p_file
                st_name(k,j) = 7;
            end
            not_TB = 1;
            if ~isempty(smp_TB_idx);
                k_exc_idx = find(smp_TB_idx == k_smp);
                if ~isempty(k_exc_idx);not_TB=0;end
            end
            if not_TB ==0
                j=j+1;
                pert = 0.0001*abs(V_As(k_exc,1));
                pert = max(pert,0.0001);
                V_As(k_exc,2) = V_As(k_exc,1) + pert;
                p_file
                st_name(k,j) = 8;
            end
            not_TA = 1;
            if ~isempty(smp_TA_idx);
                k_exc_idx = find(smp_TA_idx == k_smp);
                if ~isempty(k_exc_idx);not_TA=0;end
            end
            if not_TA ==0
                j=j+1;
                pert = 0.0001*abs(Efd(k_exc,1));
                pert = max(pert,0.0001);
                Efd(k_exc,2) = Efd(k_exc,1) + pert;
                p_file
                st_name(k,j) = 10;
            end
        end
    end
    if n_smppi ~= 0
        chk_smppi = find(smppi_idx==k_exc);
        if ~isempty(chk_smppi)~=0
            disp('disturb simple pi exciter')
            k_smppi = find(smppi_idx==k_exc);
            not_TR = 1;
            if ~isempty(smppi_TR_idx);
                k_exc_idx = find(smppi_TR_idx == k_smppi);
                if ~isempty(k_exc_idx);not_TR=0;end
            end
            if not_TR ==0
                j=j+1;
                pert = 0.0001*abs(V_TR(k_exc,1));
                pert = max(pert,0.0001);
                V_TR(k_exc,2)=V_TR(k_exc,1)+pert;
                p_file
                st_name(k,j) = 7;
            end
            j=j+1;
            pert = 0.0001*abs(V_As(k_exc,1));
            pert = max(pert,0.0001);
            V_As(k_exc,2) = V_As(k_exc,1) + pert;
            p_file
            st_name(k,j) = 8;
            j=j+1;
            pert = 0.0001*abs(Efd(k_exc,1));
            pert = max(pert,0.0001);
            Efd(k_exc,2) = Efd(k_exc,1) + pert;
            p_file
            st_name(k,j) = 9;
        end
    end

    if n_dc ~= 0
        chk_dc = find(dc_idx==k_exc);
        if ~isempty(chk_dc)
            disp('disturb dc exciter')
            k_dc = find(dc_idx==k_exc);
            not_TR = 1;
            if ~isempty(dc_TR_idx);
                k_exc_idx = find(dc_TR_idx == k_dc);
                if ~isempty(k_exc_idx);not_TR=0;end
            end
            if not_TR ==0
                j=j+1;
                pert = 0.0001*abs(V_TR(k_exc,1));
                pert = max(pert,0.0001);
                V_TR(k_exc,2)=V_TR(k_exc,1)+pert;
                p_file
                st_name(k,j) = 7;
            end
            not_TB = 1;
            if ~isempty(dc_TB_idx);
                k_exc_idx = find(dc_TB_idx == k_dc);
                if ~isempty(k_exc_idx);not_TB=0;end
            end
            if not_TB ==0
                j=j+1;
                pert = 0.0001*abs(V_As(k_exc,1));
                pert = max(pert,0.0001);
                V_As(k_exc,2) = V_As(k_exc,1) + pert;
                p_file
                st_name(k,j) = 8;
            end
            not_TA = 1;
            if ~isempty(dc_TA_idx);
                k_exc_idx = find(dc_TA_idx == k_dc);
                if ~isempty(k_exc_idx);not_TA=0;end
            end
            if not_TA ==0

                j=j+1;
                pert = 0.0001*abs(V_R(k_exc,1));
                pert = max(pert,0.0001);
                V_R(k_exc,2) = V_R(k_exc,1) + pert;
                p_file
                st_name(k,j) = 9;
            end
            not_TE = 1;
            if ~isempty(dc_TE_idx);
                k_exc_idx = find(dc_TE_idx == k_dc);
                if ~isempty(k_exc_idx);not_TE=0;end
            end
            if not_TE ==0
                j=j+1;
                pert = 0.0001*abs(Efd(k_exc,1));
                pert = max(pert,0.0001);
                Efd(k_exc,2) = Efd(k_exc,1) + pert;
                p_file
                st_name(k,j) = 10;
            end
            not_TF = 1;
            if ~isempty(dc_TF_idx);
                k_exc_idx = find(dc_TF_idx == k_dc);
                if ~isempty(k_exc_idx);not_TF=0;end
            end
            if not_TF ==0
                j=j+1;
                pert = 0.0001*abs(R_f(k_exc,1));
                pert = max(pert,0.0001);
                R_f(k_exc,2) = R_f(k_exc,1) + pert;
                p_file
                st_name(k,j) = 11;
            end
        end
    end
    if n_st3 ~= 0
        chk_st3 = find(st3_idx==k_exc);
        if ~isempty(chk_st3)~=0
            disp('disturb st3 exciter')
            k_st3 = find(st3_idx==k_exc);
            not_TR = 1;
            if ~isempty(st3_TR_idx);
                k_exc_idx = find(st3_TR_idx == k_st3);
                if ~isempty(k_exc_idx);not_TR=0;end
            end
            if not_TR == 0
                j=j+1;
                pert = 0.0001*abs(V_TR(k_exc,1));
                pert = max(pert,0.0001);
                V_TR(k_exc,2)=V_TR(k_exc,1)+pert;
                p_file
                st_name(k,j) = 7;
            end
            not_TB = 1;
            if ~isempty(st3_TB_idx);
                k_exc_idx = find(st3_TB_idx == k_st3);
                if ~isempty(k_exc_idx);not_TB=0;end
            end
            if not_TB == 0
                j=j+1;
                pert = 0.0001*abs(V_As(k_exc,1));
                pert = max(pert,0.0001);
                V_As(k_exc,2) = V_As(k_exc,1) + pert;
                p_file
                st_name(k,j) = 8;
            end
            not_TA = 1;
            if ~isempty(st3_TA_idx);
                k_exc_idx = find(st3_TA_idx == k_st3);
                if ~isempty(k_exc_idx);not_TA=0;end
            end
            if not_TA == 0
                j=j+1;
                pert = 0.0001*abs(V_R(k_exc,1));
                pert = max(pert,0.0001);
                V_R(k_exc,2) = V_R(k_exc,1) + pert;
                p_file
                st_name(k,j) = 9;
            end
        end
    end
end