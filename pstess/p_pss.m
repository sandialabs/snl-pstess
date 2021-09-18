% script file
%
% perturb the pss states
% part of svm_mgen

%-----------------------------------------------------------------------------%

k_pss = find(mac_pss == k);
if ~isempty(k_pss)
    disp('disturbing pss')

    j = j + 1;
    pert = p_ratio*abs(g.pss.pss1(k_pss,1));
    pert = max(pert,p_ratio);
    g.pss.pss1(k_pss,3) = g.pss.pss1(k_pss,1) + pert;
    g.pss.pss1(k_pss,4) = g.pss.pss1(k_pss,1) - pert;

    run('p_file');
    st_name(k,j) = 12;

    j = j + 1;
    pert = p_ratio*abs(g.pss.pss2(k_pss,1));
    pert = max(pert,p_ratio);
    g.pss.pss2(k_pss,3) = g.pss.pss2(k_pss,1) + pert;
    g.pss.pss2(k_pss,4) = g.pss.pss2(k_pss,1) - pert;

    run('p_file');
    st_name(k,j) = 13;

    if ~isempty(g.pss.pss_T4_idx)
        k_s_idx = find(g.pss.pss_T4_idx == k_pss);
    else
        k_s_idx = [];
    end

    if ~isempty(k_s_idx)
        j = j + 1;
        pert = p_ratio*abs(g.pss.pss3(k_pss,1));
        pert = max(pert,p_ratio);
        g.pss.pss3(k_pss,3) = g.pss.pss3(k_pss,1) + pert;
        g.pss.pss3(k_pss,4) = g.pss.pss3(k_pss,1) - pert;

        run('p_file');
        st_name(k,j) = 14;
    end
end

% eof
