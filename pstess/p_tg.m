% script file
%
% perturb the turbine governor states
% part of svm_mgen

%-----------------------------------------------------------------------------%

k_tg = find(mac_tg == k);
if ~isempty(k_tg)
    disp('disturbing turbine governor')
    gh = g.tg.tg_idx(k_tg);

    j = j + 1;
    pert = p_ratio*abs(g.tg.tg1(gh,1));
    pert = max(pert,p_ratio);
    g.tg.tg1(gh,3) = g.tg.tg1(gh,1) + pert;
    g.tg.tg1(gh,4) = g.tg.tg1(gh,1) - pert;

    run('p_file');
    st_name(k,j) = 21;

    j = j + 1;
    pert = p_ratio*abs(g.tg.tg2(gh,1));
    pert = max(pert,p_ratio);
    g.tg.tg2(gh,3) = g.tg.tg2(gh,1) + pert;
    g.tg.tg2(gh,4) = g.tg.tg2(gh,1) - pert;

    run('p_file');
    st_name(k,j) = 22;

    j = j + 1;
    pert = p_ratio*abs(g.tg.tg3(gh,1));
    pert = max(pert,p_ratio);
    g.tg.tg3(gh,3) = g.tg.tg3(gh,1) + pert;
    g.tg.tg3(gh,4) = g.tg.tg3(gh,1) - pert;

    run('p_file');
    st_name(k,j) = 23;
end

k_tgh = find(mac_tgh == k);
if ~isempty(k_tgh)
    disp('disturbing hydro turbine governor')
    gh = g.tg.tgh_idx(k_tgh);

    j = j + 1;
    pert = p_ratio*abs(g.tg.tg1(gh,1));
    pert = max(pert,p_ratio);
    g.tg.tg1(gh,3) = g.tg.tg1(gh,1) + pert;
    g.tg.tg1(gh,4) = g.tg.tg1(gh,1) - pert;

    run('p_file');
    st_name(k,j) = 21;

    j = j + 1;
    pert = p_ratio*abs(g.tg.tg2(gh,1));
    pert = max(pert,p_ratio);
    g.tg.tg2(gh,3) = g.tg.tg2(gh,1) + pert;
    g.tg.tg2(gh,4) = g.tg.tg2(gh,1) - pert;

    run('p_file');
    st_name(k,j) = 22;

    j = j + 1;
    pert = p_ratio*abs(g.tg.tg3(gh,1));
    pert = max(pert,p_ratio);
    g.tg.tg3(gh,3) = g.tg.tg3(gh,1) + pert;
    g.tg.tg3(gh,4) = g.tg.tg3(gh,1) - pert;

    run('p_file');
    st_name(k,j) = 23;

    j = j + 1;
    pert = p_ratio*abs(g.tg.tg4(gh,1));
    pert = max(pert,p_ratio);
    g.tg.tg4(gh,3) = g.tg.tg4(gh,1) + pert;
    g.tg.tg4(gh,4) = g.tg.tg4(gh,1) - pert;

    run('p_file');
    st_name(k,j) = 24;

    j = j + 1;
    pert = p_ratio*abs(g.tg.tg5(gh,1));
    pert = max(pert,p_ratio);
    g.tg.tg5(gh,3) = g.tg.tg5(gh,1) + pert;
    g.tg.tg5(gh,4) = g.tg.tg5(gh,1) - pert;

    run('p_file');
    st_name(k,j) = 25;
end

% eof
