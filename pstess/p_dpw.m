% script file
%
% perturb the deltaP/omega filter states
% part of svm_mgen

%-----------------------------------------------------------------------------%

k_dpw = find(mac_dpw == k);
if ~isempty(k_dpw)
    disp('disturbing deltaP/omega filter')

    j = j + 1;
    pert = p_ratio*abs(g.dpw.sdpw1(k_dpw,1));
    pert = max(pert,p_ratio);
    g.dpw.sdpw1(k_dpw,3) = g.dpw.sdpw1(k_dpw,1) + pert;
    g.dpw.sdpw1(k_dpw,4) = g.dpw.sdpw1(k_dpw,1) - pert;

    run('p_file');
    st_name(k,j) = 15;

    j = j + 1;
    pert = p_ratio*abs(g.dpw.sdpw2(k_dpw,1));
    pert = max(pert,p_ratio);
    g.dpw.sdpw2(k_dpw,3) = g.dpw.sdpw2(k_dpw,1) + pert;
    g.dpw.sdpw2(k_dpw,4) = g.dpw.sdpw2(k_dpw,1) - pert;

    run('p_file');
    st_name(k,j) = 16;

    j = j + 1;
    pert = p_ratio*abs(g.dpw.sdpw3(k_dpw,1));
    pert = max(pert,p_ratio);
    g.dpw.sdpw3(k_dpw,3) = g.dpw.sdpw3(k_dpw,1) + pert;
    g.dpw.sdpw3(k_dpw,4) = g.dpw.sdpw3(k_dpw,1) - pert;

    run('p_file');
    st_name(k,j) = 17;

    j = j + 1;
    pert = p_ratio*abs(g.dpw.sdpw4(k_dpw,1));
    pert = max(pert,p_ratio);
    g.dpw.sdpw4(k_dpw,3) = g.dpw.sdpw4(k_dpw,1) + pert;
    g.dpw.sdpw4(k_dpw,4) = g.dpw.sdpw4(k_dpw,1) - pert;

    run('p_file');
    st_name(k,j) = 18;

    j = j + 1;
    pert = p_ratio*abs(g.dpw.sdpw5(k_dpw,1));
    pert = max(pert,p_ratio);
    g.dpw.sdpw5(k_dpw,3) = g.dpw.sdpw5(k_dpw,1) + pert;
    g.dpw.sdpw5(k_dpw,4) = g.dpw.sdpw5(k_dpw,1) - pert;

    run('p_file');
    st_name(k,j) = 19;

    j = j + 1;
    pert = p_ratio*abs(g.dpw.sdpw6(k_dpw,1));
    pert = max(pert,p_ratio);
    g.dpw.sdpw6(k_dpw,3) = g.dpw.sdpw6(k_dpw,1) + pert;
    g.dpw.sdpw6(k_dpw,4) = g.dpw.sdpw6(k_dpw,1) - pert;

    run('p_file');
    st_name(k,j) = 20;
end

% eof
