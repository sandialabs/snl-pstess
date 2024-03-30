function svc_indx(svc_dc)
% Syntax: svc_indx(svc_dc)
%
% Purpose: Determines the relationship between svc and nc loads
%        - Checks for svc
%        - Determines number of SVCs
%        - Checks for user-defined damping controls

%-----------------------------------------------------------------------------%

global g;  % declaring struct of global variables

g.svc.svc_idx = [];
g.svc.svcll_idx = [];

g.svc.n_svcud = 0;
g.svc.svcud_idx = [];

if ~isempty(g.svc.svc_con)
    [g.svc.n_svc, npar] = size(g.svc.svc_con);
    g.svc.svc_idx = zeros(g.svc.n_svc,1);

    % set defaults for lead-lag
    if (npar < 9)
        g.svc.svc_con(:,8:9) = zeros(g.svc.n_svc,2);
    end

    g.svc.svcll_idx = find(g.svc.svc_con(:,9) ~= 0);
    for j = 1:g.svc.n_svc
        index = find(g.svc.svc_con(j,2) == g.ncl.load_con(:,1));
        if ~isempty(index)
            g.svc.svc_idx(j) = index;
        else
            estr = '\nsvc_indx: the svc at index %0.0f must be declared as ';
            estr = [estr, 'a non-conforming load in load_con.'];
            error(sprintf(estr,j));
        end
    end

    % check for user-defined controls
    if ~isempty(svc_dc)
        [g.svc.n_svcud,~] = size(svc_dc);
        for j = 1:g.svc.n_svcud
            g.svc.svcud_idx(j) = svc_dc{j,2};
        end
    end
end

end  % function end

% eof
