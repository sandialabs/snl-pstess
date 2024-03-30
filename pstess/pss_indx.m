function pss_indx()
% Syntax: pss_indx()
%
% Purpose: Forms indexes for the pss determines indexs for the
%          generators and exciters to which pss are connected

%-----------------------------------------------------------------------------%
% Version history
%
% Version: 1.2
% Author:  Graham Rogers
% Date:    August 1997
%
% Version: 1.0 (initial version)
% Author:  Graham Rogers
% Date:    late 1996/early 1997
%-----------------------------------------------------------------------------%

global g;  % declaring struct of global variables

lbnd = 1e-3;  % time constant lower bound

if ~isempty(g.pss.pss_con)
    g.pss.pss_idx = find((g.pss.pss_con(:,1) == 1) | (g.pss.pss_con(:,1) == 2));
    g.pss.n_pss = length(g.pss.pss_idx);
    g.pss.pss_mb_idx = g.mac.mac_int(round(g.pss.pss_con(:,2)));

    for jpss = 1:g.pss.n_pss
        g.pss.pss_exc_idx(jpss) = ...
            find(round(g.pss.pss_con(jpss,2)) == round(g.exc.exc_con(:,2)));

        if isempty(g.pss.pss_exc_idx(jpss))
            estr = '\npss_indx: exciter model needed for stabilizer at ';
            estr = [estr, 'pss index %0.0f.'];
            error(sprintf(estr,jpss));
        end
    end

    if (g.pss.n_pss ~= 0)
        g.pss.pss_T = g.pss.pss_con(g.pss.pss_idx,4);
        g.pss.pss_T2 = g.pss.pss_con(g.pss.pss_idx,6);
        g.pss.pss_T4 = g.pss.pss_con(g.pss.pss_idx,8);
        g.pss.pss_T4_idx = find(g.pss.pss_T4 >= lbnd);
        g.pss.pss_noT4_idx = find(g.pss.pss_T4 < lbnd);
        g.pss.pss_sp_idx = find(g.pss.pss_con(g.pss.pss_idx,1) == 1);
        g.pss.pss_p_idx = find(g.pss.pss_con(g.pss.pss_idx,1) == 2);
    else
        g.pss.pss_T = [];
        g.pss.pss_T2 = [];
        g.pss.pss_T4 = [];
        g.pss.pss_T4_idx = [];
        g.pss.pss_noT4_idx = [];
        g.pss.pss_sp_idx = [];
        g.pss.pss_p_idx = [];
    end
end

end  % function end

% eof
