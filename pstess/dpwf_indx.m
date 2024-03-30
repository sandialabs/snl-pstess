function dpwf_indx()
% Syntax: dpwf_indx()
%
% Purpose: Forms indexes for the dpw (deltaP/omega) filter and determines
%          indexes for the generators and pss to which the filters are
%          connected

%-----------------------------------------------------------------------------%
% Version: 1.0 (initial version)
% Author:  Graham Rogers
% Date:    July 1998
%-----------------------------------------------------------------------------%

global g;  % declaring struct of global variables

if ~isempty(g.dpw.dpw_con)
    g.dpw.n_dpw = size(g.dpw.dpw_con,1);

    err_ntz = find(g.dpw.dpw_con(:,2) > 4);
    if ~isempty(err_ntz)
        estr = '\ndpwf_indx: the maximum number of zeros is 4, ';
        estr = [estr, 'revise parameters for dpwf index %0.0f.'];
        error(sprintf(estr,err_ntz));
    end

    err_ntd = find(g.dpw.dpw_con(:,3) > 5);
    if ~isempty(err_ntd)
        estr = '\ndpwf_indx: the maximum number of poles is 5, ';
        estr = [estr, 'revise parameters for dpwf index %0.0f.'];
        error(sprintf(estr,err_ntd));
    end

    err_ntd = find(g.dpw.dpw_con(:,3) < 1);
    if ~isempty(err_ntd)
        estr = '\ndpwf_indx: the minimum number of poles is 1, ';
        estr = [estr, 'revise parameters for dpwf index %0.0f.'];
        error(sprintf(estr,err_ntd));
    end

    g.dpw.dpw_Td_idx = zeros(g.dpw.n_dpw,5);
    for kdpw = 1:g.dpw.n_dpw
        g.dpw.dpw_Td_idx(kdpw,1:g.dpw.dpw_con(kdpw,3)) = ...
            ones(1,g.dpw.dpw_con(kdpw,3));
    end

    g.dpw.dpw_Tz_idx = zeros(g.dpw.n_dpw,4);
    for kdpw = 1:g.dpw.n_dpw
        g.dpw.dpw_Tz_idx(kdpw,1:g.dpw.dpw_con(kdpw,2)) = ...
            ones(1,g.dpw.dpw_con(kdpw,2));
    end

    g.dpw.dpw_mb_idx = g.mac.mac_int(round(g.dpw.dpw_con(:,1)));
    for jdpw = 1:g.dpw.n_dpw
        g.dpw.dpw_pss_idx(jdpw) = ...
            find(round(g.dpw.dpw_con(jdpw,1)) == round(g.pss.pss_con(:,2)));

        if isempty(g.dpw.dpw_pss_idx(jdpw))
            estr = '\ndpwf_indx: there must be a pss at the same generator as ';
            estr = [estr, 'each dP/omega filter (dpwf index %0.0f).'];
            error(sprintf(estr,jdpw));
        end
    end
end

end  % function end

% eof
