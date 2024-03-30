function tcsc_indx(tcsc_dc)
% Syntax: tcsc_indx(tcsc_dc)
%
% Purpose: Determines the relationship between tcsc and nc loads
%        - Checks for tcsc
%        - Determines number of tcscs
%        - Checks for user-defined damping controls

%-----------------------------------------------------------------------------%

global g;  % declaring struct of global variables

g.tcsc.tcsc_idx = [];
g.tcsc.tcscf_idx = [];
g.tcsc.tcsct_idx = [];

g.tcsc.n_tcscud = 0;
g.tcsc.dtcscud_idx = [];

if ~isempty(g.tcsc.tcsc_con)
    g.tcsc.n_tcsc = size(g.tcsc.tcsc_con,1);
    g.tcsc.tcsc_idx = zeros(g.tcsc.n_tcsc,1);

    for j = 1:g.tcsc.n_tcsc
        index = find(g.tcsc.tcsc_con(j,2) == g.ncl.load_con(:,1));
        if ~isempty(index)
            g.tcsc.tcscf_idx(j) = index;
        else
            estr = '\ntcsc_indx: the ''from bus'' for the tcsc at index %0.0f ';
            estr = [estr, 'must be declared as a non-conforming load in load_con.'];
            error(sprintf(estr,j));
        end
        %
        index = find(g.tcsc.tcsc_con(j,3) == g.ncl.load_con(:,1));
        if ~isempty(index)
            g.tcsc.tcsct_idx(j) = index;
        else
            estr = '\ntcsc_indx: the ''to bus'' for the tcsc at index %0.0f ';
            estr = [estr, 'must be declared as a non-conforming load in load_con.'];
            error(sprintf(estr,j));
        end
    end

    % check for user-defined controls
    if ~isempty(tcsc_dc)
        g.tcsc.n_tcscud = size(tcsc_dc,1);
        for j = 1:g.tcsc.n_tcscud
            g.tcsc.dtcscud_idx(j) = tcsc_dc{j,2};
        end
    end
end

end  % function end

% eof
