function lmod_indx()
% Syntax: lmod_indx()
%
% Purpose: Determines the relationship between lmod and nc loads
%        - Checks for lmod
%        - Determines number of modulated loads

%-----------------------------------------------------------------------------%

global g;  % declaring struct of global variables

g.lmod.lmod_idx = [];

if ~isempty(g.lmod.lmod_con)
    g.lmod.n_lmod = size(g.lmod.lmod_con,1);
    g.lmod.lmod_idx = zeros(g.lmod.n_lmod,1);
    for j = 1:g.lmod.n_lmod
        index = find(g.lmod.lmod_con(j,2) == g.ncl.load_con(:,1));
        if ~isempty(index)
            g.lmod.lmod_idx(j) = index;
        else
            estr = '\nlmod_indx: the load at lmod index %0.0f must be declared as ';
            estr = [estr, 'a non-conforming load in load_con.'];
            error(sprintf(estr,j));
        end
    end
end

end  % function end

% eof
