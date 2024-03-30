function rlmod_indx()
% Syntax: rlmod_indx()
%
% Purpose: Determines the relationship between lmod and nc loads
%        - Checks for rlmod
%        - Determines number of modulated loads

%-----------------------------------------------------------------------------%

global g;  % declaring struct of global variables

g.rlmod.rlmod_idx = [];

if ~isempty(g.rlmod.rlmod_con)
    g.rlmod.n_rlmod = size(g.rlmod.rlmod_con,1);
    g.rlmod.rlmod_idx = zeros(g.rlmod.n_rlmod,1);
    for j = 1:g.rlmod.n_rlmod
        index = find(g.rlmod.rlmod_con(j,2) == g.ncl.load_con(:,1));
        if ~isempty(index)
            g.rlmod.rlmod_idx(j) = index;
        else
            estr = '\nrlmod_indx: the reactive load at rlmod index %0.0f must be ';
            estr = [estr, 'declared as a non-conforming load in load_con.'];
            error(sprintf(estr,j));
        end
    end
end

end  % function end

% eof
