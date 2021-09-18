function tg_indx()
% Syntax: tg_indx()
%
% Purpose: forms indexes from tg_con to allow different
%          governor models to be used while retaining the vector option

%-----------------------------------------------------------------------------%

global g;  % declaring struct of global variables

g.tg.n_tg = 0;
g.tg.tg_idx = [];

g.tg.n_tgh = 0;
g.tg.tgh_idx = [];

g.tg.n_tg_tot = 0;

if ~isempty(g.tg.tg_con)
    g.tg.tg_idx = find(g.tg.tg_con(:,1) == 1);
    if ~isempty(g.tg.tg_idx)
        g.tg.n_tg = length(g.tg.tg_idx);
    else
        g.tg.n_tg = 0;
    end

    g.tg.tgh_idx = find(g.tg.tg_con(:,1) == 2);
    if ~isempty(g.tg.tgh_idx)
        g.tg.n_tgh = length(g.tg.tgh_idx);
    else
        g.tg.n_tgh = 0;
    end
end

g.tg.n_tg_tot = g.tg.n_tg + g.tg.n_tgh;

end  % function end

% eof
