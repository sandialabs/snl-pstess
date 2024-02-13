function mac_indx()
% Syntax: mac_indx()
%
% Purpose: Forms indexes for the machine models to allow vector computation
%          with different generator types

%-----------------------------------------------------------------------------%

global g;  % declaring struct of global variables

[g.mac.n_mac, n_par] = size(g.mac.mac_con);

g.mac.mac_pot = zeros(g.mac.n_mac,15);
if (n_par < 22)
    g.mac.mac_con(:,22:23) = ones(g.mac.n_mac,2);
end

pqpc_idx = find((g.mac.mac_con(:,22) == 0) & (g.mac.mac_con(:,23) == 0));
if ~isempty(pqpc_idx)
    g.mac.mac_con(pqpc_idx,22:23) = ones(length(pqpc_idx),2);
end

% set up internal machine list
macmax = max(g.mac.mac_con(:,1));
g.mac.mac_int = zeros(macmax,1);
g.mac.mac_int(round(g.mac.mac_con(:,1))) = 1:g.mac.n_mac;

n_tot = g.mac.n_mac;
n_gm = g.mac.n_mac;

g.ind.n_mot = 0;
g.igen.n_ig = 0;

if ~isempty(g.ind.ind_con)
    g.ind.n_mot = length(g.ind.ind_con(:,1));
    n_tot = g.mac.n_mac + g.ind.n_mot;
    n_gm = n_tot;
    motmax = max(g.ind.ind_con(:,1));
    g.ind.ind_int = zeros(motmax,1);
    g.ind.ind_int(round(g.ind.ind_con(:,1))) = g.mac.n_mac+1:n_tot;
end

if ~isempty(g.igen.igen_con)
    g.igen.n_ig = length(g.igen.igen_con(:,1));
    n_tot = n_tot + g.igen.n_ig;
    igmax = max(g.igen.igen_con(:,1));
    g.igen.igen_int = zeros(igmax,1);
    g.igen.igen_int(round(g.igen.igen_con(:,1))) = n_gm+1:n_tot;
end

% check for types of generators
% infinite buses
g.mac.n_ib = 0;
g.mac.n_ib_em = 0;
g.mac.n_ib_tra = 0;
g.mac.n_ib_sub = 0;

% sets default to all generators not infinite buses
g.mac.not_ib_idx = (1:g.mac.n_mac)';
g.mac.mac_ib_idx = [];

if ~isempty(g.mac.ibus_con)
    if (length(g.mac.ibus_con) ~= g.mac.n_mac)
        error('mac_indx: the length of ibus_con does not match mac_con');
    end
    g.mac.mac_ib_idx = find(g.mac.ibus_con == 1);
    g.mac.not_ib_idx = find(g.mac.ibus_con == 0);
    g.mac.n_ib = length(g.mac.mac_ib_idx);
end

% em has no xd or xq
g.mac.mac_em_idx = find(g.mac.mac_con(g.mac.not_ivm_idx,6) == 0);

if ~isempty(g.mac.mac_em_idx)
    g.mac.n_em = length(g.mac.mac_em_idx);
else
    g.mac.n_em = 0;
end

% tra has no xdpp
g.mac.mac_tra_idx = find((g.mac.mac_con(g.mac.not_ivm_idx,6) ~= 0) ...
                         & (g.mac.mac_con(g.mac.not_ivm_idx,8) == 0));

if ~isempty(g.mac.mac_tra_idx)
    g.mac.n_tra = length(g.mac.mac_tra_idx);
else
    g.mac.n_tra = 0;
end

% sub has xdpp
g.mac.mac_sub_idx = find(g.mac.mac_con(g.mac.not_ivm_idx,8) ~= 0);

if ~isempty(g.mac.mac_sub_idx)
    g.mac.n_sub = length(g.mac.mac_sub_idx);
else
    g.mac.n_sub = 0;
end

% ivm indexing handled in ivm_indx.m

% finding non-ivm infinite bus generators
mac_ib_idx = g.mac.mac_ib_idx(ismember(g.mac.mac_ib_idx,g.mac.not_ivm_idx));

if (g.mac.n_ib ~= 0)
    ib_em = find(g.mac.mac_con(mac_ib_idx,6) == 0);

    if ~isempty(ib_em)
        g.mac.mac_ib_em = mac_ib_idx(ib_em);
        g.mac.n_ib_em = length(ib_em);
    end

    ib_tra = find((g.mac.mac_con(mac_ib_idx,6) ~= 0) ...
                  & (g.mac.mac_con(mac_ib_idx,8) == 0));

    if ~isempty(ib_tra)
        g.mac.mac_ib_tra = mac_ib_idx(ib_tra);
        g.mac.n_ib_tra = length(ib_tra);
    end

    ib_sub = find(g.mac.mac_con(mac_ib_idx,8) ~= 0);

    if ~isempty(ib_sub)
        g.mac.mac_ib_sub = mac_ib_idx(ib_sub);
        g.mac.n_ib_sub = length(ib_sub);
    end
end

g.mac.n_pm = g.mac.n_mac;  % number of pm_sig auxiliary inputs

end  % function end

% eof
