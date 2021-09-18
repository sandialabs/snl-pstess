function f = mac_indx
% Syntax f = mac_indx
% 9:38 am 16 December 1997
% Purpose: Forms indexes for the machine models to allow vector computation 
%          with different generator types
% f is a dummy variable  
f=0;
global mac_pot mac_con mac_int ibus_con n_mac n_em n_tra n_sub n_ib
global mac_em_idx mac_tra_idx mac_sub_idx mac_ib_idx not_ib_idx
global mac_ib_em mac_ib_tra mac_ib_sub n_ib_em n_ib_tra n_ib_sub
global ind_int ind_con n_mot
global igen_int igen_con n_ig
% insert the default values of the %percentage p and q
[n_mac, n_par] = size(mac_con);
mac_pot = zeros(n_mac,15);
if n_par < 22
   mac_con(:,22:23) = ones(n_mac,2);
end
pqpc_idx = find(mac_con(:,22)==0&mac_con(:,23)==0);
if ~isempty(pqpc_idx)
   mac_con(pqpc_idx,22:23)=ones(length(pqpc_idx),2);
end 
% set up internal machine list
macmax = max(mac_con(:,1));
mac_int = zeros(macmax,1);
mac_int(round(mac_con(:,1))) = 1:n_mac;
n_tot = n_mac; ngm = n_mac;
n_mot = 0;
n_ig = 0;
if ~isempty(ind_con)
   n_mot = length(ind_con(:,1));
   n_tot = n_mac + n_mot;
   ngm = n_tot;
   motmax= max(ind_con(:,1));
   ind_int = zeros(motmax,1);
   ind_int(round(ind_con(:,1)))=n_mac+1:n_tot;
end

if ~isempty(igen_con)
   n_ig = length(igen_con(:,1));
   n_tot = n_tot + n_ig;
   igmax= max(igen_con(:,1));
   igen_int = zeros(igmax,1);
   igen_int(round(igen_con(:,1)))=ngm+1:n_tot;
end


% check for types of generators
% infinite buses
n_ib = 0;
n_ib_em = 0;
n_ib_tra = 0;
n_ib_sub = 0;
not_ib_idx = (1:n_mac)';% sets default to all generators not infinite buses
if ~isempty(ibus_con)
   mac_ib_idx = find(ibus_con==1);
   not_ib_idx = find(ibus_con==0); 
   n_ib = length(mac_ib_idx);
end
%em has no xd or xq
mac_em_idx = find(mac_con(:,6)==0);
if ~isempty(mac_em_idx)
   n_em = length(mac_em_idx);
else
   n_em = 0;
end
%tra has no xdpp
mac_tra_idx = find((mac_con(:,6)~=0)&(mac_con(:,8)==0));
if ~isempty(mac_tra_idx)
   n_tra = length(mac_tra_idx);
else
   n_tra = 0;
end
%sub has xdpp
mac_sub_idx = find(mac_con(:,8)~=0);
if ~isempty(mac_sub_idx)
   n_sub = length(mac_sub_idx);
else
   n_sub = 0;
end
if n_ib~=0
   ib_em = find(mac_con(mac_ib_idx,6)==0);
   if ~isempty(ib_em)
      mac_ib_em = mac_ib_idx(ib_em);
      n_ib_em = length(ib_em);
   end
   ib_tra = find((mac_con(mac_ib_idx,6)~=0)&(mac_con(mac_ib_idx,8)==0));
   if ~isempty(ib_tra)
      mac_ib_tra = mac_ib_idx(ib_tra);
      n_ib_tra = length(ib_tra);
   end
   ib_sub = find(mac_con(mac_ib_idx,8)~=0);
   if ~isempty(ib_sub)
      mac_ib_sub = mac_ib_idx(ib_sub);
      ib_sub = length(mac_ib_sub);
   end
end  