function f = chq_lim(qg_max,qg_min)
%Syntax:
%       f = chq_lim(qg_max,qg_min)
% function for detecting generator vars outside limit
% sets Qg to zero if limit exceded, sets Ql to negative of limit
% sets bus_type to 3, and recalculates ang_red and volt_red
% changes  generator bus_type to type 3
% recalculates the generator index
% inputs: qg_max and qg_min are the last two clumns of the bus matrix
% outputs:f is set to zero if no limit reached, or to 1 if a limit is reached
% Version:  1.1
% Author: Graham Rogers
% Date:   May 1997
% Purpose: Addition of var limit index
% Version:  1.0
% Author:   Graham Rogers
% Date:     October 1996
%
% (c) copyright Joe Chow 1996
global Qg bus_type g_bno PQV_no PQ_no ang_red volt_red 
global Ql
global gen_chg_idx
%         gen_chg_idx indicates those generators changed to PQ buses
%         gen_cgq_idx = ones(n of bus,1) if no gen at vars limits
%                     = 0 at the corresponding bus if generator at var limit 

f = 0;
lim_flag = 0;% indicates whether limit has been reached
gen_idx = find(bus_type ==2);
qg_max_idx = find(Qg(gen_idx)>qg_max(gen_idx));
qg_min_idx = find(Qg(gen_idx)<qg_min(gen_idx));
if ~isempty(qg_max_idx)
  %some q excedes maximum
  %set Qg to zero
  Qg(gen_idx(qg_max_idx)) = zeros(length(qg_max_idx),1);
  % modify Ql
  Ql(gen_idx(qg_max_idx)) = Ql(gen_idx(qg_max_idx))...
                            - qg_max(gen_idx(qg_max_idx));
  % modify bus_type to PQ bus
  bus_type(gen_idx(qg_max_idx)) = 3*ones(length(qg_max_idx),1);
  gen_chg_idx(gen_idx(qg_max_idx)) = zeros(length(qg_max_idx),1);
  lim_flag = 1;
end
if ~isempty(qg_min_idx)
  %some q less than minimum
  %set Qg to zero
  Qg(gen_idx(qg_min_idx)) = zeros(length(qg_min_idx),1);
  % modify Ql
  Ql(gen_idx(qg_min_idx)) = Ql(gen_idx(qg_min_idx))...
                            - qg_min(gen_idx(qg_min_idx));
  % modify bus_type to PQ bus
  bus_type(gen_idx(qg_min_idx)) = 3*ones(length(qg_min_idx),1);
  gen_chg_idx(gen_idx(qg_min_idx)) = zeros(length(qg_min_idx),1);
  lim_flag = 1;
end
if lim_flag == 1
  %recalculate g_bno
  nbus = length(bus_type);
  g_bno = ones(nbus,1);
  bus_zeros=zeros(nbus,1);
  bus_index=(1:1:nbus)';
  PQV_no=find(bus_type >=2);
  PQ_no=find(bus_type==3);
  gen_index=find(bus_type==2);
  g_bno(gen_index)=bus_zeros(gen_index); 
  % construct sparse angle reduction matrix
  il = length(PQV_no);
  ii = (1:1:il)';
  ang_red = sparse(ii,PQV_no,ones(il,1),il,nbus);
  % construct sparse voltage reduction matrix
  il = length(PQ_no);
  ii = (1:1:il)';
  volt_red = sparse(ii,PQ_no,ones(il,1),il,nbus);
end
f = lim_flag;
return
