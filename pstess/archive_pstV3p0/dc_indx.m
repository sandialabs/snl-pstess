function f = dc_indx(bus,line,dci_dc,dcr_dc)
%Syntax: f = dc_indx(bus,line)
%Purpose: To form indexes for the rectifier and inverter in the dc load flow
%         and indicate the ac buses contected to the converters
%Input:   dc converter specification matrix
%         dc line matrix
%         ac bus matrix
%         ac line matrix 
%         dcr_dc - user defined damping control at rectifier cell
%         dci_dc - user defined damping control at inverter cell
%Output:  f is a dummy variable
%Modified:
%Date: January 1999
%Author: Graham Rogers
%Purpose: Addition of user defined damping controls

%Author:  Graham Rogers
%Date:    October 1996
%         (c) Copyright Joe Chow 1996 - All right reserved
%

f = 0;
global  bus_int
global  dcsp_con  dcl_con dcc_con load_con
global  r_idx  i_idx n_dcl  n_conv  ac_bus rec_ac_bus  inv_ac_bus
global  inv_ac_line  rec_ac_line ac_line dcli_idx
global  ric_idx  rpc_idx
global  no_cap_idx  cap_idx  no_ind_idx  l_no_cap  l_cap
global  ldc_idx
global  ndcr_ud ndci_ud dcrud_idx dciud_idx
% pick out ac voltages (should be LT converter transformer buses)
% check that line and cont data is consistent
n_conv = 0;
n_dcl = 0;
dcrud_idx=[];
dciud_idx=[];
ndcr_ud=0;
ndci_ud=0;
if ~isempty(dcsp_con)
   lconv = length(dcsp_con(:,1));
   lline = length(dcl_con(:,1));
   lcon = length(dcc_con(:,1));
   if lcon~=2*lline|lcon~=lconv
      nc = num2str(lconv);
      nl = num2str(lline);
      ncon = num2str(lcon); 
      disp('number of converters = ',nc)
      disp('number of lines = ',nl)
      disp('number of controls = ',ncon)
      error('dc converter and line data inconsistent')
   end
   % find index of rectifier buses
   r_idx = find(dcsp_con(:,3)==1);
   % find index of inverter buses
   i_idx = find(dcsp_con(:,3)==2);
   % find index of recitifier current control
   ric_idx = find(dcc_con(:,9)==1);
   % find index of rectifier power control
   rpc_idx = find(dcc_con(:,9)==2);
   n_dcl = lline;
   n_conv = lconv;
   inv_ac_bus = bus_int(dcsp_con(i_idx,2));
   rec_ac_bus = bus_int(dcsp_con(r_idx,2));
   ac_bus = bus_int(dcsp_con(:,2));
   inv_ac_line = zeros(lline,1);
   rec_ac_line = zeros(lline,1);
   for j = 1:lline
      acilj = find(bus_int(line(:,2))== inv_ac_bus(j));
      if isempty(acilj)
         error(' the inverter bus is not declared as a to bus')
      else
         inv_ac_line(j) = acilj;
         acilj = [];
      end
      acrlj = find(bus_int(line(:,2)) == rec_ac_bus(j));
      if isempty(acrlj)
         error(' the rectifier bus is not declared as a to bus')
      else
         rec_ac_line(j) = acrlj;
         acrlj = [];
      end
   end
   ac_line = [rec_ac_line;inv_ac_line];
   % form index of dc lines associated with the inverters
   dcli_idx = zeros(n_dcl,1);
   for k = 1:n_dcl
      dcli_idx = dcli_idx|(dcl_con(k,2)==dcsp_con(i_idx,1));
   end
   dcli_idx = find(dcli_idx~=0);
end
no_cap_idx = find(dcl_con(:,5)==0);
cap_idx = find(dcl_con(:,5)~=0);
l_no_cap = 0;
if ~isempty(no_cap_idx); 
   l_no_cap = length(no_cap_idx);
end
l_cap = n_dcl-l_no_cap;
no_ind_idx = find(dcl_con(:,4) ==0|dcl_con(:,6)==0|dcl_con(:,7)==0);

% index of converters in load_con
j = bus_int(load_con(:,1));
for k = 1: n_conv
   ldc_idx(k) = find(j==ac_bus(k));
   if isempty(ldc_idx(k))
      error('dc converter LT buses must be declared in load_con')
      % an additional load is not allowed at a converter bus
   end
end
% j(ldc_idx(k)) gives  the internal bus number of the kth converter 
% check for user defined controls
[ndcr_ud,dummy] = size(dcr_dc);
for j = 1:ndcr_ud
    if ~isempty(dcr_dc{j})
        dcrud_idx(j) = dcr_dc{j,2};
    end
    dcrud_idx(j) = 0;
end     
[ndci_ud,dummy] = size(dci_dc);
for j = 1:ndci_ud
    if ~isempty(dci_dc{j})
      dciud_idx(j) = dci_dc{j,2};
  end
  dciud_idx(j) = 0;
end     

return