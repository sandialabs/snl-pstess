function f = pss_indx;
% syntax: f = pss_indx  
% 8:03 AM 23/10/97
% Version:  1.2
% Author: Graham Rogers
% Date: August 1997
% Purpose: Forms indexes for the pss
%          determines indexs for the generators and exciters to which pss are connected 
% (c) Copyright Joe Chow/Cherry Tree Scientific Software 1997
%     All Rights Reserved
f=0;
global  pss_con pss_pot mac_con mac_int exc_con; 
global pss_idx n_pss pss_sp_idx pss_p_idx pss_mb_idx pss_exc_idx;
global pss_T  pss_T2 pss_T4 pss_T4_idx  pss_noT4_idx;
if ~isempty(pss_con)
  pss_idx = find(pss_con(:,1)==1|pss_con(:,1)==2);
  n_pss = length(pss_idx);
  pss_mb_idx = mac_int(round(pss_con(:,2)));
  for jpss = 1:n_pss
     pss_exc_idx(jpss) = find(round(pss_con(jpss,2))==round(exc_con(:,2))); 
     if isempty(pss_exc_idx(jpss));error('you must have an exciter at the same generator as a pss');end  
  end
  if n_pss~=0
     pss_T = pss_con(pss_idx,4);
     pss_T2 = pss_con(pss_idx,6);
     pss_T4 = pss_con(pss_idx,8);
     pss_T4_idx = find(pss_T4>0.001);
     pss_noT4 = find(pss_T4<0.001);
     pss_sp_idx = find(pss_con(pss_idx,1)==1);
     pss_p_idx = find(pss_con(pss_idx,1)==2);
  end
else
  n_pss = 0;
end