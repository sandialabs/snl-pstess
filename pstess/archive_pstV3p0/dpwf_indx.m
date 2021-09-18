function f = dpwf_indx;
% syntax: f = dpwf_indx  
% 6:14 pm 7/7/98
% Version:  1
% Author: Graham Rogers
% Date: July 1998
% Purpose: Forms indexes for the deltaP/w filter
%          determines indexs for the generators and pss to which the filters are connected 
% (c) Copyright Joe Chow/Cherry Tree Scientific Software 1998
%     All Rights Reserved
f=0;
global  pss_con  mac_int; 
global  dpw_con dpw_idx dpw_pss_idx n_dpw dpw_Td_idx dpw_Tz_idx dpw_mb_idx;
if ~isempty(dpw_con)
   n_dpw = length(dpw_con(:,1));
   err_ntz = find(dpw_con(:,2)>4);
   if ~isempty(err_ntz)
      error('the maximum number of zeros is 4')
   end
   err_ntd = find(dpw_con(:,3)>5);
   if ~isempty(err_ntd)
      error('the maximum number of poles is 5')
   end
   err_ntd = find(dpw_con(:,3)<1);
   if ~isempty(err_ntd)
      error('the minimum number of poles is 1')
   end

   dpw_Td_idx = zeros(n_dpw,5);
   for kdpw = 1:n_dpw
         dpw_Td_idx(kdpw,1:dpw_con(kdpw,3))= ones(1,dpw_con(kdpw,3));
   end 
   dpw_Tz_idx = zeros(n_dpw,4);
   for kdpw = 1:n_dpw
         dpw_Tz_idx(kdpw,1:dpw_con(kdpw,2))= ones(1,dpw_con(kdpw,2));
   end 

   dpw_mb_idx = mac_int(round(dpw_con(:,1)));
   for jdpw = 1:n_dpw
      dpw_pss_idx(jdpw) = find(round(dpw_con(jdpw,1))==round(pss_con(:,2))); 
      if isempty(dpw_pss_idx(jdpw));error('you must have an pss at the same generator as a dp/w filter');end  
   end
else
  n_dpw = 0;
end