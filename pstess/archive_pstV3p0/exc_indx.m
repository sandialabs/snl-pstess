function f = exc_indx
% Syntax: f = exc_indx
% 9.02am 3/7/98
% forms indexes from exc_con to allow different 
% exciter models to be used while retaining the vector option
% output is a dummy variable

% Version: 2.0
% Date:    July 1998
% Author:  Graham Rogers
% Modification: eliminate length checks in favour of isempty checks

% Version: 1.0
% Date:    June 1996
% Author:  Graham Rogers
% (c) copyright Joe Chow 1996, All rights reserved
f = 0;
%
global exc_pot exc_con n_exc 
% simple exciter
global smp_idx n_smp smppi_idx n_smppi dc_idx n_dc  dc2_idx n_dc2 st3_idx n_st3
global smp_TA smp_TA_idx smp_noTA_idx smp_TB smp_TB_idx smp_noTB_idx
global smp_TR smp_TR_idx smp_noTR_idx smppi_TR smppi_TR_idx smppi_noTR_idx
% dc exciter
global dc_TA dc_TA_idx dc_noTA_idx dc_TB dc_TB_idx dc_noTB_idx;
global dc_TE dc_TE_idx dc_noTE_idx;
global dc_TF dc_TF_idx  dc_TR dc_TR_idx dc_noTR_idx;
% st3 exciter
global st3_TA st3_TA_idx st3_noTA_idx st3_TB st3_TB_idx st3_noTB_idx;
global st3_TR st3_TR_idx st3_noTR_idx;

n_smp = 0;n_smppi =0;n_dc = 0; n_dc1 = 0.;n_dc2 = 0; n_st3 = 0; n_exc = 0;
if ~isempty(exc_con)
   %check for simple exciters 
   smp_idx = find(exc_con(:,1)== 0);
   if ~isempty(smp_idx);n_smp = length(smp_idx);end
    smppi_idx = find(exc_con(:,1)== 4);
   if ~isempty(smppi_idx);n_smppi = length(smppi_idx);end

   %check for dc exciters
   dc_idx = find((exc_con(:,1) == 1)|(exc_con(:,1)==2));
   if ~isempty(dc_idx);n_dc = length(dc_idx);end
   dc1_idx = find(exc_con(:,1)==1);
   if ~isempty(dc1_idx);n_dc1 = length(dc1_idx);end
   dc2_idx = find(exc_con(:,1)== 2);
   if ~isempty(dc2_idx);n_dc2 = length(dc2_idx);end
   %check for type 3 exciter
   st3_idx = find(exc_con(:,1) == 3);
   if ~isempty(st3_idx);n_st3 = length(st3_idx);end
   %
   %form  vectors for  time constants
   
   if n_smp ~= 0
      % TA
      smp_TA = exc_con(smp_idx,5);
      smp_TA_idx = find(smp_TA>0.001);
      smp_noTA_idx = find(smp_TA<0.001);
      % TB & TC
      smp_TB = exc_con(smp_idx,6);
      smp_TB_idx = find(smp_TB>0.001);
      smp_noTB_idx = find(smp_TB<0.001);
      % TR
      smp_TR = exc_con(smp_idx,3);
      smp_TR_idx = find(smp_TR>0.001);
      smp_noTR_idx = find(smp_TR<0.001);
   end
   if n_smppi~=0
    % TR
      smppi_TR = exc_con(smppi_idx,3);
      smppi_TR_idx = find(smppi_TR>0.001);
      smppi_noTR_idx = find(smppi_TR<0.001);
   end  
   if n_dc ~= 0
      % TA
      dc_TA = exc_con(dc_idx,5);
      dc_TA_idx = find(dc_TA >0.001);
      dc_noTA_idx = find(dc_TA<0.001);
      % TB & TC
      dc_TB = exc_con(dc_idx,6);
      dc_TB_idx = find(dc_TB >0.001);
      dc_noTB_idx = find(dc_TB<0.001);
      % TE
      dc_TE = exc_con(dc_idx,11);
      dc_TE_idx = find(dc_TE>0.001); 
      dc_noTE_idx = find(dc_TE<0.001);
      % TF
      dc_TF = exc_con(dc_idx,17);
      dc_TF_idx = find(dc_TF>0.001);
      % TR
      dc_TR = exc_con(dc_idx,3);
      dc_TR_idx = find(dc_TR>0.001);
      dc_noTR_idx = find(dc_TR<0.001);
   end
   
   if n_st3 ~=0
      % TA
      st3_TA = exc_con(st3_idx,5);
      st3_TA_idx = find(st3_TA>0.001);
      st3_noTA_idx = find(st3_TA<0.001);
      % TB & TC
      st3_TB = exc_con(st3_idx,6);
      st3_TB_idx = find(st3_TB>0.001);
      st3_noTB_idx = find(st3_TB<0.001);
      % TR
      st3_TR = exc_con(st3_idx,3);
      st3_TR_idx = find(st3_TR>0.001);
      st3_noTR_idx = find(st3_TR<0.001);
   end
   %set size of exc_pot
   n_exc = n_smp+n_smppi+n_dc+n_st3;
   exc_pot = zeros(n_exc,5);
end
