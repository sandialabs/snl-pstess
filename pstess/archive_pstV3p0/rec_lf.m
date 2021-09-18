function [alpha,mode,idc_new] = rec_lf(idc,al_min,al_max,Vdo,Rc,cm)
%Syntax: [alpha,mode,idc_new] = rec_lf(idc,al_min,al_max,Vdo,RC,cm)
% 8/12/97
%Purpose: Finds alpha for rectifier and sets rectifier taps
%         determines if mode change is necessary
%         sets new idc
%Inputs:  
%         Note: when mode(i) = 2, idc(i) should be reduced by the current margin
%         idc specified dc line current vector
%         Vdc_max vector of maximum allowed inverter dc voltages
%         Vdc_min vector of minimum allowed inverter dc voltages
%         al_min vector of minimum alpha
%         al_max vector of maximum alpha
%         Vdo vector of ideal rectifier dc voltages
%         Rc vector of rectifier commutating resistances
%         cm vector of current margins
%Output:  alpha, the rectifier firing angle in radians
%         mode - a vector giving the mode of operation of the inverter
%         mode(i) = 1 indicates that the ith inverter is operating at gamma min
%         mode(i) = 2 indicates that the ith inverter is controlling current
%         idc_new is set to specified idc when mode is 1
%         or to idc reduced by the current margin if mode is 2
%Author:  Graham Rogers
%Date:   November 1996
%        (c) Copyright Joe Chow 1996 - All rights Reserved
global  r_idx i_idx n_dcl tapr tmaxr tminr tstepr Vdc
% calculate alpha

Vdcr = Vdc(r_idx);
mode = ones(n_dcl,1);% set to default initially
calpha = (Vdcr+idc.*Rc)./Vdo;
% check that alpha is within range
max_idx = find(calpha<cos(al_max));
min_idx = find(calpha>cos(al_min));
oor_idx = find(calpha<cos(al_max)|calpha>cos(al_min));

if ~isempty(oor_idx)
   Vdo_new = zeros(n_dcl,1);
   tapn = tapr;
   tnum = (tapn - tminr)./tstepr;
   % adjust rectifier taps to bring within range
   if ~isempty(max_idx)
     Vdo_new(max_idx) = (Vdcr(max_idx)+Rc(max_idx).*idc(max_idx))...
             ./cos(al_max(max_idx));
   end
   if ~isempty(min_idx)
     Vdo_new(min_idx) = (Vdcr(min_idx)+Rc(min_idx).*idc(min_idx))...
             ./cos(al_min(min_idx));
   end
   tapn(oor_idx) = Vdo(oor_idx)./Vdo_new(oor_idx);
   % get the right tap setting
   tnum = (tapn - tminr)./tstepr;
   if ~isempty(max_idx)
     tnum(max_idx) = ceil(tnum(max_idx));% next interger up
   end
   if ~isempty(min_idx)
     tnum(min_idx) = fix(tnum(min_idx));% next integer down
   end
   tmin_idx = find(tnum<0);
   tmax_idx = find(tnum>(tmaxr-tminr)./tstepr);

   if ~isempty(tmin_idx)
      if ~isempty(max_idx)
         stop_flag = 0,
         for k = 1:length(max_idx)
            ktm = find(max_idx(k) == tmin_idx)
            if ~isempty(ktm)
               rec_num = num2str(r_idx(max_idx(ktm)));
               disp(['minimum tap setting reached at rectifier ',rec_num])
               stop_flag = 1
            end
         end
         if stp_flag==1;error('stop');end
      end
      if ~isempty(min_idx)
        if ~isempty(min_idx(tmin_idx))
          rec_num = num2str(r_idx(min_idx(tmin_idx)));
          disp('minimum tap setting reached at rectifiers ')
          disp(rec_num)
          disp('changing mode')
          n_mc = length(min_idx(tmin_idx));
          mode(i_idx(min_idx(tmin_idx))) = 2*ones(n_mc,1);
          tnum(tmin_idx) = zeros(n_mc,1);
        end
      end
   end
   if ~isempty(tmax_idx)
     if ~isempty(max_idx)  
       if ~isempty(max_idx(tmax_idx))
           rec_num = num2str(r_idx(max_idx(tmax_idx)));
           disp('maximum tap setting reached at rectifier ')
           disp(rec_num)
           error('stop')
       end
     end
     if ~isempty(min_idx)
       if ~isempty(min_idx(tmax_idx))
           rec_num = num2str(r_idx(min_idx(tmax_idx)));
           disp('maximum tap setting reached at rectifier')
           disp(rec_num)
           disp('changing mode')
           n_mc = length(min_idx(tmax_idx));
           mode(min_idx(tmax_idx)) = 2*ones(n_mc,1);
           tnum(tmax_idx) = (tmaxr(oor_idx(tmax_idx)) - tminr(oor_idx(tmax_idx)))...
                            ./tstepr(oor_idx(tmax_idx));
        end
     end
   end
   tapre = tminr + tstepr.*tnum;   
   % recalculate calpha
   Vdo= Vdo.*tapr./tapre;
   tapr = tapre;
   calpha = (Vdc(r_idx)+idc.*Rc)./Vdo;
   % apply alpha limits
   calpha = min(calpha,cos(al_min));
   calpha = max(calpha,cos(al_max));
end
% adjust idc for mode change
idc_new = idc;
mc_idx = find(mode==2);
if ~isempty(mc_idx)
    idc_new(mc_idx) = idc(mc_idx).*cm(mc_idx);
end
salpha = sqrt(ones(n_dcl,1)-calpha.*calpha);
alpha = atan2(salpha,calpha)*180/pi;
return
