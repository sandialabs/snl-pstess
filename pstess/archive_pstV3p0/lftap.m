% m-file for calculating on load tap changing
% called by Loadflow
% calculates desired tap assuming HT(from bus) voltage fixed
% and the to bus voltage is set to the voltage limit
% the tap is set to an allowable tap seting
% and the line matrix tap setting entry is changed
% if all voltages are within limits the loadflow is taken as converged
%
%Author  Graham Rogers
%Date    October 1996   

max_v_idx = find(V>=volt_max);
min_v_idx = find(V<=volt_min);
mm_chk = 0;
if (~isempty(max_v_idx)|~isempty(min_v_idx))
   mm_chk = 1;
   if ~isempty(max_v_idx)
      % change on load taps to correct voltage
      % assumes that bus to be corrected is the to bus
      for fb = 1 : length(max_v_idx)
         chk_fb = find(line(:,2) == bus(max_v_idx(fb),1));
         if ~isempty(chk_fb)
            % freeze dc taps
            if ~isempty(n_dcl)
               for kt = 1:2*n_dcl 
                  dc_chk = find(chk_fb==ac_line(kt));
                  if ~isempty(dc_chk);chk_fb(dc_chk)=[];end
               end
            end
         end
         if ~isempty(chk_fb)
            if line(chk_fb,8)~=0
               % can change tap
               % get from bus voltage
               disp('voltage high changing tap on line');disp(chk_fb)
               vm2 = volt_max(max_v_idx(fb));
               verror = V(max_v_idx(fb))-vm2;
               % voltage too high, tap must be increased
               if verror/vm2<line(chk_fb,10)
                  % increase tap by one step
                  tap = line(chk_fb,6) + line(chk_fb,10);
               else
                  tap = line(chk_fb,6) + verror/vm2;
                  tap_set = ceil( (tap-line(chk_fb,9))/line(chk_fb,10));
                  tap = tap_set*line(chk_fb,10) + line(chk_fb,9);
               end
               tap = min(line(chk_fb,8),max(tap, line(chk_fb,9)));
               % reset tap in line data
               disp('tap reset to');tap
               line(chk_fb,6) = tap;
            end
         end
      end
   end 
   if ~isempty(min_v_idx)
      % change on load taps to correct voltage
      % assumes that bus to be corrected is the to bus
      for fb = 1 : length(min_v_idx)
         chk_fb = find(line(:,2) == bus(min_v_idx(fb),1));
         if ~isempty(chk_fb)
            % freeze dc taps
            if ~isempty(n_dcl)
               for kt = 1:2*n_dcl 
                  dc_chk = find(chk_fb==ac_line(kt));
                  if ~isempty(dc_chk);chk_fb(dc_chk)=[];end
               end
            end
         end
         
         if ~isempty(chk_fb)
            if line(chk_fb,8)~=0
               % can change tap
               disp('voltage low changing tap on line');disp(chk_fb)
               vm2 = volt_min(min_v_idx(fb));
               verror = vm2 - V(min_v_idx(fb));
               % voltage too low tap must be reduced
               if verror/vm2<line(chk_fb,10)
                  % reduce tap by one increment
                  tap = line(chk_fb,6)-line(chk_fb,10);
               else
                  tap = line(chk_fb,6) - verror./vm2;
                  tap_set = fix( (tap-line(chk_fb,9))/line(chk_fb,10));
                  tap = tap_set*line(chk_fb,10) + line(chk_fb,9);
               end
               tap = min(line(chk_fb,8),max(tap, line(chk_fb,9)));
               disp('taps reset to');tap
               % reset tap in line data
               line(chk_fb,6) = tap;
            end
         end
      end
   end
end