function f = svc_indx(svc_dc)
% syntax: f = svc_indx
% 12:00 PM 8/8/97
% determines the relationship between svc and nc loads
% checks for svc
% determines number of SVCs
% checks for user defined damping controls
% f is a dummy variable
f = 0;
global svc_con load_con  n_svc  svc_idx svcll_idx % svc 
global n_dcud dcud_idx  %user defined damping controls
n_svc = 0;
svc_idx = [];
dcud_idx = [];
n_dcud = 0;
if ~isempty(svc_con)
   [n_svc npar] = size(svc_con);
   svc_idx = zeros(n_svc,1);
   % set defaults for lead lag
   if npar<9
      svc_con(:,8:9) = zeros(n_svc,2);
   end
   svcll_idx = find(svc_con(:,9)~=0);
   for j = 1:n_svc
      index = find(svc_con(j,2)==load_con(:,1));
      if ~isempty(index)
         svc_idx(j) = index;
      else
         error('you must have the svc bus declared as a non-conforming load')
      end
   end
   % check for user defined controls
   if ~isempty(svc_dc)
      [n_dcud,dummy] = size(svc_dc);
      for j = 1:n_dcud
         dcud_idx(j) = svc_dc{j,2};
      end
   end      
end

