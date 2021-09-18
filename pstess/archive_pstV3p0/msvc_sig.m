function f = msvc_sig(t,k)
% Syntax: f = msvc_sig(t,k)
% 4:39 PM 15/08/97
% defines modulation signal for svc control
global svc_sig n_svc
f=0; %dummy variable
if n_svc ~=0
  if t<=0.1
     svc_sig(:,k) = zeros(n_svc,1);
  else
     svc_sig(:,k) = zeros(n_svc,1);
  end
end
return