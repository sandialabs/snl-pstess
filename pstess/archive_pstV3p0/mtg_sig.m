function f = mtg_sig(t,k)
% Syntax: f = mtg_sig(t,k)
% 12:37 PM 7/0/98
% defines modulation signal for turbine power reference
global tg_sig n_tg n_tgh
f=0; %dummy variable
if n_tg~=0|n_tgh~=0
  tg_sig(:,k) = zeros(n_tg+n_tgh,1);
  if t<=0.0
     tg_sig(:,k) = zeros(n_tg+n_tgh,1);
  else
     tg_sig(:,k) = zeros(n_tg+n_tgh,1);
     %tg_sig(1,k) = -1.0*t;
     %tg_sig(1,k) = -0.01;
  end
end
return
