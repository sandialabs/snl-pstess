function f = mdc_sig(t,k)
% Syntax: f = mdc_sig(t,k)
% 4:40 PM 21/08/97
% defines modulation signal for dc converter control
global dc_sig  r_idx i_idx n_conv
f=0; %dummy variable
if n_conv~=0
  %if t<=0.1
    % dc_sig(:,k) = zeros(n_conv,1);
  %else
  %dc_sig(1,k) = 0.01;dc_sig(1,k+1)=0.01;
  %dc_sig(2,k) = 0.0;dc_sig(2,k+1)=0.0;
  %end
  dc_sig(:,k) = zeros(n_conv,1);
end
return