function f = mexc_sig(t,k)
% Syntax: f = mexc_sig(t,k)
% 1:20 PM 15/08/97
% defines modulation signal for exciter control
global exc_sig n_exc
f=0; %dummy variable
if n_exc~=0
%  exc_sig(:,k)=zeros(n_exc,1);
%  exc_sig(1,k)=0.1;
%end
 if t<=0
     exc_sig(:,k) = zeros(n_exc,1);
 else
    exc_sig(:,k) = zeros(n_exc,1);
    %exc_sig(1,k) = 0.05;
 end
end
return