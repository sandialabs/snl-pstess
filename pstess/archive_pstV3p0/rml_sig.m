function f = rml_sig(t,k)
% Syntax: f = rml_sig(t,k)
%5:43 PM 27/8/97
% defines modulation signal for rlmod control
global rlmod_sig n_rlmod
f=0; %dummy variable
rlmod_sig(:,k)=zeros(n_rlmod,1);
%if n_rlmod~=0
%    rlmod_sig(:,k)=zeros(n_rlmod,1);
%    rlmod_sig(:,k) = 0.1*randn(n_rlmod,1);
%  if t>=0.0
%    rlmod_sig(1,k)=.1;
%  end
%end
