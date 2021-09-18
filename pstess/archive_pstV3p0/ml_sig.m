function f = ml_sig(t,k)
% Syntax: f = ml_sig(t,k)
%4:40 PM 15/08/97
% defines modulation signal for lmod control
global lmod_sig n_lmod
f=0; %dummy variable
lmod_sig(:,k)=zeros(n_lmod,1);
%if n_lmod~=0
%    lmod_sig(:,k)=zeros(n_lmod,1);
    %lmod_sig(:,k) = 0.1*randn(n_lmod,1);
%    if t>=0.0
%       lmod_sig(1,k)=.1;
%    elseif t>0.2
%        lmod_sig(:,k)=zeros(n_lmod,1);
%    end
end