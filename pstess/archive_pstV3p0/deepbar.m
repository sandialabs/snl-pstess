function [r,x]=deepbar(rro,B,s)
% gives the equivalent rotor leakage reactance and resistance as a funtion of slip
% B is the deep bar factor
% s is the slip
% rro is the value of rotor resistance at s = 0;
% r is the effective rotor resistance
% x is the effective rotor reactance
% called by mac_ind
b = sqrt(abs(s)).*B;
r0 = rro/2;
if s~=0
	a = (1+i)*b;
	z = r0.*a.*(exp(a)+1)./(exp(a)-1);
	r = real(z);x=imag(z)./s;
else
	r = r0;
	x = (r0.*B.*B)/6;
end