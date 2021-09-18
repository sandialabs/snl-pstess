% stabilizer frequency response
% calculates the phase lag through the exciter and generator
% the inverse of this phase lag is the ideal pss phase lead
% requires output from svm_smse
% generator data must be modified to have very high H (say 300)
% Author: Graham Rogers
% Date:   September 1994
function [f,ymag,yphase] = stabf(a,b,c,tstate)
i=sqrt(-1);
for k = 0:1:180
 w = pi*i*(0.2+k*.02)*diag(ones(tstate,1));
 f(k+1)=.1+k*.01;
 x=inv(-a(1:tstate,1:tstate)+w)*b(1:tstate,1);
 y(k+1)=c(1,1:tstate)*x;
 ymag(k+1) = abs(y(k+1));
 yphase(k+1) =180.0*angle(y(k+1))/pi-180.0;
 if yphase> 180.0;yphase=yphase-360.0;end
 if yphase<-180.0;yphase=yphase+360.0;end
end
%keyboard
return
