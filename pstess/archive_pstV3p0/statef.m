% frequency response
%syntax: [f,y,ym,ya]=statef(a,b,c,d,u,fstart,fstep,fend)
%calculates frequency response from a,b,c,d matrices
% a must be square matrix of order tstate
% b must be a matrix of size tstate by ni
% c must be a matrix of size no by tstate
% u is a matrix of input disturbances of size ni by nc
% d must be a matrix of size no by ni
% fstart is the start frequency in Hz
% fstep is the frequency step in Hz
% fend is the end frequency in Hz
%
% f is the frequency vector used for plotting
% ym is the magnitude vector
% ya is the phase vector in degrees
% y is the complex response
% Author: Graham Rogers
% Date:   September 1994
% Modified September 1996
%
function [f,y,ym,ya] = statef(a,b,c,d,u,fstart,fstep,fend)
i=sqrt(-1);
tstate = length(a(:,1));
f = fstart:fstep:fend;
n_fp = length(f);
ni = length(b(1,:));
no = length(c(:,1));
nc = length(u(1,:));
if no==1&ni==1
   y=zeros(1,n_fp);
   for k = 1:1:n_fp
      w = 2*pi*i*f(k)*diag(ones(tstate,1));
      x=inv(-a+w)*b*u;
      y(k)=c*x+d*u;
      ym(k) = abs(y(k));
      ya(k) =180.0*angle(y(k))/pi;
      if ya(k)>180
         ya(k)=ya(k)-360.0;
      end
      if ya(k)<-180
         ya(k)=ya(k)+360.0;
      end
   end
else
   y=zeros(no,nc,n_fp);
   ym = y;
   ya = y;
   for k = 1:1:n_fp
      w = 2*pi*i*f(k)*diag(ones(tstate,1));
      x=inv(-a+w)*b*u;
      y(:,:,k)=c*x+d*u;
      ym(:,:,k) = abs(y(:,:,k));
      ya(:,:,k) =180.0*angle(y(:,:,k))/pi;
      for kk = 1:no
         pm = find(ya(kk,:,k)> 180.0);
         if ~isempty(pm)
            ya(kk,pm,k)=ya(kk,pm,k)-360.0*ones(size(pm));
         end
         pm = find(ya(kk,:,k)<-180.0);
         if ~isempty(pm)
            ya(kk,pm,k)=ya(kk,pm,k)+360.0*ones(size(pm));
         end
      end
   end
end  
   