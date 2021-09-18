function [Ts,Td,f,ld,pd]=sd_torque(a,c_p,c_pm,mac_state,H,fbase)
%Syntax: {Ts,Td,ra,rw,f,ld,pd]=sd_torque(a,c_p,c_pm,mac_state,H,fbase)
%Outputs
%       Ts synchronizing torque as a function of frequency
%       Td damping torque as a function of frequency
%       f frequency vector
%       ld eigenvalues of non-angle or speed terms of a (ad)
%       pd participation factors associatedwith ld
%Inputs
%       a state matrix of interconnected generator system
%       c_p electrical power output matrix
%       c_pm mechanical power output matrix
%       mac_state indicator of machine states
%       H vector of generator inertia constants
%       fbase - base system frequency

crow = find(mac_state(:,2)~=1&mac_state(:,2)~=2);
ad = a(crow,crow);
cd=c_p-c_pm;
arow = find(mac_state(:,2)==1);
wrow = arow+1;
bda = a(crow,arow);
dda = cd(:,arow);
ddw = cd(:,wrow)+ 2*a(wrow,wrow)*diag(H)/(2*pi*fbase);

n_mac = length(arow);
bdw = a(crow,wrow)/2/pi/fbase;
[ud,ld]=eig(ad);
[ld,ld_idx]=sort(diag(ld));
ud = ud(:,ld_idx);
vd = inv(ud);
pd = ud.*vd.';
[f,yda]=statef(ad,bda,cd(:,crow),zeros(n_mac),eye(n_mac),0.01,0.01,5);
[f,ydw]=statef(ad,bdw,cd(:,crow),zeros(n_mac),eye(n_mac),0.01,0.01,5);
wf = 2*pi*f;

if n_mac>1
   Ts = zeros(n_mac,n_mac,500);
   Td=Ts;
   for k = 1:length(f)
      Ts(:,:,k) = dda+real(yda(:,:,k))-imag(ydw(:,:,k))*(wf(k));
      Td(:,:,k) = ddw+real(ydw(:,:,k))+imag(yda(:,:,k))*(1/wf(k));
   end
else
   Ts = zeros(1,500); Td = Ts;
   for k = 1:length(f)
      Ts(k) = dda+real(yda(k))-imag(ydw(k))*wf(k);
      Td(k) = ddw+real(ydw(k))+imag(ydw(k))/wf(k);
   end
end

