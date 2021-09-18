 function [f]=ind_ldto(i,k)
% syntax f=ind_ldto(i,k)
% f is a dummy variable
% i is the motor number
% k is the time step
% if i is set to zero, vector computation is invoked
% Purpose
% Template for Induction Motor Load Torque Calculation as 
% a function of slip
% format for motor load data - mld_con
% 1 motor number
% 2 bus number
% 3 stiction load pu on motor base (f1)
% 4 stiction load coefficient (i1)
% 5 external load  pu on motor base(f2)
% 6 external load coefficient (i2)
% 
% load has the form
% tload = f1*slip^i1 + f2*(1-slip)^i2
% Author Graham Rogers
% Date November 1995
% (c) copyright Joe Chow 1995-1996 All rights reserved

f=0;
%set global variables
global tload slip mld_con

if i == 0  % vector computation
  numot = length(mld_con(:,2));
  tload(:,k)=mld_con(:,3).*(slip(:,k).^mld_con(:,4)) ...
             +mld_con(:,5).*(ones(numot,1)-slip(:,k)).^mld_con(:,6);
else
  tload(i,k)=mld_con(i,3)*(slip(i,k)^mld_con(i,4)) ...
             + mld_con(i,5)*(1-slip(i,k))^mld_con(i,6);
end

