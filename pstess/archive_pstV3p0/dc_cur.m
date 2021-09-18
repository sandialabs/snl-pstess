function i_ac = dc_cur(V,k,kdc)
%Syntax: i_ac = dc_cur(V,k,kdc)
%
%Purpose : calculates the ac current injected into the network
%          as a function of the equivalent HT voltage 
%Input:    V the complex equivalent HT voltage
%          k time step indicator, kdc is the dc time step indicator
%Output:   i_ac the pu current injection at the HT bus
%Called by nc_load
%Version 1.0
%Date:   March 1997
%Copyright Joe Chow 1991-1997 All Rights Reserved

global  r_idx  i_idx  dcc_pot n_conv basmva
global  i_dcr  i_dci  alpha  gamma
jay = sqrt(-1);
V0(r_idx,1) = dcc_pot(:,7).*abs(V(r_idx));
V0(i_idx,1) = dcc_pot(:,8).*abs(V(i_idx));
dc_ang(r_idx,1) = alpha(:,kdc);
dc_ang(i_idx,1) = gamma(:,kdc);
Rc(r_idx,1) = dcc_pot(:,3);
Rc(i_idx,1) = dcc_pot(:,5);
idc(r_idx,1) = i_dcr(:,kdc);
idc(i_idx,1) = i_dci(:,kdc);
Vdc = V0.*cos(dc_ang) - idc.*Rc;
cphi = Vdc./V0;
sphi = sqrt(ones(n_conv,1) - cphi.*cphi);
P = Vdc.*idc/basmva;
Q = P.*sphi./cphi;
P(i_idx) = - P(i_idx);
i_ac = (P - jay*Q)./conj(V);
return
 