function [Yrr,Yri,Yir,Yii] = dc_load(V,k,kdc)
%Syntax: [Yrr,Yri,Yir,Yii] = dc_load(V,k,kdc)
%Purpose: To calculate the non-linear Jacobian elements
%         associated with a line commutated HVDC link
% Inputs:
%         V - Equivalent HT terminal voltage
%         k - step indicator; kdc dc time step indicator
% Outputs:
%         Yrr - dir/dvr
%         Yri - dir/dvi
%         Yir - dii/dvr
%         Yii - dii/dvi
% Called by: nc_load
%Version: 1.0
%Date:    March 1997
%Author:  Graham Rogers
% Copyright (c) Joe Chow 1991-1997 All Rights Reserved

global  i_dci  i_dcr  dcc_pot  alpha  gamma  basmva  r_idx  i_idx
global  n_conv n_dcl

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
iac = (P - jay*Q)./conj(V);
ir = real(iac);
ii = imag(iac);
dV0dVr(r_idx,1) = dcc_pot(:,7).*real(V(r_idx))./abs(V(r_idx));
dV0dVr(i_idx,1) = dcc_pot(:,8).*real(V(i_idx))./abs(V(i_idx));
dV0dVi(r_idx,1) = dcc_pot(:,7).*imag(V(r_idx))./abs(V(r_idx));
dV0dVi(i_idx,1) = dcc_pot(:,8).*imag(V(i_idx))./abs(V(i_idx));
dPdVr = idc.*cos(dc_ang).*dV0dVr/basmva;
dPdVi = idc.*cos(dc_ang).*dV0dVi/basmva;
Kq = idc.*(ones(n_conv,1)-cos(dc_ang).*cphi)./sphi/basmva;
dQdVr = Kq.*dV0dVr;
dQdVi = Kq.*dV0dVi;
Vr = real(V);
Vi = imag(V);
Vmag2 = Vr.*Vr + Vi.*Vi;
Vmag4 = Vmag2.*Vmag2;
Yrr = (P + Vi.*dQdVr + Vr.*dPdVr)./Vmag2;
Yrr = Yrr - 2*(P.*Vr + Q.*Vi).*Vr./Vmag4;
Yri = (Q + Vr.*dPdVi + Vi.*dQdVi)./Vmag2;
Yri = Yri - 2*(P.*Vr + Q.*Vi).*Vi./Vmag4;
Yir = -(Q - Vi.*dPdVr + Vr.*dQdVr)./Vmag2;
Yir = Yir + 2*(Q.*Vr - P.*Vi).*Vr./Vmag4;
Yii = (P + Vi.*dPdVi - Vr.*dQdVi)./Vmag2;
Yiii = Yii + 2*(Q.*Vr - P.*Vi).*Vi./Vmag4;




