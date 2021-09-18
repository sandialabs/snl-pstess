function [Yrr,Yri,Yir,Yii] = dc_load(V,kdc)
% Syntax: [Yrr,Yri,Yir,Yii] = dc_load(V,kdc)
%
% Purpose: To calculate the non-linear Jacobian elements
%          associated with a line commutated hvdc link
%
% Inputs:  V - Equivalent HT terminal voltage
%          kdc dc time step indicator
%
% Outputs: Yrr - dir/dvr
%          Yri - dir/dvi
%          Yir - dii/dvr
%          Yii - dii/dvi
%
% Called by: nc_load

%-----------------------------------------------------------------------------%
% Version history
%
% Version: 1.0
% Date:    March 1997
% Author:  Graham Rogers
%-----------------------------------------------------------------------------%

global g;  % declaring struct of global variables

% note: Vdc is a file local variable that is not the same as g.dc.Vdc

V0(g.dc.r_idx,1) = g.dc.dcc_pot(:,7).*abs(V(g.dc.r_idx));
V0(g.dc.i_idx,1) = g.dc.dcc_pot(:,8).*abs(V(g.dc.i_idx));

dc_ang(g.dc.r_idx,1) = g.dc.alpha(:,kdc);
dc_ang(g.dc.i_idx,1) = g.dc.gamma(:,kdc);

Rc(g.dc.r_idx,1) = g.dc.dcc_pot(:,3);
Rc(g.dc.i_idx,1) = g.dc.dcc_pot(:,5);

idc(g.dc.r_idx,1) = g.dc.i_dcr(:,kdc);
idc(g.dc.i_idx,1) = g.dc.i_dci(:,kdc);

Vdc = V0.*cos(dc_ang) - idc.*Rc;
cphi = Vdc./V0;
sphi = sqrt(ones(g.dc.n_conv,1) - cphi.*cphi);

P = Vdc.*idc/g.sys.basmva;
Q = P.*sphi./cphi;
P(g.dc.i_idx) = -P(g.dc.i_idx);

iac = (P - 1j*Q)./conj(V);
ir = real(iac);
ii = imag(iac);

dV0dVr(g.dc.r_idx,1) = g.dc.dcc_pot(:,7).*real(V(g.dc.r_idx))./abs(V(g.dc.r_idx));
dV0dVr(g.dc.i_idx,1) = g.dc.dcc_pot(:,8).*real(V(g.dc.i_idx))./abs(V(g.dc.i_idx));
dV0dVi(g.dc.r_idx,1) = g.dc.dcc_pot(:,7).*imag(V(g.dc.r_idx))./abs(V(g.dc.r_idx));
dV0dVi(g.dc.i_idx,1) = g.dc.dcc_pot(:,8).*imag(V(g.dc.i_idx))./abs(V(g.dc.i_idx));

dPdVr = idc.*cos(dc_ang).*dV0dVr/g.sys.basmva;
dPdVi = idc.*cos(dc_ang).*dV0dVi/g.sys.basmva;

Kq = idc.*(ones(g.dc.n_conv,1)-cos(dc_ang).*cphi)./sphi/g.sys.basmva;
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
Yii = Yii + 2*(Q.*Vr - P.*Vi).*Vi./Vmag4;

end  % function end

% eof
