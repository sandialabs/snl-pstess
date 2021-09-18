function f=dc_vidc(k,kdc)
%syntax: f=dc_vidc(k,kdc)
% updates Vdc and i_dc assuming ac bus voltage remains constant
global bus_v 
global i_dc Vdc alpha gamma dcc_pot i_dcr  i_dci 
global r_idx i_idx ac_bus rec_ac_bus inv_ac_bus n_conv
f=0;
V0(r_idx,1) = abs(bus_v(rec_ac_bus,k)).*dcc_pot(:,7);
V0(i_idx,1) = abs(bus_v(inv_ac_bus,k)).*dcc_pot(:,8);
Vdc(r_idx,kdc) = V0(r_idx,1).*cos(alpha(:,kdc)) - i_dcr(:,kdc).*dcc_pot(:,3);
Vdc(i_idx,kdc) = V0(i_idx,1).*cos(gamma(:,kdc)) - i_dci(:,kdc).*dcc_pot(:,5);
i_dc(r_idx,kdc) = i_dcr(:,kdc);
i_dc(i_idx,kdc) = i_dci(:,kdc);
