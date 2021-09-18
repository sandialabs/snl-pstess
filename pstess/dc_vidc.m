function dc_vidc(k,kdc)
% syntax: dc_vidc(k,kdc)
% updates Vdc and i_dc assuming ac bus voltage remains constant

%-----------------------------------------------------------------------------%

global g;  % declaring struct of global variables

V0(g.dc.r_idx,1) = abs(g.bus.bus_v(g.dc.rec_ac_bus,k)).*g.dc.dcc_pot(:,7);
V0(g.dc.i_idx,1) = abs(g.bus.bus_v(g.dc.inv_ac_bus,k)).*g.dc.dcc_pot(:,8);

g.dc.Vdc(g.dc.r_idx,kdc) = V0(g.dc.r_idx,1).*cos(g.dc.alpha(:,kdc)) ...
                           - g.dc.i_dcr(:,kdc).*g.dc.dcc_pot(:,3);
g.dc.Vdc(g.dc.i_idx,kdc) = V0(g.dc.i_idx,1).*cos(g.dc.gamma(:,kdc)) ...
                           - g.dc.i_dci(:,kdc).*g.dc.dcc_pot(:,5);

g.dc.i_dc(g.dc.r_idx,kdc) = g.dc.i_dcr(:,kdc);
g.dc.i_dc(g.dc.i_idx,kdc) = g.dc.i_dci(:,kdc);

end  % function end

% eof
