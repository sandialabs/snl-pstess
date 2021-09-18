function [xr,dxr,xi,dxi] = dc_sim(k,kk,dcr,dci,dcrd_sig,dcid_sig,xr,xi,bus_sim,hdc_sol)

%-----------------------------------------------------------------------------%

global g;  % declaring struct of global variables

% predictor
kdc = 10*(k-1) + kk;
jdc = kdc + 1;

if (g.dc.n_conv ~= 0)
    if (g.dc.ndcr_ud ~= 0)
        tot_states = 0;
        for jj = 1:g.dc.ndcr_ud
            ydcrmx = dcr{jj,5};
            ydcrmn = dcr{jj,6};
            rec_num = dcr{jj,2};
            st_state = tot_states + 1;
            dcr_states = dcr{jj,7};
            tot_states = tot_states + dcr_states;
            [g.dc.dcr_dsig(rec_num,k),xr(st_state:tot_states,1),dxr(st_state:tot_states,1)] = ...
                dcr_sud(jj,kdc,2,dcr{jj,1},dcrd_sig(jj,k),ydcrmx,ydcrmn,xr(st_state:tot_states,1));
        end
    else
        dxr(1,1) = 0;
    end

    if (g.dc.ndci_ud ~= 0)
        tot_states = 0;
        for jj = 1:g.dc.ndci_ud
            ydcimx = dci{jj,5};
            ydcimn = dci{jj,6};
            inv_num = dci{jj,2};
            st_state = tot_states + 1;
            dci_states = dci{jj,7};
            tot_states = tot_states + dci_states;
            [g.dc.dci_dsig(inv_num,k),xi(st_state:tot_states,1),dxi(st_state:tot_states,1)] = ...
                dci_sud(jj,kdc,2,dci{jj,1},dcid_sig(jj,k),ydcimx,ydcimn,xi(st_state:tot_states,1));
        end
    else
        dxi(1,1) = 0;
    end

    dc_cont(0,k,kdc,bus_sim,2);
    dc_line(0,kdc,2);

end

g.dc.v_conr(:,jdc) = g.dc.v_conr(:,kdc) + hdc_sol*g.dc.dv_conr(:,kdc);
g.dc.v_coni(:,jdc) = g.dc.v_coni(:,kdc) + hdc_sol*g.dc.dv_coni(:,kdc);
g.dc.i_dcr(:,jdc) = g.dc.i_dcr(:,kdc) + hdc_sol*g.dc.di_dcr(:,kdc);
g.dc.i_dci(:,jdc) = g.dc.i_dci(:,kdc) + hdc_sol*g.dc.di_dci(:,kdc);
g.dc.v_dcc(:,jdc) = g.dc.v_dcc(:,kdc) + hdc_sol*g.dc.dv_dcc(:,kdc);
xr(:,2) = xr(:,1) + hdc_sol*dxr(:,1);
xi(:,2) = xi(:,1) + hdc_sol*dxi(:,1);

dc_cont(0,k,jdc,bus_sim,1);  % recalculate alpha and gamma
dc_vidc(k,kdc);              % update Vdc and i_dc

if (g.dc.ndcr_ud ~= 0)
    tot_states = 0;
    for jj = 1:g.dc.ndcr_ud
        ydcrmx = dcr{jj,5};
        ydcrmn = dcr{jj,6};
        rec_num = dcr{jj,2};
        st_state = tot_states + 1;
        dcr_states = dcr{jj,7};
        tot_states = tot_states + dcr_states;
        [g.dc.dcr_dsig(rec_num,k),xr(st_state:tot_states,2),dxr(st_state:tot_states,2)] = ...
            dcr_sud(jj,jdc,2,dcr{jj,1},dcrd_sig(jj,k),ydcrmx,ydcrmn,xr(st_state:tot_states,2));
    end
else
    dxr(1,2) = 0;
end

if (g.dc.ndci_ud ~= 0)
    tot_states = 0;
    for jj = 1:g.dc.ndci_ud
        ydcimx = dci{jj,5};
        ydcimn = dci{jj,6};
        inv_num = dci{jj,2};
        st_state = tot_states + 1;
        dci_states = dci{jj,7};
        tot_states = tot_states + dci_states;
        [g.dc.dci_dsig(inv_num,j),xi(st_state:tot_states,2),dxi(st_state:tot_states,2)] = ...
            dci_sud(jj,jdc,flag,dci{jj,1},dcid_sig(jj,k),ydcimx,ydcimn,xi(st_state:tot_states,2));
    end
else
    dxi(1,2) = 0;
end

dc_cont(0,k,jdc,bus_sim,2);
dc_line(0,jdc,2);

g.dc.v_conr(:,jdc) = g.dc.v_conr(:,kdc) ...
                     + 0.5*hdc_sol*(g.dc.dv_conr(:,kdc) + g.dc.dv_conr(:,jdc));
g.dc.v_coni(:,jdc) = g.dc.v_coni(:,kdc) ...
                     + 0.5*hdc_sol*(g.dc.dv_coni(:,kdc) + g.dc.dv_coni(:,jdc));
g.dc.i_dcr(:,jdc) = g.dc.i_dcr(:,kdc) ...
                    + 0.5*hdc_sol*(g.dc.di_dcr(:,kdc) + g.dc.di_dcr(:,jdc));
g.dc.i_dci(:,jdc) = g.dc.i_dci(:,kdc) ...
                    + 0.5*hdc_sol*(g.dc.di_dci(:,kdc) + g.dc.di_dci(:,jdc));
g.dc.v_dcc(:,jdc) = g.dc.v_dcc(:,kdc) ...
                    + 0.5*hdc_sol*(g.dc.dv_dcc(:,kdc) + g.dc.dv_dcc(:,jdc));

xr(:,2) = xr(:,1) + 0.5*hdc_sol*(dxr(:,1) + dxr(:,2));
xi(:,2) = xi(:,1) + 0.5*hdc_sol*(dxi(:,1) + dxi(:,2));

dc_cont(0,k,jdc,bus_sim,1);  % recalculate alpha and gamma
dc_vidc(k,jdc);              % update Vdc and i_dc

end  % function end

% eof
