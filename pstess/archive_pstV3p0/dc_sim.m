function [xr,dxr,xi,dxi] = ...
   dc_sim(k,kk,dcr,dci,xr,xi,bus_sim,hdc_sol)
global  basmva  bus_v
global  dcsp_con  dcl_con dcc_con 
global  r_idx  i_idx n_dcl  n_conv  ac_bus rec_ac_bus  inv_ac_bus
global  Vdc  i_dc  dcc_pot  alpha  gamma Vdc_ref
global  cur_ord dc_sig dcc_pot dcr_dsig dci_dsig ndcr_ud ndci_ud dcrd_sig dcid_sig
global  no_cap_idx  cap_idx  no_ind_idx  l_no_cap  l_cap

% States
%line
global i_dcr i_dci  v_dcc
global di_dcr  di_dci  dv_dcc  
%rectifier
global v_conr  
global dv_conr  
%inverter
global v_coni
global dv_coni 





% predictor
kdc = 10*(k-1)+ kk; jdc=kdc+1;
if n_conv~=0
   if ndcr_ud~=0
      tot_states=0;
      for jj = 1:ndcr_ud
         ydcrmx = dcr{jj,5};ydcrmn = dcr{jj,6};
         rec_num = dcr{jj,2};
         st_state = tot_states+1; dcr_states = dcr{jj,7}; tot_states = tot_states+dcr_states; 
         [dcr_dsig(rec_num,k),xr(st_state:tot_states,1),dxr(st_state:tot_states,1)] =...
            dcr_sud(jj,kdc,2,dcr{jj,1},dcrd_sig(jj,k),ydcrmx,ydcrmn,xr(st_state:tot_states,1));
      end
   else
      dxr(1,1)=0;
   end
   if ndci_ud~=0
      tot_states=0;
      for jj = 1:ndci_ud
         ydcimx = dci{jj,5};ydcimn = dci{jj,6};
         inv_num = dci{jj,2};
         st_state = tot_states+1; dci_states = dci{jj,7}; tot_states = tot_states+dci_states; 
         [dci_dsig(inv_num,k),xi(st_state:tot_states,1),dxi(st_state:tot_states,1)] =...
            dci_sud(jj,kdc,2,dci{jj,1},dcid_sig(jj,k),ydcimx,ydcimn,xi(st_state:tot_states,1));
      end
   else
      dxi(1,1)=0;
   end
   f = dc_cont(0,k,kdc,bus_sim,2);
   f = dc_line(0,k,kdc,bus_sim,2);
end
v_conr(:,jdc) = v_conr(:,kdc) + hdc_sol*dv_conr(:,kdc);
v_coni(:,jdc) = v_coni(:,kdc) + hdc_sol*dv_coni(:,kdc);
i_dcr(:,jdc) = i_dcr(:,kdc) + hdc_sol*di_dcr(:,kdc);
i_dci(:,jdc) = i_dci(:,kdc) + hdc_sol*di_dci(:,kdc);
v_dcc(:,jdc) = v_dcc(:,kdc) + hdc_sol*dv_dcc(:,kdc);
xr(:,2) = xr(:,1) + hdc_sol* dxr(:,1);
xi(:,2) = xi(:,1) + hdc_sol* dxi(:,1); 

f = dc_cont(0,k,jdc,bus_sim,1);% recalculate alpha and gamma
f=dc_vidc(k,kdc); % update Vdc and i_dc
if ndcr_ud~=0
   tot_states=0;
   for jj = 1:ndcr_ud
      ydcrmx = dcr{jj,5};ydcrmn = dcr{jj,6};
      rec_num = dcr{jj,2};
      st_state = tot_states+1; dcr_states = dcr{jj,7}; tot_states = tot_states+dcr_states; 
      [dcr_dsig(rec_num,k),xr(st_state:tot_states,2),dxr(st_state:tot_states,2)] =...
         dcr_sud(jj,jdc,2,dcr{jj,1},dcrd_sig(jj,k),ydcrmx,ydcrmn,xr(st_state:tot_states,2));
   end
else
   dxr(1,2)=0;
end
if ndci_ud~=0
   tot_states=0;
   for jj = 1:ndci_ud
      ydcimx = dci{jj,5};ydcimn = dci{jj,6};
      inv_num = dci{jj,2};
      st_state = tot_states+1; dci_states = dci{jj,7}; tot_states = tot_states+dci_states; 
      [dci_dsig(inv_num,j),xi(st_state:tot_states,2),dxi(st_state:tot_states,2)] =...
         dci_sud(jj,jdc,flag,dci{jj,1},dci_sig(jj,k),ydcimx,ydcimn,xi(st_state:tot_states,2));
   end
else
   dxi(1,2)=0;
end

f = dc_cont(0,k,jdc,bus_sim,2);
f = dc_line(0,k,jdc,bus_sim,2);

v_conr(:,jdc) = v_conr(:,kdc) + 0.5*hdc_sol*(dv_conr(:,kdc)+dv_conr(:,jdc));
v_coni(:,jdc) = v_coni(:,kdc) + 0.5*hdc_sol*(dv_coni(:,kdc)+dv_coni(:,jdc));
i_dcr(:,jdc) = i_dcr(:,kdc) + 0.5*hdc_sol*(di_dcr(:,kdc)+di_dcr(:,jdc));
i_dci(:,jdc) = i_dci(:,kdc) + 0.5*hdc_sol*(di_dci(:,kdc)+di_dci(:,jdc));
v_dcc(:,jdc) = v_dcc(:,kdc) + 0.5*hdc_sol*(dv_dcc(:,kdc)+dv_dcc(:,jdc));
xr(:,2) = xr(:,1) + 0.5*hdc_sol* (dxr(:,1)+dxr(:,2));
xi(:,2) = xi(:,1) + 0.5*hdc_sol* (dxi(:,1)+dxi(:,2)); 
f = dc_cont(0,k,jdc,bus_sim,1);%recalculate alpha and gamma
f=dc_vidc(k,jdc);% update Vdc and i_dc
