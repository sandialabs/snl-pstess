%perturb the deltaP/omega filter variables
k_dpw = find(mac_dpw ==k);
if ~isempty(k_dpw)
   disp(' disturb deltaP/omega filter')
   j=j+1;
   pert = 0.0001*abs(sdpw1(k_dpw,1));   
   pert = max(pert,0.0001);
   sdpw1(k_dpw,2) = sdpw1(k_dpw,1) + pert;
   p_file   % m file of perturbations
   st_name(k,j) = 15;
   j=j+1;
   pert = 0.0001*abs(sdpw2(k_dpw,1));   
   pert = max(pert,0.0001);
   sdpw2(k_dpw,2) = sdpw2(k_dpw,1) + pert;
   p_file   % m file of perturbations
   st_name(k,j) = 16;
   j=j+1;
   pert = 0.0001*abs(sdpw3(k_dpw,1));   
   pert = max(pert,0.0001);
   sdpw3(k_dpw,2) = sdpw3(k_dpw,1) + pert;
   p_file   % m file of perturbations
   st_name(k,j) = 17;
   j=j+1;
   pert = 0.0001*abs(sdpw4(k_dpw,1));   
   pert = max(pert,0.0001);
   sdpw4(k_dpw,2) = sdpw4(k_dpw,1) + pert;
   p_file   % m file of perturbations
   st_name(k,j) = 18;
   j=j+1;
   pert = 0.0001*abs(sdpw5(k_dpw,1));   
   pert = max(pert,0.0001);
   sdpw5(k_dpw,2) = sdpw5(k_dpw,1) + pert;
   p_file   % m file of perturbations
   st_name(k,j) = 19;
   j=j+1;
   pert = 0.0001*abs(sdpw6(k_dpw,1));   
   pert = max(pert,0.0001);
   sdpw6(k_dpw,2) = sdpw6(k_dpw,1) + pert;
   p_file   % m file of perturbations
   st_name(k,j) = 20;
end