%perturb the pss variables
k_pss = find(mac_pss ==k);
if ~isempty(k_pss)
      disp(' disturb pss')
      j=j+1;
      pert = 0.0001*abs(pss1(k_pss,1));   
      pert = max(pert,0.0001);
      pss1(k_pss,2) = pss1(k_pss,1) + pert;
      p_file   % m file of perturbations
      st_name(k,j) = 12;
      j = j + 1;
      pert = 0.0001*abs(pss2(k_pss,1));   
      pert = max(pert,0.0001);
      pss2(k_pss,2) = pss2(k_pss,1) + pert;
      p_file   % m file of perturbations
      st_name(k,j) = 13;
      if ~isempty(pss_T4_idx)
         k_s_idx = find(pss_T4_idx==k_pss);
      else
         k_s_idx = [];
      end
      if ~isempty(k_s_idx)
        j = j+1;
        pert = 0.0001*abs(pss3(k_pss,1));   
        pert = max(pert,0.0001);
        pss3(k_pss,2) = pss3(k_pss,1) + pert;
        p_file   % m file of perturbations
        st_name(k,j) = 14;
      end
 end