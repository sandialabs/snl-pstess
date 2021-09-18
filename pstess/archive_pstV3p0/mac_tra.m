function f = mac_tra(i,k,bus,flag)
% Syntax: f = mac_tra(i,k,bus,flag)
%
% Purpose: voltage-behind-transient-reactance generator
%          model, with vectorized computation option
%
% Input: i - generator number
%          - 0 for vectorized computation
%        k - integer time
%        bus - solved loadflow bus data
%        flag - 0 - initialization
%               1 - network interface computation
%               2 - generator dynamics computation 
%
% Output: f - dummy variable 
%
% Files:
%
% See Also: pst_var, mac_em, mac_sub

% Algorithm: 
%
% Calls:
%
% Call By:

% (c) Copyright 1991-1999 Joe H. Chow/Cherry Tree Scientific Software - All Rights Reserved

% History (in reverse chronological order)
%
% Version 2.0
% Date:   September 1997
% Saturation modified for low Eqprime

% Version:  1.0
% Author:   Joe H. Chow
% Date:     March 1991
% Modification: G.J. Rogers May 1995
% system variables
global  basmva basrad  mach_ref 
global  psi_re psi_im cur_re cur_im bus_int

% synchronous machine variables
global  mac_con mac_pot 
global  mac_ang mac_spd eqprime edprime 
global  curd curq curdg curqg fldcur
global  vex eterm theta ed eq 
global  pmech pelect qelect
global  dmac_ang dmac_spd deqprime dedprime 
global  n_tra 
global  mac_tra_idx 
global  pm_sig  
f = 0;

jay = sqrt(-1);
if n_tra~=0
 if flag == 0; % initialization         
  if i ~= 0
    % check that machine i is a transient machine
    tra = find(mac_tra_idx == i);
    if length(tra)~=0 
%       % make x'd=x'q
%       if mac_con(i,7)~=mac_con(i,12)
%          mac_con(i,12) = mac_con(i,7);
%       disp('changing xqp at generator');disp(i)
%       end
      % check Tqo'
      if mac_con(i,14)==0;mac_con(i,14)=999.0;end                                                       
      busnum = bus_int(mac_con(i,2)); % bus number 
      mac_pot(i,1) = basmva/mac_con(i,3); % scaled MVA base
      mac_pot(i,2) = 1.0; % base kv
      % extract bus information
      eterm(i,1) = bus(busnum,2);  % terminal bus voltage
      theta(busnum,1) = bus(busnum,3)*pi/180;  
                          % terminal bus angle in radians
      pelect(i,1) = bus(busnum,4)*mac_con(i,22);  
                        % electrical output power, active
      qelect(i,1) = bus(busnum,5)*mac_con(i,23);  
                        % electrical output power, reactive
      curr = sqrt(pelect(i,1)^2+qelect(i,1)^2) ...
              /eterm(i,1)*mac_pot(i,1);  % current magnitude
      phi = atan2(qelect(i,1),pelect(i,1)); 
                                        % power factor angle
      v = eterm(i,1)*exp(jay*theta(busnum,1)); 
                     % voltage in real and imaginary parts
                     % on system reference frame 
      curr = curr*exp(jay*(theta(busnum,1)-phi)); % current in real and 
                 % imaginary parts on system reference frame 
      eprime = v +(mac_con(i,5) + jay*mac_con(i,7))*curr; 
      ei = v + (mac_con(i,5) + jay*mac_con(i,11))*curr;
      mac_ang(i,1) = atan2(imag(ei),real(ei)); 
                                    % machine angle (delta)
      mac_spd(i,1) = 1; % machine speed at steady state
      rot = jay*exp(-jay*mac_ang(i,1)); % system reference frame rotation
      psi_re(i,1) = real(eprime);
      psi_im(i,1) = imag(eprime);
      eprime = eprime*rot;
      edprime(i,1) = real(eprime); 
      eqprime(i,1) = imag(eprime); 
      curr = curr*rot;
      mcurmag = abs(curr);
      pmech(i,1) = pelect(i,1)*mac_pot(i,1)...
                   + mac_con(i,5)*mcurmag*mcurmag;
      %mech power = elec power + losses on generator base
      curdg(i,1) = real(curr); curqg(i,1) = imag(curr);
      curd(i,1) = real(curr)/mac_pot(i,1); 
      curq(i,1) = imag(curr)/mac_pot(i,1);
      v = v*rot;
      ed(i,1) = real(v); eq(i,1) = imag(v);
      % compute saturation
      inv_sat = inv([0.64 0.8 1;1 1 1;1.44 1.2 1]);
      b = [0.8 1+mac_con(i,20) 1.2*(1+mac_con(i,21))];
      mac_pot(i,3) = b*inv_sat(1,:)';
      mac_pot(i,4) = b*inv_sat(2,:)';
      mac_pot(i,5) = b*inv_sat(3,:)';
      E_Isat = mac_pot(i,3)*eqprime(i,1)^2 ...
            + mac_pot(i,4)*eqprime(i,1) + mac_pot(i,5); 
      if eqprime(i,1)<0.8;E_Isat=eqprime(i,1);end
      vex(i,1) = E_Isat + (mac_con(i,6)- ...
               mac_con(i,7))*curdg(i,1);
      fldcur(i,1) = vex(i,1);  % JHC fixed typo 2011 0511, per DFK
    end
  else
    % vectorized computation
    % make xd' = xq'
%     uexp_idx = find(mac_con(mac_tra_idx,7)~=mac_con(mac_tra_idx,12));
%     if ~isempty(uexp_idx)
%        mac_con(mac_tra_idx(uexp_idx),12) = mac_con(mac_tra_idx(uexp_idx),7);
%        mtist=int2str(mac_tra_idx(uexp_idx));
%       % disp(['changing xqp at generators  ' mtist])
%        disp('changing xqp at generators') % The disp command revised by 
%        disp(mtist) % Joe Chow 12/12/2015 
%                    % incompatible dimension in concaternation 
%     end
    % make Tqo' non-zero
    notqp_idx = find(mac_con(mac_tra_idx,14)==0);
    if ~isempty(notqp_idx)
      mac_con(mac_tra_idx(notqp_idx),14) = 999.0*ones(length(notqp_idx),1);
    end
    busnum = bus_int(mac_con(mac_tra_idx,2)); % bus number 
    mac_pot(mac_tra_idx,1) = basmva*ones(n_tra,1)./mac_con(mac_tra_idx,3); 
                          % scaled MVA base
    mac_pot(mac_tra_idx,2) = ones(n_tra,1); % base kv
    % extract bus information
    eterm(mac_tra_idx,1) = bus(busnum,2);  % terminal bus voltage
    theta(busnum,1) = bus(busnum,3)*pi/180;  
                          % terminal bus angle in radians
    pelect(mac_tra_idx,1) = bus(busnum,4).*mac_con(mac_tra_idx,22);  
                        % electrical output power, active
    qelect(mac_tra_idx,1) = bus(busnum,5).*mac_con(mac_tra_idx,23);  
                        % electrical output power, reactive
    curr = sqrt(pelect(mac_tra_idx,1).^2+qelect(mac_tra_idx,1).^2) ...
            ./eterm(mac_tra_idx,1).*mac_pot(mac_tra_idx,1);  % current magnitude
    phi = atan2(qelect(mac_tra_idx,1),pelect(mac_tra_idx,1)); 
                                        % power factor angle
    v = eterm(mac_tra_idx,1).*exp(jay*theta(busnum,1)); 
                     % voltage in real and imaginary parts
                     % on system reference frame 
    curr = curr.*exp(jay*(theta(busnum,1)-phi)); % current in real and 
                                                % imaginary parts on system reference frame 
    eprime = v + (mac_con(mac_tra_idx,5)+jay*mac_con(mac_tra_idx,7)).*curr; 
    ei = v + (mac_con(mac_tra_idx,5)+jay*mac_con(mac_tra_idx,11)).*curr;
    mac_ang(mac_tra_idx,1) = atan2(imag(ei),real(ei)); 
                                    % machine angle (delta)
    mac_spd(mac_tra_idx,1) = ones(n_tra,1); 
                            % machine speed at steady state
    rot = jay*exp(-jay*mac_ang(mac_tra_idx,1)); % system reference frame rotation
    psi_re(mac_tra_idx,1)=real(eprime);
    psi_im(mac_tra_idx,1)=imag(eprime);
    eprime = eprime.*rot;
    edprime(mac_tra_idx,1) = real(eprime); 
    eqprime(mac_tra_idx,1) = imag(eprime); 
    curr = curr.*rot;
    mcurmag = abs(curr);
    pmech(mac_tra_idx,1) = pelect(mac_tra_idx,1).*mac_pot(mac_tra_idx,1)...
                            + mac_con(mac_tra_idx,5).*mcurmag.*mcurmag;
    %pmech = pelec + losses om generator base
    curdg(mac_tra_idx,1) = real(curr); 
    curqg(mac_tra_idx,1) = imag(curr);
    curd(mac_tra_idx,1) = real(curr)./mac_pot(mac_tra_idx,1); 
    curq(mac_tra_idx,1) = imag(curr)./mac_pot(mac_tra_idx,1);
    v = v.*rot;
    ed(mac_tra_idx,1) = real(v); 
    eq(mac_tra_idx,1) = imag(v);
    % compute saturation
    inv_sat = inv([0.64 0.8 1;1 1 1;1.44 1.2 1]);
    b = [0.8*ones(n_tra,1) ones(n_tra,1)+mac_con(mac_tra_idx,20) ...
           1.2*(ones(n_tra,1)+mac_con(mac_tra_idx,21))];
    mac_pot(mac_tra_idx,3) = b*inv_sat(1,:)';
    mac_pot(mac_tra_idx,4) = b*inv_sat(2,:)';
    mac_pot(mac_tra_idx,5) = b*inv_sat(3,:)';
    E_Isat = mac_pot(mac_tra_idx,3).*eqprime(mac_tra_idx,1).^2 ...
            + mac_pot(mac_tra_idx,4).*eqprime(mac_tra_idx,1)...
            + mac_pot(mac_tra_idx,5); 
    nosat_idx=find(eqprime(mac_tra_idx,1)<.8);
    if ~isempty(nosat_idx)
       E_Isat(nosat_idx)=eqprime(mac_tra_idx(nosat_idx),1);
    end
    vex(mac_tra_idx,1) = E_Isat + (mac_con(mac_tra_idx,6)- ...
               mac_con(mac_tra_idx,7)).*curdg(mac_tra_idx,1);
    fldcur(mac_tra_idx,1) = vex(mac_tra_idx,1);

  end
 % end initialization
 end
 if flag == 1 % network interface computation 
  if i ~= 0
     % check i for transient machine
     tra = find(mac_tra_idx ==i,1);
     if ~isempty(tra)
      mac_ang(i,k) = mac_ang(i,k) - mach_ref(k);   
                     % wrt machine reference
      psi_re(i,k) = sin(mac_ang(i,k))*edprime(i,k) + ...
         cos(mac_ang(i,k))*eqprime(i,k); % real part of psi
      psi_im(i,k) = -cos(mac_ang(i,k))*edprime(i,k) + ...
         sin(mac_ang(i,k))*eqprime(i,k); % imag part of psi
     end
   else
    % vectorized computation
    mac_ang(mac_tra_idx,k) = mac_ang(mac_tra_idx,k)-mach_ref(k)*ones(n_tra,1);
                     % wrt machine reference
    psi_re(mac_tra_idx,k) = sin(mac_ang(mac_tra_idx,k)).*edprime(mac_tra_idx,k) + ...
                            cos(mac_ang(mac_tra_idx,k)).*eqprime(mac_tra_idx,k); 
                            % real part of psi
    psi_im(mac_tra_idx,k) = -cos(mac_ang(mac_tra_idx,k)).*edprime(mac_tra_idx,k) + ...
                            sin(mac_ang(mac_tra_idx,k)).*eqprime(mac_tra_idx,k);
                            % imag part of psi

  end
 % end of interface calculation
 end
 
 if flag == 2 || flag == 3 % generator dynamics calculation
  if i ~= 0
    %check that i is a transient machine
    tra = find(mac_tra_idx ==i);
    if length(tra)~=0
      curd(i,k) = sin(mac_ang(i,k))*cur_re(i,k) - ...
            cos(mac_ang(i,k))*cur_im(i,k); % d-axis current
      curq(i,k) = cos(mac_ang(i,k))*cur_re(i,k) + ...
            sin(mac_ang(i,k))*cur_im(i,k); % q-axis current
      curdg(i,k) = curd(i,k)*mac_pot(i,1);
      curqg(i,k) = curq(i,k)*mac_pot(i,1);
      E_Isat = mac_pot(i,3)*eqprime(i,k)^2 ...
            + mac_pot(i,4)*eqprime(i,k) + mac_pot(i,5); 
      if eqprime(i,1)<0.8;E_Isat=eqprime(i,1);end
      dedprime(i,k) = (-edprime(i,k) + (mac_con(i,11)-...
                   mac_con(i,12))*curqg(i,k))/mac_con(i,14);
      fldcur(i,k) = E_Isat + (mac_con(i,6)-mac_con(i,7))*curdg(i,k);
      deqprime(i,k) = (vex(i,k) - fldcur(i,k))/mac_con(i,9);
      ed(i,k) = edprime(i,k) - mac_con(i,5)*curdg(i,k) + mac_con(i,7)*curqg(i,k);
      eq(i,k) = eqprime(i,k)- mac_con(i,5)*curqg(i,k) - mac_con(i,7)*curdg(i,k);
      eterm(i,k) = sqrt(ed(i,k)^2+eq(i,k)^2);
      pelect(i,k) = eq(i,k)*curq(i,k) + ed(i,k)*curd(i,k);
      qelect(i,k) = eq(i,k)*curd(i,k) - ed(i,k)*curq(i,k);
      curmag = abs(curdg(i,k) + jay*curqg(i,k));
      Te = pelect(i,k)*mac_pot(i,1) + mac_con(i,5)*curmag*curmag;
      dmac_ang(i,k) = basrad*(mac_spd(i,k)-1.);
      dmac_spd(i,k) = (pmech(i,k) + pm_sig(i,k) - Te ...   % JHC add missing pm_sig term 2011 0511, per DKF
         -mac_con(i,17)*(mac_spd(i,k)-1))/(2*mac_con(i,16)); 
    end
  else
    % vectorized computation
   
    curd(mac_tra_idx,k) = sin(mac_ang(mac_tra_idx,k)).*cur_re(mac_tra_idx,k) - ...
            cos(mac_ang(mac_tra_idx,k)).*cur_im(mac_tra_idx,k); % d-axis current
    curq(mac_tra_idx,k) = cos(mac_ang(mac_tra_idx,k)).*cur_re(mac_tra_idx,k) + ...
            sin(mac_ang(mac_tra_idx,k)).*cur_im(mac_tra_idx,k); % q-axis current
    curdg(mac_tra_idx,k) = curd(mac_tra_idx,k).*mac_pot(mac_tra_idx,1);
    curqg(mac_tra_idx,k) = curq(mac_tra_idx,k).*mac_pot(mac_tra_idx,1);
    E_Isat = mac_pot(mac_tra_idx,3).*eqprime(mac_tra_idx,k).^2 ...
            + mac_pot(mac_tra_idx,4).*eqprime(mac_tra_idx,k) + mac_pot(mac_tra_idx,5); 
    nosat_idx=find(eqprime(mac_tra_idx,1)<.8);
    if ~isempty(nosat_idx)
         E_Isat(nosat_idx)=eqprime(mac_tra_idx(nosat_idx),k);
    end
    fldcur(mac_tra_idx,k) = E_Isat + (mac_con(mac_tra_idx,6)-mac_con(mac_tra_idx,7))...
                            .*curdg(mac_tra_idx,k);
    deqprime(mac_tra_idx,k) = (vex(mac_tra_idx,k) - fldcur(mac_tra_idx,k))./mac_con(mac_tra_idx,9);
    dedprime(mac_tra_idx,k) = (-edprime(mac_tra_idx,k) +...
                              (mac_con(mac_tra_idx,11)-mac_con(mac_tra_idx,12))...
                              .*curqg(mac_tra_idx,k))./mac_con(mac_tra_idx,14);
    ed(mac_tra_idx,k) = edprime(mac_tra_idx,k) - mac_con(mac_tra_idx,5).*curdg(mac_tra_idx,k)...
                        + mac_con(mac_tra_idx,7).*curqg(mac_tra_idx,k);
    eq(mac_tra_idx,k) = eqprime(mac_tra_idx,k)- mac_con(mac_tra_idx,5).*curqg(mac_tra_idx,k)...
                        - mac_con(mac_tra_idx,7).*curdg(mac_tra_idx,k);
    eterm(mac_tra_idx,k) = sqrt(ed(mac_tra_idx,k).^2+eq(mac_tra_idx,k).^2);
    pelect(mac_tra_idx,k) = eq(mac_tra_idx,k).*curq(mac_tra_idx,k)...
                            + ed(mac_tra_idx,k).*curd(mac_tra_idx,k);
    qelect(mac_tra_idx,k) = eq(mac_tra_idx,k).*curd(mac_tra_idx,k)...
                            - ed(mac_tra_idx,k).*curq(mac_tra_idx,k);
    curmag = abs(curdg(mac_tra_idx,k) + jay*curqg(mac_tra_idx,k));
    Te = pelect(mac_tra_idx,k).*mac_pot(mac_tra_idx,1) + mac_con(mac_tra_idx,5).*curmag.*curmag;
    dmac_ang(mac_tra_idx,k) = basrad*(mac_spd(mac_tra_idx,k)-ones(n_tra,1));
    dmac_spd(mac_tra_idx,k) = (pmech(mac_tra_idx,k)+ pm_sig(mac_tra_idx,k) - Te ...
         -mac_con(mac_tra_idx,17).*(mac_spd(mac_tra_idx,k)-ones(n_tra,1)))./(2*mac_con(mac_tra_idx,16));
  end

 % end calculation of rates of change
 end
end
