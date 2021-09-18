function [f] = mac_em(i,k,bus,flag)
% Syntax: [f] = mac_em(i,k,bus,flag)
%
% Purpose: generator electromechanical model, with
%            vectorized computation option
%          state variables are: mac_ang, mac_spd
% 
% Input: i - generator number
%          - 0 for vectorized computation
%        k - integer time
%        bus - solved loadflow bus data
%        flag - 0 - initialization
%               1 - network interface computation
%               2 - generator dynamics computation and state
%                   state matrix building
%
% Output: f - dummy variable 
%
% Files:
%
% See Also:

% Algorithm: 
%
% Calls:
%
% Call By:

% (c) Copyright 1991 Joe H. Chow - All Rights Reserved

% History (in reverse chronological order)
%
% Version: 2.0
% Date:    June 1996
% Author:  Graham Rogers
% Purpose: add facility to allow different machine models in vector run
% Modification:

% Version:  1.0
% Author:   Joe H. Chow
% Date:     January 1991

% system variables
global  basmva basrad syn_ref mach_ref sys_freq
global  bus_v bus_ang psi_re psi_im cur_re cur_im bus_int

% synchronous machine variables
global  mac_con mac_pot mac_int
global  mac_ang mac_spd eqprime edprime 
global  curd curq curdg curqg fldcur
global  psidpp psiqpp vex eterm theta ed eq 
global  pmech pelect qelect
global  dmac_ang dmac_spd deqprime dedprime 
global  n_mac n_em n_tra n_sub
global  mac_em_idx mac_tra_idx mac_sub_idx
global  pm_sig n_pm 
f = 0;
if n_em ~=0
 jay = sqrt(-1);
 if flag == 0; % initialization
   if i~=0
      % non-vector calculation
      % check for em model
      em = find(mac_em_idx==i)
      if length(em)~=0
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
                                         % on generator base
        phi = atan2(qelect(i,1),pelect(i,1)); 
                                        % power factor angle
        v = eterm(i,1)*exp(jay*theta(busnum,1)); 
                     % voltage in real and imaginary parts
                     % on system reference frame 
        curr = curr*exp(jay*(theta(busnum,1)-phi)); % complex current  
                                                    % in system reference frame 
        eprime = v + jay*mac_con(i,7)*curr; 
        % voltage behind transient reactance
        ei = eprime;
        mac_ang(i,1) = atan2(imag(ei),real(ei)); 
                                    % machine angle (delta)
        mac_spd(i,1) = 1; % machine speed at steady state
        rot = jay*exp(-jay*mac_ang(i,1)); 
                          % system reference frame rotation
        psi_re(i,1) = real(eprime);
        psi_im(i,1) = imag(eprime);
        eprime = eprime*rot;
        edprime(i,1) = real(eprime); 
        eqprime(i,1) = imag(eprime); 
        curr = curr*rot;%current on Park's frame
        curdg(i,1) = real(curr); curqg(i,1) = imag(curr);
        curd(i,1) = real(curr)/mac_pot(i,1); 
        curq(i,1) = imag(curr)/mac_pot(i,1);
        % convert to system base
        v = v*rot;
        ed(i,1) = real(v); eq(i,1) = imag(v);% in Park's frame
        vex(i,1) = eqprime(i,1);
        pmech(i,1) = pelect(i,1)*mac_pot(i,1); % set input
         % mechanical power equal to electrical output power
         % since losses are zero for em model. On generator base
      end
   else
      % vectorized computation
      busnum = bus_int(mac_con(mac_em_idx,2)); % bus number 
      mac_pot(mac_em_idx,1) = basmva*ones(n_em,1)./mac_con(mac_em_idx,3); 
                          % scaled MVA base
      mac_pot(mac_em_idx,2) = 1.0*ones(n_em,1); % base kv
      % extract bus information
      eterm(mac_em_idx,1) = bus(busnum,2);  % terminal bus voltage
      theta(busnum,1) = bus(busnum,3)*pi/180;  
                          % terminal bus angle in radians
      pelect(mac_em_idx,1) = bus(busnum,4).*mac_con(mac_em_idx,22);  
                        % electrical output power, active
      qelect(mac_em_idx,1) = bus(busnum,5).*mac_con(mac_em_idx,23);  
                        % electrical output power, reactive
      curr = sqrt(pelect(mac_em_idx,1).^2+qelect(mac_em_idx,1).^2)...
            ./eterm(mac_em_idx,1).*mac_pot(mac_em_idx,1);  % current magnitude
                                                           % on generator base
      phi = atan2(qelect(mac_em_idx,1),pelect(mac_em_idx,1)); 
                                        % power factor angle
      v = eterm(mac_em_idx,1).*exp(jay*theta(busnum,1)); 
                     % voltage in real and imaginary parts
                     % on system reference frame 
      curr = curr.*exp(jay*(theta(busnum,1)-phi)); % current in real and 
                 % imaginary parts on system reference frame 
      eprime = v + jay*mac_con(mac_em_idx,7).*curr; 
      ei = eprime;
      mac_ang(mac_em_idx,1) = atan2(imag(ei),real(ei)); 
                                    % machine angle (delta)
      mac_spd(mac_em_idx,1) = ones(n_em,1); 
                            % machine speed at steady state
      rot = jay*exp(-jay*mac_ang(mac_em_idx,1)); 
                          % system reference frame rotation
      psi_re(mac_em_idx,1) = real(eprime);
      psi_im(mac_em_idx,1) = imag(eprime);
      eprime = eprime.*rot;% in Park's frame
      edprime(mac_em_idx,1) = real(eprime); 
      eqprime(mac_em_idx,1) = imag(eprime);
      curr = curr.*rot;%in Park's frame
      curdg(mac_em_idx,1) = real(curr); 
      curqg(mac_em_idx,1) = imag(curr);
      curd(mac_em_idx,1) = real(curr)./mac_pot(mac_em_idx,1); 
      curq(mac_em_idx,1) = imag(curr)./mac_pot(mac_em_idx,1);
      v = v.*rot;%in Park's frame
      ed(mac_em_idx,1) = real(v); 
      eq(mac_em_idx,1) = imag(v);
      vex(mac_em_idx,1) = eqprime(mac_em_idx,1);
      pmech(mac_em_idx,1) = pelect(mac_em_idx,1).*mac_pot(mac_em_idx,1); % set input
         % mechanical power equal to electrical output power
         % since losses are zero in em model. On generator base
   end
   %end initialization
 end

 if flag == 1 % network interface computation 
  if i ~= 0
      % check for em machine
      em = find(mac_em_idx==i);
      if length(em)~=0
        mac_ang(i,k) = mac_ang(i,k) - mach_ref(k);   
                     % wrt machine reference
        psi_re(i,k) = sin(mac_ang(i,k))*edprime(i,k) + ...
                      cos(mac_ang(i,k))*eqprime(i,k); % real part of psi
        psi_im(i,k) = -cos(mac_ang(i,k))*edprime(i,k) + ...
                      sin(mac_ang(i,k))*eqprime(i,k); % imag part of psi
      end
  else
      % vectorized computation
      mac_ang(mac_em_idx,k) = mac_ang(mac_em_idx,k)-mach_ref(k)*ones(n_em,1);
                     % wrt machine reference
      psi_re(mac_em_idx,k) = sin(mac_ang(mac_em_idx,k)).*edprime(mac_em_idx,k) + ...
         cos(mac_ang(mac_em_idx,k)).*eqprime(mac_em_idx,k); % real part of psi
      psi_im(mac_em_idx,k) = -cos(mac_ang(mac_em_idx,k)).*edprime(mac_em_idx,k) + ...
         sin(mac_ang(mac_em_idx,k)).*eqprime(mac_em_idx,k); % imag part of psi
  end
  % end interface
 end
 if flag == 2 % generator dynamics calculation
  if i ~= 0
    % check for em machine
    em = find(mac_em_idx==i);
    if length(em)~=0
      curd(i,k) = sin(mac_ang(i,k))*cur_re(i,k) - ...
            cos(mac_ang(i,k))*cur_im(i,k); % d-axis current
      curq(i,k) = cos(mac_ang(i,k))*cur_re(i,k) + ...
            sin(mac_ang(i,k))*cur_im(i,k); % q-axis current
      curdg(i,k) = curd(i,k)*mac_pot(i,1);
      curqg(i,k) = curq(i,k)*mac_pot(i,1);
      dedprime(i,k) = 0;
      deqprime(i,k) = 0;
      ed(i,k) = edprime(i,k) + mac_con(i,7)*curqg(i,k);
      eq(i,k) = eqprime(i,k) - mac_con(i,7)*curdg(i,k);
      eterm(i,k) = sqrt(ed(i,k)^2+eq(i,k)^2);
      pelect(i,k) = eq(i,k)*curq(i,k) + ed(i,k)*curd(i,k);
      qelect(i,k) = eq(i,k)*curd(i,k) - ed(i,k)*curq(i,k);
      dmac_ang(i,k) = basrad*(mac_spd(i,k)-1.);
      dmac_spd(i,k) = (pmech(i,k)-pelect(i,k)*mac_pot(i,1)... 
        -mac_con(i,17)*(mac_spd(i,k)-1))/(2.*mac_con(i,16));
    end
  else
      % vectorized computation
      curd(mac_em_idx,k) = sin(mac_ang(mac_em_idx,k)).*cur_re(mac_em_idx,k) - ...
            cos(mac_ang(mac_em_idx,k)).*cur_im(mac_em_idx,k); % d-axis current
      curq(mac_em_idx,k) = cos(mac_ang(mac_em_idx,k)).*cur_re(mac_em_idx,k) + ...
            sin(mac_ang(mac_em_idx,k)).*cur_im(mac_em_idx,k); % q-axis current
      curdg(mac_em_idx,k) = curd(mac_em_idx,k).*mac_pot(mac_em_idx,1);
      curqg(mac_em_idx,k) = curq(mac_em_idx,k).*mac_pot(mac_em_idx,1);
      dedprime(mac_em_idx,k) = zeros(n_em,1);
      deqprime(mac_em_idx,k) = zeros(n_em,1);
      ed(mac_em_idx,k) = edprime(mac_em_idx,k)...
                        + mac_con(mac_em_idx,7).*curqg(mac_em_idx,k);
      eq(mac_em_idx,k) = eqprime(mac_em_idx,k)...
                        - mac_con(mac_em_idx,7).*curdg(mac_em_idx,k);
      eterm(mac_em_idx,k) = sqrt(ed(mac_em_idx,k).^2+eq(mac_em_idx,k).^2);
      pelect(mac_em_idx,k) = eq(mac_em_idx,k).*curq(mac_em_idx,k)...
                           + ed(mac_em_idx,k).*curd(mac_em_idx,k);
      qelect(mac_em_idx,k) = eq(mac_em_idx,k).*curd(mac_em_idx,k)...
                           - ed(mac_em_idx,k).*curq(mac_em_idx,k);
      dmac_ang(mac_em_idx,k) = basrad*(mac_spd(mac_em_idx,k)-ones(n_em,1));
      dmac_spd(mac_em_idx,k) =(pmech(mac_em_idx,k)+ pm_sig(mac_em_idx,k) ...
                      -pelect(mac_em_idx,k).*mac_pot(mac_em_idx,1)...
                      -mac_con(mac_em_idx,17).*(mac_spd(mac_em_idx,k)...
                      -ones(n_em,1)))./(2*mac_con(mac_em_idx,16));
  end
  % end rate calculation
 end
end
