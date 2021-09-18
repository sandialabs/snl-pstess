function f =mac_ib(i,k,bus,flag)
% Syntax: f = mac_ib(i,k,bus,flag)
% 5:05 PM 15/08/97
% generator infinite bus model
% forms a constant internal voltage model of a synchronous generator
% takes impedance from mac_con
% infinite buses are defined in the vector ibus_con, which has zero 
% entries for generators that are not to be inf. buses and unity
% entries for generators which are to be taken as infinite buses
% uses the appropriate initialization model to determine the initial
% internal voltage

% infinite bus indexes
global mac_con mac_ib_sub mac_ib_tra mac_ib_em mac_ib_idx not_ib_idx
global n_ib n_ib_sub n_ib_tra n_ib_em
% system data
global  basmva basrad syn_ref mach_ref sys_freq
global  bus_v bus_ang psi_re psi_im cur_re cur_im bus_int
% synchronous machine variables
global  mac_con mac_int mac_pot
global  mac_ang mac_spd eqprime edprime psikd psikq
global  curd curq curdg curqg fldcur
global  psidpp psiqpp vex eterm theta ed eq 
global  pmech pelect qelect


f = 0;
% exit if no infinite buses
if n_ib~=0
   jay = sqrt(-1);
   if flag ==0
      %initialize
      if i~=0 % non-vector initialization
         
         % check that machine i is defined as an infinite bus
         ib_chk = find(mac_ib_idx==i);
         if ~isempty(ib_chk)
            mac_pot(i,1) = basmva/mac_con(i,3); 
            % scaled MVA base
            %check for em generator
            if ~isempty(mac_ib_em)
               ib_em_chk = find(mac_ib_em==i);
               if ~isempty(ib_em_chk)
                  % initialize for classical machine model
                  busnum = bus_int(mac_con(i,2)); % bus number 
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
                  eprime = v + jay*mac_con(i,7)*curr; 
                  ei = eprime;
                  mac_ang(i,1) = atan2(imag(ei),real(ei)); 
                  % machine angle (delta)
                  rot = jay*exp(-jay*mac_ang(i,1)); 
                  % system reference frame rotation
                  eprime = eprime*rot;
                  edprime(i,1) = real(eprime); 
                  eqprime(i,1) = imag(eprime); 
                  psi_re(i,1) = sin(mac_ang(i,1))*edprime(i,1) + ...
                     cos(mac_ang(i,1))*eqprime(i,1); % real part of psi
                  psi_im(i,1) = -cos(mac_ang(i,1))*edprime(i,1) + ...
                     sin(mac_ang(i,1))*eqprime(i,1); % imag part of psi
               end
            end
            % check for transient generator
            if ~isempty(mac_ib_tra)
               ib_tra_chk = find(mac_ib_tra==i);
               if ~isempty(ib_tra_chk) ~=0
                  % initialize as transient generator
                  busnum = bus_int(mac_con(i,2)); % bus number 
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
                  rot = jay*exp(-jay*mac_ang(i,1)); % system reference frame rotation
                  eprime = eprime*rot;
                  edprime(i,1) = real(eprime); 
                  eqprime(i,1) = imag(eprime); 
                  psi_re(i,1) = sin(mac_ang(i,1))*edprime(i,1) + ...
                     cos(mac_ang(i,1))*eqprime(i,1); % real part of psi
                  psi_im(i,1) = -cos(mac_ang(i,1))*edprime(i,1) + ...
                     sin(mac_ang(i,1))*eqprime(i,1); % imag part of psi
               end
            end
            % check for subtransient generator model
            if ~isempty(mac_ib_sub)
               ib_sub_chk = find(mac_ib_sub == i);
               if ~isempty(ib_sub_chk)
                  % initialize as subtransient model
                  busnum = bus_int(mac_con(i,2)); % bus number 
                  % extract bus information
                  eterm(i,1) = bus(busnum,2);  % terminal bus voltage
                  theta(busnum,1) = bus(busnum,3)*pi/180;  
                  % terminal bus angle in radians
                  pelect(i,1) = bus(busnum,4)*mac_con(i,22); % system base   
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
                  ei = v + (mac_con(i,5)+jay*mac_con(i,11))*curr;
                  mac_ang(i,1) = atan2(imag(ei),real(ei)); 
                  % machine angle (delta)
                  rot = jay*exp(-jay*mac_ang(i,1));  % system reference frame rotation
                  curr = curr*rot;
                  curdg(i,1) = real(curr); curqg(i,1) = imag(curr);
                  v = v*rot;
                  ed(i,1) = real(v); eq(i,1) = imag(v);
                  eqra = eq(i,1)+mac_con(i,5)*curqg(i,1);
                  psidpp = eqra + mac_con(i,8)*curdg(i,1);
                  psikd(i,1) = eqra + mac_con(i,4)*curdg(i,1);
                  eqprime(i,1) = eqra + mac_con(i,7)*curdg(i,1);
                  edra = -ed(i,1)-mac_con(i,5)*curdg(i,1);
                  psiqpp = edra + mac_con(i,13)*curqg(i,1);
                  psikq(i,1) = edra + mac_con(i,4)*curqg(i,1);
                  edprime(i,1) = edra + mac_con(i,12)*curqg(i,1);
                  psi_re(i,1) = sin(mac_ang(i,1)).*(-psiqpp) + ...
                     cos(mac_ang(i,1)).*psidpp; % real part of psi
                  psi_im(i,1) = -cos(mac_ang(i,1)).*(-psiqpp) + ...
                     sin(mac_ang(i,1)).*psidpp; % imag part of psi
               end
            end
         end
      else
         % vector computation
         mac_pot(mac_ib_idx,1) = basmva*ones(n_ib,1)./mac_con(mac_ib_idx,3); 
         % scaled MVA base
         % check for em models
         if n_ib_em~=0
            busnum = bus_int(mac_con(mac_ib_em,2)); % bus number 
            % extract bus information
            eterm(mac_ib_em,1) = bus(busnum,2);  % terminal bus voltage
            theta(busnum,1) = bus(busnum,3)*pi/180;  
            % terminal bus angle in radians
            pelect(mac_ib_em,1) = bus(busnum,4).*mac_con(mac_ib_em,22);  
            % electrical output power, active
            qelect(mac_ib_em,1) = bus(busnum,5).*mac_con(mac_ib_em,23);  
            % electrical output power, reactive
            curr = sqrt(pelect(mac_ib_em,1).^2+qelect(mac_ib_em,1).^2)...
               ./eterm(mac_ib_em,1).*mac_pot(mac_ib_em,1);  % current magnitude
            phi = atan2(qelect(mac_ib_em,1),pelect(mac_ib_em,1)); 
            % power factor angle
            v = eterm(mac_ib_em,1).*exp(jay*theta(busnum,1)); 
            % voltage in real and imaginary parts
            % on system reference frame 
            curr = curr.*exp(jay*(theta(busnum,1)-phi)); % current in real and 
            % imaginary parts on system reference frame 
            eprime = v + jay*mac_con(mac_ib_em,7).*curr; 
            ei = eprime;
            mac_ang(mac_ib_em,1) = atan2(imag(ei),real(ei)); 
            % machine angle (delta)
            rot = jay*exp(-jay*mac_ang(mac_ib_em,1)); 
            % system reference frame rotation
            eprime = eprime.*rot;
            edprime(mac_ib_em,1) = real(eprime); 
            eqprime(mac_ib_em,1) = imag(eprime);
            psi_re(mac_ib_em,1) = sin(mac_ang(mac_ib_em,1)).*edprime(mac_ib_em,1) + ...
               cos(mac_ang(mac_ib_em,1)).*eqprime(mac_ib_em,1); % real part of psi
            psi_im(mac_ib_em,1) = -cos(mac_ang(mac_ib_em,1)).*edprime(mac_ib_em,1) + ...
               sin(mac_ang(mac_ib_em,1)).*eqprime(mac_ib_em,1); % imag part of psi      
         end
         % check for transient models
         if n_ib_tra ~=0
            busnum = bus_int(mac_con(mac_ib_tra,2)); % bus number 
            % extract bus information
            eterm(mac_ib_tra,1) = bus(busnum,2);  % terminal bus voltage
            theta(busnum,1) = bus(busnum,3)*pi/180;  
            % terminal bus angle in radians
            pelect(mac_ib_tra,1) = bus(busnum,4).*mac_con(mac_ib_tra,22);  
            % electrical output power, active
            qelect(mac_ib_tra,1) = bus(busnum,5).*mac_con(mac_ib_tra,23);  
            % electrical output power, reactive
            curr = sqrt(pelect(mac_ib_tra,1).^2+qelect(mac_ib_tra,1).^2) ...
               ./eterm(mac_ib_tra,1).*mac_pot(mac_ib_tra,1);  % current magnitude
            phi = atan2(qelect(mac_ib_tra,1),pelect(mac_ib_tra,1)); 
            % power factor angle
            v = eterm(mac_ib_tra,1).*exp(jay*theta(busnum,1)); 
            % voltage in real and imaginary parts
            % on system reference frame 
            curr = curr.*exp(jay*(theta(busnum,1)-phi)); % current in real and 
            % imaginary parts on system reference frame 
            eprime = v + (mac_con(mac_ib_tra,5)+jay*mac_con(mac_ib_tra,7)).*curr; 
            ei = v + (mac_con(mac_ib_tra,5)+jay*mac_con(mac_ib_tra,11)).*curr;
            mac_ang(mac_ib_tra,1) = atan2(imag(ei),real(ei)); 
            % machine angle (delta)
            rot = jay*exp(-jay*mac_ang(mac_ib_tra,1)); % system reference frame rotation
            eprime = eprime.*rot;
            edprime(mac_ib_tra,1) = real(eprime); 
            eqprime(mac_ib_tra,1) = imag(eprime); 
            psi_re(mac_ib_tra,1) = sin(mac_ang(mac_ib_tra,1)).*edprime(mac_ib_tra,1) + ...
               cos(mac_ang(mac_ib_tra,1)).*eqprime(mac_ib_tra,1); 
            % real part of psi
            psi_im(mac_ib_tra,1) = -cos(mac_ang(mac_ib_tra,1)).*edprime(mac_ib_tra,1) + ...
               sin(mac_ang(mac_ib_tra,1)).*eqprime(mac_ib_tra,1);
            % imag part of psi
         end
         % check for subtransient models
         if n_ib_sub~=0
            busnum = bus_int(mac_con(mac_ib_sub,2)); % bus number 
            % extract bus information
            eterm(mac_ib_sub,1) = bus(busnum,2);  % terminal bus voltage
            theta(busnum,1) = bus(busnum,3)*pi/180;  
            % terminal bus angle in radians
            pelect(mac_ib_sub,1) = bus(busnum,4).*mac_con(mac_ib_sub,22);  
            % electrical output power, active
            qelect(mac_ib_sub,1) = bus(busnum,5).*mac_con(mac_ib_sub,23);  
            % electrical output power, reactive
            curr = sqrt(pelect(mac_ib_sub,1).^2+qelect(mac_ib_sub,1).^2) ...
               ./eterm(mac_ib_sub,1).*mac_pot(mac_ib_sub,1);  % current magnitude
            phi = atan2(qelect(mac_ib_sub,1),pelect(mac_ib_sub,1)); 
            % power factor angle
            v = eterm(mac_ib_sub,1).*exp(jay*theta(busnum,1)); 
            % voltage in real and imaginary parts
            % on system reference frame 
            curr = curr.*exp(jay*(theta(busnum,1)-phi)); % current in real and 
            % imaginary parts on system reference frame 
            ei = v + (mac_con(mac_ib_sub,5)+jay*mac_con(mac_ib_sub,11)).*curr;
            mac_ang(mac_ib_sub,1) = atan2(imag(ei),real(ei)); 
            % machine angle (delta)
            rot = jay*exp(-jay*mac_ang(mac_ib_sub,1)); % system reference frame rotation
            curr = curr.*rot;
            curdg(mac_ib_sub,1) = real(curr); 
            curqg(mac_ib_sub,1) = imag(curr);
            curd(mac_ib_sub,1) = real(curr)./mac_pot(mac_ib_sub,1); 
            curq(mac_ib_sub,1) = imag(curr)./mac_pot(mac_ib_sub,1);
            v = v.*rot;
            ed(mac_ib_sub,1) = real(v); 
            eq(mac_ib_sub,1) = imag(v);
            eqra = eq(mac_ib_sub,1)+mac_con(mac_ib_sub,5).*curqg(mac_ib_sub,1);
            psidpp = eqra + mac_con(mac_ib_sub,8).*curdg(mac_ib_sub,1);
            edra = -ed(mac_ib_sub,1)-mac_con(mac_ib_sub,5).*curdg(mac_ib_sub,1);
            psiqpp = edra + mac_con(mac_ib_sub,13).*curqg(mac_ib_sub,1);
            psi_re(mac_ib_sub,1) = sin(mac_ang(mac_ib_sub,1)).*(-psiqpp) + ...
               cos(mac_ang(mac_ib_sub,1)).*psidpp; % real part of psi
            psi_im(mac_ib_sub,1) = -cos(mac_ang(mac_ib_sub,1)).*(-psiqpp) + ...
               sin(mac_ang(mac_ib_sub,1)).*psidpp; % imag part of psi
         end
      end
      % end initialization
   end
   if flag == 1
      % network interface
      % for all infinite buses set psi to the initial value
      if i ~=0
         %test for infinite bus
         ib_chk = find(mac_ib_idx==i);
         if ~isempty(ib_chk)
            psi_re(i,k) = psi_re(i,1);
            psi_im(i,k) = psi_im(i,1);
         end
      else
         psi_re(mac_ib_idx,k) = psi_re(mac_ib_idx,1);
         psi_im(mac_ib_idx,k) = psi_im(mac_ib_idx,1);
      end
   end
   if flag ==2
      % no dynamics
   end
end





