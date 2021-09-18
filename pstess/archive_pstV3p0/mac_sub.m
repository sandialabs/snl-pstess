function f = mac_sub(i,k,bus,flag)
% Syntax: f = mac_sub(i,k,bus,flag)
%
% Purpose: voltage-behind-subtransient-reactance generator
%            model, with vectorized computation option
%          state variables are: mac_ang, mac_spd, eqprime,
%                               psikd, edprime, psikq
%
% Input: i - generator number
%          - 0 for vectorized computation
%        k - integer time
%        bus - solved loadflow bus data
%        flag - 0 - initialization
%               1 - network interface computation
%               2 - generator dynamics computation
%               3 - state matrix building
%
% Output: f - dummy variable
%
% Files:
%
% See Also: pst_var, mac_em, mac_tra

% Algorithm: R. P. Schulz, "Synchronous machine modeling".
%
% Calls:
%
% Called By:

% (c) Copyright 1991-1999 Joe H. Chow/Cherry Tree Scientific Software - All Rights Reserved

% History (in reverse chronological order)
%
% Version 2.1
% Date:   September 1997
% Saturation modified for low Eqprime
% Version:  2.0
% Date:     June 1996
% Author:   Graham Rogers
% Purpose:  change to allow multple generator model types
% Modification:

% Version:  1.0
% Author:   Joe H. Chow
% Date:     March 1991
% Modification: Correction to saturation GJR May 6 1995

% define global variables
% system variables
global  basmva basrad  mach_ref sys_freq
global  bus_int

% synchronous machine variables
global  mac_con mac_pot
global  mac_ang mac_spd eqprime edprime psikd psikq
global  psi_re psi_im cur_re cur_im
global  curd curq curdg curqg fldcur
global  psidpp psiqpp vex eterm theta ed eq
global  pmech pelect qelect
global  dmac_ang dmac_spd deqprime dedprime dpsikd dpsikq
global  n_sub mac_sub_idx
global  pm_sig
f = 0;

jay = sqrt(-1);
if n_sub~=0
    if flag == 0; % initialization
        if i ~= 0
            % check that subtransient machine
            sub = find(mac_sub_idx==i);
            if ~isempy(sub)
                % check data
                % ensure xd" = xq"
                if mac_con(i,8)~=mac_con(i,13)

                    disp(['making xdpp = xqpp for generator ',num2str(i)])
                    mac_con(i,13) = mac_con(i,8);
                end
                if mac_con(i,14) == 0;mac_con(i,14) = 999.0;end %default transient time constant
                if mac_con(i,15) == 0;mac_con(i,15) = 999.0;end %default subtransient time constant

                busnum = bus_int(mac_con(i,2)); % bus number
                mac_pot(i,1)=basmva/mac_con(i,3); % scaled MVA base
                mac_pot(i,2)=1.0; % base kv
                mac_pot(i,8)=mac_con(i,7)-mac_con(i,4);
                mac_pot(i,9)=(mac_con(i,8)-mac_con(i,4))/mac_pot(i,8);
                mac_pot(i,7)=(mac_con(i,6)-mac_con(i,7))*mac_pot(i,9);
                mac_pot(i,10)=(mac_con(i,7)-mac_con(i,8))...
                    /mac_pot(i,8);
                mac_pot(i,6)=(mac_con(i,6)-mac_con(i,7))...
                    /mac_pot(i,8)*mac_pot(i,10);
                mac_pot(i,13)=mac_con(i,12)-mac_con(i,4);
                mac_pot(i,14)=(mac_con(i,13)-mac_con(i,4))...
                    /mac_pot(i,13);
                mac_pot(i,12)=(mac_con(i,11)-mac_con(i,12))...
                    *mac_pot(i,14);
                mac_pot(i,15)=(mac_con(i,12)-mac_con(i,13))...
                    /mac_pot(i,13);
                mac_pot(i,11)=(mac_con(i,11)-mac_con(i,12))...
                    /mac_pot(i,13)*mac_pot(i,15);
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
                %on genenearor base
                phi = atan2(qelect(i,1),pelect(i,1));
                % power factor angle
                v = eterm(i,1)*exp(jay*theta(busnum,1));
                % complex voltage
                % in system reference frame
                curr = curr*exp(jay*(theta(busnum,1)-phi)); % complex current
                % in system reference frame
                ei = v + (mac_con(i,5)+jay*mac_con(i,11))*curr;
                mac_ang(i,1) = atan2(imag(ei),real(ei));
                % machine angle (delta)
                mac_spd(i,1) = 1; % machine speed at steady state
                rot = jay*exp(-jay*mac_ang(i,1));  % system reference frame rotation
                curr = curr*rot;
                curdg(i,1) = real(curr); curqg(i,1) = imag(curr);% current in Park's frame
                curd(i,1) = real(curr)/mac_pot(i,1);% current on system base
                curq(i,1) = imag(curr)/mac_pot(i,1);
                mcurmag = abs(curr); % current magnitude on machine base
                pmech(i,1) = pelect(i,1)*mac_pot(i,1) + mac_con(i,5)*(mcurmag*mcurmag);
                %pmech = pelect + losses on machine base
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
                % compute saturation
                inv_sat = inv([0.64 0.8 1;1 1 1;1.44 1.2 1]);
                b = [0.8 1+mac_con(i,20) 1.2*(1+mac_con(i,21))];
                mac_pot(i,3) = b*inv_sat(1,:)';
                mac_pot(i,4) = b*inv_sat(2,:)';
                mac_pot(i,5) = b*inv_sat(3,:)';
                E_Isat = mac_pot(i,3)*eqprime(i,1)^2 ...
                    + mac_pot(i,4)*eqprime(i,1) + mac_pot(i,5);
                if eqprime(i,1)<0.8;E_Isat=eqprime(i,1);end
                vex(i,1) = E_Isat + mac_pot(i,6)*(eqprime(i,1)-...
                    psikd(i,1))+mac_pot(i,7)*curdg(i,1);
                fldcur(i,1) = vex(i,1);
                psi_re(i,1) = sin(mac_ang(i,1)).*(-psiqpp) + ...
                    cos(mac_ang(i,1)).*psidpp; % real part of psi
                psi_im(i,1) = -cos(mac_ang(i,1)).*(-psiqpp) + ...
                    sin(mac_ang(i,1)).*psidpp; % imag part of psi
            end
        else
            % vectorized computation
            % check parameters
            uets_idx = find(mac_con(mac_sub_idx,8)~=mac_con(mac_sub_idx,13));
            if ~isempty(uets_idx)
                mac_con(mac_sub_idx(uets_idx),13)=mac_con(mac_sub_idx(uets_idx),8);
                disp('xqpp made equal to xdpp at generators  '); disp((mac_sub_idx(uets_idx))')
            end
            notp_idx = find(mac_con(mac_sub_idx,14)==0);
            if ~isempty(notp_idx)
                mac_con(mac_sub_idx(notp_idx),14) = 999.0*ones(length(notp_idx),1);
            end
            notpp_idx = find(mac_con(mac_sub_idx,15)==0);
            if ~isempty(notpp_idx)
                mac_con(mac_sub_idx(notpp_idx),15) = 999.0*ones(length(notpp_idx),1);
                % set x'q = x"q
                mac_con(mac_sub_idx(notpp_idx),12) =...
                    mac_con(mac_sub_idx(notpp_idx),13);
            end
            busnum = bus_int(mac_con(mac_sub_idx,2)); % bus number
            mac_pot(mac_sub_idx,1) = basmva*ones(n_sub,1)./mac_con(mac_sub_idx,3);
            % scaled MVA base
            mac_pot(mac_sub_idx,2) = ones(n_sub,1); % base kv
            mac_pot(mac_sub_idx,8)=mac_con(mac_sub_idx,7)-mac_con(mac_sub_idx,4);
            mac_pot(mac_sub_idx,9)=(mac_con(mac_sub_idx,8)-mac_con(mac_sub_idx,4))...
                ./mac_pot(mac_sub_idx,8);
            mac_pot(mac_sub_idx,7)=(mac_con(mac_sub_idx,6)-mac_con(mac_sub_idx,7))...
                .*mac_pot(mac_sub_idx,9);
            mac_pot(mac_sub_idx,10)=(mac_con(mac_sub_idx,7)-mac_con(mac_sub_idx,8))...
                ./mac_pot(mac_sub_idx,8);
            mac_pot(mac_sub_idx,6)=(mac_con(mac_sub_idx,6)-mac_con(mac_sub_idx,7))...
                ./mac_pot(mac_sub_idx,8).*mac_pot(mac_sub_idx,10);
            mac_pot(mac_sub_idx,13)=mac_con(mac_sub_idx,12)-mac_con(mac_sub_idx,4);
            mac_pot(mac_sub_idx,14)=(mac_con(mac_sub_idx,13)-mac_con(mac_sub_idx,4))...
                ./mac_pot(mac_sub_idx,13);
            mac_pot(mac_sub_idx,12)=(mac_con(mac_sub_idx,11)-mac_con(mac_sub_idx,12))...
                .*mac_pot(mac_sub_idx,14);
            mac_pot(mac_sub_idx,15)=(mac_con(mac_sub_idx,12)-mac_con(mac_sub_idx,13))...
                ./mac_pot(mac_sub_idx,13);
            mac_pot(mac_sub_idx,11)=(mac_con(mac_sub_idx,11)-mac_con(mac_sub_idx,12))...
                ./mac_pot(mac_sub_idx,13).*mac_pot(mac_sub_idx,15);
            % extract bus information
            eterm(mac_sub_idx,1) = bus(busnum,2);  % terminal bus voltage
            theta(busnum,1) = bus(busnum,3)*pi/180;
            % terminal bus angle in radians
            pelect(mac_sub_idx,1) = bus(busnum,4).*mac_con(mac_sub_idx,22);
            % electrical output power, active
            qelect(mac_sub_idx,1) = bus(busnum,5).*mac_con(mac_sub_idx,23);
            % electrical output power, reactive
            curr = sqrt(pelect(mac_sub_idx,1).^2+qelect(mac_sub_idx,1).^2) ...
                ./eterm(mac_sub_idx,1).*mac_pot(mac_sub_idx,1);
            % current magnitude on generator base
            phi = atan2(qelect(mac_sub_idx,1),pelect(mac_sub_idx,1));
            % power factor angle
            v = eterm(mac_sub_idx,1).*exp(jay*theta(busnum,1));
            % voltage in real and imaginary parts
            % in system reference frame
            curr = curr.*exp(jay*(theta(busnum,1)-phi));
            % complex current in system reference frame
            ei = v + (mac_con(mac_sub_idx,5)+jay*mac_con(mac_sub_idx,11)).*curr;
            % voltage behind sub-transient reactance in system frame
            mac_ang(mac_sub_idx,1) = atan2(imag(ei),real(ei));
            % machine angle (delta)
            mac_spd(mac_sub_idx,1) = ones(n_sub,1);
            % machine speed at steady state
            rot = jay*exp(-jay*mac_ang(mac_sub_idx,1));
            % system reference frame rotation to Park's frame
            curr = curr.*rot;
            % current on generator base in Park's frame
            mcurmag = abs(curr);
            pmech(mac_sub_idx,1) = pelect(mac_sub_idx,1).*mac_pot(mac_sub_idx,1)...
                + mac_con(mac_sub_idx,5).*(mcurmag.*mcurmag);
            % mechanical power = electrical power + losses on generator base
            curdg(mac_sub_idx,1) = real(curr);
            curqg(mac_sub_idx,1) = imag(curr);
            % d and q axis current on generator base
            curd(mac_sub_idx,1) = real(curr)./mac_pot(mac_sub_idx,1);
            curq(mac_sub_idx,1) = imag(curr)./mac_pot(mac_sub_idx,1);
            % d and q axis currents on system base
            v = v.*rot;% voltage in Park's frame
            ed(mac_sub_idx,1) = real(v);
            eq(mac_sub_idx,1) = imag(v);
            % d and q axis voltages in Park's frame
            eqra = eq(mac_sub_idx,1)+mac_con(mac_sub_idx,5).*curqg(mac_sub_idx,1);
            % q axis voltage behind resistance
            psidpp = eqra + mac_con(mac_sub_idx,8).*curdg(mac_sub_idx,1);
            psikd(mac_sub_idx,1) = eqra + mac_con(mac_sub_idx,4).*curdg(mac_sub_idx,1);
            eqprime(mac_sub_idx,1) = eqra + mac_con(mac_sub_idx,7).*curdg(mac_sub_idx,1);
            edra = -ed(mac_sub_idx,1)-mac_con(mac_sub_idx,5).*curdg(mac_sub_idx,1);
            psiqpp = edra + mac_con(mac_sub_idx,13).*curqg(mac_sub_idx,1);
            psikq(mac_sub_idx,1) = edra + mac_con(mac_sub_idx,4).*curqg(mac_sub_idx,1);
            edprime(mac_sub_idx,1) = edra + mac_con(mac_sub_idx,12).*curqg(mac_sub_idx,1);
            % this is the negative of Edprime in block diagram
            % compute saturation
            inv_sat = inv([0.64 0.8 1;1 1 1;1.44 1.2 1]);
            b = [0.8*ones(n_sub,1) ones(n_sub,1)+mac_con(mac_sub_idx,20)...
                1.2*(ones(n_sub,1)+mac_con(mac_sub_idx,21))];
            mac_pot(mac_sub_idx,3) = b*inv_sat(1,:)';
            mac_pot(mac_sub_idx,4) = b*inv_sat(2,:)';
            mac_pot(mac_sub_idx,5) = b*inv_sat(3,:)';
            E_Isat = mac_pot(mac_sub_idx,3).*eqprime(mac_sub_idx,1).^2 ...
                + mac_pot(mac_sub_idx,4).*eqprime(mac_sub_idx,1) + mac_pot(mac_sub_idx,5);
            nosat_idx=find(eqprime(mac_sub_idx,1)<.8);
            if ~isempty(nosat_idx)
                E_Isat(nosat_idx)=eqprime(mac_sub_idx(nosat_idx),1);
            end
            vex(mac_sub_idx,1) = E_Isat + mac_pot(mac_sub_idx,6).*(eqprime(mac_sub_idx,1)-...
                psikd(mac_sub_idx,1))+mac_pot(mac_sub_idx,7).*curdg(mac_sub_idx,1);
            fldcur(mac_sub_idx,1) = vex(mac_sub_idx,1);
            psi_re(mac_sub_idx,1) = sin(mac_ang(mac_sub_idx,1)).*(-psiqpp) + ...
                cos(mac_ang(mac_sub_idx,1)).*psidpp; % real part of psi
            psi_im(mac_sub_idx,1) = -cos(mac_ang(mac_sub_idx,1)).*(-psiqpp) + ...
                sin(mac_ang(mac_sub_idx,1)).*psidpp; % imag part of psi
            % psi is in system base and is the voltage behind xpp
        end
        %end initialization
    end
    if flag == 1 % network interface computation
        if i ~= 0
            % check for subsynchronous machine
            sub = find(mac_sub_idx == i,1);
            if ~isempty(sub)
                mac_ang(i,k) = mac_ang(i,k) - mach_ref(k);
                % wrt machine referencek
                psidpp = mac_pot(i,9)*eqprime(i,k) + ...
                    mac_pot(i,10)*psikd(i,k);
                psiqpp = mac_pot(i,14)*edprime(i,k) + ...
                    mac_pot(i,15)*psikq(i,k);
                psi_re(i,k) = sin(mac_ang(i,k))*(-psiqpp) + ...
                    cos(mac_ang(i,k))*psidpp; % real part of psi
                psi_im(i,k) = -cos(mac_ang(i,k))*(-psiqpp) + ...
                    sin(mac_ang(i,k))*psidpp; % imag part of psi
            end
        else
            % vectorized computation

            mac_ang(mac_sub_idx,k) = mac_ang(mac_sub_idx,k)-mach_ref(k)*ones(n_sub,1);
            % wrt machine reference
            psidpp = mac_pot(mac_sub_idx,9).*eqprime(mac_sub_idx,k) + ...
                mac_pot(mac_sub_idx,10).*psikd(mac_sub_idx,k);
            psiqpp = mac_pot(mac_sub_idx,14).*edprime(mac_sub_idx,k) + ...
                mac_pot(mac_sub_idx,15).*psikq(mac_sub_idx,k);
            psi_re(mac_sub_idx,k) = sin(mac_ang(mac_sub_idx,k)).*(-psiqpp) + ...
                cos(mac_ang(mac_sub_idx,k)).*psidpp; % real part of psi
            psi_im(mac_sub_idx,k) = -cos(mac_ang(mac_sub_idx,k)).*(-psiqpp) + ...
                sin(mac_ang(mac_sub_idx,k)).*psidpp; % imag part of psi
        end
        % end of interface
    end

    if flag == 2 || flag == 3 % generator dynamics calculation
        if i ~= 0
            %check for subsynchronous machine
            sub = find(mac_sub_idx==i,1);
            if ~isempty(sub)
                psiqpp = mac_pot(i,14)*edprime(i,k) + ...
                    mac_pot(i,15)*psikq(i,k);
                psidpp = mac_pot(i,9)*eqprime(i,k) + ...
                    mac_pot(i,10)*psikd(i,k);
                curd(i,k) = sin(mac_ang(i,k))*cur_re(i,k) - ...
                    cos(mac_ang(i,k))*cur_im(i,k); % d-axis current
                curq(i,k) = cos(mac_ang(i,k))*cur_re(i,k) + ...
                    sin(mac_ang(i,k))*cur_im(i,k); % q-axis current
                curdg(i,k) = curd(i,k)*mac_pot(i,1);
                curqg(i,k) = curq(i,k)*mac_pot(i,1);
                mcurmag = abs(curdg(i,k) + jay*curqg(i,k));
                E_Isat = mac_pot(i,3)*eqprime(i,k)^2 ...
                    + mac_pot(i,4)*eqprime(i,k) + mac_pot(i,5);
                if eqprime(k,1)<0.8;E_Isat=eqprime(k,1);end
                fldcur(i,k) = E_Isat + mac_pot(i,6)...
                    *(eqprime(i,k)-psikd(i,k)) + mac_pot(i,7)...
                    *curdg(i,k);
                deqprime(i,k) = (vex(i,k)-fldcur(i,k))/mac_con(i,9);
                dpsikd(i,k) = (-psikd(i,k)+eqprime(i,k)-mac_pot(i,8)...
                    *curdg(i,k))/mac_con(i,10);
                dedprime(i,k) = (-edprime(i,k) - mac_pot(i,11)...
                    *(edprime(i,k)-psikq(i,k)) - mac_pot(i,12)...
                    *curqg(i,k))/mac_con(i,14);
                dpsikq(i,k) = (edprime(i,k)-psikq(i,k)-mac_pot(i,13)...
                    *curqg(i,k))/mac_con(i,15);
                ed(i,k) = -mac_con(i,5)*curdg(i,k) - (psiqpp...
                    -mac_con(i,13)*curqg(i,k));
                eq(i,k) = -mac_con(i,5)*curqg(i,k) + (psidpp...
                    -mac_con(i,8)*curdg(i,k));
                eterm(i,k) = sqrt(ed(i,k)^2+eq(i,k)^2);
                pelect(i,k) = eq(i,k)*curq(i,k) + ed(i,k)*curd(i,k);
                qelect(i,k) = eq(i,k)*curd(i,k) - ed(i,k)*curq(i,k);
                dmac_ang(i,k) = basrad*(mac_spd(i,k)-1.);
                Te = pelect(i,k)*mac_pot(i,1) + mac_con(i,5)*mcurmag*mcurmag;
                dmac_spd(i,k) = (pmech(i,k)+pm_sig(i,k)-Te...
                    -mac_con(i,17)*(mac_spd(i,k)-1)...
                    -mac_con(i,18)*(mac_spd(i,k)-sys_freq(k)))...
                    /(2*mac_con(i,16));
            end
        else
            % vectorized computation

            psiqpp = mac_pot(mac_sub_idx,14).*edprime(mac_sub_idx,k) + ...
                mac_pot(mac_sub_idx,15).*psikq(mac_sub_idx,k);
            psidpp = mac_pot(mac_sub_idx,9).*eqprime(mac_sub_idx,k) + ...
                mac_pot(mac_sub_idx,10).*psikd(mac_sub_idx,k);
            curd(mac_sub_idx,k) = sin(mac_ang(mac_sub_idx,k)).*cur_re(mac_sub_idx,k) - ...
                cos(mac_ang(mac_sub_idx,k)).*cur_im(mac_sub_idx,k); % d-axis current
            curq(mac_sub_idx,k) = cos(mac_ang(mac_sub_idx,k)).*cur_re(mac_sub_idx,k) + ...
                sin(mac_ang(mac_sub_idx,k)).*cur_im(mac_sub_idx,k); % q-axis current
            curdg(mac_sub_idx,k) = curd(mac_sub_idx,k).*mac_pot(mac_sub_idx,1);
            curqg(mac_sub_idx,k) = curq(mac_sub_idx,k).*mac_pot(mac_sub_idx,1);
            mcurmag = abs(curdg(mac_sub_idx,k)+jay*curqg(mac_sub_idx,k));
            E_Isat = mac_pot(mac_sub_idx,3).*eqprime(mac_sub_idx,k).^2 ...
                + mac_pot(mac_sub_idx,4).*eqprime(mac_sub_idx,k) + mac_pot(mac_sub_idx,5);
            nosat_idx=find(eqprime(mac_sub_idx,1)<.8);
            if ~isempty(nosat_idx)
                E_Isat(nosat_idx)=eqprime(mac_sub_idx(nosat_idx),k);
            end
            fldcur(mac_sub_idx,k) = E_Isat + mac_pot(mac_sub_idx,6)...
                .*(eqprime(mac_sub_idx,k)-psikd(mac_sub_idx,k)) + mac_pot(mac_sub_idx,7)...
                .*curdg(mac_sub_idx,k);
            deqprime(mac_sub_idx,k) = (vex(mac_sub_idx,k)-fldcur(mac_sub_idx,k))./mac_con(mac_sub_idx,9);
            dpsikd(mac_sub_idx,k) = (-psikd(mac_sub_idx,k)+eqprime(mac_sub_idx,k)-mac_pot(mac_sub_idx,8)...
                .*curdg(mac_sub_idx,k))./mac_con(mac_sub_idx,10);
            dedprime(mac_sub_idx,k) = (-edprime(mac_sub_idx,k) - mac_pot(mac_sub_idx,11)...
                .*(edprime(mac_sub_idx,k)-psikq(mac_sub_idx,k)) - mac_pot(mac_sub_idx,12)...
                .*curqg(mac_sub_idx,k))./mac_con(mac_sub_idx,14);
            dpsikq(mac_sub_idx,k) = (edprime(mac_sub_idx,k)-psikq(mac_sub_idx,k)-mac_pot(mac_sub_idx,13)...
                .*curqg(mac_sub_idx,k))./mac_con(mac_sub_idx,15);
            ed(mac_sub_idx,k) = -mac_con(mac_sub_idx,5).*curdg(mac_sub_idx,k) - (psiqpp...
                -mac_con(mac_sub_idx,13).*curqg(mac_sub_idx,k));
            eq(mac_sub_idx,k) = -mac_con(mac_sub_idx,5).*curqg(mac_sub_idx,k) + (psidpp...
                -mac_con(mac_sub_idx,8).*curdg(mac_sub_idx,k));
            eterm(mac_sub_idx,k) = sqrt(ed(mac_sub_idx,k).^2+eq(mac_sub_idx,k).^2);
            pelect(mac_sub_idx,k) = eq(mac_sub_idx,k).*curq(mac_sub_idx,k) + ed(mac_sub_idx,k).*curd(mac_sub_idx,k);
            qelect(mac_sub_idx,k) = eq(mac_sub_idx,k).*curd(mac_sub_idx,k) - ed(mac_sub_idx,k).*curq(mac_sub_idx,k);
            dmac_ang(mac_sub_idx,k) = basrad*(mac_spd(mac_sub_idx,k)-ones(n_sub,1));
            Te = pelect(mac_sub_idx,k).*mac_pot(mac_sub_idx,1) + mac_con(mac_sub_idx,5).*mcurmag.*mcurmag;
            dmac_spd(mac_sub_idx,k) = (pmech(mac_sub_idx,k)+ pm_sig(mac_sub_idx,k)-Te...
                -mac_con(mac_sub_idx,17).*(mac_spd(mac_sub_idx,k)-ones(n_sub,1))...
                -mac_con(mac_sub_idx,18).*(mac_spd(mac_sub_idx,k)-sys_freq(k)...
                *ones(n_sub,1)))./(2*mac_con(mac_sub_idx,16));
        end
        %end rate calculation
    end
end
