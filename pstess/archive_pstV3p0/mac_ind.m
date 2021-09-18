function bus_new = mac_ind(i,k,bus,flag)
%syntax [bus_new]=mac_ind(i,k,bus,flag)
%bus_new is bus with the power and reactive power loads 
%modified to subtract the motor loads
%modification is made only when the motors are initialized
%i.e. flag = 0
%i is the motor number, 0 for vectorized computation
%k is the time step
%bus is the solved load flow bus data
%flag is 0 for initialization
%        1 for network interface
%        2 for dertermination of rates of change of states
%        3 for formation of linearized state matrix
%Purpose
% Induction Motor Model
%data format ind_con
%1 - motor number
%2 - busnumber
%3 - base MVA
%4 - rs
%5 - xs -stator leakage reactance
%6 - Xm - magnetizing reactance
%7 - rr
%8 - xr - rotor leakage reactance
%9 - H  - inertia constant motor + load in sec
%10 - r2 - double cage resistance
%11 - x2 - intercage reactance
%12 - dbf - deep bar factor
%13 - isat - current at which leakage inductance starts to saturate
%15 - fraction of bus load power taken by motor 
% if entry 15 is zero, it is assumed that the motor is to be started on the specified bus
% 
%Author Graham Rogers
%Date November 1995
%
% Version 2 added deep bar, double cage and leakage inductance saturation
%June 2002
global basmva basrad bus_int bus_v
global tload t_init p_mot q_mot vdmot vqmot  idmot iqmot ind_con ind_pot t_mot
global ind_int motbus sat_idx dbc_idx db_idx
global vdp vqp slip dvdp dvqp dslip mld_con
jay=sqrt(-1);
bus_new=bus;

if ~isempty(ind_con)
    if flag == 0;
        % initialisation
        if i == 0;
            %vector computation
            motnum=length(ind_con(:,1));
            motbus=bus_int(ind_con(:,2));
            ind_pot(:,1)=basmva./ind_con(:,3); %scaled mva base
            ind_pot(:,2)=ones(motnum,1); %base kv
            mot_vm(:,1)=bus(motbus,2); %motor terminal voltage mag
            mot_ang(:,1)=bus(motbus,3)*pi/180; %motor term voltage angle
            v=mot_vm(:,1).*exp(jay*mot_ang(:,1));
            vdmot(:,1)=real(v);
            vqmot(:,1)=imag(v);
            p_mot(:,1)=bus(motbus,6).*ind_con(:,15);%motor power demand
            %modify bus load power
            bus_new(motbus,6)=bus(motbus,6)-p_mot(:,1);
            rs=ind_con(:,4);
            xs=ind_con(:,5);
            Xm=ind_con(:,6);
            rr=ind_con(:,7);
            xr=ind_con(:,8);
            rr2 = ind_con(:,10);
            xr2 = ind_con(:,11);
            dbf = ind_con(:,12);
            isat = ind_con(:,13);
            dbc_idx = find(rr2~=0);% motors with double cage
            db_idx = find(dbf~=0);% motors with deep bars
            sat_idx = find(isat~=0);% motors with leakage inductance saturation
            ind_pot(:,3)=xs+Xm;%Xs
            ind_pot(:,4)=xr+Xm;%Xr
            ind_pot(:,5)=xs+Xm.*xr./ind_pot(:,4);%Xsp
            ind_pot(:,6)=ind_pot(:,3)-ind_pot(:,5);%(Xs-Xsp)
            ind_pot(:,7)=basrad*rr./ind_pot(:,4); %1/Tr
            
            % index of motors to be initialized for running 
            run_ind = find(ind_con(:,15)~=0);
            motrun=length(run_ind);%number of running motors
            start_ind = find(ind_con(:,15)==0);
            motstart = length(start_ind);% number of starting motors
            %assumes motor starting if power fraction zero
            % find initial slip
            slip_old=zeros(motnum,1);
            slip_new=ones(motnum,1);
            % reset ind_pot for double cage and deepbar rotor machines
            s = 0.01*ones(motnum,1);s(start_ind)=ones(motstart,1);
            if ~isempty(dbc_idx)
                [rdc,xdc]=dbcage(rr(dbc_idx),xr(dbc_idx),rr2(dbc_idx),xr2(dbc_idx),s(dbc_idx)); 
                ind_pot(dbc_idx,4) = Xm(dbc_idx)+xdc;
                ind_pot(dbc_idx,5)=xs(dbc_idx)+Xm(dbc_idx).*xdc./ind_pot(dbc_idx,4);%Xsp
                ind_pot(dbc_idx,6) = ind_pot(dbc_idx,3)-ind_pot(dbc_idx,5);
                ind_pot(dbc_idx,7)=basrad*rdc./ind_pot(dbc_idx,4); %1/Tr
            end
            if ~isempty(db_idx)
                [rdb,xdb]=deepbar(rr(db_idx),dbf(db_idx),s(db_idx));
                ind_pot(db_idx,4) = Xm(db_idx)+xr(db_idx)+xdb;
                ind_pot(db_idx,5) = xs(db_idx)+Xm(db_idx).*(xr(db_idx)+xdb)./ind_pot(db_idx,4);%Xsp
                ind_pot(db_idx,6) = ind_pot(db_idx,3)-ind_pot(db_idx,5);
                ind_pot(db_idx,7) = basrad*rdb./ind_pot(db_idx,4); %1/Tr
            end
            
            
            % Set defaults for motor starting
            imot=zeros(motnum,1);
            pem = zeros(motnum,1);
            qem = zeros(motnum,1);
            vdp(:,1)=zeros(motnum,1);
            vqp(:,1)=zeros(motnum,1);
            vp = imot;
            t_init = ones(motnum,1); % default for motor starting
            %Newton-Raphson iteration to determine initial slip for
            %running motors
            
            if motrun~=0 %check that some motors are running
                iter = 0;
                err=max(abs(slip_new-slip_old));
                while (err>=1e-8) && (iter<30)
                    iter=iter+1;
                    y=basrad.*slip_old(run_ind)./ind_pot(run_ind,7);
                    denom = ones(motrun,1)+y.*y;
                    zr=rs(run_ind)+y.*ind_pot(run_ind,6)./denom;
                    zi=ind_pot(run_ind,5)+ind_pot(run_ind,6)./denom;
                    dzr=ind_pot(run_ind,6).*(ones(motrun,1)-...
                        y.*y)./denom./denom;
                    dzi=-2*ind_pot(run_ind,6).*y./denom./denom;
                    zmod2=zr.*zr+zi.*zi;
                    dp=v(run_ind).*conj(v(run_ind)).*(dzr.*zmod2-...
                        2*zr.*(dzr.*zr+dzi.*zi));
                    dp=dp./zmod2./zmod2;
                    pem(run_ind)=v(run_ind).*conj(v(run_ind)).*zr./zmod2;
                    ynew=y-(pem(run_ind)- ...
                        p_mot(run_ind,1).*ind_pot(run_ind,1))./dp;
                    slip_new(run_ind)=ynew.*ind_pot(run_ind,7)/basrad;
                    err = max(abs(slip_new-slip_old));
                    slip_old=slip_new;
                    if ~isempty(dbc_idx)
                        [rdc,xdc]=dbcage(rr(dbc_idx),xr(dbc_idx),rr2(dbc_idx),xr2(dbc_idx),slip_new(dbc_idx)); 
                        ind_pot(dbc_idx,4) = Xm(dbc_idx)+xdc;
                        ind_pot(dbc_idx,5)=xs(dbc_idx)+Xm(dbc_idx).*xdc./ind_pot(dbc_idx,4);%Xsp
                        ind_pot(dbc_idx,6) = ind_pot(dbc_idx,3)-ind_pot(dbc_idx,5);
                        ind_pot(dbc_idx,7)=basrad*rdc./ind_pot(dbc_idx,4); %1/Tr
                    end
                    if ~isempty(db_idx)
                        [rdb,xdb]=deepbar(rr(db_idx),dbf(db_idx),slip_new(db_idx));
                        ind_pot(db_idx,4) = Xm(db_idx)+xr(db_idx)+xdb;
                        ind_pot(db_idx,5) = xs(db_idx)+Xm(db_idx).*(xr(db_idx)+xdb)./ind_pot(db_idx,4);%Xsp
                        ind_pot(db_idx,6) = ind_pot(db_idx,3)-ind_pot(db_idx,5);
                        ind_pot(db_idx,7) = basrad*rdb./ind_pot(db_idx,4); %1/Tr
                    end
                end
                if iter >=30
                    uiwait(msgbox('induction motor slip calculation failed to converge','mac_ind error','modal'))
                    return
                end
            end
            slip(:,1)=slip_new;
            ind_ldto(0,1);
            y=basrad*slip(:,1)./ind_pot(:,7);
            denom= ones(motnum,1)+y.*y;
            zr=rs+y.*ind_pot(:,6)./denom;
            zi=ind_pot(:,5)+ind_pot(:,6)./denom;
            if ~isempty(run_ind)
                imot(run_ind)=v(run_ind)./(zr(run_ind)+jay*zi(run_ind));
                sm(run_ind)=v(run_ind).*conj(imot(run_ind));
                pem(run_ind)=real(sm(run_ind));
                qem(run_ind)=imag(sm(run_ind));
                %complex initial rotor states
                vp(run_ind) = v(run_ind) - (rs(run_ind)+ jay* ind_pot(run_ind,5)).*imot(run_ind); 
                vdp(run_ind,1)=real(vp(run_ind));
                vqp(run_ind,1)=imag(vp(run_ind));
            end
            idmot(:,1)=real(imot)./ind_pot(:,1);
            iqmot(:,1)=imag(imot)./ind_pot(:,1);
            % modify qload 
            bus_new(motbus,7)=bus(motbus,7)-qem./ind_pot(:,1);
            tlm = vdp(:,k).*real(imot)+vqp(:,k).*imag(imot);
            trat = tlm./tload(:,1);
            % modify load specification to get initial load correct
            
            mld_con(run_ind,[3 5])=diag(trat(run_ind))*mld_con(run_ind,[3 5]);
        else
            error('motor by motor initialization not supported')  
        end
    end
    if flag == 1
        v = bus_v(motbus,k);
        vdmot(:,k)=real(v);
        vqmot(:,k)=imag(v);
    end
    if flag == 2
        %motor dynamics calculation
        if i == 0
            %vector calculation
            
            ind_ldto(0,k);
            idm=idmot(:,k).*ind_pot(:,1);%convert to machine base
            iqm=iqmot(:,k).*ind_pot(:,1);
            rs=ind_con(:,4);
            xs=ind_con(:,5);
            Xm=ind_con(:,6);
            rr=ind_con(:,7);
            xr=ind_con(:,8);
            rr2 = ind_con(:,10);
            xr2 = ind_con(:,11);
            dbf = ind_con(:,12);
           
            imot = abs(idm+jay*iqm);
            if ~isempty(sat_idx)
                % saturation of leakage inductance
                ism = imot(sat_idx);isat = ind_con(sat_idx,13);
                ir = jay*Xm(sat_idx).*(idm(sat_idx)+jay*iqm(sat_idx))./(rr(sat_idx)+jay*ind_pot(sat_idx,4)); 
                gs = dessat(ism,isat);
                gr = dessat(abs(ir),isat);
                xs(sat_idx) = xs(sat_idx).*(1+gs)/2;
                xr(sat_idx) = xr(sat_idx).*(1+gr)/2;
                ind_pot(sat_idx,3) = Xm(sat_idx)+xs(sat_idx);
                ind_pot(sat_idx,4) = Xm(sat_idx)+xr(sat_idx);
                ind_pot(sat_idx,5)=xs(sat_idx)+Xm(sat_idx).*xr(sat_idx)./ind_pot(sat_idx,4);%Xsp
                ind_pot(sat_idx,6) = ind_pot(sat_idx,3)-ind_pot(sat_idx,5);
                ind_pot(sat_idx,7)=basrad*rr(sat_idx)./ind_pot(sat_idx,4); %1/Tr
            end
            if ~isempty(dbc_idx)
                % reset double cage
                [rdc,xdc]=dbcage(rr(dbc_idx),xr(dbc_idx),rr2(dbc_idx),xr2(dbc_idx),slip(dbc_idx,k)); 
                ind_pot(dbc_idx,4) = Xm(dbc_idx)+xdc;
                ind_pot(dbc_idx,5)=xs(dbc_idx)+Xm(dbc_idx).*xdc./ind_pot(dbc_idx,4);%Xsp
                ind_pot(dbc_idx,6) = ind_pot(dbc_idx,3)-ind_pot(dbc_idx,5);
                ind_pot(dbc_idx,7)=basrad*rdc./ind_pot(dbc_idx,4); %1/Tr
            end
            if ~isempty(db_idx)
                % reset deepbar
                [rdb,xdb]=deepbar(rr(db_idx),dbf(db_idx),slip(db_idx,k));
                ind_pot(db_idx,4) = Xm(db_idx)+xr(db_idx)+xdb;
                ind_pot(db_idx,5) = xs(db_idx)+Xm(db_idx).*(xr(db_idx)+xdb)./ind_pot(db_idx,4);%Xsp
                ind_pot(db_idx,6) = ind_pot(db_idx,3)-ind_pot(db_idx,5);
                ind_pot(db_idx,7) = basrad*rdb./ind_pot(db_idx,4); %1/Tr
            end
            %Brereton, Lewis and Young motor model
            dvdp(:,k)=-(iqm.*ind_pot(:,6)+vdp(:,k)).*...
                ind_pot(:,7)+vqp(:,k).*slip(:,k)*basrad;
            dvqp(:,k)=(idm.*ind_pot(:,6)-vqp(:,k)).*...
                ind_pot(:,7)-vdp(:,k).*slip(:,k)*basrad;
            t_mot(:,k) = vdp(:,k).*idm+vqp(:,k).*iqm;
            dslip(:,k)=(tload(:,k)-t_mot(:,k))/2./ind_con(:,9);
        else
            error('motor by motor mode is not supported')  
        end
    end
    if flag == 3
        %linearize
        %add code later
    end
end