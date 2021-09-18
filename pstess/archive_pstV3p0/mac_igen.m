function [bus_new] = mac_igen(i,k,bus,flag)
%syntax [bus_new] = mac_igen(i,k,bus,flag)
% 6:28 PM 18/8/97
% induction generator model
% bus_new is a modified bus matrix with the induction generator
% active and reactive powers subtracted from the original loads
% at the generator bus
%i is the induction generator number, 0 for vectorized computation
%k is the time step
%bus is the solved load flow bus data
%flag is 0 for initialization
%        1 for network interface
%        2 for dertermination of rates of change of states
%        3 for formation of linearized state matrix
%Purpose
%Simple Induction Generator Model
%Single Cage
%no leakage inductance saturation
%data format igen_con
%1 - induction generator number
%2 - busnumber
%3 - base MVA
%4 - rs
%5 - xs -stator leakage reactance
%6 - Xm - magnetizing reactance
%7 - rr
%8 - xr - rotor leakage reactance
%9 - H  - inertia constant generator + turbine in sec
%15 - fraction of bus load power taken by generator
%Note: induction generator pick up power from negative load 
% 
%Purpose: Induction Generator Model
%Version: 1.0
%Author: Graham Rogers
%Date July  1997
%(c) Copyright Cherry Tree Scientific Software 1997 _ All rights reserved
global  basmva basrad bus_int 
global  tmig  pig qig vdig vqig  idig iqig igen_con igen_pot
global  igen_int igbus n_ig
global  vdpig vqpig slig dvdpig dvqpig dslig
jay=sqrt(-1);
bus_new = bus;
if ~isempty(igen_con)
   if flag == 0;
      % initialisation
      disp(' initializing')
      if i == 0;
         %vector computation
         n_ig = length(igen_con(:,1));
         igbus=bus_int(igen_con(:,2));
         igen_pot = zeros(n_ig,7);
         igen_pot(:,1)=basmva./igen_con(:,3); %scaled mva base
         igen_pot(:,2)=ones(n_ig,1); %base kv
         ig_vm(:,1)=bus(igbus,2); %ind gen terminal voltage mag
         ig_ang(:,1)=bus(igbus,3)*pi/180; %ind gen term voltage angle
         v=ig_vm(:,1).*exp(jay*ig_ang(:,1));
         vdig(:,1)=real(v);
         vqig(:,1)=imag(v);
         pig(:,1)= bus(igbus,6).*igen_con(:,15);%ind generator power 
         %modify bus load power
         bus_new(igbus,6)=bus(igbus,6)-pig(:,1);
         igen_pot(:,3)=igen_con(:,5)+igen_con(:,6);%Xs
         igen_pot(:,4)=igen_con(:,8)+igen_con(:,6);%Xr
         igen_pot(:,5)=igen_con(:,5)+igen_con(:,6).*...
            igen_con(:,8)./igen_pot(:,4);%Xsp
         igen_pot(:,6)=igen_pot(:,3)-igen_pot(:,5);%(Xs-Xsp)
         igen_pot(:,7)=basrad*igen_con(:,7)./igen_pot(:,4); %1/Tr
         rs=igen_con(:,4);
         xs=igen_con(:,5);
         Xm=igen_con(:,6);
         rr=igen_con(:,7);
         xr=igen_con(:,8);
         % find initial slip
         slip_old=zeros(n_ig,1);
         slip_new=ones(n_ig,1);
         %Newton-Raphson iteration to determine initial slip for
         %induction generators
         iter = 0;
         err=max(abs(slip_new-slip_old));
         while err>=1e-8 & iter<30
            iter=iter+1;
            y=basrad.*slip_old./igen_pot(:,7);
            denom = ones(n_ig,1)+y.*y;
            zr=rs + y.*igen_pot(:,6)./denom;
            zi=igen_pot(:,5)+igen_pot(:,6)./denom;
            dzr=igen_pot(:,6).*(ones(n_ig,1)-...
               y.*y)./denom./denom;
            dzi=-2*igen_pot(:,6).*y./denom./denom;
            zmod2=zr.*zr+zi.*zi;
            dp=v.*conj(v).*(dzr.*zmod2-...
               2*zr.*(dzr.*zr+dzi.*zi));
            dp=dp./zmod2./zmod2;
            peig =v.*conj(v).*zr./zmod2;
            ynew=y-(peig - pig(:,1).*igen_pot(:,1))./dp;
            slip_new = ynew.*igen_pot(:,7)/basrad;
            err = max(abs(slip_new-slip_old));
            slip_old=slip_new;
         end
         if iter >=30
            error('induction generator slip calculation failed to converge')
         end
         slig(:,1)=slip_new;
         y=basrad*slig(:,1)./igen_pot(:,7);
         denom= ones(n_ig,1)+y.*y;
         zr=rs+y.*igen_pot(:,6)./denom;
         zi=igen_pot(:,5)+igen_pot(:,6)./denom;
         iig =v./(zr+jay*zi);
         sig = v.*conj(iig);
         peig = real(sig);
         qeig = imag(sig);
         %complex initial rotor states
         vp = v - (rs+ jay* igen_pot(:,5)).*iig; 
         vdpig = real(vp);
         vqpig =imag(vp);
         % initial prime mover torque
         tmig(:,1) = real(vp.*conj(iig));
         idig(:,1)=real(iig)./igen_pot(:,1);
         iqig(:,1)=imag(iig)./igen_pot(:,1);
         % modify qload
         qig(:,1) = qeig./igen_pot(:,1); 
         bus_new(igbus,7)=bus(igbus,7)-qig(:,1);
      else
         % generator by generator initialization
         error('Only vector computation allowed in induction generators')
      end
   end
   if flag == 1
      %network interface
      %no interface required for induction generators
   end
   if flag == 2
      %induction generator dynamics calculation
      if i == 0
         %vector calculation    
         idigm=idig(:,k).*igen_pot(:,1);%convert to machine base
         iqigm=iqig(:,k).*igen_pot(:,1);
         %Brereton, Lewis and Young motor model
         dvdpig(:,k)=-(iqigm.*igen_pot(:,6)+vdpig(:,k)).*...
            igen_pot(:,7)+vqpig(:,k).*slig(:,k)*basrad;
         dvqpig(:,k)=(idigm.*igen_pot(:,6)-vqpig(:,k)).*...
            igen_pot(:,7)-vdpig(:,k).*slig(:,k)*basrad;	
         dslig(:,k)=(tmig(:,k)-vdpig(:,k).*...
            idigm-vqpig(:,k).*iqigm)/2./igen_con(:,9);
      else
         error(' vector computation only for induction generator')
      end
   end
   if flag == 3
      %linearize
      %add code later
   end
end