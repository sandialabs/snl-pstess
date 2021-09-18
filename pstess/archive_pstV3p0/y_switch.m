% y_switch   
% 9:45 am July 20, 1998
% script file for calculating the reduced Y matrices 
% for the various switching options
% part of s_simu
% calls red_ybus
% Author: Graham Rogers
% (c) copyright Cherry Tree Scientific Software/ J.H. Chow 1997 - 1998
% All rights reserved
% modified to add option 7, clear fault with out loss of line
% suggested by Trong Nguyen and M.A. Pai University of Illinois
% errors in line to line, line to ground and line to line to ground corrected

% create pre-fault admittance matrix
if ~isempty(load_con);nload = length(load_con(:,1));end
[Y_gprf,Y_gncprf,Y_ncgprf,Y_ncprf,V_rgprf,V_rncprf,boprf] = red_ybus(bus,line);
bus_intprf = bus_int;% store the internal bus numbers for the pre_fault system
nbus = length(bus(:,1));
% create fault on matrices
bus_f = bus;
line_f = line;
f_type = sw_con(2,6);
if f_type <4|f_type==7
   f_nearbus = sw_con(2,2);
   bus_idx = find(bus_f(:,1)==f_nearbus);
   if isempty(bus_idx)
      error('faulted bus not specified correctly')
   end
   if f_type == 0|f_type==7
      %three phase fault zero impedance to ground
      bf = 1.0/1e-7;
   elseif f_type ==1
      % line to ground fault zf = zn*z0/(zn+z0)
      xf = sw_con(2,4)*sw_con(2,5)/(sw_con(2,4) + sw_con(2,5));
      xf = max(xf,1e-7);
      bf = 1.0/xf;
   elseif f_type == 2
      % line to line to ground  zf = zn + z0
      xf = sw_con(2,4) + sw_con(2,5);
      xf = max(xf,1e-7);
      bf = 1.0/xf;
   elseif f_type==3
      %line to line  zf = zn
      xf = sw_con(2,4);
      xf = max(xf,1e-7);
      bf = 1.0/xf; 
   end
   bus_f(bus_idx(1),9) = -bf;
end
if f_type == 4
   % remove line with no fault
   f_nearbus = sw_con(2,2);
   bus_idx = find(bus_f(:,1)==f_nearbus);
   f_farbus =  sw_con(2,3);
   line_idx = find((line_f(:,1)==f_nearbus&line_f(:,2)==f_farbus)...
      |(line_f(:,2)==f_nearbus&line_f(:,1)==f_farbus));
   %choose first instance of line
   if ~isempty(line_idx)
      line_f(line_idx(1),4) = 1.0e7; %make line reactance infinite
   else
      error('can not find faulted line in line data')
   end
end
if f_type == 5
   %loss of load
   f_nearbus = sw_con(2,2);
   bus_idx = find(bus_f(:,1)==f_nearbus);
   if isempty(bus_idx)
      error('can not find faulted bus in bus data')
   end
   bus_f(bus_idx(1),6) = 0.0;
   bus_f(bus_idx(1),7) = 0.0;
end
if f_type == 6;
   %no fault  
   f_nearbus = sw_con(2,2);
   bus_idx = find(bus_f(:,1)==f_nearbus);
   if isempty(bus_idx)
      error('can not find faulted bus in bus data')
   end
end
% form fault on reduced matrices
[Y_gf,Y_gncf,Y_ncgf,Y_ncf,V_rgf,V_rncf,bof] = red_ybus(bus_f,line_f);   % fault-on
% admittance matrix
bus_intf = bus_int;
%second switching point, clear fault at near end/add new line
if f_type<4
   f_farbus = sw_con(2,3);
   line_pf1 = line;
   bus_pf1 = bus;
   line_idx = find((line_pf1(:,1)==f_nearbus & line_pf1(:,2)==f_farbus)...
      |(line_pf1(:,2)==f_nearbus & line_pf1(:,1)==f_farbus));
   if isempty(line_idx)
      fb_str = num2str(f_farbus);
      fn_str = num2str(f_nearbus)
      disp(['can not find line between ',fn_str,' & ',fb_str,' in line data'])
      error('faulted line not specified correctly')
   end 
   line_pf1(line_idx(1),4) = 1.0e7; %make faulted line reactance infinite
   new_bus = max(bus(:,1))+10; % add new faulted bus
   max_pf1b = length(bus(:,1))+1;
   bus_pf1(max_pf1b,1) = new_bus;
   bus_pf1(max_pf1b,2) = 1.0;
   bus_pf1(max_pf1b,3:7)=zeros(1,5);
   bus_pf1(max_pf1b,9) = bus_f(bus_idx(1),9);
   bus_pf1(max_pf1b,10) = 3;  
   dlpf1 = length(line_pf1(:,1))+1; % add new line
   line_pf1(dlpf1,1)=new_bus;
   line_pf1(dlpf1,2)=f_farbus;
   line_pf1(dlpf1,3:6) = line(line_idx(1),3:6);
   [Y_gpf1,Y_gncpf1,Y_ncgpf1,Y_ncpf1,V_rgpf1,V_rncpf1,bopf1]...
      = red_ybus(bus_pf1,line_pf1);  % post-fault
   % admittance matrix
   bus_intpf1 = bus_int;
elseif f_type==4|f_type==5|f_type==6;
   % fault type is 4 or 5, 6 no change in system structure
   % set post fault data to fault data (bus_pf1 = bus_f)
   bus_pf1 = bus_f;
   line_pf1 = line_f;
   Y_gpf1 = Y_gf;
   Y_gncpf1 = Y_gncf;
   Y_ncgpf1 = Y_ncgf;
   Y_ncpf1 = Y_ncf; 
   V_rgpf1 = V_rgf;
   V_rncpf1 = V_rncf;
   bus_intpf1 = bus_intf;
   bopf1 = bof;
elseif f_type == 7
   % clear fault
   bus_pf1 = bus;
   line_pf1 = line;
   Y_gpf1 = Y_gprf;
   Y_gncpf1 = Y_gncprf;
   Y_ncgpf1 = Y_ncgprf;
   Y_ncpf1 = Y_ncprf; 
   V_rgpf1 = V_rgprf;
   V_rncpf1 = V_rncprf;
   bus_intpf1 = bus_intprf;
   bopf1 = boprf;
end
%third switching point, clear fault at remote end
if f_type<4
   line_pf2 = line_pf1;
   line_pf2(dlpf1,4) = 1.0e7;%open line
   bus_pf2 = bus_pf1;
   bus_pf2(max_pf1b,9)=0.0;%remove short
   [Y_gpf2,Y_gncpf2,Y_ncgpf2,Y_ncpf2,V_rgpf2,V_rncpf2,bopf2]...
      = red_ybus(bus_pf2,line_pf2);  % post-fault
   % admittance matrix
   bus_intpf2 = bus_int;
else
   % load type = 4 or 5, 6 or 7
   % no change in system structure
   bus_pf2 = bus_pf1;
   line_pf2 = line_pf1;
   Y_gpf2 = Y_gpf1;
   Y_gncpf2 = Y_gncpf1;
   Y_ncgpf2 = Y_ncgpf1;
   Y_ncpf2 = Y_ncpf1; 
   V_rgpf2 = V_rgpf1;
   V_rncpf2 = V_rncpf1;
   bus_intpf2 = bus_intpf1;
   bopf2 = bopf1;
end