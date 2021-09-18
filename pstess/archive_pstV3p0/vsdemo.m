%vsdemo.m
%5:04 PM 15/08/97
% script file to demo vstab solutions
% performs a load flow with the load and generation increased
% from the starting value by a ratio input by the user
% modal analysis of the current loadflow Jacobian may be performed at any stage
% normally this will be close to the point at which the load flow fails to converge
% A user can choose to perform a new load flow based on either the current, previous or original voltage profile
% Version: 2.0
% Author:  Graham Rogers
% Version: 1.0
% Date:    July 1997
%Author:   Graham Rogers
%Date:     August 1994
%copyright (c) Joe Chow 1994-1997 - All rights reserved

clear
clear global
clear all
close
% global variables
pst_var 
global gen_chg_idx
global Q Ql

disp('vstab demo program')
load_bus = 3;
gen_bus=2;
swing_bus=1;
jay=sqrt(-1);
% load input data from m.file

[dfile,pathname]=uigetfile('d*.m','Select Data File');
if pathname == 0
   error(' you must select a valid data file')
else
   lfile =length(dfile);
   % strip off .m and convert to lower case
   dfile = lower(dfile(1:lfile-2));
   eval(dfile);
end
% check for valid dynamic data file
if isempty(bus)
   error(' the selected file is not a valid loadflow data file')
end
nbus = length(bus(:,1));
gen_chg_idx = ones(nbus,1);
chg_bus = find(gen_chg_idx);
g_bus_idx = find(bus(:,10)==2);
g_pq_idx = [];
gb_qload = bus(g_bus_idx,7);
bus_old = bus;
bus_sol = bus;
line_old = line;
line_sol = line;
lf_new = 2;
lf_num = 0;
pr_pre = 1;
pr_old = 1;
prat = 1;
while (lf_new >=1)
   if lf_new==1
      bus_sol = bus_old;
      line_sol = line_old;
      v_mag = [];
      pr_old = 1;
      pr_pre = 1;
      lf_num = 0;
   end
   if lf_new == 3
      bus_sol=bus_pre;
      line_sol = line_pre;
      plot_pr(lf_num) = [];
      v_mag(:,lf_num)= [];
      lf_num = lf_num - 1;
      pr_old = 1;
      prat = plot_pr(lf_num)
   end
   lf_new=2;
   while (lf_new==2)
      bus_pre = bus_sol;
      line_pre = line_sol; 
      pr_pre = prat;
      pr_text = num2str(pr_pre);
      disp(['previous power ratio = ' pr_text])
      pr_old = 1;
      lf_num = lf_num+1
      if lf_num ==1
         prat =1;
      else
         prat=input('power increase ratio >> ');
         if isempty(prat);prat = 1; end;
      end
      plot_pr(lf_num) = prat;
      bus_sol(:,4)=bus_sol(:,4)+(prat-pr_pre)*bus_sol(:,4);
      bus_sol(:,6)=bus_sol(:,6)+(prat-pr_pre)*bus_sol(:,6);
      bus_sol(chg_bus,7)=bus_sol(chg_bus,7)+(prat-pr_pre)*bus_sol(chg_bus,7);
      if ~isempty(g_pq_idx)
         bus_sol(g_bus_idx(g_pq_idx),7) = bus_sol(g_bus_idx(g_pq_idx),7) + (prat - pr_pre)*gb_qload(g_pq_idx);
      end
      gb_qload = gb_qload + (prat - pr_pre)*gb_qload;
      load_p(:,lf_num) = bus_sol(:,6);
      nbus= length(bus_sol(:,1));
      nline=length(line(:,1));
      
      [bus_sol,line_sol,line_flow] = loadflow(bus_sol,line_sol,1e-9,30, ...
         1.0,'n',2);
      chg_bus = find(gen_chg_idx);
      g_pq_idx = find(~gen_chg_idx(g_bus_idx));
      load_q(chg_bus,lf_num) = bus_sol(chg_bus,7);
      if ~isempty(g_pq_idx)
         load_q(g_bus_idx(g_pq_idx),lf_num) = gb_qload(g_pq_idx);
      end         
      gen_q(chg_bus,lf_num) = bus_sol(chg_bus,5);
      n_chg_bus = find(~gen_chg_idx);
      if ~isempty(n_chg_bus)
         gen_q(n_chg_bus,lf_num) =-( bus_sol(n_chg_bus,7)-load_q(n_chg_bus,lf_num));
      end
      gen_p(:,lf_num) = bus_sol(:,4);
      v_mag(:,lf_num) = bus_sol(:,2);
      if lf_num>1
         plot(plot_pr,v_mag)
         title('v/p curves')
         xlabel('power ratio')
         ylabel('voltage magnitude pu')
      end
      
      
      flag = 0;
      while(flag == 0)
         disp('You can examine the system data')
         disp('Type 1 to see initial bus data')
         disp('     2 to see line data')
         disp('     3 to see solved load flow bus solution')
         disp('     4 to see line flow')
         disp('     5 to see bus voltage magnitude profile')
         disp('     0 to quit')
         sel = input('enter selection ,[0]>> ');
         if isempty(sel);sel =0;end
         if sel == 1
            bus
            disp('paused: press any key to continue')
            pause 
         elseif sel == 2
            line_sol
            disp('paused: press any key to continue')
            pause
         elseif sel == 3
            bus_sol
            disp('paused: press any key to continue')
            pause
         elseif sel == 4
            line_flow
            disp('paused: press any key to continue')
            pause
         elseif sel == 5
            bar(bus_sol(:,2))
            title('bus voltage magnitude profile')
            xlabel('internal bus number')
            ylabel('voltage in pu')
            disp('paused: press any key to continue')
            pause
         elseif sel == 0
            flag = 1;
         else
            sel = 0;
            flag = 1;
         end
      end
      disp('Type 1 to do an additional loadflow')
      disp('Type 2 to go on to modal analysis')
      disp('Type 0 to quit vsdemo')
      alf=input('enter selection 0/1/2[0] >>   ');
      if isempty(alf);alf = 0;end
      if alf == 0; return;end
      if alf == 2
         lf_new=2;
         break
      else
         if alf~=1;error('invalid selection');end
         disp(' Type 1 to start from the original bus data')
         disp(' Type 2 to start load flow with current bus data')
         disp(' Type 3 to start load flow with previous bus data')
         lf_new = input('enter selection 1/2/3[2]>>  ');
         if isempty(lf_new);lf_new = 2;end
         if lf_new == 1; 
            plot_pr = [];
            break
         end
         if lf_new == 3;
            pr_pre = prat;
            break
         end
      end
   end
   if lf_new ==2
      %do modal analysis of the loadflow jacobian
      %form the sparse Y matrix
      [Y,nSW,nPV,nPQ,SB] = y_sparse(bus_sol,line_sol);
      % process bus data
      bus_no = bus_sol(:,1);
      V = bus_sol(:,2);
      ang = bus_sol(:,3)*pi/180;
      Pg = bus_sol(:,4);
      Qg = bus_sol(:,5);
      Pl = bus_sol(:,6);
      Ql = bus_sol(:,7);
      Gb = bus_sol(:,8);
      Bb = bus_sol(:,9);
      bus_type = round(bus_sol(:,10));
      sw_bno=ones(nbus,1);
      g_bno=sw_bno;
      % set up index for Jacobian calculation
      %% form PQV_no and PQ_no
      bus_zeros=zeros(nbus,1);
      bus_index=[1:1:nbus]';
      swing_index=find(bus_type==1); 
      sw_bno(swing_index)=bus_zeros(swing_index);
      PQV_no=find(bus_type >=2);
      PQ_no=find(bus_type==3);
      gen_index=find(bus_type==2);
      g_bno(gen_index)=bus_zeros(gen_index);    
      %sw_bno is a vector having ones everywhere but the swing bus locations
      %g_bno is a vector having ones everywhere but the generator bus locations
      % construct sparse angle reduction matrix
      il = length(PQV_no);
      ii = [1:1:il]';
      ang_red = sparse(ii,PQV_no,ones(il,1),il,nbus);
      % construct sparse voltage reduction matrix
      il = length(PQ_no);
      ii = [1:1:il]';
      volt_red = sparse(ii,PQ_no,ones(il,1),il,nbus);
      %form the system Jacobian
      [J11,J12,J21,J22]=form_jac(V,ang,Y,ang_red,volt_red);
      %Reduce the Jacobian to voltage terms only
      J22=J22-J21*inv(J11)*J12;
      %Find the eigenvectors and eigenvalues of the reduced inverse Jacobian
      [y,L]=eig(full(inv(J22)));
      %return to the full bus matrices for graphic output
      y=volt_red'*y*volt_red;
      L=volt_red'*diag(L);
      [my,iy]=max(abs((L)));
      ld_str = num2str(L(iy));
      disp(['the dominant eigenvalue ',ld_str])
      
      
      [lmax,ib]=max(abs(y(:,iy)));
      lm_str = num2str(lmax);
      ib_str = num2str(bus(ib,1));
      disp(['the maximum eigenvector entry is ',lm_str])
      disp(['the corresponding bus number is ',ib_str])
      eig_plot=input('do you wish to plot the eigenvector? - y/n[n]','s');
      if isempty(eig_plot); eig_plot = 'N';end
      if eig_plot ~= 'N'
         if eig_plot~='n'
            bar(abs(y(:,iy)));
            title('critical eigenvector')
            xlabel('internal bus number')
            ylabel('magnitude of eigenvector')
            
         end
      end
      disp('Type 1 to perform an additional load flow starting with the original bus data')
      disp('Type 2 to perform an additional loadflow starting with the current bus data')
      disp('Type 3 to perform an additional loadflow starting with the previous bus data')
      disp('Type 0 to exit demo')
      lf_new = input('enter your selection 0/1/2/3[0]  >>  ');
      if isempty(lf_new);lf_new=0;end
      if lf_new==0
         return
      end
   end
end

