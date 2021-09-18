function [bus_sol,line_sol,line_flow,rec_par,inv_par,line_par] = lfdcs(bus,line,dci_dc,dcr_dc)
%Syntax:
%        [bus_sol,line_sol,line_flow,rec_par,inv_par,line_par] = lfdcs(bus,line,dci_dc,dcr_dc)  
%
% Purpose:   Solves load flow with one or more dc lines
% Inputs:    bus - ac bus specification matrix
%            line - ac line specification matrix
% Outputs:   bus_sol - solved ac bus specification file
%            line_sol - solved ac line specification file
%            line_flow - calculated flow on ac lines
%            rec_par - rectifier parameters
%                    - rec_par(:,1) firing angle degrees
%                    - rec_par(:,2) dc voltage kV
%                    - rec_par(:,3) dc power MW
%                    - rec_par(:,4) equi HT bus voltage kV
%            inv_par - inverterer parameters
%                    - inv_par(:,1) extinction angle degrees
%                    - inv_par(:,2) dc voltage kV
%                    - inv_par(:,3) dc power MW
%                    - inv_par(:,4) equi HT bus voltage kV
%           line_par - dc line current kA

% Calls:     loadflow
%            dc_lf
% Called by: s_simu
% Algorithm: iterates between AC and DC solutions until 
%            converged solution is reached
% Version:   1.0
% Author:    Graham Rogers
% Date:      February 1997
%            (c) copyright Joe Chow 1991-1997  - All rights reserved
%
global basmva bus_int
global dcsp_con  dcl_con  dcc_con
global  r_idx  i_idx n_dcl  n_conv  ac_bus rec_ac_bus  inv_ac_bus
global  inv_ac_line  rec_ac_line ac_line dcli_idx
global  tap tapr tapi tmax tmin tstep tmaxr tmaxi tminr tmini tsepr tsepi
global  Vdc
disp(' load flow with HVDC')

jay = sqrt(-1);
% perform load flow iterations
errv = 0;
itermax = 40;
iter = 0;
bus_old = bus;

while (errv == 0&&iter<=itermax)
  iter = iter + 1;
  [bus_sol,line_sol,line_flow] = loadflow(bus,line,1e-6,30, ...
                                 1.0,'n',2);
  if iter==1;
     dc_indx(bus,line,dci_dc,dcr_dc);
  end

    % perform dc load flow
    [rec_par,inv_par,line_par,tap,Sr,Si] = dc_lf(bus_sol,line_sol,dci_dc,dcr_dc);

    % set taps in load flow data
    line_sol(rec_ac_line,6) = tap(r_idx);
    line_sol(inv_ac_line,6) = tap(i_idx);
    % set loads at rectifier and inverter
    bus_sol(rec_ac_bus,6) = real(Sr);
    bus_sol(inv_ac_bus,6) = real(Si);
    bus_sol(rec_ac_bus,7) = imag(Sr);
    bus_sol(inv_ac_bus,7) = imag(Si);
    Sr_old = bus(rec_ac_bus,6)+jay*bus(rec_ac_bus,7);
    Si_old = bus(inv_ac_bus,6)+jay*bus(inv_ac_bus,7);
    % check convergence
    errdc = max(abs([(Sr_old-Sr), (Si_old-Si)]));
    if errdc>1e-5
      bus = bus_sol;
      line = line_sol;
    else
      errv = 1;
     % perform ac load flow
     [bus_sol,line_sol,line_flow] = loadflow(bus,line,1e-9,30, ...
        1.0,'n',2);
     % reperform dc load flow
     [rec_par,inv_par,line_par,tap,Sr,Si] = dc_lf(bus_sol,line_sol,dci_dc,dcr_dc);
     
     % set dc quantities based on final ac load flow
     P = bus_sol(ac_bus,6);
     Q = bus_sol(ac_bus,7);
     Vac = bus_sol(ac_bus,2);
     Vang = bus_sol(ac_bus,3)*pi/180;
     i_ac = (P - jay*Q)./Vac./exp(-jay*Vang);
     Vac = Vac.*bus_sol(ac_bus,13); %convert to kV
     % convert ac currents to dc
     idceq = abs(i_ac)*pi*basmva/3/sqrt(2)...
             ./bus_sol(ac_bus,13)./dcsp_con(:,6);
     % convert ac currents to kA
     i_ac = i_ac*basmva/sqrt(3)./bus_sol(ac_bus,13);
     %calculate equivalent HT bus voltage
     xequ = sqrt(3)*dcsp_con(:,5)./dcsp_con(:,6);% eq transformer reactance
     VHT = abs(Vac.*exp(jay*Vang) + jay*xequ.*i_ac);
     Vdo = 3*sqrt(2)*VHT.*dcsp_con(:,6)/pi;% ideal dc voltages
     %calculate dc voltages
     dc_ang(r_idx,1)  = rec_par(:,1)*pi/180;
     dc_ang(i_idx,1) = inv_par(:,1)*pi/180;
     Rc = dcsp_con(:,6).*dcsp_con(:,5)*3/pi;% in series as far as dc is concerned
     Vdc = Vdo.*cos(dc_ang) - Rc.*idceq;  
   end
end
if iter >= itermax
   imstr = int2str(itermax);
   disp(['dc load flow not converged in',imstr,' iterations'])
   error('stop')
else
   % dc and ac converged
   itstr= int2str(iter);
   disp([itstr,' dc load flow iterations'])
   disp('ac/dc solution converged')
end
return