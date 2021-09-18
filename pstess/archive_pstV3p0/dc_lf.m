function [rec_par,inv_par,line_par,tap,Sr,Si] = dc_lf(bus,line,dci_dc,dcr_dc)
%Syntax [rec_par,inv_par,line_par,tap,Sr,Si] = dc_lf(bus,line)
% Dc initialization model
% Outputs: rectifier firing angle
%          inverter extinction angle
%          converter transformer tap settings
%          complex ac load at LT terminals on system base
% Input:   ac bus matrix (from data or as modified in loadflow)
%          ac line matrix
% Version 1.0
% Author Graham Rogers
% Date   October 1996
% (c) Copyright Joe Chow 1996 - All rights reserved
global  bus_int  basmva
global  dcsp_con  dcl_con dcc_con load_con
global  r_idx  i_idx n_dcl  n_conv  ac_bus rec_ac_bus  inv_ac_bus
global  inv_ac_line  rec_ac_line ac_line dcli_idx
global  tap tapr tapi tmax tmin tstep tmaxr tmaxi tminr tmini tstepr tstepi
global  Vdc
jay = sqrt(-1);
% determine dc indexes
f = dc_indx(bus,line,dci_dc,dcr_dc);
if n_conv~=0 & n_dcl~=0
  Vdc=zeros(n_conv,1);% initialize to zero vector
  dc_bus = dcsp_con(:,1);
  Vac = bus(ac_bus,2);
  Vang = bus(ac_bus,3)*pi/180;
  al_min = dcsp_con(r_idx,7)*pi/180;
  al_max = dcsp_con(r_idx,8)*pi/180;
  rdc_bus = dc_bus(r_idx);
  ga_min = dcsp_con(i_idx,7)*pi/180;
  ga_max = dcsp_con(i_idx,8)*pi/180;
  idc_bus = dc_bus(i_idx);
  % get transformer taps from load flow
  tap = line(ac_line,6);
  tapr = line(rec_ac_line,6);
  tapi = line(inv_ac_line,6);
  tmax = line(ac_line,8);
  tmaxr = line(rec_ac_line,8);
  tmaxi = line(inv_ac_line,8);
  tmin = line(ac_line,9);
  tmini = line(inv_ac_line,9);
  tminr = line(rec_ac_line,9);
  tstep = line(ac_line,10);
  tstepi = line(inv_ac_line,10);
  tstepr = line(rec_ac_line,10);
  % set Vdc limits
  Vdc_rated = dcsp_con(i_idx,4);
  % set max and min values of dc voltage at 1% from nominal
  Vdc_max = 1.01*Vdc_rated;
  Vdc_min = 0.99*Vdc_rated;
  % in load flow the rectifier and inverter are represented by 
  % P and Q loads
  % Get P and Q in MW
  P = bus(ac_bus,6)*basmva;
  Q = bus(ac_bus,7)*basmva;
  % convert ac voltages to kV line-to-line
  Vac = Vac.*bus(ac_bus,13);
  % calculate commutating voltage
  xequ = sqrt(3)*dcsp_con(:,5)./dcsp_con(:,6);% eq transformer reactance
  Rc = dcsp_con(:,6).*dcsp_con(:,5)*3/pi;% in series as far as dc is concerned
  iac = (P-jay*Q)./(Vac.*exp(-jay*Vang))/sqrt(3);
  ang_iac = angle(iac);
  % get specified idc in kA
  idc = dcl_con(dcli_idx,8)./dcsp_con(i_idx,4);%dc power/dc voltage
  % adjust ac currents to have same angle as load flow but magnitude
  % corresponding to desired dc current
  iaci = idc.*dcsp_con(i_idx,6).*exp(jay*ang_iac(i_idx))*sqrt(6)/pi;
  iacr = idc.*dcsp_con(r_idx,6).*exp(jay*ang_iac(r_idx))*sqrt(6)/pi;
  % calculate the equivalent HT bus voltage (commutating voltage)
  VHT(i_idx,1) = abs(Vac(i_idx).*exp(jay*Vang(i_idx)) + jay*xequ(i_idx).*iaci);
  VHT(r_idx,1) = abs(Vac(r_idx).*exp(jay*Vang(r_idx)) + jay*xequ(r_idx).*iacr);
  Vdo = 3*sqrt(2)*VHT.*dcsp_con(:,6)/pi; % ideal dc voltages
  rdc = dcl_con(dcli_idx,3);% dc line resistance
  cm = ones(n_dcl,1)-dcl_con(dcli_idx,9)/100;

  % Assume that inverters are operating at gamm_min
  mode = ones(n_dcl,1);
  [gamma] = inv_lf(mode,idc,Vdc_max,Vdc_min,ga_min,ga_max,Vdo(i_idx),Rc(i_idx));
  
  % calculate rectifier dc voltage
  Vdc(r_idx) = Vdc(i_idx) + idc.*rdc;

  % calculate alpha

  [alpha,mode_new,idc] = rec_lf(idc,al_min,al_max,Vdo(r_idx),Rc(r_idx),cm);
  
  % check for mode change
  mode_ch = sum(mode-mode_new);
  if mode_ch ~=0

    % mode has changed - need to recalculate the inverter conditions
    [gamma] = inv_lf(mode,idc,Vdc_max,Vdc_min,ga_min,ga_max,Vdo(i_idx),Rc(i_idx));
    
   

  end
  alpha = alpha;
  gamma = gamma;
  % recalculate Vdo based on the modified firing and extinction angles
  Vdo(r_idx) = (Vdc(r_idx)+Rc(r_idx).*idc)./cos(alpha*pi/180);
  Vdo(i_idx) = (Vdc(i_idx)+Rc(i_idx).*idc)./cos(gamma*pi/180)
  % calculate ac power factor at LT bus 
  cphi = Vdc./Vdo
  oor_idx = find(abs(cphi)>=1);
  if ~isempty(oor_idx)
    bad_num = num2str(oor_idx');
    disp(['dc_lf failed to converge at converter ',bad_num])
    error('stop')
  end
  % converter power in per unit on ac system base 
  Pr = Vdc(r_idx).*idc/basmva
  Pi = -Vdc(i_idx).*idc/basmva
  tphi = sqrt(ones(n_conv,1)-cphi.*cphi)./cphi
  Qr = Pr.*tphi(r_idx)
  Qi = -Pi.*tphi(i_idx)
  % convert to load at the converter LT bus
  iacr = idc.*dcsp_con(r_idx,6)*sqrt(6)/pi;
  %convert to per unit
  iacr = iacr*sqrt(3).*bus(rec_ac_bus,13)/basmva
  % convert xequ to per unit on ac system base
  xequ = xequ*basmva./bus(ac_bus,13)./bus(ac_bus,13)/sqrt(3)
  iaci = idc.*dcsp_con(i_idx,6)*sqrt(6)/pi;
  % convert to per unit
  iaci = iaci*sqrt(3).*bus(inv_ac_bus,13)/basmva
  % calculate LT reactive power
  Qr = Qr - xequ(r_idx).*iacr.*iacr
  Qi = Qi - xequ(i_idx).*iaci.*iaci;
  Sr = (Pr + jay*Qr)
  Si = (Pi + jay*Qi)
  rec_par = zeros(n_dcl,3);
  inv_par = rec_par;
  rec_par(:,1) = alpha;
  inv_par(:,1) = gamma;
  rec_par(:,2) = Vdc(r_idx);
  inv_par(:,2) = Vdc(i_idx);
  rec_par(:,3) = Vdc(r_idx).*idc
  inv_par(:,3) = Vdc(i_idx).*idc
  line_par = idc
  tap(i_idx) = tapi;
  tap(r_idx) = tapr
end
return
