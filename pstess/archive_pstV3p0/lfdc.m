% Purpose:   Solves load flow with one or more dc lines
% Calls:     loadflow
%            dcinit
% Algorithm: iterates between AC and DC solutions until 
%            converged solution is reached
% Version:   1.0
% Author:    Graham Rogers
% Date:      November 1996
%            (c) copyright Joe Chow 1996  - All rights reserved
%
global basmva bus_int
global  dcsp_con  dcl_con dcc_con load_con
global  r_idx  i_idx n_dcl  n_conv  con_ac_bus rec_ac_bus  inv_ac_bus
global  inv_ac_line  rec_ac_line ac_line dcli_idx
global  tap tapr tapi tmax tmin tstep tmaxr tmaxi tminr tmini tsepr tsepi
global  Vdc
disp(' load flow with HVDC')
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
   error(' the selected file is not a valid data file')
end

basmva = input('enter system base MVA - [100]');
if isempty(basmva);basmva =100;end

jay = sqrt(-1);
% perform load flow iterations
errv = 0;
itermax = 40;
iter = 0;
bus_old = bus;

while (errv == 0&iter<=itermax)
  iter = iter + 1;
%  [bus_sol,line_sol,line_flow] = loadflow(bus,line,1e-9,30, ...
%                                 1.0,'n',2);
  [bus_sol,line_sol,line_flow] = loadflow(bus,line,1e-6,30, ...
                                 1.0,'n',2); % J. Chow: convergence tolerance 
                                 % changed 09/18/2016


    % perform dc load flow
    [rec_par,inv_par,line_par,tap,Sr,Si] = dc_lf(bus_sol,line_sol,{[]},{[]});

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
    % [bus_sol,line_sol,line_flow] = loadflow(bus,line,1e-9,30, ...
    %                              1.0,'n',2);
    
    [bus_sol,line_sol,line_flow] = loadflow(bus,line,1e-6,30, ...   % change converence tolerance
                                 1.0,'n',2);                        % J. Chow, 09/18/2016
    
    %perform dc load flow
    [rec_par,inv_par,line_par,tap,Sr,Si] = dc_lf(bus_sol,line_sol,{[]},{[]});
    end
  
end
if iter >= itermax
   imstr = int2str(itermax);
   disp(['dc load flow not converged in ',imstr,' iterations'])
   error('stop')
else
   % dc and ac converged
   disp('ac/dc solution converged')
   itstr = int2str(iter);
   disp(['number of dc iterations ',itstr]) 
   flag = 0;
  while(flag == 0)
    disp('You can examine the system data')
    disp('Type 1 to see initial bus data')
    disp('     2 to see modified line data')
    disp('     3 to see solved load flow bus solution')
    disp('     4 to see line flow')
    disp('     5 to see bus voltage magnitude profile')
    disp('     6 to see bus voltage phase profile')
    disp('     7 to see the solved dc parameters')
    disp('     0 to quit')
    sel = input('enter selection >> ');
    if isempty(sel);sel=0;end
    if sel == 1
      disp('                                   Initial Bus Data')
      disp('                                      GENERATION             LOAD')
      disp('       BUS     VOLTS     ANGLE      REAL  REACTIVE      REAL  REACTIVE ')
      disp(bus_old(:,1:7))
      disp('paused: press any key to continue')
      pause
    elseif sel == 2
      disp('                               Modified Line  Data')
      disp(line_sol)
      disp('paused: press any key to continue')
      pause
    elseif sel == 3
      disp('                                    Solved Bus Data')
      disp('                                      GENERATION             LOAD')
      disp('       BUS     VOLTS     ANGLE      REAL  REACTIVE      REAL  REACTIVE ')
      disp(bus_sol(:,1:7))
      disp('paused: press any key to continue')
      pause
    elseif sel == 4
      disp('                  solved Line Flows')
      disp('                      LINE FLOWS                     ')
      disp('      LINE  FROM BUS    TO BUS      REAL  REACTIVE   ')
      disp(line_flow)
      disp('paused: press any key to continue')
      pause
    elseif sel == 5
      bar(bus_sol(:,2))
      title('bus voltage magnitude profile')
      xlabel('internal bus number')
      ylabel('voltage in pu')
      disp('paused: press any key to continue')
      pause
    elseif sel == 6
      bar(bus_sol(:,3))
      title('bus voltage phase profile')
      xlabel('internal bus number')
      ylabel('phase in degrees')
      disp('paused: press any key to continue')
      pause
    elseif sel ==7
      disp('          Rectifier')
      disp(' ')
      atxt = num2str(rec_par(:,1));
      vrectxt = num2str(rec_par(:,2));
      prectxt = num2str(rec_par(:,3));
      disp('alpha in deg')
      disp([atxt])
      disp('dc voltage in kV')
      disp([vrectxt])
      disp('Power in MW')
      disp([prectxt])
      disp('          Inverter')
      disp(' ')
      gtxt = num2str(inv_par(:,1));
      vinvtxt = num2str(inv_par(:,2));
      pinvtxt = num2str(inv_par(:,3));
      disp('gamma in degrees')
      disp([gtxt])
      disp('dc voltage in kV')
      disp([vinvtxt])
      disp('power in MW')
      disp([pinvtxt])
      idctxt = num2str(line_par);
      disp(' ')
      disp('line current in kA')
      disp([idctxt])
      disp('paused: press any key to continue')
      pause
    elseif sel == 0
      flag = 1;
    else
      error('invalid selection')
    end
  end
end
