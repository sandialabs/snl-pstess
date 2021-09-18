% Purpose:   Solves load flow with one or more dc lines
%
% Calls:     loadflow, dcinit
%
% Algorithm: iterates between AC and DC solutions until
%            converged solution is reached

%-----------------------------------------------------------------------------%
% Version history
%
% Version:   1.1
% Author:    Ryan Elliott
% Date:      July 2020
% Purpose:   Updated to accommodate batch processing mode
%
% Version:   1.0
% Author:    Graham Rogers
% Date:      November 1996
%-----------------------------------------------------------------------------%

clear all; close all; clc;

dypar = get_dypar();                     % get parameter settings
batch_mode = dypar.batch_mode;           % batch processing mode

disp('load flow with hvdc')
% input data file
if batch_mode
    [dfile,pathname] = get_path();
else
    [dfile,pathname] = uigetfile('d*.m','Select Data File');
end

% import base case data
if (pathname == 0)
    error('lfdc: you must select a valid data file.');
else
    full_dfile = fullfile(pathname, dfile);
    run(full_dfile);
    run('import_var');                   % set up global variables
end

% check for valid dynamic data file
if isempty(bus)
    error('lfdc: the selected file is not a valid data file (bus).');
end

% specify mva base for the system
if batch_mode
    g.sys.basmva = dypar.basmva;
else
    g.sys.basmva = input('Enter the system MVA base, [100]: ');
    if isempty(g.sys.basmva)
        g.sys.basmva = 100;  % default
    end
end

% perform load flow iterations
errv = 0;
itermax = 50;
iter = 0;
bus_old = bus;

while ((errv == 0) & (iter <= itermax))
    iter = iter + 1;

    % J. Chow: convergence tolerance changed 1e-9 to 1e-6 (2016)
    [bus_sol,line_sol,line_flow] = loadflow(bus,line,1e-6,50,1.0,'n',2);

    % perform dc load flow
    [rec_par,inv_par,line_par,tap,Sr,Si] = dc_lf(bus_sol,line_sol,{[]},{[]});

    % set taps in load flow data
    line_sol(g.dc.rec_ac_line,6) = tap(g.dc.r_idx);
    line_sol(g.dc.inv_ac_line,6) = tap(g.dc.i_idx);

    % set loads at rectifier and inverter
    bus_sol(g.dc.rec_ac_bus,6) = real(Sr);
    bus_sol(g.dc.inv_ac_bus,6) = real(Si);
    bus_sol(g.dc.rec_ac_bus,7) = imag(Sr);
    bus_sol(g.dc.inv_ac_bus,7) = imag(Si);
    Sr_old = bus(g.dc.rec_ac_bus,6) + 1j*bus(g.dc.rec_ac_bus,7);
    Si_old = bus(g.dc.inv_ac_bus,6) + 1j*bus(g.dc.inv_ac_bus,7);

    % check convergence
    errdc = max(abs([(Sr_old-Sr), (Si_old-Si)]));
    if (errdc > 1e-5)
        bus = bus_sol;
        line = line_sol;
    else
        errv = 1;
        % perform ac load flow
        % J. Chow: convergence tolerance changed 1e-9 to 1e-6 (2016)
        [bus_sol,line_sol,line_flow] = loadflow(bus,line,1e-6,50,1.0,'n',2);

        % perform dc load flow
        [rec_par,inv_par,line_par,tap,Sr,Si] = dc_lf(bus_sol,line_sol,{[]},{[]});
    end
end

if (iter > itermax)
    estr = '\nlfdc: dc load flow failed to converge in %0.0f iterations.';
    error(sprintf(estr,itermax));
else
    % dc and ac converged
    dstr = '\nlfdc: ac/dc solution converged after %0.0f dc load flow iterations.';
    disp(sprintf(dstr,iter));

    flag = 0;
    while (flag == 0)
        disp('You can examine the system data')
        disp('Type 1 to see initial bus data')
        disp('     2 to see modified line data')
        disp('     3 to see solved load flow bus solution')
        disp('     4 to see line flow')
        disp('     5 to see bus voltage magnitude profile')
        disp('     6 to see bus voltage phase profile')
        disp('     7 to see the solved dc parameters')
        disp('     0 to quit')

        if batch_mode
            sel = 0;      % quit without displaying solution
        else
            sel = input('Enter your selection, (0--7)[0]: ');
            if isempty(sel)
                sel = 0;  % quit without displaying solution (default)
            end
        end

        if (sel == 1)
            disp('                                   Initial Bus Data')
            disp('                                      GENERATION             LOAD')
            disp('       BUS     VOLTS     ANGLE      REAL  REACTIVE      REAL  REACTIVE ')
            disp(bus_old(:,1:7))
            disp('paused: press any key to continue')
            pause
        elseif (sel == 2)
            disp('                               Modified Line  Data')
            disp(line_sol)
            disp('paused: press any key to continue')
            pause
        elseif (sel == 3)
            disp('                                    Solved Bus Data')
            disp('                                      GENERATION             LOAD')
            disp('       BUS     VOLTS     ANGLE      REAL  REACTIVE      REAL  REACTIVE ')
            disp(bus_sol(:,1:7))
            disp('paused: press any key to continue')
            pause
        elseif (sel == 4)
            disp('                  solved Line Flows')
            disp('                      LINE FLOWS                     ')
            disp('      LINE  FROM BUS    TO BUS      REAL  REACTIVE   ')
            disp(line_flow)
            disp('paused: press any key to continue')
            pause
        elseif (sel == 5)
            bar(bus_sol(:,2))
            title('bus voltage magnitude profile')
            xlabel('internal bus number')
            ylabel('voltage in pu')
            disp('paused: press any key to continue')
            pause
        elseif (sel == 6)
            bar(bus_sol(:,3))
            title('bus voltage phase profile')
            xlabel('internal bus number')
            ylabel('phase in degrees')
            disp('paused: press any key to continue')
            pause
        elseif (sel == 7)
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
        elseif (sel == 0)
            flag = 1;
        else
            error('lfdc: invalid plotting selection.');
        end
    end
end

% eof
