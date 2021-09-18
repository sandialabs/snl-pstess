% lfdemo.m
% m file to demo loadflow solutions
% a loadflow report is placed in lf_rep.txt if requested

%-----------------------------------------------------------------------------%

clear all; close all; clc;

dypar = get_dypar();                % get parameter settings
batch_mode = dypar.batch_mode;      % batch processing mode

disp('loadflow demo program')
% input data file
if batch_mode
    [dfile,pathname] = get_path();
else
    [dfile,pathname] = uigetfile('d*.m','Select Data File');
end

% import base case data
if (pathname == 0)
    error('lfdemo: you must select a valid data file.');
else
    full_dfile = fullfile(pathname, dfile);
    run(full_dfile);
    run('import_var');              % set up global variables
end

% check for valid dynamic data file
if isempty(bus)
    error('lfdemo: the selected file is not a valid data file (bus).');
end

if batch_mode
    answer = 'y';
else
    answer = input('Would you like a load flow solution report, (y/n)[n]: ','s');
    answer = lower(answer);
    if isempty(answer)
        answer = 'y';  % default
    end
end

% write solution report
if strcmp(answer,'y')
    diary lf_rep.txt;
end

n_bus = length(bus(:,1));
g.lfac.gen_chg_idx = ones(n_bus,1);
[bus_sol,line_sol,line_flow] = loadflow(bus,line,1e-8,50,1.0,answer,2);

% close the solution report
if strcmp(answer,'y')
    diary off
    disp('lfdemo: a load flow study report is available in lf_rep.txt.')
end

flag = 0;
while (flag == 0)
    disp('You can examine the system data')
    disp('Type 1 to see initial bus data')
    disp('     2 to see modified line data')
    disp('     3 to see solved load flow bus solution')
    disp('     4 to see line flow')
    disp('     5 to see bus voltage magnitude profile')
    disp('     6 to see bus voltage phase profile')
    disp('     0 to quit')

    if batch_mode
        sel = 0;      % quit without displaying solution
    else
        sel = input('Enter your selection, (0--6)[0]: ');
        if isempty(sel);
            sel = 0;  % quit without displaying solution (default)
        end
    end

    if (sel == 1)
        disp('                                   Initial Bus Data')
        disp('                                      GENERATION             LOAD')
        disp('       BUS     VOLTS     ANGLE      REAL  REACTIVE      REAL  REACTIVE ')
        disp(bus(:,1:7))
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
    elseif (sel == 0)
        flag = 1;
    else
        error('lfdemo: invalid plotting selection.');
    end
end

% eof
