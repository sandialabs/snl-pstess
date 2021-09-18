% vsdemo.m
%
% script file to demo vstab solutions
%
% Purpose: Performs a load flow with the load and generation increased
%     from the starting value by a ratio input by the user. Modal analysis
%     of the current loadflow Jacobian may be performed at any stage. Normally
%     this will be close to the point at which the load flow fails to converge.
%     A user can choose to perform a new load flow based on either the current,
%     previous or original voltage profile

%-----------------------------------------------------------------------------%
% Version history
%
% Version: 2.1
% Author:  Ryan Elliott
% Date:    July 2020
% Purpose: Updated to accommodate batch processing mode
%
% Version: 2.0
% Author:  Graham Rogers
% Date:    July 1997
%
% Version: 1.0
% Author:  Graham Rogers
% Date:    August 1994
% Note:    Initial version
%-----------------------------------------------------------------------------%

clear all; close all; clc;

dypar = get_dypar();                % get parameter settings
batch_mode = dypar.batch_mode;      % batch processing mode

load_bus = 3;
gen_bus = 2;
swing_bus = 1;

disp('vstab demo program')
% load input data from m.file
if batch_mode
    [dfile,pathname] = get_path();
else
    [dfile,pathname] = uigetfile('d*.m','Select Data File');
end

% import base case data
if (pathname == 0)
    error('vsdemo: you must select a valid data file.');
else
    full_dfile = fullfile(pathname, dfile);
    run(full_dfile);
    run('import_var');              % set up global variables
end

% check for valid dynamic data file
if isempty(bus)
    error('vsdemo: the selected file is not a valid load flow data file (bus).');
end

n_bus = length(bus(:,1));
g.lfac.gen_chg_idx = ones(n_bus,1);
chg_bus = find(g.lfac.gen_chg_idx);
g_bus_idx = find(bus(:,10) == 2);
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

while (lf_new >= 1)
    if (lf_new == 1)
        bus_sol = bus_old;
        line_sol = line_old;
        v_mag = [];
        pr_old = 1;
        pr_pre = 1;
        lf_num = 0;
    end

    if (lf_new == 3)
        bus_sol = bus_pre;
        line_sol = line_pre;
        plot_pr(lf_num) = [];
        v_mag(:,lf_num) =  [];
        lf_num = lf_num - 1;
        pr_old = 1;
        prat = plot_pr(lf_num)
    end

    lf_new = 2;
    while (lf_new == 2)
        bus_pre = bus_sol;
        line_pre = line_sol;
        pr_pre = prat;
        pr_text = num2str(pr_pre);
        disp(['previous power ratio = ' pr_text])
        pr_old = 1;

        lf_num = lf_num + 1;
        if (lf_num == 1)
            prat = 1;
        else
            if batch_mode
                prat = 1;
            else
                prat = input('Enter the power increase ratio, [1]: ');
                if isempty(prat)
                    prat = 1;  % default
                end
            end
        end

        plot_pr(lf_num) = prat;
        bus_sol(:,4) = bus_sol(:,4) + (prat - pr_pre)*bus_sol(:,4);
        bus_sol(:,6) = bus_sol(:,6) + (prat - pr_pre)*bus_sol(:,6);
        bus_sol(chg_bus,7) = bus_sol(chg_bus,7) ...
                             + (prat - pr_pre)*bus_sol(chg_bus,7);

        if ~isempty(g_pq_idx)
            bus_sol(g_bus_idx(g_pq_idx),7) = bus_sol(g_bus_idx(g_pq_idx),7) ...
                                             + (prat - pr_pre)*gb_qload(g_pq_idx);
        end

        gb_qload = gb_qload + (prat - pr_pre)*gb_qload;
        load_p(:,lf_num) = bus_sol(:,6);
        n_bus = length(bus_sol(:,1));
        n_line = length(line(:,1));

        [bus_sol,line_sol,line_flow] = loadflow(bus_sol,line_sol,1e-8,50,1.0,'n',2);

        chg_bus = find(g.lfac.gen_chg_idx);
        g_pq_idx = find(~g.lfac.gen_chg_idx(g_bus_idx));

        load_q(chg_bus,lf_num) = bus_sol(chg_bus,7);
        if ~isempty(g_pq_idx)
            load_q(g_bus_idx(g_pq_idx),lf_num) = gb_qload(g_pq_idx);
        end

        gen_q(chg_bus,lf_num) = bus_sol(chg_bus,5);

        n_chg_bus = find(~g.lfac.gen_chg_idx);
        if ~isempty(n_chg_bus)
            gen_q(n_chg_bus,lf_num) = -(bus_sol(n_chg_bus,7)-load_q(n_chg_bus,lf_num));
        end

        gen_p(:,lf_num) = bus_sol(:,4);
        v_mag(:,lf_num) = bus_sol(:,2);

        if (lf_num > 1)
            plot(plot_pr,v_mag)
            title('v/p curves')
            xlabel('power ratio')
            ylabel('voltage magnitude pu')
        end

        flag = 0;
        while (flag == 0)
            disp('You can examine the system data')
            disp('Type 1 to see initial bus data')
            disp('     2 to see line data')
            disp('     3 to see solved load flow bus solution')
            disp('     4 to see line flow')
            disp('     5 to see bus voltage magnitude profile')
            disp('     0 to quit')

            if batch_mode
                sel = 0;
            else
                sel = input('Enter your selection, (0--5)[0]: ');
                if isempty(sel)
                    sel = 0;  % default
                end
            end

            if (sel == 1)
                bus
                disp('paused: press any key to continue')
                pause
            elseif (sel == 2)
                line_sol
                disp('paused: press any key to continue')
                pause
            elseif (sel == 3)
                bus_sol
                disp('paused: press any key to continue')
                pause
            elseif (sel == 4)
                line_flow
                disp('paused: press any key to continue')
                pause
            elseif (sel == 5)
                bar(bus_sol(:,2))
                title('bus voltage magnitude profile')
                xlabel('internal bus number')
                ylabel('voltage in pu')
                disp('paused: press any key to continue')
                pause
            elseif (sel == 0)
                flag = 1;
            else
                sel = 0;
                flag = 1;
            end
        end

        disp('Type 1 to do an additional load flow')
        disp('Type 2 to go on to modal analysis')
        disp('Type 0 to quit vsdemo')

        if batch_mode
            alf = 0;
        else
            alf = input('Enter your selection, (0--2)[0]: ');
            if isempty(alf)
                alf = 0;  % default
            end
        end

        if (alf == 0)
            return;
        elseif (alf == 2)
            lf_new = 2;
            break
        else
            if (alf ~= 1)
                error('vsdemo: invalid plotting selection.');
            end

            disp('Type 1 to start from the original bus data')
            disp('Type 2 to start load flow with current bus data')
            disp('Type 3 to start load flow with previous bus data')

            if batch_mode
                lf_new = 2;
            else
                lf_new = input('Enter your selection, (1--3)[2]: ');
                if isempty(lf_new)
                    lf_new = 2;  % default
                end
            end

            if (lf_new == 1)
                plot_pr = [];
                break
            end

            if (lf_new == 3)
                pr_pre = prat;
                break
            end
        end
    end

    if (lf_new == 2)
        % do modal analysis of the load flow jacobian
        % form the sparse Y matrix
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
        g.lfac.bus_type = round(bus_sol(:,10));
        sw_bno = ones(n_bus,1);
        g_bno = sw_bno;

        % set up index for Jacobian calculation
        % form PQV_no and PQ_no
        bus_zeros = zeros(n_bus,1);
        bus_index = [1:1:n_bus]';
        swing_index = find(g.lfac.bus_type == 1);
        sw_bno(swing_index) = bus_zeros(swing_index);

        g.lfac.PQV_no = find(g.lfac.bus_type >= 2);
        g.lfac.PQ_no = find(g.lfac.bus_type == 3);

        gen_index = find(g.lfac.bus_type == 2);
        g_bno(gen_index) = bus_zeros(gen_index);

        % sw_bno is a vector having ones everywhere but the swing bus locations
        % g_bno is a vector having ones everywhere but the generator bus locations
        % construct sparse angle reduction matrix
        il = length(g.lfac.PQV_no);
        ii = [1:1:il]';
        g.lfac.ang_red = sparse(ii,g.lfac.PQV_no,ones(il,1),il,n_bus);

        % construct sparse voltage reduction matrix
        il = length(g.lfac.PQ_no);
        ii = [1:1:il]';
        g.lfac.volt_red = sparse(ii,g.lfac.PQ_no,ones(il,1),il,n_bus);

        % form the system Jacobian
        [J11,J12,J21,J22] = form_jac(V,ang,Y);

        % reduce the Jacobian to voltage terms only
        J22 = J22-J21*inv(J11)*J12;

        % find the eigenvectors and eigenvalues of the reduced inverse Jacobian
        [y,L] = eig(full(inv(J22)));

        % return to the full bus matrices for graphic output
        y = g.lfac.volt_red'*y*g.lfac.volt_red;
        L = g.lfac.volt_red'*diag(L);
        [my,iy] = max(abs((L)));
        ld_str = num2str(L(iy));
        disp(['the dominant eigenvalue ',ld_str])

        [lmax,ib] = max(abs(y(:,iy)));
        lm_str = num2str(lmax);
        ib_str = num2str(bus(ib,1));
        disp(['the maximum eigenvector entry is ',lm_str])
        disp(['the corresponding bus number is ',ib_str])

        if batch_mode
            eig_plot = 'n';
        else
            eig_plot = input('Do you want to plot the eigenvector, (y/n)[n]: ','s');
            eig_plot = lower(eig_plot);
            if isempty(eig_plot)
                eig_plot = 'n';  % default
            end
        end

        if strcmp(eig_plot,'y')
            bar(abs(y(:,iy)));
            title('critical eigenvector')
            xlabel('internal bus number')
            ylabel('magnitude of eigenvector')
        end

        disp('Type 1 to perform an additional load flow starting with the original bus data')
        disp('Type 2 to perform an additional load flow starting with the current bus data')
        disp('Type 3 to perform an additional load flow starting with the previous bus data')
        disp('Type 0 to exit demo')

        if batch_mode
            lf_new = 0;
        else
            lf_new = input('Enter your selection, (0--3)[0]: ');
            if isempty(lf_new)
                lf_new = 0;  % default
            end
        end

        if (lf_new == 0)
            return
        end
    end
end

% eof
