function [bus_sol,line_sol,line_flow,rec_par,inv_par,line_par] = lfdcs(bus,line,dci_dc,dcr_dc)
% Syntax: [bus_sol,line_sol,line_flow,rec_par,inv_par,line_par] = ...
%            lfdcs(bus,line,dci_dc,dcr_dc)
%
% Purpose:   Solves load flow with one or more dc lines
% Inputs:    bus - ac bus specification matrix
%            line - ac line specification matrix
% Outputs:   bus_sol - solved ac bus specification file
%            line_sol - solved ac line specification file
%            line_flow - calculated flow on ac lines
%            rec_par  - rectifier parameters
%                     - rec_par(:,1) firing angle degrees
%                     - rec_par(:,2) dc voltage kV
%                     - rec_par(:,3) dc power MW
%                     - rec_par(:,4) equi HT bus voltage kV
%            inv_par  - inverterer parameters
%                     - inv_par(:,1) extinction angle degrees
%                     - inv_par(:,2) dc voltage kV
%                     - inv_par(:,3) dc power MW
%                     - inv_par(:,4) equi HT bus voltage kV
%            line_par - dc line current kA
%
% Calls:     loadflow, dc_lf
%
% Called by: s_simu
%
% Algorithm: iterates between AC and DC solutions until
%            converged solution is reached

%-----------------------------------------------------------------------------%
% Version history
%
% Version:   1.0
% Author:    Graham Rogers
% Date:      February 1997
%-----------------------------------------------------------------------------%

global g;  % declaring struct of global variables

disp('load flow with hvdc')

% perform load flow iterations
errv = 0;
itermax = 50;
iter = 0;
bus_old = bus;

while ((errv == 0) && (iter <= itermax))
    iter = iter + 1;

    [bus_sol,line_sol,line_flow] = loadflow(bus,line,1e-6,50,1.0,'n',2);

    if (iter == 1)
        dc_indx(line,dci_dc,dcr_dc);
    end

    % perform dc load flow
    [rec_par,inv_par,line_par,tap,Sr,Si] = dc_lf(bus_sol,line_sol,dci_dc,dcr_dc);

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
        [bus_sol,line_sol,line_flow] = loadflow(bus,line,1e-6,50,1.0,'n',2);

        % reperform dc load flow
        [rec_par,inv_par,line_par,tap,Sr,Si] = dc_lf(bus_sol,line_sol,dci_dc,dcr_dc);

        % set dc quantities based on final ac load flow
        P = bus_sol(g.dc.ac_bus,6);
        Q = bus_sol(g.dc.ac_bus,7);
        Vac = bus_sol(g.dc.ac_bus,2);
        Vang = bus_sol(g.dc.ac_bus,3)*pi/180;
        i_ac = (P - 1j*Q)./Vac./exp(-1j*Vang);
        Vac = Vac.*bus_sol(g.dc.ac_bus,13);  % convert to kV

        % convert ac currents to dc
        idceq = abs(i_ac)*pi*g.sys.basmva/3/sqrt(2) ...
                ./bus_sol(g.dc.ac_bus,13)./g.dc.dcsp_con(:,6);

        % convert ac currents to kA
        i_ac = i_ac*g.sys.basmva/sqrt(3)./bus_sol(g.dc.ac_bus,13);

        % xequ -- eq transformer reactance
        % VHT -- equivalent HT bus voltage
        % Vdo -- ideal dc voltages
        xequ = sqrt(3)*g.dc.dcsp_con(:,5)./g.dc.dcsp_con(:,6);
        VHT = abs(Vac.*exp(1j*Vang) + 1j*xequ.*i_ac);
        Vdo = 3*sqrt(2)*VHT.*g.dc.dcsp_con(:,6)/pi;

        % calculate dc voltages
        dc_ang(g.dc.r_idx,1) = rec_par(:,1)*pi/180;
        dc_ang(g.dc.i_idx,1) = inv_par(:,1)*pi/180;

        % Rc is in series as far as dc is concerned
        Rc = g.dc.dcsp_con(:,6).*g.dc.dcsp_con(:,5)*3/pi;
        g.dc.Vdc = Vdo.*cos(dc_ang) - Rc.*idceq;
    end
end

if (iter > itermax)
    estr = '\nlfdcs: dc load flow failed to converge in %0.0f iterations.';
    error(sprintf(estr,itermax));
else
    % dc and ac converged
    dstr = '\nlfdcs: solution converged with %0.0f dc load flow iterations.';
    disp(sprintf(dstr,iter));
end

end  % function end

% eof
