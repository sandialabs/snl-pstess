function [bus_sol,line_sol,line_flow,Jac] = loadflow(bus,line,tol,itermax,acc,display,flag)
% Syntax:  [bus_sol,line_sol,line_flow,Jac] =
%          loadflow(bus,line,tol,itermax,acc,display,flag) OR
%          [bus_sol,line_sol,line_flow] =
%          loadflow(bus,line,tol,itermax,acc,display,flag)
%
% Purpose: Solve the load-flow equations of power systems
%          modified to eliminate do loops and improve the use
%          sparse matices
%
% Input:   bus       - bus data
%          line      - line data
%          tol       - tolerance for convergence
%          itermax   - maximum number of iterations
%          acc       - acceleration factor
%          display   - 'y', generate load-flow study report
%                       else, no load-flow study report
%          flag      - 1, form new Jacobian every iteration
%                      2, form new Jacobian every other iteration
%
% Output:  bus_sol   - bus solution (see report for the solution format)
%          line_sol  - modified line matrix
%          line_flow - line flow solution (see report)
%          Jac       - Jacobian matrix
%
% Algorithm: Newton-Raphson method using the polar form of the
%            equations for P (real power) and Q (reactive power).
%
% Calls: y_sparse, calc, form_jac, chq_lim

%-----------------------------------------------------------------------------%
% Version history
%
% Version: 2.2
% Author:  Graham Rogers
% Date:    November 1997
% Purpose: Modification to correct generator var error on output
%
% Version: 2.1
% Author:  Graham Rogers
% Date:    October 1996
% Purpose: To add generator var limits and on-load tap changers
%
% Version: 2.0
% Author:  Graham Rogers
% Date:    March 1994
%
% Version: 1.0
% Authors: Kwok W. Cheung, Joe H. Chow
% Date:    March 1991
%-----------------------------------------------------------------------------%

global g;  % declaring struct of global variables

tt = clock;     % start the total time clock
load_bus = 3;
gen_bus = 2;
swing_bus = 1;

% default Jacobian update flag
if (nargin < 7)
    flag = 1;
elseif ((flag ~= 1) && (flag ~= 2))
    error('loadflow: flag not recognized.');
end

% default display option
if (nargin < 6)
    display = 'n';
end

% default acceleration factor
if (nargin < 5)
    acc = 1.0;
end

% default maximum number of iterations
if (nargin < 4)
    itermax = 50;
end

% default solution tolerance
if (nargin < 3)
    tol = 1e-8;
end

% either bus or line not supplied, cannot compute
if (nargin < 2)
    error('loadflow: bus and/or line not supplied.');
end

[n_line, n_lc] = size(line);  % number of lines and no of line col
[n_bus, n_col] = size(bus);   % number of buses and number of col

% set default bus data
if (n_col < 15)
    % set generator var limits
    if (n_col < 12)
        bus(:,11) = 9999*ones(n_bus,1);
        bus(:,12) = -9999*ones(n_bus,1);
    end

    % populating voltage rating, if necessary
    if (n_col < 13)
        bus(:,13) = ones(n_bus,1);
    end

    bus(:,14) = 1.5*ones(n_bus,1);
    bus(:,15) = 0.5*ones(n_bus,1);
    volt_min = bus(:,15);
    volt_max = bus(:,14);
else
    volt_min = bus(:,15);
    volt_max = bus(:,14);
end

no_vmin_idx = find(volt_min == 0);
if ~isempty(no_vmin_idx)
    volt_min(no_vmin_idx) = 0.5*ones(length(no_vmin_idx),1);
end

no_vmax_idx = find(volt_max == 0);
if ~isempty(no_vmax_idx)
    volt_max(no_vmax_idx) = 1.5*ones(length(no_vmax_idx),1);
end

no_mxv = find(bus(:,11) == 0);
if ~isempty(no_mxv)
    bus(no_mxv,11) = 9999*ones(length(no_mxv),1);
end

no_mnv = find(bus(:,12) == 0);
if ~isempty(no_mnv)
    bus(no_mnv,12) = -9999*ones(length(no_mnv),1);
end

no_vrate = find(bus(:,13) == 0);
if ~isempty(no_vrate)
    bus(no_vrate,13) = ones(length(no_vrate),1);
end

tap_it = 0;
tap_it_max = 10;
no_taps = 0;

% line data defaults, sets all tap ranges to zero - this fixes taps
if (n_lc < 10)
    line(:,8:10) = zeros(n_line,3);   % Joe Chow 12/15/2015 - to prevent
                                      % deletion of phase shift angle
    % line(:,7:10) = zeros(n_line,4);
    no_taps = 1;                      % disable tap changing
end

% outer loop for on-load tap changers
mm_chk = 1;
while ((tap_it <= tap_it_max) && mm_chk)
    tap_it = tap_it + 1;

    % build admittance matrix Y
    [Y,nSW,nPV,nPQ,SB] = y_sparse(bus,line);

    % process bus data
    bus_no = bus(:,1);
    V = bus(:,2);
    ang = bus(:,3)*pi/180;
    Pg = bus(:,4);
    Qg = bus(:,5);
    Pl = bus(:,6);
    Ql = bus(:,7);
    Gb = bus(:,8);
    Bb = bus(:,9);
    g.lfac.bus_type = round(bus(:,10));
    qg_max = bus(:,11);
    qg_min = bus(:,12);
    sw_bno = ones(n_bus,1);
    g_bno = sw_bno;

    % set up index for Jacobian calculation, form PQV_no and PQ_no
    % sw_bno is a vector having ones everywhere but the swing bus locations
    % g_bno is a vector having ones everywhere but the generator bus locations
    bus_zeros = zeros(n_bus,1);
    swing_index = find(g.lfac.bus_type == 1);
    sw_bno(swing_index) = bus_zeros(swing_index);

    g.lfac.PQV_no = find(g.lfac.bus_type >= 2);
    g.lfac.PQ_no = find(g.lfac.bus_type == 3);

    gen_index = find(g.lfac.bus_type == 2);
    g_bno(gen_index) = bus_zeros(gen_index);

    % construct sparse angle reduction matrix
    il = length(g.lfac.PQV_no);
    ii = (1:1:il)';
    g.lfac.ang_red = sparse(ii,g.lfac.PQV_no,ones(il,1),il,n_bus);

    % construct sparse voltage reduction matrix
    il = length(g.lfac.PQ_no);
    ii = (1:1:il)';
    g.lfac.volt_red = sparse(ii,g.lfac.PQ_no,ones(il,1),il,n_bus);

    iter = 0;  % initialize iteration counter

    % calculate the power mismatch and check convergence
    [delP,delQ,P,Q,conv_flag] = calc(V,ang,Y,Pg,Qg,Pl,Ql,sw_bno,g_bno,tol);

    % start iteration process for main Newton_Raphson solution
    st = clock;  % start the iteration time clock
    while ((conv_flag == 1) && (iter <= itermax))
        iter = iter + 1;

        % form the Jacobian matrix
        Jac_iter = form_jac(V,ang,Y);

        % reduced real and reactive power mismatch vectors
        red_delP = g.lfac.ang_red*delP;
        red_delQ = g.lfac.volt_red*delQ;

        % solve for voltage magnitude and phase angle increments
        temp = Jac_iter\[red_delP; red_delQ];

        % expand solution vectors to all buses
        delAng = g.lfac.ang_red'*temp(1:length(g.lfac.PQV_no),:);

        delV_beg = length(g.lfac.PQV_no) + 1;
        delV_end = length(g.lfac.PQV_no) + length(g.lfac.PQ_no);
        delV = g.lfac.volt_red'*temp(delV_beg:delV_end,:);

        % update voltage magnitude and phase angle
        V = V + acc*delV;
        V = max(V,volt_min);  % voltage higher than minimum
        V = min(V,volt_max);  % voltage lower than maximum
        ang = ang + acc*delAng;

        % calculate the power mismatch and check convergence
        [delP,delQ,P,Q,conv_flag] = calc(V,ang,Y,Pg,Qg,Pl,Ql,sw_bno,g_bno,tol);

        % check if Qg is outside limits
        gen_index = find(g.lfac.bus_type == 2);
        Qg(gen_index) = Q(gen_index) + Ql(gen_index);

        [Qg,Ql,g_bno,lim_flag] = chq_lim(Qg,Ql,g_bno,qg_max,qg_min);
        if (lim_flag == 1)
            warning(sprintf('\nloadflow: Qg at var limit.'));
        end
    end

    if (iter > itermax)
        wstr = '\nloadflow: inner ac load flow failed to converge after %0.0f ';
        wstr = [wstr, 'iterations at tap iteration %0.0f.'];
        warning(sprintf(wstr,[itermax,tap_it]));
    else
        disp(sprintf('loadflow: inner load flow iterations = %0.0f',iter));
    end

    if (no_taps == 0)
        run('lftap');
    else
        mm_chk = 0;
    end
end

if (tap_it > tap_it_max)
    wstr = '\nloadflow: tap iteration failed to converge after %0.0f iterations.';
    warning(sprintf(wstr,tap_it_max));
else
    disp(sprintf('loadflow: tap iterations = %0.0f',tap_it));
end

ste = clock;  % end the iteration time clock

vmx_idx = find(V == volt_max);
if ~isempty(vmx_idx)
    wstr = '\nloadflow: the voltage at bus %0.0f is at the maximum limit.';
    warning(sprintf(wstr,bus(vmx_idx,1)));
end

vmn_idx = find(V == volt_min);
if ~isempty(vmn_idx)
    wstr = '\nloadflow: the voltage at bus %0.0f is at the minimum limit.';
    warning(sprintf(wstr,bus(vmn_idx,1)));
end

gen_index = find(g.lfac.bus_type == 2);
load_index = find(g.lfac.bus_type == 3);

Pg(gen_index) = P(gen_index) + Pl(gen_index);
Qg(gen_index) = Q(gen_index) + Ql(gen_index);

gend_idx = find((bus(:,10) == 2) & (g.lfac.bus_type ~= 2));
if ~isempty(gend_idx)
    wstr = '\nloadflow: the generator at bus %0.0f is at its var limits, ';
    wstr = [wstr, 'Qg = %0.2f.'];
    warning(sprintf(wstr,[bus(gend_idx,1),Qg(gend_idx)].'));
    %
    Qlg = Ql(gend_idx) - bus(gend_idx,7); % the generator var part of the load
    Qg(gend_idx) = Qg(gend_idx) - Qlg;    % restore the generator vars
    Ql(gend_idx) = bus(gend_idx,7);       % restore the original load vars
end

Pl(load_index) = Pg(load_index) - P(load_index);
Ql(load_index) = Qg(load_index) - Q(load_index);

Pg(SB) = P(SB) + Pl(SB);
Qg(SB) = Q(SB) + Ql(SB);
VV = V.*exp(1j*ang);                      % solution voltage

% calculate the line flows and power losses
tap_index = find(abs(line(:,6)) > 0);
tap_ratio = ones(n_line,1);
tap_ratio(tap_index) = line(tap_index,6);
phase_shift(:,1) = line(:,7);
tps = tap_ratio.*exp(1j*phase_shift*pi/180);
from_bus = line(:,1);
from_int = g.bus.bus_int(round(from_bus));
to_bus = line(:,2);
to_int = g.bus.bus_int(round(to_bus));
r = line(:,3);
rx = line(:,4);
chrg = line(:,5);
z = r + 1j*rx;
y = ones(n_line,1)./z;

MW_s = VV(from_int).*conj((VV(from_int) - tps.*VV(to_int)).*y ...
                          + VV(from_int).*(1j*chrg/2))./(tps.*conj(tps));

P_s = real(MW_s);     % active power sent out by from_bus to to_bus
Q_s = imag(MW_s);     % reactive power sent out by from_bus to to_bus

MW_r = VV(to_int).*conj((VV(to_int) - VV(from_int)./tps).*y ...
                        + VV(to_int).*(1j*chrg/2));

P_r = real(MW_r);     % active power received by to_bus from from_bus
Q_r = imag(MW_r);     % reactive power received by to_bus from from_bus

iline = (1:1:n_line)';
line_ffrom = [iline, from_bus, to_bus, P_s, Q_s];
line_fto = [iline, to_bus, from_bus, P_r, Q_r];

P_loss = sum(P_s) + sum(P_r);
Q_loss = sum(Q_s) + sum(Q_r);

bus_sol = [bus_no, V, ang*180/pi, Pg, Qg, Pl, Ql, Gb, Bb, ...
           g.lfac.bus_type, qg_max, qg_min, bus(:,13), volt_max, volt_min];

line_sol = line;
line_flow(1:n_line,:) = [iline, from_bus, to_bus, P_s, Q_s];
line_flow(1+n_line:2*n_line,:) = [iline, to_bus, from_bus, P_r, Q_r];

% halt execution upon non-convergence
if (conv_flag == 1)
    error('loadflow: ac load flow failed to converge.');
end

% display results
if strcmp(display,'y')
    clc
    disp('                             LOAD-FLOW STUDY')
    disp('                    REPORT OF POWER FLOW CALCULATIONS ')
    disp(' ')
    disp(date)
    fprintf('SWING BUS                  : BUS %g \n', SB)
    fprintf('NUMBER OF ITERATIONS       : %g \n', iter)
    fprintf('SOLUTION TIME              : %g sec.\n', etime(ste,st))
    fprintf('TOTAL TIME                 : %g sec.\n', etime(clock,tt))
    fprintf('TOTAL REAL POWER LOSSES    : %g.\n', P_loss)
    fprintf('TOTAL REACTIVE POWER LOSSES: %g.\n\n', Q_loss)

    if (conv_flag == 0)
        disp('                                      GENERATION             LOAD')
        disp('       BUS     VOLTS     ANGLE      REAL  REACTIVE      REAL  REACTIVE ')
        disp(bus_sol(:,1:7))

        disp('                      LINE FLOWS                     ')
        disp('      LINE  FROM BUS    TO BUS      REAL  REACTIVE   ')
        disp(line_ffrom)
        disp(line_fto)
    end
end

if (iter > itermax)
    wstr = '\nloadflow: solution did not converge in %g iterations.';
    warning(sprintf(wstr,itermax));
end

% outputting Jacobian matrix for voltage stability studies
if (nargout > 3)
    Jac = Jac_iter;
end

end  % function end

% eof
