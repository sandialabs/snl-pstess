% A 3-machine 9-bus system from Chow's book pp.70
% data3m9b.m

% bus data format
% bus: number, voltage(pu), angle(degree), p_gen(pu), q_gen(pu),
%      p_load(pu), q_load(pu),G shunt,B shunt, bus_type
%      bus_type - 1, swing bus
%               - 2, generator bus (PV bus)
%               - 3, load bus (PQ bus)

bus = [ 1 1.04    0.00   0.00   0.00  0.00  0.00  0.00  0.00 1;
	2 1.02533 0.00   1.63   0.00  0.00  0.00  0.00  0.00 2;
	3 1.02536 0.00   0.85   0.00  0.00  0.00  0.00  0.00 2;
	4 1.00    0.00   0.00   0.00  0.00  0.00  0.00  0.00 3;
	5 1.00    0.00   0.00   0.00  0.90  0.30  0.00  0.00 3;
	6 1.00    0.00   0.00   0.00  0.00  0.00  0.00  0.00 3;
	7 1.00    0.00   0.00   0.00  1.00  0.35  0.00  0.00 3;
	8 1.00    0.00   0.00   0.00  0.00  0.00  0.00  0.00 3;
	9 1.00    0.00   0.00   0.00  1.25  0.50  0.00  0.00 3];

% line data format
% line: from bus, to bus, resistance(pu), reactance(pu),
%       line charging(pu), tap ratio

line = [1 4 0.0    0.0576 0.     1. 0. ;
	4 5 0.017  0.092  0.158  1. 0. ;
	5 6 0.039  0.17   0.358  1. 0. ;
	3 6 0.0    0.0586 0.     1. 0. ;
	6 7 0.0119 0.1008 0.209  1. 0. ;
	7 8 0.0085 0.072  0.149  1. 0. ;
	8 2 0.0    0.0625 0.     1. 0. ;
	8 9 0.032  0.161  0.306  1. 0. ;
	9 4 0.01   0.085  0.176  1. 0. ];


% Machine data format
% Machine data format
%       1. machine number,
%       2. bus number,
%       3. base mva,
%       4. leakage reactance x_l(pu),
%       5. resistance r_a(pu),
%       6. d-axis sychronous reactance x_d(pu),
%       7. d-axis transient reactance x'_d(pu),
%       8. d-axis subtransient reactance x"_d(pu),
%       9. d-axis open-circuit time constant T'_do(sec),
%      10. d-axis open-circuit subtransient time constant
%                T"_do(sec),
%      11. q-axis sychronous reactance x_q(pu),
%      12. q-axis transient reactance x'_q(pu),
%      13. q-axis subtransient reactance x"_q(pu),
%      14. q-axis open-circuit time constant T'_qo(sec),
%      15. q-axis open circuit subtransient time constant
%                T"_qo(sec),
%      16. inertia constant H(sec),
%      17. damping coefficient d_o(pu),
%      18. dampling coefficient d_1(pu),
%      19. bus number
%
% note: all the following machines use electro-mechanical model
mac_con = [ ...
1 1 100 0.000  0.000  0.  0.0608  0 0  0 0 0 0 0 0  13.64  9.6 0  1;
2 2 100 0.000  0.000  0.  0.1198  0 0  0 0 0 0 0 0   6.4   2.5 0  2;
3 3 100 0.000  0.000  0.  0.1813  0 0  0 0 0 0 0 0   3.01  1.0 0  3];
% induction motor data
% 1. Motor Number
% 2. Bus Number
% 3. Motor MVA Base
% 4. rs pu
% 5. xs pu - stator leakage reactance
% 6. Xm pu - magnetizing reactance
% 7. rr pu 
% 8. xr pu - rotor leakage reactance
% 9. H  s  - motor plus load inertia constant
% 15. fraction of bus power drawn by motor ( if zero motor statrts at t=0)
ind_con = [ ...
  1  7  25. .001 .01 3 .009 .01  .7 0 0 0 0 0 0
  2  9  25. .001 .01 3 .009 .01  .7 0 0 0 0 0 .15
];
% Motor Load Data 
% format for motor load data - mld_con
% 1 motor number
% 2 bus number
% 3 stiction load pu on motor base (f1)
% 4 stiction load coefficient (i1)
% 5 external load  pu on motor base(f2)
% 6 external load coefficient (i2)
% 
% load has the form
% tload = f1*slip^i1 + f2*(1-slip)^i2
mld_con = [ ...
1  7  .1  1  .7  2
2  9  .1  1  .7  2
];

%Switching file defines the simulation control
% row 1 col1  simulation start time (s) (cols 2 to 6 zeros)
%       col7  initial time step (s)
% row 2 col1  fault application time (s)
%       col2  bus number at which fault is applied
%       col3  bus number defining far end of faulted line
%       col4  zero sequence impedance in pu on system base
%       col5  negative sequence impedance in pu on system base
%       col6  type of fault  - 0 three phase
%                            - 1 line to ground
%                            - 2 line-to-line to ground
%                            - 3 line-to-line
%                            - 4 loss of line with no fault
%                            - 5 loss of load at bus
%       col7  time step for fault period (s)
% row 3 col1  near end fault clearing time (s) (cols 2 to 6 zeros)
%       col7  time step for second part of fault (s)
% row 4 col1  far end fault clearing time (s) (cols 2 to 6 zeros)
%       col7  time step for fault cleared simulation (s)
% row 5 col1  time to change step length (s)
%       col7  time step (s)
%
%
%
% row n col1 finishing time (s)  (n indicates that intermediate rows may be inserted)

sw_con = [...
0    0    0    0    0    0    0.0025;%sets intitial time step
2.0  4    5  0    0    0    0.0025; %apply three phase fault at bus 4, on line 4-5
2.05  0    0    0    0    0    0.0025; %clear fault at bus 4
2.10 0    0    0    0    0    0.0025; %clear remote end
4.0 0    0    0    0    0    0.01; % increase time step 
5.0  0    0    0    0    0    0.01; % increase time step
6.0  0    0    0    0    0    0]; % end simulation

