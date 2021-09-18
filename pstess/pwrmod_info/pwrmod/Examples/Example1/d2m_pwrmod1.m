% 2-machine PST simulation with pwrmod control on buses 2 and 3. 
% D. Trudnowski, Sp. 2015

bus = [ ...
% num volt  angle p_gen q_gen p_load q_load G_shunt B_shunt type q_max q_min v_rated v_max v_min
  1   1.0   0.0   0.12  0     0      0      0       0       1    100   -100  20      1.5   0.5;
  2   0.98  0.0   0.005 0.001 0      0      0       0       2    1     -1    20      1.5   0.5; %Solar gen
  3   0.97  0.0   0.05  0.01  0      0       0      0       2    1     -1    20      1.5   0.5; %Solar gen
  4   1.0   0.0   0.15  0     0      0      0       0       2    100   -100  20      1.5   0.5
  5   1.0   0.0   0     0     0.275  0.2    0       0       3    0     0     20      1.5   0.5;
  6   1.0   0.0   0     0     0      0      0       0       3    0     0     20      1.5   0.5];

line = [ ...
% bus bus r    x    y    tapratio tapphase tapmax tapmin tapsize
  4   3   0.0  0.1  0.0  1.0      0        0      0      0; %transformer
  1   6   0.0  0.75 0.0  1.0      0        0      0      0; %line
  5   6   0.0  0.75 0.0  1.0      0        0      0      0; %line
  5   3   0.0  1.5  0.0  1.0      0        0      0      0; %line
  2   3   0.0  1.0  0.0  1.0      0        0      0      0; %line
  1   2   0.0  5.0  0.0  1.0      0        0      0      0]; %line

%Both generator parameters are from the example machine in chap 4 of
%Anderson and Fouad with a salient pole assumption.
%This model is a sub-transient model.
%From eqn (4.289), T"_qo=T'_qo if it is a 2-axis transient model.
mac_con = [ ...
% num bus base  x_l  r_a    x_d x'_d   x"_d  T'_do T"_do x_q   x'_q  x"_q  T'_qo T"_qo H      d_0  
   1   1  100   0.15 0.0011 1.7 0.245  0.185 5.9   0.03  1.64  1.64  0.185 0.082 0.082 2.37   0   0   1;
   2   4  100   0.15 0.0011 1.7 0.245  0.185 5.9   0.03  1.64  1.64  0.185 0.082 0.082 2.37   0   0   4];
  
%Exciter
%From p. 1137 of Kundur
exc_con = [...
%  type mach  T_R   K_A   T_A   T_B   T_C   Vmax  Vmin
    0    1    0     212   0.01  12    1     7.5   -6.0; %Fast static exciter, with TGR
    0    2    0     212   0.01  12    1     7.5   -6.0]; %Fast static exciter, with TGR

%PSS
%Designed using the "large-inertia, infinite-bus" method in Roger's book and
%Kundur's book.
%type gen# K  Tw T1  T2   T3  T4   max min
% pss_con = [ ...
%   1   1    80 10 0.4 0.04 0.4 0.06 0.1 -0.1;
%   1   2    80 10 0.4 0.04 0.4 0.06 0.1 -0.1];

tg_con = [...
%     mach wf 1/R     Tmax   Ts     Tc     T3     T4     T5
   1   1   1  20.0    1.00   0.40   75.0   10.0  -2.4    1.2];%hydro

% non-conforming load
% col 1       bus number
% col 2       fraction const active power load
% col 3       fraction const reactive power load
% col 4       fraction const active current load
% col 5       fraction const reactive current load
load_con = [...
%bus Pcont Qconst P_Icont Q_Iconst
  2   1     1     0       0;  %Modulation bus
  3   1     1     0       0;  %Modulation bus
  5   0     0     0.1     0]; %Load bus

% Power modulation data (sets up modulation of real and reac power injection)
% col 1       bus number
% col 2       real-power time constant (Tp)
% col 3       max real-power modulation (pu on system base)
% col 4       min real-power modulation (pu on system base)
% col 5       reac-power time constant (Tq)
% col 6       max reac-power modulation (pu on system base)
% col 7       min reac-power modulation (pu on system base)
% NOTE: This creates b_pwrmod_p and b_pwrmod_q for the linear analysis.
%       b_pwrmod_*(:,j) corresponds to row j of pwrmod_con.
pwrmod_con=[...
%bus  T    Pmax Pmin  Tq   Qmax  Qmin
  2   0.05 1    -1    0.05  1     -1; 
  3   0.05 1    -1    0.05  1     -1];

% Switching
%Switching file defines the simulation control
% row 1 col1  simulation start time (s) (cols 2 to 6 zeros)
%     col7  initial time step (s)
% row 2 col1  fault application time (s)
%     col2  bus number at which fault is applied
%     col3  bus number defining far end of faulted line
%     col4  zero sequence impedance in pu on system base
%     col5  negative sequence impedance in pu on system base
%     col6  type of fault  - 0 three phase
%                  - 1 line to ground
%                  - 2 line-to-line to ground
%                  - 3 line-to-line
%                  - 4 loss of line with no fault
%                  - 5 loss of load at bus
%     col7  time step for fault period (s)
% row 3 col1  near end fault clearing time (s) (cols 2 to 6 zeros)
%     col7  time step for second part of fault (s)
% row 4 col1  far end fault clearing time (s) (cols 2 to 6 zeros)
%     col7  time step for fault cleared simulation (s)
% row 5 col1  time to change step length (s)
%     col7  time step (s)
sw_con = [...
0        0    0    0    0    0    1/120;%sets intitial time step
5.0      6    1    0    0    5    1/120; %no fault
5+1/60   0    0    0    0    0    1/120; %clear near end of fault
5+2/60   0    0    0    0    0    1/120; %clear far end of fault
5.1      0    0    0    0    0    1/120; % increase time step
10       0    0    0    0    0    1/120]; % end simulation