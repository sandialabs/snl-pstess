%-----------------------------------------------------------------------------%
% Two Area Test Case
% sub transient generators with static exciters, turbine/governors
% 50% constant current active loads
% load modulation
% with power system stabilizers
%-----------------------------------------------------------------------------%

disp('Two-area test case with subtransient generator models')
disp('Static exciters')
disp('turbine/governors')
% bus data format
% bus:
% col1 number
% col2 voltage magnitude(pu)
% col3 voltage angle(degree)
% col4 p_gen(pu)
% col5 q_gen(pu),
% col6 p_load(pu)
% col7 q_load(pu)
% col8 G shunt(pu)
% col9 B shunt(pu)
% col10 bus_type
%       bus_type - 1, swing bus
%                - 2, generator bus (PV bus)
%                - 3, load bus (PQ bus)
% col11 q_gen_max(pu)
% col12 q_gen_min(pu)
% col13 v_rated (kV)
% col14 v_max  pu
% col15 v_min  pu

% 171004 - Ryan's note: Changed the power flow due to new inertia distribution
% Bus 2 is now the slack bus.

bus = [...
%   num  volt   angle      p_gen    q_gen p_load     q_load G_shunt B_shunt type q_max q_min v_rated v_max v_min
    1  1.0100   -7.8436   1.6800     0.2340       0        0  0.00   0.00       2  1.33  -2.5   22.0  1.10  0.90;
    2  1.0100         0  12.3343+1.5 1.9995       0        0  0.00   0.00       1  9.73  -2.5   22.0  1.10  0.90;
    3  0.9978  -19.9001        0        0         0        0  0.00   3.70       3  0.0    0.0  230.0  1.50  0.50;
    4  1.0174  -22.5865        0        0    9.7600   1.0000  0.00   0.00       3  0.0    0.0  115.0  1.05  0.95;
   10  1.0065   -9.4251        0        0         0        0  0.00   0.00       3  0.0    0.0  230.0  1.50  0.50;
   11  1.0100  -32.8705   1.6800     0.3071       0        0  0.00   0.00       2  1.33  -2.5   22.0  1.10  0.90;
   12  1.0100  -25.0255  12.3000-1.5 2.1789       0        0  0.00   0.00       2  9.73  -2.5   22.0  1.10  0.90;
   13  0.9923  -44.9639        0        0         0        0  0.00   6.10       3  0.0    0.0  230.0  1.50  0.50;
   14  1.0156  -49.8303        0        0   17.6500   1.0000  0.00   0.00       3  0.0    0.0  115.0  1.05  0.95;
   20  0.9980  -11.7916        0        0         0        0  0.00   0.00       3  0.0    0.0  230.0  1.50  0.50;
  101  1.0177  -32.6699        0        0         0        0  0.00   1.40       3  0.0    0.0  230.0  1.50  0.50;
  110  1.0053  -34.4538        0        0         0        0  0.00   0.00       3  0.0    0.0  230.0  1.50  0.50;
  120  0.9950  -36.8200        0        0         0        0  0.00   0.00       3  0.0    0.0  230.0  1.50  0.50;
  3101 0.9978  -19.9001        0        0         0        0  0.00   0.00       3  0.0    0.0  230.0  1.50  0.50];

% line data format
% line:
%      col1     from bus
%      col2     to bus
%      col3     resistance(pu)
%      col4     reactance(pu)
%      col5     line charging(pu)
%      col6     tap ratio
%      col7     tap phase
%      col8     tapmax
%      col9     tapmin
%      col10    tapsize

line = [...
% f_bus t_bus   r           x      b    tapratio tapphase tapmax tapmin tapsize
    1   10      0.0     0.0167   0.00    1.0     0.0        0.0     0.0  0.0;
    2   20      0.0     0.0167   0.00    1.0     0.0        0.0     0.0  0.0;
    3    4      0.0     0.005    0.00    0.975   0.0        1.2     0.8  0.00625;
    3   20      0.001   0.0100   0.0175  1.0     0.0        0.0     0.0  0.0;
    3   101     0.011   0.110    0.1925  1.0     0.0        0.0     0.0  0.0;
    3   3101    0.000   0.001    0.0000  1.0     0.0        0.0     0.0  0.0;
   3101 101     0.011   0.109    0.1925  1.0     0.0        0.0     0.0  0.0;
    10  20      0.0025  0.025    0.0437  1.0     0.0        0.0     0.0  0.0;
    11  110     0.0     0.0167   0.0     1.0     0.0        0.0     0.0  0.0;
    12  120     0.0     0.0167   0.0     1.0     0.0        0.0     0.0  0.0;
    13  101     0.011   0.11     0.1925  1.0     0.0        0.0     0.0  0.0;
    13  101     0.011   0.11     0.1925  1.0     0.0        0.0     0.0  0.0;
    13   14     0.0     0.005    0.00    0.9688  0.0        1.2     0.8  0.00625;
    13  120     0.001   0.01     0.0175  1.0     0.0        0.0     0.0  0.0;
    110 120     0.0025  0.025    0.0437  1.0     0.0        0.0     0.0  0.0];

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
%      18. damping coefficient d_1(pu),
%      19. bus number
%
% note: all the following machines use sub-transient model
% 171004 - Ryan's note: Redistributed the system inertia (via the MVA base)
% Reduced the inertia constant from 6.5 to 5.75

my_H = 5.75;
my_A = 0.05;
my_B = (0.5 - my_A);
my_A = 4*my_A;
my_B = 4*my_B;
my_MVA = 1000;

mac_con = [...
% num bus  mva          xl    ra      xd   xdp   xdpp  tdop  tdopp  xq   xqp   xqpp  tqop  tqopp h     do dl bus  s1    s2
  1    1   my_MVA*my_A  0.20  0.0025  1.8  0.30  0.25  8.00  0.03   1.7  0.55  0.24  0.4   0.05  my_H  0  0   1   0.05  0.3;
  2    2   my_MVA*my_B  0.20  0.0025  1.8  0.30  0.25  8.00  0.03   1.7  0.55  0.24  0.4   0.05  my_H  0  0   2   0.05  0.3;
  3   11   my_MVA*my_A  0.20  0.0025  1.8  0.30  0.25  8.00  0.03   1.7  0.55  0.24  0.4   0.05  my_H  0  0  11   0.05  0.3;
  4   12   my_MVA*my_B  0.20  0.0025  1.8  0.30  0.25  8.00  0.03   1.7  0.55  0.24  0.4   0.05  my_H  0  0  12   0.05  0.3];

clear('my_H','my_A','my_B','my_MVA');

% 171005 - Ryan's note: Changed exciter gain from 200 to 100
exc_con = [...
    0  1  0.02  100.0  0.05  0  0  5.0  -5.0  0  0  0  0  0  0  0  0  0  0  0;
    0  2  0.02  100.0  0.05  0  0  5.0  -5.0  0  0  0  0  0  0  0  0  0  0  0;
    0  3  0.02  100.0  0.05  0  0  5.0  -5.0  0  0  0  0  0  0  0  0  0  0  0;
    0  4  0.02  100.0  0.05  0  0  5.0  -5.0  0  0  0  0  0  0  0  0  0  0  0];

% power system stabilizer model
%       col1    type 1 speed input; 2 power input
%       col2    generator number
%       col3    pssgain*washout time constant
%       col4    washout time constant
%       col5    first lead time constant
%       col6    first lag time constant
%       col7    second lead time constant
%       col8    second lag time constant
%       col9    maximum output limit
%       col10   minimum output limit

pss_con = [...
    1  1  100  10 0.05 0.01 0.05 0.01 0.2 -0.05;
    1  2  100  10 0.05 0.01 0.05 0.01 0.2 -0.05;
    1  3  100  10 0.05 0.01 0.05 0.01 0.2 -0.05;
    1  4  100  10 0.05 0.01 0.05 0.01 0.2 -0.05];

% governor model
% tg_con matrix format
%column        data                     unit
%  1    turbine model number (=1)
%  2    machine number
%  3    speed set point   wf            pu
%  4    steady state gain 1/R           pu
%  5    maximum power order  Tmax       pu on generator base
%  6    servo time constant   Ts        sec
%  7    governor time constant  Tc      sec
%  8    transient gain time constant T3 sec
%  9    HP section time constant   T4   sec
% 10    reheater time constant    T5    sec

% 171005 - Ryan's note: Changed droop constant from 0.04 to 0.05
my_Rg = 0.05;     % droop constant
my_Kg = 1/my_Rg;

tg_con = [...
1  1  1  my_Kg  1.0  0.1  0.5  0.0  1.25  5.0;
1  2  1  my_Kg  1.0  0.1  0.5  0.0  1.25  5.0;
1  3  1  my_Kg  1.0  0.1  0.5  0.0  1.25  5.0;
1  4  1  my_Kg  1.0  0.1  0.5  0.0  1.25  5.0];

% load model (zip parameters)
%       col1    bus number
%       col2    proportion of constant active power load
%       col3    proportion of constant reactive power load
%       col4    proportion of constant active current load
%       col5    proportion of constant reactive current load

%           bus Pcont Qconst P_Icont Q_Iconst
load_con = [  4     0  0        0.5  0.0;  % 50% constant I/Z real power
             14     0  0        0.5  0.0;  % 50% constant I/Z real power
             10     0  0        1.0  1.0;  % constant current (ess)
             20     0  0        1.0  1.0;  % constant current (ess)
            110     0  0        1.0  1.0;  % constant current (ess)
            120     0  0        1.0  1.0]; % constant current (ess)

% ess_con matrix format
% col   data                                        units
%  1    energy storage system/converter number      integer
%  2    bus number                                  integer
%  3    voltage transducer time constant            sec
%  4    use Pade approximation flag (0 = bypass)    binary
%  5    voltage magnitude time delay (Pade)         sec
%  6    power capacity                              MVA
%  7    energy capacity                             MWh
%  8    converter voltage for current limit         pu
%  9    complex power modulation priority mode      integer
%           1 = active; 2 = reactive; 3 = propor.   --
% 10    converter allowable power factor            pf
% 11    initial state of charge [0,1]               fraction of cap.
% 12    minimum state of charge [0,1]               fraction of cap.
% 13    maximum state of charge [0,1]               fraction of cap.
% 14    converter interface time constant           sec
% 15    active current ramp rate limit              pu/sec on ess base
% 16    reactive current ramp rate limit            pu/sec on ess base
% 17    LVPL breakpoint current (1.22)              pu on ess base
% 18    LVPL zero crossing voltage (0.5)            pu
% 19    LVPL breakpoint voltage (0.9)               pu
% 20    charge/discharge only indicator             integer
%           1 = charge only; 2 = discharge only     --
% 21    charge/discharge efficiency [0,1]           pct as fraction

my_paf = 0;      % Pade flag (0=no delay)
my_Tdv = 0.020;  % voltage magnitude time delay (Pade)
Ppercent = 10;

tmp_ec1 = [...
%  no  bus  Tv    paf     Td      P   E   Vr   Pr  pf
    1,  10, 0.02, my_paf, my_Tdv, Ppercent/100*mac_con(1,3), Ppercent/100*mac_con(1,3)*0.5, 0.93, 1, 0.0;
    2,  20, 0.02, my_paf, my_Tdv, Ppercent/100*mac_con(2,3), Ppercent/100*mac_con(2,3)*0.5, 0.93, 1, 0.0;
    3, 110, 0.02, my_paf, my_Tdv, Ppercent/100*mac_con(3,3), Ppercent/100*mac_con(3,3)*0.5, 0.93, 1, 0.0;
    4, 120, 0.02, my_paf, my_Tdv, Ppercent/100*mac_con(4,3), Ppercent/100*mac_con(4,3)*0.5, 0.93, 1, 0.0];

tmp_ec2 = [
%   Ei   Emn  Emx  Tg    rrp  rrq ilvpl1 zerox brkpt cdi eta
    0.5  0.3  0.7  0.02  30   30   1.22  0.40  0.90   0  0.92;
    0.5  0.3  0.7  0.02  30   30   1.22  0.40  0.90   0  0.92;
    0.5  0.3  0.7  0.02  30   30   1.22  0.40  0.90   0  0.92;
    0.5  0.3  0.7  0.02  30   30   1.22  0.40  0.90   0  0.92];

ess_con = [tmp_ec1, tmp_ec2];

% ess_con = [];

clear('my_paf','my_Tdv','tmp_ec1','tmp_ec2');

% Switching file defines the simulation control
% row 1 col1  simulation start time (s) (cols 2 to 6 zeros)
%       col7  initial time step (s)
% row 2 col1  fault application time (s)
%       col2  bus number at which fault is applied
%       col3  bus number defining far end of faulted line
%       col4  zero sequence impedance in pu on system base
%       col5  negative sequence impedance in pu on system base
%       col6  type of fault - 0 three phase
%                           - 1 line to ground
%                           - 2 line-to-line to ground
%                           - 3 line-to-line
%                           - 4 loss of line with no fault
%                           - 5 loss of load at bus
%                           - 6 no action
%                           - 7 three phase fault without loss of line
%                           - 8 three phase fault with nonzero impedance
%       col7  time step for fault period (s)
%       col8  shunt conductance (pu)
% row 3 col1  near end fault clearing time (s) (cols 2 to 6 zeros)
%       col7  time step for second part of fault (s)
% row 4 col1  far end fault clearing time (s) (cols 2 to 6 zeros)
%       col7  time step for fault cleared simulation (s)
% row 5 col1  time to change step length (s)
%       col7  time step (s)
%
% row n col1 finishing time (s)  (n indicates that intermediate rows may be inserted)

my_Ts = 1/(8*60);
sw_con = [...
    0.0           0    0    0    0    0    my_Ts;   % sets intitial time step
    1.0        3101    3    0    0    0    my_Ts;   % no fault
    1.0+2.5/60    0    0    0    0    0    my_Ts;   % clear near end, original value = 14
    1.0+3/60      0    0    0    0    0    my_Ts;   % clear remote end, original value = 15
    3.8           0    0    0    0    0    my_Ts];  % end simulation (the system goes unstable at t>2.7s with no control)

% eof
