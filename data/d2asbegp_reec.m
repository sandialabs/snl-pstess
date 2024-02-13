% Two Area Test Case
% sub transient generators with static exciters, turbine/governors
% 50% constant current active loads
% load modulation
% with power system stabilizers

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
% col11 q_gen_max (pu)
% col12 q_gen_min (pu)
% col13 v_rated (kV)
% col14 v_max (pu)
% col15 v_min (pu)

% 230118 - Ryan's note: Changed the power flow due to new capacity/inertia distribution
bus = [...
    1  1.0100   -7.8436   1.6800    0.2340         0        0  0.00  0.00 2  1.2   -2.0  22.0   1.1   0.9;
    2  1.0100         0  12.3343    1.9995         0        0  0.00  0.00 1  8.75  -2.0  22.0   1.1   0.9;
    3  0.9978  -19.9001        0         0         0        0  0.00  3.33 3  0.0    0.0  230.0  1.5   0.5;
    4  1.0174  -22.5865        0         0    9.7600   1.0000  0.00  0.00 3  0.0    0.0  115.0  1.05  0.95;
   10  1.0065   -9.4251        0         0         0        0  0.00  0.00 3  0.0    0.0  230.0  1.5   0.5;
   11  1.0100  -32.8705   1.6800    0.3071         0        0  0.00  0.00 2  1.2   -2.0  22.0   1.1   0.9;
   12  1.0100  -25.0255  12.3000    2.1789         0        0  0.00  0.00 2  8.75  -2.0  22.0   1.1   0.9;
   13  0.9923  -44.9639        0         0         0        0  0.00  4.85 3  0.0    0.0  230.0  1.5   0.5;
   14  1.0156  -49.8303        0         0   17.6500   1.0000  0.00  0.00 3  0.0    0.0  115.0  1.05  0.95;
   20  0.9980  -11.7916        0         0         0        0  0.00  0.00 3  0.0    0.0  230.0  1.5   0.5;
  101  1.0177  -32.6699        0         0         0        0  0.00  1.27 3  0.0    0.0  230.0  1.5   0.5;
  110  1.0053  -34.4538        0         0         0        0  0.00  0.00 3  0.0    0.0  230.0  1.5   0.5;
  120  0.9950  -36.8200        0         0         0        0  0.00  0.00 3  0.0    0.0  230.0  1.5   0.5;
];

% bus(:,[4,6,7,9]) = 1.15*bus(:,[4,6,7,9]);  % increasing the loading while maintaining the power factor

% Original power flow
% bus = [...
%    1  1.01     18.5     7.2615 1.274  0.00  0.00  0.00  0.00 1  5.0  -2.0  22.0   1.1  .9;
%    2  1.01      7.725   7.00   1.948  0.00  0.00  0.00  0.00 2  5.0  -2.0  22.0   1.1  .9;
%    3  0.9791   -7.441   0.00   0.00   0.00  0.00  0.00  3.00 3  0.0   0.0  230.0  1.5  .5;
%    4  0.998   -10.232   0.00   0.00   9.76  1.00  0.00  0.00 3  0.0   0.0  115.0  1.05 .95;
%   10  0.9962   11.578   0.00   0.00   0.00  0.00  0.00  0.00 3  0.0   0.0  230.0  1.5  .5;
%   11  1.01     -9.0124  7.00   1.143  0.00  0.00  0.00  0.00 2  5.0  -2.0  22.0   1.1  .9;
%   12  1.01    -19.12    7.00   1.784  0.00  0.00  0.00  0.00 2  5.0  -2.0  22.0   1.1  .9;
%   13  0.9863  -34.063   0.00   0.00   0.00  0.00  0.00  5.00 3  0.0   0.0  230.0  1.5  .5;
%   14  1.0062  -39.02    0.00   0.00   17.65 1.00  0.00  0.00 3  0.0   0.0  115.0  1.05 .95;
%   20  0.9847    0.9748  0.00   0.00   0.00  0.00  0.00  0.00 3  0.0   0.0  230.0  1.5  .5;
%   101 1.0003  -21.054   0.00   0.00   0.00  0.00  0.00  1.27 3  0.0   0.0  230.0  1.5  .5;
%   110 0.9980  -15.69    0.00   0.00   0.00  0.00  0.00  0.00 3  0.0   0.0  230.0  1.5  .5;
%   120 0.9877  -25.87    0.00   0.00   0.00  0.00  0.00  0.00 3  0.0   0.0  230.0  1.5  .5;
% ];

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
  1   10  0.0     0.0167  0.00    1.0     0.0  0.0  0.0  0.0;
  2   20  0.0     0.0167  0.00    1.0     0.0  0.0  0.0  0.0;
  3    4  0.0     0.005   0.00    0.975   0.0  1.2  0.8  0.00625;
  3   20  0.001   0.0100  0.0175  1.0     0.0  0.0  0.0  0.0;
  3  101  0.011   0.110   0.1925  1.0     0.0  0.0  0.0  0.0;
  3  101  0.011   0.110   0.1925  1.0     0.0  0.0  0.0  0.0;
 10   20  0.0025  0.025   0.0437  1.0     0.0  0.0  0.0  0.0;
 11  110  0.0     0.0167  0.0     1.0     0.0  0.0  0.0  0.0;
 12  120  0.0     0.0167  0.0     1.0     0.0  0.0  0.0  0.0;
 13  101  0.011   0.11    0.1925  1.0     0.0  0.0  0.0  0.0;
 13  101  0.011   0.11    0.1925  1.0     0.0  0.0  0.0  0.0;
 13   14  0.0     0.005   0.00    0.9688  0.0  1.2  0.8  0.00625;
 13  120  0.001   0.01    0.0175  1.0     0.0  0.0  0.0  0.0;
110  120  0.0025  0.025   0.0437  1.0     0.0  0.0  0.0  0.0];

% line(:,4) = 1.05*line(:,4);  % increase line reactance to make the system less stable

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
%      20. s(1.0)
%      21. s(1.2)
%
% note: all the following machines use sub-transient model

% 230118 - Ryan's note: Redistributed the system inertia (via the MVA base)
%                       Reduced the inertia constant from 6.5 to 5.75
%                       Gen 1 becomes a PV plant

my_H = 5.75;
my_A = 0.05;
my_B = (0.5 - my_A);
my_MVA = 4*1000;
mac_con = [...
% num bus  mva          xl    ra      xd   xdp   xdpp  tdop  tdopp  xq   xqp   xqpp  tqop  tqopp h     do dl bus  s1    s2
% 1    1   my_MVA*my_A  0.20  0.0025  1.8  0.30  0.25  8.00  0.03   1.7  0.55  0.25  0.4   0.05  my_H  0  0   1   0.05  0.3;
  2    2   my_MVA*my_B  0.20  0.0025  1.8  0.30  0.25  8.00  0.03   1.7  0.55  0.25  0.4   0.05  my_H  0  0   2   0.05  0.3;
% 3   11   my_MVA*my_A  0.20  0.0025  1.8  0.30  0.25  8.00  0.03   1.7  0.55  0.25  0.4   0.05  my_H  0  0  11   0.05  0.3;
  4   12   my_MVA*my_B  0.20  0.0025  1.8  0.30  0.25  8.00  0.03   1.7  0.55  0.25  0.4   0.05  my_H  0  0  12   0.05  0.3];

%     mva
% 1   200 (changed to reec/ess)
% 2  1800 (first synchronous machine)
% 3   200
% 4  1800

% % original mac_con data
% mac_con = [...
% % num bus  mva          xl    ra      xd   xdp   xdpp  tdop  tdopp  xq   xqp   xqpp  tqop  tqopp h     do dl bus  s1    s2
%   1    1   my_MVA*my_A  0.20  0.0025  1.8  0.30  0.25  8.00  0.03   1.7  0.55  0.24  0.4   0.05  my_H  0  0   1   0.05  0.3;
%   2    2   my_MVA*my_B  0.20  0.0025  1.8  0.30  0.25  8.00  0.03   1.7  0.55  0.24  0.4   0.05  my_H  0  0   2   0.05  0.3;
%   3   11   my_MVA*my_A  0.20  0.0025  1.8  0.30  0.25  8.00  0.03   1.7  0.55  0.24  0.4   0.05  my_H  0  0  11   0.05  0.3;
%   4   12   my_MVA*my_B  0.20  0.0025  1.8  0.30  0.25  8.00  0.03   1.7  0.55  0.24  0.4   0.05  my_H  0  0  12   0.05  0.3];

% exc_con matrix format
% column  data
%      1  exciter type (=3 for exc_st3, IEEE type ST3)
%      2  machine number
%      3  transducer filter time constant (T_R - sec)
%      4  voltage regulator gain (K_A)
%      5  voltage regulator time constant (T_A - sec)
%      6  transient gain reduction time constnat (T_B - sec) -- denominator
%      7  transient gain reduction time constnat (T_C - sec) -- numerator
%      8  max voltage regulator output (V_Rmax - pu)
%      9  min voltage regulator output (V_Rmin - pu)
%     10  max internal signal (VImax - pu)
%     11  min internal signal (VImin - pu)
%     12  first state regulator gain (KJ)
%     13  potential circuit gain coef (KP)
%     14  potential circuit phase angle (qP - degrees)
%     15  current circuit gain coef (KI)
%     16  potential source reactance (XL - pu)
%     17  rectifier loading factor (KC)
%     18  max field voltage (Efdmax - pu)
%     19  inner loop feedback constant (KG)
%     20  max innerloop voltage feedback (VGmax - pu)

% 171005 - Ryan's note: Changed exciter gain from 200 to 100
exc_con = [...
%  0  1  0.02  100.0  0.05  0  0  5.0 -5.0  0  0  0  0  0  0  0  0  0  0  0;
   0  2  0.02  100.0  0.05  0  0  5.0 -5.0  0  0  0  0  0  0  0  0  0  0  0;
%  0  3  0.02  100.0  0.05  0  0  5.0 -5.0  0  0  0  0  0  0  0  0  0  0  0;
   0  4  0.02  100.0  0.05  0  0  5.0 -5.0  0  0  0  0  0  0  0  0  0  0  0];

% simple exciter model, type 0; there are three exciter models
% exc_con = [...
% 0 1 0.01 200.0  0.05     0     0    5.0  -5.0...
%     0    0      0     0     0    0    0      0      0    0   0;
% 0 2 0.01 200.0  0.05     0     0    5.0  -5.0...
%     0    0      0     0     0    0    0      0      0    0   0;
% 0 3 0.01 200.0  0.05     0     0    5.0  -5.0...
%     0    0      0     0     0    0    0      0      0    0   0;
% 0 4 0.01 200.0  0.05     0     0    5.0  -5.0...
%     0    0      0     0     0    0    0      0      0    0   0];

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

% pss model parameters
pss_con = [...
% 1  1  100  10 0.05 0.01 0.05 0.01 0.2 -0.05;
  1  2  100  10 0.05 0.01 0.05 0.01 0.2 -0.05;
% 1  3  100  10 0.05 0.01 0.05 0.01 0.2 -0.05;
  1  4  100  10 0.05 0.01 0.05 0.01 0.2 -0.05];

% pss_con = [];

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
my_Rg = 0.05;     % droop constant (K=20)
my_Kg = 1/my_Rg;

tg_con = [...
% 1  1  1  my_Kg  1.0  0.5  0.01  0.0  1.25  5.0;
  1  2  1  my_Kg  1.0  0.5  0.01  0.0  1.25  5.0;
% 1  3  1  my_Kg  1.0  0.5  0.01  0.0  1.25  5.0;
  1  4  1  my_Kg  1.0  0.5  0.01  0.0  1.25  5.0];

clear('my_Rg','my_Kg');

% original governor models
% tg_con = [...
% 1  1  1  25.0  1.0  0.1  0.5 0.0 1.25 5.0;
% 1  2  1  25.0  1.0  0.1  0.5 0.0 1.25 5.0;
% 1  3  1  25.0  1.0  0.1  0.5 0.0 1.25 5.0;
% 1  4  1  25.0  1.0  0.1  0.5 0.0 1.25 5.0;
% ];

% Induction motor data
% 1. Motor Number
% 2. Bus Number
% 3. Motor MVA Base
% 4. rs pu
% 5. xs pu - stator leakage reactance
% 6. Xm pu - magnetizing reactance
% 7. rr pu
% 8. xr pu - rotor leakage reactance
% 9. H  s  - motor plus load inertia constant
% 10. rr1 pu - second cage resistance
% 11. xr1 pu - intercage reactance
% 12. dbf    - deepbar factor
% 13. isat pu - saturation current
% 15. fraction of bus power drawn by motor ( if zero motor statrts at t=0)

ind_con = [];
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
mld_con = [];

%       col1    bus number
%       col2    proportion of constant active power load
%       col3    proportion of constant reactive power load
%       col4    proportion of constant active current load
%       col5    proportion of constant reactive current load
load_con = [  1  0  0  1.0  1.0;   % constant current (reec/ess)
             11  0  0  1.0  1.0;   % constant current (reec/ess)
              4  0  0  1.0  0.0;   % 50% constant I/Z real power
             14  0  0  1.0  0.0];  % 50% constant I/Z real power

disp('100% constant current real power load')

%disp('load modulation')
%active and reactive load modulation enabled
%       col1    load number
%       col2    bus number
%       col3    MVA rating
%       col4    maximum output limit pu
%       col4    minimum output limit pu
%       col6    Time constant

lmod_con = [...
  1    4  100  1 -1  1  0.05;
  2   14  100  1 -1  1  0.05];

% lmod_con = [];

rlmod_con = [...
  1    4  100  1 -1  1  0.05;
  2   14  100  1 -1  1  0.05];

% rlmod_con = [];

% reec_con matrix format
%
% col   data                                                units
% 1     MVA base (mvab), inherited from ess                 MVA
% 2     bus number                                          integer
% 3     voltage dip logic upper limit (vdiph)               pu
% 4     voltage dip logic lower limit (vdipl)               pu
% 5     transducer time constant (trv)                      sec
% 6     lower deadband in voltage error (dbd1)              pu
% 7     upper deadband in voltage error (dbd2)              pu
% 8     reactive current injection gain (kqv)               pu
% 9     maximum limit of reactive current injection (iqh1)  pu
% 10    minimum limit of reactive current injection (iql1)  pu
% 11    value at which iqcmd is held for "thldq" seconds    pu
%           following a voltage dip (iqfrz)
% 12    post voltage dip hold time for iqcmd (thldq)        sec
% 13    post voltage dip hold time for ipcmd (thldp)        sec
% 14    active power transducer time constant (tp)          sec
% 15    reactive power transducer time constant (tq)        sec
% 16    reactive power voltage maximum limit (qvmax)        pu
% 17    reactive power voltage minimum limit (qvmin)        pu
% 18    voltage control maximum limit (vmax)                pu
% 19    voltage control minimum limit (vmin)                pu
% 20    anti-windup flag for Q, V control (=1 anti-windup)  binary
% 21    proportional gain (kqp)                             pu
% 22    integral gain (kqi)                                 pu
% 23    proportional gain (kvp)                             pu
% 24    integral gain (kvi)                                 pu
% 25    time constant (tiq)                                 sec
% 26    up ramp rate on power reference (dpmax)             pu/sec
% 27    down ramp rate on power reference (dpmin)           pu/sec
% 28    maximum power reference (pmax)                      pu
% 29    minimum power reference (pmin)                      pu
% 30    time constant (tpord)                               pu
% 31    voltage compensation flag (vcmpflag)                integer
%           1 = current compensation; 0 = reactive droop
% 32    power factor flag (pfflag)                          integer
%           1 = power factor control; 0 = Q control
% 33    voltage control flag (vflag)                        integer
%           1 = Q control; 0 = Voltage control
% 34    reactive power control flag (qflag)                 integer
%           1 = Voltage/Q control;
%           0 = Constant power factor or Q control
% 35    power reference flag (pflag) [UNUSED]               integer
%           1 = Pref*speed; 0 = Pref
% 36    current compensation resistance on gen. base (rc)   pu
% 37    current compensation reactance on gen. base (xc)    pu
% 38    filter time constant for voltage measurement (tr1)  sec
% 39    reactive current compensation gain (kc)             pu
% 40    maximum converter blocking voltage (vblkh)          pu
% 41    minimum converter blocking voltage (vblkl)          pu
% 42    voltage blocking delay (tblkdelay)                  sec
% 43    time constant for active current command filter     sec
% 44    time constant for reactive current command filter   sec

%          mva bus vdiph vdipl trv   dbd1 dbd2 kqv iqh1  iql1
tmp_rc1 = [ 0   1   99    -99  0.02  0.0  0.0  0.1 1.05 -1.05;
            0  11   99    -99  0.02  0.0  0.0  0.1 1.05 -1.05];

%         iqfrz thldq thdlp tp   tq   qvmax  qvmin vmax vmin qvwindup
tmp_rc2 = [0.15  0.0  0.0   0.05 0.05 0.436 -0.436 1.10 0.90  1;
           0.15  0.0  0.0   0.05 0.05 0.436 -0.436 1.10 0.90  1];

%         kqp   kqi  kvp  kvi   tiq   dpmax dpmin pmax  pmin  tpord
tmp_rc3 = [0.0  0.5  0.0  40.0  0.02  99    -99   1.12  0.04  0.02;
           0.0  0.5  0.0  40.0  0.02  99    -99   1.12  0.04  0.02];

%         ------flags-------
%         vcmp pf  v   q   p  rc   xc    tr1   kc    vblkh vblkl tblkdelay tipc  tiqc
tmp_rc4 = [1   0   1   1   0  0.0  0.01  0.03  0.01  1.2   0.1   0.05      0.0   0.0;
           1   0   1   1   0  0.0  0.01  0.03  0.01  1.2   0.1   0.05      0.0   0.0];

reec_con = [tmp_rc1, tmp_rc2, tmp_rc3, tmp_rc4];

% reec_con = [];  % uncommenting this line turns off reec controllers

clear('tmp_rc1','tmp_rc2','tmp_rc3','tmp_rc4');

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
% 22    scaling on limit ipmin (ke) [0,1]           pu

% note: since E=0, this models a PV plant and
%       ess_con(,11:13) are ignored

my_paf = 0;      % bypass Pade approximation
my_Vr = 1/20.0;  % converter voltage for current limit

tmp_ec1 = [...
%  no  bus  Tv    paf Td   P            E  Vr    Pr  pf
    1    1  0.02   0  0.0  my_MVA*my_A  0  my_Vr  2  0.0;
    2   11  0.02   0  0.0  my_MVA*my_A  0  my_Vr  2  0.0];

tmp_ec2 = [...
%   Ei   Emn  Emx  Tg    rrp  rrq ilvpl1 zerox brkpt cdi  eta  ke
    0.5  0.0  1.0  0.02  300  300  -1.22  0.40  0.90   0  0.92 0.0;
    0.5  0.0  1.0  0.02  300  300  -1.22  0.40  0.90   0  0.92 0.0];

ess_con = [tmp_ec1, tmp_ec2];

% ess_con = [];  % uncomment this line to turn off ess models
clear('my_H','my_A','my_B','my_MVA');
clear('my_paf','my_Vr','tmp_ec1','tmp_ec2');

%----------------------------------------------------------------------------%
% Monitored lines
%
% When conducting an eigenanalysis, initializing lmon_con causes the
% following variables to be created: c_pf1, c_pf2, c_qf1, c_qf2, c_ilif,
% c_ilit, c_ilmf, c_ilmt, c_ilrf, c_ilrt. The ith row of c_pf1 and c_qf1
% correspond to line(lmon_con(i),:) of the line matrix.

lmon_con = 1:size(line,1);  % all lines

% Switching file defines the simulation control
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
%                            - 6 no action
%                            - 7 three-phase fault without loss of line
%                            - 8 three-phase resistive shunt insertion
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

my_case = 'do_nothing';

switch my_case
    case 'do_nothing'
        sw_con = [...
            0.0        0    0    0    0    0    my_Ts;   % sets initial time step
            1.0       11  110    0    0    6    my_Ts;   % do nothing
            1.5        0    0    0    0    0    my_Ts;   % clear near end
            1.5+1/60   0    0    0    0    0    my_Ts;   % clear remote end
           20.0        0    0    0    0    0    my_Ts];  % end simulation
    case 'gen3_drop'
        sw_con = [...
            0.0        0    0    0    0    0    my_Ts;   % sets initial time step
            5.0       11  110    0    0    4    my_Ts;   % gen trip
            5.5        0    0    0    0    0    my_Ts;   % clear near end
            5.5+1/60   0    0    0    0    0    my_Ts;   % clear remote end
           40.0        0    0    0    0    0    my_Ts];  % end simulation
    case 'brake_insertion'
        sw_con = [...
            0.0        0    0    0    0    0    my_Ts;   % sets initial time step
            5.0       20    0 1/0.50  0    8    my_Ts;   % brake insertion
            5.5        0    0    0    0    0    my_Ts;   % clear near end
            5.5+1/60   0    0    0    0    0    my_Ts;   % clear remote end
           40.0        0    0    0    0    0    my_Ts];  % end simulation
    case '3phase_fault'
        % 3-phase fault at bus 3 (far end 101)
        sw_con = [...
            0.0        0    0    0    0    0    my_Ts;   % sets initial time step
            5.0        3  101  1e-1   0    0    my_Ts;   % 3-ph fault
            5.0+6/60   0    0    0    0    0    my_Ts;   % clear near end
            5.0+7/60   0    0    0    0    0    my_Ts;   % clear remote end
           40.0        0    0    0    0    0    my_Ts];  % end simulation
    otherwise
        error('d2asbegp_reec: contingency not defined.\n')
end

% ibus_con = [0 1 1];  % sets infinite buses

clear('my_Ts','my_case');

% eof
