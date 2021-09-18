% Two Area Test Case
% with dc link added


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
%               - 2, generator bus (PV bus)
%               - 3, load bus (PQ bus)
% col11 q_gen_max(pu)
% col12 q_gen_min(pu)
% col13 v_rated (kV)
% col14 v_max  pu
% col15 v_min  pu

bus = [ 1  1.03     18.5  7.50   1.61  0.00  0.00  0.00  0.00 1  99.0  -99.0  22.0   1.1  .9;
	    2  1.01     8.80  7.50   5.76  0.00  0.00  0.00  0.00 2  99.0  -99.0  22.0   1.1  .9;
	    3  1.0     -6.1   0.00   0.00  1.00  0.00  0.00  6.00 2  99.0  -99.0  500.0  1.5  .5;
        4  0.95    -10    0.00   0.00  9.76  1.00  0.00  0.00 3  0.0     0.0  115.0  1.05 .95;
%       5  1.0     -10    0.00   0.00  10.9  3.6   0.00  0.00 3  99.0  -99.0  110.0  1.2  .5;  
        5  1.0     -10    0.00   0.00  10.9  3.6   0.00  0.00 3  99.0  -99.0  115.0  1.2  .5;  
	    10 1.0103   12.1  0.00   0.00  0.00  0.00  0.00  0.00 3  0.0     0.0  230.0  1.5  .5;
	    11 1.03    -6.8   7.16   1.49  0.00  0.00  0.00  0.00 2  99.0  -99.0  22.0   1.1  .9;
	    12 1.01    -16.9  7.00   1.39  0.00  0.00  0.00  0.00 2  99.0  -99.0  22.0   1.1  .9;
	    13 1.0     -31.8  0.00   0.00  0.00  0.00  0.00  0.00 2  99.0  -99.0  500.0  1.5  .5;
        14 0.95    -38    0.00   0.00  17.67 1.00  0.00  0.00 3  0.0     0.0  115.0  1.05 .95; 
%	    15 1.0     -14    0.00   0.00 -10.1  6.09  0.00  6.00 3  99.0  -99.0  110.0  1.2  .5; 
	    15 1.0     -14    0.00   0.00 -10.1  6.09  0.00  6.00 3  99.0  -99.0  115.0  1.2  .5; 
        20 0.9876   2.1   0.00   0.00  0.00  0.00  0.00  0.00 3  0.0     0.0  230.0  1.5  .5;
       101 1.05    -19.3  0.00   2.00  0.00  0.00  0.00  0.00 3  99.0  -99.0  500.0  1.5  .5;
       110 1.0125  -13.4  0.00   0.00  0.00  0.00  0.00  0.00 3  0.0     0.0  230.0  1.5  .5;
       120 0.9938  -23.6  0.00   0.00  0.00  0.00  0.00  0.00 3  0.0     0.0  230.0  1.5  .5 ];


% line data format
% line: from bus, to bus, resistance(pu), reactance(pu),
%       line charging(pu), tap ratio, tap phase, tapmax, tapmin, tapsize

line = [...
1    10  0.0     0.0167   0.00    1.0  0. 0.  0.  0.;
2    20  0.0     0.0167   0.00    1.0  0. 0.  0.  0.;
3     4  0.0     0.005    0.00    1.0  0. 1.5 0.5 0.05;
3     5  0.0     0.001    0.00    1.0  0. 5.0 0.5 0.005;
3    20  0.001   0.0100   0.0175  1.0  0. 0.  0.  0.;
3   101  0.011   0.110    0.1925  1.0  0. 0.  0.  0.;
3   101  0.011   0.110    0.1925  1.0  0. 0.  0.  0.;
10   20  0.0025  0.025    0.0437  1.0  0. 0.  0.  0.;
11  110  0.0     0.0167   0.0     1.0  0. 0.  0.  0.;
12  120  0.0     0.0167   0.0     1.0  0. 0.  0.  0.;
13   14  0.0     0.005    0.00    1.0  0. 1.5 0.5 0.05;
13   15  0.0     0.001    0.00    1.0  0. 5.0 0.5 0.005;     
13  101  0.011   0.11     0.1925  1.0  0. 0.  0.  0.;
13  101  0.011   0.11     0.1925  1.0  0. 0.  0.  0.;
13  120  0.001   0.01     0.0175  1.0  0. 0.  0.  0.;
110 120  0.0025  0.025    0.0437  1.0  0. 0.  0.  0.];

% col 1  hvdc converter number
% col 2  LT bus number in load flow data
% col 3  converter type
%        1 - rectifier
%        2 - inverter
% col 4  dc rated voltage (kV)
% col 5  commutating reactance ohms/bridge
% col 6  number of bridges in series
% col 7  rectifier - alpha min
%        inverter  - gamma min
% col 8  rectifier - alpha max
%        inverter  - gamma max
      
dcsp_con = [...
1  5  1  500  3    4   5   30;
2  15 2  500  3    4  18   25];
%1  5  1  500  6    4   5   30;
%2  15 2  500  6    4  18   25];

% col 1   rectifier number
% col 2   inverter number
% col 3   dc line resistance ohms
% col 4   dc line inductance ( milli H)
% col 5   dc line capacitance (micro F)
% col 6   rectifier smoothing inductance (milli H)
% col 7   inverter smoothing inductance (milli H)
% col 8   dc line rating (MW)
% col 9   current margin for inverter current control %


dcl_con = [...
1  2  20  100.0  0  1000.0  1000.0  1000  15];

% col 1 converter number
% col 2 proportional gain
% col 3 integral gain
% col 4 output gain
% col 5 max integral limit
% col 6 min integral limit
% col 7 max output limit
% col 8 min output limit
% col 9  control type 
%       - 1 rectifier current control
%       - 2 rectifier power control
%       - 0 for inverter controls
% note: the order of the converters must be the same as that in dcsp_con
dcc_con = [...
1  1.0  1.0  1.0  1.5 0  90  5   1;
2  1.0  1.0  0.1  15  0  90  15  0];


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
% note: all the following machines use sub-transient model
mac_con = [ ...

1 1 900 0.200  0.00    1.8  0.30  0.25 8.00  0.03...
                       1.7  0.55  0.24 0.4   0.05...
                       6.5  0  0  1;
2 2 900 0.200  0.00    1.8  0.30  0.25 8.00  0.03...
                       1.7  0.55  0.25 0.4   0   ...
                       6.5  0  0  2;
3 11 900 0.200  0.00   1.8  0.30  0.25 8.00  0.03...
                       1.7  0.55  0.24 0.4   0.05...
                       6.5  0  0  11;
4 12 900 0.200  0.00   1.8  0.30  0.25 8.00  0.03...
                       1.7  0.55  0.25 0.4   0.05...
                       6.5  0  0  12];

%exciter models all smpexc - additional zeros provided to match other exciter models


exc_con = [...
0 1 0.01 200.0  0     0     0    5.0  -5.0...
    0    0      0     0     0    0    0      0      0    0   0;
0 2 0.01 200.0  0     0     0    5.0  -5.0...
    0    0      0     0     0    0    0      0      0    0   0;
0 3 0.01 200.0  0     0     0    5.0  -5.0...
    0    0      0     0     0    0    0      0      0    0   0;
0 4 0.01 200.0  0     0     0    5.0  -5.0...
    0    0      0     0     0    0    0      0      0    0   0];


%pss on all exciters
pss_con = [...
1 1 300.0  20.0  0.06  0.04  0.08  0.04  0.2  -0.05;
1 2 300.0  20.0  0.06  0.04  0.08  0.04  0.2  -0.05;
1 3 300.0  20.0  0.06  0.04  0.08  0.04  0.2  -0.05;
1 4 300.0  20.0  0.06  0.04  0.08  0.04  0.2  -0.05];

% governor model
% tg_con matrix format
%column	       data			unit
%  1	turbine model number (=1)	
%  2	machine number	
%  3	speed set point   wf		pu
%  4	steady state gain 1/R		pu
%  5	maximum power order  Tmax	pu on generator base
%  6	servo time constant   Ts	sec
%  7	governor time constant  Tc	sec
%  8	transient gain time constant T3	sec
%  9	HP section time constant   T4	sec
% 10	reheater time constant    T5	sec

tg_con = [...
1  1  1  25.0  1.0  0.1  0.5 0.0 1.25 5.0;
1  2  1  25.0  1.0  0.1  0.5 0.0 1.25 5.0;
1  3  1  25.0  1.0  0.1  0.5 0.0 1.25 5.0;
1  4  1  25.0  1.0  0.1  0.5 0.0 1.25 5.0];



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
0    0    0    0    0    0    0.005;%sets intitial time step
0.1  3    101  0.4    0.4    1    0.001; %line to ground fault bus 3 line 3-101 
0.2 0    0    0    0    0    0.005; %clear fault at bus 3
0.25 0    0    0    0    0    0.005;%clear remote end
4.0  0    0    0    0    0    0]; % end simulation

% non-conforming load declaration of dc LT buses
% col 1           bus number
% col 2           fraction const active power load
% col 3           fraction const reactive power load
% col 4           fraction const active current load
% col 5           fraction const reactive current load
load_con = [...
             5  0  0  0  0;
            15  0  0  0  0];
