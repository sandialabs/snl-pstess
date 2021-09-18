% IEEE 10-generator, 39-bus system
% Also known as the 10-machine New England Power System
%
% Based on dynamic data from, ``IEEE PES Task Force on Benchmark Systems for Stability Controls,'' Hiskens.
% Case first presented in, ``A Practical Method for the Direct Analysis of Transient,'' Athay et al.
% Additional information provided in, ``Energy function analysis for power system stability,'' M.A. Pai.
%
% Prepared by: Ryan Elliott, Sandia National Laboratories, rtellio@sandia.gov
% Date: July 10, 2020

gov_flag = true;  % (true=turbine governor models, false=no governors)

%-------------------------------------------------------------------------------------------------------------%
% power flow data

% bus format
%    1  bus   : number
%    2  volt  : voltage magnitude (pu)
%    3  ang   : voltage angle (deg)
%    4  pgen  : active power generation (pu)
%    5  qgen  : reactive power generation (pu)
%    6  pload : active power load (pu)
%    7  qload : reactive power load (pu)
%    8  gsh   : shunt conductance (pu)
%    9  bsh   : shunt susceptance (pu)
%   10  type  : bus type (3=load, 2=gen, 1=swing)
%   11  qgmax : maximum reactive power generation (pu)
%   12  qgmin : minimum reactive power generation (pu)
%   13  vrate : rated voltage (kV)
%   14  vmax  : maximum allowable voltage at bus (pu)
%   15  vmin  : minimum allowable voltage at bus (pu)
%
% note: in power flow load buses (type 3) are treated as p-q buses,
%       gen buses (type 2) are treated as p-v buses, and the swing
%       bus is treated as p-delta

disp('10-machine New England power system data')
% bus  volt  ang  pgen  qgen  pload  qload  gshunt  bshunt  bus_type  qgmax  qgmin  vrated  vmax  vmin
bus = [
   1  1.0474  -8.44   0.0000  0.0000   0.0000   0.0000  0.0  0.0  3  0.0  0.0  345.0  1.5  0.5;  % new scotland
   2  1.0487  -5.75   0.0000  0.0000   0.0000   0.0000  0.0  0.0  3  0.0  0.0  345.0  1.5  0.5;  % northfield
   3  1.0302  -8.60   0.0000  0.0000   3.2200   0.0240  0.0  0.0  3  0.0  0.0  345.0  1.5  0.5;  % ludlow
   4  1.0039  -9.61   0.0000  0.0000   5.0000   1.8400  0.0  0.0  3  0.0  0.0  345.0  1.5  0.5;  % manchester
   5  1.0053  -8.61   0.0000  0.0000   0.0000   0.0000  0.0  0.0  3  0.0  0.0  345.0  1.5  0.5;  % scovill rock
   6  1.0077  -7.95   0.0000  0.0000   0.0000   0.0000  0.0  0.0  3  0.0  0.0  345.0  1.5  0.5;  % haddam neck
   7  0.9970 -10.12   0.0000  0.0000   2.3380   0.8400  0.0  0.0  3  0.0  0.0  345.0  1.5  0.5;  % black pond
   8  0.9960 -10.62   0.0000  0.0000   5.2200   1.7600  0.0  0.0  3  0.0  0.0  345.0  1.5  0.5;  % southington
   9  1.0282 -10.32   0.0000  0.0000   0.0000   0.0000  0.0  0.0  3  0.0  0.0  345.0  1.5  0.5;  % pleasant valley
  10  1.0172  -5.43   0.0000  0.0000   0.0000   0.0000  0.0  0.0  3  0.0  0.0  345.0  1.5  0.5;  % millstone
  11  1.0127  -6.28   0.0000  0.0000   0.0000   0.0000  0.0  0.0  3  0.0  0.0  345.0  1.5  0.5;  % montville
  12  1.0002  -6.24   0.0000  0.0000   0.0750   0.8800  0.0  0.0  3  0.0  0.0  115.0  1.5  0.5;  % montville
  13  1.0143  -6.10   0.0000  0.0000   0.0000   0.0000  0.0  0.0  3  0.0  0.0  345.0  1.5  0.5;  % montville
  14  1.0117  -7.66   0.0000  0.0000   0.0000   0.0000  0.0  0.0  3  0.0  0.0  345.0  1.5  0.5;  % card street
  15  1.0154  -7.74   0.0000  0.0000   3.2000   1.5300  0.0  0.0  3  0.0  0.0  345.0  1.5  0.5;  % sherman road
  16  1.0318  -6.19   0.0000  0.0000   3.2900   0.3230  0.0  0.0  3  0.0  0.0  345.0  1.5  0.5;  % medway
  17  1.0336  -7.30   0.0000  0.0000   0.0000   0.0000  0.0  0.0  3  0.0  0.0  345.0  1.5  0.5;  % millbury
  18  1.0309  -8.22   0.0000  0.0000   1.5800   0.3000  0.0  0.0  3  0.0  0.0  345.0  1.5  0.5;  % carpenter hill
  19  1.0499  -1.02   0.0000  0.0000   0.0000   0.0000  0.0  0.0  3  0.0  0.0  345.0  1.5  0.5;  % brayton point
  20  0.9912  -2.01   0.0000  0.0000   6.2800   1.0300  0.0  0.0  3  0.0  0.0  115.0  1.5  0.5;  % brayton point
  21  1.0318  -3.78   0.0000  0.0000   2.7400   1.1500  0.0  0.0  3  0.0  0.0  345.0  1.5  0.5;  % bridgewater
  22  1.0498   0.67   0.0000  0.0000   0.0000   0.0000  0.0  0.0  3  0.0  0.0  345.0  1.5  0.5;  % pilgrim
  23  1.0448   0.47   0.0000  0.0000   2.4750   0.8460  0.0  0.0  3  0.0  0.0  345.0  1.5  0.5;  % canal
  24  1.0373  -6.07   0.0000  0.0000   3.0860  -0.9220  0.0  0.0  3  0.0  0.0  345.0  1.5  0.5;  % walpole
  25  1.0576  -4.36   0.0000  0.0000   2.2400   0.4720  0.0  0.0  3  0.0  0.0  345.0  1.5  0.5;  % vermont yankee
  26  1.0521  -5.53   0.0000  0.0000   1.3900   0.1700  0.0  0.0  3  0.0  0.0  345.0  1.5  0.5;  % scobie pond
  27  1.0377  -7.50   0.0000  0.0000   2.8100   0.7550  0.0  0.0  3  0.0  0.0  345.0  1.5  0.5;  % sandy pond
  28  1.0501  -2.01   0.0000  0.0000   2.0600   0.2760  0.0  0.0  3  0.0  0.0  345.0  1.5  0.5;  % surowiec (pownal)
  29  1.0499   0.74   0.0000  1.0000   2.8350   0.2690  0.0  0.0  2  0.0  0.0  345.0  1.5  0.5;  % maine yankee
  30  1.0475  -3.33   2.5000  1.4616   0.0000   0.0000  0.0  0.0  2  3.0 -2.0   13.8  1.5  0.5;  % gen 10 // northfield
  31  0.9820   0.00   5.2081  1.9825   0.0920   0.0460  0.0  0.0  1  2.5 -1.5   13.8  1.5  0.5;  % gen 2 // connecticut yankee (swing bus)
  32  0.9831   2.57   6.5000  2.0514   0.0000   0.0000  0.0  0.0  2  3.0 -2.0   13.8  1.5  0.5;  % gen 3 // millstone
  33  0.9972   4.19   6.3200  1.0991   0.0000   0.0000  0.0  0.0  2  3.0 -2.0   13.8  1.5  0.5;  % gen 4 // brayton point
  34  1.0123   3.17   5.0800  1.6576   0.0000   0.0000  0.0  0.0  2  2.5 -1.5   13.8  1.5  0.5;  % gen 5 // brayton point
  35  1.0493   5.63   6.5000  2.1241   0.0000   0.0000  0.0  0.0  2  3.0 -2.0   13.8  1.5  0.5;  % gen 6 // pilgrim
  36  1.0635   8.32   5.6000  1.0117   0.0000   0.0000  0.0  0.0  2  2.5 -1.5   13.8  1.5  0.5;  % gen 7 // canal
  37  1.0278   2.42   5.4000  0.0044   0.0000   0.0000  0.0  0.0  2  2.5 -1.5   13.8  1.5  0.5;  % gen 8 // vermont yankee
  38  1.0265   7.81   8.3000  0.2284   0.0000   0.0000  0.0  0.0  2  4.0 -2.5   13.8  1.5  0.5;  % gen 9 // maine yankee
  39  1.0300 -10.05  10.0000  0.8828  11.0400   2.5000  0.0  0.0  2  4.5 -3.0   13.8  1.5  0.5]; % gen 1 // external (aggregate)

% line format
%    1  f_bus    : from bus number
%    2  t_bus    : to bus number
%    3  r        : line resistance (pu)
%    4  x        : line reactance (pu)
%    5  b        : line-charging susceptance (pu)
%    6  tap_mag  : tap ratio
%    7  tap_ang  : tap phase (for phase shifters)
%    8  tap_max  : tap maximum (pu)
%    9  tap_min  : tap minimum (pu)
%   10  tap_step : tap step size (pu)

% from_bus  to_bus  r  x  b  tap_mag  tap_ang  tap_max  tap_min  tap_step
line = [
   1   2  0.00350  0.04110  0.69870  0.0000  0.00  0.00  0.00  0.00;
   1  39  0.00100  0.02500  0.75000  0.0000  0.00  0.00  0.00  0.00;
   2   3  0.00130  0.01510  0.25720  0.0000  0.00  0.00  0.00  0.00;
   2  25  0.00700  0.00860  0.14600  0.0000  0.00  0.00  0.00  0.00;
   3   4  0.00130  0.02130  0.22140  0.0000  0.00  0.00  0.00  0.00;
   3  18  0.00110  0.01330  0.21380  0.0000  0.00  0.00  0.00  0.00;
   4   5  0.00080  0.01280  0.13420  0.0000  0.00  0.00  0.00  0.00;
   4  14  0.00080  0.01290  0.13820  0.0000  0.00  0.00  0.00  0.00;
   5   6  0.00020  0.00260  0.04340  0.0000  0.00  0.00  0.00  0.00;
   5   8  0.00080  0.01120  0.14760  0.0000  0.00  0.00  0.00  0.00;
   6   7  0.00060  0.00920  0.11300  0.0000  0.00  0.00  0.00  0.00;
   6  11  0.00070  0.00820  0.13890  0.0000  0.00  0.00  0.00  0.00;
   7   8  0.00040  0.00460  0.07800  0.0000  0.00  0.00  0.00  0.00;
   8   9  0.00230  0.03630  0.38040  0.0000  0.00  0.00  0.00  0.00;
   9  39  0.00100  0.02500  1.20000  0.0000  0.00  0.00  0.00  0.00;
  10  11  0.00040  0.00430  0.07290  0.0000  0.00  0.00  0.00  0.00;
  10  13  0.00040  0.00430  0.07290  0.0000  0.00  0.00  0.00  0.00;
  13  14  0.00090  0.01010  0.17230  0.0000  0.00  0.00  0.00  0.00;
  14  15  0.00180  0.02170  0.36600  0.0000  0.00  0.00  0.00  0.00;
  15  16  0.00090  0.00940  0.17100  0.0000  0.00  0.00  0.00  0.00;
  16  17  0.00070  0.00890  0.13420  0.0000  0.00  0.00  0.00  0.00;
  16  19  0.00160  0.01950  0.30400  0.0000  0.00  0.00  0.00  0.00;
  16  21  0.00080  0.01350  0.25480  0.0000  0.00  0.00  0.00  0.00;
  16  24  0.00030  0.00590  0.06800  0.0000  0.00  0.00  0.00  0.00;
  17  18  0.00070  0.00820  0.13190  0.0000  0.00  0.00  0.00  0.00;
  17  27  0.00130  0.01730  0.32160  0.0000  0.00  0.00  0.00  0.00;
  21  22  0.00080  0.01400  0.25650  0.0000  0.00  0.00  0.00  0.00;
  22  23  0.00060  0.00960  0.18460  0.0000  0.00  0.00  0.00  0.00;
  23  24  0.00220  0.03500  0.36100  0.0000  0.00  0.00  0.00  0.00;
  25  26  0.00320  0.03230  0.51300  0.0000  0.00  0.00  0.00  0.00;
  26  27  0.00140  0.01470  0.23960  0.0000  0.00  0.00  0.00  0.00;
  26  28  0.00430  0.04740  0.78020  0.0000  0.00  0.00  0.00  0.00;
  26  29  0.00570  0.06250  1.02900  0.0000  0.00  0.00  0.00  0.00;
  28  29  0.00140  0.01510  0.24900  0.0000  0.00  0.00  0.00  0.00;
  12  11  0.00160  0.04350  0.00000  1.0060  0.00  1.20  0.80  0.05;  % xfrmr // millstone
  12  13  0.00160  0.04350  0.00000  1.0060  0.00  1.20  0.80  0.05;  % xfrmr // millstone
  19  20  0.00070  0.01380  0.00000  1.0600  0.00  1.20  0.80  0.05;  % xfrmr // brayton point
   2  30  0.00000  0.01810  0.00000  1.0250  0.00  1.20  0.80  0.05;  % step-up xfrmr // northfield
   6  31  0.00000  0.02500  0.00000  1.0700  0.00  1.20  0.80  0.05;  % step-up xfrmr // connecticut yankee (swing bus)
  10  32  0.00000  0.02000  0.00000  1.0700  0.00  1.20  0.80  0.05;  % step-up xfrmr // millstone
  19  33  0.00070  0.01420  0.00000  1.0700  0.00  1.20  0.80  0.05;  % step-up xfrmr // brayton point
  20  34  0.00090  0.01800  0.00000  1.0090  0.00  1.20  0.80  0.05;  % step-up xfrmr // brayton point
  22  35  0.00000  0.01430  0.00000  1.0250  0.00  1.20  0.80  0.05;  % step-up xfrmr // pilgrim
  23  36  0.00050  0.02720  0.00000  1.0000  0.00  1.20  0.80  0.05;  % step-up xfrmr // canal
  25  37  0.00060  0.02320  0.00000  1.0250  0.00  1.20  0.80  0.05;  % step-up xfrmr // vermont yankee
  29  38  0.00080  0.01560  0.00000  1.0250  0.00  1.20  0.80  0.05]; % step-up xfrmr // maine yankee

% note: unit 1 (bus 39) has no step-up transformer and is modeled as
%       directly connected at high voltage because it represents a large
%       aggregation (primarily generation from new york state).

%-------------------------------------------------------------------------------------------------------------%
% dynamic data

% mac_con format
%    1  gen    : machine number (may be different from bus number)
%    2  bus    : bus number
%    3  mva    : machine base mva (mva)
%    4  xl     : leakage reactance x_l (pu)
%    5  ra     : resistance r_a (pu)
%    6  xd     : d-axis sychronous reactance (pu)
%    7  xdp    : d-axis transient reactance (pu)
%    8  xdpp   : d-axis subtransient reactance (pu)
%    9  tdop   : d-axis open-circuit time constant (s)
%   10  tdopp  : d-axis open-circuit subtransient time constant (s)
%   11  xq     : q-axis sychronous reactance (pu)
%   12  xqp    : q-axis transient reactance (pu)
%   13  xqpp   : q-axis subtransient reactance (pu)
%   14  tqop   : q-axis open-circuit time constant (s)
%   15  tqopp  : q-axis open-circuit subtransient time constant (s)
%   16  h      : inertia constant (s)
%   17  do     : damping coefficient (pu)
%   18  dl     : damping coefficient relative to pmech (pu)
%   19  bus    : bus number
%   20  s(1.0) : saturation factor at 1.0 pu flux
%   21  s(1.2) : saturation factor at 1.2 pu flux
%
% note 1: all per-unit parameters above are spec'd on machine base unless otherwise noted
% note 2: pst requires that xqpp = xdpp
% note 3: if either tqop = 0 or tqopp = 0, the time constant is bypassed

% gen bus mva  xl     ra   xd      xdp     xdpp    tdop  tdopp  xq      xqp     xqpp    tqop  tqopp  h      do   dl  bus  s(1.0) s(1.2)
mac_con = [
  1  39 5000  0.1500  0.0  1.0000  0.3000  0.1500   7.00  0.03  0.9500  0.4000  0.1500  0.70  0.05  10.000  0.0  0.0  39  0.050  0.30;  % external
  2  31  560  0.1960  0.0  1.6520  0.3903  0.1952   6.56  0.03  1.5792  0.9520  0.1952  1.50  0.05   5.411  0.0  0.0  31  0.050  0.30;  % connecticut yankee
  3  32  870  0.2645  0.0  2.1707  0.4620  0.2310   5.70  0.03  2.0619  0.7621  0.2310  1.50  0.05   4.115  0.0  0.0  32  0.050  0.30;  % millstone
  4  33  950  0.2803  0.0  2.4890  0.4142  0.2071   5.69  0.03  2.4510  1.5770  0.2071  1.50  0.05   3.011  0.0  0.0  33  0.050  0.30;  % brayton point
  5  34  540  0.2916  0.0  3.6180  0.7128  0.3564   5.40  0.03  3.3480  0.8964  0.3564  0.44  0.05   4.815  0.0  0.0  34  0.050  0.30;  % brayton point
  6  35  690  0.1546  0.0  1.7526  0.3450  0.1725   7.30  0.03  1.6629  0.5617  0.1725  0.40  0.05   5.043  0.0  0.0  35  0.050  0.30;  % pilgrim
  7  36  800  0.2576  0.0  2.3600  0.3920  0.1960   5.66  0.03  2.3360  1.4880  0.1960  1.50  0.05   3.300  0.0  0.0  36  0.050  0.30;  % canal
  8  37  620  0.1736  0.0  1.7980  0.3534  0.1767   6.70  0.03  1.7360  0.5648  0.1767  0.41  0.05   3.919  0.0  0.0  37  0.050  0.30;  % vermont yankee
  9  38  880  0.2622  0.0  1.8533  0.5016  0.2508   4.79  0.03  1.8040  0.5166  0.2508  1.96  0.05   3.920  0.0  0.0  38  0.050  0.30;  % maine yankee
 10  30 1200  0.1500  0.0  1.2000  0.3720  0.1860  10.20  0.03  0.8280  0.0960  0.1860  0.00  0.05   3.500  0.0  0.0  30  0.050  0.30]; % northfield

% tg_con format
%    1  type  : turbine governor type (1=thermal model, 2=hydro model)
%    2  num   : machine number
%    3  wf    : speed set point (pu)
%    4  1/R   : steady-state gain (reciprocal of droop constant, pu)
%    5  t_max : maximum power order (pu on machine base)
%    6  ts    : servo time constant (s)
%    7  tc    : governor time constant (s)
%    8  t3    : transient gain time constant (s)
%    9  t4    : high-pressure section time constant (s)
%   10  t5    : reheater time constant (s)

% type  num  wf    1/R  t_max  ts    tc    t3    t4   t5
tg_con = [
   1   1  1.0  0.50*20.0  1.0  0.10  10.0  3.0   0.0  0.01;  % mix, 50 pct baseloaded
   1   2  1.0  0.00*20.0  1.0  0.10  10.0  3.0   0.0  0.01;  % steam (nuclear), no gov
   1   3  1.0  0.00*20.0  1.0  0.10  10.0  3.0   0.0  0.01;  % steam (nuclear), no gov
   1   4  1.0       20.0  1.0  0.10  10.0  3.0   0.0  0.01;  % steam (coal)
   1   5  1.0       20.0  1.0  0.10  10.0  3.0   0.0  0.01;  % steam (coal)
   1   6  1.0  0.00*20.0  1.0  0.10  10.0  3.0   0.0  0.01;  % steam (nuclear), no gov
   1   7  1.0       20.0  1.0  0.50  10.0  4.0   0.0  1.00;  % gas turbine
   1   8  1.0  0.00*20.0  1.0  0.10  10.0  3.0   0.0  0.01;  % steam (nuclear), no gov
   1   9  1.0  0.00*20.0  1.0  0.10  10.0  3.0   0.0  0.01;  % steam (nuclear), no gov
   1  10  1.0       20.0  1.0  0.50  15.0  1.0  -1.0  0.50;  % slow hydro
];

% check for the inclusion of turbine governor data
if ~gov_flag
    tg_con(:,4) = 0.0;  % baseload all units
end

% exc_con format
%    1  type   : exciter type (0=smpexc, 1=dc1, 2=dc2, 3=st3, 4=smppi)
%    2  num    : machine number
%    3  tr     : input filter time constant (s)
%    4  ka     : voltage regulator gain (pu)
%    5  ta     : voltage regulator time constant (s)
%    6  tb     : voltage regulator time constant (s)
%    7  tc     : voltage regulator time constant (s)
%    8  vr_max : maximum voltage regulator output (pu)
%    9  vr_min : minimum voltage regulator output (pu)
%
% note: the meaning of the remaining columns varies depending on the
%       specific exciter model being employed (see the pst documentation
%       for details).

% type  num  tr  ka   ta     tb    tc  vr_max  vr_min  vi_max  vi_min  kj  kp   qp   ki   xl   kc  efd_max kg  vg_max
exc_con = [
   3   1  0.0  200.0  0.02  10.0  1.0  5.000  -5.000  0.100  -0.100  1.0  1.0  0.0  0.0  0.0  0.0  5.000  0.0  5.000;  % external
   3   2  0.0  200.0  0.02  10.0  1.0  5.000  -5.000  0.100  -0.100  1.0  1.0  0.0  0.0  0.0  0.0  5.000  0.0  5.000;  % connecticut yankee
   3   3  0.0  200.0  0.02  10.0  1.0  5.000  -5.000  0.100  -0.100  1.0  1.0  0.0  0.0  0.0  0.0  5.000  0.0  5.000;  % millstone
   3   4  0.0  200.0  0.02  10.0  1.0  5.000  -5.000  0.100  -0.100  1.0  1.0  0.0  0.0  0.0  0.0  5.000  0.0  5.000;  % brayton point
   3   5  0.0  200.0  0.02  10.0  1.0  5.000  -5.000  0.100  -0.100  1.0  1.0  0.0  0.0  0.0  0.0  5.000  0.0  5.000;  % brayton point
   3   6  0.0  200.0  0.02  10.0  1.0  5.000  -5.000  0.100  -0.100  1.0  1.0  0.0  0.0  0.0  0.0  5.000  0.0  5.000;  % pilgrim
   3   7  0.0  200.0  0.02  10.0  1.0  5.000  -5.000  0.100  -0.100  1.0  1.0  0.0  0.0  0.0  0.0  5.000  0.0  5.000;  % canal
   3   8  0.0  200.0  0.02  10.0  1.0  5.000  -5.000  0.100  -0.100  1.0  1.0  0.0  0.0  0.0  0.0  5.000  0.0  5.000;  % vermont yankee
   3   9  0.0  200.0  0.02  10.0  1.0  5.000  -5.000  0.100  -0.100  1.0  1.0  0.0  0.0  0.0  0.0  5.000  0.0  5.000;  % maine yankee
   3  10  0.0  200.0  0.02  10.0  1.0  5.000  -5.000  0.100  -0.100  1.0  1.0  0.0  0.0  0.0  0.0  5.000  0.0  5.000]; % northfield

% note: the exciter data is spec'd in the exc_st3 format (ieee type st3)

% pss_con format
%    1  type   : pss type (1=speed input, 2=power input)
%    2  bus    : machine number
%    3  gain   : the product of the gain and the washout time constant (pu*s)
%    4  tw     : washout filter time constant (s)
%    5  tn1    : lead time constant (s)
%    6  tn2    : lag time constant (s)
%    7  td1    : lead time constant (s)
%    8  td2    : lag time constant (s)
%    9  vs_max : maximum pss output (pu)
%   10  vs_min : minimum pss output (pu)

% type  num  k*tw  tw  tn1   td1   tn2   td2  vs_max  vs_min
pss_con = [
   1   1   10.0  10.0  5.00  0.60  3.00  0.50  0.200  -0.200;
   1   2    5.0  10.0  5.00  0.40  1.00  0.10  0.200  -0.200;
   1   3    5.0  10.0  3.00  0.20  2.00  0.20  0.200  -0.200;
   1   4   20.0  10.0  1.00  0.10  1.00  0.30  0.200  -0.200;
   1   5   10.0  10.0  1.50  0.20  1.00  0.10  0.200  -0.200;
   1   6   40.0  10.0  0.50  0.10  0.50  0.05  0.200  -0.200;
   1   7   75.0  10.0  0.20  0.02  0.50  0.10  0.200  -0.200;
   1   8   20.0  10.0  1.00  0.20  1.00  0.10  0.200  -0.200;
   1   9   20.0  10.0  1.00  0.50  2.00  0.10  0.200  -0.200;
   1  10   10.0  10.0  1.00  0.05  3.00  0.50  0.200  -0.200];

%-------------------------------------------------------------------------------------------------------------%
% simulation control data

% sw_con format
%    row 1 col1 : simulation start time (s) (cols 2 to 6 zeros)
%          col7 : initial time step (s)
%    row 2 col1 : fault application time (s)
%          col2 : bus number at which fault is applied
%          col3 : bus number defining far end of faulted line
%          col4 : zero sequence impedance in pu on system base
%          col5 : negative sequence impedance in pu on system base
%          col6 : event type - 0 three-phase fault
%               :            - 1 line-to-ground fault
%               :            - 2 line-to-line to ground fault
%               :            - 3 line-to-line fault
%               :            - 4 loss of line with no fault
%               :            - 5 loss of load at bus
%               :            - 6 no fault (do nothing)
%               :            - 7 three-phase fault with no loss of line
%          col7 : time step for fault period (s)
%    row 3 col1 : near end fault clearing time (s) (cols 2 to 6 zeros)
%          col7 : time step for second part of fault (s)
%    row 4 col1 : far end fault clearing time (s) (cols 2 to 6 zeros)
%          col7 : time step for fault cleared simulation (s)
%    row 5 col1 : time to change step length (s)
%          col7 : time step (s)

ne_tstep = 0.002;  % 2 ms, or approximately 1/8th of a cycle
sw_con = [
  0.00   0    0    0    0    0    ne_tstep;  % sets intitial time step
  1.00  16   17    0    0    7    ne_tstep;  % apply three-phase fault at bus 16
  1.10   0    0    0    0    0    ne_tstep;  % clear fault at bus 16
  1.13   0    0    0    0    0    ne_tstep;  % clear remote end
  1.50   0    0    0    0    0    ne_tstep;  % increase time step
  2.00   0    0    0    0    0    ne_tstep;  % increase time step
 10.00   0    0    0    0    0    ne_tstep]; % end simulation

% eof
