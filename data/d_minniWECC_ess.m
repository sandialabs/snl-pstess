% d_minniWECC_base.m
% Based on d_minniWECC_V3C_C3_6.m
% PST data file
%
% MinniWECC system, version 3C, case 3, number 6.  Same as
% d_minniWECC_V3C_C3_5.m but added PSS units as a variable and Alberta
% connection as a variable.  Alberta is disconnected by removing gen 34,
% removing the load on bus 120, and reducing the line impedances to 1%
% (this allows me to maintain the same bus numbers).
%
% MT_P reduces the flow by reducing the power for Colstrip.  Alberta
% import is decreased by decreasing the power for Alberta generation.
% For both cases, the power is equally increased at gen 7 thru 11.  COI
% loading is reduced decreasing bus 43 Pload by COI_P, gen 7 is reduced
% as it is the swing bus.
%
% System bases:  100 MVA, 60 Hz.
%
% Written by Dan Trudnowski, 2011
% Modified by Ryan Elliott, 2018
%
% Case C
% Same settings as d_minniWECC_V3C_C3_5

AlbCon = 1;       % If == 1, Alberta is connected; otherwise 0
COI_P = 20;       % COI_P = amount to reduce COI by
MT_P = 5;         % MT_P = amount to reduce Colstrip by
Alberta_P = 7.5;  % Alberta_P = amount to reduce Alberta import by

% pssGen = [1;13;14];                   % units Dan T. picked
% pssGen = [1;2;3;7;8;10;13;16;24;32];  % hydro units
pssGen = (1:1:34).';                    % all units
% pssGen = 2;                           % one unit for transfer function analysis
% pssGen = [];                          % no PSS

pssType = 1*ones(length(pssGen),1);     % delta-omega type

% pssType() = 1;                        % delta-omega type
% pssType() = 2;                        % delta-p/omega type

%-----------------------------------------------------------------------------%
% Operating parameters
% PDC = DC line real-power flows

PDC = [
    28.5;  % PDCI
    19     % Intermountain
];

F = (MT_P + Alberta_P)/5;  % Amount to add to gen 7 thru 11

if ~(AlbCon==0 || AlbCon==1)
    error('AlbCon must be 0 or 1');
end

% TCSC (not used)
YTCSCf = 1000;           % TCSC fixed series capacitive admittance
YTCSCcMax = 0.2*YTCSCf;  % TCSC controllable admittance limit
TCSC_Gen = [10;21];      % Generators to use for speed feedback
TCSC_Gain = (0/0.1)*60;  % TCSC gain in admittance/speed (pu/pu)

%-----------------------------------------------------------------------------%
% Bus data
% bus data format
% bus:
% col1 number
% col2 voltage magnitude(pu)
% col3 voltage angle(degree)
% col4 p_gen(pu)
% col5 q_gen(pu)
% col6 p_load(pu)
% col7 q_load(pu)
% col8 G shunt(pu)
% col9 B shunt(pu)
% col10 bus_type
%     bus_type
%     - 1, swing bus
%     - 2, generator bus (PV bus)
%     - 3, load bus (PQ bus)
% col11 q_gen_max(pu)
% col12 q_gen_min(pu)
% col13 v_rated (kV)
% col14 v_max  pu
% col15 v_min  pu

bus = [...
% bus    V          Angle       P_gen         Q_gen     P_load    Q_load          G       B      type QgenMax   QgenMin      V_rate vmax vmin
    1    1.025000   27.543403   22.0000-2     0.3847         0         0         0         0    2    45        -45           20     1.5  0.5;  % Gen 1 bus
    2    1.028803   24.644353         0            0   -0.0000    0.0000         0    2.0000    3    0.0        0.0         500     1.5  0.5;  % HV bus for Gen 1, BC
    3    1.025000   17.150142   19.2500-2     0.3251         0         0         0         0    2    37        -37           20     1.5  0.5;  % Gen 2 bus
    4    1.030283   14.113397         0            0    0.0000   -0.0000         0    2.0000    3    0.0        0.0         500     1.5  0.5;  % HV bus for Gen 2, BC
    5    1.025000   35.329953   29.7500-2     3.9890         0         0         0         0    2    37        -37           20     1.5  0.5;  % Gen 3 bus
    6    1.016987   30.377037         0            0    0.0000   -0.0000         0    2.0000    3    0.0        0.0         500     1.5  0.5;  % HV bus for Gen 3, BC
    7    1.038539   -1.196955         0            0   -0.0000   -0.0000         0   28.0000    3    0.0        0.0         500     1.5  0.5;  % HV bus for bus 8 load, BC
    8    1.023717   -4.162308         0            0   44.0000   11.0000         0         0    3    0.0        0.0         110     1.5  0.5;  % Load bus
    9    1.025000   -4.312305   12.7500-2     1.4852         0         0         0         0    2    20        -20           20     1.5  0.5;  % Gen 4 bus
   10    1.028327   -7.820619         0            0   -0.0000   -0.0000         0   19.5000    3    0.0        0.0         500     1.5  0.5;  % HV bus for Gen 4 and bus 11 load, North
   11    1.001939  -12.831969         0            0   54.0000   13.5000         0         0    3    0.0        0.0         110     1.5  0.5;  % Load bus
   12    1.025000    3.216846   27.2000-10    4.6084         0         0         0         0    2    40        -40           20     1.5  0.5;  % Gen 5 bus
   13    1.022422    0.394608         0            0    0.0000    0.0000         0         0    3    0.0        0.0         500     1.5  0.5;  % HV bus for Gen 5, North
   14    1.025000   -3.719713    8.5000       1.7107         0         0         0         0    2    12        -12           20     1.5  0.5;  % Gen 6 bus
   15    1.018393   -8.390425         0            0    0.0000    0.0000         0   14.0000-3  3    0.0        0.0         500     1.5  0.5;  % HV bus for Gen 6 and bus 16 load, North
   16    0.991652  -13.503317         0            0   36.0000    9.0000         0         0    3    0.0        0.0         110     1.5  0.5;  % Load bus
   17    1.025000    0.0        33.4224+F     3.1060         0         0         0         0    1    83        -83           20     1.5  0.5;  % Gen 7 bus (swing bus)
   18    1.025485   -0.370200         0            0   -0.0000   -0.0000         0         0    3    0.0        0.0         500     1.5  0.5;  % HV bus for Gen 7, North
   19    1.025000    8.342673   25.5000+F     2.1979         0         0         0         0    2    36        -36           20     1.5  0.5;  % Gen 8 bus
   20    1.029723    3.269518         0            0    0.0000    0.0000         0    9.2500    3    0.0        0.0         500     1.5  0.5;  % HV bus for Gen 8, 9, and bus 21 load, North
   21    1.003385   -1.727688         0            0    9.0000    2.2500         0         0    3    0.0        0.0         110     1.5  0.5;  % Load bus
   22    1.025000    9.784483    9.7750+F     0.8428         0         0         0         0    2    12.3      -12.3         20     1.5  0.5;  % Gen 9 bus
   23    1.025000    3.962714   32.3000+F     2.3076         0         0         0         0    2    50        -30           20     1.5  0.5;  % Gen 10 bus
   24    1.030654   -0.571788         0            0     PDC(1)   0.0000         0   11.1250-3  3    0.0        0.0         500     1.5  0.5;  % HV bus for Gen 10, 11, and bus 26 load, North
   25    1.025000    4.587579   21.2500+F     1.5300         0         0         0         0    2    30        -30           20     1.5  0.5;  % Gen 11 bus
   26    1.004348   -5.559672         0            0    4.5000    1.1250         0         0    3    0.0        0.0         110     1.5  0.5;  % Load bus
   27    1.025000   11.399068    2.1250       0.5795         0         0         0         0    2     3        -3            20     1.5  0.5;  % Gen 12 bus
   28    1.009307    6.686211         0            0    0.0000    0.0000         0    7.5375    3    0.0        0.0         500     1.5  0.5;  % HV bus for Gen 12 and bus 29 load, North
   29    0.984853    1.491441         0            0   15.7500    3.5000         0         0    3    0.0        0.0         110     1.5  0.5;  % Load bus
   30    1.025000   32.674940   12.7500       0.6454         0         0         0         0    2    15        -15           20     1.5  0.5;  % Gen 13 bus
   31    1.028110   28.048443         0            0   -0.0000    0.0000         0         0    3    0.0        0.0         500     1.5  0.5;  % HV bus for Gen 13, North
   32    1.025000   45.931198   19.2950-MT_P  0.7683         0         0         0         0    2    23        -23           20     1.5  0.5;  % Gen 14 bus
   33    1.021492   41.846406         0            0    0.0000   -0.0000         0    3.0000-3  3    0.0        0.0         500     1.5  0.5;  % HV bus for Gen 14, North
   34    1.025000   -7.645986    4.6750       0.3233         0         0         0         0    2    6.5       -6.5          20     1.5  0.5;  % Gen 15 bus
   35    1.042121  -12.280612         0            0   -0.0000   -0.0000         0    6.2500    3    0.0        0.0         500     1.5  0.5;  % HV bus for Gen 15 and bus 36 load, North
   36    1.016209  -17.155758         0            0    9.0000    2.2500         0         0    3    0.0        0.0         110     1.5  0.5;  % Load bus
   37    1.025000  -15.577639   25.0000       2.1043         0         0         0         0    2    30        -30           20     1.5  0.5;  % Gen 16 bus
   38    1.031759  -17.242934         0            0   -0.0000   -0.0000         0    4.0000-4  3    0.0        0.0         500     1.5  0.5;  % HV bus for Gen 16, nWest
   39    1.025000  -23.010601   20.0000       2.2455         0         0         0         0    2    25        -25           20     1.5  0.5;  % Gen 17 bus
   40    1.028149  -21.187695         0            0   -0.0000   -0.0000         0   12.0000-3  3    0.0        0.0         500     1.5  0.5;  % HV bus, nWest
   41    1.025000  -24.309083   50.0000       4.5989         0         0         0         0    2    52        -52           20     1.5  0.5;  % Gen 18 bus
   42    1.010189  -27.502985         0            0    0.0000    0.0000         0    7.5000-3  3    0.0        0.0         500     1.5  0.5;  % HV bus for Gen 18, nWest
   43    0.987479  -31.638462         0            0  162.4000-COI_P-10  40.6    0         0    3    0.0        0.0         110     1.5  0.5;  % Load bus
   44    1.006861  -28.347986         0            0   -0.0000   -0.0000         0   60.6000-22 3    0.0        0.0         500     1.5  0.5;  % HV bus for bus 43 load and gen 17, nWest
   45    1.025000  -14.896937   26.1000       0.7803         0         0         0         0    2    35        -35           20     1.5  0.5;  % Gen 19 bus
   46    1.025150  -19.782279         0            0    0.0000    0.0000         0         0    3    0.0        0.0         500     1.5  0.5;  % HV bus for Gen 19, sWest
   47    1.025000   -7.755116   19.1400       1.5863         0         0         0         0    2    26        -26           20     1.5  0.5;  % Gen 20 bus
   48    1.025000  -20.207229  110.4900       0.3951         0         0         0         0    2    178       -178          20     1.5  0.5;  % Gen 21 bus
   49    1.026301  -24.267669         0            0    -PDC(1)   0.0000         0   52.6250    3    0.0        0.0         500     1.5  0.5;  % HV bus for Gen 21 and load bus 50, sWest
   50    1.001570  -29.010970         0            0  110.5000   27.6250         0         0    3    0.0        0.0         110     1.5  0.5;  % Load bus
   51    1.025000  -33.645733   10.4400       0.4226         0         0         0         0    2    25        -25           20     1.5  0.5;  % Gen 22 bus
   52    1.023229  -36.384361         0            0    0.0000   -0.0000         0    2.0000    3    0.0        0.0         230     1.5  0.5;  % HV bus for Gen 22, sWest
   53    1.025000  -47.821352    4.3500       1.4722         0         0         0         0    2    35        -35           20     1.5  0.5;  % Gen 23 bus
   54    1.020325  -48.638462         0            0   -0.0000    0.0000         0   19.7125    3    0.0        0.0         500     1.5  0.5;  % HV bus for Gen 23 and load bus 55, sWest
   55    0.995397  -53.439226         0            0   34.8500    8.7125         0         0    3    0.0        0.0         110     1.5  0.5;  % Load bus
   56    0.997811  -33.658851         0            0   68.0000   17.0000         0         0    3    0.0        0.0         110     1.5  0.5;  % Load bus
   57    1.022661  -28.880647         0            0    0.0000    0.0000         0   35.0000    3    0.0        0.0         500     1.5  0.5;  % HV bus for Gen 24, 25, load bus 56, sWest
   58    1.025000  -23.523031   26.1000       1.6470         0         0         0         0    2    32        -32           20     1.5  0.5;  % Gen 24 bus
   59    1.025000  -23.936208   39.1500       2.3825         0         0         0         0    2    52        -52           20     1.5  0.5;  % Gen 25 bus
   60    1.025000  -14.311610   42.3000       4.3692         0         0         0         0    2    47        -47           20     1.5  0.5;  % Gen 26 bus
   61    1.019748  -20.242319         0            0   -0.0000   -0.0000         0    5.0000    3    0.0        0.0         500     1.5  0.5;  % HV bus for Gen 26, East
   62    1.025000  -20.708165  100.8000      11.4939         0         0         0         0    2    113       -113          20     1.5  0.5;  % Gen 27 bus
   63    1.018073  -26.595632         0            0   -0.0000   -0.0000         0   40.1500    3    0.0        0.0         500     1.5  0.5;  % HV bus for Gen 27 and load bus 64, East
   64    0.991319  -31.711866         0            0  120.6000   30.1500         0         0    3    0.0        0.0         110     1.5  0.5;  % Load bus
   65    1.025000  -13.317594   67.5000      -3.3380         0         0         0         0    2    100       -100          20     1.5  0.5;  % Gen 28 bus
   66    1.022786  -17.748900         0            0   -0.0000   -0.0000         0   10.0000-10 3    0.0        0.0         500     1.5  0.5;  % HV bus for Gen 28, East
   67    1.025184  -31.024276         0            0    0.0000    0.0000         0   14.7500    3    0.0        0.0         500     1.5  0.5;  % HV bus for load bus 95, East
   68    1.025000  -10.694976   81.0000       5.5317         0         0         0         0    2    105       -105          20     1.5  0.5;  % Gen 29 bus
   69    1.024907  -15.750558         0            0    PDC(2)    0.0000         0   22.6250    3    0.0        0.0         345     1.5  0.5;  % HV bus for Gen 29 and load bus 70, East
   70    0.998398  -20.796578         0            0   58.5000   14.6250         0         0    3    0.0        0.0         110     1.5  0.5;  % Load bus
   71    1.025000  -11.271883   85.5000      15.5046         0         0         0         0    2    97        -97           20     1.5  0.5;  % Gen 30 bus
   72    1.011887  -17.125082         0            0   -0.0000    0.0000         0   20.2500    3    0.0        0.0         345     1.5  0.5;  % HV bus for Gen 30 and load bus 73, East
   73    0.984908  -22.306219         0            0   81.0000   20.2500         0         0    3    0.0        0.0         110     1.5  0.5;  % Load bus
   74    1.025000   10.665854   18.7000       0.1526         0         0         0         0    2    24        -24           20     1.5  0.5;  % Gen 31 bus
   75    1.023584    5.553005         0            0    0.0000    0.0000         0    2.0000-2  3    0.0        0.0         345     1.5  0.5;  % HV bus for Gen 31, North
   76    1.025000  -19.618753   15.3000       3.2617         0         0         0         0    2    23        -23           20     1.5  0.5;  % Gen 32 bus
   77    1.017847  -24.206789         0            0    0.0000    0.0000         0   11.7500    3    0.0        0.0         230     1.5  0.5;  % (MODIFIED DUE TO EXCESS QGEN AT BUS 76)
   78    0.991086  -29.325313         0            0   27.0000    6.7500         0         0    3    0.0        0.0         110     1.5  0.5;  % (MODIFIED DUE TO EXCESS QGEN AT BUS 76)
   79    1.027607   10.159352         0            0   -0.0000   -0.0000         0    3.0000    3    0.0        0.0         500     1.5  0.5;  % HV node, BC
   80    1.025519   10.485472         0            0   -0.0000    0.0000         0    2.0000    3    0.0        0.0         500     1.5  0.5;  % HV node, BC
   81    1.032711    4.565500         0            0    0.0000   -0.0000         0    1.0000    3    0.0        0.0         500     1.5  0.5;  %
   82    1.038442   -4.197903         0            0         0         0         0    7.0000    3    0.0        0.0         500     1.5  0.5;  %
   83    1.039326   24.713728         0            0   -0.0000   -0.0000         0    3.5000    3    0.0        0.0         500     1.5  0.5;  %
   84    1.042926   19.531768         0            0   -0.0000   -0.0000         0    7.0000    3    0.0        0.0         500     1.5  0.5;  %
   85    1.028456   28.068485         0            0   -0.0000    0.0000         0         0    3    0.0        0.0         500     1.5  0.5;  %
   86    1.045737   -8.757775         0            0    0.0000    0.0000         0    0.5000    3    0.0        0.0         500     1.5  0.5;  % (MODIFIED DUE TO EXCESSIVELY HIGH VOLTAGE)
   87    1.045490  -11.048114         0            0   -0.0000   -0.0000         0    1.0000    3    0.0        0.0         500     1.5  0.5;  %
   88    1.042117  -10.333129         0            0   -0.0000   -0.0000         0    1.5000    3    0.0        0.0         500     1.5  0.5;  %
   89    1.040304  -13.257499         0            0    0.0000    0.0000         0    5.0000    3    0.0        0.0         500     1.5  0.5;  %
   90    1.024546  -19.359512         0            0   -0.0000    0.0000         0    1.0000    3    0.0        0.0         500     1.5  0.5;  %
   91    1.024940  -22.026234         0            0   -0.0000    0.0000         0         0    3    0.0        0.0         500     1.5  0.5;  %
   92    1.007111  -29.900670         0            0   -0.0000    0.0000         0         0    3    0.0        0.0         500     1.5  0.5;  %
   93    1.020164  -27.740709         0            0   -0.0000   -0.0000         0         0    3    0.0        0.0         500     1.5  0.5;  %
   94    1.010306  -26.863679         0            0   -0.0000    0.0000         0         0    3    0.0        0.0         345     1.5  0.5;  % HV node, x-former, sWest
   95    0.998686  -36.067290         0            0   27.0000    6.7500         0         0    3    0.0        0.0         110     1.5  0.5;  % Load bus to HV bus 67
   96    1.013893  -28.454219         0            0   -0.0000    0.0000         0    2.0000    3    0.0        0.0         345     1.5  0.5;  % HV node, x-former, East
   97    1.017390  -22.162030         0            0    0.0000   -0.0000         0         0    3    0.0        0.0         500     1.5  0.5;  %
   98    1.017902  -23.314433         0            0   -0.0000   -0.0000         0         0    3    0.0        0.0         500     1.5  0.5;  %
   99    1.020432  -20.091802         0            0    0.0000    0.0000         0         0    3    0.0        0.0         500     1.5  0.5;  %
  100    1.023005  -12.903233         0            0   -0.0000   -0.0000         0         0    3    0.0        0.0         345     1.5  0.5;  % HV node, North-East
  101    1.023636  -10.052369         0            0   -0.0000    0.0000         0    2.0000    3    0.0        0.0         345     1.5  0.5;  % HV node, North
  102    1.029336  -12.591782         0            0   -0.0000   -0.0000         0    3.0000    3    0.0        0.0         230     1.5  0.5;  % HV node, x-former, North
  103    1.039724   21.979581         0            0    0.0000   -0.0000         0    0.5000    3    0.0        0.0         230     1.5  0.5;  % HV node, x-former, North
  104    1.026476  -11.951378         0            0   -0.0000   -0.0000         0         0    3    0.0        0.0         345     1.5  0.5;  % HV node, x-former, North
  105    1.037317  -11.467126         0            0         0         0         0    1.0000    3    0.0        0.0         500     1.5  0.5;  % HV node, x-former, North
  106    1.018798   39.263278         0            0    0.0000    0.0000         0         0    3    0.0        0.0         500     1.5  0.5;  % HV node, Alberta
  107    1.025000  -12.382405   10.0000       1.1583         0         0         0         0    2    12        -12           20     1.5  0.5;  % Gen 33 bus
  108    1.023023  -17.854752         0            0   -0.0000    0.0000         0    5.5000    3    0.0        0.0         345     1.5  0.5;  % HV bus for Gen 33 and load bus 109, nWest
  109    1.003196  -21.765738         0            0   14.0000    3.5000         0         0    3    0.0        0.0         110     1.5  0.5;  % Load bus
  110    1.045046  -13.854199         0            0   -0.0000    0.0000         0    0.0000    3    0.0        0.0         345     1.5  0.5;  % HV node, x-former, nWest (MODIFIED DUE TO EXCESSIVELY HIGH VOLTAGE)
  111    1.026358    2.204998         0            0   -0.0000    0.0000         0         0    3    0.0        0.0         230     1.5  0.5;  % HV node, x-former, North
  112    0.998502  -44.296629         0            0  131.7500   32.9375         0         0    3    0.0        0.0         110     1.5  0.5;  % Load bus
  113    1.023331  -39.524883         0            0   -PDC(2)    0.0000         0   60.9375    3    0.0        0.0         500     1.5  0.5;  % HV bus for load bus 112, sWest
  114    1.016149  -23.982093         0            0   -0.0000   -0.0000         0    2.0000    3    0.0        0.0         500     1.5  0.5;  % HV node, nWest
  115    1.021067  -12.596982         0            0    0.0000    0.0000         0         0    3    0.0        0.0         500     1.5  0.5;  % HV bus for Gen 20, sWest
  116    1.018586  -34.204512         0            0    0.0000    0.0000         0         0    3    0.0        0.0         500     1.5  0.5;  % HV node, sWest
  117    1.019085   39.953112         0            0   -0.0000    0.0000         0         0    3    0.0        0.0         500     1.5  0.5;  % HV node, Alberta
  118    1.025000   44.923045   [92.0000-Alberta_P   17.3938]*AlbCon   0    0    0         0    2    130       -130          20     1.5  0.5;  % Gen 34 bus
  119    1.019519   40.642458         0            0   -0.0000   -0.0000         0  23.2500*AlbCon 3 0.0        0.0         500     1.5  0.5;  % HV bus for Gen 34 and load bus 120, Alberta
  120    0.994853   35.655219         0            0   [93.0000-11   23.2500-4]*AlbCon  0  0    3    0.0        0.0         110     1.5  0.5;  % Load bus
  121    1.040304  -13.257499         0            0    0.0000    0.0000         0         0    3    0.0        0.0         500     1.5  0.5;  % TCSC bus 1
  122    1.044771  -14.652687         0            0    0.0000    0.0000         0         0    3    0.0        0.0         500     1.5  0.5;  % TCSC bus 2
  123    1.025000   27.543403         0            0    0.0000    0.0000         0         0    3    0.0        0.0          20     1.5  0.5;  % ESS bus for gen at bus 1
  124    1.025000   17.150142         0            0    0.0000    0.0000         0         0    3    0.0        0.0          20     1.5  0.5;  % ESS bus for gen at bus 3
  125    1.025000   35.329953         0            0    0.0000    0.0000         0         0    3    0.0        0.0          20     1.5  0.5;  % ESS bus for gen at bus 5
  126    1.025000   -4.312305         0            0    0.0000    0.0000         0         0    3    0.0        0.0          20     1.5  0.5;  % ESS bus for gen at bus 9
  127    1.025000    3.216846         0            0    0.0000    0.0000         0         0    3    0.0        0.0          20     1.5  0.5;  % ESS bus for gen at bus 12
  128    1.025000   -3.719713         0            0    0.0000    0.0000         0         0    3    0.0        0.0          20     1.5  0.5;  % ESS bus for gen at bus 14
  129    1.025000    0.000000         0            0    0.0000    0.0000         0         0    3    0.0        0.0          20     1.5  0.5;  % ESS bus for gen at bus 17
  130    1.025000    8.342673         0            0    0.0000    0.0000         0         0    3    0.0        0.0          20     1.5  0.5;  % ESS bus for gen at bus 19
  131    1.025000    9.784483         0            0    0.0000    0.0000         0         0    3    0.0        0.0          20     1.5  0.5;  % ESS bus for gen at bus 22
  132    1.025000    3.962714         0            0    0.0000    0.0000         0         0    3    0.0        0.0          20     1.5  0.5;  % ESS bus for gen at bus 23
  133    1.025000    4.587579         0            0    0.0000    0.0000         0         0    3    0.0        0.0          20     1.5  0.5;  % ESS bus for gen at bus 25
  134    1.025000   11.399068         0            0    0.0000    0.0000         0         0    3    0.0        0.0          20     1.5  0.5;  % ESS bus for gen at bus 27
  135    1.025000   32.674940         0            0    0.0000    0.0000         0         0    3    0.0        0.0          20     1.5  0.5;  % ESS bus for gen at bus 30
  136    1.025000   45.931198         0            0    0.0000    0.0000         0         0    3    0.0        0.0          20     1.5  0.5;  % ESS bus for gen at bus 32
  137    1.025000   -7.645986         0            0    0.0000    0.0000         0         0    3    0.0        0.0          20     1.5  0.5;  % ESS bus for gen at bus 34
  138    1.025000  -15.577639         0            0    0.0000    0.0000         0         0    3    0.0        0.0          20     1.5  0.5;  % ESS bus for gen at bus 37
  139    1.025000  -23.010601         0            0    0.0000    0.0000         0         0    3    0.0        0.0          20     1.5  0.5;  % ESS bus for gen at bus 39
  140    1.025000  -24.309083         0            0    0.0000    0.0000         0         0    3    0.0        0.0          20     1.5  0.5;  % ESS bus for gen at bus 41
  141    1.025000  -14.896937         0            0    0.0000    0.0000         0         0    3    0.0        0.0          20     1.5  0.5;  % ESS bus for gen at bus 45
  142    1.025000   -7.755116         0            0    0.0000    0.0000         0         0    3    0.0        0.0          20     1.5  0.5;  % ESS bus for gen at bus 47
  143    1.025000  -20.207229         0            0    0.0000    0.0000         0         0    3    0.0        0.0          20     1.5  0.5;  % ESS bus for gen at bus 48
  144    1.025000  -33.645733         0            0    0.0000    0.0000         0         0    3    0.0        0.0          20     1.5  0.5;  % ESS bus for gen at bus 51
  145    1.025000  -47.821352         0            0    0.0000    0.0000         0         0    3    0.0        0.0          20     1.5  0.5;  % ESS bus for gen at bus 53
  146    1.025000  -23.523031         0            0    0.0000    0.0000         0         0    3    0.0        0.0          20     1.5  0.5;  % ESS bus for gen at bus 58
  147    1.025000  -23.936208         0            0    0.0000    0.0000         0         0    3    0.0        0.0          20     1.5  0.5;  % ESS bus for gen at bus 59
  148    1.025000  -14.311610         0            0    0.0000    0.0000         0         0    3    0.0        0.0          20     1.5  0.5;  % ESS bus for gen at bus 60
  149    1.025000  -20.708165         0            0    0.0000    0.0000         0         0    3    0.0        0.0          20     1.5  0.5;  % ESS bus for gen at bus 62
  150    1.025000  -13.317594         0            0    0.0000    0.0000         0         0    3    0.0        0.0          20     1.5  0.5;  % ESS bus for gen at bus 65
  151    1.025000  -10.694976         0            0    0.0000    0.0000         0         0    3    0.0        0.0          20     1.5  0.5;  % ESS bus for gen at bus 68
  152    1.025000  -11.271883         0            0    0.0000    0.0000         0         0    3    0.0        0.0          20     1.5  0.5;  % ESS bus for gen at bus 71
  153    1.025000   10.665854         0            0    0.0000    0.0000         0         0    3    0.0        0.0          20     1.5  0.5;  % ESS bus for gen at bus 74
  154    1.025000  -19.618753         0            0    0.0000    0.0000         0         0    3    0.0        0.0          20     1.5  0.5;  % ESS bus for gen at bus 76
  155    1.025000  -12.382405         0            0    0.0000    0.0000         0         0    3    0.0        0.0          20     1.5  0.5;  % ESS bus for gen at bus 107
  156    1.025000   44.923045         0            0    0.0000    0.0000         0         0    3    0.0        0.0          20     1.5  0.5]; % ESS bus for gen at bus 118
%  29    0.9786   -2.6100         0         0   15.7500    3.9375         0         0    3    0.0        0.0         110    1.5  0.5;  % Load bus (MODIFIED DUE TO MODERATELY LOW VOLTAGE)
%  77    1.0109  -52.1521         0         0    0.0000    0.0000         0   11.7500    3    0.0        0.0         230    1.5  0.5;  % HV bus for Gen 32 and load bus 78, North
%  78    0.9839  -57.3435         0         0   27.0000    6.7500         0         0    3    0.0        0.0         110    1.5  0.5;  % Load bus
%  83    1.0328   19.8541         0         0   -0.0000   -0.0000         0    7.0000    3    0.0        0.0         500    1.5  0.5;  % (MODIFIED DUE TO EXCESSIVELY HIGH VOLTAGE)
%  86    1.0229  -32.5634         0         0    0.0000    0.0000         0    2.0000    3    0.0        0.0         500    1.5  0.5;  % (MODIFIED DUE TO EXCESSIVELY HIGH VOLTAGE)
%  87    1.0216  -37.1888         0         0   -0.0000   -0.0000         0    3.0000    3    0.0        0.0         500    1.5  0.5;  % (MODIFIED DUE TO EXCESSIVELY HIGH VOLTAGE)
%  88    1.0154  -35.8480         0         0   -0.0000   -0.0000         0    3.0000    3    0.0        0.0         500    1.5  0.5;  % (MODIFIED DUE TO EXCESSIVELY HIGH VOLTAGE)
% 103    1.0203   15.6828         0         0    0.0000   -0.0000         0    1.0000    3    0.0        0.0         230    1.5  0.5;  % HV node, x-former, North (MODIFIED DUE TO EXCESSIVELY HIGH VOLTAGE)
% 110    1.0355  -42.4121         0         0   -0.0000    0.0000         0    2.0000    3    0.0        0.0         345    1.5  0.5;  % HV node, x-former, nWest (MODIFIED DUE TO EXCESSIVELY HIGH VOLTAGE)

%-----------------------------------------------------------------------------%
% Line data

AlbCon2 = 0.01;
if ( AlbCon )
    AlbCon2 = 1;
end

line = [
%         Fbus         Tbus           R            X            Y
            1            2            0    0.0026667            0;  % Gen 1 stepup x-former
            3            4            0    0.0032432            0;  % Gen 2 stepup x-former
            5            6            0    0.0032432            0;  % Gen 3 stepup x-former
            9           10            0        0.006            0;  % Gen 4 stepup x-former
           12           13            0        0.003            0;  % Gen 5 stepup x-former
           14           15            0         0.01            0;  % Gen 6 stepup x-former
           17           18            0    0.0014458            0;  % Gen 7 stepup x-former
           19           20            0    0.0033333            0;  % Gen 8 stepup x-former
           22           20            0    0.0097561            0;  % Gen 9 stepup x-former
           23           24            0       0.0024            0;  % Gen 10 stepup x-former
           25           24            0        0.004            0;  % Gen 11 stepup x-former
           27           28            0         0.04            0;  % Gen 12 stepup x-former
           30           31            0    0.0066667            0;  % Gen 13 stepup x-former
           32           33            0    0.0052174            0;  % Gen 14 stepup x-former
           34           35            0     0.018462            0;  % Gen 15 stepup x-former
           37           38            0        0.004            0;  % Gen 16 stepup x-former 1 (see below for x-former 2)
           39           44            0       0.0048            0;  % Gen 17 stepup x-former
           41           42            0    0.0011538            0;  % Gen 18 stepup x-former
           45           46            0    0.0034286            0;  % Gen 19 stepup x-former
           47          115            0    0.0046154            0;  % Gen 20 stepup x-former
           48           49            0   0.00067416            0;  % Gen 21 stepup x-former
           51           52            0       0.0048            0;  % Gen 22 stepup x-former
           53           54            0    0.0034286            0;  % Gen 23 stepup x-former
           58           57            0      0.00375            0;  % Gen 24 stepup x-former
           59           57            0    0.0023077            0;  % Gen 25 stepup x-former
           60           61            0    0.0025532            0;  % Gen 26 stepup x-former
           62           63            0    0.0010619            0;  % Gen 27 stepup x-former
           65           66            0       0.0012            0;  % Gen 28 stepup x-former
           68           69            0    0.0011429            0;  % Gen 29 stepup x-former
           71           72            0    0.0012371            0;  % Gen 30 stepup x-former
           74           75            0        0.005            0;  % Gen 31 stepup x-former
           76           77            0    0.0054545            0;  % Gen 32 stepup x-former
          107          108            0         0.01            0;  % Gen 33 stepup x-former
          118          119            0   0.00092308            0;  % Gen 34 stepup x-former
            8            7            0      0.00125            0;  % Load bus 8 x-former
           11           10            0    0.0016667            0;  % Load bus 11 x-former
           16           15            0       0.0025            0;  % Load bus 16 x-former
           21           20            0         0.01            0;  % Load bus 21 x-former
           26           24            0         0.02            0;  % Load bus 26 x-former
           29           28            0    0.0057143            0;  % Load bus 29 x-former
           36           35            0         0.01            0;  % Load bus 36 x-former
           43           44            0   0.00043103            0;  % Load bus 43 x-former
           50           49            0   0.00076923            0;  % Load bus 50 x-former
           55           54            0     0.002439            0;  % Load bus 55 x-former
           56           57            0      0.00125            0;  % Load bus 56 x-former
           64           63            0   0.00074627            0;  % Load bus 64 x-former
           95           67            0    0.0033333            0;  % Load bus 95 x-former
           70           69            0    0.0015385            0;  % Load bus 70 x-former
           73           72            0    0.0011111            0;  % Load bus 73 x-former
           78           77            0    0.0033333            0;  % Load bus 78 x-former
          109          108            0        0.005            0;  % Load bus 109 x-former
          112          113            0   0.00064516            0;  % Load bus 112 x-former
          120          119            0    0.0010753            0;  % Load bus 120 x-former
          123            2            0    0.0026667            0;  % ESS bus x-former 123 (gen bus 1)
          124            4            0    0.0032432            0;  % ESS bus x-former 124 (gen bus 3)
          125            6            0    0.0032432            0;  % ESS bus x-former 125 (gen bus 5)
          126           10            0        0.006            0;  % ESS bus x-former 126 (gen bus 9)
          127           13            0        0.003            0;  % ESS bus x-former 127 (gen bus 12)
          128           15            0         0.01            0;  % ESS bus x-former 128 (gen bus 14)
          129           18            0    0.0014458            0;  % ESS bus x-former 129 (gen bus 17)
          130           20            0    0.0033333            0;  % ESS bus x-former 130 (gen bus 19)
          131           20            0    0.0097561            0;  % ESS bus x-former 131 (gen bus 22)
          132           24            0       0.0024            0;  % ESS bus x-former 132 (gen bus 23)
          133           24            0        0.004            0;  % ESS bus x-former 133 (gen bus 25)
          134           28            0         0.04            0;  % ESS bus x-former 134 (gen bus 27)
          135           31            0    0.0066667            0;  % ESS bus x-former 135 (gen bus 30)
          136           33            0    0.0052174            0;  % ESS bus x-former 136 (gen bus 32)
          137           35            0     0.018462            0;  % ESS bus x-former 137 (gen bus 34)
          138           38            0        0.004            0;  % ESS bus x-former 138 (gen bus 37)
          139           44            0       0.0048            0;  % ESS bus x-former 139 (gen bus 39)
          140           42            0    0.0011538            0;  % ESS bus x-former 140 (gen bus 41)
          141           46            0    0.0034286            0;  % ESS bus x-former 141 (gen bus 45)
          142          115            0    0.0046154            0;  % ESS bus x-former 142 (gen bus 47)
          143           49            0   0.00067416            0;  % ESS bus x-former 143 (gen bus 48)
          144           52            0       0.0048            0;  % ESS bus x-former 144 (gen bus 51)
          145           54            0    0.0034286            0;  % ESS bus x-former 145 (gen bus 53)
          146           57            0      0.00375            0;  % ESS bus x-former 146 (gen bus 58)
          147           57            0    0.0023077            0;  % ESS bus x-former 147 (gen bus 59)
          148           61            0    0.0025532            0;  % ESS bus x-former 148 (gen bus 60)
          149           63            0    0.0010619            0;  % ESS bus x-former 149 (gen bus 62)
          150           66            0       0.0012            0;  % ESS bus x-former 150 (gen bus 65)
          151           69            0    0.0011429            0;  % ESS bus x-former 151 (gen bus 68)
          152           72            0    0.0012371            0;  % ESS bus x-former 152 (gen bus 71)
          153           75            0        0.005            0;  % ESS bus x-former 153 (gen bus 74)
          154           77            0    0.0054545            0;  % ESS bus x-former 154 (gen bus 76)
          155          108            0         0.01            0;  % ESS bus x-former 155 (gen bus 107)
          156          119            0   0.00092308            0;  % ESS bus x-former 156 (gen bus 118)
            2           80     0.001679     0.012916            0;
           80           79       0.0014       0.0095            0;
           80           81      0.00145       0.0117            0;
           81            7      0.00145       0.0117            0;
           80            7       0.0033       0.0221            0;
            4           79     0.000309     0.004232            0;
            6           79      0.00265       0.0366            0;
            6           79       0.0028       0.0396            0;
           79            7      0.00125      0.00825            0;
           79            7       0.0025       0.0184            0;
            6           28     0.003125       0.0375            0;  % BC-North (includes x-former)
            7           10     0.000746     0.012954            0;
           18           10       0.0013       0.0251            0;
           18           82       0.0004       0.0063            0;
           82           10     0.000203     0.003377            0;
           10           13      0.00075        0.018            0;
           13           15      0.00075        0.018            0;
           82           24     0.000343     0.007409            0;
           18           28      0.00064      0.01456            0;
           33           83       0.0035       0.0438            0;
           33           85       0.0028      0.03504            0;
           83           85       0.0007      0.00876            0;
           83           84       0.0012       0.0222            0;
           83           84       0.0012       0.0222            0;
           31           84       0.0007       0.0124            0;
           84           20      0.00215      0.03135            0;
           84           28       0.0012       0.0202            0;
           18           20       0.0011       0.0207            0;
           20           24            0        0.002            0;
           24           15      0.00025        0.006            0;
           24           86       0.0007       0.0112            0;
           24           86       0.0011       0.0211            0;
           24           86       0.0011       0.0211            0;
           86           87     0.001075     0.004175            0;
           87          121     0.001075     0.004175            0;  % Added TCSC
           86           88       0.0012       0.0034            0;
           86          121       0.0022       0.0084            0;  % Added TCSC
           88          121       0.0008       0.0076            0;  % Added TCSC
           88          105       0.0039        0.024            0;
           89           35       0.0017       0.0261            0;
           15           35       0.0014       0.0136            0;
           77          102        0.005         0.02            0;
          105          104            0         0.01            0;  % x-former, North
           83          103            0         0.01            0;  % x-former, North
          103          102         0.03         0.12            0;
          102          104            0        0.002            0;  % x-former, North
          104          101       0.0025     0.008333            0;
          104          101       0.0025     0.008333            0;
          101           75       0.0025        0.015            0;
           20          111            0         0.01            0;  % x-former, North
          111           77         0.02         0.24            0;
          101           69        0.002         0.02            0;
          101          100        0.001         0.01            0;
          100           69        0.001         0.01            0;
          104          108      0.01275     0.041667            0;
           89           38        0.001        0.009            0;  % Malin--Round Mt. Ck 1
           89           38       0.0011       0.0089            0;  % Malin--Round MT. Ck 2
           89           90       0.0011       0.0123            0;  % Capt. Jack--Olinda
           89          110            0         0.01            0;  % x-former, North-West
          110          108         0.02     0.066667            0;
           38           40       0.0011       0.0059            0;
           38           40       0.0011       0.0059            0;
           40           44       0.0014       0.0117            0;
           40           44       0.0015       0.0102            0;
           90           44       0.0021       0.0124            0;
           90           44        0.002        0.012            0;
           42           44            0       0.0003            0;
           44          114       0.0007       0.0159            0;
           44          114       0.0006       0.0141            0;
           44          114       0.0015       0.0307            0;
          114           46       0.0015       0.0145            0;
          114           99       0.0008        0.009            0;
           99          115       0.0008       0.0202            0;
           99           46       0.0007       0.0055            0;
          115           46      0.00045      0.01055            0;
           46           49       0.0011       0.0082            0;
           46           49       0.0011        0.008            0;
           46           91      0.00055      0.00335            0;
           91           49      0.00055      0.00335            0;
           49          113       0.0002      0.00565            0;
           49           52        0.002         0.02            0;  % West (includes x-former)
           52           54        0.002         0.02            0;  % West (includes x-former)
          113           57       0.0008       0.0126            0;
          113           57       0.0012       0.0226            0;
          113           57       0.0017       0.0128            0;
          113          116      0.00115      0.01715            0;
          116           57      0.00115      0.01715            0;
           57           94            0         0.01            0;  % West-East (includes x-former)
           57           66       0.0021       0.0143            0;
           57           98     0.000512     0.003589            0;
           98           66     0.000512     0.003589            0;
           57           63       0.0085        0.095            0;
           57           93       0.0009       0.0082            0;
           93           63       0.0009       0.0082            0;
          113           92    0.0010875     0.016819            0;
           92           61    0.0010875     0.016819            0;
           54           61       0.0021     0.024473            0;
           61           63            0         0.01            0;
           94           69        0.015         0.05            0;
           63           66       0.0018       0.0275            0;
           63           97      0.00135        0.016            0;
           97           66      0.00135        0.016            0;
           66           67       0.0009      0.02085            0;
           67           96            0        0.003            0;  % x-former, East
           63           96      0.00375       0.0125            0;
           96           69       0.0075        0.025            0;
           96           72      0.01125       0.0375            0;
           72           69     0.014062     0.046875            0;
          108           69         0.02     0.066667            0;
            6          106       0.0046*AlbCon   0.0641*AlbCon2        0;  % original 0.0641*AlbCon2
          106          119            0*AlbCon     0.02*AlbCon2        0;  % original 0.02*AlbCon2
          106          117            0*AlbCon     0.01*AlbCon2        0;  % original 0.01*AlbCon2
          117          119            0*AlbCon     0.01*AlbCon2        0;  % original 0.01*AlbCon2
           37           90            0        0.004            0;  % Gen 16 stepup x-former 2
           52          113       0.0002      0.00565            0;
           63           66        0.004        0.052            0;
          121          122            0         1/YTCSCf        0;  % Inductor to make TCSC net impedance zero
          122           89            0        -1/YTCSCf        0]; % TCSC fixed capacitance

%----------------------------------------------------------------------------%
% Machine data format
%   Column
%     1. machine number (may be different from bus number),
%     2. bus number,
%     3. machine base mva,
%     4. leakage reactance x_l(pu),
%     5. resistance r_a(pu),
%     6. d-axis sychronous reactance x_d(pu),
%     7. d-axis transient reactance x'_d(pu),
%     8. d-axis subtransient reactance x"_d(pu),
%     9. d-axis open-circuit time constant T'_do(sec),
%    10. d-axis open-circuit subtransient time constant T"_do(sec),
%    11. q-axis sychronous reactance x_q(pu),
%    12. q-axis transient reactance x'_q(pu),
%    13. q-axis subtransient reactance x"_q(pu),
%        NOTE:  PST requires that x"_q = x"_d
%    14. q-axis open-circuit time constant T'_qo(sec),
%        if T'_q0=0, PST sets it to 999 (i.e., removes from model),
%    15. q-axis open circuit subtransient time constant T"_qo(sec),
%        if T"_q0=0, PST sets it to 999,
%    16. inertia constant H(sec),
%    17. damping coefficient d_o(pu),
%    18. damping coefficient d_1(pu), (relative to pmech)
%    19. bus number
%    20. s(1.0)
%    21. s(1.2)
%
% John Undrill reccomended the x'_q, T'_qo, and T"_qo values for the
% salient machines.
d_0S = 0.0;  % Salient damping
d_0R = 0.0;  % Round rotor damping
d_1 = 1;     % Pmech damping (Not used in new mac_sub)

mac_con = [
% 1   2   3      4x    5y   6y   7y    8y    9y     10y    11y   12y   13y   14y   15y    16y   17    18   19  20x    21x
% num bus base   x_l   r_a  x_d  x'_d  x"_d  T'_do  T"_do  x_q   x'_q  x"_q  T'_qo T"_qo  H     d_0   d_1  bus s(1.0) s(1.2)
   1   1  4500   0.17  0.0  1.2  0.3   0.22  6.0    0.025  0.7   0.23  0.22  0.06  0.04   5.0   d_0S  d_1   1  0.05   0.3;   % Hydro, salient pole
   2   3  3700   0.17  0.0  1.2  0.3   0.22  6.0    0.025  0.7   0.23  0.22  0.06  0.04   5.0   d_0S  d_1   3  0.05   0.3;   % Hydro, salient pole
   3   5  3700   0.17  0.0  1.2  0.3   0.22  6.0    0.025  0.7   0.23  0.22  0.06  0.04   5.0   d_0S  d_1   5  0.05   0.3;   % Hydro, salient pole
   4   9  2000   0.17  0.0  2.0  0.2   0.18  7.5    0.025  1.9   0.7   0.18  0.5   0.06   6.5   d_0R  d_1   9  0.05   0.3;   % Gas Turbine, round rotor
   5  12  4000   0.17  0.0  2.0  0.2   0.18  7.5    0.025  1.9   0.7   0.18  0.5   0.06   6.5   d_0R  d_1  12  0.05   0.3;   % Steam, round rotor
   6  14  1200   0.17  0.0  2.0  0.2   0.18  7.5    0.025  1.9   0.7   0.18  0.5   0.06   6.5   d_0R  d_1  14  0.05   0.3;   % Gas Turbine, round rotor
   7  17  8300   0.17  0.0  1.2  0.3   0.22  6.0    0.025  0.7   0.23  0.22  0.06  0.04   5.0   d_0S  d_1  17  0.05   0.3;   % Hydro, salient pole
   8  19  3600   0.17  0.0  1.2  0.3   0.22  6.0    0.025  0.7   0.23  0.22  0.06  0.04   5.0   d_0S  d_1  19  0.05   0.3;   % Hydro, salient pole
   9  22  1230   0.17  0.0  2.0  0.2   0.18  7.5    0.025  1.9   0.7   0.18  0.5   0.06   6.5   d_0R  d_1  22  0.05   0.3;   % Steam, round rotor
  10  23  5000   0.17  0.0  1.2  0.3   0.22  6.0    0.025  0.7   0.23  0.22  0.06  0.04   5.0   d_0S  d_1  23  0.05   0.3;   % Hydro, salient pole
  11  25  3000   0.17  0.0  2.0  0.2   0.18  7.5    0.025  1.9   0.7   0.18  0.5   0.06   6.5   d_0R  d_1  25  0.05   0.3;   % Gas Turbine, round rotor
  12  27   300   0.17  0.0  2.0  0.2   0.18  7.5    0.025  1.9   0.7   0.18  0.5   0.06   6.5   d_0R  d_1  27  0.05   0.3;   % Gas Turbine, round rotor
  13  30  1800   0.17  0.0  1.2  0.3   0.22  6.0    0.025  0.7   0.23  0.22  0.06  0.04   5.0   d_0S  d_1  30  0.05   0.3;   % Hydro, salient pole
  14  32  2300   0.17  0.0  2.0  0.2   0.18  7.5    0.025  1.9   0.7   0.18  0.5   0.06   6.5   d_0R  d_1  32  0.05   0.3;   % Steam, round rotor
  15  34   650   0.17  0.0  2.0  0.2   0.18  7.5    0.025  1.9   0.7   0.18  0.5   0.06   6.5   d_0R  d_1  34  0.05   0.3;   % Gas Turbine, round rotor
  16  37  3000   0.17  0.0  1.2  0.3   0.22  6.0    0.025  0.7   0.23  0.22  0.06  0.04   5.0   d_0S  d_1  37  0.05   0.3;   % Hydro, salient pole
  17  39  2500   0.17  0.0  2.0  0.2   0.18  7.5    0.025  1.9   0.7   0.18  0.5   0.06   6.5   d_0R  d_1  39  0.05   0.3;   % Gas Turbine, round rotor
  18  41  5200   0.17  0.0  2.0  0.2   0.18  7.5    0.025  1.9   0.7   0.18  0.5   0.06   6.5   d_0R  d_1  41  0.05   0.3;   % Gas Turbine, round rotor
  19  45  3500   0.17  0.0  2.0  0.2   0.18  7.5    0.025  1.9   0.7   0.18  0.5   0.06   6.5   d_0R  d_1  45  0.05   0.3;   % Gas Turbine, round rotor
  20  47  2600   0.17  0.0  2.0  0.2   0.18  7.5    0.025  1.9   0.7   0.18  0.5   0.06   6.5   d_0R  d_1  47  0.05   0.3;   % Steam, round rotor
  21  48 17800   0.17  0.0  2.0  0.2   0.18  7.5    0.025  1.9   0.7   0.18  0.5   0.06   6.5   d_0R  d_1  48  0.05   0.3;   % Steam, round rotor? LA
  22  51  2500   0.17  0.0  2.0  0.2   0.18  7.5    0.025  1.9   0.7   0.18  0.5   0.06   6.5   d_0R  d_1  51  0.05   0.3;   % Steam, round rotor
  23  53  3500   0.17  0.0  2.0  0.2   0.18  7.5    0.025  1.9   0.7   0.18  0.5   0.06   6.5   d_0R  d_1  53  0.05   0.3;   % Steam, round rotor? SD
  24  58  3200   0.17  0.0  1.2  0.3   0.22  6.0    0.025  0.7   0.23  0.22  0.06  0.04   5.0   d_0S  d_1  58  0.05   0.3;   % Hydro, salient pole
  25  59  5200   0.17  0.0  2.0  0.2   0.18  7.5    0.025  1.9   0.7   0.18  0.5   0.06   6.5   d_0R  d_1  59  0.05   0.3;   % Gas Turbine, round rotor
  26  60  4700   0.17  0.0  2.0  0.2   0.18  7.5    0.025  1.9   0.7   0.18  0.5   0.06   6.5   d_0R  d_1  60  0.05   0.3;   % Steam, round rotor
  27  62 11300   0.17  0.0  2.0  0.2   0.18  7.5    0.025  1.9   0.7   0.18  0.5   0.06   6.5   d_0R  d_1  62  0.05   0.3;   % Gas Turbine, round rotor
  28  65 10000   0.17  0.0  2.0  0.2   0.18  7.5    0.025  1.9   0.7   0.18  0.5   0.06   6.5   d_0R  d_1  65  0.05   0.3;   % Steam, round rotor
  29  68 10500   0.17  0.0  2.0  0.2   0.18  7.5    0.025  1.9   0.7   0.18  0.5   0.06   6.5   d_0R  d_1  68  0.05   0.3;   % Steam, round rotor
  30  71  9700   0.17  0.0  2.0  0.2   0.18  7.5    0.025  1.9   0.7   0.18  0.5   0.06   6.5   d_0R  d_1  71  0.05   0.3;   % Steam, round rotor
  31  74  2400   0.17  0.0  2.0  0.2   0.18  7.5    0.025  1.9   0.7   0.18  0.5   0.06   6.5   d_0R  d_1  74  0.05   0.3;   % Steam, round rotor
  32  76  2200   0.17  0.0  1.2  0.3   0.22  6.0    0.025  0.7   0.23  0.22  0.06  0.04   5.0   d_0S  d_1  76  0.05   0.3;   % Hydro, salient pole
  33 107  1200   0.17  0.0  2.0  0.2   0.18  7.5    0.025  1.9   0.7   0.18  0.5   0.06   6.5   d_0R  d_1 107  0.05   0.3];  % Steam, round rotor? Reno

if AlbCon
  mac_con = [mac_con;
  34 118 13000   0.17  0.0  2.0  0.2   0.18  7.5    0.025  1.9   0.7   0.18  0.5   0.06   6.5   d_0R  d_1 118  0.05   0.3];  % Steam, round rotor
end

%-----------------------------------------------------------------------------%
% Exciter models
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

my_MW_Ka = 200;  % K_A = 200 (original)

exc_con = [...
%   1   2   3     4          5     6     7    8     9     10    11    12  13  14  15  16  17 18?    19  20
%  type num T_R   K_A        T_A   T_B   T_C  Vrmax Vrmin VImax VImin KJ  KP  qP  KI  XL  KC Efdmax KG  VGmax
    3    1  0.0   my_MW_Ka   0.02  10.0  1.0  5.0   -5.0  0.1   -0.1  1.0 1   0   0   0   0  5.0    0  5.0;
    3    2  0.0   my_MW_Ka   0.02  10.0  1.0  5.0   -5.0  0.1   -0.1  1.0 1   0   0   0   0  5.0    0  5.0;
    3    3  0.0   my_MW_Ka   0.02  10.0  1.0  5.0   -5.0  0.1   -0.1  1.0 1   0   0   0   0  5.0    0  5.0;
    3    4  0.0   my_MW_Ka   0.02  10.0  1.0  5.0   -5.0  0.1   -0.1  1.0 1   0   0   0   0  5.0    0  5.0;
    3    5  0.0   my_MW_Ka   0.02  10.0  1.0  5.0   -5.0  0.1   -0.1  1.0 1   0   0   0   0  5.0    0  5.0;
    3    6  0.0   my_MW_Ka   0.02  10.0  1.0  5.0   -5.0  0.1   -0.1  1.0 1   0   0   0   0  5.0    0  5.0;
    3    7  0.0   my_MW_Ka   0.02  10.0  1.0  5.0   -5.0  0.1   -0.1  1.0 1   0   0   0   0  5.0    0  5.0;
    3    8  0.0   my_MW_Ka   0.02  10.0  1.0  5.0   -5.0  0.1   -0.1  1.0 1   0   0   0   0  5.0    0  5.0;
    3    9  0.0   my_MW_Ka   0.02  10.0  1.0  5.0   -5.0  0.1   -0.1  1.0 1   0   0   0   0  5.0    0  5.0;
    3   10  0.0   my_MW_Ka   0.02  10.0  1.0  5.0   -5.0  0.1   -0.1  1.0 1   0   0   0   0  5.0    0  5.0;
    3   11  0.0   my_MW_Ka   0.02  10.0  1.0  5.0   -5.0  0.1   -0.1  1.0 1   0   0   0   0  5.0    0  5.0;
    3   12  0.0   my_MW_Ka   0.02  10.0  1.0  5.0   -5.0  0.1   -0.1  1.0 1   0   0   0   0  5.0    0  5.0;
    3   13  0.0   my_MW_Ka   0.02  10.0  1.0  5.0   -5.0  0.1   -0.1  1.0 1   0   0   0   0  5.0    0  5.0;
    3   14  0.0   my_MW_Ka   0.02  10.0  1.0  5.0   -5.0  0.1   -0.1  1.0 1   0   0   0   0  5.0    0  5.0;
    3   15  0.0   my_MW_Ka   0.02  10.0  1.0  5.0   -5.0  0.1   -0.1  1.0 1   0   0   0   0  5.0    0  5.0;
    3   16  0.0   my_MW_Ka   0.02  10.0  1.0  5.0   -5.0  0.1   -0.1  1.0 1   0   0   0   0  5.0    0  5.0;
    3   17  0.0   my_MW_Ka   0.02  10.0  1.0  5.0   -5.0  0.1   -0.1  1.0 1   0   0   0   0  5.0    0  5.0;
    3   18  0.0   my_MW_Ka   0.02  10.0  1.0  5.0   -5.0  0.1   -0.1  1.0 1   0   0   0   0  5.0    0  5.0;
    3   19  0.0   my_MW_Ka   0.02  10.0  1.0  5.0   -5.0  0.1   -0.1  1.0 1   0   0   0   0  5.0    0  5.0;
    3   20  0.0   my_MW_Ka   0.02  10.0  1.0  5.0   -5.0  0.1   -0.1  1.0 1   0   0   0   0  5.0    0  5.0;
    3   21  0.0   my_MW_Ka   0.02  10.0  1.0  5.0   -5.0  0.1   -0.1  1.0 1   0   0   0   0  5.0    0  5.0;
    3   22  0.0   my_MW_Ka   0.02  10.0  1.0  5.0   -5.0  0.1   -0.1  1.0 1   0   0   0   0  5.0    0  5.0;
    3   23  0.0   my_MW_Ka   0.02  10.0  1.0  5.0   -5.0  0.1   -0.1  1.0 1   0   0   0   0  5.0    0  5.0;
    3   24  0.0   my_MW_Ka   0.02  10.0  1.0  5.0   -5.0  0.1   -0.1  1.0 1   0   0   0   0  5.0    0  5.0;
    3   25  0.0   my_MW_Ka   0.02  10.0  1.0  5.0   -5.0  0.1   -0.1  1.0 1   0   0   0   0  5.0    0  5.0;
    3   26  0.0   my_MW_Ka   0.02  10.0  1.0  5.0   -5.0  0.1   -0.1  1.0 1   0   0   0   0  5.0    0  5.0;
    3   27  0.0   my_MW_Ka   0.02  10.0  1.0  5.0   -5.0  0.1   -0.1  1.0 1   0   0   0   0  5.0    0  5.0;
    3   28  0.0   my_MW_Ka   0.02  10.0  1.0  5.0   -5.0  0.1   -0.1  1.0 1   0   0   0   0  5.0    0  5.0;
    3   29  0.0   my_MW_Ka   0.02  10.0  1.0  5.0   -5.0  0.1   -0.1  1.0 1   0   0   0   0  5.0    0  5.0;
    3   30  0.0   my_MW_Ka   0.02  10.0  1.0  5.0   -5.0  0.1   -0.1  1.0 1   0   0   0   0  5.0    0  5.0;
    3   31  0.0   my_MW_Ka   0.02  10.0  1.0  5.0   -5.0  0.1   -0.1  1.0 1   0   0   0   0  5.0    0  5.0;
    3   32  0.0   my_MW_Ka   0.02  10.0  1.0  5.0   -5.0  0.1   -0.1  1.0 1   0   0   0   0  5.0    0  5.0;
    3   33  0.0   my_MW_Ka   0.02  10.0  1.0  5.0   -5.0  0.1   -0.1  1.0 1   0   0   0   0  5.0    0  5.0
%   3   X   0.0  7.04   0.4   6.67  1.0 7.57    0    0.2   -0.2  200 4.365 20 4.83 0.091 1.096 6.53 1 6.53  % Typical data from PST manual
];

if AlbCon
exc_con = [exc_con;
    3   34  0.0   my_MW_Ka   0.02  10.0  1.0  5.0   -5.0  0.1   -0.1  1.0 1   0   0   0   0  5.0    0  5.0];
end

clear('my_MW_Ka');

%-----------------------------------------------------------------------------%
% PSS model
% pss_con matrix format
% column  data (units)
%      1  Type input (1=spd, 2=power)
%      2  machine number
%      3  gain K
%      4  Washout time const Tw (sec)
%      5  1st lead time const T1 (sec)
%      6  1st lag time const T2 (sec)
%      7  2nd lead time const T3 (sec)
%      8  2nd lag time const T4 (sec)
%      9  max output limit (pu)
%     10  min output limit (pu)

my_Tw = 10;
my_Kp = 5.0;
for i_pss = 1:length(pssGen)
    pss_con(i_pss,:) = [...
    %   type            gen#           K*Tw         Tw     T1    T2    T3   T4    max   min
        pssType(i_pss), pssGen(i_pss), my_Kp*my_Tw, my_Tw, 0.25, 0.04, 0.2, 0.03, 0.1, -0.1];
end

% pss_con = [];

clear('i_pss','my_Kp','my_Tw');

%-----------------------------------------------------------------------------%
% Governor models
% tg_con matrix format
% column  data (units)
%      1  turbine model number (=1)
%      2  machine number
%      3  speed set point, wf (pu)
%      4  steady state gain, 1/R (pu)
%      5  maximum power order, Tmax (pu on generator base)
%      6  servo time constant, Ts (sec)
%      7  governor time constant, Tc (sec)
%      8  transient gain time constant, T3 (sec)
%      9  HP section time constant, T4 (sec)
%     10  reheater time constant, T5 (sec)

tg_con = [...
%     num  wf 1/R        Tmax  Ts    Tc   T3     T4     T5
1     1   1  20         1.0   0.50  15   1.0    -1.0   0.5;    % Slow Hydro
1     2   1  20         1.0   0.50  15   1.0    -1.0   0.5;    % Slow Hydro
1     3   1  20         1.0   0.50  5    1.0    -1.0   0.5;    % Fast Hydro
%1    4   1  20         1.0   0.50  10   4.0    0      1.0;    % Gas Turbine - NO GOV
%1    5   1  20         1.0   0.50  10   3.0    0      0.05;   % Steam - NO GOV
%1    6   1  20         1.0   0.50  10   4.0    0      1.0;    % Gas Turbine - NO GOV
1     7   1  20         1.0   0.50  15   1.0    -1.0   0.5;    % Slow Hydro
1     8   1  20         1.0   0.50  15   1.0    -1.0   0.5;    % Slow Hydro
%1    9   1  20         1.0   0.10  10   3.0    0      0.01;   % Steam - NO GOV
1    10   1  20         1.0   0.50  15   1.0    -1.0   0.5;    % Slow Hydro
%1   11   1  20         1.0   0.50  10   4.0    0      1.0;    % Gas Turbine - NO GOV
%1   12   1  20         1.0   0.50  10   4.0    0      1.0;    % Gas Turbine - NO GOV
1    13   1  20         1.0   0.50  5    1.0    -1.0   0.5;    % Fast Hydro
1    14   1  20         1.0   0.10  10   3.0    0      0.01;   % Steam
%1   15   1  20         1.0   0.50  10   4.0    0      1.0;    % Gas Turbine - NO GOV
1    16   1  20         1.0   0.50  5    1.0    -1.0   0.5;    % Fast Hydro
%1   17   1  20         1.0   0.10  10   3.0    0      0.01;   % Steam - NO GOV
1    18   1  20*0.25    1.0   0.50  10   4.0    0      1.0;    % Gas Turbine - 75% BASE LOADED
1    19   1  20*0.25    1.0   0.50  10   4.0    0      1.0;    % Gas Turbine - 75% BASE LOADED
%1   20   1  20         1.0   0.10  10   3.0    0      0.01;   % Steam - NO GOV
1    21   1  20         1.0   0.10  10   3.0    0      0.01;   % Steam
%1   22   1  20         1.0   0.10  10   3.0    0      0.01;   % Steam - NO GOV
1    23   1  20         1.0   0.10  10   3.0    0      0.01;   % Steam
1    24   1  20         1.0   0.50  15   1.0    -1.0   0.5;    % Slow Hydro
1    25   1  20*0.25    1.0   0.50  10   4.0    0      1.0;    % Gas Turbine - 75% BASE LOADED
%1   26   1  20         1.0   0.10  10   3.0    0      0.01;   % Steam - NO GOV
1    27   1  20*0.25    1.0   0.50  10   4.0    0      1.0;    % Gas Turbine - 75% BASE LOADED
1    28   1  20*0.2     1.0   0.10  10   3.0    0      0.01;   % Steam - 80% BASE LOADED
%1   29   1  20         1.0   0.10  10   3.0    0      0.01;   % Steam - NO GOV
1    30   1  20*0.25    1.0   0.10  10   3.0    0      0.01;   % Steam - 75% BASE LOADED
%1   31   1  20         1.0   0.10  10   3.0    0      0.01;   % Steam - NO GOV
1    32   1  20         1.0   0.50  5    1.0    -1.0   0.5;    % Fast Hydro
1    33   1  20         1.0   0.10  10   3.0    0      0.01];  % Steam

if AlbCon
    tg_con = [tg_con;
1    34   1  20*0.25    1.0   0.10  10   3.0    0      0.01];  % Steam - 75% BASE LOADED
end

%-----------------------------------------------------------------------------%
% Non-conforming loads
% col 1       bus number
% col 2       fraction const active power load
% col 3       fraction const reactive power load
% col 4       fraction const active current load
% col 5       fraction const reactive current load
load_con = [...
% bus   Pc   Qc   P_Ic Q_Ic
   8    0    0    1    0;   % Load
  11    0    0    1    0;   % Load
  16    0    0    1    0;   % Load
  21    0    0    1    0;   % Load
  24    0    0    1    0;   % PDCI North
  26    0    0    1    0;   % Load
  29    0    0    1    0;   % Load
  36    0    0    1    0;   % Load
  43    0    0    1    0;   % Load
  49    0    0    1    0;   % PDCI South
  50    0    0    1    0;   % Load
  55    0    0    1    0;   % Load
  56    0    0    1    0;   % Load
  64    0    0    1    0;   % Load
  69    0    0    1    0;   % DC North
  70    0    0    1    0;   % Load
  73    0    0    1    0;   % Load
  78    0    0    1    0;   % Load
  95    0    0    1    0;   % Load
 109    0    0    1    0;   % Load
 112    0    0    1    0;   % Load
 113    0    0    1    0;   % DC South
 120    0    0    1    0;   % Load
   2    0    0    1    0;   % BC modulation
   4    0    0    1    0;   % BC modulation
   6    0    0    1    0;   % BC modulation
   7    0    0    1    0;   % BC modulation
  79    0    0    1    0;   % BC modulation
  80    0    0    1    0;   % BC modulation
  10    0    0    1    0;   % North modulation
  13    0    0    1    0;   % North modulation
  15    0    0    1    0;   % North modulation
  18    0    0    1    0;   % North modulation
  20    0    0    1    0;   % North modulation
  28    0    0    1    0;   % North modulation
  31    0    0    1    0;   % North modulation
  33    0    0    1    0;   % North modulation
  35    0    0    1    0;   % North modulation
  75    0    0    1    0;   % North modulation
  77    0    0    1    0;   % North modulation
  83    0    0    1    0;   % North modulation
  88    0    0    1    0;   % North modulation
  89    0    0    1    0;   % North modulation
 101    0    0    1    0;   % North modulation
  38    0    0    1    0;   % nWest modulation
  40    0    0    1    0;   % nWest modulation
  42    0    0    1    0;   % nWest modulation
  44    0    0    1    0;   % nWest modulation
  90    0    0    1    0;   % nWest modulation
 108    0    0    1    0;   % nWest modulation
 114    0    0    1    0;   % nWest modulation
  46    0    0    1    0;   % sWest modulation
  52    0    0    1    0;   % sWest modulation
  54    0    0    1    0;   % sWest modulation
  57    0    0    1    0;   % sWest modulation
  99    0    0    1    0;   % sWest modulation
 115    0    0    1    0;   % sWest modulation
  61    0    0    1    0;   % East modulation
  63    0    0    1    0;   % East modulation
  66    0    0    1    0;   % East modulation
  67    0    0    1    0;   % East modulation
  72    0    0    1    0;   % East modulation
 106    0    0    1    0;   % Alberta modulation
 119    0    0    1    0;   % Alberta modulation
 121    0    0    0    0;   % TCSC bus
 122    0    0    0    0;   % TCSC bus
 123    0    0    1    1;   % ESS bus for gen at bus 1
 124    0    0    1    1;   % ESS bus for gen at bus 3
 125    0    0    1    1;   % ESS bus for gen at bus 5
 126    0    0    1    1;   % ESS bus for gen at bus 9
 127    0    0    1    1;   % ESS bus for gen at bus 12
 128    0    0    1    1;   % ESS bus for gen at bus 14
 129    0    0    1    1;   % ESS bus for gen at bus 17
 130    0    0    1    1;   % ESS bus for gen at bus 19
 131    0    0    1    1;   % ESS bus for gen at bus 22
 132    0    0    1    1;   % ESS bus for gen at bus 23
 133    0    0    1    1;   % ESS bus for gen at bus 25
 134    0    0    1    1;   % ESS bus for gen at bus 27
 135    0    0    1    1;   % ESS bus for gen at bus 30
 136    0    0    1    1;   % ESS bus for gen at bus 32
 137    0    0    1    1;   % ESS bus for gen at bus 34
 138    0    0    1    1;   % ESS bus for gen at bus 37
 139    0    0    1    1;   % ESS bus for gen at bus 39
 140    0    0    1    1;   % ESS bus for gen at bus 41
 141    0    0    1    1;   % ESS bus for gen at bus 45
 142    0    0    1    1;   % ESS bus for gen at bus 47
 143    0    0    1    1;   % ESS bus for gen at bus 48
 144    0    0    1    1;   % ESS bus for gen at bus 51
 145    0    0    1    1;   % ESS bus for gen at bus 53
 146    0    0    1    1;   % ESS bus for gen at bus 58
 147    0    0    1    1;   % ESS bus for gen at bus 59
 148    0    0    1    1;   % ESS bus for gen at bus 60
 149    0    0    1    1;   % ESS bus for gen at bus 62
 150    0    0    1    1;   % ESS bus for gen at bus 65
 151    0    0    1    1;   % ESS bus for gen at bus 68
 152    0    0    1    1;   % ESS bus for gen at bus 71
 153    0    0    1    1;   % ESS bus for gen at bus 74
 154    0    0    1    1;   % ESS bus for gen at bus 76
 155    0    0    1    1;   % ESS bus for gen at bus 107
 156    0    0    1    1];  % ESS bus for gen at bus 118

%-----------------------------------------------------------------------------%
% Load modulation data (sets up modulation of real part of load)
% col 1       load modulation number (Index for b_lmod and lmod_sig)
% col 2       bus number
% col 3       modulation base MVA (MVA)
% col 4       max conductance (pu)
% col 5       min conductance (pu)
% col 6       regulator gain (K)
% col 7       regulator time constant (TR)
% NOTE: This creates b_lmod for the linear analysis.
lmod_con = [...
%num bus  MVA  Max  Min  K  TR
  1   2   100  100 -100  1  0.01;  % BC
  2   4   100  100 -100  1  0.01;  % BC
  3   6   100  100 -100  1  0.01;  % BC
  4   7   100  100 -100  1  0.01;  % BC
  5  79   100  100 -100  1  0.01;  % BC
  6  80   100  100 -100  1  0.01;  % BC
  7  10   100  100 -100  1  0.01;  % North
  8  13   100  100 -100  1  0.01;  % North
  9  15   100  100 -100  1  0.01;  % North
 10  18   100  100 -100  1  0.01;  % North
 11  20   100  100 -100  1  0.01;  % North
 12  24   100  100 -100  1  0.01;  % North PDCI
 13  28   100  100 -100  1  0.01;  % North
 14  31   100  100 -100  1  0.01;  % North
 15  33   100  100 -100  1  0.01;  % North
 16  35   100  100 -100  1  0.01;  % North
 17  75   100  100 -100  1  0.01;  % North
 18  77   100  100 -100  1  0.01;  % North
 19  83   100  100 -100  1  0.01;  % North
 20  88   100  100 -100  1  0.01;  % North
 21  89   100  100 -100  1  0.01;  % North
 22 101   100  100 -100  1  0.01;  % North
 23  38   100  100 -100  1  0.01;  % nWest
 24  40   100  100 -100  1  0.01;  % nWest
 25  42   100  100 -100  1  0.01;  % nWest
 26  44   100  100 -100  1  0.01;  % nWest
 27  90   100  100 -100  1  0.01;  % nWest
 28 108   100  100 -100  1  0.01;  % nWest
 29 114   100  100 -100  1  0.01;  % nWest
 30  46   100  100 -100  1  0.01;  % sWest
 31  49   100  100 -100  1  0.01;  % sWest PDCI
 32  52   100  100 -100  1  0.01;  % sWest
 33  54   100  100 -100  1  0.01;  % sWest
 34  57   100  100 -100  1  0.01;  % sWest
 35  99   100  100 -100  1  0.01;  % sWest
 36 113   100  100 -100  1  0.01;  % sWest DC
 37 115   100  100 -100  1  0.01;  % sWest
 38  61   100  100 -100  1  0.01;  % East
 39  63   100  100 -100  1  0.01;  % East
 40  66   100  100 -100  1  0.01;  % East
 41  67   100  100 -100  1  0.01;  % East
 42  69   100  100 -100  1  0.01;  % East DC
 43  72   100  100 -100  1  0.01;  % East
 44 106   100  100 -100  1  0.01;  % Alberta
 45 119   100  100 -100  1  0.01;  % Alberta
 46   8   100  100 -100  bus(8,6)   0.01;   % Load
 47  11   100  100 -100  bus(11,6)  0.01;   % Load
 48  16   100  100 -100  bus(16,6)  0.01;   % Load
 49  21   100  100 -100  bus(21,6)  0.01;   % Load
 50  26   100  100 -100  bus(26,6)  0.01;   % Load
 51  29   100  100 -100  bus(29,6)  0.01;   % Load
 52  36   100  100 -100  bus(36,6)  0.01;   % Load
 53  43   100  100 -100  bus(43,6)  0.01;   % Load
 54  50   100  100 -100  bus(50,6)  0.01;   % Load
 55  55   100  100 -100  bus(55,6)  0.01;   % Load
 56  56   100  100 -100  bus(56,6)  0.01;   % Load
 57  64   100  100 -100  bus(64,6)  0.01;   % Load
 58  70   100  100 -100  bus(70,6)  0.01;   % Load
 59  73   100  100 -100  bus(73,6)  0.01;   % Load
 60  78   100  100 -100  bus(78,6)  0.01;   % Load
 61  95   100  100 -100  bus(95,6)  0.01;   % Load
 62 109   100  100 -100  bus(109,6) 0.01;   % Load
 63 112   100  100 -100  bus(112,6) 0.01;   % Load
 64 120   100  100 -100  bus(120,6) 0.01;   % Load
 65 123   100  100 -100  1  0.01;  % ESS bus for gen at bus 1
 66 124   100  100 -100  1  0.01;  % ESS bus for gen at bus 3
 67 125   100  100 -100  1  0.01;  % ESS bus for gen at bus 5
 68 126   100  100 -100  1  0.01;  % ESS bus for gen at bus 9
 69 127   100  100 -100  1  0.01;  % ESS bus for gen at bus 12
 70 128   100  100 -100  1  0.01;  % ESS bus for gen at bus 14
 71 129   100  100 -100  1  0.01;  % ESS bus for gen at bus 17
 72 130   100  100 -100  1  0.01;  % ESS bus for gen at bus 19
 73 131   100  100 -100  1  0.01;  % ESS bus for gen at bus 22
 74 132   100  100 -100  1  0.01;  % ESS bus for gen at bus 23
 75 133   100  100 -100  1  0.01;  % ESS bus for gen at bus 25
 76 134   100  100 -100  1  0.01;  % ESS bus for gen at bus 27
 77 135   100  100 -100  1  0.01;  % ESS bus for gen at bus 30
 78 136   100  100 -100  1  0.01;  % ESS bus for gen at bus 32
 79 137   100  100 -100  1  0.01;  % ESS bus for gen at bus 34
 80 138   100  100 -100  1  0.01;  % ESS bus for gen at bus 37
 81 139   100  100 -100  1  0.01;  % ESS bus for gen at bus 39
 82 140   100  100 -100  1  0.01;  % ESS bus for gen at bus 41
 83 141   100  100 -100  1  0.01;  % ESS bus for gen at bus 45
 84 142   100  100 -100  1  0.01;  % ESS bus for gen at bus 47
 85 143   100  100 -100  1  0.01;  % ESS bus for gen at bus 48
 86 144   100  100 -100  1  0.01;  % ESS bus for gen at bus 51
 87 145   100  100 -100  1  0.01;  % ESS bus for gen at bus 53
 88 146   100  100 -100  1  0.01;  % ESS bus for gen at bus 58
 89 147   100  100 -100  1  0.01;  % ESS bus for gen at bus 59
 90 148   100  100 -100  1  0.01;  % ESS bus for gen at bus 60
 91 149   100  100 -100  1  0.01;  % ESS bus for gen at bus 62
 92 150   100  100 -100  1  0.01;  % ESS bus for gen at bus 65
 93 151   100  100 -100  1  0.01;  % ESS bus for gen at bus 68
 94 152   100  100 -100  1  0.01;  % ESS bus for gen at bus 71
 95 153   100  100 -100  1  0.01;  % ESS bus for gen at bus 74
 96 154   100  100 -100  1  0.01;  % ESS bus for gen at bus 76
 97 155   100  100 -100  1  0.01;  % ESS bus for gen at bus 107
 98 156   100  100 -100  1  0.01]; % ESS bus for gen at bus 118

%-----------------------------------------------------------------------------%
% Reactive load modulation matrix.
% col 1       load modulation number (Index for b_rlmod and rlmod_sig)
% col 2       bus number
% col 3       modulation base MVA (MVA)
% col 4       max conductance (pu)
% col 5       min conductance (pu)
% col 6       regulator gain (K)
% col 7       regulator time constant (TR)
% NOTE: This creates b_rlmod for the linear analysis.

rlmod_con = [...
%num bus  MVA  Max  Min  K  TR
  1   2   100  100 -100  1  0.01;  % BC
  2   4   100  100 -100  1  0.01;  % BC
  3   6   100  100 -100  1  0.01;  % BC
  4   7   100  100 -100  1  0.01;  % BC
  5  79   100  100 -100  1  0.01;  % BC
  6  80   100  100 -100  1  0.01;  % BC
  7  10   100  100 -100  1  0.01;  % North
  8  13   100  100 -100  1  0.01;  % North
  9  15   100  100 -100  1  0.01;  % North
 10  18   100  100 -100  1  0.01;  % North
 11  20   100  100 -100  1  0.01;  % North
 12  24   100  100 -100  1  0.01;  % North PDCI
 13  28   100  100 -100  1  0.01;  % North
 14  31   100  100 -100  1  0.01;  % North
 15  33   100  100 -100  1  0.01;  % North
 16  35   100  100 -100  1  0.01;  % North
 17  75   100  100 -100  1  0.01;  % North
 18  77   100  100 -100  1  0.01;  % North
 19  83   100  100 -100  1  0.01;  % North
 20  88   100  100 -100  1  0.01;  % North
 21  89   100  100 -100  1  0.01;  % North
 22 101   100  100 -100  1  0.01;  % North
 23  38   100  100 -100  1  0.01;  % nWest
 24  40   100  100 -100  1  0.01;  % nWest
 25  42   100  100 -100  1  0.01;  % nWest
 26  44   100  100 -100  1  0.01;  % nWest
 27  90   100  100 -100  1  0.01;  % nWest
 28 108   100  100 -100  1  0.01;  % nWest
 29 114   100  100 -100  1  0.01;  % nWest
 30  46   100  100 -100  1  0.01;  % sWest
 31  49   100  100 -100  1  0.01;  % sWest PDCI
 32  52   100  100 -100  1  0.01;  % sWest
 33  54   100  100 -100  1  0.01;  % sWest
 34  57   100  100 -100  1  0.01;  % sWest
 35  99   100  100 -100  1  0.01;  % sWest
 36 113   100  100 -100  1  0.01;  % sWest DC
 37 115   100  100 -100  1  0.01;  % sWest
 38  61   100  100 -100  1  0.01;  % East
 39  63   100  100 -100  1  0.01;  % East
 40  66   100  100 -100  1  0.01;  % East
 41  67   100  100 -100  1  0.01;  % East
 42  69   100  100 -100  1  0.01;  % East DC
 43  72   100  100 -100  1  0.01;  % East
 44 106   100  100 -100  1  0.01;  % Alberta
 45 119   100  100 -100  1  0.01;  % Alberta
 46   8   100  100 -100  bus(8,7)    0.01;   % Load
 47  11   100  100 -100  bus(11,7)   0.01;   % Load
 48  16   100  100 -100  bus(16,7)   0.01;   % Load
 49  21   100  100 -100  bus(21,7)   0.01;   % Load
 50  26   100  100 -100  bus(26,7)   0.01;   % Load
 51  29   100  100 -100  bus(29,7)   0.01;   % Load
 52  36   100  100 -100  bus(36,7)   0.01;   % Load
 53  43   100  100 -100  bus(43,7)   0.01;   % Load
 54  50   100  100 -100  bus(50,7)   0.01;   % Load
 55  55   100  100 -100  bus(55,7)   0.01;   % Load
 56  56   100  100 -100  bus(56,7)   0.01;   % Load
 57  64   100  100 -100  bus(64,7)   0.01;   % Load
 58  70   100  100 -100  bus(70,7)   0.01;   % Load
 59  73   100  100 -100  bus(73,7)   0.01;   % Load
 60  78   100  100 -100  bus(78,7)   0.01;   % Load
 61  95   100  100 -100  bus(95,7)   0.01;   % Load
 62 109   100  100 -100  bus(109,7)  0.01;   % Load
 63 112   100  100 -100  bus(112,7)  0.01;   % Load
 64 120   100  100 -100  bus(120,7)  0.01;   % Load
 65 123   100  100 -100  1  0.01;  % ESS bus for gen at bus 1
 66 124   100  100 -100  1  0.01;  % ESS bus for gen at bus 3
 67 125   100  100 -100  1  0.01;  % ESS bus for gen at bus 5
 68 126   100  100 -100  1  0.01;  % ESS bus for gen at bus 9
 69 127   100  100 -100  1  0.01;  % ESS bus for gen at bus 12
 70 128   100  100 -100  1  0.01;  % ESS bus for gen at bus 14
 71 129   100  100 -100  1  0.01;  % ESS bus for gen at bus 17
 72 130   100  100 -100  1  0.01;  % ESS bus for gen at bus 19
 73 131   100  100 -100  1  0.01;  % ESS bus for gen at bus 22
 74 132   100  100 -100  1  0.01;  % ESS bus for gen at bus 23
 75 133   100  100 -100  1  0.01;  % ESS bus for gen at bus 25
 76 134   100  100 -100  1  0.01;  % ESS bus for gen at bus 27
 77 135   100  100 -100  1  0.01;  % ESS bus for gen at bus 30
 78 136   100  100 -100  1  0.01;  % ESS bus for gen at bus 32
 79 137   100  100 -100  1  0.01;  % ESS bus for gen at bus 34
 80 138   100  100 -100  1  0.01;  % ESS bus for gen at bus 37
 81 139   100  100 -100  1  0.01;  % ESS bus for gen at bus 39
 82 140   100  100 -100  1  0.01;  % ESS bus for gen at bus 41
 83 141   100  100 -100  1  0.01;  % ESS bus for gen at bus 45
 84 142   100  100 -100  1  0.01;  % ESS bus for gen at bus 47
 85 143   100  100 -100  1  0.01;  % ESS bus for gen at bus 48
 86 144   100  100 -100  1  0.01;  % ESS bus for gen at bus 51
 87 145   100  100 -100  1  0.01;  % ESS bus for gen at bus 53
 88 146   100  100 -100  1  0.01;  % ESS bus for gen at bus 58
 89 147   100  100 -100  1  0.01;  % ESS bus for gen at bus 59
 90 148   100  100 -100  1  0.01;  % ESS bus for gen at bus 60
 91 149   100  100 -100  1  0.01;  % ESS bus for gen at bus 62
 92 150   100  100 -100  1  0.01;  % ESS bus for gen at bus 65
 93 151   100  100 -100  1  0.01;  % ESS bus for gen at bus 68
 94 152   100  100 -100  1  0.01;  % ESS bus for gen at bus 71
 95 153   100  100 -100  1  0.01;  % ESS bus for gen at bus 74
 96 154   100  100 -100  1  0.01;  % ESS bus for gen at bus 76
 97 155   100  100 -100  1  0.01;  % ESS bus for gen at bus 107
 98 156   100  100 -100  1  0.01]; % ESS bus for gen at bus 118

%----------------------------------------------------------------------------%
% Monitored lines
%
% When conducting an eigenanalysis, initializing lmon_con causes the
% following variables to be created: c_pf1, c_pf2, c_qf1, c_qf2, c_ilif,
% c_ilit, c_ilmf, c_ilmt, c_ilrf, c_ilrt. The ith row of c_pf1 and c_qf1
% correspond to line(lmon_con(i),:) of the line matrix.

lmon_con = [1:size(line,1)];  % All lines

%----------------------------------------------------------------------------%
% Linear time-varying synchronizing controllers
%
% lsc_con matrix format
%    col   data                                        units
%     1    LTV synchronizing controller index number   integer
%     2    bus number                                  integer
%     3    angle transducer time constant              sec
%     4    use Pade approximation flag (0 = bypass)    binary
%     5    remote signal time delay (Pade)             sec
%     6    local signal time delay (Pade)              sec
%     7    remote angle setpoint rate limit            rad per step
%     8    remote angle setpoint time constant         sec
%     9    local angle setpoint rate limit             rad per step
%    10    local angle setpoint time constant          sec
%    11    local LTV tuning weight (alpha 1)           pu
%    12    local LTV highpass numerator 1              sec
%    13    local LTV highpass denominator 1            sec
%    14    local LTV highpass numerator 2              sec
%    15    local LTV highpass denominator 2            sec
%    16    local LTV lead-lag numerator 1              sec
%    17    local LTV lead-lag denominator 1            sec
%    18    local LTV lead-lag numerator 2              sec
%    19    local LTV lead-lag denominator 2            sec
%    20    local LTV modulation command lower bound    pu on lsc base
%    21    local LTV modulation command upper bound    pu on lsc base
%    22    center LTI tuning weight (alpha 2)          pu
%    23    center LTI highpass numerator 1             sec
%    24    center LTI highpass denominator 1           sec
%    25    center LTI highpass numerator 2             sec
%    26    center LTI highpass denominator 2           sec
%    27    center LTI lead-lag numerator 1             sec
%    28    center LTI lead-lag denominator 1           sec
%    29    center LTI lead-lag numerator 2             sec
%    30    center LTI lead-lag denominator 2           sec
%    31    center LTI modulation command lower bound   pu on lsc base
%    32    center LTI modulation command upper bound   pu on lsc base
%    33    transient stability control gain            pu
%    34    lowpass filter time constant                sec
%    35    total modulation command lower bound        pu on lsc base
%    36    total modulation command upper bound        pu on lsc base

my_Kts = 6.0;   % transient stability gain
my_a1 = 1.00;   % alpha 1 tuning parameter (local vs wide-area)
my_a2 = 0.0015; % alpha 2 tuning parameter (wide-area vs constant reference)

my_paf = 0;     % use pade approximation flag (=0 no pade, =1 use pade)
my_Tdc = 0.00;  % remote time delay for Pade approximation
my_Tdl = 0.00;  % local time delay for Pade approximation

my_Tw1 = 1/(2*pi*0.1);
my_Tw2 = 10;
[my_bw1,my_aw1] = butter(1,(1/(2*pi*my_Tw1))*2*pi,'high','s');  % highpass filters
[my_bw2,my_aw2] = butter(2,(1/(2*pi*my_Tw2))*2*pi,'high','s');

% the filter must be a second-order section
if (length(my_aw1) < 3)
    my_bw1 = [my_bw1, 0.0];
    my_aw1 = [my_aw1, 0.0];
elseif (length(my_aw2) < 3)
    my_bw2 = [my_bw2, 0.0];
    my_aw2 = [my_aw2, 0.0];
end

[my_q1,my_bw1] = deconv(my_bw1,my_aw1);
[my_q2,my_bw2] = deconv(my_bw2,my_aw2);

% checking for incompatible filter specifications
my_flag = false;
if (my_q1 ~= 1 || my_q2 ~= 1)
    my_flag = true;
elseif ((length(my_bw1) ~= 3) || (length(my_aw1) ~= 3))
    my_flag = true;
elseif ((length(my_bw2) ~= 3) || (length(my_aw2) ~= 3))
    my_flag = true;
end

if my_flag
    error('lsc: washout filters specified incorrectly');
end

my_eb = [123:1:156];     % buses with storage
my_neb = length(my_eb);

tmp_lc1 = zeros(my_neb,10);
tmp_lc2 = zeros(my_neb,11);
tmp_lc3 = zeros(my_neb,11);
tmp_lc4 = zeros(my_neb,4);

for i_lsc = 1:my_neb
    tmp_lc1(i_lsc,:) = [...
    %   no     bus           Tv    paf     Tdc     Tdl     rrc  Tsc  rrl  Tsl
        i_lsc, my_eb(i_lsc), 0.02, my_paf, my_Tdc, my_Tdl, 0.0, 0.0, 0.0, 0.0];

    tmp_lc2(i_lsc,:) = [...
    %   a1     Tn1wl      Td1wl      Tn2wl      Td2wl      Tn1cl   Td1cl   Tn2cl   Td2cl   lbl ubl
        my_a1, my_bw1(2), my_aw1(2), my_bw1(3), my_aw1(3), 0.0000, 0.0000, 0.0000, 0.0000, -99, 99];

    tmp_lc3(i_lsc,:) = [...
    %   a2     Tn1wc      Td1wc      Tn2wc      Td2wc      Tn1cc   Td1cc   Tn2cc   Td2cc   lbc ubc
        my_a2, my_bw2(2), my_aw2(2), my_bw2(3), my_aw2(3), 7.2306, 5.6051, 7.2306, 5.6051, -99, 99];

    tmp_lc4(i_lsc,:) = [...
    %   K       Tlp   lbt ubt
        my_Kts, 0.02, -99, 99];
end

lsc_con = [tmp_lc1, tmp_lc2, tmp_lc3, tmp_lc4];
% lsc_con = [];

clear('i_lsc','my_Kts','my_a1','my_a2','my_Tdl','my_Tdc');
clear('my_Tw1','my_Tw2','my_aw1','my_aw2','my_bw1','my_bw2');
clear('my_q1','my_q2','my_flag','tmp_lc1','tmp_lc2','tmp_lc3','tmp_lc4');

%----------------------------------------------------------------------------%
% Energy storage systems
%
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

my_Tdv = 0.00;  % voltage magnitude time delay (Pade)

tmp_ec1 = zeros(my_neb,5);
tmp_ec2 = zeros(my_neb,5);
tmp_ec3 = zeros(my_neb,11);

for i_ess = 1:my_neb
    tmp_ec1(i_ess,:) = [...
    %   no     bus           Tv    paf     Td
        i_ess, my_eb(i_ess), 0.02, my_paf, my_Tdv];

    tmp_ec2(i_ess,:) = [...
    %   P                       E                       Vr   Pr  pf
        0.015*mac_con(i_ess,3), 0.030*mac_con(i_ess,3), 0.95, 2, 0.0];

    tmp_ec3(i_ess,:) = [
    %   Ei   Emn  Emx  Tg    rrp  rrq ilvpl1 zerox brkpt cdi eta
        0.5  0.3  0.7  0.02  30   30   1.22  0.40  0.90   0  0.92];
end

ess_con = [tmp_ec1, tmp_ec2, tmp_ec3];
% ess_con = [];

clear('i_ess','my_eb','my_neb','my_paf','my_Tdv','tmp_ec1','tmp_ec2','tmp_ec3');

%----------------------------------------------------------------------------%
% Switching sequence
%
% Switching file defines the simulation control
% row 1 col1  simulation start time (s) (cols 2 to 6 zeros)
%       col7  initial time step (s)
% row 2 col1  fault application time (s)
%       col2  bus number at which fault is applied
%       col3  bus number defining far end of faulted line
%       col4  zero sequence impedance in pu on system base
%       col5  negative sequence impedance in pu on system base
%       col6  type of fault
%             - 0 three phase
%             - 1 line to ground
%             - 2 line-to-line to ground
%             - 3 line-to-line
%             - 4 loss of line with no fault
%             - 5 loss of load at bus
%             - 6 no action
%             - 7 three phase fault without loss of line
%             - 8 three phase fault with nonzero impedance
%       col7  time step for fault period (s)
%       col8  shunt reactance (pu)
% row 3 col1  near end fault clearing time (s) (cols 2 to 6 zeros)
%       col7  time step for second part of fault (s)
% row 4 col1  far end fault clearing time (s) (cols 2 to 6 zeros)
%       col7  time step for fault cleared simulation (s)
% row 5 col1  time to change step length (s)
%       col7  time step (s)

my_Ts = 1/(8*60);  % pick an integer fraction of the cycle length

% Do nothing
sw_con = [...
    0.0        0    0    0    0    0    my_Ts;   % sets intitial time step
    1.0       18   17    0    0    6    my_Ts;   % brake insertion
    1.5        0    0    0    0    0    my_Ts;   % clear near end
    1.5+1/60   0    0    0    0    0    my_Ts;   % clear remote end
   20.0        0    0    0    0    0    my_Ts];  % end simulation

% SLG fault near Gen 2 in British Columbia (high-voltage side of step-up xfmr)
% 190303 - Ryan's note: Not cooperating!
% sw_con = [...
%     0.0        0    0    0    0    0    my_Ts;   % sets intitial time step
%     1.0        4   79    0    0    1    my_Ts;   % apply fault (three phase)
%     1.0+3/60   0    0    0    0    0    my_Ts;   % clear near end of fault
%     1.0+4/60   0    0    0    0    0    my_Ts;   % clear far end of fault
%     1.5        0    0    0    0    0    my_Ts;   % adjust time step
%    20.0        0    0    0    0    0    my_Ts];  % end simulation

% 3-ph fault on COI without line trip
% sw_con = [...
%     0.0        0    0    0    0    0    my_Ts;   % sets intitial time step
%     1.0       89   90    0    0    0    my_Ts;   % apply fault (three phase)
%     1.1        0    0    0    0    0    my_Ts;   % clear near end of fault
%     1.1+1/60   0    0    0    0    0    my_Ts;   % clear far end of fault
%     1.5        0    0    0    0    0    my_Ts;   % adjust time step
%    20.0        0    0    0    0    0    my_Ts];  % end simulation

% Brake insertion
% sw_con = [...
%     0.0        0    0    0    0    0    my_Ts;   % sets intitial time step
%     1.0       18   17  1/14   0    8    my_Ts;   % brake insertion
%     1.5        0    0    0    0    0    my_Ts;   % clear near end
%     1.5+1/60   0    0    0    0    0    my_Ts;   % clear remote end
%    20.0        0    0    0    0    0    my_Ts];  % end simulation

% Palo Verde drop (gen 26)
% sw_con = [...
%     0.0        0    0    0    0    0    my_Ts;   % sets intitial time step
%     1.0       60   61    0    0    4    my_Ts;   % gen trip (no fault)
%     1.0+1/60   0    0    0    0    0    my_Ts;   % adjust time step
%    40.0        0    0    0    0    0    my_Ts];  % end simulation

% Colstrip drop (gen 14)
% sw_con = [...
%     0.0        0    0    0    0    0    my_Ts;   % sets intitial time step
%     1.0       32   33    0    0    4    my_Ts;   % gen trip (no fault)
%     1.0+1/60   0    0    0    0    0    my_Ts;   % adjust time step
%    40.0        0    0    0    0    0    my_Ts];  % end simulation

% Line 6-106 (Alberta)
% sw_con = [...
%     0.0        0    0    0    0    0    my_Ts;   % sets intitial time step
%     1.0        6  106    0    0    0    my_Ts;   % apply fault (three phase)
%     1.1        0    0    0    0    0    my_Ts;   % clear near end of fault
%     1.1+1/60   0    0    0    0    0    my_Ts;   % clear far end of fault
%     1.5        0    0    0    0    0    my_Ts;   % adjust time step
%    20.0        0    0    0    0    0    my_Ts];  % end simulation

clear('my_Ts');

% eof
