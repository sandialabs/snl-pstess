% npcc 48 machine system data
%      Joe Chow, 4/92
% loadflow data
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
bus = [ ...
      1 1.0113  2.63    0.     0.     0.     0.  0. 0.     3  0  0;
   2 1.0070  1.86    0.     0.     0.     0.  0. 0.     3  0  0;
   3 1.0473  1.79    0.     0.     9.0   88.0 0. 0.     3  0  0;
   4 1.0079  1.86    0.     0.     0.     0.  0. 0.     3  0  0;
   5 1.0040  0.04    0.     0.     0.     0.  0. 0.     3  0  0;
   6 1.0073 -0.59    0.     0.   320.0  153.0 0. 0.     3  0  0;
   7 1.0235  0.76    0.     0.   329.0   32.0 0. 0.     3  0  0;
   8 1.0229 -0.54    0.     0.     0.     0.  0. 0.     3  0  0;
   9 1.0186 -1.48    0.     0.   158.0   30.0 0. 0.     3  0  0;
   10 1.0421  5.40    0.     0.     0.     0.  0. 0.     3  0  0;
   11 0.9842  3.92    0.     0.   680.0  103.0 0. 0.     3  0  0;
   12 1.0252  3.43    0.     0.   274.0  115.0 0. 0.     3  0  0;
   13 1.0458  8.15    0.     0.     0.     0.  0. 0.     3  0  0;
   14 1.0400  7.86    0.     0.   248.0   85.0 0. 0.     3  0  0;
   15 1.0294  0.93    0.     0.   309.0  -92.0 0. 0.     3  0  0;
   16 1.0400  1.74    0.     0.   224.0   47.0 0. 0.     3  0  0;
   17 1.0402  0.64    0.     0.   139.0   17.0 0. 0.     3  0  0;
   18 1.0264 -1.09    0.     0.   281.0   76.0 0. 0.     3  0  0;
   19 1.0411  3.83    0.     0.   206.0   28.0 0. 0.     3  0  0;
   20 1.0418  6.51    0.     0.   284.0   27.0 0. 0.     3  0  0;
   21 0.9800 10.70  650.0  216.5    0.     0.  0. 0.     2  999  -999;
   22 0.9900 10.69  632.0  109.8    0.     0.  0. 0.     2  999  -999;
   23 1.0066  9.12  503.0  172.0    0.     0.  0. 0.     2  999  -999;
   24 1.0500 13.51  700.0  250.6    0.     0.  0. 0.     2  999  -999;
   25 1.0600 15.78  560.0  106.1    0.     0.  0. 0.     2  999  -999;
   26 1.0200 13.41  800.0   30.9    0.     0.  0. 0.     2  999  -999;
   27 1.0200  8.68  540.0   42.2    0.     0.  0. 0.     2  999  -999;
   28 1.0513  0.48    0.     0.     0.     0.  0. 0.     3  0 0;
   29 1.0257  0.48    0.     0.     0.     0.  0. 0.     3  0 0;
   30 1.0152 -1.84    0.     0.   322.0    2.0 0. 0.     3  0  0;
   31 0.9951 -1.93    0.     0.   500.0  184.0 0. 0.     3  0  0;
   32 0.9996 -0.33    0.     0.     0.     0.  0. 0.     3  0  0;
   33 1.0023  0.38    0.     0.     0.     0.  0. 0.     3  0  0;
   34 0.9922 -1.61    0.     0.   234.0   84.0 0. 0.     3  0  0;
   35 0.9916 -2.01    0.     0.   522.0  177.0 0. 0.     3  0  0;
   36 0.9800  8.92  554.0  215.3    9.0    5.0 0. 0.     2  999  -999;
   37 1.0375  1.28    0.     0.     0.     0.  0. 0.     3  0  0;
   38 1.0360  2.33    0.     0.     0.     0.  0. 0.     3  0  0;
   39 1.0365  0.76    0.     0.     0.     0.  0. 0.     3  0  0;
   40 1.0029 -3.14    0.     0.     0.     0.  0. 0.     3  0  0;
   41 1.0353 -7.00    0.     0.   700.0  150.0 0. 0.2500 3  0  0;
   42 1.0400 -6.11  370.0   30.8  177.0   17.0 0. 0.     2  999  -999;
   43 1.0454  4.32    0.     0.     0.     0.  0. 0.     3  0  0;
   44 1.0418  2.55    0.     0.     0.     0.  0. 0.     3  0  0;
   45 1.0409 -0.28    0.     0.   256.0   42.0 0. 0.2500 3  0  0;
   46 1.0500  2.49  154.0   67.5  177.0   25.0 0. 0.     2  999  -999;
   47 1.0300  8.20  400.0   54.2  528.0  160.0 0. 0.6000 2  999  -999;
   48 1.0500 11.57  400.0   85.3    0.     0.  0. 0.     2  999  -999;
   49 0.9889 -4.54    0.     0.   208.0   47.0 0. 0.3000 3  0  0;
   50 1.0500  8.65  600.0   75.6    0.     0.  0. 0.     2  999  -999;
   51 1.0400  2.71  141.0   89.8  731.0  100.0 0. 0.2500 2  999  -999;
   52 1.0448 15.28    0.     0.     0.     0.  0. 0.     3  0  0;
   53 1.0300 13.57  833.0  164.6  923.0  250.0 0. 0.     2  999  -999;
   54 1.0408 23.53 1115.0 -150.0   20.0    0.  0. 0.     2  999  -999;
   55 1.0400 22.75 1095.0  197.1  939.0  134.0 0. 0.0200 2  999  -999;
   56 1.0300 23.41  410.0  171.0  172.0   50.0 0. 0.0600 2  999  -999;
   57 1.0300 22.35  298.0  -54.4  240.0   -6.0 0. 0.     2  999  -999;
   58 1.0024 21.06    0.     0.    84.0   40.0 0. 0.0200 3  0  0;
   59 1.0252 18.79    0.     0.   705.0  248.0 0. 0.5000 3  0  0;
   60 1.0400 27.21  410.0   34.5    0.     0.  0. 0.     2  999  -999;
   61 1.0400 23.57  194.0   72.8  237.0   12.0 0. 0.0400 2  999  -999;
   62 1.0085 20.07    0.     0.     0.     0.  0. 0.     3  0  0;
   63 0.9926 11.42    0.     0.     0.     0.  0. 0.     3  0  0;
   64 0.9990  8.35    0.     0.   143.0   25.0 0. 0.     3  0  0;
   65 1.0400 10.34   91.0   63.8   91.0    0.  0. 0.     2  999  -999;
   66 0.9888 10.08    0.     0.     0.     0.  0. 0.     3  0  0;
   67 0.9837  7.72    0.     0.     0.     0.  0. 0.     3  0  0;
   68 0.9900  4.10   66.0   47.8  211.0   41.0 0. 0.     2  999  -999;
   69 1.0073  5.70    0.     0.     0.     0.  0. 0.     3  0  0;
   70 0.9909  4.73    0.     0.     0.     0.  0. 0.     3  0  0;
   71 1.0000  2.48  134.0   13.5  371.0   53.0 0. 0.     2  999 -999;
   72 1.0500  1.85  150.0   40.3    0.0    0.  0. 0.     2  999 -999;
   73 1.0305  0.15    0.     0.     0.     0.  0. 0.     3  0  0;
   74 1.0269  0.07    0.     0.   102.0    0.  0. 0.     3  0  0;
   75 1.0414  1.55    0.     0.     0.     0.  0. 0.     3  0  0;
   76 1.0484  0.88    0.     0.     0.     0.  0. 0.     3  0  0;
   77 1.0355  0.44    0.     0.     0.     0.  0. 0.     3  0  0;
   78 1.0200  0.    632.1   71.9 2000.0  600.0 0. 0.     1  0  0;
   79 1.0500  8.17 1000.0  257.3    0.     0.  0. 0.     2  999  -999;
   80 1.0232  0.63  720.0   -0.0  700.0   25.0 0. 0.     2  999  -999;
   81 1.0192 -3.02    0.     0.     0.     0.  0. 0.     3  0  0;
   82 1.0500  7.57  460.0  125.1    0.     0.  0. 0.     2  999  -999;
   83 1.0213 22.86    0.     0.    80.0   16.0 0. 0.     3  0  0;
   84 1.0248 23.43    0.     0.     0.     0.  0. 0.     3  0  0;
   85 1.0455 34.34    0.     0.    18.0    0.  0. 0.     3  0  0;
   86 1.0000 39.69 1650.0  487.0   50.0   25.0 0. 0.     2  999  -999;
   87 1.0200 32.97    0.     0.     0.     0.  0. 0.     3  0  0;
   88 1.0619 31.21    0.     0.   110.0   55.0 0. 0.     3  0  0;
   89 1.0331 29.02    0.     0.   100.0   70.0 0. 0.     3  0  0;
   90 1.0284 28.94    0.     0.   115.0   45.0 0. 0.     3  0  0;
   91 1.0400 10.61  930.0  947.6 1650.0  500.0 0. 0.     2  999  -999;
   92 1.0400 17.61  450.0  -75.2   10.0  320.0 0. 0.     2  999  -999;
   93 1.0709  7.63    0.     0.   125.0    0.  0. 0.     3  0  0;
   94 1.0554  7.94    0.     0.    95.0    0.  0. 0.     3  0  0;
   95 1.1597 12.58  125.0   13.0   25.0  100.0 0. 0.     3  0  0;
   96 1.1604 13.53   60.0   25.0   20.0    0.  0. 0.     3  0  0;
   97 1.0400 13.16  600.0  319.1  910.0   50.0 0. 0.     2  999  -999;
   98 1.0400 28.67  585.0   48.8    0.0    0.0 0. 0.     2  999  -999;
   99 1.0973 28.03    0.     0.     0.     0.  0. 0.     3  0  0;
   100 1.0830 23.07    0.     0.     0.     0.  0. 0.     3  0  0;
   101 1.050  21.91 1200.0  399.2  350.0  130.0 0. 0.     2  999  -999;
   102 1.0472 22.64    0.     0.     0.     0.  0. 0.     3  0  0;
   103 1.0413 22.52    0.     0.     0.     0.  0. 0.     3  0  0;
   104 1.0567 12.38    0.     0.   200.0   90.0 0. 0.     3  0  0;
   105 1.0559 23.12    0.     0.   325.0    0.  0. 0.     3  0  0;
   106 1.0553 13.80    0.     0.     0.   -12.0 0. 0.     3  0  0;
   107 1.0666 14.87    0.     0.     0.   -19.0 0. 0.     3  0  0;
   108 1.0630 13.01    0.     0.     0.   -13.0 0. 0.     3  0  0;
   109 1.0630 12.99    0.     0.     0.   -13.0 0. 0.     3  0  0;
   110 1.0617 11.71    0.     0.   718.0  125.0 0. 0.     3  0  0;
   111 1.0420  9.97    0.     0.   399.0   45.0 0. 0.     3  0  0;
   112 1.0610 28.34    0.     0.     0.     0.  0. 0.     3  0  0;
   113 1.0229 22.34    0.     0.   206.0    0.  0. 0.     3  0  0;
   114 1.0242 29.52    0.     0.     0.    90.0 0. 0.     3  0  0;
   115 1.0200 32.86  550.0  108.0  100.0    0.  0. 0.     2  999  -999;
   116 1.0128 26.63    0.     0.   200.0    0.  0. 0.     3  0  0;
   117 1.0087 30.39    0.     0.   400.0   97.0 0. 0.     3  0  0;
   118 1.0023 31.38    0.     0.   280.0   86.0 0. 0.     3  0  0;
   119 1.0000 29.78  550.0   94.7  750.0  200.0 0. 0.     2  999  -999;
   120 1.0200 26.28  700.0   43.6  946.0  -40.0 0. 0.     2  999  -999;
   121 1.0400 31.90  145.0  221.3    0.0    0.0 0. 0.     2  999  -999;
   122 1.0500 32.78  155.0   90.7    0.0    0.0 0. 0.     2  999  -999;
   123 1.0100 37.20  655.0   37.0   63.0   10.0 0. 0.     2  999  -999;
   124 1.0386  3.19    0.    22.0  370.0   71.0 0. 0.     3  0  0;
   125 1.0191 -3.14    0.     0.   320.0    5.0 0. 0.     3  0  0;
   126 1.0182 -5.11    0.     0.   450.0    0.  0. 0.     3  0  0;
   127 1.0422  7.55    0.   150.0  200.0    0.  0. 0.     3  0  0;
   128 1.0348  3.86    0.   162.0  970.0   60.0 0. 0.     3  0  0;
   129 1.0227 10.60    0.    50.0  610.0    0.  0. 0.     3  0  0;
   130 1.0200 13.89  490.0 -365.8  180.0    0.  0. 0.     2  999  -999;
   131 1.0248 15.30    0.   900.0 1160.0  468.0 0. 0.     3  0  0;
   132 1.0351 21.03    0.   337.0  210.0    0.  0. 0.     3  0  0;
   133 1.0200 37.34 1700.0  -75.2    0.     0.  0. 0.     2  999  -999;
   134 1.0300 40.69  880.0  103.2    0.     0.  0. 0.     2  999  -999;
   135 1.0200 38.48 2330.0  -42.8  170.0    0.  0. 0.     2  999  -999;
   136 0.9972 29.02    0.     0.   270.0  223.0 0. 0.     3  0  0;
   137 1.0400 28.67   70.0  180.0    0.     0.  0. 0.     2  999  -999;
   138 0.9500  6.87    0.     0.   110.0   64.0 0. 0.     3  0  0;
   139 1.0100 31.91  115.0   23.5    0.     0.  0. 0.     2  999  -999;
   140 1.0412 27.93    0.     0.     0.     0.  0. 0.     3  0  0];
bus(:,4:7) = bus(:,4:7)/100;   % MVA to pu conversion

%   line data
%   format: col 1 -- from bus
%           col 2 -- to bus
%           col 3 -- R
%           col 4 -- X
%           col 5 -- B
%           col 6 -- tap ratio
%           col 7 -- phase shifter angle
line = [...
      1   21  0.        0.020000  0.        1.070000    0;
   2    1  0.000400  0.004300  0.070000  0.          0;
   2   33  0.000700  0.008200  0.140000  0.          0;
   3    2  0.001600  0.043500  0.        1.060000    0;
   3    4  0.001600  0.043500  0.        1.060000    0;
   4    1  0.000400  0.004300  0.070000  0.          0;
   5    4  0.000900  0.010100  0.170000  0.          0;
   5   31  0.000800  0.012900  0.140000  0.          0;
   6    5  0.001800  0.021700  0.370000  0.          0;
   7    6  0.000900  0.009400  0.170000  0.          0;
   8    7  0.000700  0.008900  0.130000  0.          0;
   9    8  0.000700  0.008200  0.130000  0.          0;
   9   30  0.001100  0.013300  0.210000  0.          0;
   10    7  0.001600  0.019500  0.300000  0.          0;
   10   11  0.000700  0.013800  0.        1.060000    0;
   10   22  0.000700  0.014200  0.        1.070000    0;
   11   23  0.000900  0.018000  0.        1.009000    0;
   12    7  0.000800  0.013500  0.250000  0.          0;
   13   12  0.000800  0.014000  0.260000  0.          0;
   13   24  0.        0.014300  0.        1.025000    0;
   14   13  0.000600  0.009600  0.180000  0.          0;
   14   25  0.000500  0.027200  0.        1.000000    0;
   15    7  0.000300  0.005900  0.070000  0.          0;
   15   14  0.002200  0.035000  0.360000  0.          0;
   16   27  0.000600  0.023200  0.        1.025000    0;
   16   29  0.007000  0.008600  0.150000  0.          0;
   17   16  0.003200  0.032300  0.530000  0.          0;
   18    8  0.001300  0.017300  0.320000  0.          0;
   18   17  0.001400  0.014700  0.240000  0.          0;
   19   17  0.004300  0.047400  0.780000  0.          0;
   20   17  0.005700  0.062500  1.030000  0.          0;
   20   19  0.001400  0.015100  0.250000  0.          0;
   20   26  0.000800  0.015600  0.        1.025000    0;
   28   29  0.        0.018100  0.        1.025000    0;
   30   29  0.001300  0.015100  0.260000  0.          0;
   31   30  0.001300  0.021300  0.220000  0.          0;
   32   31  0.000800  0.012800  0.130000  0.          0;
   33   32  0.000200  0.002600  0.040000  0.          0;
   33   36  0.        0.025000  0.        1.070000    0;
   34   33  0.000600  0.009200  0.110000  0.          0;
   35   32  0.000800  0.011200  0.150000  0.          0;
   35   34  0.000400  0.004600  0.080000  0.          0;
   37   29  0.003500  0.041100  0.700000  0.          0;
   38   37  0.001600  0.016300  0.250000  0.          0;
   39   37  0.000800  0.007400  0.480000  0.          0;
   40   37  0.000900  0.031500  0.060000  0.          0;
   40   41  0.        0.016100  0.        0.950000    0;
   42   41  0.001800  0.008900  0.070000  0.          0;
   43   37  0.001300  0.018800  1.310000  0.          0;
   44   40  0.008000  0.053000  0.400000  0.          0;
   44   43  0.        0.011000  0.        0.          0;
   45   41  0.077600  0.224100  0.130000  0.          0;
   45   44  0.        0.020500  0.        0.          0;
   46   45  0.123500  0.354800  0.160000  0.          0;
   47   46  0.063200  0.121000  0.050000  0.          0;
   48   44  0.015000  0.105000  0.790000  0.          0;
   48   47  0.        0.033000  0.        0.          0;
   48  100  0.001000  0.159500  0.010000  0.970000    0;
   49   46  0.164000  0.578000  0.060000  0.          0;
   49   48  0.015600  0.153600  0.200000  1.014000    0;
   50   43  0.002500  0.026800  0.400000  0.          0;
   50   43  0.002500  0.026800  0.400000  0.          0;
   51   45  0.060500  0.189600  0.080000  0.          0;
   51   50  0.        0.020700  0.        0.          0;
   52   50  0.002000  0.022000  1.280000  0.          0;
   53   51  0.084400  0.216100  0.220000  0.          0;
   53   52  0.001400  0.038800  0.        0.          0;
   54   52  0.002000  0.023300  1.120000  1.003000    0;
   54   55  0.        0.013000  0.        1.010000    0;
   55   53  0.030000  0.107000  0.300000  0.          0;
   56   54  0.010000  0.008000  0.070000  0.          0;
   56   57  0.        0.142000  0.        0.980000    0;
   57   55  0.007400  0.040000  0.020000  0.          0;
   58   56  0.002000  0.015000  0.110000  0.          0;
   58   59  0.        0.009300  0.        0.960000    0;
   59   55  0.016000  0.081700  0.040000  0.          0;
   59   57  0.014000  0.068000  0.030000  0.          0;
   60   58  0.003500  0.034400  0.270000  0.          0;
   60   61  0.        0.064000  0.        0.990000    0;
   61   59  0.068000  0.143000  0.060000  0.          0;
   62   54  0.008000  0.069000  0.130000  0.          0;
   62   58  0.002000  0.018000  0.030000  0.          0;
   63   62  0.012900  0.080300  0.150000  0.          0;
   64   63  0.        0.046000  0.        0.          0;
   65   53  0.035000  0.080000  0.050000  0.          0;
   65   64  0.035000  0.080000  0.050000  0.          0;
   67   63  0.012000  0.093000  0.180000  0.          0;
   67   66  0.        0.023100  0.        0.          0;
   68   64  0.098800  0.231000  0.130000  0.          0;
   68   67  0.        0.049000  0.        0.          0;
   69   38  0.004300  0.048500  0.790000  0.          0;
   69   66  0.001800  0.027400  0.270000  0.          0;
   70   67  0.007000  0.061900  0.120000  0.          0;
   71   45  0.175100  0.557500  0.120000  0.          0;
   71   68  0.120300  0.204800  0.070000  0.          0;
   71   69  0.        0.037500  0.        0.          0;
   71   70  0.        0.049000  0.        0.          0;
   72   51  0.0051000  0.0232000  0.30000  0.          0;
   72   71  0.0053000  0.0135000  0.10000  0.          0;
   73   35  0.002300  0.036300  0.380000  0.          0;
   73   39  0.001900  0.018300  0.290000  0.          0;
   73   39  0.001900  0.018300  0.290000  0.          0;
   74   73  0.002200  0.019600  0.340000  0.          0;
   74   73  0.002200  0.019600  0.340000  0.          0;
   75   76  0.000100  0.007400  0.        0.990000    0;
   77   74  0.000300  0.004300  0.080000  0.          0;
   77   76  0.003000  0.006800  1.370000  0.          0;
   78   74  0.000500  0.004500  0.320000  0.          0;
   78   74  0.000500  0.004500  0.320000  0.          0;
   79   78  0.000300  0.015300  0.        0.          0;
   80   78  0.013500  0.058200  0.500000  0.          0;
   81   78  0.000500  0.027600  0.        0.          0;
   82   78  0.000500  0.030800  0.        0.          0;
   83  112  0.018500  0.153700  0.        0.          0;
   83  113  0.003700  0.020200  0.        0.          0;
   84   83  0.005200  0.017400  0.        1.003000    0;
   85   86  0.        0.005700  0.        1.110000    0;
   85  105  0.010000  0.109000  0.180000  0.          0;
   85  105  0.010000  0.109000  0.180000  0.          0;
   85  112  0.006900  0.059500  0.100000  0.          0;
   85  112  0.006900  0.059500  0.100000  0.          0;
   87   85  0.000200  0.024600  0.        0.939000    0;
   88   85  0.001600  0.020200  0.030000  0.          0;
   88   85  0.001600  0.020200  0.030000  0.          0;
   88  105  0.011200  0.099300  0.170000  0.          0;
   88  105  0.011200  0.099300  0.170000  0.          0;
   89   88  0.000500  0.020500  0.        0.          0;
   89  105  0.089400  0.336300  0.        0.          0;
   89  113  0.145600  0.411500  0.        0.          0;
   90   89  0.001300  0.004900  0.        0.          0;
   91  110  0.001700  0.011700  0.        0.          0;
   92   91  0.005800  0.092800  0.        0.          0;
   93   91  0.042100  0.293200  0.        0.          0;
   94   93  0.004200  0.028300  0.        0.          0;
   94  111  0.004900  0.033600  0.        0.          0;
   95   91  0.026100  0.184400  0.        0.          0;
   95   93  0.018700  0.128400  0.        0.          0;
   96   95  0.013600  0.095500  0.        0.          0;
   97   92  0.000800  0.030300  0.        0.          0;
   97   96  0.016000  0.076700  0.        0.          0;
   98   91  0.012000  0.085400  0.        0.          0;
   99   98  0.000300  0.009200  0.        0.          0;
   100   99  0.001000  0.070200  0.        0.          0;
   101  102  0.000700  0.021100  0.        1.045000    0;
   101  103  0.000600  0.015100  0.        1.056000    0;
   102   54  0.001600  0.024800  0.060000  0.          0;
   103   54  0.001600  0.024800  0.060000  0.          0;
   104  101  0.007700  0.077300  0.        0.          0;
   105  101  0.018800  0.149300  0.        0.          0;
   106  104  0.001800  0.018000  0.        0.          0;
   106  105  0.016500  0.115600  0.        0.          0;
   107  101  0.007600  0.075100  0.        0.          0;
   107  105  0.011500  0.110300  0.        0.          0;
   108  101  0.007800  0.077200  0.        0.          0;
   109  101  0.007800  0.078300  0.        0.          0;
   110  101  0.048200  0.263600  0.        0.          0;
   110  104  0.000600  0.006200  0.        0.          0;
   110  107  0.002100  0.018600  0.        0.          0;
   110  108  0.002600  0.019400  0.        0.          0;
   110  109  0.002600  0.019400  0.        0.          0;
   111  105  0.013100  0.080900  0.130000  0.          0;
   111  108  0.005900  0.057700  0.        0.          0;
   111  109  0.005900  0.057700  0.        0.          0;
   112  105  0.008000  0.049700  0.300000  0.          0;
   113  105  0.382400  0.947900  0.        0.          0;
   113  112  0.010100  0.087300  0.        0.          0;
   114   90  0.002400  0.015000  0.        0.          0;
   115   87  0.000100  0.001800  0.020000  0.          0;
   116   84  0.001500  0.104500  0.        0.          0;
   116  115  0.073300  0.979900  0.        0.          0;
   117  114  0.066300  0.369600  0.        0.          0;
   117  115  0.001000  0.011400  0.        0.          0;
   117  116  0.047400  0.263100  0.        0.          0;
   118  115  0.165700  0.155100  0.        0.          0;
   118  116  0.011700  0.074400  0.        0.          0;
   118  117  0.002000  0.019900  0.        0.          0;
   119  114  0.087600  0.254000  0.        0.          0;
   119  115  0.003300  0.039200  0.        0.          0;
   119  117  0.001700  0.018900  0.        0.          0;
   120  116  0.230700  0.704100  0.        0.          0;
   120  117  0.006700  0.079700  0.        0.          0;
   120  118  0.005000  0.062900  0.        0.          0;
   120  119  0.106500  0.349800  0.        0.          0;
   121  114  0.006300  0.080200  0.        0.          0;
   121  115  0.022000  0.298400  0.        0.          0;
   121  116  0.095100  0.457600  0.        0.          0;
   121  117  0.015000  0.086900  0.        0.          0;
   121  118  0.139100  0.479800  0.        0.          0;
   121  119  0.163600  0.580400  0.        0.          0;
   122  116  0.003300  0.128700  0.        0.          0;
   122  117  0.038300  0.088600  0.        0.          0;
   122  118  0.095000  0.507100  0.        0.          0;
   122  121  0.125000  0.469900  0.        0.          0;
   123  118  0.001400  0.017400  0.        0.          0;
   124   75  0.000700  0.017500  1.390000  0.          0;
   125   81  0.000100  0.001100  0.        0.          0;
   125  124  0.002500  0.073000  0.        0.          0;
   126  124  0.        0.083900  0.        0.          0;
   126  125  0.        0.041100  0.        0.          0;
   127  124  0.000400  0.010500  0.720000  0.          0;
   127  126  0.        0.253200  0.        0.          0;
   128  124  0.000500  0.011700  0.840000  0.          0;
   128  125  0.004900  0.263400  0.        0.          0;
   128  126  0.        0.274100  0.        0.          0;
   128  127  0.        0.133900  0.        0.          0;
   129  128  0.        0.138100  0.        0.          0;
   130  128  0.000700  0.017500  1.260000  0.          0;
   130  129  0.000400  0.008500  0.600000  0.          0;
   131  130  0.000200  0.003800  0.540000  0.          0;
   132  127  0.000900  0.022100  1.620000  0.          0;
   132  130  0.000700  0.017300  1.250000  0.          0;
   133  131  0.002000  0.039000  2.770000  0.          0;
   133  132  0.001200  0.029300  2.090000  0.          0;
   134   66  0.0007600  0.0114100  1.160000  0.          0;
   134  132  0.0023500  0.0438800  0.        0.          0;
   135  115  0.        0.630900  0.        0.          0;
   135  122  0.057300  0.758100  0.        0.          0;
   135  131  0.005700  0.069900  0.        0.          0;
   135  132  0.001200  0.028800  2.060000  0.          0;
   135  133  0.000300  0.006300  0.450000  0.          0;
   135  134  0.004600  0.062500  0.        0.          0;
   136  134  0.004600  0.078200  0.790000  0.          0;
   137  136  0.000800  0.023900  0.        1.000000    0;
   138   67  0.006300  0.049100  0.090000  0.          0;
   138  126  0.        0.680500  0.        0.          0;
   138  132  0.033100  0.422400  0.        0.          0;
   138  137  0.255300  0.718000  0.        0.          0;
   139  115  0.        0.124300  0.        0.          0;
   139  122  0.        0.153200  0.        0.          0;
   139  131  0.005900  0.105200  0.        0.          0;
   139  135  0.017700  0.066100  0.        0.          0;
   139  136  0.        0.205400  0.        0.          0;
   140   60  0.003900  0.036300  0.070000  0.          0 ];

%  generator data

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
%      20. saturation factor S(1.0)
%      21. saturation factor S(1.2)
%      22. active power fraction
%      23. reactive power fraction
% note: all machines use transient reactance model
%   machine inertia data and damping coefficients need to be 
%     divided by base mvc
%#  b# MVA   xl  ra  xd    x'd  x"d T'do T"do 
%                    xq    xq'  xq" T'qo T"qo
%                    H     d0    d1 b#   S1 S1.2 frP  frQ

mac_con = [ 
   1  36 100 0.046  0 0 0.081  0   0  0 ... 
      0 0      0   0  0 ...
      30.3     0   0  36    0  0   1.00 1.00;
   2  21 100 0.029  0 0 0.048  0   0  0 ...
      0 0      0   0  0 ...
      34.8     0   0  21    0  0   1.00 1.00;
   3  22 100 0.0295 0 0 0.0436 0   0  0 ...
      0 0      0   0  0 ...
      28.6     0   0  22    0  0   1.00 1.00;
   4  23 100 0.0333 0 0 0.0815 0   0  0 ...
      0 0      0   0  0 ...
      7.34     0   0  23    0  0   0.55 0.55;
   5  23 100 0.0745 0 0 0.182  0   0  0 ...
      0 0      0   0  0 ...
      18.6     0   0  23    0  0   0.45 0.45;
   6  24 100 0.029  0 0 0.048  0   0  0 ...
      0 0      0   0  0 ...
      34.8     0   0  24    0  0   1.00 1.00;
   7  25 100 0.0322 0 0 0.048  0   0  0 ...
      0 0      0   0  0 ...
      26.4     0   0  25    0  0   1.00 1.00;
   8  26 100 0.0313 0 0 0.058  0   0  0 ...
      0 0      0   0  0 ...
      34.2     0   0  26    0  0   1.00 1.00;
   9  27 100 0.030  0 0 0.055  0   0  0 ...
      0 0      0   0  0 ...
      24.30    0   0  27    0  0   1.00 1.00;
   10  42 100 0.0445 0 0 0.0595 0   0  0 ...
      0 0      0   0  0 ...
      18.86    0   0  42    0  0   1.00 1.00;
   11  47 100 0.0782 0 0  0.110 0   0  0 ...
      0  0.110 0   0  0 ...
      15.17    0   0  47    0  0   1.00 1.00;
   12  48 100 0.0782 0 0  0.110 0   0  0 ... 
      0  0     0   0  0 ...
      15.17    0   0  48    0  0   1.00 1.00;
   13  50 100 0.059  0 0  0.078 0   0  0 ...
      0  0     0   0  0 ...
      34.44    0   0  50    0  0   1.00 1.00;
   14  51 100 0.0413 0 0 0.0625 0   0  0 ...
      0 0      0   0  0 ...
      21.46    0   0  51    0  0   1.00 1.00;
   15  53 100 0.0000 0 0 0.070  0   0  0 ...
      0 0      0   0  0 ...
      37.0      0   0  53    0  0   1.00 1.00;
   16  54 100 0.050  0 0 0.0566 0   0  0 ...
      0 0      0   0  0 ...
      51.24    0   0  54    0  0   0.50 0.50;
   17  54 100 0.050 0  0 0.0566 0   0  0 ...
      0 0.0566 0   0  0 ...
      51.24    0   0  54    0  0   0.50 0.50;
   18  55 100 0.0302 0 0 0.0364 0   0  0 ...
      0 0      0   0  0 ...
      59.78    0   0  55    0  0   1.00 1.00;
   19  56 100 0.0515 0 0 0.0715 0   0  0 ...
      0 0.0715 0   0  0 ...
      21.940   0   0  56    0  0   1.00 1.00;
   20  57 100 0.0464 0 0 0.0699 0   0  0 ...
      0 0      0   0  0 ...
      17.71    0   0  57    0  0   1.00 1.00;
   21  60 100 0.0523 0 0 0.0722 0   0  0 ...
      0 0      0   0  0 ...
      21.74    0   0  60    0  0   1.00 1.00;
   22  61 100 0.0947 0 0 0.1248 0   0  0 ...
      0 0.1248 0   0  0 ...
      8.320    0   0  61    0  0   1.00 1.00;
   23  65 100 0.0000 0 0      0.1220 0 0  0 ...
      0      0      0 0  0 ...
      10.71  0 0 65    0  0   1.00 1.00;
   24  68 100 0.0000 0 0      0.2686 0 0.00 0 ...
      0      0 0 0 0 ...
      3.77   0  0 68    0  0   1.00 1.00;
   25  71 100 0.0000 0 0      0.1900 0 0.00 0 ...
      0      0      0 0    0 ...
      8.50   0  0 71    0  0   1.00 1.00;
   26  72 100 0.0000 0 0      0.1185 0 0.00 0 ...
      0      0 0 0 0 ...
      11.77  0  0 72    0  0   1.00 1.00;
   27  78 100 0.0000 0 0      0.0001 0 0.00 0 ...
      0      0 0 0 0 ...
      1000.0 0 0 78    0  0   1.00 1.00;
   28  79 100 0.0280 0 0 0.032  0  0  0  ...
      0 0      0  0  0 ...
      48.0     0  0  79    0  0   1.00 1.00;
   29  80 100 0.0250 0 0 0.047  0  0  0 ...
      0 0      0  0  0 ...
      23.8     0  0  80    0  0   1.00 1.00;
   30  82 100 0.0500 0 0 0.0620 0  0  0 ...
      0 0      0  0  0 ...
      19.600   0  0  82    0  0   1.00 1.00;
   31 101 100 0.0170 0 0 0.0260 0  0  0 ...
      0 0      0  0  0 ...
      55.000   0  0  101   0  0   1.00 1.00;
   32  86 100 0.0100 0 0 0.020  0  0  0 ...
      0 0      0  0  0 ... 
      79.000   0  0  86    0  0   1.00 1.00;
   33  91 100 0.0000 0 0      0.0150 0 0.00 0 ...
      0      0 0 0 0 ...
      98.7   0   0 91    0  0   1.00 1.00;
   34  92 100 0.0000 0 0      0.1000 0 0.00 0 ...
      0      0 0 0 0 ...
      27.0   0  0 92    0   0   1.00 1.00;
   35  97 100 0.0000 0 0      0.0650 0 0.00 0 ...
      0      0 0 0 0 ...
      36.0   0  0 97    0  0   1.00 1.00;
   36  98 100 0.0150 0 0 0.0220 0 0 0 ...
      0 0      0 0 0 ...
      72.0     0  0  98    0  0   1.00 1.00;
   37 115 100 0.0000 0 0      0.0300 0 0.00 0 ...
      0      0 0 0 0 ...
      46.0   0  0 115   0  0   1.00 1.00;
   38 119 100 0.0000 0 0      0.0357 0 0.00 0 ...
      0      0 0 0 0 ...
      30.0   0  0 119   0  0   1.00 1.00;
   39 120 100 0.0000 0 0      0.0001 0 0.00 0 ...
      0      0 0 0 0 ...
      1000   0  0 120   0  0   1.00 1.00;
   40 121 100 0.0000 0 0      0.0250 0 0.00 0 ...
      0      0 0 0 0 ...
      66.0   0  0 121   0  0   1.00 1.00;
   41 122 100 0.0000 0 0      0.0080 0 0.00 0 ...
      0      0 0 0 0 ...
      190.0  0 0 122   0  0   1.00 1.00;
   42 123 100 0.0000 0 0      0.0230 0 0.00 0 ...
      0      0 0 0 0 ...
      33.0   0  0 123   0  0   1.00 1.00;
   43 130 100 0.0000 0 0      0.0306 0 0.00 0 ...
      0      0 0 0 0 ...
      24.0   0  0 130   0  0   1.00 1.00;
   44 133 100 0.0000 0 0      0.0001 0 0.00 0 ...
      0      0 0 0 0 ...
      1000   0  0 133   0  0   1.00 1.00;
   45 134 100 0.0000 0 0      0.0170 0 0.00 0 ...
      0      0 0 0 0 ...
      44.0   0  0 134   0  0   1.00 1.00;
   46 135 100 0.0000 0 0      0.0064 0 0.00 0 ...
      0      0 0 0 0 ...
      115.0  0 0 135   0  0   1.00 1.00;
   47 137 100 0.0000 0 0      0.2100 0 0.00 0 ...
      0      0 0 0 0 ...
      3.5    0  0 137   0  0   1.00 1.00;
   48 139 100 0.0000 0 0      0.0001 0 0.00 0 ...
      0      0 0 0 0 ...
      1000   0  0 139   0  0   1.00 1.00];


%ibus_con=zeros(48,1);
%ibus_con([27 39 44 48])=ones(4,1);



