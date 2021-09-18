% Version: $Id: convert.pld,v 1.4 1996/11/18 18:44:00 jservice Exp $
% Command-Line: convert.pld bruce1a.raw
% Date: Mon Nov 18 16:26:58 1996
%
% PTI data file header:
%  0  100.0
% BRUCE UNIT 3 FULL LOAD REJECTION
% TITLE OF STUDY : EXISTING GOVERNOR
%
% PTI generator and bus data variables:
% $bus_number, $v_bus_pu, $v_ang_deg, $p_gen/$mva_base, $q_gen/$mva_base,   
%p_load_MW/$mva_base, $q_load_MVAR/$mva_base, $g_shunt_MW, $b_shunt_MVAR,   
%m_bus_id{$bus_code}, $q_max/$mva_base, $q_min/$mva_base, $v_base_kv,   
%vmaxpu_def, $vminpu_def
%
bus = [    1 1.05   1.796  8      3.75  0     0     0 0 2 3.9 3.75    18.5 1.1 0.9 ;
           2 1      0      0      0     7.449 1.241 0 0 1 10  -10    500   1.1 0.9 ;
           3 1.05  -3.302  0      0     0     0     0 0 3 0    0     500   1.1 0.9 ;
           4 1.045 -4.068  0      0     0     0     0 0 3 0    0     500   1.1 0.9 ;
           5 1.05  -3.302  0      0     0     0     0 0 3 0    0     500   1.1 0.9 ;
           6 1.007 -2.18   0      0     0.258 0.132 0 0 3 0    0      13.8 1.1 0.9 ;
           7 1.041  0.9821 0      0     0     0     0 0 3 0    0      18.5 1.1 0.9 ;
           8 1.007 -2.18   0      0     0.258 0.132 0 0 3 0    0      13.8 1.1 0.9 ;
           9 1.045 -4.068  0      0     0     0     0 0 3 0    0     230   1.1 0.9 ];
%
% PTI branch data variables:
% $bus_from, $bus_to, $r_branch_pu, $x_branch_pu, $b_charge_admit_pu,   
%t_o_n_ratio, $t_phase_shift, $tap_max, $tap_min, $tap_incr,   
%i_rate_a/$i_base
line = [2 4 0.00028 0.02305 0.01     0      0 0   0   0       0 ;
        3 5 0.0001  0.001   0.00205  0      0 0   0   0       0 ;
        4 5 0.0001  0.002   0.00205  0      0 0   0   0       0 ;
        1 5 0.0002  0.0137  0        0.9625 0 1.1 0.9 0.00625 0 ;
        4 9 0.0002  0.0206  0        1      0 1.1 0.9 0.00625 0 ;
        6 7 0.0097  0.2291  0        1      0 1.1 0.9 0.00625 0 ;
        8 7 0.0097  0.2291  0        1      0 1.1 0.9 0.00625 0 ] ;
