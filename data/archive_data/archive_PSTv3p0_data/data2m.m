% Two machine system


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
%  the voltage on bus 2 is adjusted for zero Q at machine 1
bus = [ 1 1.05    0.0   8.0   2.0  0.0  0.0  0.0  0.0  2  999  -999;
        2 1.08103 0.0   8.0   0.0  0.0  0.0  0.0  0.0  1  999  -999;
        3 1.0     0     0     0    16.0 10.0    0.0  0.0  3    0     0];

% line data format
% line: from bus, to bus, resistance(pu), reactance(pu),
%       line charging(pu), tap ratio, phase shifter angle

line = [1  3  0.0  0.01    0.  1. 0.;
        2  3  0.0  0.03999 0.  1. 0.;
        2  3  0.0  0.03999 0.  1. 0.];

% Machine data format
% machine:  1. machine number 
%           2. bus number 
%           3. base mva
%           4. leakage reactance x_l(pu)
%           5. resistance r_a(pu)
%           6. d-axis sychronous reactance x_d(pu)
%           7. d-axis transient reactance x'_d(pu)
%           8. d-axis subtransient reactance x"_d(pu)
%           9. d-axis open-circuit time constant T'_do(sec)
%          10. d-axis open-circuit subtransient time
%                constant T"_do(sec)
%          11. q-axis sychronous reactance x_q(pu)
%          12. q-axis transient reactance x'_q(pu)
%          13. q-axis subtransient reactance x"_q(pu)
%          14. q-axis open-circuit time constant T'_qo(sec)
%          15. q-axis open circuit subtransient time
%                constant T"_qo(sec)
%          16. inertia constant H(sec)
%          17. damping coefficient d_o(pu)
%          18. dampling coefficient d_1(pu)
%          19. bus number
%          20. S(1.0) - saturation factor
%          21. S(1.2) - saturation factor
%  mac_tra models, 
mac_con = [
1 1 991    0.15 0  2.0    0.245 0 5.0   0.0 ...
                   1.91   0.245 0 0.061 0.0...
                   2.875 2.875   0   1    0     0;
2 2 991    0.15 0  2.0    0.245 0 5.0   0.0 ...
                   1.91   0.245 0 0.061 0.0...
                   2.875 2.875   0   1    0     0];



