function [Jac11,Jac12,Jac21,Jac22] = form_jac(V,ang,Y)
% Syntax:  [Jac] = form_jac(V,ang,Y)
%          [Jac11,Jac12,Jac21,Jac22] = form_jac(V,ang,Y)
%
% Purpose: form the Jacobian matrix using sparse matrix techniques
%
% Input:   V        - magnitude of bus voltage
%          ang      - angle(rad) of bus voltage
%          Y        - admittance matrix
%          ang_red  - matrix to eliminate swing bus voltage magnitude and
%                     angle entries
%          volt_red - matrix to eliminate generator bus voltage magnitude
%                     entries
% Output:  Jac      - jacobian matrix
%          Jac11,Jac12,Jac21,Jac22 - submatrices of jacobian matrix
%
% Called By: vsdemo, loadflow

%-----------------------------------------------------------------------------%
% Version history
%
% Version: 2.0
% Author:  Graham Rogers
% Date:    March 1994
% Purpose: eliminated do loops to improve speed
%
% Version: 1.0
% Author:  Kwok W. Cheung, Joe H. Chow
% Date:    March 1991
%-----------------------------------------------------------------------------%

global g;  % declaring struct of global variables

exp_ang = exp(1j*ang);

% voltage in rectangular coordinates
V_rect = V.*exp_ang;
CV_rect = conj(V_rect);
Y_con = conj(Y);

% vector of conjugate currents
i_c = Y_con*CV_rect;

% complex power vector
S = V_rect.*i_c;
S = sparse(diag(S));

Vdia = sparse(diag(V_rect));
CVdia = conj(Vdia);
Vmag = sparse(diag(abs(V)));

S1 = Vdia*Y_con*CVdia;
t1 = ((S+S1)/Vmag)*g.lfac.volt_red';
t2 = (S-S1)*g.lfac.ang_red';

J11 = -g.lfac.ang_red*imag(t2);
J12 = g.lfac.ang_red*real(t1);
J21 = g.lfac.volt_red*real(t2);
J22 = g.lfac.volt_red*imag(t1);

if (nargout == 1)
    Jac11 = [J11 J12;
             J21 J22];
else
    Jac11 = J11;
    Jac12 = J12;
    Jac21 = J21;
    Jac22 = J22;
end

end  % function end

% eof
