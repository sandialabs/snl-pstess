function [Jac11,Jac12,Jac21,Jac22]=form_jac(V,ang,Y,ang_red,volt_red)
% Syntax:  [Jac] = form_jac(V,ang,Y,ang_red,volt_red)
%          [Jac11,Jac12,Jac21,Jac22] = form_jac(V,ang,Y,...
%                                      ang_red,volt_red)
%
% Purpose: form the Jacobian matrix using sparse matrix techniques
%
% Input:   V        - magnitude of bus voltage
%          ang      - angle(rad) of bus voltage
%          Y        - admittance matrix
%          ang_red  - matrix to eliminate swing bus  voltage magnitude and angle 
%                     entries
%          volt_red - matrix to eliminate generator bus voltage magnitude
%                     entries
% Output:  Jac      - jacobian matrix
%          Jac11,Jac12,Jac21,Jac22 - submatrices of 
%                                      jacobian matrix  
% See also:   
%
% Calls:
%
% Called By:   vsdemo loadflow

% (c) Copyright 1991-1996 Joe H. Chow - All Rights Reserved
%
% History (in reverse chronological order)
% Version:   2.0
% Author:    Graham Rogers
% Date:      March 1994
% Purpose:   eliminated do loops to improve speed
% Version:   1.0
% Author:    Kwok W. Cheung, Joe H. Chow
% Date:      March 1991
%
% ***********************************************************
jay = sqrt(-1);
exp_ang = exp(jay*ang);
% Voltage rectangular coordinates
V_rect = V.*exp_ang;
CV_rect=conj(V_rect);
Y_con = conj(Y);
%vector of conjugate currents
i_c=Y_con*CV_rect;
% complex power vector
S=V_rect.*i_c;
S=sparse(diag(S));
Vdia=sparse(diag(V_rect));
CVdia=conj(Vdia);
Vmag=sparse(diag(abs(V)));
S1=Vdia*Y_con*CVdia;
t1=((S+S1)/Vmag)*volt_red';
t2=(S-S1)*ang_red';
J11=-ang_red*imag(t2);
J12=ang_red*real(t1);
J21=volt_red*real(t2);
J22=volt_red*imag(t1);
if nargout > 3
   Jac11 = J11; clear J11
   Jac12 = J12; clear J12
   Jac21 = J21; clear J21
   Jac22 = J22; clear J22
else
   Jac11 = [J11 J12;
	       J21 J22];
end

