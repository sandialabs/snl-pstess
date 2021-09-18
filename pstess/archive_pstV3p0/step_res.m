function [res,t] = step_res(a,b,c,d,v_in,t_max)
% Calculation of Step Response from State Matrix 
% 4:46 PM 15/08/97
% Syntax: [res, t] = step_res(a,b,c,d,v_in,t_max)
% Inputs: 
%         a,b,c,d - state matrices
%         v_in the magnitude of the input disturbance
%         t_max the maximum value of the response in sec
% Outputs: res  - step response matrix
%           t   - time vector
% Use plot(t,res) to view the response
% Author: Graham Rogers
% Date:   August 1997
% (c) Copyright Joe Chow/Cherry Tree Scientific Software 1997
% All Rights Reserved

% check input compatibility
n_in = length(v_in);
[mb nb] = size(b);
[ma na] = size(a);
[mc,nc] = size(c);
[md,nd] = size(d);

if n_in~=nb
   error('the length of the input vector must be the same as the number of columns in b')
end
if ma~=na
   error('the a matrix must be square')
end
if mb~=ma
   error('the b matrix and a matrix must have the same number of rows')
end
if nc~=na
   error('the c matrix must have the same number as columns as the a matrix')
end
if md~=mc
   error('the d matrix must have the same number of rows as the c matrix')
end
if nd~=nb
   error('the d matrix must have the same number of columns as the b matrix')
end

% find the eigenvalues of a
l = eig(a);


% determine the required time step
t_stepf = 0.01;
t_stepr = 0.01;
comp_idx = find(abs(imag(l))>0);
if ~isempty(comp_idx)
   freq = abs(imag(l(comp_idx)))/2/pi;
   freq_max = max(freq);
   t_stepf = 1/10/freq_max;
end

real_idx = find(abs(imag(l))==0);
if ~isempty(real_idx)
   alpha_max = max(abs(real(l(real_idx))));
   t_stepr = 1/5/alpha_max;
end
t_step = min(t_stepf,t_stepr);
if t_step == 0
   num_steps = 500;
else
   num_steps = round(t_max/t_step);
end
t_step = t_max/num_steps;

% define time vector

t = 0:t_step:t_max;
b = inv(a)*b*v_in;
a = a*t_step;
a = expm(a);
x = zeros(ma,num_steps+1);
a1 = eye(na) - a;
for k = 2:num_steps+1
   x(:,k) = a*x(:,k-1) - a1*b;
end
res = c*x + d*v_in*ones(size(t));
plot(t,res)
return