% switching point generation m file
% syntax [t_switch k_switch] = switch(h)
function [t_switch k_switch] = switch(h) 
  t_switch(1) = 0;     % all time in second+s, start time
  t_switch(2) = 0.5;  % time to apply fault
  t_switch(3) = 0.55; % time to clear fault, 3 cycles
  t_switch(4) = 3.0;  % end time
  if h~=0;stepsize = h;else;stepsize = 0.01;end
  %step to apply fault
  k_switch(1) = round((t_switch(2)-t_switch(1))/stepsize)+1;
  %step to clear fault
  k_switch(2) = round((t_switch(3)-t_switch(1))/stepsize)+1;
  %last step
  k_switch(3) = round((t_switch(4)-t_switch(1))/stepsize)+1;
 
return