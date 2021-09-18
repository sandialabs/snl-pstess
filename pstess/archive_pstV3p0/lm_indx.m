function f = lm_indx
% syntax: f = lm_indx
% 5:02 PM 15/08/97
% determines the relationship between lmod and nc loads
% checks for lmod
% determines number of modulated loads
% f is a dummy variable
f = 0;
global lmod_con load_con  n_lmod  lmod_idx
n_lmod = 0;
lmod_idx = [];
if ~isempty(lmod_con)
    n_lmod = length(lmod_con(:,1));
    lmod_idx = zeros(n_lmod,1);
    for j = 1:n_lmod
       index = find(lmod_con(j,2)==load_con(:,1));
       if ~isempty(index)
          lmod_idx(j) = index;
       else
          error('you must have the load modulation bus declared as a non-conforming load')
       end
    end
end
       
    