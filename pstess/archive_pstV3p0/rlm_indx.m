function f = rlm_indx
% syntax: f = rlm_indx
% 5:28 PM 27/8/97
% determines the relationship between lmod and nc loads
% checks for rlmod
% determines number of modulated loads
% f is a dummy variable
f = 0;
global rlmod_con load_con  n_rlmod  rlmod_idx
n_rlmod = 0;
rlmod_idx = [];
if ~isempty(rlmod_con)
    n_rlmod = length(rlmod_con(:,1));
    rlmod_idx = zeros(n_rlmod,1);
    for j = 1:n_rlmod
       index = find(rlmod_con(j,2)==load_con(:,1));
       if ~isempty(index)
          rlmod_idx(j) = index;
       else
          error('you must have the reactive load modulation bus declared as a non-conforming load')
       end
    end
end
       
    