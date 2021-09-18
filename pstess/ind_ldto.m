function ind_ldto(i,k)
% Syntax: ind_ldto(i,k)
%
% Purpose: Template for induction motor load torque calculation as
%          a function of slip
%
% Inpute:  i is the motor number
%          k is the time step
%          if i is set to zero, vector computation is invoked
%
% Data format for mld_con
%          1 motor number
%          2 bus number
%          3 stiction load pu on motor base (f1)
%          4 stiction load coefficient (i1)
%          5 external load  pu on motor base(f2)
%          6 external load coefficient (i2)
%
%          load has the form
%          tload = f1*slip^i1 + f2*(1-slip)^i2

%-----------------------------------------------------------------------------%
% Version history
%
% Author: Graham Rogers
% Date:   November 1995
%-----------------------------------------------------------------------------%

global g;  % declaring struct of global variables

if (i == 0)  % vector computation
    numot = length(g.ind.mld_con(:,2));
    g.ind.tload(:,k) = ...
        g.ind.mld_con(:,3).*(g.ind.slip(:,k).^g.ind.mld_con(:,4)) ...
        + g.ind.mld_con(:,5).*(ones(numot,1) - g.ind.slip(:,k)).^g.ind.mld_con(:,6);
else
    g.ind.tload(i,k) = ...
        g.ind.mld_con(i,3)*(g.ind.slip(i,k)^g.ind.mld_con(i,4)) ...
        + g.ind.mld_con(i,5)*(1 - g.ind.slip(i,k))^g.ind.mld_con(i,6);
end

end  % function end

% eof
