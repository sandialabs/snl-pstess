function Y = small20(X,small)
%SMALL20 Y = small20(X,small)
%        Zeros out small numbers in an array.
%Inputs:  X - an array with some small elements that
%              are to be replaced by zero
%        small - max absolute value of numbers to be  
%              set to zero
%Output: Y - X array with small elements set to zero

%%%%%%%%%%%%%%%%%%%%% small20.m %%%%%%%%%%%%%%%%%%%%%
%       Feedback Control Problems with MATLAB 
%           and the Control System Toolbox
%      D. K. Frederick and J. H. Chow, Nov. 94
%----------------------------------------------------

Y = X;
vec = find(abs(X) <= abs(small) & X ~= 0);
Y(vec) = 0*vec;

%%%%%%%%%%%%%%%%%% end of small20.m %%%%%%%%%%%%%%%%%
