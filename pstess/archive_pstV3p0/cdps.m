function p = cdps(dir)
if nargin==0
    p = 'C:\Program Files\MATLAB\R2006b\PSTV2\';
else
    p = ['C:\Program Files\MATLAB\R2006b\PSTV2\' dir];
end
cd(p)