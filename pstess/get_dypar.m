function [ dypar ] = get_dypar( )
% Syntax: dypar = get_dypar()
%
% Purpose: Returns a struct of PSTess parameter settings. This is a
%          functional global variable that controls the behavior of
%          PSTess. These variables are used in many different places
%          throughout the code, so this function allows you to set
%          them in one place rather than in each file.
%
% Outputs: dypar - a struct of PSTess application parameter settings

%-----------------------------------------------------------------------------%
% Version history
%
% Version: 1.0
% Author:  Ryan T. Elliott
% Date:    July 2020
% Note:    Initial version
%-----------------------------------------------------------------------------%

dypar.batch_mode = true;  % may be either true or false
dypar.basmva = 100;       % system MVA base for batch mode
dypar.sys_freq = 60;      % nominal frequency (Hz) for batch mode

end  % function end

% eof
