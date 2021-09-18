function p = cdps(dir)
% Syntax: p = cdps(dir)
%
% Purpose: changes the current directory in a convenient way for use
%          with PSTess. This function works best if it's in a folder
%          included in your MATLAB path.
%
% Input:   dir - the directory under the PSTess parent folder that you
%                want to switch to
%
% Output:  the full name of the current directory

%-----------------------------------------------------------------------------%
% Version history
%
% Author:  Ryan Elliott
% Date:    July 2020
% Purpose: Updated to accommodate batch processing mode
%
% Version: 1.0
% Author:  Graham Rogers and/or Joe Chow
% Note:    Initial version
%-----------------------------------------------------------------------------%

dypar = get_dypar();                % get parameter settings
batch_mode = dypar.batch_mode;      % batch processing mode

% get base path information
if batch_mode
    [~,pathname] = get_path();
    slash_char = strfind(pathname,'\');
    pathname = pathname(1:slash_char(end));
else
    pathname = uigetdir('C:\Users\','Select base directory');
end

% specify new directory
if (nargin == 0)
    p = pathname;                   % no argument given
else
    p = [pathname, '\', dir];       % directory from input
end

cd(p);

fprintf('\nChanging current directory to:\n\n%s\n\n', p);

end  % function end

% eof
