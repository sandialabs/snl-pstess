function set_path(data_file,path_name)
% Syntax: set_path(data_file) OR
%         set_path(data_file,path_name)
%
% Purpose: Changes the data file and path name returned by get_path().
%          This function works by rewriting the text file that defines
%          get_path().
%
% Input:   data_file - name of the desired base case file
%          path_name - full path of the directory where data_file resides
%
% Note:    If path_name is not provided, set_path() assumes you do not
%          wish to change it, and only the data_file will be updated.

%-----------------------------------------------------------------------------%
% Here is the format of get_path.m
%
% function [ data_file, path_name ] = get_path( )
% % Syntax: [data_file,path_name] = get_path()
% %
% % Purpose: Returns the absolute path corresponding to the data folder and
% %          current base case file. This is a simple function that accepts
% %          no arguments. Its sole purpose is to store and return the
% %          absolute path leading to the PSTess base case data files.
% %
% % Output:  data_file - the name of the current base case file
% %          path_name - the full name of the directory where data_file resides
%
% %-----------------------------------------------------------------------------%
%
% path_name = 'C:\Users\joebloggs\Documents\';
%
% data_file = 'd_minniWECC_ess.m';
%
% end  % function end
%
% % eof
%
%-----------------------------------------------------------------------------%

%-----------------------------------------------------------------------------%
% Version history
%
% Version: 1.0
% Author:  Ryan Elliott
% Note:    Initial version
%-----------------------------------------------------------------------------%

if (nargin < 2)
    [~, path_name] = get_path();
end

str_cell{1,:} = 'function [ data_file, path_name ] = get_path( )';
str_cell{2,:} = '% Syntax: [data_file,path_name] = get_path()';
str_cell{3,:} = '%';
str_cell{4,:} = '% Purpose: Returns the absolute path corresponding to the data folder and';
str_cell{5,:} = '%          current base case file. This is a simple function that accepts';
str_cell{6,:} = '%          no arguments. Its sole purpose is to store and return the';
str_cell{7,:} = '%          absolute path leading to the PSTess base case data files.';
str_cell{8,:} = '%';
str_cell{9,:} = '% Output:  data_file - the name of the current base case file';
str_cell{10,:} = '%          path_name - the full name of the directory where data_file resides';
str_cell{11,:} = '';
str_cell{12,:} = '%-----------------------------------------------------------------------------%';
str_cell{13,:} = '';
str_cell{14,:} = ['path_name = ''', path_name, ''';'];
str_cell{15,:} = '';
str_cell{16,:} = ['data_file = ''', data_file, ''';'];
str_cell{17,:} = '';
str_cell{18,:} = 'end  % function end';
str_cell{19,:} = '';
str_cell{20,:} = '% eof';

fid = fopen('get_path.m','w');
fprintf(fid,'%s\r\n',str_cell{:});
fclose(fid);

end  % function end

% eof
