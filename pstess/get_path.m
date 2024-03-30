function [ data_file, path_name ] = get_path( )
% Syntax: [data_file,path_name] = get_path()
%
% Purpose: Returns the absolute path corresponding to the data folder and
%          current base case file. This is a simple function that accepts
%          no arguments. Its sole purpose is to store and return the
%          absolute path leading to the PSTess base case data files.
%
% Output:  data_file - the name of the current base case file
%          path_name - the full name of the directory where data_file resides

%-----------------------------------------------------------------------------%

path_name = pwd;
path_name = [path_name(1:end-6),'data'];

data_file = 'd2asbegp_ess.m';

end  % function end

% eof
