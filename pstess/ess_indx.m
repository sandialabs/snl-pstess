function ess_indx( )
% Syntax: ess_indx()
%
% Purpose: Determines the indexes in load_con corresponding to ess
%          locations, i.e., maps the relationship between ess and nc_load
%        - Checks for ess
%        - Checks for user-defined energy storage controls

%-----------------------------------------------------------------------------%
% Version history
%
% Version: 1.0
% Author:  Ryan T. Elliott
% Date:    August 2019
% Note:    Initial version
%-----------------------------------------------------------------------------%

global g;  % declaring struct of global variables

g.ess.n_ess = 0;
g.ess.ess_idx = [];

g.ess.n_essud = 0;
g.ess.dessud_idx = [];

if ~isempty(g.ess.ess_con)
    g.ess.n_ess = size(g.ess.ess_con,1);
    g.ess.ess_idx = zeros(g.ess.n_ess,1);

    for j = 1:g.ess.n_ess
        index = find(g.ess.ess_con(j,2) == g.ncl.load_con(:,1));
        if ~isempty(index)
            g.ess.ess_idx(j) = index;           % index relative to load_con
        else
            estr = '\ness_indx: the ess at bus %0.0f must be declared ';
            estr = [estr, 'as a non-conforming load in load_con.'];
            error(sprintf(estr,g.ess.ess_con(j,2)));
        end
    end

    % check for user-defined controls
    if ~isempty(g.ess.essud_con)
        g.ess.n_essud = size(g.ess.essud_con,1);
        for j = 1:g.ess.n_essud
            ess_idx = find(g.ess.essud_con(j,2) == g.ess.ess_con(:,2));
            if ~isempty(ess_idx)
                g.ess.dessud_idx(j) = ess_idx;  % index relative to ess_con
            else
                estr = '\ness_indx: the ess_sud at bus %0.0f must be connected to ';
                estr = [estr, 'an ess model.'];
                error(sprintf(estr,g.ess.essud_con(j,2)));
            end
        end
    end
end  % if-statement end

end  % function end

% eof
