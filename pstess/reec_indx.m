function reec_indx( )
% Syntax: reec_indx()
%
% Purpose: Determines the indexes in ess_con corresponding to reec
%          controllers, i.e., maps the relationship between reec and ess

%-----------------------------------------------------------------------------%
% Version history
%
% Version: 1.0
% Author:  Ryan T. Elliott
% Date:    November 2022
% Note:    Initial version
%-----------------------------------------------------------------------------%

global g;  % declaring struct of global variables

g.reec.ess_idx = [];

if ~isempty(g.reec.reec_con)
    g.reec.n_reec = size(g.reec.reec_con,1);
    g.reec.ess_idx = zeros(g.reec.n_reec,1);

    for j = 1:g.reec.n_reec
        index = find(g.reec.reec_con(j,2) == g.ess.ess_con(:,2));
        if ~isempty(index)
            g.reec.ess_idx(j) = index;           % index relative to ess_con
        else
            estr = '\nreec_indx: the reec at bus %0.0f must be connected ';
            estr = [estr, 'to an ess listed in ess_con.'];
            error(sprintf(estr,g.reec.reec_con(j,2)));
        end
    end
end

end  % function end

% eof
