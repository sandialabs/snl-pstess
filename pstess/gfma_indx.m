function gfma_indx( )
% Syntax: gfma_indx()
%
% Purpose: Determines the indexes in ivm_con and mac_con corresponding
%          to gfma controllers, i.e., maps the relationship between
%          gfma and ivm models

%-----------------------------------------------------------------------------%
% Version history
%
% Version: 1.0
% Author:  Ryan T. Elliott
% Date:    October 2022
% Note:    Initial version
%-----------------------------------------------------------------------------%

global g;  % declaring struct of global variables

g.gfma.ivm_idx = [];
g.gfma.mac_idx = [];

if ~isempty(g.gfma.gfma_con)
    g.gfma.n_gfma = size(g.gfma.gfma_con,1);

    for j = 1:g.gfma.n_gfma
        ivm_idx = find(g.gfma.gfma_con(j,2) == g.ivm.ivm_con(:,2));
        if ~isempty(ivm_idx)
            g.gfma.ivm_idx(j) = ivm_idx;  % index relative to ivm_con
            if (g.ivm.ivm_con(ivm_idx,3) ~= 1)
                estr = '\ngfma_indx: the ivm controller type is incorrectly ';
                estr = [estr, 'specified at bus %0.0f.'];
                error(sprintf(estr,g.gfma.gfma_con(j,2)));
            elseif (g.ivm.ivm_con(ivm_idx,7) > 0)
                estr = '\ngfma_indx: the ivm commanded angle time constant ';
                estr = [estr, 'must be set to zero at bus %0.0f.'];
                error(sprintf(estr,g.gfma.gfma_con(j,2)));
            elseif (g.ivm.ivm_con(ivm_idx,8) > 0)
                estr = '\ngfma_indx: the ivm commanded voltage magnitude time ';
                estr = [estr, 'constant must be set to zero at bus %0.0f.'];
                error(sprintf(estr,g.gfma.gfma_con(j,2)));
            end
        else
            estr = '\ngfma_indx: the gfma controller at bus %0.0f must be ';
            estr = [estr, 'connected to an ivm model.'];
            error(sprintf(estr,g.gfma.gfma_con(j,2)));
        end

        mac_idx = find(g.gfma.gfma_con(j,2) == g.mac.mac_con(:,2));
        if ~isempty(mac_idx)
            g.gfma.mac_idx(j) = mac_idx;  % index relative to mac_con
        else
            estr = '\ngfma_indx: the gfma controller at bus %0.0f must ';
            estr = [estr, 'correspond to a machine model.'];
            error(sprintf(estr,g.gfma.gfma_con(j,2)));
        end
    end
end

end  % function end

% eof
