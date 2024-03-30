function ivm_indx(varargin)
% Syntax: ivm_indx()
%
% Purpose: Determines the indexes in mac_con corresponding to ivm
%          generators, i.e., maps the relationship between ivm and mac_con
%        - Checks for ivm
%        - Checks for user-defined ivm controls

%-----------------------------------------------------------------------------%
% Version history
%
% Version: 1.0
% Author:  Ryan T. Elliott
% Date:    October 2022
% Note:    Initial version
%-----------------------------------------------------------------------------%

global g;  % declaring struct of global variables

import_flag = false;
if (nargin > 0)
    import_flag = varargin{1};
end

g.ivm.divmud_idx = [];
g.ivm.divmud_mac_idx = [];

if (~isempty(g.ivm.ivm_con) && import_flag)
    % importing ivm_con into mac_con

    g.mac.n_ivm = size(g.ivm.ivm_con,1);
    g.mac.mac_ivm_idx = ...
        (g.mac.n_not_ivm+1:g.mac.n_not_ivm+g.mac.n_ivm);

    if ~isempty(g.mac.mac_con)
        ivm_mac_con = zeros(g.mac.n_ivm,size(g.mac.mac_con,2));
        mac_num_max = max(g.mac.mac_con(:,1));
    else
        ivm_mac_con = zeros(g.mac.n_ivm,21);
        mac_num_max = 0;
    end

    for ii = 1:size(ivm_mac_con,1)
        ivm_mac_con(ii,1) = mac_num_max + ii;      % machine number
        ivm_mac_con(ii,2) = g.ivm.ivm_con(ii,2);   % bus number
        ivm_mac_con(ii,3) = g.ivm.ivm_con(ii,4);   % mva base
        ivm_mac_con(ii,5) = g.ivm.ivm_con(ii,5);   % r_a, source resistance
        ivm_mac_con(ii,7) = g.ivm.ivm_con(ii,6);   % x'_d, source reactance
        ivm_mac_con(ii,9) = g.ivm.ivm_con(ii,7);   % Td, ang time constant
        ivm_mac_con(ii,10) = g.ivm.ivm_con(ii,8);  % Tv, mag time constant
        ivm_mac_con(ii,16) = g.ivm.ivm_con(ii,9);  % H, inertia time constant
        ivm_mac_con(ii,19) = g.ivm.ivm_con(ii,2);  % bus number
    end

    g.mac.mac_con = [g.mac.mac_con; ivm_mac_con];
elseif (~isempty(g.ivm.ivm_con) && ~import_flag)
    g.ivm.n_ivm = size(g.ivm.ivm_con,1);
    g.mac.n_ib_ivm = 0;                            % ivm infinite buses

    if (g.mac.n_ib ~= 0)
        ib_ivm = ismember(g.mac.mac_ib_idx,g.mac.mac_ivm_idx);

        if ~isempty(ib_ivm)
            g.mac.mac_ib_ivm = g.mac.mac_ib_idx(ib_ivm);
            g.mac.n_ib_ivm = length(g.mac.mac_ib_ivm);
        end
    end

    for j = 1:g.ivm.n_ivm
        if (g.ivm.ivm_con(j,6) < 1e-12)
            estr = '\nivm_indx: the ivm at bus %0.0f must have a nonzero ';
            estr = [estr, 'coupling reactance.'];
            error(sprintf(estr,g.ivm.ivm_con(j,2)));
        end
    end

    % check for user-defined controls
    if ~isempty(g.ivm.ivmud_con)
        g.ivm.n_ivmud = size(g.ivm.ivmud_con,1);
        for j = 1:g.ivm.n_ivmud
            ivm_idx = find(g.ivm.ivmud_con(j,2) == g.ivm.ivm_con(:,2));
            if ~isempty(ivm_idx)
                g.ivm.divmud_idx(j) = ivm_idx;      % index relative to ivm_con
                if (g.ivm.ivm_con(ivm_idx,3) ~= 0)
                    estr = '\nivm_indx: the ivm controller type is incorrectly ';
                    estr = [estr, 'specified at bus %0.0f.'];
                    error(sprintf(estr,g.ivm.ivmud_con(j,2)));
                end
            else
                estr = '\nivm_indx: the ivm_sud at bus %0.0f must be connected to ';
                estr = [estr, 'an ivm model.'];
                error(sprintf(estr,g.ivm.ivmud_con(j,2)));
            end

            mac_idx = find(g.ivm.ivmud_con(j,2) == g.mac.mac_con(:,2));
            if ~isempty(mac_idx)
                g.ivm.divmud_mac_idx(j) = mac_idx;  % index relative to mac_con
            else
                estr = '\nivm_indx: the ivm_sud at bus %0.0f must correspond to ';
                estr = [estr, 'a machine model.'];
                error(sprintf(estr,g.ivm.ivmud_con(j,2)));
            end
        end
    end
end  % if-statement end

end  % function end

% eof
