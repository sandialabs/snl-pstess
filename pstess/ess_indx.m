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

g.ess.ess_idx = [];
g.ess.dessud_idx = [];

if ~isempty(g.ess.ess_con)
    [g.ess.n_ess,n_par] = size(g.ess.ess_con);
    g.ess.ess_idx = zeros(g.ess.n_ess,1);

    if (n_par < 22)                             % check for ke parameters
        g.ess.ess_con(:,22) = ones(g.ess.n_ess,1);
    end

    vdl_default.vp = [-1,2];                    % default vdl curves (from PSLF)
    vdl_default.ip = [1.10,1.10];
    vdl_default.vq = [-1,2];
    vdl_default.iq = [1.45,1.45];

    if ~isfield(g.ess,'vdl')
        g.ess.vdl = cell(g.ess.n_ess,1);
    else
        if (numel(g.ess.vdl) ~= g.ess.n_ess)
            error('ess_index: incorrect number of VDL curves specified.');
        end
    end

    for j = 1:g.ess.n_ess
        if isempty(g.ess.vdl{j})
            g.ess.vdl{j} = vdl_default;
        else
            vdl = g.ess.vdl{j};
            if ~(isfield(vdl,'vp') & isfield(vdl,'ip') ...
                 & isfield(vdl,'vq') & isfield(vdl,'iq'))
                estr = '\ness_indx: VDL curves missing required fields ';
                estr = [estr, 'for the ess at bus %0.0f.'];
                error(sprintf(estr,g.ess.ess_con(j,2)));
            end

            if ~(length(vdl.vp) == length(vdl.ip) ...
                 & length(vdl.vq) == length(vdl.iq))
                estr = '\ness_indx: VDL curve dimensions do not match ';
                estr = [estr, 'for the ess at bus %0.0f.'];
                error(sprintf(estr,g.ess.ess_con(j,2)));
            end

            if any(diff(vdl.vp) <= 0) | any(diff(vdl.vq) <= 0)
                estr = '\ness_indx: VDL curve voltages must be in ascending ';
                estr = [estr, 'order for the ess at bus %0.0f.'];
                error(sprintf(estr,g.ess.ess_con(j,2)));
            end

            if any(vdl.ip < 0)
                estr = '\ness_indx: active current VDL curve contains negative ';
                estr = [estr, 'value(s) at bus %0.0f.'];
                error(sprintf(estr,g.ess.ess_con(j,2)));
            end
        end
    end

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
