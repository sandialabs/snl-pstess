function dc_indx(line,dci_dc,dcr_dc)
% Syntax: dc_indx(line,dci_dc,dcr_dc)
%
% Purpose: To form indexes for the rectifier and inverter in the dc load flow
%          and indicate the ac buses contected to the converters
%
% Input:   ac line matrix
%          dcr_dc - user defined damping control at rectifier cell
%          dci_dc - user defined damping control at inverter cell

%-----------------------------------------------------------------------------%
% Version history
%
% Date:    January 1999
% Author:  Graham Rogers
% Purpose: Addition of user defined damping controls
%
% Author:  Graham Rogers
% Date:    October 1996
%-----------------------------------------------------------------------------%

global g;  % declaring struct of global variables

% pick out ac voltages (should be LT converter transformer buses)
% check that line and cont data is consistent

g.dc.dcrud_idx = [];
g.dc.dciud_idx = [];
g.dc.ndcr_ud = 0;
g.dc.ndci_ud = 0;

if ~isempty(g.dc.dcsp_con)
    lconv = size(g.dc.dcsp_con,1);
    lline = size(g.dc.dcl_con,1);
    lcon = size(g.dc.dcc_con,1);
    if (lcon ~= 2*lline || lcon ~= lconv)
        estr = '\ndc_indx: dc converter and line data are inconsistent.';
        estr = [estr, '\n  lconv: number of converters = %0.0f'];
        estr = [estr, '\n  lline: number of lines = %0.0f'];
        estr = [estr, '\n   lcon: number of controls = %0.0f'];
        error(sprintf(estr,[lconv,lline,lcon]));
    end

    % r_idx -- index of rectifier buses
    % i_idx -- index of inverter buses
    g.dc.r_idx = find(g.dc.dcsp_con(:,3) == 1);
    g.dc.i_idx = find(g.dc.dcsp_con(:,3) == 2);

    % ric_idx -- index of rectifier current control
    % rpc_idx -- index of rectifier power control
    g.dc.ric_idx = find(g.dc.dcc_con(:,9) == 1);
    g.dc.rpc_idx = find(g.dc.dcc_con(:,9) == 2);

    g.dc.n_dcl = lline;
    g.dc.n_conv = lconv;

    g.dc.inv_ac_bus = g.bus.bus_int(g.dc.dcsp_con(g.dc.i_idx,2));
    g.dc.rec_ac_bus = g.bus.bus_int(g.dc.dcsp_con(g.dc.r_idx,2));
    g.dc.ac_bus = g.bus.bus_int(g.dc.dcsp_con(:,2));

    g.dc.inv_ac_line = zeros(lline,1);
    g.dc.rec_ac_line = zeros(lline,1);

    for j = 1:lline
        acilj = find(g.bus.bus_int(line(:,2)) == g.dc.inv_ac_bus(j));
        if isempty(acilj)
            estr = '\ndc_indx: the inverter bus is not declared as a ''to bus'' ';
            estr = [estr, 'at dc line index %0.0f.'];
            error(sprintf(estr,j));
        else
            g.dc.inv_ac_line(j) = acilj;
            acilj = [];
        end

        acrlj = find(g.bus.bus_int(line(:,2)) == g.dc.rec_ac_bus(j));
        if isempty(acrlj)
            estr = '\ndc_indx: the rectifier bus is not declared as a ''to bus'' ';
            estr = [estr, 'at dc line index %0.0f.'];
            error(sprintf(estr,j));
        else
            g.dc.rec_ac_line(j) = acrlj;
            acrlj = [];
        end
    end

    g.dc.ac_line = [g.dc.rec_ac_line; g.dc.inv_ac_line];

    % form index of dc lines associated with the inverters
    g.dc.dcli_idx = zeros(g.dc.n_dcl,1);
    for k = 1:g.dc.n_dcl
        g.dc.dcli_idx = (g.dc.dcli_idx ...
                         | (g.dc.dcl_con(k,2) == g.dc.dcsp_con(g.dc.i_idx,1)));
    end

    g.dc.dcli_idx = find(g.dc.dcli_idx ~= 0);
end

g.dc.no_cap_idx = find(g.dc.dcl_con(:,5) == 0);
g.dc.cap_idx = find(g.dc.dcl_con(:,5) ~= 0);
g.dc.l_no_cap = 0;

if ~isempty(g.dc.no_cap_idx);
    g.dc.l_no_cap = length(g.dc.no_cap_idx);
end

g.dc.l_cap = g.dc.n_dcl - g.dc.l_no_cap;
g.dc.no_ind_idx = find(g.dc.dcl_con(:,4) == 0 ...
                       | g.dc.dcl_con(:,6) == 0 | g.dc.dcl_con(:,7) == 0);

% index of converters in load_con
j = g.bus.bus_int(g.ncl.load_con(:,1));
for k = 1:g.dc.n_conv
    g.dc.ldc_idx(k) = find(j == g.dc.ac_bus(k));
    if isempty(g.dc.ldc_idx(k))
        estr = 'dc_indx: dc converter LT buses must be declared in load_con, ';
        estr = [estr, 'check converter index %0.0f.'];
        error(sprintf(estr,k));
        % an additional load is not allowed at a converter bus
    end
end
% j(g.dc.ldc_idx(k)) gives the internal bus number of the kth converter

% check for user defined controls
[g.dc.ndcr_ud,~] = size(dcr_dc);
for j = 1:g.dc.ndcr_ud
    if ~isempty(dcr_dc{j})
        g.dc.dcrud_idx(j) = dcr_dc{j,2};
    end
    g.dc.dcrud_idx(j) = 0;
end

[g.dc.ndci_ud,~] = size(dci_dc);
for j = 1:g.dc.ndci_ud
    if ~isempty(dci_dc{j})
        g.dc.dciud_idx(j) = dci_dc{j,2};
    end
    g.dc.dciud_idx(j) = 0;
end

end  % function end

% eof
