function pwrmod_indx(bus)
% Syntax: pwrmod_indx
%
% Purpose: Determines the relationship between pwrmod and nc loads
%        - Checks for pwrmod
%        - Determines number of modulated buses

%-----------------------------------------------------------------------------%
% Version history
%
% Version: 1.0 (initial version)
% Date:    2015
% Author:  D. Trudnowski
%-----------------------------------------------------------------------------%

global g;  % declaring struct of global variables

g.pwr.pwrmod_idx = [];

if ~isempty(g.pwr.pwrmod_con)
    if isempty(g.ncl.load_con)
        estr = 'pwrmod_indx: pwrmod buses must be declared as ';
        estr = [estr, 'non-conforming loads in load_con.'];
        error(estr);
    end

    g.pwr.n_pwrmod = size(g.pwr.pwrmod_con,1);
    g.pwr.pwrmod_idx = zeros(g.pwr.n_pwrmod,1);
    for j = 1:g.pwr.n_pwrmod
        index = find(g.pwr.pwrmod_con(j,1)==g.ncl.load_con(:,1));
        if ~isempty(index)
            if (abs(sum(g.ncl.load_con(index,2:end)')-2) > 1e-8)
                estr = '\npwrmod_indx: the pwrmod device at index %0.0f must be ';
                estr = [estr, 'declared as either 100%% constant power or '];
                estr = [estr, '100%% constant current.'];
                error(sprintf(estr,j));
            end

            if ((g.ncl.load_con(index,2) == 1 && g.ncl.load_con(index,3) == 1) ...
                || (g.ncl.load_con(index,4) == 1 && g.ncl.load_con(index,5) == 1));
                g.pwr.pwrmod_idx(j) = index;
            else
                estr = '\npwrmod_indx: the pwrmod device at index %0.0f must be ';
                estr = [estr, 'declared as either 100%% constant power or '];
                estr = [estr, '100%% constant current.'];
                error(sprintf(estr,j));
            end

            kk = g.bus.bus_int(g.pwr.pwrmod_con(:,1));
            if any(abs(bus(kk,10)-2))
                estr = '\npwrmod_indx: the pwrmod device at index %0.0f must be ';
                estr = [estr, 'declared as a type 2 (PV) bus in power flow (bus).'];
                error(sprintf(estr,j));
            end

            if ((max(abs(bus(kk,6))) > 1e-10) || (max(abs(bus(kk,7))) > 1e-10))
                estr = '\npwrmod_indx: the pwrmod device at index %0.0f must have ';
                estr = [estr, 'zero load upon initialization.'];
                error(sprintf(estr,j));
            end
        else
            estr = '\npwrmod_indx: the pwrmod device at index %0.0f must be ';
            estr = [estr, 'declared as a non-conforming load in load_con.'];
            error(sprintf(estr,j));
        end
    end
end

end  % function end

% eof
