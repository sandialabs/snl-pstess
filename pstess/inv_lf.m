function [gamma] = inv_lf(mode,idc,Vdc_max,Vdc_min,ga_min,ga_max,Vdo,Rc)
% Syntax: [gamma] = inv_lf(mode,idc,Vdc_max,Vdc_min,ga_min,ga_max,Vdo,RC)
%
% Purpose: Finds gamma for inverter and sets inverter taps
%
% Input:   mode - a vector giving the mode of operation of the inverter
%          mode(i) = 1 indicates that the ith inverter is operating at gamma min
%          mode(i) = 2 indicates that the ith inverter is controlling current
%          Note: when mode(i) = 2, idc(i) should be reduced by the current margin
%          idc specified dc line current vector
%          Vdc_max vector of maximum allowed inverter dc voltages
%          Vdc_min vector of minimum allowed inverter dc voltages
%          ga_min vector of minimum gamma
%          ga_max vector of maximum gamma
%          Vdo vector of ideal inverter dc voltages
%          Rc vector of inverter commutating resistances
%
% Output:  gamma, the inverter extinction angle in radians

%-----------------------------------------------------------------------------%
% Version history
%
% Author:  Graham Rogers
% Date:    November 1996
%-----------------------------------------------------------------------------%

global g;  % declaring struct of global variables

% note: the output gamma is in degrees so it is not the same as g.dc.gamma,
%       which is in radians

Vdo_new = Vdo;
mode_idx = find(mode == 1);
if ~isempty(mode_idx)
    % inverter at constant ga_min
    % calculate inverter dc voltage
    ga1 = ga_min(mode_idx);
    Rc1 = Rc(mode_idx);
    idc1 = idc(mode_idx);
    Vdo1 = Vdo(mode_idx);
    Vdc1 = Vdo1.*cos(ga1) - Rc1.*idc1;
    g.dc.Vdc(g.dc.i_idx) = Vdc1;

    % check that Vdc is within range
    max_idx = find(Vdc1 > Vdc_max(mode_idx));
    min_idx = find(Vdc1 < Vdc_min(mode_idx));
    oor_idx = find(Vdc1 > Vdc_max(mode_idx) | Vdc1 < Vdc_min(mode_idx));
    if ~isempty(max_idx)
        for k = 1:length(max_idx)
            mo_idx(k) = find(oor_idx == max_idx(k));
        end
    end

    if ~isempty(min_idx)
        for k = 1:length(min_idx)
            mno_idx(k) = find(oor_idx == min_idx(k));
        end
    end

    % reset taps if voltage out of limits
    t_idx = mode_idx(oor_idx);
    if ~isempty(oor_idx)
        % reset inverter taps to bring dc voltage back to limit
        if ~isempty(max_idx)
            Vdo_new(max_idx,1) = (Vdc_max(mode_idx(max_idx)) ...
                                  + Rc1(max_idx).*idc1(max_idx))./cos(ga1(max_idx));
        end

        if ~isempty(min_idx)
            Vdo_new(min_idx,1) = (Vdc_min(mode_idx(min_idx)) ...
                                  + Rc1(min_idx).*idc1(min_idx))./cos(ga1(min_idx));
        end

        % useful debugging variables
        % Vdo1(oor_idx)
        % Vdo_new(oor_idx)
        % g.dc.tapi(t_idx)

        tapn = g.dc.tapi(t_idx).*Vdo1(oor_idx)./Vdo_new(oor_idx);

        % get the right tap setting
        tnum = (tapn - g.dc.tmini(t_idx))./g.dc.tstepi(t_idx);

        % set tap step to the nearest integer higher for the maximum limit
        if ~isempty(max_idx)
            tnum(mo_idx) = ceil(tnum(mo_idx));
        end

        % set tap step to the nearest integer lower for the minimum limit
        if ~isempty(min_idx)
            tnum(mno_idx) = fix(tnum(mno_idx));
        end

        tapin = g.dc.tmini(t_idx) + g.dc.tstepi(t_idx).*tnum;
        tmax_idx = find(tapin > g.dc.tmaxi(t_idx));
        if ~isempty(tmax_idx)
            inv_num = g.dc.i_idx(t_idx(tmax_idx));
            estr = '\ninv_lf: maximum tap setting reached at ';
            estr = [estr, 'converter index %0.0f.'];
            error(sprintf(estr,inv_num));
        end

        tmin_idx = find(tnum == 0);
        if ~isempty(tmin_idx)
            inv_num = g.dc.i_idx(t_idx(tmin_idx));
            estr = '\ninv_lf: minimum tap setting reached at ';
            estr = [estr, 'converter index %0.0f.'];
            error(sprintf(estr,inv_num));
        end

        % recalculate Vdo
        Vdo1(t_idx) = Vdo1(t_idx).*g.dc.tapi(t_idx)./tapin;
        g.dc.tapi(t_idx) = tapin;
        Vdc1 = Vdo1.*cos(ga1)-Rc1.*idc1;
        g.dc.Vdc(g.dc.i_idx) = Vdc1;
    end

    % set gamma at gamma min
    gamma(mode_idx,1) = ga1*180/pi;
end

mode_idx = find(mode == 2);
if ~isempty(mode_idx)
    % inverter controlling current
    % idc must be set in the calling program to the lower current in this mode
    % dc voltage is specified at the rectifier
    % determine gamma and inverter taps required to maintain idc
    Vdc2 = g.dc.Vdc(g.dc.i_idx(mode_idx));  % inverter dc voltage
    Rc2 = Rc(mode_idx);
    idc2 = idc(mode_idx);
    Vdo2 = Vdo(mode_idx);
    gamn2 = ga_min(mode_idx);
    gamx2 = ga_max(mode_idx);
    cgamma = (Vdc2 + Rc2.*idc2)./Vdo2;

    % find gamma out of range
    min_idx = find(cgamma > cos(gamn2));
    max_idx = find(cgamma < cos(gamx2));
    oor_idx = find(cgamma > cos(gamn2) | cgamma < cos(gamx2));
    t_idx = mode_idx(oor_idx);
    if ~isempty(max_idx)
        for k = 1:length(max_idx)
            mo_idx(k) = find(oor_idx == max_idx(k));
        end
    end

    if ~isempty(min_idx)
        for k = 1:length(min_idx)
            mno_idx(k) = find(oor_idx == min_idx(k));
        end
    end

    if ~isempty(oor_idx)
        if ~isempty(min_idx)
            cgamma(min_idx) = cos(gamn2(min_idx));
        end

        if ~isempty(max_idx)
            cgamma(max_idx) = cos(gamx2(max_idx));
        end

        % calculate the required Vdo
        Vdo_new = (Vdc2(oor_idx) + Rc2(oor_idx).*idc2(oor_idx))./cgamma(oor_idx);

        % calculate ideal tap setting
        tapn = g.dc.tapi(t_idx).*Vdo2(oor_idx)./Vdo_new;

        % get the right tap setting
        tnum = (tapn - g.dc.tmini(t_idx))./g.dc.tstepi(t_idx);
        if ~isempty(max_idx)  % set to next whole number up
            tnum(mo_idx) = ceil(tnum(mo_idx));
        end

        if ~isempty(min_idx)  % set to next whole number down
            tnum(mno_idx) = fix(tnum(mno_idx));
        end

        tapin = g.dc.tmini(t_idx) + g.dc.tstepi(t_idx).*tnum;
        tmax_idx = find(tapin > g.dc.tmaxi(t_idx));
        if ~isempty(tmax_idx)
            inv_num = g.dc.i_idx(t_idx(tmax_idx));
            estr = '\ninv_lf: maximum tap setting reached at ';
            estr = [estr, 'converter index %0.0f.'];
            error(sprintf(estr,inv_num));
        end

        tmin_idx = find(tnum == 0);
        if ~isempty(tmin_idx)
            inv_num = g.dc.i_idx(t_idx(tmin_idx));
            estr = '\ninv_lf: minimum tap setting reached at ';
            estr = [estr, 'converter index %0.0f.'];
            error(sprintf(estr,inv_num));
        end

        % recalculate Vdo
        Vdo2 = Vdo2.*g.dc.tapi(t_idx)./tapin;
        g.dc.tapi(t_idx) = tapin;
        Vdc2 = Vdo2.*cgamma-Rc2.*idc2;
        g.dc.Vdc(g.dc.i_idx) = Vdc2;

        sgamma = sqrt(ones(length(mode_idx),1)-cgamma.*cgamma);
        gamma = atan2(sgamma,cgamma)*180/pi;
    end
end

end  % function end

% eof
