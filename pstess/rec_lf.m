function [alpha,mode,idc_new] = rec_lf(idc,al_min,al_max,Vdo,Rc,cm)
% Syntax: [alpha,mode,idc_new] = rec_lf(idc,al_min,al_max,Vdo,RC,cm)
%
% Purpose: Finds alpha for rectifier and sets rectifier taps
%          determines if mode change is necessary
%          sets new idc
%
% Note:    When mode(i) = 2, idc(i) should be reduced by the current margin
%
% Inputs:  idc specified dc line current vector
%          al_min vector of minimum alpha
%          al_max vector of maximum alpha
%          Vdo vector of ideal rectifier dc voltages
%          Rc vector of rectifier commutating resistances
%          cm vector of current margins
%
% Output:  alpha, the rectifier firing angle in radians
%          mode - a vector giving the mode of operation of the inverter
%          mode(i) = 1 indicates that the ith inverter is operating at gamma min
%          mode(i) = 2 indicates that the ith inverter is controlling current
%          idc_new is set to specified idc when mode is 1
%          or to idc reduced by the current margin if mode is 2

%-----------------------------------------------------------------------------%
% Version history
%
% Version: 1.0 (initial version)
% Author:  Graham Rogers
% Date:    November 1996
%-----------------------------------------------------------------------------%

global g;  % declaring struct of global variables

% note: the output alpha is in degrees so it is not the same as g.dc.alpha,
%       which is in radians

% calculate alpha
Vdcr = g.dc.Vdc(g.dc.r_idx);
mode = ones(g.dc.n_dcl,1);     % set to default initially
calpha = (Vdcr+idc.*Rc)./Vdo;

% check that alpha is within range
max_idx = find(calpha < cos(al_max));
min_idx = find(calpha > cos(al_min));
oor_idx = find((calpha < cos(al_max)) | (calpha > cos(al_min)));

if ~isempty(oor_idx)
    Vdo_new = zeros(g.dc.n_dcl,1);
    tapn = g.dc.tapr;
    tnum = (tapn - g.dc.tminr)./g.dc.tstepr;

    % adjust rectifier taps to bring within range
    if ~isempty(max_idx)
        Vdo_new(max_idx) = (Vdcr(max_idx) + Rc(max_idx).*idc(max_idx)) ...
                           ./cos(al_max(max_idx));
    end

    if ~isempty(min_idx)
        Vdo_new(min_idx) = (Vdcr(min_idx) + Rc(min_idx).*idc(min_idx)) ...
                           ./cos(al_min(min_idx));
    end

    tapn(oor_idx) = Vdo(oor_idx)./Vdo_new(oor_idx);

    % get the right tap setting
    tnum = (tapn - g.dc.tminr)./g.dc.tstepr;
    if ~isempty(max_idx)
        tnum(max_idx) = ceil(tnum(max_idx));  % next interger up
    end
    %
    if ~isempty(min_idx)
        tnum(min_idx) = fix(tnum(min_idx));   % next integer down
    end

    tmin_idx = find(tnum < 0);
    tmax_idx = find(tnum > (g.dc.tmaxr-g.dc.tminr)./g.dc.tstepr);

    if ~isempty(tmin_idx)
        if ~isempty(max_idx)
            for k = 1:length(max_idx)
                ktm = find(max_idx(k) == tmin_idx)
                if ~isempty(ktm)
                    rec_num = g.dc.r_idx(max_idx(ktm));
                    estr = '\nrec_lf: minimum tap setting reached at ';
                    estr = [estr, 'converter index %0.0f.'];
                    error(sprintf(estr,rec_num));
                end
            end
        end
        %
        if ~isempty(min_idx)
            if ~isempty(min_idx(tmin_idx))
                rec_num = g.dc.r_idx(min_idx(tmin_idx));
                wstr = '\nrec_lf: minimum tap setting reached at ';
                wstr = [wstr, 'converter index %0.0f, changing mode.'];
                warning(sprintf(wstr,rec_num));
                %
                n_mc = length(min_idx(tmin_idx));
                mode(g.dc.i_idx(min_idx(tmin_idx))) = 2*ones(n_mc,1);
                tnum(tmin_idx) = zeros(n_mc,1);
            end
        end
    end
    %
    if ~isempty(tmax_idx)
        if ~isempty(max_idx)
            if ~isempty(max_idx(tmax_idx))
                rec_num = g.dc.r_idx(max_idx(tmax_idx));
                estr = '\nrec_lf: maximum tap setting reached at ';
                estr = [estr, 'converter index %0.0f.'];
                error(sprintf(estr,rec_num));
            end
        end
        %
        if ~isempty(min_idx)
            if ~isempty(min_idx(tmax_idx))
                rec_num = g.dc.r_idx(min_idx(tmax_idx));
                wstr = '\nrec_lf: maximum tap setting reached at ';
                wstr = [wstr, 'converter index %0.0f, changing mode.'];
                warning(sprintf(wstr,rec_num));
                %
                n_mc = length(min_idx(tmax_idx));
                mode(min_idx(tmax_idx)) = 2*ones(n_mc,1);
                tnum(tmax_idx) = (g.dc.tmaxr(oor_idx(tmax_idx)) ...
                                  - g.dc.tminr(oor_idx(tmax_idx))) ...
                                 ./g.dc.tstepr(oor_idx(tmax_idx));
            end
        end
    end

    tapre = g.dc.tminr + g.dc.tstepr.*tnum;

    % recalculate calpha
    Vdo = Vdo.*g.dc.tapr./tapre;
    g.dc.tapr = tapre;
    calpha = (g.dc.Vdc(g.dc.r_idx)+idc.*Rc)./Vdo;

    % apply alpha limits
    calpha = min(calpha,cos(al_min));
    calpha = max(calpha,cos(al_max));
end

% adjust idc for mode change
idc_new = idc;
mc_idx = find(mode == 2);
if ~isempty(mc_idx)
    idc_new(mc_idx) = idc(mc_idx).*cm(mc_idx);
end

salpha = sqrt(ones(g.dc.n_dcl,1)-calpha.*calpha);
alpha = atan2(salpha,calpha)*180/pi;

end  % function end

% eof
