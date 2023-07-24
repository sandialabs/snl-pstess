function dc_cont(i,k,kdc,bus,flag)
% Syntax: dc_cont(i,k,kdc,bus,flag)
%
% Purpose: models hvdc pole controls
%
% Input: i - 0 vector computation only for hvdc controls
%        k - integer time for overall simulation
%        kdc - integer time for dc
%        bus - solved loadflow bus data
%        flag - 0 - initialization
%               1 - network interface computation
%               2 - hvdc dynamics computation
%               3 - state matrix building

%-----------------------------------------------------------------------------%
% Version history
%
% Version: 1.0
% Date:    April 1997
% Author:  Graham Rogers
%-----------------------------------------------------------------------------%

global g;  % declaring struct of global variables

% check that dc controls are defined
if ~isempty(g.dc.dcc_con)
    if (flag == 0)
        % initialize controls
        if (i == 0)
            % vector computation
            % rectifier controls
            % calculate current order
            g.dc.cur_ord(g.dc.r_idx,1) = g.dc.i_dc(g.dc.r_idx,1);
            g.dc.cur_ord(g.dc.i_idx,1) = g.dc.i_dc(g.dc.i_idx,1) ...
                                         .*(ones(g.dc.n_dcl,1) ...
                                            - g.dc.dcl_con(:,9)/100.0);

            g.dc.Vdc_ref = g.dc.Vdc(g.dc.i_idx,1);
            g.dc.v_conr(:,1) = g.dc.alpha(:,1)./g.dc.dcc_con(g.dc.r_idx,4);

            vrmax = find(g.dc.v_conr(:,1) > g.dc.dcc_con(g.dc.r_idx,5));
            vrmin = find(g.dc.v_conr(:,1) < g.dc.dcc_con(g.dc.r_idx,6));

            % check rectifier integrator limit violation
            if ~isempty(vrmax)
                estr = '\ndc_cont: v_conr is greater than the maximum limit at ';
                estr = [estr, 'rectifier index %0.0f.'];
                error(sprintf(estr,vrmax));
            end

            if ~isempty(vrmin)
                estr = '\ndc_cont: v_conr is less than the minimum limit at ';
                estr = [estr, 'rectifier index %0.0f.'];
                error(sprintf(estr,vrmin));
            end

            g.dc.v_coni(:,1) = g.dc.gamma(:,1)./g.dc.dcc_con(g.dc.i_idx,4);
            vimax = find(g.dc.v_coni(:,1) > g.dc.dcc_con(g.dc.i_idx,5));
            vimin = find(g.dc.v_coni(:,1) < g.dc.dcc_con(g.dc.i_idx,6));

            % check rectifier integrator limit violation
            if ~isempty(vimax)
                estr = '\ndc_cont: v_coni is greater than the maximum limit at ';
                estr = [estr, 'inverter index %0.0f.'];
                error(sprintf(estr,vimax));
            end

            if ~isempty(vimin)
                estr = '\ndc_cont: v_coni is less than the minimum limit at ';
                estr = [estr, 'inverter index %0.0f.'];
                error(sprintf(estr,vimin));
            end

            % dcc_pot(:,1)  -- gamma reference value
            % dcc_pot(:,2)  -- rectifier xeqpu
            % dcc_pot(:,3)  -- rectifier Rc
            % dcc_pot(:,4)  -- inverter xeqpu
            % dcc_pot(:,5)  -- inverter Rc
            % dcc_pot(:,6)  -- multiplier ideal rectifier dc voltage
            % dcc_pot(:,7)  -- multiplier ideal inverter dc voltage
            % dcc_pot(:,8)  -- multiplier acpu to dc amps
            % dcc_pot(:,9)  -- rectifier (perhaps rated power, unused)
            % dcc_pot(:,10) -- inverter (perhaps rated power, unused)
            g.dc.dcc_pot(:,1) = g.dc.gamma(:,1);
            g.dc.dcc_pot(:,2) = g.dc.dcsp_con(g.dc.r_idx,5) ...
                                ./g.dc.dcsp_con(g.dc.r_idx,6);
            g.dc.dcc_pot(:,2) = g.dc.dcc_pot(:,2)*g.sys.basmva ...
                                ./bus(g.dc.rec_ac_bus,13)./bus(g.dc.rec_ac_bus,13);
            g.dc.dcc_pot(:,3) = g.dc.dcsp_con(g.dc.r_idx,6) ...
                                .*g.dc.dcsp_con(g.dc.r_idx,5)*3/pi;
            g.dc.dcc_pot(:,4) = g.dc.dcsp_con(g.dc.i_idx,5) ...
                                ./g.dc.dcsp_con(g.dc.i_idx,6);
            g.dc.dcc_pot(:,4) = g.dc.dcc_pot(:,4)*g.sys.basmva ...
                                ./bus(g.dc.inv_ac_bus,13)./bus(g.dc.inv_ac_bus,13);
            g.dc.dcc_pot(:,5) = g.dc.dcsp_con(g.dc.i_idx,6) ...
                                .*g.dc.dcsp_con(g.dc.i_idx,5)*3/pi;
            g.dc.dcc_pot(:,6) = g.dc.dcc_pot(:,5);
            g.dc.dcc_pot(:,7) = 3*sqrt(2).*g.dc.dcsp_con(g.dc.r_idx,6) ...
                                .*bus(g.dc.rec_ac_bus,13)/pi;
            g.dc.dcc_pot(:,8) = 3*sqrt(2).*g.dc.dcsp_con(g.dc.i_idx,6) ...
                                .*bus(g.dc.inv_ac_bus,13)/pi;
            % g.dc.dcc_pot(:,9) = pi*g.sys.basmva/3/sqrt(2) ...
            %                     ./bus(g.dc.rec_ac_bus,13) ...
            %                     ./g.dc.dcsp_con(g.dc.r_idx,6);
            % g.dc.dcc_pot(:,10) = pi*g.sys.basmva/3/sqrt(2) ...
            %                      ./bus(g.dc.inv_ac_bus,13) ...
            %                      ./g.dc.dcsp_con(g.dc.i_idx,6);

            dc_dsig(:,1) = zeros(g.dc.n_conv,1);  % init. damping control signals
            if (g.dc.dc_sig(:,1) ~= zeros(g.dc.n_conv,1))
                % reset initial values of alpha and gamma
                g.dc.alpha(:,1) = ...
                    ((-g.dc.cur_ord(g.dc.r_idx,1) ...
                      + g.dc.dc_sig(g.dc.r_idx,1) + g.dc.dcr_dsig(:,1) ...
                      + g.dc.i_dcr(:,1)).*g.dc.dcc_con(g.dc.r_idx,2) ...
                     + g.dc.v_conr(:,1)).*g.dc.dcc_con(g.dc.r_idx,4);

                g.dc.gamma(:,1) = ...
                    ((g.dc.Vdc(g.dc.i_idx,1) - g.dc.Vdc_ref) ...
                     ./g.dc.Vdc_ref.*g.dc.dcc_con(g.dc.i_idx,2) ...
                     + g.dc.v_coni(:,1)).*g.dc.dcc_con(g.dc.i_idx,4);
            end
        else
            error('dc_cont: initialization must be vectorized.');
        end

    end  % end of initialization

    if (flag == 1)
        % network interface
        if (i ~= 0)
            error('dc_cont: network interface calculation must be vectorized.');
        else
            % vector computation
            % i_dc, v_dcc, and the control states are fixed

            % determine firing and extinction angles
            g.dc.alpha(:,kdc) = ...
                ((-g.dc.cur_ord(g.dc.r_idx,k) ...
                  + g.dc.dc_sig(g.dc.r_idx,k) + g.dc.dcr_dsig(:,k) ...
                  + g.dc.i_dcr(:,kdc)).*g.dc.dcc_con(g.dc.r_idx,2) ...
                 + g.dc.v_conr(:,kdc)).*g.dc.dcc_con(g.dc.r_idx,4);

            % check for alpha limits
            g.dc.alpha(:,kdc) = max(g.dc.alpha(:,kdc), ...
                                    g.dc.dcc_con(g.dc.r_idx,8)*pi/180);
            g.dc.alpha(:,kdc) = min(g.dc.alpha(:,kdc), ...
                                    g.dc.dcc_con(g.dc.r_idx,7)*pi/180);

            g.dc.gamma(:,kdc) = ((g.dc.Vdc(g.dc.i_idx,kdc) - g.dc.Vdc_ref) ...
                                 ./g.dc.Vdc_ref.*g.dc.dcc_con(g.dc.i_idx,2) ...
                                 + g.dc.v_coni(:,kdc)).*g.dc.dcc_con(g.dc.i_idx,4);

            cur_error = g.dc.i_dci(:,kdc) - g.dc.cur_ord(g.dc.i_idx,k);
            ce_idx = find(cur_error < 0);
            if ~isempty(ce_idx)
                g.dc.gamma(ce_idx,kdc) = ...
                    g.dc.gamma(ce_idx,kdc) ...
                    + cur_error(ce_idx).*g.dc.dcc_con(g.dc.i_idx(ce_idx),2) ...
                      .*g.dc.dcc_con(g.dc.i_idx(ce_idx),4);
            end

            % check gamma limits
            g.dc.gamma(:,kdc) = max(g.dc.gamma(:,kdc), ...
                                    g.dc.dcc_con(g.dc.i_idx,8)*pi/180);
            g.dc.gamma(:,kdc) = min(g.dc.gamma(:,kdc), ...
                                    g.dc.dcc_con(g.dc.i_idx,7)*pi/180);
        end
    end

    if (flag == 2)
        % calculate rates of change of states
        if (i == 0)  % vectorized computation
            % rectifier
            g.dc.dv_conr(:,kdc) = ...
                (-g.dc.cur_ord(g.dc.r_idx,k) ...
                 + g.dc.i_dcr(:,kdc) + g.dc.dc_sig(g.dc.r_idx,k) ...
                 + g.dc.dcr_dsig(:,k)).*g.dc.dcc_con(g.dc.r_idx,3);

            % check for state limits
            recmx = find(g.dc.v_conr(:,kdc) > g.dc.dcc_con(g.dc.r_idx,5));
            if ~isempty(recmx)
                g.dc.v_conr(recmx,kdc) = g.dc.dcc_con(g.dc.r_idx(recmx),5);
                recdmx = find(g.dc.dv_conr(recmx,k) > 0);
                if ~isempty(recdmx)
                    g.dc.dv_conr(recmx(recdmx),kdc) = zeros(length(recdmx),1);
                end
            end

            recmn = find(g.dc.v_conr(:,kdc) < g.dc.dcc_con(g.dc.r_idx,6));
            if ~isempty(recmn)
                g.dc.v_conr(recmn,kdc) = g.dc.dcc_con(g.dc.r_idx(recmn),6);
                recdmn = find(g.dc.dv_conr(recmn,kdc) < 0);
                if ~isempty(recdmn)
                    g.dc.dv_conr(recmn(recdmn),kdc) = zeros(length(recdmn),1);
                end
            end

            % inverter
            cur_err = g.dc.cur_ord(g.dc.i_idx,k) - g.dc.i_dci(:,kdc);
            n_ce = find(cur_err< 0);
            if ~isempty(n_ce)
                cur_err(n_ce) = zeros(length(n_ce),1);
            end

            inv_err = (g.dc.Vdc(g.dc.i_idx,kdc) - g.dc.Vdc_ref) ...
                      ./g.dc.Vdc_ref - cur_err;

            g.dc.dv_coni(:,kdc) = (inv_err + g.dc.dc_sig(g.dc.i_idx,k) ...
                                   + g.dc.dci_dsig(:,k)).*g.dc.dcc_con(g.dc.i_idx,3);

            % check state limits
            invmx = find(g.dc.v_coni(:,kdc) > g.dc.dcc_con(g.dc.i_idx,5));
            if ~isempty(invmx)
                g.dc.v_coni(invmx,kdc) = g.dc.dcc_con(g.dc.i_idx(invmx),5);
                invdmx = find(g.dc.dv_coni(invmx,k) > 0);
                if ~isempty(invdmx)
                    g.dc.dv_coni(recmx(invdmx),kdc) = zeros(length(invdmx),1);
                end
            end

            invmn = find(g.dc.v_coni(:,kdc) < g.dc.dcc_con(g.dc.i_idx,6));
            if ~isempty(invmn)
                g.dc.v_coni(invmn,kdc) = g.dc.dcc_con(g.dc.i_idx(invmn),6);
                invdmn = find(g.dc.dv_coni(invmn,kdc) < 0);
                if ~isempty(invdmn)
                    g.dc.dv_conr(invmn(invdmn),kdc) = zeros(length(invdmn),1);
                end
            end
        else
            error('dc_cont: dynamics calculation must be vectorized.');
        end
    end
end

end  % function end

% eof
