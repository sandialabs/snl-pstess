%-----------------------------------------------------------------------------%
% lsc model verification script
%
% Purpose: This is a script for verifying that the differential equations
%          encoded in the lsc model match those in the block diagram. It
%          takes the transfer function input to each `block' and passes it
%          through an equivalent discrete-time filter. The worst-case error
%          for each state should be no greater than a few percent, but
%          ultimately the eye test is the most reliable. To compute the
%          error percentage, the maximum absolute error is divided by the
%          maximum state deviation. Due to numerical issues, this is more
%          consistently useful than dividing the error by the value of the
%          state at the sample where the maximum error occurred.
%
% Ryan Elliott, rtellio@sandia.gov
% December 2020
%-----------------------------------------------------------------------------%

%-----------------------------------------------------------------------------%
% lsc1   transducer for remote angle signal
% lsc2   transducer for local angle signal
% lsc3   Pade approx. for remote angle signal
% lsc4   Pade approx. for local angle signal
% lsc5   remote angle signal setpoint tracking filter
% lsc6   local angle signal setpoint tracking filter
% lsc7   local LTV highpass filter state 1
% lsc8   local LTV highpass filter state 2
% lsc9   local LTV lead-lag stage 1
% lsc10  local LTV lead-lag stage 2
% lsc11  center LTI highpass filter state 1
% lsc12  center LTI highpass filter state 2
% lsc13  center LTI lead-lag stage 1
% lsc14  center LTI lead-lag stage 2
% lsc15  lowpass filter
%-----------------------------------------------------------------------------%

clear all; close all; clc;

load('./mat/ess_example1_nonlinear_sim.mat');  % you may choose any valid mat file
% load('./mat/ess_example2_nonlinear_sim.mat');

Fs = round(1/median(diff(t)));
lbnd = 1e-3;
t_ax = [0,0.6*t(end)];  % time window for plotting (axes only)
y_stretch = 0.08;       % stretch the y-axis to make the plots easier to read

%-----------------------------------------------------------------------------%
% lsc1 -- transducer for remote angle signal

vlsc1 = zeros(size(g.lsc.lsc1));
max_pct_err = 0;

for ii = 1:size(g.lsc.lsc_con,1)
    if (g.lsc.lsc_con(ii,3) > lbnd)
        b = [0,1];
        a = [g.lsc.lsc_con(ii,3),1];

        [bd,ad] = bilinear(b,a,Fs);
    else
        bd = 1;
        ad = 1;
    end

    theta_coi_est = mean(g.lsc.theta_sensor,1);
    U = theta_coi_est(1)*ones(3,1);
    Y = theta_coi_est(1)*ones(3,1);
    Z = filtic(bd,ad,Y,U);
    vlsc1(ii,:) = filter(bd,ad,theta_coi_est,Z);

    tmp_err = norm(vlsc1(ii,:)-g.lsc.lsc1(ii,:),'inf');
    tmp_pct_err = tmp_err;

    pct_den = max(abs(vlsc1(ii,:)-vlsc1(ii,1)));
    if (pct_den > 0)
        tmp_pct_err = 100*tmp_pct_err/pct_den;
    end

    if (tmp_pct_err >= max_pct_err)
        lsc1_err = tmp_err;
        max_pct_err = tmp_pct_err;
    end
end

fprintf('\n lsc1_err: %0.2e, %0.2f%%',lsc1_err,max_pct_err);

fig{1} = figure;
ax{1}(1) = subplot(2,1,1,'parent',fig{1});
ax{1}(2) = subplot(2,1,2,'parent',fig{1});
hold(ax{1}(1),'on');
plot(ax{1}(1),t,vlsc1);
plot(ax{1}(1),t,g.lsc.lsc1,'--');
plot(ax{1}(2),t,vlsc1-g.lsc.lsc1);

for ii = 1:2
    v = axis(ax{1}(ii));
    y_span = v(4)-v(3);
    alpha = y_stretch*y_span;
    axis(ax{1}(ii),[t_ax(1),t_ax(2),v(3)-alpha,v(4)+alpha]);
end

xlabel(ax{1}(2), 'Time (s)');
ylabel(ax{1}(1), 'lsc1 (pu)');
ylabel(ax{1}(2), 'error (pu)');

%-----------------------------------------------------------------------------%
% lsc2 -- transducer for local angle signal

vlsc2 = zeros(size(g.lsc.lsc2));
max_pct_err = 0;

for ii = 1:size(g.lsc.lsc_con,1)
    if (g.lsc.lsc_con(ii,3) > lbnd)
        b = [0,1];
        a = [g.lsc.lsc_con(ii,3),1];

        [bd,ad] = bilinear(b,a,Fs);
    else
        bd = 1;
        ad = 1;
    end

    U = g.lsc.theta_sensor(g.lsc.lsc_sensor_idx(ii),1)*ones(3,1);
    Y = g.lsc.theta_sensor(g.lsc.lsc_sensor_idx(ii),1)*ones(3,1);
    Z = filtic(bd,ad,Y,U);
    vlsc2(ii,:) = filter(bd,ad,g.lsc.theta_sensor(g.lsc.lsc_sensor_idx(ii),:),Z);

    tmp_err = norm(vlsc2(ii,:)-g.lsc.lsc2(ii,:),'inf');
    tmp_pct_err = tmp_err;

    pct_den = max(abs(vlsc2(ii,:)-vlsc2(ii,1)));
    if (pct_den > 0)
        tmp_pct_err = 100*tmp_pct_err/pct_den;
    end

    if (tmp_pct_err >= max_pct_err)
        lsc2_err = tmp_err;
        max_pct_err = tmp_pct_err;
    end
end

fprintf('\n lsc2_err: %0.2e, %0.2f%%',lsc2_err,max_pct_err);

fig{2} = figure;
ax{2}(1) = subplot(2,1,1,'parent',fig{2});
ax{2}(2) = subplot(2,1,2,'parent',fig{2});
hold(ax{2}(1),'on');
plot(ax{2}(1),t,vlsc2);
plot(ax{2}(1),t,g.lsc.lsc2,'--');
plot(ax{2}(2),t,vlsc2-g.lsc.lsc2);

for ii = 1:2
    v = axis(ax{2}(ii));
    y_span = v(4)-v(3);
    alpha = y_stretch*y_span;
    axis(ax{2}(ii),[t_ax(1),t_ax(2),v(3)-alpha,v(4)+alpha]);
end

xlabel(ax{2}(2), 'Time (s)');
ylabel(ax{2}(1), 'lsc2 (pu)');
ylabel(ax{2}(2), 'error (pu)');

%-----------------------------------------------------------------------------%
% lsc3 -- Pade approx. for remote angle signal

vlsc3 = zeros(size(g.lsc.lsc3));
max_pct_err = 0;

for ii = 1:size(g.lsc.lsc_con,1)
    if (g.lsc.lsc_con(ii,5)/2 > lbnd)
        b = [-g.lsc.lsc_con(ii,5)/2,1];
        a = [g.lsc.lsc_con(ii,5)/2,1];

        [bd,ad] = bilinear(b,a,Fs);
    else
        bd = 1;
        ad = 1;
    end

    U = g.lsc.lsc1(ii,1)*ones(3,1);
    Y = g.lsc.lsc1(ii,1)*ones(3,1);
    Z = filtic(bd,ad,Y,U);
    vlsc3(ii,:) = filter(bd,ad,g.lsc.lsc1(ii,:),Z);

    tmp_err = norm(vlsc3(ii,:)-g.lsc.theta_coi_pade(ii,:),'inf');
    tmp_pct_err = tmp_err;

    pct_den = max(abs(vlsc3(ii,:)-vlsc3(ii,1)));
    if (pct_den > 0)
        tmp_pct_err = 100*tmp_pct_err/pct_den;
    end

    if (tmp_pct_err >= max_pct_err)
        lsc3_err = tmp_err;
        max_pct_err = tmp_pct_err;
    end
end

fprintf('\n lsc3_err: %0.2e, %0.2f%%',lsc3_err,max_pct_err);

fig{3} = figure;
ax{3}(1) = subplot(2,1,1,'parent',fig{3});
ax{3}(2) = subplot(2,1,2,'parent',fig{3});
hold(ax{3}(1),'on');
plot(ax{3}(1),t,vlsc3);
plot(ax{3}(1),t,g.lsc.theta_coi_pade,'--');
plot(ax{3}(2),t,vlsc3-g.lsc.theta_coi_pade);

for ii = 1:2
    v = axis(ax{3}(ii));
    y_span = v(4)-v(3);
    alpha = y_stretch*y_span;
    axis(ax{3}(ii),[t_ax(1),t_ax(2),v(3)-alpha,v(4)+alpha]);
end

xlabel(ax{3}(2), 'Time (s)');
ylabel(ax{3}(1), 'lsc3 (pu)');
ylabel(ax{3}(2), 'error (pu)');

%-----------------------------------------------------------------------------%
% lsc4 -- Pade approx. for local angle signal

vlsc4 = zeros(size(g.lsc.lsc4));
max_pct_err = 0;

for ii = 1:size(g.lsc.lsc_con,1)
    if (g.lsc.lsc_con(ii,6)/2 > lbnd)
        b = [-g.lsc.lsc_con(ii,6)/2,1];
        a = [g.lsc.lsc_con(ii,6)/2,1];

        [bd,ad] = bilinear(b,a,Fs);
    else
        bd = 1;
        ad = 1;
    end

    U = g.lsc.lsc2(ii,1)*ones(3,1);
    Y = g.lsc.lsc2(ii,1)*ones(3,1);
    Z = filtic(bd,ad,Y,U);
    vlsc4(ii,:) = filter(bd,ad,g.lsc.lsc2(ii,:),Z);

    tmp_err = norm(vlsc4(ii,:)-g.lsc.theta_lsc_pade(ii,:),'inf');
    tmp_pct_err = tmp_err;

    pct_den = max(abs(vlsc4(ii,:)-vlsc4(ii,1)));
    if (pct_den > 0)
        tmp_pct_err = 100*tmp_pct_err/pct_den;
    end

    if (tmp_pct_err >= max_pct_err)
        lsc4_err = tmp_err;
        max_pct_err = tmp_pct_err;
    end
end

fprintf('\n lsc4_err: %0.2e, %0.2f%%',lsc4_err,max_pct_err);

fig{4} = figure;
ax{4}(1) = subplot(2,1,1,'parent',fig{4});
ax{4}(2) = subplot(2,1,2,'parent',fig{4});
hold(ax{4}(1),'on');
plot(ax{4}(1),t,vlsc4);
plot(ax{4}(1),t,g.lsc.theta_lsc_pade,'--');
plot(ax{4}(2),t,vlsc4-g.lsc.theta_lsc_pade);

for ii = 1:2
    v = axis(ax{4}(ii));
    y_span = v(4)-v(3);
    alpha = y_stretch*y_span;
    axis(ax{4}(ii),[t_ax(1),t_ax(2),v(3)-alpha,v(4)+alpha]);
end

xlabel(ax{4}(2), 'Time (s)');
ylabel(ax{4}(1), 'lsc4 (pu)');
ylabel(ax{4}(2), 'error (pu)');

%-----------------------------------------------------------------------------%
% lsc5 -- remote angle signal setpoint tracking filter

vlsc5 = zeros(size(g.lsc.lsc5));
max_pct_err = 0;

for ii = 1:size(g.lsc.lsc_con,1)
    if (g.lsc.lsc_con(ii,8) > lbnd)
        b = [0,1];
        a = [g.lsc.lsc_con(ii,8),1];

        [bd,ad] = bilinear(b,a,Fs);
    else
        bd = 1;
        ad = 1;
    end

    u_coi_set(ii,:) = g.lsc.dlsc5(ii,:).*max(g.lsc.lsc_con(ii,8),lbnd) ...
                      + g.lsc.lsc5(ii,:);

    U = u_coi_set(ii,1)*ones(3,1);
    Y = u_coi_set(ii,1)*ones(3,1);
    Z = filtic(bd,ad,Y,U);
    vlsc5(ii,:) = filter(bd,ad,u_coi_set(ii,:),Z);

    tmp_err = norm(vlsc5(ii,:)-g.lsc.lsc5(ii,:),'inf');
    tmp_pct_err = tmp_err;

    pct_den = max(abs(vlsc5(ii,:)-vlsc5(ii,1)));
    if (pct_den > 0)
        tmp_pct_err = 100*tmp_pct_err/pct_den;
    end

    if (tmp_pct_err >= max_pct_err)
        lsc5_err = tmp_err;
        max_pct_err = tmp_pct_err;
    end
end

fprintf('\n lsc5_err: %0.2e, %0.2f%%',lsc5_err,max_pct_err);

fig{5} = figure;
ax{5}(1) = subplot(2,1,1,'parent',fig{5});
ax{5}(2) = subplot(2,1,2,'parent',fig{5});
hold(ax{5}(1),'on');
plot(ax{5}(1),t,vlsc5);
plot(ax{5}(1),t,g.lsc.lsc5,'--');
plot(ax{5}(2),t,vlsc5-g.lsc.lsc5);

for ii = 1:2
    v = axis(ax{5}(ii));
    y_span = v(4)-v(3);
    alpha = y_stretch*y_span;
    axis(ax{5}(ii),[t_ax(1),t_ax(2),v(3)-alpha,v(4)+alpha]);
end

xlabel(ax{5}(2), 'Time (s)');
ylabel(ax{5}(1), 'lsc5 (pu)');
ylabel(ax{5}(2), 'error (pu)');

%-----------------------------------------------------------------------------%
% lsc6 -- local angle signal setpoint tracking filter

vlsc6 = zeros(size(g.lsc.lsc6));
max_pct_err = 0;

for ii = 1:size(g.lsc.lsc_con,1)
    if (g.lsc.lsc_con(ii,10) > lbnd)
        b = [0,1];
        a = [g.lsc.lsc_con(ii,10),1];

        [bd,ad] = bilinear(b,a,Fs);
    else
        bd = 1;
        ad = 1;
    end

    u_local_set(ii,:) = g.lsc.dlsc6(ii,:).*max(g.lsc.lsc_con(ii,10),lbnd) ...
                        + g.lsc.lsc6(ii,:);

    U = u_local_set(ii,1)*ones(3,1);
    Y = u_local_set(ii,1)*ones(3,1);
    Z = filtic(bd,ad,Y,U);
    vlsc6(ii,:) = filter(bd,ad,u_local_set(ii,:),Z);

    tmp_err = norm(vlsc6(ii,:)-g.lsc.lsc6(ii,:),'inf');
    tmp_pct_err = tmp_err;

    pct_den = max(abs(vlsc6(ii,:)-vlsc6(ii,1)));
    if (pct_den > 0)
        tmp_pct_err = 100*tmp_pct_err/pct_den;
    end

    if (tmp_pct_err >= max_pct_err)
        lsc6_err = tmp_err;
        max_pct_err = tmp_pct_err;
    end
end

fprintf('\n lsc6_err: %0.2e, %0.2f%%',lsc6_err,max_pct_err);

fig{6} = figure;
ax{6}(1) = subplot(2,1,1,'parent',fig{6});
ax{6}(2) = subplot(2,1,2,'parent',fig{6});
hold(ax{6}(1),'on');
plot(ax{6}(1),t,vlsc6);
plot(ax{6}(1),t,g.lsc.lsc6,'--');
plot(ax{6}(2),t,vlsc6-g.lsc.lsc6);

for ii = 1:2
    v = axis(ax{6}(ii));
    y_span = v(4)-v(3);
    alpha = y_stretch*y_span;
    axis(ax{6}(ii),[t_ax(1),t_ax(2),v(3)-alpha,v(4)+alpha]);
end

xlabel(ax{6}(2), 'Time (s)');
ylabel(ax{6}(1), 'lsc6 (pu)');
ylabel(ax{6}(2), 'error (pu)');

%-----------------------------------------------------------------------------%
% lsc7 -- local LTV highpass filter state 1
% lsc8 -- local LTV highpass filter state 2

vlsc7 = zeros(size(g.lsc.lsc1));
max_pct_err = 0;

for ii = 1:size(g.lsc.lsc_con,1)
    b = [0,g.lsc.lsc_con(ii,12),g.lsc.lsc_con(ii,14)];
    a = [1,g.lsc.lsc_con(ii,13),g.lsc.lsc_con(ii,15)];

    if ((g.lsc.lsc_con(ii,13) == 0) && (g.lsc.lsc_con(ii,15) == 0))
        error('main_lsc_verify: the washout filters may not be bypassed.');
    end

    b = conv(a,1) + b;
    [bd,ad] = bilinear(b,a,Fs);

    u_a1hp(ii,:) = g.lsc.lsc6(ii,:) - g.lsc.lsc5(ii,:) ...
                   - g.lsc.lsc2(ii,:) + g.lsc.lsc1(ii,:);

    if (g.lsc.lsc_con(ii,4) ~= 0)  % check Pade flag
        u_a1hp(ii,:) = g.lsc.lsc6(ii,:) - g.lsc.lsc5(ii,:) ...
                       - g.lsc.theta_lsc_pade(ii,:) ...
                       + g.lsc.theta_coi_pade(ii,:);
    end

    U = u_a1hp(ii,1)*ones(3,1);
    Y = u_a1hp(ii,1)*ones(3,1);
    Z = filtic(bd,ad,Y,U);
    vlsc7(ii,:) = filter(bd,ad,u_a1hp(ii,:),Z);

    tmp_err = norm(vlsc7(ii,:)-(g.lsc.lsc7(ii,:)+u_a1hp(ii,:)),'inf');
    tmp_pct_err = tmp_err;

    pct_den = max(abs(vlsc7(ii,:)-vlsc7(ii,1)));
    if (pct_den > 0)
        tmp_pct_err = 100*tmp_pct_err/pct_den;
    end

    if (tmp_pct_err >= max_pct_err)
        lsc7_err = tmp_err;
        max_pct_err = tmp_pct_err;
    end
end

fprintf('\n lsc7_err: %0.2e, %0.2f%%',lsc7_err,max_pct_err);
fprintf('\n lsc8_err: %0.2e, %0.2f%%',lsc7_err,max_pct_err);

fig{7} = figure;
ax{7}(1) = subplot(2,1,1,'parent',fig{7});
ax{7}(2) = subplot(2,1,2,'parent',fig{7});
hold(ax{7}(1),'on');
plot(ax{7}(1),t,vlsc7);
plot(ax{7}(1),t,g.lsc.lsc7 + u_a1hp,'--');
plot(ax{7}(2),t,vlsc7 - (g.lsc.lsc7+u_a1hp));

for ii = 1:2
    v = axis(ax{7}(ii));
    y_span = v(4)-v(3);
    alpha = y_stretch*y_span;
    axis(ax{7}(ii),[t_ax(1),t_ax(2),v(3)-alpha,v(4)+alpha]);
end

xlabel(ax{7}(2), 'Time (s)');
ylabel(ax{7}(1), 'lsc7 (pu)');
ylabel(ax{7}(2), 'error (pu)');

%-----------------------------------------------------------------------------%
% lsc9  -- local LTV lead-lag stage 1
% lsc10 -- local LTV lead-lag stage 2

vlsc10 = zeros(size(g.lsc.lsc1));
max_pct_err = 0;

for ii = 1:size(g.lsc.lsc_con,1)
    if ((g.lsc.lsc_con(ii,17) == 0) && (g.lsc.lsc_con(ii,19) == 0))
        bd = 1;
        ad = 1;
    else
        b = conv([g.lsc.lsc_con(ii,16),1],[g.lsc.lsc_con(ii,18),1]);
        a = conv([g.lsc.lsc_con(ii,17),1],[g.lsc.lsc_con(ii,19),1]);

        [bd,ad] = bilinear(b,a,Fs);
    end

    u_a1c1(ii,:) = g.lsc.lsc_con(ii,11).*(u_a1hp(ii,:) + g.lsc.lsc7(ii,:));

    tmp_a1c1(ii,:) = g.lsc.lsc9(ii,:) ...
                     .*(1 - g.lsc.lsc_con(ii,16)./max(g.lsc.lsc_con(ii,17),lbnd));

    u_a1c2(ii,:) = ...
        tmp_a1c1(ii,:) ...
        + u_a1c1(ii,:).*(g.lsc.lsc_con(ii,16)./max(g.lsc.lsc_con(ii,17),lbnd));

    tmp_a1c2(ii,:) = g.lsc.lsc10(ii,:) ...
                     .*(1 - g.lsc.lsc_con(ii,18)./max(g.lsc.lsc_con(ii,19),lbnd));

    y_a1c2(ii,:) = ...
        tmp_a1c2(ii,:) ...
        + u_a1c2(ii,:).*(g.lsc.lsc_con(ii,18)./max(g.lsc.lsc_con(ii,19),lbnd));

    U = u_a1c1(1)*ones(3,1);
    Y = u_a1c1(1)*ones(3,1);
    Z = filtic(bd,ad,Y,U);
    vlsc10(ii,:) = filter(bd,ad,u_a1c1(ii,:),Z);

    tmp_err = norm(vlsc10(ii,:)-y_a1c2(ii,:),'inf');
    tmp_pct_err = tmp_err;

    pct_den = max(abs(vlsc10(ii,:)-vlsc10(ii,1)));
    if (pct_den > 0)
        tmp_pct_err = 100*tmp_pct_err/pct_den;
    end

    if (tmp_pct_err >= max_pct_err)
        lsc10_err = tmp_err;
        max_pct_err = tmp_pct_err;
    end
end

fprintf('\n lsc9_err: %0.2e, %0.2f%%',lsc10_err,max_pct_err);
fprintf('\nlsc10_err: %0.2e, %0.2f%%',lsc10_err,max_pct_err);

fig{10} = figure;
ax{10}(1) = subplot(2,1,1,'parent',fig{10});
ax{10}(2) = subplot(2,1,2,'parent',fig{10});
hold(ax{10}(1),'on');
plot(ax{10}(1),t,vlsc10);
plot(ax{10}(1),t,y_a1c2,'--');
plot(ax{10}(2),t,vlsc10 - y_a1c2);

for ii = 1:2
    v = axis(ax{10}(ii));
    y_span = v(4)-v(3);
    alpha = y_stretch*y_span;
    axis(ax{10}(ii),[t_ax(1),t_ax(2),v(3)-alpha,v(4)+alpha]);
end

xlabel(ax{10}(2), 'Time (s)');
ylabel(ax{10}(1), 'lsc10 (pu)');
ylabel(ax{10}(2), 'error (pu)');

%-----------------------------------------------------------------------------%
% lsc11 -- center LTI highpass filter state 1
% lsc12 -- center LTI highpass filter state 2

vlsc11 = zeros(size(g.lsc.lsc1));
max_pct_err = 0;

for ii = 1:size(g.lsc.lsc_con,1)
    b = [0,g.lsc.lsc_con(ii,23),g.lsc.lsc_con(ii,25)];
    a = [1,g.lsc.lsc_con(ii,24),g.lsc.lsc_con(ii,26)];

    if ((g.lsc.lsc_con(ii,24) == 0) && (g.lsc.lsc_con(ii,26) == 0))
        error('main_lsc_verify: the washout filters may not be bypassed.');
    end

    b = conv(a,1) + b;
    [bd,ad] = bilinear(b,a,Fs);

    u_a2hp(ii,:) = g.lsc.lsc5(ii,:) - g.lsc.lsc1(ii,:);  % highpass input

    if (g.lsc.lsc_con(ii,4) ~= 0)                        % check Pade flag
        u_a2hp(ii,:) = g.lsc.lsc5(ii,:) - g.lsc.theta_coi_pade(ii,:);
    end

    U = u_a2hp(1)*ones(3,1);
    Y = u_a2hp(1)*ones(3,1);
    Z = filtic(bd,ad,Y,U);
    vlsc11(ii,:) = filter(bd,ad,u_a2hp(ii,:),Z);

    tmp_err = norm(vlsc11(ii,:)-(g.lsc.lsc11(ii,:)+u_a2hp(ii,:)),'inf');
    tmp_pct_err = tmp_err;

    pct_den = max(abs(vlsc11(ii,:)-vlsc11(ii,1)));
    if (pct_den > 0)
        tmp_pct_err = 100*tmp_pct_err/pct_den;
    end

    if (tmp_pct_err >= max_pct_err)
        lsc11_err = tmp_err;
        max_pct_err = tmp_pct_err;
    end
end

fprintf('\nlsc11_err: %0.2e, %0.2f%%',lsc11_err,max_pct_err);
fprintf('\nlsc12_err: %0.2e, %0.2f%%',lsc11_err,max_pct_err);

fig{11} = figure;
ax{11}(1) = subplot(2,1,1,'parent',fig{11});
ax{11}(2) = subplot(2,1,2,'parent',fig{11});
hold(ax{11}(1),'on');
plot(ax{11}(1),t,vlsc11);
plot(ax{11}(1),t,g.lsc.lsc11 + u_a2hp,'--');
plot(ax{11}(2),t,vlsc11 - (g.lsc.lsc11 + u_a2hp));

for ii = 1:2
    v = axis(ax{11}(ii));
    y_span = v(4)-v(3);
    alpha = y_stretch*y_span;
    axis(ax{11}(ii),[t_ax(1),t_ax(2),v(3)-alpha,v(4)+alpha]);
end

xlabel(ax{11}(2), 'Time (s)');
ylabel(ax{11}(1), 'lsc11 (pu)');
ylabel(ax{11}(2), 'error (pu)');

%-----------------------------------------------------------------------------%
% lsc13 -- center LTI lead-lag stage 1
% lsc14 -- center LTI lead-lag stage 2

vlsc14 = zeros(size(g.lsc.lsc1));
max_pct_err = 0;

for ii = 1:size(g.lsc.lsc_con,1)
    if ((g.lsc.lsc_con(ii,28) == 0) && (g.lsc.lsc_con(ii,30) == 0))
        bd = 1;
        ad = 1;
    else
        b = conv([g.lsc.lsc_con(ii,27),1],[g.lsc.lsc_con(ii,29),1]);
        a = conv([g.lsc.lsc_con(ii,28),1],[g.lsc.lsc_con(ii,30),1]);

        [bd,ad] = bilinear(b,a,Fs);
    end

    u_a2c1(ii,:) = g.lsc.lsc_con(ii,22).*(u_a2hp(ii,:) + g.lsc.lsc11(ii,:));

    tmp_a2c1(ii,:) = g.lsc.lsc13(ii,:) ...
                     .*(1 - g.lsc.lsc_con(ii,27)./max(g.lsc.lsc_con(ii,28),lbnd));

    u_a2c2(ii,:) = ...
        tmp_a2c1(ii,:) ...
        + u_a2c1(ii,:).*(g.lsc.lsc_con(ii,27)./max(g.lsc.lsc_con(ii,28),lbnd));

    tmp_a2c2(ii,:) = g.lsc.lsc14(ii,:) ...
                     .*(1 - g.lsc.lsc_con(ii,29)./max(g.lsc.lsc_con(ii,30),lbnd));

    y_a2c2(ii,:) = ...
        tmp_a2c2(ii,:) ...
        + u_a2c2(ii,:).*(g.lsc.lsc_con(ii,29)./max(g.lsc.lsc_con(ii,30),lbnd));

    U = u_a2c1(1)*ones(3,1);
    Y = u_a2c1(1)*ones(3,1);
    Z = filtic(bd,ad,Y,U);
    vlsc14(ii,:) = filter(bd,ad,u_a2c1(ii,:),Z);

    tmp_err = norm(vlsc14(ii,:)-y_a2c2(ii,:),'inf');
    tmp_pct_err = tmp_err;

    pct_den = max(abs(vlsc14(ii,:)-vlsc14(ii,1)));
    if (pct_den > 0)
        tmp_pct_err = 100*tmp_pct_err/pct_den;
    end

    if (tmp_pct_err >= max_pct_err)
        lsc14_err = tmp_err;
        max_pct_err = tmp_pct_err;
    end
end

fprintf('\nlsc13_err: %0.2e, %0.2f%%',lsc14_err,max_pct_err);
fprintf('\nlsc14_err: %0.2e, %0.2f%%',lsc14_err,max_pct_err);

fig{14} = figure;
ax{14}(1) = subplot(2,1,1,'parent',fig{14});
ax{14}(2) = subplot(2,1,2,'parent',fig{14});
hold(ax{14}(1),'on');
plot(ax{14}(1),t,vlsc14);
plot(ax{14}(1),t,y_a2c2,'--');
plot(ax{14}(2),t,vlsc14 - y_a2c2);

for ii = 1:2
    v = axis(ax{14}(ii));
    y_span = v(4)-v(3);
    alpha = y_stretch*y_span;
    axis(ax{14}(ii),[t_ax(1),t_ax(2),v(3)-alpha,v(4)+alpha]);
end

xlabel(ax{14}(2), 'Time (s)');
ylabel(ax{14}(1), 'lsc14 (pu)');
ylabel(ax{14}(2), 'error (pu)');

%-----------------------------------------------------------------------------%
% lsc15 -- lowpass filter

vlsc15 = zeros(size(g.lsc.lsc15));
max_pct_err = 0;

for ii = 1:size(g.lsc.lsc_con,1)
    if (g.lsc.lsc_con(ii,34) > lbnd)
        b = [0,1];
        a = [g.lsc.lsc_con(ii,34),1];

        [bd,ad] = bilinear(b,a,Fs);
    else
        bd = 1;
        ad = 1;
    end

    u_lp(ii,:) = g.lsc.dlsc15(ii,:).*max(g.lsc.lsc_con(ii,34),lbnd) ...
                 + g.lsc.lsc15(ii,:);

    U = u_lp(ii,1)*ones(3,1);
    Y = u_lp(ii,1)*ones(3,1);
    Z = filtic(bd,ad,Y,U);
    vlsc15(ii,:) = filter(bd,ad,u_lp(ii,:),Z);

    tmp_err = norm(vlsc15(ii,:)-g.lsc.lsc15(ii,:),'inf');
    tmp_pct_err = tmp_err;

    pct_den = max(abs(vlsc15(ii,:)-vlsc15(ii,1)));
    if (pct_den > 0)
        tmp_pct_err = 100*tmp_pct_err/pct_den;
    end

    if (tmp_pct_err >= max_pct_err)
        lsc15_err = tmp_err;
        max_pct_err = tmp_pct_err;
    end
end

fprintf('\nlsc15_err: %0.2e, %0.2f%%',lsc15_err,max_pct_err);

fig{15} = figure;
ax{15}(1) = subplot(2,1,1,'parent',fig{15});
ax{15}(2) = subplot(2,1,2,'parent',fig{15});
hold(ax{15}(1),'on');
plot(ax{15}(1),t,vlsc15);
plot(ax{15}(1),t,g.lsc.lsc15,'--');
plot(ax{15}(2),t,vlsc15-g.lsc.lsc15);

for ii = 1:2
    v = axis(ax{15}(ii));
    y_span = v(4)-v(3);
    alpha = y_stretch*y_span;
    axis(ax{15}(ii),[t_ax(1),t_ax(2),v(3)-alpha,v(4)+alpha]);
end

xlabel(ax{15}(2), 'Time (s)');
ylabel(ax{15}(1), 'lsc15 (pu)');
ylabel(ax{15}(2), 'error (pu)');

fprintf('\n\n');

% eof
