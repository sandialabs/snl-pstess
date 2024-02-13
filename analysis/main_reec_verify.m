%-----------------------------------------------------------------------------%
% reec model verification script
%
% Purpose: This is a script for verifying that the differential equations
%          encoded in the reec model match those in the block diagram. It
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
% January 2023
%-----------------------------------------------------------------------------%

%-----------------------------------------------------------------------------%
% reec1    terminal voltage transducer
% reec2    active power transducer
% reec3    reactive power transducer
% reec4    reactive PI loop integrator
% reec5    voltage PI loop integrator
% reec6    reactive current order filter (PI bypass)
% reec7    active power order filter
% reec8    voltage compensation filter (vcmpflag)
%-----------------------------------------------------------------------------%

clear all; close all; clc;

load('./mat/reec_example1_nonlinear_sim.mat');  % you may choose any valid mat file
% load('./mat/reec_example2_nonlinear_sim.mat');

Fs = round(1/median(diff(t)));
lbnd = 1e-3;
t_ax = [0,0.6*t(end)];  % time window for plotting (axes only)
y_stretch = 0.08;       % stretch the y-axis to make the plots easier to read
busnum = g.bus.bus_int(g.reec.reec_con(:,2));

%-----------------------------------------------------------------------------%
% reec1 -- terminal voltage transducer

vreec1 = zeros(size(g.reec.reec1));
max_pct_err = 0;

for ii = 1:size(g.reec.reec_con,1)
    if (g.reec.reec_con(ii,5) > lbnd)
        b = [0,1];
        a = [g.reec.reec_con(ii,5),1];

        [bd,ad] = bilinear(b,a,Fs);
    else
        bd = 1;
        ad = 1;
    end

    vt_raw = abs(g.bus.bus_v(busnum(ii),:));
    U = vt_raw(1)*ones(3,1);
    Y = vt_raw(1)*ones(3,1);
    Z = filtic(bd,ad,Y,U);
    vreec1(ii,:) = filter(bd,ad,vt_raw,Z);

    tmp_err = norm(vreec1(ii,:)-g.reec.reec1(ii,:),'inf');
    tmp_pct_err = tmp_err;

    pct_den = max(abs(vreec1(ii,:)-vreec1(ii,1)));
    if (pct_den > 0)
        tmp_pct_err = 100*tmp_pct_err/pct_den;
    end

    if (tmp_pct_err >= max_pct_err)
        reec1_err = tmp_err;
        max_pct_err = tmp_pct_err;
    end
end

fprintf('\n reec1_err: %0.2e, %0.2f%%',reec1_err,max_pct_err);

fig{1} = figure;
ax{1}(1) = subplot(2,1,1,'parent',fig{1});
ax{1}(2) = subplot(2,1,2,'parent',fig{1});
hold(ax{1}(1),'on');
plot(ax{1}(1),t,vreec1);
plot(ax{1}(1),t,g.reec.reec1,'--');
plot(ax{1}(2),t,vreec1-g.reec.reec1);

for ii = 1:2
    v = axis(ax{1}(ii));
    y_span = v(4)-v(3);
    alpha = y_stretch*y_span;
    axis(ax{1}(ii),[t_ax(1),t_ax(2),v(3)-alpha,v(4)+alpha]);
end

xlabel(ax{1}(2), 'Time (s)');
ylabel(ax{1}(1), 'reec1 (pu)');
ylabel(ax{1}(2), 'error (pu)');

%-----------------------------------------------------------------------------%
% reec2 -- active power transducer

vreec2 = zeros(size(g.reec.reec2));
max_pct_err = 0;

for ii = 1:size(g.reec.reec_con,1)
    if (g.reec.reec_con(ii,14) > lbnd)
        b = [0,1];
        a = [g.reec.reec_con(ii,14),1];

        [bd,ad] = bilinear(b,a,Fs);
    else
        bd = 1;
        ad = 1;
    end

    pgen = real(g.ess.ess_sinj(g.reec.ess_idx(ii),:)).*g.reec.reec_pot(ii,1);
    U = pgen(1)*ones(3,1);
    Y = pgen(1)*ones(3,1);
    Z = filtic(bd,ad,Y,U);
    vreec2(ii,:) = filter(bd,ad,pgen,Z);

    tmp_err = norm(vreec2(ii,:)-g.reec.reec2(ii,:),'inf');
    tmp_pct_err = tmp_err;

    pct_den = max(abs(vreec2(ii,:)-vreec2(ii,1)));
    if (pct_den > 0)
        tmp_pct_err = 100*tmp_pct_err/pct_den;
    end

    if (tmp_pct_err >= max_pct_err)
        reec2_err = tmp_err;
        max_pct_err = tmp_pct_err;
    end
end

fprintf('\n reec2_err: %0.2e, %0.2f%%',reec2_err,max_pct_err);

fig{2} = figure;
ax{2}(1) = subplot(2,1,1,'parent',fig{2});
ax{2}(2) = subplot(2,1,2,'parent',fig{2});
hold(ax{2}(1),'on');
plot(ax{2}(1),t,vreec2);
plot(ax{2}(1),t,g.reec.reec2,'--');
plot(ax{2}(2),t,vreec2-g.reec.reec2);

for ii = 1:2
    v = axis(ax{2}(ii));
    y_span = v(4)-v(3);
    alpha = y_stretch*y_span;
    axis(ax{2}(ii),[t_ax(1),t_ax(2),v(3)-alpha,v(4)+alpha]);
end

xlabel(ax{2}(2), 'Time (s)');
ylabel(ax{2}(1), 'reec2 (pu)');
ylabel(ax{2}(2), 'error (pu)');

%-----------------------------------------------------------------------------%
% reec3 -- reactive power transducer

vreec3 = zeros(size(g.reec.reec3));
max_pct_err = 0;

for ii = 1:size(g.reec.reec_con,1)
    if (g.reec.reec_con(ii,15) > lbnd)
        b = [0,1];
        a = [g.reec.reec_con(ii,15),1];

        [bd,ad] = bilinear(b,a,Fs);
    else
        bd = 1;
        ad = 1;
    end

    qgen = imag(g.ess.ess_sinj(g.reec.ess_idx(ii),:)).*g.reec.reec_pot(ii,1);
    U = qgen(1)*ones(3,1);
    Y = qgen(1)*ones(3,1);
    Z = filtic(bd,ad,Y,U);
    vreec3(ii,:) = filter(bd,ad,qgen,Z);

    tmp_err = norm(vreec3(ii,:)-g.reec.reec3(ii,:),'inf');
    tmp_pct_err = tmp_err;

    pct_den = max(abs(vreec3(ii,:)-vreec3(ii,1)));
    if (pct_den > 0)
        tmp_pct_err = 100*tmp_pct_err/pct_den;
    end

    if (tmp_pct_err >= max_pct_err)
        reec3_err = tmp_err;
        max_pct_err = tmp_pct_err;
    end
end

fprintf('\n reec3_err: %0.2e, %0.2f%%',reec3_err,max_pct_err);

fig{3} = figure;
ax{3}(1) = subplot(2,1,1,'parent',fig{3});
ax{3}(2) = subplot(2,1,2,'parent',fig{3});
hold(ax{3}(1),'on');
plot(ax{3}(1),t,vreec3);
plot(ax{3}(1),t,g.reec.reec3,'--');
plot(ax{3}(2),t,vreec3-g.reec.reec3);

for ii = 1:2
    v = axis(ax{3}(ii));
    y_span = v(4)-v(3);
    alpha = y_stretch*y_span;
    axis(ax{3}(ii),[t_ax(1),t_ax(2),v(3)-alpha,v(4)+alpha]);
end

xlabel(ax{3}(2), 'Time (s)');
ylabel(ax{3}(1), 'reec3 (pu)');
ylabel(ax{3}(2), 'error (pu)');

%-----------------------------------------------------------------------------%
% reec4 -- reactive PI loop integrator

vreec4 = zeros(size(g.reec.reec4));
max_pct_err = 0;

for ii = 1:size(g.reec.reec_con,1)
    b = [0,1].*g.reec.reec_con(ii,22);
    a = [1,0];

    [bd,ad] = bilinear(b,a,Fs);

    qlim = g.reec.qref(ii,:);

    pfflag = g.reec.reec_con(ii,32);
    if (pfflag == 1)
        qlim = g.reec.reec2(ii,:).*tan(g.reec.pfaref(ii));
    end

    qlim = min(qlim,g.reec.reec_con(ii,16));
    qlim = max(qlim,g.reec.reec_con(ii,17));
    u_pi = qlim - g.reec.reec3(ii,:) + imag(g.reec.reec_sig(ii,:));

    U = u_pi(1)*ones(3,1);
    Y = g.reec.reec4(ii,1)*ones(3,1);
    Z = filtic(bd,ad,Y,U);
    vreec4(ii,:) = filter(bd,ad,u_pi,Z);

    tmp_err = norm(vreec4(ii,:)-g.reec.reec4(ii,:),'inf');
    tmp_pct_err = tmp_err;

    pct_den = max(abs(vreec4(ii,:)-vreec4(ii,1)));
    if (pct_den > 0)
        tmp_pct_err = 100*tmp_pct_err/pct_den;
    end

    if (tmp_pct_err >= max_pct_err)
        reec4_err = tmp_err;
        max_pct_err = tmp_pct_err;
    end
end

fprintf('\n reec4_err: %0.2e, %0.2f%%',reec4_err,max_pct_err);

fig{4} = figure;
ax{4}(1) = subplot(2,1,1,'parent',fig{4});
ax{4}(2) = subplot(2,1,2,'parent',fig{4});
hold(ax{4}(1),'on');
plot(ax{4}(1),t,vreec4);
plot(ax{4}(1),t,g.reec.reec4,'--');
plot(ax{4}(2),t,vreec4-g.reec.reec4);

for ii = 1:2
    v = axis(ax{4}(ii));
    y_span = v(4)-v(3);
    alpha = y_stretch*y_span;
    axis(ax{4}(ii),[t_ax(1),t_ax(2),v(3)-alpha,v(4)+alpha]);
end

xlabel(ax{4}(2), 'Time (s)');
ylabel(ax{4}(1), 'reec4 (pu)');
ylabel(ax{4}(2), 'error (pu)');

%-----------------------------------------------------------------------------%
% reec5 -- reactive PI loop integrator

vreec5 = zeros(size(g.reec.reec5));
max_pct_err = 0;

for ii = 1:size(g.reec.reec_con,1)
    b = [0,1].*g.reec.reec_con(ii,24);
    a = [1,0];

    [bd,ad] = bilinear(b,a,Fs);

    u_pi = g.reec.dreec5(ii,:)./g.reec.reec_con(ii,24);

    U = u_pi(1)*ones(3,1);
    Y = g.reec.reec5(ii,1)*ones(3,1);
    Z = filtic(bd,ad,Y,U);
    vreec5(ii,:) = filter(bd,ad,u_pi,Z);

    tmp_err = norm(vreec5(ii,:)-g.reec.reec5(ii,:),'inf');
    tmp_pct_err = tmp_err;

    pct_den = max(abs(vreec5(ii,:)-vreec5(ii,1)));
    if (pct_den > 0)
        tmp_pct_err = 100*tmp_pct_err/pct_den;
    end

    if (tmp_pct_err >= max_pct_err)
        reec5_err = tmp_err;
        max_pct_err = tmp_pct_err;
    end
end

fprintf('\n reec5_err: %0.2e, %0.2f%%',reec5_err,max_pct_err);

fig{5} = figure;
ax{5}(1) = subplot(2,1,1,'parent',fig{5});
ax{5}(2) = subplot(2,1,2,'parent',fig{5});
hold(ax{5}(1),'on');
plot(ax{5}(1),t,vreec5);
plot(ax{5}(1),t,g.reec.reec5,'--');
plot(ax{5}(2),t,vreec5-g.reec.reec5);

for ii = 1:2
    v = axis(ax{5}(ii));
    y_span = v(4)-v(3);
    alpha = y_stretch*y_span;
    axis(ax{5}(ii),[t_ax(1),t_ax(2),v(3)-alpha,v(4)+alpha]);
end

xlabel(ax{5}(2), 'Time (s)');
ylabel(ax{5}(1), 'reec5 (pu)');
ylabel(ax{5}(2), 'error (pu)');

%-----------------------------------------------------------------------------%
% reec6 -- reactive current order filter (PI bypass)

vreec6 = zeros(size(g.reec.reec6));
max_pct_err = 0;

for ii = 1:size(g.reec.reec_con,1)
    if (g.reec.reec_con(ii,25) > lbnd)
        b = [0,1];
        a = [g.reec.reec_con(ii,25),1];

        [bd,ad] = bilinear(b,a,Fs);
    else
        bd = 1;
        ad = 1;
    end

    qlim = g.reec.qref(ii,:);

    pfflag = g.reec.reec_con(ii,32);
    if any(pfflag)
        qlim(pfflag) = g.reec.reec2(pfflag,:).*tan(g.reec.pfaref(pfflag));
    end

    qlim = min(qlim,g.reec.reec_con(ii,16));
    qlim = max(qlim,g.reec.reec_con(ii,17));

    u_filt = qlim./max(g.reec.reec1(ii,:),0.01);

    U = u_filt(1)*ones(3,1);
    Y = u_filt(1)*ones(3,1);
    Z = filtic(bd,ad,Y,U);
    vreec6(ii,:) = filter(bd,ad,u_filt,Z);

    tmp_err = norm(vreec6(ii,:)-g.reec.reec6(ii,:),'inf');
    tmp_pct_err = tmp_err;

    pct_den = max(abs(vreec6(ii,:)-vreec6(ii,1)));
    if (pct_den > 0)
        tmp_pct_err = 100*tmp_pct_err/pct_den;
    end

    if (tmp_pct_err >= max_pct_err)
        reec6_err = tmp_err;
        max_pct_err = tmp_pct_err;
    end
end

fprintf('\n reec6_err: %0.2e, %0.2f%%',reec6_err,max_pct_err);

fig{6} = figure;
ax{6}(1) = subplot(2,1,1,'parent',fig{6});
ax{6}(2) = subplot(2,1,2,'parent',fig{6});
hold(ax{6}(1),'on');
plot(ax{6}(1),t,vreec6);
plot(ax{6}(1),t,g.reec.reec6,'--');
plot(ax{6}(2),t,vreec6-g.reec.reec6);

for ii = 1:2
    v = axis(ax{6}(ii));
    y_span = v(4)-v(3);
    alpha = y_stretch*y_span;
    axis(ax{6}(ii),[t_ax(1),t_ax(2),v(3)-alpha,v(4)+alpha]);
end

xlabel(ax{6}(2), 'Time (s)');
ylabel(ax{6}(1), 'reec6 (pu)');
ylabel(ax{6}(2), 'error (pu)');

%-----------------------------------------------------------------------------%
% reec7 -- reactive current order filter (PI bypass)

vreec7 = zeros(size(g.reec.reec7));
max_pct_err = 0;

for ii = 1:size(g.reec.reec_con,1)
    if (g.reec.reec_con(ii,30) > lbnd)
        b = [0,1];
        a = [g.reec.reec_con(ii,30),1];

        [bd,ad] = bilinear(b,a,Fs);
    else
        bd = 1;
        ad = 1;
    end

    U = g.reec.pref(ii,1)*ones(3,1);
    Y = g.reec.pref(ii,1)*ones(3,1);
    Z = filtic(bd,ad,Y,U);
    vreec7(ii,:) = filter(bd,ad,g.reec.pref(ii,:),Z);

    tmp_err = norm(vreec7(ii,:)-g.reec.reec7(ii,:),'inf');
    tmp_pct_err = tmp_err;

    pct_den = max(abs(vreec7(ii,:)-vreec7(ii,1)));
    if (pct_den > 0)
        if (pct_den > 1e-10)
            tmp_pct_err = 100*tmp_pct_err/pct_den;
        else
            tmp_pct_err = 0.01;
        end
    end

    if (tmp_pct_err >= max_pct_err)
        reec7_err = tmp_err;
        max_pct_err = tmp_pct_err;
    end
end

fprintf('\n reec7_err: %0.2e, %0.2f%%',reec7_err,max_pct_err);

fig{7} = figure;
ax{7}(1) = subplot(2,1,1,'parent',fig{7});
ax{7}(2) = subplot(2,1,2,'parent',fig{7});
hold(ax{7}(1),'on');
plot(ax{7}(1),t,vreec7);
plot(ax{7}(1),t,g.reec.reec7,'--');
plot(ax{7}(2),t,vreec7-g.reec.reec7);

for ii = 1:2
    v = axis(ax{7}(ii));
    y_span = v(4)-v(3);
    alpha = y_stretch*y_span;
    axis(ax{7}(ii),[t_ax(1),t_ax(2),v(3)-alpha,v(4)+alpha]);
end

xlabel(ax{7}(2), 'Time (s)');
ylabel(ax{7}(1), 'reec7 (pu)');
ylabel(ax{7}(2), 'error (pu)');

%-----------------------------------------------------------------------------%
% reec8 -- voltage compensation filter (vcmpflag)

vreec8 = zeros(size(g.reec.reec8));
max_pct_err = 0;

for ii = 1:size(g.reec.reec_con,1)
    if (g.reec.reec_con(ii,38) > lbnd)
        b = [0,1];
        a = [g.reec.reec_con(ii,38),1];

        [bd,ad] = bilinear(b,a,Fs);
    else
        bd = 1;
        ad = 1;
    end

    u_filt = g.reec.dreec8(ii,:).*max(g.reec.reec_con(ii,38),lbnd) ...
             + g.reec.reec8(ii,:);

    U = u_filt(1)*ones(3,1);
    Y = u_filt(1)*ones(3,1);
    Z = filtic(bd,ad,Y,U);
    vreec8(ii,:) = filter(bd,ad,u_filt,Z);

    tmp_err = norm(vreec8(ii,:)-g.reec.reec8(ii,:),'inf');
    tmp_pct_err = tmp_err;

    pct_den = max(abs(vreec8(ii,:)-vreec8(ii,1)));
    if (pct_den > 0)
        tmp_pct_err = 100*tmp_pct_err/pct_den;
    end

    if (tmp_pct_err >= max_pct_err)
        reec8_err = tmp_err;
        max_pct_err = tmp_pct_err;
    end
end

fprintf('\n reec8_err: %0.2e, %0.2f%%',reec8_err,max_pct_err);

fig{8} = figure;
ax{8}(1) = subplot(2,1,1,'parent',fig{8});
ax{8}(2) = subplot(2,1,2,'parent',fig{8});
hold(ax{8}(1),'on');
plot(ax{8}(1),t,vreec8);
plot(ax{8}(1),t,g.reec.reec8,'--');
plot(ax{8}(2),t,vreec8-g.reec.reec8);

for ii = 1:2
    v = axis(ax{8}(ii));
    y_span = v(4)-v(3);
    alpha = y_stretch*y_span;
    axis(ax{8}(ii),[t_ax(1),t_ax(2),v(3)-alpha,v(4)+alpha]);
end

xlabel(ax{8}(2), 'Time (s)');
ylabel(ax{8}(1), 'reec8 (pu)');
ylabel(ax{8}(2), 'error (pu)');

%-----------------------------------------------------------------------------%
% reec9 -- active current command filter (ipcmd)

vreec9 = zeros(size(g.reec.reec9));
max_pct_err = 0;

for ii = 1:size(g.reec.reec_con,1)
    if (g.reec.reec_con(ii,43) > lbnd)
        b = [0,1];
        a = [g.reec.reec_con(ii,43),1];

        [bd,ad] = bilinear(b,a,Fs);
    else
        bd = 1;
        ad = 1;
    end

    % ipcmd
    u_filt = (g.reec.reec7(ii,:) + g.reec.paux(ii,:)) ...
             ./max(g.reec.reec1(ii,:),0.01);

    U = u_filt(1)*ones(3,1);
    Y = u_filt(1)*ones(3,1);
    Z = filtic(bd,ad,Y,U);
    vreec9(ii,:) = filter(bd,ad,u_filt,Z);

    tmp_err = norm(vreec9(ii,:)-g.reec.reec9(ii,:),'inf');
    tmp_pct_err = tmp_err;

    pct_den = max(abs(vreec9(ii,:)-vreec9(ii,1)));
    if (pct_den > 0)
        tmp_pct_err = 100*tmp_pct_err/pct_den;
    end

    if (tmp_pct_err >= max_pct_err)
        reec9_err = tmp_err;
        max_pct_err = tmp_pct_err;
    end
end

fprintf('\n reec9_err: %0.2e, %0.2f%%',reec9_err,max_pct_err);

fig{9} = figure;
ax{9}(1) = subplot(2,1,1,'parent',fig{9});
ax{9}(2) = subplot(2,1,2,'parent',fig{9});
hold(ax{9}(1),'on');
plot(ax{9}(1),t,vreec9);
plot(ax{9}(1),t,g.reec.reec9,'--');
plot(ax{9}(2),t,vreec9-g.reec.reec9);

for ii = 1:2
    v = axis(ax{9}(ii));
    y_span = v(4)-v(3);
    alpha = y_stretch*y_span;
    axis(ax{9}(ii),[t_ax(1),t_ax(2),v(3)-alpha,v(4)+alpha]);
end

xlabel(ax{9}(2), 'Time (s)');
ylabel(ax{9}(1), 'reec9 (pu)');
ylabel(ax{9}(2), 'error (pu)');

%-----------------------------------------------------------------------------%
% reec10 -- reactive current command filter (iqcmd)

vreec10 = zeros(size(g.reec.reec10));
max_pct_err = 0;

for ii = 1:size(g.reec.reec_con,1)
    if (g.reec.reec_con(ii,44) > lbnd)
        b = [0,1];
        a = [g.reec.reec_con(ii,44),1];

        [bd,ad] = bilinear(b,a,Fs);
    else
        bd = 1;
        ad = 1;
    end

    % local voltage control loop
    verr = g.reec.vref0(ii) - g.reec.reec1(ii,:) + real(g.reec.reec_sig(ii,:));

    if (verr >= g.reec.reec_con(ii,6)) & (verr <= g.reec.reec_con(ii,7))
        verr = 0.0;
    end

    if (verr > 0)
        verr = verr - g.reec.reec_con(ii,7);
    elseif (verr < 0)
        verr = verr - g.reec.reec_con(ii,6);
    end

    iqv = g.reec.reec_con(ii,8).*verr;
    iqinj = min(iqv,g.reec.reec_con(ii,9));
    iqinj = max(iqinj,g.reec.reec_con(ii,10));

    % cascaded PI loops
    qlim = g.reec.qref(ii,:);

    if (g.reec.reec_con(ii,32) == 1)
        qlim = g.reec.reec2(ii,:).*tan(g.reec.pfaref(ii));
    end

    qlim = min(qlim,g.reec.reec_con(ii,16));
    qlim = max(qlim,g.reec.reec_con(ii,17));

    qerr = qlim - g.reec.reec3(ii,:) + imag(g.reec.reec_sig(ii,:));

    y_piq = g.reec.reec_con(ii,21).*qerr + g.reec.reec4(ii,:);

    if (y_piq > g.reec.reec_con(ii,18))
        y_piq = g.reec.reec_con(ii,18);
    end

    if (y_piq < g.reec.reec_con(ii,19));
        y_piq = g.reec.reec_con(ii,19);
    end

    vlim = qlim + g.reec.vref1(ii);

    if (g.reec.reec_con(ii,33) == 1)
        vlim = y_piq;
    end

    vlim = min(vlim,g.reec.reec_con(ii,18));
    vlim = max(vlim,g.reec.reec_con(ii,19));

    y_piv = g.reec.reec_con(ii,23).*(vlim - g.reec.reec8(ii,:)) ...
            + g.reec.reec5(ii,:);

    if (y_piv > g.reec.iqmax(ii))
        y_piv = g.reec.iqmax(ii);
    end

    if (y_piv < g.reec.iqmin(ii))
        y_piv = g.reec.iqmin(ii);
    end

    % summing the results
    iq2 = g.reec.reec6(ii,:);

    if (g.reec.reec_con(ii,34) == 1)
        iq2 = y_piv;
    end

    u_filt = iqinj + iq2;

    U = u_filt(1)*ones(3,1);
    Y = u_filt(1)*ones(3,1);
    Z = filtic(bd,ad,Y,U);
    vreec10(ii,:) = filter(bd,ad,u_filt,Z);

    tmp_err = norm(vreec10(ii,:)-g.reec.reec10(ii,:),'inf');
    tmp_pct_err = tmp_err;

    pct_den = max(abs(vreec10(ii,:)-vreec10(ii,1)));
    if (pct_den > 0)
        tmp_pct_err = 100*tmp_pct_err/pct_den;
    end

    if (tmp_pct_err >= max_pct_err)
        reec10_err = tmp_err;
        max_pct_err = tmp_pct_err;
    end
end

fprintf('\n reec10_err: %0.2e, %0.2f%%',reec10_err,max_pct_err);

fig{10} = figure;
ax{10}(1) = subplot(2,1,1,'parent',fig{10});
ax{10}(2) = subplot(2,1,2,'parent',fig{10});
hold(ax{10}(1),'on');
plot(ax{10}(1),t,vreec10);
plot(ax{10}(1),t,g.reec.reec10,'--');
plot(ax{10}(2),t,vreec10-g.reec.reec10);

for ii = 1:2
    v = axis(ax{10}(ii));
    y_span = v(4)-v(3);
    alpha = y_stretch*y_span;
    axis(ax{10}(ii),[t_ax(1),t_ax(2),v(3)-alpha,v(4)+alpha]);
end

xlabel(ax{10}(2), 'Time (s)');
ylabel(ax{10}(1), 'reec10 (pu)');
ylabel(ax{10}(2), 'error (pu)');

%-----------------------------------------------------------------------------%
% wrapping up

fprintf('\n\n');

% eof
