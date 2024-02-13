%-----------------------------------------------------------------------------%
% gfma model verification script
%
% Purpose: This is a script for verifying that the differential equations
%          encoded in the gfma model match those in the block diagram. It
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
% October 2022
%-----------------------------------------------------------------------------%

%-----------------------------------------------------------------------------%
% gfma1    commanded voltage angle (integrator)
% gfma2    max P overload mitigation control state
% gfma3    min P overload mitigation control state
% gfma4    commanded voltage magnitude (first-order)
% gfma5    voltage regulation integral control state
% gfma6    max Q overload mitigation control state
% gfma7    min Q overload mitigation control state
% gfma8    inverter active power transducer
% gfma9    inverter reactive power transducer
% gfma10   inverter voltage transducer
%-----------------------------------------------------------------------------%

clear all; close all; clc;

load('./mat/gfma_example1_nonlinear_sim.mat');  % you may choose any valid mat file

Fs = round(1/median(diff(t)));
lbnd = 1e-4;
t_ax = [0,0.6*t(end)];  % time window for plotting (axes only)
y_stretch = 0.08;       % stretch the y-axis to make the plots easier to read

%-----------------------------------------------------------------------------%
% gfma1 -- commanded voltage angle (integrator)

vgfma1 = zeros(size(g.gfma.gfma1));
max_pct_err = -1;

for ii = 1:size(g.gfma.gfma_con,1)
    b = [0,1];                        % integrator
    a = [1,0];

    [bd,ad] = bilinear(b,a,Fs);

    U = g.gfma.dgfma1(ii,1)*ones(3,1);
    Y = g.gfma.gfma1(ii,1)*ones(3,1);
    Z = filtic(bd,ad,Y,U);
    vgfma1(ii,:) = filter(bd,ad,g.gfma.dgfma1(ii,:),Z);

    tmp_err = norm(vgfma1(ii,:)-g.gfma.gfma1(ii,:),'inf');
    tmp_pct_err = tmp_err;

    pct_den = max(abs(vgfma1(ii,:)-vgfma1(ii,1)));
    if (pct_den > 0)
        tmp_pct_err = 100*tmp_pct_err/pct_den;
    end

    if (tmp_pct_err > max_pct_err)
        gfma1_err = tmp_err;
        max_pct_err = tmp_pct_err;
    end
end

fprintf('\ngfma1_err: %0.2e, %0.2f%%',gfma1_err,max_pct_err);

fig{1} = figure;
ax{1}(1) = subplot(2,1,1,'parent',fig{1});
ax{1}(2) = subplot(2,1,2,'parent',fig{1});
hold(ax{1}(1),'on');
plot(ax{1}(1),t,vgfma1);
plot(ax{1}(1),t,g.gfma.gfma1,'--');
plot(ax{1}(2),t,vgfma1-g.gfma.gfma1);

for ii = 1:2
    v = axis(ax{1}(ii));
    y_span = v(4)-v(3);
    alpha = y_stretch*y_span;
    axis(ax{1}(ii),[t_ax(1),t_ax(2),v(3)-alpha,v(4)+alpha]);
end

xlabel(ax{1}(2), 'Time (s)');
ylabel(ax{1}(1), 'gfma1 (pu)');
ylabel(ax{1}(2), 'error (pu)');

%-----------------------------------------------------------------------------%
% gfma2 -- max P overload mitigation control state

vgfma2 = zeros(size(g.gfma.gfma1));
max_pct_err = -1;

for ii = 1:size(g.gfma.gfma_con,1)
    b = [0,1];                        % integrator
    a = [1,0];

    [bd,ad] = bilinear(b,a,Fs);

    U = g.gfma.dgfma2(ii,1)*ones(3,1);
    Y = g.gfma.gfma2(ii,1)*ones(3,1);
    Z = filtic(bd,ad,Y,U);
    vgfma2(ii,:) = filter(bd,ad,g.gfma.dgfma2(ii,:),Z);
    vgfma2(ii,:) = min(vgfma2(ii,:),0);

    tmp_err = norm(vgfma2(ii,:)-g.gfma.gfma2(ii,:),'inf');
    tmp_pct_err = tmp_err;

    pct_den = max(abs(vgfma2(ii,:)-vgfma2(ii,1)));
    if (pct_den > 0)
        tmp_pct_err = 100*tmp_pct_err/pct_den;
    end

    if (tmp_pct_err > max_pct_err)
        gfma2_err = tmp_err;
        max_pct_err = tmp_pct_err;
    end
end

fprintf('\ngfma2_err: %0.2e, %0.2f%%',gfma2_err,max_pct_err);

fig{2} = figure;
ax{2}(1) = subplot(2,1,1,'parent',fig{2});
ax{2}(2) = subplot(2,1,2,'parent',fig{2});
hold(ax{2}(1),'on');
plot(ax{2}(1),t,vgfma2);
plot(ax{2}(1),t,g.gfma.gfma2,'--');
plot(ax{2}(2),t,vgfma2-g.gfma.gfma2);

for ii = 1:2
    v = axis(ax{2}(ii));
    y_span = v(4)-v(3);
    alpha = y_stretch*y_span;
    axis(ax{2}(ii),[t_ax(1),t_ax(2),v(3)-alpha,v(4)+alpha]);
end

xlabel(ax{2}(2), 'Time (s)');
ylabel(ax{2}(1), 'gfma2 (pu)');
ylabel(ax{2}(2), 'error (pu)');

%-----------------------------------------------------------------------------%
% gfma3 -- min P overload mitigation control state

vgfma3 = zeros(size(g.gfma.gfma1));
max_pct_err = -1;

for ii = 1:size(g.gfma.gfma_con,1)
    b = [0,1];                        % integrator
    a = [1,0];

    [bd,ad] = bilinear(b,a,Fs);

    U = g.gfma.dgfma3(ii,1)*ones(3,1);
    Y = g.gfma.gfma3(ii,1)*ones(3,1);
    Z = filtic(bd,ad,Y,U);
    vgfma3(ii,:) = filter(bd,ad,g.gfma.dgfma3(ii,:),Z);
    vgfma3(ii,:) = max(vgfma3(ii,:),0);

    tmp_err = norm(vgfma3(ii,:)-g.gfma.gfma3(ii,:),'inf');
    tmp_pct_err = tmp_err;

    pct_den = max(abs(vgfma3(ii,:)-vgfma3(ii,1)));
    if (pct_den > 0)
        tmp_pct_err = 100*tmp_pct_err/pct_den;
    end

    if (tmp_pct_err > max_pct_err)
        gfma3_err = tmp_err;
        max_pct_err = tmp_pct_err;
    end
end

fprintf('\ngfma3_err: %0.2e, %0.2f%%',gfma3_err,max_pct_err);

fig{3} = figure;
ax{3}(1) = subplot(2,1,1,'parent',fig{3});
ax{3}(2) = subplot(2,1,2,'parent',fig{3});
hold(ax{3}(1),'on');
plot(ax{3}(1),t,vgfma3);
plot(ax{3}(1),t,g.gfma.gfma3,'--');
plot(ax{3}(2),t,vgfma3-g.gfma.gfma3);

for ii = 1:2
    v = axis(ax{3}(ii));
    y_span = v(4)-v(3);
    alpha = y_stretch*y_span;
    axis(ax{3}(ii),[t_ax(1),t_ax(2),v(3)-alpha,v(4)+alpha]);
end

xlabel(ax{3}(2), 'Time (s)');
ylabel(ax{3}(1), 'gfma3 (pu)');
ylabel(ax{3}(2), 'error (pu)');

%-----------------------------------------------------------------------------%
% gfma4 -- commanded voltage magnitude (first-order)

vgfma4 = zeros(size(g.gfma.gfma1));
max_pct_err = -1;

for ii = 1:size(g.gfma.gfma_con,1)
    if (g.gfma.gfma_con(ii,14) > lbnd)
        b = [0,1];                        % lowpass filter
        a = [g.gfma.gfma_con(ii,14),1];

        [bd,ad] = bilinear(b,a,Fs);
    else
        bd = 1;
        ad = 1;
    end

    Q1 = g.gfma.gfma_con(ii,9).*(g.gfma.qset(ii,:) - g.gfma.gfma9(ii,:) ...
                                + imag(g.gfma.gfma_sig(ii,:)));

    Q2 = g.gfma.gfma_con(ii,15).*(g.gfma.gfma_con(ii,10) - g.gfma.gfma9(ii,:)) ...
         + g.gfma.gfma6(ii,:);

    Q2 = min(Q2,0);

    Q3 = g.gfma.gfma_con(ii,15).*(g.gfma.gfma_con(ii,11) - g.gfma.gfma9(ii,:)) ...
         + g.gfma.gfma7(ii,:);

    Q3 = max(Q3,0);

    vref = g.gfma.vset(ii,:) + Q1 + Q2 + Q3;

    e_cmd = vref;
    if (g.gfma.gfma_con(ii,25) == 1)
        e_cmd = g.gfma.gfma_con(ii,12) ...
                .*(vref - g.gfma.gfma10(ii,:)) + g.gfma.gfma5(ii,:);
    end

    U = e_cmd(1)*ones(3,1);
    Y = g.gfma.gfma4(ii,1)*ones(3,1);
    Z = filtic(bd,ad,Y,U);
    vgfma4(ii,:) = filter(bd,ad,e_cmd,Z);
    vgfma4(ii,:) = ...
        max(min(vgfma4(ii,:),g.gfma.gfma_con(ii,18)),g.gfma.gfma_con(ii,19));

    tmp_err = norm(vgfma4(ii,:)-g.gfma.gfma4(ii,:),'inf');
    tmp_pct_err = tmp_err;

    pct_den = max(abs(vgfma4(ii,:)-vgfma4(ii,1)));
    if (pct_den > 0)
        tmp_pct_err = 100*tmp_pct_err/pct_den;
    end

    if (tmp_pct_err > max_pct_err)
        gfma4_err = tmp_err;
        max_pct_err = tmp_pct_err;
    end
end

fprintf('\ngfma4_err: %0.2e, %0.2f%%',gfma4_err,max_pct_err);

fig{4} = figure;
ax{4}(1) = subplot(2,1,1,'parent',fig{4});
ax{4}(2) = subplot(2,1,2,'parent',fig{4});
hold(ax{4}(1),'on');
plot(ax{4}(1),t,vgfma4);
plot(ax{4}(1),t,g.gfma.gfma4,'--');
plot(ax{4}(2),t,vgfma4-g.gfma.gfma4);

for ii = 1:2
    v = axis(ax{4}(ii));
    y_span = v(4)-v(3);
    alpha = y_stretch*y_span;
    axis(ax{4}(ii),[t_ax(1),t_ax(2),v(3)-alpha,v(4)+alpha]);
end

xlabel(ax{4}(2), 'Time (s)');
ylabel(ax{4}(1), 'gfma4 (pu)');
ylabel(ax{4}(2), 'error (pu)');

%-----------------------------------------------------------------------------%
% gfma5 -- voltage regulation integral control state

vgfma5 = zeros(size(g.gfma.gfma1));
max_pct_err = -1;

for ii = 1:size(g.gfma.gfma_con,1)
    b = [0,1];                        % integrator
    a = [1,0];

    [bd,ad] = bilinear(b,a,Fs);

    Q1 = g.gfma.gfma_con(ii,9).*(g.gfma.qset(ii,:) - g.gfma.gfma9(ii,:) ...
                                + imag(g.gfma.gfma_sig(ii,:)));

    Q2 = g.gfma.gfma_con(ii,15).*(g.gfma.gfma_con(ii,10) - g.gfma.gfma9(ii,:)) ...
         + g.gfma.gfma6(ii,:);

    Q2 = min(Q2,0);

    Q3 = g.gfma.gfma_con(ii,15).*(g.gfma.gfma_con(ii,11) - g.gfma.gfma9(ii,:)) ...
         + g.gfma.gfma7(ii,:);

    Q3 = max(Q3,0);

    vref = g.gfma.vset(ii,:) + Q1 + Q2 + Q3;

    e_cmd = vref;
    if (g.gfma.gfma_con(ii,25) == 1)
        e_cmd = g.gfma.gfma_con(ii,12) ...
                .*(vref - g.gfma.gfma10(ii,:)) + g.gfma.gfma5(ii,:);
    end

    u_int = g.gfma.gfma_con(ii,13).*(vref - g.gfma.gfma10(ii,:));

    U = u_int(1)*ones(3,1);
    Y = g.gfma.gfma5(ii,1)*ones(3,1);
    Z = filtic(bd,ad,Y,U);
    vgfma5(ii,:) = filter(bd,ad,u_int,Z);
    vgfma5(ii,:) = ...
        max(min(vgfma5(ii,:),g.gfma.gfma_con(ii,18)),g.gfma.gfma_con(ii,19));

    tmp_err = norm(vgfma5(ii,:)-g.gfma.gfma5(ii,:),'inf');
    tmp_pct_err = tmp_err;

    pct_den = max(abs(vgfma5(ii,:)-vgfma5(ii,1)));
    if (pct_den > 0)
        tmp_pct_err = 100*tmp_pct_err/pct_den;
    end

    if (tmp_pct_err > max_pct_err)
        gfma5_err = tmp_err;
        max_pct_err = tmp_pct_err;
    end
end

fprintf('\ngfma5_err: %0.2e, %0.2f%%',gfma5_err,max_pct_err);

fig{5} = figure;
ax{5}(1) = subplot(2,1,1,'parent',fig{5});
ax{5}(2) = subplot(2,1,2,'parent',fig{5});
hold(ax{5}(1),'on');
plot(ax{5}(1),t,vgfma5);
plot(ax{5}(1),t,g.gfma.gfma5,'--');
plot(ax{5}(2),t,vgfma5-g.gfma.gfma5);

for ii = 1:2
    v = axis(ax{5}(ii));
    y_span = v(4)-v(3);
    alpha = y_stretch*y_span;
    axis(ax{5}(ii),[t_ax(1),t_ax(2),v(3)-alpha,v(4)+alpha]);
end

xlabel(ax{5}(2), 'Time (s)');
ylabel(ax{5}(1), 'gfma5 (pu)');
ylabel(ax{5}(2), 'error (pu)');

%-----------------------------------------------------------------------------%
% gfma6 -- max Q overload mitigation control state

vgfma6 = zeros(size(g.gfma.gfma1));
max_pct_err = -1;

for ii = 1:size(g.gfma.gfma_con,1)
    b = [0,1];                        % integrator
    a = [1,0];

    [bd,ad] = bilinear(b,a,Fs);

    U = g.gfma.dgfma6(ii,1)*ones(3,1);
    Y = g.gfma.gfma6(ii,1)*ones(3,1);
    Z = filtic(bd,ad,Y,U);
    vgfma6(ii,:) = filter(bd,ad,g.gfma.dgfma6(ii,:),Z);
    vgfma6(ii,:) = min(vgfma6(ii,:),0);

    tmp_err = norm(vgfma6(ii,:)-g.gfma.gfma6(ii,:),'inf');
    tmp_pct_err = tmp_err;

    pct_den = max(abs(vgfma6(ii,:)-vgfma6(ii,1)));
    if (pct_den > 0)
        tmp_pct_err = 100*tmp_pct_err/pct_den;
    end

    if (tmp_pct_err > max_pct_err)
        gfma6_err = tmp_err;
        max_pct_err = tmp_pct_err;
    end
end

fprintf('\ngfma6_err: %0.2e, %0.2f%%',gfma6_err,max_pct_err);

fig{6} = figure;
ax{6}(1) = subplot(2,1,1,'parent',fig{6});
ax{6}(2) = subplot(2,1,2,'parent',fig{6});
hold(ax{6}(1),'on');
plot(ax{6}(1),t,vgfma6);
plot(ax{6}(1),t,g.gfma.gfma6,'--');
plot(ax{6}(2),t,vgfma6-g.gfma.gfma6);

for ii = 1:2
    v = axis(ax{6}(ii));
    y_span = v(4)-v(3);
    alpha = y_stretch*y_span;
    axis(ax{6}(ii),[t_ax(1),t_ax(2),v(3)-alpha,v(4)+alpha]);
end

xlabel(ax{6}(2), 'Time (s)');
ylabel(ax{6}(1), 'gfma6 (pu)');
ylabel(ax{6}(2), 'error (pu)');

%-----------------------------------------------------------------------------%
% gfma7 -- min Q overload mitigation control state

vgfma7 = zeros(size(g.gfma.gfma1));
max_pct_err = -1;

for ii = 1:size(g.gfma.gfma_con,1)
    b = [0,1];                        % integrator
    a = [1,0];

    [bd,ad] = bilinear(b,a,Fs);

    U = g.gfma.dgfma7(ii,1)*ones(3,1);
    Y = g.gfma.gfma7(ii,1)*ones(3,1);
    Z = filtic(bd,ad,Y,U);
    vgfma7(ii,:) = filter(bd,ad,g.gfma.dgfma7(ii,:),Z);
    vgfma7(ii,:) = max(vgfma7(ii,:),0);

    tmp_err = norm(vgfma7(ii,:)-g.gfma.gfma7(ii,:),'inf');
    tmp_pct_err = tmp_err;

    pct_den = max(abs(vgfma7(ii,:)-vgfma7(ii,1)));
    if (pct_den > 0)
        tmp_pct_err = 100*tmp_pct_err/pct_den;
    end

    if (tmp_pct_err > max_pct_err)
        gfma7_err = tmp_err;
        max_pct_err = tmp_pct_err;
    end
end

fprintf('\ngfma7_err: %0.2e, %0.2f%%',gfma7_err,max_pct_err);

fig{7} = figure;
ax{7}(1) = subplot(2,1,1,'parent',fig{7});
ax{7}(2) = subplot(2,1,2,'parent',fig{7});
hold(ax{7}(1),'on');
plot(ax{7}(1),t,vgfma7);
plot(ax{7}(1),t,g.gfma.gfma7,'--');
plot(ax{7}(2),t,vgfma7-g.gfma.gfma7);

for ii = 1:2
    v = axis(ax{7}(ii));
    y_span = v(4)-v(3);
    alpha = y_stretch*y_span;
    axis(ax{7}(ii),[t_ax(1),t_ax(2),v(3)-alpha,v(4)+alpha]);
end

xlabel(ax{7}(2), 'Time (s)');
ylabel(ax{7}(1), 'gfma7 (pu)');
ylabel(ax{7}(2), 'error (pu)');

%-----------------------------------------------------------------------------%
% gfma8 -- inverter active power transducer

vgfma8 = zeros(size(g.gfma.gfma1));
max_pct_err = -1;

for ii = 1:size(g.gfma.gfma_con,1)
    if (g.gfma.gfma_con(ii,21) > lbnd)
        b = [0,1];                        % lowpass filter
        a = [g.gfma.gfma_con(ii,21),1];

        [bd,ad] = bilinear(b,a,Fs);
    else
        bd = 1;
        ad = 1;
    end

    P_inv = g.mac.pelect(g.gfma.mac_idx(ii),:) ...
            .*g.mac.mac_pot(g.gfma.mac_idx(ii),1);

    U = P_inv(1)*ones(3,1);
    Y = g.gfma.gfma8(ii,1)*ones(3,1);
    Z = filtic(bd,ad,Y,U);
    vgfma8(ii,:) = filter(bd,ad,P_inv,Z);

    tmp_err = norm(vgfma8(ii,:)-g.gfma.gfma8(ii,:),'inf');
    tmp_pct_err = tmp_err;

    pct_den = max(abs(vgfma8(ii,:)-vgfma8(ii,1)));
    if (pct_den > 0)
        tmp_pct_err = 100*tmp_pct_err/pct_den;
    end

    if (tmp_pct_err > max_pct_err)
        gfma8_err = tmp_err;
        max_pct_err = tmp_pct_err;
    end
end

fprintf('\ngfma8_err: %0.2e, %0.2f%%',gfma8_err,max_pct_err);

fig{8} = figure;
ax{8}(1) = subplot(2,1,1,'parent',fig{8});
ax{8}(2) = subplot(2,1,2,'parent',fig{8});
hold(ax{8}(1),'on');
plot(ax{8}(1),t,vgfma8);
plot(ax{8}(1),t,g.gfma.gfma8,'--');
plot(ax{8}(2),t,vgfma8-g.gfma.gfma8);

for ii = 1:2
    v = axis(ax{8}(ii));
    y_span = v(4)-v(3);
    alpha = y_stretch*y_span;
    axis(ax{8}(ii),[t_ax(1),t_ax(2),v(3)-alpha,v(4)+alpha]);
end

xlabel(ax{8}(2), 'Time (s)');
ylabel(ax{8}(1), 'gfma8 (pu)');
ylabel(ax{8}(2), 'error (pu)');

%-----------------------------------------------------------------------------%
% gfma9 -- inverter reactive power transducer

vgfma9 = zeros(size(g.gfma.gfma1));
max_pct_err = -1;

for ii = 1:size(g.gfma.gfma_con,1)
    if (g.gfma.gfma_con(ii,22) > lbnd)
        b = [0,1];                        % lowpass filter
        a = [g.gfma.gfma_con(ii,22),1];

        [bd,ad] = bilinear(b,a,Fs);
    else
        bd = 1;
        ad = 1;
    end

    Q_inv = g.mac.qelect(g.gfma.mac_idx(ii),:) ...
            .*g.mac.mac_pot(g.gfma.mac_idx(ii),1);

    U = Q_inv(1)*ones(3,1);
    Y = g.gfma.gfma9(ii,1)*ones(3,1);
    Z = filtic(bd,ad,Y,U);
    vgfma9(ii,:) = filter(bd,ad,Q_inv,Z);

    tmp_err = norm(vgfma9(ii,:)-g.gfma.gfma9(ii,:),'inf');
    tmp_pct_err = tmp_err;

    pct_den = max(abs(vgfma9(ii,:)-vgfma9(ii,1)));
    if (pct_den > 0)
        tmp_pct_err = 100*tmp_pct_err/pct_den;
    end

    if (tmp_pct_err > max_pct_err)
        gfma9_err = tmp_err;
        max_pct_err = tmp_pct_err;
    end
end

fprintf('\ngfma9_err: %0.2e, %0.2f%%',gfma9_err,max_pct_err);

fig{9} = figure;
ax{9}(1) = subplot(2,1,1,'parent',fig{9});
ax{9}(2) = subplot(2,1,2,'parent',fig{9});
hold(ax{9}(1),'on');
plot(ax{9}(1),t,vgfma9);
plot(ax{9}(1),t,g.gfma.gfma9,'--');
plot(ax{9}(2),t,vgfma9-g.gfma.gfma9);

for ii = 1:2
    v = axis(ax{9}(ii));
    y_span = v(4)-v(3);
    alpha = y_stretch*y_span;
    axis(ax{9}(ii),[t_ax(1),t_ax(2),v(3)-alpha,v(4)+alpha]);
end

xlabel(ax{9}(2), 'Time (s)');
ylabel(ax{9}(1), 'gfma9 (pu)');
ylabel(ax{9}(2), 'error (pu)');

%-----------------------------------------------------------------------------%
% gfma10 -- inverter voltage transducer

vgfma10 = zeros(size(g.gfma.gfma1));
max_pct_err = -1;

busnum = g.bus.bus_int(g.gfma.gfma_con(:,2));
for ii = 1:size(g.gfma.gfma_con,1)
    if (g.gfma.gfma_con(ii,23) > lbnd)
        b = [0,1];                        % lowpass filter
        a = [g.gfma.gfma_con(ii,23),1];

        [bd,ad] = bilinear(b,a,Fs);
    else
        bd = 1;
        ad = 1;
    end

    Vt = abs(g.bus.bus_v(busnum(ii),:));

    U = real(Vt(1))*ones(3,1);
    Y = g.gfma.gfma10(ii,1)*ones(3,1);
    Z = filtic(bd,ad,Y,U);
    vgfma10(ii,:) = filter(bd,ad,Vt,Z);

    tmp_err = norm(vgfma10(ii,:)-g.gfma.gfma10(ii,:),'inf');
    tmp_pct_err = tmp_err;

    pct_den = max(abs(vgfma10(ii,:)-vgfma10(ii,1)));
    if (pct_den > 0)
        tmp_pct_err = 100*tmp_pct_err/pct_den;
    end

    if (tmp_pct_err > max_pct_err)
        gfma10_err = tmp_err;
        max_pct_err = tmp_pct_err;
    end
end

fprintf('\ngfma10_err: %0.2e, %0.2f%%',gfma10_err,max_pct_err);

fig{10} = figure;
ax{10}(1) = subplot(2,1,1,'parent',fig{10});
ax{10}(2) = subplot(2,1,2,'parent',fig{10});
hold(ax{10}(1),'on');
plot(ax{10}(1),t,vgfma10);
plot(ax{10}(1),t,g.gfma.gfma10,'--');
plot(ax{10}(2),t,vgfma10-g.gfma.gfma10);

for ii = 1:2
    v = axis(ax{10}(ii));
    y_span = v(4)-v(3);
    alpha = y_stretch*y_span;
    axis(ax{10}(ii),[t_ax(1),t_ax(2),v(3)-alpha,v(4)+alpha]);
end

xlabel(ax{10}(2), 'Time (s)');
ylabel(ax{10}(1), 'gfma10 (pu)');
ylabel(ax{10}(2), 'error (pu)');

%-----------------------------------------------------------------------------%
% eqprime -- ivm voltage magnitude

vivm1 = zeros(size(g.mac.eqprime(g.gfma.mac_idx,:)));
max_pct_err = -1;

for ii = 1:size(g.ivm.ivm_con,1)
    if (g.mac.mac_con(g.gfma.mac_idx(ii),10) > lbnd)
        b = [0,1];
        a = [g.mac.mac_con(g.gfma.mac_idx(ii),10),1];

        [bd,ad] = bilinear(b,a,Fs);
    else
        bd = 1;
        ad = 1;
    end

    e_ord = g.mac.vex(g.gfma.mac_idx(ii),:);

    % e_ord = g.mac.vex(g.gfma.mac_idx(ii),:) ...
    %         + abs(g.ivm.ivm_sig(g.gfma.ivm_idx(ii),:));

    U = e_ord(1)*ones(3,1);
    Y = g.mac.eqprime(g.gfma.mac_idx(ii),1)*ones(3,1);
    Z = filtic(bd,ad,Y,U);
    vivm1(ii,:) = filter(bd,ad,e_ord,Z);

    tmp_err = norm(vivm1(ii,:)-g.mac.eqprime(g.gfma.mac_idx(ii),:),'inf');
    tmp_pct_err = tmp_err;

    pct_den = max(abs(vivm1(ii,:)-vivm1(ii,1)));
    if (pct_den > 0)
        tmp_pct_err = 100*tmp_pct_err/pct_den;
    end

    if (tmp_pct_err > max_pct_err)
        ivm1_err = tmp_err;
        max_pct_err = tmp_pct_err;
    end
end

fprintf('\n\nivm1_err: %0.2e, %0.2f%%',ivm1_err,max_pct_err);

fig{11} = figure;
ax{11}(1) = subplot(2,1,1,'parent',fig{11});
ax{11}(2) = subplot(2,1,2,'parent',fig{11});
hold(ax{11}(1),'on');
plot(ax{11}(1),t,vivm1);
plot(ax{11}(1),t,g.mac.eqprime(g.gfma.mac_idx,:),'--');
plot(ax{11}(2),t,vivm1-g.mac.eqprime(g.gfma.mac_idx,:));

for ii = 1:2
    v = axis(ax{11}(ii));
    y_span = v(4)-v(3);
    alpha = y_stretch*y_span;
    axis(ax{11}(ii),[t_ax(1),t_ax(2),v(3)-alpha,v(4)+alpha]);
end

xlabel(ax{11}(2), 'Time (s)');
ylabel(ax{11}(1), 'ivm eqprime (pu)');
ylabel(ax{11}(2), 'error (pu)');

%-----------------------------------------------------------------------------%
% mac_ang -- ivm voltage angle

vivm2 = zeros(size(g.mac.mac_ang(g.gfma.mac_idx,:)));
max_pct_err = -1;

for ii = 1:size(g.ivm.ivm_con,1)
    if (g.mac.mac_con(g.gfma.mac_idx(ii),9) > lbnd)
        b = [0,1];
        a = [g.mac.mac_con(g.gfma.mac_idx(ii),9),1];

        [bd,ad] = bilinear(b,a,Fs);
    else
        bd = 1;
        ad = 1;
    end

    d_ord = g.mac.fldcur(g.gfma.mac_idx(ii),:);

    % d_ord = g.mac.fldcur(g.gfma.mac_idx(ii),:) ...
    %         + angle(g.ivm.ivm_sig(g.gfma.ivm_idx(ii),:));

    U = d_ord(1)*ones(3,1);
    Y = g.mac.mac_ang(g.gfma.mac_idx(ii),1)*ones(3,1);
    Z = filtic(bd,ad,Y,U);
    vivm2(ii,:) = filter(bd,ad,d_ord,Z);

    tmp_err = norm(vivm2(ii,:)-g.mac.mac_ang(g.gfma.mac_idx(ii),:),'inf');
    tmp_pct_err = tmp_err;

    pct_den = max(abs(vivm2(ii,:)-vivm2(ii,1)));
    if (pct_den > 0)
        tmp_pct_err = 100*tmp_pct_err/pct_den;
    end

    if (tmp_pct_err > max_pct_err)
        ivm2_err = tmp_err;
        max_pct_err = tmp_pct_err;
    end
end

fprintf('\nivm2_err: %0.2e, %0.2f%%',ivm2_err,max_pct_err);

fig{12} = figure;
ax{12}(1) = subplot(2,1,1,'parent',fig{12});
ax{12}(2) = subplot(2,1,2,'parent',fig{12});
hold(ax{12}(1),'on');
plot(ax{12}(1),t,vivm2);
plot(ax{12}(1),t,g.mac.mac_ang(g.gfma.mac_idx,:),'--');
plot(ax{12}(2),t,vivm2-g.mac.mac_ang(g.gfma.mac_idx,:));

for ii = 1:2
    v = axis(ax{12}(ii));
    y_span = v(4)-v(3);
    alpha = y_stretch*y_span;
    axis(ax{12}(ii),[t_ax(1),t_ax(2),v(3)-alpha,v(4)+alpha]);
end

xlabel(ax{12}(2), 'Time (s)');
ylabel(ax{12}(1), 'ivm mac\_ang (pu)');
ylabel(ax{12}(2), 'error (pu)');

fprintf('\n\n');

% eof
