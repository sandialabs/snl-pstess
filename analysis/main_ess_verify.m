%-----------------------------------------------------------------------------%
% ess model verification script
%
% Purpose: This is a script for verifying that the differential equations
%          encoded in the ess model match those in the block diagram. It
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
% ess1  voltage transducer
% ess2  Pade approximation
% ess3  active current converter interface
% ess4  reactive current converter interface
% ess5  ess SOC integrator
%-----------------------------------------------------------------------------%

clear all; close all; clc;

load('./mat/ess_example1_nonlinear_sim.mat');  % you may choose any valid mat file
% load('./mat/ess_example2_nonlinear_sim.mat');

Fs = round(1/median(diff(t)));
lbnd = 1e-3;
t_ax = [0,0.6*t(end)];  % time window for plotting (axes only)
y_stretch = 0.08;       % stretch the y-axis to make the plots easier to read

%-----------------------------------------------------------------------------%
% ess1 -- voltage transducer

vess1 = zeros(size(g.ess.ess1));
max_pct_err = -1;

for ii = 1:size(g.ess.ess_con,1)
    if (g.ess.ess_con(ii,3) > lbnd)
        b = [0,1];                        % lowpass filter
        a = [g.ess.ess_con(ii,3),1];

        [bd,ad] = bilinear(b,a,Fs);
    else
        bd = 1;
        ad = 1;
    end

    U = g.ess.ess_vmag(ii,1)*ones(3,1);
    Y = g.ess.ess_vmag(ii,1)*ones(3,1);
    Z = filtic(bd,ad,Y,U);
    vess1(ii,:) = filter(bd,ad,g.ess.ess_vmag(ii,:),Z);

    tmp_err = norm(vess1(ii,:)-g.ess.ess1(ii,:),'inf');
    tmp_pct_err = tmp_err;

    pct_den = max(abs(vess1(ii,:)-vess1(ii,1)));
    if (pct_den > 0)
        tmp_pct_err = 100*tmp_pct_err/pct_den;
    end

    if (tmp_pct_err > max_pct_err)
        ess1_err = tmp_err;
        max_pct_err = tmp_pct_err;
    end
end

fprintf('\ness1_err: %0.2e, %0.2f%%',ess1_err,max_pct_err);

fig{1} = figure;
ax{1}(1) = subplot(2,1,1,'parent',fig{1});
ax{1}(2) = subplot(2,1,2,'parent',fig{1});
hold(ax{1}(1),'on');
plot(ax{1}(1),t,vess1);
plot(ax{1}(1),t,g.ess.ess1,'--');
plot(ax{1}(2),t,vess1-g.ess.ess1);

for ii = 1:2
    v = axis(ax{1}(ii));
    y_span = v(4)-v(3);
    alpha = y_stretch*y_span;
    axis(ax{1}(ii),[t_ax(1),t_ax(2),v(3)-alpha,v(4)+alpha]);
end

xlabel(ax{1}(2), 'Time (s)');
ylabel(ax{1}(1), 'ess1 (pu)');
ylabel(ax{1}(2), 'error (pu)');

%-----------------------------------------------------------------------------%
% ess2 -- local voltage magnitude Pade approximation

vess2 = zeros(size(g.ess.ess1));
max_pct_err = -1;

for ii = 1:size(g.ess.ess_con,1)
    if (g.ess.ess_con(ii,5)/2 >= lbnd)
        b = [-g.ess.ess_con(ii,5)/2,1];   % Pade approximation
        a = [g.ess.ess_con(ii,5)/2,1];

        [bd,ad] = bilinear(b,a,Fs);
    else
        bd = 1;
        ad = 1;
    end

    U = g.ess.ess1(ii,1)*ones(3,1);
    Y = g.ess.ess1(ii,1)*ones(3,1);
    Z = filtic(bd,ad,Y,U);
    vess2(ii,:) = filter(bd,ad,g.ess.ess1(ii,:),Z);

    tmp_err = norm(vess2(ii,:)-g.ess.ess_vmag_pade(ii,:),'inf');
    tmp_pct_err = tmp_err;

    pct_den = max(abs(vess2(ii,:)-vess2(ii,1)));
    if (pct_den > 0)
        tmp_pct_err = 100*tmp_pct_err/pct_den;
    end

    if (tmp_pct_err > max_pct_err)
        ess2_err = tmp_err;
        max_pct_err = tmp_pct_err;
    end
end

fprintf('\ness2_err: %0.2e, %0.2f%%',ess2_err,max_pct_err);

fig{2} = figure;
ax{2}(1) = subplot(2,1,1,'parent',fig{2});
ax{2}(2) = subplot(2,1,2,'parent',fig{2});
hold(ax{2}(1),'on');
plot(ax{2}(1),t,vess2);
plot(ax{2}(1),t,g.ess.ess_vmag_pade,'--');
plot(ax{2}(2),t,vess2-g.ess.ess_vmag_pade);

for ii = 1:2
    v = axis(ax{2}(ii));
    y_span = v(4)-v(3);
    alpha = y_stretch*y_span;
    axis(ax{2}(ii),[t_ax(1),t_ax(2),v(3)-alpha,v(4)+alpha]);
end

xlabel(ax{2}(2), 'Time (s)');
ylabel(ax{2}(1), 'ess2 (pu)');
ylabel(ax{2}(2), 'error (pu)');

%-----------------------------------------------------------------------------%
% ess3 -- active current converter interface

vess3 = zeros(size(g.ess.ess1));
max_pct_err = -1;

for ii = 1:size(g.ess.ess_con,1)
    if (g.ess.ess_con(ii,14) > lbnd)
        b = [0,1];                        % lowpass filter
        a = [g.ess.ess_con(ii,14),1];

        [bd,ad] = bilinear(b,a,Fs);
    else
        bd = 1;
        ad = 1;
    end

    ip_ord(ii,:) = g.ess.dess3(ii,:).*max(g.ess.ess_con(ii,14),lbnd) ...
                   + g.ess.ess3(ii,:);

    U = ip_ord(ii,1)*ones(3,1);
    Y = ip_ord(ii,1)*ones(3,1);
    Z = filtic(bd,ad,Y,U);
    vess3(ii,:) = filter(bd,ad,ip_ord(ii,:),Z);

    tmp_err = norm(vess3(ii,:)-g.ess.ess3(ii,:),'inf');
    tmp_pct_err = tmp_err;

    pct_den = max(abs(vess3(ii,:)-vess3(ii,1)));
    if (pct_den > 0)
        tmp_pct_err = 100*tmp_pct_err/pct_den;
    end

    if (tmp_pct_err > max_pct_err)
        ess3_err = tmp_err;
        max_pct_err = tmp_pct_err;
    end
end

fprintf('\ness3_err: %0.2e, %0.2f%%',ess3_err,max_pct_err);

fig{3} = figure;
ax{3}(1) = subplot(2,1,1,'parent',fig{3});
ax{3}(2) = subplot(2,1,2,'parent',fig{3});
hold(ax{3}(1),'on');
plot(ax{3}(1),t,vess3);
plot(ax{3}(1),t,g.ess.ess3,'--');
plot(ax{3}(2),t,vess3-g.ess.ess3);

for ii = 1:2
    v = axis(ax{3}(ii));
    y_span = v(4)-v(3);
    alpha = y_stretch*y_span;
    axis(ax{3}(ii),[t_ax(1),t_ax(2),v(3)-alpha,v(4)+alpha]);
end

xlabel(ax{3}(2), 'Time (s)');
ylabel(ax{3}(1), 'ess3 (pu)');
ylabel(ax{3}(2), 'error (pu)');

%-----------------------------------------------------------------------------%
% ess4 -- reactive current converter interface

vess4 = zeros(size(g.ess.ess1));
max_pct_err = -1;

for ii = 1:size(g.ess.ess_con,1)
    if (g.ess.ess_con(ii,14) > lbnd)
        b = [0,1];                        % lowpass filter
        a = [g.ess.ess_con(ii,14),1];

        [bd,ad] = bilinear(b,a,Fs);
    else
        bd = 1;
        ad = 1;
    end

    iq_ord(ii,:) = g.ess.dess4(ii,:).*max(g.ess.ess_con(ii,14),lbnd) ...
                   + g.ess.ess4(ii,:);

    U = iq_ord(ii,1)*ones(3,1);
    Y = iq_ord(ii,1)*ones(3,1);
    Z = filtic(bd,ad,Y,U);
    vess4(ii,:) = filter(bd,ad,iq_ord(ii,:),Z);

    tmp_err = norm(vess4(ii,:)-g.ess.ess4(ii,:),'inf');
    tmp_pct_err = tmp_err;

    pct_den = max(abs(vess4(ii,:)-vess4(ii,1)));
    if (pct_den > 0)
        tmp_pct_err = 100*tmp_pct_err/pct_den;
    end

    if (tmp_pct_err > max_pct_err)
        ess4_err = tmp_err;
        max_pct_err = tmp_pct_err;
    end
end

fprintf('\ness4_err: %0.2e, %0.2f%%',ess4_err,max_pct_err);

fig{4} = figure;
ax{4}(1) = subplot(2,1,1,'parent',fig{4});
ax{4}(2) = subplot(2,1,2,'parent',fig{4});
hold(ax{4}(1),'on');
plot(ax{4}(1),t,vess4);
plot(ax{4}(1),t,g.ess.ess4,'--');
plot(ax{4}(2),t,vess4-g.ess.ess4);

for ii = 1:2
    v = axis(ax{4}(ii));
    y_span = v(4)-v(3);
    alpha = y_stretch*y_span;
    axis(ax{4}(ii),[t_ax(1),t_ax(2),v(3)-alpha,v(4)+alpha]);
end

xlabel(ax{4}(2), 'Time (s)');
ylabel(ax{4}(1), 'ess4 (pu)');
ylabel(ax{4}(2), 'error (pu)');

%-----------------------------------------------------------------------------%
% ess5 -- ess SOC integrator

vess5 = zeros(size(g.ess.ess1));
max_pct_err = -1;

for ii = 1:size(g.ess.ess_con,1)
    b = [0,1];                            % integrator
    a = [1,0];

    [bd,ad] = bilinear(b,a,Fs);

    p_ch = -min(real(g.ess.ess_sinj(ii,:)),0);
    p_dis = max(real(g.ess.ess_sinj(ii,:)),0);

    if (g.ess.ess_con(ii,21) < lbnd)
        error('main_ess_verify: the efficiency must be nonzero.');
    end

    tmp_ch = g.ess.ess_pot(ii,2).*(p_ch.*g.ess.ess_con(ii,21) ...
                                   - p_dis./g.ess.ess_con(ii,21));

    U = 0*ones(3,1);
    Y = 0*ones(3,1);
    Z = filtic(bd,ad,Y,U);
    vess5(ii,:) = filter(bd,ad,tmp_ch) + g.ess.ess_soc(ii,1);

    tmp_err = norm(vess5(ii,:)-g.ess.ess_soc(ii,:),'inf');
    tmp_pct_err = tmp_err;

    pct_den = max(abs(vess5(ii,:)-vess5(ii,1)));
    if (pct_den > 0)
        tmp_pct_err = 100*tmp_pct_err/pct_den;
    end

    if (tmp_pct_err > max_pct_err)
        ess5_err = tmp_err;
        max_pct_err = tmp_pct_err;
    end
end

fprintf('\ness5_err: %0.2e, %0.2f%%',ess5_err,max_pct_err);

fig{5} = figure;
ax{5}(1) = subplot(2,1,1,'parent',fig{5});
ax{5}(2) = subplot(2,1,2,'parent',fig{5});
hold(ax{5}(1),'on');
plot(ax{5}(1),t,vess5);
plot(ax{5}(1),t,g.ess.ess_soc,'--');
plot(ax{5}(2),t,vess5-g.ess.ess_soc);

for ii = 1:2
    v = axis(ax{5}(ii));
    y_span = v(4)-v(3);
    alpha = y_stretch*y_span;
    axis(ax{5}(ii),[t_ax(1),t_ax(2),v(3)-alpha,v(4)+alpha]);
end

xlabel(ax{5}(2), 'Time (s)');
ylabel(ax{5}(1), 'ess5 (pu)');
ylabel(ax{5}(2), 'error (pu)');

fprintf('\n\n');

% eof
