% reec active modulation example
%
% Purpose: Example of how to use reec.m to modulate the control 
%          references of a grid-follinwg inverter. This example
%          compares the output of the nonlinear and linearized models
%          for a case where the perturbation comes from reec_sig
%          (set in mreec_sig).
%
% Data file = d2asbegp_reec.m
%
% Note: You must run this example from the `analysis' folder. Before
%       running, rename the original /pstess/mreec_sig.m file to
%       mreec_sig_original.m. Then move mreec_sig_example1.m into /pstess
%       and rename it mreec_sig.m. Make sure the working directory is
%       configured in get_path.m. Before running double check that sw_con
%       in d2asbegp_reec.m corresponds to ftype=6 (do nothing).
%-----------------------------------------------------------------------------%

clear all; close all; clc;
run('plot_options');

%% Run nonlinear simulation and store results
cd('../pstess/');
set_path('d2asbegp_reec.m');  % set the working directory in get_path

run('s_simu');  % run PSTess
bus_v = g.bus.bus_v;
t_step = 1/round(1/median(diff(t)));
t_end = round(t(end));

% save t and bus_v results
save('../analysis/mat/reec_example1_nonlinear_sim','t','bus_v','g');

%% Build linear model, simulate, and store results

clearvars -except t_step t_end; close all; clc;

run('svm_mgen');  % conduct linearization
save('../analysis/mat/reec_example1_linear','a_mat','b_reec_v','b_reec_q','c_v','c_ang');

%% Simulate linear model

tL = 0:t_step:t_end;

u1 = zeros(size(b_reec_q,2),length(tL));
u2 = zeros(size(b_reec_v,2),length(tL));
x = zeros(size(a_mat,1),length(tL));
dx = zeros(size(x));

y_v = zeros(size(c_v,1),length(tL));
y_ang = zeros(size(c_ang,1),length(tL));
y_qinv = zeros(size(c_reec_q,1),length(tL));
y_vloc = zeros(size(c_reec_v,1),length(tL));

u1_mask = (tL >= 1 & tL < 1.5);
u1(1,u1_mask) = -3e-3;

u1_mask = (tL >= 4.5 & tL < 5.0);
u2(2,u1_mask) = 5e-3;

for ii = 1:length(tL)
    dx(:,ii) = a_mat*x(:,ii) + b_reec_q*u1(:,ii) + b_reec_v*u2(:,ii);
    y_v(:,ii) = c_v*x(:,ii);
    y_ang(:,ii) = c_ang*x(:,ii);
    y_qinv(:,ii) = c_reec_q*x(:,ii);
    y_vloc(:,ii) = c_reec_v*x(:,ii);

    if (ii < length(tL))
        x(:,ii+1) = x(:,ii) + dx(:,ii)*t_step;
    end
end

load('../analysis/mat/reec_example1_nonlinear_sim','t','bus_v');
y_v = y_v + abs(bus_v(:,1))*ones(1,length(tL));
y_ang = y_ang + angle(bus_v(:,1))*ones(1,length(tL));
y_qinv = y_qinv + g.reec.reec3(:,1)*ones(1,length(tL));
y_vloc = y_vloc + g.reec.reec1(:,1)*ones(1,length(tL));
bus_vL = transpose(y_v.*exp(1j*y_ang)).';
save('../analysis/mat/reec_example1_linear_sim','tL','bus_vL','u1','u2','y_qinv','y_vloc');

%% Plot Nonlinear vs Linearized
clear all; close all; clc;
cd('../analysis');
load('./mat/reec_example1_nonlinear_sim','t','bus_v','g');
load('./mat/reec_example1_linear_sim','tL','bus_vL','u1','u2','y_qinv','y_vloc');

nb1 = 5;   % first bus to plot (bus 10 in Kundur system)
nb2 = 12;  % second bus to plot (bus 110 in Kundur system)

% tustin-bessel derivative filter ~10 Hz corner
Fs = round(1/median(diff(t)));
t_step = 1/Fs;

[b,a] = besself(2,2*pi*10);
b = [b,0]/2/pi;
a = [0,a];

[bd,ad] = bilinear(b,a,Fs);

bus_f = zeros(size(bus_v));
bus_fL = zeros(size(bus_vL));
for ii = 1:size(bus_v,1)
    Zi = filtic(bd,ad,zeros(3,1),angle(bus_v(ii,1))*ones(3,1));
    bus_f(ii,:) = 1e3*filter(bd,ad,angle(bus_v(ii,:)),Zi);
    bus_fL(ii,:) = 1e3*filter(bd,ad,angle(bus_vL(ii,:)),Zi);
end

fig{1} = figure;
ax{1}(1) = subplot(2,1,1,'parent',fig{1});
ax{1}(2) = subplot(2,1,2,'parent',fig{1});
hold(ax{1}(1),'on');
hold(ax{1}(2),'on');

plot(ax{1}(1),t,abs(bus_v(nb1,:)),'b',tL,abs(bus_vL(nb1,:)),'m:','linewidth',2);
plot(ax{1}(2),t,abs(bus_v(nb2,:)),'b',tL,abs(bus_vL(nb2,:)),'m:','linewidth',2);
ylabel(ax{1}(1),'$|V_{1}|$ (pu)','interpreter','latex');
ylabel(ax{1}(2),'$|V_{3}|$ (pu)','interpreter','latex');
xlabel(ax{1}(2),'Time (s)','interpreter','latex');

legend(ax{1}(1),'nonlinear','linearized','location','best');
legend(ax{1}(2),'nonlinear','linearized','location','best');

fig{2} = figure;
ax{2}(1) = subplot(2,1,1,'parent',fig{2});
ax{2}(2) = subplot(2,1,2,'parent',fig{2});
hold(ax{2}(1),'on');
hold(ax{2}(2),'on');

plot(ax{2}(1),t,bus_f(nb1,:),'b',tL,bus_fL(nb1,:),'m:','linewidth',2);
plot(ax{2}(2),t,bus_f(nb2,:),'b',tL,bus_fL(nb2,:),'m:','linewidth',2);

ylabel(ax{2}(1),'$\Delta f_{1}$ (mHz)','interpreter','latex');
ylabel(ax{2}(2),'$\Delta f_{3}$ (mHz)','interpreter','latex');
xlabel(ax{2}(2),'Time (s)','interpreter','latex');

legend(ax{2}(1),'nonlinear','linearized','location','best');
legend(ax{2}(2),'nonlinear','linearized','location','best');

fig{3} = figure;
ax{3}(1) = subplot(2,1,1,'parent',fig{3});
ax{3}(2) = subplot(2,1,2,'parent',fig{3});
hold(ax{3}(1),'on');
hold(ax{3}(2),'on');

plot(ax{3}(1),t,g.reec.reec3(1,:),'b',tL,y_qinv(1,:),'m:','linewidth',2);
plot(ax{3}(2),t,g.reec.reec3(2,:),'b',tL,y_qinv(2,:),'m:','linewidth',2);

ylabel(ax{3}(1),'$Q_{1}$ (pu)','interpreter','latex');
ylabel(ax{3}(2),'$Q_{3}$ (pu)','interpreter','latex');
xlabel(ax{3}(2),'Time (s)','interpreter','latex');

legend(ax{3}(1),'nonlinear','linearized','location','best');
legend(ax{3}(2),'nonlinear','linearized','location','best');

for ii = 1:numel(fig)
    for jj = 1:numel(ax{ii})
        v = axis(ax{ii}(jj));
        ystretch = 0.10*(v(4) - v(3));
        axis(ax{ii}(jj),[v(1), 10, v(3)-ystretch, v(4)+ystretch]);
    end
end

print(fig{1},'-dpng','-r600','./fig/reec_example1_voltage_mag.png');
print(fig{2},'-dpng','-r600','./fig/reec_example1_bus_freq.png');
print(fig{3},'-dpng','-r600','./fig/reec_example1_pwr_cmd.png');

% eof
