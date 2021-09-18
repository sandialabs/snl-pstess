% Example of how to use pwrmod_dyn.m to modulate the power into a bus.
% Data file = d2m_pwrmod1.m

%% Add pst path to MATLAB
% generate relative path generically
folderDepth = 6; % depth of current directory from main PST directory
pstVer = 'pstV3p1'; %pstSETO';
pathParts = strsplit(pwd, filesep);
PSTpath = pathParts(1);

for pNdx = 2:max(size(pathParts))-folderDepth
    PSTpath = [char(PSTpath), filesep, char(pathParts(pNdx))];
end
PSTpath = [char(PSTpath), filesep, pstVer, filesep];

addpath(PSTpath)
save PSTpath.mat PSTpath
clear folderDepth pathParts pNdx PSTpath
%% Run nonlinear simulation and store results
clear all; close all; clc
load PSTpath
delete([PSTpath 'DataFile.m']); copyfile('d2m_pwrmod1.m',[PSTpath 'DataFile.m']); %System data file
delete([PSTpath 'pwrmod_dyn.m']); copyfile('pwrmod_dyn_Example1.m',[PSTpath 'pwrmod_dyn.m']); %Modulation file
s_simu_Batch %Run PST
save('Exmaple1_NonlinearSim','t','bus_v'); %Save t and bus_v results

%% Build linear model, simulate, and store results
%Build linear model
clear all; close all; clc
load PSTpath
delete([PSTpath 'DataFile.m']); copyfile('d2m_pwrmod1.m',[PSTpath 'DataFile.m']); %System data file
delete([PSTpath 'pwrmod_dyn.m']); copyfile('pwrmod_dyn_Example1.m',[PSTpath 'pwrmod_dyn.m']); %Modulation file
svm_mgen_Batch %Conduct linearization
%%
%Simulate linear model
Gv = ss(a_mat,b_pwrmod_p,c_v,zeros(6,2));
Ga = ss(a_mat,b_pwrmod_p,c_ang,zeros(6,2));
tL = [0:1/120:10]';
u = [-0.0001*(stepfun(tL,1)-stepfun(tL,1.5)) 0.0002*(stepfun(tL,4)-stepfun(tL,4.5))]; 
load('Exmaple1_NonlinearSim','t','bus_v');
v = lsim(Gv,u,tL);
v = v + ones(size(v,1),1)*abs(bus_v(1:6,1))';
a = lsim(Ga,u,tL);
a = a + ones(size(a,1),1)*angle(bus_v(1:6,1))';
bus_vL = transpose(v.*exp(1i*a));
tL = tL';
save('Exmaple1_LinearSim','tL','bus_vL'); %Save linear results

%% Plot
clear all; close all; clc
load('Exmaple1_NonlinearSim','t','bus_v');
load('Exmaple1_LinearSim','tL','bus_vL');
figure
subplot(411)
nb = 2; %Bus to plot
plot(t,abs(bus_v(nb,:)),'k',tL,abs(bus_vL(nb,:)),'r');
ylabel(['bus ' num2str(nb) ' V (abs)'])
f = 1e3*angle(bus_v(nb,2:end)./bus_v(nb,1:end-1))./(2*pi*diff(t)); f = [f f(end)];
fL = 1e3*angle(bus_vL(nb,2:end)./bus_vL(nb,1:end-1))./(2*pi*diff(tL)); fL = [fL fL(end)];
subplot(412)
plot(t,f,'k',tL,fL,'r')
ylabel(['bus ' num2str(nb) ' (mHz)'])
subplot(413)
nb = 3;
nb = 2; %Bus to plot
plot(t,abs(bus_v(nb,:)),'k',tL,abs(bus_vL(nb,:)),'r');
ylabel(['bus ' num2str(nb) ' V (abs)'])
f = 1e3*angle(bus_v(nb,2:end)./bus_v(nb,1:end-1))./(2*pi*diff(t)); f = [f f(end)];
fL = 1e3*angle(bus_vL(nb,2:end)./bus_vL(nb,1:end-1))./(2*pi*diff(tL)); fL = [fL fL(end)];
subplot(414)
plot(t,f,'k',tL,fL,'r')
ylabel(['bus ' num2str(nb) ' (mHz)'])

set(gcf,'Position',[360 202 560 720]);