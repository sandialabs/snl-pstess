% m file for pss design
% syntax: [tw,t1,t2,t3,t4] = pss_des(a,b,c,d)
% tw is the washout time constant
% t1 is the first lead time constant
% t2 is the first lag time constant
% t3 is the second lead time constant
% t4 is the second lag time constant
% requires a,b,c,d as output from svm_mgen
% calls statef and pss_phse
function[tw,t1,t2,t3,t4] = pss_des(a,b,c,d)
% input the frequency response range
fstart = input('enter the start frequency (Hz) [0.1]');
if isempty(fstart); fstart = 0.1;end
fstep = input('enter the frequency step (Hz) [0.01]');
if isempty(fstep); fstep = 0.01;end
fend = input('enter the end frequency (Hz) [2.0]');
if isempty(fend); fend = 2.0;end
% calculate the ideal frequency response
do_stab = 1;
% calculate the ideal pss magnitude and phase
[f,y,ymag,yphase]=statef(a,b,c,d,1,fstart,fstep,fend);
while(do_stab ==1),
  %calculate the pss phase lead
  [phase,tw,t1,t2,t3,t4] = pss_phse(f);
  %plot the ideal and pss phase
  plot(f,-yphase,'r',f,phase,'b');
  title('pss phase(blue) and ideal phase(red)');
  xlabel('frequency (Hz)');
  ylabel('phase (degrees)');
  more_plots = input('Do you wish to try another pss design: Y/N[Y]','s');
  if isempty(more_plots);more_plots='Y';end;
  if more_plots == 'y';more_plots = 'Y';end;
  if more_plots ~='Y';do_stab = 0; end;
end
return
