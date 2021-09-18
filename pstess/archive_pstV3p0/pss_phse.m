function [phase,tw,t1,t2,t3,t4] = pss_phse(f)
tw=input('input the washout time constant in secs:[5]');
if isempty(tw);tw=5;end
t1=input('the first lead time constant in secs:[.2]');
if isempty(t1);t1=.2;end
t2=input('the first lag time constant in secs:[.02]');
if isempty(t2);t2=.02;end
t3=input('the second lead time constant in secs:[.2]');
if isempty(t3);t3=.2;end
t4=input('the second lag time constant in secs:[.02]');
if isempty(t4);t4=.02;end
i=sqrt(-1);
w=i*2*pi*f;
%keyboard
phase = 0.5*pi*ones(size(f))-angle(ones(size(f))+w*tw);
phase = phase + angle(ones(size(f))+w*t1)-angle(ones(size(f))+w*t2);
phase = phase + angle(ones(size(f))+w*t3)-angle(ones(size(f))+w*t4);
phase =phase*180/pi;
return

