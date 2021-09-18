function [phase,tw,t1,t2,t3,t4] = pss_phse(f)
% Syntax: [phase,tw,t1,t2,t3,t4] = pss_phse(f)

%-----------------------------------------------------------------------------%

dypar = get_dypar();                % get parameter settings
batch_mode = dypar.batch_mode;      % batch processing mode

if batch_mode
    tw = 10.0;
    t1 = 0.2;
    t2 = 0.02;
    t3 = 0.2;
    t4 = 0.02;
else
    tw = input('Enter the washout time constant in sec, [10]: ');
    if isempty(tw)
        tw = 10.0;
    end
    %
    t1 = input('Enter the first lead time constant in sec, [0.2]: ');
    if isempty(t1)
        t1 = 0.2;
    end
    %
    t2 = input('Enter the first lag time constant in sec, [0.02]: ');
    if isempty(t2)
        t2 = 0.02;
    end
    %
    t3 = input('Enter the second lead time constant in sec, [0.2]: ');
    if isempty(t3)
        t3 = 0.2;
    end
    %
    t4 = input('Enter the second lag time constant in sec, [0.02]: ');
    if isempty(t4)
        t4 = 0.02;
    end
end

jw = 1j*2*pi*f;

phase = 0.5*pi*ones(size(f)) - angle(ones(size(f)) + jw*tw);
phase = phase + angle(ones(size(f)) + jw*t1) - angle(ones(size(f)) + jw*t2);
phase = phase + angle(ones(size(f)) + jw*t3) - angle(ones(size(f)) + jw*t4);
phase = phase*180/pi;

end  % function end

% eof
