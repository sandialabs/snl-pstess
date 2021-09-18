function [tw,t1,t2,t3,t4] = pss_des(a,b,c,d)
% Syntax: [tw,t1,t2,t3,t4] = pss_des(a,b,c,d)
%
% Purpose: function for designing power system stabilizers
%
% tw is the washout time constant
% t1 is the first lead time constant
% t2 is the first lag time constant
% t3 is the second lead time constant
% t4 is the second lag time constant
%
% requires a,b,c,d as output from svm_mgen
% calls state_f and pss_phse

%-----------------------------------------------------------------------------%
% Version history
%
% Version: 1.1
% Author:  Ryan Elliott
% Date:    July 2020
% Purpose: Updated to accommodate batch processing mode
%
% Version: 1.0
% Author:  Graham Rogers and/or Joe Chow
% Note:    Initial version
%-----------------------------------------------------------------------------%

dypar = get_dypar();                % get parameter settings
batch_mode = dypar.batch_mode;      % batch processing mode

% input the frequency response range
if batch_mode
    fstart = 0.1;  % defaults
    fstep = 0.01;
    fend = 2.0;
else
    fstart = input('Enter the start frequency in Hz, [0.1]: ');
    if isempty(fstart)
        fstart = 0.1;
    end
    %
    fstep = input('Enter the frequency step in Hz, [0.01]: ');
    if isempty(fstep)
        fstep = 0.01;
    end
    %
    fend = input('Enter the end frequency in Hz, [2.0]: ');
    if isempty(fend)
        fend = 2.0;
    end
end

% calculate the ideal frequency response
do_stab = 1;

% calculate the ideal pss magnitude and phase
[f,y,ymag,yphase] = state_f(a,b,c,d,1,fstart,fstep,fend);
while (do_stab == 1)
    % calculate the pss phase lead
    [phase,tw,t1,t2,t3,t4] = pss_phse(f);

    % plot the ideal and pss phase
    plot(f,-yphase,'r',f,phase,'b');
    title('pss phase(blue) and ideal phase(red)');
    xlabel('frequency (Hz)');
    ylabel('phase (degrees)');
    % set(gca,'xscale','log');

    if batch_mode
        more_plots = 'n';
    else
        more_plots = input('Would you like to try another design, (y/n)[y]: ','s');
        more_plots = lower(more_plots);
        if isempty(more_plots)
            more_plots = 'n';  % default
        end
    end

    if ~strcmp(more_plots,'y')
        do_stab = 0;
    end
end

end  % function end

% eof
