function [tw,t1,t2,t3,t4] = stab_d(a,b,c,tstate)
% Syntax: [tw,t1,t2,t3,t4] = stab_d(a,b,c,tstate)
%
% Purpose: function for power system stabilizer design
%
% Note:    requires system data from svm_smse

%-----------------------------------------------------------------------------%
% Version history
%
% Purpose: Added batch processing option to bypass user query
% Date:    July 2020
% Author:  Ryan Elliott
%
% Author:  Graham Rogers
% Date:    September 1994
%-----------------------------------------------------------------------------%

dypar = get_dypar();                % get parameter settings
batch_mode = dypar.batch_mode;      % batch processing mode

[f,yphase] = stab_f(a,b,c,tstate);

if batch_mode
    st_des = 'y';
else
    st_des = input('Would you like to design a stabilizer, (y/n)[y]: ','s')
    st_des = lower(st_des);
    if isempty(st_des)
        st_des = 'y';  % default
    end
end

if strcmp(st_des,'y')
    if batch_mode
        stab_name = 'null_gen';
    else
        stab_name = input('What is the name of the generator: ','s')
        stab_name = lower(stab_name);
        if isempty(stab_name)
            stab_name = 'null_gen';  % default
        end
    end

    point = length(f);
    ideal_phase = -yphase;

    if batch_mode
        st_type = 's';
    else
        st_type = input('Power or speed input stabilizer, (p/s)[s]: ','s');
        st_type = lower(st_type);
        if isempty(st_type)
            st_type = 's';  % default
        end
    end

    if ~strcmp(st_type,'s')
        new_stab = 'y';
        while strcmp(new_stab,'y')
            % power input stabilizer
            [phase,tw,t1,t2,t3,t4] = pss_phse(f);
            phase = phase + 90.*ones(1,point);
            subplot(111), plot(f,phase,f,ideal_phase);
            grid
            tit = ['Ideal and ',stab_name,' pss phase'];
            title(tit)

            if batch_mode
                new_stab = 'n';
            else
                new_stab = input('Would you like to try another design, (y/n)[n]: ','s');
                new_stab = lower(new_stab);
                if isempty(new_stab)
                    new_stab = 'n';  % default
                end
            end
        end
    else
        new_stab = 'y';
        while strcmp(new_stab,'y')
            % speed input stabilizer
            [phase,tw,t1,t2,t3,t4] = pss_phse(f);
            subplot(111), plot(f,phase,f,ideal_phase);
            grid
            tit = ['Ideal and ',stab_name,' phase'];
            title(tit)

            if batch_mode
                new_stab = 'n';
            else
                new_stab = input('Would you like to try another design, (y/n)[n]: ','s');
                new_stab = lower(new_stab);
                if isempty(new_stab)
                    new_stab = 'n';  % default
                end
            end
        end
    end
end

end  % function end

% eof
