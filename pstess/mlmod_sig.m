function mlmod_sig(t,k)
% Syntax: mlmod_sig(t,k)
%
% Purpose: Defines modulation signal for lmod control

%-----------------------------------------------------------------------------%

global g;  % declaring struct of global variables

if (g.lmod.n_lmod ~= 0)
    g.lmod.lmod_sig(:,k) = zeros(g.lmod.n_lmod,1);

    % if (t > 1.0)
    %     g.lmod.lmod_sig(1,k) = -0.10;
    % end
end

end  % function end

% eof
