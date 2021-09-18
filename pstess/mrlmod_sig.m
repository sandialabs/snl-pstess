function mrlmod_sig(t,k)
% Syntax: mrlmod_sig(t,k)
%
% Purpose: Defines modulation signal for rlmod control

%-----------------------------------------------------------------------------%

global g;  % declaring struct of global variables

if (g.rlmod.n_rlmod ~= 0)
    g.rlmod.rlmod_sig(:,k) = zeros(g.rlmod.n_rlmod,1);

    % if (t > 1.0)
    %     g.rlmod.rlmod_sig(1,k) = -0.10;
    % end
end

end  % function end

% eof
