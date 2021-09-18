function mess_sig(t,k)
% Syntax: mess_sig(t,k)
%
% Purpose: Defines modulation signal for ess control

%-----------------------------------------------------------------------------%

global g;  % declaring struct of global variables

if (g.ess.n_ess ~= 0)
    g.ess.ess_sig(:,k) = zeros(g.ess.n_ess,1);

    % if (t > 1.0 && t < 1.5)
    %     g.ess.ess_sig(1,k) = g.ess.ess_sig(1,1) - 2e-3;
    %     g.ess.ess_sig(2,k) = g.ess.ess_sig(2,1) + 4e-3;
    % end
end

end  % function end

% eof
