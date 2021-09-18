function mpm_sig(t,k)
% Syntax: mpm_sig(t,k)
%
% Purpose: Defines modulation signal for generator mechanical power

%-----------------------------------------------------------------------------%

global g;  % declaring struct of global variables

if (g.mac.n_pm ~= 0)
    g.mac.pm_sig(:,k) = zeros(g.mac.n_pm,1);

    % if (t > 1.0)
    %     g.mac.pm_sig(1,k) = 0.01;
    %     g.mac.pm_sig(2,k) = -0.01;
    % end
end

end  % function end

% eof
