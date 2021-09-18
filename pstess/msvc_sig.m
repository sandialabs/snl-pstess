function msvc_sig(t,k)
% Syntax: msvc_sig(t,k)
%
% Purpose: Defines modulation signal for svc control

%-----------------------------------------------------------------------------%

global g;  % declaring struct of global variables

if (g.svc.n_svc ~= 0)
    g.svc.svc_sig(:,k) = zeros(g.svc.n_svc,1);

    % if (t > 1.0)
    %     g.svc.svc_sig(1,k) = 0.02;
    % end
end

end  % function end

% eof
