function mexc_sig(t,k)
% Syntax: mexc_sig(t,k)
%
% Purpose: Defines modulation signal for exciter control

%-----------------------------------------------------------------------------%

global g;  % declaring struct of global variables

if (g.exc.n_exc ~= 0)
    g.exc.exc_sig(:,k) = zeros(g.exc.n_exc,1);

    % if (t > 1.0)
    %     g.exc.exc_sig(1,k) = 0.01;
    % end
end

end  % function end

% eof
