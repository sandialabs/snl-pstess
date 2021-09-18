function mtcsc_sig(t,k)
% Syntax: mtcsc_sig(t,k)
%
% Purpose: Defines modulation signal for tcsc control

%-----------------------------------------------------------------------------%

global g;  % declaring struct of global variables

if (g.tcsc.n_tcsc ~= 0)
    g.tcsc.tcsc_sig(:,k) = zeros(g.tcsc.n_tcsc,1);

    % if (t > 1.0)
    %     g.tcsc.tcsc_sig(1,k) = 0.02;
    % end
end

end  % function end

% eof
