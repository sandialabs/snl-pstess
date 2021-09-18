function mtg_sig(t,k)
% Syntax: mtg_sig(t,k)
%
% Purpose: Defines modulation signal for turbine power reference

%-----------------------------------------------------------------------------%

global g;  % declaring struct of global variables

if (g.tg.n_tg_tot ~= 0)
    g.tg.tg_sig(:,k) = zeros(g.tg.n_tg_tot,1);

    % if (t > 1.0)
    %     g.tg.tg_sig(1,k) = -0.01;
    %     % g.tg.tg_sig(1,k) = -1.0*t;
    % end
end

end  % function end

% eof
