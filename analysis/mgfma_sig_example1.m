function mgfma_sig(t,k)
% Syntax: mgfma_sig(t,k)
%
% Purpose: Defines modulation signal for gfma control

%-----------------------------------------------------------------------------%

global g;  % declaring struct of global variables

if (g.gfma.n_gfma ~= 0)
    g.gfma.gfma_sig(:,k) = zeros(g.gfma.n_gfma,1);

    if (t >= 1.0 && t < 1.5)
        g.gfma.gfma_sig(1,k) = g.gfma.gfma_sig(1,1) + 2e-3;
    elseif (t >= 4.5 && t < 5.0)
        g.gfma.gfma_sig(2,k) = g.gfma.gfma_sig(2,1) - 1j*5e-3;
    end
end

end  % function end

% eof
