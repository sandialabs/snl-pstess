function mreec_sig(t,k)
% Syntax: mreec_sig(t,k)
%
% Purpose: Defines modulation signal for reec control
%          reec_sig, real part -- local voltage reference perturbation (pu)
%          reec_sig, imag part -- reactive power reference perturbation (pu)

%-----------------------------------------------------------------------------%

global g;  % declaring struct of global variables

if (g.reec.n_reec ~= 0)
    g.reec.reec_sig(:,k) = zeros(g.reec.n_reec,1);

    % if (t >= 1.0 && t < 1.5)
    %     g.reec.reec_sig(1,k) = g.reec.reec_sig(1,1) - 1j*3e-3;
    % elseif (t >= 4.5 && t < 5.0)
    %     g.reec.reec_sig(2,k) = g.reec.reec_sig(2,1) + 5e-3;
    % end
end

end  % function end

% eof
