function mdc_sig(t,k)
% Syntax: mdc_sig(t,k)
%
% Purpose: Defines modulation signal for dc converter control

%-----------------------------------------------------------------------------%

global g;  % declaring struct of global variables

if (g.dc.n_conv ~= 0)
    g.dc.dc_sig(:,k) = zeros(g.dc.n_conv,1);

    % if (t > 1.0)
    %     g.dc.dc_sig(1,k) = 0.01;
    %     g.dc.dc_sig(1,k+1) = 0.01;
    %     g.dc.dc_sig(2,k) = 0.0;
    %     g.dc.dc_sig(2,k+1) = 0.0;
    % end
end

end  % function end

% eof
