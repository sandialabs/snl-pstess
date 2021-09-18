function [f,ymag,yphase] = stab_f(a,b,c,tstate)
% Syntax: [f,ymag,yphase] = stab_f(a,b,c,tstate)
%
% Purpose: stabilizer frequency response
%
% Method:  calculates the phase lag through the exciter and generator
%          the inverse of this phase lag is the ideal pss phase lead
%          requires output from svm_smse
%
% Note:    generator data must be modified to have very high H (say 300)

%-----------------------------------------------------------------------------%
% Version history
%
% Author: Graham Rogers
% Date:   September 1994
%-----------------------------------------------------------------------------%

for k = 0:1:180
    jw = 1j*pi*(0.2+k*0.02)*diag(ones(tstate,1));
    f(k+1) = 0.1 + k*0.01;
    x = inv(-a(1:tstate,1:tstate)+jw)*b(1:tstate,1);
    y(k+1) = c(1,1:tstate)*x;
    ymag(k+1) = abs(y(k+1));

    yphase(k+1) = 180.0*angle(y(k+1))/pi - 180.0;
    if (yphase > 180.0)
        yphase = yphase - 360.0;
    end
    %
    if (yphase < -180.0)
        yphase = yphase + 360.0;
    end
end

end  % function end

% eof
