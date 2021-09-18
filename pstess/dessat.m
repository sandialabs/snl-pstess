function gsat = dessat(a,isat)
% vectorized describing funtion for saturation
% 1 if a<=s
% called by mac_ind

%-----------------------------------------------------------------------------%

b_idx = find(a>=isat);
gsat = ones(length(a),1);
if ~isempty(b_idx)
    wt = atan2(isat(b_idx),sqrt(a(b_idx).*a(b_idx)-isat(b_idx).*isat(b_idx)));
    gsat(b_idx) = (2/pi)*(wt + 0.5*sin(2*wt));
end

end  % function end

% eof
