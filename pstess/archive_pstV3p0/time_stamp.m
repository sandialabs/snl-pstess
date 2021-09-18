function tstr = time_stamp
y = clock;
tstr = [num2str(y(4)) ':' num2str(y(5))];