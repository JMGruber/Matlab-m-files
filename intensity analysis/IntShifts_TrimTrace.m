function tlive = IntShifts_TrimTrace(trace,thrliveI)

tlive=0;        % survival time (in units of time resolution)
for t=length(trace):-1:3
   % if mean(trace(t-2:t,2))>thrliveI*timeres  % mean of 3 consecutive points
    if trace(t,2)>thrliveI*trace(1,1)
        if (trace(t-1,2))&&(trace(t-2,2))>thrliveI*trace(1,1)  % all 3 consecutive points > intensity threshold
           tlive=t;
           break
        end
    end
end