function SM = IntShifts_testSM(trace,maxbg,tlive,thrliveI) % boolean

SM=true;
if tlive==length(trace)
    Q=find(trace(4:tlive-1,2)<maxbg*trace(1,1))+3;
    if isempty(Q)
        SM=false;
    else
        for j=1:length(Q)
           if j<length(Q)
               if max(trace(Q(j)-2,2),max(trace(Q(j)-1,2),max(trace(Q(j)+1,2),trace(Q(j)+2,2))))>thrliveI*trace(1,1) %1 particle (un)quenches in one step
                  break
               end
           else
               if max(trace(Q(j)-2,2),max(trace(Q(j)-1,2),trace(Q(j)+1,2)))>thrliveI*trace(1,1)  %1 particle (un)quenches in one step
                   break
               end
           end
        end
    end
else
    Q=find(trace(4:tlive+1,2)<maxbg*trace(1,1))+3;
    if isempty(Q)
        SM=false;
    else
        for j=1:length(Q)
           if Q(j)<tlive
                if max(trace(Q(j)-2,2),max(trace(Q(j)-1,2),max(trace(Q(j)+1,2),trace(Q(j)+2,2))))>thrliveI*trace(1,1) %1 particle (un)quenches in one step
                   break
                end
           end
        end
    end
end
if (SM)&&(j==length(Q))
   if ~((Q(j)>tlive)&&(max(trace(Q(j)-2,2),max(trace(Q(j)-1,2)))>thrliveI*trace(1,1)))
       SM=false;
   end
end