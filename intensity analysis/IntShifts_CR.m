function trace = IntShifts_CR(trace,thr2c,tlive);

% Remove excessively large intensities due to cosmic rays

largeints = find(trace(1:tlive,2)>thr2c*2);
for j=1:length(largeints)
    if largeints(j)==1
        newint=(trace(largeints(j)+1,2)+trace(largeints(j)+2,2))/2;
    elseif largeints(j)==tlive
        newint=(trace(largeints(j)-1,2)+trace(largeints(j)-2,2))/2;
    else
        newint=(trace(largeints(j)-1,2)+trace(largeints(j)+1,2))/2;
    end
    trace(largeints(j),2)=newint;
end