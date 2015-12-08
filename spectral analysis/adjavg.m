function avgdata=adjavg(data,n)

% This function averages adjacent data points to smoothen the data
% n = number of points to average
newdata=[];
for i=1:length(data)
    j=round(n/2)-1;
    while (i<=j)||(i>length(data)-j-1)
        j=j-1;
    end
    if j<=0
        newdata(i)=data(i);
    else
        newdata(i)=mean(data(i-j:i-j+n-1));
    end;
end
avgdata=newdata';