function newdata=subtrCR(data,threshold)

% Remove outliers, usually due to Cosmic Rays, bmo interpolation 
% data = 1D array of intensity values, not full matrix!
% threshold = intensity/pixel for an outlier

if ~isempty(find(data>threshold,1))
    px=find(data>threshold);
    if (px(1)<4)||(px(length(px))>length(data)-4)
        newval=zeros(length(px),1);
    else
        startslope=mean(data(px(1)-4:px(1)-1));
        endslope=mean(data(px(length(px))+1:px(length(px))+4));
        slope=(endslope-startslope)/length(px);
        newval=startslope+slope*(1:length(px));
    end
    data(px)=newval;
end
newdata=data;