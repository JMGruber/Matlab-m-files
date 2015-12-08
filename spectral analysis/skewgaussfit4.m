function fitdata = skewgaussfit4(data,startq,fntype,lb,ub)  

% startq = 1x3 array with starting values of FWHM, Offset, and Skewness
%          Leave empty if the values should be determined bmo a Gaussian fit.

% fntype = type of skewed Gaussian used
% fntype=3: Danielis's function
% fntype=2: Mikas's function (more flexible)
% fntype=1: Composite function

% lb = lower bound of fitting parameters = 1x4 array [Ampl_lb FWHM_lb Offset_lb Skew_lb]
% ub = upper bound of fitting parameters = 1x4 array [Ampl_ub FWHM_ub Offset_ub Skew_ub]

if isempty(startq)
    paraest=gaussfit(data);
    start=[paraest(1) paraest(2) paraest(3) 0.01];
else
    start=[max(data(:,2)),startq];
end

switch fntype
    case 1
        fitdata = lsqcurvefit(@skewgaussian1,start,data(:,1),data(:,2),lb,ub);
    case 2
        fitdata = lsqcurvefit(@skewgaussian2,start,data(:,1),data(:,2),lb,ub);
    case 3
        fitdata = lsqcurvefit(@skewgaussian3,start,data(:,1),data(:,2),lb,ub);
end

Step=mean(diff(data(:,1)))/30; 
GaussX=min(data(:,1)):Step:max(data(:,1));  
plot(data(:,1),data(:,2),'-r');  
hold on;  
%plot(GaussX,skewgaussian2(GaussX,fitdata(1),fitdata(2),fitdata(3),fitdata(4))+fitdata(5),'b');  
plot(GaussX,skewgaussian1([fitdata(1),fitdata(2),fitdata(3),fitdata(4)],GaussX),'b');  
title(['SkewGAUSSFIT:  Width: ', num2str(fitdata(2)), '    Center: ',...
    num2str(fitdata(3)),'   Skewness:  ',num2str(fitdata(4)) ]);  
xlabel('X-axis');  
ylabel('Intensity');  
grid; 
hold off 
axis tight; 