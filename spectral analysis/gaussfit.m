	function fitdata = gaussfit(data,FWHMq,Offsetq,BaseLine)  

% function fitdata = gaussfit(data,FWHMq,Offsetq,BaseLine);  
% Gaussfit requires the use of gaussian.m, mingauss.m and of course  
%     the file gaussfit.m.  
% Usage:    gaussfit(data,FWHMq,Offsetq,BaseLineq);  
% Where: data is a xy-matrix, FWHMq is the esimated FWHM  
% 	   Offsetq is the esimated offset of the peak (X-axis)  
%	   BaseLineq is the estimated baseline offset (y-axis)  
%  
% Output:   DataVector(1,4)  
% Where: DataVector(1) is the fitted Amplitude, A  
% 	   DataVector(2) is the fitted FHWM   
%	   DataVector(3) is the fitted Offset  
% 	   DataVector(4) is the fitted baseline   
% 
%	Data: 10/4/99 
%	Written by DSL, modified by TPJK
%
%	The input parameters will be guestimated from data if
%	desired....  just don't put them in...
%

%data=spec; FWHMq=20;Offsetq=680;BaseLine=0;

warning off;
ShowDataVar=1;

data=sortrows(data);		% Get into ascending form...
  
[Aq,Index]=max(data(:,2));  
Option=[0;5e-8;5e-8]; 
Option(14)=20000; 

if exist('BaseLine','var')==0 
	Sizer=max(size(data));
	BaseLine=mean(data([1,2,Sizer-1,Sizer],2));
end


if exist('Offsetq','var')==0 
	Offsetq=data(Index,1);
end


if exist('FWHMq','var')==0 
	TEMPDATA=data;
%	TEMPDATA(:,2)=normal(data(:,2));
	TEMPDATA(:,2)=(data(:,2))./max(data(:,2));
    
	%assuming the data is somewhat resolved....
	A=find(( TEMPDATA(:,2) > 0.3) & (TEMPDATA(:,2) < 0.7));
	A1=A(find(A<Index));
	A2=A(find(A>Index));

	if (isempty(A1) || isempty(A2))       
		FWHMq=20;
	else
		FWHMq=TEMPDATA(round(mean(A2)),1)-TEMPDATA(round(mean(A1)),1);
	end
end

fitdata=fminsearch('mingauss',[Aq,FWHMq,Offsetq,BaseLine],[],data);  
%fitdata=fminsearch(mingauss([Aq,FWHMq,Offsetq,BaseLine],data),Option,[],data);  

% if ShowDataVar==1
%     Step=mean(diff(data(:,1)))/30; 
%     GaussX=min(data(:,1)):Step:max(data(:,1));  
%     plot(data(:,1),data(:,2),'r-');  
%     hold on;  
%     plot(GaussX,gaussian(GaussX,fitdata(1),fitdata(2),fitdata(3))+fitdata(4),'b');  
%     title(['GAUSSFIT:  Width: ', num2str(fitdata(2)), '    Center: ', num2str(fitdata(3))])  
%     xlabel('X-axis');  
%     ylabel('Intensity');  
%     grid; 
%     hold off 
%     axis tight; 
% end

return 
 
 
%************************************************************************************************************************************ 
function MinimizingFunction = mingauss(Point,data)  
 
MinimizingFunction=sum((data(:,2)-gaussian(data(:,1),Point(1),Point(2),Point(3))-Point(4)).^2); 
 
return 
