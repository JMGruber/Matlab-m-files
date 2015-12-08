function fitdata = dblgaussfit(data,FWHMq1,Offsetq1,FWHMq2,Offsetq2,lb,ub)

wav=data(:,1);
x=data(:,2);
midwavx=(Offsetq1+Offsetq2)/2;
midwav=find(wav>midwavx,1);
[A1,Index1]=max(x(1:midwav));
[A2,Index2]=max(x(midwav:length(x)));
% if Index1==midwav
%     A1=A1*0.7;
% end
% if Index2==midwav
%     A2=A2*0.7;
% end

start=[A1,FWHMq1,Offsetq1,A2,FWHMq2,Offsetq2];

%fitdata = nlinfit(data(:,1),data(:,2),@dblgaussian2,[A1,FWHMq1,Offsetq1,A2,FWHMq2,Offsetq2]);
fitdata =lsqcurvefit(@dblgaussian2,start,data(:,1),data(:,2),lb,ub);
fitdata=real(fitdata);

    Step=mean(diff(data(:,1)))/30; 
    GaussX=min(data(:,1)):Step:max(data(:,1));  
    plot(data(:,1),data(:,2),'-ro');  
    hold on;  
    plot(GaussX,dblgaussian(GaussX,fitdata(1),fitdata(2),fitdata(3),fitdata(4),fitdata(5),fitdata(6)),'b');  
    plot(GaussX,gaussian(GaussX,fitdata(1),fitdata(2),fitdata(3)),':k');
    plot(GaussX,gaussian(GaussX,fitdata(4),fitdata(5),fitdata(6)),':k');
    title(['DOUBLEGAUSSFIT:  Width1: ', num2str(fitdata(2)), '    Center1: ', num2str(fitdata(3)),'  Width2: ', num2str(fitdata(5)), '    Center2: ', num2str(fitdata(6))])  
    xlabel('Wavelength');  
    ylabel('Intensity');  
    grid; 
    hold off 
    axis tight; 

return