function curve = skewgaussian1(bb,x)
A=bb(1); FWHM=bb(2); Offset=bb(3); Skew=bb(4);
%
%	Composite skewed Gaussian
%
%


Sizer=size(x);

if Sizer(1)==1
   x=x';
elseif Sizer(2)==1
   %do nothing and continue
else
   disp('x is not a vector');
	return;   
end

Rightside=find(x>Offset);
Leftside=find(x<=Offset);

if isempty(Rightside)==0
	curveR(:,1)=x(Rightside);
end
if isempty(Leftside)==0
   curveL(:,1)=x(Leftside);
end

%if (Skew < 1e-3 | Skew > -1e-3)
%   Skew=0;
%end



if Skew>0 && isempty(Rightside)==0 && isempty(Leftside)==0
	sigma=FWHM*sqrt(log(2))*exp(-Skew/2);
    %xplus=sigma/sqrt(log(2))*(exp(Skew/2-1))/(1-exp(-Skew));    %skewness enhances a lot
    %FWHMGauss=2*(FWHM-xplus);                                   %when these 2 lines are removed!
    g=sqrt(log(2))*(1-exp(-abs(Skew)))/sigma;

    curveR(:,2)=A.*exp((-4*log(2)/Skew.^2).*(log(1+g*(x(Rightside)-Offset))).^2);
%	curveL(:,2)=A.*exp((-4*log(2)/FWHMGauss.^2).*(x(Leftside)-Offset).^2);
	curveL(:,2)=A.*exp((-4*log(2)/FWHM.^2).*(x(Leftside)-Offset).^2);
   
    curve=([curveL' curveR'])';
	curve=sortrows(curve);

	curve(:,1)=[];

elseif Skew<0  && isempty(Rightside)==0 && isempty(Leftside)==0
    sigma=FWHM*sqrt(log(2))*exp(Skew/2);
    %xplus=sigma/sqrt(log(2))*(exp(-Skew/2-1))/(1-exp(Skew));    %skewness enhances a lot
    %FWHMGauss=2*(FWHM-xplus);                                   %when these 2 lines are removed!
    g=sqrt(log(2))*(1-exp(-abs(Skew)))/sigma;
    
    curveL(:,2)=A.*exp((-4*log(2)/Skew.^2).*(log(1+g*(-x(Leftside)+Offset))).^2);
%	curveR(:,2)=A.*exp((-4*log(2)/FWHMGauss.^2).*(x(Rightside)-Offset).^2);
    curveR(:,2)=A.*exp((-4*log(2)/FWHM.^2).*(x(Rightside)-Offset).^2);
    
    curve=([curveR' curveL'])';
	curve=sortrows(curve);

	curve(:,1)=[];
else
   curve=gaussian(x,A,FWHM,Offset);
end   
