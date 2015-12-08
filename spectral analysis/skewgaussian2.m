%function curve = skewgaussian2(x,A,FWHM,Offset,Skew)
function curve = skewgaussian2(bb,x)
A=bb(1); FWHM=bb(2); Offset=bb(3); Skew=bb(4);
%
%Usage: 	curve = skewgaussian([A,FWHM,Offset,Skew],x);
%

if Skew==0
    curve=gaussian(x,A,FWHM,Offset);
else
    sigma=FWHM*sqrt(log(2))*exp(-Skew/2);
    g=sqrt(log(2))*(1-exp(-Skew))/sigma;
    curve=real(A.*exp((-4*log(2)/Skew.^2).*(log(1+g*(x-Offset))).^2));
    if Skew>0
        curve2=curve(x>Offset-1/g);
        curve=[zeros(1,length(curve)-length(curve2)) curve2']';
    else
        curve2=curve(x<Offset-1/g);
        curve=[curve2' zeros(1,length(curve)-length(curve2))]';
    end
end
