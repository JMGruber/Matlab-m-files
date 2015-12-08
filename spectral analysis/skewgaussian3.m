%function curve = skewgaussian3(x,A,FWHM,Offset,Skew)
function curve = skewgaussian3(bb,x)
A=bb(1); FWHM=bb(2); Offset=bb(3); Skew=bb(4); %BaseLine=bb(5);

if Skew==0 
    curve=gaussian(x,A,FWHM,Offset);
else
    sigma=Skew*FWHM/sinh(Skew);
    curve=real(A.*exp(-log(2)/Skew.^2.*(log(1+2.*Skew./sigma.*(x-Offset))).^2));
    if Skew>0
        curve2=curve(x>Offset-sigma/2/Skew);
        curve=[zeros(1,length(curve)-length(curve2)) curve2']';
    else
        curve2=curve(x<Offset-sigma/2/Skew);
        curve=[curve2' zeros(1,length(curve)-length(curve2))]';
    end
%    curve=real(A.*exp(-log(2)/Skew.^2.*(log(1+2.*Skew./sigma.*(x-Offset))).^2))+BaseLine;
    %curve=A.*exp((-4*log(2)/Skew.^2).*(log(1+g*(x-Offset))).^2);
end   
