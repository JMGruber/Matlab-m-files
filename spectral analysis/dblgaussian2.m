%function curve = dblgaussian(x,A1,FWHM1,Offset1,A2,FWHM2,Offset2)
function curve = dblgaussian2(bb,x)
A1=bb(1); FWHM1=bb(2); Offset1=bb(3); A2=bb(4); FWHM2=bb(5); Offset2=bb(6);
%  Usage:  curve = dblgaussian2([A1,FWHM1,Offset1,A2,FWHM2,Offset2],x);

curve=A1.*exp((-4*log(2)/FWHM1.^2).*(x-Offset1).^2)+A2.*exp((-4*log(2)/FWHM2.^2).*(x-Offset2).^2);