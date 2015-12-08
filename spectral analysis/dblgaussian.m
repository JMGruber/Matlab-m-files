function curve = dblgaussian(x,A1,FWHM1,Offset1,A2,FWHM2,Offset2)
%  Usage:  curve = dblgaussian(x,A1,FWHM1,Offset1,A2,FWHM2,Offset2);
%     x is a vector of x data values
%     A1 is the amplitude of the first gaussian
%     FWHM1 is the full width at half maximum of the first gaussian
%     Offset1 is the center of the first gaussian
%     A2 is the amplitude of the second gaussian
%     FWHM2 is the full width at half maximum of the second gaussian
%     Offset2 is the center of the second gaussian
%     curve is a vector of y values corresponding 
%         to the provided x values.

curve=A1.*exp((-4*log(2)/FWHM1.^2).*(x-Offset1).^2)+A2.*exp((-4*log(2)/FWHM2.^2).*(x-Offset2).^2);