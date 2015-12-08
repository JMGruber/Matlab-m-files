function curve = gaussian(x,A,FWHM,Offset)
%  Usage:  curve = gaussian(x,A,FWHM,Offset);
%     x is a vector of x data values
%     A is the amplitude of the gaussian
%     FWHM is the full width at half maximum
%     Offset is the center of the gaussian
%     curve is a vector of y values corresponding 
%         to the provided x values.

curve=A.*exp((-4*log(2)/FWHM.^2).*(x-Offset).^2);


