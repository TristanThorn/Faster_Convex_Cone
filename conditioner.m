function signal2 = conditioner(signal,rad,downsampl,laser) 

rad2 = sqrt(rad^2 + laser^2); 

blur = cos(2*pi*[-2*rad:2*rad]'/(4*rad))+1; 
blur = blur/sum(blur); 

signal2 = convn(signal,blur,'same'); 
signal2 = resample(signal2,1,downsampl); 