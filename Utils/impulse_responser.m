% Generates signal for an impulse response corresponding to the derivative of a given order
% The impulse response is normalized based on the square root of the squared-sum of coefficients
function signal2 = impulse_responser(signal,order,steps,norm) 
signal2 = signal; 

if order==0
    coeff = 1;
end
if order==1
    coeff = [1 1];
end
if order==2
    coeff = [1 2 1];
end
if order==3
    coeff = [1 3 3 1];
end
if order==4
    coeff = [1 4 6 4 1];
end
if order==5
    coeff = [1 5 10 10 5 1];
end
if order==6
    coeff = [1 6 15 20 15 6 1];
end

for i=1:order 
    signal = circshift(signal,[steps 0]); signal(1:steps,:) = 0; 
    signal2 = signal2 + (-1)^i*coeff(i+1)*signal; 
end 

signal2 = signal2 * sqrt((order+1)/sum(coeff.^2)); % normalization 
if norm
    signal2 = signal2 * sqrt(pi);
end



