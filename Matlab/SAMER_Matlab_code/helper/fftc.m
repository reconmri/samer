function data = fftc( data, num )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    data = fftshift(fft(fftshift(data,num),[], num),num);

end

