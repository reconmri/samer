function data = ifftc( data, num )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    data = ifftshift(ifft(ifftshift(data,num),[], num),num);

end

