% apply2DSense.m

% These scripts and data provide example code for the method described in:
% Polak et al. 2021, "Scout Accelerated Motion Estimation and Correction (SAMER)"


% This function evaluates the SENSE model (forward and hermitean operation, no motion operators). 


function [ res, tflag ] = apply2DSense( in, params, sens_white,  mask, tflag )

    % precomupte the fft-scaling
    fft_norm = sqrt(params.img_size(1)*params.img_size(2));


    % hermitean operation is requested
    if strcmp(tflag,'transp')

        % prepare the data structures
        ksp = reshape(in, [params.img_size, params.num_chan]);    

        %Mask, inverse Fourier transform and multiply with conjugate
        %coil sensitivity, sum across the channels
        img = sum(conj(sens_white) .*ifftc(ifftc(mask.*ksp,1),2),3);  

        % scale using pre-computed FFT factor 
        res = img(:) * fft_norm ;


       
    % forward operation is requested
    else

        % prepare the data structures
        img = reshape(in, [params.img_size]);   
        
        % Multiply with the coil sensitivity, Fourier transform and mask
        ksp = mask.*fftc(fftc(sens_white .*img,1),2);   
        
         % scale using pre-computed FFT factor 
        res = ksp(:) / fft_norm;


    end

end

