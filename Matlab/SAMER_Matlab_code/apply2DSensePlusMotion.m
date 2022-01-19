% apply2DSensePlusMotion.m

% These scripts and data provide example code for the method described in:
% Polak et al. 2021, "Scout Accelerated Motion Estimation and Correction (SAMER)"


% This function evaluates the SENSE+motion model (forward and hermitean operation). 
% For motion simulations with added noise, gaussian random noise can be added on demand.


function [ res, tflag ] = apply2DSensePlusMotion( in, params, sens_white,  mask, motion_params, addNoise, tflag )


    % precomupte the fft-scaling
    fft_norm = sqrt(params.img_size(1)*params.img_size(2));

    % hermitean operation is requested
    if strcmp(tflag,'transp')

        % prepare the data structures
        img = zeros(params.img_size);
        ksp = reshape(in, [params.img_size, params.num_chan]);

        % evaluate the hermitean SENSE + motion model for each shot
       for shot = [1:size(mask,3)]

            %Mask, inverse Fourier transform and multiply with conjugate
            %coil sensitivity, sum across the channels
            temp = sum(conj(sens_white) .*ifftc(ifftc(mask(:,:,shot).*ksp,1),2),3);

            % perform inverse translation and rotation
            img = img + imrotate(imtranslate(real(temp),-motion_params(shot,1:2),'bilinear'),-motion_params(shot,3), 'bilinear','crop') + ...
                + 1i*imrotate(imtranslate(imag(temp),-motion_params(shot,1:2),'bilinear'),-motion_params(shot,3), 'bilinear','crop');

       end

       % scale the output using the pre-computed FFT scale factor
       res = img(:) * fft_norm ;


       
    % forward operation is requested
    else
        
        % pre-compute random gaussian noise for each channel image
        noise = 0;
        if addNoise == 1
            noise = 0.03.*rand(size( sens_white));
        end
        
        
        % prepare the data-structures
        img = reshape(in, [params.img_size]);
        ksp = zeros([params.img_size, params.num_chan]);

        % evaluate the SENSE + motion model for each shot
        for shot = [1:size(mask,3)]

            % translate and rotate the current image estimate
            img_warp = imtranslate(imrotate(real(img),motion_params(shot, 3), 'bilinear','crop'),motion_params(shot, 1:2),'bilinear') + ...
                + 1i*imtranslate(imrotate(imag(img),motion_params(shot, 3), 'bilinear','crop'),motion_params(shot,1:2),'bilinear');

            % multiply with coil sensitivity, Fourier transform and mask
            % the data. If requested add random gaussian noise to each coil channel
            % image (Note, this option is only available during the simulation
            % of corrupted k-space data)
            ksp = ksp + mask(:,:,shot).*fftc(fftc(noise+sens_white .*img_warp,1),2);

        end

        % scale the output using the pre-computed FFT scale factor
        res = ksp(:) / fft_norm ;


    end

end

