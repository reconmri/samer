% computeSensePlusMotionDataConsistency.m

% These scripts and data provide example code for the method described in:
% Polak et al. 2021, "Scout Accelerated Motion Estimation and Correction (SAMER)"


% This function evalutes the SENSE plus motion forward model for a given
% input image and motion parameters and computes the data consistency error
% with respect to the acquired/simualted k-space data


function [obj]  = computeSensePlusMotionDataConsistency( x, img, kspace_data, sens, shot_mask)
    
    % pre-compute fft-scaling
    fft_norm = sqrt(size(img,1)*size(img,2));

    % perform translation and rotation 
    img_warp = imtranslate(imrotate(real(img),x(3), 'bilinear','crop'),x(1:2),'bilinear') + ...
                  + 1i*imtranslate(imrotate(imag(img),x(3), 'bilinear','crop'),x(1:2),'bilinear');
    
    % Multiply with the coil sensitivity and Fourier transform
    forwardModelEval = fftc(fftc(sens .*img_warp,1),2)/fft_norm;

    % Mask the data for a given shot
    forwardModelEval_masked = shot_mask.* forwardModelEval;
    kspace_data_masked = shot_mask.* kspace_data ;  
    
    % compute the data consistency error
    data_consistency_error = norm(forwardModelEval_masked(:)-kspace_data_masked(:)) / norm(kspace_data_masked(:));
    
    % return the data consistency error
    obj = double(data_consistency_error);
     
end



