% samer_main_script.m

% These scripts and data provide example code for the method described in:
% Polak et al. 2021, "Scout Accelerated Motion Estimation and Correction (SAMER)"


% This script demonstrates the SAMER motion estimation and reconstruction approach using simulated motion data. 
% In the provided examples the SAMER optimization is performed on different reordering schemes and scout scans and 
% the use of coil compression is explored. Ultimately, the SAMER motion estimation is evaluated on noisy k-space data.


%% preparation

set(0,'DefaultFigureWindowStyle','docked')
addpath('helper')


% load data
load data/img
load data/motion_gt
load data/shot_mask
load data/fminunc_options


% number of shots
N_shots = size(shot_mask_lin,3);

% plot the reference image used for the motion simulation
mosaic(img_gt,1,1,1, 'Reference image',[0,0.8])

% plot the motion parameters used for the motion simulation
subplot_title = ["Tx [voxel]","Ty [voxel]","Rz [°]"];
x_label = "shot";
y_label = ["Tx /mm", "Ty /mm", "Rz /°"] ;
figure(2), clf(), 
for idx = [1:3]
    subplot(3,1,idx), plot(motion_gt(:,idx).'), title(subplot_title(idx)), xlabel (x_label), ylabel(y_label(idx)), legend('ground truth'), grid on,axis([1,24,-3,3])
end
sgtitle('Ground truth motion')



%% Motion simulation with linear reordering

% plot the shot mask for linear reordering at R=2x2 acceleration
mosaic(sum(shot_mask_lin.*permute([1:N_shots],[1,3,2]),3),1,1,10,'Shot mask for linear reordering (R=2x2)'), colormap jet


% simulate the motion corrupted k-space data with linear reordering
params = [];
params.num_chan = 16;
params.img_size = [256,192];

kspace_data = reshape(apply2DSensePlusMotion(img_gt(:), params,sens_svd, shot_mask_lin, motion_gt , 0, 'ntransp'),[params.img_size,params.num_chan]);


% start motion estimation using R=2x2 accelerated low-resolution scout scan
motion_est = zeros(N_shots,3);

% load scout (R=2x2, 1x4x4 mm³)
img_in = img_scoutR4;

% plot low-resolution scout
mosaic(img_scoutR4,1,1,11,'Low-resolution scout (1x4x4 mm³, R=4)',[0,0.8])

% choose number of coil channels used for motion estimation
svd_coil_channels = 16;

% start per shot motion optimization
for shot = [1:N_shots]
    
    display(['shot: ', num2str(shot),'/',num2str(N_shots)]);

    % initialize motion parameters from zeros
    x0 = [0,0,0];
    
    % optimize over data-consistency error to estimate motion parameters
    [mot_params, loss_org] = fminunc(@(x)computeSensePlusMotionDataConsistency( x, img_in, kspace_data(:,:,1:svd_coil_channels), sens_svd(:,:,1:svd_coil_channels), shot_mask_lin(:,:,shot)), x0, options);

    % save the estimated motion parameters
    motion_est(shot,:) = mot_params;

end


% plot the estimated motion parameters
subplot_title = ["Tx [voxel]","Ty [voxel]","Rz [°]"];
x_label = "shot";
y_label = ["Tx /mm", "Ty /mm", "Rz /°"] ;

figure(12), clf(), 
for idx = [1:3]
    const_offset = mean(motion_gt(:,idx).'-motion_est(:,idx) .',2);   % de-trend motion parameters -> find constant offset due to intra-scan motion between scout and imaging scan
    motion_est(:,idx) = motion_est(:,idx)+const_offset;   
    subplot(3,1,idx), plot(motion_gt(:,idx).'), hold on, plot(motion_est(:,idx) .'),  xlabel (x_label), ylabel(y_label(idx)), legend('ground truth', 'estimated'), title(subplot_title(idx)), grid on, axis([1,24,-3,3])
end
sgtitle('Estimated motion using linear reordering and R=4 low-res scout')




%% Motion simulation with linear + checkered reordering (R=2x2 scout)

% plot shot mask for linear + checkered sampling
mosaic(sum(shot_mask_lincheck.*permute([1:N_shots],[1,3,2]),3),1,1,20,'Shot mask for linear+checkered reordering (R=4)'), colormap jet


% simulate the motion corrupted k-space data with linear reordering
params = [];
params.num_chan = 16;
params.img_size = [256,192];

kspace_data = reshape(apply2DSensePlusMotion(img_gt(:), params,sens_svd, shot_mask_lincheck, motion_gt , 0, 'ntransp'),[params.img_size,params.num_chan]);


% start the motion estimation for linear + checkered sampling
motion_est = zeros(N_shots,3);

% chose R=2x2 low-resolution scout scan
img_in = img_scoutR4;

% use 16 coil channels for motion estimation
svd_coil_channels = 16;


% start the motion optimization
for shot = [1:N_shots]
    
    display(['shot: ', num2str(shot),'/',num2str(N_shots)]);

    x0 = [0,0,0];
    [mot_params, loss_org] = fminunc(@(x)computeSensePlusMotionDataConsistency( x, img_in, kspace_data(:,:,1:svd_coil_channels), sens_svd(:,:,1:svd_coil_channels), shot_mask_lincheck_opt(:,:,shot)), x0, options);

    motion_est(shot,:) = mot_params;

end


% Plot the estimated motion parameters
subplot_title = ["Tx [voxel]","Ty [voxel]","Rz [°]"];
x_label = "shot";
y_label = ["Tx /mm", "Ty /mm", "Rz /°"] ;

figure(21), clf(), 
for idx = [1:3]
    const_offset = mean(motion_gt(:,idx).'-motion_est(:,idx) .',2);   % intra-scan motion between scout and imaging scan
    motion_est(:,idx) = motion_est(:,idx)+const_offset;
    subplot(3,1,idx), plot(motion_gt(:,idx).', 'LineWidth', 1.5), hold on, plot(motion_est(:,idx) .', 'LineWidth', 1.5),  xlabel (x_label), ylabel(y_label(idx)), legend('ground truth', 'estimated'), title(subplot_title(idx)), grid on, axis([1,24,-3,3])
end

sgtitle('Estimated motion using linear+checkered reordering and R=4 low-res scout')



%% Motion simulation with linear + checkered reordering (R=12 scout)

motion_est = zeros(N_shots,3);

% choose R=12 accelerated scout
img_in = img_scoutR12;

% plot the scout scan
mosaic(img_scoutR12,1,1,22,'Highly-accelerated low-resolution scout (1x4x4 mm³, R=12)',[0,0.8])

% choose 16 coil channels for motion optimization
svd_coil_channels = 16;


% start the motion optimization
for shot = [1:N_shots]
    
    display(['shot: ', num2str(shot),'/',num2str(N_shots)]);
   
    x0 = [0,0,0];
    [mot_params, loss_org] = fminunc(@(x)computeSensePlusMotionDataConsistency( x, img_in, kspace_data(:,:,1:svd_coil_channels), sens_svd(:,:,1:svd_coil_channels), shot_mask_lincheck_opt(:,:,shot)), x0, options);

    motion_est(shot,:) = mot_params;

end


% plot the estimated motion parameters

subplot_title = ["Tx [voxel]","Ty [voxel]","Rz [°]"];
x_label = "shot";
y_label = ["Tx /mm", "Ty /mm", "Rz /°"] ;

figure(23), clf(), 
for idx = [1:3]
    
    const_offset = mean(motion_gt(:,idx).'-motion_est(:,idx) .',2);  
    motion_est(:,idx) = motion_est(:,idx)+const_offset;
    subplot(3,1,idx), plot(motion_gt(:,idx).', 'LineWidth', 1.5), hold on, plot(motion_est(:,idx) .', 'LineWidth', 1.5),  xlabel (x_label), ylabel(y_label(idx)), legend('ground truth', 'estimated'), title(subplot_title(idx)), grid on, axis([1,24,-3,3])

end

sgtitle('Estimated motion using linear+checkered reordering and R=12 low-res scout')


%% Motion simulation with linear + checkered reordering and coil compression (R=12 scout)

motion_est = zeros(N_shots,3);

% choose R=12 accelerated low-resolution scout scan
img_in = img_scoutR12;

% use only 4 coil channels for motion estimation
svd_coil_channels = 4;


% start motion optimization
for shot = [1:N_shots]
    
    display(['shot: ', num2str(shot),'/',num2str(N_shots)]);
    
    x0 = [0,0,0];
    [mot_params, loss_org] = fminunc(@(x)computeSensePlusMotionDataConsistency( x, img_in, kspace_data(:,:,1:svd_coil_channels), sens_svd(:,:,1:svd_coil_channels), shot_mask_lincheck_opt(:,:,shot)), x0, options);

    motion_est(shot,:) = mot_params;

end


% plot the motion parameters
subplot_title = ["Tx [voxel]","Ty [voxel]","Rz [°]"];
x_label = "shot";
y_label = ["Tx /mm", "Ty /mm", "Rz /°"] ;

figure(24), clf(), 
for idx = [1:3]
    const_offset = mean(motion_gt(:,idx).'-motion_est(:,idx) .',2);  
    motion_est(:,idx) = motion_est(:,idx)+const_offset;
    subplot(3,1,idx), plot(motion_gt(:,idx).', 'LineWidth', 1.5), hold on, plot(motion_est(:,idx).', 'LineWidth', 1.5),  xlabel (x_label), ylabel(y_label(idx)), legend('ground truth', 'estimated'), title(subplot_title(idx)), grid on, axis([1,24,-3,3])
end

sgtitle('Estimated motion using linear+checkered reordering, N_{ch}=4 coil compression and R=12 low-res scout')



%% Perform the image reconstruction using conjugate gradient optimization

lsqr_tol = 1e-2;
lsqr_iter=5;


% SENSE reconstruction ( no motion correction)
[res, flag, relres, iter] = lsqr(@apply2DSense,kspace_data(:), lsqr_tol, lsqr_iter, [], [], [], params,sens_svd, sum(shot_mask_lincheck,3) );    
img_no_moco= reshape(res, params.img_size);

% plot the reconstructed SENSE image
mosaic(img_no_moco,1,1,25, 'SENSE image reconstruction (no motion correction)',[0,0.8])


% SAMER reconstruction ( no motion correction)
[res, flag, relres, iter] = lsqr(@apply2DSensePlusMotion,kspace_data(:), lsqr_tol, lsqr_iter, [], [], [], params,sens_svd, shot_mask_lincheck,motion_est,0 );    
img_moco= reshape(res, params.img_size);

% plot the reconstructed SAMER image
mosaic(img_moco,1,1,26, 'SAMER image reconstruction using estimated motion parameters',[0,0.8])



%% Motion simulation with noisy k-space data (linear + checkered reordering)


% simulate motion corrupted k-space data and added noise 
params = [];
params.num_chan = 16;
params.img_size = [256,192];

kspace_data_with_noise = reshape(apply2DSensePlusMotion(img_gt(:), params,sens_svd, shot_mask_lincheck, motion_gt , 1, 'ntransp'),[params.img_size,params.num_chan]);


% start the motion estimation
motion_est = zeros(N_shots,3);

% use R=12 low-resolution scout scan
img_in = img_scoutR12 ;

% use only 4 coil channels for the motion estimation
svd_coil_channels = 4;


% start the motion optimization
for shot = [1:N_shots]
    
    display(['shot: ', num2str(shot),'/',num2str(N_shots)]);

    x0 = [0,0,0];
    [mot_params, loss_org] = fminunc(@(x)computeSensePlusMotionDataConsistency( x, img_in, kspace_data_with_noise(:,:,1:svd_coil_channels), sens_svd(:,:,1:svd_coil_channels), shot_mask_lincheck_opt(:,:,shot)), x0, options);

    motion_est(shot,:) = mot_params;

end


% plot the motion parameters
subplot_title = ["Tx [voxel]","Ty [voxel]","Rz [°]"];
x_label = "shot";
y_label = ["Tx /mm", "Ty /mm", "Rz /°"] ;

figure(31), clf(), 
for idx = [1:3]
     const_offset = mean(motion_gt(:,idx).'-motion_est(:,idx) .',2);  
    motion_est(:,idx) = motion_est(:,idx)+const_offset;
    subplot(3,1,idx), plot(motion_gt(:,idx).', 'LineWidth', 1.5), hold on, plot(motion_est(:,idx) .', 'LineWidth', 1.5),  xlabel (x_label), ylabel(y_label(idx)), legend('ground truth', 'estimated'), title(subplot_title(idx)), grid on, axis([1,24,-3,3])
end

sgtitle('Estimated motion on noisy data using linear+checkered reordering, N_{ch}=4 coil compression and R=12 low-res scout')



%% Perform the image reconstruction on noisy k-space data

lsqr_tol = 1e-2;
lsqr_iter=5;

% SENSE reconstruction (no motion correction)
[res, flag, relres, iter] = lsqr(@apply2DSense,kspace_data_with_noise(:), lsqr_tol, lsqr_iter, [], [], [], params,sens_svd, sum(shot_mask_lincheck,3) );    
img_no_moco= reshape(res, params.img_size);

% plot SENSE reconstruction
mosaic(img_no_moco,1,1,32, 'SENSE image reconstruction with noisy data (no motion correction)',[0,0.8])


% SAMER reconstruction (no motion correction)
[res, flag, relres, iter] = lsqr(@apply2DSensePlusMotion,kspace_data_with_noise(:), lsqr_tol, lsqr_iter, [], [], [], params,sens_svd, shot_mask_lincheck,motion_est,0 );    
img_moco= reshape(res, params.img_size);

% plot SAMER reconstruction
mosaic(img_moco,1,1,33, 'SAMER image reconstruction with noisy data using estimated motion parameters',[0,0.8])

