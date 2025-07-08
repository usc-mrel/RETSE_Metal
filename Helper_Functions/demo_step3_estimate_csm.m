% demo_step3_estimate_csm.m

%% Start a stopwatch timer
start_time = tic;

%% Read a .json file
tstart = tic; fprintf('%s: Reading a .json file: %s... ', datetime, json_file);
fid = fopen(json_file); 
json_txt = fread(fid, [1 inf], 'char=>char'); 
fclose(fid); 
json = jsondecode(json_txt);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% Define the full path of a filename
%--------------------------------------------------------------------------
if ispc
    siemens_twix_file  = strrep(json.siemens_twix_file, '/', '\');
    ismrmrd_data_file  = strrep(json.ismrmrd_data_file, '/', '\');
    ismrmrd_noise_file = strrep(json.ismrmrd_noise_file, '/', '\');
    output_path        = strrep(json.output_path, '/', '\');
else
    siemens_twix_file  = json.siemens_twix_file;
    ismrmrd_data_file  = json.ismrmrd_data_file;
    ismrmrd_noise_file = json.ismrmrd_noise_file;
    output_path        = json.output_path;
end

%--------------------------------------------------------------------------
% Define the BART directory
%--------------------------------------------------------------------------
bart_path = json.bart_path;

%--------------------------------------------------------------------------
% Reconstruction parameters
%--------------------------------------------------------------------------
lambda                 = json.recon_parameters.lambda;                 % l2 regularization parameter
maxiter                = json.recon_parameters.maxiter;                % PCG maximum iteration 
slice_type             = json.recon_parameters.slice_type;             % type of an excitation slice: "curved" vs "flat"
cal_size               = json.recon_parameters.cal_size.';             % size of calibration region
gnl_correction_flag    = json.recon_parameters.gnl_correction_flag;    % 1=yes, 0=no

%--------------------------------------------------------------------------
% Number of slices
%--------------------------------------------------------------------------
if isfield(json, 'nr_slices')
    nr_slices = json.nr_slices;
else
    nr_slices = 1;
end

%% Make an output path
mkdir(output_path);

%% Set up BART commands
%--------------------------------------------------------------------------
% Define a BART command
%--------------------------------------------------------------------------
if ispc
    command_prefix = 'wsl';
else
    command_prefix = '';
end
bart_command = sprintf('%s %s/bart', command_prefix, bart_path);

%--------------------------------------------------------------------------
% Translate from a Windows path to a WSL path 
%--------------------------------------------------------------------------
if ispc
    bart_output_path = strrep(output_path, '\', '/');
    bart_output_path = sprintf('/mnt/%s/%s/', lower(bart_output_path(1)), bart_output_path(4:end));
else
    bart_output_path = sprintf('%s/', output_path);
end

%% Read an ISMRMRD file (k-space data)
tstart = tic; fprintf('%s: Reading an ISMRMRD file: %s... ', datetime, ismrmrd_data_file);
if exist(ismrmrd_data_file, 'file')
    dset = ismrmrd.Dataset(ismrmrd_data_file, 'dataset');
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
else
    error('File %s does not exist.  Please generate it.' , ismrmrd_data_file);
end

%% Get imaging parameters from an XML header
header = ismrmrd.xml.deserialize(dset.readxml);

%--------------------------------------------------------------------------
% Encoding Space (Nkx, Nky, Nkz)
%--------------------------------------------------------------------------
encoded_fov(1) = header.encoding.encodedSpace.fieldOfView_mm.x * 1e-3; % [m] RO
encoded_fov(2) = header.encoding.encodedSpace.fieldOfView_mm.y * 1e-3; % [m] PE
encoded_fov(3) = header.encoding.encodedSpace.fieldOfView_mm.z * 1e-3; % [m] SL

Nkx = header.encoding.encodedSpace.matrixSize.x; % number of readout samples in k-space
Nky = header.encoding.encodedSpace.matrixSize.y; % number of phase-encoding steps in k-space
Nkz = header.encoding.encodedSpace.matrixSize.z; % number of slice-encoding steps in k-space

encoded_resolution = encoded_fov ./ [Nkx Nky Nkz]; % [m]

%--------------------------------------------------------------------------
% Number of receive channels
%--------------------------------------------------------------------------
Nc = header.acquisitionSystemInformation.receiverChannels;

%% Calculate the number of total slices
nr_recons = nr_slices;

%% Perform image reconstruction per slice
for idx = 1:nr_recons

    %% Get the slice number
    slice_number = ind2sub(nr_slices, idx);

    %% Read a .cfl file
    %----------------------------------------------------------------------
    % ksp_cal_cartesian (Nkx x Nky x Nkz x Nc)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path, sprintf('ksp_slc%d', slice_number));
    tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
    ksp = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % mask (Nkx x Nky x Nkz)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path, sprintf('mask_slc%d', slice_number));
    tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
    mask = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % circle_mask (Nkx x Nky x Nkz)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path, sprintf('circle_mask_slc%d_%s', slice_number, slice_type));
    tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
    circle_mask = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % x (Nkx x Nky x Nkz)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path, sprintf('x_slc%d_%s', slice_number, slice_type));
    tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
    x = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % y (Nkx x Nky x Nkz)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path, sprintf('y_slc%d_%s', slice_number, slice_type));
    tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
    y = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % z (Nkx x Nky x Nkz)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path, sprintf('z_slc%d_%s', slice_number, slice_type));
    tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
    z = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % dx (Nkx x Nky x Nkz)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path, sprintf('dx_slc%d_%s', slice_number, slice_type));
    tstart = tic; fprintf('%s: Writing a .cfl file: %s... ', datetime, cfl_file);
    dx = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % dy (Nkx x Nky x Nkz)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path, sprintf('dy_slc%d_%s', slice_number, slice_type));
    tstart = tic; fprintf('%s: Writing a .cfl file: %s... ', datetime, cfl_file);
    dy = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % dz (Nkx x Nky x Nkz)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path, sprintf('dz_slc%d_%s', slice_number, slice_type));
    tstart = tic; fprintf('%s: Writing a .cfl file: %s... ', datetime, cfl_file);
    dz = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % u (Nkx x Nky x Nkz)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path, sprintf('u_slc%d_%s', slice_number, slice_type));
    tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
    u = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % v (Nkx x Nky x Nkz)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path, sprintf('v_slc%d_%s', slice_number, slice_type));
    tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
    v = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % w (Nkx x Nky x Nkz)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path, sprintf('w_slc%d_%s', slice_number, slice_type));
    tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
    w = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % du (Nkx x Nky x Nkz)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path, sprintf('du_slc%d_%s', slice_number, slice_type));
    tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
    du = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % dv (Nkx x Nky x Nkz)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path, sprintf('dv_slc%d_%s', slice_number, slice_type));
    tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
    dv = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % dw (Nkx x Nky x Nkz)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path, sprintf('dw_slc%d_%s', slice_number, slice_type));
    tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
    dw = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %% Fix a sampling mask!
    mask(:,(mask(1,:) > 0)) = 1;

    %% Calculate parameters for Type-1 and Type-2 NUFFTs
    if gnl_correction_flag
        p1 = double(u + du) / encoded_fov(1) * (2 * pi); % RO [-pi,pi]
        p2 = double(v + dv) / encoded_fov(2) * (2 * pi); % PE [-pi,pi]
    else
        p1 = double(u) / encoded_fov(1) * (2 * pi); % RO [-pi,pi]
        p2 = double(v) / encoded_fov(2) * (2 * pi); % PE [-pi,pi]
    end

    support_mask = zeros(Nkx, Nky, Nkz, 'single');
    support_mask((abs(p1) < pi) & (abs(p2) < pi) & (circle_mask > 0)) = 1;

    p1 = p1(support_mask > 0);
    p2 = p2(support_mask > 0);

    nj = length(p1);
    eps = 1e-12;
    iflag = -1;

    %% Calculate "fake" coil sensitivity maps (Nkx x Nky x Nkz)
    sens = complex(ones(Nkx, Nky, Nkz, 'single'));

    %% Type-1 NUFFT based SENSE operators
    Ah = @(x) type1_sense_adjoint(x, sens, mask, nj, p1, p2, iflag, eps, support_mask);
    AhA = @(x) type1_sense_normal(x, sens, mask, nj, p1, p2, iflag, eps, support_mask, lambda);

    %% Type-1 NUFFT based MaxGIRF reconstruction
    cimg = complex(zeros(Nkx, Nky, Nkz, Nc, 'single'));
    for c = 1:Nc
        clear type1_sense_normal;
        tstart = tic; fprintf('%s:(SLC=%2d/%2d) Performing Type-1 NUFFT based SENSE reconstruction (c=%2d/%2d):\n', datetime, slice_number, nr_slices, c, Nc);
        [img, flag, relres, iter, resvec] = pcg(@(x) AhA(x), Ah(ksp(:,:,1,c)), 1e-4, maxiter);
        cimg(:,:,:,c) = reshape(img, [Nkx Nky]);
        fprintf('%s: done! (%6.4f/%6.4f sec)\n', datetime, toc(tstart), toc(start_time));
    end

    %% Write a .cfl file
    %----------------------------------------------------------------------
    % cimg (Nkx x Nky x Nkz x Nc)
    %----------------------------------------------------------------------
    cimg_filename = sprintf('cimg_cal_slc%d_%s_gnl%d_i%d_l%4.2f', slice_number, slice_type, gnl_correction_flag, maxiter, lambda);
    cfl_file = fullfile(output_path, cimg_filename);
    tstart = tic; fprintf('%s: Writing a .cfl file: %s... ', datetime, cfl_file);
    writecfl(cfl_file, cimg);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %% Calculate gridded k-space
    %----------------------------------------------------------------------
    % Siemens: k-space <=> image space
    % BART:                image space <=> k-space
    %----------------------------------------------------------------------
    tstart = tic; fprintf('%s:(SLC=%2d/%2d) Applying forward FFT to move from image space to k-space... ', datetime, slice_number, nr_slices);
    kgrid = cimg;
    for dim = 1:3
        kgrid = 1 / sqrt(size(kgrid,dim)) * fftshift(fft(ifftshift(kgrid, dim), [], dim), dim);
    end
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %% Write a .cfl file
    %----------------------------------------------------------------------
    % kgrid (Nkx x Nky x Nkz x Nc)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path, sprintf('kgrid_slc%d_%s_gnl%d', slice_number, slice_type, gnl_correction_flag));
    tstart = tic; fprintf('%s: Writing a .cfl file: %s... ', datetime, cfl_file);
    writecfl(cfl_file, kgrid);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %% Estimate coil sensitivity using ESPIRiT calibration
    %----------------------------------------------------------------------
    % ESPIRiT calibration
    %----------------------------------------------------------------------
    % Usage: ecalib [-t f] [-c f] [-k d:d:d] [-r d:d:d] [-m d] [-S] [-W] [-I] [-1]
    %               [-P] [-v f] [-a] [-d d] <kspace> <sensitivities> [<ev-maps>]
    %
    % Estimate coil sensitivities using ESPIRiT calibration.
    % Optionally outputs the eigenvalue maps.
    %
    % -t threshold     This determines the size of the null-space
    % -c crop_value    Crop the sensitivities if the eigenvalue is smaller than {crop_value}
    % -k ksize         kernel size
    % -r cal_size      Limits the size of the calibration region
    % -m maps          Number of maps to compute
    % -S               create maps with smooth transitions (Soft-SENSE)
    % -W               soft-weighting of the singular vectors
    % -I               intensity correction
    % -1               perform only first part of the calibration
    % -P               Do not rotate the phase with respect to the first principal component
    % -v variance      Variance of noise in data
    % -a               Automatically pick thresholds
    % -d level         Debug level
    % -h               help
    %----------------------------------------------------------------------
    kgrid_file   = strcat(bart_output_path, sprintf('kgrid_slc%d_%s_gnl%d', slice_number, slice_type, gnl_correction_flag));
    sens_file    = strcat(bart_output_path, sprintf('sens_slc%d_%s_gnl%d', slice_number, slice_type, gnl_correction_flag));
    ev_maps_file = strcat(bart_output_path, sprintf('ev_maps_slc%d_%s_gnl%d', slice_number, slice_type, gnl_correction_flag));
    command = sprintf('%s ecalib -t 0.001 -c 0 -k6:6:1 -r%d:%d:%d -m 1 -d5 %s %s %s', bart_command, cal_size(1), cal_size(2), cal_size(3), kgrid_file, sens_file, ev_maps_file);
    tstart = tic; fprintf('%s:(SLC=%2d/%2d)[BART] Estimating coil sensitivities using ESPIRiT calibration:\n%s\n', datetime, slice_number, nr_slices, command);
    [status_ecalib,result_ecalib] = system(command);
    fprintf('%s: done! (%6.4f/%6.4f sec)\n', datetime, toc(tstart), toc(start_time));

    %% Read a .cfl file
    %----------------------------------------------------------------------
    % sens (Nkx x Nky x Nkz x Nc)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path, sprintf('sens_slc%d_%s_gnl%d', slice_number, slice_type, gnl_correction_flag));
    tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
    sens = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %% Read a .cfl file
    %----------------------------------------------------------------------
    % ev_maps (Nkx x Nky x Nkz)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path, sprintf('ev_maps_slc%d_%s_gnl%d', slice_number, slice_type, gnl_correction_flag));
    tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
    ev_maps = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %% Display coil images
    nr_rows = 3;
    nr_cols = 2;
    cimg_montage = complex(zeros(Nkx * nr_rows, Nky * nr_cols, 'single'));
    sens_montage = complex(zeros(Nkx * nr_rows, Nky * nr_cols, 'single'));

    count = 0;
    for idx1 = 1:nr_rows
        for idx2 = 1:nr_cols
            idx1_range = (1:Nkx).' + (idx1 - 1) * Nkx;
            idx2_range = (1:Nky).' + (idx2 - 1) * Nky;
            count = count + 1;
            cimg_montage(idx1_range,idx2_range) = cimg(:,:,1,count);
            sens_montage(idx1_range,idx2_range) = sens(:,:,1,count);
            if count >= Nc
                break;
            end
        end
    end

    figure('Color', 'k', 'Position', [1 5 1239 973]);
    imagesc(abs(cimg_montage));
    axis image off;
    colormap(gray(256));
    caxis([0 5]);
    title({'Magnitude of coil images (calibration)', ...
           sprintf('Type-1 NUFFT based SENSE, SLC = %d, %s slice', slice_number, slice_type), ...
           sprintf('GNL = %d', gnl_correction_flag)}, 'Color', 'w', 'Interpreter', 'latex', 'FontWeight', 'normal', 'FontSize', 16);
    fig_filename = sprintf('cimg_cal_slc%d_%s_gnl%d_i%d_l%4.2f_mag', slice_number, slice_type, gnl_correction_flag, maxiter, lambda);
    export_fig(fullfile(output_path, fig_filename), '-r300', '-tif', '-c[0, 200, 310, 540]'); % [top,right,bottom,left]
    close gcf;

    figure('Color', 'k', 'Position', [1 5 1239 973]);
    imagesc(angle(cimg_montage) * 180 / pi);
    axis image off;
    caxis([-180 180]);
    colormap(hsv(256));
    hc = colorbar;
    set(hc, 'Color', 'w', 'FontSize', 14, 'Position', [0.8916 0.1530 0.0124 0.7185], 'TickLabelInterpreter', 'latex');
    title(hc, '[deg]', 'Color', 'w', 'Interpreter', 'latex', 'Position', [11.76135 529.225375 0]);
    title({'Phase of coil images (calibration)', ...
           sprintf('Type-1 NUFFT based SENSE, SLC = %d, %s slice', slice_number, slice_type), ...
           sprintf('GNL = %d', gnl_correction_flag)}, 'Color', 'w', 'Interpreter', 'latex', 'FontWeight', 'normal', 'FontSize', 16);
    fig_filename = sprintf('cimg_cal_slc%d_%s_gnl%d_i%d_l%4.2f_phase', slice_number, slice_type, gnl_correction_flag, maxiter, lambda);
    export_fig(fullfile(output_path, fig_filename), '-r300', '-tif', '-c[0, 200, 310, 540]');
    close gcf;

    figure('Color', 'k', 'Position', [1 5 1239 973]);
    imagesc(abs(sens_montage));
    axis image off;
    colormap(gray(256));
    caxis([0 1]);
    title({'Magnitude of coil sensitivity maps', ...
           sprintf('ESPIRiT, SLC = %d, %s slice', slice_number, slice_type), ...
           sprintf('GNL = %d', gnl_correction_flag)}, 'Color', 'w', 'Interpreter', 'latex', 'FontWeight', 'normal', 'FontSize', 16);
    fig_filename = sprintf('sens_cal_slc%d_%s_gnl%d_i%d_l%4.2f_mag', slice_number, slice_type, gnl_correction_flag, maxiter, lambda);
    export_fig(fullfile(output_path, fig_filename), '-r300', '-tif', '-c[0, 200, 310, 540]'); % [top,right,bottom,left]
    close gcf;

    figure('Color', 'k', 'Position', [1 5 1239 973]);
    imagesc(angle(sens_montage) * 180 / pi);
    axis image off;
    caxis([-180 180]);
    colormap(hsv(256));
    hc = colorbar;
    set(hc, 'Color', 'w', 'FontSize', 14, 'Position', [0.8916 0.1530 0.0124 0.7185], 'TickLabelInterpreter', 'latex');
    title(hc, '[deg]', 'Color', 'w', 'Interpreter', 'latex', 'Position', [11.76135 529.225375 0]);
    title({'Phase of coil sensitivity maps', ...
           sprintf('ESPIRiT, SLC = %d, %s slice', slice_number, slice_type), ...
           sprintf('GNL = %d', gnl_correction_flag)}, 'Color', 'w', 'Interpreter', 'latex', 'FontWeight', 'normal', 'FontSize', 16);
    fig_filename = sprintf('sens_cal_slc%d_%s_gnl%d_i%d_l%4.2f_phase', slice_number, slice_type, gnl_correction_flag, maxiter, lambda);
    export_fig(fullfile(output_path, fig_filename), '-r300', '-tif', '-c[0, 200, 310, 540]');
    close gcf;
end
