% demo_step4_cartesian_recon.mjson_file

%% Clean slate
%close all; clear all; clc;

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
    offres_path         = strrep(json.offres_path, '/', '\');
else
    siemens_twix_file  = json.siemens_twix_file;
    ismrmrd_data_file  = json.ismrmrd_data_file;
    ismrmrd_noise_file = json.ismrmrd_noise_file;
    output_path        = json.output_path;
    offres_path         = json.offres_path;
end

%--------------------------------------------------------------------------
% Reconstruction parameters
%--------------------------------------------------------------------------
%Lmax                  = json.recon_parameters.Lmax;                  % maximum rank of the SVD approximation of a higher-order encoding matrix
%L                     = json.recon_parameters.L;                     % rank of the SVD approximation of a higher-order encoding matrix
lambda                = json.recon_parameters.lambda;                % l2 regularization parameter
tol                   = json.recon_parameters.tol;                   % PCG tolerance
maxiter               = json.recon_parameters.maxiter;               % PCG maximum iteration 
slice_type            = json.recon_parameters.slice_type;            % type of an excitation slice: "curved" vs "flat"
readout_gridding_flag = json.recon_parameters.readout_gridding_flag; % 1=yes, 0=no
phase_correction_flag = json.recon_parameters.phase_correction_flag; % 1=yes, 0=no
gnl_correction_flag   = json.recon_parameters.gnl_correction_flag;   % 1=yes, 0=no
conc_correction_flag  = json.recon_parameters.conc_correction_flag;  % 1=yes, 0=no
offres_correction_flag = json.recon_parameters.offres_correction_flag; % 1=yes, 0=no

%--------------------------------------------------------------------------
% Number of slices
%--------------------------------------------------------------------------
if isfield(json, 'nr_slices')
    nr_slices = json.nr_slices;
else
    nr_slices = 1;
end

%--------------------------------------------------------------------------
% Number of repetitions
%--------------------------------------------------------------------------
if isfield(json, 'nr_repetitions')
    nr_repetitions = json.nr_repetitions;
else
    nr_repetitions = 1;
end

%--------------------------------------------------------------------------
% main_orientation (SAGITTAL/CORONAL/TRANSVERSAL = 0/1/2)
%--------------------------------------------------------------------------
if isfield(json, 'main_orientation')
    main_orientation = json.main_orientation;
else
    main_orientation = 2;
end

%--------------------------------------------------------------------------
% offres_sign
%--------------------------------------------------------------------------
if isfield(json, 'offres_sign')
    offres_sign = json.offres_sign;
else
    offres_sign = 1;
end

%% Make an output path
mkdir(output_path);

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
nr_recons = nr_slices * nr_repetitions;

%% Perform image reconstruction per slice
for idx = 1:nr_recons

    %% Get information about the current slice
    [slice_number, repetition_number] = ind2sub([nr_slices nr_repetitions], idx);

    %% Read a .cfl file
    %----------------------------------------------------------------------
    % ksp_cartesian (Nkx x Nky x Nkz x Nc)
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
    % sens (Nkx x Nky x Nkz x Nc)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path, sprintf('sens_slc%d_%s_gnl%d', slice_number, slice_type, gnl_correction_flag));
    tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
    sens = readcfl(cfl_file);
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
    
    %----------------------------------------------------------------------
    % displacement (Nkx x Nky x Nkz) saving the actual Siemens slice number
    %----------------------------------------------------------------------

    cfl_file = fullfile(offres_path, sprintf('displacement_slc%d_%s', slice_number, slice_type));
    tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
    displacement = readcfl(cfl_file);
    displacement = offres_sign * displacement;
    displacement = displacement;
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %% Fix a sampling mask!
    mask(:,(mask(1,:) > 0)) = 1;
    mask(1,(mask(end,:) > 0)) = 1;

    %% Calculate parameters for Type-1 and Type-2 NUFFTs
    if gnl_correction_flag && ~offres_correction_flag
        p1 = double(u + du) / encoded_fov(1) * (2 * pi); % RO [-pi,pi]
        p2 = double(v + dv) / encoded_fov(2) * (2 * pi); % PE [-pi,pi]
    elseif gnl_correction_flag && offres_correction_flag
        p1 = double(u + du + displacement) / encoded_fov(1) * (2 * pi); % RO [-pi,pi]
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

    %% Cartesian MaxGIRF operators
    Ah = @(x) type1_sense_adjoint(x, sens, mask, nj, p1, p2, iflag, eps, support_mask);
    AhA = @(x) type1_sense_normal(x, sens, mask, nj, p1, p2, iflag, eps, support_mask, lambda);

    %% Cartesian MaxGIRF reconstruction
    clear cartesian_maxgirf_normal;
    tstart = tic; fprintf('%s:(SLC=%2d/%2d)(REP=%2d/%2d) Performing Cartesian MaxGIRF reconstruction:\n', datetime, slice_number, nr_slices, repetition_number, nr_repetitions);
    [img, flag, relres, iter, resvec] = pcg(@(x) AhA(x), Ah(ksp), 1e-4, maxiter);
    img = reshape(img, [Nkx Nky]);
    fprintf('%s: done! (%6.4f/%6.4f sec)\n', datetime, toc(tstart), toc(start_time));

    %% Write a .cfl file
    %----------------------------------------------------------------------
    % img (Nkx x Nky x Nkz)
    %----------------------------------------------------------------------
    img_filename = sprintf('img_type1_slc%d_rep%d_%s_gridding%d_phc%d_conc%d_gnl%d_offres%d_i%d_l%4.2f', slice_number, repetition_number, slice_type, readout_gridding_flag, phase_correction_flag, conc_correction_flag, gnl_correction_flag, offres_correction_flag, maxiter, lambda);
    cfl_file = fullfile(output_path, img_filename);
    tstart = tic; fprintf('%s: Writing a .cfl file: %s... ', datetime, cfl_file);
    writecfl(cfl_file, img);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % support_mask (Nkx x Nky x Nkz)
    %----------------------------------------------------------------------
    mask_filename = sprintf('support_mask_slc%d_%s', slice_number, slice_type);
    cfl_file = fullfile(output_path, mask_filename);
    tstart = tic; fprintf('%s: Writing a .cfl file: %s... ', datetime, cfl_file);
    writecfl(cfl_file, support_mask);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %% Display the magnitude of an image
    FontSize = 14;

    c1 = floor(Nkx/2) + 1;
    c2 = floor(Nky/2) + 1;
    c3 = floor(Nkz/2) + 1;

    % FLASH, EPI
    zlimits = [-0.27 0.27];
    xlimits = [-0.15 0.15];

    if main_orientation == 2 % TRANSVERSAL=2
        Position = [680 476 1136 502];
        slice_direction = 'z';
        slice_offset = z(c1,c2) * 1e3;
        ax = 0;
        el = 90;
    elseif main_orientation == 0 % SAGITTAL=0
        Position = [1038 25 778 953];
        slice_direction = 'x';
        slice_offset = x(c1,c2) * 1e3;
        ax = -90;
        el = 0;
    elseif main_orientation == 1 % CORONAL=1
        Position = [680 28 1029 950];
        slice_direction = 'y';
        slice_offset = y(c1,c2) * 1e3;
        ax = 0;
        el = 0;
    end
    figure('Color', 'w', 'Position', Position);
    surf(x * 1e3, y * 1e3, z * 1e3, abs(img), 'EdgeColor', 'none');
    colormap(gray(256));
    axis image;
    xlim(xlimits * 1e3);
    zlim(zlimits * 1e3);
    if repetition_number > 1
        caxis([0 12]);
    else
        caxis([0 25]);
    end
    title({'Cartesian Type-1 based SENSE Recon.', sprintf('SLC = %d, REP = %d, %s slice, %s = %4.1f mm', slice_number, repetition_number, slice_type, slice_direction, slice_offset), ...
           sprintf('GNC/offres = %d/%d', gnl_correction_flag, offres_correction_flag), ...
           sprintf('max. iterations = %d, \\lambda = %4.2f', maxiter, lambda)}, ...
           'FontWeight', 'normal', 'FontSize', FontSize);
    xlabel('x [mm]', 'Interpreter', 'latex', 'FontSize', FontSize);
    ylabel('y [mm]', 'Interpreter', 'latex', 'FontSize', FontSize);
    zlabel('z [mm]', 'Interpreter', 'latex', 'FontSize', FontSize);
    set(gca, 'TickLabelInterpreter', 'latex');
    view(ax,el);
    set(gca, 'ZDir', 'reverse');
    drawnow;
    fig_filename = sprintf('img_type1_slc%d_rep%d_%s_gridding%d_phc%d_conc%d_gnl%d_offres%d_i%d_l%4.2f_mag', slice_number, repetition_number, slice_type, readout_gridding_flag, phase_correction_flag, conc_correction_flag, gnl_correction_flag, offres_correction_flag, maxiter, lambda);
    export_fig(fullfile(output_path, fig_filename), '-r300', '-tif');
    close gcf;

    %% Display the phase of an image
    figure('Color', 'w', 'Position', Position);
    surf(x * 1e3, y * 1e3, z * 1e3, angle(img) * 180 / pi, 'EdgeColor', 'none');
    colormap(hsv(256));
    axis image;
    xlim(xlimits * 1e3);
    zlim(zlimits * 1e3); 
    caxis([-180 180]);
    title({'Cartesian MaxGIRF', sprintf('SLC = %d, REP = %d, %s slice, %s = %4.1f mm', slice_number, repetition_number, slice_type, slice_direction, slice_offset), ...
           sprintf('gridding/PHC/CFC/GNC/offres = %d/%d/%d/%d/%d', readout_gridding_flag, phase_correction_flag, conc_correction_flag, gnl_correction_flag, offres_correction_flag), ...
           sprintf('max. iterations = %d, \\lambda = %4.2f', maxiter, lambda)}, ...
           'FontWeight', 'normal', 'FontSize', FontSize);
    xlabel('x [mm]', 'Interpreter', 'latex', 'FontSize', FontSize);
    ylabel('y [mm]', 'Interpreter', 'latex', 'FontSize', FontSize);
    zlabel('z [mm]', 'Interpreter', 'latex', 'FontSize', FontSize);
    set(gca, 'TickLabelInterpreter', 'latex');
    view(ax,el);
    set(gca, 'ZDir', 'reverse');
    drawnow;
    fig_filename = sprintf('img_type1_slc%d_rep%d_%s_gridding%d_phc%d_conc%d_gnl%d_offres%d_i%d_l%4.2f_phase', slice_number, repetition_number, slice_type, readout_gridding_flag, phase_correction_flag, conc_correction_flag, gnl_correction_flag, offres_correction_flag, maxiter, lambda);
    export_fig(fullfile(output_path, fig_filename), '-r300', '-tif');
    close gcf;

    %% Display a support mask
%     figure('Color', 'w', 'Position', Position);
%     surf(x * 1e3, y * 1e3, z * 1e3, support_mask, 'EdgeColor', 'none');
%     colormap(gray(256));
%     axis image;
%     xlim(xlimits * 1e3);
%     zlim(zlimits * 1e3);
%     caxis([0 1]);
%     title({'Cartesian MaxGIRF (mask)', sprintf('SLC = %d, REP = %d, %s slice, %s = %4.1f mm', slice_number, repetition_number, slice_type, slice_direction, slice_offset), ...
%            sprintf('gridding/PHC/CFC/GNC/offres = %d/%d/%d/%d/%d', readout_gridding_flag, phase_correction_flag, conc_correction_flag, gnl_correction_flag, offres_correction_flag), ...
%            sprintf('max. iterations = %d, \\lambda = %4.2f', maxiter, lambda)}, ...
%            'FontWeight', 'normal', 'FontSize', FontSize);
%     xlabel('x [mm]', 'Interpreter', 'latex', 'FontSize', FontSize);
%     ylabel('y [mm]', 'Interpreter', 'latex', 'FontSize', FontSize);
%     zlabel('z [mm]', 'Interpreter', 'latex', 'FontSize', FontSize);
%     set(gca, 'TickLabelInterpreter', 'latex', 'Color', 'k');
%     view(ax,el);
%     set(gca, 'ZDir', 'reverse');
%     drawnow;
%     export_fig(fullfile(output_path, mask_filename), '-r300', '-tif');
%     close gcf;
end
