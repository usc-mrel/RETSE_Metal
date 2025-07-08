% demo_step2_prepare_ksp.m

%% Set source directories
grad_file_path = fullfile('./Helper_Functions/src/GradientCoils');

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
% Define the BART path
%--------------------------------------------------------------------------
bart_path = json.bart_path;

%--------------------------------------------------------------------------
% Reconstruction parameters
%--------------------------------------------------------------------------
slice_type               = json.recon_parameters.slice_type;               % type of an excitation slice: "curved" vs "flat"
remove_oversampling_flag = json.recon_parameters.remove_oversampling_flag; % 1=yes, 0=no

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
% measurement information
%--------------------------------------------------------------------------
patient_position = header.measurementInformation.patientPosition;

gamma = 4257.59 * (1e4 * 2 * pi); % gyromagnetic ratio for 1H [rad/sec/T]
B0 = header.experimentalConditions.H1resonanceFrequency_Hz * (2 * pi / gamma); % [Hz] * [2pi rad/cycle] / [rad/sec/T] => [T]

%--------------------------------------------------------------------------
% sequence parameters
%--------------------------------------------------------------------------
TE = header.sequenceParameters.TE * 1e-3; % [msec] * [sec/1e3 msec] => * 1e-3 [sec]

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

kspace_encoding_step_1_center = header.encoding.encodingLimits.kspace_encoding_step_1.center;
kspace_encoding_step_2_center = header.encoding.encodingLimits.kspace_encoding_step_2.center;

kspace_encoding_step_1_maximum = header.encoding.encodingLimits.kspace_encoding_step_1.maximum;
kspace_encoding_step_1_minimum = header.encoding.encodingLimits.kspace_encoding_step_1.minimum;

kspace_encoding_step_2_maximum = header.encoding.encodingLimits.kspace_encoding_step_2.maximum;
kspace_encoding_step_2_minimum = header.encoding.encodingLimits.kspace_encoding_step_2.minimum;

acceleration_factor1 = header.encoding.parallelImaging.accelerationFactor.kspace_encoding_step_1;
acceleration_factor2 = header.encoding.parallelImaging.accelerationFactor.kspace_encoding_step_2;

%--------------------------------------------------------------------------
% Recon Space (Nx, Ny, Nz)
%--------------------------------------------------------------------------
recon_fov(1) = header.encoding.reconSpace.fieldOfView_mm.x * 1e-3; % [m] RO
recon_fov(2) = header.encoding.reconSpace.fieldOfView_mm.y * 1e-3; % [m] PE
recon_fov(3) = header.encoding.reconSpace.fieldOfView_mm.z * 1e-3; % [m] SL

Nx = header.encoding.reconSpace.matrixSize.x; % number of samples in image space (RO)
Ny = header.encoding.reconSpace.matrixSize.y; % number of samples in image space (PE)
Nz = header.encoding.reconSpace.matrixSize.z; % number of samples in image space (SL)

recon_resolution = recon_fov ./ [Nx Ny Nz]; % [m]

%--------------------------------------------------------------------------
% Readout oversampling factor
%--------------------------------------------------------------------------
readout_os_factor = encoded_fov(1) / recon_fov(1); % readout oversampling factor

%--------------------------------------------------------------------------
% Number of receive channels
%--------------------------------------------------------------------------
Nc = header.acquisitionSystemInformation.receiverChannels;

%% Write a .cfl file
%--------------------------------------------------------------------------
% TE (1 x 1)
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path, 'TE');
tstart = tic; fprintf('%s: Writing a .cfl file: %s... ', datetime, cfl_file);
writecfl(cfl_file, TE);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%% Read acquisitions
tstart = tic; fprintf('%s: Reading acquisitions... ', datetime);
raw_data = dset.readAcquisition(); % read all the acquisitions
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% Get data type
%--------------------------------------------------------------------------
acq_is_noise_measurement                = raw_data.head.flagIsSet('ACQ_IS_NOISE_MEASUREMENT');
acq_is_parallel_calibration             = raw_data.head.flagIsSet('ACQ_IS_PARALLEL_CALIBRATION');
acq_is_parallel_calibration_and_imaging = raw_data.head.flagIsSet('ACQ_IS_PARALLEL_CALIBRATION_AND_IMAGING');
acq_is_reverse                          = raw_data.head.flagIsSet('ACQ_IS_REVERSE');
acq_is_navigation_data                  = raw_data.head.flagIsSet('ACQ_IS_NAVIGATION_DATA');
acq_is_phasecorr_data                   = raw_data.head.flagIsSet('ACQ_IS_PHASECORR_DATA');
acq_is_hpfeedback_data                  = raw_data.head.flagIsSet('ACQ_IS_HPFEEDBACK_DATA');
acq_is_dummyscan_data                   = raw_data.head.flagIsSet('ACQ_IS_DUMMYSCAN_DATA');
acq_is_rtfeedback_data                  = raw_data.head.flagIsSet('ACQ_IS_RTFEEDBACK_DATA');
acq_is_surfacecoilcorrectionscan_data   = raw_data.head.flagIsSet('ACQ_IS_SURFACECOILCORRECTIONSCAN_DATA');

%% Get imaging data
img_data = raw_data.select(find(~acq_is_noise_measurement & ~acq_is_navigation_data & ~acq_is_parallel_calibration & ~acq_is_phasecorr_data));

%% Display data type
figure('Color', 'w', 'Position', [7 258 1809 734]);
subplot(2,5,1);
plot(acq_is_noise_measurement);
title('NOISE\_MEASUREMENT', 'FontWeight', 'normal');
grid on;

subplot(2,5,2);
plot(acq_is_parallel_calibration);
title('PARALLEL\_CALIBRATION', 'FontWeight', 'normal');
grid on;

subplot(2,5,3);
plot(acq_is_parallel_calibration_and_imaging);
title('PARALLEL\_CALIBRATION\_AND\_IMAGING', 'FontWeight', 'normal');
grid on;

subplot(2,5,4);
plot(acq_is_reverse);
title('REVERSE', 'FontWeight', 'normal');
grid on;

subplot(2,5,5);
plot(acq_is_navigation_data);
title('NAVIGATION\_DATA', 'FontWeight', 'normal');
grid on;

subplot(2,5,6);
plot(acq_is_phasecorr_data);
title('PHASECORR\_DATA', 'FontWeight', 'normal');
grid on;

subplot(2,5,7);
plot(acq_is_hpfeedback_data);
title('HPFEEDBACK\_DATA', 'FontWeight', 'normal');
grid on;

subplot(2,5,8);
plot(acq_is_dummyscan_data);
title('DUMMYSCAN\_DATA', 'FontWeight', 'normal');
grid on;

subplot(2,5,9);
plot(acq_is_rtfeedback_data);
title('RTFEEDBACK\_DATA', 'FontWeight', 'normal');
grid on;

subplot(2,5,10);
plot(acq_is_surfacecoilcorrectionscan_data);
title('SURFACECOILCORRECTIONSCAN\_DATA', 'FontWeight', 'normal');
grid on;

[filepath,json_filename,ext] = fileparts(json_file);

export_fig(fullfile(output_path, sprintf('data_type_%s', json_filename)), '-r300', '-tif');
close gcf;

fprintf('sum(acq_is_noise_measurement)                = %d\n', sum(acq_is_noise_measurement));
fprintf('sum(acq_is_parallel_calibration)             = %d\n', sum(acq_is_parallel_calibration));
fprintf('sum(acq_is_parallel_calibration_and_imaging) = %d\n', sum(acq_is_parallel_calibration_and_imaging));
fprintf('sum(acq_is_reverse)                          = %d\n', sum(acq_is_reverse));
fprintf('sum(acq_is_navigation_data)                  = %d\n', sum(acq_is_navigation_data));
fprintf('sum(acq_is_phasecorr_data)                   = %d\n', sum(acq_is_phasecorr_data));
fprintf('sum(acq_is_hpfeedback_data)                  = %d\n', sum(acq_is_hpfeedback_data));
fprintf('sum(acq_is_dummyscan_data)                   = %d\n', sum(acq_is_dummyscan_data));
fprintf('sum(acq_is_rtfeedback_data)                  = %d\n', sum(acq_is_rtfeedback_data));
fprintf('sum(acq_is_surfacecoilcorrectionscan_data)   = %d\n', sum(acq_is_surfacecoilcorrectionscan_data));

%% Parse an ISMRMRD header
%--------------------------------------------------------------------------
% ISMRMRD header
%--------------------------------------------------------------------------
% uint16_t version;                                    /**< First unsigned int indicates the version */
% uint64_t flags;                                      /**< bit field with flags */
% uint32_t measurement_uid;                            /**< Unique ID for the measurement */
% uint32_t scan_counter;                               /**< Current acquisition number in the measurement */
% uint32_t acquisition_time_stamp;                     /**< Acquisition clock */
% uint32_t physiology_time_stamp[ISMRMRD_PHYS_STAMPS]; /**< Physiology time stamps, e.g. ecg, breating, etc. */
% uint16_t number_of_samples;                          /**< Number of samples acquired */
% uint16_t available_channels;                         /**< Available coils */
% uint16_t active_channels;                            /**< Active coils on current acquisiton */
% uint64_t channel_mask[ISMRMRD_CHANNEL_MASKS];        /**< Mask to indicate which channels are active. Support for 1024 channels */
% uint16_t discard_pre;                                /**< Samples to be discarded at the beginning of  acquisition */
% uint16_t discard_post;                               /**< Samples to be discarded at the end of acquisition */
% uint16_t center_sample;                              /**< Sample at the center of k-space */
% uint16_t encoding_space_ref;                         /**< Reference to an encoding space, typically only one per acquisition */
% uint16_t trajectory_dimensions;                      /**< Indicates the dimensionality of the trajectory vector (0 means no trajectory) */
% float sample_time_us;                                /**< Time between samples in micro seconds, sampling BW */
% float position[3];                                   /**< Three-dimensional spatial offsets from isocenter */
% float read_dir[3];                                   /**< Directional cosines of the readout/frequency encoding */
% float phase_dir[3];                                  /**< Directional cosines of the phase */
% float slice_dir[3];                                  /**< Directional cosines of the slice direction */
% float patient_table_position[3];                     /**< Patient table off-center */
% ISMRMRD_EncodingCounters idx;                        /**< Encoding loop counters, see above */
% int32_t user_int[ISMRMRD_USER_INTS];                 /**< Free user parameters */
% float user_float[ISMRMRD_USER_FLOATS];               /**< Free user parameters */
%--------------------------------------------------------------------------
% Where EncodingCounters are defined as:
% uint16_t kspace_encode_step_1;    /**< e.g. phase encoding line number */
% uint16_t kspace_encode_step_2;    /**< e.g. partition encoding number */
% uint16_t average;                 /**< e.g. signal average number */
% uint16_t slice;                   /**< e.g. imaging slice number */
% uint16_t contrast;                /**< e.g. echo number in multi-echo */
% uint16_t phase;                   /**< e.g. cardiac phase number */
% uint16_t repetition;              /**< e.g. dynamic number for dynamic scanning */
% uint16_t set;                     /**< e.g. flow encoding set */
% uint16_t segment;                 /**< e.g. segment number for segmented acquisition */
% uint16_t user[ISMRMRD_USER_INTS]; /**< Free user parameters */
%--------------------------------------------------------------------------
number_of_samples           = double(max(img_data.head.number_of_samples));
discard_pre                 = double(max(img_data.head.discard_pre));
discard_post                = double(max(img_data.head.discard_post));
center_sample               = double(max(img_data.head.center_sample));
nr_channels                 = double(max(img_data.head.active_channels));
nr_phase_encoding_steps     = double(max(img_data.head.idx.kspace_encode_step_1)) + 1;
nr_partition_encoding_steps = double(max(img_data.head.idx.kspace_encode_step_2)) + 1;
nr_averages                 = double(max(img_data.head.idx.average)) + 1;
nr_slices                   = double(max(img_data.head.idx.slice)) + 1;
nr_contrasts                = double(max(img_data.head.idx.contrast)) + 1;
nr_phases                   = double(max(img_data.head.idx.phase)) + 1;
nr_repetitions              = double(max(img_data.head.idx.repetition)) + 1;
nr_sets                     = double(max(img_data.head.idx.set)) + 1;
nr_segments                 = double(max(img_data.head.idx.segment)) + 1;
nr_samples                  = number_of_samples - discard_pre - discard_post;

%--------------------------------------------------------------------------
% Get the dimensionality of the trajectory vector (0 means no trajectory)
%--------------------------------------------------------------------------
trajectory_dimensions = double(max(img_data.head.trajectory_dimensions));

%% Display an ISMRMRD header
fprintf('========================= ISMRMRD header =========================\n');
fprintf('encoded_fov        = %8.4f %8.4f %8.4f [mm]\n', encoded_fov(1) * 1e3, encoded_fov(2) * 1e3, encoded_fov(3) * 1e3);
fprintf('Nkx Nky Nkz        = %d      %d        %d\n', Nkx, Nky, Nkz);
fprintf('encoded_resolution = %8.4f %8.4f %8.4f [mm]\n', encoded_resolution(1) * 1e3, encoded_resolution(2) * 1e3, encoded_resolution(3) * 1e3);
fprintf('------------------------------------------------------------------\n');
fprintf('recon_fov          = %8.4f %8.4f %8.4f [mm]\n', recon_fov(1) * 1e3, recon_fov(2) * 1e3, recon_fov(3) * 1e3);
fprintf('Nx Ny Nz           = %d      %d        %d\n', Nx, Ny, Nz);
fprintf('recon_resolution   = %8.4f %8.4f %8.4f [mm]\n', recon_resolution(1) * 1e3, recon_resolution(2) * 1e3, recon_resolution(3) * 1e3);
fprintf('------------------------------------------------------------------\n');
fprintf('trajectory                  = %s\n', header.encoding.trajectory);
fprintf('number_of_samples           = %d\n', number_of_samples);
fprintf('discard_pre                 = %d\n', discard_pre);
fprintf('discard_post                = %d\n', discard_post);
fprintf('center_sample               = %d\n', center_sample);
fprintf('nr_channels                 = %d\n', nr_channels);
fprintf('nr_phase_encoding_steps     = %d\n', nr_phase_encoding_steps); % EPI factor
fprintf('nr_partition_encoding_steps = %d\n', nr_partition_encoding_steps);
fprintf('nr_averages                 = %d\n', nr_averages);
fprintf('nr_slices                   = %d\n', nr_slices);
fprintf('nr_contrasts                = %d\n', nr_contrasts);
fprintf('nr_phases                   = %d\n', nr_phases);
fprintf('nr_repetitions              = %d\n', nr_repetitions);
fprintf('nr_sets                     = %d\n', nr_sets);
fprintf('nr_segments                 = %d\n', nr_segments);
fprintf('==================================================================\n');

%% Define parameters for image reconstruction
%--------------------------------------------------------------------------
% Set the total number of voxels in image space
%--------------------------------------------------------------------------
N = Nkx * Nky * Nkz;

%% Calculate the receiver noise matrix (Nc x Nc)
[Psi,inv_L] = calculate_receiver_noise_matrix(ismrmrd_noise_file);

%% Read a Siemens .dat file
fprintf('%s: Reading a Siemens .dat file: %s\n', datetime, siemens_twix_file);
twix = mapVBVD(siemens_twix_file);

%% Get the name of a gradient set
if isfield(twix{1}.hdr.Meas, 'tGradientCoil')
    tGradientCoil = twix{1}.hdr.Meas.tGradientCoil;
elseif isfield(twix{2}.hdr.Meas, 'tGradientCoil')
    tGradientCoil = twix{2}.hdr.Meas.tGradientCoil;
end

%% Reduce the TWIX dataset
if length(twix) > 1
    twix = twix{end};
end

%% Get a slice normal vector from Siemens TWIX format
%--------------------------------------------------------------------------
% dNormalSag: Sagittal component of a slice normal vector (in the PCS)
%--------------------------------------------------------------------------
if isfield(twix.hdr.MeasYaps.sSliceArray.asSlice{1}.sNormal, 'dSag')
    dNormalSag = twix.hdr.MeasYaps.sSliceArray.asSlice{1}.sNormal.dSag;
else
    dNormalSag = 0;
end

%--------------------------------------------------------------------------
% dNormalCor: Coronal component of a slice normal vector (in the PCS)
%--------------------------------------------------------------------------
if isfield(twix.hdr.MeasYaps.sSliceArray.asSlice{1}.sNormal, 'dCor')
    dNormalCor = twix.hdr.MeasYaps.sSliceArray.asSlice{1}.sNormal.dCor;
else
    dNormalCor = 0;
end

%--------------------------------------------------------------------------
% dNormalTra: Transverse component of a slice normal vector (in the PCS)
%--------------------------------------------------------------------------
if isfield(twix.hdr.MeasYaps.sSliceArray.asSlice{1}.sNormal, 'dTra')
    dNormalTra = twix.hdr.MeasYaps.sSliceArray.asSlice{1}.sNormal.dTra;
else
    dNormalTra = 0;
end

%--------------------------------------------------------------------------
% dRotAngle: Slice rotation angle ("swap Fre/Pha")
%--------------------------------------------------------------------------
if isfield(twix.hdr.MeasYaps.sSliceArray.asSlice{1}, 'dInPlaneRot')
    dRotAngle = twix.hdr.MeasYaps.sSliceArray.asSlice{1}.dInPlaneRot; % [rad]
else
    dRotAngle = 0; % [rad]
end

%% Calculate a transformation matrix from the RCS to the GCS [r,c,s] <=> [PE,RO,SL]
R_rcs2gcs = [0    1    0 ; % [PE]   [0 1 0] * [r]
             1    0    0 ; % [RO] = [1 0 0] * [c]
             0    0    1]; % [SL]   [0 0 1] * [s]

%% Calculate a rotation matrix from the GCS to the PCS
[R_gcs2pcs, phase_sign, read_sign, main_orientation] = siemens_calculate_transform_gcs_to_pcs(dNormalSag, dNormalCor, dNormalTra, dRotAngle);

%% Calculate a rotation matrix from the PCS to the DCS
R_pcs2dcs = siemens_calculate_transform_pcs_to_dcs(patient_position);

%% Calculate a rotation matrix from the GCS to the DCS
R_gcs2dcs = R_pcs2dcs * R_gcs2pcs;

%% Read a Siemens .grad file
grad_filename = sprintf('coeff_%s.grad', tGradientCoil); % Gradient Coil daVinci (r0=0.25m, s.u.), Siemens 0.55T Aera
grad_file = fullfile(grad_file_path, grad_filename);
grad = read_siemens_grad_file(grad_file);

%% Parse a Siemens .grad file
grad_info = grad.info;
R0        = grad.R0;      % radius of a gradient coil [m]
Gref      = grad.lnorm;   % reference gradient strength [mT/m]
alpha_z   = grad.alpha_z; % spherical harmonic coefficients (cosine) for the z gradient coil
alpha_x   = grad.alpha_x; % spherical harmonic coefficients (cosine) for the x gradient coil
alpha_y   = grad.alpha_y; % spherical harmonic coefficients (cosine) for the y gradient coil
beta_z    = grad.beta_z;  % spherical harmonic coefficients (sine)   for the z gradient coil
beta_x    = grad.beta_x;  % spherical harmonic coefficients (sine)   for the x gradient coil
beta_y    = grad.beta_y;  % spherical harmonic coefficients (sine)   for the y gradient coil

%% Calculate the number of total slices
nr_recons = nr_slices;

%% Prepare k-space data per slice
for idx = 1:nr_recons

    %% Get information about the current slice
    slice_number = ind2sub(nr_slices, idx);

    %% Define variables
    %----------------------------------------------------------------------
    % Get a list of imaging acquisitions for the current slice
    %----------------------------------------------------------------------
    img_acq_list = find((img_data.head.idx.slice == (slice_number - 1)));

    %----------------------------------------------------------------------
    % Get a list of segments (imaging)
    %----------------------------------------------------------------------
    img_segment_list = double(img_data.head.idx.segment(img_acq_list));

    %% Get "Cartesian-format" "imaging" k-space data (Nkx x Nky x Nkz x Nc)
    %----------------------------------------------------------------------
    % ( 1)not used   x ( 2)readout dimension x ( 3)number of TRs  x ( 4)COIL_DIM
    % ( 5)MAPS_DIM   x ( 6)TE_DIM            x ( 7)COEFF_DIM      x ( 8)COEFF2_DIM
    % ( 9)ITER_DIM   x (10)CSHIFT_DIM        x (11)TIME_DIM       x (12)TIME2_DIM
    % (13)LEVEL_DIM  x (14)SLICE_DIM         x (15)AVG_DIM        x (16)BATCH_DIM
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    % Calculate the index range (1st index) of a readout acquisition
    %----------------------------------------------------------------------
    kx_range = (discard_pre+1:number_of_samples-discard_post).' - (center_sample + 1) + floor(Nkx/2) + 1;

    %----------------------------------------------------------------------
    % Preallocate k-space data (Nkx x Nky x Nkz x Nc)
    %----------------------------------------------------------------------
    ksp = complex(zeros(Nkx, Nky, Nkz, Nc, 'single'));

    %----------------------------------------------------------------------
    % Sort one acquisition at a time
    %----------------------------------------------------------------------
    nr_acquisitions = length(img_acq_list);

    for i = 1:nr_acquisitions
        tstart = tic; fprintf('%s:(SLC=%2d/%2d) Reading "imaging" k-space data (%3d/%3d)... ', datetime, slice_number, nr_slices, i, nr_acquisitions);

        %------------------------------------------------------------------
        % Get the phase encoding line number
        %------------------------------------------------------------------
        kspace_encode_step_1 = double(img_data.head.idx.kspace_encode_step_1(img_acq_list(i)));
        phase_encoding_line_number = kspace_encode_step_1 - (kspace_encoding_step_1_maximum - kspace_encoding_step_1_center + 1) + floor(Nky/2) + 1;

        %------------------------------------------------------------------
        % Get the partition encoding number
        %------------------------------------------------------------------
        kspace_encode_step_2 = double(img_data.head.idx.kspace_encode_step_2(img_acq_list(i)));
        partition_encoding_number = kspace_encode_step_2 - (kspace_encoding_step_2_maximum - kspace_encoding_step_2_center + 1) + floor(Nkz/2) + 1;
        if partition_encoding_number == 0, partition_encoding_number = partition_encoding_number + 1; end % For 2D imaging

        %------------------------------------------------------------------
        % Prewhiten k-space data
        %------------------------------------------------------------------
        readout = img_data.data{img_acq_list(i)}; % number_of_samples x Nc
        readout = (inv_L * readout.').';

        %------------------------------------------------------------------
        % Accumulate k-space
        %------------------------------------------------------------------
        ksp(kx_range,phase_encoding_line_number,partition_encoding_number,:) = ksp(kx_range,phase_encoding_line_number,partition_encoding_number,:) + reshape(readout, [Nkx 1 1 Nc]);
        fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
    end

    %----------------------------------------------------------------------
    % Divide k-space by the number of averages
    %----------------------------------------------------------------------
    ksp = ksp / nr_averages;

    %% Flip the RO direction of k-space data (Nkx x Nky x Nkz x Nc)
    tstart = tic; fprintf('%s:(SLC=%2d/%2d) Flipping the RO direction of k-space data... ', datetime, slice_number, nr_slices);
    if read_sign == -1
        ksp = flip(ksp,1);
    end
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %% Flip the PE direction of k-space data (Nkx x Nky x Nkz x Nc)
    tstart = tic; fprintf('%s:(SLC=%2d/%2d) Flipping the PE direction of k-space data... ', datetime, slice_number, nr_slices);
    if phase_sign == -1
        ksp = flip(ksp,2);
    end
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %% Define filenames
    ksp_filename = sprintf('ksp_slc%d', slice_number);

    %% Write a .cfl file
    %----------------------------------------------------------------------
    % ksp (Nkx x Nky x Nkz x Nc)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path, ksp_filename);
    tstart = tic; fprintf('%s: Writing a .cfl file: %s... ', datetime, cfl_file);
    writecfl(cfl_file, ksp);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %% Create a sampling mask
    clf_file = fullfile(output_path, sprintf('mask_slc%d', slice_number)); % Nkx x Nky x Nkz
    if ~exist(strcat(clf_file, '.cfl'), 'file')
        ksp_file = strcat(bart_output_path, ksp_filename); % Nkx x Nky x Nkz x Nc
        pat_file = strcat(bart_output_path, sprintf('mask_slc%d', slice_number)); % Nkx x Nky x Nkz
        command = sprintf('%s pattern %s %s', bart_command, ksp_file, pat_file);
        tstart = tic; fprintf('%s:[BART] Creating a sampling mask:\n%s\n', datetime, command);
        [status_pattern,result_pattern] = system(command);
        fprintf('%s: done! (%6.4f/%6.4f sec)\n', datetime, toc(tstart), toc(start_time));
    end

    %% Calculate coil images
    cimg = 1 / sqrt(Nkx * Nky) * fftshift(fft2(ifftshift(ksp)));

    %% Display coil images
    if 0
    nr_rows = 3;
    nr_cols = 7;
    cimg_montage = complex(zeros(Nkx * nr_rows, Nky * nr_cols, 'single'));

    count = 0;
    for idx1 = 1:nr_rows
        for idx2 = 1:nr_cols
            idx1_range = (1:Nkx).' + (idx1 - 1) * Nkx;
            idx2_range = (1:Nky).' + (idx2 - 1) * Nky;
            count = count + 1;
            cimg_montage(idx1_range,idx2_range) = cimg(:,:,1,count);
            if count >= Nc
                break;
            end
        end
    end

    figure('Color', 'k', 'Position', [1 5 1239 973]);
    imagesc(abs(cimg_montage));
    axis image off;
    colormap(gray(256));
    clim([0 5]);
    title({sprintf('FFT recon, SLC = %d', slice_number)}, 'Color', 'w', 'Interpreter', 'latex', 'FontWeight', 'normal', 'FontSize', 16);
    cimg_filename = sprintf('cimg_img_fft_slc%d_mag', slice_number);
    export_fig(fullfile(output_path, cimg_filename), '-r300', '-tif', '-c[50, 400, 310, 540]'); % [top,right,bottom,left]
    close gcf;

    figure('Color', 'k', 'Position', [1 5 1239 973]);
    imagesc(angle(cimg_montage) * 180 / pi);
    axis image off;
    colormap(hsv(256));
    clim([-180 180]);
    title({sprintf('FFT recon, SLC = %d', slice_number)}, 'Color', 'w', 'Interpreter', 'latex', 'FontWeight', 'normal', 'FontSize', 16);
    cimg_filename = sprintf('cimg_img_fft_slc%d_phase', slice_number);
    export_fig(fullfile(output_path, cimg_filename), '-r300', '-tif', '-c[50, 400, 310, 540]'); % [top,right,bottom,left]
    close gcf;
    end

    %% Calculate the actual slice number for Siemens interleaved multislice imaging
    %----------------------------------------------------------------------
    % slice_number: acquisition slice number
    % actual_slice_number is used to retrieve slice information from TWIX format
    % Q: For coronal? acq_slice_nr = slice_nr;
    %----------------------------------------------------------------------
    if nr_slices > 1 % multi-slice
        if mod(nr_slices,2) == 0 % even
            offset1 = 0;
            offset2 = 1;
        else % odd
            offset1 = 1;
            offset2 = 0;
        end
        if slice_number <= ceil(nr_slices / 2)
            actual_slice_number = 2 * slice_number - offset1;
        else
            actual_slice_number = 2 * (slice_number - ceil(nr_slices / 2)) - offset2;
        end
    else
        actual_slice_number = slice_number;
    end

    %% Get a slice offset in the PCS from Siemens TWIX format
    if isfield(twix.hdr.MeasYaps.sSliceArray.asSlice{actual_slice_number}, 'sPosition')
        if isfield(twix.hdr.MeasYaps.sSliceArray.asSlice{actual_slice_number}.sPosition, 'dSag')
            sag_offset_twix = twix.hdr.MeasYaps.sSliceArray.asSlice{actual_slice_number}.sPosition.dSag; % [mm]
        else
            sag_offset_twix = 0; % [mm]
        end
        if isfield(twix.hdr.MeasYaps.sSliceArray.asSlice{actual_slice_number}.sPosition, 'dCor')
            cor_offset_twix = twix.hdr.MeasYaps.sSliceArray.asSlice{actual_slice_number}.sPosition.dCor; % [mm]
        else
            cor_offset_twix = 0; % [mm]
        end
        if isfield(twix.hdr.MeasYaps.sSliceArray.asSlice{actual_slice_number}.sPosition, 'dTra')
            tra_offset_twix = twix.hdr.MeasYaps.sSliceArray.asSlice{actual_slice_number}.sPosition.dTra; % [mm]
        else
            tra_offset_twix = 0; % [mm]
        end
    else
        sag_offset_twix = 0; % [mm]
        cor_offset_twix = 0; % [mm]
        tra_offset_twix = 0; % [mm]
    end

    %% Get a slice offset of a stack in the PCS from Siemens TWIX format
    pcs_offset = [sag_offset_twix; cor_offset_twix; tra_offset_twix] * 1e-3; % [mm] * [m/1e3mm] => [m]

    %% Get a slice offset in the PCS from ISMRMRD format
    sag_offset_ismrmrd = double(img_data.head.position(1,img_acq_list(1))); % [mm]
    cor_offset_ismrmrd = double(img_data.head.position(2,img_acq_list(1))); % [mm]
    tra_offset_ismrmrd = double(img_data.head.position(3,img_acq_list(1))); % [mm]
    pcs_offset_ismrmrd = [sag_offset_ismrmrd; cor_offset_ismrmrd; tra_offset_ismrmrd] * 1e-3; % [mm] * [m/1e3mm] => [m]

    %% Get a rotation matrix from the GCS to the PCS (ISMRMRD format)
    phase_dir = double(img_data.head.phase_dir(:,img_acq_list(1)));
    read_dir  = double(img_data.head.read_dir(:,img_acq_list(1)));
    slice_dir = double(img_data.head.slice_dir(:,img_acq_list(1)));
    R_gcs2pcs_ismrmrd = [phase_dir read_dir slice_dir];

    %% Calculate a slice offset in the DCS [m]
    dcs_offset = R_pcs2dcs * pcs_offset; % 3 x 1

    %% Display slice information
    fprintf('======================= SLICE INFORMATION ========================\n');
    fprintf('main_orientation = %d (SAGITTAL/CORONAL/TRANSVERSAL = 0/1/2)\n', main_orientation);
    fprintf('slice_nr = %d, actual_slice_nr = %d\n', slice_number, actual_slice_number);
    fprintf('dNormalSag = %+g \ndNormalCor = %+g \ndNormalTra = %+g \ndRotAngle = %g [rad]\n', dNormalSag, dNormalCor, dNormalTra, dRotAngle);
    fprintf('phase_sign = %+g, read_sign = %+g\n', phase_sign, read_sign);
    fprintf('---------------------- From Siemens TWIX format ------------------\n');
    fprintf('                   [sag]   %10.5f [mm]\n', sag_offset_twix);
    fprintf('slice offset(PCS): [cor] = %10.5f [mm]\n', cor_offset_twix);
    fprintf('                   [tra]   %10.5f [mm]\n', tra_offset_twix);
    fprintf('---------------------- From ISMRMRD format -----------------------\n');
    fprintf('                   [sag]   %10.5f [mm]\n', sag_offset_ismrmrd);
    fprintf('slice offset(PCS): [cor] = %10.5f [mm]\n', cor_offset_ismrmrd);
    fprintf('                   [tra]   %10.5f [mm]\n', tra_offset_ismrmrd);
    fprintf('-------------- From Siemens TWIX format (w/ gradient flips) ------\n');
    fprintf('                   [sag]   [%10.5f %10.5f %10.5f][PE]\n', R_gcs2pcs(1,1), R_gcs2pcs(1,2), R_gcs2pcs(1,3));
    fprintf('R_gcs2pcs        : [cor] = [%10.5f %10.5f %10.5f][RO]\n', R_gcs2pcs(2,1), R_gcs2pcs(2,2), R_gcs2pcs(2,3));
    fprintf('                   [tra]   [%10.5f %10.5f %10.5f][SL]\n', R_gcs2pcs(3,1), R_gcs2pcs(3,2), R_gcs2pcs(3,3));
    fprintf('-------------- From ISMRMRD format (w/o gradient flips) ----------\n');
    fprintf('                   [sag]   [%10.5f %10.5f %10.5f][PE]\n', R_gcs2pcs_ismrmrd(1,1), R_gcs2pcs_ismrmrd(1,2), R_gcs2pcs_ismrmrd(1,3));
    fprintf('R_gcs2pcs        : [cor] = [%10.5f %10.5f %10.5f][RO]\n', R_gcs2pcs_ismrmrd(2,1), R_gcs2pcs_ismrmrd(2,2), R_gcs2pcs_ismrmrd(2,3));
    fprintf('                   [tra]   [%10.5f %10.5f %10.5f][SL]\n', R_gcs2pcs_ismrmrd(3,1), R_gcs2pcs_ismrmrd(3,2), R_gcs2pcs_ismrmrd(3,3));
    fprintf('------------------------------------------------------------------\n');
    fprintf('                   [ x ]   [%10.5f %10.5f %10.5f][sag]\n', R_pcs2dcs(1,1), R_pcs2dcs(1,2), R_pcs2dcs(1,3));
    fprintf('R_pcs2dcs        : [ y ] = [%10.5f %10.5f %10.5f][cor]\n', R_pcs2dcs(2,1), R_pcs2dcs(2,2), R_pcs2dcs(2,3));
    fprintf('                   [ z ]   [%10.5f %10.5f %10.5f][tra]\n', R_pcs2dcs(3,1), R_pcs2dcs(3,2), R_pcs2dcs(3,3));
    fprintf('------------------------------------------------------------------\n');
    fprintf('                   [ x ]   [%10.5f %10.5f %10.5f][PE]\n', R_gcs2dcs(1,1), R_gcs2dcs(1,2), R_gcs2dcs(1,3));
    fprintf('R_gcs2dcs        : [ y ] = [%10.5f %10.5f %10.5f][RO]\n', R_gcs2dcs(2,1), R_gcs2dcs(2,2), R_gcs2dcs(2,3));
    fprintf('                   [ z ]   [%10.5f %10.5f %10.5f][SL]\n', R_gcs2dcs(3,1), R_gcs2dcs(3,2), R_gcs2dcs(3,3));
    fprintf('==================================================================\n');

    %% Read a .cfl file
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

    %% Calculate a circular mask or a cube mask
    %----------------------------------------------------------------------
    % The correction is limited to a 2 * r0 cube enclosing the specified FOV.
    % Pixels outside the 2 * r0 FOV are set to zero for safety reasons.
    %----------------------------------------------------------------------
    circle_mask = zeros(Nkx, Nky, Nkz, 'single');
    %circle_mask((x.^2 + y.^2 + z.^2) < R0^2) = 1;
    circle_mask((abs(x) < R0) & (abs(y) < R0) & (abs(z) < R0)) = 1;

    %% Convert Cartesian coordinates to spherical coordinates
    tstart = tic; fprintf('%s:(SLC=%2d/%2d) Converting Cartesian coordinates (x,y,z) to spherical coordinates (r,theta,phi)... ', datetime, slice_number, nr_slices);
    radius = reshape(sqrt(x.^2 + y.^2 + z.^2), [N 1]); % N x 1 [m]
    theta = reshape(atan2(sqrt(x.^2 + y.^2), z), [N 1]); % same as acos(z./r)
    phi = reshape(atan2(y, x), [N 1]); % azimuth
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %% Calculate the z component of the nonlinear magnetic field generated by the x gradient coil [mT]
    tstart = tic; fprintf('%s:(SLC=%2d/%2d) Calculating the z component of the nonlinear magnetic field generated by the x gradient coil... ', datetime, slice_number, nr_slices);
    Bz_xgrad_nonlinear = siemens_B(radius, theta, phi, R0, alpha_x, beta_x); % N x 1 [mT]
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %% Calculate the z component of the nonlinear magnetic field generated by the y gradient coil [mT]
    tstart = tic; fprintf('%s:(SLC=%2d/%2d) Calculating the z component of the nonlinear magnetic field generated by the y gradient coil... ', datetime, slice_number, nr_slices);
    Bz_ygrad_nonlinear = siemens_B(radius, theta, phi, R0, alpha_y, beta_y); % N x 1 [mT]
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %% Calculate the z component of the nonlinear magnetic field generated by the z gradient coil [mT]
    tstart = tic; fprintf('%s:(SLC=%2d/%2d) Calculating the z component of the nonlinear magnetic field generated by the z gradient coil... ', datetime, slice_number, nr_slices);
    Bz_zgrad_nonlinear = siemens_B(radius, theta, phi, R0, alpha_z, beta_z); % N x 1 [mT]
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %% Calculate displacement fields in the DCS [m] (3 x N)
    % [?] / [mT/m] = [m] => [?] = [m] * [mT/m] => [mT]
    dx = reshape(Bz_xgrad_nonlinear / Gref, [Nkx Nky Nkz]); % [mT] / [mT/m] => [m]
    dy = reshape(Bz_ygrad_nonlinear / Gref, [Nkx Nky Nkz]); % [mT] / [mT/m] => [m]
    dz = reshape(Bz_zgrad_nonlinear / Gref, [Nkx Nky Nkz]); % [mT] / [mT/m] => [m]
    clear Bz_xgrad_nonlinear Bz_ygrad_nonlinear Bz_zgrad_nonlinear;

    %% Calculate displacement fields in the GCS [m] (3 x N) [v,u,w] = [PE,RO,SL]
    dr_gcs = R_gcs2dcs.' * cat(2, dx(:), dy(:), dz(:)).'; % 3 x N
    du = reshape(dr_gcs(2,:), [Nkx Nky Nkz]); % Nkx x Nky x Nkz [m] (RO)
    dv = reshape(dr_gcs(1,:), [Nkx Nky Nkz]); % Nkx x Nky x Nkz [m] (PE)
    dw = reshape(dr_gcs(3,:), [Nkx Nky Nkz]); % Nkx x Nky x Nkz [m] (SL)

    %% Write a .cfl file
    %----------------------------------------------------------------------
    % circle_mask (Nkx x Nky x Nkz)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path, sprintf('circle_mask_slc%d_%s', slice_number, slice_type));
    tstart = tic; fprintf('%s: Writing a .cfl file: %s... ', datetime, cfl_file);
    writecfl(cfl_file, circle_mask);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % dx (Nkx x Nky x Nkz)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path, sprintf('dx_slc%d_%s', slice_number, slice_type));
    tstart = tic; fprintf('%s: Writing a .cfl file: %s... ', datetime, cfl_file);
    writecfl(cfl_file, dx);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % dy (Nkx x Nky x Nkz)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path, sprintf('dy_slc%d_%s', slice_number, slice_type));
    tstart = tic; fprintf('%s: Writing a .cfl file: %s... ', datetime, cfl_file);
    writecfl(cfl_file, dy);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % dz (Nkx x Nky x Nkz)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path, sprintf('dz_slc%d_%s', slice_number, slice_type));
    tstart = tic; fprintf('%s: Writing a .cfl file: %s... ', datetime, cfl_file);
    writecfl(cfl_file, dz);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % du (Nkx x Nky x Nkz)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path, sprintf('du_slc%d_%s', slice_number, slice_type));
    tstart = tic; fprintf('%s: Writing a .cfl file: %s... ', datetime, cfl_file);
    writecfl(cfl_file, du);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % dv (Nkx x Nky x Nkz)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path, sprintf('dv_slc%d_%s', slice_number, slice_type));
    tstart = tic; fprintf('%s: Writing a .cfl file: %s... ', datetime, cfl_file);
    writecfl(cfl_file, dv);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % dw (Nkx x Nky x Nkz)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path, sprintf('dw_slc%d_%s', slice_number, slice_type));
    tstart = tic; fprintf('%s: Writing a .cfl file: %s... ', datetime, cfl_file);
    writecfl(cfl_file, dw);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
end
