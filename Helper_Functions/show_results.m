%% Clean slate
clear all; clc;
start_time = tic;
%% Read a .json file
json_files = auto_find_json_files();
nr_json_files = length(json_files);
for json_number = 1:nr_json_files
    json_file = json_files{json_number};
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
        output_path        = strrep(json.output_path, '/', '\');

    else
        output_path        = json.output_path;
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

    %% Write a .cfl file
    %----------------------------------------------------------------------
    % img (Nkx x Nky x Nkz)
    %----------------------------------------------------------------------
    %% Calculate the number of total slices
    nr_recons = nr_slices * nr_repetitions;

    %% Read two directions image per slice
    for idx = 1:nr_recons
        %% Get information about the current slice
        [slice_number, repetition_number] = ind2sub([nr_slices nr_repetitions], idx);
        
        % Load uncorrected image
        img_filename = sprintf('img_type1_slc%d_rep%d_%s_gridding%d_phc%d_conc%d_gnl%d_offres%d_i%d_l%4.2f', slice_number, repetition_number, slice_type, readout_gridding_flag, phase_correction_flag, conc_correction_flag, gnl_correction_flag, 0, maxiter, lambda);
        cfl_file = fullfile(output_path, img_filename);
        tstart = tic; fprintf('Reading a .cfl file: %s... ', cfl_file);
        img_recon_temp(:,:,idx,2*json_number-1) = readcfl(cfl_file);
        fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

        % Load correcgted image from one direction
        img_filename = sprintf('img_type1_slc%d_rep%d_%s_gridding%d_phc%d_conc%d_gnl%d_offres%d_i%d_l%4.2f', slice_number, repetition_number, slice_type, readout_gridding_flag, phase_correction_flag, conc_correction_flag, gnl_correction_flag, 1, maxiter, lambda);
        cfl_file = fullfile(output_path, img_filename);
        tstart = tic; fprintf('Reading a .cfl file: %s... ', cfl_file);
        img_recon_temp(:,:,idx,2*json_number) = readcfl(cfl_file);
        fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
    end
end

%% Read the corrected recon image
json_file = auto_find_multi_json();

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
        output_path        = strrep(json.output_path, '/', '\');

    else
        output_path        = json.output_path;
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

    %% Write a .cfl file
    %----------------------------------------------------------------------
    % img (Nkx x Nky x Nkz)
    %----------------------------------------------------------------------
    %% Calculate the number of total slices
    nr_recons = nr_slices * nr_repetitions;

    %% Read two directions image per slice
    for idx = 1:nr_recons
        %% Get information about the current slice
        [slice_number, repetition_number] = ind2sub([nr_slices nr_repetitions], idx);

        img_filename = sprintf('img_type1_slc%d_rep%d_%s_gridding%d_phc%d_conc%d_gnl%d_offres%d_i%d_l%4.2f', slice_number, repetition_number, slice_type, readout_gridding_flag, phase_correction_flag, conc_correction_flag, gnl_correction_flag, 1, maxiter, lambda);
        cfl_file = fullfile(output_path, img_filename);
        tstart = tic; fprintf('Reading a .cfl file: %s... ', cfl_file);
        img_recon_temp(:,:,idx,5) = readcfl(cfl_file);
        fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
    end

%% Re-order to correct slice
for nimg = 1:size(img_recon_temp,4)
    img_recon(:,:,:,nimg) =reorder_slice_3rd(img_recon_temp(:,:,:,nimg), 'int', nr_slices);
end
