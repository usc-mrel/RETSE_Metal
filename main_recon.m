%% main_recon
% Main pipeline to correct gradient non-linearity and off-resonance distortion

%% Step 0: Clean slate
restoredefaultpath;
close all; clear; clc

%% Set source directories
ismrmrd_path   = ''; % Set your ismrmrd path
bart_path = '';   % Set your bart path
functions_path = './Helper_Functions';

%% Add source directories to search path
addpath(genpath(ismrmrd_path));
addpath(genpath(functions_path));

%% Step 1: Generate and Define JSON files (only GNL correction is on)
auto_generate_json( ...
    'offres_correction_flag', [0 0], ...
    'offres_sign', [1,-1], ...
    'nr_slices', 20,...
    'bart_path', bart_path);

json_files = auto_find_json_files();
nr_json_files = length(json_files);

%% Step 2: Convert TWIX to ISMRMRD (skip if already exists)
auto_convert_twix_to_ismrmrd();

%% Step 3: Perform GNL-only reconstruction on each dataset
for json_number = 1:nr_json_files
    json_file = json_files{json_number};
    
    demo_step1_calculate_voxel_coordinates;
    demo_step2_prepare_ksp;
    demo_step3_estimate_csm;
    demo_step4_cartesian_gnl_recon;
end

%% Step 4: Estimate off-resonance displacement map (Saved in cfl file)
opts = struct();
opts.epsilon = 5e-6;
opts.maxIter = 800;
opts.stag_tol = 6.5e-4;
opts.stag_lim = 100;
opts.stepsize = 5e-5;
opts.xres = 0.875; % readout resolution [mm]
opts.rBW = 150;    % readout bandwidth [Hz/Px]
opts.spacing = [4,4,1];

calculate_offres(opts);

%% Step 5: Re-generate JSON with off-resonance correction enabled
auto_generate_json( ...
    'offres_correction_flag', [1 1], ...
    'offres_sign', [1,-1], ...
    'nr_slices', 20,...
    'bart_path', bart_path);

%% Step 6: Apply GNL + Off-resonance correction on each dataset
json_files = auto_find_json_files();
nr_json_files = length(json_files);
for json_number = 1:2
    json_file = json_files{json_number};
    demo_cartesian_recon;
end

%% Step 7: Perform joint dual-encoding reconstruction (PEC-SENSE)
demo_cartesian_dual_pecsense_recon;

%% Show results 
% img_recon: 5 images (Dir.1; corrected from Dir.1; Dir.2; corrected from Dir. 2; Correctino image from Phase-constrained reconstruciton)
% Show (Dir.1,Dir.2 and Corrected images)
show_results;

[N1,N2,N3] = size(img_recon,[1,2,3]);
rowcrop = 240;
colcrop = 150;
idx1 = (-floor(rowcrop/2):ceil(rowcrop/2)-1).' + floor(N1/2) + 1;
idx2 = (-floor(colcrop/2):ceil(colcrop/2)-1).' + floor(N2/2) + 1;
figure;

imagesc(abs(cat(2,img_recon(idx1,idx2,12,1), img_recon(idx1,idx2,12,3), img_recon(idx1,idx2,12,5))),[-1.5,15.5]);
axis image off
colormap(gray)