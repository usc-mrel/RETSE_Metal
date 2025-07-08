%% calculate_offres.m
% Combined script to read GNL-corrected .cfl per slie and estimate
% off-resonance displacmenet along readout direciton

function calculate_offres(opts)
addpath(genpath('./Helper_Functions'));

%% Step 1: Locate and read .cfl slices from /result folders
root_path = pwd;
result_root = fullfile(root_path, 'result');

% List result subfolders, excluding hidden and dualrecon_pcs
folder_info = dir(result_root);
all_folders = folder_info([folder_info.isdir] & ~startsWith({folder_info.name}, '.') & ...
    ~strcmp({folder_info.name}, 'dualrecon_pcs'));
result_folders = {all_folders.name};

% Require at least 2 valid folders
if numel(result_folders) < 2
    error('At least two valid result folders are required, found %d.', numel(result_folders));
end

% Only use the first two non-dualrecon_pcs folders
fprintf('Using folders for fplus and fminus:\n 1. %s\n 2. %s\n', result_folders{1}, result_folders{2});

gnl_volumes = cell(1, 2);

for idx_folder = 1:2
    cfl_folder = fullfile(result_root, result_folders{idx_folder});
    cfl_files = dir(fullfile(cfl_folder, 'img_type1_slc*_gnl1_offres0*.cfl'));

    assert(~isempty(cfl_files), 'No matching .cfl files found in %s', cfl_folder);

    % Extract and sort slice numbers
    slice_numbers = zeros(1, numel(cfl_files));
    for k = 1:numel(cfl_files)
        tokens = regexp(cfl_files(k).name, 'slc(\d+)', 'tokens');
        assert(~isempty(tokens), 'Cannot parse slice number from %s', cfl_files(k).name);
        slice_numbers(k) = str2double(tokens{1}{1});
    end
    [~, sort_idx] = sort(slice_numbers);
    cfl_files = cfl_files(sort_idx);

    % Load all slices
    for slice_idx = 1:numel(cfl_files)
        [~, prefix, ~] = fileparts(cfl_files(slice_idx).name);
        slice_img = readcfl(fullfile(cfl_folder, prefix));
        if slice_idx == 1
            [Nx, Ny] = size(slice_img);
            Ns = numel(cfl_files);
            vol = zeros(Nx, Ny, Ns, 'single');
        end
        vol(:,:,slice_idx) = single(slice_img);
    end

    gnl_volumes{idx_folder} = abs(vol);  % Use abs() as required for estimation
    fprintf('Loaded %s with size [%d %d %d]\n', result_folders{idx_folder}, size(vol));
end

fplus = abs(reorder_slice_3rd(gnl_volumes{1},'int',numel(cfl_files))); % Convert to nominal slice order
fminus = abs(reorder_slice_3rd(gnl_volumes{2},'int',numel(cfl_files)));

%% Step 3: offfreq B0 Estimation
fprintf('Estimating B0 field...\n');
tic;
[Tx, cost] = calc_bfield(fplus, fminus, opts);
toc;
Tx = gather(Tx);

%% Step 4: Invert distortion (voxel-wise kernel, you can comment here)
% fprintf('Solving for distortion-corrected image...\n');
[Nkx, Nky, nr_slices] = size(Tx);
recon = zeros(Nkx, Nky, nr_slices, 'like', fplus);
x = 1:Nkx;

tic;
parfor iy = 1:Nky
    for iz = 1:nr_slices
        t1 = x.' - Tx(:,iy,iz);
        kplus = sinc(t1.' - x.');

        t2 = x.' + Tx(:,iy,iz);
        kminus = sinc(t2.' - x.');

        ktot = [kplus; kminus];
        recon(:,iy,iz) = ktot \ [fplus(:,iy,iz); fminus(:,iy,iz)];
    end
end
toc;

save(sprintf('recon_%d',opts.rBW),'recon');
fprintf(' Finished Off-Resonance Displacement Estimation and Image-Doamin Correction.\n');

%% Step 5: Writ into cfl files
b0_output_paths = fullfile(pwd, 'offres');

if ~exist(b0_output_paths, 'dir')
    mkdir(b0_output_paths);
end
%% Save Tx into cfl files
start_time = tic;
% Coronal only!
%reorient = @(x) rot90(x,1); % counterclockwise rotation
reorient = @(x) x;

for slice_number = 1:nr_slices

    %-------------------------------------------------    %------------------------------------------------------------------
    % Calculate the actual slice number for Siemens interleaved multislice imaging
    %------------------------------------------------------------------
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

    Tx_save = complex(zeros(Nkx, Nky, 'single'));
    Tx_save = reorient(Tx(:,:,actual_slice_number));

    %------------------------------------------------------------------
    % Write as a .cfl file-----------------
    slice_type = 'flat';
    cfl_file = fullfile(b0_output_paths, sprintf('fieldmap_slc%d_%s', slice_number, slice_type));
    tstart = tic; fprintf('%s: Writing a .cfl file: %s... ', datetime, cfl_file);
    writecfl(cfl_file, Tx_save / 150); % [Hz]
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %------------------------------------------------------------------
    % Write as a .cfl file
    %------------------------------------------------------------------
    cfl_file = fullfile(b0_output_paths, sprintf('displacement_slc%d_%s', slice_number, slice_type));
    tstart = tic; fprintf('%s: Writing a .cfl file: %s... ', datetime, cfl_file);
    writecfl(cfl_file, Tx_save * opts.xres * 1e-3); % [Hz] / [Hz/mm] * [m/1e3mm] => [m]
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
end
end