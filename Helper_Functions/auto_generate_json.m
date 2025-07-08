function auto_generate_json(varargin)
% AUTO_GENERATE_offres_JSON Generate JSON files for MRI recon using offres fieldmap.
%
% Example usage:
%    auto_generate_offres_json('offres_correction_flag', [1 1], 'offres_sign', [-1 1], 'nr_slices', 20)

%% Parse input arguments
p = inputParser;
addParameter(p, 'offres_correction_flag', []);
addParameter(p, 'offres_sign', []);
addParameter(p, 'nr_slices', 20);
addParameter(p, 'bart_path', []);
parse(p, varargin{:});

offres_correction_flag = p.Results.offres_correction_flag;
offres_sign = p.Results.offres_sign;
nr_slices = p.Results.nr_slices;

%% Define root paths
root_path = pwd;
raw_root = fullfile(root_path, 'raw');
bart_path = p.Results.bart_path;
offres_path = fullfile(root_path, 'offres');

%% Find subfolders under /raw/
folder_info = dir(raw_root);
folders = {folder_info([folder_info.isdir] & ~startsWith({folder_info.name}, '.')).name};

if isempty(offres_correction_flag) || isempty(offres_sign)
    error('Both offres_correction_flag and offres_sign must be provided.');
end

if numel(offres_correction_flag) ~= numel(folders) || numel(offres_sign) ~= numel(folders)
    error('offres_correction_flag and offres_sign must match the number of folders under /raw/');
end

%% Prepare multi.json
multi_json = struct();

%% Loop through folders to generate JSON
for idx = 1:numel(folders)
    folder = folders{idx};
    raw_folder = fullfile(raw_root, folder);

    % Find .dat file
    dat_files = dir(fullfile(raw_folder, '*.dat'));
    if isempty(dat_files)
        warning('No .dat file found in %s, skipping.', raw_folder);
        continue;
    end
    [~, prefix, ~] = fileparts(dat_files(1).name);

    result_folder = fullfile(root_path, 'result', prefix);

    json_struct = struct( ...
        'siemens_twix_file', fullfile(raw_folder, [prefix, '.dat']), ...
        'ismrmrd_data_file', fullfile(raw_folder, [prefix, '.h5']), ...
        'ismrmrd_noise_file', fullfile(raw_folder, ['noise_', prefix, '.h5']), ...
        'output_path', result_folder, ...
        'sens_path', result_folder, ...
        'bart_path', bart_path, ...
        'offres_path', offres_path, ...
        'recon_parameters', struct( ...
            'lambda', 0, ...
            'tol', 1e-5, ...
            'maxiter', 6, ...
            'slice_type', 'flat', ...
            'cal_size', [24, 24, 1], ...
            'remove_oversampling_flag', 0, ...
            'readout_gridding_flag', 0, ...
            'phase_correction_flag', 0, ...
            'gnl_correction_flag', 1, ...
            'conc_correction_flag', 0, ...
            'offres_correction_flag', offres_correction_flag(idx) ...
        ), ...
        'nr_slices', nr_slices, ...
        'main_orientation', 1, ...
        'offres_sign', offres_sign(idx) ...
    );

    % Save individual JSON
    json_text = jsonencode(json_struct, 'PrettyPrint', true);
    json_filename = fullfile(raw_folder, [prefix, '.json']);
    fid = fopen(json_filename, 'w');
    if fid == -1
        error('Cannot create JSON file: %s', json_filename);
    end
    fwrite(fid, json_text, 'char');
    fclose(fid);

    fprintf('Saved JSON: %s (offres_correction_flag = %d, offres_sign = %d)\n', json_filename, offres_correction_flag(idx), offres_sign(idx));

    % Fill multi.json
    multi_json.(['siemens_twix_file_', num2str(idx)]) = fullfile(raw_folder, [prefix, '.dat']);
    multi_json.(['ismrmrd_data_file_', num2str(idx)]) = fullfile(raw_folder, [prefix, '.h5']);
    multi_json.(['ismrmrd_noise_file_', num2str(idx)]) = fullfile(raw_folder, ['noise_', prefix, '.h5']);
    multi_json.(['input_path_', num2str(idx)]) = result_folder;
    multi_json.(['sens_path_', num2str(idx)]) = result_folder;
    multi_json.(['offres_sign_', num2str(idx)]) = offres_sign(idx);
end

% Final fields for multi.json
multi_json.output_path = fullfile(root_path, 'result', 'dualrecon_pcs');
multi_json.bart_path = bart_path;
multi_json.offres_path = offres_path;
multi_json.recon_parameters = struct( ...
    'lambda', 0, ...
    'tol', 1e-5, ...
    'maxiter', 6, ...
    'slice_type', 'flat', ...
    'cal_size', [24, 24, 1], ...
    'remove_oversampling_flag', 0, ...
    'readout_gridding_flag', 0, ...
    'phase_correction_flag', 0, ...
    'gnl_correction_flag', 1, ...
    'conc_correction_flag', 0, ...
    'offres_correction_flag', 1);

multi_json.nr_slices = nr_slices;
multi_json.main_orientation = 1;

% Save multi.json
multi_json_text = jsonencode(multi_json, 'PrettyPrint', true);
multi_json_filename = fullfile(raw_root, 'multi.json');
fid = fopen(multi_json_filename, 'w');
if fid == -1
    error('Cannot create multi.json');
end
fwrite(fid, multi_json_text, 'char');
fclose(fid);

fprintf('Saved multi JSON: %s\n', multi_json_filename);

end