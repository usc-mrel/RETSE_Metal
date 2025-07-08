function auto_convert_twix_to_ismrmrd()
% AUTO_CONVERT_TWIX_TO_ISMRMRD Convert .dat files under /raw/*/ to ISMRMRD .h5 files.
%
% Skips conversion if .h5 files already exist.

%% Add /usr/local/bin to PATH (for siemens_to_ismrmrd)
setenv('PATH', [getenv('PATH') ':/usr/local/bin']);

%% Define raw root directory
raw_root = fullfile(pwd, 'raw');

%% Find all subfolders
folder_info = dir(raw_root);
folders = {folder_info([folder_info.isdir] & ~startsWith({folder_info.name}, '.')).name};

%% Start conversion
start_time = tic;

for idx_folder = 1:numel(folders)
    folder = folders{idx_folder};
    folder_path = fullfile(raw_root, folder);
    fprintf('Processing folder: %s\n', folder_path);
    
    % Find .dat files inside this folder
    dat_files = dir(fullfile(folder_path, '*.dat'));
    
    if isempty(dat_files)
        warning('No .dat files found in %s, skipping...', folder_path);
        continue;
    end
    
    for idx_dat = 1:numel(dat_files)
        dat_filename = dat_files(idx_dat).name;
        [~, prefix, ~] = fileparts(dat_filename); % get file name without extension
        
        noise_h5_file = fullfile(folder_path, ['noise_', prefix, '.h5']);
        data_h5_file = fullfile(folder_path, [prefix, '.h5']);
        
        % Check if both .h5 files already exist
        if isfile(noise_h5_file) && isfile(data_h5_file)
            fprintf('Both noise and data .h5 files already exist for %s, skipping conversion.\n', prefix);
            continue; % Skip to next .dat
        end
        
        % --- Noise conversion (if missing) ---
        if ~isfile(noise_h5_file)
            linux_command1 = sprintf('siemens_to_ismrmrd -f "%s" -z 1 -o "%s"', ...
                fullfile(folder_path, dat_filename), noise_h5_file);
            
            fprintf('Running noise conversion: %s\n', linux_command1);
            tic;
            [status1, result1] = system(linux_command1);
            if status1 ~= 0
                warning('Noise conversion failed: %s\n%s', linux_command1, result1);
            else
                fprintf('Noise conversion done (%.2f sec)\n', toc);
            end
        else
            fprintf('Noise .h5 file already exists: %s\n', noise_h5_file);
        end
        
        % --- K-space conversion (if missing) ---
        if ~isfile(data_h5_file)
            linux_command2 = sprintf('siemens_to_ismrmrd -f "%s" -z 2 -o "%s"', ...
                fullfile(folder_path, dat_filename), data_h5_file);
            
            fprintf('Running k-space conversion: %s\n', linux_command2);
            tic;
            [status2, result2] = system(linux_command2);
            if status2 ~= 0
                warning('K-space conversion failed: %s\n%s', linux_command2, result2);
            else
                fprintf('K-space conversion done (%.2f sec)\n', toc);
            end
        else
            fprintf('K-space .h5 file already exists: %s\n', data_h5_file);
        end
        
    end
end

fprintf('All TWIX to ISMRMRD conversion completed in %.2f seconds.\n', toc(start_time));

end
