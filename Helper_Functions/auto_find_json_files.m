function json_files = auto_find_json_files()
% AUTO_FIND_JSON_FILES Automatically find individual JSON files under /raw/
% 
% OUTPUT:
%   json_files - a cell array of full paths to each dataset's JSON file
%
% This function:
% - Scans the /raw/ folder
% - Finds subfolders (e.g., PE_FH, PE_HF)
% - Finds .json file in each (ignores multi.json)
% - Returns full paths as a cell array

%% Define raw path relative to current folder
raw_root = fullfile(pwd, 'raw');

%% Find all subfolders
folder_info = dir(raw_root);
folders = {folder_info([folder_info.isdir] & ~startsWith({folder_info.name}, '.')).name};

%% Initialize output
json_files = {};

%% Loop through each folder
for idx = 1:numel(folders)
    folder = folders{idx};
    raw_folder = fullfile(raw_root, folder);
    
    % Find .json file inside folder
    json_list = dir(fullfile(raw_folder, '*.json'));
    
    % Exclude multi.json if it exists
    json_list = json_list(~contains({json_list.name}, 'multi.json'));
    
    if isempty(json_list)
        warning('No JSON file found in %s', raw_folder);
        continue;
    elseif numel(json_list) > 1
        warning('Multiple JSON files found in %s, using the first one.', raw_folder);
    end
    
    % Save full path
    json_files{idx} = fullfile(raw_folder, json_list(1).name);
end

%% Display detected json files
fprintf('Detected %d individual JSON files:\n', numel(json_files));
disp(json_files);

end
