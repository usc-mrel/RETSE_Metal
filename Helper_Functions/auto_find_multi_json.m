function multi_json_file = auto_find_multi_json()
% AUTO_FIND_MULTI_JSON Automatically find multi.json file under raw/ folder
%
% OUTPUT:
%   multi_json_file - full path to the multi.json file

% Assume pwd is at main project level, so raw/ folder exists
raw_folder = fullfile(pwd, 'raw');

if ~exist(raw_folder, 'dir')
    error('Cannot find raw/ folder in current path: %s', raw_folder);
end

% Search for multi.json
multi_json_info = dir(fullfile(raw_folder, 'multi.json'));

if isempty(multi_json_info)
    error('No multi.json file found under raw/.');
elseif length(multi_json_info) > 1
    error('More than one multi.json found under raw/. Please check.');
end

multi_json_file = fullfile(raw_folder, multi_json_info.name);

fprintf('%s: Found multi.json at %s\n', datetime, multi_json_file);
end
