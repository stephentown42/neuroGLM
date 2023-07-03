function settings = load_settings(file_path, file_name)
% function settings = load_settings(file_path, file_name)
%
% Load settings for basis functions from file

% JSON
if contains(file_name,'.json')
    json_txt = fileread(fullfile( file_path, file_name));
    settings = jsondecode(json_txt);
    return
end

% Add support for additional file formats below