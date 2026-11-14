% 1. Define the base path and the names of your two data folders
base_path = '/home/khushbu-lab1/Tanmat/';
folder_names = {'Dizygotic twins mat (2)/mat', 'Monozygotic Twins mat/mat'}; % <-- REPLACE 'Another_Data_Folder_Name'

all_subjects = []; % Initialize the final matrix

% 2. Loop through each data folder
for f = 1:length(folder_names)
    
    current_folder = fullfile(base_path, folder_names{f});
    fprintf('Processing files in: %s\n', current_folder); % Display current folder
    
    % Get all .mat files in the current folder
    files = dir(fullfile(current_folder, '*.mat'));
    
    % 3. Loop through each file in the current folder
    for i = 1:length(files)
        fname = files(i).name;
        
        % Skip hidden/system files (start with "._")
        if startsWith(fname, '._')
            continue
        end
        
        full_file_path = fullfile(current_folder, fname);
        
        % Load file into structure 'd'
        d = load(full_file_path);
        
        % Check for the confirmed variable name 'ts_data'
        if isfield(d, 'ts_data')
            % Extract your matrix
            ts = d.ts_data;
            
            % Concatenate
            all_subjects = [all_subjects; ts];
        else
            warning('File %s does not contain the variable "ts_data". Skipping.', fname);
            % (Optional: You could add logic here to check for an alternate variable name)
        end
    end
end

% 4. Save the final concatenated data
save(fullfile(base_path, 'all_subjects_twins.mat'), 'all_subjects');
fprintf('Data concatenation complete. Saved to all_subjects_twins.mat.\n');