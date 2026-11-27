function step05_concatenate_subjects(inputDirs, outputFile)
% STEP05_CONCATENATE_SUBJECTS  Stack every subject's time series into one matrix.
%
%   step05_concatenate_subjects(inputDirs, outputFile)
%
%   inputDirs   cell array of folders containing per-subject .mat files, each with
%               a variable 'ts_data' (time points x ROIs). Folders are processed in
%               the order given, files alphabetically within each folder.
%   outputFile  full path of the .mat file to write; it contains 'all_subjects'.
%
%   Example:
%     step05_concatenate_subjects({'data/dizygotic', 'data/monozygotic'}, ...
%                                 'data/all_subjects_twins.mat');

arguments
    inputDirs cell
    outputFile (1,:) char
end

all_subjects = [];

for f = 1:numel(inputDirs)
    current_folder = inputDirs{f};
    fprintf('Processing files in: %s\n', current_folder);

    files = dir(fullfile(current_folder, '*.mat'));

    for i = 1:numel(files)
        fname = files(i).name;

        % Skip hidden/system files (start with "._")
        if startsWith(fname, '._')
            continue
        end

        d = load(fullfile(current_folder, fname));

        if isfield(d, 'ts_data')
            all_subjects = [all_subjects; d.ts_data]; %#ok<AGROW>
        else
            warning('File %s does not contain the variable "ts_data". Skipping.', fname);
        end
    end
end

outDir = fileparts(outputFile);
if ~isempty(outDir) && ~exist(outDir, 'dir'), mkdir(outDir); end
save(outputFile, 'all_subjects');
fprintf('Data concatenation complete. Saved to %s.\n', outputFile);
end
