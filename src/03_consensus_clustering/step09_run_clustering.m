function step09_run_clustering(inputFile, outputFile)
% STEP09_RUN_CLUSTERING  Consensus clustering of the concatenated 64-ROI data.
%
%   step09_run_clustering(inputFile, outputFile)
%
%   inputFile   .mat file containing 'all_subjects' (nSubjects*nTimepoints x ROIs),
%               e.g. the Step 5 output after the Step 4 ROI selection
%   outputFile  .mat file to write; contains the 'Consensus' struct
%
% Requires the TbCAPs toolbox (CAP_ConsensusClustering_efficient) on the MATLAB path.
%
%   Example:
%     step09_run_clustering('data/all_subjects_twins64.mat', ...
%                           'results/ConsensusResults_downsampled10.mat');

arguments
    inputFile (1,:) char
    outputFile (1,:) char
end

%% --- Parameters ---
nSubjects = 240;
nTimepoints = 125;   % reduced from the full length to avoid out-of-memory errors during clustering
K_range = 2:20;      % number of clusters to test
replicates = double(20);

%% --- Load the data ---
% Load only 'all_subjects' to prevent a memory-doubling crash.
temp_data = load(inputFile, 'all_subjects');
big_matrix = temp_data.all_subjects;
clear temp_data
disp('Data loaded.');

%% --- Split into subjects ---
X_double = cell(1, nSubjects);   % explicit double-precision input
for s = 1:nSubjects
    idx_start = (s-1) * nTimepoints + 1;
    idx_end = s * nTimepoints;
    X_double{s} = double(big_matrix(idx_start:idx_end, :));
end
disp('Subjects loaded');

%% --- Run consensus clustering ---
disp('Running consensus clustering...');
Consensus = CAP_ConsensusClustering_efficient(X_double, K_range, 'items', 0.8, replicates, 'correlation');

%% --- Save results ---
outDir = fileparts(outputFile);
if ~isempty(outDir) && ~exist(outDir, 'dir'), mkdir(outDir); end
save(outputFile, 'Consensus');
disp('Consensus clustering completed and results saved.');
end
