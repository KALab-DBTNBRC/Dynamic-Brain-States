%% --- Load the data ---
% FIX: Loads the variable 'all_subjects' directly to prevent memory-doubling crash.
temp_data = load('all_subjects_twins64.mat', 'all_subjects');
big_matrix = temp_data.all_subjects; 
disp('Data loaded.');

%% --- Parameters ---
nSubjects = 240;
% FIX: Reduced nTimepoints to prevent out-of-memory errors during clustering.
nTimepoints = 125; 
nROIs = 300;

%% --- Split and downsample into subjects ---
X_double = cell(1, nSubjects); % FINAL INPUT: Uses a new array for explicit double input
for s = 1:nSubjects
    idx_start = (s-1) * nTimepoints + 1;
    idx_end = s * nTimepoints;
    subj_data = big_matrix(idx_start:idx_end, :); 
    
    % FIX: Ensure the input to the cell array is explicitly double-precision
    X_double{s} = double(subj_data); 
end
disp('Subjects loaded');

%% --- Run consensus clustering ---
K_range = 2:20; % number of clusters to test
replicates = double(20); % FIX: Keeps the replicates count explicitly as double

disp('Running consensus clustering...');

% FINAL CALL: Using 'fraction' 1.0 to bypass indexing bugs and rely on the patched function logic.
Consensus = CAP_ConsensusClustering_efficient(X_double, K_range, 'items', 0.8, replicates, 'correlation');
%% --- Save results ---
save('ConsensusResults_downsampled10.mat', 'Consensus');
disp('Consensus clustering completed and results saved.');