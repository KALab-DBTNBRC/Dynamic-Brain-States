function step06_consensus_metrics(prototypesFile, outputDir, doNull)
% STEP06_CONSENSUS_METRICS  Corrected consensus clustering for choosing K
% (Supplementary Fig. S2A-D, Table Sx).
%
%   step06_consensus_metrics(prototypesFile, outputDir)
%   step06_consensus_metrics(prototypesFile, outputDir, true)   % ROI-permuted null
%
%   prototypesFile  .mat file containing X_prototypes (1 x 240 cell, each 200 x 64)
%   outputDir       folder for the .mat, .csv and .png outputs (created if missing)
%   doNull          false (default) = real data. true = null reference in which each
%                   ROI is shuffled independently across prototypes.
%
% WHAT WAS WRONG BEFORE
%   CAP_ConsensusClustering(_efficient) (TbCAPs) expects each subject's data as
%   ROIs x frames. X_prototypes{s} is 200 prototypes x 64 ROIs, so the toolbox
%   clustered the 64 ROI columns of each subject (240 x 64 = 15,360 items)
%   instead of the 48,000 prototypes. This script clusters the prototypes.
%
% WHY NOT JUST TRANSPOSE AND CALL THE TOOLBOX
%   A 48,000 x 48,000 consensus matrix is 9.2 GB (single) per K, i.e. ~175 GB
%   for K = 2..20. This script never stores it: it keeps the cluster labels of
%   every resample (48,000 x 20 int8 per K) and computes every metric EXACTLY
%   over all 1.15e9 prototype pairs, block by block.
%
% PROCEDURE (same as the toolbox call in step09_run_clustering.m)
%   items = 48,000 prototypes, features = 64 ROIs, K = 2..20,
%   20 resamples per K, 80% of items drawn without replacement per resample,
%   kmeans 'correlation' (Start 'plus', 5 replicates per resample). 5 replicates
%   reproduced the saved workspace run (SD of consensus within 0.7-8% at every K
%   tested) far better than 1 replicate (up to 26% off at K = 2-4).
%
% METRICS (consensus index of a pair = times clustered together / times co-drawn)
%   PAC      = fraction of pairs with 0.1 < consensus < 0.9  (Senbabaoglu et al. 2014)
%   CDF area = ConsensusClusterPlus definition: 100-bin histogram of the consensus
%              index, area = sum over bins of cumulative fraction x bin width
%   DeltaK   = A(2) for K = 2; (A(K) - A(K-1)) / A(K-1) for K > 2 (Monti et al. 2003,
%              as implemented in ConsensusClusterPlus)
%
% RUNTIME / MEMORY: roughly 1-2 h and ~2 GB RAM on a lab workstation.
% Requires: Statistics and Machine Learning Toolbox. MATLAB R2019b+ (arguments blocks).

arguments
    prototypesFile (1,:) char
    outputDir (1,:) char
    doNull (1,1) logical = false
end

if ~exist(outputDir, 'dir'), mkdir(outputDir); end
rng(4 + doNull, 'twister');                % reproducible

%% ---- Load prototypes --------------------------------------------------
S = load(prototypesFile, 'X_prototypes');
X_prototypes = S.X_prototypes; clear S;

X = cell2mat(X_prototypes(:));             % (240*200) x 64 : items = prototypes
assert(isequal(size(X), [48000 64]), 'Expected 48,000 prototypes x 64 ROIs, got %d x %d', size(X,1), size(X,2));
if doNull
    for j = 1:size(X,2), X(:,j) = X(randperm(size(X,1)), j); end
    tag = 'null';
else
    tag = 'real';
end

%% ---- Parameters -------------------------------------------------------
K_range = 2:20;
N       = 20;                              % resamples per K
Pcc     = 0.8;                             % fraction of items per resample
n       = size(X,1);
m       = floor(Pcc*n);
blk     = 1000;                            % row block size for pair counting
nSub    = 2000;                            % items shown in consensus-matrix figures
subIdx  = sort(randperm(n, nSub));
K_show  = 2:6;                             % store subset consensus matrices for these K

nK       = numel(K_range);
L_all    = zeros(n, N, nK, 'int8');        % cluster label per item per resample (0 = not drawn)
T_all    = zeros(N+1, N+1, nK);            % T(a+1,b+1): #pairs co-clustered a times, co-drawn b times
PAC      = zeros(1,nK);  CDF_area = zeros(1,nK);  DeltaK = zeros(1,nK);
meanC    = zeros(1,nK);  sdC      = zeros(1,nK);
F_edges  = zeros(100,nK);                  % CDF at 0.01,0.02,...,1.00
cdf_x    = cell(1,nK);   cdf_F    = cell(1,nK);
subC     = zeros(nSub, nSub, numel(K_show), 'single');
subOrder = zeros(nSub, numel(K_show));

%% ---- Main loop --------------------------------------------------------
for kk = 1:nK
    K = K_range(kk);  t0 = tic;
    L = zeros(n, N, 'int8');
    for h = 1:N
        idx = randperm(n, m);
        lab = kmeans(X(idx,:), K, 'Distance','correlation', 'Start','plus', ...
                     'Replicates',5, 'MaxIter',100, 'EmptyAction','singleton');
        L(idx,h) = int8(lab);
    end
    L_all(:,:,kk) = L;

    % ---- exact pair-count table over all pairs i<j, block by block ----
    T = zeros(N+1, N+1);
    D = single(L > 0);
    for i0 = 1:blk:n
        rows = i0:min(i0+blk-1, n);
        cols = i0:n;
        B = D(rows,:) * D(cols,:)';                        % times co-drawn
        A = zeros(numel(rows), numel(cols), 'single');     % times co-clustered
        for r = 1:N
            li = L(rows,r);  lj = L(cols,r)';
            A = A + single((li == lj) & (li > 0) & (lj > 0));
        end
        upper = cols > rows(:);                            % keep j > i only
        lin = sub2ind([N+1 N+1], double(A(upper))+1, double(B(upper))+1);
        T = T + reshape(accumarray(lin, 1, [(N+1)^2 1]), N+1, N+1);
    end
    T_all(:,:,kk) = T;

    % ---- metrics from the table ----
    [aa, bb] = ndgrid(0:N, 0:N);
    ok  = (bb > 0) & (T > 0);
    cnt = T(ok); a = aa(ok); b = bb(ok); v = a ./ b; tot = sum(cnt);
    meanC(kk) = sum(cnt .* v) / tot;
    sdC(kk)   = sqrt(sum(cnt .* (v - meanC(kk)).^2) / tot);
    PAC(kk)   = sum(cnt(10*a > b & 10*a < 9*b)) / tot;          % 0.1 < c < 0.9 (exact)
    for e = 1:100
        F_edges(e,kk) = sum(cnt(100*a <= e*b)) / tot;            % CDF at e/100
    end
    CDF_area(kk) = sum(F_edges(:,kk)) / 100;
    [ux, ~, ic] = unique(v);
    cdf_x{kk} = ux;  cdf_F{kk} = cumsum(accumarray(ic, cnt)) / tot;

    % ---- consensus matrix on a random subset (for Fig. S2D) ----
    s = find(K_show == K, 1);
    if ~isempty(s) && ~doNull
        Ls = double(L(subIdx,:));  Ac = zeros(nSub);  Bc = zeros(nSub);
        for r = 1:N
            dr = Ls(:,r) > 0;
            Bc = Bc + double(dr & dr');
            Ac = Ac + double((Ls(:,r) == Ls(:,r)') & dr & dr');
        end
        C = Ac ./ Bc;  C(1:nSub+1:end) = 1;  C(isnan(C)) = 0;
        Dm = 1 - C;  Dm(1:nSub+1:end) = 0;  Dm = (Dm + Dm')/2;
        Z = linkage(squareform(Dm, 'tovector'), 'average');
        f = figure('Visible','off'); [~,~,ord] = dendrogram(Z, 0); close(f);
        subC(:,:,s) = single(C);  subOrder(:,s) = ord(:);
    end
    fprintf('K=%2d  PAC=%.4f  CDF area=%.4f  mean=%.4f  sd=%.4f  (%.0f s)\n', ...
            K, PAC(kk), CDF_area(kk), meanC(kk), sdC(kk), toc(t0));
end
DeltaK(1) = CDF_area(1);
DeltaK(2:end) = (CDF_area(2:end) - CDF_area(1:end-1)) ./ CDF_area(1:end-1);

%% ---- Save compact results (well under 2 GB) ----------------------------
outFile = fullfile(outputDir, sprintf('consensus_metrics_%s.mat', tag));
save(outFile, 'K_range','PAC','CDF_area','DeltaK','meanC','sdC','F_edges', ...
     'cdf_x','cdf_F','T_all','L_all','subIdx','subC','subOrder','K_show','N','Pcc','-v7');
fprintf('Saved %s\n', outFile);

%% ---- Table Sx ------------------------------------------------------------
TableSx = table(K_range(:), PAC(:), CDF_area(:), DeltaK(:), meanC(:), sdC(:), ...
    'VariableNames', {'K','PAC','CDF_area','DeltaK','MeanConsensus','SDConsensus'});
writetable(TableSx, fullfile(outputDir, sprintf('TableSx_consensus_metrics_%s.csv', tag)));
disp(TableSx);
if doNull, return; end   % figures are made from the real-data run
% Publication figures with the null overlay are made from consensus_metrics_real.mat
% and consensus_metrics_null.mat (the plotting script is not part of this repository).

%% ---- Quick-look figures (S2A-S2D) ---------------------------------------
figure('Color','w','Position',[100 100 1400 900]);
subplot(2,2,1); plot(K_range, PAC, 'k-o','LineWidth',1.5,'MarkerFaceColor','k');
xlabel('K'); ylabel('PAC (0.1-0.9)'); title('S2A  PAC'); grid on; xlim([2 20]);
subplot(2,2,2); hold on; cm = parula(nK);
for kk = 1:nK, stairs([0; cdf_x{kk}], [0; cdf_F{kk}], 'Color', cm(kk,:), 'LineWidth', 1.2); end
xlabel('Consensus index'); ylabel('CDF'); title('S2B  CDF, K = 2-20'); colormap(parula);
cb = colorbar; cb.Ticks = linspace(0,1,5); cb.TickLabels = round(linspace(2,20,5)); cb.Label.String = 'K';
subplot(2,2,3); plot(K_range, DeltaK, 'k-o','LineWidth',1.5,'MarkerFaceColor','k');
xlabel('K'); ylabel('\DeltaK (relative change in CDF area)'); title('S2C  \DeltaK'); grid on; xlim([2 20]);
saveas(gcf, fullfile(outputDir, 'FigS2_ABC_quicklook.png'));
figure('Color','w','Position',[100 100 1500 450]);
showK = [3 4 5];
for p = 1:3
    s = find(K_show == showK(p)); o = subOrder(:,s);
    subplot(1,3,p); imagesc(subC(o,o,s), [0 1]); axis square off;
    title(sprintf('S2D  K = %d', showK(p)));
end
colormap(flipud(gray)); colorbar;
saveas(gcf, fullfile(outputDir, 'FigS2D_quicklook.png'));
end
