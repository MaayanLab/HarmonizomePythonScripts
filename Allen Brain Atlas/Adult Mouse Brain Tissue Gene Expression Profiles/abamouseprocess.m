


clear all;

% get path to mapping files and add utilities to search path
currentpath = pwd;
parentpath = currentpath(1:find(currentpath=='\', 1, 'last')-1);
mappingfilespath = [parentpath '\mapping'];
utilitiespath = [parentpath '\utilities'];
searchpaths = strsplit(path, ';')';
if ~ismember(utilitiespath, searchpaths)
    addpath(utilitiespath, '-begin');
end
clear currentpath parentpath utilitiespath searchpaths;






% load connectivity matrix
gene_atb = load('input/gene_tissue_brainatlas.mat', '-mat', 'cm');



% map gene symbols to entrez gene symbols and discard rows corresponding to
% un-mapped symbols (4290 rows discarded)
[gene_atb.cm.term, gene_atb.cm.termname, gene_atb.cm.termid, gene_atb.cm.termidname, discard] = genesymlookup(removespecialchars(gene_atb.cm.term), [], 'gene', 'mouse', true, true, mappingfilespath);
gene_atb.cm = cmrowdiscard(gene_atb.cm, discard);



gene_atb.cm.matrix(gene_atb.cm.matrix == 0) = NaN;
if sum(isnan(gene_atb.cm.matrix(:))) > 0 % TRUE
    % remove rows and columns with at more than 5% missing values (1242
    % rows discarded, 271 columns discarded)
    gene_atb.cm = cmnantrim_frac(gene_atb.cm, 0.90, 1, 0.90, 1, 'column');
    gene_atb.cm = cmnantrim_frac(gene_atb.cm, 0.95, 1, 0.95, 1, 'row');
end



if sum(isnan(gene_atb.cm.matrix(:))) > 0 % TRUE
    % impute remaining missing values
    gene_atb.cm = cmnanimpute(gene_atb.cm, 'median', 'row');
end



% view column distributions
% figure(1);
% clf;
% subplot(2, 2, 1);
% hist(gene_atb.cm.matrix, max([10 gene_atb.cm.numterms/100]));
% subplot(2, 2, 2);
% hist(gene_atb.cm.matrix(:,randperm(gene_atb.cm.numentries,4)), max([10 gene_atb.cm.numterms/100]));
% subplot(2, 2, 3);
% hist(gene_atb.cm.matrix(:,randperm(gene_atb.cm.numentries,4)), max([10 gene_atb.cm.numterms/100]));
% subplot(2, 2, 4);
% hist(gene_atb.cm.matrix(:,randperm(gene_atb.cm.numentries,4)), max([10 gene_atb.cm.numterms/100]));



% quantile normalization not needed. data already normalized according to Supplemental Methods 2: Informatics Data Processing; nature05453-s02.doc; Lein et al. 2007. Genome-wide atlas of gene expression in the adult mouse brain. Nature. 445: 168-176.
if ~colsarenorm(gene_atb.cm.matrix) % TRUE
    % quantile normalize
    gene_atb.cm.matrix = quantilenormalization(gene_atb.cm.matrix);
end



% view column distributions
figure(2);
clf;
subplot(2, 2, 1);
hist(gene_atb.cm.matrix, max([10 gene_atb.cm.numterms/100]));
subplot(2, 2, 2);
hist(gene_atb.cm.matrix(:,randperm(gene_atb.cm.numentries,4)), max([10 gene_atb.cm.numterms/100]));
subplot(2, 2, 3);
hist(gene_atb.cm.matrix(:,randperm(gene_atb.cm.numentries,4)), max([10 gene_atb.cm.numterms/100]));
subplot(2, 2, 4);
hist(gene_atb.cm.matrix(:,randperm(gene_atb.cm.numentries,4)), max([10 gene_atb.cm.numterms/100]));



if (max(gene_atb.cm.matrix(:)) - min(gene_atb.cm.matrix(:))) > 30 % TRUE
    % log2 transformation
    gene_atb.cm.matrix = log2(gene_atb.cm.matrix);
end



if numel(unique(gene_atb.cm.term)) < gene_atb.cm.numterms % TRUE
    % merge measurements corresponding to the same gene (14263 unique gene
    % symbols out of 14385 total gene symbols)
    gene_atb.cm = cmrowmerge(gene_atb.cm, 'mean');
end



% view column distributions
figure(3);
clf;
subplot(2, 2, 1);
hist(gene_atb.cm.matrix, max([10 gene_atb.cm.numterms/100]));
subplot(2, 2, 2);
hist(gene_atb.cm.matrix(:,randperm(gene_atb.cm.numentries,4)), max([10 gene_atb.cm.numterms/100]));
subplot(2, 2, 3);
hist(gene_atb.cm.matrix(:,randperm(gene_atb.cm.numentries,4)), max([10 gene_atb.cm.numterms/100]));
subplot(2, 2, 4);
hist(gene_atb.cm.matrix(:,randperm(gene_atb.cm.numentries,4)), max([10 gene_atb.cm.numterms/100]));



sm = cm2sm_cosine_nocluster(gene_atb.cm);
numidentical = sum(sm.matrix > (1 - 1e-12), 2) - 1;
suspectgenes = sm.term(numidentical > 0);
clear sm;
if ~isempty(suspectgenes) % FALSE
    % remove suspect rows that have identical values to at least one other row
    discard = ismember(gene_atb.cm.term, suspectgenes);
    gene_atb.cm = cmrowdiscard(gene_atb.cm, discard);
end



sm = cm2sm_cosine_nocluster(cmtranspose(gene_atb.cm));
numidentical = sum(sm.matrix > (1 - 1e-12), 2) - 1;
suspectattributes = sm.term(numidentical > 0);
clear sm;
if ~isempty(suspectattributes) % FALSE
    % remove suspect columns that have identical values to at least one other column
    discard = ismember(gene_atb.cm.entry, suspectattributes);
    gene_atb.cm = cmcoldiscard(gene_atb.cm, discard);
end



% view distributions of row and col stats
[~] = cmrowcolstats(gene_atb.cm, true, 4);



% discard rows and cols with extremely low median or extremely low median
% absolute deviation (mad)
threshquantile = 0.001;
gene_atb.cm = cmtrim_lowmedlowmad(gene_atb.cm, threshquantile);



% view distributions of row and col stats
[~] = cmrowcolstats(gene_atb.cm, true, 5);



% cluster and view matrix
[gene_atb.cm, ~] = cmcluster(gene_atb.cm, true);



% save result
save('output/gene_attribute_matrix_imported.mat', '-struct', 'gene_atb');






clearvars -except gene_atb;

% normalize, prefer KS since it is non-parametric and therefore doesn't
% assume anything about the shape of the distribution of values per row
threshfrac = 0.05; % not relevant for normalization only
type = 'tertiary';
method = 'rows';
discardemptyvectors = false;
[~, gene_atb.cm] = cmthresh_ks(gene_atb.cm, threshfrac, type, method, discardemptyvectors);  % Davies-Bouldin Index = 
% gene_atb.cm.matrix = robustzscore(gene_atb.cm.matrix, 2);  % Davies-Bouldin Index = 
% gene_atb.cm.matrix = zscore(gene_atb.cm.matrix, 0, 2);  % Davies-Bouldin Index = 
% daviesbouldin = cm2daviesbouldin(cmtranspose(gene_atb.cm));  % can't calculate, don't have group labels for columns



% search for best threshold (lowest davies-bouldin index)
% threshfrac = [(0.01:0.01:0.04) (0.05:0.05:0.95) (0.96:0.01:0.99)]';
% type = 'tertiary';
% method = 'matrix';
% discardemptyvectors = false;
% 
% daviesbouldin = zeros([numel(threshfrac) 1]);
% 
% for i = 1:1:numel(threshfrac)
%     
%     tm = cmthresh(gene_atb.cm, threshfrac(i), type, method, discardemptyvectors);
% 
%     daviesbouldin(i) = cm2daviesbouldin(cmtranspose(tm));
%     
% end
% 
% figure(5); clf; plot(threshfrac, daviesbouldin, '-ok'); ylabel('davies-bouldin index'); xlabel('threshold');




% threshold.  note, we could equivalently just use cmthresh to get the
% thresholded matrix, but using cmthresh_ks here in case we ever want the
% final normalized matrix (the output that is currently "~")
threshfrac = 0.1;
type = 'tertiary';
method = 'matrix';
discardemptyvectors = true;
[gene_atb.cm, ~] = cmthresh_ks(gene_atb.cm, threshfrac, type, method, discardemptyvectors);  % Davies-Bouldin Index = 0.282
% daviesbouldin = cm2daviesbouldin(cmtranspose(gene_atb.cm));  % can't calculate, don't have group labels for columns



% cluster and view matrix
[gene_atb.cm, ~] = cmcluster(gene_atb.cm, true);



% view distributions of row and col stats
[~] = cmrowcolstats(gene_atb.cm, true, 6);



% convert to sparse matrix
gene_atb.cm.matrix = sparse(gene_atb.cm.matrix);



% save result
save('output/gene_attribute_matrix_prepared.mat', '-struct', 'gene_atb');






% get standardized matrix

gene_atb = load('output/gene_attribute_matrix_imported.mat', '-mat', 'cm');

threshfrac = 0.05;
type = 'tertiary';
method = 'rows';
discardemptyvectors = false;
[~, gene_atb.cm] = cmthresh_ks(gene_atb.cm, threshfrac, type, method, discardemptyvectors);

threshfrac = 0.1;
type = 'tertiary';
method = 'matrix';
discardemptyvectors = true;
[~, gene_atb.cm] = cmthresh_ks(gene_atb.cm, threshfrac, type, method, discardemptyvectors);

% cluster and view matrix
[gene_atb.cm, ~] = cmcluster(gene_atb.cm, false);

% save result
save('output/gene_attribute_matrix_standardized.mat', '-struct', 'gene_atb');


