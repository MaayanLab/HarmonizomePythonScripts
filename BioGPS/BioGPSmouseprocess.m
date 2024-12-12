


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




% initialize matrix structure
gene_atb.cm = cminit(36182, 153, [], 'GeneSym', [], 'ProbeID', [], 'GeneID', [], 'Tissue', [], [], [], [], []);


% read platform annotation
numprobes = 36118;
probeid = cell(numprobes, 1);
genesym = cell(numprobes, 1);
fid = fopen('input/gnf1m.annot2007.tsv', 'r');
currline = fgetl(fid);
for i = 1:1:numprobes
    currline = fgetl(fid);
    currcells = strsplit(currline, '\t', 'CollapseDelimiters', false);
    probeid{i} = currcells{1};
    genesym{i} = currcells{7};
end
fclose(fid);
discard = ismember(genesym, {'' 'obsoleted by Celera'});
probeid(discard) = [];
genesym(discard) = [];


% read data
fid = fopen('input/GNF1M_plus_macrophage_small.bioGPS.txt', 'r');
currline = fgetl(fid);
currcells = strsplit(currline, '\t', 'CollapseDelimiters', false);
entry = strtrim(currcells(2:end))';
for i = 1:1:gene_atb.cm.numentries
    if ~ismember([entry{i} '_1'], gene_atb.cm.entry)
        gene_atb.cm.entry{i} = [entry{i} '_1'];
    else
        gene_atb.cm.entry{i} = [entry{i} '_2'];
    end
end
for i = 1:1:gene_atb.cm.numterms
    currline = fgetl(fid);
    idx = regexp(currline, '\t');
    gene_atb.cm.termdesc{i} = currline(1:idx(1)-1);
end
fclose(fid);
gene_atb.cm.matrix = dlmread('input/GNF1M_plus_macrophage_small.bioGPS.txt', '\t', 1, 1);


% map probeids to genesyms
[o1, o2] = ismember(gene_atb.cm.termdesc, probeid);
o2(o2 == 0) = [];
gene_atb.cm.term(o1) = genesym(o2);
gene_atb.cm = cmrowdiscard(gene_atb.cm, ~o1);


% map gene identifiers to NCBI Entrez Gene Symbols and Gene IDs
[gene_atb.cm.term, gene_atb.cm.termname, gene_atb.cm.termid, gene_atb.cm.termidname, discard] = genesymlookup(gene_atb.cm.term, [], 'gene', 'both', true, true, mappingfilespath);
gene_atb.cm = cmrowdiscard(gene_atb.cm, discard);


% view column distributions
figure(1);
clf;
subplot(2, 2, 1);
hist(gene_atb.cm.matrix, max([10 gene_atb.cm.numterms/100]));
subplot(2, 2, 2);
hist(gene_atb.cm.matrix(:,randperm(gene_atb.cm.numentries,4)), max([10 gene_atb.cm.numterms/100]));
subplot(2, 2, 3);
hist(gene_atb.cm.matrix(:,randperm(gene_atb.cm.numentries,4)), max([10 gene_atb.cm.numterms/100]));
subplot(2, 2, 4);
hist(gene_atb.cm.matrix(:,randperm(gene_atb.cm.numentries,4)), max([10 gene_atb.cm.numterms/100]));


% remove rows and columns with more than 5% missing values
% gene_atb.cm.matrix(gene_atb.cm.matrix == 0) = NaN; % no zeros
% if sum(isnan(gene_atb.cm.matrix(:))) > 0 % FALSE
%     gene_atb.cm = cmnantrim_frac(gene_atb.cm, 0.6, Inf, 0.6, Inf, 'column');
%     gene_atb.cm = cmnantrim_frac(gene_atb.cm, 0.8, Inf, 0.8, Inf, 'column');
%     gene_atb.cm = cmnantrim_frac(gene_atb.cm, 0.9, Inf, 0.9, Inf, 'column');
%     gene_atb.cm = cmnantrim_frac(gene_atb.cm, 0.95, Inf, 0.95, Inf, 'column');
% end


% impute remaining missing values
% if sum(isnan(gene_atb.cm.matrix(:))) > 0 % FALSE
%     gene_atb.cm = cmnanimpute(gene_atb.cm, 'median', 'row');
% end


% log2 transform
gene_atb.cm.matrix = log2(gene_atb.cm.matrix);


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


% quantile normalize
gene_atb.cm.matrix = quantilenormalization(gene_atb.cm.matrix);


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


% merge rows corresponding to the same gene
if numel(unique(gene_atb.cm.term)) < gene_atb.cm.numterms
    gene_atb.cm = cmrowmerge(gene_atb.cm, 'mean');
end


% view column distributions
figure(4);
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
% if ~isempty(suspectgenes) % FALSE
%     % remove suspect rows that have identical values to at least one other row
%     discard = ismember(gene_atb.cm.term, suspectgenes);
%     gene_atb.cm = cmrowdiscard(gene_atb.cm, discard);
% end


sm = cm2sm_cosine_nocluster(cmtranspose(gene_atb.cm));
numidentical = sum(sm.matrix > (1 - 1e-12), 2) - 1;
suspectattributes = sm.term(numidentical > 0);
clear sm;
% if ~isempty(suspectattributes) % FALSE
%     % remove suspect columns that have identical values to at least one
%     % other column
%     discard = ismember(gene_atb.cm.entry, suspectattributes);
%     gene_atb.cm = cmcoldiscard(gene_atb.cm, discard);
% end


% average columns mapping to same attribute
gene_atb.cm.entry = cellfun(@(x) x(1:find(x == '_', 1, 'first')-1), gene_atb.cm.entry, 'UniformOutput', false);
gene_atb.cm.entryname = 'Tissue';

gene_atb.cm = cmcolmerge(gene_atb.cm, 'mean');


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


% view distributions of row and col stats
[~] = cmrowcolstats(gene_atb.cm, true, 4);


% cluster and view matrix
[gene_atb.cm, ~] = cmcluster(gene_atb.cm, false);
close force all;


% save result
save('output/gene_attribute_matrix_imported.mat', '-struct', 'gene_atb');






clearvars -except gene_atb;

% normalize, prefer KS since it is non-parametric and therefore doesn't
% assume anything about the shape of the distribution of values per row
threshfrac = 0.05; % not relevant for normalization only
type = 'tertiary';
method = 'rows';
discardemptyvectors = false;
[~, gene_atb.cm] = cmthresh_ks(gene_atb.cm, threshfrac, type, method, discardemptyvectors);
% gene_atb.cm.matrix = robustzscore(gene_atb.cm.matrix, 2);
% gene_atb.cm.matrix = zscore(gene_atb.cm.matrix, 0, 2);
% daviesbouldin = cm2daviesbouldin(cmtranspose(gene_atb.cm)); % can't measure, no group labels


% search for best threshold (lowest davies-bouldin index)
% can't measure, no group labels
%{
threshfrac = [(0.01:0.01:0.04) (0.05:0.05:0.95) (0.96:0.01:0.99)]';
type = 'tertiary';
method = 'matrix';
discardemptyvectors = false;

daviesbouldin = zeros([numel(threshfrac) 1]);

for i = 1:1:numel(threshfrac)
    
    tm = cmthresh(gene_atb.cm, threshfrac(i), type, method, discardemptyvectors);

    daviesbouldin(i) = cm2daviesbouldin(cmtranspose(tm));
    
end

figure(5); clf; plot(threshfrac, daviesbouldin, '-ok'); ylabel('davies-bouldin index'); xlabel('threshold');
%}


% threshold.  note, we could equivalently just use cmthresh to get the
% thresholded matrix, but using cmthresh_ks here in case we ever want the
% final normalized matrix
threshfrac = 0.15;
type = 'tertiary';
method = 'matrix';
discardemptyvectors = true;
[gene_atb.cm, nlpm] = cmthresh_ks(gene_atb.cm, threshfrac, type, method, discardemptyvectors);
% daviesbouldin = cm2daviesbouldin(cmtranspose(gene_atb.cm)); % can't measure, no group labels


% cluster and view matrix
[gene_atb.cm, ~] = cmcluster(gene_atb.cm, false);
close force all;


% view distributions of row and col stats
[~] = cmrowcolstats(gene_atb.cm, true, 6);


% convert to sparse matrix
gene_atb.cm.matrix = sparse(gene_atb.cm.matrix);


% save result
save('output/gene_attribute_matrix_prepared.mat', '-struct', 'gene_atb');






clearvars -except gene_atb nlpm;

gene_atb.cm = nlpm;
% daviesbouldin = cm2daviesbouldin(cmtranspose(gene_atb.cm)); % can't measure, no group labels
clear nlpm;


% cluster and view matrix
[gene_atb.cm, ~] = cmcluster(gene_atb.cm, false);
close force all;


% view distributions of row and col stats
[~] = cmrowcolstats(gene_atb.cm, true, 7);


% save result
save('output/gene_attribute_matrix_standardized.mat', '-struct', 'gene_atb');


