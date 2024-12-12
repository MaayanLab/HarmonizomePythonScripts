


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

%{

% initialize matrix structure
pos_atb.cm = cminit(117304, 26, [], 'Chr_Start_End', [], [], [], [], [], 'Tissue', [], 'TissueID', [], [], []);


% read column annotation
numcols = 26;
tissueid = cell(numcols, 1);
tissue = cell(numcols, 1);
fid = fopen('input/RRBS_DMRs_v2_column_annotation.txt', 'r');
for i = 1:1:numcols
    currline = fgetl(fid);
    currcells = strsplit(currline, '\t', 'CollapseDelimiters', false);
    tissueid{i} = currcells{1};
    tissue{i} = currcells{2};
end
fclose(fid);


% read data
fid = fopen('input/RRBS_DMRs_v2.tsv', 'r');
currline = fgetl(fid);
currcells = strsplit(currline, '\t', 'CollapseDelimiters', false);
pos_atb.cm.entrydesc = currcells(4:end)';
pos_atb.cm.entrydesc = cellfun(@(x) strrep(x, 'methylation_level_', ''), pos_atb.cm.entrydesc, 'UniformOutput', false);
for i = 1:1:pos_atb.cm.numterms
    currline = fgetl(fid);
    idx = regexp(currline, '\t');
    pos_atb.cm.term{i} = ['chr' currline(1:idx(1)-1) '_' currline(idx(1)+1:idx(2)-1) '_' currline(idx(2)+1:idx(3)-1)];
end
fclose(fid);
pos_atb.cm.matrix = dlmread('input/RRBS_DMRs_v2.tsv', '\t', 1, 3);


% map tissues to tissue ids
[o1, o2] = ismember(pos_atb.cm.entrydesc, tissueid);
o2(o2 == 0) = [];
pos_atb.cm.entry(o1) = tissue(o2);
pos_atb.cm = cmcoldiscard(pos_atb.cm, ~o1);


% view column distributions
figure(1);
clf;
subplot(2, 2, 1);
hist(pos_atb.cm.matrix, max([10 pos_atb.cm.numterms/100]));
subplot(2, 2, 2);
hist(pos_atb.cm.matrix(:,randperm(pos_atb.cm.numentries,4)), max([10 pos_atb.cm.numterms/100]));
subplot(2, 2, 3);
hist(pos_atb.cm.matrix(:,randperm(pos_atb.cm.numentries,4)), max([10 pos_atb.cm.numterms/100]));
subplot(2, 2, 4);
hist(pos_atb.cm.matrix(:,randperm(pos_atb.cm.numentries,4)), max([10 pos_atb.cm.numterms/100]));


% remove rows and columns with more than 5% missing values
if sum(isnan(pos_atb.cm.matrix(:))) > 0 % TRUE
    pos_atb.cm = cmnantrim_frac(pos_atb.cm, 0.6, Inf, 0.6, Inf, 'column');
    pos_atb.cm = cmnantrim_frac(pos_atb.cm, 0.8, Inf, 0.8, Inf, 'column');
    pos_atb.cm = cmnantrim_frac(pos_atb.cm, 0.9, Inf, 0.9, Inf, 'column');
    pos_atb.cm = cmnantrim_frac(pos_atb.cm, 0.95, Inf, 0.95, Inf, 'column');
end


% impute remaining missing values
if sum(isnan(pos_atb.cm.matrix(:))) > 0 % TRUE
    pos_atb.cm = cmnanimpute(pos_atb.cm, 'median', 'row');
end


% view column distributions
figure(2);
clf;
subplot(2, 2, 1);
hist(pos_atb.cm.matrix, max([10 pos_atb.cm.numterms/100]));
subplot(2, 2, 2);
hist(pos_atb.cm.matrix(:,randperm(pos_atb.cm.numentries,4)), max([10 pos_atb.cm.numterms/100]));
subplot(2, 2, 3);
hist(pos_atb.cm.matrix(:,randperm(pos_atb.cm.numentries,4)), max([10 pos_atb.cm.numterms/100]));
subplot(2, 2, 4);
hist(pos_atb.cm.matrix(:,randperm(pos_atb.cm.numentries,4)), max([10 pos_atb.cm.numterms/100]));


% quantile normalize
pos_atb.cm.matrix = quantilenormalization(pos_atb.cm.matrix);


% view column distributions
figure(3);
clf;
subplot(2, 2, 1);
hist(pos_atb.cm.matrix, max([10 pos_atb.cm.numterms/100]));
subplot(2, 2, 2);
hist(pos_atb.cm.matrix(:,randperm(pos_atb.cm.numentries,4)), max([10 pos_atb.cm.numterms/100]));
subplot(2, 2, 3);
hist(pos_atb.cm.matrix(:,randperm(pos_atb.cm.numentries,4)), max([10 pos_atb.cm.numterms/100]));
subplot(2, 2, 4);
hist(pos_atb.cm.matrix(:,randperm(pos_atb.cm.numentries,4)), max([10 pos_atb.cm.numterms/100]));


% merge rows corresponding to the same chromosome position
if numel(unique(pos_atb.cm.term)) < pos_atb.cm.numterms
    pos_atb.cm = cmrowmerge(pos_atb.cm, 'mean');
end


% view column distributions
figure(4);
clf;
subplot(2, 2, 1);
hist(pos_atb.cm.matrix, max([10 pos_atb.cm.numterms/100]));
subplot(2, 2, 2);
hist(pos_atb.cm.matrix(:,randperm(pos_atb.cm.numentries,4)), max([10 pos_atb.cm.numterms/100]));
subplot(2, 2, 3);
hist(pos_atb.cm.matrix(:,randperm(pos_atb.cm.numentries,4)), max([10 pos_atb.cm.numterms/100]));
subplot(2, 2, 4);
hist(pos_atb.cm.matrix(:,randperm(pos_atb.cm.numentries,4)), max([10 pos_atb.cm.numterms/100]));


% sm = cm2sm_cosine_nocluster(pos_atb.cm);
% numidentical = sum(sm.matrix > (1 - 1e-12), 2) - 1;
% suspectgenes = sm.term(numidentical > 0);
% clear sm;
% if ~isempty(suspectgenes) % CAN'T DO TOO MANY TERMS
%     % remove suspect rows that have identical values to at least one other row
%     discard = ismember(pos_atb.cm.term, suspectgenes);
%     pos_atb.cm = cmrowdiscard(pos_atb.cm, discard);
% end


sm = cm2sm_cosine_nocluster(cmtranspose(pos_atb.cm));
numidentical = sum(sm.matrix > (1 - 1e-12), 2) - 1;
suspectattributes = sm.term(numidentical > 0);
clear sm;
% if ~isempty(suspectattributes) % FALSE
%     % remove suspect columns that have identical values to at least one
%     % other column
%     discard = ismember(pos_atb.cm.entry, suspectattributes);
%     pos_atb.cm = cmcoldiscard(pos_atb.cm, discard);
% end


% view distributions of row and col stats
[~] = cmrowcolstats(pos_atb.cm, true, 5);


% cluster and view matrix
% [pos_atb.cm, ~] = cmcluster(pos_atb.cm, false);
% close force all;


% save result
save('output/gene_attribute_matrix_imported.mat', '-struct', 'pos_atb');
%}


pos_atb = load('output/gene_attribute_matrix_imported.mat', '-mat');


% map genomic loci to genes
loci = load('C:\Users\Andrew\Dropbox\GeneInfo\gene_loci_GRCh37_hg19.mat', '-mat', 'NumGenes', 'GeneSym', 'GeneID', 'ChrNum', 'StartPos', 'EndPos');
loci.Centroid = (loci.StartPos + loci.EndPos)/2;
loci.Len = loci.EndPos - loci.StartPos;
dmrs.Label = pos_atb.cm.term;
dmrs.NumLabels = pos_atb.cm.numterms;
dmrs.ChrNum = cell(dmrs.NumLabels, 1);
dmrs.StartPos = zeros(dmrs.NumLabels, 1);
dmrs.EndPos = zeros(dmrs.NumLabels, 1);
for i = 1:1:dmrs.NumLabels
    idx = find(dmrs.Label{i}=='_');
    dmrs.ChrNum{i} = dmrs.Label{i}(1:idx(1)-1);
    dmrs.StartPos(i) = str2double(dmrs.Label{i}(idx(1)+1:idx(2)-1));
    dmrs.EndPos(i) = str2double(dmrs.Label{i}(idx(2)+1:end));
end
dmrs.Centroid = (dmrs.StartPos + dmrs.EndPos)/2;
dmrs.Len = dmrs.EndPos - dmrs.StartPos;

loci_dmrs.lists = listsinit(loci.NumGenes, loci.GeneSym, 'GeneSym', [], [], loci.GeneID, 'GeneID', [], 'Chr_Start_End', [], [], [], [], true, []);

tic;
for i = 1:1:loci_dmrs.lists.numterms
    
    distance = sqrt((dmrs.Centroid - loci.Centroid(i)).^2)./strcmp(dmrs.ChrNum, loci.ChrNum{i});
    overlap = (dmrs.Len/2 + loci.Len(i)/2) - distance;
    overlapfraction = overlap/loci.Len(i);
    overlapfraction(overlapfraction > 1) = 1;
    overlapfraction(overlapfraction < 0) = 0;
    
    keep = overlapfraction > 0;
    entries = dmrs.Label(keep);
    weights = overlapfraction(keep);
    [loci_dmrs.lists.weights{i}, si] = sort(weights, 'descend');
    loci_dmrs.lists.entries{i} = entries(si);
    loci_dmrs.lists.numentries(i) = numel(loci_dmrs.lists.entries{i});
    
end
toc;

loci_dmrs.edges = lists2edges(loci_dmrs.lists);
gene_pos.edges = edgesunique_mean(loci_dmrs.edges);
gene_pos.lists = edges2lists(gene_pos.edges);


% map gene identifiers to NCBI Entrez Gene Symbols and Gene IDs
[gene_pos.lists.term, gene_pos.lists.termname, gene_pos.lists.termid, gene_pos.lists.termidname, discard] = genesymlookup(gene_pos.lists.term, [], 'gene', 'human', false, false, mappingfilespath);
gene_pos.lists = listsdiscard(gene_pos.lists, discard);


atb_pos.cm = cmtranspose(pos_atb.cm);

gene_atb.cm = cminit(gene_pos.lists.numterms, atb_pos.cm.numterms, gene_pos.lists.term, gene_pos.lists.termname, [], [], gene_pos.lists.termid, gene_pos.lists.termidname, atb_pos.cm.term, atb_pos.cm.termname, atb_pos.cm.termdesc, atb_pos.cm.termdescname, [], [], []);

tic;
for i = 1:1:gene_pos.lists.numterms
    
%     for j = 1:1:atb_pos.cm.numterms
%         
%         hit = ~isnan(atb_pos.cm.matrix(j,:));
%         pos = atb_pos.cm.entry(hit);
%         val = atb_pos.cm.matrix(j,hit)';
%         [o1, o2] = ismember(pos, gene_pos.lists.entries{i});
%         o2(o2==0) = [];
%         wgt = zeros(numel(pos), 1);
%         wgt(o1) = gene_pos.lists.weights{i}(o2);
%         
%         gene_atb.cm.matrix(i,j) = wgt'*val/sum(wgt);
%         
%     end
    
    [o1, o2] = ismember(atb_pos.cm.entry, gene_pos.lists.entries{i});
    o2(o2==0) = [];
    weights = zeros(1, atb_pos.cm.numentries);
    weights(o1) = gene_pos.lists.weights{i}(o2);
    weights = weights/sum(weights);
    
    gene_atb.cm.matrix(i,:) = weights*atb_pos.cm.matrix';
    
end
toc;

% view column distributions
figure(6);
clf;
subplot(2, 2, 1);
hist(gene_atb.cm.matrix, max([10 gene_atb.cm.numterms/100]));
subplot(2, 2, 2);
hist(gene_atb.cm.matrix(:,randperm(gene_atb.cm.numentries,4)), max([10 gene_atb.cm.numterms/100]));
subplot(2, 2, 3);
hist(gene_atb.cm.matrix(:,randperm(gene_atb.cm.numentries,4)), max([10 gene_atb.cm.numterms/100]));
subplot(2, 2, 4);
hist(gene_atb.cm.matrix(:,randperm(gene_atb.cm.numentries,4)), max([10 gene_atb.cm.numterms/100]));


% quantile normalize % YES OR NO?  WHAT ABOUT ON pos_atb?
gene_atb.cm.matrix = quantilenormalization(gene_atb.cm.matrix);


% view column distributions
figure(7);
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
[~] = cmrowcolstats(gene_atb.cm, true, 8);


% cluster and view matrix
[gene_atb.cm, ~] = cmcluster(gene_atb.cm, false);
% close force all;


% save result
save('output/gene_attribute_matrix_imported.mat', '-struct', 'gene_atb');






clearvars -except gene_atb;

% normalize, prefer KS since it is non-parametric and therefore doesn't
% assume anything about the shape of the distribution of values per row
threshfrac = 0.05; % not relevant for normalization only
type = 'tertiary';
method = 'rows';
support = [0-1e-6 1+1e-6];
discardemptyvectors = false;
[~, gene_atb.cm] = cmthresh_ks(gene_atb.cm, threshfrac, type, method, discardemptyvectors, support);
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
% support = [-1-1e-6 1+1e-6]; % [0 1] if type was binary in previous function call. [-1 1] if type was tertiary in previous function call.
support = 'unbounded';
discardemptyvectors = true;
[gene_atb.cm, nlpm] = cmthresh_ks(gene_atb.cm, threshfrac, type, method, discardemptyvectors, support);
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
% close force all;


% view distributions of row and col stats
[~] = cmrowcolstats(gene_atb.cm, true, 7);


% save result
save('output/gene_attribute_matrix_standardized.mat', '-struct', 'gene_atb');


