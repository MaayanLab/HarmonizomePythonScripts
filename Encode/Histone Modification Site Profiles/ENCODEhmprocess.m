


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

% initialize lists structure
atb_gene.lists = listsinit(2113, [], 'GeneSym_CellType_Species_Replicate', [], 'GeneSym', [], [], [], 'GeneSym', [], [], [], [], true, []);
condition = cell(atb_gene.lists.numterms, 1);
peaktype = cell(atb_gene.lists.numterms, 1);
halfmaxdistance = 2000;
sigma = sqrt(-halfmaxdistance^2/2/log(0.5));

% read data
fid = fopen('input/dataset_20160408_original.txt', 'r');
for i = 1:1:atb_gene.lists.numterms
    currline = fgetl(fid);
    currcells = strsplitbyadr(currline, '\t');
    atb_gene.lists.termdesc{i} = currcells{1};
    atb_gene.lists.term{i} = [currcells{1} '_' currcells{2} '_' currcells{3} '_' currcells{4}];
    condition{i} = [currcells{1} '_' currcells{2} '_' currcells{3}];
    peaktype{i} = currcells{5};
    atb_gene.lists.entries{i} = currcells(6:2:end-1)';
    distance = str2double(currcells(7:2:end)');
    atb_gene.lists.weights{i} = exp(-distance.^2/2/sigma.^2);
    atb_gene.lists.numentries(i) = numel(atb_gene.lists.entries{i});
end
fclose(fid);
discard = atb_gene.lists.numentries == 0 | cellfun(@isempty, atb_gene.lists.term) | ismember(atb_gene.lists.term, {'-666' '' '-'}) | cellfun(@isempty, atb_gene.lists.termdesc) | ismember(atb_gene.lists.termdesc, {'-666' '' '-'});
atb_gene.lists = listsdiscard(atb_gene.lists, discard);
condition(discard) = [];
peaktype(discard) = [];
clear fid ans currcells currline discard i distance;

% save unmodified data
save('output/attribute_gene_lists_unmodified.mat', '-struct', 'atb_gene');
save('output/extra_vars_unmodified.mat', '-mat', 'peaktype', 'halfmaxdistance', 'sigma', 'condition');
clear halfmaxdistance sigma peaktype condition;

% convert to matrix format
gene_atb.cm = cmtranspose(lists2cm(atb_gene.lists));
clear atb_gene;

% map gene symbols to entrez gene symbols and discard rows corresponding to un-mapped symbols
% 5708 unmapped out of 30376 (19%), many non-coding genes
[gene_atb.cm.term, gene_atb.cm.termname, gene_atb.cm.termid, gene_atb.cm.termidname, discard] = genesymlookup(gene_atb.cm.term, [], 'gene', 'both', true, true, mappingfilespath);
gene_atb.cm = cmrowdiscard(gene_atb.cm, discard);
disp([num2str(sum(discard)) ' unmapped out of ' num2str(numel(discard)) ' (' num2str(round(100*sum(discard)/numel(discard))) '%)']);
clear discard;

% merge rows corresponding to the same gene
if numel(unique(gene_atb.cm.term)) < gene_atb.cm.numterms
    gene_atb.cm = cmrowmerge(gene_atb.cm, 'max');
end

% merge cols corresponding to replicates of the same condition
gene_atb.cm.entry = cellfun(@(x) x(1:find(x=='_', 1, 'last')-1), gene_atb.cm.entry, 'uniformoutput', false);
gene_atb.cm.entryname = 'GeneSym_CellType_Species';
if numel(unique(gene_atb.cm.entry)) < gene_atb.cm.numentries
    gene_atb.cm = cmcolmerge(gene_atb.cm, 'mean'); % note using mean here because we want to average replicates
end

% remove empty rows and cols
gene_atb.cm = cmtrim(gene_atb.cm, 1, Inf, 1, Inf);

% save cleaned data
gene_atb.cm.matrix = sparse(gene_atb.cm.matrix);
save('output/gene_attribute_matrix_cleaned.mat', '-struct', 'gene_atb');
gene_atb.cm.matrix = full(gene_atb.cm.matrix);

% standardize
gene_atb.cm.matrix(gene_atb.cm.matrix==0) = NaN;
figure(3); clf; hist(nanmax(gene_atb.cm.matrix, [], 2), 100); xlabel('rowmax');
figure(4); clf; hist(nanmax(gene_atb.cm.matrix, [], 1), 100); xlabel('colmax');
gene_atb.cm.matrix(isnan(gene_atb.cm.matrix)) = 0;
tftype = 'raw'; % values in entire matrix are already on the same scale, already normalized
gene_atb.cm = cm_tfidf_standardization(gene_atb.cm, tftype);
figure(1); clf; hist(gene_atb.cm.matrix(gene_atb.cm.matrix~=0), 100);
support = 'unbounded';
fignum = 2;
gene_atb.cm = cm_ksdensity_standardization_sparsematrix(gene_atb.cm, support, fignum);
clear tftype support fignum;

% separate transcription factors from histone modifications
histonemodifications = {'H3K27ac';'H3K27me3';'H3K36me3';'H3K4me1';'H3K4me2';'H3K4me3';'H3K79me2';'H3K79me3';'H3K9ac';'H3K9me1';'H3K9me3';'H3ac';'H4K20me1'};
ishm = ismember(gene_atb.cm.entrydesc, histonemodifications);
gene_atb.cm = cmtrim(cmcoldiscard(gene_atb.cm, ~ishm), 1, Inf, 1, Inf);
clear histonemodifications ishm;

% cluster and view standardized matrix
gene_atb.cm = cmcluster_nogram(gene_atb.cm, true);

% save standardized data
gene_atb.cm.matrix = sparse(gene_atb.cm.matrix);
save('output/gene_attribute_matrix_standardized.mat', '-struct', 'gene_atb');
gene_atb.cm.matrix = full(gene_atb.cm.matrix);
writecm('output/gene_attribute_matrix_standardized', gene_atb.cm);
gzip('output/gene_attribute_matrix_standardized.txt');
delete('output/gene_attribute_matrix_standardized.txt');
genes = gene_atb.cm.term;
atbs = gene_atb.cm.entry;

% threshold
gene_atb.cm.matrix = (gene_atb.cm.matrix > 0.5) - (gene_atb.cm.matrix < 0); % keep only edges with above median signal

% save thresholded data
gene_atb.cm.matrix = sparse(gene_atb.cm.matrix);
save('output/gene_attribute_matrix_thresholded.mat', '-struct', 'gene_atb');
gene_atb.cm.matrix = full(gene_atb.cm.matrix);
writecm('output/gene_attribute_matrix_thresholded', gene_atb.cm);
gzip('output/gene_attribute_matrix_thresholded.txt');
delete('output/gene_attribute_matrix_thresholded.txt');
clear gene_atb;

% align cleaned matrix with standardized matrix and save
load('output/gene_attribute_matrix_cleaned.mat', '-mat', 'cm');
cm = conmatmap(cm, genes, atbs);
save('output/gene_attribute_matrix_cleaned.mat', '-mat', 'cm');
writecm('output/gene_attribute_matrix_cleaned', cm);
gzip('output/gene_attribute_matrix_cleaned.txt');
delete('output/gene_attribute_matrix_cleaned.txt');
clear cm genes atbs mappingfilespath;


