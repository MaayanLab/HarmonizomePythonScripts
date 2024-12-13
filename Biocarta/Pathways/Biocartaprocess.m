


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
atb_gene.lists = listsinit(254, [], 'Pathway', [], 'Biocarta ID', [], [], [], 'GeneID', [], [], [], 'GeneID', false, []);

% read data
fid = fopen('input/dataset_20141217_original.txt', 'r');
for i = 1:1:atb_gene.lists.numterms
    currline = fgetl(fid);
    currcells = strsplitbyadr(currline, '\t');
    atb_gene.lists.term{i} = currcells{1};
    atb_gene.lists.termdesc{i} = currcells{2};
    atb_gene.lists.entries{i} = currcells(4:end)';
    atb_gene.lists.entryids{i} = str2double(currcells(4:end)');
    atb_gene.lists.numentries(i) = numel(atb_gene.lists.entries{i});
end
fclose(fid);
discard = atb_gene.lists.numentries == 0;
atb_gene.lists = listsdiscard(atb_gene.lists, discard);
clear fid ans currcells currline discard i;

% get new pathway names, ids and urls
nump = 314;
pname = cell(nump, 1);
purl = cell(nump, 1);
pid = cell(nump, 1);
fid = fopen('input/pathwaylist_20160404_parsed.txt', 'r');
currline = fgetl(fid);
for i = 1:1:nump
    currline = fgetl(fid);
    currcells = strsplitbyadr(currline, '\t');
    pname{i} = currcells{1};
    pid{i} = currcells{2};
    purl{i} = currcells{3};
end
fclose(fid);
clear fid currline i currcells ans;

% update pathways
[o1, o2] = ismember(removespecialchars(lower(atb_gene.lists.term)), removespecialchars(lower(pname)));
o2(o2 == 0) = [];
atb_gene.lists.term(o1) = pname(o2);
atb_gene.lists.termdesc(o1) = pid(o2);
atb_gene.lists = listsdiscard(atb_gene.lists, ~o1);
clear nump pname purl pid o1 o2;

% convert to matrix format
gene_atb.cm = cmtranspose(edges2cm(edgesunique(lists2edges(atb_gene.lists))));
clear atb_gene;

% map gene symbols to entrez gene symbols and discard rows corresponding to un-mapped symbols
% 74 unmapped out of 1478 (5%)
[gene_atb.cm.term, gene_atb.cm.termname, gene_atb.cm.termid, gene_atb.cm.termidname, discard] = genesymlookup([], gene_atb.cm.termid, 'gene', 'human', false, false, mappingfilespath);
gene_atb.cm = cmrowdiscard(gene_atb.cm, discard);
disp([num2str(sum(discard)) ' unmapped out of ' num2str(numel(discard)) ' (' num2str(round(100*sum(discard)/numel(discard))) '%)']);
clear discard;

% merge rows corresponding to the same gene
if numel(unique(gene_atb.cm.term)) < gene_atb.cm.numterms
    gene_atb.cm = cmrowmerge(gene_atb.cm, 'max');
end

% remove empty rows and cols
gene_atb.cm = cmtrim(gene_atb.cm, 1, Inf, 1, Inf);

% save cleaned data
gene_atb.cm.matrix = sparse(gene_atb.cm.matrix);
save('output/gene_attribute_matrix_cleaned.mat', '-struct', 'gene_atb');
gene_atb.cm.matrix = full(gene_atb.cm.matrix);

% standardize
tftype = 'raw'; % doesn't matter for binary data
gene_atb.cm = cm_tfidf_standardization(gene_atb.cm, tftype);
figure(1); clf; hist(gene_atb.cm.matrix(gene_atb.cm.matrix~=0), 100);
support = 'unbounded';
fignum = 2;
gene_atb.cm = cm_ksdensity_standardization_sparsematrix(gene_atb.cm, support, fignum);
clear tftype support fignum;

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
gene_atb.cm.matrix = (gene_atb.cm.matrix > 0) - (gene_atb.cm.matrix < 0);

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


