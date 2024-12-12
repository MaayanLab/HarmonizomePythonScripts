


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
gene_atb.cm = cminit(26374, 668, [], 'GeneID', [], [], [], 'GeneID', [], 'CellLine', [], 'Tissue', [], [], []);


% read cell line tissue info
cellline = cell(676, 1);
tissue = cell(676, 1);
fid = fopen('input/cellline_tissue_list.txt', 'r');
for i = 1:1:676
    currline = fgetl(fid);
    currcells = strsplit(currline, '\t');
    cellline{i} = currcells{1};
    tissue{i} = currcells{2};
end
fclose(fid);


% read data
fid = fopen('input/Data4_genecnvs/140331_cellLineCopyNumber_NA2NAN.txt', 'r');

currline = fgetl(fid);

currcells = strsplit(currline, '\t', 'CollapseDelimiters', false);

gene_atb.cm.entry = upper(currcells(2:end))';

[o1, o2] = ismember(gene_atb.cm.entry, cellline);
o2(o2 == 0) = [];
gene_atb.cm.entrydesc(o1) = tissue(o2);

fclose(fid);

gene_atb.cm.matrix = dlmread('input/Data4_genecnvs/140331_cellLineCopyNumber_NA2NAN.txt', '\t', 1, 0);

gene_atb.cm.termid = gene_atb.cm.matrix(:,1);

gene_atb.cm.matrix(:,1) = [];

gene_atb.cm.term = cellfun(@num2str, num2cell(gene_atb.cm.termid), 'UniformOutput', false);


% map gene identifiers to NCBI Entrez Gene Symbols and Gene IDs
[gene_atb.cm.term, gene_atb.cm.termname, gene_atb.cm.termid, gene_atb.cm.termidname, discard] = genesymlookup([], gene_atb.cm.termid, 'gene', 'human', false, false, mappingfilespath);
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
if sum(isnan(gene_atb.cm.matrix(:))) > 0
    gene_atb.cm = cmnantrim_frac(gene_atb.cm, 0.6, Inf, 0.6, Inf, 'column');
    gene_atb.cm = cmnantrim_frac(gene_atb.cm, 0.8, Inf, 0.8, Inf, 'column');
    gene_atb.cm = cmnantrim_frac(gene_atb.cm, 0.9, Inf, 0.9, Inf, 'column');
    gene_atb.cm = cmnantrim_frac(gene_atb.cm, 0.95, Inf, 0.95, Inf, 'column');
end


% impute remaining missing values
if sum(isnan(gene_atb.cm.matrix(:))) > 0
    gene_atb.cm = cmnanimpute(gene_atb.cm, 'median', 'row');
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


% view distributions of row and col stats
[~] = cmrowcolstats(gene_atb.cm, true, 5);


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





%{

% map cell lines and tissues
ontologies = {'mouse_adult_gross_anatomy_ontology_nodes' 'uber_anatomy_ontology_nodes' 'experimental_factor_ontology_nodes' 'cell_ontology_nodes' 'cell_line_ontology_nodes' 'brenda_tissue_enzyme_source_ontology_nodes'}';

folder = {'ImportedData' 'PreparedData' 'StandardizedData'}';
attribute = 'cellline';
filesuffix = 'klijn_cnv_20150415';

for fi = 1:1:numel(folder)
    
    gene_atb = load(['../' folder{fi} '/gene_' attribute '_' filesuffix '.mat'], '-mat', 'cm');

    if strcmpi(gene_atb.cm.entryname, 'cellline') || strcmpi(gene_atb.cm.entryname, 'cell line') || strcmpi(gene_atb.cm.entryname, 'tissue')

        [gene_atb.cm.entry, gene_atb.cm.entryname, isunmapped] = celltissuelabelmapper(gene_atb.cm.entry, gene_atb.cm.entryname, ontologies);

        writeunmappedlabels(gene_atb.cm.entry, gene_atb.cm.entryname, isunmapped, ['../UnmappedAttributes/' folder{fi} '_entry_LABELNAME_' filesuffix]);

    end

    if strcmpi(gene_atb.cm.entrydescname, 'cellline') || strcmpi(gene_atb.cm.entrydescname, 'cell line') || strcmpi(gene_atb.cm.entrydescname, 'tissue')

        [gene_atb.cm.entrydesc, gene_atb.cm.entrydescname, isunmapped] = celltissuelabelmapper(gene_atb.cm.entrydesc, gene_atb.cm.entrydescname, ontologies);
        
        writeunmappedlabels(gene_atb.cm.entrydesc, gene_atb.cm.entrydescname, isunmapped, ['../UnmappedAttributes/' folder{fi} '_entrydesc_LABELNAME_' filesuffix]);

    end
    
    save(['../' folder{fi} '/gene_' attribute '_' filesuffix '.mat'], '-struct', 'gene_atb');
    
end
%}

