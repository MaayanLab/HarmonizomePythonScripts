


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
% initialize edges structure
gene_atb.edges = edgesinit(1319440, [], 'GeneSym', [], 'Ensemble Acc', [], [], [], 'Tissue', [], 'Cell Type', [], [], true, []);

level = {'Not detected' 'Low' 'Medium' 'High'};

% read file
fid = fopen('input/proteinihc_normal_tissue.csv', 'r');

currline = fgetl(fid);

for i = 1:1:gene_atb.edges.numedges
    
    currline = fgetl(fid);
    
    currline(1) = [];
    currline(end) = [];
    
    currcells = strsplit(currline, '","', 'CollapseDelimiters', false);
    
    gene_atb.edges.source{i} = currcells{1};
    
    gene_atb.edges.sourcedesc{i} = currcells{1};
    
    gene_atb.edges.target{i} = lower(currcells{2});
    
    gene_atb.edges.targetdesc{i} = lower(currcells{3});
    
    gene_atb.edges.weight(i) = find(strcmp(level, currcells{4})) - 1;
    
end

fclose(fid);

save('input/gene_tissue_hpa_proteinihc_unmodified.mat', '-struct', 'gene_atb');
%}
gene_atb = load('input/gene_tissue_hpa_proteinihc_unmodified.mat', '-mat', 'edges');

gene_atb.edges.target(strcmp(gene_atb.edges.target, 'endometrium 1')) = {'endometrium'};
gene_atb.edges.target(strcmp(gene_atb.edges.target, 'endometrium 2')) = {'endometrium'};
gene_atb.edges.target(strcmp(gene_atb.edges.target, 'skin 1')) = {'skin'};
gene_atb.edges.target(strcmp(gene_atb.edges.target, 'skin 2')) = {'skin'};
gene_atb.edges.target(strcmp(gene_atb.edges.target, 'soft tissue 1')) = {'soft tissue'};
gene_atb.edges.target(strcmp(gene_atb.edges.target, 'soft tissue 2')) = {'soft tissue'};
gene_atb.edges.target(strcmp(gene_atb.edges.target, 'stomach 1')) = {'stomach'};
gene_atb.edges.target(strcmp(gene_atb.edges.target, 'stomach 2')) = {'stomach'};

e = sourcetargetcat(gene_atb.edges.target, gene_atb.edges.targetdesc);
ue = unique(e);
ut = unique(gene_atb.edges.target);

gene_atb.edges = rmfield(gene_atb.edges, {'targetdesc' 'targetdescname'});

gene_atb.edges = edgesunique_mean(gene_atb.edges); % this is super super slow
gene_atb.cm = edges2cm_nanfill(gene_atb.edges);

% gene_atb.cm = edges2cm_nanfill_mean(gene_atb.edges); % faster?

gene_atb = rmfield(gene_atb, 'edges');

fid = fopen('input/proteinihc_normal_tissue_ensemblegeneids.txt', 'w');

for i = 1:1:gene_atb.cm.numterms
    
    fprintf(fid, '%s\r\n', gene_atb.cm.term{i});
    
end

fclose(fid);

% get conversion from biomart

fid = fopen('input/proteinihc_normal_tissue_ensemblegeneids2entrezgeneids.txt', 'r');

currline = fgetl(fid);

ensgid = cell(100000, 1);
geneid = zeros([100000 1]);

i = 0;

while ~feof(fid)
    
    currline = fgetl(fid);
    
    currcells = strsplit(currline, '\t', 'CollapseDelimiters', false);
    
    if numel(currcells{1}) > 0 && numel(currcells{2}) > 0
        
        i = i + 1;
        
        ensgid{i} = currcells{1};
        
        geneid(i) = str2double(currcells{2});
        
    end
    
end

fclose(fid);

if i < 100000
    ensgid(i+1:end) = [];
    geneid(i+1:end) = [];
end

[ensgid, ui, ~] = unique(ensgid, 'stable');  % need to resolve this!
geneid = geneid(ui);

[o1, o2] = ismember(gene_atb.cm.termdesc, ensgid);
o2(o2 == 0) = [];
gene_atb.cm.termid = zeros([gene_atb.cm.numterms 1]);
gene_atb.cm.termidname = 'GeneID';
gene_atb.cm.termid(o1) = geneid(o2);
gene_atb.cm = cmrowdiscard(gene_atb.cm, ~o1);



% map gene symbols to entrez gene symbols and discard rows corresponding to
% un-mapped symbols (3 rows discarded, 15991 x 44 remaining)
[gene_atb.cm.term, gene_atb.cm.termname, gene_atb.cm.termid, gene_atb.cm.termidname, discard] = genesymlookup([], gene_atb.cm.termid, 'gene', 'human', false, false, mappingfilespath);
gene_atb.cm = cmrowdiscard(gene_atb.cm, discard);



if sum(isnan(gene_atb.cm.matrix(:))) > 0 % TRUE
    % remove rows and columns with more than 5% missing values (94
    % rows discarded, 0 columns discarded. 15897 x 44 remaining)
    gene_atb.cm = cmnantrim_frac(gene_atb.cm, 0.95, Inf, 0.95, Inf, 'column');
end



if sum(isnan(gene_atb.cm.matrix(:))) > 0 % TRUE
    % impute remaining missing values
    gene_atb.cm = cmnanimpute(gene_atb.cm, 'median', 'row');
end



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



% no need to quantile normalize or log2 this type of data, I think



if numel(unique(gene_atb.cm.term)) < gene_atb.cm.numterms % TRUE
    % merge measurements corresponding to the same gene (15802 unique gene
    % symbols out of 15897 total gene symbols, 15802 x 44 remaining)
    gene_atb.cm = cmrowmerge(gene_atb.cm, 'mean');
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



sm = cm2sm_cosine_nocluster(gene_atb.cm);
numidentical = sum(sm.matrix > (1 - 1e-12), 2) - 1;
suspectgenes = sm.term(numidentical > 0);
clear sm;
% if ~isempty(suspectgenes) % TRUE, but can't think of a good rationale for doing this on categorical data
%     % remove suspect rows that have identical values to at least one other row
%     discard = ismember(gene_atb.cm.term, suspectgenes);
%     gene_atb.cm = cmrowdiscard(gene_atb.cm, discard);
% end



sm = cm2sm_cosine_nocluster(cmtranspose(gene_atb.cm));
numidentical = sum(sm.matrix > (1 - 1e-12), 2) - 1;
suspectattributes = sm.term(numidentical > 0);
clear sm;
if ~isempty(suspectattributes) % FALSE
    % remove suspect columns that have identical values to at least one
    % other column
    discard = ismember(gene_atb.cm.entry, suspectattributes);
    gene_atb.cm = cmcoldiscard(gene_atb.cm, discard);
end



% view distributions of row and col stats
[~] = cmrowcolstats(gene_atb.cm, true, 3);



% % discard rows and cols with extremely low median or extremely low median
% % absolute deviation (mad) (discarded X row and Y col, # x # remaining)
% % DOES THIS MAKE ANY SENSE FOR THIS DATA?
% threshquantile = 0.001;
% gene_atb.cm = cmtrim_lowmedlowmad(gene_atb.cm, threshquantile);



% discard rows and cols with all zero values (not detected), 15788 x 44
% remaining
gene_atb.cm = cmtrim(gene_atb.cm, 1, Inf, 1, Inf);



% view distributions of row and col stats
[~] = cmrowcolstats(gene_atb.cm, true, 4);



% cluster and view matrix
[gene_atb.cm, ~] = cmcluster(gene_atb.cm, true);



% convert to sparse matrix
gene_atb.cm.matrix = sparse(gene_atb.cm.matrix);



% save result
save('output/gene_attribute_matrix_imported.mat', '-struct', 'gene_atb');



% NEED TO WORK THROUGH THIS PART

clearvars -except gene_atb;

gene_atb.cm.matrix = full(gene_atb.cm.matrix);

% normalize, prefer KS since it is non-parametric and therefore doesn't
% assume anything about the shape of the distribution of values per row
threshfrac = 0.05; % not relevant for normalization only
type = 'tertiary';
method = 'rows';
discardemptyvectors = false;
[~, gene_atb.cm] = cmthresh_ks(gene_atb.cm, threshfrac, type, method, discardemptyvectors);
% gene_atb.cm.matrix = robustzscore(gene_atb.cm.matrix, 2);
% gene_atb.cm.matrix = zscore(gene_atb.cm.matrix, 0, 2);
% daviesbouldin = cm2daviesbouldin(cmtranspose(gene_atb.cm));  % can't compute, no column labels



% search for best threshold (lowest davies-bouldin index)
% can't compute, no column labels
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

figure(6); clf; plot(threshfrac, daviesbouldin, '-ok'); ylabel('davies-bouldin index'); xlabel('threshold');
%}



% threshold.  note, we could equivalently just use cmthresh to get the
% thresholded matrix, but using cmthresh_ks here in case we ever want the
% final normalized matrix (the output that is currently "~")
threshfrac = 0.2;
type = 'tertiary';
method = 'matrix';
discardemptyvectors = true;
[gene_atb.cm, ~] = cmthresh_ks(gene_atb.cm, threshfrac, type, method, discardemptyvectors);  % Davies-Bouldin Index = 0.0848
% daviesbouldin = cm2daviesbouldin(cmtranspose(gene_atb.cm));



% cluster and view matrix
[gene_atb.cm, ~] = cmcluster(gene_atb.cm, true);



% view distributions of row and col stats
[~] = cmrowcolstats(gene_atb.cm, true, 7);



% convert to sparse matrix
gene_atb.cm.matrix = sparse(gene_atb.cm.matrix);



% save result
save('output/gene_attribute_matrix_prepared.mat', '-struct', 'gene_atb');
%}






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
if exist('output/gene_attribute_matrix_standardized.mat', 'file') == 0
    save('output/gene_attribute_matrix_standardized.mat', '-struct', 'gene_atb');
else
    error('file already exists?!');
end


