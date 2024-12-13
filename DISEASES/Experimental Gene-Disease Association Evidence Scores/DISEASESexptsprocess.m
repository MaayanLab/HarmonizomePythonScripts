


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




% initialize edges structure
gene_atb.edges = edgesinit(36650, [], 'GeneSym', [], 'Ensemble Acc', [], [], [], 'Disease', [], 'Disease Ontology ID', [], [], true, []);



% read file
fid = fopen('input/human_disease_experiments_full_20150202.tsv', 'r');

for i = 1:1:gene_atb.edges.numedges
    
    currline = fgetl(fid);
    
    currcells = strsplitbyadr(currline, '\t');
    
    gene_atb.edges.source{i} = currcells{2};
    
    gene_atb.edges.sourcedesc{i} = currcells{1};
    
    gene_atb.edges.target{i} = lower(currcells{4});
    
    gene_atb.edges.targetdesc{i} = currcells{3};
    
    e = find(currcells{6} == '=', 1, 'last');
    
    gene_atb.edges.weight(i) = -log10(str2double(currcells{6}(e+2:end)));  % -log10(p-value)
    
end

fclose(fid);

gene_atb.edges.weight(isinf(gene_atb.edges.weight)) = max(gene_atb.edges.weight(~isinf(gene_atb.edges.weight))) + 1;


% convert to connectivity matrix
gene_atb.cm = edges2cm(gene_atb.edges);
gene_atb = rmfield(gene_atb, 'edges');


% map gene symbols to entrez gene symbols and discard rows corresponding to
% un-mapped symbols
[gene_atb.cm.term, gene_atb.cm.termname, gene_atb.cm.termid, gene_atb.cm.termidname, discard] = genesymlookup(gene_atb.cm.term, [], 'gene', 'human', true, false, mappingfilespath);
gene_atb.cm = cmrowdiscard(gene_atb.cm, discard);


if numel(unique(gene_atb.cm.term)) < gene_atb.cm.numterms
    % merge measurements corresponding to the same gene
    gene_atb.cm = cmrowmerge(gene_atb.cm, 'max');
end


% map diseases to disease ontology and discard un-mapped diseases
atb_atb = load('output/gene_attribute_matrix_imported.mat', '-mat');
id_alt = load('output/gene_attribute_matrix_imported.mat', '-mat');

[o1, o2] = ismember(gene_atb.cm.entrydesc, id_alt.edges.targetdesc);
o2(o2 == 0) = [];

[o3, o4] = ismember(gene_atb.cm.entrydesc, atb_atb.cm.entrydesc);
o4(o4 == 0) = [];

gene_atb.cm.entry(o1) = id_alt.edges.source(o2);
gene_atb.cm.entrydesc(o1) = id_alt.edges.sourcedesc(o2);
gene_atb.cm.entry(o3) = atb_atb.cm.entry(o4);
gene_atb.cm.entrydesc(o3) = atb_atb.cm.entrydesc(o4);

discard = ~o1 & ~o3;
gene_atb.cm = cmcoldiscard(gene_atb.cm, discard);


if numel(unique(gene_atb.cm.entry)) < gene_atb.cm.numentries
    % merge measurements corresponding to the same disease
    gene_atb.cm = cmcolmerge(gene_atb.cm, 'max');
end


% remove empty rows and cols
gene_atb.cm = cmtrim(gene_atb.cm, 1, Inf, 1, Inf);


% make sure all appropriate connections to parent diseases are made
% gene_atb.em = cm2em_qryavg(gene_atb.cm, atb_atb.cm, false, []);
% discard = ismember(gene_atb.em.entry, gene_atb.cm.entry);
% gene_atb.em = cmcoldiscard(gene_atb.em, discard);
% if gene_atb.em.numentries > 0 && gene_atb.em.numterms == gene_atb.cm.numterms && all(strcmp(gene_atb.em.term, gene_atb.cm.term))
%     gene_atb.cm = cmtranspose(cmvertcat_prealigned(cmtranspose(gene_atb.cm), cmtranspose(gene_atb.em)));
%     gene_atb = rmfield(gene_atb, 'em');
% else
%     error('problem!');
% end
% gene_atb.cm = cmtrim(gene_atb.cm, 1, Inf, 1, Inf);
gene_atb.em = cm2em_qryavg(gene_atb.cm, atb_atb.cm, false, []);
gene_atb.cm = conmatmap(gene_atb.cm, gene_atb.em.term, gene_atb.em.entry);
gene_atb.em.matrix(gene_atb.cm.matrix~=0) = gene_atb.cm.matrix(gene_atb.cm.matrix~=0);
gene_atb.cm = gene_atb.em;
gene_atb = rmfield(gene_atb, 'em');
gene_atb.cm = cmtrim(gene_atb.cm, 1, Inf, 1, Inf);


% cluster and view matrix
[gene_atb.cm, ~] = cmcluster(gene_atb.cm, false);


% save result
if sum(sum(gene_atb.cm.matrix == 0)) > 1/3*numel(gene_atb.cm.matrix)
    gene_atb.cm.matrix = sparse(gene_atb.cm.matrix);
end
save('output/gene_attribute_matrix_cleaned.mat', '-struct', 'gene_atb');
if issparse(gene_atb.cm.matrix)
    gene_atb.cm.matrix = full(gene_atb.cm.matrix);
end
savepath = 'output/';
mkdir(savepath);
writecm([savepath 'gene_attribute_matrix_cleaned'], gene_atb.cm);
gzip([savepath 'gene_attribute_matrix_cleaned.txt']);
delete([savepath 'gene_attribute_matrix_cleaned.txt']);






% get standardized matrix
clearvars -except savepath gene_atb;

% threshfrac = 0.1;
% type = 'binary';
% method = 'matrix';
% discardemptyvectors = true;
% [~, gene_atb.cm] = cmthresh_ks(gene_atb.cm, threshfrac, type, method, discardemptyvectors);
gene_atb.cm = cm_standardize_ignorezeros(gene_atb.cm);

% cluster and view matrix
[gene_atb.cm, ~] = cmcluster(gene_atb.cm, false);

% save result
if sum(sum(gene_atb.cm.matrix == 0)) > 1/3*numel(gene_atb.cm.matrix)
    gene_atb.cm.matrix = sparse(gene_atb.cm.matrix);
end
save('output/gene_attribute_matrix_standardized.mat', '-struct', 'gene_atb');
if issparse(gene_atb.cm.matrix)
    gene_atb.cm.matrix = full(gene_atb.cm.matrix);
end
writecm([savepath 'gene_attribute_matrix_standardized'], gene_atb.cm);
gzip([savepath 'gene_attribute_matrix_standardized.txt']);
delete([savepath 'gene_attribute_matrix_standardized.txt']);


