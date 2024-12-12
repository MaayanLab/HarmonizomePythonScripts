


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
gene_atb.edges = edgesinit(2261, [], 'GeneSym', [], [], [], [], [], 'GeneSym', [], [], [], [], true, []);


% read data
fid = fopen('input/ppi-weighted-highconfidence.txt', 'r');

for i = 1:1:gene_atb.edges.numedges
    
    currline = fgetl(fid);

    currcells = strsplitbyadr(currline, '\t');
    
    gene_atb.edges.source{i} = currcells{1};
    
    gene_atb.edges.target{i} = currcells{2};
    
    gene_atb.edges.weight(i) = str2double(currcells{3});
    
end

fclose(fid);


% map gene identifiers to NCBI Entrez Gene Symbols and Gene IDs
[gene_atb.edges.source, gene_atb.edges.sourcename, gene_atb.edges.sourceid, gene_atb.edges.sourceidname, discard] = genesymlookup(gene_atb.edges.source, [], 'gene', 'human', true, false, mappingfilespath);
gene_atb.edges = edgesdiscard(gene_atb.edges, discard);


% map attribute identifiers to NCBI Entrez Gene Symbols and Gene IDs
[gene_atb.edges.target, gene_atb.edges.targetname, gene_atb.edges.targetid, gene_atb.edges.targetidname, discard] = genesymlookup(gene_atb.edges.target, [], 'gene', 'human', true, false, mappingfilespath);
gene_atb.edges = edgesdiscard(gene_atb.edges, discard);


% remove duplicate edges
gene_atb.edges = edgesunique(gene_atb.edges);


% get edges from source to target and from target to source
gene_atb.edges = edges2symedges(gene_atb.edges);


% convert to matrix format (attribute table)
gene_atb.cm = edges2cm(gene_atb.edges);
gene_atb = rmfield(gene_atb, 'edges');


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


