


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
atb_gene.lists = listsinit(269, [], 'GeneSym_PMID_Cell or Tissue_Perturbation_Organism_GEO Platform_GEO Accession_Direction', [], 'GeneSym', [], [], [], 'GeneSym', [], [], [], [], true, []);


% read data
fid = fopen('input/TF-LOF_Expression_from_GEO.txt', 'r');

for i = 1:1:atb_gene.lists.numterms
    
    currline = fgetl(fid);

    currcells = strsplitbyadr(currline, '\t');
    
    atb_gene.lists.term{i} = currcells{1};
    
    ui = regexp(currcells{1}, '_');
    atb_gene.lists.termdesc{i} = currcells{1}(1:ui(1)-1);
    
    gv = cellfun(@(x) {x(1:regexp(x, ',')-1) str2double(x(regexp(x, ',')+1:end))}, currcells(3:end-1)', 'UniformOutput', false);
    atb_gene.lists.entries{i} = cellfun(@(x) x{1}, gv, 'UniformOutput', false);
    atb_gene.lists.weights{i} = cellfun(@(x) x{2}, gv);

    atb_gene.lists.numentries(i) = numel(atb_gene.lists.entries{i});
    
    if strcmp(currcells{1}(ui(end)+1:end), 'up')
        atb_gene.lists.weights{i} = abs(atb_gene.lists.weights{i});
    elseif strcmp(currcells{1}(ui(end)+1:end), 'down')
        atb_gene.lists.weights{i} = -1*abs(atb_gene.lists.weights{i});
    else
        atb_gene.lists.weights{i} = -666*ones(atb_gene.lists.numentries(i), 1);
    end
    
end

fclose(fid);


% merge up and down
gene_atb.edges = sourcetargetswap(lists2edges(atb_gene.lists));
clear atb_gene;
gene_atb.edges.target = cellfun(@(x) x(1:find(x=='_', 1, 'last')-1), gene_atb.edges.target, 'UniformOutput', false);
gene_atb.edges.targetname = 'GeneSym_PMID_Cell or Tissue_Perturbation_Organism_GEO Platform_GEO Accession';
discard = gene_atb.edges.weight == -666 | gene_atb.edges.weight == 0 | isnan(gene_atb.edges.weight);
gene_atb.edges = edgesdiscard(gene_atb.edges, discard);


% map gene symbols to entrez gene symbols and discard edges corresponding to un-mapped symbols
[gene_atb.edges.source, gene_atb.edges.sourcename, gene_atb.edges.sourceid, gene_atb.edges.sourceidname, discard] = genesymlookup(gene_atb.edges.source, [], 'gene', 'both', true, true, mappingfilespath);
gene_atb.edges = edgesdiscard(gene_atb.edges, discard);


% remove duplicate edges
gene_atb.edges = edgesunique(gene_atb.edges);


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
gene_atb.cm = cm_standardize_ignorezeros_signed(gene_atb.cm);

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


