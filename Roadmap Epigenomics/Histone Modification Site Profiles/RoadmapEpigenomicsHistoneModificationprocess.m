


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
% initialize lists structure
atb_gene.lists = listsinit(383, [], 'Histone Modification_Cell Type', [], [], [], [], [], 'GeneSym', [], [], [], [], true, []);


% read data
fid = fopen('input/Epigenomics_Roadmap_HM_ChIP-seq_fromEnrichr.txt', 'r');

for i = 1:1:atb_gene.lists.numterms
    
    currline = fgetl(fid);

    currcells = strsplitbyadr(currline, '\t');
    
    si = regexp(currcells{1}, ' ');
    hm = currcells{1}(1:si(1)-1);
    ct = currcells{1}(si(1)+1:end);
    atb_gene.lists.term{i} = [hm '_' ct];
    
    gv = cellfun(@(x) {x(1:regexp(x, ',')-1) str2double(x(regexp(x, ',')+1:end))}, currcells(3:end-1)', 'UniformOutput', false);
    atb_gene.lists.entries{i} = cellfun(@(x) x{1}, gv, 'UniformOutput', false);
    atb_gene.lists.weights{i} = cellfun(@(x) x{2}, gv);
    
    atb_gene.lists.numentries(i) = numel(atb_gene.lists.entries{i});
    
end

fclose(fid);


% convert to edges format
gene_atb.edges = sourcetargetswap(lists2edges(atb_gene.lists));
clear atb_gene;
discard = isnan(gene_atb.edges.weight) | gene_atb.edges.weight == 0;
gene_atb.edges = edgesdiscard(gene_atb.edges, discard);
save('input/gene_hm_roadmapepigenomics_nomerge_unprocessed_20150730.mat', '-struct', 'gene_atb');
%}

% load data
gene_atb = load('input/gene_hm_roadmapepigenomics_nomerge_unprocessed_20150730.mat', '-mat', 'edges');


% map gene symbols to entrez gene symbols and discard edges corresponding to un-mapped symbols
[gene_atb.edges.source, gene_atb.edges.sourcename, gene_atb.edges.sourceid, gene_atb.edges.sourceidname, discard] = genesymlookup(gene_atb.edges.source, [], 'gene', 'human', true, false, mappingfilespath);
gene_atb.edges = edgesdiscard(gene_atb.edges, discard);


% remove duplicate edges
gene_atb.edges = edgesunique(gene_atb.edges);


% discard edges with weight below threshold
% not initially sparse (20%) so use threshold to improve sparsity
% similar to sparsity of conserved matrix (~2.25%)
% at min weight 2^0=1, have sparsity 11% (warning: throwing away nearly half the values!)
% hist(log2(gene_atb.edges.weight), 1000);
minweight = 2^0;
discard = gene_atb.edges.weight < minweight;
gene_atb.edges = edgesdiscard(gene_atb.edges, discard);
gene_atb.edges = rmfield(gene_atb.edges, 'weight');


% convert to matrix format (attribute table)
gene_atb.cm = edges2cm(gene_atb.edges);
gene_atb = rmfield(gene_atb, 'edges');


% view distributions of row and col stats
[~] = cmrowcolstats(gene_atb.cm, true, 1);


% cluster and view matrix
[gene_atb.cm, ~] = cmcluster(gene_atb.cm, false);
close force all;


% convert to sparse matrix
gene_atb.cm.matrix = sparse(gene_atb.cm.matrix);


% save result
save('output/gene_attribute_matrix_imported.mat', '-struct', 'gene_atb');






clearvars -except gene_atb;


% convert to full matrix
gene_atb.cm.matrix = full(gene_atb.cm.matrix);


% remove rows and cols with too many connections
threshfrac = 1/2;
gene_atb.cm = cmtrim_frac(gene_atb.cm, 0, threshfrac, 0, threshfrac, 'column');
gene_atb.cm = cmtrim(gene_atb.cm, 1, Inf, 1, Inf);


% cluster and view matrix
% [gene_atb.cm, ~] = cmcluster(gene_atb.cm, true);
% close force all;


% view distributions of row and col stats
% [~] = cmrowcolstats(gene_atb.cm, true, 2);


% convert to sparse matrix
gene_atb.cm.matrix = sparse(gene_atb.cm.matrix);


% save result
save('output/gene_attribute_matrix_prepared.mat', '-struct', 'gene_atb');


