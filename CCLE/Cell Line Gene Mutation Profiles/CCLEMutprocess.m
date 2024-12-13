


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




gene_cl.edges = edgesinit(122566, [], 'GeneSym', [], [], [], 'GeneID', [], 'CellLine', [], 'Tissue', [], [], false, []);

fid = fopen('input/CCLE_hybrid_capture1650_hg19_NoCommonSNPs_NoNeutralVariants_CDS_2012.05.07.maf', 'r');

currline = fgetl(fid);

tic;
for i = 1:1:gene_cl.edges.numedges
    
    currline = upper(fgetl(fid));
    
    currcells = strsplit(currline, '\t', 'CollapseDelimiters', false);
    
    gene_cl.edges.source{i} = currcells{1};
    gene_cl.edges.sourceid(i) = str2double(currcells{2});
    
    subcells = strsplit(currcells{16}, '_');
    
    gene_cl.edges.target{i} = subcells{1};
    gene_cl.edges.targetdesc{i} = lower(strjoin(subcells(2:end), '_'));
    
end

fclose(fid);

gene_cl.cm = edges2cm(gene_cl.edges);

gene_cl = rmfield(gene_cl, 'edges');

% map gene symbols to entrez gene symbols and discard rows corresponding to un-mapped symbols
[gene_cl.cm.term, gene_cl.cm.termname, gene_cl.cm.termid, gene_cl.cm.termidname, discard] = genesymlookup([], gene_cl.cm.termid, 'gene', 'human', false, false, mappingfilespath);
gene_cl.cm = cmrowdiscard(gene_cl.cm, discard);

% cluster and view data
cgo = cm2clustergram(gene_cl.cm, 'none', 'all', 'cosine', 'average');

close force all;

gene_cl.cm = conmatmap(gene_cl.cm, cgo.RowLabels, cgo.ColumnLabels');

HeatMap(gene_cl.cm.matrix, 'Colormap', redbluecmap);

% save result
save('output/gene_attribute_matrix_imported.mat', '-struct', 'gene_cl');


