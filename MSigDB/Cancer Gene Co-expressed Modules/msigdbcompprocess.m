


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




% initialize lists structure (atb=cancer gene, gene=coexpressed gene)
atb_gene.lists = listsinit(427, [], 'GeneSym_ExpressionDatasetID', [], 'GeneSym', [], [], [], 'GeneSym', [], [], [], [], false, []);


% read data
fid = fopen('input/cancerneighborhoods_c4.cgn.v5.0.symbols.gmt', 'r');

for i = 1:1:atb_gene.lists.numterms
    
    currline = fgetl(fid);

    currcells = strsplitbyadr(currline, '\t');
    
    ui = regexp(currcells{1}, '_');
    genesym = currcells{1}(ui+1:end);
    dataset = currcells{2}(1:ui-1);
    
    atb_gene.lists.term{i} = [genesym '_' dataset];
    atb_gene.lists.termdesc{i} = genesym;
    atb_gene.lists.entries{i} = currcells(3:end)';
    atb_gene.lists.numentries(i) = numel(atb_gene.lists.entries{i});
    
end

fclose(fid);


% convert to matrix format (attribute table)
gene_atb.cm = cmtranspose(lists2cm(atb_gene.lists));
clear atb_gene;


% map gene identifiers to NCBI Entrez Gene Symbols and Gene IDs
[gene_atb.cm.term, gene_atb.cm.termname, gene_atb.cm.termid, gene_atb.cm.termidname, discard] = genesymlookup(gene_atb.cm.term, [], 'gene', 'human', true, false, mappingfilespath);
gene_atb.cm = cmrowdiscard(gene_atb.cm, discard);


% merge rows corresponding to the same gene
if numel(unique(gene_atb.cm.term)) < gene_atb.cm.numterms
    gene_atb.cm = cmrowmerge(gene_atb.cm, 'max');
end


% map attribute identifiers to NCBI Entrez Gene Symbols and Gene IDs
gene_atb.cm.entry = gene_atb.cm.entrydesc;
gene_atb.cm.entryname = gene_atb.cm.entrydescname;
gene_atb.cm = rmfield(gene_atb.cm, {'entrydesc' 'entrydescname'});
[gene_atb.cm.entry, gene_atb.cm.entryname, gene_atb.cm.entryid, gene_atb.cm.entryidname, discard] = genesymlookup(gene_atb.cm.entry, [], 'gene', 'human', true, false, mappingfilespath);
gene_atb.cm = cmcoldiscard(gene_atb.cm, discard);


% merge cols corresponding to the same gene
if numel(unique(gene_atb.cm.entry)) < gene_atb.cm.numentries
    gene_atb.cm = cmcolmerge(gene_atb.cm, 'max');
end


% remove empty rows and cols
gene_atb.cm = cmtrim(gene_atb.cm, 1, Inf, 1, Inf);


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


