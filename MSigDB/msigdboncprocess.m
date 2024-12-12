


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
atb_gene.lists = listsinit(181, [], 'GeneSym_Perturbation or Condition_Direction', [], 'GeneSym', [], [], [], 'GeneSym', [], [], [], [], true, []);


% read data
fid = fopen('input/oncogenicsigs_c6.all.v5.0.symbols.gmt', 'r');

for i = 1:1:atb_gene.lists.numterms
    
    currline = fgetl(fid);

    currcells = strsplitbyadr(currline, '\t');
    
    atb_gene.lists.term{i} = currcells{1};
    
    di = regexp(currcells{1}, '.');
    ui = regexp(currcells{1}, '_');
    if numel(di) == 0
        mi = ui;
    elseif numel(ui) == 0
        mi = di;
    else
        mi = min(di(1), ui(1));
    end
    atb_gene.lists.termdesc{i} = currcells{1}(1:mi);
    
    atb_gene.lists.entries{i} = currcells(3:end)';
    atb_gene.lists.numentries(i) = numel(atb_gene.lists.entries{i});
    
    if strcmp(currcells{1}(end-1:end), 'UP')
        atb_gene.lists.weights{i} = ones(atb_gene.lists.numentries(i), 1);
    elseif strcmp(currcells{1}(end-1:end), 'DN')
        atb_gene.lists.weights{i} = -1*ones(atb_gene.lists.numentries(i), 1);
    else
        atb_gene.lists.weights{i} = -666*ones(atb_gene.lists.numentries(i), 1);
    end
    
end

fclose(fid);


% merge up and down
atb_gene.edges = lists2edges(atb_gene.lists);
atb_gene.edges.source = cellfun(@(x) x(1:find(x=='_', 1, 'last')-1), atb_gene.edges.source, 'UniformOutput', false);
atb_gene.edges.sourcename = 'GeneSym_Perturbation or Condition';
discard = atb_gene.edges.weight == -666;
atb_gene.edges = edgesdiscard(atb_gene.edges, discard);


% convert to matrix format (attribute table)
gene_atb.cm = cmtranspose(edges2cm(atb_gene.edges));
clear atb_gene;


% map gene identifiers to NCBI Entrez Gene Symbols and Gene IDs
[gene_atb.cm.term, gene_atb.cm.termname, gene_atb.cm.termid, gene_atb.cm.termidname, discard] = genesymlookup(gene_atb.cm.term, [], 'gene', 'human', true, false, mappingfilespath);
gene_atb.cm = cmrowdiscard(gene_atb.cm, discard);


% merge rows corresponding to the same gene
if numel(unique(gene_atb.cm.term)) < gene_atb.cm.numterms
    gene_atb.cm = cmrowmerge(gene_atb.cm, 'max');
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


