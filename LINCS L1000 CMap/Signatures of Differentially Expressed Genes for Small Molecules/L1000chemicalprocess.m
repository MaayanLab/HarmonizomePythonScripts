


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
atb_gene.uplists = listsinit(33166, [], 'Perturbation ID_Perturbagen_Cell Line_Time_Time Unit_Dose_Dose Unit', [], [], [], [], [], 'GeneSym', [], [], [], [], false, []);
atb_gene.dnlists = listsinit(33166, [], 'Perturbation ID_Perturbagen_Cell Line_Time_Time Unit_Dose_Dose Unit', [], [], [], [], [], 'GeneSym', [], [], [], [], false, []);


% read data
fid = fopen('input/lincs_l1000cds2_cpcd-gse70138_up.gmt', 'r');

for i = 1:1:atb_gene.uplists.numterms
    
    currline = strrep(fgetl(fid), '"', '');

    currcells = strsplitbyadr(currline, '\t');
    
    atb_gene.uplists.term{i} = [currcells{2} '_' currcells{1}];
%     atb_gene.uplists.termdesc{i} = currcells{2};
    atb_gene.uplists.entries{i} = currcells(3:end)';
    atb_gene.uplists.numentries(i) = numel(atb_gene.uplists.entries{i});
    
end

fclose(fid);

fid = fopen('input/lincs_l1000cds2_cpcd-gse70138_dn.gmt', 'r');

for i = 1:1:atb_gene.dnlists.numterms
    
    currline = strrep(fgetl(fid), '"', '');

    currcells = strsplitbyadr(currline, '\t');
    
    atb_gene.dnlists.term{i} = [currcells{2} '_' currcells{1}];
%     atb_gene.dnlists.termdesc{i} = currcells{2};
    atb_gene.dnlists.entries{i} = currcells(3:end)';
    atb_gene.dnlists.numentries(i) = numel(atb_gene.dnlists.entries{i});
    
end

fclose(fid);


% merge up and down
atb_gene.upedges = lists2edges(atb_gene.uplists);
atb_gene.upedges.weight = ones(atb_gene.upedges.numedges, 1);
atb_gene.dnedges = lists2edges(atb_gene.dnlists);
atb_gene.dnedges.weight = -1*ones(atb_gene.dnedges.numedges, 1);
atb_gene.edges = edgesvertcat(atb_gene.upedges, atb_gene.dnedges);


% convert to matrix format (attribute table)
gene_atb.cm = cmtranspose(edges2cm(atb_gene.edges));
clear atb_gene;


% map gene identifiers to NCBI Entrez Gene Symbols and Gene IDs
[gene_atb.cm.term, gene_atb.cm.termname, gene_atb.cm.termid, gene_atb.cm.termidname, discard] = genesymlookup(gene_atb.cm.term, [], 'gene', 'human', true, false, mappingfilespath);
gene_atb.cm = cmrowdiscard(gene_atb.cm, discard);


% merge rows corresponding to the same gene
if numel(unique(gene_atb.cm.term)) < gene_atb.cm.numterms
    gene_atb.cm = cmrowmerge(gene_atb.cm, 'mean');
end
gene_atb.cm.matrix(gene_atb.cm.matrix > 0) = 1;
gene_atb.cm.matrix(gene_atb.cm.matrix < 0) = -1;


% remove empty rows and cols
gene_atb.cm = cmtrim(gene_atb.cm, 1, Inf, 1, Inf);


% view distributions of row and col stats
% [~] = cmrowcolstats(gene_atb.cm, true, 1);


% cluster and view matrix
% [gene_atb.cm, ~] = cmcluster(gene_atb.cm, false);
% close force all;


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


