


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





% load data (going to use signatures thresholded at 300 up and 300 down)
v = load('input/virus_degene_resource_geo.mat', '-mat', 'ud');
e = load('input/virus_degene_resource_geo_ebola.mat', '-mat', 'ud');
ud.lists = listsvertcat(v.ud.lists, e.ud.lists);
clear v e;

% discard geo accessions that didn't profile enough genes
figure(1); clf; hist(ud.lists.numentries, 100);
badgse = ud.lists.termid(ud.lists.numentries < 7500); % 11300);
discard = ismember(ud.lists.termid, badgse);
atb_gene.lists = listsdiscard(ud.lists, discard);
clear ud;


% append tissue and geo accessions to attribute labels
term = atb_gene.lists.term;
termname = atb_gene.lists.termname;
for i = 1:1:atb_gene.lists.numterms
    atb_gene.lists.term{i} = [atb_gene.lists.term{i} '_' atb_gene.lists.termdesc{i} '_GSE' num2str(atb_gene.lists.termid(i))];
end
atb_gene.lists.termname = 'Virus Perturbation_PMID_GEO Accession';
atb_gene.lists.termdesc = term;
atb_gene.lists.termdescname = termname;


% convert to edges format
gene_atb.edges = sourcetargetswap(lists2edges(atb_gene.lists));
clear atb_gene;
discard = isnan(gene_atb.edges.weight) | gene_atb.edges.weight == 0;
gene_atb.edges = edgesdiscard(gene_atb.edges, discard);


% map gene symbols to entrez gene symbols and discard edges corresponding to un-mapped symbols
[gene_atb.edges.source, gene_atb.edges.sourcename, gene_atb.edges.sourceid, gene_atb.edges.sourceidname, discard] = genesymlookup(gene_atb.edges.source, [], 'gene', 'human', true, false, mappingfilespath);
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


