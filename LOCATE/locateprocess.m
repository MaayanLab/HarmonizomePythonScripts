


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
gene_atb.edges = edgesinit(40961, [], 'HGNC ID', [], 'UniprotAcc', [], 'GeneID', [], 'Cellular Component', [], 'GO ID', [], [], false, []);


% read data
fid = fopen('input/parse_results_fixed.txt', 'r');

currline = fgetl(fid);

for i = 1:1:gene_atb.edges.numedges
    
    currline = fgetl(fid);
    
    currcells = strsplitbyadr(currline, '\t');
    
    gene_atb.edges.source{i} = ['HGNC:' currcells{3}];
    
    gene_atb.edges.sourcedesc{i} = currcells{4};
    
    gene_atb.edges.sourceid(i) = str2double(currcells{2});
    
    gene_atb.edges.target{i} = lower(currcells{7});
    
    gene_atb.edges.targetdesc{i} = currcells{6};
    
end

fclose(fid);

discard = (strcmp(gene_atb.edges.source, '-666') & strcmp(gene_atb.edges.sourcedesc, '-666') & gene_atb.edges.sourceid == -666) | strcmp(gene_atb.edges.targetdesc, '-666');
gene_atb.edges = edgesdiscard(gene_atb.edges, discard);
save('input/gene_cc_locate_curated_unprocessed_20150730.mat', '-struct', 'gene_atb');
%}

% load data
% gene_atb = load('input/gene_cc_locate_curated_unprocessed_20150730.mat', '-mat', 'edges');


% map gene identifiers to NCBI Entrez Gene Symbols and Gene IDs
uniprot_entrez = load('C:/Users/Andrew/Dropbox/GeneInfo/uniprot_entrez_resource_human_withisoforms.mat', '-mat');
gene_atb.edges.sourcedesc = cellfun(@(x) strsplitbyadr(x, '-'), gene_atb.edges.sourcedesc, 'uniformoutput', false);
gene_atb.edges.sourcedesc = cellfun(@(x) x{1}, gene_atb.edges.sourcedesc, 'uniformoutput', false);
[o1, o2] = ismember(gene_atb.edges.sourcedesc, uniprot_entrez.edges.source);
o2(o2 == 0) = [];

hgnc = load('C:/Users/Andrew/Dropbox/GeneInfo/hgnc_id_symbol.mat', '-mat');
[o3, o4] = ismember(gene_atb.edges.source, hgnc.termdesc);
o4(o4 == 0) = [];

entrez = load('C:\Users\Andrew\Dropbox\GeneInfo\gene_synonym_resource_human.mat', '-mat');
[o5, o6] = ismember(gene_atb.edges.sourceid, entrez.lists.termid);
o6(o6 == 0) = [];

gene_atb.edges.source = repmat({'-666'}, gene_atb.edges.numedges, 1);
gene_atb.edges.sourcename = 'GeneSym';
gene_atb.edges.sourceid = -666*ones([gene_atb.edges.numedges 1]);
gene_atb.edges.sourceidname = 'GeneID';
gene_atb.edges = rmfield(gene_atb.edges, {'sourcedesc' 'sourcedescname'});
gene_atb.edges.source(o1) = uniprot_entrez.edges.target(o2);
gene_atb.edges.source(o3) = hgnc.term(o4);
gene_atb.edges.source(o5) = entrez.lists.term(o6);
discard = ~(o1 | o3 | o5);
gene_atb.edges = edgesdiscard(gene_atb.edges, discard);

clear uniprot_entrez hgnc entrez;

[gene_atb.edges.source, gene_atb.edges.sourcename, gene_atb.edges.sourceid, gene_atb.edges.sourceidname, discard] = genesymlookup(gene_atb.edges.source, [], 'gene', 'human', true, false, mappingfilespath);
gene_atb.edges = edgesdiscard(gene_atb.edges, discard);


% map attributes to ontology (CHECK THIS!)
atbid_altid = load('output/gene_attribute_matrix_imported.mat', '-mat');

[o1, o2] = ismember(gene_atb.edges.targetdesc, atbid_altid.edges.targetdesc);
o2(o2 == 0) = [];

gene_atb.edges.target(o1) = atbid_altid.edges.source(o2);
gene_atb.edges.targetdesc(o1) = atbid_altid.edges.sourcedesc(o2);

atb_atb = load('output/gene_attribute_matrix_imported.mat', '-mat');

[o1, o2] = ismember(gene_atb.edges.targetdesc, atb_atb.cm.termdesc);
o2(o2 == 0) = [];

gene_atb.edges.targetname = atb_atb.cm.termname;
gene_atb.edges.target(o1) = atb_atb.cm.term(o2);
gene_atb.edges.targetdesc(o1) = atb_atb.cm.termdesc(o2);

gene_atb.edges = edgesdiscard(gene_atb.edges, ~o1);


% remove duplicate edges
gene_atb.edges = edgesunique(gene_atb.edges);


% % discard edges with weight below threshold
% % not initially sparse (49%)
% % can get sparsity down to 15% using threshold=30
% hist(gene_atb.edges.weight, 100);
% minweight = 30;
% discard = gene_atb.edges.weight < minweight;
% gene_atb.edges = edgesdiscard(gene_atb.edges, discard);
% gene_atb.edges = rmfield(gene_atb.edges, 'weight');


% convert to matrix format (attribute table)
gene_atb.cm = edges2cm(gene_atb.edges);
gene_atb = rmfield(gene_atb, 'edges');


% view distributions of row and col stats
[~] = cmrowcolstats(gene_atb.cm, true, 1);


% cluster and view matrix
[gene_atb.cm, ~] = cmcluster(gene_atb.cm, false);
% close force all;


% convert to sparse matrix
gene_atb.cm.matrix = sparse(gene_atb.cm.matrix);


% save result
save('output/gene_attribute_matrix_imported.mat', '-struct', 'gene_atb');






% ensure propagation of gene associations up ontology to ancestor terms
gene_atb.cm = cm2em_cosine(gene_atb.cm, atb_atb.cm, false, []);
gene_atb.cm.matrix = double(gene_atb.cm.matrix > 0);
gene_atb.cm = cmtrim(gene_atb.cm, 1, Inf, 1, Inf);


% view distributions of row and col stats
[~] = cmrowcolstats(gene_atb.cm, true, 2);


% cluster and view matrix
[gene_atb.cm, ~] = cmcluster(gene_atb.cm, false);
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
[gene_atb.cm, ~] = cmcluster(gene_atb.cm, true);
% close force all;


% view distributions of row and col stats
[~] = cmrowcolstats(gene_atb.cm, true, 3);


% convert to sparse matrix
gene_atb.cm.matrix = sparse(gene_atb.cm.matrix);


% save result
save('output/gene_attribute_matrix_prepared.mat', '-struct', 'gene_atb');


