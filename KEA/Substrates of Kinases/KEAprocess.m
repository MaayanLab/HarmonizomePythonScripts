


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




% initialize lists structure (term=kinase, entries=substrate)
lists = listsinit(428, [], 'GeneSym', [], [], [], [], [], 'GeneSym', [], [], [], [], false, []);


% read data
fid = fopen('input/KEA_2015_from_enrichr.txt', 'r');

for i = 1:1:lists.numterms
    
    currline = fgetl(fid);
    
    currcells = strsplit(currline, '\t', 'CollapseDelimiters', false);

    lists.term{i} = currcells{1};
    
    entries = currcells(3:end)';
    entries = cellfun(@(x) x(1:find(x == ',', 1, 'first')-1), entries, 'UniformOutput', false);
    
    discard = strcmp(entries, '');

    lists.entries{i} = entries(~discard);
    lists.numentries(i) = numel(lists.entries{i});
    
end

fclose(fid);

discard = lists.numentries == 0;
lists = listsdiscard(lists, discard);

edges = lists2edges(lists);


% swap data so source=substrate and target=kinase
gene_atb.edges = sourcetargetswap(edges);
clearvars -except gene_atb;


% map gene identifiers to NCBI Entrez Gene Symbols and Gene IDs
[gene_atb.edges.source, gene_atb.edges.sourcename, gene_atb.edges.sourceid, gene_atb.edges.sourceidname, discard] = genesymlookup(gene_atb.edges.source, [], 'gene', 'both', true, true, mappingfilespath);
gene_atb.edges = edgesdiscard(gene_atb.edges, discard);


% map atb identifiers to NCBI Entrez Gene Symbols and Gene IDs
[gene_atb.edges.target, gene_atb.edges.targetname, gene_atb.edges.targetid, gene_atb.edges.targetidname, discard] = genesymlookup(gene_atb.edges.target, [], 'gene', 'both', true, true, mappingfilespath);
gene_atb.edges = edgesdiscard(gene_atb.edges, discard);


% discard duplicate edges
gene_atb.edges = edgesunique(gene_atb.edges);


% convert to matrix format (attribute table)
gene_atb.cm = edges2cm(gene_atb.edges);
gene_atb = rmfield(gene_atb, 'edges');


% view distributions of row and col stats
[~] = cmrowcolstats(gene_atb.cm, true, 1);


% cluster and view matrix
[gene_atb.cm, ~] = cmcluster(gene_atb.cm, true);
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


