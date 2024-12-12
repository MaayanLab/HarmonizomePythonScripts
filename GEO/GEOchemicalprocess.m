


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
atb_gene.lists = listsinit(700, [], 'Drug Perturbation_Organism_GEO Platform_GEO Accession_Direction', [], 'Drug', [], [], [], 'GeneSym', [], [], [], [], true, []);


% read data
fid = fopen('input/Drug_Perturbations_from_GEO.txt', 'r');

for i = 1:1:atb_gene.lists.numterms
    
    currline = fgetl(fid);

    currcells = strsplitbyadr(currline, '\t');
    
    currcells{1} = strrep(currcells{1}, '_chdir', '');
    
    atb_gene.lists.term{i} = currcells{1};
    
    ui = regexp(currcells{1}, '_');
    atb_gene.lists.termdesc{i} = currcells{1}(1:ui(1)-1);
    
    atb_gene.lists.entries{i} = currcells(3:end-1)';
    atb_gene.lists.numentries(i) = numel(atb_gene.lists.entries{i});
    
    if strcmp(currcells{1}(ui(end)+1:end), 'up')
        atb_gene.lists.weights{i} = ones(atb_gene.lists.numentries(i), 1);
    elseif strcmp(currcells{1}(ui(end)+1:end), 'down')
        atb_gene.lists.weights{i} = -1*ones(atb_gene.lists.numentries(i), 1);
    else
        atb_gene.lists.weights{i} = -666*ones(atb_gene.lists.numentries(i), 1);
    end
    
end

fclose(fid);


% merge up and down
gene_atb.edges = sourcetargetswap(lists2edges(atb_gene.lists));
clear atb_gene;
gene_atb.edges.target = cellfun(@(x) x(1:find(x=='_', 1, 'last')-1), gene_atb.edges.target, 'UniformOutput', false);
gene_atb.edges.targetname = 'Drug Perturbation_Organism_GEO Platform_GEO Accession';
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


