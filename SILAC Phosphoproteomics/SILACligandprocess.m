


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
atb_gene.lists = listsinit(84, [], 'Direction Perturbation_Cell Line (Organism) [PMID]', [], [], [], [], [], 'GeneSym', [], [], [], [], true, []);


% read data
fid = fopen('input/silac_experiment_protein_level_set_library_20150730.tsv', 'r');

for i = 1:1:atb_gene.lists.numterms
    
    currline = fgetl(fid);

    currcells = strsplitbyadr(currline, '\t');
    
    atb_gene.lists.term{i} = currcells{1};
%     atb_gene.lists.termdesc{i} = currcells{2};
    atb_gene.lists.entries{i} = currcells(3:end)';
    atb_gene.lists.numentries(i) = numel(atb_gene.lists.entries{i});
    
    if currcells{1}(1) == 'u'
        atb_gene.lists.weights{i} = ones(atb_gene.lists.numentries(i), 1);
    else
        atb_gene.lists.weights{i} = -1*ones(atb_gene.lists.numentries(i), 1);
    end
    
end

fclose(fid);


% merge up and down
atb_gene.edges = lists2edges(atb_gene.lists);
atb_gene.edges.source = cellfun(@(x) x(find(x==' ', 1, 'first')+1:end), atb_gene.edges.source, 'UniformOutput', false);
atb_gene.edges.sourcename = 'Perturbation_Cell Line (Organism) [PMID]';


% convert to matrix format (attribute table)
gene_atb.cm = cmtranspose(edges2cm(atb_gene.edges));
clear atb_gene;


% write attributes to file
% fid = fopen('input/attribute_list_20150730.txt', 'w');
% for i = 1:1:gene_atb.cm.numentries
%     fprintf(fid, '%s\r\n', gene_atb.cm.entry{i});
% end
% fclose(fid);


% read after assigning attribute categories
atb = cell(gene_atb.cm.numentries, 1);
atbcat = cell(gene_atb.cm.numentries, 1);
fid = fopen('input/attribute_list_20150730.txt', 'r');
for i = 1:1:gene_atb.cm.numentries
    currline = fgetl(fid);
    currcells = strsplitbyadr(currline, '\t');
    atb{i} = currcells{1};
    atbcat{i} = currcells{2};
end
fclose(fid);


% map attribute categories
gene_atb.cm.entrydesc = repmat({'-666'}, gene_atb.cm.numentries, 1);
gene_atb.cm.entrydescname = 'Perturbation Type';
[o1, o2] = ismember(gene_atb.cm.entry, atb);
o2(o2 == 0) = [];
gene_atb.cm.entrydesc(o1) = atbcat(o2);


% map gene identifiers to NCBI Entrez Gene Symbols and Gene IDs
[gene_atb.cm.term, gene_atb.cm.termname, gene_atb.cm.termid, gene_atb.cm.termidname, discard] = genesymlookup(gene_atb.cm.term, [], 'gene', 'both', true, true, mappingfilespath);
gene_atb.cm = cmrowdiscard(gene_atb.cm, discard);


% merge rows corresponding to the same gene
if numel(unique(gene_atb.cm.term)) < gene_atb.cm.numterms
    gene_atb.cm = cmrowmerge(gene_atb.cm, 'max');
end


% remove empty rows and cols
gene_atb.cm = cmtrim(gene_atb.cm, 1, Inf, 1, Inf);


% view distributions of row and col stats
% [~] = cmrowcolstats(gene_atb.cm, true, 1);


% cluster and view matrix
[gene_atb.cm, ~] = cmcluster(gene_atb.cm, false);
close force all;


% convert to sparse matrix
gene_atb.cm.matrix = sparse(gene_atb.cm.matrix);


% save result
save('output/gene_attribute_matrix_imported.mat', '-struct', 'gene_atb');
discard = ~strcmp(gene_atb.cm.entrydesc, 'drug perturbation');
gene_drug.cm = cmcoldiscard(gene_atb.cm, discard);
gene_drug.cm = cmtrim(gene_drug.cm, 1, Inf, 1, Inf);
save('output/gene_attribute_matrix_imported.mat', '-struct', 'gene_drug');
discard = ~strcmp(gene_atb.cm.entrydesc, 'gene perturbation');
gene_gene.cm = cmcoldiscard(gene_atb.cm, discard);
gene_gene.cm = cmtrim(gene_gene.cm, 1, Inf, 1, Inf);
save('output/gene_attribute_matrix_imported.mat', '-struct', 'gene_gene');
discard = ~strcmp(gene_atb.cm.entrydesc, 'ligand (protein) perturbation');
gene_ligand.cm = cmcoldiscard(gene_atb.cm, discard);
gene_ligand.cm = cmtrim(gene_ligand.cm, 1, Inf, 1, Inf);
save('output/gene_attribute_matrix_imported.mat', '-struct', 'gene_ligand');




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
discard = ~strcmp(gene_atb.cm.entrydesc, 'drug perturbation');
gene_drug.cm = cmcoldiscard(gene_atb.cm, discard);
gene_drug.cm = cmtrim(gene_drug.cm, 1, Inf, 1, Inf);
save('output/gene_attribute_matrix_prepared.mat', '-struct', 'gene_drug');
discard = ~strcmp(gene_atb.cm.entrydesc, 'gene perturbation');
gene_gene.cm = cmcoldiscard(gene_atb.cm, discard);
gene_gene.cm = cmtrim(gene_gene.cm, 1, Inf, 1, Inf);
save('output/gene_attribute_matrix_prepared.mat', '-struct', 'gene_gene');
discard = ~strcmp(gene_atb.cm.entrydesc, 'ligand (protein) perturbation');
gene_ligand.cm = cmcoldiscard(gene_atb.cm, discard);
gene_ligand.cm = cmtrim(gene_ligand.cm, 1, Inf, 1, Inf);
save('output/gene_attribute_matrix_prepared.mat', '-struct', 'gene_ligand');


