


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





% load data (see get_crowdkd_cm.m for how this file was generated)
gene_atb = load('input/gene_perturbation_crowd_20150204.mat', '-mat', 'cm');


% append to attribute labels
for i = 1:1:gene_atb.cm.numentries
    gene_atb.cm.entry{i} = [gene_atb.cm.entry{i} '_' gene_atb.cm.entrydesc{i}];
    gene_atb.cm.entrydesc{i} = gene_atb.cm.entry{i}(1:find(gene_atb.cm.entry{i}=='_', 1, 'first')-1);
end
gene_atb.cm.entryname = 'GeneSym_Perturbation_GEO Accession_ID_Organism_Cell or Tissue Type';
gene_atb.cm.entrydescname = 'GeneSym';


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


