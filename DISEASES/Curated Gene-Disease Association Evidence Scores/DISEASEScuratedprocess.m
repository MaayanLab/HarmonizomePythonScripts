


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
gene_atb.edges = edgesinit(23477, [], 'GeneSym', [], 'Ensemble Acc', [], [], [], 'Disease', [], 'Disease Ontology ID', [], [], false, []);



% read file
fid = fopen('input/human_disease_knowledge_full_20150202.tsv', 'r');

for i = 1:1:gene_atb.edges.numedges
    
    currline = fgetl(fid);
    
    currcells = strsplit(currline, '\t', 'CollapseDelimiters', false);
    
    gene_atb.edges.source{i} = currcells{2};
    
    gene_atb.edges.sourcedesc{i} = currcells{1};
    
    gene_atb.edges.target{i} = lower(currcells{4});
    
    gene_atb.edges.targetdesc{i} = currcells{3};
    
end

fclose(fid);



% convert to connectivity matrix
gene_atb.cm = edges2cm(gene_atb.edges);
gene_atb = rmfield(gene_atb, 'edges');



% map gene symbols to entrez gene symbols and discard rows corresponding to
% un-mapped symbols (34 rows discarded, 2263 x 770 remaining)
[gene_atb.cm.term, gene_atb.cm.termname, gene_atb.cm.termid, gene_atb.cm.termidname, discard] = genesymlookup(gene_atb.cm.term, [], 'gene', 'human', true, false, mappingfilespath);
gene_atb.cm = cmrowdiscard(gene_atb.cm, discard);



if numel(unique(gene_atb.cm.term)) < gene_atb.cm.numterms % TRUE
    % merge measurements corresponding to the same gene (2252 unique gene
    % symbols out of 2263 total gene symbols, 2252 x 770 remaining)
    gene_atb.cm = cmrowmerge(gene_atb.cm, 'max');
end



% map diseases to disease ontology and discard un-mapped diseases (2 cols
% discarded, 2252 x 768 remaining)
atb_atb = load('output/gene_attribute_matrix_imported.mat', '-mat');
id_alt = load('output/gene_attribute_matrix_imported.mat', '-mat');

[o1, o2] = ismember(gene_atb.cm.entrydesc, id_alt.edges.targetdesc);
o2(o2 == 0) = [];

[o3, o4] = ismember(gene_atb.cm.entrydesc, atb_atb.cm.entrydesc);
o4(o4 == 0) = [];

gene_atb.cm.entry(o1) = id_alt.edges.source(o2);
gene_atb.cm.entrydesc(o1) = id_alt.edges.sourcedesc(o2);
gene_atb.cm.entry(o3) = atb_atb.cm.entry(o4);
gene_atb.cm.entrydesc(o3) = atb_atb.cm.entrydesc(o4);

discard = ~o1 & ~o3;
gene_atb.cm = cmcoldiscard(gene_atb.cm, discard);



if numel(unique(gene_atb.cm.entry)) < gene_atb.cm.numentries % TRUE
    % merge measurements corresponding to the same disease (764 unique gene
    % symbols out of 768 total gene symbols, 2252 x 764 remaining)
    gene_atb.cm = cmcolmerge(gene_atb.cm, 'max');
end



% remove empty rows and cols, 2252 x 759 remaining
gene_atb.cm = cmtrim(gene_atb.cm, 1, Inf, 1, Inf);



sm = cm2sm_cosine_nocluster(gene_atb.cm);
numidentical = sum(sm.matrix > (1 - 1e-12), 2) - 1;
suspectgenes = sm.term(numidentical > 0);
clear sm;
% if ~isempty(suspectgenes) % TRUE, but can't think of a good rationale for doing this on curated data
%     % remove suspect rows that have identical values to at least one other row
%     discard = ismember(gene_atb.cm.term, suspectgenes);
%     gene_atb.cm = cmrowdiscard(gene_atb.cm, discard);
% end



sm = cm2sm_cosine_nocluster(cmtranspose(gene_atb.cm));
numidentical = sum(sm.matrix > (1 - 1e-12), 2) - 1;
suspectattributes = sm.term(numidentical > 0);
clear sm;
% if ~isempty(suspectattributes) % TRUE, but can't think of a good rationale for doing this on curated data
%     % remove suspect columns that have identical values to at least one
%     % other column
%     discard = ismember(gene_atb.cm.entry, suspectattributes);
%     gene_atb.cm = cmcoldiscard(gene_atb.cm, discard);
% end



% make sure all appropriate connections to parent diseases are made, 2252 x
% 770 remaining
gene_atb.cm = cm2em_cosine(gene_atb.cm, atb_atb.cm, false, []);
gene_atb.cm.matrix = double(gene_atb.cm.matrix > 0);
gene_atb.cm = cmtrim(gene_atb.cm, 1, Inf, 1, Inf);



% view distributions of row and col stats
[~] = cmrowcolstats(gene_atb.cm, true, 4);



% cluster and view matrix
[gene_atb.cm, ~] = cmcluster(gene_atb.cm, true);



% convert to sparse matrix
gene_atb.cm.matrix = sparse(gene_atb.cm.matrix);



% save result
save('output/gene_attribute_matrix_imported.mat', '-struct', 'gene_atb');






clearvars -except gene_atb;

% convert to full matrix
gene_atb.cm.matrix = full(gene_atb.cm.matrix);



% remove rows and cols with too many connections, 2252 x 767 remaining
threshfrac = 1/2;
gene_atb.cm = cmtrim_frac(gene_atb.cm, 0, threshfrac, 0, threshfrac, 'column');
gene_atb.cm = cmtrim(gene_atb.cm, 1, Inf, 1, Inf);



% cluster and view matrix
% [gene_atb.cm, ~] = cmcluster(gene_atb.cm, true);



% view distributions of row and col stats
[~] = cmrowcolstats(gene_atb.cm, true, 7);



% convert to sparse matrix
gene_atb.cm.matrix = sparse(gene_atb.cm.matrix);



% save result
save('output/gene_attribute_matrix_prepared.mat', '-struct', 'gene_atb');


