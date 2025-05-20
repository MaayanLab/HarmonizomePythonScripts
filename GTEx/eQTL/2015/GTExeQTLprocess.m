


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
% get gene_snp, then get snp_snp, which is based on genomic distance, then
% relate genes to eachother based on closeness of their snps



% get eQTL
%{
gene_atb.edges = edgesinit(10000, [], 'Ensemble Acc', [], [], [], [], [], 'rsID', [], [], [], [], true, []);
% gene_atb.edges = edgesinit(11000000, [], 'Ensemble Acc', [], [], [], [], [], 'rsID', [], [], [], [], true, []);

fid = fopen('input/Multi_tissue_eQTL_GTEx_Pilot_Phase_datasets/res_final_amean_com_genes_com_snps.txt', 'r');
% fid = fopen('input/Multi_tissue_eQTL_GTEx_Pilot_Phase_datasets/res_final_amean_com_genes_com_snps_all.txt', 'r');

currline = fgetl(fid);

i = 0;

while ~feof(fid)
    
    i = i + 1;
    
    currline = fgetl(fid);
    
    currcells = strsplit(currline, '\t', 'CollapseDelimiters', false);
    
    gene_atb.edges.source{i} = currcells{1};
    
    gene_atb.edges.target{i} = currcells{2};
    
    gene_atb.edges.weight(i) = mean(str2double(currcells(3:end)));
    
end

fclose(fid);

if i < gene_atb.edges.numedges
    discard = false([gene_atb.edges.numedges 1]);
    discard(i+1:end) = true;
    gene_atb.edges = edgesdiscard(gene_atb.edges, discard);
end

save('input/gene_bestsnp_gtex_unmodified.mat', '-struct', 'gene_atb');
%}
gene_atb = load('input/gene_bestsnp_gtex_unmodified.mat', '-mat', 'edges');



% get SNP info
%{
snp.numsnps = 6820471;
snp.term = cell(snp.numsnps, 1);
snp.termname = 'rsid';
snp.termdesc = cell(snp.numsnps, 1);
snp.termdescname = 'chromosome';
snp.termid = zeros([snp.numsnps 1]);
snp.termidname = 'position';

fid = fopen('input/GTEx_genot_imputed_variants_info4_maf05_CR95_CHR_POSb37_ID_REF_ALT.txt', 'r');

currline = fgetl(fid);

for i = 1:1:snp.numsnps
    
    currline = fgetl(fid);
    
    currcells = strsplit(currline, ' ', 'CollapseDelimiters', false);
    
    snp.term{i} = currcells{3};
    
    snp.termdesc{i} = ['chr' currcells{1}];
    
    snp.termid(i) = str2double(currcells{2});
    
end

fclose(fid);

save('input/snp_list_gtex_unmodified.mat', '-struct', 'snp');
%}
snp = load('input/snp_list_gtex_unmodified.mat', '-mat');



% map SNP info (chromosome and position) to SNPs
[o1, o2] = ismember(gene_atb.edges.target, snp.term);
o2(o2 == 0) = [];

gene_atb.edges.targetdesc = repmat({'-666'}, gene_atb.edges.numedges, 1);
gene_atb.edges.targetid = -666*ones([gene_atb.edges.numedges 1]);

gene_atb.edges.targetdescname = snp.termdescname;
gene_atb.edges.targetdesc(o1) = snp.termdesc(o2);
gene_atb.edges.targetidname = snp.termidname;
gene_atb.edges.targetid(o1) = snp.termid(o2);

clear snp;



% map gene symbols to ensemble accession numbers
load('input/gene_tissuesample_gtex_unmodified.mat', '-mat', 'cm');

ens_qry = cellfun(@(x) x(1:find(x == '.', 1, 'first')-1), gene_atb.edges.source, 'UniformOutput', false);
numel(unique(ens_qry)) == numel(ens_qry) % make sure removing decimal doesn't affect mapping

ens_ref = cellfun(@(x) x(1:find(x == '.', 1, 'first')-1), cm.termdesc, 'UniformOutput', false);
numel(unique(ens_ref)) == numel(ens_ref) % make sure removing decimal doesn't affect mapping

[o1, o2] = ismember(ens_qry, ens_ref);
o2(o2 == 0) =[];

gene_atb.edges.sourcedesc = gene_atb.edges.source;
gene_atb.edges.sourcedescname = gene_atb.edges.sourcename;

gene_atb.edges.sourcename = 'GeneSym';
gene_atb.edges.source(o1) = cm.term(o2);

missed = gene_atb.edges.source(~o1);
gene_atb.edges = edgesdiscard(gene_atb.edges, ~o1);

clear cm ens_sqry ens_ref o1 o2;



% convert to connectivity matrix, 9644 x 9397 (very few genes share SNPs!)
gene_atb.cm = edges2cm(gene_atb.edges);
gene_atb = rmfield(gene_atb, 'edges');



% map gene symbols to entrez gene symbols and discard rows corresponding to
% un-mapped symbols (1740 rows discarded, 7904 x 9397 remaining)
[gene_atb.cm.term, gene_atb.cm.termname, gene_atb.cm.termid, gene_atb.cm.termidname, discard] = genesymlookup(gene_atb.cm.term, [], 'gene', 'human', true, false, mappingfilespath);
missed = gene_atb.cm.term(discard);
gene_atb.cm = cmrowdiscard(gene_atb.cm, discard);



if numel(unique(gene_atb.cm.term)) < gene_atb.cm.numterms % TRUE
    % merge measurements corresponding to the same gene (7898 unique gene
    % symbols out of 7904 total gene symbols, 7898 x 9397 remaining)
    gene_atb.cm = cmrowmerge(gene_atb.cm, 'max');
end




% remove empty rows and cols, 7898 x 7815 remaining
gene_atb.cm = cmtrim(gene_atb.cm, 1, Inf, 1, Inf);



sm = cm2sm_cosine_nocluster(gene_atb.cm);
numidentical = sum(sm.matrix > (1 - 1e-12), 2) - 1;
suspectgenes = sm.term(numidentical > 0);
clear sm;
% if ~isempty(suspectgenes) % TRUE, but can't think of a good rationale for doing this on this data
%     % remove suspect rows that have identical values to at least one other row
%     discard = ismember(gene_atb.cm.term, suspectgenes);
%     gene_atb.cm = cmrowdiscard(gene_atb.cm, discard);
% end



sm = cm2sm_cosine_nocluster(cmtranspose(gene_atb.cm));
numidentical = sum(sm.matrix > (1 - 1e-12), 2) - 1;
suspectattributes = sm.term(numidentical > 0);
clear sm;
% if ~isempty(suspectattributes) % TRUE, but can't think of a good rationale for doing this on this data
%     % remove suspect columns that have identical values to at least one
%     % other column
%     discard = ismember(gene_atb.cm.entry, suspectattributes);
%     gene_atb.cm = cmcoldiscard(gene_atb.cm, discard);
% end




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



% remove rows and cols with too many connections, 7898 x 7815 remaining
threshfrac = 1/4;
gene_atb.cm = cmtrim_frac(gene_atb.cm, 0, threshfrac, 0, threshfrac, 'column');
gene_atb.cm = cmtrim(gene_atb.cm, 1, Inf, 1, Inf);



% cluster and view matrix
[gene_atb.cm, ~] = cmcluster(gene_atb.cm, true);



% view distributions of row and col stats
[~] = cmrowcolstats(gene_atb.cm, true, 7);



% convert to sparse matrix
gene_atb.cm.matrix = sparse(gene_atb.cm.matrix);



% save result
save('output/gene_attribute_matrix_prepared.mat', '-struct', 'gene_atb');



atb_atb.cm = cminit(gene_atb.cm.numentries, gene_atb.cm.numentries, gene_atb.cm.entry, gene_atb.cm.entryname, gene_atb.cm.entrydesc, gene_atb.cm.entrydescname, gene_atb.cm.entryid, gene_atb.cm.entryidname, gene_atb.cm.entry, gene_atb.cm.entryname, gene_atb.cm.entrydesc, gene_atb.cm.entrydescname, gene_atb.cm.entryid, gene_atb.cm.entryidname, []);

for i = 1:1:atb_atb.cm.numterms
        
    atb_atb.cm.matrix(i,:) = (abs(atb_atb.cm.entryid - atb_atb.cm.termid(i))./strcmp(atb_atb.cm.entrydesc, atb_atb.cm.termdesc{i}))';
    
end

figure(101);
clf;
hist(atb_atb.cm.matrix(~isinf(atb_atb.cm.matrix)), 100);

mu = mean(atb_atb.cm.matrix(~isinf(atb_atb.cm.matrix)));
sigma = std(atb_atb.cm.matrix(~isinf(atb_atb.cm.matrix)));

% !!! probably have to define a characteristic distance for each chromosome!!!

figure(102);
clf;
hist(exp(-atb_atb.cm.matrix(~isinf(atb_atb.cm.matrix)).^2/2/sigma^2), 100);

figure(103);
clf;
hist(exp(-atb_atb.cm.matrix(~isinf(atb_atb.cm.matrix))/mu), 100);

atb_atb.sm = rmfield(atb_atb.cm, {'entry' 'entryname' 'entrydesc' 'entrydescname' 'entryid' 'entryidname' 'numentries'});
atb_atb.sm.matrix = exp(-atb_atb.cm.matrix.^2/2/(sigma)^2);
atb_atb.sm = smcluster(atb_atb.sm);
HeatMap(atb_atb.sm.matrix, 'Colormap', redbluecmap);


