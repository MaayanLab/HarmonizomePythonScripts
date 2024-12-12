


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
atb_gene.lists = listsinit(353, [], 'GeneSym_PMID_CellType_Organism', [], 'GeneSym', [], [], [], 'GeneSym', [], [], [], [], false, []);

% read data
fid = fopen('input/dataset_20141013_original.gmt', 'r');
for i = 1:1:atb_gene.lists.numterms
    currline = fgetl(fid);
    currcells = cellfun(@strtrim, strsplitbyadr(currline, '\t'), 'uniformoutput', false);
    subcells = cellfun(@strtrim, strsplitbyadr(currcells{1}, '-'), 'uniformoutput', false);
    for j = 1:1:numel(subcells)
        match = regexp(subcells{j}, '(?<pmid>\d+)', 'names');
        if ~isempty(match) && numel(match.pmid) > 4 && numel(match.pmid) < 12 && numel(match.pmid) == numel(subcells{j})
            pmididx = j;
            break
        end
    end
    genesym = upper(strjoin(subcells(1:pmididx-1), '-'));
    pmid = subcells{pmididx};
    organism = lower(subcells{end});
    celltype = strjoin(subcells(pmididx+1:end-1), '-');
    atb_gene.lists.termdesc{i} = genesym;
    atb_gene.lists.term{i} = [genesym '_' pmid '_' celltype '_' organism];
    atb_gene.lists.entries{i} = currcells(3:end)';
    atb_gene.lists.numentries(i) = numel(atb_gene.lists.entries{i});
end
fclose(fid);
discard = atb_gene.lists.numentries == 0 | cellfun(@isempty, atb_gene.lists.term) | ismember(atb_gene.lists.term, {'' '-' '-666'}) | cellfun(@isempty, atb_gene.lists.termdesc) | ismember(atb_gene.lists.termdesc, {'' '-' '-666'});
atb_gene.lists = listsdiscard(atb_gene.lists, discard);
clear fid ans currcells currline discard i celltype genesym j match organism pmid pmididx subcells;

% convert to matrix format
gene_atb.cm = cmtranspose(edges2cm(edgesunique(lists2edges(atb_gene.lists))));
clear atb_gene;

% map gene identifiers to entrez gene symbols and discard rows corresponding to un-mapped symbols
% 18138 unmapped out of 47722 (38%) lots of garbage symbols from other species. why???
[gene_atb.cm.term, gene_atb.cm.termname, gene_atb.cm.termid, gene_atb.cm.termidname, discard] = genesymlookup(gene_atb.cm.term, [], 'gene', 'both', true, true, mappingfilespath);
gene_atb.cm = cmrowdiscard(gene_atb.cm, discard);
disp([num2str(sum(discard)) ' unmapped out of ' num2str(numel(discard)) ' (' num2str(round(100*sum(discard)/numel(discard))) '%)']);
clear discard;

% map attribute identifiers to entrez gene symbols and discard rows corresponding to un-mapped symbols
% 7 unmapped out of 353 (2%)
gene_atb.cm.entry = gene_atb.cm.entrydesc;
gene_atb.cm.entryname = gene_atb.cm.entrydescname;
gene_atb.cm = rmfield(gene_atb.cm, {'entrydesc' 'entrydescname'});
[gene_atb.cm.entry, gene_atb.cm.entryname, gene_atb.cm.entryid, gene_atb.cm.entryidname, discard] = genesymlookup(gene_atb.cm.entry, [], 'gene', 'both', true, true, mappingfilespath);
gene_atb.cm = cmcoldiscard(gene_atb.cm, discard);
disp([num2str(sum(discard)) ' unmapped out of ' num2str(numel(discard)) ' (' num2str(round(100*sum(discard)/numel(discard))) '%)']);
clear discard;

% merge rows corresponding to the same gene
if numel(unique(gene_atb.cm.term)) < gene_atb.cm.numterms
    gene_atb.cm = cmrowmerge(gene_atb.cm, 'max');
end

% merge cols corresponding to the same gene
if numel(unique(gene_atb.cm.entry)) < gene_atb.cm.numentries
    gene_atb.cm = cmcolmerge(gene_atb.cm, 'mean'); % note using mean here
end

% remove empty rows and cols
gene_atb.cm = cmtrim(gene_atb.cm, 1, Inf, 1, Inf);

% get columnwise thresholds
targetsetsize = 750;
numt = zeros(gene_atb.cm.numentries, 1);
tval = zeros(gene_atb.cm.numentries, 1000);
numg = zeros(gene_atb.cm.numentries, 1000);
numc = zeros(gene_atb.cm.numentries, 1);
for i = 1:1:gene_atb.cm.numentries
    t = unique(gene_atb.cm.matrix(:,i));
    numt(i) = numel(t);
    tval(i,1:numt(i)) = t';
    for j = 1:1:numt(i)
        numg(i,j) = sum(gene_atb.cm.matrix(:,i) > t(j));
    end
    d2 = (numg(i,:) - targetsetsize).^2;
    numc(i) = mean(find(d2 == min(d2)));
end
[unumt, ui, ri] = unique(numt);
unumc = zeros(numel(unumt), 1);
for i = 1:1:numel(unumt)
    unumc(i) = nanmedian(numc(ri==i));
end
figure(3);
clf;
plot(unumt-1, unumc, '-ok');
xlabel('number of profiles');
ylabel(['median consensus threshold to get around ' num2str(targetsetsize) ' genes']);
axis([0 10 0 10]);
x = unumt(~isnan(unumc))-1;
y = unumc(~isnan(unumc));
y(x<3) = [];
x(x<3) = [];
numprofiles = numt-1;
p = polyfit(log10(x), log10(y), 1);
numconsensus = round(10.^(polyval(p, log10(numprofiles))));
% numconsensus = round(10.^(spline(log10(x),log10(y),log10(numprofiles))));
% % % numconsensus = round(10.^(csaps(log10(x),log10(y),1,log10(numprofiles))));
% sp = spaps(log10(x), log10(y), 0);
% numconsensus = round(10.^(fnval(sp, log10(numprofiles))));
numconsensus(numconsensus<1) = 1;
figure(3);
hold on;
plot(numprofiles, numconsensus, 'sr');
hold off;
axis([1 10 1 10]);
threshold = zeros(gene_atb.cm.numentries, 1);
for i = 1:1:gene_atb.cm.numentries
    threshold(i) = tval(i,numconsensus(i));
    gene_atb.cm.matrix(gene_atb.cm.matrix(:,i) <= threshold(i),i) = 0;
end
discard = numprofiles < 3;
gene_atb.cm = cmtrim(cmcoldiscard(gene_atb.cm, discard), 1, Inf, 1, Inf);
figure(4); clf; hist(sum(gene_atb.cm.matrix ~= 0, 1)', 10); xlabel('targets'); ylabel('TFs');
figure(5); clf; hist(log10(sum(gene_atb.cm.matrix ~= 0, 1))', 10); xlabel('targets'); ylabel('TFs');

% save cleaned data
gene_atb.cm.matrix = sparse(gene_atb.cm.matrix);
save('output/gene_attribute_matrix_cleaned.mat', '-struct', 'gene_atb');
gene_atb.cm.matrix = full(gene_atb.cm.matrix);

% standardize
tftype = 'raw'; % no need to normalize because all values in matrix are already on the same scale, ranging from 1 (gene is a hit in all profiles) to 0 (gene is a hit in no profiles)
gene_atb.cm = cm_tfidf_standardization(gene_atb.cm, tftype);
figure(1); clf; hist(gene_atb.cm.matrix(gene_atb.cm.matrix~=0), 100);
support = 'unbounded';
fignum = 2;
gene_atb.cm = cm_ksdensity_standardization_sparsematrix(gene_atb.cm, support, fignum);
clear tftype support fignum;

% cluster and view standardized matrix
gene_atb.cm = cmcluster_nogram(gene_atb.cm, true);

% save standardized data
gene_atb.cm.matrix = sparse(gene_atb.cm.matrix);
save('output/gene_attribute_matrix_standardized.mat', '-struct', 'gene_atb');
gene_atb.cm.matrix = full(gene_atb.cm.matrix);
writecm('output/gene_attribute_matrix_standardized', gene_atb.cm);
gzip('output/gene_attribute_matrix_standardized.txt');
delete('output/gene_attribute_matrix_standardized.txt');
genes = gene_atb.cm.term;
atbs = gene_atb.cm.entry;

% threshold
gene_atb.cm.matrix = (gene_atb.cm.matrix > 0) - (gene_atb.cm.matrix < 0);

% save thresholded data
gene_atb.cm.matrix = sparse(gene_atb.cm.matrix);
save('output/gene_attribute_matrix_thresholded.mat', '-struct', 'gene_atb');
gene_atb.cm.matrix = full(gene_atb.cm.matrix);
writecm('output/gene_attribute_matrix_thresholded', gene_atb.cm);
gzip('output/gene_attribute_matrix_thresholded.txt');
delete('output/gene_attribute_matrix_thresholded.txt');
clear gene_atb;

% align cleaned matrix with standardized matrix and save
load('output/gene_attribute_matrix_cleaned.mat', '-mat', 'cm');
cm = conmatmap(cm, genes, atbs);
save('output/gene_attribute_matrix_cleaned.mat', '-mat', 'cm');
writecm('output/gene_attribute_matrix_cleaned', cm);
gzip('output/gene_attribute_matrix_cleaned.txt');
delete('output/gene_attribute_matrix_cleaned.txt');
clear cm genes atbs mappingfilespath;


