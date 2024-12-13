


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



%{
% contents = dir('input/ExpressionMatrices');
% cancer = {contents(3:end).name}';
% numcancers = numel(cancer);
fid = fopen('input/ExpressionMatrices/Cancer_Acronyms.txt', 'r');
cancer.numterms = 25;
cancer.term = cell(cancer.numterms, 1);
cancer.termname = 'Cancer Name';
cancer.termdesc = cell(cancer.numterms, 1);
cancer.termdescname = 'Cancer Acronym';
for i = 1:1:cancer.numterms
    currline = fgetl(fid);
    currcells = strsplitbyadr(currline, '\t');
    cancer.term{i} = currcells{2};
    cancer.termdesc{i} = currcells{1};
end
fclose(fid);

genes = cell(cancer.numterms, 1);
geneids = cell(cancer.numterms, 1);
samples = cell(cancer.numterms, 1);
numgenes = zeros(cancer.numterms, 1);
numsamples = zeros(cancer.numterms, 1);
for i = 1:1:cancer.numterms
    fid = fopen(['input/ExpressionMatrices/' cancer.termdesc{i} '/datamatrix_genes_normalized.txt'], 'r');
    currline = fgetl(fid);
    currcells = strsplitbyadr(currline, '\t');
    samples{i} = currcells(3:end)';
    numsamples(i) = numel(samples{i});
    genes{i} = cell(100000, 1);
    geneids{i} = zeros(100000, 1);
    j = 0;
    while ~feof(fid)
        j = j + 1;
        currline = fgetl(fid);
        si = regexp(currline, '\t');
        genes{i}{j} = currline(1:si(1)-1);
        geneids{i}(j) = str2double(currline(si(1)+1:si(2)-1));
    end
    fclose(fid);
    if j < 100000
        genes{i}(j+1:end) = [];
        geneids{i}(j+1:end) = [];
    end
    numgenes(i) = numel(genes{i});
    gene_atb.cm(i) = cminit(numgenes(i), numsamples(i), genes{i}, 'GeneSym', [], [], geneids{i}, 'GeneID', samples{i}, 'TCGA Barcode', repmat({[cancer.term{i} ' (' cancer.termdesc{i} ')']}, numsamples(i), 1), 'Cancer Name (Cancer Acronym)', [], [], []);
    gene_atb.cm(i).matrix = dlmread(['input/ExpressionMatrices/' cancer.termdesc{i} '/datamatrix_genes_normalized.txt'], '\t', 1, 2);
end

genesmatch = false(cancer.numterms-1, 1);
for i = 2:1:cancer.numterms
    genesmatch(i-1) = all(strcmp(genes{1}, genes{i}));
end
genesmatch = all(genesmatch);

atb_gene.cm = cmtranspose(gene_atb.cm(1));
for i = 2:1:cancer.numterms
    atb_gene.cm = cmvertcat(atb_gene.cm, cmtranspose(gene_atb.cm(i)));
end
gene_atb.cm = cmtranspose(atb_gene.cm);
clear atb_gene;

sample_type = cell(gene_atb.cm.numentries, 1);
for i = 1:1:gene_atb.cm.numentries
    barcodeparts = strsplitbyadr(gene_atb.cm.entry{i}, '-');
    sample_type{i} = barcodeparts{4}(1:2);
end
[u_sample_type, ui, ri] = unique(sample_type);
num_sample_types = numel(u_sample_type);
c_sample_type = zeros(num_sample_types, 1);
for i = 1:1:num_sample_types
    c_sample_type(i) = sum(ri==i);
end
% 01 = primary solid tumor, count = 6587
% 02 = recurrent solid tumor, count = 34
% 03 = primary blood derived cancer - peripheral blood, count = 173
% 05 = additional - new primary, count = 4
% 06 = metastatic, count = 314
% 11 = solid tissue normal, count = 549

% keep only primary cancer samples (01 or 03)
discard = ~ismember(sample_type, {'01'; '03'});
gene_atb.cm = cmcoldiscard(gene_atb.cm, discard);

save('input/gene_tissuesample_tcga_primarycancersamples_unmodified_20160106.mat', '-struct', 'gene_atb');
%}

% load data
gene_atb = load('input/gene_tissuesample_tcga_primarycancersamples_unmodified_20160106.mat', '-mat', 'cm');


% map gene symbols to entrez gene symbols and discard rows corresponding to un-mapped symbols
[gene_atb.cm.term, gene_atb.cm.termname, gene_atb.cm.termid, gene_atb.cm.termidname, discard] = genesymlookup([], gene_atb.cm.termid, 'gene', 'human', false, false, mappingfilespath);
gene_atb.cm = cmrowdiscard(gene_atb.cm, discard);
clear discard;


% handle NaNs (there are none)
if sum(isnan(gene_atb.cm.matrix(:))) > 0 % FALSE
    % remove rows and columns with more than 5% missing values
    gene_atb.cm = cmnantrim_frac(gene_atb.cm, 0.95, Inf, 0.95, Inf, 'column');
end

rowhasnansorzeros = sum(isnan(gene_atb.cm.matrix) | (gene_atb.cm.matrix==0), 2) > 0;
if sum(isnan(gene_atb.cm.matrix(:))) > 0 % FALSE
    % impute remaining missing values
    gene_atb.cm = cmnanimpute(gene_atb.cm, 'median', 'row');
end


% % view column distributions
% figure(1);
% clf;
% subplot(2, 2, 1);
% hist(log2(gene_atb.cm.matrix(~rowhasnansorzeros,:)), max([10 gene_atb.cm.numterms/100]));
% subplot(2, 2, 2);
% hist(log2(gene_atb.cm.matrix(~rowhasnansorzeros,randperm(gene_atb.cm.numentries,4))), max([10 gene_atb.cm.numterms/100]));
% subplot(2, 2, 3);
% hist(log2(gene_atb.cm.matrix(~rowhasnansorzeros,randperm(gene_atb.cm.numentries,4))), max([10 gene_atb.cm.numterms/100]));
% subplot(2, 2, 4);
% hist(log2(gene_atb.cm.matrix(~rowhasnansorzeros,randperm(gene_atb.cm.numentries,4))), max([10 gene_atb.cm.numterms/100]));


% quantify extent to which signatures for same condition have high similarity
atb_gene.sm = cm2sm_cosine_nocluster(cmtranspose(gene_atb.cm));
subject_id = cell(atb_gene.sm.numterms, 1);
for i = 1:1:atb_gene.sm.numterms
    di = find(atb_gene.sm.term{i} == '-');
    subject_id{i} = atb_gene.sm.term{i}(1:di(3)-1);
end
[i, j] = find(tril(true(size(atb_gene.sm.matrix)),-1));
k = sub2ind(size(atb_gene.sm.matrix), i, j);
X = atb_gene.sm.matrix(k);
Y = double(strcmp(atb_gene.sm.termdesc(i), atb_gene.sm.termdesc(j)));
discard = strcmp(subject_id(i), subject_id(j));
Y(discard) = [];
X(discard) = [];
[stats_genesmapped, fpr, tpr] = classifierperformance(Y, X, [], []);
fpr = [0; fpr; 1];
tpr = [0; tpr; 1];
[fpr, ui] = unique(fpr);
tpr = tpr(ui);
fpr_x = (0:0.01:1)';
tpr_genesmapped = interp1(fpr, tpr, fpr_x);
figure(11);
clf;
subplot(2,3,1);
plot(fpr_x, tpr_genesmapped, '-k');
title('genesmapped');
axis([0 1 0 1]);
axis square;
clear atb_gene subject_id i j k discard fpr tpr ui X Y;


% % find lower detection limit
% detlim = zeros([1000 gene_atb.cm.numentries]);
% posmin = zeros([1 gene_atb.cm.numentries]);
% pos25 = zeros([1 gene_atb.cm.numentries]);
% pos50 = zeros([1 gene_atb.cm.numentries]);
% pos75 = zeros([1 gene_atb.cm.numentries]);
% posmax = zeros([1 gene_atb.cm.numentries]);
% 
% for i = 1:1:gene_atb.cm.numentries
%     x = gene_atb.cm.matrix(gene_atb.cm.matrix(:,i) > 0,i);
%     x = sort(x);
%     detlim(:,i) = x(1:1000);
%     posmin(i) = x(1);
%     pos25(i) = x(round(0.25*numel(x)));
%     pos50(i) = x(round(0.50*numel(x)));
%     pos75(i) = x(round(0.75*numel(x)));
%     posmax(i) = x(end);
% end
% 
% [pos50, si] = sort(pos50);
% pos25 = pos25(si);
% posmin = posmin(si);
% pos75 = pos75(si);
% posmax = posmax(si);
% 
% figure(2);
% clf;
% subplot(1,2,1);
% hist(log2(detlim), 100);
% subplot(1,2,2);
% plot(1:1:gene_atb.cm.numentries, log2(posmin), '-ok', 1:1:gene_atb.cm.numentries, log2(pos25), '-ok', 1:1:gene_atb.cm.numentries, log2(pos50), '-ok', 1:1:gene_atb.cm.numentries, log2(pos75), '-ok', 1:1:gene_atb.cm.numentries, log2(posmax), '-ok');


% handle zeros (expression level below detection)
gene_atb.cm.matrix(gene_atb.cm.matrix == 0) = NaN;
if sum(isnan(gene_atb.cm.matrix(:))) > 0 % TRUE
    % remove rows and columns with more than 5% undetected values
    frac = (0.05:0.05:0.95)';
    for i = 1:1:numel(frac)
        gene_atb.cm = cmnantrim_frac(gene_atb.cm, frac(i), Inf, frac(i), Inf, 'column');
    end
end

rowhaszeros = sum(isnan(gene_atb.cm.matrix), 2) > 0;
if sum(isnan(gene_atb.cm.matrix(:))) > 0 % TRUE
    % impute remaining undetected values
%     gene_atb.cm.matrix(isnan(gene_atb.cm.matrix)) = mean(posmin)/2;
    gene_atb.cm = cmnanimpute(gene_atb.cm, 'median', 'row');
end
clear detlim frac pos25 pos50 pos75 posmax posmin si x;


% % view column distributions
% figure(3);
% clf;
% subplot(2, 2, 1);
% hist(log2(gene_atb.cm.matrix), max([10 gene_atb.cm.numterms/100]));
% subplot(2, 2, 2);
% hist(log2(gene_atb.cm.matrix(:,randperm(gene_atb.cm.numentries,4))), max([10 gene_atb.cm.numterms/100]));
% subplot(2, 2, 3);
% hist(log2(gene_atb.cm.matrix(:,randperm(gene_atb.cm.numentries,4))), max([10 gene_atb.cm.numterms/100]));
% subplot(2, 2, 4);
% hist(log2(gene_atb.cm.matrix(:,randperm(gene_atb.cm.numentries,4))), max([10 gene_atb.cm.numterms/100]));


% quantify extent to which signatures for same condition have high similarity
atb_gene.sm = cm2sm_cosine_nocluster(cmtranspose(gene_atb.cm));
subject_id = cell(atb_gene.sm.numterms, 1);
for i = 1:1:atb_gene.sm.numterms
    di = find(atb_gene.sm.term{i} == '-');
    subject_id{i} = atb_gene.sm.term{i}(1:di(3)-1);
end
[i, j] = find(tril(true(size(atb_gene.sm.matrix)),-1));
k = sub2ind(size(atb_gene.sm.matrix), i, j);
X = atb_gene.sm.matrix(k);
Y = double(strcmp(atb_gene.sm.termdesc(i), atb_gene.sm.termdesc(j)));
discard = strcmp(subject_id(i), subject_id(j));
Y(discard) = [];
X(discard) = [];
[stats_zeroshandled, fpr, tpr] = classifierperformance(Y, X, [], []);
fpr = [0; fpr; 1];
tpr = [0; tpr; 1];
[fpr, ui] = unique(fpr);
tpr = tpr(ui);
fpr_x = (0:0.01:1)';
tpr_zeroshandled = interp1(fpr, tpr, fpr_x);
figure(11);
subplot(2,3,2);
plot(fpr_x, tpr_zeroshandled, '-k');
title('zeroshandled');
axis([0 1 0 1]);
axis square;
clear atb_gene subject_id i j k discard fpr tpr ui X Y;


% log2 transform
gene_atb.cm.matrix = log2(gene_atb.cm.matrix);


% % view column distributions
% figure(4);
% clf;
% subplot(2, 2, 1);
% hist(gene_atb.cm.matrix, max([10 gene_atb.cm.numterms/100]));
% subplot(2, 2, 2);
% hist(gene_atb.cm.matrix(:,randperm(gene_atb.cm.numentries,4)), max([10 gene_atb.cm.numterms/100]));
% subplot(2, 2, 3);
% hist(gene_atb.cm.matrix(:,randperm(gene_atb.cm.numentries,4)), max([10 gene_atb.cm.numterms/100]));
% subplot(2, 2, 4);
% hist(gene_atb.cm.matrix(:,randperm(gene_atb.cm.numentries,4)), max([10 gene_atb.cm.numterms/100]));


% quantify extent to which signatures for same condition have high similarity
atb_gene.sm = cm2sm_cosine_nocluster(cmtranspose(gene_atb.cm));
subject_id = cell(atb_gene.sm.numterms, 1);
for i = 1:1:atb_gene.sm.numterms
    di = find(atb_gene.sm.term{i} == '-');
    subject_id{i} = atb_gene.sm.term{i}(1:di(3)-1);
end
[i, j] = find(tril(true(size(atb_gene.sm.matrix)),-1));
k = sub2ind(size(atb_gene.sm.matrix), i, j);
X = atb_gene.sm.matrix(k);
Y = double(strcmp(atb_gene.sm.termdesc(i), atb_gene.sm.termdesc(j)));
discard = strcmp(subject_id(i), subject_id(j));
Y(discard) = [];
X(discard) = [];
[stats_log2, fpr, tpr] = classifierperformance(Y, X, [], []);
fpr = [0; fpr; 1];
tpr = [0; tpr; 1];
[fpr, ui] = unique(fpr);
tpr = tpr(ui);
fpr_x = (0:0.01:1)';
tpr_log2 = interp1(fpr, tpr, fpr_x);
figure(11);
subplot(2,3,3);
plot(fpr_x, tpr_log2, '-k');
title('log2');
axis([0 1 0 1]);
axis square;
clear atb_gene subject_id i j k discard fpr tpr ui X Y;


% quantile normalization
gene_atb.cm.matrix = quantilenormalization(gene_atb.cm.matrix);


% % view column distributions
% figure(5);
% clf;
% subplot(2, 2, 1);
% hist(gene_atb.cm.matrix, max([10 gene_atb.cm.numterms/100]));
% subplot(2, 2, 2);
% hist(gene_atb.cm.matrix(:,randperm(gene_atb.cm.numentries,4)), max([10 gene_atb.cm.numterms/100]));
% subplot(2, 2, 3);
% hist(gene_atb.cm.matrix(:,randperm(gene_atb.cm.numentries,4)), max([10 gene_atb.cm.numterms/100]));
% subplot(2, 2, 4);
% hist(gene_atb.cm.matrix(:,randperm(gene_atb.cm.numentries,4)), max([10 gene_atb.cm.numterms/100]));


% quantify extent to which signatures for same condition have high similarity
atb_gene.sm = cm2sm_cosine_nocluster(cmtranspose(gene_atb.cm));
subject_id = cell(atb_gene.sm.numterms, 1);
for i = 1:1:atb_gene.sm.numterms
    di = find(atb_gene.sm.term{i} == '-');
    subject_id{i} = atb_gene.sm.term{i}(1:di(3)-1);
end
[i, j] = find(tril(true(size(atb_gene.sm.matrix)),-1));
k = sub2ind(size(atb_gene.sm.matrix), i, j);
X = atb_gene.sm.matrix(k);
Y = double(strcmp(atb_gene.sm.termdesc(i), atb_gene.sm.termdesc(j)));
discard = strcmp(subject_id(i), subject_id(j));
Y(discard) = [];
X(discard) = [];
[stats_qn, fpr, tpr] = classifierperformance(Y, X, [], []);
fpr = [0; fpr; 1];
tpr = [0; tpr; 1];
[fpr, ui] = unique(fpr);
tpr = tpr(ui);
fpr_x = (0:0.01:1)';
tpr_qn = interp1(fpr, tpr, fpr_x);
figure(11);
subplot(2,3,4);
plot(fpr_x, tpr_qn, '-k');
title('qn');
axis([0 1 0 1]);
axis square;
clear atb_gene subject_id i j k discard fpr tpr ui X Y;


% merge measurements corresponding to the same gene
if numel(unique(gene_atb.cm.term)) < gene_atb.cm.numterms % FALSE
    gene_atb.cm = cmrowmerge(gene_atb.cm, 'mean');
end


% % view column distributions
% figure(6);
% clf;
% subplot(2, 2, 1);
% hist(gene_atb.cm.matrix, max([10 gene_atb.cm.numterms/100]));
% subplot(2, 2, 2);
% hist(gene_atb.cm.matrix(:,randperm(gene_atb.cm.numentries,4)), max([10 gene_atb.cm.numterms/100]));
% subplot(2, 2, 3);
% hist(gene_atb.cm.matrix(:,randperm(gene_atb.cm.numentries,4)), max([10 gene_atb.cm.numterms/100]));
% subplot(2, 2, 4);
% hist(gene_atb.cm.matrix(:,randperm(gene_atb.cm.numentries,4)), max([10 gene_atb.cm.numterms/100]));


% % quantify extent to which signatures for same condition have high similarity
% atb_gene.sm = cm2sm_cosine_nocluster(cmtranspose(gene_atb.cm));
% subject_id = cell(atb_gene.sm.numterms, 1);
% for i = 1:1:atb_gene.sm.numterms
%     di = find(atb_gene.sm.term{i} == '-');
%     subject_id{i} = atb_gene.sm.term{i}(1:di(3)-1);
% end
% [i, j] = find(tril(true(size(atb_gene.sm.matrix)),-1));
% k = sub2ind(size(atb_gene.sm.matrix), i, j);
% X = atb_gene.sm.matrix(k);
% Y = double(strcmp(atb_gene.sm.termdesc(i), atb_gene.sm.termdesc(j)));
% discard = strcmp(subject_id(i), subject_id(j));
% Y(discard) = [];
% X(discard) = [];
% [stats_mergegenes, fpr, tpr] = classifierperformance(Y, X, [], []);
% fpr = [0; fpr; 1];
% tpr = [0; tpr; 1];
% [fpr, ui] = unique(fpr);
% tpr = tpr(ui);
% fpr_x = (0:0.01:1)';
% tpr_mergegenes = interp1(fpr, tpr, fpr_x);
% figure(11);
% subplot(2,3,5);
% plot(fpr_x, tpr_mergegenes, '-k');
% title('mergegenes');
% axis([0 1 0 1]);
% axis square;
% clear atb_gene subject_id i j k discard fpr tpr ui X Y;


% cluster and view matrix
gene_atb.cm = cmcluster(gene_atb.cm, false);


% save result
save('output/gene_attribute_matrix_imported.mat', '-struct', 'gene_atb');




clearvars -except gene_atb;

% normalize, prefer KS since it is non-parametric and therefore doesn't
% assume anything about the shape of the distribution of values per row
threshfrac = 0.05; % not relevant for normalization only
type = 'tertiary';
method = 'rows';
discardemptyvectors = false;
[~, gene_atb.cm] = cmthresh_ks(gene_atb.cm, threshfrac, type, method, discardemptyvectors);


% standardize and threshold
threshfrac = 0.1;
type = 'tertiary';
method = 'matrix';
discardemptyvectors = true;
[tm, gene_atb.cm] = cmthresh_ks(gene_atb.cm, threshfrac, type, method, discardemptyvectors);


% cluster and view matrix
gene_atb.cm = cmcluster(gene_atb.cm, false);


% quantify extent to which signatures for same condition have high similarity
atb_gene.sm = cm2sm_cosine_nocluster(cmtranspose(gene_atb.cm));
subject_id = cell(atb_gene.sm.numterms, 1);
for i = 1:1:atb_gene.sm.numterms
    di = find(atb_gene.sm.term{i} == '-');
    subject_id{i} = atb_gene.sm.term{i}(1:di(3)-1);
end
[i, j] = find(tril(true(size(atb_gene.sm.matrix)),-1));
k = sub2ind(size(atb_gene.sm.matrix), i, j);
X = atb_gene.sm.matrix(k);
Y = double(strcmp(atb_gene.sm.termdesc(i), atb_gene.sm.termdesc(j)));
discard = strcmp(subject_id(i), subject_id(j));
Y(discard) = [];
X(discard) = [];
[stats_std, fpr, tpr] = classifierperformance(Y, X, [], []);
fpr = [0; fpr; 1];
tpr = [0; tpr; 1];
[fpr, ui] = unique(fpr);
tpr = tpr(ui);
fpr_x = (0:0.01:1)';
tpr_std = interp1(fpr, tpr, fpr_x);
figure(12);
clf;
subplot(1,2,1);
plot(fpr_x, tpr_std, '-k');
title('std');
axis([0 1 0 1]);
axis square;
clear atb_gene subject_id i j k discard fpr tpr ui X Y;


% save standardized matrix
save('output/gene_attribute_matrix_standardized.mat', '-struct', 'gene_atb');


% cluster thresholded matrix
gene_atb.cm = conmatmap(tm, gene_atb.cm.term, gene_atb.cm.entry);
HeatMap(gene_atb.cm.matrix, 'colormap', redbluecmap);
clear tm;


% quantify extent to which signatures for same condition have high similarity
atb_gene.sm = cm2sm_cosine_nocluster(cmtranspose(gene_atb.cm));
subject_id = cell(atb_gene.sm.numterms, 1);
for i = 1:1:atb_gene.sm.numterms
    di = find(atb_gene.sm.term{i} == '-');
    subject_id{i} = atb_gene.sm.term{i}(1:di(3)-1);
end
[i, j] = find(tril(true(size(atb_gene.sm.matrix)),-1));
k = sub2ind(size(atb_gene.sm.matrix), i, j);
X = atb_gene.sm.matrix(k);
Y = double(strcmp(atb_gene.sm.termdesc(i), atb_gene.sm.termdesc(j)));
discard = strcmp(subject_id(i), subject_id(j));
Y(discard) = [];
X(discard) = [];
[stats_thr, fpr, tpr] = classifierperformance(Y, X, [], []);
fpr = [0; fpr; 1];
tpr = [0; tpr; 1];
[fpr, ui] = unique(fpr);
tpr = tpr(ui);
fpr_x = (0:0.01:1)';
tpr_thr = interp1(fpr, tpr, fpr_x);
figure(12);
subplot(1,2,2);
plot(fpr_x, tpr_thr, '-k');
title('thr');
axis([0 1 0 1]);
axis square;
clear atb_gene subject_id i j k discard fpr tpr ui X Y;


% convert to sparse matrix
gene_atb.cm.matrix = sparse(gene_atb.cm.matrix);


% save thresholded matrix
save('output/gene_attribute_matrix_prepared.mat', '-struct', 'gene_atb');


