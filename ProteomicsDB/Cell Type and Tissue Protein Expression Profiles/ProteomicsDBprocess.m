


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
% initialize edges structure
gene_atb1.edges = edgesinit(2261302 + 7374, [], 'GeneSym', [], 'UniprotAcc', [], [], [], 'Tissue/Cell Line', [], 'Tissue', [], 'Brenda Tissue Ontology BTO:', true, []);



% read file part 1
fid = fopen('input/protein_expression_mslev1_scpsel1_calc0_part1_20150129.txt', 'r');

currline = fgetl(fid);

for i = 1:1:gene_atb1.edges.numedges
    
    currline = fgetl(fid);
    
    currcells = strsplit(currline, '\t', 'CollapseDelimiters', false);
    
    gene_atb1.edges.source{i} = currcells{1};
    
    gene_atb1.edges.sourcedesc{i} = currcells{3};
    
    gene_atb1.edges.target{i} = currcells{4};
    
    gene_atb1.edges.targetdesc{i} = currcells{5};
    
    gene_atb1.edges.targetid(i) = str2double(currcells{6}(5:end));
    
    gene_atb1.edges.weight(i) = str2double(currcells{7});
    
end

fclose(fid);

% read file part 2
fid = fopen('input/protein_expression_mslev1_scpsel1_calc0_part2_20150129.txt', 'r');

currline = fgetl(fid);

for j = i+1:1:gene_atb1.edges.numedges
    
    currline = fgetl(fid);
    
    currcells = strsplit(currline, '\t', 'CollapseDelimiters', false);
    
    gene_atb1.edges.source{j} = currcells{1};
    
    gene_atb1.edges.sourcedesc{j} = currcells{3};
    
    gene_atb1.edges.target{j} = currcells{4};
    
    gene_atb1.edges.targetdesc{j} = currcells{5};
    
    gene_atb1.edges.targetid(j) = str2double(currcells{6}(5:end));
    
    gene_atb1.edges.weight(j) = str2double(currcells{7});
    
end

fclose(fid);

% discard bad values
discard = gene_atb1.edges.weight <= 0 | isnan(gene_atb1.edges.weight); % 14 negative values, no NaN values
gene_atb1.edges = edgesdiscard(gene_atb1.edges, discard);

% get mean expression for each gene-tissue pair
gene_atb1.uedgesmean = edgesunique_mean(gene_atb1.edges);

% save imported data
save('input/protein_expression_mslev1_scpsel1_calc0_20150129.mat', '-struct', 'gene_atb1');



% initialize edges structure
gene_atb2.edges = edgesinit(437875, [], 'GeneSym', [], 'UniprotAcc', [], [], [], 'Tissue/Cell Line', [], 'Tissue', [], 'Brenda Tissue Ontology BTO:', true, []);

% read file
fid = fopen('input/protein_expression_mslev1_scpsel2_calc0_20150129.txt', 'r');

currline = fgetl(fid);

for i = 1:1:gene_atb2.edges.numedges
    
    currline = fgetl(fid);
    
    currcells = strsplit(currline, '\t', 'CollapseDelimiters', false);
    
    gene_atb2.edges.source{i} = currcells{1};
    
    gene_atb2.edges.sourcedesc{i} = currcells{3};
    
    gene_atb2.edges.target{i} = currcells{4};
    
    gene_atb2.edges.targetdesc{i} = currcells{5};
    
    gene_atb2.edges.targetid(i) = str2double(currcells{6}(5:end));
    
    gene_atb2.edges.weight(i) = str2double(currcells{7});
    
end

fclose(fid);

% discard bad values
discard = gene_atb2.edges.weight <= 0 | isnan(gene_atb2.edges.weight); % 0 negative values, no NaN values
gene_atb2.edges = edgesdiscard(gene_atb2.edges, discard);

% get mean expression for each gene-tissue pair
gene_atb2.uedgesmean = edgesunique_mean(gene_atb2.edges);

% save imported data
save('input/protein_expression_mslev1_scpsel2_calc0_20150129.mat', '-struct', 'gene_atb2');



% initialize edges structure
gene_atb3.edges = edgesinit(336300, [], 'GeneSym', [], 'UniprotAcc', [], [], [], 'Tissue/Cell Line', [], 'Tissue', [], 'Brenda Tissue Ontology BTO:', true, []);

% read file
fid = fopen('input/protein_expression_mslev2_scpsel1_calc0_20150129.txt', 'r');

currline = fgetl(fid);

for i = 1:1:gene_atb3.edges.numedges
    
    currline = fgetl(fid);
    
    currcells = strsplit(currline, '\t', 'CollapseDelimiters', false);
    
    gene_atb3.edges.source{i} = currcells{1};
    
    gene_atb3.edges.sourcedesc{i} = currcells{3};
    
    gene_atb3.edges.target{i} = currcells{4};
    
    gene_atb3.edges.targetdesc{i} = currcells{5};
    
    gene_atb3.edges.targetid(i) = str2double(currcells{6}(5:end));
    
    gene_atb3.edges.weight(i) = str2double(currcells{7});
    
end

fclose(fid);

% discard bad values
discard = gene_atb3.edges.weight <= 0 | isnan(gene_atb3.edges.weight); % 1 negative values, no NaN values
gene_atb3.edges = edgesdiscard(gene_atb3.edges, discard);

% get mean expression for each gene-tissue pair
gene_atb3.uedgesmean = edgesunique_mean(gene_atb3.edges);

% save imported data
save('input/protein_expression_mslev2_scpsel1_calc0_20150129.mat', '-struct', 'gene_atb3');



% initialize edges structure
gene_atb4.edges = edgesinit(181475, [], 'GeneSym', [], 'UniprotAcc', [], [], [], 'Tissue/Cell Line', [], 'Tissue', [], 'Brenda Tissue Ontology BTO:', true, []);

% read file
fid = fopen('input/protein_expression_mslev2_scpsel2_calc0_20150129.txt', 'r');

currline = fgetl(fid);

for i = 1:1:gene_atb4.edges.numedges
    
    currline = fgetl(fid);
    
    currcells = strsplit(currline, '\t', 'CollapseDelimiters', false);
    
    gene_atb4.edges.source{i} = currcells{1};
    
    gene_atb4.edges.sourcedesc{i} = currcells{3};
    
    gene_atb4.edges.target{i} = currcells{4};
    
    gene_atb4.edges.targetdesc{i} = currcells{5};
    
    gene_atb4.edges.targetid(i) = str2double(currcells{6}(5:end));
    
    gene_atb4.edges.weight(i) = str2double(currcells{7});
    
end

fclose(fid);

% discard bad values
discard = gene_atb4.edges.weight <= 0 | isnan(gene_atb4.edges.weight); % no negative values, no NaN values
gene_atb4.edges = edgesdiscard(gene_atb4.edges, discard);

% get mean expression for each gene-tissue pair
gene_atb4.uedgesmean = edgesunique_mean(gene_atb4.edges);

% save imported data
save('input/protein_expression_mslev2_scpsel2_calc0_20150129.mat', '-struct', 'gene_atb4');



% diagnostics
figure(1);
clf;
subplot(2,2,1);
hist(gene_atb1.uedgesmean.weight, 100);
set(gca, 'PlotBoxAspectRatio', [1 1 1]);
subplot(2,2,2);
hist(gene_atb2.uedgesmean.weight, 100);
set(gca, 'PlotBoxAspectRatio', [1 1 1]);
subplot(2,2,3);
hist(gene_atb3.uedgesmean.weight, 100);
set(gca, 'PlotBoxAspectRatio', [1 1 1]);
subplot(2,2,4);
hist(gene_atb4.uedgesmean.weight, 100);
set(gca, 'PlotBoxAspectRatio', [1 1 1]);

e1 = sourcetargetcat(gene_atb1.uedgesmean.source, gene_atb1.uedgesmean.target);
e2 = sourcetargetcat(gene_atb2.uedgesmean.source, gene_atb2.uedgesmean.target);
e3 = sourcetargetcat(gene_atb3.uedgesmean.source, gene_atb3.uedgesmean.target);
e4 = sourcetargetcat(gene_atb4.uedgesmean.source, gene_atb4.uedgesmean.target);
e12 = intersect(e1, e2);
e13 = intersect(e1, e3);
e14 = intersect(e1, e4);
e23 = intersect(e2, e3);
e24 = intersect(e2, e4);
e34 = intersect(e3, e4);
n12 = numel(e12);
n13 = numel(e13);
n14 = numel(e14);
n23 = numel(e23);
n24 = numel(e24);
n34 = numel(e34);
v12 = zeros([n12 2]);
v13 = zeros([n13 2]);
v14 = zeros([n14 2]);
v23 = zeros([n23 2]);
v24 = zeros([n24 2]);
v34 = zeros([n34 2]);

[o1, o2] = ismember(e12, e1);
o2(o2 == 0) = [];
v12(o1,1) = gene_atb1.uedgesmean.weight(o2);
[o1, o2] = ismember(e12, e2);
o2(o2 == 0) = [];
v12(o1,2) = gene_atb2.uedgesmean.weight(o2);

[o1, o2] = ismember(e13, e1);
o2(o2 == 0) = [];
v13(o1,1) = gene_atb1.uedgesmean.weight(o2);
[o1, o2] = ismember(e13, e3);
o2(o2 == 0) = [];
v13(o1,2) = gene_atb3.uedgesmean.weight(o2);

[o1, o2] = ismember(e14, e1);
o2(o2 == 0) = [];
v14(o1,1) = gene_atb1.uedgesmean.weight(o2);
[o1, o2] = ismember(e14, e4);
o2(o2 == 0) = [];
v14(o1,2) = gene_atb4.uedgesmean.weight(o2);

[o1, o2] = ismember(e23, e2);
o2(o2 == 0) = [];
v23(o1,1) = gene_atb2.uedgesmean.weight(o2);
[o1, o2] = ismember(e23, e3);
o2(o2 == 0) = [];
v23(o1,2) = gene_atb3.uedgesmean.weight(o2);

[o1, o2] = ismember(e24, e2);
o2(o2 == 0) = [];
v24(o1,1) = gene_atb2.uedgesmean.weight(o2);
[o1, o2] = ismember(e24, e4);
o2(o2 == 0) = [];
v24(o1,2) = gene_atb4.uedgesmean.weight(o2);

[o1, o2] = ismember(e34, e3);
o2(o2 == 0) = [];
v34(o1,1) = gene_atb3.uedgesmean.weight(o2);
[o1, o2] = ismember(e34, e4);
o2(o2 == 0) = [];
v34(o1,2) = gene_atb4.uedgesmean.weight(o2);

figure(2);
clf;
subplot(2,3,1);
plot(v12(:,1), v12(:,2), '.k', [0 10], [0 10], '-r');
axis([0 10 0 10]);
xlabel('1');
ylabel('2');
set(gca, 'PlotBoxAspectRatio', [1 1 1]);
subplot(2,3,2);
plot(v13(:,1), v13(:,2), '.k', [0 10], [0 10], '-r');
axis([0 10 0 10]);
xlabel('1');
ylabel('3');
set(gca, 'PlotBoxAspectRatio', [1 1 1]);
subplot(2,3,3);
plot(v14(:,1), v14(:,2), '.k', [0 10], [0 10], '-r');
axis([0 10 0 10]);
xlabel('1');
ylabel('4');
set(gca, 'PlotBoxAspectRatio', [1 1 1]);
subplot(2,3,4);
plot(v23(:,1), v23(:,2), '.k', [0 10], [0 10], '-r');
axis([0 10 0 10]);
xlabel('2');
ylabel('3');
set(gca, 'PlotBoxAspectRatio', [1 1 1]);
subplot(2,3,5);
plot(v24(:,1), v24(:,2), '.k', [0 10], [0 10], '-r');
axis([0 10 0 10]);
xlabel('2');
ylabel('4');
set(gca, 'PlotBoxAspectRatio', [1 1 1]);
subplot(2,3,6);
plot(v34(:,1), v34(:,2), '.k', [0 10], [0 10], '-r');
axis([0 10 0 10]);
xlabel('3');
ylabel('4');
set(gca, 'PlotBoxAspectRatio', [1 1 1]);



% transform values (I don't think this is necessary)
%{
b12 = robustfit(v12(:,2), v12(:,1));
v12(:,3) = b12(1) + b12(2)*v12(:,2);
% gene_atb2.edges.weight = b12(1) + b12(2)*gene_atb2.edges.weight;

b13 = robustfit(v13(:,2), v13(:,1));
v13(:,3) = b13(1) + b13(2)*v13(:,2);
% gene_atb3.edges.weight = b13(1) + b13(2)*gene_atb3.edges.weight;

b14 = robustfit(v14(:,2), v14(:,1));
v14(:,3) = b14(1) + b14(2)*v14(:,2);
% gene_atb4.edges.weight = b14(1) + b14(2)*gene_atb4.edges.weight;

figure(3);
clf;
subplot(1,3,1);
plot(v12(:,2), v12(:,1), '.k', [0 10], [0 10], '-r', v12(:,2), v12(:,3), '-b');
axis([0 10 0 10]);
xlabel('2');
ylabel('1');
set(gca, 'PlotBoxAspectRatio', [1 1 1]);
subplot(1,3,2);
plot(v13(:,2), v13(:,1), '.k', [0 10], [0 10], '-r', v13(:,2), v13(:,3), '-b');
axis([0 10 0 10]);
xlabel('3');
ylabel('1');
set(gca, 'PlotBoxAspectRatio', [1 1 1]);
subplot(1,3,3);
plot(v14(:,2), v14(:,1), '.k', [0 10], [0 10], '-r', v14(:,2), v14(:,3), '-b');
axis([0 10 0 10]);
xlabel('4');
ylabel('1');
set(gca, 'PlotBoxAspectRatio', [1 1 1]);

figure(4);
clf;
subplot(2,3,1);
plot(v12(:,1), v12(:,3), '.k', [0 10], [0 10], '-r');
axis([0 10 0 10]);
xlabel('1');
ylabel('2');
set(gca, 'PlotBoxAspectRatio', [1 1 1]);
subplot(2,3,2);
plot(v13(:,1), v13(:,3), '.k', [0 10], [0 10], '-r');
axis([0 10 0 10]);
xlabel('1');
ylabel('3');
set(gca, 'PlotBoxAspectRatio', [1 1 1]);
subplot(2,3,3);
plot(v14(:,1), v14(:,3), '.k', [0 10], [0 10], '-r');
axis([0 10 0 10]);
xlabel('1');
ylabel('4');
set(gca, 'PlotBoxAspectRatio', [1 1 1]);
subplot(2,3,4);
plot(v23(:,1), v23(:,3), '.k', [0 10], [0 10], '-r');
axis([0 10 0 10]);
xlabel('2');
ylabel('3');
set(gca, 'PlotBoxAspectRatio', [1 1 1]);
subplot(2,3,5);
plot(v24(:,1), v24(:,3), '.k', [0 10], [0 10], '-r');
axis([0 10 0 10]);
xlabel('2');
ylabel('4');
set(gca, 'PlotBoxAspectRatio', [1 1 1]);
subplot(2,3,6);
plot(v34(:,1), v34(:,3), '.k', [0 10], [0 10], '-r');
axis([0 10 0 10]);
xlabel('3');
ylabel('4');
set(gca, 'PlotBoxAspectRatio', [1 1 1]);
%}



% average all data
clearvars -except gene_atb1 gene_atb2 gene_atb3 gene_atb4;

gene_atb.edges = edgesvertcat(gene_atb4.edges, gene_atb3.edges);
clear gene_atb4 gene_atb3;
gene_atb.edges = edgesvertcat(gene_atb.edges, gene_atb2.edges);
clear gene_atb2;
gene_atb.edges = edgesvertcat(gene_atb.edges, gene_atb1.edges);
clear gene_atb1;

gene_atb.edges = edgesunique_mean(gene_atb.edges);



% convert to connectivity matrix (all zeros are "no data from database",
% but if no data was in database, does that mean "not probed" or "probed
% and not detected"?)
gene_atb.cm = edges2cm(gene_atb.edges);
gene_atb = rmfield(gene_atb, 'edges');



% map gene symbols to entrez gene symbols and discard rows corresponding to
% un-mapped symbols ( rows discarded)
[gene_atb.cm.term, gene_atb.cm.termname, gene_atb.cm.termid, gene_atb.cm.termidname, discard] = genesymlookup(gene_atb.cm.term, [], 'gene', 'human', true, false, mappingfilespath);
gene_atb.cm = cmrowdiscard(gene_atb.cm, discard);



save('input/gene_tissue_proteomicsdb_unprocessed_20150129.mat', '-struct', 'gene_atb');
%}
gene_atb = load('input/gene_tissue_proteomicsdb_unprocessed_20150129.mat', '-mat', 'cm');



gene_atb.cm.matrix(gene_atb.cm.matrix == 0) = NaN;
if sum(isnan(gene_atb.cm.matrix(:))) > 0 % TRUE
    % remove rows and columns with at more than 5% missing values (14940
    % rows discarded, 157 columns discarded. 2801 x 62 remaining)
    l = (0.05:0.01:0.95)';
    p1 = 1.2;
    s1 = 0.95;
    p2 = 1/1.2;
    s2 = 0.95;
    c1 = l.^p1/l(end)^p1*s1;
    c2 = l.^p2/l(end)^p2*s2;
    for i = 1:1:numel(c1)
        gene_atb.cm = cmnantrim_frac(gene_atb.cm, c1(i), Inf, c2(i), Inf, 'column');
    end
end



if sum(isnan(gene_atb.cm.matrix(:))) > 0 % TRUE
    % impute remaining missing values
    gene_atb.cm = cmnanimpute(gene_atb.cm, 'median', 'row');
end



% view column distributions
% figure(1);
% clf;
% subplot(2, 2, 1);
% hist(gene_atb.cm.matrix, max([10 gene_atb.cm.numterms/100]));
% subplot(2, 2, 2);
% hist(gene_atb.cm.matrix(:,randperm(gene_atb.cm.numentries,4)), max([10 gene_atb.cm.numterms/100]));
% subplot(2, 2, 3);
% hist(gene_atb.cm.matrix(:,randperm(gene_atb.cm.numentries,4)), max([10 gene_atb.cm.numterms/100]));
% subplot(2, 2, 4);
% hist(gene_atb.cm.matrix(:,randperm(gene_atb.cm.numentries,4)), max([10 gene_atb.cm.numterms/100]));



% quantile normalization not needed. data already normalized.
% if ~colsarenorm(gene_atb.cm.matrix) % TRUE
%     % quantile normalize
%     gene_atb.cm.matrix = quantilenormalization(gene_atb.cm.matrix);
% end



% view column distributions
figure(2);
clf;
subplot(2, 2, 1);
hist(gene_atb.cm.matrix, max([10 gene_atb.cm.numterms/100]));
subplot(2, 2, 2);
hist(gene_atb.cm.matrix(:,randperm(gene_atb.cm.numentries,4)), max([10 gene_atb.cm.numterms/100]));
subplot(2, 2, 3);
hist(gene_atb.cm.matrix(:,randperm(gene_atb.cm.numentries,4)), max([10 gene_atb.cm.numterms/100]));
subplot(2, 2, 4);
hist(gene_atb.cm.matrix(:,randperm(gene_atb.cm.numentries,4)), max([10 gene_atb.cm.numterms/100]));



if (max(gene_atb.cm.matrix(:)) - min(gene_atb.cm.matrix(:))) > 30 % FALSE
    % log2 transformation
    gene_atb.cm.matrix = log2(gene_atb.cm.matrix);
end



if numel(unique(gene_atb.cm.term)) < gene_atb.cm.numterms % TRUE
    % merge measurements corresponding to the same gene (2779 unique gene
    % symbols out of 2802 total gene symbols)
    gene_atb.cm = cmrowmerge(gene_atb.cm, 'mean');
end



% view column distributions
figure(3);
clf;
subplot(2, 2, 1);
hist(gene_atb.cm.matrix, max([10 gene_atb.cm.numterms/100]));
subplot(2, 2, 2);
hist(gene_atb.cm.matrix(:,randperm(gene_atb.cm.numentries,4)), max([10 gene_atb.cm.numterms/100]));
subplot(2, 2, 3);
hist(gene_atb.cm.matrix(:,randperm(gene_atb.cm.numentries,4)), max([10 gene_atb.cm.numterms/100]));
subplot(2, 2, 4);
hist(gene_atb.cm.matrix(:,randperm(gene_atb.cm.numentries,4)), max([10 gene_atb.cm.numterms/100]));



sm = cm2sm_cosine_nocluster(gene_atb.cm);
numidentical = sum(sm.matrix > (1 - 1e-12), 2) - 1;
suspectgenes = sm.term(numidentical > 0);
clear sm;
if ~isempty(suspectgenes) % FALSE
    % remove suspect rows that have identical values to at least one other row
    discard = ismember(gene_atb.cm.term, suspectgenes);
    gene_atb.cm = cmrowdiscard(gene_atb.cm, discard);
end



sm = cm2sm_cosine_nocluster(cmtranspose(gene_atb.cm));
numidentical = sum(sm.matrix > (1 - 1e-12), 2) - 1;
suspectattributes = sm.term(numidentical > 0);
clear sm;
if ~isempty(suspectattributes) % TRUE
    % remove suspect columns that have identical values to at least one
    % other column (8 columns discarded, 54 remaining)
    discard = ismember(gene_atb.cm.entry, suspectattributes);
    gene_atb.cm = cmcoldiscard(gene_atb.cm, discard);
end



% view distributions of row and col stats
[~] = cmrowcolstats(gene_atb.cm, true, 4);



% discard rows and cols with extremely low median or extremely low median
% absolute deviation (mad) (discarded 1 row and 1 col, 2776 x 53 remaining)
threshquantile = 0.001;
gene_atb.cm = cmtrim_lowmedlowmad(gene_atb.cm, threshquantile);



% view distributions of row and col stats
[~] = cmrowcolstats(gene_atb.cm, true, 5);



% cluster and view matrix
[gene_atb.cm, ~] = cmcluster(gene_atb.cm, true);



% save result
save('output/gene_attribute_matrix_imported.mat', '-struct', 'gene_atb');






clearvars -except gene_atb;

% normalize, prefer KS since it is non-parametric and therefore doesn't
% assume anything about the shape of the distribution of values per row
threshfrac = 0.05; % not relevant for normalization only
type = 'tertiary';
method = 'rows';
discardemptyvectors = false;
[~, gene_atb.cm] = cmthresh_ks(gene_atb.cm, threshfrac, type, method, discardemptyvectors);  % Davies-Bouldin Index = 0.0831 vs 0.0183 originally
% gene_atb.cm.matrix = robustzscore(gene_atb.cm.matrix, 2);  % Davies-Bouldin Index = 0.0846 vs 0.0183 originally
% gene_atb.cm.matrix = zscore(gene_atb.cm.matrix, 0, 2);  % Davies-Bouldin Index = 0.0842 vs 0.0183 originally
% daviesbouldin = cm2daviesbouldin(cmtranspose(gene_atb.cm));



% search for best threshold (lowest davies-bouldin index)
threshfrac = [(0.01:0.01:0.04) (0.05:0.05:0.95) (0.96:0.01:0.99)]';
type = 'tertiary';
method = 'matrix';
discardemptyvectors = false;

daviesbouldin = zeros([numel(threshfrac) 1]);

for i = 1:1:numel(threshfrac)
    
    tm = cmthresh(gene_atb.cm, threshfrac(i), type, method, discardemptyvectors);

    daviesbouldin(i) = cm2daviesbouldin(cmtranspose(tm));
    
end

figure(6); clf; plot(threshfrac, daviesbouldin, '-ok'); ylabel('davies-bouldin index'); xlabel('threshold');




% threshold.  note, we could equivalently just use cmthresh to get the
% thresholded matrix, but using cmthresh_ks here in case we ever want the
% final normalized matrix (the output that is currently "~")
threshfrac = 0.15;
type = 'tertiary';
method = 'matrix';
discardemptyvectors = true;
[gene_atb.cm, ~] = cmthresh_ks(gene_atb.cm, threshfrac, type, method, discardemptyvectors);  % Davies-Bouldin Index = 0.0856
% daviesbouldin = cm2daviesbouldin(cmtranspose(gene_atb.cm));



% cluster and view matrix
[gene_atb.cm, ~] = cmcluster(gene_atb.cm, true);



% view distributions of row and col stats
[~] = cmrowcolstats(gene_atb.cm, true, 7);



% convert to sparse matrix
gene_atb.cm.matrix = sparse(gene_atb.cm.matrix);



% save result
save('output/gene_attribute_matrix_prepared.mat', '-struct', 'gene_atb');






% get standardized matrix

gene_atb = load('output/gene_attribute_matrix_imported.mat', '-mat', 'cm');

threshfrac = 0.05;
type = 'tertiary';
method = 'rows';
discardemptyvectors = false;
[~, gene_atb.cm] = cmthresh_ks(gene_atb.cm, threshfrac, type, method, discardemptyvectors);

threshfrac = 0.1;
type = 'tertiary';
method = 'matrix';
discardemptyvectors = true;
[~, gene_atb.cm] = cmthresh_ks(gene_atb.cm, threshfrac, type, method, discardemptyvectors);

% cluster and view matrix
[gene_atb.cm, ~] = cmcluster(gene_atb.cm, false);

% save result
if exist('output/gene_attribute_matrix_standardized.mat', 'file') == 0
    save('output/gene_attribute_matrix_standardized.mat', '-struct', 'gene_atb');
else
    error('file already exists?!');
end


