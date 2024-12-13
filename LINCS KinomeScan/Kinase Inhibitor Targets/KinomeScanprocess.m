


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
gene_atb.edges = edgesinit(100000, [], 'GeneSym', [], [], [], [], [], 'Chemical', [], [], [], [], true, []);


% read data
dircontents = dir('input/Files20150120');
filename = {dircontents(:).name}';
filename([dircontents(:).isdir]) = [];

i = 0;

for j = 1:1:numel(filename)
    
    fid = fopen(['input/Files20150120/' filename{j}], 'r');
    
    currline = fgetl(fid);
    
    currcells = strsplit(currline, ',', 'CollapseDelimiters', false);
    
    if strcmp(currcells{9}, 'percentControl')
    
        while ~feof(fid);

            currline = fgetl(fid);

            qo = strfind(currline, ',"');
            qc = strfind(currline, '",');

            if numel(qo) > 0 && numel(qc) > 0

                if numel(qo) == numel(qc)

                    for k = 1:1:numel(qo)

                        prestr = currline(1:qo(k)+1);
                        substr = currline(qo(k)+2:qc(k)-1);
                        poststr = currline(qc(k):end);
                        currline = [prestr strrep(substr, ',', '|') poststr];

                    end

                    currline = strrep(strrep(regexprep(currline, ',', '\t'), '|', ','), '"', '');

                    currcells = strsplit(currline, '\t', 'CollapseDelimiters', false);

                else

                    disp('problem!');

                    currline = strrep(currline, '"', '');

                    currcells = strsplit(currline, ',', 'CollapseDelimiters', false);

                end

            else

                currline = strrep(currline, '"', '');

                currcells = strsplit(currline, ',', 'CollapseDelimiters', false);

            end

            chemical = currcells{2};
            gene = currcells{7};
            weight = str2double(currcells{9});  % percent control

            if numel(chemical) > 0 && numel(gene) > 0  && numel(weight) > 0 && ~isnan(weight)

                i = i + 1;

                gene_atb.edges.source{i} = gene;

                gene_atb.edges.target{i} = chemical;

                gene_atb.edges.weight(i) = 1 - weight/100;  % convert to inhibition fraction (1 = complete inhibition)

            end

        end
    
    end
    
    fclose(fid);
    
end

if i < gene_atb.edges.numedges
    discard = false([gene_atb.edges.numedges 1]);
    discard(i+1:end) = true;
    gene_atb.edges = edgesdiscard(gene_atb.edges, discard);
end


% map gene identifiers to NCBI Entrez Gene Symbols and Gene IDs
[gene_atb.edges.source, gene_atb.edges.sourcename, gene_atb.edges.sourceid, gene_atb.edges.sourceidname, discard] = genesymlookup(gene_atb.edges.source, [], 'kinase20150306', 'both', true, true, mappingfilespath);
gene_atb.edges = edgesdiscard(gene_atb.edges, discard);


% remove duplicate edges
gene_atb.edges = edgesunique_mean(gene_atb.edges);


% convert to matrix format (attribute table)
gene_atb.cm = edges2cm_nanfill(gene_atb.edges);
gene_atb = rmfield(gene_atb, 'edges');


% view column distributions
figure(1);
clf;
subplot(2, 2, 1);
hist(gene_atb.cm.matrix, max([10 gene_atb.cm.numterms/100]));
subplot(2, 2, 2);
hist(gene_atb.cm.matrix(:,randperm(gene_atb.cm.numentries,4)), max([10 gene_atb.cm.numterms/100]));
subplot(2, 2, 3);
hist(gene_atb.cm.matrix(:,randperm(gene_atb.cm.numentries,4)), max([10 gene_atb.cm.numterms/100]));
subplot(2, 2, 4);
hist(gene_atb.cm.matrix(:,randperm(gene_atb.cm.numentries,4)), max([10 gene_atb.cm.numterms/100]));


% remove rows and columns with more than 5% missing values
if sum(isnan(gene_atb.cm.matrix(:))) > 0
    gene_atb.cm = cmnantrim_frac(gene_atb.cm, 0.6, Inf, 0.6, Inf, 'column');
    gene_atb.cm = cmnantrim_frac(gene_atb.cm, 0.8, Inf, 0.8, Inf, 'column');
    gene_atb.cm = cmnantrim_frac(gene_atb.cm, 0.9, Inf, 0.9, Inf, 'column');
    gene_atb.cm = cmnantrim_frac(gene_atb.cm, 0.95, Inf, 0.95, Inf, 'column');
end


% impute remaining missing values
if sum(isnan(gene_atb.cm.matrix(:))) > 0
    gene_atb.cm = cmnanimpute(gene_atb.cm, 'median', 'row');
end


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


% view distributions of row and col stats
[~] = cmrowcolstats(gene_atb.cm, true, 1);


% cluster and view matrix
[gene_atb.cm, ~] = cmcluster(gene_atb.cm, false);
% close force all;


% save result
save('output/gene_attribute_matrix_imported.mat', '-struct', 'gene_atb');






clearvars -except gene_atb;


% threshold values at 90% inhibition
gene_atb.cm.matrix = double(gene_atb.cm.matrix >= 0.9);
gene_atb.cm = cmtrim(gene_atb.cm, 1, Inf, 1, Inf);


% cluster and view matrix
[gene_atb.cm, ~] = cmcluster(gene_atb.cm, true);
% close force all;


% view distributions of row and col stats
[~] = cmrowcolstats(gene_atb.cm, true, 2);


% convert to sparse matrix
gene_atb.cm.matrix = sparse(gene_atb.cm.matrix);


% save result
save('output/gene_attribute_matrix_prepared.mat', '-struct', 'gene_atb');







% get cleaned matrix

gene_atb = load('output/gene_attribute_matrix_imported.mat', '-mat');
gene_atb.cm.matrix = full(gene_atb.cm.matrix);

if sum(sum(gene_atb.cm.matrix == 0)) > 1/3*numel(gene_atb.cm.matrix)
    gene_atb.cm.matrix = sparse(gene_atb.cm.matrix);
end
save('output/gene_attribute_matrix_cleaned.mat', '-struct', 'gene_atb');
if issparse(gene_atb.cm.matrix)
    gene_atb.cm.matrix = full(gene_atb.cm.matrix);
end
savepath = 'output/';
mkdir(savepath);
writecm([savepath 'gene_attribute_matrix_cleaned'], gene_atb.cm);
gzip([savepath 'gene_attribute_matrix_cleaned.txt']);
delete([savepath 'gene_attribute_matrix_cleaned.txt']);






% get standardized matrix
clearvars -except savepath gene_atb;

% threshfrac = 0.1;
% type = 'binary';
% method = 'matrix';
% discardemptyvectors = true;
% [~, gene_atb.cm] = cmthresh_ks(gene_atb.cm, threshfrac, type, method, discardemptyvectors);
gene_atb.cm = cm_standardize_ignorezeros(gene_atb.cm);

% cluster and view matrix
[gene_atb.cm, ~] = cmcluster(gene_atb.cm, false);

% save result
if sum(sum(gene_atb.cm.matrix == 0)) > 1/3*numel(gene_atb.cm.matrix)
    gene_atb.cm.matrix = sparse(gene_atb.cm.matrix);
end
save('output/gene_attribute_matrix_standardized.mat', '-struct', 'gene_atb');
if issparse(gene_atb.cm.matrix)
    gene_atb.cm.matrix = full(gene_atb.cm.matrix);
end
writecm([savepath 'gene_attribute_matrix_standardized'], gene_atb.cm);
gzip([savepath 'gene_attribute_matrix_standardized.txt']);
delete([savepath 'gene_attribute_matrix_standardized.txt']);


