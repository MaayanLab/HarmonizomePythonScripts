


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




gene_atb.edges = edgesinit(100000, [], 'GeneSym', [], [], [], [], [], 'Chemical_CellLine', [], [], [], [], true, []);

filename = {'20093' '20094' '20095' '20097' '20098' '20099' '20100' '20101' '20102' '20103' '20104' '20105' '20106' '20204' '20206' '20207' '20209' '20210'};

i = 0;

for j = 1:1:numel(filename)
    
    fid = fopen(['input/Files20150120/' filename{j} '.txt'], 'r');
    
    currline = fgetl(fid);
    
    while ~feof(fid);
        
        currline = strrep(fgetl(fid), '"', '');
        
        currcells = strsplit(currline, '\t', 'CollapseDelimiters', false);
        
        chemical = currcells{2};
        cellline = currcells{8};
        genes = strsplit(currcells{9}, ',');
        weight = str2double(currcells{13});
        
        if numel(chemical) > 0 && numel(cellline) > 0  && numel(weight) > 0 && ~isnan(weight)
            
            for k = 1:1:numel(genes)
                
                subcells = strsplit(strrep(strrep(genes{k}, '_', ' '), '(', ' '));
                
                gene = subcells{1};
                
                if numel(gene) > 0
                    
                    i = i + 1;
                    
                    gene_atb.edges.source{i} = gene;
                    
                    gene_atb.edges.target{i} = [chemical '_' cellline];
                    
                    gene_atb.edges.weight(i) = weight/100;
                    
                end
                
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
gene_atb.cm = edges2cm(gene_atb.edges);
gene_atb = rmfield(gene_atb, 'edges');


% view distributions of row and col stats
[~] = cmrowcolstats(gene_atb.cm, true, 1);


% cluster and view matrix
[gene_atb.cm, ~] = cmcluster(gene_atb.cm, false);
% close force all;


% convert to sparse matrix
gene_atb.cm.matrix = sparse(gene_atb.cm.matrix);


% save result
save('../ImportedData/gene_chemical-cellline_kinativ_20150120.mat', '-struct', 'gene_atb');






clearvars -except gene_atb;


% convert to full matrix
gene_atb.cm.matrix = full(gene_atb.cm.matrix);


% % remove rows and cols with too many connections
% threshfrac = 1/2;
% gene_atb.cm = cmtrim_frac(gene_atb.cm, 0, threshfrac, 0, threshfrac, 'column');
% gene_atb.cm = cmtrim(gene_atb.cm, 1, Inf, 1, Inf);


% threshold values
posweight = sort(gene_atb.cm.matrix(gene_atb.cm.matrix > 0), 'descend');
ub = posweight(round(0.05*numel(posweight)));  % 0.868
% ub = 0.8;

negweight = sort(gene_atb.cm.matrix(gene_atb.cm.matrix < 0), 'ascend');
lb = negweight(round(0.05*numel(negweight)));  % -0.428
% lb = -0.8;

gene_atb.cm.matrix = (gene_atb.cm.matrix >= ub) - (gene_atb.cm.matrix <= lb);

gene_atb.cm = cmtrim(gene_atb.cm, 1, Inf, 1, Inf);


% cluster and view matrix
[gene_atb.cm, ~] = cmcluster(gene_atb.cm, true);
% close force all;


% view distributions of row and col stats
[~] = cmrowcolstats(gene_atb.cm, true, 2);


% convert to sparse matrix
gene_atb.cm.matrix = sparse(gene_atb.cm.matrix);


% save result
save('../PreparedData/gene_chemical-cellline_kinativ_20150120.mat', '-struct', 'gene_atb');







% get cleaned matrix

gene_atb = load('../ImportedData/gene_chemical-cellline_kinativ_20150120.mat', '-mat');
gene_atb.cm.matrix = full(gene_atb.cm.matrix);

if sum(sum(gene_atb.cm.matrix == 0)) > 1/3*numel(gene_atb.cm.matrix)
    gene_atb.cm.matrix = sparse(gene_atb.cm.matrix);
end
save('../CleanedData/gene_chemical-cellline_kinativ_20150120.mat', '-struct', 'gene_atb');
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

threshfrac = 0.1;
type = 'tertiary';
method = 'matrix';
discardemptyvectors = true;
[~, gene_atb.cm] = cmthresh_ks(gene_atb.cm, threshfrac, type, method, discardemptyvectors);
% gene_atb.cm = cm_standardize_ignorezeros(gene_atb.cm);

% cluster and view matrix
[gene_atb.cm, ~] = cmcluster(gene_atb.cm, false);

% save result
if sum(sum(gene_atb.cm.matrix == 0)) > 1/3*numel(gene_atb.cm.matrix)
    gene_atb.cm.matrix = sparse(gene_atb.cm.matrix);
end
save('../StandardizedData/gene_chemical-cellline_kinativ_20150120.mat', '-struct', 'gene_atb');
if issparse(gene_atb.cm.matrix)
    gene_atb.cm.matrix = full(gene_atb.cm.matrix);
end
writecm([savepath 'gene_attribute_matrix_standardized'], gene_atb.cm);
gzip([savepath 'gene_attribute_matrix_standardized.txt']);
delete([savepath 'gene_attribute_matrix_standardized.txt']);


