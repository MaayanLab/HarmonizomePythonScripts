


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




gene_atb.lists = listsinit(12589, [], 'GeneSym', [], [], [], 'GeneID', [], 'Disease', [], 'UMLS CUI', [], [], false, []);

fid = fopen('input/GeneID-Disease.txt', 'r');

for i = 1:1:4
    
    currline = fgetl(fid);
    
end

for i = 1:1:gene_atb.lists.numterms
    
    currline = strrep(strrep(fgetl(fid), '(', '|'), ')', '');
    
    currcells = strsplit(currline, '\t', 'CollapseDelimiters', false);
    
    b = find(currcells{1} == '|');
    
    gene_atb.lists.term{i} = currcells{1}(1:b-1);
    
    gene_atb.lists.termid(i) = str2double(currcells{1}(b+1:end));
    
    if numel(currcells) > 1
    
        gene_atb.lists.entries{i} = repmat({'-666'}, numel(currcells)-1, 1);
        gene_atb.lists.entrydescs{i} = repmat({'-666'}, numel(currcells)-1, 1);

        for j = 2:1:numel(currcells)

            b = find(currcells{j} == '|');

            gene_atb.lists.entries{i}{j-1} = currcells{j}(1:b-1);

            gene_atb.lists.entrydescs{i}{j-1} = currcells{j}(b+1:end);

        end

        gene_atb.lists.numentries(i) = numel(gene_atb.lists.entries{i});    
    
    else
        
        gene_atb.lists.numentries(i) = 0;
        
    end
    
end

fclose(fid);

discard = gene_atb.lists.numentries == 0;
gene_atb.lists = listsdiscard(gene_atb.lists, discard);

gene_atb.edges = lists2edges(gene_atb.lists);

gene_atb.edges = edgesunique(gene_atb.edges);


% map gene identifiers to NCBI Entrez Gene Symbols and Gene IDs
[gene_atb.edges.source, gene_atb.edges.sourcename, gene_atb.edges.sourceid, gene_atb.edges.sourceidname, discard] = genesymlookup([], gene_atb.edges.sourceid, 'gene', 'human', false, false, mappingfilespath);
gene_atb.edges = edgesdiscard(gene_atb.edges, discard);


% convert to matrix format (attribute table)
gene_atb.cm = edges2cm(gene_atb.edges);
gene_atb = rmfield(gene_atb, {'edges' 'lists'});


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


% view distributions of row and col stats
[~] = cmrowcolstats(gene_atb.cm, true, 2);


% convert to sparse matrix
gene_atb.cm.matrix = sparse(gene_atb.cm.matrix);


% save result
save('output/gene_attribute_matrix_prepared.mat', '-struct', 'gene_atb');


