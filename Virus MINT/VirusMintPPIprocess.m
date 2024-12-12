


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




% initialize edges structure (gene=human genes, atb=viral genes)
gene_atb.edges = edgesinit(2707, [], 'UniprotACC', [], 'Organism', [], 'TaxID',  [], 'UniprotACC', [], 'Organism', [], 'TaxID', false, []);


% read data
fid = fopen('input/2012-10-26-mint-viruses-binary.mitab26.txt', 'r');

currline = fgetl(fid);

for i = 1:1:gene_atb.edges.numedges
    
    currline = strrep(fgetl(fid), '"', '');
    
    currcells = strsplitbyadr(currline, '\t');
    
    gene_atb.edges.source{i} = currcells{1}(11:end);
    
    gene_atb.edges.target{i} = currcells{2}(11:end);
    
    taxinfo = currcells{10}(7:end);
    
    p = regexp(taxinfo, '(');
    
    gene_atb.edges.sourceid(i) = str2double(taxinfo(1:p-1));
    
    gene_atb.edges.sourcedesc{i} = taxinfo(p+1:end-1);
    
    taxinfo = currcells{11}(7:end);
    
    p = regexp(taxinfo, '(');
    
    gene_atb.edges.targetid(i) = str2double(taxinfo(1:p-1));
    
    gene_atb.edges.targetdesc{i} = taxinfo(p+1:end-1);
    
end

fclose(fid);


% arrange data so human genes are source and viral genes are target
gene_atb.edges = edges2symedges(gene_atb.edges);
discard = gene_atb.edges.targetid == 9606 | gene_atb.edges.sourceid ~= 9606;
gene_atb.edges = edgesdiscard(gene_atb.edges, discard);


% map gene identifiers to NCBI Entrez Gene Symbols and Gene IDs
uniprot_entrez = load('C:/Users/Andrew/Dropbox/GeneInfo/uniprot_entrez_resource_human_withisoforms.mat', '-mat');
gene_atb.edges.sourcedesc = gene_atb.edges.source;
gene_atb.edges.sourcedescname = gene_atb.edges.sourcename;
gene_atb.edges.source = repmat({'-666'}, gene_atb.edges.numedges, 1);
gene_atb.edges.sourcename = 'GeneSym';
gene_atb.edges.sourceid = -666*ones([gene_atb.edges.numedges 1]);
gene_atb.edges.sourceidname = 'GeneID';
gene_atb.edges.sourcedesc = cellfun(@(x) strsplitbyadr(x, '-'), gene_atb.edges.sourcedesc, 'uniformoutput', false);
gene_atb.edges.sourcedesc = cellfun(@(x) x{1}, gene_atb.edges.sourcedesc, 'uniformoutput', false);
[o1, o2] = ismember(gene_atb.edges.sourcedesc, uniprot_entrez.edges.source);
o2(o2 == 0) = [];
gene_atb.edges.source(o1) = uniprot_entrez.edges.target(o2);
gene_atb.edges.sourceid(o1) = uniprot_entrez.edges.targetid(o2);
gene_atb.edges = edgesdiscard(gene_atb.edges, ~o1);
clear uniprot_entrez;


% remove duplicate edges
gene_atb.edges = edgesunique(gene_atb.edges);


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
save('output/gene_attribute_matrix_imported.mat', '-struct', 'gene_atb');






clearvars -except gene_atb;


% convert to full matrix
gene_atb.cm.matrix = full(gene_atb.cm.matrix);


% remove rows and cols with too many connections
threshfrac = 1/2;
gene_atb.cm = cmtrim_frac(gene_atb.cm, 0, threshfrac, 0, threshfrac, 'column');
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





