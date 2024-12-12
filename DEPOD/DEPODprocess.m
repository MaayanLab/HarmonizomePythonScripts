


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




% initialize edges structure (atb=phosphatase, gene=substrate)
gene_atb.edges = edgesinit(1095, [], 'GeneSym', [], 'UniprotACC', [], 'GeneID', [], 'GeneSym', [], 'UniprotACC', [], 'GeneID', false, []);
discard = true(gene_atb.edges.numedges, 1);


% read data
fid = fopen('input/DEPOD_201405_human_phosphatase-substrate.mitab', 'r');

currline = fgetl(fid);

for i = 1:1:gene_atb.edges.numedges
    
    currline = fgetl(fid);
    
    si = regexp(currline, '\t');
    
    source = currline(si(1)+1:si(2)-1);
    target = currline(1:si(1)-1);
    
    if numel(strfind(source, 'uniprotkb:')) > 0 && numel(strfind(target, 'uniprotkb:')) > 0
        
        discard(i) = false;
        
        gene_atb.edges.sourcedesc{i} = source(11:end);
        
        gene_atb.edges.targetdesc{i} = target(11:end);
        
    end
    
end

fclose(fid);

gene_atb.edges = edgesdiscard(gene_atb.edges, discard);


% map gene identifiers to NCBI Entrez Gene Symbols and Gene IDs
uniprot_entrez = load('C:/Users/Andrew/Dropbox/GeneInfo/uniprot_entrez_resource_human_withisoforms.mat', '-mat');
for i = 1:1:gene_atb.edges.numedges
    si = regexp(gene_atb.edges.sourcedesc{i}, '-');
    if numel(si) > 0
        gene_atb.edges.sourcedesc{i} = gene_atb.edges.sourcedesc{i}(1:si(1)-1);
    end
end
gene_atb.edges.sourceid = -666*ones([gene_atb.edges.numedges 1]);
gene_atb.edges.sourceidname = 'GeneID';
gene_atb.edges.source = repmat({'-666'}, gene_atb.edges.numedges, 1);
gene_atb.edges.sourcename = 'GeneSym';
[o1, o2] = ismember(gene_atb.edges.sourcedesc, uniprot_entrez.edges.source);
o2(o2 == 0) = [];
gene_atb.edges.source(o1) = uniprot_entrez.edges.target(o2);
gene_atb.edges.sourceid(o1) = uniprot_entrez.edges.targetid(o2);
gene_atb.edges = edgesdiscard(gene_atb.edges, ~o1);


% map attribute identifiers to NCBI Entrez Gene Symbols and Gene IDs
uniprot_entrez = load('C:/Users/Andrew/Dropbox/GeneInfo/uniprot_entrez_resource_human_withisoforms.mat', '-mat');
for i = 1:1:gene_atb.edges.numedges
    si = regexp(gene_atb.edges.targetdesc{i}, '-');
    if numel(si) > 0
        gene_atb.edges.targetdesc{i} = gene_atb.edges.targetdesc{i}(1:si(1)-1);
    end
end
gene_atb.edges.targetid = -666*ones([gene_atb.edges.numedges 1]);
gene_atb.edges.targetidname = 'GeneID';
gene_atb.edges.target = repmat({'-666'}, gene_atb.edges.numedges, 1);
gene_atb.edges.targetname = 'GeneSym';
[o1, o2] = ismember(gene_atb.edges.targetdesc, uniprot_entrez.edges.source);
o2(o2 == 0) = [];
gene_atb.edges.target(o1) = uniprot_entrez.edges.target(o2);
gene_atb.edges.targetid(o1) = uniprot_entrez.edges.targetid(o2);
gene_atb.edges = edgesdiscard(gene_atb.edges, ~o1);


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
% [gene_atb.cm, ~] = cmcluster(gene_atb.cm, true);


% view distributions of row and col stats
[~] = cmrowcolstats(gene_atb.cm, true, 2);


% convert to sparse matrix
gene_atb.cm.matrix = sparse(gene_atb.cm.matrix);


% save result
save('output/gene_attribute_matrix_prepared.mat', '-struct', 'gene_atb');


