


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
gene_site.edges = edgesinit(853, [], 'GeneSym', [], 'UniprotAcc', [], 'GeneID', [], 'Site', [], 'Disease', [], 'SiteGroupID', false, []);

discard = false([gene_site.edges.numedges 1]);

fid = fopen('input/Disease-associated_sites_20150324_fixed.txt', 'r');

currline = fgetl(fid);
currline = fgetl(fid);
currline = fgetl(fid);
currline = fgetl(fid);

for i = 1:1:gene_site.edges.numedges
    
    currline = fgetl(fid);
    
    currcells = strsplit(currline, '\t', 'CollapseDelimiters', false);
    
    if strcmpi(currcells{9}, 'human') && numel(currcells{6}) > 0
        
        gene_site.edges.source{i} = currcells{6};
        
        gene_site.edges.sourcedesc{i} = currcells{4};
        
        gene_site.edges.sourceid(i) = str2double(currcells{5});
        
        gene_site.edges.target{i} = currcells{11};
        
        gene_site.edges.targetdesc{i} = currcells{1};
        
        gene_site.edges.targetid(i) = str2double(currcells{10});
        
    else
        
        discard(i) = true;
        
    end
    
end

fclose(fid);

gene_site.edges = edgesdiscard(gene_site.edges, discard);

save('input/gene_site_human_disease_unmodified_20150324.mat', '-struct', 'gene_site');
%}

% load data
gene_atb = load('input/gene_site_human_disease_unmodified_20150324.mat', '-mat', 'edges');


% map gene identifiers to NCBI Entrez Gene Symbols and Gene IDs
[gene_atb.edges.source, gene_atb.edges.sourcename, gene_atb.edges.sourceid, gene_atb.edges.sourceidname, discard] = genesymlookup(gene_atb.edges.source, [], 'gene', 'both', true, false, mappingfilespath);
gene_atb.edges = edgesdiscard(gene_atb.edges, discard);


% save intermediate result
save('output/gene_attribute_matrix_imported.mat', '-struct', 'gene_atb');
substrate = struct;
[substrate.term, ui, ri] = unique(gene_atb.edges.source);
substrate.termname = gene_atb.edges.sourcename;
substrate.termid = gene_atb.edges.sourceid(ui);
substrate.termidname = gene_atb.edges.sourceidname;
substrate.numterms = numel(substrate.term);
save('output/gene_attribute_matrix_imported.mat', '-struct', 'substrate');


% swap data so target=disease and targetdesc=site
target = gene_atb.edges.target;
targetname = gene_atb.edges.targetname;
gene_atb.edges.target = gene_atb.edges.targetdesc;
gene_atb.edges.targetname = gene_atb.edges.targetdescname;
gene_atb.edges.targetdesc = target;
gene_atb.edges.targetdescname = targetname;


% discard duplicate edges
gene_atb.edges = edgesunique(gene_atb.edges);


% convert to matrix format (attribute table)
gene_atb.cm = edges2cm(gene_atb.edges);
gene_atb = rmfield(gene_atb, 'edges');


% view distributions of row and col stats
[~] = cmrowcolstats(gene_atb.cm, true, 1);


% cluster and view matrix
[gene_atb.cm, ~] = cmcluster(gene_atb.cm, true);
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
% close force all;


% view distributions of row and col stats
% [~] = cmrowcolstats(gene_atb.cm, true, 2);


% convert to sparse matrix
gene_atb.cm.matrix = sparse(gene_atb.cm.matrix);


% save result
save('output/gene_attribute_matrix_prepared.mat', '-struct', 'gene_atb');


