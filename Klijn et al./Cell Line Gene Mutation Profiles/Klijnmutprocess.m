


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
gene_atb.edges = edgesinit(280535, [], 'GeneSym', [], 'Chromosome_Position_Context', [], 'GeneID', [], 'Cell Line', [], 'Tissue', [], [], false, []);
context = repmat({'-666'}, gene_atb.edges.numedges, 1);


% read data
fid = fopen('input/Data3_genemutations/140331_cellLineMutations.txt', 'r');

currline = fgetl(fid);

for i = 1:1:gene_atb.edges.numedges
    
    currline = fgetl(fid);
    
    currcells = strsplit(currline, '\t', 'CollapseDelimiters', false);
    
    gene_atb.edges.source{i} = currcells{27};

    gene_atb.edges.sourceid(i) = str2double(currcells{5});
    
    gene_atb.edges.sourcedesc{i} = [currcells{16} '_' currcells{17} '_' currcells{8}];
    
    context{i} = currcells{8};
    
    gene_atb.edges.target{i} = upper(currcells{66});
    
    gene_atb.edges.targetdesc{i} = lower(currcells{58});
    
end

fclose(fid);


% write cell line tissue info
% [ucl, ui, ~] = unique(gene_atb.edges.target);
% ut = gene_atb.edges.targetdesc(ui);
% fid = fopen('input/cellline_tissue_list.txt', 'w');
% for i = 1:1:numel(ucl);
%     fprintf(fid, '%s\t%s\r\n', ucl{i}, ut{i});
% end
% fclose(fid);
% clear ucl ui ut;


% map gene identifiers to NCBI Entrez Gene Symbols and Gene IDs
[gene_atb.edges.source, gene_atb.edges.sourcename, gene_atb.edges.sourceid, gene_atb.edges.sourceidname, discard] = genesymlookup([], gene_atb.edges.sourceid, 'gene', 'human', false, false, mappingfilespath);
gene_atb.edges = edgesdiscard(gene_atb.edges, discard);


% remove duplicate edges
gene_atb.edges = edgesunique(gene_atb.edges);


% convert to matrix format (attribute table)
gene_atb.cm = edges2cm(gene_atb.edges);
gene_atb = rmfield(gene_atb, 'edges');


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







% map cell lines and tissues
ontologies = {'mouse_adult_gross_anatomy_ontology_nodes' 'uber_anatomy_ontology_nodes' 'experimental_factor_ontology_nodes' 'cell_ontology_nodes' 'cell_line_ontology_nodes' 'brenda_tissue_enzyme_source_ontology_nodes'}';

folder = {'ImportedData' 'PreparedData'}';
attribute = 'cellline';
filesuffix = 'klijn_mut_20150415';

for fi = 1:1:numel(folder)
    
    gene_atb = load(['../' folder{fi} '/gene_' attribute '_' filesuffix '.mat'], '-mat', 'cm');

    if strcmpi(gene_atb.cm.entryname, 'cellline') || strcmpi(gene_atb.cm.entryname, 'cell line') || strcmpi(gene_atb.cm.entryname, 'tissue')

        [gene_atb.cm.entry, gene_atb.cm.entryname, isunmapped] = celltissuelabelmapper(gene_atb.cm.entry, gene_atb.cm.entryname, ontologies);

        writeunmappedlabels(gene_atb.cm.entry, gene_atb.cm.entryname, isunmapped, ['../UnmappedAttributes/' folder{fi} '_entry_LABELNAME_' filesuffix]);

    end

    if strcmpi(gene_atb.cm.entrydescname, 'cellline') || strcmpi(gene_atb.cm.entrydescname, 'cell line') || strcmpi(gene_atb.cm.entrydescname, 'tissue')

        [gene_atb.cm.entrydesc, gene_atb.cm.entrydescname, isunmapped] = celltissuelabelmapper(gene_atb.cm.entrydesc, gene_atb.cm.entrydescname, ontologies);
        
        writeunmappedlabels(gene_atb.cm.entrydesc, gene_atb.cm.entrydescname, isunmapped, ['../UnmappedAttributes/' folder{fi} '_entrydesc_LABELNAME_' filesuffix]);

    end
    
    save(['../' folder{fi} '/gene_' attribute '_' filesuffix '.mat'], '-struct', 'gene_atb');
    
end


