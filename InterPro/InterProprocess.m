


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





% initialize attribute structure
domain.numterms = 100000;
domain.term = cell(domain.numterms, 1);
domain.termname = 'Domain';
domain.termdesc = cell(domain.numterms, 1);
domain.termdescname = 'DomainID';


% read attribute list
fid = fopen('input/interpro.xml', 'r');

i = 0;

while ~feof(fid)
    
    currline = fgetl(fid);
    
    if numel(currline) > 13 && strcmp(currline(1:13), '<interpro id=')
        
        i = i + 1;
        
        q = find(currline == '"');
        
        domain.termdesc{i} = currline(q(1)+1:q(2)-1);
        
        gotname = false;
        
        while ~gotname
            
            currline = fgetl(fid);
            
            t1 = strfind(currline, '<name>');
            
            if ~isempty(t1)
                
                t2 = strfind(currline, '</name>');
                
                domain.term{i} = currline(t1+6:t2-1);
                
                gotname = true;
                
            end
            
        end
        
        disp(['i = ' num2str(i)]);
        
    end
    
end

fclose(fid);

if i < domain.numterms
    domain.term(i+1:end) = [];
    domain.termdesc(i+1:end) = [];
end
domain.numterms = numel(domain.term);

[domain.termdesc, ui, ~] = unique(domain.termdesc);
domain.term = domain.term(ui);
domain.numterms = numel(domain.term);


% initialize lists structure
gene_atb.lists = listsinit(20201, [], 'GeneSym', [], 'UniprotAcc', [], 'GeneID', [], 'DomainID', [], [], [], [], false, []);

discard = false([gene_atb.lists.numterms 1]);


% read data
fid = fopen('input/uniprot_human_domains_20150324.tab', 'r');

currline = fgetl(fid);

for i = 1:1:gene_atb.lists.numterms
    
    currline = upper(fgetl(fid));
    
    currcells = strsplit(currline, '\t', 'CollapseDelimiters', false);
    
    if numel(currcells{3}) > 0 && numel(currcells{4}) > 0 && numel(currcells{5}) > 3

        gene_atb.lists.term{i} = currcells{3};
        gene_atb.lists.termdesc{i} = currcells{1};
        gene_atb.lists.termid(i) = str2double(currcells{4}(1:find(currcells{4} == ';', 1, 'first')-1));
        
        subcells = strsplit(currcells{5}, ';')';
        disc = strcmp(subcells, '');
        gene_atb.lists.entries{i} = subcells(~disc);
        gene_atb.lists.numentries(i) = numel(gene_atb.lists.entries{i});
    
    else
        
        discard(i) = true;
        
    end
    
end

fclose(fid);

gene_atb.lists = listsdiscard(gene_atb.lists, discard);

gene_atb.edges = lists2edges(gene_atb.lists);

gene_atb.lists = edges2lists(gene_atb.edges);


% convert to matrix format (attribute table)
gene_atb.cm = lists2cm(gene_atb.lists);
gene_atb = rmfield(gene_atb, {'lists' 'edges'});


% map gene identifiers to NCBI Entrez Gene Symbols and Gene IDs
[gene_atb.cm.term, gene_atb.cm.termname, gene_atb.cm.termid, gene_atb.cm.termidname, discard] = genesymlookup([], gene_atb.cm.termid, 'gene', 'human', false, false, mappingfilespath);
gene_atb.cm = cmrowdiscard(gene_atb.cm, discard);


% merge rows corresponding to the same gene
if numel(unique(gene_atb.cm.term)) < gene_atb.cm.numterms
    gene_atb.cm = cmrowmerge(gene_atb.cm, 'max');
end


% map attributes to attribute list
gene_atb.cm.entrydesc = gene_atb.cm.entry;
gene_atb.cm.entrydescname = gene_atb.cm.entryname;
[o1, o2] = ismember(gene_atb.cm.entrydesc, domain.termdesc);
o2(o2 == 0) = [];
gene_atb.cm.entry(o1) = domain.term(o2);
gene_atb.cm.entryname = domain.termname;
gene_atb.cm = cmcoldiscard(gene_atb.cm, ~o1);


% merge columns corresponding to the same attribute
if numel(unique(gene_atb.cm.entry)) < gene_atb.cm.numentries
    gene_atb.cm = cmcolmerge(gene_atb.cm, 'max');
end


% remove empty rows and cols
gene_atb.cm = cmtrim(gene_atb.cm, 1, Inf, 1, Inf);


% view distributions of row and col stats
[~] = cmrowcolstats(gene_atb.cm, true, 1);


% cluster and view matrix
[gene_atb.cm, ~] = cmcluster(gene_atb.cm, true);


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


