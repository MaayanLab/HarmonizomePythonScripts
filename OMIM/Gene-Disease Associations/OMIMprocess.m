


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



% get omim to entrez gene mapping

omim_entrez.edges = edgesinit(100000, [], 'OMIMID', [], [], [], [], [], 'GeneSym', [], [], [], 'GeneID', false, []);

fid = fopen('input/mim2gene_20150115.txt', 'r');

currline = fgetl(fid);

i = 0;

while ~feof(fid)
    
    currline = fgetl(fid);
    
    currcells = strsplit(currline, '\t', 'CollapseDelimiters', false);
    
    if strcmp(currcells{2}, 'gene')
        
        i = i + 1;
        
        omim_entrez.edges.source{i} = ['OMIM:' currcells{1}];
        
        omim_entrez.edges.target{i} = currcells{4};
        
        omim_entrez.edges.targetid(i) = str2double(currcells{3});
        
    end
    
end

fclose(fid);

if i < omim_entrez.edges.numedges
%     discard = false([omim_entrez.edges.numedges 1]);
    discard = isnan(omim_entrez.edges.targetid);
    discard(i+1:end) = true;
    omim_entrez.edges = edgesdiscard(omim_entrez.edges, discard);
end

[omim_entrez.edges.target, omim_entrez.edges.targetname, omim_entrez.edges.targetid, omim_entrez.edges.targetidname, discard] = genesymlookup([], omim_entrez.edges.targetid, 'gene', 'human', false, false, mappingfilespath);
omim_entrez.edges = edgesdiscard(omim_entrez.edges, discard);

save('input/omim_entrez.mat', '-struct', 'omim_entrez');



% get disease genes

disease_gene.edges = edgesinit(6778, [], 'Disease', [], 'OMIMID', [], 'PhenotypeClass', [], 'GeneSymbol', [], 'OMIMID', [], 'GeneID', false, []);

mimregexp = '\d\d\d\d\d\d';

phnregexp = '\(\d\)';

fid = fopen('input/morbidmap_20150115.txt', 'r');

for i = 1:1:disease_gene.edges.numedges
    
    currline = fgetl(fid);
    
    currcells = strsplit(currline, '|', 'CollapseDelimiters', false);
    
    subcells = strsplit(currcells{1}, ', ');
    
    [b, e] = regexp(subcells{end}, phnregexp);
    
    if ~isempty(b)
        
        disease_gene.edges.sourceid(i) = str2double(subcells{end}(b+1:e-1));
        
        subcells{end} = subcells{end}(1:b-2);
        
    end
    
    [b, e] = regexp(subcells{end}, mimregexp);
    
    if ~isempty(b)
        
        disease_gene.edges.sourcedesc{i} = ['OMIM:' subcells{end}(b:e)];
        disease_gene.edges.source{i} = lower(strjoin(subcells(1:end-1), ', '));
        
    else
        
        disease_gene.edges.source{i} = lower(strjoin(subcells, ', '));
        
    end
    
    subcells = strsplit(currcells{2}, ', ');
    disease_gene.edges.target{i} = subcells{1};
    disease_gene.edges.targetdesc{i} = ['OMIM:' currcells{3}];
    
end

fclose(fid);

[disease_gene.edges.target, disease_gene.edges.targetname, disease_gene.edges.targetid, disease_gene.edges.targetidname, discard] = genesymlookup(disease_gene.edges.target, [], 'gene', 'human', false, false, mappingfilespath);

omim_entrez = load('input/omim_entrez.mat', '-mat');

[o1, o2] = ismember(disease_gene.edges.targetdesc, omim_entrez.edges.source);
o2(o2 == 0) = [];
disease_gene.edges.target(o1) = omim_entrez.edges.target(o2);
disease_gene.edges.targetid(o1) = omim_entrez.edges.targetid(o2);

discard = discard & ~o1;
disease_gene.edges = edgesdiscard(disease_gene.edges, discard);

disease_gene.cm = edges2cm(disease_gene.edges);

gene_disease.cm = cmtranspose(disease_gene.cm);

figure(1);
clf;
subplot(1,2,1);
hist(sum(gene_disease.cm.matrix, 1), 100);
ylabel('number of columns');
xlabel('column sum');
subplot(1,2,2);
hist(sum(gene_disease.cm.matrix, 2), 30);
ylabel('number of rows');
xlabel('row sum');

cgo = cm2clustergram(gene_disease.cm, 'none', 'all', 'cosine', 'average');

gene_disease.cm = conmatmap(gene_disease.cm, cgo.RowLabels, cgo.ColumnLabels');

HeatMap(gene_disease.cm.matrix, 'Colormap', redbluecmap);

gene_disease.cm.matrix = sparse(gene_disease.cm.matrix);

save('output/gene_attribute_matrix_imported.mat', '-struct', 'gene_disease');



% get disease list and fill in some missing omim ids
clearvars -except gene_disease;

disease.numterms = 100000;
disease.term = cell(disease.numterms, 1);
disease.termname = 'Disease';
disease.termdesc = cell(disease.numterms, 1);
disease.termdescname = 'OMIMID';

mimregexp = '\d\d\d\d\d\d';

fid = fopen('input/omim_20150115.txt', 'r');

i = 0;
linecount = 0;
tic;

while ~feof(fid)
    
    currline = fgetl(fid);
    linecount = linecount + 1;
    
    if numel(currline) >= 10 && strcmp(currline(1:10), '*FIELD* TI')
        
        currline = fgetl(fid);
        
        titleline = ' ';
        
        while isempty(strfind(currline, '*FIELD* TX'))
            
            titleline = [titleline ' ' currline];
            
            currline = fgetl(fid);
            
        end
        
        i = i + 1;
        
        [b, e] = regexp(titleline, mimregexp);
        
        disease.termdesc{i} = ['OMIM:' titleline(b:e)];
        
        sc = find(titleline == ';', 1, 'first');
        
        if ~isempty(sc)
        
            disease.term{i} = lower(titleline(e+2:sc-1));
            
        else
            
            disease.term{i} = lower(titleline(e+2:end));
            
        end
        
        while disease.term{i}(1) == ' '
            
            disease.term{i}(1) = [];
            
        end
        
    end
    
    if mod(linecount, 400000) == 0
        disp(['Finished line ' num2str(linecount) ' of 4181041.']);
        toc;
        disp(' ');
        tic;
    end
    
end

disp(['Finished line ' num2str(linecount) ' of 4181041.']);
toc;
disp('');

fclose(fid);

if i < disease.numterms
    disease.term(i+1:end) = [];
    disease.termdesc(i+1:end) = [];
end
disease.numterms = i;

missingid = strcmp(gene_disease.cm.entrydesc, '-666');

entrydesc = gene_disease.cm.entrydesc;
[o1, o2] = ismember(removespecialchars(gene_disease.cm.entry), removespecialchars(disease.term));
o2(o2 == 0) = [];
entrydesc(o1) = disease.termdesc(o2);
gene_disease.cm.entrydesc(missingid) = entrydesc(missingid);

save('output/gene_attribute_matrix_imported.mat', '-struct', 'gene_disease');

missingid = strcmp(gene_disease.cm.entrydesc, '-666');

disease.term = [gene_disease.cm.entry(~missingid); disease.term];
disease.termdesc = [gene_disease.cm.entrydesc(~missingid); disease.termdesc];
disease.numterms = numel(disease.term);

[disease.termdesc, ui, ~] = unique(disease.termdesc, 'stable');
disease.term = disease.term(ui);
disease.numterms = numel(disease.term);

disease.term = [disease.term; gene_disease.cm.entry(missingid)];
disease.termdesc = [disease.termdesc; gene_disease.cm.entrydesc(missingid)];
disease.numterms = numel(disease.term);

[disease.term, ui, ~] = unique(disease.term, 'stable');
disease.termdesc = disease.termdesc(ui);
disease.numterms = numel(disease.term);

save('output/gene_attribute_matrix_imported.mat', '-struct', 'disease');


