


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




gene_locus = load('C:\Users\Andrew\Dropbox\GeneInfo\gene_loci_GRCh37_hg19.mat', '-mat');

numsites = 4474877;

gene_tf.edges = edgesinit(numsites, [], 'GeneSym', [], [], [], 'GeneID', [], 'TF Name', [], 'ZScore', [], 'GeneID', true, []);

fid = fopen('input/HUMAN_hg19_BBLS_1_00_FDR_0_10.bed', 'r');
k = 0;
tic;
for i = 1:1:numsites
    
    currline = fgetl(fid);
    
    currcells = strsplit(currline, '\t', 'CollapseDelimiters', false);
    
    chr = currcells{1};
    
    startpos = str2double(currcells{2});
    stoppos = str2double(currcells{3});
    meanpos = (startpos + stoppos)/2;
    
    distances = abs((gene_locus.StartPos - meanpos)./double(strcmp(gene_locus.ChrNum, chr)));

    hidx = find(distances == min(distances));
    
    mindist = distances(hidx(1));
    
    if mindist < 2000
        
        for j = 1:1:numel(hidx)
            
            k = k + 1;
            
            gene_tf.edges.source{k} = gene_locus.GeneSym{hidx(j)};
            gene_tf.edges.sourceid(k) = gene_locus.GeneID(hidx(j));
            gene_tf.edges.target{k} = currcells{4}(find(currcells{4} == '=', 1, 'last')+1:end);
            gene_tf.edges.targetdesc{k} = currcells{5};
            gene_tf.edges.weight(k) = -mindist/2000 + 1;
            
        end
        
    end
    
    if mod(i, 100000) == 0
        disp(['Finished ' num2str(i) ' of ' num2str(numsites) '.']);
        toc;
        disp(' ');
        tic;
    end
    
end
disp(['Finished ' num2str(i) ' of ' num2str(numsites) '.']);
toc;
disp(' ');
fclose(fid);

if k < gene_tf.edges.numedges
    discard = false([gene_tf.edges.numedges 1]);
    discard(k+1:end) = true;
    gene_tf.edges = edgesdiscard(gene_tf.edges, discard);
end

gene_tf.edges = edgesunique(gene_tf.edges);

gene_tf.cm = edges2cm(gene_tf.edges);

gene_tf = rmfield(gene_tf, 'edges');

[gene_tf.cm.term, gene_tf.cm.termname, gene_tf.cm.termid, gene_tf.cm.termidname, discard] = genesymlookup(gene_tf.cm.term, [], 'gene', 'human', false, false, mappingfilespath);
gene_tf.cm = cmrowdiscard(gene_tf.cm, discard);
gene_tf.cm = cmrowmerge(gene_tf.cm, 'max');



[gene_tf.cm.entrydesc, gene_tf.cm.entrydescname, gene_tf.cm.entryid, gene_tf.cm.entryidname, discard] = genesymlookup(gene_tf.cm.entry, [], 'gene', 'human', true, false, mappingfilespath);

nummissed = 142;
tfname = cell(nummissed, 1);
genesym = cell(nummissed, 1);
geneid = zeros([nummissed 1]);

fid = fopen('input/missedTFs.txt', 'r');

currline = fgetl(fid);

for i = 1:1:nummissed
    
    currline = fgetl(fid);
    
    currcells = strsplit(currline, '\t', 'CollapseDelimiters', false);
    
    tfname{i} = currcells{1};
    
    genesym{i} = currcells{2};
    
    geneid(i) = str2double(currcells{3});
    
end

fclose(fid);

geneid(isnan(geneid)) = -666;

[o1, o2] = ismember(gene_tf.cm.entry, tfname);
o2(o2 == 0) = [];

gene_tf.cm.entrydesc(o1) = genesym(o2);
gene_tf.cm.entryid(o1) = geneid(o2);

discard = discard & ~o1;
gene_tf.cm = cmcoldiscard(gene_tf.cm, discard);
gene_tf.cm = cmcolmerge(gene_tf.cm, 'max');



figure(1);
clf;
subplot(1,2,1);
hist(sum(gene_tf.cm.matrix, 1), 100);
ylabel('number of columns');
xlabel('column sum');
subplot(1,2,2);
hist(sum(gene_tf.cm.matrix, 2), 30);
ylabel('number of rows');
xlabel('row sum');

[gene_tf.cm, cgo] = cmcluster(gene_tf.cm, true);

gene_tf.cm.matrix = sparse(gene_tf.cm.matrix);

save('output/gene_attribute_matrix_imported.mat', '-struct', 'gene_tf');

gene_tf.cm.matrix = full(gene_tf.cm.matrix);

gene_tf.cm.matrix = double(gene_tf.cm.matrix > 0);

gene_tf.cm.matrix = sparse(gene_tf.cm.matrix);

save('output/gene_attribute_matrix_imported.mat', '-struct', 'gene_tf');

temp = gene_tf.cm.entry;
tempname = gene_tf.cm.entryname;
gene_tf.cm.entry = gene_tf.cm.entrydesc;
gene_tf.cm.entryname = gene_tf.cm.entrydescname;
gene_tf.cm.entrydesc = temp;
gene_tf.cm.entrydescname = tempname;
discard = gene_tf.cm.entryid == -666;
gene_tf.cm = cmcoldiscard(gene_tf.cm, discard);
gene_tf.cm = cmcolmerge(gene_tf.cm, 'max');
gene_tf.cm = cmtrim(gene_tf.cm, 1, Inf, 1, Inf);

[gene_tf.cm, cgo] = cmcluster(gene_tf.cm, true);

gene_tf.cm.matrix = sparse(gene_tf.cm.matrix);

save('output/gene_attribute_matrix_imported.mat', '-struct', 'gene_tf');


