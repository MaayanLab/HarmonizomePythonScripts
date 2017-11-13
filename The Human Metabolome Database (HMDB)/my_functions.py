import sys, datetime, os
import numpy as np
import pandas as pd
import scipy.stats.mstats as stat
import scipy.spatial.distance as dist


def merge(inputDF, axis, method):

    if axis == 'column' or axis == 'Column' and method == 'mean':
        inputDF = inputDF.groupby(inputDF.columns, axis=1).mean()
        return(inputDF)

    elif axis == 'row' or axis == 'Row' and method == 'mean':
        inputDF = inputDF.groupby(level=0, axis=0).mean()
        return(inputDF)

def zscore(inputDF, axis):
    if axis == 'row':
        for i,index in enumerate(inputDF.index):

            progressPercent = ((i+1)/len(inputDF.index))*100

            sys.stdout.write("Progress: %d%%  %d Out of %d   \r" % (progressPercent, (i+1), len(inputDF.index)))
            sys.stdout.flush()

            mean = inputDF.ix[index].mean()
            std = inputDF.ix[index].std()
            inputDF.ix[index] = inputDF.ix[index].apply(lambda x: ((x-mean)/std))
    elif axis == 'column':
        for i,col in enumerate(inputDF.columns):

            progressPercent = ((i+1)/len(inputDF.columns))*100

            sys.stdout.write("Progress: %d%%  %d Out of %d   \r" % (progressPercent, (i+1), len(inputDF.columns)))
            sys.stdout.flush()

            mean = inputDF[col].mean()
            std = inputDF[col].std()
            inputDF[col] = inputDF[col].apply(lambda x: ((x-mean)/std))


def log2(inputDF):
    DF = inputDF.copy()
    DF = DF.apply(lambda x: np.log2(x+1))
    return(DF)


def quantileNormalize(inputDF):
    df = inputDF.copy()
    #compute rank
    dic = {}
    for i,col in enumerate(df):

        progressPercent = ((i+1)/len(df.columns))*100

        sys.stdout.write("Step 1/2 progress: %d%%  %d Out of %d   \r" % (progressPercent, (i+1), len(df.columns)))
        sys.stdout.flush()

        dic.update({col : sorted(df[col])})

    sorted_df = pd.DataFrame(dic)
    rank = sorted_df.mean(axis = 1).tolist()
    #sort
    for i,col in enumerate(df):

        progressPercent = ((i+1)/len(df.columns))*100

        sys.stdout.write("Step 2/2 progress: %d%%  %d Out of %d   \r" % (progressPercent, (i+1), len(df.columns)))
        sys.stdout.flush()

        t = np.searchsorted(np.sort(df[col]), df[col])
        df[col] = [rank[i] for i in t]
    return df

mappingDF = pd.read_csv('/Users/moshesilverstein/Documents/Harmonizome/mappingFile_2017.txt', sep='\t', header=None, index_col=0)

def mapgenesymbols(inputDF):

    inputDF.reset_index(inplace=True)


    lst1 = []

    for i, index in enumerate(inputDF.index):

        progressPercent = ((i+1)/len(inputDF.index))*100

        sys.stdout.write("Progeres: %d%%  %d Out of %d   \r" % (progressPercent, (i+1), len(inputDF.index)))
        sys.stdout.flush()

        if inputDF.ix[index, inputDF.columns[0]] in mappingDF.index:
            lst1.append(mappingDF.ix[inputDF.ix[index, inputDF.columns[0]], 1])
        else:
            lst1.append(np.nan)


    inputDF[inputDF.columns[0]] = lst1


    inputDF.dropna(inplace=True, subset=[inputDF.columns[0]])
    inputDF.set_index(inputDF.columns[0], inplace=True)

def createTertiaryMarix(inputDF):
    df = inputDF.copy()
    for i,col in enumerate(inputDF.columns):

        progressPercent = ((i+1)/len(inputDF.columns))*100

        sys.stdout.write("Progeres: %d%%  %d Out of %d   \r" % (progressPercent, (i+1), len(inputDF.columns)))
        sys.stdout.flush()

        lst = []
        lst = inputDF[col].values.tolist()
        lst.sort()

        lst = lst[::-1]
        upCutoff = lst[int(0.1*len(lst))]
        index = inputDF[inputDF[col] > upCutoff].index
        df.ix[index, col] = 1

        zeroIndex = inputDF[inputDF[col] < upCutoff].index
        df.ix[zeroIndex, col] = 0

        lst = lst[::-1]
        downCutoff = lst[int(0.1*len(lst))]
        index = inputDF[inputDF[col] < downCutoff].index
        df.ix[index, col] = -1


    return(df)


def createUpGeneSetLib(inputDF, path, name, details=None):

    filenameGMT = name+'_%s.gmt'% str(datetime.date.today())[0:7].replace('-', '_')

    if os.path.isfile(path+filenameGMT):
        os.remove(path+filenameGMT)

    for i,col in enumerate(inputDF.columns):

        progressPercent = ((i+1)/len(inputDF.columns))*100

        sys.stdout.write("Progeres: %d%%  %d Out of %d   \r" % (progressPercent, (i+1), len(inputDF.columns)))
        sys.stdout.flush()

        index = inputDF[inputDF[col] == 1].index

        lst = index.values.tolist()

        if len(lst) > 5 and len(lst) <= 2000:

            lst.insert(0, col)
            if details:
                lst.insert(1, details[i])
            else:
                lst.insert(1, 'NA')
            lst = ['{0}\t'.format(elem) for elem in lst] # add tabs between terms in the lst
            lst.insert(len(lst), '\n') # add a newline char at the end of each lst

            with open(path+filenameGMT, 'a') as the_file:
                the_file.writelines(lst)

def createDownGeneSetLib(inputDF, path, name, details=None):

    filenameGMT = name+'_%s.gmt'% str(datetime.date.today())[0:7].replace('-', '_')

    if os.path.isfile(path+filenameGMT):
        os.remove(path+filenameGMT)

    for i,col in enumerate(inputDF.columns):

        progressPercent = ((i+1)/len(inputDF.columns))*100

        sys.stdout.write("Progeres: %d%%  %d Out of %d   \r" % (progressPercent, (i+1), len(inputDF.columns)))
        sys.stdout.flush()

        index = inputDF[inputDF[col] == -1].index

        lst = index.values.tolist()

        lst.insert(0, col)
        if details:
            lst.insert(1, details[i])
        else:
            lst.insert(1, 'NA')
        lst = ['{0}\t'.format(elem) for elem in lst] # add tabs between terms in the lst
        lst.insert(len(lst), '\n') # add a newline char at the end of each lst

        with open(path+filenameGMT, 'a') as the_file:
            the_file.writelines(lst)


def createUpAttributeSetLib(inputDF, path, name):

    inputDF = inputDF.T

    filenameGMT = name+'_%s.gmt'% str(datetime.date.today())[0:7].replace('-', '_')

    if os.path.isfile(path+filenameGMT):
        os.remove(path+filenameGMT)

    for i,col in enumerate(inputDF.columns):

        progressPercent = ((i+1)/len(inputDF.columns))*100

        sys.stdout.write("Progeres: %d%%  %d Out of %d   \r" % (progressPercent, (i+1), len(inputDF.columns)))
        sys.stdout.flush()

        index = inputDF[inputDF[col] == 1].index

        lst = index.values.tolist()

        if len(lst) > 5:
            lst.insert(0, col)
            lst.insert(1, 'NA')
            lst = ['{0}\t'.format(elem) for elem in lst] # add tabs between terms in the lst
            lst.insert(len(lst), '\n') # add a newline char at the end of each lst

            with open(path+filenameGMT, 'a') as the_file:
                the_file.writelines(lst)

    inputDF = inputDF.T

def createDownAttributeSetLib(inputDF, path, name):

    inputDF = inputDF.T

    filenameGMT = name+'_%s.gmt'% str(datetime.date.today())[0:7].replace('-', '_')

    if os.path.isfile(path+filenameGMT):
        os.remove(path+filenameGMT)

    for i,col in enumerate(inputDF.columns):

        progressPercent = ((i+1)/len(inputDF.columns))*100

        sys.stdout.write("Progeres: %d%%  %d Out of %d   \r" % (progressPercent, (i+1), len(inputDF.columns)))
        sys.stdout.flush()

        index = inputDF[inputDF[col] == -1].index

        lst = index.values.tolist()
        lst.insert(0, col)
        lst.insert(1, 'NA')
        lst = ['{0}\t'.format(elem) for elem in lst] # add tabs between terms in the lst
        lst.insert(len(lst), '\n') # add a newline char at the end of each lst

        with open(path+filenameGMT, 'a') as the_file:
            the_file.writelines(lst)

    inputDF = inputDF.T

def createSimilarityMatrix(inputDF, metric):
    similarity_matrix = dist.pdist(inputDF, metric)
    similarity_matrix = dist.squareform(similarity_matrix)
    similarity_df = pd.DataFrame(data=similarity_matrix[0:,0:], index=inputDF.index, columns=inputDF.index)
    similarity_df = similarity_df.applymap(lambda x: 1-x)
    return(similarity_df)

getGeneIDs = pd.read_csv('/Users/moshesilverstein/Documents/Harmonizome/GeneSymbolAndIDS_2017.txt', sep='\t', index_col=0)

def createGeneList(inputDf):


    gene_list = pd.DataFrame(columns=['GeneSym', 'GeneID'])

    gene_list['GeneSym'] = inputDf.index

    for i,index in enumerate(gene_list.index):

        progressPercent = ((i+1)/len(gene_list.index))*100

        sys.stdout.write("Progeres: %d%%  %d Out of %d   \r" % (progressPercent, (i+1), len(gene_list.index)))
        sys.stdout.flush()

        if gene_list.ix[index, 'GeneSym'] in getGeneIDs.index:
            gene_list.ix[index, 'GeneID'] = getGeneIDs.ix[gene_list.ix[index, 'GeneSym'], 'Entrez Gene ID(supplied by NCBI)']


    return(gene_list)

def createAttributeList(inputDF):

    attribute_list = pd.DataFrame(columns=['Attributes'])

    attribute_list['Attributes'] = inputDF.columns

    return(attribute_list)

def createGeneAttributeEdgeList(inputDF, genelist, path, name):

    count = 0

    filenameGMT = name+'_%s.tsv'% str(datetime.date.today())[0:7].replace('-', '_')

    if os.path.isfile(path+filenameGMT):
        os.remove(path+filenameGMT)

    lst = ['Attribute', 'Gene', 'GeneID', 'Weight']
    lst = ['{0}\t'.format(elem) for elem in lst]
    lst.insert(len(lst), '\n')

    with open(path+filenameGMT, 'a') as the_file:
        the_file.writelines(lst)

    # df = pd.DataFrame(columns=['Attribute', 'Gene', 'GeneID', 'Weight'])
    temp = pd.DataFrame(columns=['Attribute', 'Gene', 'GeneID', 'Weight'])

    for i,col in enumerate(inputDF.columns):

        progressPercent = ((i+1)/len(inputDF.columns))*100

        sys.stdout.write("Progeres: %d%%  %d Out of %d   \r" % (progressPercent, (i+1), len(inputDF.columns)))
        sys.stdout.flush()

        temp['Gene'] = inputDF[col].index
        temp['Weight'] = inputDF[col].values.tolist()
        temp['Attribute'] = [col]*len(temp['Gene'])
        temp['GeneID'] = genelist['GeneID']

        with open(path+filenameGMT, 'a') as the_file:
                 temp.to_csv(the_file, header=False, index=False, sep='\t')

        count += temp[temp['Weight'] != 0].shape[0]

        # for index in temp.index:
        #     lst = [temp.ix[index, 'Attribute'], temp.ix[index, 'Gene'], str(temp.ix[index, 'GeneID']), temp.ix[index, 'Weight']]
        #     lst = ['{0}\t'.format(elem) for elem in lst]
        #     lst.insert(len(lst), '\n')
        #
        #     with open(path+filenameGMT, 'a') as the_file:
        #         the_file.writelines(lst)
    print('\n\n The number of statisticaly relevent gene-attribute associations is: %d' %count)


def createBinaryMatix(inputDF, ppi=False):

    if ppi:

        genes = list(set(inputDF.iloc[:,0].unique().tolist()+inputDF.iloc[:,1].unique().tolist()))

        matrix = pd.DataFrame(index=genes, columns=genes, data=0)

        for i, gene in enumerate(genes):

            progressPercent = ((i+1)/len(genes))*100

            sys.stdout.write("Progeres: %d%%  %d Out of %d   \r" % (progressPercent, (i+1), len(genes)))
            sys.stdout.flush()

            lst = inputDF[inputDF.iloc[:,0] == gene].iloc[:,1].tolist()
            lst += inputDF[inputDF.iloc[:,1] == gene].iloc[:,0].tolist()
            lst = set(lst)
            lst.discard(gene)
            lst = list(lst)

            matrix.ix[gene, lst] = 1

        return(matrix)

    else:
        genes = list(set(inputDF.iloc[:,0].unique().tolist()))

        attributes = list(set(inputDF.iloc[:,1].unique().tolist()))

        matrix = pd.DataFrame(index=genes, columns=attributes, data=0)

        for i, gene in enumerate(genes):

            progressPercent = ((i+1)/len(genes))*100

            sys.stdout.write("Progeres: %d%%  %d Out of %d   \r" % (progressPercent, (i+1), len(genes)))
            sys.stdout.flush()

            lst = inputDF[inputDF.iloc[:,0] == gene].iloc[:,1].tolist()
            lst = set(lst)
            lst.discard(gene)
            lst = list(lst)

            matrix.ix[gene, lst] = 1

        return(matrix)
