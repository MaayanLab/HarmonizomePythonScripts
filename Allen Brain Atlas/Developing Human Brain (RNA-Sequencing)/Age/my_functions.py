import sys, datetime, os
import numpy as np
import pandas as pd
import scipy.stats.mstats as stat
import scipy.spatial.distance as dist
from statsmodels.distributions.empirical_distribution import ECDF
from scipy import stats
import zipfile


def merge(inputDF, axis, method):

    if axis == 'column' or axis == 'Column' and method == 'mean':
        inputDF = inputDF.groupby(inputDF.columns, axis=1).mean()
        return(inputDF)

    elif axis == 'row' or axis == 'Row' and method == 'mean':
        inputDF = inputDF.groupby(level=0, axis=0).mean()
        return(inputDF)

def zscore(inputDF, axis):
    if axis == 'row':
        inputDF.reset_index(inplace=True)
        genes = inputDF.iloc[:, 0].values.tolist()
        inputDF.drop(inputDF.columns[0], axis=1, inplace=True)
        # newMatrloc = stats.zscore(inputDF, axis=0, ddof=1)
        # newDF = pd.DataFrame(data=newMatrloc, index=inputDF.index, columns=inputDF.columns)

        # for i,index in enumerate(inputDF.index):
        #
        #     progressPercent = ((i+1)/len(inputDF.index))*100
        #
        #     sys.stdout.write("Progress: %d%%  %d Out of %d   \r" % (progressPercent, (i+1), len(inputDF.index)))
        #     sys.stdout.flush()
        #
        #     mean = inputDF.loc[index].mean()
        #     std = inputDF.loc[index].std()
        #     inputDF.loc[index] = inputDF.loc[index].apply(lambda x: ((x-mean)/std))

        #Modefied Z Score
        for i,index in enumerate(inputDF.index):

            progressPercent = ((i+1)/len(inputDF.index))*100

            sys.stdout.write("Progress: %d%%  %d Out of %d   \r" % (progressPercent, (i+1), len(inputDF.index)))
            sys.stdout.flush()

            median_y = np.median(inputDF.iloc[index])
            median_absolute_deviation_y = np.median([np.abs(y - median_y) for y in inputDF.iloc[index]])
            mean_absolute_deviation_y = np.mean([np.abs(y - median_y) for y in inputDF.iloc[index]])

            if median_absolute_deviation_y !=0:
                modified_z_scores = [0.6745 * (y - median_y) / median_absolute_deviation_y
                                     for y in inputDF.iloc[index]]
            else:
                modified_z_scores = [(y - median_y) / ( 1.253314 * mean_absolute_deviation_y )
                                     for y in inputDF.iloc[index]]
            inputDF.iloc[index] = modified_z_scores
            inputDF.index = genes

    elif axis == 'column':
        for i,col in enumerate(inputDF.columns):

            progressPercent = ((i+1)/len(inputDF.columns))*100

            sys.stdout.write("Progress: %d%%  %d Out of %d   \r" % (progressPercent, (i+1), len(inputDF.columns)))
            sys.stdout.flush()

            mean = inputDF[col].mean()
            std = inputDF[col].std()
            inputDF[col] = inputDF[col].apply(lambda x: ((x-mean)/std))

    # return(newDF)


def log2(inputDF):
    DF = inputDF.copy()
    DF = DF.apply(lambda x: np.log2(x+1))
    return(DF)


def quantileNormalize(inputDF):
    df = inputDF.copy()

    attributes = df.columns.values.tolist()

    df.columns = np.arange(0, len(attributes))

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

    df.columns = attributes
    return df

mappingDF = pd.read_csv('/Users/moshesilverstein/Documents/Harmonizome/mappingFile_2017.txt', sep='\t', header=None, index_col=0)

def mapgenesymbols(inputDF):

    inputDF.reset_index(inplace=True)


    lst1 = []

    for i, index in enumerate(inputDF.index):

        progressPercent = ((i+1)/len(inputDF.index))*100

        sys.stdout.write("Progeres: %d%%  %d Out of %d   \r" % (progressPercent, (i+1), len(inputDF.index)))
        sys.stdout.flush()

        if inputDF.loc[index, inputDF.columns[0]] in mappingDF.index:
            lst1.append(mappingDF.loc[inputDF.loc[index, inputDF.columns[0]], 1])
        else:
            lst1.append(np.nan)


    inputDF[inputDF.columns[0]] = lst1


    inputDF.dropna(inplace=True, subset=[inputDF.columns[0]])
    inputDF.set_index(inputDF.columns[0], inplace=True)

def createTertiaryMatrix(inputDF):

    # percent = 0.10

    df = inputDF.copy()

    # def mapter2(x):
    #
    #     if x > up:
    #         return 1
    #     elif x < down:
    #         return -1
    #     else:
    #         return 0

    # def mapter(x):
    #
    #     if x > up:
    #         return 1
    #     elif x < down:
    #         return -1
    #     else:
    #         return 0

    # def mapter2(x):
    #     if x >= 0.95:
    #         return x
    #     elif x <= -0.95:
    #         return x
    #     else:
    #         return 0

    def mapter(x):
        if x >= 0.95:
            return 1
        elif x <= -0.95:
            return -1
        else:
            return 0

    # for i,col in enumerate(inputDF.columns):
    #
    #     progressPercent = ((i+1)/len(inputDF.columns))*100
    #
    #     sys.stdout.write("Progeres: %d%%  %d Out of %d   \r" % (progressPercent, (i+1), len(inputDF.columns)))
    #     sys.stdout.flush()
    #
    #
    #     values = df[col].values.flatten().tolist()
    #
    #     values.sort(reverse=True)
    #
    #     up = values[int(np.floor(percent*len(values)))]
    #
    #     values = values[::-1]
    #
    #     down = values[int(np.floor(percent*len(values)))]
    #
    #     df[col] = df[col].apply(mapter2)



    # values = df.values.flatten().tolist()

    # values.sort(reverse=True)

    # up = values[int(np.floor(percent*len(values)))]

    # values = values[::-1]

    # down = values[int(np.floor(percent*len(values)))]

    df = df.applymap(mapter)

    # df = df.applymap(mapter2)


    #
    #     lst = []
    #     lst = inputDF[col].values.tolist()
    #     lst.sort()
    #
    #     lst = lst[::-1]
    #     # upCutoff = lst[int(0.1*len(lst))]
    #     upCutoff = lst[499]
    #     index = inputDF[inputDF[col] > upCutoff].index
    #     df.loc[index, col] = 1
    #
    #     zeroIndex = inputDF[inputDF[col] < upCutoff].index
    #     df.loc[zeroIndex, col] = 0
    #
    #     lst = lst[::-1]
    #     # downCutoff = lst[int(0.1*len(lst))]
    #     downCutoff = lst[499]
    #     index = inputDF[inputDF[col] < downCutoff].index
    #     df.loc[index, col] = -1


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
    similarity_df.index.name = ''
    similarity_df.columns.name = ''
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

        if gene_list.loc[index, 'GeneSym'] in getGeneIDs.index:
            gene_list.loc[index, 'GeneID'] = getGeneIDs.loc[gene_list.loc[index, 'GeneSym'], 'Entrez Gene ID(supplied by NCBI)']


    return(gene_list)

def createAttributeList(inputDF, metaData=pd.DataFrame()):

    if not metaData.empty:

        cols = metaData.columns.tolist()

        cols.insert(0, 'Attributes')

        attribute_list = pd.DataFrame(columns=cols)

        attribute_list['Attributes'] = inputDF.columns

        attribute_list.set_index('Attributes', inplace=True)

        for i,attribute in enumerate(attribute_list.index):

            progressPercent = ((i+1)/len(attribute_list.index))*100

            sys.stdout.write("Progeres: %d%%  %d Out of %d   \r" % (progressPercent, (i+1), len(attribute_list.index)))
            sys.stdout.flush()

            for col in attribute_list.columns:
                attribute_list.loc[attribute, col] = metaData.loc[attribute, col]

    else:
        attribute_list = pd.DataFrame(columns=['Attributes'])
        attribute_list['Attributes'] = inputDF.columns
        attribute_list.set_index('Attributes', inplace=True)

    return(attribute_list)

def createGeneAttributeEdgeList(inputDF, attributelist, genelist, path, name):

    count = 0

    filenameGMT = name+'_%s.tsv'% str(datetime.date.today())[0:7].replace('-', '_')

    if os.path.isfile(path+filenameGMT):
        os.remove(path+filenameGMT)


    col = np.append(genelist.columns.values, attributelist.columns.values)
    col = col.flatten().tolist()
    col.insert(2, 'Attribute')
    col.append('Weight')

    temp = pd.DataFrame(columns=col)
    # col = ['Attribute', 'Gene', 'GeneID', 'Weight']
    col = ['{0}\t'.format(elem) for elem in col]

    col.insert(len(col), '\n')

    with open(path+filenameGMT, 'a') as the_file:
        the_file.writelines(col)

    # df = pd.DataFrame(columns=['Attribute', 'Gene', 'GeneID', 'Weight'])

    for i,col in enumerate(inputDF.columns):

        progressPercent = ((i+1)/len(inputDF.columns))*100

        sys.stdout.write("Progeres: %d%%  %d Out of %d   \r" % (progressPercent, (i+1), len(inputDF.columns)))
        sys.stdout.flush()

        temp['GeneSym'] = inputDF[col].index
        temp['GeneID'] = genelist['GeneID']
        temp['Attribute'] = [col]*len(temp['GeneSym'])
        for col2 in attributelist.columns:
            temp[col2] = [attributelist.loc[col, col2]]*len(temp['GeneSym'])
        temp['Weight'] = inputDF[col].values.tolist()

        with open(path+filenameGMT, 'a') as the_file:
                 temp.to_csv(the_file, header=False, index=False, sep='\t')

        count += temp[temp['Weight'] >= 0.95].shape[0]
        count += temp[temp['Weight'] <= -0.95].shape[0]

        # for index in temp.index:
        #     lst = [temp.loc[index, 'Attribute'], temp.loc[index, 'Gene'], str(temp.loc[index, 'GeneID']), temp.loc[index, 'Weight']]
        #     lst = ['{0}\t'.format(elem) for elem in lst]
        #     lst.insert(len(lst), '\n')
        #
        #     with open(path+filenameGMT, 'a') as the_file:
        #         the_file.writelines(lst)
    print('\n\n The number of statisticaly relevent gene-attribute associations is: %d' %count)


def createBinaryMatrix(inputDF, ppi=False):

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

            matrloc.loc[gene, lst] = 1

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

            matrix.loc[gene, lst] = 1

        return(matrix)

def createStandardizedMatrix(inputDF):
    df = inputDF.copy()

    # def mapter2(x):
    #     if x >= 0.95:
    #         return x
    #     elif x <= -0.95:
    #         return x
    #     else:
    #         return 0

    for i,index in enumerate(inputDF.index):

        progressPercent = ((i+1)/len(inputDF.index))*100

        sys.stdout.write("Progeres: %d%%  %d Out of %d   \r" % (progressPercent, (i+1), len(inputDF.index)))
        sys.stdout.flush()

        ourECDF = ECDF(df.loc[index])

        df.loc[index] = ourECDF(df.loc[index])

        mean = df.loc[index].values.flatten().mean()

        df.loc[index] = df.loc[index].apply(lambda x: 2*(x-mean))

        # df.loc[index] = df.loc[index].apply(mapter2)

    values = df.values.flatten()
    ourECDF = ECDF(values)

    ourECDF = ourECDF(values).reshape((len(df.index), len(df.columns)))

    newDF = pd.DataFrame(data = ourECDF, index=df.index, columns= df.columns)

    newDF.index.name = ''

    mean = newDF.values.flatten().mean()

    newDF = newDF.applymap(lambda x: 2*(x-mean))

    # for i,col in enumerate(inputDF.columns):
    #
    #     progressPercent = ((i+1)/len(inputDF.columns))*100
    #
    #     sys.stdout.write("Progeres: %d%%  %d Out of %d   \r" % (progressPercent, (i+1), len(inputDF.columns)))
    #     sys.stdout.flush()
    #
    #     positiveAssociation = np.abs(df[df[col] > 0][col].values.tolist())
    #     positiveAssociationIndex = df[df[col] > 0][col].index
    #     positiveECDF = ECDF(positiveAssociation)
    #
    #     for j,value in enumerate(positiveAssociation):
    #         df.loc[positiveAssociationIndex[j], col] = positiveECDF(value)
    #
    #     negativeAssociation = np.abs(df[df[col] < 0][col].values.tolist())
    #     negativeAssociationIndex = df[df[col] < 0][col].index
    #     negativeECDF = ECDF(negativeAssociation)
    #
    #     for k,value in enumerate(negativeAssociation):
    #         df.loc[negativeAssociationIndex[k], col] = -negativeECDF(value)

    return(newDF)

def removeAndImpute(inputDF):

    df = inputDF.copy()

    #mean = np.mean(df.values.flatten())

    df.reset_index(inplace=True)

    df.replace(0.0, np.nan, inplace=True)

    df.dropna(thresh=int(0.05*df.shape[0]), axis=1, inplace=True)
    df.dropna(thresh=int(0.05*df.shape[1]), axis=0, inplace=True)

    genes = df.iloc[:, 0].values.tolist()

    df.drop(df.columns[0], axis=1, inplace=True)

    #df.replace(np.nan, mean, inplace=True)

    df = df.T.fillna(df.mean(axis=1)).T

    df.index = genes

    return(df)

# def writeToZipFile(inputDF, path, name):
#
#     filename= name+'_%s'% str(datetime.date.today())[0:7].replace('-', '_')
#
#     if os.path.isfile(path+filename):
#         os.remove(path+filename)
#
#
#     for i, index in enumerate(inputDF.index):
#
#             progressPercent = ((i+1)/len(inputDF.columns))*100
#
#             sys.stdout.write("Progeres: %d%%  %d Out of %d   \r" % (progressPercent, (i+1), len(inputDF.columns)))
#             sys.stdout.flush()
#
#
#             lst = inputDF.iloc[i, :].values.tolist()
#             lst.insert(0, index)
#             lst = ['{0}\t'.format(elem) for elem in lst] # add tabs between terms in the lst
#             lst.insert(len(lst), '\n') # add a newline char at the end of each lst
#
#             with open(path+filename, 'a') as the_file:
#                 the_file.writelines(lst)
#
#     with ZipFile(path+filename+'.zip', 'w') as myzip:
#         myzip.write(path+filename)
