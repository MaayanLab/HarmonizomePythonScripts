# Creates the dictionaries and converts them to GMT files

import pandas as pd
from collections import defaultdict

df = pd.DataFrame()
print('Running...')

kinome_df = pd.read_excel('http://lincs.hms.harvard.edu/wordpress/wp-content/uploads/2013/11/HMS-LINCS_KinomeScan_Datasets_2018-01-18.xlsx')
kinome_df['sm_hms_id'] = kinome_df['sm_hms_id'].str.replace(r'\D', '') #remove HMSL before the ID

# retrieve the CSV file for the small molecules, 2-183
for csv_name in range(0, 182):  
    hms_id_string = str(kinome_df.iloc[csv_name]['sm_hms_id'])
    # Eliminate the small molecules that are not in the LINCS database (which produce random molecules and kinases)
    if (int(hms_id_string) < 10520): 
        url = 'http://lincs.hms.harvard.edu/db/datasets/20000/results?small+molecules={0}&output_type=.csv'.format(hms_id_string)
        data = pd.read_csv(url)
        data_df = pd.DataFrame(data)
        df = df.append(data_df.loc[data_df['Binding Class'] <= 3]) # Only append the rows with target affinity <= 3

# create dictionary of dictionaries, one corresponding to each target affinity level. Set is used to eliminate duplicates 
levels = {'1': defaultdict(set), '2': defaultdict(set), '3': defaultdict(set)}

# kinases are keys and small molecules are values 
for x in range(0, df.shape[0] - 1):
    kinase_name = df.iloc[x]['HUGO Gene Symbol']
    sm_name = df.iloc[x]['HMSL Small Mol Name']
    target_aff = df.iloc[x]['Binding Class']
    
    levels[str(target_aff)][kinase_name].add(sm_name)
    
fw1 = open('level1.gmt', 'w')
for k, v in levels['1'].items():
    print(str(k) + '\t\t', '\t'.join(v), file = fw1)
fw1.close()

fw2 = open('level2.gmt', 'w')
for k, v in levels['2'].items():
    print(str(k) + '\t\t', '\t'.join(v), file = fw2)
fw2.close()

fw3 = open('level3.gmt', 'w')
for k, v in levels['3'].items():
    print(str(k) + '\t\t', '\t'.join(v), file = fw3)
fw3.close()

print('Done.')
