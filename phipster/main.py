import json
import requests
from os import listdir
from os.path import isfile, join
from collections import defaultdict


def download_interactions():
    # Download viral protein - human protein interactions for each virus
    # Save by virus id
    url = 'http://ec2-18-215-31-69.compute-1.amazonaws.com:8000'
    viri = json.load(open('virus.json', 'r'))

    for index, virus in enumerate(viri[580:]):
        if index % 10 == 0:
            print('{0}/1001'.format(index+580))
        virus_interact_url = '{0}/home/ppi/virus-lrfilter/{1}/100?format=json'.format(url, virus['id'])
        virus_interact_json = requests.get(virus_interact_url).json()
        with open('interactions/{0}.json'.format(virus['id']), 'w') as interact_file:
            json.dump(virus_interact_json, interact_file)
    return None


def extract_id_name(json):
    id_name_dict = dict()
    for rec in json:
        id_name_dict[rec["id"]] = rec["name"]
    return id_name_dict


def main():
    # http://ec2-18-215-31-69.compute-1.amazonaws.com:8000/home/hp/?format=json
    hp_id_name = extract_id_name(json.load(open('human_proteins.json', 'r')))
    # http://ec2-18-215-31-69.compute-1.amazonaws.com:8000/home/vp/?format=json
    vp_id_name = extract_id_name(json.load(open('viral_proteins.json', 'r')))
    # http://ec2-18-215-31-69.compute-1.amazonaws.com:8000/home/virus/?format=json
    virus_id_name = extract_id_name(json.load(open('virus.json', 'r')))

    interaction_files = [f for f in listdir('interactions') if isfile(join('interactions', f))]
    gmt = defaultdict(lambda: defaultdict(list))

    # Iterate over viral protein - human protein interaction JSONs
    # Make a nested dictionary: virus[viral protein][gene1, gene2, ..., geneN]
    for interaction_file in interaction_files:
        interaction = json.load(open('interactions/{0}'.format(interaction_file), 'r'))
        virus_id = interaction_file.split('.')[0]
        for vp_interaction in interaction:
            vp = vp_id_name[vp_interaction["viralprotein_id"]]
            hp = hp_id_name[vp_interaction["humanprotein_id"]]
            virus = virus_id_name[int(virus_id)]
            gmt[virus][vp].append(hp)

    # Flatten the dictionary
    gmt_plain = []
    for virus in sorted(gmt.keys()):
        for vp in gmt[virus]:
            if len(gmt[virus][vp]) >= 5:
                term = '{0} {1}\t\t{2}'.format(virus, vp, '\t'.join(sorted(gmt[virus][vp])))
                gmt_plain.append(term)

    with open('phipster.gmt', 'w') as phipster:
        phipster.write('\n'.join(gmt_plain))

    return None


if __name__ == '__main__':
    main()