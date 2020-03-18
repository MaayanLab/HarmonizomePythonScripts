import json
import requests
from os import listdir
from os.path import isfile, join
from collections import defaultdict, Counter


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


def check_dups(proteins):
    k_low = [k.lower() for k in proteins.keys()]
    if len(k_low) - len(set(k_low)) > 0:
        return [item for item, count in Counter(k_low).items() if count > 1]
    return None


def dup_term(dup, jk):
    return jk[jk.index(dup.lower)]


def main():
    # http://ec2-18-215-31-69.compute-1.amazonaws.com:8000/home/hp/?format=json
    hp_id_name = extract_id_name(json.load(open('human_proteins.json', 'r')))
    # http://ec2-18-215-31-69.compute-1.amazonaws.com:8000/home/vp/?format=json
    vp_id_name = extract_id_name(json.load(open('viral_proteins.json', 'r')))
    # http://ec2-18-215-31-69.compute-1.amazonaws.com:8000/home/virus/?format=json
    virus_id_name = extract_id_name(json.load(open('virus.json', 'r')))

    interaction_files = [f for f in listdir('interactions') if isfile(join('interactions', f))]
    gmt = defaultdict(lambda: defaultdict(set))

    # Iterate over viral protein - human protein interaction JSONs
    # Make a nested dictionary: virus[viral protein][gene1, gene2, ..., geneN]
    for interaction_file in interaction_files:
        interaction = json.load(open('interactions/{0}'.format(interaction_file), 'r'))
        virus_id = interaction_file.split('.')[0]
        for vp_interaction in interaction:
            vp = vp_id_name[vp_interaction["viralprotein_id"]]
            hp = hp_id_name[vp_interaction["humanprotein_id"]]
            virus = virus_id_name[int(virus_id)]
            # Case-insensitive check for duplicate terms. If found, use the first term name
            gk = list(gmt[virus].keys())
            if gk:
                gk_low = [k.lower() for k in gk]
                if vp.lower() in gk_low:
                    vp = gk[gk_low.index(vp.lower())]
            gmt[virus][vp].add(hp)

    # Flatten the dictionary
    gmt_plain = []
    for virus in sorted(gmt.keys()):
        for vp in gmt[virus]:
            if len(gmt[virus][vp]) >= 5:
                term = '{0} {1}\t\t{2}'.format(virus, vp, '\t'.join(sorted(list(gmt[virus][vp]))))
                gmt_plain.append(term)

    with open('phipster.gmt', 'w') as phipster:
        phipster.write('\n'.join(gmt_plain))

    return None


if __name__ == '__main__':
    main()