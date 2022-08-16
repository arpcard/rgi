import os, json, argparse
from argparse import RawTextHelpFormatter
from app.settings import APP_NAME, SOFTWARE_VERSION

def check_for_all_classifications(classtype, class_dict):
    dc = 0
    rm = 0
    gf = 0
    for key in class_dict:
        if class_dict[key][classtype] == "Drug Class":
            dc += 1
        if class_dict[key][classtype] == "Resistance Mechanism":
            rm += 1
        if class_dict[key][classtype] == "AMR Gene Family":
            gf += 1
    if dc < 1:
        print("Error: could not find a Drug Class classification")
    if rm < 1:
        print("Error: could not find a Resistance Mechanism classification")
        no_rm[rgi_data[orf][ordered[0]]["type_match"]].append(ordered[0])
        print(no_rm)
    if gf < 1:
        print("Error: could not find an AMR Gene Family classification")

def read_file(j):
    with open(os.path.join("{}".format(j))) as f:
        j = json.load(f)
    return j

def identify_snps(j):
    snps = {}
    for orf, hsp in j.items():
        if isinstance(hsp, dict):
            try:
                snps[orf] = {}
                for key in hsp:
                    if hsp[key]['model_id'] in snps[orf]:
                        snps[orf][hsp[key]['model_id']].append(hsp[key]['snp']['original']+str(hsp[key]['snp']['position'])+hsp[key]['snp']['change']+':'+hsp[key]['model_name'])
                    else:
                        snps[orf][hsp[key]['model_id']] = [hsp[key]['snp']['original']+str(hsp[key]['snp']['position'])+hsp[key]['snp']['change']+':'+hsp[key]['model_name']]
            except KeyError:
                snps[orf] = {}
    return snps

def main(j):
    criteria = ['Perfect', 'Strict', 'Loose']
    dc = {'Perfect': {}, 'Strict': {}, 'Loose': {}}
    rm = {'Perfect': {}, 'Strict': {}, 'Loose': {}}
    gf = {'Perfect': {}, 'Strict': {}, 'Loose': {}}
    genes = {'Perfect': {}, 'Strict': {}, 'Loose': {}}
    no_rm = {'Perfect': [], 'Strict': [], 'Loose': []}

    for orf, hsp in j.items():
        if isinstance(hsp, dict):
            best_hsp = max(hsp.keys(), key=(lambda key: hsp[key]['bit_score']))
            if j[orf][best_hsp]['evalue'] > 1e-10:
                continue
            if 'ARO_category' in j[orf][best_hsp]:
                check_for_all_classifications('category_aro_class_name', j[orf][best_hsp]['ARO_category'])
                for key in j[orf][best_hsp]['ARO_category']:
                    for c in criteria:
                        if j[orf][best_hsp]['type_match'] == c:

                            if j[orf][best_hsp]['model_name'] not in genes[c]:
                                genes[c][j[orf][best_hsp]['model_name']] = []
                            for i in genes[c]:
                                if j[orf][best_hsp]['model_name'] == i:
                                    genes[c][i].append({orf: best_hsp})

                            if "category_aro_class_name" in j[orf][best_hsp]["ARO_category"][key]:
                                if j[orf][best_hsp]["ARO_category"][key]["category_aro_class_name"] == "Drug Class":
                                    if j[orf][best_hsp]["ARO_category"][key]["category_aro_name"] not in dc[c]:
                                        dc[c][j[orf][best_hsp]["ARO_category"][key]["category_aro_name"]] = []
                                    for i in dc[c]:
                                        if j[orf][best_hsp]['ARO_category'][key]['category_aro_name'] == i:
                                            dc[c][i].append({orf: best_hsp})

                                elif j[orf][best_hsp]["ARO_category"][key]["category_aro_class_name"] == "Resistance Mechanism":
                                    if j[orf][best_hsp]["ARO_category"][key]["category_aro_name"] not in rm[c]:
                                        rm[c][j[orf][best_hsp]["ARO_category"][key]["category_aro_name"]] = []
                                    for i in rm[c]:
                                        if j[orf][best_hsp]["ARO_category"][key]["category_aro_name"] == i:
                                            rm[c][i].append({orf: best_hsp})
                                    if not not no_rm:
                                        if len(no_rm[j[orf][best_hsp]["type_match"]]) > 0:
                                            rm[j[orf][best_hsp]["type_match"]]["unclassified"] = no_rm[j[orf][best_hsp]["type_match"]]

                                elif j[orf][best_hsp]["ARO_category"][key]["category_aro_class_name"] == "AMR Gene Family":
                                    if j[orf][best_hsp]["ARO_category"][key]["category_aro_name"] not in gf[c]:
                                        gf[c][j[orf][best_hsp]["ARO_category"][key]["category_aro_name"]] = []
                                    for i in gf[c]:
                                        if j[orf][best_hsp]["ARO_category"][key]["category_aro_name"] == i:
                                            gf[c][i].append({orf: best_hsp})
    return dc, rm, gf, genes

def make_json(m,j,f,t,s):
    finalgene = {'name': 'ARO', 'children': []}
    finaldrugclass = {'name': 'Drug Classes', 'children': []}
    finalresistmech = {'name': 'Resistance Mechanisms', 'children': []}
    finalgenefam = {'name': 'AMR Gene Family', 'children': []}
    counting_dict = {'Perfect': 0, 'Strict': 0, 'Loose': 0}

    for dc in m:
        #print(dc)
        for crit in dc:
            if bool(dc[crit]) is False:
                if not any(crit in item.values() for item in finaldrugclass['children']):
                    finaldrugclass['children'].append({'name': crit, 'children': []})
                    finalresistmech['children'].append({'name': crit, 'children': []})
                    finalgenefam['children'].append({'name': crit, 'children': []})
                    finalgene['children'].append({'name': crit, 'children': []})
                else:
                    pass
            elif crit == 'Loose' and f is False:
                cdict = {'name': crit, 'children': []}
                if dc == m[0]:
                    finaldrugclass['children'].append(cdict)
                elif dc == m[1]:
                    finalresistmech['children'].append(cdict)
                elif dc == m[2]:
                    finalgenefam['children'].append(cdict)
                elif dc == m[3]:
                    finalgene['children'].append(cdict)
                else:
                    print("Something has gone horribly wrong")
                continue
            else:
                cdict = {'name': crit, 'children': []}
                for i in dc[crit]:
                    temp = {}
                    temp = {'name': i, 'children': []}
                    innerdict = {}
                    for hit in dc[crit][i]:
                        #print(hit)
                        for orf in j:
                            for hsp in j[orf]:
                                if orf in hit.keys() and hsp in hit.values() and j[orf][hsp]['type_match'] == crit:
                                    innerdict[j[orf][hsp]['model_id']+'__gene_header__'+orf] = []
                                    ndict = {}
                                    ndict['ORF_ID'] = orf.split()[0]
                                    if t == 'protein':
                                        ndict['CONTIG'] = ''
                                        ndict['START'] = ''
                                        ndict['STOP'] = ''
                                        ndict['orf_orientation'] = ''
                                        ndict['Predicted_DNA'] = ''
                                    else:
                                        ndict['CONTIG'] = orf.split()[8]
                                        ndict['START'] = j[orf][hsp]['orf_start']
                                        ndict['STOP'] = j[orf][hsp]['orf_end']
                                        ndict['orf_orientation'] = j[orf][hsp]['orf_strand']
                                        ndict['Predicted_DNA'] = j[orf][hsp]['orf_dna_sequence']
                                    ndict['CUT_OFF'] = j[orf][hsp]['type_match']
                                    ndict['Best_Hit_evalue'] = j[orf][hsp]['evalue']
                                    ndict['ARO_name'] = j[orf][hsp]['model_name']
                                    ndict['Best_Identities'] = j[orf][hsp]['perc_identity']
                                    ndict['ARO'] = j[orf][hsp]['ARO_accession']
                                    ndict['ARO_name'] = j[orf][hsp]['model_name']
                                    ndict['model_type'] = j[orf][hsp]['model_type']
                                    ndict['Best_Hit_ARO_classes'] = {'Drug Class': [], 'Resistance Mechanism': [],
                                                                    'AMR Gene Family': [], 'Antibiotic': [],
                                                                    'Adjuvant': [], 'Efflux Component': [],
                                                                    'Efflux Regulator': [], 'Antibiotic+Adjuvant': []}
                                    ndict['pass_bitscore'] = j[orf][hsp]['pass_bitscore']
                                    ndict['bit_score'] = j[orf][hsp]['bit_score']
                                    ndict['Predicted_Protein'] = j[orf][hsp]['orf_prot_sequence']
                                    ndict['CARD_Protein_Sequence'] = j[orf][hsp]['sequence_from_broadstreet']
                                    ndict['match'] = j[orf][hsp]['match']
                                    ndict['query'] = j[orf][hsp]['query']
                                    ndict['sbjct'] = j[orf][hsp]['sequence_from_db']
                                    ndict['LABEL'] = orf
                                    ndict['ID'] = list(hit.values())[0]
                                    ndict['Model_ID'] = j[orf][hsp]['model_id']

                                    if orf in s and j[orf][hsp]['model_type'] != 'protein homolog model':
                                        for i in s[orf]:
                                            if i == j[orf][hsp]['model_id']:
                                                ndict['SNPs_in_Best_Hit_ARO'] = ', '.join(x.split(':')[0] for x in s[orf][i])
                                            else:
                                                ndict['Other_SNPs'] = ', '.join(x+':'+i for x in s[orf][i])
                                    else:
                                        ndict['SNPs_in_Best_Hit_ARO'] = "n/a"
                                        ndict['Other_SNPs'] = "n/a"

                                    if 'SNPs_in_Best_Hit_ARO' not in ndict:
                                        ndict['SNPs_in_Best_Hit_ARO'] = "n/a"
                                    if 'Other_SNPs' not in ndict:
                                        ndict['Other_SNPs'] = "n/a"

                                    for key in j[orf][hsp]['ARO_category']:
                                        ndict['Best_Hit_ARO_classes'][j[orf][hsp]['ARO_category'][key]['category_aro_class_name']].append(j[orf][hsp]['ARO_category'][key]['category_aro_name'])

                                    innerdict[j[orf][hsp]['model_id']+'__gene_header__'+orf].append(ndict)
                                    #if i == 'cephalosporin':
                                    #    print(innerdict,'\n')
                    else:
                        temp['children'].append(innerdict)
                        cdict['children'].append(temp)
                else:
                    if dc == m[0]:
                        finaldrugclass['children'].append(cdict)
                    elif dc == m[1]:
                        finalresistmech['children'].append(cdict)
                    elif dc == m[2]:
                        finalgenefam['children'].append(cdict)
                    elif dc == m[3]:
                        finalgene['children'].append(cdict)
                    else:
                        print('Something has gone horribly wrong')
                        #print(cdict,'\n')
                                    # if dc[crit][i].index(hit) == 0:
                                    #     cdict['children'].append(temp)
                                    #     print(dc[crit][i])
                                    # else:
                                    #     pass

    else:
        return finaldrugclass, finalresistmech, finalgenefam, finalgene

def create_parser():
    parser = argparse.ArgumentParser(prog="rgi parser",description='{} - {} - Parser \n\nCreates categorical .json files RGI wheel visualization. An input .json file containing the RGI results must be input.'.format(APP_NAME,SOFTWARE_VERSION), formatter_class=RawTextHelpFormatter)
    parser.add_argument('-i', '--input', dest="input", required=True, help="RGI results in a .json file")
    parser.add_argument('-o', '--output', dest="output", default="RGIResults", help="Name/identifier for the output categorical .json files")
    parser.add_argument('--include_loose', dest="loose", action='store_true', help="Include loose hits in addition to strict and perfect hits")
    parser.add_argument('-t', '--type', dest="type", default='contig', help="type of input sequence: contig, protein or read")
    return parser

def write_output(results, counter, output):
    with open(os.path.join(os.getcwd(), output + '-tree-dc.json'), 'w') as o1:
        json.dump(results[0], o1)

    with open(os.path.join(os.getcwd(), output + '-tree-rm.json'), 'w') as o2:
        json.dump(results[1], o2)

    with open(os.path.join(os.getcwd(), output + '-tree-gf.json'), 'w') as o3:
        json.dump(results[2], o3)

    with open(os.path.join(os.getcwd(), output + '-count-hits.json'), 'w') as o4:
        json.dump(counter, o4)

    with open(os.path.join(os.getcwd(), output + '-tree-og.json'), 'w') as o5:
        json.dump(results[3], o5)

def calc_number_of_hits(m,f):
    perfect_list = []
    strict_list = []
    loose_list = []

    for c in m:
        for key,value in c['Perfect'].items():
            for i in value:
                for k,v in i.items():
                    if k not in perfect_list:
                        perfect_list.append(k)
        for key,value in c['Strict'].items():
            for i in value:
                for k,v in i.items():
                    if k not in strict_list:
                        strict_list.append(k)
        if f is True:
            for key,value in c['Loose'].items():
                for i in value:
                    for k,v in i.items():
                        if k not in loose_list:
                            loose_list.append(k)

    count_json = {'Perfect': len(perfect_list), 'Strict': len(strict_list), 'Loose': len(loose_list)}
    return count_json

def api_main(args):
    j = read_file(args.input)
    m = main(j)
    count = calc_number_of_hits(m,args.loose)
    snps = identify_snps(j)
    res = make_json(m,j,args.loose,args.type, snps)
    write_output(res, count, args.output)

def run():
    parser = create_parser()
    args = parser.parse_args()
    api_main(args)

if __name__ == '__main__':
    run()
