# The purpose of this script is for finding the cocited pmids in batch (by using PMC api)
# The purpose the codes in "Cochrane_PMCnet_cocited.ipynb" is for demonstration and testing

import pandas as pd
import numpy as np
import glob
import json
from lxml import etree as ET
import requests
import time
import re


def divide_chunks(l, n):
    # looping till length l
    for i in range(0, len(l), n):
        yield l[i:i + n]

# result = list(divide_chunks( target_for_breaking_down , 5000))

df0 = pd.read_csv('/Users/iwishsomeday/Desktop/2019_SummerProject/Cochrane/outfiles/IncludedStudies_AllVersions_withPMID_v2.tsv', sep='\t', dtype=str)
#df0.info()

dates = pd.read_csv('/Users/iwishsomeday/Desktop/2019_SummerProject/Cochrane/outfiles/IncludedStudies_AllVersions_Dates.tsv', sep='\t')
dates = dates[['accession_number','LastSearchYR']]
#dates.info()

df = df0.merge(dates, how='left')
#df.info()

target = pd.read_csv('/Users/iwishsomeday/Desktop/2019_SummerProject/Cochrane/outfiles/IncludedStudies_AllVersions_withPMID_HasSeednNew.csv')
#target.info()

df['target'] = df.accession_number.isin(target.accession_number.tolist())
dft = df.loc[df.target == True]
dft.accession_number.nunique()
#dft.info()

mask_seed = (dft.review_version == 'old') | (dft.review_version == 'both')
df_seed = dft.loc[mask_seed]
#df_seed.info()

df_updated = dft.loc[dft.review_version == 'updated']
#df_updated.info()

SR2seed_df = df_seed.dropna().groupby('accession_number')['pmid'].apply(list).to_frame().reset_index()
SR2seed_df['accession_number'] = SR2seed_df['accession_number'].astype(str)
SR2seed_df = SR2seed_df.rename(columns={'pmid':'seed_pmid'})
#SR2seed_df.info()

SR2updated_df = df_updated.dropna().groupby('accession_number')['pmid'].apply(list).to_frame().reset_index()
SR2updated_df['accession_number'] = SR2updated_df['accession_number'].astype(str)
SR2updated_df = SR2updated_df.rename(columns={'pmid':'updated_pmid'})
#SR2updated_df.info()

review_pairs_df = SR2seed_df.merge(SR2updated_df)
review_pairs_df = review_pairs_df.merge(df_seed[['accession_number','LastSearchYR','pubyr_review_old', 'NumberOfStudies', 'NumberOfStudies_old']].drop_duplicates())
review_pairs_df = review_pairs_df.replace('NONE', np.nan)
review_pairs_df.info()

review_pairs = review_pairs_df.values.tolist()
#print(review_pairs[0])

pout = '/Users/iwishsomeday/Desktop/2019_SummerProject/Cochrane/outfiles/Cochrance_cocited_counts.tsv'
pout_net = '/Users/iwishsomeday/Desktop/2019_SummerProject/Cochrane/outfiles/Cochrance_cocited_network.tsv'

'''
with open(pout, 'a') as fout:
    fout.write('accession_number' + '\t' + 'N_Seed' + '\t' + 'N_SeedHasPMID' + '\t' 'N_SeedNoPMID' + '\t' + 'N_Cocited' + '\t' + 'N_CocitedInTimewindow' + '\t' + 'N_StudiesFoundInCocited' + "\t" + 'N_StudiesNotFoundInCocited' + '\n')
'''
'''
with open(pout_net, 'a') as fout:
    fout.write('accession_number' + '\t' + 'seed_pmid' + '\t' + 'citing_pmid' + '\t' + 'co_cited_pmid' + '\t' + 'included' + '\n')
'''

# half = pd.read_csv('/Volumes/tkh_mypassport/CochraneXML/Cochrance_cocited_counts.tsv', sep='\t')
# half.info()
# print(half.tail(2))
# #
# slice_df = review_pairs_df.loc[review_pairs_df.index > 979]
# print(slice_df.head())
#
# review_pairs = slice_df.values.tolist()
# print(len(review_pairs))


'''
Get the seed ready
'''

for review_pair in review_pairs:
    accession_number = review_pair[0]
    seed = list(set(review_pair[1]))
    updated_pmid = list(set(review_pair[2]))
    lastsearch = review_pair[3]
    pubyr_review_old = review_pair[4]
    N_included_new = review_pair[5]
    N_included_old = int(review_pair[6])
    len_seed = len(seed)
    diff = N_included_old - len_seed

    #print(accession_number, lastsearch, pubyr_review_old, N_included_old, seed, updated_pmid)

    ### Find the articles which cite seeds

    pmids = '&id=' + '&id='.join(seed)
    # pmids = '&id=18613223'
    # print(pmids)

    pmid_citing_seed_pair = []
    pmid_citing_seed = []

    get_pmid = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=pubmed&linkname=pubmed_pubmed_citedin" + pmids + "&retmode=json"
    response = requests.get(get_pmid)
    get_result = response.content.decode("utf-8")
    json_result = json.loads(get_result)
    # print(json_result)
    citpair = json_result['linksets']
    # print(citpair)
    c = 0
    for pair in citpair:
        seed_pmid = pair['ids'][0]

        #citing_pmids = pair['linksetdbs']
        #if len(citing_pmids) > 0:
        try:
            citing_pmids = pair['linksetdbs']
            for citing in citing_pmids:
                citing_pmids = citing['links']
                pmid_citing_seed_pair.append([seed_pmid, citing_pmids])
                for citing_pmid in citing_pmids:
                    if citing_pmid not in pmid_citing_seed:
                        pmid_citing_seed.append(citing_pmid)
        #else:
        except:
            citing_pmids = ['NONE']
            c += 1
            pmid_citing_seed_pair.append([seed_pmid, citing_pmids])
    time.sleep(1)

    ## Number of pmids which cited the seeds
    N_pmid_citing_seed = len(pmid_citing_seed)
    ## Number of seeds not being cited in PMC
    N_seed_NotCiedInPMC = c
    ## Number of seeds being cited in PMC
    N_seed_citedInPMC = len(seed) - c

    pmid_citing_seed_pair2 = [['seed_pmid', 'citing_pmid']]
    for pcsp in pmid_citing_seed_pair:
        seed_pmid = pcsp[0]
        for pid in pcsp[1]:
            citing_pmid = pid
            pmid_citing_seed_pair2.append([seed_pmid, citing_pmid])
    pmid_citing_seed_pair_df = pd.DataFrame(pmid_citing_seed_pair2[1:], columns=pmid_citing_seed_pair2[0])

    #pmid_citing_seed_pair_df.info()
    #pmid_citing_seed_pair_df.head()


    ### Find the co-cited articles

    PMID_chunks = list(divide_chunks(pmid_citing_seed, 300))

    cited_pairs = []

    for chunk in PMID_chunks:
        pmids = '&id=' + '&id='.join(chunk)
        # print(pmids)

        get_pmid = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=pubmed&linkname=pubmed_pubmed_refs" + pmids + "&retmode=json"
        response = requests.get(get_pmid)
        get_result = response.content.decode("utf-8")
        json_result = json.loads(get_result)
        citpair = json_result['linksets']
        # print(citpair)
        for pair in citpair:
            citing_pmid = pair['ids'][0]
            try:
                cited_set = pair['linksetdbs']
                for cs in cited_set:
                    cited_pmids = cs['links']
                    cited_pairs.append([citing_pmid, cited_pmids])

            except:
                cited_pairs.append([citing_pmid, ['NoSeed']])
        time.sleep(1)

    #print(len(cited_pairs))

    cited_pairs2 = [['citing_pmid', 'cocited_pmid']]
    for cp in cited_pairs:
        citing_pmid = cp[0]
        for cocited in cp[1]:
            cited_pairs2.append([citing_pmid, cocited])
    cited_pairs_df = pd.DataFrame(cited_pairs2[1:], columns=cited_pairs2[0])
    #cited_pairs_df.info()

    cited = []

    for c in cited_pairs:
        cited_pmids = c[1]
        for cited_id in cited_pmids:
            cited.append(cited_id)

    cited = list(set(cited))
    co_cited = list(set(cited) - set(seed))
    #print(len(cited))
    #print(len(co_cited))
    #print(co_cited[0:3])

    ## Number of cocited PMID
    N_cocited = len(co_cited)

    ### get pubyr of the co-cited articles

    cocited_chunks = list(divide_chunks(co_cited, 450))
    #print(len(cocited_chunks))

    cocited_inrange = []
    errors = []

    for cocited in cocited_chunks:
        pmid = ','.join(cocited)
        # print(pmid)
        get_pmid = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&id=' + str(
            pmid) + '&version=2.0'

        response = requests.get(get_pmid)
        get_result = response.content
        # print(get_result)

        tree = ET.fromstring(get_result)
        docsum = tree.xpath('.//DocumentSummary')
        for doc in docsum:
            cocited_pmid = doc.attrib['uid']

            # if len(doc.xpath('./PubDate')) < 1:
            #     pass
            # else:
            #     pubdate = doc.xpath('./PubDate')[0].text
            #     pubyr = re.findall(r'[12]\d\d\d', pubdate)[0]
            #     pubyr = int(pubyr)
            #     if int(lastsearch) >= pubyr >= int(pubyr_review_old):
            #         cocited_inrange.append([cocited_pmid, pubyr])

            try:
                pubdate = doc.xpath('./PubDate')[0].text
                pubyr = re.findall(r'[12]\d\d\d', pubdate)[0]
                pubyr = int(pubyr)
                if int(lastsearch) >= pubyr >= int(pubyr_review_old):
                    cocited_inrange.append([cocited_pmid, pubyr])
            except:
                errors.append(cocited_pmid)

        time.sleep(1)

    cocited_inrange_df = pd.DataFrame(cocited_inrange, columns=['cocited_pmid', 'cocited_pubyr'])
    #cocited_inrange_df.info()

    ### Find how many articles being cited in the updated review are in the co-cited network

    cocited_inrange_pmid = []
    for cid in cocited_inrange:
        cocited_inrange_pmid.append(cid[0])

    ## Number of cocited pmid fit the time window
    N_time_fit = len(cocited_inrange_pmid)

    NotFoundInCocited = list(set(updated_pmid) - set(cocited_inrange_pmid))
    FoundInCocited = list(set(updated_pmid) - set(NotFoundInCocited))

    #print(accession_number, N_included_old, len_seed, diff, N_cocited, N_time_fit, len(FoundInCocited), len(NotFoundInCocited))
    '''
    with open(pout, 'a') as fout:
        fout.write(accession_number + '\t' + str(N_included_old) + '\t' + str(len_seed) + '\t' + str(diff) + '\t' + str(N_cocited) + '\t' + str(N_time_fit) + '\t' + str(len(FoundInCocited)) + '\t' + str(len(NotFoundInCocited)) + '\n')
    '''

    ### Merge dadaframe and get the network

    # print(pmid_citing_seed_pair_df.info())
    # print(cited_pairs_df.info())
    # print(cocited_inrange_df.info())

    M1 = pmid_citing_seed_pair_df.merge(cited_pairs_df)
    M2 = M1.merge(cocited_inrange_df).drop_duplicates()
    M2['included'] = M2.cocited_pmid.isin(FoundInCocited)
    M2.insert(loc=0, column='accession_number', value=accession_number)
    M2 = M2[['accession_number', 'seed_pmid', 'citing_pmid', 'cocited_pmid', 'included']]
    M2 = M2.sort_values(['seed_pmid', 'citing_pmid', 'cocited_pmid'])
    #M2.info()

    M2out = M2.values.tolist()
    for mo in M2out:
        '''
        with open(pout_net, 'a') as fout:
            fout.write(mo[0] + '\t' + str(mo[1]) + '\t' + str(mo[2]) + '\t' + str(mo[3]) + '\t' + str(mo[4]) + '\n')
        '''









