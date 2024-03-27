#import libraries
import requests
import json
import re


# PART 1
#function 1 for Uniprot
def get_uniprot(ids: list):
  accessions = ','.join(ids)
  http_args = {'params': {'accessions': accessions}}
  return requests.get("https://rest.uniprot.org/uniprotkb/accessions", **http_args)

#testing for function 1 (Uniprot)
my_ids = ['P11473', 'P13053']
resp = get_uniprot(my_ids)
print("Response code for get_uniprot function: ", resp)

#function 2 for Uniprot
def parse_response_uniprot(resp: dict):
    resp = resp.json()
    resp = resp["results"]
    output = {}
    for val in resp:
        acc = val['primaryAccession']
        species = val['organism']['scientificName']
        gene = val['genes']
        seq = val['sequence']
        output[acc] = {'organism':species, 'geneInfo':gene, 'sequenceInfo':seq, 'type':'protein', "========================": "========================"}

    return output

#testing of function 2 for Uniprot
res = parse_response_uniprot(resp)
print("Printing output for provided Uniprot IDs : ", res)

#function 1 for Ensembl
def get_ensembl(ids: list):
  server = "https://rest.ensembl.org"
  ext = "/lookup/id"
  id_s = {"ids": ids}
  headers={ "Content-Type" : "application/json", "Accept" : "application/json"}
  return requests.post(server+ext, headers=headers, data=json.dumps(id_s))

#testing for function 1 (Ensembl)
my_ids = ['ENSMUSG00000041147', 'ENSG00000139618']
resp = get_ensembl(my_ids)
print("Response code for get_ensembl function: ", resp)

# function 2 for Ensembl
def parse_response_ensembl(resp: dict):
    resp = resp.json()
    output = {}
    for key, val in resp.items():
        species = val["species"]
        gene = val["display_name"]
        description = val["description"]
        feature = val["biotype"]
        transcript = val["canonical_transcript"]
        output[key] = {"organism": species, "gene": gene, "description": description, "feature": feature, "transcript name": transcript, "============": "============"}

    return output

#testing of function 2 (Ensembl)
res = parse_response_ensembl(resp)
print("Printing output for provided Ensembl IDs : ", res)

#PART 2
#function for identifying the source by ID using regex and outputing corresponding parsed result
def regex_fun(ids: list):
  if re.fullmatch('[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}', ids[0]):
    resp = get_uniprot(my_ids)
    res = parse_response_uniprot(resp)

  elif re.fullmatch('ENS[A-Z]{0}[0-9]{11}|ENS[A-Z]{3,4}[0-9]{11}', ids[0]):
    resp = get_ensembl(my_ids)
    res = parse_response_ensembl(resp)

  return res

#testing the regex_fun function (with Ensembl IDs)
my_ids = ['ENSMUSG00000041147', 'ENSG00000139618']
result = regex_fun(my_ids)
print("Ensembl IDs were provided. Printing output: ", result)

#testing the regex_fun function (with Uniprot IDs)
my_ids = ['P11473', 'P13053']
result = regex_fun(my_ids)
print("Uniprot IDs were provided. Printing output: ", result)