#import libraries
import requests
import json
import re
import subprocess
from Bio import SeqIO
import json


# load the testing data
file = "test_file.fasta"
ext = "fasta"

#Seqkit part
def seqkit_parser(filename):
  seqkit = subprocess.run(("seqkit", "stats", filename, "-a"),
                            capture_output=True,
                            text=True)
  if seqkit.stderr!="":
    print("Something went wrong:", seqkit.stderr)

  else:
    fasta_sk = seqkit.stdout.strip().split('\n')
    fasta_stats_sk = dict(zip(fasta_sk[0].split()[1:], fasta_sk[1].split()[1:]))

    return fasta_stats_sk
  
#parse using seqkit
res = seqkit_parser(file)
print("Output of seqkit_parser: ", res)

#get type of sequences
type_f = list(seqkit_parser(file).values())[1]
print("The type of sequence is: ", type_f)

#Biopython part
def biopython_parser(filename):
  type_f = list(seqkit_parser(filename).values())[1]

  biopy = SeqIO.parse(filename, ext)
  ids= []

  if type_f == "Protein":
    for seq in biopy:
      s_match = re.fullmatch('[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}', seq.description)
      ids.append(s_match)

  elif type_f == "DNA":
    for seq in biopy:
      s_match = re.fullmatch('ENS([A-Z]{0})(E|FM|G|GT|P|R|T)[0-9]{11}|ENS[A-Z]{3,4}(E|FM|G|GT|P|R|T)[0-9]{11}', seq.description)
      ids.append(s_match)

  return ids

#the code breaks at extracting IDs:
res = biopython_parser(file)
print(res) #[None, None]

#Part from HW2_1

#function 1 for Uniprot
def get_uniprot(ids: list):
  accessions = ','.join(ids)
  http_args = {'params': {'accessions': accessions}}
  return requests.get("https://rest.uniprot.org/uniprotkb/accessions", **http_args)

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

#function 1 for Ensembl
def get_ensembl(ids: list):
  server = "https://rest.ensembl.org"
  ext = "/lookup/id"
  id_s = {"ids": ids}
  headers={ "Content-Type" : "application/json", "Accept" : "application/json"}
  return requests.post(server+ext, headers=headers, data=json.dumps(id_s))

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

#PART 2
#function for identifying the source by ID using regex and outputing corresponding parsed result
def regex_fun(ids: list):
  if re.fullmatch('[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}', ids[0]):
    resp = get_uniprot(my_ids)
    res = parse_response_uniprot(resp)

  elif re.fullmatch('ENS([A-Z]{0})(E|FM|G|GT|P|R|T)[0-9]{11}|ENS[A-Z]{3,4}(E|FM|G|GT|P|R|T)[0-9]{11}', ids[0]):
    resp = get_ensembl(my_ids)
    res = parse_response_ensembl(resp)

  return res


