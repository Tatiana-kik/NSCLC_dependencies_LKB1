"""
This module works as single file application.
It makes preparation for Phase2 - the list of files with expression data
from TCGA-LUAD.

As result of execution this file we will have file (phase2/request_tsv.txt)
for HTTP request.

To download files you can use curl command with parameters:

  curl --output phase2/expression_files.tar.gz \
       --remote-name --remote-header-name --request POST \
       --header "Content-Type: application/json" \
       --data @phase2/request_tsv.txt "https://api.gdc.cancer.gov/data"
"""

import requests
import pickle
import json
import os

# Phase II: Getting the list of files with expression data from TCGA-LUAD

# Preparing request to get TCGA-LUAD expression data files with case ids and
# primary site "lung" from TCGA "https://api.gdc.cancer.gov/files" endpoint to
# download.
fields = [
    "file_name",
    "files.file_id",
    "cases.case_id",
    "cases.disease_type"
    ]

fields = ",".join(fields)

files_endpt = "https://api.gdc.cancer.gov/files"

filters = {
    "op": "and",
    "content":[
        {
        "op": "in",
        "content":{
            "field": "cases.project.primary_site",
            "value": ["Lung"]
            }
        },
        {
        "op": "in",
        "content":{
            "field": "files.experimental_strategy",
            "value": ["RNA-Seq"]
            }
        },
        {
        "op": "in",
        "content":{
            "field": "files.data_format",
            "value": ["TSV"]
            }
        },
        {
        "op": "in",
        "content": {
            "field": "cases.project.project_id",
            "value": ["TCGA-LUAD"]
            }
        }
    ]
}

params = {
    "filters": filters,
    "fields": fields,
    "format": "TSV",
    "size": "3000"
}

# Requesting the list of files and case ids
response = requests.post(files_endpt,
                         headers={"Content-Type": "application/json"},
                         json=params)

print(response.content.decode("utf-8"))
files = response.content.decode("utf-8").split('\n')
file_ids = []
file_dict = {}
disease = []

for file in files:
    if file.split('\t')[0] != 'cases.0.case_id' and len(file)>0:
        if file.split('\t')[2].split('.')[2] == 'augmented_star_gene_counts':
            # Combining the list of files to download
            file_ids.append(file.split('\t')[3].split('\r')[0])
            # Combining the dictionary with case ids and file names
            file_dict[file.split('\t')[2]] = file.split('\t')[0]

print(file_ids)
print(file_dict)
print(len(file_ids))

# create the directory
try:
    os.makedirs('phase2')
except FileExistsError:
    # directory already exists
    pass

# Dumping case ids dictionary for further data split
path = os.path.join('phase2', 'file_dict.p')
pickle.dump(file_dict, open(path, 'wb'))
print(f'\nCase-ids dictionary dumped at {path}.')

params_files = {
    "ids": file_ids
}

# Saving file ids dictionary to download with CURL
path = os.path.join('phase2', 'request_tsv.txt')
with open(path, 'w') as f:
    f.write(json.dumps(params_files))
    print(f'\nFile-ids dictionary saved at {path}.\n')

# prepare directory to extract files after downloading archive by cURL
os.mkdir(os.path.join('phase2', 'TCGA'))
