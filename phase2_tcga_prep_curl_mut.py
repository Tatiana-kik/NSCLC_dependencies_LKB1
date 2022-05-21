"""
This module works as single file application.
It makes part of Phase2:  prepare request file phase2/request_maf.txt
                          to download MAF files by curl.

After execution of this file you can download files by curl:

  curl --output phase2/mutation_files.tar.gz \
       --remote-name --remote-header-name --request POST \
       --header "Content-Type: application/json" \
       --data @phase2/request_maf.txt "https://api.gdc.cancer.gov/data"
"""

import requests
import json
import pickle
import os

# Phase II: Getting the list of files with mutation data from TCGA-LUAD

# Preparing request to get TCGA-LUAD mutation data files with case ids and
# primary site "lung" from TCGA "https://api.gdc.cancer.gov/files" endpoint
# to download
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
            "value": ["WXS"]
            }
        },
        {
        "op": "in",
        "content":{
            "field": "files.data_format",
            "value": ["MAF"]
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

# A POST is used, so the filter parameters can be passed directly as
# a Dict object.
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
for file in files:
    if file.split('\t')[0] != 'cases.0.case_id' and len(file)>0:
        if file.split('\t')[2].split('.')[2] == 'aliquot_ensemble_masked':
            # Combining the list of files to download
            file_ids.append(file.split('\t')[3].split('\r')[0])

print(file_ids)

print(len(file_ids))
params_files = {
    "ids": file_ids
}

# Saving file ids dictionary to download with CURL
path = os.path.join('phase2', 'request_maf.txt')
with open(path, 'w') as f:
    f.write(json.dumps(params_files))
    print(f'\nFile-ids dictionary saved at {path}.\n')

# prepare directory to extract files after downloading archive by cURL
os.mkdir(os.path.join('phase2', 'TCGA_MAF'))
