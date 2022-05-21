"""
This module works as single file application.
It downloads required into 'phase1/csv_data' folder.
"""

import wget
import os


# required files
files = [{'name':'CRISPR_gene_effect.csv',
          'link':'https://ndownloader.figshare.com/files/34008491'},
         {'name':'CRISPR_gene_dependency.csv',
          'link':'https://ndownloader.figshare.com/files/34008485'},
         {'name':'CCLE_expression.csv',
          'link':'https://ndownloader.figshare.com/files/34008404'},
         {'name':'CCLE_mutations.csv',
          'link':'https://ndownloader.figshare.com/files/34008434'},
         {'name':'sample_info.csv',
          'link':'https://ndownloader.figshare.com/files/34008503'}]


def download_csv_files():
    """
    Make http request for remote file data.
    """

    # create the directory
    try:
        os.makedirs(os.path.join('phase1', 'csv_data'))
    except FileExistsError:
        # directory already exists
        pass

    # download required files
    for f in files:
        path = os.path.join('phase1', 'csv_data', f['name'])
        print(path)
        if os.path.exists(path):
            os.remove(path)  # if exist -> remove it
        wget.download(f['link'], path)
        print('')


if __name__ == '__main__':
    download_csv_files()
