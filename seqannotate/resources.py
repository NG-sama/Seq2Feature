import os
import subprocess
import sys
from datetime import date
from importlib.resources import files
from tempfile import NamedTemporaryFile
from io import StringIO

import pandas as pd
import yaml
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

"""
Global varibales
"""
rawdatabase_loc = "/workspaces/universal/Seq2Feature/seqannotate/data/data/BLAST_dbs" #Tweak location of rawdatabase. Primary database can be found in pLAnnotate repo [https://github.com/mmcguffi/pLannotate.git]


ROOT_DIR = os.path.dirname(os.path.abspath(__file__))

valid_genbank_exts = [".gbk", ".gb", ".gbf", ".gbff"]
valid_fasta_exts = [".fa", ".fasta"]
MAX_PLAS_SIZE = 50000

DF_COLS = [
    "sseqid",
    "qstart",
    "qend",
    "sstart",
    "send",
    "sframe",
    "score",
    "evalue",
    "qseq",
    "length",
    "slen",
    "pident",
    "qlen",
    "db",
    "Feature",
    "Description",
    "Type",
    "priority",
    "percmatch",
    "abs percmatch",
    "pi_permatch",
    "wiggle",
    "wstart",
    "wend",
    "kind",
    "qstart_dup",
    "qend_dup",
    "fragment",
]


def get_resource(group, name):
    return str(files(__package__) / f"data/{group}/{name}")

"""
def get_image(name):
    return get_resource("images", name)


def get_template(name):
    return get_resource("templates", name)


def get_example_fastas():
    return get_resource("fastas", "")

"""
def get_yaml_path():
    return get_resource("data", "databases.yml")

def get_details(name):
    return get_resource("data", name)

def get_yaml(yaml_file_loc):
    # file_name = get_resource("data", "databases.yml")
    with open(yaml_file_loc, "r") as f:
        dbs = yaml.load(f, Loader=yaml.SafeLoader)

    # collapes list
    for db in dbs.keys():

        blast_database_loc = dbs[db]["location"]
        if blast_database_loc == "Default":
            blast_database_loc = str(rawdatabase_loc) # Refer to global variable for location of rawdatabase
        try:
            parameters = " ".join(dbs[db]["parameters"])
        except KeyError:
            parameters = ""
        dbs[db]["parameters"] = parameters

        db_loc = os.path.join(blast_database_loc, db)
        # dbs[db]['name'] = db
        dbs[db]["db_loc"] = db_loc
        # db_list.append(dbs[db])
    return dbs
def parse_fasta_string(fasta_string):
    """
    Parses a FASTA string and returns a list of SeqRecord objects.
    """
    fasta_io = StringIO(fasta_string)
    return list(SeqIO.parse(fasta_io, "fasta"))

def process_row(fasta_data):
    """
    Processes a single FASTA data string and returns a dictionary with record attributes.
    """
    records = parse_fasta_string(fasta_data)
    if records:
        record = records[0]  # Assuming each row contains a single FASTA record
        return {
            'id': record.id,
            'name': record.name,
            'description': record.description,
            'dbxrefs': ','.join(record.dbxrefs) if record.dbxrefs else pd.NA,
            'seq': str(record.seq)
        }
    return {
        'id': pd.NA,
        'name': pd.NA,
        'description': pd.NA,
        'dbxrefs': pd.NA,
        'seq': pd.NA
    }

def get_csv_to_gene(pd_file, seq_col="Full Gene Seq"):
    """
    Processes a dataframe containing rows of FASTA data.
    Extracts FASTA descriptions and sequences and adds them to new columns in the dataframe.
    
    :param pd_file: DataFrame containing FASTA data in a specified column
    :param seq_col: Name of the column containing FASTA data
    :return: DataFrame with added columns for FASTA description and sequence
    """
    # Create a copy of the DataFrame to avoid issues with views
    pd_file = pd_file.copy()

    # Initialize new columns
    fasta_cols = ['id', 'name', 'description', 'dbxrefs', 'seq']
    pd_file[fasta_cols] = pd.NA

    # Process each row using apply
    def apply_process_row(row):
        return process_row(row[seq_col])

    results = pd_file.apply(apply_process_row, axis=1, result_type='expand')
    
    # Assign results back to the DataFrame
    pd_file[fasta_cols] = results

    return pd_file

"""
Code to download and setup the database as described in pLannotate
"""
def download_databases():

    # dynamic version number for the databases
    # this is locked at minor version bumps
    # need to upload a new database into github every minor update
    # patch number bumps just refer to the X.X.0 version
    #db_loc = f"https://github.com/mmcguffi/pLannotate/releases/download/v{plannotate_version.rsplit('.',1)[0]}.0/BLAST_dbs.tar.gz"
    db_loc = "https://github.com/barricklab/pLannotate/releases/download/v1.1.0/BLAST_dbs.tar.gz"

    # subprocess.call(["wget", "-P", f"{ROOT_DIR}/data/", db_loc])
    subprocess.call(["curl", "-L", "-o", f"{ROOT_DIR}/data/BLAST_dbs.tar.gz", db_loc])

    # check if download was successful
    if not os.path.exists(f"{ROOT_DIR}/data/BLAST_dbs.tar.gz"):
        print("Error downloading databases. Please try again or contact the developer.")
        sys.exit()

    print("Download complete.")
    print()

    print("Extracting...")
    subprocess.call(
        ["tar", "-xzf", f"{ROOT_DIR}/data/BLAST_dbs.tar.gz", "-C", f"{ROOT_DIR}/data/"] #If you already have the downloaded .tar file, change the sequence here or unzip it directly from bash.
    )
    print("Extraction complete.")
    print()

    print("Removing archive...")
    subprocess.call(["rm", f"{ROOT_DIR}/data/BLAST_dbs.tar.gz"])
    print("Removal complete.")
    print()

    print("Done.")
    print()