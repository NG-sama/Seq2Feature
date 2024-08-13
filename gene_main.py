from seqannotate import main, resources as rsc
import pandas as pd
from tqdm import tqdm

"""
List of global variables defined below
"""
save_loc = f"/workspaces/universal/Seq2Feature/seqannotate/data/gene_data/annotated files" #Location of where your annotated files are saved
read_loc = f"/workspaces/universal/Seq2Feature/seqannotate/data/gene_data/gene_data.csv" #Location of the file that is to be read and annotated.
def process_rows(row):
    """
    Runs row-wise operations to annotate each row
    """
    annotate_main = main.annotate(row['seq'])
    file_name = row['Catalog Number']
    annotate_main.to_csv(f"{save_loc}"+f"/{file_name}.csv", index = False)


def loop_primary():
    gene_data = pd.read_csv(read_loc)
    gene_list = rsc.get_csv_to_gene(gene_data[:3]) #Change the code here depending on the CSV file
    # Wrap df.apply with tqdm to track progress
    tqdm.pandas()  # Register pandas with tqdm
    gene_list.progress_apply(process_rows, axis=1)
    

if __name__ == "__main__":
    loop_primary()