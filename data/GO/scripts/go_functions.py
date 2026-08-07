import requests
import gzip
import shutil
import os
import pandas as pd

def download_files(file_name, output_dir, overwrite=False, download_obo=False):
    url = f"https://current.geneontology.org/annotations/{file_name}.gaf.gz"
    gz_file = os.path.join(output_dir, f"{file_name}.gaf.gz")
    output_file = os.path.join(output_dir, f"{file_name}.gaf")

    if os.path.exists(gz_file) and not overwrite:
        print(f"{os.path.normpath(gz_file)} already exists. Skipping download.")
    else:
        response = requests.get(url, stream=True)
        if response.status_code == 200:
            with open(gz_file, 'wb') as f:
                for chunk in response.iter_content(chunk_size=1024):
                    if chunk:
                        f.write(chunk)
            print(f"GO annotations downloaded successfully and saved to {os.path.normpath(gz_file)}")
        else:
            print(f"Failed to download GO annotations. Status code: {response.status_code}")
            return

    if os.path.exists(output_file) and not overwrite:
        print(f"{os.path.normpath(output_file)} already exists. Skipping unzip.")
    else:
        with gzip.open(gz_file, 'rb') as f_in:
            with open(output_file, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        print(f"File {os.path.normpath(gz_file)} unzipped successfully to {os.path.normpath(output_file)}")

    if download_obo:
        # Download the go-basic.obo file
        obo_url = "https://purl.obolibrary.org/obo/go/go-basic.obo"
        obo_file = os.path.join(output_dir, "go-basic.obo")

        response = requests.get(obo_url, stream=True)
        if response.status_code == 200:
            with open(obo_file, 'wb') as f:
                for chunk in response.iter_content(chunk_size=1024):
                    if chunk:
                        f.write(chunk)
            print(f"GO basic ontology downloaded successfully and saved to {os.path.normpath(obo_file)}")
        else:
            print(f"Failed to download GO basic ontology. Status code: {response.status_code}")
            return None

    return output_file

def parse_obo_file(obo_file):
    if not os.path.exists(obo_file):
        print(f"File {obo_file} does not exist.")
        return None

    go_terms = {}
    obsolete_terms = {}
    with open(obo_file, 'r') as f:
        term = None
        is_obsolete = False
        for line in f:
            line = line.strip()
            if line == "[Term]":
                if term:
                    if is_obsolete:
                        obsolete_terms[term['id']] = term
                    else:
                        go_terms[term['id']] = term
                term = {}
                is_obsolete = False
            elif line.startswith("id: "):
                term['id'] = line.split(": ")[1]
            elif line.startswith("name: "):
                term['name'] = line.split(": ")[1]
            elif line.startswith("namespace: "):
                term['namespace'] = line.split(": ")[1]
            elif line.startswith("def: "):
                term['def'] = line.split(": ")[1]
            elif line.startswith("is_a: "):
                if 'is_a' not in term:
                    term['is_a'] = []
                term['is_a'].append(line.split(" ! ")[0].split(": ")[1])
            elif line.startswith("is_obsolete: true"):
                is_obsolete = True
            elif line.startswith("replaced_by: "):
                term['replaced_by'] = line.split(": ")[1]
            elif line.startswith("consider: "):
                if 'consider' not in term:
                    term['consider'] = []
                term['consider'].append(line.split(": ")[1])
        if term:
            if is_obsolete:
                obsolete_terms[term['id']] = term
            else:
                go_terms[term['id']] = term

    print(f"Successfully parsed {len(go_terms)} GO terms and {len(obsolete_terms)} obsolete terms from {obo_file}")
    return go_terms, obsolete_terms

def read_gaf_to_dataframe(gaf_file):
    # Read the .gaf file into a pandas DataFrame
    column_names = [
        "DB", "DB_Object_ID", "DB_Object_Symbol", "Qualifier", "GO_ID", "DB_Reference",
        "Evidence_Code", "With_From", "Aspect", "DB_Object_Name", "DB_Object_Synonym",
        "DB_Object_Type", "Taxon", "Date", "Assigned_By", "Annotation_Extension", "Gene_Product_Form_ID"
    ]
    dtype = {
        "DB": str, "DB_Object_ID": str, "DB_Object_Symbol": str, "Qualifier": str, "GO_ID": str,
        "DB_Reference": str, "Evidence_Code": str, "With_From": str, "Aspect": str, "DB_Object_Name": str,
        "DB_Object_Synonym": str, "DB_Object_Type": str, "Taxon": str, "Date": str, "Assigned_By": str,
        "Annotation_Extension": str, "Gene_Product_Form_ID": str
    }
    df = pd.read_csv(gaf_file, sep='\t', comment='!', names=column_names, header=None, dtype=dtype)
    return df

def create_annotations_csv(overwrite=False):
    # Ensure the directory exists
    output_dir = os.path.normpath(os.path.join(os.path.dirname(__file__), ".."))
    os.makedirs(output_dir, exist_ok=True)

    obo_file = os.path.normpath(os.path.join(output_dir, "go-basic.obo"))  # Specify the path to the go-basic.obo file
    gaf_file = os.path.normpath(os.path.join(output_dir, "goa_human.gaf"))  # Specify the path to the .gaf file
    csv_file = os.path.normpath(os.path.join(output_dir, 'GO_clean.csv'))  # Specify the path to the output CSV file

    # Check if required files exist
    if not os.path.exists(obo_file):
        print(f"Required file {obo_file} does not exist. Aborting.")
        return
    if not os.path.exists(gaf_file):
        print(f"Required file {gaf_file} does not exist. Aborting.")
        return

    go_terms, obsolete_terms = parse_obo_file(obo_file)

    df = read_gaf_to_dataframe(gaf_file)
    
    # Number of rows before removing duplicates
    num_rows_before = len(df)
    
    # Drop duplicates
    df = df.drop_duplicates()
    
    # Number of rows after removing duplicates
    num_rows_after = len(df)
    
    # Check for duplicates
    num_duplicates = num_rows_before - num_rows_after
    if num_duplicates > 0:
        print(f"There were {num_duplicates} duplicate rows in the {gaf_file}.")
    else:
        print("There were no duplicate rows in the DataFrame.")
    
    print(f"Number of rows before removing duplicates: {num_rows_before}")
    print(f"Number of rows after removing duplicates: {num_rows_after}")

    # Example: Map annotations to GO terms
    if go_terms:
        df['GO_Term_Name'] = df['GO_ID'].map(lambda go_id: go_terms[go_id]['name'] if go_id in go_terms else None)
        df['GO_Term_Namespace'] = df['GO_ID'].map(lambda go_id: go_terms[go_id]['namespace'] if go_id in go_terms else None)
        
        # Create the desired CSV files for each namespace
        namespaces = ['biological_process', 'molecular_function', 'cellular_component']
        for namespace in namespaces:
            namespace_df = df[df['GO_Term_Namespace'] == namespace]
            csv_file = os.path.normpath(os.path.join(output_dir, f'GO_clean_{namespace}.csv'))
            
            if os.path.exists(csv_file) and not overwrite:
                print(f"CSV file '{csv_file}' already exists. Skipping creation.")
            else:
                annotation_df = namespace_df[['GO_Term_Name', 'DB_Object_Symbol']].dropna().rename(columns={'GO_Term_Name': 'annotation', 'DB_Object_Symbol': 'SYMBOL'})
                annotation_df.to_csv(csv_file, index=False)
                print(f"CSV file '{csv_file}' created successfully.")