from go_functions import download_files
import os

if __name__ == "__main__":
    file_name = "goa_human"  # Change this to the desired file name to download
    overwrite = False  # Set to True if you want to overwrite existing files
    download_obo = False  # Set to True if you want to download the go-basic.obo file

    # Ensure the directory exists
    output_dir = os.path.normpath(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
    os.makedirs(output_dir, exist_ok=True)

    gaf_file = download_files(file_name, output_dir, overwrite=overwrite, download_obo=download_obo)