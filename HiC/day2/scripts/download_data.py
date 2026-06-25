import os
import requests
import subprocess
import sys
from tqdm import tqdm
import cooler
import urllib3
from hic2cool.hic2cool_utils import hic2cool_convert

# URLs for Hi-C maps
URLS = [
    "https://genedev.bionet.nsc.ru/ftp/by_User/DashaPanchenko/OMICS_course_spring_2026/day1/data/MoPh7_enr_v2.hic",
    "https://genedev.bionet.nsc.ru/ftp/by_User/DashaPanchenko/OMICS_course_spring_2026/day1/data/MoPh14_enr_v2.hic",
    "https://genedev.bionet.nsc.ru/ftp/by_User/DashaPanchenko/OMICS_course_spring_2026/day1/data/MoPh15_enr_v2.hic",
    "https://genedev.bionet.nsc.ru/ftp/by_User/DashaPanchenko/OMICS_course_spring_2026/day1/data/MoPh11_enr_v2.hic"
]

DATA_DIR = "data"
NPROC = 1
VERIFY_SSL = False
RESOLUTIONS_TO_BALANCE = [10_000, 100_000, 1_000_000]

if not VERIFY_SSL:
    urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

def download_file(url, dest_folder):
    """Downloads a file from a URL to a output folder"""
    os.makedirs(dest_folder, exist_ok=True)
    filename = url.split("/")[-1]
    filepath = os.path.join(dest_folder, filename)

    if os.path.exists(filepath):
        print(f"{filename} already exists. Skipping download.")
        return filepath

    print(f"Downloading {filename}...")
    response = requests.get(url, stream=True, verify=VERIFY_SSL)
    response.raise_for_status()

    total_size = int(response.headers.get('content-length', 0))
    block_size = 1024 * 1024

    with open(filepath, 'wb') as file, tqdm(
        desc=filename,
        total=total_size,
        unit='iB',
        unit_scale=True,
        unit_divisor=1024,
    ) as bar:
        for data in response.iter_content(block_size):
            bar.update(len(data))
            file.write(data)
            
    return filepath

def balance_mcool(mcool_path, nproc=NPROC):
    """
    Balance all resolutions inside an .mcool file.
    Balancing writes weights into bins/weight for each resolution.
    """

    print(f"Balancing {mcool_path}...")

    resolutions = cooler.fileops.list_coolers(mcool_path)
    selected_resolutions = {
        f"/resolutions/{resolution}"
        for resolution in RESOLUTIONS_TO_BALANCE
    }

    for resolution_group in resolutions:
        if resolution_group not in selected_resolutions:
            continue

        uri = f"{mcool_path}::{resolution_group}"

        if "weight" in cooler.Cooler(uri).bins().columns:
            print(f"  {uri} already balanced. Skipping.")
            continue

        print(f"  Balancing {uri}")

        cmd = [
            sys.executable,
            "-m", "cooler",
            "balance",
            "--ignore-diags", "2",
            "--mad-max", "5",
            "--min-nnz", "10",
            "--nproc", str(nproc),
            uri
        ]

        subprocess.run(cmd, check=True)

if __name__ == "__main__":
    print("Starting data download...")
    for url in URLS:
        hic_path = download_file(url, DATA_DIR)

        cool_path = hic_path.replace(".hic", ".mcool")

        if not os.path.exists(cool_path):
            print(f"Converting {hic_path} to {cool_path}...")

            # resolution=0 creates a multi-resolution .mcool
            hic2cool_convert(
                hic_path,
                cool_path,
                resolution=0,
                nproc=NPROC
            )

            balance_mcool(cool_path, nproc=NPROC)

        else:
            print(f"{cool_path} already exists. Skipping conversion.")
            balance_mcool(cool_path, nproc=NPROC)
            
            
    print("All files downloaded and converted successfully!")
