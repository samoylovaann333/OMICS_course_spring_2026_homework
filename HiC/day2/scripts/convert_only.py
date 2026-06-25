import os
import subprocess
import sys
import cooler
from hic2cool.hic2cool_utils import hic2cool_convert

DATA_DIR = "data"
NPROC = 1
RESOLUTIONS_TO_BALANCE = [10_000, 100_000, 1_000_000]

SAMPLES = ["MoPh7_enr_v2", "MoPh11_enr_v2", "MoPh14_enr_v2", "MoPh15_enr_v2"]

def balance_mcool(mcool_path, nproc=NPROC):
    print(f"Balancing {mcool_path}...")
    resolutions = cooler.fileops.list_coolers(mcool_path)
    selected_resolutions = {f"/resolutions/{r}" for r in RESOLUTIONS_TO_BALANCE}
    for resolution_group in resolutions:
        if resolution_group not in selected_resolutions:
            continue
        uri = f"{mcool_path}::{resolution_group}"
        if "weight" in cooler.Cooler(uri).bins().columns:
            print(f"  {uri} already balanced. Skipping.")
            continue
        print(f"  Balancing {uri}")
        cmd = [sys.executable, "-m", "cooler", "balance",
               "--ignore-diags", "2", "--mad-max", "5", "--min-nnz", "10",
               "--nproc", str(nproc), uri]
        subprocess.run(cmd, check=True)

for sample in SAMPLES:
    hic_path = os.path.join(DATA_DIR, f"{sample}.hic")
    cool_path = os.path.join(DATA_DIR, f"{sample}.mcool")

    if not os.path.exists(hic_path):
        print(f"MISSING: {hic_path}")
        continue

    if not os.path.exists(cool_path):
        print(f"Converting {hic_path} to {cool_path}...")
        hic2cool_convert(hic_path, cool_path, resolution=0, nproc=NPROC)

    balance_mcool(cool_path, nproc=NPROC)

print("All conversions and balancing done!")
