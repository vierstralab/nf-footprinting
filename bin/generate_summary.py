import sys

from footprint_tools.stats.differential_bayesian.api import load_data_results

if __name__ == "__main__":
    dhs_id, npz_path, output_path = sys.argv[1:4]
    table = load_data_results(npz_path).footprints.region_table()
    table.insert(0, "dhs_id", dhs_id)
    table.to_csv(output_path, sep="\t", index=False)
