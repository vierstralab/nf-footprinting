from genome_tools import df_to_genomic_intervals

from genome_tools.plotting.modular_plot.api import DataBundle
from genome_tools.plotting.modular_plot.loaders.footprint import FootprintsDataLoader

from differential_test_api import DifferentialFitLoader, WindowedScoreLoader, DifferentialModelConfig

import sys
import pandas as pd



def extract_data_for_dhs_interval(interval, sample_data):
    data = DataBundle(interval=interval)

    data.groups_data = sample_data.loc[:, 'extended_annotation'].copy()
    data.grouping_column = 'extended_annotation'

    data = FootprintsDataLoader()._load(
        data,
        footprints_metadata=sample_data,
        calc_posteriors=False
    )
    print(data.obs.shape)
    data = DifferentialFitLoader()._load(
        data,
        keep_sig2_loglik=True,
        config=DifferentialModelConfig(
            mu_grid_params=(-5.0, 5.0, 201),
            sig2_grid_params=(1e-3, 10.0, 81)
        )
    )
    print('diff loaded')
    data = WindowedScoreLoader()._load(data)

    return data


if __name__ == "__main__":
    dhs_id = sys.argv[1]  # 'chunk2499_1915_8'
    dhs_index = pd.read_table(sys.argv[2]).set_index('dhs_id')

    sample_data = pd.read_table(sys.argv[3]).set_index('sample_id')

    # for dhs_id in human_anndata.var.query('num_samples >= 4000').index[4:5]:
    dhs_region = df_to_genomic_intervals(
        dhs_index.loc[[dhs_id]].reset_index(),
        extra_columns=['dhs_id']
    )[0]

    data = extract_data_for_dhs_interval(
        dhs_region,
        sample_data
    )

    data.differential.to_npz(path=sys.argv[4])
    # dhs_region = GenomicInterval("chr15", 74995300, 74995600)



