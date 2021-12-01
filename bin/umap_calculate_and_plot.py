#!/usr/bin/env python

__date__ = '2020-03-13'
__version__ = '0.0.1'

import argparse
import os
import random
import numpy as np
import pandas as pd
import scanpy as sc
import itertools
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib import gridspec
import warnings

# Silence NumbaPerformanceWarning in umap. See below:
# https://github.com/lmcinnes/umap/issues/252
try:
    from numba.errors import NumbaPerformanceWarning
    warnings.filterwarnings('ignore', category=NumbaPerformanceWarning)
except ImportError:
    pass

# To error out when 'FloatingPointError: divide by zero encountered in power'
# from UMAP. This is usually caused by a poor choice of min_dist and spread
# parameters.
np.seterr(all='raise')

# Set seed for reproducibility
seed_value = 0
# 0. Set `PYTHONHASHSEED` environment variable at a fixed value
# os.environ['PYTHONHASHSEED']=str(seed_value)
# 1. Set `python` built-in pseudo-random generator at a fixed value
random.seed(seed_value)
# 2. Set `numpy` pseudo-random generator at a fixed value
np.random.seed(seed_value)

# Nice large palette.
COLORS_LARGE_PALLETE = [
    '#0F4A9C', '#3F84AA', '#C9EBFB', '#8DB5CE', '#C594BF', '#DFCDE4',
    '#B51D8D', '#6f347a', '#683612', '#B3793B', '#357A6F', '#989898',
    '#CE778D', '#7F6874', '#E09D37', '#FACB12', '#2B6823', '#A0CC47',
    '#77783C', '#EF4E22', '#AF1F26'
]


# This function is based off of scanpy:
# https://github.com/theislab/scanpy/blob/master/scanpy/plotting/_tools/scatterplots.py
def panel_grid(hspace, wspace, ncols, num_panels):
    """Init plot."""
    n_panels_x = min(ncols, num_panels)
    n_panels_y = np.ceil(num_panels / n_panels_x).astype(int)
    if wspace is None:
        #  try to set a wspace that is not too large or too small given the
        #  current figure size
        wspace = 0.75 / rcParams['figure.figsize'][0] + 0.02
    # each panel will have the size of rcParams['figure.figsize']
    fig = plt.figure(
        figsize=(
            n_panels_x * rcParams['figure.figsize'][0] * (1 + wspace),
            n_panels_y * rcParams['figure.figsize'][1],
        )
    )
    left = 0.2 / n_panels_x
    bottom = 0.13 / n_panels_y
    gs = gridspec.GridSpec(
        nrows=n_panels_y,
        ncols=n_panels_x,
        left=left,
        right=1 - (n_panels_x - 1) * left - 0.01 / n_panels_x,
        bottom=bottom,
        top=1 - (n_panels_y - 1) * bottom - 0.1 / n_panels_y,
        hspace=hspace,
        wspace=wspace
    )
    return fig, gs


def save_plot(
    adata,
    list__umap_keys,
    out_file_base,
    color_var,
    colors_quantitative=True,
    colors_large_palette=COLORS_LARGE_PALLETE,
    drop_legend=-1
):
    """Save a plot."""
    fig, grid = panel_grid(
        hspace=0.40,  # NOTE: if increase paramters in title, change this
        wspace=None,
        ncols=4,
        num_panels=len(list__umap_keys)
    )
    i__ax = 0
    for i__umap in list__umap_keys:

        # Set the facet title.
        plt_title = ''
        for parameter, value in list__umap_keys[i__umap].items():
            plt_title = '{}{}={}\n'.format(
                plt_title,
                parameter,
                value
            )
        plt_title = plt_title.rstrip()

        # Get the proper axis for this plot.
        ax = plt.subplot(grid[i__ax])

        # We could avoid this line by setting basis=i__umap, but then not
        # consistent axis labels (e.g., X_umap__n_neighbors_151).
        adata.obsm['X_umap'] = adata.obsm[i__umap]

        # If not colors_quantitative, then boot up categorical plots.
        legend_loc = 'right margin'
        color_palette = 'viridis'
        if colors_quantitative is False:
            # Cast to category - required for booleans.
            adata.obs[color_var] = adata.obs[color_var].astype('category')
            n_categories = len(adata.obs[color_var].cat.categories)
            color_palette = None
            if n_categories <= len(plt.get_cmap('Dark2').colors):
                color_palette = 'Dark2'
            elif n_categories <= len(colors_large_palette):
                color_palette = colors_large_palette
            if drop_legend >= 0 and n_categories >= drop_legend:
                legend_loc = None

        # For some reason, there is an error if return_fig = False, but
        # not when show = False.
        # sc.pl.embedding(
        #     basis='X_umap',
        if color_var != 'embedding_density':
            sc.pl.umap(
                adata=adata,
                color=color_var,
                palette=color_palette,
                alpha=0.4,
                title=plt_title,
                legend_loc=legend_loc,
                ax=ax,
                show=False
            )
        else:
            # NOTE: i__umap looks something like X_umap__n_neighbors=15...
            adata.obs['umap_density'] = adata.obs[
                '{}__density'.format(i__umap).replace('X_', '')
            ]
            adata.uns['umap_density_params'] = adata.uns[
                '{}__density_params'.format(i__umap).replace('X_', '')
            ]
            sc.pl.embedding_density(
                adata=adata,
                basis='umap',
                alpha=0.4,
                title=plt_title,
                ax=ax,
                show=False
            )
            del adata.obs['umap_density']
            del adata.uns['umap_density_params']

        del adata.obsm['X_umap']
        i__ax += 1

    fig.savefig(
        '{}-{}.png'.format(out_file_base, color_var),
        #dpi=300,
        bbox_inches='tight'
    )


def main():
    """Run CLI."""
    parser = argparse.ArgumentParser(
        description="""
            Read AnnData object and PCs file. Generates UMAP.
            """
    )

    parser.add_argument(
        '-v', '--version',
        action='version',
        version='%(prog)s {version}'.format(version=__version__)
    )

    parser.add_argument(
        '-h5', '--h5_anndata',
        action='store',
        dest='h5',
        required=True,
        help='H5 AnnData file.'
    )

    parser.add_argument(
        '-pc', '--tsv_pcs',
        action='store',
        dest='pc',
        default='',
        help='Tab-delimited file of PCs for each cell. First column is\
            cell_barcode. Subsequent columns are PCs. If "", uses pca\
            slot in AnnData.\
            (default: "")'
    )

    parser.add_argument(
        '-cq', '--colors_quantitative',
        action='store',
        dest='cq',
        default='',
        help='Comma seperated list of quantitative variable names for colors.\
            (default: "")'
    )

    parser.add_argument(
        '-cc', '--colors_categorical',
        action='store',
        dest='cc',
        default='',
        help='Comma seperated list of categorical variable names for colors.\
            (default: "")'
    )

    parser.add_argument(
        '-npc', '--number_pcs',
        action='store',
        dest='npc',
        default=0,
        type=int,
        help='Number of PCs to use.\
            (default: maximum number in tsv_pcs file)'
    )

    parser.add_argument(
        '-nn', '--n_neighbors',
        action='store',
        dest='n_neighbors',
        default='15',
        type=str,
        help='Number of neighbors for sc.pp.neighbors call\
            (default: %(default)s)'
    )

    parser.add_argument(
        '-uinit', '--umap_init',
        action='store',
        dest='umap_init',
        default='X_pca',
        help='How to initialize the low dimensional embedding.\
            Valid options: any key for adata.obsm,\
            ’paga’: positions from paga(),\
            ’spectral’: use a spectral embedding of the graph,\
            ’random’: assign initial embedding positions at random.\
            (default: X_pca, the slot where tsv_pcs is stored if provided)'
    )

    parser.add_argument(
        '-umd', '--umap_min_dist',
        action='store',
        dest='umap_min_dist',
        default='0.5',
        type=str,
        help='The effective minimum distance between embedded points. Smaller\
            values will result in a more clustered/clumped embedding where\
            nearby points on the manifold are drawn closer together, while\
            larger values will result on a more even dispersal of points.\
            The value should be set relative to the spread value, which\
            determines the scale at which embedded points will be spread out.\
            (default: %(default)s)'
    )

    parser.add_argument(
        '-us', '--umap_spread',
        action='store',
        dest='umap_spread',
        default='1.0',
        type=str,
        help='The minimum distance apart that points are allowed to be in the\
            low dimensional representation (effective scale of embedded points\
            ). In combination with min_dist this determines how\
            clustered/clumped the embedded points are.\
            (default: %(default)s)'
    )

    parser.add_argument(
        '-dln', '--drop_legend_n',
        action='store',
        dest='drop_legend',
        default=-1,
        type=int,
        help='Drop the legend for categorical colors with >= drop_legend_n\
            categories. If drop_legend_n < 0, then no legend drops are\
            performed.\
            (default: %(default)s)'
    )

    parser.add_argument(
        '--force_recalculate_neighbors',
        action='store_true',
        dest='calculate_neighbors',
        default=False,
        help='Calculate neighbor graph even if it already exists in the\
            AnnData (which it my do so if you already ran BBKNN).\
            (default: %(default)s)'
    )

    parser.add_argument(
        '-ncpu', '--number_cpu',
        action='store',
        dest='ncpu',
        default=4,
        type=int,
        help='Number of CPUs to use.\
            (default: %(default)s)'
    )

    parser.add_argument(
        '-of', '--output_file',
        action='store',
        dest='of',
        default='',
        help='Basename of output files, assuming output in current working \
            directory.\
            (default: <h5_anndata>-<tsv_pcs>-umap)'
    )

    options = parser.parse_args()

    # Fixed settings.
    verbose = True

    # Scanpy settings
    sc.settings.figdir = os.getcwd()  # figure output directory to match base.
    sc.settings.n_jobs = options.ncpu  # number CPUs
    # sc.settings.max_memory = 500  # in Gb
    #sc.set_figure_params(dpi_save=300)

    # Load the AnnData file.
    adata = sc.read_h5ad(filename=options.h5)

    # Load the PCs.
    if options.pc == '':
        df_pca = pd.DataFrame(
            data=adata.obsm['X_pca'],
            index=adata.obs.index,
            columns=[
                'PC{}'.format(x) for x in
                range(1, adata.obsm['X_pca'].shape[1]+1)
            ]
        )
    else:
        df_pca = pd.read_csv(options.pc, sep='\t', index_col='cell_barcode')
    # df_pca = pd.read_csv(
    #     'adata-pcs-harmony.tsv.gz',
    #     sep='\t',
    #     index_col='cell_barcode'
    # )

    # Check that nPCs is valid.
    n_pcs = options.npc
    if 'neighbors' in adata.uns and not options.calculate_neighbors:
        # If we are using the pre-calculated neighbors use the PCs from that
        n_pcs = adata.uns['neighbors']['params']['n_pcs']
    if n_pcs == 0:
        n_pcs = len(df_pca.columns)
    elif n_pcs > len(df_pca.columns):
        raise Exception(
            '--number_pcs ({}) is > than n_pcs in --tsv_pcs ({}).'.format(
                n_pcs,
                len(df_pca.columns)
            )
        )
    if verbose:
        print('Using {} PCs.'.format(n_pcs))
    # Subset number of PCs to be exactly nPCs - here we assume PCs are ordered.
    print('Subetting PCs - we assume they are ordered by column index.')
    df_pca = df_pca.iloc[:, range(0, n_pcs)]
    print('PC columns:\t{}'.format(np.array_str(df_pca.columns)))

    # Add the reduced dimensions to the AnnData object.
    # NOTE: We need to do this for BBKNN in the case were we init with X_pca
    adata.obsm['X_pca__umap'] = df_pca.loc[adata.obs.index, :].values.copy()

    # Get the init position for UMAP
    umap_init = options.umap_init
    if umap_init == 'X_pca':
        umap_init = 'X_pca__umap'

    # Get the out file base.
    out_file_base = options.of
    if out_file_base == '':
        out_file_base = '{}-{}-umap'.format(
            os.path.basename(options.h5.rstrip('h5ad').rstrip('.')),
            os.path.basename(options.pc.rstrip('tsv.gz').rstrip('.'))
        )
    # Append the parameters to the output file.
    out_file_base = '{},number_pcs={}'.format(
        out_file_base,
        n_pcs
    )
    out_file_base = '{},umap_init={}'.format(
        out_file_base,
        options.umap_init
    )

    # Parse the color variables.
    colors_quantitative = []
    if options.cq != '':
        colors_quantitative = options.cq.split(',')

    colors_categorical = []
    if options.cc != '':
        colors_categorical = options.cc.split(',')

    if len(colors_quantitative) == 0 and len(colors_categorical) == 0:
        raise Exception('Specify a color value.')

    # Add colors_large_palette to adata.uns.
    # adata.uns["annotation_colors"] = COLORS_LARGE_PALLETE

    # Parse the neighbors iterations.
    list__n_neighbors = []
    if options.n_neighbors != '':
        list__n_neighbors = list(map(int, options.n_neighbors.split(',')))

    # Parse the min_dist iterations.
    list__min_dist = []
    if options.umap_min_dist != '':
        list__min_dist = list(map(float, options.umap_min_dist.split(',')))

    # Parse the neighbors iterations.
    list__spread = []
    if options.umap_spread != '':
        list__spread = list(map(float, options.umap_spread.split(',')))

    # Update the out base if only one of any iteration.
    if len(list__n_neighbors) == 1:
        out_file_base = '{},n_neighbors={}'.format(
            out_file_base,
            list__n_neighbors[0]
        )
    if len(list__min_dist) == 1:
        out_file_base = '{},umap_min_dist={}'.format(
            out_file_base,
            list__min_dist[0]
        )
    if len(list__spread) == 1:
        out_file_base = '{},umap_spread={}'.format(
            out_file_base,
            list__spread[0]
        )

    # Loop over all combinations of the different paramters we want to analyse.
    list__umap_keys = {}
    for i__n_neighbors, i__min_dist, i__spread in itertools.product(
        list__n_neighbors, list__min_dist, list__spread
    ):
        # Check input parameters
        if not (2 <= i__n_neighbors <= 100):
            # Recommended in parameter documentation:
            # https://umap-learn.readthedocs.io/en/latest/api.html
            warnings.warn(
                'WARNING: it is suggested to set n_neighbors to a {}'.format(
                    'value between 2-100.'
                )
            )
        if not (0.0 <= i__min_dist <= 1.0):
            # Recommended here: https://github.com/lmcinnes/umap/issues/249
            warnings.warn(
                'WARNING: it is suggested to set umap_min_dist to a {}'.format(
                    'value between 0-1.'
                )
            )
        if not (0.0 <= i__spread <= 3.0):
            # Recommendation based on single cell experience.
            warnings.warn(
                'WARNING: it is suggested to set umap_spread to a {}'.format(
                    'value between 0-3.'
                )
            )

        # Set the plot label.
        plt__label = 'n_neighbors={}'.format(i__n_neighbors)
        plt__label = '{},umap_min_dist={}'.format(
            plt__label,
            str(i__min_dist).replace('.', 'pt')
        )
        plt__label = '{},umap_spread={}'.format(
            plt__label,
            str(i__spread).replace('.', 'pt')
        )

        # Calculate neighbors for on the specified PCs.
        # By default saved to adata.uns['neighbors']
        #
        # First, however, check to see if adata.uns['neighbors'] already exists
        # ...and unless the user tells us not to, use that slot, not calculate
        # neighbors. This default behaviour is to accommodate the instance when
        # bbknn has been run on the data.
        if 'neighbors' not in adata.uns or options.calculate_neighbors:
            sc.pp.neighbors(
                adata,
                use_rep='X_pca__umap',
                n_pcs=n_pcs,
                n_neighbors=i__n_neighbors,  # Scanpy default = 15
                copy=False,
                random_state=0
            )
        else:
            warnings.warn(
                'WARNING: found neighbors slot in adata.uns. {}'.format(
                    'Not calculating neighbors (ignoring n_neighbors).'
                )
            )
            # If we are using the pre-calculated neighbors drop npcs note.
            # if 'n_pcs' in adata.uns['neighbors']['params']:
            # n_pcs = adata.uns['neighbors']['params']['n_pcs']
            i__n_neighbors = adata.uns['neighbors']['params']['n_neighbors']

        # Save the parameters to a dict
        list__umap_keys['X_umap__{}'.format(plt__label)] = {
            'n_neighbors': i__n_neighbors,
            'umap_min_dist': i__min_dist,
            'umap_spread': i__spread
        }

        adata.uns['neighbors__{}'.format(plt__label)] = adata.uns['neighbors']

        # TODO: add paga
        # # If init with paga, plot paga first - NOTE we can only do this if
        # if options.umap_init == 'paga' and 'paga' not in adata.uns:
        #     print(
        #         'Trying to call sc.tl.paga.',
        #         'NOTE: requires one to have clustered the data.'
        #     )
        #     sc.tl.paga(
        #         adata,
        #         use_rna_velocity=False,
        #         copy=False
        #     )

        # UMAP
        # Saved to adata.uns['umap'] and adata.obsm['X_umap']
        # NOTE: If umap_init == X_pca, then X_umap will have an equal number
        #       of n_components to X_pca (n_components is overridden).
        sc.tl.umap(
            adata,
            n_components=2,
            min_dist=i__min_dist,  # Scanpy default = 0.05
            spread=i__spread,  # Scanpy default = 1.0
            init_pos=umap_init,  # Scanpy default = spectral
            # For some reason cannot access neighbors key slot, thus we
            # must keep uns['neighbors'] until we have run this.
            # neighbors_key='neighbors__{}'.format(plt__label),
            copy=False,
            random_state=0
        )

        if 'embedding_density' in colors_quantitative:
            sc.tl.embedding_density(
                adata,
                basis='umap'
            )
            # Rename density estimates
            adata.obs[
                'umap__{}__density'.format(plt__label)
            ] = adata.obs.pop('umap_density')
            adata.uns[
                'umap__{}__density_params'.format(plt__label)
            ] = adata.uns.pop('umap_density_params')

        # Rename UMAP
        adata.uns[
            'umap__{}__params'.format(plt__label)
        ] = adata.uns.pop('umap')
        adata.obsm[
            'X_umap__{}'.format(plt__label)
        ] = adata.obsm.pop('X_umap')

        # Delete key that we no longer need since already copied and we have
        # run umap.
        del adata.uns['neighbors']

    # NOTE: If the color var is a gene, you should color by ln(CPM+1).
    #       By default these sc.pl.umap uses the .raw attribute of AnnData
    #       if present which is assumed to be ln(CPM+1).

    # For each color to plot, loop over the different iterations.
    for color_var in colors_quantitative:
        save_plot(
            adata=adata,
            list__umap_keys=list__umap_keys,
            out_file_base=out_file_base,
            color_var=color_var,
            colors_quantitative=True,
            drop_legend=options.drop_legend
        )
    for color_var in colors_categorical:
        save_plot(
            adata=adata,
            list__umap_keys=list__umap_keys,
            out_file_base=out_file_base,
            color_var=color_var,
            colors_quantitative=False,
            drop_legend=options.drop_legend
        )

    adata.write(
        '{}.h5ad'.format('test'),
        compression='gzip'
    )
    # In some ocassions, you might still observe disconnected clusters and
    # similar connectivity violations. They can usually be remedied by running:
    # sc.tl.paga(adata)
    # From below, remove `plot=False` if you want to see the coarse-grained
    # graph
    # sc.pl.paga(adata, plot=False)
    # sc.tl.umap(adata, init_pos='paga')


if __name__ == '__main__':
    main()
