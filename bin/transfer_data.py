#!/usr/bin/env python

__author__ = 'Matiss Ozols'
__date__ = '2021-07-28'
__version__ = '0.0.1'

import os
from shutil import copyfile,copytree
from os import listdir
import glob
import argparse



def main_data_colection(pipeline='',name='',directory='',out_directory='',extras=None):
    name_dir=f"{out_directory}/handover/{name}"
    try:
        os.mkdir(f'{name_dir}')
    except:
        print("dir exists")


    if (pipeline=='Cellbender'):
        print('prepearing Cellbender folder')
        try:
            os.mkdir(f'{name_dir}/Cellbender')
        except:
            print('dire exists')
        Folders = listdir(f'{directory}/outputs')
        print(f'coppying Celbender to {name_dir}/Cellbender')
        for folder in Folders:
            if (folder == 'qc_cluster_input_files'):
                continue
            
            copyfile(f'{directory}/outputs/{folder}/cellbender-epochs_250__learnrt_1pt0Eneg7__zdim_100__zlayer_500__lowcount_10/plots/cellbender_results-cellbender_FPR_0pt01_filtered-ambient_signature-scatter_genenames.png', f'{name_dir}/Cellbender/{folder}_ambient_signature-scatter_genenames.png')
            copyfile(f'{directory}/outputs/{folder}/cellbender-epochs_250__learnrt_1pt0Eneg7__zdim_100__zlayer_500__lowcount_10/plots/cellbender_results-cellbender_FPR_0pt01_filtered-abs_count_difference-boxplot.png', f'{name_dir}/Cellbender/{folder}_ount_difference-boxplot.png')
            copyfile(f'{directory}/outputs/{folder}/cellbender-epochs_250__learnrt_1pt0Eneg7__zdim_100__zlayer_500__lowcount_10/plots/cellbender.pdf', f'{name_dir}/Cellbender/cellbender_{folder}.pdf')
            copyfile(f'{directory}/outputs/{folder}/compare_cellranger/fpr_0.1/boxplots_cellranger_vs_cellbender.png', f'{name_dir}/Cellbender/{folder}_boxplots_cellranger_vs_cellbender.png')
            copyfile(f'{directory}/outputs/{folder}/compare_cellranger/fpr_0.1/barcode_vs_total_counts.png', f'{name_dir}/Cellbender/{folder}_barcode_vs_total_counts.png')
            copyfile(f'{directory}/outputs/{folder}/compare_cellranger/fpr_0.1/boxplot_topgenes_cellranger_vs_cellbender.png', f'{name_dir}/Cellbender/{folder}_boxplot_topgenes_cellranger_vs_cellbender.png')

    elif (pipeline=='Fetch'):
        # 'Note that the names for the future projects may be different - have to be handled on the Nextflow modules'
        print('prepearing fetch folder')
        try:
            os.mkdir(f'{name_dir}/Fetch Pipeline')
        except:
            print('exists')
        copyfile(f'{directory}/results/Submission_Data_Pilot_UKB.file_metadata.tsv', f'{name_dir}/Fetch Pipeline/Submission_Data_Pilot_UKB.file_metadata.tsv')
        print(f'coppying Fetch to: {name_dir}/Fetch Pipeline')
        Folders = listdir(f'{directory}/results/iget_study_cellranger')
        for folder in Folders:
            Folders2 = listdir(f'{directory}/results/iget_study_cellranger/{folder}')
            for folder2 in Folders2:
                print(folder2)
                Folders3 = glob.glob(f'{directory}/results/iget_study_cellranger/{folder}/{folder2}/cellranger*')
                copyfile(f'{Folders3[0]}/web_summary.html', f'{name_dir}/Fetch Pipeline/html_{folder2}.html')

    elif (pipeline=='Deconvolution'):
        print('prepearing Deconvolution folder')
        try:
            os.mkdir(f'{name_dir}/Deconvolution')
        except:
            print('dire exists')
        Folders = listdir(f'{directory}/results/split_donor_h5ad')
        for folder in Folders:
            copyfile(f'{directory}/results/split_donor_h5ad/{folder}/Vireo_plots.pdf', f'{name_dir}/Deconvolution/Vireo_plots_{folder}.pdf')
        print(f'prepearing fetch folder {name_dir}/Deconvolution')

    elif (pipeline=='QC'):
        print('prepearing QC folder')
        try:
            os.mkdir(f'{name_dir}/QC metrics')
        except:
            print('dire exists')
        try: 
            os.mkdir(f'{name_dir}/Clustering')
        except:
            print('dire exists')
        try:
            os.mkdir(f'{name_dir}/Cell-type assignment')
        except:
            print('dire exists')

        #copy the QC mertic plots    
        copyfile(f'{directory}/plots/adata-cell_desity.png', f'{name_dir}/QC metrics/adata-cell_desity.png')
        copyfile(f'{directory}/plots/adata-cell_filtered_per_experiment-n_cells_before_after.png', f'{name_dir}/QC metrics/adata-cell_filtered_per_experiment-n_cells_before_after.png')
        copyfile(f'{directory}/plots/scatterplot-sex_sample_swap_check.png', f'{name_dir}/QC metrics/scatterplot-sex_sample_swap_check.png')
        copyfile(f'{directory}/plots/adata-outlier_cells.png', f'{name_dir}/QC metrics/adata-outlier_cells.png')
        copyfile(f'{directory}/plots/plot_ecdf-x_log10.var=total_counts.color=experiment_id-adata.png', f'{name_dir}/QC metrics/plot_ecdf-x_log10.var=total_counts.color=experiment_id-adata.png')
        
        copytree(f'{directory}/handover/minimal_dataset_summary', f'{name_dir}/Summary')
        
        # copy the clustering plots
        Folders3 = glob.glob(f'{directory}/normalize=*')
        for fold1 in Folders3:
            subfolder_1=fold1.split('/')[-1]
            Folders4 = glob.glob(f'{directory}/{subfolder_1}/reduced_dims*')
            for fold2 in Folders4:
                subfolder_2=fold2.split('/')[-1]
                try:
                    copyfile(f'{directory}/{subfolder_1}/{subfolder_2}/plots/umap-total_counts.png', f'{name_dir}/Clustering/umap-total_counts.png')
                    copyfile(f'{directory}/{subfolder_1}/plots/adata-normalized_pca-knee-variance-spline=None-knee_raw.png', f'{name_dir}/Clustering/adata-normalized_pca-knee-variance-spline=None-knee_raw.png')
                    copyfile(f'{directory}/{subfolder_1}/{subfolder_2}/plots/umap-experiment_id.png', f'{name_dir}/Clustering/umap-experiment_id.png')
                    copyfile(f'{directory}/{subfolder_1}/{subfolder_2}/plots/umap-n_cells.png', f'{name_dir}/Clustering/umap-n_cells.png')
                    break
                except:
                    print('not the correct folder')

        # copy Azimuth assignements
        Folders = listdir(f'{directory}/azimuth')
        for file in Folders:
            copyfile(f'{directory}/azimuth/{file}', f'{name_dir}/Cell-type assignment/{file}')



    else:
        print('This pipeline has not currently been integrated within Web-design. Please inform Sanger HGI team if you wish this to be added.')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="""
            Collect important graphs for web visualisations
            """
    )

    parser.add_argument(
        '-v', '--version',
        action='version',
        version='%(prog)s {version}'.format(version=__version__)
    )

    parser.add_argument(
        '-p', '--pipeline',
        action='store',
        dest='pipeline',
        required=True,
        help='Which pipeline is running this job.'
    )

    parser.add_argument(
        '-n', '--name',
        action='store',
        dest='name',
        required=True,
        help='The project name.'
    )

    parser.add_argument(
        '-d', '--directory',
        action='store',
        dest='directory',
        required=True,
        help='The nextflow running directory.'
    )

    parser.add_argument(
        '-od', '--out_directory',
        action='store',
        dest='out_directory',
        required=True,
        help='The nextflow output directory.'
    )

    parser.add_argument(
        '-e', '--extras',
        action='store',
        dest='extras',
        required=False,
        help='Sometimes, nextflow has generated a dynamic paths, hence these can be parsed from this extra argument'
    )

    options = parser.parse_args()
    pipeline = options.pipeline
    name = options.name
    directory = options.directory
    out_directory=options.out_directory
    extras = options.extras
    main_data_colection(pipeline=pipeline,name=f"{name}",directory=directory,out_directory=out_directory,extras=extras)
