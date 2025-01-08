import scanpy as sc
import os
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import argparse

import warnings
warnings.filterwarnings('ignore')

def extract_t_stat(adata, file_name, use_hvg=False, write_t_stat=False, write_hvg_t_stat_mat=False):

    if use_hvg:
        hvg_num = 6000
        sc.pp.highly_variable_genes(adata, n_top_genes=hvg_num)
        adata = adata[:, adata.var['highly_variable']]
    else:
        hvg_num = adata.shape[1]

    sc.tl.rank_genes_groups(adata, 'cell_type', method='t-test')

    prediction_folder = 'prediction'
    try:
        os.mkdir(prediction_folder)
        print('Folder %s is created' % (prediction_folder))
    except FileExistsError:
        print('Folder %s already exists' % (prediction_folder))

    hvg_t_stat_df = sc.get.rank_genes_groups_df(adata, group=None)
    
    if write_t_stat:
        hvg_t_stat_df.to_csv('prediction/t_stat.%s.csv' % (file_name.split('.')[0]), header=True, index=False)
    
    type_num = len(hvg_t_stat_df['group'].unique())
    print('Number of cell types in query\'s expression profile is %d' % (type_num))

    hvg_t_stat_rows = []
    for i in range(type_num):
    
        hvg_t_stat_df_subset = hvg_t_stat_df.loc[i*hvg_num:((i+1)*hvg_num-1), :].copy()
        type = hvg_t_stat_df_subset.loc[:, 'group'].copy().drop_duplicates().values
    
        if len(type) != 1:
            print('Error occurs in hvg t-stat file.')
            break
        else:
            type = type[0]

        hvg_t_stat_row = hvg_t_stat_df_subset[['names', 'scores']].copy().set_index('names').T
        hvg_t_stat_row.index = [type]
        
        hvg_t_stat_rows.append(hvg_t_stat_row)

    query_hvg_t_stat_mat_df = pd.concat(hvg_t_stat_rows)
    if write_hvg_t_stat_mat:
        query_hvg_t_stat_mat_df.to_csv('prediction/hvg_t_stat_mat.%s.csv' % (file_name.split('.')[0]), header=True, index=True)

    return query_hvg_t_stat_mat_df

# Layer 1: T cells or not
def corr_t_cells(hvg_t_stat_mat_df, query_exp_df, layer, alpha, beta):

    print(query_exp_df.columns)
    print(hvg_t_stat_mat_df.columns)
    intersect_genes = list(set(query_exp_df.columns) & set(hvg_t_stat_mat_df.columns))
    print(intersect_genes)
    print('The number of intersect genes between query and reference is %d' % (len(intersect_genes)))
    
    df_concat = pd.concat([hvg_t_stat_mat_df, query_exp_df], join='inner')
    
    query_clusters = query_exp_df.index
    reference_types = hvg_t_stat_mat_df.index

    corr_dfs = []
    predictions = []
    for query_cluster in query_clusters:
        print(query_cluster)
        query_exp = df_concat.loc[query_cluster, :]
        corr_list = []
        
        for reference_pmid_ct in reference_types:
            reference_exp = df_concat.loc[reference_pmid_ct, :]
            corr = query_exp.corr(reference_exp)
            corr_list.append(corr)
        
        corr_df = pd.DataFrame({layer: reference_types,
                                'correlation': corr_list})
        corr_df['query_name'] = query_cluster
        corr_df = corr_df.sort_values(by='correlation', ascending=False)
        corr_dfs.append(corr_df)

        # The below criteria can be changed (addition or deletion) by users.
        prediction = corr_df[layer].values[0] if (corr_df['correlation'].values[0] - corr_df['correlation'].values[1] >= alpha) & (corr_df['correlation'].values[1] < beta) else 'N'
        predictions.append(prediction)
        
    corr_df = pd.concat(corr_dfs)

    prediction_df = pd.DataFrame({'query_name': query_clusters, 'predicted_t_cells_or_not': predictions})

    return corr_df, prediction_df

def anno_t_cells(adata, query_exp_df, file_name, reference_path, alpha, beta, layer='t_cells_or_not', write_adata=False):

    hvg_t_stat_mat_df = pd.read_csv(reference_path, header=0, index_col=0)

    corr_df, prediction_df = corr_t_cells(query_exp_df=query_exp_df, 
                                        hvg_t_stat_mat_df=hvg_t_stat_mat_df,
                                        layer=layer,
                                        alpha=alpha,
                                        beta=beta)

    corr_df.to_csv('prediction/corr.%s.%s.csv' % (layer, file_name.split('.')[0]), header=True, index=False)
    #prediction_df.to_csv('prediction/prediction.%s.%s.csv' % (layer, file_name.split('.')[0]), header=True, index=False)

    predicted_type_column = 'predicted_' + layer
    adata.obs[predicted_type_column] = 'null'
    for i in range(prediction_df.shape[0]):
        adata.obs.loc[adata.obs['cell_type']==prediction_df.loc[i, 'query_name'], predicted_type_column] = prediction_df.loc[i, predicted_type_column]
    
    if write_adata:
        adata.write_h5ad('prediction/' + file_name.replace('h5ad', predicted_type_column + '.h5ad'))

    return(adata)

# Layer 2: pmid_ct    
def corr_pmid_ct(hvg_t_stat_mat_df, query_exp_df, layer, query_adata, info_df):

    df_concat = pd.concat([hvg_t_stat_mat_df, query_exp_df], join='inner')

    query_obs_df = query_adata.obs[['cell_type', 'predicted_t_cells_or_not']].drop_duplicates().copy()

    query_clusters = query_exp_df.index

    reference_pmid_cts_t_cells = info_df.loc[info_df['compartment']=='T cells', 'pmid_ct'].values
    reference_pmid_cts_not_t_cells = info_df.loc[info_df['compartment']!='T cells', 'pmid_ct'].values

    corr_dfs = []
    predictions = []
    predicted_t_cells_or_not_list = []
    for query_cluster in query_clusters:

        query_exp = df_concat.loc[query_cluster, :]

        query_t_cells = query_obs_df.loc[query_obs_df['cell_type']==query_cluster, 'predicted_t_cells_or_not'].values[0]
        predicted_t_cells_or_not_list.append(query_t_cells)

        if query_t_cells == 'Y':
            references = reference_pmid_cts_t_cells
        elif query_t_cells == 'N':
            references = reference_pmid_cts_not_t_cells
        else:
            print('t cells error')

        corr_list = []
        for reference_pmid_ct in references:
            reference_exp = df_concat.loc[reference_pmid_ct, :]
            corr = query_exp.corr(reference_exp)
            corr_list.append(corr)
        
        corr_df = pd.DataFrame({layer: references,
                                'correlation': corr_list})
        corr_df['query_name'] = query_cluster
        corr_df = corr_df.sort_values(by='correlation', ascending=False)
        print(corr_df)
        corr_dfs.append(corr_df)

        prediction = corr_df[layer].values[0]
        predictions.append(prediction)
        
    corr_df = pd.concat(corr_dfs)

    prediction_df = pd.DataFrame({'query_name': query_clusters,
                                  'predicted_t_cells_or_not': predicted_t_cells_or_not_list,
                                  'predicted_pmid_ct': predictions})

    return corr_df, prediction_df

def anno_pmid_ct(adata, query_exp_df, file_name, reference_path, layer='pmid_ct', write_adata=False):

    info_df = pd.read_csv('references/cluster_info_summary.csv', header=0)
    info_df = info_df[['compartment', 'cell_class', 'pmid_ct', 'pd_cc_cl_tf']]
    info_df.drop_duplicates(inplace=True)
    info_df.reset_index(drop=True, inplace=True)

    hvg_t_stat_mat_df = pd.read_csv(reference_path, header=0, index_col=0)

    corr_df, prediction_df = corr_pmid_ct(query_exp_df=query_exp_df, 
                                          hvg_t_stat_mat_df=hvg_t_stat_mat_df,
                                          layer=layer,
                                          query_adata=adata,
                                          info_df=info_df)

    corr_df.to_csv('prediction/corr.%s.%s.csv' % (layer, file_name.split('.')[0]), header=True, index=False)
    #prediction_df.to_csv('prediction/prediction.%s.%s.csv' % (layer, file_name.split('.')[0]), header=True, index=False)

    predicted_type_column = 'predicted_' + layer
    adata.obs[predicted_type_column] = 'null'
    for i in range(prediction_df.shape[0]):
        adata.obs.loc[adata.obs['cell_type']==prediction_df.loc[i, 'query_name'], predicted_type_column] = prediction_df.loc[i, predicted_type_column]
    
    if write_adata:
        adata.write_h5ad('prediction/' + file_name.replace('h5ad', predicted_type_column + '.h5ad'))

    return(adata)

# Layer 3: pd_cc_cl_tf
def corr_pd_cc_cl_tf(hvg_t_stat_mat_df, query_exp_df, layer, query_adata, info_df):

    df_concat = pd.concat([hvg_t_stat_mat_df, query_exp_df], join='inner')

    query_obs_df = query_adata.obs[['cell_type', 'predicted_t_cells_or_not', 'predicted_pmid_ct']].drop_duplicates().copy()

    query_clusters = query_exp_df.index

    corr_dfs = []
    corrs_for_predictions = []
    predictions = []
    predicted_pmid_cts = []
    predicted_t_cells_or_not_list = []
    for query_cluster in query_clusters:
        print(query_cluster)

        query_exp = df_concat.loc[query_cluster, :]

        predicted_t_cells_or_not = query_obs_df.loc[query_obs_df['cell_type']==query_cluster, 'predicted_t_cells_or_not'].values[0]
        predicted_t_cells_or_not_list.append(predicted_t_cells_or_not)

        predicted_pmid_ct = query_obs_df.loc[query_obs_df['cell_type']==query_cluster, 'predicted_pmid_ct'].values[0]
        predicted_pmid_cts.append(predicted_pmid_ct)

        references = info_df.loc[info_df['pmid_ct']==predicted_pmid_ct, 'pd_cc_cl_tf'].values
        
        print('reference pd_cc_cl_tf num %d' % (len(references)))

        corr_list = []
        for reference_pmid_ct in references:
            reference_exp = df_concat.loc[reference_pmid_ct, :]
            corr = query_exp.corr(reference_exp)
            corr_list.append(corr)
        
        corr_df = pd.DataFrame({layer: references,
                                'correlation': corr_list})
        corr_df['query_name'] = query_cluster
        corr_df = corr_df.sort_values(by='correlation', ascending=False)
        corr_dfs.append(corr_df)

        prediction = corr_df[layer].values[0] if corr_df['correlation'].values[0] > 0.2 else 'others' # use 0.2 as threshold
        predictions.append(prediction)
        corr_for_prediction = corr_df['correlation'].values[0]
        corrs_for_predictions.append(corr_for_prediction)
        print('%s predicted cell type is %s with correlation of %s' % (query_cluster, prediction, str(corr_for_prediction)))
        
    corr_df = pd.concat(corr_dfs)

    prediction_df = pd.DataFrame({'query_name': query_clusters,
                                  'predicted_t_cells_or_not': predicted_t_cells_or_not_list,
                                  'predicted_pmid_ct': predicted_pmid_cts,
                                  'predicted_pd_cc_cl_tf': predictions,
                                  'correlation': corrs_for_predictions})

    return corr_df, prediction_df

def anno_pd_cc_cl_tf(adata, query_exp_df, file_name, reference_path, gamma, layer='pd_cc_cl_tf', write_adata=False):

    info_df = pd.read_csv('references/cluster_info_summary.csv', header=0)
    info_df = info_df[['compartment', 'cell_class', 'pmid_ct', 'pd_cc_cl_tf']]
    info_df.drop_duplicates(inplace=True)
    info_df.reset_index(drop=True, inplace=True)

    hvg_t_stat_mat_df = pd.read_csv(reference_path, header=0, index_col=0)

    corr_df, prediction_df = corr_pd_cc_cl_tf(query_exp_df=query_exp_df, 
                                          hvg_t_stat_mat_df=hvg_t_stat_mat_df,
                                          layer=layer,
                                          query_adata=adata,
                                          info_df=info_df)

    corr_df.to_csv('prediction/corr.%s.%s.csv' % (layer, file_name.split('.')[0]), header=True, index=False)
    prediction_df.to_csv('prediction/prediction.%s.%s.csv' % (layer, file_name.split('.')[0]), header=True, index=False)

    predicted_type_column = 'predicted_' + layer
    adata.obs[predicted_type_column] = 'null'
    for i in range(prediction_df.shape[0]):
        adata.obs.loc[adata.obs['cell_type']==prediction_df.loc[i, 'query_name'], predicted_type_column] = prediction_df.loc[i, predicted_type_column]
    
    adata = reanno_low_corr(adata, query_exp_df, file_name, gamma=gamma)

    if write_adata:
        adata.write_h5ad('prediction/' + file_name.replace('h5ad', predicted_type_column + '.h5ad'))

    return(adata)

# Reannotation
def reanno_low_corr(adata, query_exp_df, file_name, gamma):

    prediction_df = pd.read_csv('prediction/prediction.%s.%s.csv' % ('pd_cc_cl_tf', file_name.split('.')[0]), header=0)
    query_reanno_indices = (prediction_df['predicted_t_cells_or_not'] == 'Y') & (prediction_df['correlation'] < gamma)
    
    if query_reanno_indices.sum() > 0:
    
        print('############################################')
        print('Reannotation begins')
    
        query_reanno_list = prediction_df.loc[query_reanno_indices, 'query_name'].values
        prediction_df = prediction_df[~query_reanno_indices].copy()
        reanno_prediction_dfs = []
        for query_reanno in query_reanno_list:
            print('Query cluster %s needs to be reannotated' % (query_reanno))
            adata_reanno = adata[adata.obs['cell_type']==query_reanno].copy()
            adata_reanno.obs['predicted_t_cells_or_not'] = 'N'
    
            query_exp_df_reanno = query_exp_df.loc[[query_reanno]]
            reference_path = 'references/hvg_t_stat_mat.pmid_ct.csv'
            adata_reanno = anno_pmid_ct(adata_reanno, query_exp_df_reanno, 'reanno_' + file_name, reference_path)

            reference_path = 'references/hvg_t_stat_mat.pd_cc_cl_tf.csv'
            adata_reanno = anno_pd_cc_cl_tf(adata_reanno, query_exp_df_reanno, 'reanno_' + file_name, reference_path, gamma=gamma, write_adata=False)

            reanno_adata_indices = (adata.obs['cell_type'] == query_reanno)
            adata.obs.loc[reanno_adata_indices, 'predicted_t_cells_or_not'] = adata_reanno.obs['predicted_t_cells_or_not']
            adata.obs.loc[reanno_adata_indices, 'predicted_pmid_ct'] = adata_reanno.obs['predicted_pmid_ct']
            adata.obs.loc[reanno_adata_indices, 'predicted_pd_cc_cl_tf'] = adata_reanno.obs['predicted_pd_cc_cl_tf']

            reanno_prediction_df = pd.read_csv('prediction/prediction.%s.%s.csv' % ('pd_cc_cl_tf', 'reanno_' + file_name.split('.')[0]), header=0)
            reanno_prediction_dfs.append(reanno_prediction_df)
        
        reanno_prediction_dfs.append(prediction_df)    
        final_prediction_df = pd.concat(reanno_prediction_dfs)
        final_prediction_df.to_csv('prediction/prediction.%s.%s.csv' % ('pd_cc_cl_tf', file_name.split('.')[0]), header=True, index=False)

    return adata

# Perform three-layered cell-type inference
def panda(file_path, cell_type_column, alpha=0.1, beta=0, gamma=0.2, t_cell_only=False, write_adata=False, plot=True, preprocess=False):

    adata = sc.read_h5ad(file_path)
    file_name = file_path.split('/')[-1]
    adata.var_names_make_unique()

    if cell_type_column != 'cell_type':
        adata.obs['cell_type'] = adata.obs[cell_type_column]
    
    if t_cell_only:
        query_exp_df = extract_t_stat(adata, file_name, use_hvg=True)
        adata.obs['predicted_t_cells_or_not'] = 'Y'
    else:
        query_exp_df = extract_t_stat(adata, file_name)
        reference_path = 'references/hvg_t_stat_mat.t_cells_or_not.csv'
        adata = anno_t_cells(adata, query_exp_df, file_name, reference_path, alpha=alpha, beta=beta, write_adata=write_adata)

    reference_path = 'references/hvg_t_stat_mat.pmid_ct.csv'
    adata = anno_pmid_ct(adata, query_exp_df, file_name, reference_path, write_adata=write_adata)

    reference_path = 'references/hvg_t_stat_mat.pd_cc_cl_tf.csv'
    adata = anno_pd_cc_cl_tf(adata, query_exp_df, file_name, reference_path, gamma=gamma, write_adata=write_adata)

    if plot:
        plot_umap(adata, file_name, preprocess=preprocess)

    return adata

# Plot UMAP visualization of the original and inferred clusters in query data
def plot_umap(adata, file_name, predicted_type='predicted_pd_cc_cl_tf', preprocess=False):
    
    if preprocess:
        sc.tl.pca(adata)
        sc.pp.neighbors(adata)
        sc.tl.umap(adata)

        adata.write('prediction/' + file_name.replace('h5ad', 'umap.h5ad'))

    fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(12, 6))
    sc.pl.umap(adata, color=['cell_type'], legend_loc='on data',ax=axs[0], show=False) 
    axs[0].set_title('%s: original cell type' % (file_name.split('.')[0])) 

    sc.pl.umap(adata, color=[predicted_type], legend_loc='on data', ax=axs[1], show=False)
    axs[1].set_title('%s: predicted cell type' % (file_name.split('.')[0])) 

    plt.tight_layout()
    plt.show()
    fig.savefig('prediction/umap.%s.png' % (file_name.split('.')[0]))

def main():
    parser = argparse.ArgumentParser(description="Use PANDA the three-layer annotation tool to infer cell types of blood scRNA-seq datasets.")

    parser.add_argument('file_path', type=str, help="The file path of the query dataset")
    parser.add_argument('cell_type_column', type=str, help="The name of the cell type column in the query dataset")
    
    # optional arguments
    parser.add_argument('--alpha', type=float, default=0.1, help="The 1st correlation threshold for Layer 1 inference (T cells or not)")
    parser.add_argument('--beta', type=float, default=0, help="The 2nd correlation threshold for Layer 1 inference (T cells or not)")
    parser.add_argument('--gamma', type=float, default=0.2, help="The correlation threshold for reannotation")
    
    parser.add_argument('-t', '--t_cell_only', action='store_true', help="The correlation threshold for reannotation")
    
    parser.add_argument('-w', '--write_adata', action='store_true', help="The correlation threshold for reannotation")
    parser.add_argument('-p', '--plot', action='store_true', help="The correlation threshold for reannotation")
    parser.add_argument('-pp', '--preprocess', action='store_true', help="The correlation threshold for reannotation")

    args = parser.parse_args()

    panda(file_path=args.file_path, 
          cell_type_column=args.cell_type_column, 
          alpha=args.alpha, 
          beta=args.beta, 
          gamma=args.gamma, 
          t_cell_only=args.t_cell_only, 
          write_adata=args.write_adata, 
          plot=args.plot, 
          preprocess=args.preprocess)

if __name__ == "__main__":
    main()