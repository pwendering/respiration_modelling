#!/usr/bin/python3
from matplotlib import pyplot as plt
import seaborn as sb
import pandas as pd
import numpy as np
from pathlib import Path


class TissueGeneExpression:

    def __init__(self, expr_file, tissue_dict_file):
        self.expr_file = expr_file
        self.tissue_expr_file = 'gene_expr_per_tissue.csv'
        self.geneIDs = [
            ["At1g59900", "At5g50850", "At1g24180"],  # EC 1.2.4.1
            ["At1g30120", "At1g01090", "At2g34590"],  # EC 1.2.4.1
            ["At4g35260", "At2g17130", "At4g35650", "At5g03290", "At3g09810"],  # EC 1.1.1.41, 1.1.1.286
            ["At1g65930"],  # EC 1.1.1.42
            ["At5g14590"],  # EC 1.1.1.42
            ["At3g55410", "At5g65750", "At4g26910", "At5g55070"],  # E1 and E2 subunit
            ["At1g64190"],  # EC 1.1.1.44 (cyt/chl)
            ["At3g02360"],  # EC 1.1.1.44 (cyt/perox) https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4854672/#s9title
            ["At5g41670"],  # EC 1.1.1.44 (cyt/chl)
            ["At2g19900"],  # EC 1.1.1.40 (cytosol) NADP-ME1
            ["At5g11670"],  # EC 1.1.1.40 (cytosol) NADP-ME2
            ["At5g25880"],  # EC 1.1.1.40 (cytosol) NADP-ME3
            ["At1g79750"],  # EC 1.1.1.40 (chloroplast) NADP-ME4
            ["At2g13560"],  # EC 1.1.1.39 (mitochondrion) NAD-ME1
            ["At4g00570"],  # EC 1.1.1.39 (mitochondrion) NAD-ME2
        ]
        self.enzymes = [
            "PDH (mit)",
            "PDH (chl)",
            "IDH (mit)",
            "IDH (cyt)",
            "IDH (mit/chl)",
            "2OGDH (mit)",
            "PGD1 (cyt/chl)",
            "PGD2 (cyt/per)",
            "PGD3 (cyt/chl)",
            "NADP-ME1 (cyt)",
            "NADP-ME2 (cyt)",
            "NADP-ME3 (cyt)",
            "NADP-ME4 (chl)",
            "NAD-ME1 (mit)",
            "NAD-ME2 (mit)",
        ]

        self.tissue_dict = self.get_tissue_dict(tissue_dict_file)

    def get_tissue_dict(self, tissue_dict_file):
        tissue_dict = dict()
        with open(tissue_dict_file, 'r') as tdfile:
            tdfile.readline()
            for line in tdfile:
                split = line.split("\t")
                tissue_dict[split[1].strip()] = split[0]
        return tissue_dict

    def create_expr_df(self):
        df_orig = pd.read_csv(self.expr_file, sep=";")

        print("Creating file for tissue-specific expression")

        tissues_uniq = list(set(self.tissue_dict.keys()))

        # create dataframe with columns: gene_id, tpm, tissue, pathway
        gene_ids = []
        tpm = []
        tissues = []
        enzymes = []

        # loop over gene IDs
        for i in range(0, df_orig.shape[0]):

            g_id = df_orig.iloc[i]['AGI_code'].lower()

            # loop over enzymes
            for j in range(0, len(self.enzymes)):
                if g_id in [x.lower() for x in self.geneIDs[j]]:
                    # loop over tissues
                    for k in range(0, len(tissues_uniq)):
                        t_idx = df_orig.columns.str.fullmatch("TPM_" + tissues_uniq[k])
                        ibaq_val = df_orig.iloc[i, t_idx][0]

                        if not np.isnan(ibaq_val):
                            gene_ids.append(g_id)
                            tpm.append(ibaq_val)
                            enzymes.append(self.enzymes[j])
                            tissues.append(self.tissue_dict[tissues_uniq[k]])

        df_new = pd.DataFrame.from_dict({'gene_id': gene_ids,
                                         'tpm': tpm,
                                         'tissue': tissues,
                                         'enzyme': enzymes})

        df_new.to_csv(self.tissue_expr_file, index=False, quoting=False)

        return df_new

    def get_expr_matrix(self, df, enzymes, tissues):
        mat = np.empty((len(enzymes), len(tissues)))
        for i in range(0, len(enzymes)):
            for j in range(0, len(tissues)):
                enz_bool = list(df['enzyme'] == enzymes[i])
                ts_bool = list(df['tissue'] == tissues[j])

                row_idx = [enz_bool[idx] and ts_bool[idx] for idx in range(0, len(enz_bool))]
                if np.any(row_idx):
                    values = df.iloc[row_idx]['tpm']
                    mat[i, j] = np.min(values)
                else:
                    mat[i, j] = 0  # np.nan
        return mat

    def plot_expr_heatmap(self):

        # read file with tissue-specific expression
        p = Path(self.expr_file)
        if not p.exists():
            print("Expression data file not found")
            exit(1)
        else:
            expr_df = pd.read_csv(self.tissue_expr_file)

        expr_mat = self.get_expr_matrix(expr_df, self.enzymes, list(self.tissue_dict.values()))
        df_plot = pd.DataFrame(expr_mat, columns=list(self.tissue_dict.values()), index=self.enzymes)

        plt.rcParams['font.family'] = "Arial"
        plt.rcParams['font.size'] = 20
        kws = dict(cbar_kws=dict(orientation='horizontal', ticks=[-3, -2, -1, 0, 1, 2, 3, 4]))
        ax = sb.clustermap(df_plot, cmap="crest", figsize=(15, 7), z_score=0,
                           row_cluster=False, col_cluster=True, **kws)

        x0, _y0, _w, _h = ax.cbar_pos
        ax.ax_cbar.set_position([0.9, 0.15, ax.ax_row_dendrogram.get_position().width, 0.02])
        ax.ax_cbar.set_title("row Z-score tpm")
        ax.ax_cbar.tick_params(axis='x', length=10)
        for spine in ax.ax_cbar.spines:
            ax.ax_cbar.spines[spine].set_color('crimson')
            ax.ax_cbar.spines[spine].set_linewidth(2)

        ax.ax_row_dendrogram.set_visible(False)
        ax.ax_col_dendrogram.set_visible(False)
        ax.figure.tight_layout()
        plt.savefig("tissue_gene_expression_z_score", dpi=300, bbox_inches='tight')


def cm2inch(cm):
    return cm/2.54
