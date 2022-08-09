#!/usr/bin/python3
import re
import requests
from matplotlib import pyplot as plt
import seaborn as sb
from zipfile import ZipFile
import pandas as pd
from xml.dom import minidom
import numpy as np
import csv
from pathlib import Path


class TissueGeneExpression:

    def __init__(self, expr_dir, org_name="Arabidopsis thaliana", org_id_plantcyc="ARA"):
        self.expr_dir = expr_dir
        self.gene_file_name = '../../pw_genes.csv'
        self.tissue_expr_file = '../../gene_expr_per_tissue.csv'
        self.pathways_uniprot = [
            "pentose phosphate pathway",
            "tricarboxylic acid cycle",
            "glycolysis",
            "photorespiration"
        ]
        self.pathways_plantcyc = [
            "PENTOSE-P-PWY",
            "PWY-5690",  # TCA cycle
            "PWY-1042",  # glycolysis
            "PWY-3781",  # aerobic respiration (cytochrome c)
            "PWY-4302"  # aerobic respiration III (alternative oxidase pathway)
        ]
        self.pathway_names_plantcyc = [
            "pentose phosphate pathway",
            "tca cycle",  # TCA cycle
            "glycolysis",  # glycolysis
            "aerobic respiration (cytochrome c)",
            "aerobic respiration (alternative oxidase)"
        ]
        self.org_name = org_name
        self.org_id_plantcyc = org_id_plantcyc
        self.tissue_dict = {
            "leaf": ["leaf", "rosette leaf [0-9]+", "rosette", "adult vascular leaf"],
            "flower": ["receptacle", "flower", "carpel", "pollen"],
            "fruit": ["fruit"],
            "floral bud": ["floral bud", "flower bud"],
            "mixed leaf and apical meristem": ["mixed leaf and apical meristem"],
            "seed": ["seed", "endosperm"],
            "root": ["root"]
        }

    def find_genes_from_pathway_uniprot(self, pw, org_name):
        if org_name == "":
            org_name = self.org_name

        base = "https://rest.uniprot.org/uniprotkb/search?query=%28reviewed%3Atrue%29%2BAND%2B%28organism_name%3A"
        fmt = "tsv"

        url = base + \
              org_name.replace(" ", "%20") + \
              "%29%2BAND%2B%28cc_pathway%3A" + pw.replace(" ", "%20") + \
              "%29&fields=gene_oln" + \
              "&format=" + fmt

        response = requests.get(url).text.split("\n")

        if len(response) > 1:
            genes = response[1:]
            print(pw + ": " + str(len(genes)) + " genes found in UniProt")
        else:
            print(pw + ": No genes found in UniProt")
            genes = None

        return genes

    def find_genes_from_pathway_plantcyc(self, pw, org_id):

        url = "https://pmn.plantcyc.org/" + org_id + "/pathway-genes?object=" + pw

        response = requests.get(url).text.split("\n")

        if len(response) > 1:
            genes = list(set([s.split("\t")[0] for s in response[3:] if s != ""]))
            print(pw + ": " + str(len(genes)) + " genes found in PlantCyc")
        else:
            print(pw + ": No genes found in PlantCyc")
            genes = None

        return genes

    def find_genes_for_pathways(self, pathways=None, org=None, database="plantcyc"):

        if org is None and database == "plantcyc":
            org = self.org_id_plantcyc
        elif org is None:
            org = self.org_name

        if pathways is None:
            pathways = self.pathways_plantcyc

        genes = []
        for p in pathways:
            if database == "plantcyc":
                genes.append(self.find_genes_from_pathway_plantcyc(p, org))
            else:
                genes.append(self.find_genes_from_pathway_uniprot(p, org))
        return genes

    def create_expr_df(self, pathways=None, org=None, database="plantcyc"):

        print("Creating file for tissue-specific expression")

        if pathways is None:
            pathways = self.pathways_plantcyc
        print(pathways)
        allowed_labels = [x for xs in list(self.tissue_dict.values()) for x in xs]

        # find genes associated with pathways of interest
        print("Retrieving genes from pathways")
        p = Path(self.gene_file_name)
        genes_read = False
        if p.exists():
            print("==> loading from file")
            with open(self.gene_file_name, 'r') as gene_csv_file:
                reader = csv.reader(gene_csv_file)
                genes = []
                for row in reader:
                    genes.append(row)
                if len(genes) == len(pathways):
                    genes_read = True

        if not genes_read:
            print("==> Retrieving via request from " + database)
            genes = self.find_genes_for_pathways(pathways, org, database=database)

            with open(self.gene_file_name, 'w', newline='') as gene_csv_file:
                writer = csv.writer(gene_csv_file)
                writer.writerows(genes)

        print("Open zip with gene expression datasets")
        with ZipFile(self.expr_dir, 'r') as zipObj:
            tpms_files = [s for s in zipObj.namelist() if "tpms" in s]
            conf_files = [s for s in zipObj.namelist() if "configuration" in s]

            # loop over studies and save TPM values for genes of interest
            gene_ids = []
            tpms = []
            tpm_sds = []
            tissues = []
            categories = []
            pwys = []

            for i in range(0, len(tpms_files)):

                print(re.findall(r'^[^/]+', conf_files[i]))
                c_file = zipObj.open(conf_files[i])

                # find tissue name to abbreviation translation
                xml = minidom.parseString(c_file.read().decode())
                assay_groups = xml.getElementsByTagName('assay_group')

                ids = [s.attributes['id'].value for s in assay_groups]

                labels = [s.attributes['label'].value for s in assay_groups]
                curr_tissues = [s.split(";")[-1].strip() for s in labels]
                if curr_tissues[0].startswith("E-"):
                    curr_tissues = [s.split(";")[-2].strip() for s in labels]

                # create dictionary for sample IDs and tissue names
                label_dict = {}
                for j in range(0, len(ids)):
                    if np.any([True for s in allowed_labels if bool(re.search(r'^' + s + r'$', curr_tissues[j]))]):
                        label_dict[ids[j]] = curr_tissues[j]

                with zipObj.open(tpms_files[i]) as t_file:
                    h_line = t_file.readline().decode()
                    headers = h_line.split("\t")
                    for line in t_file:
                        for pw_idx in range(0, len(pathways)):
                            gene_id = [s for s in genes[pw_idx] if line.decode().lower().find(s.lower()) >= 0]

                            if gene_id is not None and gene_id != []:
                                row_values = line.decode().split("\t")
                                for col_idx in range(0, len(headers)):
                                    if headers[col_idx] in label_dict.keys():
                                        gene_ids.append(gene_id[0])
                                        curr_tpm_val = [float(s) for s in row_values[col_idx].split(",")]
                                        tpms.append(np.mean(curr_tpm_val))
                                        tpm_sds.append(np.std(curr_tpm_val).round(3))
                                        curr_tissue = label_dict[headers[col_idx]]
                                        tissues.append(curr_tissue)
                                        cat_found = False
                                        for k, val in self.tissue_dict.items():
                                            if not cat_found:
                                                for v in val:
                                                    if bool(re.search(r'^' + v + r'$', curr_tissue)):
                                                        categories.append(k)
                                                        cat_found = True
                                                        break
                                        if not cat_found:
                                            print(curr_tissue)
                                            print(self.tissue_dict.items())
                                            exit(0)
                                        pwys.append(pathways[pw_idx])
                print([len(gene_ids), len(tpms), len(tpm_sds), len(tissues), len(categories), len(pwys)])
            df = pd.DataFrame.from_dict({'gene_id': gene_ids,
                                         'tpm': tpms,
                                         'tpm_sd': tpm_sds,
                                         'tissue': tissues,
                                         'category': categories,
                                         'pathway': pwys})
            df.to_csv(self.tissue_expr_file, index=False, quoting=False)

    def plot_diff_expr(self, database="plantcyc"):

        if database == "plantcyc":
            pw_names = self.pathway_names_plantcyc
        else:
            pw_names = self.pathways_uniprot

        # read file with tissue-specific expression
        p = Path(self.tissue_expr_file)
        if not p.exists():
            print("Expression data file not found")
            exit(1)
        else:
            expr_df = pd.read_csv(self.tissue_expr_file)

        gt_zero_idx = [x > 0.0 for x in expr_df['tpm']]
        expr_df = expr_df.iloc[gt_zero_idx, ]
        expr_df['tpm'] = [np.log10(x) for x in expr_df['tpm']]

        plt.rcParams['font.family'] = "Arial"

        ax = sb.violinplot(data=expr_df,
                           x='pathway',
                           y='tpm',
                           hue='category',
                           palette='pastel',
                           hue_order=['root', 'leaf', 'mixed leaf and apical meristem', 'floral bud',
                                      'flower', 'fruit', 'seed'])
        ax.set_xticklabels(pw_names, rotation=10, fontsize=14)
        ax.set_ylabel(r'$\mathrm{log_{10}\ tpm}$', fontsize=14)
        ax.set_xlabel('')
        ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.1),
                  ncol=7, fancybox=True, shadow=True)
        for pos in ['right', 'top']:
            ax.spines[pos].set_visible(False)
        plt.subplots_adjust(left=0.04, bottom=0.1, right=0.98, top=0.9, wspace=0.2, hspace=0.2)
        plt.tight_layout()

        f = plt.gcf()
        f.set_size_inches(cm2inch(40), cm2inch(8))
        f.savefig('tissue_expression.png', dpi=300, bbox_inches='tight')
        # plt.show()


def cm2inch(cm):
    return cm/2.54
