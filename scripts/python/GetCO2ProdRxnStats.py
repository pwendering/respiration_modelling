#!/usr/bin/python3
import pandas as pd
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot
from collections import Counter, OrderedDict
import seaborn as sns
import re
from itertools import compress
from wordcloud import WordCloud, STOPWORDS

class GetCO2ProdRxnStats:

    def __init__(self, co2_prod_file, param_file):
        self.co2_prod = str(co2_prod_file)
        self.param_file = str(param_file)

        self.df = self.read_data_frame(self.co2_prod)
        self.param_df = self.read_data_frame(self.param_file)

        self.ec_level_1 = get_ec_level_1()
        self.ec_level_2 = get_ec_level_2()

    def read_data_frame(self, file_name):
        df = pd.read_csv(file_name, sep='\t')
        print(df.head())
        return df

    def plot_ec_distribution(self):

        viridiplantae_idx = np.where([s == 1 for s in self.df['is_viridiplantae']])
        ec_all = self.df.loc[:, 'ec_number']
        ec_viridiplantae = self.df.loc[viridiplantae_idx, 'ec_number']

        re_level_1 = r'^[0-9]+\.'
        re_level_2 = r'^[0-9]+\.[0-9]+\.'

        uniq_ec_level_1 = sorted(list(set([re.findall(re_level_1, s)[0] for s in ec_viridiplantae])))
        uniq_ec_level_2 = sorted(list(set([re.findall(re_level_2, s)[0] for s in ec_viridiplantae])))

        print(self.df.tail())

        n_per_ec_level_1 = np.zeros([1, len(uniq_ec_level_1)])
        n_per_ec_viridiplantae_level_1 = np.zeros([1, len(uniq_ec_level_1)])
        for i in range(0, len(uniq_ec_level_1)):
            ec_in_class = [s.startswith(uniq_ec_level_1[i]) for s in ec_all]
            n_per_ec_level_1[0, i] = np.size(np.where(ec_in_class))
            ec_in_class = [s.startswith(uniq_ec_level_1[i]) for s in ec_viridiplantae]
            n_per_ec_viridiplantae_level_1[0, i] = np.size(np.where(ec_in_class))

        n_per_ec_level_2 = np.zeros([1, len(uniq_ec_level_2)])
        n_per_ec_viridiplantae_level_2 = np.zeros([1, len(uniq_ec_level_2)])
        for i in range(0, len(uniq_ec_level_2)):
            ec_in_class = [s.startswith(uniq_ec_level_2[i]) for s in ec_all]
            n_per_ec_level_2[0, i] = np.size(np.where(ec_in_class))
            ec_in_class = [s.startswith(uniq_ec_level_2[i]) for s in ec_viridiplantae]
            n_per_ec_viridiplantae_level_2[0, i] = np.size(np.where(ec_in_class))

        y_tick_level_1 = [self.ec_level_1[ec] for ec in uniq_ec_level_1]
        y_tick_level_2 = [self.ec_level_2[ec] for ec in uniq_ec_level_2]

        x_tick_pos_l1 = np.arange(len(y_tick_level_1))
        x_tick_pos_l2 = np.arange(len(y_tick_level_2))

        print(n_per_ec_level_1)

        pyplot.rcParams['font.family'] = "Arial"
        pyplot.rcParams["figure.figsize"] = [7.50, 3.50]

        fig, (ax1, ax2) = pyplot.subplots(1, 2)
        for pos in ['right', 'top', 'left']:
            ax1.spines[pos].set_visible(False)
            ax2.spines[pos].set_visible(False)

        # plot level 1
        markerline, stemlines, baseline = ax1.stem(x_tick_pos_l1, n_per_ec_level_1[0],
                                                   label='all',
                                                   orientation='horizontal')
        pyplot.setp(markerline, 'color', (0.2, 0.4, 0.6, 1))
        pyplot.setp(stemlines, 'color', (0.2, 0.4, 0.6, 1))
        pyplot.setp(baseline, 'color', (0.6, 0.6, 0.6, 0.6))

        markerline, stemlines, baseline = ax1.stem(x_tick_pos_l1, n_per_ec_viridiplantae_level_1[0],
                                                   label='Viridiplantae',
                                                   orientation='horizontal')
        pyplot.setp(markerline, 'color', (0.4, 0.1, 0.2, 1))
        pyplot.setp(stemlines, 'color', (0.4, 0.1, 0.2, 1))
        pyplot.setp(baseline, 'color', (0.6, 0.6, 0.6, 0.6))

        ax1.set_yticks(x_tick_pos_l1, y_tick_level_1)
        ax1.legend(frameon=False)
        ax1.tick_params(axis="x", length=5)
        ax1.tick_params(axis="y", length=0)
        ax1.invert_yaxis()

        # plot level 2
        markerline, stemlines, baseline = ax2.stem(x_tick_pos_l2, n_per_ec_level_2[0],
                                                   label='all',
                                                   orientation='horizontal')
        pyplot.setp(markerline, 'color', (0.2, 0.4, 0.6, 1))
        pyplot.setp(stemlines, 'color', (0.2, 0.4, 0.6, 1))
        pyplot.setp(baseline, 'color', (0.6, 0.6, 0.6, 0.6))

        markerline, stemlines, baseline = ax2.stem(x_tick_pos_l2, n_per_ec_viridiplantae_level_2[0],
                                                   label='Viridiplantae',
                                                   orientation='horizontal')
        pyplot.setp(markerline, 'color', (0.4, 0.1, 0.2, 1))
        pyplot.setp(stemlines, 'color', (0.4, 0.1, 0.2, 1))
        pyplot.setp(baseline, 'color', (0.6, 0.6, 0.6, 0.6))

        ax2.set_yticks(x_tick_pos_l2, y_tick_level_2)
        ax2.tick_params(axis="x", length=5)
        ax2.tick_params(axis="y", length=0)
        ax2.yaxis.tick_right()
        ax2.yaxis.set_ticks_position('both')
        ax2.invert_xaxis()
        ax2.invert_yaxis()

        pyplot.subplots_adjust(left=0.15, bottom=0.11, right=0.5, top=0.9, wspace=0.2, hspace=0.2)
        pyplot.figtext(0.3, 0.01, 'Number')
        pyplot.figtext(0.01, 0.92, 'a', fontsize=14, fontweight='bold')
        pyplot.figtext(0.34, 0.92, 'b', fontsize=14, fontweight='bold')

        #f = pyplot.gcf()
        #f.set_size_inches(cm2inch(30), cm2inch(10))
        #f.savefig('ec_distribution.png', dpi=300)
        pyplot.show()

    def write_inhibitors(self):
        viridiplantae_idx = [s == 1 for s in self.df['is_viridiplantae']]
        df_viridiplantae = self.df.loc[viridiplantae_idx, ['ec_number', 'inhibitor']].copy()
        inh_per_ec_class_viridiplantae = []
        top_ten_inh_per_class = []
        for i in range(1, 8):
            ec_class_idx = [str(s).startswith(str(i)) for s in df_viridiplantae['ec_number']]
            tmp_inhibitors = df_viridiplantae.loc[ec_class_idx, 'inhibitor']
            tmp_inhibitors_all = tmp_inhibitors.str.split('|').explode()
            keep_idx = [s not in ['more'] for s in tmp_inhibitors_all]
            tmp_inhibitors_all = tmp_inhibitors_all[keep_idx].dropna()
            tmp_inhibitors_uniq = tmp_inhibitors_all.unique()

            # count the number of occurrences for unique inhibitors
            c = Counter(tmp_inhibitors_all)
            top_ten_inh_per_class.append(c.most_common(10))
            inh_per_ec_class_viridiplantae.append(tmp_inhibitors_uniq)

        # write results to file
        c_labels = ["EC" + str(i) for i in range(1, 8)]
        inh_per_class_df = pd.DataFrame(inh_per_ec_class_viridiplantae, c_labels).T
        inh_per_class_df.to_csv('inhibitors_per_ec_class.tsv', sep="\t")

        # write top ten inhibitors to file
        top_ten_single_string = []
        for i in range(0, 7):
            data = [n for (s, n) in top_ten_inh_per_class[i]]
            labels = [s for (s, n) in top_ten_inh_per_class[i]]
            top_ten_single_string.append([labels[i] + " (" + str(data[i]) + ")" for i in range(0, len(labels))])

        top_ten_inh_df = pd.DataFrame(top_ten_single_string, c_labels).T
        top_ten_inh_df.to_csv('top_ten_inhibitors_per_ec_class.tsv', sep="\t")

    def write_activators(self):
        viridiplantae_idx = [s == 1 for s in self.df['is_viridiplantae']]
        df_viridiplantae = self.df.loc[viridiplantae_idx, ['ec_number', 'activator']].copy()
        act_per_ec_class_viridiplantae = []
        top_ten_act_per_class = []
        for i in range(1, 8):
            ec_class_idx = [str(s).startswith(str(i)) for s in df_viridiplantae['ec_number']]
            tmp_activators = df_viridiplantae.loc[ec_class_idx, 'activator']
            tmp_activators_all = tmp_activators.str.split('|').explode()
            keep_idx = [s not in ['more'] for s in tmp_activators_all]
            tmp_activators_all = tmp_activators_all[keep_idx].dropna()
            tmp_activators_uniq = tmp_activators_all.unique()

            # count the number of occurrences for unique activators
            c = Counter(tmp_activators_all)
            top_ten_act_per_class.append(c.most_common(10))
            act_per_ec_class_viridiplantae.append(tmp_activators_uniq)

        # write results to file
        c_labels = ["EC" + str(i) for i in range(1, 8)]
        act_per_class_df = pd.DataFrame(act_per_ec_class_viridiplantae, c_labels).T
        act_per_class_df.to_csv('activators_per_ec_class.tsv', sep="\t")

        # write top ten activators to file
        top_ten_single_string = []
        for i in range(0, 7):
            data = [n for (s, n) in top_ten_act_per_class[i]]
            labels = [s for (s, n) in top_ten_act_per_class[i]]
            top_ten_single_string.append([labels[i] + " (" + str(data[i]) + ")" for i in range(0, len(labels))])

        top_ten_act_df = pd.DataFrame(top_ten_single_string, c_labels).T
        top_ten_act_df.to_csv('top_ten_activators_per_ec_class.tsv', sep="\t")

    def plot_kinetic_param(self, ec_level=1, param_type="kcat"):

        if ec_level == 1:
            ec_pattern = r'^[0-9]+\.'
        elif ec_level == 2:
            ec_pattern = r'^[0-9]+\.[0-9]+\.'

        ec_class = []
        is_vir = []
        values = []

        for i in range(0, len(self.param_df.index)):
            ec = self.param_df.loc[i, 'ec_number']
            ec_full = re.search(r'([0-9]+\.){3}[0-9]+', ec)
            is_ec = pd.Series.any(self.df['ec_number'] == ec)
            is_sub = re.search(r'(^|\|)CO2($|\|)', self.param_df.loc[i, 'substrate']) != None
            is_valid_value = self.param_df.loc[i, param_type] > 0
            if is_ec and ec_full and is_valid_value and not is_sub:
                ec_class.append(re.findall(ec_pattern, self.param_df.loc[i, 'ec_number'])[0])
                values.append(np.log10(self.param_df.loc[i, param_type]))
                if self.param_df.loc[i, 'lineage'].find("Viridiplantae") > 0:
                    is_vir.append("Viridiplantae")
                else:
                    is_vir.append("not Viridiplantae")

        param_df_filter = pd.DataFrame({
            param_type: values,
            'ec_class': ec_class,
            'is_vir': is_vir})

        uniq_ec_classes = sorted(param_df_filter['ec_class'].unique())
        if ec_level == 1:
            x_tick_labels = [self.ec_level_1[ec] for ec in uniq_ec_classes]
        else:
            x_tick_labels = [self.ec_level_2[ec] for ec in uniq_ec_classes]

        pyplot.rcParams['font.family'] = "Arial"
        pyplot.rcParams['font.style'] = "normal"

        ax = sns.violinplot(data=param_df_filter, x="ec_class", y=param_type,
                            hue="is_vir",
                            split=True,
                            scale='area',
                            palette='pastel')
        ax.set_xticks(range(0, len(x_tick_labels)))
        ax.set_xticklabels(labels=x_tick_labels)

        if ec_level == 1:
            ax.set_xticklabels(ax.get_xticklabels(),rotation = 30)

        if param_type == "kcat":
            ax.set_ylabel(r'$\mathrm{\log_{10}\ k_{cat}\ [s^{-1}]}$', fontsize=14)
        elif param_type == "Km":
            ax.set_ylabel(r'$\mathrm{\log_{10}\ K_M\ [M]}$', fontsize=14)

        ax.set_xlabel('EC class', fontsize=14)
        ax.get_legend().remove()

        pyplot.subplots_adjust(left=0.13, bottom=0.3, right=0.5, top=0.762, wspace=0.2, hspace=0.2)

        f = pyplot.gcf()
        f.set_size_inches(cm2inch(13), cm2inch(10))
        if param_type == "kcat":
            f.savefig('kcat_distribution.png', dpi=300)
        elif param_type == "Km":
            f.savefig('km_distribution.png', dpi=300)
        else:
            pyplot.show()

    def plot_location(self):
        pyplot.rcParams['font.family'] = "Arial"

        viridiplantae_idx = [s == 1 for s in self.df['is_viridiplantae']]
        df_viridiplantae = self.df.loc[viridiplantae_idx, ['ec_number', 'location']].copy()
        df_viridiplantae = df_viridiplantae.dropna()
        all_locations = [str(s).split('|') for s in df_viridiplantae['location']]
        all_locations = [sub_item for item in all_locations for sub_item in item]
        all_locations = [re.sub('\[.*\]:','',s) for s in all_locations]

        comps = ['Cytoplasm', 'Membrane', 'Mitochondrion', 'Endoplasmic reticulum',
                 'Nucleus', 'Peroxisome', 'Secreted',
                 'Chloroplast']
        all_locations_replaced = all_locations.copy()
        for i in range(0,len(all_locations_replaced)):
            found = False
            for j in range(0,len(comps)):
                if not found and (comps[j].lower() in all_locations_replaced[i].lower()):
                    all_locations_replaced[i] = comps[j]
                    found = True

        c_all = Counter(all_locations)

        c_all_rep = Counter(all_locations_replaced)
        c_all_rep = OrderedDict(c_all_rep.most_common())
        print([(i, round(c_all_rep[i] / len(all_locations_replaced) * 100.0, 2)) for i in c_all_rep])
        val = list(c_all_rep.values())
        lab = list(c_all_rep.keys())

        colors = sns.color_palette('pastel')[0:len(val)]
        ax = sns.barplot([val[i] for i in range(0, len(lab)) if val[i] > 1],
                         [lab[i] for i in range(0, len(lab)) if val[i] > 1],
                         palette=colors)
        ax.set_xlabel("Number of Enzymes")

        pyplot.savefig("barplot_loc_all.svg", format="svg", bbox_inches='tight')
        pyplot.close()

        membrane_comps = []
        for i in range(0, len(all_locations)):
            if all_locations_replaced[i].lower().find('membrane') >= 0:
                membrane_comps.append(all_locations[i])

        c_mem = Counter(membrane_comps)
        c_mem = OrderedDict(c_mem.most_common())
        val = list(c_mem.values())
        lab = list(c_mem.keys())
        colors = sns.color_palette('pastel')[0:len(val)]
        ax = sns.barplot([val[i] for i in range(0, len(lab)) if val[i] > 1],
                         [lab[i] for i in range(0, len(lab)) if val[i] > 1],
                         palette=colors)
        ax.set_xlabel("Number of Enzymes")

        pyplot.savefig("barplot_loc_membrane.svg", format="svg", bbox_inches='tight')
        pyplot.close()

def get_ec_level_1():
    return {
        '1.': 'Oxidoreductases (1.-.-.-)',
        '2.': 'Transferases (2.-.-.-)',
        '3.': 'Hydrolases (3.-.-.-)',
        '4.': 'Lyases (4.-.-.-)',
        '5.': 'Isomerases (5.-.-.-)',
        '6.': 'Ligases (6.-.-.-)',
        '7.': 'Translocases (7.-.-.-)'
    }
def get_ec_level_2():
    return {
        '1.1.':  'Acting on the CH-OH group of donors (1.1.-.-)',
        '1.13.': 'Acting on single donors with incorporation of molecular oxygen (oxygenases) (1.13.-.-)',
        '1.14.': 'Acting on paired donors, with incorporation or reduction of molecular oxygen (1.14.-.-)',
        '1.17.': 'Acting on CH or CH(2) groups (1.17.-.-)',
        '1.2.':  'Acting on the aldehyde or oxo group of donors (1.2.-.-)',
        '1.3.':  'Acting on the CH-CH group of donors (1.3.-.-)',
        '1.4.':  'Acting on the CH-NH(2) group of donors (1.4.-.-)',
        '1.5.':  'Acting on the CH-NH group of donors (1.5.-.-)',
        '2.2.':  'Transferring aldehyde or ketonic groups (2.2.-.-)',
        '2.3.':  'Acyltransferases (2.3.-.-)',
        '2.4.':  'Glycosyltransferases (2.4.-.-)',
        '2.5.':  'Transferring alkyl or aryl groups, other than methyl groups (2.5.-.-)',
        '2.7.':  'Transferring phosphorus-containing groups (2.7.-.-)',
        '3.1.':  'Acting on ester bonds (3.1.-.-)',
        '3.5.':  'Acting on carbon-nitrogen bonds, other than peptide bonds (3.5.-.-)',
        '4.1.':  'Carbon-carbon lyases (4.1.-.-)',
        '4.2.':  'Carbon-oxygen lyases (4.2.-.-)'
    }

def cm2inch(cm):
    return cm/2.54

def get_kcats_for_class(df, ec_class):
    ec_class_idx = [str(s).startswith(ec_class) for s in df['ec_number']]
    kcats = df.loc[ec_class_idx, 'kcat']
    kcats = kcats.str.split('|').explode()
    keep_idx = [float(k) > 0 for k in kcats]
    kcats = kcats[keep_idx].dropna()
    kcats = [float(k) for k in kcats]
    return kcats

