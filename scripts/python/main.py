import sys
from GetCO2ProdRxnStats import GetCO2ProdRxnStats
from RetrieveDataBrenda import RetrieveDataBrenda
from RetrieveDataSabio import RetrieveDataSabio
from CombineBrendaSabio import CombineBrendaSabio
from TissueProteinExpression import TissueProteinExpression
from TissueGeneExpression import TissueGeneExpression

def main(co2_prod_file):

    # retrieve enzyme data from SABIO-RK
    # rds = RetrieveDataSabio()
    # rds.retrieve_kinetic_params()

    # Generate combined data for KM values from BRENDA and SABIO-RK
    # rdb = RetrieveDataBrenda()
    # rdb.retrieve_km()
    # cb = CombineBrendaSabio('../../data/brenda_km.tsv', '../../data/sabio_params.tsv', 'Km')
    # cb.add_lineagues()

    # Generate combined data for kcat values from BRENDA and SABIO-RK
    # rdb = RetrieveDataBrenda()
    # rdb.retrieve_kcats()
    # cb = CombineBrendaSabio('../../data/brenda_kcats.tsv', '../../data/sabio_params.tsv', 'kcat')
    # cb.add_lineagues()

    # plot related data for CO2 producing reactions
    stats = GetCO2ProdRxnStats(co2_prod_file, 'brenda_sabio_kcat_comb_lineage.tsv')
    stats.plot_ec_distribution()
    stats.plot_location()
    stats.write_inhibitors()
    stats.write_activators()
    stats.plot_kinetic_param(ec_level=1, param_type="kcat")
    stats.plot_kinetic_param(ec_level=2, param_type="kcat")

    # KM plot
    stats = GetCO2ProdRxnStats(co2_prod_file, '../../data/brenda_sabio_Km_comb_lineage.tsv')
    stats.plot_kinetic_param(ec_level=1, param_type="Km")
    stats.plot_kinetic_param(ec_level=2, param_type="Km")

    # Tissue-specific gene expression
    # gene expression
    teg = TissueGeneExpression('../../data/gene_tissue_atlas.csv', '../../data/tissue_name_dict.csv')
    teg.create_expr_df()
    teg.plot_expr_heatmap()

    # Tissue-specific protein expression
    tep = TissueProteinExpression('../../data/prot_tissue_atlas.csv', '../../data/tissue_name_dict.csv')
    tep.create_expr_df()
    tep.plot_diff_expr()
    tep.plot_expr_heatmap()

if __name__ == '__main__':
    nargin = len(sys.argv)
    rxn_data_file = "../../data/full-enzyme-data-table.tsv"
    main(rxn_data_file)

