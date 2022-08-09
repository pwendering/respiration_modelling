#!/usr/bin/bash
########################################################################
# Classify CO2-producing reactions known to occur in photoautotrophic 
# organisms.
#
# INPUT
#
# SmartTable downloaded from BioCyc that contains all CO2-
# producing reactions in the database.
# Matched columns are
# * Matches (all CARBON-DIOXIDE)
# * EC number
# * MetaNetX ID
# * Spontaneous?
# * Pathway
# * BioCyc reaction ID
# * Data source
# * Regulated by
# * Requirements
#
# OUTPUT
# * rxn_table.txt   (list of reactions that produce CO2 with associated 
#                   properties EC number, MetaNetX ID, isSpontaneous,
#                   and Pathway)
# * uniq_ec_rxn.txt (list of unique EC numbers with associated reactions
#                   and pathways
if [[ $# == 0 ]];
then
  filename="../../data/co2-producing-reactions-biocyc.txt"
else
  filename=$1
fi

if test -f "$filename"
then
  lnr=$(cat $filename | wc -l)
  echo "Parsing BioCyc smarttable with $(($lnr-1)) lines" 
else
  echo "Input file $filename does not exist"
  exit 1
fi

if test -f "uniq_ec_rxn.txt"
then
  rm uniq_ec_rxn.txt
fi

awk -F"\t" -v OFS="\t" '{
  mnx=match($4,/MNXR[0-9]+/);
  gsub("EC-","",$3);
  gsub(/[\ ]+\/\/[\ ]+/,"|",$3);
  gsub(/[\ ]+\/\/[\ ]+/,"|",$6);
  print $1,$3,substr($4,mnx,RLENGTH),sub(/T/,1,$5),$6
}' <(sed '1d' $filename) > rxn_table.txt

# find EC numbers that are contained in the BRENDA file but not in the BioCyc file
ec_numbers_biocyc=$(cut -f 2 rxn_table.txt | grep -Eo "([0-9]+\.){2}[0-9]+?(\.[0-9]+)" | sort | uniq)
ec_numbers_brenda=$(grep -vwFf \
  <(cut -f 2 rxn_table.txt | grep -Eo "([0-9]+\.){2}[0-9]+?(\.[0-9]+)" | sort | uniq) \
  <(grep -i product ../../data/brenda_co2_prod.csv | cut -f 2 | sort | uniq))

# find unique EC numbers and associated reactions
for ec in $(echo $ec_numbers_biocyc $ec_numbers_brenda)
do
  # reaction IDs
  rxn=$(awk -v ec=$ec '$2 ~ ec {print $1}' rxn_table.txt | tr "\n" "|" | sed 's/\|$//')
  # pathway IDs
  pwy=$(awk -v ec=$ec '$2 ~ ec {print $5}' rxn_table.txt | tr "\n" "|" | sed -e 's/|\{2,\}/|/g' -e 's/|$//' -e 's/^|//')
  
  # lineages of associated organism in UniProt contains term "Viridiplantae"?
  url=$(printf '%s' "https://rest.uniprot.org/uniprotkb/search?" \
    "query=ec:$ec+AND+taxonomy_id:33090&" \
    "fields=lineage,cc_subcellular_location&" \
    "size=500&" \
    "format=tsv"
    )
    
  curl -s ${url} | sed '1d' > tmp_lineage.txt
  
  isphot=$(cat tmp_lineage.txt | wc -l)

  if [[ $isphot > 0 ]]
  then
    photBool=1
  else
    photBool=0
  fi
  
  # subcellular location
  loc=$(awk -F"\t" '/SUBCELLULAR LOCATION/ {\
    gsub(/SUBCELLULAR LOCATION: /,"",$2);\
    gsub(/\{[^;]*\}/,"",$2);\
    gsub(/[\ \.]*$/,"",$2);\
    print $2\
    }' tmp_lineage.txt |\
    tr ";" "\n" |\
    sed -e 's/^\ //' -e 's/\ $//' |\
    sort |\
    uniq |\
    tr "\n" "| " |\
    sed -e 's/^|//' -e 's/|$//')

  echo -e "$ec\t$rxn\t$pwy\t$photBool\t$loc" | tee -a uniq_ec_rxn.txt
  
  rm tmp_lineage.txt

done

# associate EC numbers with activators, inhibitors and turnover numbers
echo "Searching for activators, inhibitors and turnover numbers in BRENDA"
../python/getInhActKcatBrenda.py <(cut -f 1 uniq_ec_rxn.txt) enzyme-data.tsv

# concatenate files
paste uniq_ec_rxn.txt <(cut -f2- enzyme-data.tsv | sed '1d') > ../../data/full-enzyme-data-table.tsv

