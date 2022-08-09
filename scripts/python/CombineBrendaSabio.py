#!/usr/bin/python3
import requests
from pathlib import Path
import csv
import re


class CombineBrendaSabio:

    def __init__(self, brenda_file, sabio_file, param_type):
        self.brenda_file = brenda_file
        self.sabio_file = sabio_file
        self.param_type = param_type
        self.comb_file = "../../brenda_sabio_" + param_type + "_comb.tsv"
        self.param_type = param_type

    def get_lineage(self, org_name):
        base = 'https://rest.uniprot.org/uniprotkb/search?compressed=false'
        query = org_name.replace(" ", "%20")
        field = 'lineage'
        fmt = 'tsv'
        url = base + '&' + 'format=' + fmt + '&' + 'query=' + query + '&' + 'fields=' + field
        response = requests.get(url).text.split("\n")
        if len(response) > 1:
            lineage = re.sub(r' \([a-z ]+\)', '', response[1])
            return lineage.replace(",", ";")
        else:
            print("No lineage found for " + org_name)
            return org_name

    def combine_files(self):

        if self.param_type == "kcat":
            unit = 's^(-1)'
        elif self.param_type == "Km":
            unit = 'M'
        else:
            unit = ''

        comb_file = open(self.comb_file, "a+")
        writer = csv.writer(comb_file, delimiter="\t", lineterminator='\n')

        # adjust format SABIO-RK file
        with open(self.sabio_file, 'r') as sab_file:
            for line in sab_file:
                if line.find("mutant") < 0:
                    l_split = line.replace(";", "|").strip().split("\t")
                    is_param_type = bool(l_split[3] == self.param_type)
                    is_unit = bool(l_split[8] == unit)
                    is_param_num = bool(re.search(r'^[E0-9\.\-]+$', l_split[5]))
                    if is_param_type and is_unit and is_param_num:
                        writer.writerow([l_split[i] for i in [0, 1, 2, 5]])

        # adjust format BRENDA file
        with open(self.brenda_file, 'r') as brenda_file:
            brenda_file.readline()
            for line in brenda_file:
                if line.find("mutant") < 0:
                    l_split = line.strip().split("\t")
                    writer.writerow([l_split[i] for i in [0, 1, 2, 3]])

        comb_file.close()

    def add_lineagues(self):

        path = Path(self.comb_file)

        if not path.is_file():
            print("Combined file not does exist yet, combining files...")
            self.combine_files()
            print("Done.")

        outfile = self.comb_file.replace('.', '_lineage.')

        # initialize lineage dictionary
        lineage_dict = dict()

        # open output file
        out_file = open(outfile, "w", newline='')
        writer = csv.writer(out_file, delimiter="\t")

        writer.writerow(["ec_number", "substrate", "lineage", self.param_type])

        # go through combined file
        with open(self.comb_file, "r") as cf:
            for ln in cf:
                line_split = ln.strip().split("\t")
                org_name = line_split[2]
                if org_name in lineage_dict:
                    lineage = lineage_dict[org_name]
                else:
                    lineage = self.get_lineage(org_name)
                    lineage_dict[org_name] = lineage

                if lineage != "":
                    line_split[2] = lineage

                writer.writerow(line_split)

        out_file.close()
