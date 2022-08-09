#!/usr/bin/python3
# retrieve activators, inhibitors, and turnover values from BRENDA for a 
# list of EC numbers
import sys
from zeep import Client
import hashlib
import csv
import os.path


def process_brenda_output(entry, key, unique=True):
	if unique:
		return('|'.join(set([s[key] for s in entry])))
	else:
		return('|'.join([s[key] for s in entry]))


def setup_client():
	# try three times to get the client
	success = False
	client = []
	n=0
	while n<3 and not success:
		try:
			wsdl = "https://www.brenda-enzymes.org/soap/brenda_zeep.wsdl"
			client = Client(wsdl)
			success = True
		except:
			n += 1
	return success, client


def assert_non_empty(item):
	return item is not None and item != [] and item != ""


def assert_invalid_entry(comment):
	return any(s in comment for s in ['recombinant', 'mutant'])


def main(ec_file, outfile_name, username, password):

	success, client = setup_client()
	if not success:
		sys.exit("Client could not be set up, exiting")

	print("write output to "+outfile_name+"\n")

	# open output file
	if os.path.exists(outfile_name):
		write_header = False
		with open(outfile_name) as f:
			for line in f:
				pass
			last_ec = line.split('\t')[0]
		fout = open(outfile_name, 'a', newline='')

	else:
		fout = open(outfile_name, 'w', newline='')
		last_ec = ""
		write_header = False

	writer = csv.writer(fout, delimiter='\t')

	# write the header line
	if write_header:
		writer.writerow(["ec_number", "inhibitor", "activator", "turnover_number", "substrate"])

	# open input file
	fin = open(ec_file, 'r', newline='')

	for ec in fin:
		
		ec_query = ec.strip()
		
		if last_ec != "" and ec_query != last_ec:
			continue
		elif ec_query == last_ec:
			last_ec = ""
			continue

		row = [ec_query + '\t']

		# inhibitors
		parameters = (username, password, "ecNumber*"+ec_query, "inhibitor*",
		"commentary*", "organism*", "ligandStructureId*", "literature*")
		result_string = client.service.getInhibitors(*parameters)
		if result_string != "" and result_string != []:
			inhibitors = process_brenda_output(result_string, "inhibitor")
		else:
			inhibitors = ""

		# activators
		parameters = (username, password, "ecNumber*"+ec_query, "activatingCompound*",
		"commentary*", "organism*", "ligandStructureId*", "literature*")
		result_string = client.service.getActivatingCompound(*parameters)
		if result_string != "" and result_string != []:
			activators = process_brenda_output(result_string, "activatingCompound")
		else:
			activators = ""

		# turnover numbers
		parameters = (username, password, "ecNumber*"+ec_query, "turnoverNumber*",
		"turnoverNumberMaximum*", "substrate*", "commentary*",
		"organism*", "ligandStructureId*", "literature*")
		result_string = client.service.getTurnoverNumber(*parameters)

		if result_string != "" and result_string != []:
			# remove entries that are associated with mutant/recombinant enzymes
			comments = [entry["commentary"] for entry in result_string]
			to_remove = []
			for i in range(0, len(comments)):
				if assert_non_empty(comments[i]) and assert_invalid_entry(comments[i]):
					to_remove.append(i)

			to_remove = sorted(to_remove, reverse=True)
			for i in to_remove:
				result_string.pop(i)

		kcats = process_brenda_output(result_string, "turnoverNumber", unique=False)
		substrates = process_brenda_output(result_string, "substrate", unique=False)

		# append row to file
		row = [ec_query, inhibitors, activators, kcats, substrates]
		writer.writerow(row)
		
	fout.close()
	fin.close()


if __name__ == "__main__":
	narg = len(sys.argv)
	if narg < 2:
		print("USAGE: getKMBRENDA.py ec_number_file outFileName(opt) username(opt) password(opt)")
	else:

		ec_number_file = sys.argv[1]
		
		if narg < 3:
			outFileName = "enzyme-data.tsv"
		else:
			outFileName = sys.argv[2]

		if narg < 4:
			username = ""
			password = ""
		else:
			username = sys.argv[3]
			password = sys.argv[4]
	
		password = hashlib.sha256(password.encode("utf-8")).hexdigest()
	
		main(ec_number_file, outFileName, username, password)
