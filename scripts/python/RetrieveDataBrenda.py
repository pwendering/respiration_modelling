#!/usr/bin/python3
# retrieve kcat values from BRENDA
from zeep import Client
import hashlib
import csv
import time


class RetrieveDataBrenda:
	def __init__(self):
		self.brenda_kcat_file = '../../brenda_kcats.tsv'
		self.brenda_km_file = '../../brenda_km.tsv'
		self.username = ""
		self.password = hashlib.sha256("".encode("utf-8")).hexdigest()

	def set_up_client(self):
		# try three times to get the client
		success = False
		n = 0
		client = None
		while n < 3 and not success:
			try:
				wsdl = "https://www.brenda-enzymes.org/soap/brenda_zeep.wsdl"
				client = Client(wsdl)
				success = True
				print("Client set up (" + str(n + 1) + " attempt(s)).")
			except:
				n += 1

		return client

	def process_brenda_output(self, entry, keys):
		line = []
		for key in keys:
			line.append(entry[key])
		return line

	def append_entry_to_file(self, csvwriter, result, keys):
		for entry in result:
			line = self.process_brenda_output(entry, keys)
			csvwriter.writerow(line)

	def get_ec_with_turnover_number(self, client):
		parameters = (self.username, self.password)
		ec_numbers = client.service.getEcNumbersFromTurnoverNumber(*parameters)
		return ec_numbers

	def get_ec_with_km_value(self, client):
		parameters = (self.username, self.password)
		ec_numbers = client.service.getEcNumbersFromKmValue(*parameters)
		return ec_numbers

	def retrieve_kcats(self):

		keys = ["ecNumber", "substrate", "organism", "turnoverNumber", "commentary"]

		# open output file
		fout = open(self.brenda_kcat_file, "w", newline='')
		writer = csv.writer(fout, delimiter="\t")

		# write the header line
		writer.writerow(["ec_number", "substrate", "organism", "kcat", "commentary"])

		client = self.set_up_client()

		if client is None:
			print("Client could not be set up")
			exit(1)

		ec_numbers = self.get_ec_with_turnover_number(client)

		for i in range(0, len(ec_numbers)):
			parameters = (self.username, self.password, "ecNumber*" + ec_numbers[i],
						  "organism*", "turnoverNumberMaximum*", "substrate*", "commentary*",
						  "turnoverNumber*", "ligandStructureId*", "literature*")
			success = False
			n = 0
			while n < 10 and not success:
				# send the request
				try:
					result_string = client.service.getTurnoverNumber(*parameters)
				except:
					time.sleep(5)
					n += 1
					print("Request attempt " + str(n) + " failed")

				if result_string != "":
					success = True

			# write the output to file
			self.append_entry_to_file(writer, result_string, keys)

			if i % 10 == 0 and i > 0:
				time.sleep(0.5)
			if i % 100 == 0 and i > 0:
				print("done with " + str(i))
				time.sleep(5)

		fout.close()

	def retrieve_km(self):

		keys = ["ecNumber", "substrate", "organism", "kmValue", "commentary"]

		# open output file
		fout = open(self.brenda_km_file, "w", newline='')
		writer = csv.writer(fout, delimiter="\t")

		# write the header line
		writer.writerow(["ec_number", "substrate", "organism", "Km", "commentary"])

		client = self.set_up_client()

		ec_numbers = self.get_ec_with_km_value(client)

		for i in range(0, len(ec_numbers)):
			parameters = (self.username, self.password, "ecNumber*" + ec_numbers[i],
						  "kmValue*", "kmValueMaximum*","substrate*", "commentary*",
						  "organism*", "ligandStructureId*", "literature*")
			success = False
			n = 0
			while n < 10 and not success:
				# send the request
				try:
					result_string = client.service.getKmValue(*parameters)
				except:
					time.sleep(5)
					n += 1
					print("Request attempt " + str(n) + " failed")

				if result_string != "":
					success = True

			# write the output to file
			self.append_entry_to_file(writer, result_string, keys)

			if i % 10 == 0 and i > 0:
				time.sleep(0.5)
			if i % 100 == 0 and i > 0:
				print("done with " + str(i))
				time.sleep(5)

		fout.close()
