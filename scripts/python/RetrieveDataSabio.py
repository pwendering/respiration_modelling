#!/usr/bin/python3
# retrieve kinetic parameters from SABIO-RK (http://sabiork.h-its.org/)

import requests
import math


class RetrieveDataSabio:

	def __init__(self):
		self.sabio_outfile = '../../sabio_params.tsv'

	def get_org_entries(self):
		""" obtain all entry IDs for a given organism """
		# URL
		ENTRYID_QUERY_URL = 'http://sabiork.h-its.org/sabioRestWebServices/searchKineticLaws/entryIDs'

		query = {'format': 'txt', 'q': '*'}

		# make GET request
		request = requests.get(ENTRYID_QUERY_URL, params=query)
		request.raise_for_status()  # raise if 404 error

		# extract entry IDs
		entry_ids = [int(x) for x in request.text.strip().split('\n')]
		print('> %d matching entries found.' % len(entry_ids))

		return entry_ids

	def get_params(self, entry_ids):
		""" obtain kinetic parameters for each given entry ID """

		# URL
		PARAM_QUERY_URL = 'http://sabiork.h-its.org/entry/exportToExcelCustomizable'

		# construct request
		data_field = {'entryIDs[]': entry_ids}
		query = {'format': 'tsv',
				 'fields[]': ['ECNumber', 'Substrate', 'Organism', 'Parameter', 'Temperature', 'EnzymeType']}

		# make POST request
		request = requests.post(PARAM_QUERY_URL, params=query, data=data_field)
		request.raise_for_status()

		return request.text

	def write_output(self, req_output, outfile_name):
		"""
		write the output of the parameter request to a tab-separated file
		"""
		f = open(outfile_name, "a+", encoding="utf-8")
		f.write(req_output)
		f.close()

	def retrieve_kinetic_params(self):

		entry_ids = self.get_org_entries()

		# separate data into bins of 1000 for faster processing
		bin_size = 1000
		round_val = round(math.log(bin_size, 10))
		i = 1
		while i*bin_size <= round(len(entry_ids), -round_val):
			bin_start = (i - 1) * bin_size
			bin_end = i * bin_size
			print("Processing IDs " + str(bin_start) + "--" + str(bin_end))
			tmp_entry_ids = entry_ids[bin_start:bin_end]
			response_text = self.get_params(tmp_entry_ids).encode("utf-8", errors='replace')

			self.write_output(response_text.decode("utf-8"), self.sabio_outfile)
			i += 1
