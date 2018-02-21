#!/usr/bin/env python

""" MultiQC module to parse output from KAT """
import logging
from collections import OrderedDict

from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import table

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
	def __init__(self):
		# Initialise the parent object
		super(MultiqcModule, self).__init__(name='KAT', anchor='kat',
											href="https://github.com/TGAC/KAT",
											info="is an toolkit for analysing sequencing data via its k-mer spectra.")

		# Find and load any KAT dist analysis reports
		self.kat_data = dict()
		for c_file in self.find_log_files('kat'):
			s_name = self.clean_s_name(c_file['s_name'][:-4], c_file['root'])
			content = c_file['f'].splitlines()
			self.kat_data[s_name] = self.parse_kat_report(content)

			# Filter to strip out ignored sample names
		self.kat_data = self.ignore_samples(self.kat_data)

		if len(self.kat_data) == 0:
			raise UserWarning

		log.info("Found {} reports".format(len(self.kat_data)))

		# Write parsed report data to a file
		self.write_data_file(self.kat_data, 'multiqc_kat')

		headers = OrderedDict()
		headers["kmer_peaks"] = {
			'title': '# of Kmer Peaks',
			'description': "Number of peaks identified in the K-mer spectra",
			'format': '{:,.0f}',
			'custom_style': 'width=10%'
		}
		headers["gc_peaks"] = {
			'title': '# of GC Peaks',
			'description': "Number of peaks identified in the GC distribution",
			'format': '{:,.0f}',
			'custom_style': 'width=10%'
		}
		headers["est_genome_size"] = {
			'title': 'Est. genome Size',
			'description': "Estimated Genome Size based on K-mer spectra",
			'format': '{:,.0f}',
			'custom_style': 'width=20%'
		}
		headers["mean_kmer_freq"] = {
			'title': 'Mean K-mer Freq.',
			'description': "Mean K-mer Frequency, provides an estimate of sequencing coverage",
			'format': '{:,.0f}',
			'custom_style': 'width=20%',
			'suffix': 'x'
		}

		kat_config = {
			'namespace': 'KAT',
			'scale': 'RdYlGn',
			'save_file': True,
			'raw_data_fn': 'mqc_kat_table'
		}

		# Basic Stats Table
		self.add_section(
			name='KAT Distribution Analysis',
			anchor='kat-first',
			description='Table showing k-mer coverage distributions and if available GC distributions',
			helptext="This table can give a quick idea of potential contaminants that can be identified via unexpected numbers of k-mer or gc peaks in the data",
			plot=table.plot(self.kat_data, headers, kat_config)
		)


	def parse_kat_report(self, content):
		table_data = {}
		i = 0
		while i < len(content):
			line = content[i].strip()
			if line.startswith('K-mer frequency spectra statistics'):
				table_data['kmer_peaks'] = int(content[i + 3].split(':')[1].strip())
				table_data['mean_kmer_freq'] = int(content[i + 6].split(':')[1].strip()[:-1])
				i += 7
			if line.startswith('Estimated genome size'):
				table_data['est_genome_size'] = int(float(content[i].split(':')[1].strip()[:-4]) * 1000000.0)
			if line.startswith('GC distribution statistics'):
				table_data['gc_peaks'] = int(content[i + 3].split(':')[1].strip())
				i += 4
			i+=1

		return table_data
