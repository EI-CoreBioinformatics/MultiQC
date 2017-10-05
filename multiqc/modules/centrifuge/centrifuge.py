from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import table, bargraph
from enum import Enum, unique
import logging
log = logging.getLogger(__name__)

from collections import OrderedDict

@unique
class TaxRank(Enum):
	ROOT = 0
	SUPERKINGDOM = 1
	KINGDOM = 2
	PHYLUM = 3
	CLASS = 4
	ORDER = 5
	FAMILY = 6
	GENUS = 7
	SPECIES = 8


mapping = {
    '0': 'first',
    '1': 'second',
    '2': 'third',
    '3': 'fourth',
    '4': 'fifth'
}


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='Centrifuge', anchor='centrifuge',
        href="https://ccb.jhu.edu/software/centrifuge/index.shtml",
        info="Summarises centrifuge reports from multiple samples")
        
        self.table_data = {}
        self.cls_lineage_data = {}
        self.els_lineage_data = {}
        self.kingdom_data = {}
        self.rank_data={}
        self.top5_data = {}

        for c_file in self.find_log_files('centrifuge'):
            s_name = self.clean_s_name(c_file['s_name'][:-11], c_file['root'])

            content = c_file['f'].splitlines()

            self.table_data[s_name], self.cls_lineage_data[s_name], self.els_lineage_data[s_name], self.kingdom_data[s_name], self.rank_data[s_name], self.top5_data[s_name] = self.parse_cf_reports(content)
            self.add_data_source(c_file, s_name)

        self.table_data = self.ignore_samples(self.table_data)
        self.cls_lineage_data = self.ignore_samples(self.cls_lineage_data)
        self.els_lineage_data = self.ignore_samples(self.els_lineage_data)
        self.kingdom_data = self.ignore_samples(self.kingdom_data)
        self.rank_data = self.ignore_samples(self.rank_data)
        self.top5_data = self.ignore_samples(self.top5_data)

        if len(self.table_data) == 0:
            log.debug("Could not find any data in {}".format(config.analysis_dir))
            raise UserWarning

        log.info("Found {} reports".format(len(self.table_data)))

        headers = OrderedDict()
        headers["total_hits"] = {
            'title': 'Total Hits',
            'description': "Total hits to DB",
            'format': '{:,.0f}'
        }
        headers["actual_taxon"] = {
            'title': 'Classified Taxon',
            'description': "Classified taxon"
        }
        headers["actual_taxon_rank"] = {
            'title': 'Cls Rank',
            'description': "Classified taxon's rank",
        }
        headers["actual_taxon_rank_id"] = {
            'title': 'Cls Rank Level',
            'description': "Classified taxon's rank level",
            'scale': 'RdYlGn',
            'min': 0,
            'max': 8,
            'format': '{:,.0f}'
        }
        headers["expected_taxon"] = {
            'title': 'Expected Taxon',
            'description': "Expected taxon"
        }
        headers["expected_taxon_rank"] = {
            'title': 'Exp Rank',
            'description': "Expected taxon's rank"
        }
        headers["expected_taxon_rank_id"] = {
            'title': 'Exp Rank Level',
            'description': "Expected taxon's rank level",
            'scale': 'RdYlGn',
            'min': 0,
            'max': 8,
            'format': '{:,.0f}'
        }
        headers["mrca_taxon"] = {
            'title': 'MRCA Taxon',
            'description': "Most recent common ancestor taxon"
        }
        headers["mrca_taxon_rank"] = {
            'title': 'MRCA Rank',
            'description': "Most recent common ancestor taxon's rank"
        }
        headers["mrca_taxon_rank_id"] = {
            'title': 'MRCA Rank Level',
            'description': "Most recent common ancestor taxon's rank level",
            'scale': 'RdYlGn',
            'min': 0,
            'max': 8,
            'format': '{:,.0f}'
        }
        tc_config = {
            'namespace': 'Centrifuge',
            'scale': 'RdYlGn',
            'save_file': True,
            'raw_data_fn': 'multiqc_centrifuge_taxatable'
        }

        self.add_section(
            name='Taxonomic classification',
            anchor='centrifuge-first',
            description='Table showing classified taxon for each sample.  The classfied taxon is typically the first category (starting from species level) to contain > 50% of total hits.  If sample has an expected taxon provided in samplesheet then this is also displayed, along with the most recent common ancestor of the classified and expected taxa.',
            helptext="Help?",
            plot=table.plot(self.table_data, headers, tc_config)
        )


        t5_headers = OrderedDict()
        t5_headers['first_name'] = {
            'title': "First"
        }
        t5_headers['first_count'] = {
            'title': "Hits %",
            'scale': 'RdYlGn',
            'min': 0,
            'max': 100.0,
            'format': '{:,.1f}'
        }

        t5_headers['second_name'] = {
            'title': "Second"
        }
        t5_headers['second_count'] = {
            'title': "Hits %",
            'scale': 'RdYlGn',
            'min': 0,
            'max': 100.0,
            'format': '{:,.1f}'
        }
        t5_headers['third_name'] = {
            'title': "Third"
        }
        t5_headers['third_count'] = {
            'title': "Hits %",
            'scale': 'RdYlGn',
            'min': 0,
            'max': 100.0,
            'format': '{:,.1f}'
        }
        t5_headers['fourth_name'] = {
            'title': "Fourth"
        }
        t5_headers['fourth_count'] = {
            'title': "Hits %",
            'scale': 'RdYlGn',
            'min': 0,
            'max': 100.0,
            'format': '{:,.1f}'
        }
        t5_headers['fifth_name'] = {
            'title': "Fifth"
        }
        t5_headers['fifth_count'] = {
            'title': "Hits %",
            'scale': 'RdYlGn',
            'min': 0.0,
            'max': 100.0,
            'format': '{:,.1f}'
        }
        t5_config = {
            'namespace': 'Centrifuge',
            'scale': 'RdYlGn',
            'save_file': True,
            'raw_data_fn': 'multiqc_centrifuge_top5'
        }
        self.add_section(
            name='Top 5 Species',
            anchor='centrifuge-second',
            description='Table showing the top 5 largest species for each sample..',
            helptext="Help?",
            plot=table.plot(self.top5_data, t5_headers, t5_config)
        )

        k_cats = OrderedDict()
        k_cats['bacteria'] = {
            'name': "Bacteria",
            'color': 'orange'
        }
        k_cats['archaea'] = {
            'name': "Archaea",
            'color': 'black'
        }
        k_cats['viruses'] = {
            'name': "Viruses",
            'color': 'blue'
        }
        k_cats['fungi'] = {
            'name': 'Fungi',
            'color': 'gray'
        }
        k_cats['viridiplantae'] = {
            'name': 'Viridiplantae',
            'color': 'green'
        }
        k_cats['metazoa'] = {
            'name': 'Metazoa',
            'color': 'red'
        }
        k_config = {
            'id': 'kingdom',
            'cpswitch': True,
            'cpswitch_c_active': False,
            'title': 'Hits in each kingdom',
            'ymin': 0

        }

        self.add_section(
            name='Kingdom classification',
            anchor='centrifuge-third',
            description='Bar plot showing hits falling into each Kingdom (or Super Kingdom).',
            helptext="Help?",
            plot=bargraph.plot(self.kingdom_data, k_cats, k_config)
        )


        cats = OrderedDict()
        cats['species'] = {
            'name': 'Species',
            'color': 'darkgreen'
        }
        cats['genus'] = {
            'name': 'Genus',
            'color': 'green'
        }
        cats['family'] = {
            'name': 'Family',
            'color': 'lightgreen'
        }
        cats['order'] = {
            'name': 'Order',
            'color': 'yellow'
        }
        cats['class'] = {
            'name': 'Class',
            'color': 'orange'
        }
        cats['phylum'] = {
            'name': 'Phylum',
            'color': 'pink'
        }
        cats['kingdom'] = {
            'name': 'Kingdom',
            'color': 'red'
        }
        cats['superkingdom'] = {
            'name': 'Domain',
            'color': 'darkred'
        }
        cats['root'] = {
            'name': 'Root',
            'color': 'gray'
        }
        lineage_config = {
            'id': 'cls_lineage',
            'cpswitch': True,
            'cpswitch_c_active': False,
            'title': 'Hits at each rank in taxons lineage',
            'ymin': 0

        }

        self.add_section(
            name='Rank hits',
            anchor='centrifuge-fourth',
            description='Bar plot showing the rank at which centrifuge hits were classified.',
            helptext="Help?",
            plot=bargraph.plot(self.rank_data, cats, lineage_config)
        )

        self.add_section(
            name='Classified taxons lineage',
            anchor='centrifuge-fifth',
            description='Bar plot showing the classified taxons lineage (largest taxa at each taxonomic rank incorperating the classified taxon).',
            helptext="Help?",
            plot=bargraph.plot(self.cls_lineage_data, cats, lineage_config)
        )

        self.add_section(
            name='Expected taxons lineage',
            anchor='centrifuge-sixth',
            description='Bar plot showing the expected taxons lineage (largest taxa at each taxonomic rank incoperating the expected taxon).',
            helptext="Help?",
            plot=bargraph.plot(self.els_lineage_data, cats, lineage_config)
        )
        



        #self.write_data_file(self.cf_data, 'multiqc_centrifuge')

    def parse_taxon(self, taxon_str, rank=None):

        parts = taxon_str.strip().split('(')
        name = parts[0][:-1].strip().lower()
        name = name[0].capitalize() + name[1:]

        if rank:
            id = parts[1][3:-1]
            rank_name = rank.name.lower()
            rank_id = rank.value
            return name, id, rank_name, rank_id
        else:
            parts2 = parts[1].split(';')
            id_part = parts2[0]
            id = id_part.split(':')[1].strip()
            rank_part = parts2[1]
            rank = rank_part.split(':')[1][:-1].strip()
            rank = rank.strip()
            rank_id = TaxRank[rank.upper()].value
            return name, id, rank, rank_id



    def parse_cf_reports(self, content):

        table_data = {}

        exp_present = False

        for line in content:

            if line.startswith("Classified Taxon:"):
                name, id, rank, rank_id = self.parse_taxon(line[18:].strip())
                table_data["actual_taxon"] = name + " (" + id + ")"
                table_data["actual_taxon_rank"] = rank
                table_data["actual_taxon_rank_id"] = rank_id
            elif line.startswith("Expected Taxon:"):
                name, id, rank, rank_id = self.parse_taxon(line[16:].strip())
                table_data["expected_taxon"] = name + " (" + id + ")"
                table_data["expected_taxon_rank"] = rank
                table_data["expected_taxon_rank_id"] = rank_id
            elif line.startswith("MRCA Taxon:"):
                name, id, rank, rank_id = self.parse_taxon(line[11:].strip())
                table_data["mrca_taxon"] = name + " (" + id + ")"
                table_data["mrca_taxon_rank"] = rank
                table_data["mrca_taxon_rank_id"] = rank_id
            elif line.startswith("Expected Taxon ID provided"):
                exp_present = True
            elif line.startswith("Total hits"):
                table_data["total_hits"] = int(line.split(':')[1].strip())

        cls_lineage_data = self.extract_lineage(content, "Classified Taxon's Lineage")
        kingdom_data = self.extract_kingdom(content)
        rank_data = self.extract_lineage(content, "Rank Perc")
        t5_data = self.extract_t5(content, "Top 5 ")

        els_lineage_data = {}
        if exp_present:
            els_lineage_data = self.extract_lineage(content, "Expected Taxon's Lineage")
        else:
            for i, r in enumerate(TaxRank):
                els_lineage_data[r.name.lower()] = 100.0 if i == 0 else 0.0

        return table_data, cls_lineage_data, els_lineage_data, kingdom_data, rank_data, t5_data



    def extract_t5(self, content, header):

        desc_str = ""
        perc_str = ""
        count_str = ""
        in_section = 0
        t5_data = {}

        for l in content:
            line = l.strip()

            if in_section == 1:
                desc_str = line
                in_section = 2
            elif in_section == 2:
                perc_str = line
                in_section = 3
            elif in_section == 3:
                count_str = line
                in_section = 0

            if line.startswith(header):
                in_section = 1

        desc_parts = desc_str.split('\t')
        perc_parts = perc_str.split('\t')
        count_parts = count_str.split('\t')

        if len(desc_parts) == 0 or len(perc_parts) == 0 or len(count_parts) == 0:
            raise ValueError("Error: Could not find lineage with header: " + header + ".  Invalid centrifuge summary file.")

        if len(desc_parts) != len(perc_parts) or len(desc_parts) != len(count_parts) or len(desc_parts) != 5:
            raise ValueError("Error: Lineage with header: " + header + "; has different number of header and value columns.")

        for i in range(len(desc_parts)):
            h = desc_parts[i].strip()
            #print(h)
            name, id, rank, rank_id = self.parse_taxon(h, rank=TaxRank.SPECIES)
            t5_data[mapping[str(i)] + "_name"] = h
            t5_data[mapping[str(i)] + "_count"] = float(perc_parts[i])

        return t5_data


    def extract_lineage(self, content, header):

        desc_str = ""
        perc_str = ""
        count_str = ""
        in_section = 0
        lineage_data = {}

        for l in content:
            line = l.strip()

            if in_section == 1:
                desc_str = line
                in_section = 2
            elif in_section == 2:
                perc_str = line
                in_section = 3
            elif in_section == 3:
                count_str = line
                in_section = 0

            if line.startswith(header):
                in_section = 1

        desc_parts = desc_str.split('\t')
        perc_parts = perc_str.split('\t')
        count_parts = count_str.split('\t')

        if len(desc_parts) == 0 or len(perc_parts) == 0 or len(count_parts) == 0:
            raise ValueError("Error: Could not find lineage with header: " + header + ".  Invalid centrifuge summary file.")

        if len(desc_parts) != len(perc_parts) or len(desc_parts) != len(count_parts):
            raise ValueError("Error: Lineage with header: " + header + "; has different number of header and value columns.")

        for i in range(len(desc_parts)):
            h = desc_parts[i].strip()
            #print(h)
            name, id, rank, rank_id = self.parse_taxon(h)
            lineage_data[TaxRank[rank.upper()].name.lower()] = int(count_parts[i]) if i == 0 else (int(count_parts[i]) - int(count_parts[i - 1]))

        return lineage_data


    def extract_kingdom(self, content):

        desc_str = ""
        perc_str = ""
        count_str = ""
        in_section = 0

        for l in content:
            line = l.strip()

            if in_section == 1:
                desc_str = line
                in_section = 2
            elif in_section == 2:
                perc_str = line
                in_section = 3
            elif in_section == 3:
                count_str = line
                in_section = 0

            if line.startswith("Kingdom Perc"):
                in_section = 1

        desc_parts = desc_str.split('\t')
        perc_parts = perc_str.split('\t')
        count_parts = count_str.split('\t')
        kingdom_data = {}

        if len(desc_parts) == 0 or len(perc_parts) == 0 or len(count_parts) == 0:
            raise ValueError("Error: Could not find kingdom percentages.  Invalid centrifuge summary file.")

        if len(perc_parts) != len(desc_parts) or len(desc_parts) != len(count_parts):
            raise ValueError("Error: Kingdom percentages has different number of header and percentage columns.")

        for i in range(len(desc_parts)):
            h = desc_parts[i].strip()
            #print(h)
            name, id, rank, rank_id = self.parse_taxon(h, rank=TaxRank.KINGDOM)
            kingdom_data[name.lower()] = int(count_parts[i])

        return kingdom_data