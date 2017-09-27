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
        for f in self.find_log_files('centrifuge', filehandles=True):
            s_name = self.clean_s_name(f['s_name'][:-11], f['root'])
            self.table_data[s_name], self.cls_lineage_data[s_name], self.els_lineage_data[s_name], self.kingdom_data[s_name] = self.parse_cf_reports(f)
            self.add_data_source(f, s_name)

        self.table_data = self.ignore_samples(self.table_data)
        self.cls_lineage_data = self.ignore_samples(self.cls_lineage_data)
        self.els_lineage_data = self.ignore_samples(self.els_lineage_data)
        self.kingdom_data = self.ignore_samples(self.kingdom_data)

        if len(self.table_data) == 0:
            log.debug("Could not find any data in {}".format(config.analysis_dir))
            raise UserWarning

        log.info("Found {} reports".format(len(self.table_data)))

        headers = OrderedDict()
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
        config = {
            'namespace': 'Centrifuge',
            'scale': 'RdYlGn',
            'save_file': True,
            'raw_data_fn': 'multiqc_centrifuge_taxatable'
        }

        self.add_section(
            name='Taxonomic classification',
            anchor='centrifuge-first',
            description='Table showing classified taxon for each sample.  If expected taxon is  provided in samplesheet then this is also displayed, along with the most recent common ancestor of the classified and expected taxa.',
            helptext="Help?",
            plot=table.plot(self.table_data, headers, config)
        )

        k_cats = OrderedDict()
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
            'cpswitch': False,
            'cpswitch_c_active': False,
            'title': 'Proportion of hits in each kingdom',
            'ylab': "% hits",
            'ymax': 100,
            'ymin': 0

        }

        self.add_section(
            name='Kingdom percentages',
            anchor='centrifuge-second',
            description='Bar plot showing the proportion of hits in each kingdom.',
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
            'cpswitch': False,
            'cpswitch_c_active': False,
            'title': 'Proportion of hits in rank for taxons lineage',
            'ylab': "% hits",
            'ymax': 100,
            'ymin': 0

        }

        self.add_section(
            name='Classified taxons lineage',
            anchor='centrifuge-third',
            description='Bar plot showing the classified taxons lineage (largest taxa at each taxonomic rank).',
            helptext="Help?",
            plot=bargraph.plot(self.cls_lineage_data, cats, lineage_config)
        )

        self.add_section(
            name='Expected taxons lineage',
            anchor='centrifuge-fourth',
            description='Bar plot showing the expected taxons lineage (largest taxa at each taxonomic rank).',
            helptext="Help?",
            plot=bargraph.plot(self.els_lineage_data, cats, lineage_config)
        )
        



        #self.write_data_file(self.cf_data, 'multiqc_centrifuge')

    def parse_taxon(self, taxon_str, kingdom=False):

        parts = taxon_str.strip().split('(')
        name = parts[0][:-1].strip().lower()

        if kingdom:
            id = parts[1][3:-1]
            rank = "kingdom"
            rank_id = TaxRank.KINGDOM.value
            return name, id, rank, rank_id
        else:
            parts2 = parts[1].split(';')
            id_part = parts2[0]
            id = id_part.split(':')[1].strip()
            rank_part = parts2[1]
            rank = rank_part.split(':')[1][:-1].strip()
            rank = rank.strip()
            rank_id = TaxRank[rank.upper()].value
            return name, id, rank, rank_id



    def parse_cf_reports(self, f):

        table_data = {}
        cls_lineage_data = {}
        els_lineage_data = None

        at_start = True
        in_ctl = 0
        in_etl = 0
        in_k = 0

        ctl_header = ""
        ctl_perc = ""
        etl_header = ""
        etl_perc = ""
        k_header = ""
        k_perc = ""

        i = 0
        for l in f['f']:
            line = l.strip()

            if line.startswith("Expected Taxon ID provided"):
                at_start = False
                continue

            if in_ctl == 1:
                ctl_header = line
                in_ctl = 2
            elif in_ctl == 2:
                ctl_perc = line
                in_ctl = 0

            if in_etl == 1:
                etl_header = line
                in_etl = 2
            elif in_etl == 2:
                etl_perc = line
                in_etl = 0

            if in_k == 1:
                k_header = line
                in_k = 2
            elif in_k == 2:
                k_perc = line
                in_k = 0

            if at_start:
                if line.startswith("Classified Taxon:"):
                    name, id, rank, rank_id = self.parse_taxon(line[18:].strip())
                    table_data["actual_taxon"] = name + " (" + id + ")"
                    table_data["actual_taxon_rank"] = rank
                    table_data["actual_taxon_rank_id"] = rank_id
                elif line.startswith("Classified Taxon's Lineage"):
                    in_ctl = 1
                elif line.startswith("Kingdom Perc"):
                    in_k = 1



            else:
                if line.startswith("Expected Taxon:"):
                    name, id, rank, rank_id = self.parse_taxon(line[16:].strip())
                    table_data["expected_taxon"] = name + " (" + id + ")"
                    table_data["expected_taxon_rank"] = rank
                    table_data["expected_taxon_rank_id"] = rank_id

                elif line.startswith("MRCA Taxon:"):
                    name, id, rank, rank_id = self.parse_taxon(line[11:].strip())
                    table_data["mrca_taxon"] = name + " (" + id + ")"
                    table_data["mrca_taxon_rank"] = rank
                    table_data["mrca_taxon_rank_id"] = rank_id
                elif line.startswith("Expected Taxon's Lineage"):
                    in_etl = 1

            i+=1

        k_head_parts = k_header.split('\t')
        k_perc_parts = k_perc.split('\t')
        kingdom_data = {}

        if len(k_head_parts) == 0 or len(k_perc_parts) == 0:
            raise ValueError("Error: Could not find kingdom percentages.  Invalid centrifuge summary file.")

        if len(k_perc_parts) != len(k_head_parts):
            raise ValueError("Error: Kingdom percentages has different number of header and percentage columns.")

        for i in range(len(k_head_parts)):
            h = k_head_parts[i]
            name, id, rank, rank_id = self.parse_taxon(h.strip(), kingdom=True)
            kingdom_data[name] = float(k_perc_parts[i])

        ctl_head_parts = ctl_header.split('\t')
        ctl_perc_parts = ctl_perc.split('\t')

        if len(ctl_head_parts) == 0 or len(ctl_perc_parts) == 0:
            raise ValueError("Error: Could not find classified taxon's lineage.  Invalid centrifuge summary file.")

        if len(ctl_perc_parts) != len(ctl_head_parts):
            raise ValueError("Error: Classified taxon's lineage has different number of header and percentage columns.")

        for i in range(len(ctl_head_parts)):
            h = ctl_head_parts[i]
            name, id, rank, rank_id = self.parse_taxon(h.strip())
            cls_lineage_data[TaxRank[rank.upper()].name.lower()] = float(ctl_perc_parts[i]) if i == 0 else (float(ctl_perc_parts[i]) - float(ctl_perc_parts[i-1]))

        els_lineage_data = {}

        if etl_header and etl_header != "":
            etl_head_parts = etl_header.split('\t')
            etl_perc_parts = etl_perc.split('\t')
            if len(etl_head_parts) == 0 or len(etl_perc_parts) == 0:
                raise ValueError("Error: Could not find expected taxon's lineage.  Invalid centrifuge summary file.")

            if len(etl_perc_parts) != len(etl_head_parts):
                raise ValueError(
                    "Error: Expected taxon's lineage has different number of header and percentage columns.")

            for i in range(len(etl_head_parts)):
                h = etl_head_parts[i]
                name, id, rank, rank_id = self.parse_taxon(h.strip())
                els_lineage_data[TaxRank[rank.upper()].name.lower()] = float(etl_perc_parts[i]) if i == 0 else (float(etl_perc_parts[i]) - float(etl_perc_parts[i - 1]))
        else:
            for i in range(len(ctl_head_parts)):
                h = ctl_head_parts[i]
                name, id, rank, rank_id = self.parse_taxon(h.strip())
                els_lineage_data[TaxRank[rank.upper()].name.lower()] = 100.0 if i == 8 else 0.0

        #print(kingdom_data)

        return table_data, cls_lineage_data, els_lineage_data, kingdom_data

        
