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
	SUBGENUS = 7.5
	SPECIES = 8
	SUBSPECIES = 8.5
	NO_RANK = 9

	@classmethod
	def get(cls, item, alt=cls.ROOT):
		return cls.__members__.get(item, alt).value

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
        info=""" is a very rapid and memory-efficient system for the classification of DNA sequences.  The system uses a 
    novel indexing scheme based on the Burrows-Wheeler transform (BWT) and the Ferragina-Manzini (FM) index, optimized 
    specifically for the metagenomic classification problem.  The version of centrifuge used here was modified by EI to accomodate
    efficient processing of multiple samples against the NCBI NT Database.  Centrifuge reports were parsed and interpretted 
    using unique read counts per taxon to produce the summary statistics displayed below.""")
        
        self.table_data = {}
        self.cls_lineage_data = {}
        self.els_lineage_data = {}
        self.kingdom_data = {}
        self.rank_data={}
        self.top5_species_data = {}
        self.top5_class_data = {}

        for c_file in self.find_log_files('centrifuge'):
            s_name = self.clean_s_name(c_file['s_name'][:-11], c_file['root'])

            content = c_file['f'].splitlines()

            self.table_data[s_name], self.cls_lineage_data[s_name], self.els_lineage_data[s_name], self.kingdom_data[s_name], self.rank_data[s_name], self.top5_species_data[s_name], self.top5_class_data[s_name] = self.parse_cf_reports(content)
            self.add_data_source(c_file, s_name)

        self.table_data = self.ignore_samples(self.table_data)
        self.cls_lineage_data = self.ignore_samples(self.cls_lineage_data)
        self.els_lineage_data = self.ignore_samples(self.els_lineage_data)
        self.kingdom_data = self.ignore_samples(self.kingdom_data)
        self.rank_data = self.ignore_samples(self.rank_data)
        self.top5_species_data = self.ignore_samples(self.top5_species_data)
        self.top5_class_data = self.ignore_samples(self.top5_class_data)

        if len(self.table_data) == 0:
            log.debug("Could not find any data in {}".format(config.analysis_dir))
            raise UserWarning

        log.info("Found {} reports".format(len(self.table_data)))

        headers = OrderedDict()
        headers["total_hits"] = {
            'title': 'Total Hits',
            'description': "Total hits to DB",
            'format': '{:,.0f}',
            'custom_style': 'width=10%'
        }
        headers["actual_taxon"] = {
            'title': 'Classified Taxon',
            'description': "Classified taxon"
        }
        headers["actual_taxon_rank"] = {
            'title': 'Rank',
            'description': "Classified taxon's rank",
            'custom_style': 'width=6%'
        }
        headers["actual_taxon_rank_id"] = {
            'title': 'Level',
            'description': "Classified taxon's rank level (0=Root -> 8=Species)",
            'scale': 'RdYlGn',
            'min': 0,
            'max': 8,
            'format': '{:,.0f}',
            'custom_style': 'width=3%'
        }
        headers["expected_taxon"] = {
            'title': 'Expected Taxon',
            'description': "Expected taxon"
        }
        headers["expected_taxon_rank"] = {
            'title': 'Rank',
            'description': "Expected taxon's rank",
            'custom_style': 'width=6%'
        }
        headers["expected_taxon_rank_id"] = {
            'title': 'Level',
            'description': "Expected taxon's rank level (0=Root -> 8=Species)",
            'scale': 'RdYlGn',
            'min': 0,
            'max': 8,
            'format': '{:,.0f}',
            'custom_style': 'width=3%'
        }
        headers["mrca_taxon"] = {
            'title': 'MRCA Taxon',
            'description': "Most recent common ancestor taxon"
        }
        headers["mrca_taxon_rank"] = {
            'title': 'Rank',
            'description': "Most recent common ancestor taxon's rank",
            'custom_style': 'width=6%'
        }
        headers["mrca_taxon_rank_id"] = {
            'title': 'Level',
            'description': "Most recent common ancestor taxon's rank level (0=Root -> 8=Species)",
            'scale': 'RdYlGn',
            'min': 0,
            'max': 8,
            'format': '{:,.0f}',
            'custom_style': 'width=3%'
        }
        tc_config = {
            'namespace': 'Centrifuge',
            'scale': 'RdYlGn',
            'save_file': True,
            'raw_data_fn': 'mqc_centrifuge_taxclas'
        }

        self.add_section(
            name='Taxonomic Classification',
            anchor='centrifuge-first',
            description='Table showing classified taxon for each sample.  The classfied taxon is typically the first category '
                        '(starting from species level) to contain > 50% of total hits.  If sample has an expected taxon provided '
                        'in samplesheet then this is also displayed, along with the most recent common ancestor of the classified and expected taxa.',
            helptext="This table is intended to give a quick view of what makes up the majority of the sample, and if an expected taxon id was specified,"
                     "whether that matches up with what was detected.",
            plot=table.plot(self.table_data, headers, tc_config)
        )


        t5_headers = OrderedDict()
        t5_headers['first_name'] = {
            'title': "1st Name",
            'description': "Taxonomic name and id at 1st position",
            'format': '{:50}'
        }
        t5_headers['first_count'] = {
            'title': "1st %",
            'description': "Percentage of centrifuge hits allocated to 1st position",
            'scale': 'RdYlGn',
            'min': 0,
            'max': 100.0,
            'format': '{:,.1f}',
            'custom_style': 'width=5%'
        }

        t5_headers['second_name'] = {
            'title': "2nd Name",
            'description': "Taxonomic name and id at 2nd position",
            'format': '{:50}'
        }
        t5_headers['second_count'] = {
            'title': "2nd %",
            'description': "Percentage of centrifuge hits allocated to 2nd position",
            'scale': 'RdYlGn',
            'min': 0,
            'max': 100.0,
            'format': '{:,.1f}',
            'custom_style': 'width=5%'
        }
        t5_headers['third_name'] = {
            'title': "3rd Name",
            'description': "Taxonomic name and id at 3rd position",
            'format': '{:50}'
        }
        t5_headers['third_count'] = {
            'title': "3rd %",
            'description': "Percentage of centrifuge hits allocated to 3rd position",
            'scale': 'RdYlGn',
            'min': 0,
            'max': 100.0,
            'format': '{:,.1f}',
            'custom_style': 'width=5%'
        }
        t5_headers['fourth_name'] = {
            'title': "4th Name",
            'description': "Taxonomic name and id at 4th position",
            'format': '{:50}',
            'hidden': True
        }
        t5_headers['fourth_count'] = {
            'title': "4th %",
            'description': "Percentage of centrifuge hits allocated to 4th position",
            'scale': 'RdYlGn',
            'min': 0,
            'max': 100.0,
            'format': '{:,.1f}',
            'custom_style': 'width=5%',
            'hidden': True
        }
        t5_headers['fifth_name'] = {
            'title': "5th Name",
            'description': "Taxonomic name and id at 5th position",
            'format': '{:50}',
            'hidden': True
        }
        t5_headers['fifth_count'] = {
            'title': "5th %",
            'description': "Percentage of centrifuge hits allocated to 5th position",
            'scale': 'RdYlGn',
            'min': 0.0,
            'max': 100.0,
            'format': '{:,.1f}',
            'custom_style': 'width=5%',
            'hidden': True
        }
        t5_species_config = {
            'namespace': 'Centrifuge',
            'scale': 'RdYlGn',
            'save_file': True,
            'raw_data_fn': 'mqc_centrifuge_top5_species'
        }
        self.add_section(
            name='Top Species',
            anchor='centrifuge-second',
            description='This table shows the most abundant species for each sample. Top 3 shown by default although up to top 5 available through configure columns button.',
            helptext="The intention of this table is to show a low-level view of what is present in each sample.",
            plot=table.plot(self.top5_species_data, t5_headers, t5_species_config)
        )
        t5_classes_config = {
            'namespace': 'Centrifuge',
            'scale': 'RdYlGn',
            'save_file': True,
            'raw_data_fn': 'mqc_centrifuge_top5_classes'
        }
        self.add_section(
            name='Top Classes',
            anchor='centrifuge-third',
            description='This table shows the most abundant taxonomic classes for each sample. Top 3 shown by default although up to top 5 available through configure columns button.',
            helptext='Occasionally, groupings from lower ranks, such as "Order" or "Family" may be promoted to class level '
                     'and registered here if most abundant species has no defined class in it\'s taxonomic lineage.'
                     'The intention of this table is to give a mid-level taxonomic view of what is present in the sample',
            plot=table.plot(self.top5_class_data, t5_headers, t5_classes_config)
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
            'id': 'centrifuge_kingdom',
            'cpswitch': True,
            'cpswitch_c_active': False,
            'title': 'Hits in each kingdom',
            'ymin': 0

        }

        self.add_section(
            name='Kingdom Classification',
            anchor='centrifuge-fourth',
            description='Bar plot showing hits falling into each Kingdom (or Super Kingdom).',
            helptext="Strictly speaking this plot does not necessarily contain all Kingdom-level taxonomic groupings, "
                     "also it contains some groupings at super kingdom level.  The definitely list shown here is: 'Bacteria', 'Archaea', 'Viruses', 'Fungi', 'Viridiplantae' and 'Metazoa'."
                     "Anything else is not shown.  The intention is to give a high-level view of what is in the sample.",
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
        rank_config = {
            'id': 'centrifuge_ranks',
            'cpswitch': True,
            'cpswitch_c_active': False,
            'title': 'Hits at each rank in taxons lineage',
            'ymin': 0
        }

        self.add_section(
            name='Rank Hits',
            anchor='centrifuge-fifth',
            description='Bar plot showing the rank at which centrifuge hits were classified.',
            helptext="Most of the time we should expect most hits to occur at the species level, hence the plot should look"
                     "mostly green.  However, on occasion centrifuge may create hits at higher taxonomic level if there is"
                     "a high level of ambiguity in the alignments.  If that's the case here then it warrents further investigation.",
            plot=bargraph.plot(self.rank_data, cats, rank_config)
        )

        cls_lineage_config = {
            'id': 'centrifuge_classified_lineage',
            'cpswitch': True,
            'cpswitch_c_active': False,
            'title': 'Hits at each rank in taxons lineage',
            'ymin': 0
        }

        self.add_section(
            name='Classified Taxon Lineage',
            anchor='centrifuge-sixth',
            description='Bar plot showing the classified taxons lineage (largest taxa at each taxonomic rank incorperating the classified taxon).',
            helptext="This plot is showing the proportions of hits that occur along a taxonomic lineage that includes the"
                     "classified taxon id.  If this id is above species level then we take the largest taxa at each level"
                     "down to species.  If the sequenceing was of a particular organism that is well captured in the database (typically NCBI NT, unless otherwise directed),"
                     "then we would expect the classified taxon id to be at species level and the species (green) proportion of this"
                     "plot to by greater than 50%.  If this was a metagenomic sample, or a non-model organism was sequenced, or"
                     "if the sample has heavy contamination, then we might expect the middle to upper ranks to contain larger"
                     "proportions of the bar graph.",
            plot=bargraph.plot(self.cls_lineage_data, cats, cls_lineage_config)
        )

        exp_lineage_config = {
            'id': 'centrifuge_expected_lineage',
            'cpswitch': True,
            'cpswitch_c_active': False,
            'title': 'Hits at each rank in taxons lineage',
            'ymin': 0
        }

        self.add_section(
            name='Expected Taxon Lineage',
            anchor='centrifuge-seventh',
            description='Bar plot showing the expected taxons lineage (largest taxa at each taxonomic rank incoperating the expected taxon).',
            helptext="See help for classified taxon's lineage for more details but to summarise this shows the lineage containing"
                     "the expected taxa if present.  If it wasn't provided then this plot will be all gray.  If the expected"
                     "taxa was a close match to the classfied taxa this should look very similar to the bar chart above.",
            plot=bargraph.plot(self.els_lineage_data, cats, exp_lineage_config)
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
            rank_id = TaxRank.get(rank.upper().replace(" ", "_")).value # TaxRank[rank.upper()].value if rank != '?' else 0
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
        t5_species_data = self.extract_t5(content, "Top 5 species", TaxRank.SPECIES)
        t5_class_data = self.extract_t5(content, "Top 5 class", TaxRank.CLASS)

        els_lineage_data = {}
        if exp_present:
            els_lineage_data = self.extract_lineage(content, "Expected Taxon's Lineage")
        else:
            for i, r in enumerate(TaxRank):
                els_lineage_data[r.name.lower()] = 100.0 if i == 0 else 0.0

        return table_data, cls_lineage_data, els_lineage_data, kingdom_data, rank_data, t5_species_data, t5_class_data



    def extract_t5(self, content, header, in_rank):

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

        if len(desc_parts) != len(perc_parts) or len(desc_parts) != len(count_parts):
            raise ValueError("Error: Lineage with header: " + header + "; has different number of header and value columns.")

        for i in range(len(desc_parts)):
            h = desc_parts[i].strip()
            #print(h)
            #name, id, rank, rank_id = self.parse_taxon(h, rank=in_rank)
            t5_data[mapping[str(i)] + "_name"] = h
            t5_data[mapping[str(i)] + "_count"] = float(perc_parts[i])

        if (len(desc_parts) < 5):
            for i in range(len(desc_parts),5):
                t5_data[mapping[str(i)] + "_name"] = "N/A"
                t5_data[mapping[str(i)] + "_count"] = 0.0

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
            #lineage_data[TaxRank[rank.upper()].name.lower()] = int(count_parts[i]) if i == 0 else (int(count_parts[i]) - int(count_parts[i - 1]))
            lineage_data[TaxRank.get(rank.upper()).name.lower()] = int(count_parts[i]) if i == 0 else (int(count_parts[i]) - int(count_parts[i - 1]))

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
