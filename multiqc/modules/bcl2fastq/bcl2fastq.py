from multiqc.modules.base_module import BaseMultiqcModule
import logging
import os
import json
from collections import OrderedDict
from multiqc import config
from multiqc.plots import bargraph, table

log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='bcl2fastq', anchor='bcl2fastq',
        href="https://support.illumina.com/sequencing/sequencing_software/bcl2fastq-conversion-software.html",
        info="can be used to both demultiplex data and convert BCL files to FASTQ file formats for downstream analysis.")

        # Gather data from all json files
        self.bcl2fastq_data = dict()
        for myfile in self.find_log_files('bcl2fastq'):
            self.parse_file_as_json(myfile)

        # Collect counts by lane and sample (+source_files)
        self.bcl2fastq_bylane = dict()
        self.bcl2fastq_undetermined = dict()
        self.bcl2fastq_bysample = dict()
        self.bcl2fastq_bysample_lane = dict()
        self.source_files = dict()
        self.split_data_by_lane_and_sample()

        # Filter to strip out ignored sample names
        self.bcl2fastq_bylane = self.ignore_samples(self.bcl2fastq_bylane)
        self.bcl2fastq_undetermined = self.ignore_samples(self.bcl2fastq_undetermined)
        self.bcl2fastq_bysample = self.ignore_samples(self.bcl2fastq_bysample)
        self.bcl2fastq_bysample_lane = self.ignore_samples(self.bcl2fastq_bysample_lane)

        # Return with Warning if no files are found
        if len(self.bcl2fastq_bylane) == 0 and len(self.bcl2fastq_bysample) == 0:
            raise UserWarning

        # Print source files
        for s in self.source_files.keys():
            self.add_data_source(s_name=s, source=",".join(list(set(self.source_files[s]))), module='bcl2fastq', section='bcl2fastq-bysample')

        # Add sample counts to general stats table
        self.add_general_stats()
        self.write_data_file({str(k): self.bcl2fastq_bylane[k] for k in self.bcl2fastq_bylane.keys()}, 'multiqc_bcl2fastq_bylane')
        self.write_data_file(self.bcl2fastq_bysample, 'multiqc_bcl2fastq_bysample')

        # Add section for summary stats per flow cell
        self.add_section (
            name = 'Lane Statistics',
            anchor = 'bcl2fastq-lanestats',
            description = 'Statistics about each lane for each flowcell',
            plot = self.lane_stats_table()
        )

        self.add_section(
            name = 'Top Undetermined Barcodes',
            anchor='bcl2fastq-undetermined',
            description='The undetermined barcodes from each lane with the highest counts',
            plot = self.undetermined_table()
        )

        # Add section for counts by lane
        cats = OrderedDict()
        cats["perfect"] = {'name': 'Perfect Index Reads'}
        cats["imperfect"] = {'name': 'Mismatched Index Reads'}
        cats["undetermined"] = {'name': 'Undetermined Reads'}
        self.add_section (
            name = 'Clusters by lane',
            anchor = 'bcl2fastq-bylane',
            description = 'Number of reads per lane (with number of perfect index reads).',
            helptext = """Perfect index reads are those that do not have a single mismatch.
                All samples of a lane are combined. Undetermined reads are treated as a third category.""",
            plot = bargraph.plot(
                self.get_bar_data_from_counts(self.bcl2fastq_bylane),
                cats,
                {
                    'id': 'bcl2fastq_lane_counts',
                    'title': 'bcl2fastq: Clusters by lane',
                    'hide_zero_cats': False
                }
            )
        )

        # Add section for counts by sample
        # get cats for per-lane tab
        lcats = set()
        for s_name in self.bcl2fastq_bysample_lane:
            lcats.update(self.bcl2fastq_bysample_lane[s_name].keys())
        lcats = sorted(list(lcats))
        self.add_section (
            name = 'Clusters by sample',
            anchor = 'bcl2fastq-bysample',
            description = 'Number of reads per sample.',
            helptext = """Perfect index reads are those that do not have a single mismatch.
                All samples are aggregated across lanes combinned. Undetermined reads are ignored.
                Undetermined reads are treated as a separate sample.""",
            plot = bargraph.plot(
                [
                    self.get_bar_data_from_counts(self.bcl2fastq_bysample),
                    self.bcl2fastq_bysample_lane
                ],
                [cats, lcats],
                {
                    'id': 'bcl2fastq_sample_counts',
                    'title': 'bcl2fastq: Clusters by sample',
                    'hide_zero_cats': False,
                    'ylab': 'Number of clusters',
                    'data_labels': ['Index mismatches', 'Counts per lane']
                }
            )
        )


    def parse_file_as_json(self, myfile):
        try:
            content = json.loads(myfile["f"])
        except ValueError:
            log.warn('Could not parse file as json: {}'.format(myfile["fn"]))
            return
        runId = content["RunId"]
        if not runId in self.bcl2fastq_data:
            self.bcl2fastq_data[runId] = dict()
        run_data = self.bcl2fastq_data[runId]
        for conversionResult in content.get("ConversionResults", []):
            lane = 'L{}'.format(conversionResult["LaneNumber"])
            if lane in run_data:
                log.debug("Duplicate runId/lane combination found! Overwriting: {}".format(self.prepend_runid(runId, lane)))
            run_data[lane] = {
                "total": 0,
                "total_yield": 0,
                "perfectIndex": 0,
                "samples": dict(),
                "yieldQ30": 0,
                "qscore_sum": 0,
                "unknown_barcodes": dict()
            }
            for demuxResult in conversionResult.get("DemuxResults", []):
                sample = demuxResult["SampleName"]
                if sample in run_data[lane]["samples"]:
                    log.debug("Duplicate runId/lane/sample combination found! Overwriting: {}, {}".format(self.prepend_runid(runId, lane),sample))
                run_data[lane]["samples"][sample] = {
                    "total": 0,
                    "total_yield": 0,
                    "perfectIndex": 0,
                    "barcode": "",
                    "filename": os.path.join(myfile['root'],myfile["fn"]),
                    "yieldQ30": 0,
                    "qscore_sum": 0
                }
                run_data[lane]["total"] += demuxResult["NumberReads"]
                run_data[lane]["total_yield"] += demuxResult["Yield"]
                run_data[lane]["samples"][sample]["total"] += demuxResult["NumberReads"]
                run_data[lane]["samples"][sample]["total_yield"] += demuxResult["Yield"]
                run_data[lane]["samples"][sample]["barcode"] = ""
                if "IndexMetrics" in demuxResult:
                    for indexMetric in demuxResult["IndexMetrics"]:
                        run_data[lane]["samples"][sample]["barcode"] = indexMetric["IndexSequence"]
                        run_data[lane]["perfectIndex"] += indexMetric["MismatchCounts"]["0"]
                        run_data[lane]["samples"][sample]["perfectIndex"] += indexMetric["MismatchCounts"]["0"]
                for readMetric in demuxResult.get("ReadMetrics", []):
                    run_data[lane]["yieldQ30"] += readMetric["YieldQ30"]
                    run_data[lane]["qscore_sum"] += readMetric["QualityScoreSum"]
                    run_data[lane]["samples"][sample]["yieldQ30"] += readMetric["YieldQ30"]
                    run_data[lane]["samples"][sample]["qscore_sum"] += readMetric["QualityScoreSum"]
            undeterminedYieldQ30 = 0
            undeterminedQscoreSum = 0
            if "Undetermined" in conversionResult:
                for readMetric in conversionResult["Undetermined"]["ReadMetrics"]:
                    undeterminedYieldQ30 += readMetric["YieldQ30"]
                    undeterminedQscoreSum += readMetric["QualityScoreSum"]
                run_data[lane]["samples"]["undetermined"] = {
                    "total": conversionResult["Undetermined"]["NumberReads"],
                    "total_yield": conversionResult["Undetermined"]["Yield"],
                    "perfectIndex": 0,
                    "yieldQ30": undeterminedYieldQ30,
                    "qscore_sum": undeterminedQscoreSum,
                    "barcode": ""
                }

        for unknownBarcodes in content.get("UnknownBarcodes", []):
            lane = 'L{}'.format(unknownBarcodes["Lane"])
            lane_unknown_barcodes = unknownBarcodes.get("Barcodes", [])
            nb_to_count = min(len(lane_unknown_barcodes), 5)
            run_data[lane]["unknown_barcodes"] = {
                "1": {"barcode":"N/A", "count": 0},
                "2": {"barcode":"N/A", "count": 0},
                "3": {"barcode":"N/A", "count": 0},
                "4": {"barcode":"N/A", "count": 0},
                "5": {"barcode":"N/A", "count": 0}
            }
            for i, barcode in enumerate(lane_unknown_barcodes):
                if i == nb_to_count:
                    break
                run_data[lane]["unknown_barcodes"][str(i+1)]["barcode"] = barcode
                run_data[lane]["unknown_barcodes"][str(i+1)]["count"] = int(lane_unknown_barcodes[barcode])


        # Calculate Percents and averages
        for lane in run_data:
            run_data[lane]["percent_Q30"] = safe_div(float(run_data[lane]["yieldQ30"]), float(run_data[lane]["total_yield"])) * 100.0
            run_data[lane]["percent_perfectIndex"] = safe_div(float(run_data[lane]["perfectIndex"]), float(run_data[lane]["total"])) * 100.0
            run_data[lane]["mean_qscore"] = safe_div(float(run_data[lane]["qscore_sum"]), float(run_data[lane]["total_yield"]))
            for sample, d in run_data[lane]["samples"].items():
                tyd = float(d["total_yield"])
                td = float(d["total"])
                run_data[lane]["samples"][sample]["percent_Q30"] = safe_div(float(d["yieldQ30"]), float(d["total_yield"])) * 100.0
                run_data[lane]["samples"][sample]["percent_perfectIndex"] = safe_div(float(d["perfectIndex"]), float(d["total"])) * 100.0
                run_data[lane]["samples"][sample]["mean_qscore"] = safe_div(float(d["qscore_sum"]), float(d["total_yield"]))

    def split_data_by_lane_and_sample(self):
        for runId in self.bcl2fastq_data.keys():
            for lane in self.bcl2fastq_data[runId].keys():
                uniqLaneName = self.prepend_runid(runId, lane)
                self.bcl2fastq_bylane[uniqLaneName] = {
                    "total": self.bcl2fastq_data[runId][lane]["total"],
                    "total_yield": self.bcl2fastq_data[runId][lane]["total_yield"],
                    "perfectIndex": self.bcl2fastq_data[runId][lane]["perfectIndex"],
                    "undetermined": self.bcl2fastq_data[runId][lane]["samples"].get("undetermined", {}).get("total", "NA"),
                    "yieldQ30": self.bcl2fastq_data[runId][lane]["yieldQ30"],
                    "qscore_sum": self.bcl2fastq_data[runId][lane]["qscore_sum"],
                    "percent_Q30": self.bcl2fastq_data[runId][lane]["percent_Q30"],
                    "percent_perfectIndex": self.bcl2fastq_data[runId][lane]["percent_perfectIndex"],
                    "mean_qscore": self.bcl2fastq_data[runId][lane]["mean_qscore"]
                }
                self.bcl2fastq_undetermined[uniqLaneName] = {
                    "first_barcode": self.bcl2fastq_data[runId][lane]["unknown_barcodes"]['1']['barcode'],
                    "first_count": self.bcl2fastq_data[runId][lane]["unknown_barcodes"]['1']['count'],
                    "second_barcode": self.bcl2fastq_data[runId][lane]["unknown_barcodes"]['2']['barcode'],
                    "second_count": self.bcl2fastq_data[runId][lane]["unknown_barcodes"]['2']['count'],
                    "third_barcode": self.bcl2fastq_data[runId][lane]["unknown_barcodes"]['3']['barcode'],
                    "third_count": self.bcl2fastq_data[runId][lane]["unknown_barcodes"]['3']['count'],
                    "fourth_barcode": self.bcl2fastq_data[runId][lane]["unknown_barcodes"]['4']['barcode'],
                    "fourth_count": self.bcl2fastq_data[runId][lane]["unknown_barcodes"]['4']['count'],
                    "fifth_barcode": self.bcl2fastq_data[runId][lane]["unknown_barcodes"]['5']['barcode'],
                    "fifth_count": self.bcl2fastq_data[runId][lane]["unknown_barcodes"]['5']['count']
                }
                for sample in self.bcl2fastq_data[runId][lane]["samples"].keys():
                    if not sample in self.bcl2fastq_bysample:
                        self.bcl2fastq_bysample[sample] = {
                            "total": 0,
                            "total_yield": 0,
                            "perfectIndex": 0,
                            "yieldQ30": 0,
                            "qscore_sum": 0,
                            "barcode": ""
                        }
                    if not sample in self.bcl2fastq_bysample_lane:
                        self.bcl2fastq_bysample_lane[sample] = dict()
                    self.bcl2fastq_bysample_lane[sample][lane] = self.bcl2fastq_data[runId][lane]["samples"][sample]["total"]
                    self.bcl2fastq_bysample[sample]["total"] += self.bcl2fastq_data[runId][lane]["samples"][sample]["total"]
                    self.bcl2fastq_bysample[sample]["total_yield"] += self.bcl2fastq_data[runId][lane]["samples"][sample]["total_yield"]
                    self.bcl2fastq_bysample[sample]["perfectIndex"] += self.bcl2fastq_data[runId][lane]["samples"][sample]["perfectIndex"]
                    self.bcl2fastq_bysample[sample]["yieldQ30"] += self.bcl2fastq_data[runId][lane]["samples"][sample]["yieldQ30"]
                    self.bcl2fastq_bysample[sample]["qscore_sum"] += self.bcl2fastq_data[runId][lane]["samples"][sample]["qscore_sum"]
                    self.bcl2fastq_bysample[sample]["barcode"] = self.bcl2fastq_data[runId][lane]["samples"][sample]["barcode"]
                    self.bcl2fastq_bysample[sample]["percent_Q30"] = safe_div(float(self.bcl2fastq_bysample[sample]["yieldQ30"]) / float(self.bcl2fastq_bysample[sample]["total_yield"])) * 100.0
                    self.bcl2fastq_bysample[sample]["percent_perfectIndex"] = safe_div(float(self.bcl2fastq_bysample[sample]["perfectIndex"]) / float(self.bcl2fastq_bysample[sample]["total"])) * 100.0
                    self.bcl2fastq_bysample[sample]["mean_qscore"] = safe_div(float(self.bcl2fastq_bysample[sample]["qscore_sum"]) / float(self.bcl2fastq_bysample[sample]["total_yield"]))
                    if sample != "undetermined":
                        if not sample in self.source_files:
                            self.source_files[sample] = []
                        self.source_files[sample].append(self.bcl2fastq_data[runId][lane]["samples"][sample]["filename"])

    def add_general_stats(self):
        data = {
            key: {
                "yieldQ30": self.bcl2fastq_bysample[key]["yieldQ30"],
                "total": self.bcl2fastq_bysample[key]["total"],
                "perfectPercent": '{0:.1f}'.format(
                    float( 100.0 * safe_div(self.bcl2fastq_bysample[key]["perfectIndex"], self.bcl2fastq_bysample[key]["total"] ))),
                "barcode": self.bcl2fastq_bysample[key]["barcode"]
            } for key in self.bcl2fastq_bysample.keys()
        }
        headers = OrderedDict()
        headers['total'] = {
            'title': 'Clusters',
            'description': 'Total number of reads for this sample as determined by bcl2fastq demultiplexing',
            'scale': 'Blues',
            'format': '{:,.0f}'
        }
        headers['yieldQ30'] = {
            'title': '{} Yield &ge; Q30'.format(config.base_count_prefix),
            'description': 'Number of bases with a Phred score of 30 or higher ({})'.format(config.base_count_desc),
            'scale': 'Greens',
            'shared_key': 'base_count'
        }
        headers['perfectPercent'] = {
            'title': '% Perfect Index',
            'description': 'Percent of reads with perfect index (0 mismatches)',
            'max': 100,
            'min': 0,
            'scale': 'RdYlGn',
            'suffix': '%'
        }
        headers['barcode'] = {
            'title': 'Barcode',
            'description': 'The combined index sequence (i7 + i5) if present'
        }
        self.general_stats_addcols(data, headers)

    def lane_stats_table(self):
        """ Return a table with overview stats for each bcl2fastq lane for a single flow cell """
        headers = OrderedDict()
        headers['total_yield'] = {
            'title': '{} Total Yield'.format(config.base_count_prefix),
            'description': 'Number of bases ({})'.format(config.base_count_desc),
            'scale': 'Greens',
            'shared_key': 'base_count'
        }
        headers['total'] = {
            'title': 'Total Clusters',
            'description': 'Total number of clusters for this lane',
            'scale': 'Blues',
            'format': '{:,.0f}'
        }
        headers['percent_Q30'] = {
            'title': '% bases &ge; Q30',
            'description': 'Percentage of bases with greater than or equal to Q30 quality score',
            'suffix': '%',
            'max': 100,
            'min': 0,
            'scale': 'RdYlGn'
        }
        headers['mean_qscore'] = {
            'title': 'Mean Quality',
            'description': 'Average phred qualty score',
            'min': 0,
            'scale': 'Spectral'
        }
        headers['percent_perfectIndex'] = {
            'title': '% Perfect Index',
            'description': 'Percent of reads with perfect index (0 mismatches)',
            'max': 100,
            'min': 0,
            'scale': 'RdYlGn',
            'suffix': '%'
        }
        table_config = {
            'namespace': 'bcl2fastq',
            'id': 'bcl2fastq-lane-stats-table',
            'table_title': 'bcl2fastq Lane Statistics',
            'col1_header': 'Run ID - Lane',
            'no_beeswarm': True
        }
        return table.plot(self.bcl2fastq_bylane, headers, table_config)

    def undetermined_table(self):
        """ Return a table with overview stats for each bcl2fastq lane for a single flow cell """
        headers = OrderedDict()
        headers['first_barcode'] = {
            'title': '1st Barcode',
            'description': 'Barcode with highest count',
            'format': '{:50}'
        }
        headers['first_count'] = {
            'title': "1st count",
            'description': "Number of reads associated with top undetermined barcode",
            'scale': 'RdYlGn',
            'format': '{:,.0f}',
            'min': 0,
            'custom_style': 'width=5%'
        }
        headers['second_barcode'] = {
            'title': '2nd Barcode',
            'description': 'Barcode with second highest count',
            'format': '{:50}'
        }
        headers['second_count'] = {
            'title': "2nd count",
            'description': "Number of reads associated with second top undetermined barcode",
            'scale': 'RdYlGn',
            'format': '{:,.0f}',
            'min': 0,
            'custom_style': 'width=5%'
        }
        headers['third_barcode'] = {
            'title': '3rd Barcode',
            'description': 'Barcode with third highest count',
            'format': '{:50}'
        }
        headers['third_count'] = {
            'title': "3rd count",
            'description': "Number of reads associated with third top undetermined barcode",
            'scale': 'RdYlGn',
            'format': '{:,.0f}',
            'min': 0,
            'custom_style': 'width=5%'
        }
        headers['fourth_barcode'] = {
            'title': '4th Barcode',
            'description': 'Barcode with fourth highest count',
            'format': '{:50}',
            'hidden': True
        }
        headers['fourth_count'] = {
            'title': "4th count",
            'description': "Number of reads associated with fourth top undetermined barcode",
            'scale': 'RdYlGn',
            'format': '{:,.0f}',
            'min': 0,
            'custom_style': 'width=5%',
            'hidden': True
        }
        headers['fifth_barcode'] = {
            'title': '5th Barcode',
            'description': 'Barcode with fifth highest count',
            'format': '{:50}',
            'hidden': True
        }
        headers['fifth_count'] = {
            'title': "5th count",
            'description': "Number of reads associated with fifth top undetermined barcode",
            'scale': 'RdYlGn',
            'format': '{:,.0f}',
            'min': 0,
            'custom_style': 'width=5%',
            'hidden': True
        }

        table_config = {
            'namespace': 'bcl2fastq',
            'id': 'bcl2fastq-undetermined-table',
            'table_title': 'bcl2fastq Undetermined Barcodes',
            'col1_header': 'Run ID - Lane',
            'no_beeswarm': True
        }
        return table.plot(self.bcl2fastq_undetermined, headers, table_config)

    def prepend_runid(self, runId, rest):
        return str(runId)+" - "+str(rest)

    def get_bar_data_from_counts(self, counts):
        bar_data = {}
        for key, value in counts.items():
            bar_data[key] = {
                "perfect": value["perfectIndex"],
                "imperfect": value["total"] - value["perfectIndex"],
            }
            if "undetermined" in value:
                bar_data[key]["undetermined"] = value["undetermined"]
        return bar_data

def safe_div(x,y):
    if y == 0:
        return 0
    return x / y