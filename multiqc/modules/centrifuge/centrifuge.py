from multiqc.modules.base_module import BaseMultiqcModule
import logging
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='Centrifuge', anchor='centrifuge',
        href="https://ccb.jhu.edu/software/centrifuge/index.shtml",
        info="Summarises centrifuge reports from multiple samples")
        
        self.cf_data = dict()
        for f in self.find_log_files('centrifuge', filehandles=True):
            self.cf_data[f['s_name']] = self.parse_cf_reports(f)
            self.add_data_source(f)
        
        self.cf_data = self.ignore_samples(self.cf_data)

        if len(self.cf_data) == 0:
            log.debug("Could not find any data in {}".format(config.analysis_dir))
            raise UserWarning

        log.info("Found {} reports".format(len(self.cf_data)))

        self.write_data_file(self.cf_data, 'multiqc_centrifuge')

    def parse_cf_reports(self, f):
        i = 0
        for l in f['f']:
            i += 1
        return i

        
