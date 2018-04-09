from ceci import PipelineStage
from descformats import MetacalCatalog, HDFFile, YamlFile
from desc.wlmassmap.mocks.extract_footprint import extract_footprint
from desc.wlmassmap.mocks.mock_observation import mock_observation

class extractFootprintPipe(PipelineStage):
    name = 'extractFootprintPipe'
    inputs = []
    outputs = [('truth_catalog', HDFFile)]
    config_options = {'catalog':'protoDC2',
                      'ra_range':[0.,5.],
                      'dec_range':[0.,5.]}

    def run(self):
        config = self.read_config()
        config['output_filename'] = self.get_output('truth_catalog')
        config['footprint'] = {'type':'patch', 'ra_range':config['ra_range'],
                               'dec_range':config['dec_range']}
        extract_footprint(config)

class mockShearMeasurementPipe(PipelineStage):
    name = 'mockShearMeasurementPipe'
    inputs = [('truth_catalog', HDFFile)]
    outputs = [('shear_catalog', MetacalCatalog)]
    config_options = {'reduced_shear':True,
                      'sigma_noise': float,
                      'delta_gamma': 0.01}

    def run(self):
        """
        Runs the code
        """
        config = self.read_config()
        config['input_filename'] = self.get_input('truth_catalog')
        config['output_filename'] = self.get_output('shear_catalog')
        config['shape_noise'] = {'type':'Gaussian', 'sigma':config['sigma_noise']}
        config['format'] = {'type':'metacal', 'R':[[1,0],[1,0]],
                            'delta_gamma':config['delta_gamma']}
        mock_observation(config)

if __name__ == '__main__':
    cls = PipelineStage.main()
