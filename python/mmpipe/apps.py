from ceci import PipelineStage
from descformats import HDFFile, YamlFile, FitsFile
from descformats.tx import MetacalCatalog, TomographyCatalog
from desc.wlmassmap.mocks.extract_footprint import extract_footprint
from desc.wlmassmap.mocks.mock_observation import mock_observation
from desc.wlmassmap.shear_map import shear_map
from desc.wlmassmap.convergence_map import convergence_map
from collections import defaultdict

class extractFootprintPipe(PipelineStage):
    name = 'extractFootprintPipe'
    inputs = []
    outputs = [('truth_catalog', HDFFile)]
    config_options = {'catalog':'protoDC2',
                      'ra_range':[0.,5.],
                      'dec_range':[0.,5.]}

    def run(self):
        config = self.read_config(defaultdict(lambda :None))
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
        config = self.read_config(defaultdict(lambda :None))
        config['input_filename'] = self.get_input('truth_catalog')
        config['output_filename'] = self.get_output('shear_catalog')
        config['shape_noise'] = {'type':'Gaussian', 'sigma':config['sigma_noise']}
        config['format'] = {'type':'metacal', 'R':[[1,0],[0,1]],
                            'delta_gamma':config['delta_gamma']}
        mock_observation(config)

class shearMapPipe(PipelineStage):
    name = 'shearMapPipe'
    inputs = [('shear_catalog', MetacalCatalog)]
    outputs = [('shear_map', FitsFile)]
    config_options = {'center_ra':float,
                      'center_dec':float,
                      'pixel_size':1.,
                      'nx':300,
                      'ny':300}

    def run(self):
        config = self.read_config(defaultdict(lambda :None))
        config['input_filename'] = self.get_input('shear_catalog')
        config['output_filename'] = self.get_output('shear_map')
        config['projection'] = {'type':'gnomonic',
                                'center_ra':config['center_ra'],
                                'center_dec':config['center_dec'],
                                'pixel_size':config['pixel_size'],
                                'nx':config['nx'],
                                'ny':config['ny']}
        shear_map(config)


class convergenceMapPipe(PipelineStage):
    name = 'convergenceMapPipe'
    inputs = [('shear_map', FitsFile)]
    outputs = [('converenge_map', FitsFile)]
    config_options = {'smoothing':1.,
                      'zero_padding': 128}

    def run(self):
        config = self.read_config(defaultdict(lambda :None))
        config['input_filename'] = self.get_input('shear_map')
        config['output_filename'] = self.get_output('converenge_map')
        config['algorithm'] = {'name':'flat_ks',
                                'smoothing':config['smoothing'],
                                'zero_padding':config['zero_padding']}
        convergence_map(config)

if __name__ == '__main__':
    cls = PipelineStage.main()
