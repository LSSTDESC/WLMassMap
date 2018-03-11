from pipette import PipelineStage
from pipette.types import HDFFile
from .dtypes import ShearCatFile
from desc.wlmassmap.mocks import extract_footprint, mock_observation

class extractFootprintPipe(PipelineStage):
    name = 'extractFootprintPipe'
    inputs = []
    outputs = [('truth_catalog', HDFFile)]

    config_options = {'gcr_catalog':'proto-dc2_v2.0',
                      'footprint':{
                        'type':'patch',
                        'ra_range': [0,5],
                        'dec_range': [0,5]}}

    def run(self):
        config = self.read_config()
        config['output_filename'] = self.get_output('truth_catalog')
        extract_footprint(config)


class mockShearMeasurementPipe(PipelineStage):
    name = 'mockShearMeasurementPipe'
    inputs = [('truth_catalog', HDFFile)]
    outputs = [('shear_catalog', ShearCatFile)]
    config_options = {'reduced_shear':True,
                      'shape_noise':{
                            'type':'Gaussian',
                            'sigma': 0.13},
                      'format':{
                            'type':'metacal',
                            'R':[[1,0],[0,1]],
                            'delta_gamma': 0.01}}

    def run(self):
        """
        Runs the code
        """
        config = self.read_config()
        config['input_filename'] = self.get_input('truth_catalog')
        config['output_filename'] = self.get_output('shear_catalog')
        mock_observation(config)
