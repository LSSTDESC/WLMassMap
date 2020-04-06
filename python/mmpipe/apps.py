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

class TXConvergenceMaps(PipelineStage):
    name = 'TXConvergenceMaps'
    inputs = [
        ('diagnostic_maps', DiaognisticMaps)
    ]
    outputs = [
        ('mass_maps', DiagnosticMaps)
    ]
    config_options = {
        'smoothing': 5.,
        'algorithm': 'spherical_ks',
        'flip_g1': False,
        'flip_g2': False,
        'zero_padding': 128,
        'lmax': 4096
    }

    def run(self):
        from desc.wlmassmap.kaiser_squires import flat_KS_map, healpix_KS_map

        config = self.config
        shear_maps, weights = self.load_maps()

        if config['algorithm'] == 'flat_ks':
            kappa_e, kappa_b = flat_KS_map(gmap)

        elif config['algorithm'] == 'healpix_ks':
            kappa_e, kappa_b = healpix_KS_map(gmap, lmax=config['lmax'])
        else:
            raise NotImplementedError


    def load_maps():
        import healpy

        # Load the various input maps and their metadata
        map_file = self.open_input('diagnostic_maps', wrapper=True)
        pix_info = map_file.read_map_info('mask')
        area = map_file.file['maps'].attrs["area"]

        nbin_source = map_file.file['maps'].attrs['nbin_source']

        # Choose pixelization and read mask and systematics maps
        pixel_scheme = choose_pixelization(**pix_info)

        # Load the mask. It should automatically be the same shape as the
        # others, based on how it was originally generated.
        # We remove any pixels that are at or below our threshold (default=0)
        mask = map_file.read_map('mask')
        mask[np.isnan(mask)] = 0
        mask[mask==healpy.UNSEEN] = 0

        # Load all the maps in.
        # TODO: make it possible to just do a subset of these
        ngal_maps = [map_file.read_map(f'ngal_{b}') for b in range(nbin_lens)]
        g1_maps = [map_file.read_map(f'g1_{b}') for b in range(nbin_source)]
        g2_maps = [map_file.read_map(f'g2_{b}') for b in range(nbin_source)]
        lensing_weights = [map_file.read_map(f'lensing_weight_{b}') for b in range(nbin_source)]

        # Mask any pixels which have the healpix bad value
        for (g1, g2, lw) in zip(g1_maps, g2_maps, lensing_weights):
            lw[g1 == healpy.UNSEEN] = 0
            lw[g2 == healpy.UNSEEN] = 0

        # When running on the CosmoDC2 mocks I've found I need to flip
        # both g1 and g2 in order to get both positive galaxy-galaxy lensing
        # and shear-shear spectra.
        if self.config['flip_g1']:
            for g1 in g1_maps:
                w = np.where(g1!=healpy.UNSEEN)
                g1[w]*=-1

        if self.config['flip_g2']:
            for g2 in g2_maps:
                w = np.where(g2!=healpy.UNSEEN)
                g2[w]*=-1

        return [np.stack([g1,g2], axis=0) for (g1,g2) in zip(g1_maps, g2_maps)], lensing_weights

if __name__ == '__main__':
    cls = PipelineStage.main()
