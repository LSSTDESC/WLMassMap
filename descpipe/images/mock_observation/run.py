#!/usr/bin/env python3

"""
This pipeline step is trivial enough that we just code it up
in this run program rather than calling another script.
"""
# Make sure this file is executable.
import descpipe

class Stage(descpipe.Stage):
    name = "mock_observation"
    config = {
        "config":"config.yaml"
    }

    inputs = {"ground_truth": "hdf5"}

    outputs = {
        "shape_catalog": "fits"
    }


    def run(self):
        # Imports must be in here
        import yaml
        from desc.wlmassmap.mocks.mock_observation import mock_observation

        config_file = self.get_config_path("config")
        #Configuration options
        config = yaml.load(open(config_file))['mock_observation']

        # Overwriting configuration with pipeline output
        config['input_filename'] = self.get_input_path("ground_truth")
        config['output_filename'] = self.get_output_path("shape_catalog")

        # Execute the code
        mock_observation(config)


# Always end with this
if __name__ == '__main__':
    Stage.main()
