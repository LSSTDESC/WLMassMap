#!/usr/bin/env python3

"""
This pipeline step is trivial enough that we just code it up
in this run program rather than calling another script.
"""
# Make sure this file is executable.
import descpipe

class Stage(descpipe.Stage):
    name = "extract_footprint"
    config = {
        "config":"config.yaml"
    }

    inputs = {}

    outputs = {
        "ground_truth": "hdf5",
    }


    def run(self):
        # Imports must be in here
        import yaml
        from desc.wlmassmap.mocks.extract_footprint import extract_footprint

        config_file = self.get_config_path("config")
        output_file = self.get_output_path("ground_truth")

        #Configuration options
        config = yaml.load(open(config_file))['extract_footprint']

        # Overwriting configuration with pipeline output
        config['output_filename'] = output_file
        print(config)
        # Execute the code
        extract_footprint(config)


# Always end with this
if __name__ == '__main__':
    Stage.main()
