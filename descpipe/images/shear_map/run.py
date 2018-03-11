#!/usr/bin/env python3

"""
This pipeline step is trivial enough that we just code it up
in this run program rather than calling another script.
"""
# Make sure this file is executable.
import descpipe

class Stage(descpipe.Stage):
    name = "shear_map"
    config = {
        "config":"config.yaml"
    }

    inputs = {"shape_catalog":"fits"}

    outputs = {
        "shear_map": "fits",
    }


    def run(self):
        # Imports must be in here
        import yaml
        from desc.wlmassmap.shear_map import shear_map

        config_file = self.get_config_path("config")

        #Configuration options
        config = yaml.load(open(config_file))['shear_map']

        # Overwriting configuration with pipeline output
        config['input_filename'] = self.get_input_path("shape_catalog")
        config['output_filename'] = self.get_output_path("shear_map")

        # Execute the code
        shear_map(config)


# Always end with this
if __name__ == '__main__':
    Stage.main()
