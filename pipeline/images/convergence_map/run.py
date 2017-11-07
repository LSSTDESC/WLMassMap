#!/usr/bin/env python3

"""
This pipeline step is trivial enough that we just code it up
in this run program rather than calling another script.
"""
# Make sure this file is executable.
import descpipe

class Stage(descpipe.Stage):
    name = "convergence_map"
    config = {
        "config":"config.yaml"
    }

    inputs = {"shear_map":"fits"}

    outputs = {"convergence_map": "fits"}


    def run(self):
        # Imports must be in here
        import yaml
        from desc.wlmassmap.convergence_map import convergence_map

        config_file = self.get_config_path("config")

        #Configuration options
        config = yaml.load(open(config_file))['convergence_map']

        # Overwriting configuration with pipeline output
        config['input_filename'] = self.get_input_path("shear_map")
        config['output_filename'] = self.get_output_path("convergence_map")
        # Execute the code
        convergence_map(config)


# Always end with this
if __name__ == '__main__':
    Stage.main()
