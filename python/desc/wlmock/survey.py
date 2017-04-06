import h5py

from .subfield import subfield

class survey(object):
    """
    Class representing the HSC survey, defining the list of subfields,
    size of survey, etc...
    """

    def __init__(self, subfield_dict):
        """
        Initialise the survey class.
        Subfield_dict is a dictionary of subfield name and path to the
        data (for now, HSC randoms)
        """
        self.subfields = {}

        # Runs through the list of subfields
        for f in subfield_dict:
            # If a string is provided
            if isinstance(subfield_dict[f], basestring):
                self.subfields[f] = subfield(f,
                                             randoms_filename=subfield_dict[f])
            else:
                self.subfields[f] = subfield(f,
                                             dset=subfield_dict[f])

    def write(self, filename):
        """
        Saves the details of the surrvey to HDF5
        """
        # Opening output file
        f = h5py.File(filename, 'w')

        # Saving each subfield in own group
        for s in self.subfields:
            self.subfields[s].write(f)

        # closing
        f.close()

    @classmethod
    def read(cls, filename):
        """
        Reads  in a sruvey from a file
        """
        F = h5py.File(filename, 'r')

        subfield_dict = {}

        # Gp through the file and losd each dataset
        for f in F:
            subfield_dict[f] = F[f]

        return cls(subfield_dict)
