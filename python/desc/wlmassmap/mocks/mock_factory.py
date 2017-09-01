from abc import ABCMeta, abstractmethod
from astropy.extern import six

@six.add_metaclass(ABCMeta)
class MockFactory(object):
    """
    Abstract base class for mock factories. Implementations of this class should

    """

    __baseFields = ['galaxy_id',                        # Unique galaxy id
                    'ra', 'ra_true',                    # Magnified and
                    'dec', 'dec_true',                  # non Magnified coordinates
                    'redshift',                         # Cosmological + motion
                    'ellipticity_1', 'ellipticity_2',
                    'ellipticity_1_true', 'ellipticity_2_true',
                    'shear_1', 'shear_2', 'convergence',
                    'reduced_shear_1', 'reduced_shear_2']

    def __init__(self, center=None):
        """
        Mock factory initialisation

        Parameters
        ----------
            center: tuple of float
                (Ra, Dec) center of the survey
        """
        self.center = center

    @abstractmethod
    def _populate(self, **kwargs):
        """
        Generate galaxy positions and shapes, needs to be implemented for
        specific providers
        """
        raise NotImplementedError("All subclasses of MockFactory"
        " must include a _populate method")

    def generate(self, **kwargs):
        """
        Generates a shape catalog, following the parameters of the model

        Returns
        -------
        cat: Table
            Realisation of the shape catalog.
        """

        # Sample galaxy positions and shapes
        cat = self._populate()
