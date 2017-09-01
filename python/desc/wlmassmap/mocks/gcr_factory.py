
from GCR import load_catalog

from .mock_factory import MockFactory

class GCRFactory(MockFactory):
    """
    Implementation of mock catalog factory based on the Generic Catalog Reader
    """

    def __init__(self, catalog_path, catalog_name, **kwargs):
        """

        Parameters
        ----------
        catalog_path: string
            Path to the catalog directory

        catalog_name: string
            Name of the
        """
        super(this.__class__, this).__init__(**kwargs)

        # Loads the catalog
        self.catalog = load_catalog(os.path.join(catalog_path, catalog_name))

    def _populate(self):
        """
        Read galaxy catalog using the GCR
        """

        cat = self.catalog.get_quantities(['galaxy_id'])
