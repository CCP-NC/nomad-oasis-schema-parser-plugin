import os

from nomad.datamodel.metainfo.eln import Entity
from nomad.metainfo import Quantity


class MetadataVoilaNotebook(Entity):
    notebook_file = Quantity(
        type=str,
        a_eln=dict(component='FileEditQuantity'),
        a_browser=dict(adaptor='RawFileAdaptor'),
    )

    tags = Quantity(
        type=str,
        shape=['*'],
        description='Add a tag that can be used for search.',
        a_eln=dict(component='StringEditQuantity'),
    )

    def normalize(self, archive, logger):
        super().normalize(archive, logger)

        if self.notebook_file and os.path.splitext(self.notebook_file)[-1] != '.ipynb':
            logger.error('Please upload a jupyter notebook file (.ipynb).')