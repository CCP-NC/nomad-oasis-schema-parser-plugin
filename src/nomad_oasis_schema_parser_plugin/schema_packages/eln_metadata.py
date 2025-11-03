from nomad.datamodel.data import (
    ArchiveSection,
    EntryData,
)
from nomad.datamodel.metainfo.annotations import (
    ELNAnnotation,
    ELNComponentEnum,
)
from nomad.metainfo import (
    Quantity,
    SchemaPackage,
    Section,
    SubSection,
)

m_package = SchemaPackage()


class ORCIDInputELN(ArchiveSection):
    """Simple ORCID input section for ELN"""
    
    orcid_id = Quantity(
        type=str,
        description='ORCID identifier for the author',
        a_eln=ELNAnnotation(
            component=ELNComponentEnum.StringEditQuantity,
            label='ORCID ID',
        ),
    )


class CCPNCMetadataELN(EntryData):
    """
    ELN entry for CCPNC metadata input.
    Currently supports ORCID only.
    """
    
    m_def = Section(
        label='CCPNC Metadata Entry',
    )
    
    orcid = SubSection(
        section_def=ORCIDInputELN,
        description='Author ORCID information',
    )

m_package.__init_metainfo__()
