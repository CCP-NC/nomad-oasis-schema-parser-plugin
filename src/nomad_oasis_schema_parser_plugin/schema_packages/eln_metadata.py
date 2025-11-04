from nomad.datamodel.data import (
    ArchiveSection,
    EntryData,
)
from nomad.datamodel.metainfo.annotations import (
    ELNAnnotation,
    ELNComponentEnum,
)
from nomad.metainfo import (
    MEnum,
    Quantity,
    SchemaPackage,
    Section,
)

m_package = SchemaPackage()


class CCPNCMetadataELN(EntryData):
    """
    ELN entry for CCPNC metadata input.
    All fields are at the top level for easy data entry.
    """
    
    m_def = Section(
        label='CCPNC Metadata Entry',
    )
    
    # ORCID field
    orcid_id = Quantity(
        type=str,
        description='ORCID identifier for the author',
        a_eln=ELNAnnotation(
            component=ELNComponentEnum.StringEditQuantity,
            label='Author\'s ORCID ID',
        ),
    )
    
    # Material Properties field
    chemical_name = Quantity(
        type=str,
        description='Free-text chemical name assigned by users',
        a_eln=ELNAnnotation(
            component=ELNComponentEnum.StringEditQuantity,
            label='Chemical Name',
        ),
    )

    # Publication license field - DROP-DOWN
    license = Quantity(
        type=MEnum([
            'PDDL v1.0',
            'ODC BY v1.0',
            'CC BY 4.0',
        ]),
        description='License under which the record is released.',
        a_eln=ELNAnnotation(
            component=ELNComponentEnum.EnumEditQuantity,
            label='License',
        ),
    )

    # Publication DOI field
    doi = Quantity(
        type=str,
        description='Digital Object Identifier if results are part of a publication.',
        a_eln=ELNAnnotation(
            component=ELNComponentEnum.StringEditQuantity,
            label='Publication DOI',
        ),
    )

    # External Database Reference fields
    # Database name - DROP-DOWN
    external_database_name = Quantity(
        type=MEnum([
            'CSD',
            'ICSD',
            'COD',
            'Other',
            'N/A',
        ]),
        description='External database name where additional information on the material exists',
        a_eln=ELNAnnotation(
            component=ELNComponentEnum.EnumEditQuantity,
            label='External Database Name',
        ),
    )

    # Other database name
    external_database_name_other = Quantity(
        type=str,
        description='External database name where additional information on the material exists.',
        a_eln=ELNAnnotation(
            component=ELNComponentEnum.StringEditQuantity,
            label='Other External Database Name',
        ),
    )

    # Database reference code
    external_database_reference_code = Quantity(
        type=str,
        description='Specific database code pointing to the material or a polymorphic form of the material.',
        a_eln=ELNAnnotation(
            component=ELNComponentEnum.StringEditQuantity,
            label='External Database Reference Code',
        ),
    )

    # Free text metadata fields
    structural_descriptor_notes = Quantity(
        type=str,
        description='Additional notes specific to the polymorphic forms of the material.',
        a_eln=ELNAnnotation(
            component=ELNComponentEnum.StringEditQuantity,
            label='Additional relevant descriptors (e.g.: phase, conformers, etc.)',
        ),
    )

    uploader_author_notes = Quantity(
        type=str,
        description='Additional metadata that authors want to indicate about the computation.',
        a_eln=ELNAnnotation(
            component=ELNComponentEnum.StringEditQuantity,
            label='Author\'s Notes',
        ),
    )

m_package.__init_metainfo__()
