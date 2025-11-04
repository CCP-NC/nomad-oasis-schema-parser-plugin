from typing import (
    TYPE_CHECKING,
)

if TYPE_CHECKING:
    pass

import numpy as np
from nomad.config import config
from nomad.datamodel.data import ArchiveSection, EntryData
from nomad.datamodel.metainfo.annotations import ELNAnnotation
from nomad.metainfo import JSON, Quantity, Reference, SchemaPackage, SubSection
from nomad_simulations.schema_packages.general import Simulation

from nomad_oasis_schema_parser_plugin.schema_packages.eln_metadata import (
    CCPNCMetadataELN,
)
from nomad_oasis_schema_parser_plugin.schema_packages.metadata_voila import (
    MetadataVoilaNotebook,
)

configuration = config.get_plugin_entry_point(
    'nomad_oasis_schema_parser_plugin.schema_packages:ccpnc_schema_entry_point'
)

m_package = SchemaPackage()

class MaterialProperties(ArchiveSection):
    # Note from @JosePizarro3: note we have all these information somewhere else in the 
    # `nomad_simulations` schema. Nevertheless, if you feel it is better to keep these 
    # quantities here for clarity, it is totally fine for me.
    chemical_name = Quantity(
        type=str,
        description="""
        Free-text chemical name assigned by users.
        """,
    )

    chemical_name_tokens = Quantity(
        type=str,
        shape=['*'],
        description="""
        Free-text chemical name, but tokenised to take individual words in the name to 
        assist in wildcard searches.
        """,
    )

    formula = Quantity(
        # type=[(str, int)],  # better to use JSON type
        type=JSON,
        shape=['*'],
        description="""
        Dictionary containing the species (chemical symbol of an element in the 
        material) as keys and number of atoms of that element in the material as their 
        value.
        """,
    )

    stoichiometry = Quantity(
        type=JSON,
        shape=['*'],
        description="""
        Reduced proportion of materials details.
        """,
    )

    elements_ratios = Quantity(
        type=np.float64,
        shape=['*'],
        description="""
        Ratio of constituent elements (each element is a number between 0 and 1).
        """,
    )

class PublicationRecord(ArchiveSection):
    doi = Quantity(
        type=str,
        description="""
        Digital Object Identifier if results are part of a publication.
        """,
    )


class ORCID(ArchiveSection):
    orcid_id = Quantity(
        type=str,
        description="""
        Dictionary containing the ORCID IDs of the author (keys) and uploader (values) 
        profiles.
        """,
    )


class CCPNCRecord(ArchiveSection):
    visible = Quantity(
        type=bool,
        description="""
        A boolean value that indicates if the record is to be hidden or available to be 
        returned when searched.
        """,
    )

    immutable_id = Quantity(
        type=str,
        description="""
        7 digit unique record identifier.
        """,
    )

    license = Quantity(
        type=str,
        description="""
        License under which the record is released.
        """,
    )


class ExternalDatabaseReference(ArchiveSection):
    external_database_name = Quantity(
        type=str,
        description="""
        External database name where additional information on the material exists
        """,
    )

    external_database_name_other = Quantity(
        type=str,
        description="""
        External database name where additional information on the material exists.
        Use this field if 'Other' is selected in the 'external_database_name' field.
        """,
    )

    external_database_reference_code = Quantity(
        type=str,
        description="""
        Specific database code pointing to the material or a polymorphic form of the 
        material.
        """,
    )


class FreeTextMetadata(ArchiveSection):
    uploader_author_notes = Quantity(
        type=str,
        description="""
        Additional metadata that authors want to indicate about the computation.
        """,
    )

    structural_descriptor_notes = Quantity(
        type=str,
        description="""
        Additional notes specific to the polymorphic forms of the material.
        """,
    )

class CCPNC_VoilaNotebook(MetadataVoilaNotebook, EntryData):
    # m_def = Section(a_eln=dict(hide=['lab_id']))

    def normalize(self, archive, logger):
        super().normalize(archive, logger)


class CCPNCMetadata(ArchiveSection):
    material_properties = SubSection(section_def=MaterialProperties)
    orcid = SubSection(section_def=ORCID)
    ccpnc_record = SubSection(section_def=CCPNCRecord)
    external_database_reference = SubSection(section_def=ExternalDatabaseReference)
    free_text_metadata = SubSection(section_def=FreeTextMetadata)
    publication_record = SubSection(section_def=PublicationRecord)

# Define the CCPNCSimulation class holding CCP-NC specific metadata
class CCPNCSimulation(Simulation):
    ccpnc_metadata = SubSection(section_def=CCPNCMetadata)

    # Add reference to the metadata ELN entry
    metadata_eln_reference = Quantity(
        type=Reference(CCPNCMetadataELN.m_def),
        description='Reference to the external metadata ELN entry',
    )

m_package.__init_metainfo__()
