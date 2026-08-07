from typing import (
    TYPE_CHECKING,
    NamedTuple,
)

if TYPE_CHECKING:
    pass

import numpy as np
from ase.data import atomic_names, chemical_symbols
from nomad.config import config
from nomad.datamodel.data import ArchiveSection
from nomad.metainfo import JSON, Quantity, Reference, SchemaPackage, Section, SubSection
from nomad.metainfo.elasticsearch_extension import Elasticsearch
from nomad_simulations.schema_packages.general import Simulation
from nomad_simulations.schema_packages.model_method import DFT

from nomad_oasis_schema_parser_plugin.schema_packages.eln_metadata import (
    CCPNCMetadataELN,
)

configuration = config.get_plugin_entry_point(
    'nomad_oasis_schema_parser_plugin.schema_packages:ccpnc_schema_entry_point'
)

m_package = SchemaPackage()


class MaterialProperties(ArchiveSection):
    chemical_name = Quantity(
        type=str,
        description="""
        Free-text chemical name assigned by users. Each distinct name
        (including differences in case or formatting) appears as its own
        checkbox option below — near-duplicate names are not merged.
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
        7-digit unique record identifier, zero-padded on the left. To search
        for a specific record, enter all 7 digits, e.g. record 1 must be
        entered as '0000001', not '1'. Search requires an exact match.
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
        External database where additional information on the material exists.
        Permitted values, shown as checkboxes below: 'csd' (Cambridge
        Structural Database), 'icsd' (Inorganic Crystal Structure Database),
        'cod' (Crystallography Open Database), 'other', or 'n/a'.
        """,
    )

    external_database_name_other = Quantity(
        type=str,
        description="""
        External database name where additional information on the material exists.
        Use this field if 'other' is selected in the 'external_database_name' field.
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


class CCPNCMetadata(ArchiveSection):
    material_properties = SubSection(section_def=MaterialProperties)
    orcid = SubSection(section_def=ORCID)
    ccpnc_record = SubSection(section_def=CCPNCRecord)
    external_database_reference = SubSection(section_def=ExternalDatabaseReference)
    free_text_metadata = SubSection(section_def=FreeTextMetadata)
    publication_record = SubSection(section_def=PublicationRecord)

_SITE_LABEL_DESCRIPTION = (
    'Atomic site label as it appears in the source magres file, e.g. \'H_1\'.'
)

# `m_def.more.label_quantity` (used by the GUI archive browser to label list
# items) and `key_quantity` (the backend equivalent, set on the SubSections
# below) both need to point at `site_label` for site labels to display correctly

# New section for element-resolved magnetic shielding isotropy values
class ElementIsotropyEntry(ArchiveSection):
    m_def = Section(label_quantity='site_label')

    element = Quantity(type=str, description="Element symbol, e.g. 'H', 'C', 'O'.")
    site_label = Quantity(type=str, description=_SITE_LABEL_DESCRIPTION)
    isotropy = Quantity(
        type=float, 
        unit='ppm',
        a_eln=dict(defaultDisplayUnit='ppm'),
        description='Shielding isotropy value for an atomic site.'
    )


# New section for element-resolved electric field gradient Vzz values
class ElementVzzEntry(ArchiveSection):
    m_def = Section(label_quantity='site_label')

    element = Quantity(type=str, description="Element symbol, e.g. 'H', 'C', 'O'.")
    site_label = Quantity(type=str, description=_SITE_LABEL_DESCRIPTION)
    Vzz = Quantity(
        type=float, 
        unit='a_u_efg',
        description='Electric field gradient Vzz value for an atomic site.'
    )


class IsotropyEntry(ArchiveSection):
    m_def = Section(label_quantity='site_label')

    site_label = Quantity(type=str, description=_SITE_LABEL_DESCRIPTION)
    isotropy = Quantity(
        type=float, 
        unit='ppm',
        a_eln=dict(defaultDisplayUnit='ppm'),
        description='Shielding isotropy value for any element'
    )


class VzzEntry(ArchiveSection):
    m_def = Section(label_quantity='site_label')

    site_label = Quantity(type=str, description=_SITE_LABEL_DESCRIPTION)
    Vzz = Quantity(
        type=float,
        unit='a_u_efg',
        description='Electric field gradient Vzz value for any element',
        a_elasticsearch=Elasticsearch(  # Add this annotation
            metrics={'min': 'min', 'max': 'max'}  # Enables min/max aggregations
        ),
    )


# Elements for which per-element magnetic shielding / electric field gradient
# subsections are exposed below. Kept as a single list so the per-element
# entry classes and the SubSection declarations that use them stay in sync.
_NMR_ELEMENT_SYMBOLS = [
    'Al', 'B', 'Ba', 'Bi', 'Br', 'C', 'Cd', 'Cl', 'Cr', 'Cs', 'Cu', 'F', 'Fe',
    'Ga', 'H', 'Hf', 'I', 'In', 'La', 'Li', 'Mg', 'N', 'Na', 'O', 'P', 'S',
    'Sb', 'Sc', 'Se', 'Si', 'Sn', 'Sr', 'Ta', 'Te', 'Ti', 'V', 'Y', 'Zn', 'Zr',
]


def _element_name(symbol: str) -> str:
    return atomic_names[chemical_symbols.index(symbol)]


def _make_isotropy_entry_class(symbol: str) -> type:
    """Create an IsotropyEntry subclass whose `isotropy` quantity carries an
    element-specific description, so the search filter tooltip for e.g. the
    Hydrogen isotropy list reads differently from the Carbon one."""
    element_name = _element_name(symbol)
    return type(
        f'{symbol}IsotropyEntry',
        (IsotropyEntry,),
        {
            '__module__': __name__,
            # A subclass's auto-generated `m_def` does not inherit the
            # parent's `label_quantity`, so it must be set again here.
            'm_def': Section(label_quantity='site_label'),
            'isotropy': Quantity(
                type=float,
                unit='ppm',
                a_eln=dict(defaultDisplayUnit='ppm'),
                description=(
                    f'Isotropic magnetic shielding value for {element_name} '
                    f'({symbol}) atomic sites.'
                ),
            ),
        },
    )


def _make_vzz_entry_class(symbol: str) -> type:
    """Create a VzzEntry subclass whose `Vzz` quantity carries an
    element-specific description, analogous to `_make_isotropy_entry_class`."""
    element_name = _element_name(symbol)
    return type(
        f'{symbol}VzzEntry',
        (VzzEntry,),
        {
            '__module__': __name__,
            # A subclass's auto-generated `m_def` does not inherit the
            # parent's `label_quantity`, so it must be set again here.
            'm_def': Section(label_quantity='site_label'),
            'Vzz': Quantity(
                type=float,
                unit='a_u_efg',
                description=(
                    f'Electric field gradient Vzz value for {element_name} '
                    f'({symbol}) atomic sites.'
                ),
                a_elasticsearch=Elasticsearch(
                    metrics={'min': 'min', 'max': 'max'}
                ),
            ),
        },
    )


ISOTROPY_ENTRY_CLASSES = {
    symbol: _make_isotropy_entry_class(symbol) for symbol in _NMR_ELEMENT_SYMBOLS
}
VZZ_ENTRY_CLASSES = {
    symbol: _make_vzz_entry_class(symbol) for symbol in _NMR_ELEMENT_SYMBOLS
}


class _GenericEntry(NamedTuple):
    """The combined (all-elements) list shown alongside the per-element ones."""

    attr_name: str
    section_def: type
    description: str


class _PerElementEntries(NamedTuple):
    """Per-element list config. `description_template` takes `{symbol}`."""

    attr_suffix: str
    classes: dict
    description_template: str


def _make_element_resolved_class(
    class_name: str, generic: _GenericEntry, per_element: _PerElementEntries
) -> type:
    """Build an element-resolved container class (magnetic shielding or EFG)
    with one SubSection per element in `per_element.classes`, instead of
    hand-writing one attribute per element.
    """
    namespace = {
        '__module__': __name__,
        generic.attr_name: SubSection(
            section_def=generic.section_def,
            repeats=True,
            description=generic.description,
            key_quantity='site_label',
        ),
    }
    for symbol, entry_class in per_element.classes.items():
        namespace[f'{symbol}{per_element.attr_suffix}'] = SubSection(
            section_def=entry_class,
            repeats=True,
            description=per_element.description_template.format(symbol=symbol),
            key_quantity='site_label',
        )
    return type(class_name, (ArchiveSection,), namespace)


ElementResolvedMagneticShielding = _make_element_resolved_class(
    'ElementResolvedMagneticShielding',
    _GenericEntry(
        'element_isotropy_list',
        ElementIsotropyEntry,
        'List of element/isotropy entries.',
    ),
    _PerElementEntries(
        '_isotropy_list',
        ISOTROPY_ENTRY_CLASSES,
        'List of shielding isotropy values for {symbol}',
    ),
)

ElementResolvedElectricFieldGradient = _make_element_resolved_class(
    'ElementResolvedElectricFieldGradient',
    _GenericEntry(
        'element_vzz_list', ElementVzzEntry, 'List of element/vzz entries.'
    ),
    _PerElementEntries(
        '_vzz_list', VZZ_ENTRY_CLASSES, 'List of EFG Vzz values for {symbol}'
    ),
)


class ElementResolvedNMRSearch(ArchiveSection):
    element_resolved_magnetic_shielding = SubSection(
        section_def=ElementResolvedMagneticShielding
    )

    element_resolved_electric_field_gradient = SubSection(
        section_def=ElementResolvedElectricFieldGradient
    )


# Define the CCPNCSimulation class holding CCP-NC specific metadata
class CCPNCSimulation(Simulation):
    ccpnc_metadata = SubSection(section_def=CCPNCMetadata)
    model_method = SubSection(sub_section=DFT.m_def, repeats=True)

    # New subsection for element-resolved magnetic shielding
    element_resolved_nmr_search = SubSection(section_def=ElementResolvedNMRSearch)

    # Add reference to the metadata ELN entry
    metadata_eln_reference = Quantity(
        type=Reference(CCPNCMetadataELN.m_def),
        description='Reference to the external metadata ELN entry',
    )


m_package.__init_metainfo__()
