from typing import (
    TYPE_CHECKING,
)

if TYPE_CHECKING:
    pass

import numpy as np
from ase.data import atomic_names, chemical_symbols
from nomad.config import config
from nomad.datamodel.data import ArchiveSection
from nomad.metainfo import JSON, Quantity, Reference, SchemaPackage, SubSection
from nomad.metainfo.elasticsearch_extension import Elasticsearch
from nomad_simulations.schema_packages.general import Program, Simulation
from nomad_simulations.schema_packages.model_method import DFT

from nomad_oasis_schema_parser_plugin.schema_packages.eln_metadata import (
    CCPNCMetadataELN,
)

configuration = config.get_plugin_entry_point(
    'nomad_oasis_schema_parser_plugin.schema_packages:ccpnc_schema_entry_point'
)

m_package = SchemaPackage()


class CCPNCProgram(Program):
    """Program section with descriptions tailored to the CCPNC/magres NMR
    calculation context. The base `Program.name`/`version` descriptions from
    nomad_simulations ('The name of the program.') are too generic to guide
    users browsing this filter menu.
    """

    name = Quantity(
        type=str,
        description="""
        Name of the DFT/simulation code used to run the NMR calculation that
        produced this magres file, e.g. 'CASTEP' or 'QuantumESPRESSO'.
        """,
    )

    version = Quantity(
        type=str,
        description="""
        Version string of the simulation code, as reported in the magres
        file's calculation metadata (e.g. '19.1').
        """,
    )


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


# New section for element-resolved magnetic shielding isotropy values
class ElementIsotropyEntry(ArchiveSection):
    element = Quantity(type=str, description="Element symbol, e.g. 'H', 'C', 'O'.")
    isotropy = Quantity(
        type=float, 
        unit='ppm',
        a_eln=dict(defaultDisplayUnit='ppm'),
        description='Shielding isotropy value for an atomic site.'
    )


# New section for element-resolved electric field gradient Vzz values
class ElementVzzEntry(ArchiveSection):
    element = Quantity(type=str, description="Element symbol, e.g. 'H', 'C', 'O'.")
    Vzz = Quantity(
        type=float, 
        unit='a_u_efg',
        description='Electric field gradient Vzz value for an atomic site.'
    )


class IsotropyEntry(ArchiveSection):
    isotropy = Quantity(
        type=float, 
        unit='ppm',
        a_eln=dict(defaultDisplayUnit='ppm'),
        description='Shielding isotropy value for any element'
    )


class VzzEntry(ArchiveSection):
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


class ElementResolvedMagneticShielding(ArchiveSection):
    element_isotropy_list = SubSection(
        section_def=ElementIsotropyEntry,
        repeats=True,
        description='List of element/isotropy entries.',
    )
    Al_isotropy_list = SubSection(
        section_def=ISOTROPY_ENTRY_CLASSES['Al'],
        repeats=True,
        description='List of shielding isotropy values for Al',
    )
    B_isotropy_list = SubSection(
        section_def=ISOTROPY_ENTRY_CLASSES['B'],
        repeats=True,
        description='List of shielding isotropy values for B',
    )
    Ba_isotropy_list = SubSection(
        section_def=ISOTROPY_ENTRY_CLASSES['Ba'],
        repeats=True,
        description='List of shielding isotropy values for Ba',
    )
    Bi_isotropy_list = SubSection(
        section_def=ISOTROPY_ENTRY_CLASSES['Bi'],
        repeats=True,
        description='List of shielding isotropy values for Bi',
    )
    Br_isotropy_list = SubSection(
        section_def=ISOTROPY_ENTRY_CLASSES['Br'],
        repeats=True,
        description='List of shielding isotropy values for Br',
    )
    C_isotropy_list = SubSection(
        section_def=ISOTROPY_ENTRY_CLASSES['C'],
        repeats=True,
        description='List of shielding isotropy values for C',
    )
    Cd_isotropy_list = SubSection(
        section_def=ISOTROPY_ENTRY_CLASSES['Cd'],
        repeats=True,
        description='List of shielding isotropy values for Cd',
    )
    Cl_isotropy_list = SubSection(
        section_def=ISOTROPY_ENTRY_CLASSES['Cl'],
        repeats=True,
        description='List of shielding isotropy values for Cl',
    )
    Cr_isotropy_list = SubSection(
        section_def=ISOTROPY_ENTRY_CLASSES['Cr'],
        repeats=True,
        description='List of shielding isotropy values for Cr',
    )
    Cs_isotropy_list = SubSection(
        section_def=ISOTROPY_ENTRY_CLASSES['Cs'],
        repeats=True,
        description='List of shielding isotropy values for Cs',
    )
    Cu_isotropy_list = SubSection(
        section_def=ISOTROPY_ENTRY_CLASSES['Cu'],
        repeats=True,
        description='List of shielding isotropy values for Cu',
    )
    F_isotropy_list = SubSection(
        section_def=ISOTROPY_ENTRY_CLASSES['F'],
        repeats=True,
        description='List of shielding isotropy values for F',
    )
    Fe_isotropy_list = SubSection(
        section_def=ISOTROPY_ENTRY_CLASSES['Fe'],
        repeats=True,
        description='List of shielding isotropy values for Fe',
    )
    Ga_isotropy_list = SubSection(
        section_def=ISOTROPY_ENTRY_CLASSES['Ga'],
        repeats=True,
        description='List of shielding isotropy values for Ga',
    )
    H_isotropy_list = SubSection(
        section_def=ISOTROPY_ENTRY_CLASSES['H'],
        repeats=True,
        description='List of shielding isotropy values for H',
    )
    Hf_isotropy_list = SubSection(
        section_def=ISOTROPY_ENTRY_CLASSES['Hf'],
        repeats=True,
        description='List of shielding isotropy values for Hf',
    )
    I_isotropy_list = SubSection(
        section_def=ISOTROPY_ENTRY_CLASSES['I'],
        repeats=True,
        description='List of shielding isotropy values for I',
    )
    In_isotropy_list = SubSection(
        section_def=ISOTROPY_ENTRY_CLASSES['In'],
        repeats=True,
        description='List of shielding isotropy values for In',
    )
    La_isotropy_list = SubSection(
        section_def=ISOTROPY_ENTRY_CLASSES['La'],
        repeats=True,
        description='List of shielding isotropy values for La',
    )
    Li_isotropy_list = SubSection(
        section_def=ISOTROPY_ENTRY_CLASSES['Li'],
        repeats=True,
        description='List of shielding isotropy values for Li',
    )
    Mg_isotropy_list = SubSection(
        section_def=ISOTROPY_ENTRY_CLASSES['Mg'],
        repeats=True,
        description='List of shielding isotropy values for Mg',
    )
    N_isotropy_list = SubSection(
        section_def=ISOTROPY_ENTRY_CLASSES['N'],
        repeats=True,
        description='List of shielding isotropy values for N',
    )
    Na_isotropy_list = SubSection(
        section_def=ISOTROPY_ENTRY_CLASSES['Na'],
        repeats=True,
        description='List of shielding isotropy values for Na',
    )
    O_isotropy_list = SubSection(
        section_def=ISOTROPY_ENTRY_CLASSES['O'],
        repeats=True,
        description='List of shielding isotropy values for O',
    )
    P_isotropy_list = SubSection(
        section_def=ISOTROPY_ENTRY_CLASSES['P'],
        repeats=True,
        description='List of shielding isotropy values for P',
    )
    S_isotropy_list = SubSection(
        section_def=ISOTROPY_ENTRY_CLASSES['S'],
        repeats=True,
        description='List of shielding isotropy values for S',
    )
    Sb_isotropy_list = SubSection(
        section_def=ISOTROPY_ENTRY_CLASSES['Sb'],
        repeats=True,
        description='List of shielding isotropy values for Sb',
    )
    Sc_isotropy_list = SubSection(
        section_def=ISOTROPY_ENTRY_CLASSES['Sc'],
        repeats=True,
        description='List of shielding isotropy values for Sc',
    )
    Se_isotropy_list = SubSection(
        section_def=ISOTROPY_ENTRY_CLASSES['Se'],
        repeats=True,
        description='List of shielding isotropy values for Se',
    )
    Si_isotropy_list = SubSection(
        section_def=ISOTROPY_ENTRY_CLASSES['Si'],
        repeats=True,
        description='List of shielding isotropy values for Si',
    )
    Sn_isotropy_list = SubSection(
        section_def=ISOTROPY_ENTRY_CLASSES['Sn'],
        repeats=True,
        description='List of shielding isotropy values for Sn',
    )
    Sr_isotropy_list = SubSection(
        section_def=ISOTROPY_ENTRY_CLASSES['Sr'],
        repeats=True,
        description='List of shielding isotropy values for Sr',
    )
    Ta_isotropy_list = SubSection(
        section_def=ISOTROPY_ENTRY_CLASSES['Ta'],
        repeats=True,
        description='List of shielding isotropy values for Ta',
    )
    Te_isotropy_list = SubSection(
        section_def=ISOTROPY_ENTRY_CLASSES['Te'],
        repeats=True,
        description='List of shielding isotropy values for Te',
    )
    Ti_isotropy_list = SubSection(
        section_def=ISOTROPY_ENTRY_CLASSES['Ti'],
        repeats=True,
        description='List of shielding isotropy values for Ti',
    )
    V_isotropy_list = SubSection(
        section_def=ISOTROPY_ENTRY_CLASSES['V'],
        repeats=True,
        description='List of shielding isotropy values for V',
    )
    Y_isotropy_list = SubSection(
        section_def=ISOTROPY_ENTRY_CLASSES['Y'],
        repeats=True,
        description='List of shielding isotropy values for Y',
    )
    Zn_isotropy_list = SubSection(
        section_def=ISOTROPY_ENTRY_CLASSES['Zn'],
        repeats=True,
        description='List of shielding isotropy values for Zn',
    )
    Zr_isotropy_list = SubSection(
        section_def=ISOTROPY_ENTRY_CLASSES['Zr'],
        repeats=True,
        description='List of shielding isotropy values for Zr',
    )


class ElementResolvedElectricFieldGradient(ArchiveSection):
    element_vzz_list = SubSection(
        section_def=ElementVzzEntry,
        repeats=True,
        description='List of element/vzz entries.',
    )
    Al_vzz_list = SubSection(
        section_def=VZZ_ENTRY_CLASSES['Al'],
        repeats=True,
        description='List of EFG Vzz values for Al',
    )
    B_vzz_list = SubSection(
        section_def=VZZ_ENTRY_CLASSES['B'],
        repeats=True,
        description='List of EFG Vzz values for B',
    )
    Ba_vzz_list = SubSection(
        section_def=VZZ_ENTRY_CLASSES['Ba'],
        repeats=True,
        description='List of EFG Vzz values for Ba',
    )
    Bi_vzz_list = SubSection(
        section_def=VZZ_ENTRY_CLASSES['Bi'],
        repeats=True,
        description='List of EFG Vzz values for Bi',
    )
    Br_vzz_list = SubSection(
        section_def=VZZ_ENTRY_CLASSES['Br'],
        repeats=True,
        description='List of EFG Vzz values for Br',
    )
    C_vzz_list = SubSection(
        section_def=VZZ_ENTRY_CLASSES['C'],
        repeats=True,
        description='List of EFG Vzz values for C',
    )
    Cd_vzz_list = SubSection(
        section_def=VZZ_ENTRY_CLASSES['Cd'],
        repeats=True,
        description='List of EFG Vzz values for Cd',
    )
    Cl_vzz_list = SubSection(
        section_def=VZZ_ENTRY_CLASSES['Cl'],
        repeats=True,
        description='List of EFG Vzz values for Cl',
    )
    Cr_vzz_list = SubSection(
        section_def=VZZ_ENTRY_CLASSES['Cr'],
        repeats=True,
        description='List of EFG Vzz values for Cr',
    )
    Cs_vzz_list = SubSection(
        section_def=VZZ_ENTRY_CLASSES['Cs'],
        repeats=True,
        description='List of EFG Vzz values for Cs',
    )
    Cu_vzz_list = SubSection(
        section_def=VZZ_ENTRY_CLASSES['Cu'],
        repeats=True,
        description='List of EFG Vzz values for Cu',
    )
    F_vzz_list = SubSection(
        section_def=VZZ_ENTRY_CLASSES['F'],
        repeats=True,
        description='List of EFG Vzz values for F',
    )
    Fe_vzz_list = SubSection(
        section_def=VZZ_ENTRY_CLASSES['Fe'],
        repeats=True,
        description='List of EFG Vzz values for Fe',
    )
    Ga_vzz_list = SubSection(
        section_def=VZZ_ENTRY_CLASSES['Ga'],
        repeats=True,
        description='List of EFG Vzz values for Ga',
    )
    H_vzz_list = SubSection(
        section_def=VZZ_ENTRY_CLASSES['H'],
        repeats=True,
        description='List of EFG Vzz values for H',
    )
    Hf_vzz_list = SubSection(
        section_def=VZZ_ENTRY_CLASSES['Hf'],
        repeats=True,
        description='List of EFG Vzz values for Hf',
    )
    I_vzz_list = SubSection(
        section_def=VZZ_ENTRY_CLASSES['I'],
        repeats=True,
        description='List of EFG Vzz values for I',
    )
    In_vzz_list = SubSection(
        section_def=VZZ_ENTRY_CLASSES['In'],
        repeats=True,
        description='List of EFG Vzz values for In',
    )
    La_vzz_list = SubSection(
        section_def=VZZ_ENTRY_CLASSES['La'],
        repeats=True,
        description='List of EFG Vzz values for La',
    )
    Li_vzz_list = SubSection(
        section_def=VZZ_ENTRY_CLASSES['Li'],
        repeats=True,
        description='List of EFG Vzz values for Li',
    )
    Mg_vzz_list = SubSection(
        section_def=VZZ_ENTRY_CLASSES['Mg'],
        repeats=True,
        description='List of EFG Vzz values for Mg',
    )
    N_vzz_list = SubSection(
        section_def=VZZ_ENTRY_CLASSES['N'],
        repeats=True,
        description='List of EFG Vzz values for N',
    )
    Na_vzz_list = SubSection(
        section_def=VZZ_ENTRY_CLASSES['Na'],
        repeats=True,
        description='List of EFG Vzz values for Na',
    )
    O_vzz_list = SubSection(
        section_def=VZZ_ENTRY_CLASSES['O'],
        repeats=True,
        description='List of EFG Vzz values for O',
    )
    P_vzz_list = SubSection(
        section_def=VZZ_ENTRY_CLASSES['P'],
        repeats=True,
        description='List of EFG Vzz values for P',
    )
    S_vzz_list = SubSection(
        section_def=VZZ_ENTRY_CLASSES['S'],
        repeats=True,
        description='List of EFG Vzz values for S',
    )
    Sb_vzz_list = SubSection(
        section_def=VZZ_ENTRY_CLASSES['Sb'],
        repeats=True,
        description='List of EFG Vzz values for Sb',
    )
    Sc_vzz_list = SubSection(
        section_def=VZZ_ENTRY_CLASSES['Sc'],
        repeats=True,
        description='List of EFG Vzz values for Sc',
    )
    Se_vzz_list = SubSection(
        section_def=VZZ_ENTRY_CLASSES['Se'],
        repeats=True,
        description='List of EFG Vzz values for Se',
    )
    Si_vzz_list = SubSection(
        section_def=VZZ_ENTRY_CLASSES['Si'],
        repeats=True,
        description='List of EFG Vzz values for Si',
    )
    Sn_vzz_list = SubSection(
        section_def=VZZ_ENTRY_CLASSES['Sn'],
        repeats=True,
        description='List of EFG Vzz values for Sn',
    )
    Sr_vzz_list = SubSection(
        section_def=VZZ_ENTRY_CLASSES['Sr'],
        repeats=True,
        description='List of EFG Vzz values for Sr',
    )
    Ta_vzz_list = SubSection(
        section_def=VZZ_ENTRY_CLASSES['Ta'],
        repeats=True,
        description='List of EFG Vzz values for Ta',
    )
    Te_vzz_list = SubSection(
        section_def=VZZ_ENTRY_CLASSES['Te'],
        repeats=True,
        description='List of EFG Vzz values for Te',
    )
    Ti_vzz_list = SubSection(
        section_def=VZZ_ENTRY_CLASSES['Ti'],
        repeats=True,
        description='List of EFG Vzz values for Ti',
    )
    V_vzz_list = SubSection(
        section_def=VZZ_ENTRY_CLASSES['V'],
        repeats=True,
        description='List of EFG Vzz values for V',
    )
    Y_vzz_list = SubSection(
        section_def=VZZ_ENTRY_CLASSES['Y'],
        repeats=True,
        description='List of EFG Vzz values for Y',
    )
    Zn_vzz_list = SubSection(
        section_def=VZZ_ENTRY_CLASSES['Zn'],
        repeats=True,
        description='List of EFG Vzz values for Zn',
    )
    Zr_vzz_list = SubSection(
        section_def=VZZ_ENTRY_CLASSES['Zr'],
        repeats=True,
        description='List of EFG Vzz values for Zr',
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
    program = SubSection(sub_section=CCPNCProgram.m_def, repeats=False)

    # New subsection for element-resolved magnetic shielding
    element_resolved_nmr_search = SubSection(section_def=ElementResolvedNMRSearch)

    # Add reference to the metadata ELN entry
    metadata_eln_reference = Quantity(
        type=Reference(CCPNCMetadataELN.m_def),
        description='Reference to the external metadata ELN entry',
    )


m_package.__init_metainfo__()
