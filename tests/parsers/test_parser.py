import os

# Import first to register the plugin normalizer entry points, avoiding
# a circular import when CCPNCNormalizer is imported directly below.
import nomad.normalizing  # noqa: F401
import structlog
from nomad.datamodel import EntryArchive

from nomad_oasis_schema_parser_plugin.normalizers.normalizer import CCPNCNormalizer
from nomad_oasis_schema_parser_plugin.parsers.parser import CCPNCMagresParser

DATA_DIR = os.path.join(os.path.dirname(__file__), '..', 'data')


def get_logger():
    return structlog.get_logger('test')


def test_parser_read_magres_json():
    """Test parser can read magres and json and populate metadata correctly."""
    parser = CCPNCMagresParser()
    archive = EntryArchive()
    logger = get_logger()

    # Test parsing when accompanying json file is found
    magres_path = os.path.join(DATA_DIR, 'qe_generated.magres')
    parser.parse(magres_path, archive, logger)
    # Check some key metadata fields
    assert archive.data.ccpnc_metadata.material_properties.chemical_name == 'Kaolinite'
    assert (
        archive.data.ccpnc_metadata.publication_record.doi
        == '10.1016/j.clay.2018.12.013'
    )
    # Raw source data has 'pddl'; normalize_license() maps it to the
    # canonical ELN dropdown label.
    assert archive.data.ccpnc_metadata.ccpnc_record.license == 'PDDL v1.0'
    assert (
        archive.data.ccpnc_metadata.external_database_reference.external_database_reference_code
        == '1896953'
    )

    # Check number of atoms, elements, or other quantities
    sim = archive.data
    mat = sim.ccpnc_metadata.material_properties
    TOLERANCE = 1e-6
    # For Kaolinite, check elements and ratios
    assert set(mat.formula[i]['species'] for i in range(len(mat.formula))) == {
        'Al',
        'H',
        'O',
        'Si',
    }
    assert abs(sum(mat.elements_ratios) - 1.0) < TOLERANCE


def test_parser_read_magres_csv():
    """Test parser can read magres and csv and populate metadata correctly."""
    parser = CCPNCMagresParser()
    archive = EntryArchive()
    logger = get_logger()

    # Test parsing when accompanying csv file is found (json not found)
    magres_path = os.path.join(DATA_DIR, 'multi_upload_file1.magres')
    parser.parse(magres_path, archive, logger)
    # Check some key metadata fields
    assert archive.data.ccpnc_metadata.material_properties.chemical_name.startswith(
        '5-Hydroxy'
    )
    assert (
        archive.data.ccpnc_metadata.publication_record.doi == '10.1002/anie.201908914'
    )
    # Raw source data has 'pddl'; normalize_license() maps it to the
    # canonical ELN dropdown label.
    assert archive.data.ccpnc_metadata.ccpnc_record.license == 'PDDL v1.0'
    assert (
        archive.data.ccpnc_metadata.external_database_reference.external_database_reference_code
        == 'BINMEQ05'
    )


def test_parser_handles_qe_generated_magres_with_missing_metadata():
    """
    Test that the parser correctly extracts program name and version from
    a QE-GIPAW magres file with minimal/incomplete calculation metadata.
    """
    DATA_DIR = os.path.join(os.path.dirname(__file__), '..', 'data')
    magres_path = os.path.join(DATA_DIR, 'qe_generated.magres')

    parser = CCPNCMagresParser()
    archive = EntryArchive()
    logger = get_logger()

    parser.parse(magres_path, archive, logger)

    # The parser should set the program name and version correctly
    program = archive.data.program
    assert program.name == 'Quantum ESPRESSO'
    # For this test file, version should be extracted from calc_code ("QE-GIPAW 5.x")
    assert program.version == '5.x'


def test_normalizer_sets_site_labels():
    """
    Site labels (e.g. 'H_1') are copied from the matched atom's `label` onto
    the element-resolved isotropy/Vzz entries during normalization. Check
    they survive both the combined list and the per-element grouping.
    """
    parser = CCPNCMagresParser()
    archive = EntryArchive()
    logger = get_logger()

    magres_path = os.path.join(DATA_DIR, 'multi_upload_file1.magres')
    parser.parse(magres_path, archive, logger)

    normalizer = CCPNCNormalizer()
    normalizer._populate_element_resolved_magnetic_shielding(archive, logger)
    normalizer._populate_element_resolved_electric_field_gradient(archive, logger)

    element_resolved = archive.data.element_resolved_nmr_search

    magnetic_shielding = element_resolved.element_resolved_magnetic_shielding
    combined_labels = [
        entry.site_label for entry in magnetic_shielding.element_isotropy_list[:3]
    ]
    per_element_labels = [
        entry.site_label for entry in magnetic_shielding.H_isotropy_list[:3]
    ]
    assert combined_labels == ['H_1', 'H_2', 'H_3']
    assert per_element_labels == ['H_1', 'H_2', 'H_3']

    electric_field_gradient = element_resolved.element_resolved_electric_field_gradient
    efg_labels = [entry.site_label for entry in electric_field_gradient.H_vzz_list[:3]]
    assert efg_labels == ['H_1', 'H_2', 'H_3']
