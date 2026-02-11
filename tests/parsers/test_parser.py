import os

import structlog
from nomad.datamodel import EntryArchive

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
    assert archive.data.ccpnc_metadata.ccpnc_record.license == 'pddl'
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
    assert archive.data.ccpnc_metadata.ccpnc_record.license == 'pddl'
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
