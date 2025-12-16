import csv
import json
import os

import numpy as np

from nomad_oasis_schema_parser_plugin.schema_packages.schema_package import (
    ORCID,
    CCPNCMetadata,
    CCPNCRecord,
    ExternalDatabaseReference,
    MaterialProperties,
    PublicationRecord,
)

DATA_DIR = os.path.join(os.path.dirname(__file__), '..', 'data')
EXPECTED_STOICHIOMETRY_COUNT = 21

def test_material_properties_fields():
    mat = MaterialProperties()
    mat.chemical_name = "Ethanol"
    mat.elements_ratios = np.array(
        [0.2222222222222222,0.6666666666666666,0.1111111111111111])
    assert mat.chemical_name == "Ethanol"
    assert np.allclose(
        mat.elements_ratios, [0.22222222, 0.66666667, 0.11111111], rtol=1e-6)

def test_ccpnc_metadata_from_json():
    json_path = os.path.join(DATA_DIR, "upload_single_file.json")
    with open(json_path) as f:
        metadata = json.load(f)

    ccpnc_metadata = CCPNCMetadata()
    # Populate fields as the parser would
    ccpnc_metadata.material_properties = MaterialProperties()
    ccpnc_metadata.material_properties.chemical_name = metadata["chemname"]
    ccpnc_metadata.material_properties.elements_ratios = np.array(
        metadata["elements_ratios"])
    ccpnc_metadata.material_properties.formula = metadata["formula"]
    ccpnc_metadata.material_properties.stoichiometry = metadata["stochiometry"]

    ccpnc_metadata.orcid = ORCID()
    ccpnc_metadata.orcid.orcid_id = metadata["ORCID"]

    ccpnc_metadata.ccpnc_record = CCPNCRecord()
    ccpnc_metadata.ccpnc_record.immutable_id = metadata["immutable_id"]

    ccpnc_metadata.external_database_reference = ExternalDatabaseReference()
    ccpnc_metadata.external_database_reference.external_database_name = (
        metadata["version_metadata"]["extref_type"]
    )
    ccpnc_metadata.external_database_reference.external_database_reference_code = (
        metadata["version_metadata"]["extref_code"]
    )

    ccpnc_metadata.publication_record = PublicationRecord()
    ccpnc_metadata.publication_record.doi = metadata["version_metadata"]["doi"]

    # Assertions
    assert ccpnc_metadata.material_properties.chemical_name.startswith("5-Hydroxy")
    assert ccpnc_metadata.orcid.orcid_id == "0000-0002-3377-6759"
    assert ccpnc_metadata.ccpnc_record.immutable_id == "0005631"
    assert ccpnc_metadata.external_database_reference.external_database_name == "csd"
    assert ccpnc_metadata.material_properties.formula[0]["species"] == "C"
    assert (ccpnc_metadata.material_properties.stoichiometry[0]["n"] == 
            EXPECTED_STOICHIOMETRY_COUNT)

def test_ccpnc_metadata_from_csv():
    csv_path = os.path.join(DATA_DIR, "metadata_info.csv")
    expected = {
        "multi_upload_file1.magres": {
            "doi": "10.1002/anie.201908914",
            "license": "pddl",
            "extref_type": "csd",
            "extref_code": "BINMEQ05",
            "chemname_prefix": "5-Hydroxy",
        },
        "multi_upload_file2.magres": {
            "doi": "10.1002/anie.201908914",
            "license": "pddl",
            "extref_type": "csd",
            "extref_code": "BINMEQ05",
            "chemname_prefix": "5-Hydroxy",
        },
    }
    found = set()
    with open(csv_path, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            fname = row["filename"]
            if fname in expected:
                exp = expected[fname]
                # MaterialProperties
                mat = MaterialProperties()
                mat.chemical_name = row["chemname"]
                assert mat.chemical_name.startswith(exp["chemname_prefix"])
                # DOI
                assert row["doi"] == exp["doi"]
                # License
                ccpnc_record = CCPNCRecord()
                ccpnc_record.license = row["license"]
                assert ccpnc_record.license == exp["license"]
                # ExternalDatabaseReference
                ext_db_ref = ExternalDatabaseReference()
                ext_db_ref.external_database_name = row["extref_type"]
                ext_db_ref.external_database_reference_code = row["extref_code"]
                assert ext_db_ref.external_database_name == exp["extref_type"]
                assert ext_db_ref.external_database_reference_code == exp["extref_code"]
                found.add(fname)
    for fname in expected:
        assert fname in found, f"Metadata for {fname} not found in CSV."
