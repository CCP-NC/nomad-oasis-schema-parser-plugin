import json
import os
from typing import (
    TYPE_CHECKING,
)

if TYPE_CHECKING:
    from nomad.datamodel.datamodel import (
        EntryArchive,
    )
    from structlog.stdlib import (
        BoundLogger,
    )

from nomad.config import config
from nomad.datamodel import EntryArchive
from nomad.datamodel.metainfo.workflow import Workflow
# utility function used to get auxiliary files next to the `mainfile`
from nomad_oasis_schema_parser_plugin.parsers.utils.utils import get_files
from nomad_oasis_schema_parser_plugin.schema_packages.schema_package import (
    ORCID,
    CCPNCMetadata,
    CCPNCRecord,
    ExternalDatabaseReference,
    FreeTextMetadata,
    MaterialProperties,
)
from nomad_oasis_schema_parser_plugin.schema_packages.schema_package import (
    CCPNCSimulation as Simulation,
)

configuration = config.get_plugin_entry_point(
    'nomad_oasis_schema_parser_plugin.parsers:ccpnc_parser_entry_point'
)


class CCPNCMagresParser(MagresParser):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
    def test_file_discovery(self, logger: "BoundLogger") -> None:
        """Test method to check file discovery"""
        logger.info(f"Testing file discovery for mainfile: {self.mainfile}")
        logger.info(f"Maindir: {self.maindir}")
        logger.info(f"Basename: {self.basename}")
        
        # List all files in the directory
        import os
        try:
            files_in_dir = os.listdir(self.maindir)
            logger.info(f"Files in directory: {files_in_dir}")
            
            # Check for JSON files specifically
            json_files = [f for f in files_in_dir if f.endswith('.json')]
            logger.info(f"JSON files found: {json_files}")
            
            # Test the get_files function
            found_files = get_files(
                pattern="*.json", filepath=self.mainfile, stripname=self.basename
            )
            logger.info(f"get_files result for *.json: {found_files}")
            
            found_mrd_files = get_files(
                pattern="MRD*.json", filepath=self.mainfile, stripname=self.basename
            )
            logger.info(f"get_files result for MRD*.json: {found_mrd_files}")
            
        except Exception as e:
            logger.error(f"Error during file discovery test: {e}")

    def parse_json_file(
        self, filepath: str, logger: "BoundLogger"
    ) -> CCPNCMetadata | None:
        """Parse the JSON file and extract relevant information."""
        logger.info(
            f"Looking for JSON files with pattern 'MRD*.json' near mainfile: "
            f"{filepath}"
        )
        magres_json_file = get_files(
            pattern="MRD*.json", filepath=filepath, stripname=self.basename
        )
        if not magres_json_file:
            logger.warning("No JSON file found.")
            return None
        
        json_file_path = magres_json_file[0]
        logger.info(f"Found JSON file: {json_file_path}")

        try:
            with open(json_file_path) as f:
                magres_json_data = json.load(f)
            logger.info(f"Successfully loaded JSON data from {json_file_path}")
            logger.debug(f"JSON data keys: {list(magres_json_data.keys())}")
        except (OSError, json.JSONDecodeError) as e:
            logger.error(f"Failed to read or parse JSON file {json_file_path}: {e}")
            return None
    
        # Create metadata objects
        ccpnc_metadata = CCPNCMetadata()
        material_properties = MaterialProperties()
        orcid = ORCID()
        ccpnc_record = CCPNCRecord()
        external_database_reference = ExternalDatabaseReference()
        free_text_metadata = FreeTextMetadata()

        # Parse material properties
        material_properties.chemical_name = magres_json_data.get("chemname", "")
        material_properties.formula = magres_json_data.get("formula", "")
        material_properties.stoichiometry = magres_json_data.get("stochiometry", "")
        material_properties.elements_ratios = magres_json_data.get(
            "elements_ratios", ""
            )
        logger.debug(f"Extracted chemical_name: {material_properties.chemical_name}")

        # material_properties.chemical_name_tokens =
        orcid.orcid_id = magres_json_data.get("ORCID", "")
        logger.debug(f"Extracted ORCID: {orcid.orcid_id}")

        # ccpnc_record.visible =
        ccpnc_record.immutable_id = magres_json_data.get("immutable_id", "")
        logger.debug(f"Extracted immutable_id: {ccpnc_record.immutable_id}")

        # Parse version metadata
        version_metadata = magres_json_data.get("version_metadata", {})
        external_database_reference.external_database_name = version_metadata.get(
            "extref_type", ""
            )
        external_database_reference.external_database_reference_code = (
            version_metadata.get("extref_code", "")
        )
        free_text_metadata.uploader_author_notes = version_metadata.get("notes", "")
        free_text_metadata.structural_descriptor_notes = version_metadata.get(
            "chemform", ""
            )

        # Assemble the metadata
        ccpnc_metadata.material_properties = material_properties
        ccpnc_metadata.orcid = orcid
        ccpnc_metadata.ccpnc_record = ccpnc_record
        ccpnc_metadata.external_database_reference = external_database_reference
        ccpnc_metadata.free_text_metadata = free_text_metadata

        logger.info("Successfully created CCPNCMetadata object")
        return ccpnc_metadata
    def parse(
        self,
        filepath: str,
        archive: 'EntryArchive',
        logger: 'BoundLogger',
        child_archives: dict[str, 'EntryArchive'] = None,
    ) -> None:
        self.mainfile = filepath
        self.maindir = os.path.dirname(self.mainfile)
        self.basename = os.path.basename(self.mainfile)
        self.archive = archive
        
        logger.info(f'CCPNCMagresParser.parse starting for file: {filepath}')
        logger.info(f'Configuration parameter: {configuration.parameter}')

        # Initialize the magres file parser (from parent class)
        self.init_parser(logger=logger)
        self._check_units_magres(logger=logger)
        
        # Adding this temporarily for debugging
        self.test_file_discovery(logger)

        # Create workflow
        archive.workflow2 = Workflow(name='CCPNC Magres Processing')

        # Adding Simulation to data
        simulation = Simulation()
        logger.info("Created Simulation object")
        # archive.data = simulation

        # Parse magres file structure and calculation parameters (from parent class)
        calculation_params = self.magres_file_parser.get('calculation', {})
        if calculation_params.get('code', '') != 'CASTEP':
            logger.warning(
                'Non-CASTEP NMR simulations may not be fully supported.'
            )
        
        # Set program information
        simulation.program = Program(
            name=calculation_params.get('code', 'Unknown'),
            version=calculation_params.get('code_version', ''),
        )
        logger.info(
            f"Set program: {simulation.program.name} v{simulation.program.version}"
            )
        # Parse JSON file and extract metadata
        ccpnc_metadata = self.parse_json_file(filepath=self.mainfile, logger=logger)
        if ccpnc_metadata:
            simulation.ccpnc_metadata = ccpnc_metadata
            logger.info("Successfully assigned CCPNC metadata to simulation")
        else:
            logger.warning("No CCPNC metadata could be extracted")

        archive.data = simulation
        logger.info("Successfully assigned simulation to archive.data")
