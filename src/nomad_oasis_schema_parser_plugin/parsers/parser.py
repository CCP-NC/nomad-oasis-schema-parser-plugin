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
from nomad_nmr_schema.schema_packages.schema_package import (
    ElectricFieldGradient,
    IndirectSpinSpinCoupling,
    IndirectSpinSpinCouplingFermiContact,
    IndirectSpinSpinCouplingOrbitalDiamagnetic,
    IndirectSpinSpinCouplingOrbitalParamagnetic,
    IndirectSpinSpinCouplingSpinDipolar,
    MagneticShielding,
    MagneticSusceptibility,
)
from nomad_nmr_schema.schema_packages.schema_package import (
    Outputs as NMROutputs,
)

# Import the original magres parser and NMR schema components
from nomad_parser_magres.parsers.parser import MagresParser

# from nomad_parser_magres.parsers.parser import MagresParser
from nomad_simulations.schema_packages.atoms_state import AtomsState
from nomad_simulations.schema_packages.general import Program

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
        # Override the schema classes to use NMR schema for magnetic shielding
        self.mag_shielding = MagneticShielding
        self.e_field_gradient_class = ElectricFieldGradient
        self.indirect_spin_spin_couplings_class = IndirectSpinSpinCoupling
        self.indirect_spin_spin_couplings_fc_class = (
            IndirectSpinSpinCouplingFermiContact
        )
        self.indirect_spin_spin_couplings_orbital_d_class = (
            IndirectSpinSpinCouplingOrbitalDiamagnetic
        )
        self.indirect_spin_spin_couplings_orbital_p_class = (
            IndirectSpinSpinCouplingOrbitalParamagnetic
        )
        self.indirect_spin_spin_couplings_spin_class = (
            IndirectSpinSpinCouplingSpinDipolar
        )
        self.mag_susceptibility_class = MagneticSusceptibility
        self.magres_outputs_class = NMROutputs

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

    def parse_outputs_with_nmr_schema(
        self,
        simulation: Simulation,
        logger: "BoundLogger",
    ) -> NMROutputs | None:
        """
        Parse the NMR outputs section using the NMR schema, focusing on magnetic 
        shielding.
        This method combines the original magres parsing with the NMR schema.
        """
        logger.info("Starting to parse NMR outputs with magnetic shielding")
        
        # Initial check on simulation.model_system
        if simulation.model_system is None or len(simulation.model_system) == 0:
            logger.warning(
                'Could not find the ModelSystem that the outputs reference to.'
            )
            return None
            
        model_system = simulation.model_system[-1]
        if not model_system.cell or len(model_system.cell) == 0:
            logger.warning('Could not find the cell sub-section.')
            return None
            
        if (
            not hasattr(model_system, 'particle_states') 
            or not model_system.particle_states
        ):
            logger.warning('Could not find the particle_states list.')
            return None

        # Create outputs object with references
        outputs = NMROutputs(
            model_method_ref=(
            simulation.model_method[-1] if simulation.model_method else None
            ),
            model_system_ref=model_system,
        )
        
        # Check if [magres][/magres] was correctly parsed
        magres_data = self.magres_file_parser.get('magres')
        if not magres_data:
            logger.warning('Could not find [magres] data block in magres file.')
            return None

        logger.info("Found magres data block, proceeding to parse ALL NMR quantities")

        cell = model_system.cell[-1]
        atom_state_class = AtomsState
        
        # Parse magnetic shieldings
        logger.info("Parsing magnetic shieldings...")
        ms = self.parse_magnetic_shieldings(
            magres_data=magres_data,
            cell=cell,
            atom_state_class=atom_state_class,
            model_system=model_system,
            logger=logger,
        )
        if len(ms) > 0:
            outputs.magnetic_shieldings = ms
            logger.info(f"Successfully parsed {len(ms)} magnetic shielding tensors")
        else:
            logger.info("No magnetic shielding data found")

        # Parse electric field gradients
        logger.info("Parsing electric field gradients...")
        efg = self.parse_electric_field_gradients(
            magres_data=magres_data,
            cell=cell,
            atom_state_class=atom_state_class,
            model_system=model_system,
            logger=logger,
        )
        if len(efg) > 0:
            outputs.electric_field_gradients = efg
            logger.info(
                f"Successfully parsed {len(efg)} electric field gradient tensors"
            )
        else:
            logger.info("No electric field gradient data found")

        # Parse indirect spin-spin couplings (total)
        logger.info("Parsing indirect spin-spin couplings (total)...")
        isc = self.parse_indirect_spin_spin_couplings(
            magres_data=magres_data,
            cell=cell,
            atom_state_class=atom_state_class,
            model_system=model_system,
            logger=logger,
        )
        if len(isc) > 0:
            # Filter out None values that might exist from the original parser
            isc_filtered = [coupling for coupling in isc if coupling is not None]
            outputs.indirect_spin_spin_couplings = isc_filtered
            logger.info(
                f"Successfully parsed {len(isc_filtered)} total indirect spin-spin "
                f"couplings"
            )
        else:
            logger.info("No total indirect spin-spin coupling data found")

        # Parse Fermi contact contribution
        logger.info("Parsing Fermi contact spin-spin couplings...")
        isc_fc = self.parse_indirect_spin_spin_couplings_fc(
            magres_data=magres_data,
            cell=cell,
            atom_state_class=atom_state_class,
            model_system=model_system,
            logger=logger,
        )
        if len(isc_fc) > 0:
            isc_fc_filtered = [coupling for coupling in isc_fc if coupling is not None]
            outputs.indirect_spin_spin_couplings_fermi_contact = isc_fc_filtered
            logger.info(
                f"Successfully parsed {len(isc_fc_filtered)} Fermi contact couplings"
            )
        else:
            logger.info("No Fermi contact coupling data found")

        # Parse orbital diamagnetic contribution
        logger.info("Parsing orbital diamagnetic spin-spin couplings...")
        isc_orbital_d = self.parse_indirect_spin_spin_couplings_orbital_d(
            magres_data=magres_data,
            cell=cell,
            atom_state_class=atom_state_class,
            model_system=model_system,
            logger=logger,
        )
        if len(isc_orbital_d) > 0:
            isc_orbital_d_filtered = [
                coupling for coupling in isc_orbital_d if coupling is not None
            ]
            outputs.indirect_spin_spin_couplings_orbital_d = isc_orbital_d_filtered
            logger.info(
                f"Successfully parsed {len(isc_orbital_d_filtered)} orbital "
                f"diamagnetic couplings"
            )
        else:
            logger.info("No orbital diamagnetic coupling data found")

        # Parse orbital paramagnetic contribution
        logger.info("Parsing orbital paramagnetic spin-spin couplings...")
        isc_orbital_p = self.parse_indirect_spin_spin_couplings_orbital_p(
            magres_data=magres_data,
            cell=cell,
            atom_state_class=atom_state_class,
            model_system=model_system,
            logger=logger,
        )
        if len(isc_orbital_p) > 0:
            isc_orbital_p_filtered = [
                coupling for coupling in isc_orbital_p if coupling is not None
            ]
            outputs.indirect_spin_spin_couplings_orbital_p = isc_orbital_p_filtered
            logger.info(
                f"Successfully parsed {len(isc_orbital_p_filtered)} orbital "
                f"paramagnetic couplings"
            )
        else:
            logger.info("No orbital paramagnetic coupling data found")

        # Parse spin dipolar contribution
        logger.info("Parsing spin dipolar spin-spin couplings...")
        isc_spin = self.parse_indirect_spin_spin_couplings_spin(
            magres_data=magres_data,
            cell=cell,
            atom_state_class=atom_state_class,
            model_system=model_system,
            logger=logger,
        )
        if len(isc_spin) > 0:
            isc_spin_filtered = [
                coupling for coupling in isc_spin if coupling is not None
            ]
            outputs.indirect_spin_spin_couplings_spin_dipolar = isc_spin_filtered
            logger.info(
                f"Successfully parsed {len(isc_spin_filtered)} spin dipolar couplings"
            )
        else:
            logger.info("No spin dipolar coupling data found")

        # Parse magnetic susceptibilities
        logger.info("Parsing magnetic susceptibilities...")
        mag_sus = self.parse_magnetic_susceptibilities(
            magres_data=magres_data, 
            logger=logger
        )
        if len(mag_sus) > 0:
            outputs.magnetic_susceptibilities = mag_sus
            logger.info(
                f"Successfully parsed {len(mag_sus)} magnetic susceptibility tensors"
            )
        else:
            logger.info("No magnetic susceptibility data found")

        logger.info("Completed parsing all available NMR quantities")

        return outputs

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

        # Parse model system (from parent class)
        model_system = self.parse_model_system(logger=logger)
        if model_system is not None:
            simulation.model_system.append(model_system)
            logger.info(
                f"Successfully parsed model system with "
                f"{len(model_system.particle_states)} atoms"
            )
        else:
            logger.warning("Could not parse model system from magres file")

        # Parse model method (from parent class)
        model_method = self.parse_model_method(calculation_params=calculation_params)
        simulation.model_method.append(model_method)
        logger.info("Successfully parsed model method")

        # Parse NMR outputs with magnetic shielding
        outputs = self.parse_outputs_with_nmr_schema(
            simulation=simulation,
            logger=logger,
        )
        if outputs is not None:
            simulation.outputs.append(outputs)
            logger.info("Successfully parsed NMR outputs")

            # Log summary of parsed quantities
            summary = []
            if hasattr(outputs, 'magnetic_shieldings') and outputs.magnetic_shieldings:
                summary.append(
                    f"{len(outputs.magnetic_shieldings)} magnetic shieldings"
                )
            if hasattr(outputs, 'electric_field_gradients') and \
                    outputs.electric_field_gradients:
                summary.append(
                    f"{len(outputs.electric_field_gradients)} electric field gradients"
                )
            if hasattr(outputs, 'indirect_spin_spin_couplings') and \
                    outputs.indirect_spin_spin_couplings:
                summary.append(
                    f"{len(outputs.indirect_spin_spin_couplings)} total spin-spin "
                    "couplings"
                )
            if hasattr(outputs, 'indirect_spin_spin_couplings_fermi_contact') and \
                    outputs.indirect_spin_spin_couplings_fermi_contact:
                summary.append(
                    f"{len(outputs.indirect_spin_spin_couplings_fermi_contact)} Fermi "
                    "contact couplings"
                )
            if hasattr(outputs, 'indirect_spin_spin_couplings_orbital_d') and \
                    outputs.indirect_spin_spin_couplings_orbital_d:
                summary.append(
                    f"{len(outputs.indirect_spin_spin_couplings_orbital_d)} orbital "
                    "diamagnetic couplings"
                )
            if hasattr(outputs, 'indirect_spin_spin_couplings_orbital_p') and \
                    outputs.indirect_spin_spin_couplings_orbital_p:
                summary.append(
                    f"{len(outputs.indirect_spin_spin_couplings_orbital_p)} orbital "
                    "paramagnetic couplings"
                )
            if hasattr(outputs, 'indirect_spin_spin_couplings_spin_dipolar') and \
                    outputs.indirect_spin_spin_couplings_spin_dipolar:
                summary.append(
                    f"{len(outputs.indirect_spin_spin_couplings_spin_dipolar)} spin "
                    "dipolar couplings"
                )
            if hasattr(outputs, 'magnetic_susceptibilities') and \
                    outputs.magnetic_susceptibilities:
                summary.append(
                    f"{len(outputs.magnetic_susceptibilities)} "
                    "magnetic susceptibilities"
                )
            
            if summary:
                logger.info(f"Parsed NMR quantities summary: {', '.join(summary)}")
        else:
            logger.warning("Could not parse NMR outputs")

        # Parse JSON file and extract metadata
        ccpnc_metadata = self.parse_json_file(filepath=self.mainfile, logger=logger)
        if ccpnc_metadata:
            simulation.ccpnc_metadata = ccpnc_metadata
            logger.info("Successfully assigned CCPNC metadata to simulation")
        else:
            logger.warning("No CCPNC metadata could be extracted")

        archive.data = simulation
        logger.info("Successfully assigned simulation to archive.data")
