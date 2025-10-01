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
    PublicationRecord,
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

    def parse_json_file(
        self, filepath: str, logger: "BoundLogger"
    ) -> CCPNCMetadata | None:
        """Parse the JSON file and extract relevant information with exact 
        filename matching."""
        # Extract the base filename without extension
        base_filename = os.path.splitext(self.basename)[0]
        expected_json_filename = f"{base_filename}.json"
    
        logger.info(
            f"Looking for exact matching JSON file '{expected_json_filename}' "
            f"for magres file: {filepath}"
        )
    
        # First try exact filename matching
        magres_json_file = get_files(
            pattern=expected_json_filename, 
            filepath=filepath, 
            stripname=self.basename
        )

        if not magres_json_file:
            logger.warning("No JSON file found.")
            logger.warning(
                "If you have added json file, please name it the same as its magres "
                "file."
            )
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
        publication_record = PublicationRecord()

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
        ccpnc_record.license = version_metadata.get("license", "")
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

        # Parse publication record
        publication_record.doi = version_metadata.get("doi", "")

        # Assemble the metadata
        ccpnc_metadata.material_properties = material_properties
        ccpnc_metadata.orcid = orcid
        ccpnc_metadata.ccpnc_record = ccpnc_record
        ccpnc_metadata.external_database_reference = external_database_reference
        ccpnc_metadata.free_text_metadata = free_text_metadata
        ccpnc_metadata.publication_record = publication_record

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
        # Initial validation
        if not self._validate_simulation_data(simulation, logger):
            return None
        
        # Get validated data
        model_system = simulation.model_system[-1]
        magres_data = self.magres_file_parser.get('magres')
        if not magres_data:
            logger.warning('Could not find [magres] data block in magres file.')
            return None

        # Create outputs object with references
        outputs = NMROutputs(
            model_method_ref=(
                simulation.model_method[-1] if simulation.model_method else None
            ),
            model_system_ref=model_system,
        )
        
        # Parse all NMR quantities
        self._parse_all_nmr_quantities(outputs, magres_data, model_system, logger)
        
        return outputs

    def _validate_simulation_data(
        self, 
        simulation: Simulation, 
        logger: "BoundLogger"
    ) -> bool:
        """Validate that simulation has the required model_system data."""
        if simulation.model_system is None or len(simulation.model_system) == 0:
            logger.warning(
                'Could not find the ModelSystem that the outputs reference to.'
            )
            return False
        return True

    def _parse_all_nmr_quantities(
        self,
        outputs: NMROutputs,
        magres_data: dict,
        model_system,
        logger: "BoundLogger",
    ) -> None:
        """Parse all NMR quantities and assign them to outputs."""
        cell = model_system.cell[-1]
        atom_state_class = AtomsState

        # Create parser context
        parser_context = {
            'magres_data': magres_data,
            'cell': cell,
            'atom_state_class': atom_state_class,
            'model_system': model_system,
            'logger': logger,
        }
        
        # Define parsing configurations
        nmr_parsers = [
            {
                'method': 'parse_magnetic_shieldings',
                'output_attr': 'magnetic_shieldings',
                'log_name': 'magnetic shielding',
                'requires_filtering': False,
            },
            {
                'method': 'parse_electric_field_gradients',
                'output_attr': 'electric_field_gradients',
                'log_name': 'electric field gradient',
                'requires_filtering': False,
            },
            {
                'method': 'parse_indirect_spin_spin_couplings',
                'output_attr': 'indirect_spin_spin_couplings',
                'log_name': 'total indirect spin-spin coupling',
                'requires_filtering': True,
            },
            {
                'method': 'parse_indirect_spin_spin_couplings_fc',
                'output_attr': 'indirect_spin_spin_couplings_fermi_contact',
                'log_name': 'Fermi contact coupling',
                'requires_filtering': True,
            },
            {
                'method': 'parse_indirect_spin_spin_couplings_orbital_d',
                'output_attr': 'indirect_spin_spin_couplings_orbital_d',
                'log_name': 'orbital diamagnetic coupling',
                'requires_filtering': True,
            },
            {
                'method': 'parse_indirect_spin_spin_couplings_orbital_p',
                'output_attr': 'indirect_spin_spin_couplings_orbital_p',
                'log_name': 'orbital paramagnetic coupling',
                'requires_filtering': True,
            },
            {
                'method': 'parse_indirect_spin_spin_couplings_spin',
                'output_attr': 'indirect_spin_spin_couplings_spin_dipolar',
                'log_name': 'spin dipolar coupling',
                'requires_filtering': True,
            },
        ]
        
        # Parse standard NMR quantities
        for parser_config in nmr_parsers:
            self._parse_single_nmr_quantity(outputs, parser_context, parser_config)
        
        # Parse magnetic susceptibilities (different signature)
        self._parse_magnetic_susceptibilities(outputs, magres_data, logger)

    def _parse_single_nmr_quantity(
        self,
        outputs: NMROutputs,
        parser_context: dict,
        config: dict,
    ) -> None:
        """Parse a single NMR quantity based on configuration."""
        method = getattr(self, config['method'])
        
        # Call the parser method
        results = method(
            magres_data=parser_context['magres_data'],
            cell=parser_context['cell'],
            atom_state_class=parser_context['atom_state_class'],
            model_system=parser_context['model_system'],
            logger=parser_context['logger'],
        )
        
        # Filter results if needed
        if config['requires_filtering'] and results:
            results = [item for item in results if item is not None]
        
        # Assign results to outputs
        if results and len(results) > 0:
            setattr(outputs, config['output_attr'], results)
        else:
            parser_context['logger'].info(f"No {config['log_name']} data found")

    def _parse_magnetic_susceptibilities(
        self,
        outputs: NMROutputs,
        magres_data: dict,
        logger: "BoundLogger",
    ) -> None:
        """Parse magnetic susceptibilities (different method signature)."""
        mag_sus = self.parse_magnetic_susceptibilities(
            magres_data=magres_data, 
            logger=logger
        )
        
        if len(mag_sus) > 0:
            outputs.magnetic_susceptibilities = mag_sus
        else:
            logger.info("No magnetic susceptibility data found")

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
  
        # Initialize the magres file parser (from parent class)
        self.init_parser(logger=logger)
        self._check_units_magres(logger=logger)

        # Create workflow
        archive.workflow2 = Workflow(name='CCPNC Magres Processing')

        # Adding Simulation to data
        simulation = Simulation()

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
