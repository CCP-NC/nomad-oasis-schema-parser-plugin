import csv
import json
import os
from typing import (
    TYPE_CHECKING,
)

import ase.data

if TYPE_CHECKING:
    from nomad.datamodel.datamodel import (
        EntryArchive,
    )
    from structlog.stdlib import (
        BoundLogger,
    )

import numpy as np
from nomad.config import config
from nomad.datamodel import EntryArchive
from nomad.datamodel.metainfo.workflow import Workflow
from nomad.units import ureg
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
from nomad_nmr_schema.schema_packages.tensor_utils import NMRTensor, TensorConvention

# Import the original magres parser and NMR schema components
from nomad_parser_magres.parsers.parser import MagresParser

# from nomad_parser_magres.parsers.parser import MagresParser
from nomad_simulations.schema_packages.atoms_state import AtomsState
from nomad_simulations.schema_packages.general import Program
from runschema.run import Program as RunSchemaProgram
from runschema.run import Run as OldRun
from runschema.system import Atoms as RunSchemaAtoms
from runschema.system import System as RunSchemaSystem

# utility function used to get auxiliary files next to the `mainfile`
from nomad_oasis_schema_parser_plugin.parsers.utils.utils import (
    create_archive,
    get_files,
)
from nomad_oasis_schema_parser_plugin.schema_packages.eln_metadata import (
    CCPNCMetadataELN,
)
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

    def parse_csv_metadata(
        self,
        filepath: str,
        target_filename: str,
        logger: "BoundLogger"
    ) -> dict | None:
        """Parse CSV file to extract metadata for a specific magres file.
    
        Args:
            filepath: Path to the magres file (used to locate CSV)
            target_filename: The filename to look for in CSV (e.g., 'ethanol.magres')
            logger: Logger instance
        
        Returns:
            Dictionary with metadata matching JSON structure, or None if not found
        """
        # Look for CSV file
        csv_files = get_files(
            pattern='metadata_info.csv',
            filepath=filepath,
            stripname=self.basename,
            deep=True
        )

        if not csv_files:
            logger.info("No CSV metadata file found")
            return None

        csv_file_path = csv_files[0]

        try:
            with open(csv_file_path, encoding='utf-8') as f:
                csv_reader = csv.DictReader(f)
                
                # Find the row matching magres filename
                for row in csv_reader:
                    if row.get('filename', '').strip() == target_filename:
                        # Helper function to convert empty strings to None
                        def clean_value(value):
                            """Convert empty/whitespace strings to None"""
                            if value is None:
                                return None
                            cleaned = value.strip()
                            return None if cleaned == '' else cleaned
                        
                        # Construct metadata dict matching JSON structure
                        metadata_dict = {
                            'chemname': clean_value(row.get('chemname')),
                            'type': 'magres',
                            'version_metadata': {
                                'license': clean_value(row.get('license')),
                                'doi': clean_value(row.get('doi')),
                                'extref_type': clean_value(row.get('extref_type')),
                                'extref_code': clean_value(row.get('extref_code')),
                                'extref_other': clean_value(row.get('extref_other')),
                                'chemform': clean_value(row.get('chemform')),
                                'notes': clean_value(row.get('notes')),
                            }
                        }
                        
                        return metadata_dict
                
                logger.warning(f"No entry found for {target_filename} in CSV file")
                return None
                
        except Exception as e:
            logger.error(f"Failed to read or parse CSV file {csv_file_path}: {e}")
            return None

    def populate_metadata_from_dict(
        self,
        metadata_dict: dict,
        calculation_params: dict,
        logger: "BoundLogger"
    ) -> CCPNCMetadata:
        """Populate CCPNCMetadata from a dictionary (from JSON or CSV).
        
        Args:
            metadata_dict: Dictionary containing metadata
            calculation_params: Calculation parameters from magres file
            logger: Logger instance
            
        Returns:
            CCPNCMetadata object
        """
        def get_value_or_none(data, key, default=None):
            """Get value from dict, converting empty strings to None"""
            value = data.get(key, default)
            if isinstance(value, str) and value.strip() == '':
                return None
            return value
        
        ccpnc_metadata = CCPNCMetadata()
        material_properties = MaterialProperties()
        orcid = ORCID()
        ccpnc_record = CCPNCRecord()
        external_database_reference = ExternalDatabaseReference()
        free_text_metadata = FreeTextMetadata()
        publication_record = PublicationRecord()
        
        # Parse material properties
        material_properties.chemical_name = get_value_or_none(metadata_dict, "chemname")
        material_properties.formula = get_value_or_none(metadata_dict, "formula")
        material_properties.stoichiometry = get_value_or_none(
            metadata_dict, "stochiometry"
        )
        material_properties.elements_ratios = get_value_or_none(
            metadata_dict, "elements_ratios"
        )

        # Parse ORCID
        orcid.orcid_id = get_value_or_none(metadata_dict, "ORCID")

        # Parse CCPNC record
        ccpnc_record.immutable_id = get_value_or_none(metadata_dict, "immutable_id")
        
        # Parse version metadata
        version_metadata = metadata_dict.get("version_metadata", {})
        ccpnc_record.license = get_value_or_none(version_metadata, "license")
        external_database_reference.external_database_name = get_value_or_none(
            version_metadata, "extref_type"
        )
        external_database_reference.external_database_reference_code = (
            get_value_or_none(version_metadata, "extref_code")
        )
        free_text_metadata.uploader_author_notes = get_value_or_none(
            version_metadata, "notes"
        )
        free_text_metadata.structural_descriptor_notes = get_value_or_none(
            version_metadata, "chemform"
        )
        
        # Parse publication record
        publication_record.doi = get_value_or_none(
            version_metadata, "doi"
        )
        
        # Add magres_calc from calculation_params if not in metadata_dict
        if 'magres_calc' not in version_metadata and calculation_params:
            version_metadata['magres_calc'] = {
                'calc_code': calculation_params.get('code', ''),
                'calc_code_version': calculation_params.get('code_version', ''),
                'calc_xcfunctional': calculation_params.get('functional', ''),
            }
        
        # Assemble the metadata
        ccpnc_metadata.material_properties = material_properties
        ccpnc_metadata.orcid = orcid
        ccpnc_metadata.ccpnc_record = ccpnc_record
        ccpnc_metadata.external_database_reference = external_database_reference
        ccpnc_metadata.free_text_metadata = free_text_metadata
        ccpnc_metadata.publication_record = publication_record
        
        logger.info("Successfully created CCPNCMetadata object")
        return ccpnc_metadata

    def create_metadata_eln(
        self,
        archive: 'EntryArchive',
        logger: "BoundLogger",
    ) -> None:
        """
        Create a separate metadata.archive.json ELN entry if no metadata files exist.
        """
        try:
            # Create the ELN with initialized CCPNC metadata subsection
            eln_entry = CCPNCMetadataELN()

        except Exception as e:
            logger.error(f"Metadata object creation failed: {e}")
            import traceback
            logger.error(traceback.format_exc())
            return None

        file_name = 'metadata.archive.json'

        try:
            reference = create_archive(
                entity=eln_entry,
                archive=archive,
                file_name=file_name,
                overwrite=False,
            )
            return reference
        except Exception as e:
            logger.error(f"Failed to create metadata ELN: {e}")
            import traceback
            logger.error(traceback.format_exc())
            return None

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
            logger.warning(
                "No JSON file found. If you have added a JSON file, please name it the "
                "same as its magres file."
            )
            return None
        
        json_file_path = magres_json_file[0]
        logger.info(f"Found JSON file: {json_file_path}")

        try:
            with open(json_file_path) as f:
                magres_json_data = json.load(f)

            # Use the common populate method
            return self.populate_metadata_from_dict(
                metadata_dict=magres_json_data,
                calculation_params=None,  # JSON already has all data
                logger=logger
            )

        except (OSError, json.JSONDecodeError) as e:
            logger.error(f"Failed to read or parse JSON file {json_file_path}: {e}")
            return None

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
            logger.error('Could not find [magres] data block in magres file.')
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
            logger.error(
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

    def parse_system_old(
        self,
        logger: 'BoundLogger',
        sec_run: OldRun):
        """
        Testing old run for preparing atoms for crystal structure viewing
        """
        sec_atoms = RunSchemaAtoms()
        sec_system = RunSchemaSystem()

        atoms_old = self.magres_file_parser.get('atoms', [])
        if not atoms_old:
            logger.error("Parse error - No atoms found in atoms object")
            return None

        # Store lattice_vectors and periodic boundary conditions
        lattice_vectors_old = np.reshape(np.array(atoms_old.get('lattice', [])), (3, 3))
        sec_atoms.lattice_vectors = lattice_vectors_old * ureg.angstrom
        pbc = (
            [True, True, True] 
            if lattice_vectors_old is not None 
            else [False, False, False]
        )
        sec_atoms.periodic = pbc

        # Storing atom positions and labels
        atoms_list = atoms_old.get('atom', [])
        if len(atoms_list) == 0:
            logger.error("No atom lists found in atoms object")
            return None
        atom_labels = []
        atom_positions = []
        for atom in atoms_list:
            atom_labels.append(atom[0])
            atom_positions.append(atom[3:])  # Ensure only x,y,z are taken
        sec_atoms.labels = atom_labels
        sec_atoms.positions = atom_positions * ureg.angstrom

        # Add species (atomic numbers) based on labels
        try:
            sec_atoms.species = [
                ase.data.atomic_numbers.get(label, 0) for label in atom_labels]
        except Exception as e:
            logger.error(f"Failed to assign species (atomic numbers) to atoms: {e}")

        sec_system.atoms = sec_atoms
        sec_run.system.append(sec_system)

    def _prepare_sample_and_sec_run(
        self,
        simulation,
        archive,
        logger,
        program_name,
        program_version):
        # Create a sample class to populate results using normalizer
        class CCPNCSample:
            def __init__(self):
                self.elements = []
                self.chemical_formula = None
                self.name = None

        # Create ccpnc sample object without pre-populating elements
        ccpnc_sample = CCPNCSample()

        # Store sample data for normalizer
        archive._ccpnc_sample = ccpnc_sample

        # Prepare sec_run for system parsing
        sec_run = OldRun()
        self.parse_system_old(logger=logger, sec_run=sec_run)
        sec_run.program = RunSchemaProgram(
            name=program_name,
            version=program_version,
        )

        archive._ccpnc_sec_run = sec_run
        archive.run.append(sec_run)

    def _prepare_measurement_and_params(self, archive, calculation_params):
        class CCPNCMeasurement:
            def __init__(self):
                self.method_abbreviation = 'NMR'
                self.sample = []

            def m_xpath(self, path):
                """Mock m_xpath method that returns None for any path"""
                return None

        # Create ccpnc measurement object
        ccpnc_measurement = CCPNCMeasurement()

        # Store measurement data and calculation params for normalizer
        archive._ccpnc_measurement = ccpnc_measurement
        archive._ccpnc_calculation_params = calculation_params

    def _parse_and_attach_metadata(
        self,
        simulation,
        calculation_params,
        archive,
        logger):
        ccpnc_metadata = None
        metadata_source = None

        json_metadata = self.parse_json_file(filepath=self.mainfile, logger=logger)
        if json_metadata:
            ccpnc_metadata = json_metadata
            metadata_source = 'json'
        else:
            metadata_dict = self.parse_csv_metadata(
                filepath=self.mainfile,
                target_filename=self.basename,
                logger=logger
            )
            if metadata_dict:
                ccpnc_metadata = self.populate_metadata_from_dict(
                    metadata_dict=metadata_dict,
                    calculation_params=calculation_params,
                    logger=logger
                )
                metadata_source = 'csv'

        if ccpnc_metadata:
            simulation.ccpnc_metadata = ccpnc_metadata

        if metadata_source is None:
            metadata_reference = self.create_metadata_eln(
                archive=archive, logger=logger
            )
            if metadata_reference:
                simulation.metadata_eln_reference = metadata_reference
    
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
        code = calculation_params.get('code', '')

        # If code is a list as in QE-GIPAW cases, join it into a string and update dict
        if isinstance(code, list):
            code = ' '.join(str(c) for c in code)
            calculation_params['code'] = code  # Update the dict with the fixed value

        # Validate supported codes
        supported_codes = ['CASTEP', 'QE']
        is_supported = any(supported in code for supported in supported_codes)

        if not is_supported:
            logger.error(
                'Only CASTEP and QE-GIPAW based NMR simulations are currently supported'
                'by the CCPNC magres parser. Found calc_code: "%s"', code
            )
            return

        # Add XC functional mappings to calculation_params for the normalizer
        calculation_params['_xc_functional_type_map'] = self._xc_functional_type_map
        calculation_params['_xc_functional_map'] = self._xc_functional_map

        # Parse program information
        # Note: Older QE-GIPAW generated magres files may have limited metadata in the
        # calculation block, (e.g., calc_code_version='git'). The parser attempts to 
        # extract version from calc_code field when necessary.
        program_name, program_version = self._parse_program_info(
            calculation_params, logger
        )
        simulation.program = Program(
            name=program_name,
            version=program_version,
        )
        # Parse model system (from parent class)
        model_system = self.parse_model_system(logger=logger)
        if model_system is not None:
            simulation.model_system.append(model_system)
            self._prepare_sample_and_sec_run(
                simulation, archive, logger, program_name, program_version)   
        else:
            logger.error("Could not parse model system from magres file")

        # Parse model method (from parent class)
        model_method = self.parse_model_method(calculation_params=calculation_params)
        simulation.model_method.append(model_method)
        self._prepare_measurement_and_params(archive, calculation_params)

        # Parse NMR outputs with magnetic shielding
        outputs = self.parse_outputs_with_nmr_schema(
            simulation=simulation,
            logger=logger,
        )
        if outputs is not None:
            simulation.outputs.append(outputs)
        else:
            logger.error("Could not parse NMR outputs")

        self._parse_and_attach_metadata(simulation, calculation_params, archive, logger)

        archive.data = simulation
