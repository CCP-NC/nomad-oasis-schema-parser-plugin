from typing import TYPE_CHECKING

from nomad import atomutils
from nomad.atomutils import Formula
from nomad.datamodel import EntryArchive
from nomad.datamodel.results import DFT as OldModelDFT
from nomad.datamodel.results import Relation, System
from nomad.datamodel.results import Simulation as OldModelSimulation
from nomad.normalizing import Normalizer
from nomad.normalizing.common import (
    ase_atoms_from_nomad_atoms,
    cell_from_ase_atoms,
)
from nomad.normalizing.results import ResultsNormalizer

if TYPE_CHECKING:
    pass


from nomad_oasis_schema_parser_plugin.schema_packages.schema_package import (
    ORCID,
    CCPNCMetadata,
    CCPNCRecord,
    CCPNCSimulation,
    ExternalDatabaseReference,
    FreeTextMetadata,
    MaterialProperties,
    PublicationRecord,
)


class CCPNCNormalizer(Normalizer):
    """Normalizer for the CCPNC custom parser."""

    normalizer_level = 1

    def normalize(
        self,
        archive: EntryArchive,
        logger=None,
    ) -> None:

        if logger is not None:
            self.logger = logger.bind(normalizer=self.__class__.__name__)

        # Create normalizer instance
        results_normalizer = ResultsNormalizer()
        results_normalizer.entry_archive = archive
        results_normalizer.logger = self.logger
        
        if hasattr(archive, '_ccpnc_sample') and hasattr(archive, '_ccpnc_measurement'):
            # Normalize measurement using the results normalizer
            results_normalizer.normalize_measurement(archive._ccpnc_measurement)

            # Populate the simulation.dft section if we have calculation params
            if hasattr(archive, '_ccpnc_calculation_params'):
                self._populate_simulation_dft(
                    archive, archive._ccpnc_calculation_params
                )

            # Populate topology using the run section atoms reference
            if hasattr(archive, '_ccpnc_sec_run'):
                self._populate_topology(
                    archive, 
                    archive._ccpnc_sec_run.system[0].atoms,
                    logger=self.logger
                )

        # METADATA SYNCHRONIZATION - Only for magres entries (CCPNCSimulation)
        if isinstance(archive.data, CCPNCSimulation):
            self.logger.info("Normalizer: Processing CCPNCSimulation")
            self._synchronize_metadata_from_eln(archive)

    def _populate_topology(
        self, archive: EntryArchive, atoms_data, logger=None
    ) -> None:
        """Populate the topology section with atoms_ref and cell information."""
        
        # Check if topology already exists
        existing_topology = archive.m_xpath('results.material.topology')
        if existing_topology:
            self.logger.info("Topology already exists, skipping population")
            return
        
        try:
            # Get masses if available
            masses = atomutils.get_masses_from_computational_model(
                archive, repr_system=None
            )

            # Create the topology list
            topology: list[System] = []

            # Create ASE compatible Atoms object from atoms_data
            ase_atoms = ase_atoms_from_nomad_atoms(atoms_data)

            # Extract composition data from atoms
            elements = list(set(ase_atoms.get_chemical_symbols()))
            n_atoms = len(ase_atoms)
        
            # Create chemical formula object using NOMAD's Formula utility
            formula_obj = Formula(ase_atoms.get_chemical_formula())
            
            # Create the original/root system with composition data
            original_system = System(
                method='parser',
                label='original',
                description='A representative system from the CCPNC calculation.',
                system_relation=Relation(type='root'),
                atoms_ref=atoms_data,
                # Add composition information
                elements=elements,
                chemical_formula_hill=formula_obj.format('hill'),
                chemical_formula_reduced=formula_obj.format('reduced'),
                chemical_formula_iupac=formula_obj.format('iupac'),
                chemical_formula_descriptive=formula_obj.format('descriptive'),
                chemical_formula_anonymous=formula_obj.format('anonymous'),
                n_atoms=n_atoms,
            )
            
            # Update cell information for structure viewer
            original_system.cell = cell_from_ase_atoms(
                ase_atoms,
                masses=masses,
                atom_labels=getattr(atoms_data, 'labels', None)
            )
            index = len(topology)
            original_system.system_id = f'results/material/topology/{index}'
            topology.append(original_system)
            archive.results.material.topology = topology
            
        except Exception as e:
            self.logger.error(
                'Failed to populate topology',
                exc_info=e,
                error=str(e)
            )

    def _populate_simulation_dft(
        self, archive: EntryArchive, calculation_params: dict
    ) -> None:
        """Populate the simulation.dft section with CCPNC-specific data."""

        method = archive.results.method

        # Create simulation section (this was missing!)
        method.simulation = OldModelSimulation()

        # Get program information and populate simulation section
        program_name = calculation_params.get('code', 'Unknown')
        program_version = calculation_params.get('code_version', 'Unknown')

        # Add debugging for program version
        method.simulation.program_name = program_name
        method.simulation.program_version = program_version

        # Create DFT section
        method.simulation.dft = OldModelDFT()

        # Extract and remove the XC functional mappings
        xc_functional_type_map = calculation_params.pop('_xc_functional_type_map', {})
        xc_functional_map = calculation_params.pop('_xc_functional_map', {})

        # Extract XC functional information using the passed mappings
        xc_functional_raw = calculation_params.get('xcfunctional', 'LDA')
        method.simulation.dft.xc_functional_type = xc_functional_type_map.get(
            xc_functional_raw, 'GGA'
        )
        method.simulation.dft.xc_functional_names = xc_functional_map.get(
            xc_functional_raw, []
        )

    def _synchronize_metadata_from_eln(self, archive: EntryArchive) -> None:
        """
        Synchronize metadata from metadata.archive.json.
        Maps ELN level fields to hierarchical CCPNCMetadata structure.
    
        This runs when the user clicks 'Reprocess' on the magres entry.
        """
        # Helper function
        def update_field(parent_obj, parent_class, field_name, new_value, display_name):
            """Helper to update a field with logging"""
            if not new_value:  # Skip empty strings and None
                return parent_obj
            
            # Initialise parent if needed
            if not parent_obj:
                parent_obj = parent_class()
            
            old_value = getattr(parent_obj, field_name, None)
            if old_value != new_value:
                setattr(parent_obj, field_name, new_value)
            
            return parent_obj
        
        try:
            from nomad.datamodel.context import ServerContext
            if not isinstance(archive.m_context, ServerContext):
                self.logger.warning("Not in ServerContext")
                return
            
            metadata_file = "metadata.archive.json"
            
            # Read the metadata file
            import json
            with archive.m_context.raw_file(metadata_file, "r") as f:
                metadata_dict = json.load(f)
            
            eln_data = metadata_dict['data']
            
            # Extract ELN DATA
            simulation = archive.data
            
            # Initialize ccpnc_metadata if needed
            if not simulation.ccpnc_metadata:
                simulation.ccpnc_metadata = CCPNCMetadata()
            
            # Material Properties mapping
            simulation.ccpnc_metadata.material_properties = update_field(
                simulation.ccpnc_metadata.material_properties,
                MaterialProperties,
                'chemical_name',
                eln_data.get('chemical_name'),
                'Chemical Name'
            )
            
            # ORCID mapping
            simulation.ccpnc_metadata.orcid = update_field(
                simulation.ccpnc_metadata.orcid,
                ORCID,
                'orcid_id',
                eln_data.get('orcid_id'),
                'ORCID ID'
            )

            # License mapping
            simulation.ccpnc_metadata.ccpnc_record = update_field(
                simulation.ccpnc_metadata.ccpnc_record,
                CCPNCRecord,
                'license',
                eln_data.get('license'),
                'License'
            )

            # Publication DOI mapping
            simulation.ccpnc_metadata.publication_record = update_field(
                simulation.ccpnc_metadata.publication_record,
                PublicationRecord,
                'doi',
                eln_data.get('doi'),
                'Publication DOI'
            )

            # External Database Reference mapping
            simulation.ccpnc_metadata.external_database_reference = update_field(
                simulation.ccpnc_metadata.external_database_reference,
                ExternalDatabaseReference,
                'external_database_name',
                eln_data.get('external_database_name'),
                'External Database Name'
            )

            simulation.ccpnc_metadata.external_database_reference = update_field(
                simulation.ccpnc_metadata.external_database_reference,
                ExternalDatabaseReference,
                'external_database_name_other',
                eln_data.get('external_database_name_other'),
                'Other External Database Name'
            )

            simulation.ccpnc_metadata.external_database_reference = update_field(
                simulation.ccpnc_metadata.external_database_reference,
                ExternalDatabaseReference,
                'external_database_reference_code',
                eln_data.get('external_database_reference_code'),
                'External Database Reference Code'
            )

            # Free Text Metadata mapping
            simulation.ccpnc_metadata.free_text_metadata = update_field(
                simulation.ccpnc_metadata.free_text_metadata,
                FreeTextMetadata,
                'structural_descriptor_notes',
                eln_data.get('structural_descriptor_notes'),
                'Additonal structural descriptors'
            )

            simulation.ccpnc_metadata.free_text_metadata = update_field(
                simulation.ccpnc_metadata.free_text_metadata,
                FreeTextMetadata,
                'uploader_author_notes',
                eln_data.get('uploader_author_notes'),
                'Author\'s Notes'
            )
            
        except Exception as e:
            self.logger.warning(f"Error: {e}")
            import traceback
            self.logger.warning(traceback.format_exc())
