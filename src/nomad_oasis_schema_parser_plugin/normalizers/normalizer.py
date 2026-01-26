import json
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
    ExternalDatabaseReference,
    FreeTextMetadata,
    MaterialProperties,
    PublicationRecord,
)


class CCPNCNormalizer(Normalizer):
    """Normalizer for the CCPNC custom parser."""

    normalizer_level = 3

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

            # Populate topology using the run section atoms reference
            if hasattr(archive, '_ccpnc_sec_run'):
                self._populate_topology(
                    archive, 
                    archive._ccpnc_sec_run.system[0].atoms,
                    logger=self.logger
                )

            # Populate the simulation.dft section if we have calculation params
            if hasattr(archive, '_ccpnc_calculation_params'):
                self._populate_simulation_dft(
                    archive, archive._ccpnc_calculation_params
                )
    
        if (archive.results and archive.results.method and 
            archive.results.method.simulation and 
            archive.results.method.simulation.dft):
            dft = archive.results.method.simulation.dft
            self.logger.info({
                "event": "CCPNCNormalizer FINAL DFT VALUES",
                "xc_functional_type": getattr(dft, 'xc_functional_type', 'MISSING'),
                "xc_functional_names": getattr(dft, 'xc_functional_names', 'MISSING'),
                "jacobs_ladder": getattr(dft, 'jacobs_ladder', 'MISSING'),
                "normalizer": "CCPNCNormalizer"
            })

        # List of fields to check and populate from topology[0] if missing or empty
        fields = [
            "elements",
            "chemical_formula_descriptive",
            "chemical_formula_reduced",
            "chemical_formula_hill",
            "chemical_formula_iupac",
            "chemical_formula_anonymous",]
        if (
            hasattr(archive, "results")
            and hasattr(archive.results, "material")
            and hasattr(archive.results.material, "topology")
            and archive.results.material.topology
        ):
            topo_buffer = archive.results.material.topology[0]
            for field in fields:
                value = getattr(archive.results.material, field, None)
                topo_value = getattr(topo_buffer, field, None)
                if not value and topo_value:
                    setattr(archive.results.material, field, topo_value)
        
        # Populate element-resolved magnetic shielding from normalized outputs
        self._populate_element_resolved_magnetic_shielding(archive, logger)
        # Only try to sync from ELN if no metadata exists yet
        if not (
            hasattr(archive.data, 'ccpnc_metadata') 
            and archive.data.ccpnc_metadata is not None
        ):
            logger.info("No metadata found, attempting to synchronize from ELN")
            self._synchronize_metadata_from_eln(archive, logger)
        else:
            logger.info("Metadata already present, skipping ELN synchronization")

    def _process_magnetic_shielding_entry_normalized(
        self, i, ms, particle_states_ref, logger, ElementIsotropyEntry
    ):
        """
        Process a single magnetic shielding entry using pre-computed isotropy values.
        
        Args:
            i: Index of the entry
            ms: MagneticShielding object with normalized isotropy
            particle_states_ref: List of particle states for matching entity_ref
            logger: Logger instance
            ElementIsotropyEntry: Class for creating isotropy entries
            
        Returns:
            ElementIsotropyEntry or None if entry is invalid
        """
        entity_ref = getattr(ms, 'entity_ref', None)
        isotropy = getattr(ms, 'isotropy', None)
        
        if entity_ref is None:
            logger.warning(f"Skipping MS entry {i}: missing entity_ref.")
            return None
            
        if isotropy is None:
            logger.warning(f"Skipping MS entry {i}: isotropy not computed.")
            return None
        
        atom = next((ps for ps in particle_states_ref if entity_ref is ps), None)
        if atom is None:
            logger.warning(
                f"Skipping MS entry {i}: could not match entity_ref to any "
                "particle_state."
            )
            return None
        
        chemical_symbol = getattr(atom, 'chemical_symbol', None)
        if chemical_symbol is None:
            logger.warning(f"Skipping MS entry {i}: could not resolve chemical_symbol.")
            return None
        
        entry = ElementIsotropyEntry()
        entry.element = chemical_symbol
        entry.isotropy = isotropy
        return entry

    def _group_and_set_isotropies(
        self,
        element_isotropy_list,
        ms_section,
        IsotropyEntry
    ):
        """
        Group isotropy values by element and set element-specific isotropy lists.
        
        Args:
            element_isotropy_list: List of ElementIsotropyEntry objects
            ms_section: ElementResolvedMagneticShielding section to populate
            IsotropyEntry: Class for creating isotropy entries
        """
        element_groups = {}
        for entry in element_isotropy_list:
            element = entry.element
            if element not in element_groups:
                element_groups[element] = []
            element_groups[element].append(entry.isotropy)
        
        for element, isotropies in element_groups.items():
            attr_name = f"{element}_isotropy_list"
            if hasattr(ms_section, attr_name):
                setattr(
                    ms_section,
                    attr_name,
                    [IsotropyEntry(isotropy=iso) for iso in isotropies]
                )

    def _populate_element_resolved_magnetic_shielding(
        self, 
        archive: EntryArchive, 
        logger
    ) -> None:
        """
        Populate the element_resolved_magnetic_shielding section using 
        pre-computed isotropy values from normalized MagneticShielding objects.
        """
        # Validate archive structure
        if not hasattr(archive.data, 'outputs') or len(archive.data.outputs) == 0:
            logger.info("No outputs found for element-resolved magnetic shielding.")
            return
        
        if not hasattr(archive.data, 'model_system') or len(archive.data.model_system) == 0:
            logger.info("No model_system found for element-resolved magnetic shielding.")
            return

        outputs_ref = archive.data.outputs[0]
        model_system_ref = archive.data.model_system[0]
        particle_states_ref = getattr(model_system_ref, 'particle_states', None)
        
        if not particle_states_ref:
            logger.warning("No particle_states found in model_system.")
            return

        ms_list = getattr(outputs_ref, 'magnetic_shieldings', None)
        if not ms_list or len(ms_list) == 0:
            logger.info("No magnetic_shieldings found in outputs.")
            return

        # Import schema classes
        from nomad_oasis_schema_parser_plugin.schema_packages.schema_package import (
            ElementIsotropyEntry,
            ElementResolvedMagneticShielding,
            ElementResolvedNMRSearch,
            IsotropyEntry,
        )

        # Ensure all magnetic shielding objects are normalized
        for ms in ms_list:
            if hasattr(ms, 'normalize') and not hasattr(ms, '_normalized'):
                ms.normalize(archive, logger)
                ms._normalized = True

        # Process each magnetic shielding entry
        element_isotropy_list = [
            self._process_magnetic_shielding_entry_normalized(
                i, ms, particle_states_ref, logger, ElementIsotropyEntry
            )
            for i, ms in enumerate(ms_list)
        ]
        element_isotropy_list = [
            entry for entry in element_isotropy_list if entry is not None
        ]

        if not element_isotropy_list:
            logger.warning("No valid magnetic shielding entries to populate.")
            return

        # Create element-resolved sections
        ms_section = ElementResolvedMagneticShielding()
        self._group_and_set_isotropies(element_isotropy_list, ms_section, IsotropyEntry)
        ms_section.element_isotropy_list = element_isotropy_list

        # Get or create element_resolved_nmr_search section
        if hasattr(archive.data, 'element_resolved_nmr_search') and archive.data.element_resolved_nmr_search:
            element_section = archive.data.element_resolved_nmr_search
        else:
            element_section = ElementResolvedNMRSearch()
            archive.data.element_resolved_nmr_search = element_section

        element_section.element_resolved_magnetic_shielding = ms_section
        
        logger.info(
            f"Successfully populated element-resolved magnetic shielding with "
            f"{len(element_isotropy_list)} entries."
        )
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

        # Ensure results.method exists
        if not hasattr(archive, 'results') or not archive.results:
            self.logger.warning(
                "No results section found, cannot populate simulation.dft"
            )
            return
            
        if not hasattr(archive.results, 'method') or not archive.results.method:
            self.logger.warning(
                "No method section found, cannot populate simulation.dft"
            )
            return
        
        method = archive.results.method
        # Create simulation section
        if not hasattr(method, 'simulation') or not method.simulation:
            method.simulation = OldModelSimulation()
            self.logger.info("Created new simulation section")
        else:
            self.logger.info("Using existing simulation section")

        # Get program information and populate simulation section
        program_name = 'Unknown'
        program_version = 'Unknown'
        if (hasattr(archive, 'run') and archive.run and 
            hasattr(archive.run[0], 'program')):
            run_program = archive.run[0].program
            program_name = getattr(run_program, 'name', 'Unknown')
            program_version = getattr(run_program, 'version', 'Unknown')

        # Set program name and version
        method.simulation.program_name = program_name
        method.simulation.program_version = program_version

        # Create DFT section
        if not hasattr(method.simulation, 'dft') or not method.simulation.dft:
            method.simulation.dft = OldModelDFT()
            self.logger.info("Created new DFT section")
        else:
            self.logger.info("Using existing DFT section")

        # Extract and remove the XC functional mappings
        xc_functional_type_map = calculation_params.pop('_xc_functional_type_map', {})
        xc_functional_map = calculation_params.pop('_xc_functional_map', {})

        # Extract XC functional information using the passed mappings
        xc_functional_raw = calculation_params.get('xcfunctional', 'LDA')
        # Get mapped values
        xc_functional_type = xc_functional_type_map.get(xc_functional_raw, 'GGA')
        xc_functional_names = xc_functional_map.get(xc_functional_raw, [])
        
        # Set the values - set both jacobs_ladder AND xc_functional_type
        method.simulation.dft.jacobs_ladder = xc_functional_type
        method.simulation.dft.xc_functional_type = xc_functional_type
        method.simulation.dft.xc_functional_names = xc_functional_names

    def _synchronize_metadata_from_eln(
        self, 
        archive: EntryArchive,
        logger=None
        ) -> bool:
        """
        Synchronize metadata from ELN entry if it exists and no metadata is already 
        present.
        Maps ELN level fields to hierarchical CCPNCMetadata structure.

        Returns:
            bool: True if metadata was successfully synchronized, False otherwise
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

        # Check if metadata already exists (from JSON or CSV parsing)
        if (
            hasattr(archive.data, 'ccpnc_metadata')
            and archive.data.ccpnc_metadata is not None
        ):
            logger.info(
                "Metadata already exists from JSON/CSV parsing, "
                "skipping ELN synchronization"
            )
            return True

        # Check if ELN reference exists
        if (
            not hasattr(archive.data, 'metadata_eln_reference')
            or archive.data.metadata_eln_reference is None
        ):
            logger.info("No ELN reference found, skipping metadata synchronization")
            return False
    
        metadata_file = "metadata.archive.json"
        
        try:
            with archive.m_context.raw_file(metadata_file, "r") as f:
                metadata_data = json.load(f)

            # Extract the ELN data and convert to CCPNCMetadata
            eln_data = metadata_data.get("data", {})
            # Check if there's any ELN data at all
            if not eln_data:
                self.logger.warning("No ELN data found in metadata file")
                return False

            # Check if any relevant fields exist in the ELN data
            relevant_fields = [
                'chemical_name', 'orcid_id', 'license', 'doi', 
                'external_database_name', 'external_database_name_other',
                'external_database_reference_code', 'structural_descriptor_notes',
                'uploader_author_notes'
            ]
        
            has_relevant_data = any(eln_data.get(field) for field in relevant_fields)
            if not has_relevant_data:
                self.logger.warning("No relevant CCPNC metadata fields found in ELN")
                return False
            
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

            self.logger.info("Successfully synchronized metadata from ELN")
            return True

        except KeyError:
            self.logger.error(
                "No ELN metadata file found, which is expected when metadata "
                "comes from JSON/CSV"
            )
        except Exception as e:
            self.logger.error(f"Error: {e}")
            import traceback
            self.logger.error(traceback.format_exc())
        return False