from typing import TYPE_CHECKING

from nomad.datamodel import EntryArchive
from nomad.datamodel.results import DFT as OldModelDFT
from nomad.datamodel.results import Simulation as OldModelSimulation
from nomad.normalizing import Normalizer
from nomad.normalizing.results import ResultsNormalizer

if TYPE_CHECKING:
    pass

from nomad_oasis_schema_parser_plugin.schema_packages.eln_metadata import (
    CCPNCMetadataELN,
)
from nomad_oasis_schema_parser_plugin.schema_packages.schema_package import (
    ORCID,
    CCPNCMetadata,
    CCPNCSimulation,
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
            # Normalize sample using the results normalizer
            results_normalizer.normalize_sample(archive._ccpnc_sample)

            # Normalize measurement using the results normalizer
            results_normalizer.normalize_measurement(archive._ccpnc_measurement)

            # Populate the simulation.dft section if we have calculation params
            if hasattr(archive, '_ccpnc_calculation_params'):
                self._populate_simulation_dft(
                    archive, archive._ccpnc_calculation_params
                )

        # === METADATA SYNCHRONIZATION ===
        # Only for magres entries (CCPNCSimulation)
        if isinstance(archive.data, CCPNCSimulation):
            self.logger.warning("🔵 Normalizer: Processing CCPNCSimulation")
            self._synchronize_metadata_from_eln(archive)

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
        Synchronize ORCID metadata from metadata.archive.json.
        
        This runs when the user clicks 'Reprocess' on the magres entry.
        """
        self.logger.warning("="*60)
        self.logger.warning("📂 METADATA SYNCHRONIZATION (ORCID only)")
        self.logger.warning("="*60)
        
        try:
            from nomad.datamodel.context import ServerContext
            if not isinstance(archive.m_context, ServerContext):
                self.logger.warning("⚠️ Not in ServerContext")
                return
            
            metadata_file = "metadata.archive.json"
            
            if not archive.m_context.raw_path_exists(metadata_file):
                self.logger.warning(f"⚠️ {metadata_file} not found")
                self.logger.warning("="*60)
                return
            
            self.logger.warning(f"✅ Found {metadata_file}")
            
            # Read the file
            with archive.m_context.raw_file(metadata_file, "r") as f:
                import json
                metadata_dict = json.load(f)
            
            if 'data' not in metadata_dict:
                self.logger.warning("⚠️ No 'data' section")
                self.logger.warning("="*60)
                return
            
            data_dict = metadata_dict['data']
            
            # Verify schema
            if 'm_def' in data_dict and 'CCPNCMetadataELN' not in data_dict['m_def']:
                self.logger.warning("⚠️ Wrong schema type")
                self.logger.warning("="*60)
                return
            
            # === SYNC ORCID ===
            self.logger.warning("🔄 Syncing ORCID to simulation...")
            
            simulation = archive.data
            
            # Initialize ccpnc_metadata if needed
            if not simulation.ccpnc_metadata:
                simulation.ccpnc_metadata = CCPNCMetadata()
                self.logger.warning("   Created CCPNCMetadata section")
            
            # Extract ORCID from nested structure
            if 'orcid' in data_dict and isinstance(data_dict['orcid'], dict):
                orcid_data = data_dict['orcid']
                if 'orcid_id' in orcid_data and orcid_data['orcid_id']:
                    # Initialize ORCID if needed
                    if not simulation.ccpnc_metadata.orcid:
                        simulation.ccpnc_metadata.orcid = ORCID()
                    
                    old_value = simulation.ccpnc_metadata.orcid.orcid_id
                    new_value = orcid_data['orcid_id']
                    
                    if old_value != new_value:
                        simulation.ccpnc_metadata.orcid.orcid_id = new_value
                        self.logger.warning(f"   ✅ ORCID: {old_value or '(empty)'} → {new_value}")
                    else:
                        self.logger.warning(f"   ℹ️ ORCID unchanged: {old_value}")
                else:
                    self.logger.warning("   ⚠️ ORCID data found but orcid_id is empty")
            else:
                self.logger.warning("   ⚠️ No ORCID data in metadata file")
            
            self.logger.warning("="*60)
            
        except Exception as e:
            self.logger.warning(f"❌ Error: {e}")
            import traceback
            self.logger.warning(traceback.format_exc())
            self.logger.warning("="*60)
