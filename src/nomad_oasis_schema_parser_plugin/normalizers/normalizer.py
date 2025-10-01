from typing import TYPE_CHECKING

from nomad.datamodel import EntryArchive
from nomad.datamodel.results import DFT as OldModelDFT
from nomad.datamodel.results import Simulation as OldModelSimulation
from nomad.normalizing import Normalizer
from nomad.normalizing.results import ResultsNormalizer

if TYPE_CHECKING:
    pass


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

        else:
            self.logger.warning("No CCPNC-specific data found for normalization")

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

