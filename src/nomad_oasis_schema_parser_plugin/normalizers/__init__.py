from nomad.config.models.plugins import NormalizerEntryPoint


class CCPNCNormalizerEntryPoint(NormalizerEntryPoint):

    def load(self):
        from nomad_oasis_schema_parser_plugin.normalizers.normalizer import (
            CCPNCNormalizer,
        )

        return CCPNCNormalizer(**self.dict())


ccpnc_normalizer_entry_point = CCPNCNormalizerEntryPoint(
    name = ' CCPNC Normalizer',
    description = 'CCPNC custom normalizer.',
)