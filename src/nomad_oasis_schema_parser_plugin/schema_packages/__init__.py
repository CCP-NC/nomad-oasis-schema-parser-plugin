from nomad.config.models.plugins import SchemaPackageEntryPoint
from pydantic import Field


class CCPNCSchemaPackageEntryPoint(SchemaPackageEntryPoint):
    parameter: int = Field(0, description='Custom configuration parameter')

    def load(self):
        from nomad_oasis_schema_parser_plugin.schema_packages.schema_package import (
            m_package,
        )

        return m_package


ccpnc_schema_entry_point = CCPNCSchemaPackageEntryPoint(
    name='CCPNCSchemaPackage',
    description='CCPNC schema package entry point configuration.',
)
