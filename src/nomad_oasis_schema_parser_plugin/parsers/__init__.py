from nomad.config.models.plugins import ParserEntryPoint
from pydantic import Field


class CCPNCParserEntryPoint(ParserEntryPoint):
    parameter: int = Field(0, description='Custom configuration parameter')

    def load(self):
        from nomad_oasis_schema_parser_plugin.parsers.parser import CCPNCMagresParser

        return CCPNCMagresParser(**self.dict())


ccpnc_parser_entry_point = CCPNCParserEntryPoint(
    name='CCPNCParserEntryPoint',
    description='CCPNC parser entry point configuration.',
    level=1,
    parser_as_interface=False,  # in order to use `child_archives` and auto workflows
    mainfile_name_re=r'.*\.magres',
)
