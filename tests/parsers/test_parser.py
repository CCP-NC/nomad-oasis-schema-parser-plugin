import logging

from nomad.datamodel import EntryArchive

from nomad_oasis_schema_parser_plugin.parsers.parser import CCPNCMagresParser


def test_parse_file():
    parser = CCPNCMagresParser()
    archive = EntryArchive()
    # parser.parse('tests/data/example.out', archive, logging.getLogger())

    # assert archive.workflow2.name == 'test'
