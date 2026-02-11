from nomad_oasis_schema_parser_plugin.schema_packages.schema_package import (
    ORCID,
    CCPNCMetadata,
    CCPNCRecord,
    ExternalDatabaseReference,
    FreeTextMetadata,
    MaterialProperties,
    PublicationRecord,
)


def create_ccpnc_metadata_from_dict(metadata_dict, calculation_params, logger):
    def get_value_or_none(data, key, default=None):
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
    material_properties.chemical_name = get_value_or_none(metadata_dict, 'chemname')
    material_properties.formula = get_value_or_none(metadata_dict, 'formula')
    material_properties.stoichiometry = get_value_or_none(metadata_dict, 'stochiometry')
    material_properties.elements_ratios = get_value_or_none(
        metadata_dict, 'elements_ratios'
    )

    # Parse ORCID
    orcid.orcid_id = get_value_or_none(metadata_dict, 'ORCID')

    # Parse CCPNC record
    ccpnc_record.immutable_id = get_value_or_none(metadata_dict, 'immutable_id')

    # Parse version metadata
    version_metadata = metadata_dict.get('version_metadata', {})
    ccpnc_record.license = get_value_or_none(version_metadata, 'license')
    external_database_reference.external_database_name = get_value_or_none(
        version_metadata, 'extref_type'
    )
    external_database_reference.external_database_reference_code = get_value_or_none(
        version_metadata, 'extref_code'
    )
    free_text_metadata.uploader_author_notes = get_value_or_none(
        version_metadata, 'notes'
    )
    free_text_metadata.structural_descriptor_notes = get_value_or_none(
        version_metadata, 'chemform'
    )

    # Parse publication record
    publication_record.doi = get_value_or_none(version_metadata, 'doi')

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

    logger.info('Successfully created CCPNCMetadata object')
    return ccpnc_metadata
