import re

from nomad_oasis_schema_parser_plugin.schema_packages.schema_package import (
    ORCID,
    CCPNCMetadata,
    CCPNCRecord,
    ExternalDatabaseReference,
    FreeTextMetadata,
    MaterialProperties,
    PublicationRecord,
)


def get_value_or_none(data, key, default=None):
    """Extract value from dict, returning None for empty strings."""
    value = data.get(key, default)
    if isinstance(value, str) and value.strip() == '':
        return None
    return value


def tokenize_name(name):
    """Tokenize chemical name with position-aware wildcards for search.
    
    Args:
        name: Chemical name string to tokenize
        
    Returns:
        List of tokens including original name and wildcard variants
        based on token position (start/middle/end)
    """
    sep = re.compile(r'[0-9,\'\(\)\[\]\s\-]+')
    tokens = [tk.lower() for tk in sep.split(name) if len(tk) > 1]
    result_tokens = set()
    name_lower = name.lower()
    # Always add the original name (lowercased)
    result_tokens.add(name_lower)
    if not tokens:
        return list(result_tokens)
    for tk in tokens:
        # Find all positions of token in the original name
        for match in re.finditer(re.escape(tk), name_lower):
            start, end = match.start(), match.end()
            # Check if token is at the very start of the string
            at_start = start == 0
            # Check if token is at the very end of the string
            at_end = end == len(name_lower)
            
            # Always add the base token
            result_tokens.add(tk)
            
            # Add appropriate wildcard variants based on position
            if at_start and at_end:
                # Standalone token, no wildcards needed
                pass
            elif at_start:
                # Token at start, add token* (something after)
                result_tokens.add(f'{tk}*')
            elif at_end:
                # Token at end, add *token (something before)
                result_tokens.add(f'*{tk}')
            else:
                # Token in the middle, add *token* (something before and after)
                result_tokens.add(f'*{tk}*')
    return list(result_tokens)


def create_ccpnc_metadata_from_dict(metadata_dict, calculation_params, logger):

    ccpnc_metadata = CCPNCMetadata()
    material_properties = MaterialProperties()
    orcid = ORCID()
    ccpnc_record = CCPNCRecord()
    external_database_reference = ExternalDatabaseReference()
    free_text_metadata = FreeTextMetadata()
    publication_record = PublicationRecord()

    # Parse material properties
    material_properties.chemical_name = get_value_or_none(metadata_dict, 'chemname')
    # Tokenize chemical name for search
    chemname_val = material_properties.chemical_name
    if chemname_val:
        material_properties.chemical_name_tokens = tokenize_name(chemname_val)
    else:
        material_properties.chemical_name_tokens = []
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
