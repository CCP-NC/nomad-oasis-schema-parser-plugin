import os
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from nomad.datamodel.data import ArchiveSection
    from nomad.datamodel.datamodel import EntryArchive

from glob import glob


def get_files(pattern: str, filepath: str, stripname: str = '', deep: bool = True):
    """Get files following the `pattern` with respect to the file `stripname` (usually
    this being the mainfile of the given parser) up to / down from the `filepath` 
    (`deep=True` going down, `deep=False` up)

    Args:
        pattern (str): targeted pattern to be found
        filepath (str): filepath to start the search
        stripname (str, optional): name with respect to which do the search. Defaults 
        to ''.
        deep (bool, optional): boolean setting the path in the folders to scan (down or 
        up). Defaults to down=True.

    Returns:
        list: List of found files.
    """
    for _ in range(10):
        filenames = glob(f'{os.path.dirname(filepath)}/{pattern}')
        pattern = os.path.join('**' if deep else '..', pattern)
        if filenames:
            break

    if len(filenames) > 1:
        # filter files that match
        suffix = os.path.basename(filepath).strip(stripname)
        matches = [f for f in filenames if suffix in f]
        filenames = matches if matches else filenames

    filenames = [f for f in filenames if os.access(f, os.F_OK)]
    return filenames

def get_reference(upload_id: str, entry_id: str) -> str:
    """Create a reference path to another entry in the same upload."""
    return f'../uploads/{upload_id}/archive/{entry_id}#data'


def get_entry_id_from_file_name(file_name: str, archive: 'EntryArchive') -> str:
    """Generate entry ID from filename using NOMAD's hash function."""
    from nomad.utils import hash
    return hash(archive.metadata.upload_id, file_name)

def create_archive(
    entity: 'ArchiveSection',
    archive: 'EntryArchive',
    file_name: str,
    overwrite: bool = False,
) -> str:
    """
    Create a child archive file (.archive.json) that will be processed as a separate 
    entry.

    Args:
        entity: ArchiveSection instance to be written as the data section
        archive: The parent archive
        file_name: Name for the .archive.json file (e.g., 'metadata.archive.json')
        overwrite: Whether to overwrite existing file

    Returns:
        Reference path to the created archive
    """
    if overwrite or not archive.m_context.raw_path_exists(file_name):
        with archive.m_context.update_entry(
            file_name, write=True, process=True
        ) as entry:
            entry['data'] = entity.m_to_dict(with_root_def=True)

    return get_reference(
        archive.metadata.upload_id, 
        get_entry_id_from_file_name(file_name, archive)
    )
