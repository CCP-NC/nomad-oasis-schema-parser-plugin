import os
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from nomad.datamodel.data import ArchiveSection
    from nomad.datamodel.datamodel import EntryArchive

from glob import glob


def get_files(
    pattern: str,
    filepath: str,
    deep: bool = True,
    search_dir: str = None,
):
    """Get files following the `pattern` with respect to the file `stripname` (usually
    this being the mainfile of the given parser) up to / down from the `filepath`
    (`deep=True` going down, `deep=False` up)

    Args:
        pattern (str): targeted pattern to be found
        filepath (str): filepath to start the search (ignored when search_dir is given)
        deep (bool, optional): boolean setting the path in the folders to scan (down or
        up). Defaults to down=True.
        search_dir (str, optional): explicit root directory to start the search from,
        overriding os.path.dirname(filepath). Pass the upload raw root to ensure the
        search covers the whole upload regardless of mainfile nesting depth.

    Returns:
        list: List of found files.
    """
    base_dir = search_dir if search_dir is not None else os.path.dirname(filepath)
    for _ in range(10):
        filenames = glob(f'{base_dir}/{pattern}', recursive=True)
        if filenames:
            break
        pattern = os.path.join('**' if deep else '..', pattern)

    if len(filenames) > 1:
        raise ValueError(
            f'Ambiguous upload structure: found {len(filenames)} files matching '
            f'"{pattern}" when at most one is expected.\n'
            f'Found: {filenames}\n'
            f'Please ensure only one such file exists relative to the mainfile.'
        )

    filenames = [f for f in filenames if os.access(f, os.F_OK)]
    return filenames


def find_all_files(pattern: str, search_root: str) -> list:
    """Find all files matching `pattern` recursively from `search_root`.

    Unlike get_files, this collects every match and never raises on multiple results.
    Intended for census queries such as counting how many .magres files exist in an upload.

    Args:
        pattern: glob pattern, e.g. '*.magres'
        search_root: root directory to search from; should be the upload 'raw' root to
                     avoid picking up files from neighbouring uploads

    Returns:
        List of matching absolute file paths (may be empty).
    """
    return glob(f'{search_root}/**/{pattern}', recursive=True)


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
        archive.metadata.upload_id, get_entry_id_from_file_name(file_name, archive)
    )
