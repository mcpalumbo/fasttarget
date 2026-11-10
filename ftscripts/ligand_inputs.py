"""Standalone ligand-source configuration and local database resolution."""
import json
import math
import re
from pathlib import Path


def validate_options(source, options):
    defaults = dict(database=None, include_human=False, min_identity=90,
                    min_coverage=70, strict_errors=False)
    if source == 'chembl':
        defaults.update(direct_only=False, min_confidence_score=None,
                        require_pchembl=False, exclude_problematic=True,
                        include_species_context=True)
    unknown = set(options) - set(defaults) - {'enabled'}
    if unknown:
        raise ValueError(f'Unknown {source} options: {sorted(unknown)}')
    result = dict(defaults, **options)
    for key, default in defaults.items():
        value = result[key]
        if isinstance(default, bool) and not isinstance(value, bool):
            raise ValueError(f'{key} must be boolean')
    for key in ('min_identity', 'min_coverage'):
        value = result[key]
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or not 0 <= value <= 100:
            raise ValueError(f'{key} must be between 0 and 100')
    score = result.get('min_confidence_score')
    if score is not None and (type(score) is not int or not 0 <= score <= 9):
        raise ValueError('min_confidence_score must be an integer between 0 and 9, or null')
    database = result['database']
    if database is not None and (not isinstance(database, str) or not database.strip()):
        raise ValueError('database must be a non-empty path or null')
    return result


def local_path(value, databases_path):
    path = Path(value).expanduser()
    return path if path.is_absolute() else Path(databases_path) / path


def resolve_database(source, value, databases_path):
    if value:
        path = local_path(value, databases_path)
        if not path.is_file():
            raise FileNotFoundError(path)
        return path
    root = Path(databases_path) / source
    candidates = []
    if source == 'chembl':
        for manifest in root.glob('*/metadata.json'):
            info = json.loads(manifest.read_text())
            if info.get('status') == 'ready' and info.get('database'):
                path = manifest.parent / info['database']
                if path.is_file():
                    candidates.append(path)
    else:
        candidates = list(root.glob('*/*.sqlite'))
        if not candidates:
            candidates = list(root.glob('*/*.dmp'))
    if not candidates:
        raise FileNotFoundError(f'No local {source} database in {root}; run databases.py --download {source}')
    return max(candidates, key=lambda path: (tuple(int(n) for n in re.findall(r'\d+', str(path.relative_to(root)))), str(path)))
