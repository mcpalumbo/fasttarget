import os
import json
import sys
import pandas as pd
import logging
import math
import tempfile
import hashlib
from datetime import datetime


def file_to_list(file_path):
    
    """
    Read a text file and return a list with each line as an element.

    :param file_path: The name of the text file to read from.

    :return: List with each line as an element.
    
    """
    
    my_list = []
    
    # Returning None on failure pushed the error downstream, where it surfaced as
    # an AttributeError or TypeError on None far from the missing file. Raise here
    # so the message names the file.
    try:
        with open(file_path, 'r') as file:
            my_list = [line.strip() for line in file]
        return my_list
    except FileNotFoundError:
        logging.error(f'File not found: {file_path}')
        raise
    except Exception as e:
        logging.exception(f"An unexpected error occurred while reading the file: {e}")
        raise

def list_to_file(file_path, my_list):

    """
    Write a list to a .txt file. Each element as a new line.

    :param file_path: Path of .txt file.
    :param my_list: List to covert.

    """

    with open(file_path, 'w') as file:
        for item in my_list:
            file.write(str(item) + '\n')

def file_check(file_path):

    """
    Check if a file exists and is not empty.

    :param file_path: Path of the file.

    :return: True if the file exists and it is not empty, False otherwise.

    """

    file_res = False

    if file_path:

        # Check if the file exists
        if not os.path.isfile(file_path):
            file_res = False
        else:
            # Check if the file is empty
            if os.path.getsize(file_path) > 0:
                file_res = True
            else:
                file_res = False

    return file_res

# Completion markers
#
# file_check answers "is there usable content here?", which is the right question
# for most resume gates. It cannot answer "did this step already run?", because an
# empty output is ambiguous: a search that legitimately found no hits and a search
# that crashed after creating the file both leave zero bytes.
#
# Steps whose empty output is a valid result write a marker once the tool exits
# successfully. A resume can then skip them without having to guess, while an
# empty output with no marker is still treated as work to redo.

def completion_marker(output_path):
    """
    Path of the completion marker for an output file.

    :param output_path: Path of the output file the step produces.
    :return: Path of the marker file.
    """
    return f'{output_path}.complete'


def mark_complete(output_path):
    """
    Record that the step producing output_path finished successfully.

    :param output_path: Path of the output file the step produces.
    """
    marker = completion_marker(output_path)
    try:
        with open(marker, 'w') as handle:
            handle.write(f'{datetime.now().isoformat()}\n')
    except OSError as error:
        # The marker is an optimization: failing to write it only costs a redo.
        logging.warning(f'Could not write completion marker {marker}: {error}')


def step_completed(output_path):
    """
    Check whether a step already ran to completion, even with empty output.

    :param output_path: Path of the output file the step produces.
    :return: True when the output has content, or when it exists alongside a
        completion marker. False when the output is missing, or empty with no
        marker (which may be the remains of a crashed run).
    """
    if file_check(output_path):
        return True
    return os.path.isfile(output_path) and os.path.isfile(completion_marker(output_path))

# Parameter manifests
#
# Resume gates ask only whether an output file exists. That is correct when the
# output depends on nothing but its inputs, but most steps are tuned by config
# values: raising core.min_identity and re-running reuses a table built at the
# old identity, with nothing in the output saying so.
#
# A step that depends on parameters writes a manifest next to its output
# recording the values actually used. The gate then checks the manifest as well
# as the file, so a changed parameter is a cache miss instead of a silent reuse.

MANIFEST_SUFFIX = '.manifest.json'


def manifest_path(output_path):
    """
    Path of the parameter manifest for an output file.

    :param output_path: Path of the output file the step produces.
    :return: Path of the manifest file.
    """
    return f'{output_path}{MANIFEST_SUFFIX}'


def _normalize_manifest_value(value):
    """
    Reduce a parameter to a form that compares equal across runs.

    Values travel through YAML, Nextflow string interpolation and JSON, so the
    same setting can arrive as 95, 95.0 or '95'. Comparing the raw values would
    report a mismatch for a parameter that did not change.

    :param value: Parameter value.
    :return: Normalized value.
    """
    if value is None:
        return None
    if isinstance(value, bool):
        # Tagged, because numbers become floats below and True == 1.0 in Python:
        # a flag turned on would compare equal to a threshold of 1.
        return f'<bool:{str(value).lower()}>'
    if isinstance(value, (int, float)):
        return float(value)
    if isinstance(value, str):
        try:
            return float(value.strip())
        except ValueError:
            return value.strip()
    if isinstance(value, (list, tuple)):
        return [_normalize_manifest_value(item) for item in value]
    if isinstance(value, dict):
        return {str(key): _normalize_manifest_value(item) for key, item in sorted(value.items())}
    return str(value)


def _normalize_parameters(parameters):
    """
    Normalize every value of a parameter mapping.

    :param parameters: Mapping of parameter name to value.
    :return: Mapping with normalized values.
    """
    return {str(key): _normalize_manifest_value(value) for key, value in (parameters or {}).items()}


def write_manifest(output_path, parameters, extra=None, hash_output=True):
    """
    Record the parameters a step used to produce its output.

    :param output_path: Path of the output file the step produced.
    :param parameters: Mapping of parameter name to the value actually used.
    :param extra: Optional mapping of additional context to store alongside.
    :param hash_output: Whether to record the output's SHA-256 and size, so a
        later run can tell a complete output from a partially written one. Pass
        False only for an output too large to be worth re-reading.
    """
    manifest = {
        'output': os.path.basename(output_path),
        'written': datetime.now().isoformat(),
        'parameters': _normalize_parameters(parameters),
    }

    # The parameters say how the output was meant to be produced; this says what
    # was actually written. Without it a write interrupted halfway left a
    # truncated file with a valid manifest, which the next run reused: a table cut
    # from twenty rows to two was a cache hit.
    if hash_output:
        try:
            manifest['output_sha256'] = sha256_file(output_path)
            manifest['output_bytes'] = os.path.getsize(output_path)
        except OSError as error:
            logging.warning(f'Could not checksum {output_path} for its manifest: {error}')

    if extra:
        manifest['extra'] = extra
    try:
        atomic_write_json(manifest, manifest_path(output_path))
    except OSError as error:
        # Like the completion marker, the manifest is an optimization: failing to
        # write it costs a recompute, it does not invalidate the output.
        logging.warning(f'Could not write manifest for {output_path}: {error}')


def read_manifest(output_path):
    """
    Read the parameter manifest of an output file.

    :param output_path: Path of the output file the step produces.
    :return: Manifest dictionary, or None when missing or unreadable.
    """
    try:
        with open(manifest_path(output_path), 'r', encoding='utf-8') as handle:
            manifest = json.load(handle)
    except (OSError, json.JSONDecodeError):
        return None
    return manifest if isinstance(manifest, dict) else None


def _display_manifest_value(value):
    """
    Render a normalized value the way it was written in the config.

    :param value: Normalized parameter value.
    :return: String for a message: the tag of a boolean and the trailing '.0' of a
        whole number are internal to the comparison and would only confuse.
    """
    if value == '<bool:true>':
        return 'True'
    if value == '<bool:false>':
        return 'False'
    if isinstance(value, float) and value.is_integer():
        return str(int(value))
    return str(value)


def manifest_differences(output_path, parameters):
    """
    Parameters an existing output was demonstrably built with other values for.

    Only recorded parameters are compared. A parameter the manifest does not
    mention is unknown rather than different: manifests written by an earlier
    version of FastTarget do not carry fields that were added later, and treating
    those as changed would stop a run whose parameters are in fact untouched.
    Use manifest_unrecorded to find them.

    :param output_path: Path of the output file the step produces.
    :param parameters: Mapping of parameter name to the value for this run.
    :return: List of '<name> <recorded> -> <current>' strings, empty when every
        recorded parameter matches.
    """
    expected = _normalize_parameters(parameters)
    recorded = _normalize_parameters((read_manifest(output_path) or {}).get('parameters'))

    differences = []
    for name in sorted(expected):
        if name in recorded and recorded[name] != expected[name]:
            differences.append(
                f'{name} {_display_manifest_value(recorded[name])} '
                f'-> {_display_manifest_value(expected[name])}'
            )
    return differences


def manifest_output_changed(output_path):
    """
    Check whether an output differs from what its manifest recorded.

    :param output_path: Path of the output file the step produces.
    :return: True when the manifest holds a checksum that no longer matches, so
        the file was truncated, overwritten or edited since. False when it
        matches, and False when no checksum was recorded: an output from an
        earlier version is unverifiable, not known to be wrong.
    """
    manifest = read_manifest(output_path) or {}
    recorded = manifest.get('output_sha256')
    if not recorded:
        return False

    size = manifest.get('output_bytes')
    try:
        # Comparing the size first avoids re-reading a large file in the common
        # case where a truncated write also changed its length.
        if size is not None and os.path.getsize(output_path) != size:
            return True
        return sha256_file(output_path) != recorded
    except OSError:
        return True


def manifest_unrecorded(output_path, parameters):
    """
    Parameters the manifest of an existing output says nothing about.

    :param output_path: Path of the output file the step produces.
    :param parameters: Mapping of parameter name to the value for this run.
    :return: Sorted list of parameter names absent from the manifest. Every name
        when there is no manifest at all.
    """
    expected = _normalize_parameters(parameters)
    recorded = _normalize_parameters((read_manifest(output_path) or {}).get('parameters'))
    return sorted(name for name in expected if name not in recorded)


def manifest_matches(output_path, parameters, label=None):
    """
    Check whether an existing output was built with the given parameters.

    :param output_path: Path of the output file the step produces.
    :param parameters: Mapping of parameter name to the value for this run.
    :param label: Optional name used in the log line, defaults to the file name.
    :return: True when the output exists and its manifest matches. False when the
        output is missing, has no manifest, or was built with other values.
    """
    if not file_check(output_path):
        return False

    # A cheap output is rebuilt for an unknown parameter as well as a changed one:
    # recomputing costs little, and it replaces a result of unknown provenance with
    # one that carries its parameters.
    reasons = manifest_differences(output_path, parameters)
    reasons += [f'{name} unrecorded' for name in manifest_unrecorded(output_path, parameters)]
    if manifest_output_changed(output_path):
        reasons.append('the output no longer matches its recorded checksum')
    if not reasons:
        return True

    # print only: stdout is captured into the log (see ftscripts.logger), so also
    # calling logging.info here would write the same line to the file twice.
    name = label or os.path.basename(output_path)
    print(f'Rebuilding {name}: {"; ".join(reasons)}')
    return False


def json_to_dict (file_path):

    """
    Read a .json file to a dict.

    :param file_path: Path of .json file.

    :return: Dictionary representing JSON object.

    """

    if os.path.exists(file_path):
        # Read dictionary from the file_path given
        with open(file_path, 'r') as file:
            loaded_dict = json.load(file)
    else:
        logging.error(f"The file '{file_path}' not found.")
        loaded_dict = None
    return loaded_dict

def jsonl_to_dict(file_path):
    """
    Read a .jsonl file to a list of dicts.

    :param file_path: Path of .jsonl file.

    :return: List of dictionaries representing JSON objects.
    """

    if os.path.exists(file_path):
        loaded_dicts = []
        with open(file_path, 'r') as file:
            for line in file:
                loaded_dicts.append(json.loads(line.strip()))
    else:
        logging.error(f"The file '{file_path}' not found.")
        loaded_dicts = None
    return loaded_dicts

def dict_to_json(output_path, file_name, my_dict):

    """
    Write a .json file from a a dictionary.

    :param output_path:  Path to save the .json file.
    :param file_name: Name of the .json file.
    :param my_dict: Dictionary to save.

    """

    def _sanitize_for_json(obj):
        if isinstance(obj, dict):
            return {k: _sanitize_for_json(v) for k, v in obj.items()}
        if isinstance(obj, list):
            return [_sanitize_for_json(v) for v in obj]
        if isinstance(obj, tuple):
            return [_sanitize_for_json(v) for v in obj]
        if hasattr(obj, "item"):
            return _sanitize_for_json(obj.item())
        if isinstance(obj, float) and (math.isnan(obj) or math.isinf(obj)):
            return None
        return obj

    full_path = os.path.join(output_path, file_name)

    # A missing directory used to print a line and return, so the data was
    # dropped and the caller went on believing it had been written. Every caller
    # writes into the organism output tree, which the pipeline owns, so creating
    # the directory is what they all want and removes the silent-loss path.
    if not os.path.isdir(output_path):
        logging.info(f"Creating '{output_path}' to save '{file_name}'.")
        os.makedirs(output_path, exist_ok=True)

    with open(full_path, 'w') as file:
        json.dump(_sanitize_for_json(my_dict), file, allow_nan=False)
    logging.debug(f"File '{full_path}' saved.")

def sha256_file(file_path):
    """
    Calculates the SHA-256 checksum of a file.

    :param file_path: Path to the input file.
    :return: Hexadecimal SHA-256 checksum.
    """
    sha256 = hashlib.sha256()
    with open(file_path, "rb") as input_file:
        for chunk in iter(lambda: input_file.read(1024 * 1024), b""):
            sha256.update(chunk)
    return sha256.hexdigest()

def fileset_fingerprint(paths):
    """
    Short digest identifying a set of input files.

    Recording the file names themselves would make a manifest unreadable for a
    few thousand genomes, and a mismatch message unusable. The digest is reported
    alongside a count, which is what makes it meaningful to a reader.

    :param paths: Iterable of file paths. Only the base names are used, so the set
        survives being moved or staged into a work directory.
    :return: First 12 hex characters of the SHA-256 of the sorted base names.
    """
    names = sorted(os.path.basename(path) for path in paths)
    digest = hashlib.sha256('\n'.join(names).encode('utf-8')).hexdigest()
    return digest[:12]


def create_organism_subfolders(output_path, organism_name):

    """
    Create subfolders for an organism. 
    
    :param output_path: Output path where the organism folder will be created.
    :param organism_names: Name of the organism.
    """
    organism_dir = os.path.join(output_path, organism_name)

    if not os.path.exists(organism_dir):
        os.makedirs(organism_dir, exist_ok=True)
        print(f'Created directory: {organism_dir}')
    
    list_dir = ['offtarget', 'metabolism', 'structures', 'essentiality', 'conservation', 'metadata', 'genome', 'localization']
    # Subfolders
    for dir in list_dir:
        dir_path = os.path.join(organism_dir, dir)
        if not os.path.exists(dir_path):
            os.makedirs(dir_path, exist_ok=True)
            print(f"Created directory: {dir_path}")

def read_blast_output(file_path, len=False):
    """
    Read BLASTP output file and return a pandas DataFrame.

    :param file_path: File path of blast output file.
    :len: If True, the length of the query and subject sequences will be added to the DataFrame.

    :return: Pandas DataFrame with blast output.
    """

    if os.path.exists(file_path):
        blast_output_df = pd.read_csv(file_path, sep='\t', header=None)

        # Columns used for blastp
        blast_columns = [
        "qseqid",   # query or source (gene) sequence id
        "sseqid",   # subject or target (reference genome) sequence id
        "pident",   # percentage of identical positions
        "length",   # alignment length (sequence overlap)
        "mismatch", # number of mismatches
        "gapopen",  # number of gap openings
        "qstart",   # start of alignment in query
        "qend",     # end of alignment in query
        "sstart",   # start of alignment in subject
        "send",     # end of alignment in subject
        "evalue",   # expect value
        "bitscore", # bit score
        "qcovhsp",  # Query Coverage hsp
        "qcovs"     # Query Coverage full
        ]

        if len: #add qlen and slen to columns
            blast_columns.extend(['qlen', 'slen'])

        blast_output_df.columns = blast_columns

        return blast_output_df

    else:
        # Previously returned None, so callers failed later on None.iterrows()
        # instead of being told which file was missing.
        raise FileNotFoundError(f'BLAST output file not found: {file_path}')


def atomic_write_dataframe_parquet(dataframe, output_path):
    """
    Writes a dataframe atomically in Parquet format.

    :param dataframe: Pandas dataframe to write.
    :param output_path: Path to the output Parquet file.
    """
    output_dir = os.path.dirname(output_path)
    temporary_path = None
    try:
        with tempfile.NamedTemporaryFile(
            dir=output_dir,
            prefix=f".{os.path.basename(output_path)}.",
            suffix=".tmp",
            delete=False,
        ) as temporary_file:
            temporary_path = temporary_file.name
        dataframe.to_parquet(temporary_path, index=False)
        os.replace(temporary_path, output_path)
    finally:
        if temporary_path and os.path.exists(temporary_path):
            os.remove(temporary_path)


def atomic_write_json(data, output_path):
    """
    Writes data atomically in JSON format.

    :param data: Data to serialize as JSON.
    :param output_path: Path to the output JSON file.
    """
    output_dir = os.path.dirname(output_path)
    temporary_path = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            dir=output_dir,
            prefix=f".{os.path.basename(output_path)}.",
            suffix=".tmp",
            delete=False,
        ) as temporary_file:
            temporary_path = temporary_file.name
            json.dump(data, temporary_file, indent=2, sort_keys=True)
            temporary_file.write("\n")
            temporary_file.flush()
            os.fsync(temporary_file.fileno())
        os.replace(temporary_path, output_path)
    finally:
        if temporary_path and os.path.exists(temporary_path):
            os.remove(temporary_path)
