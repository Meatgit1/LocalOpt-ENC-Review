"""Read and bind this run's files. Changed files are checked at their next boundary."""

from pathlib import Path
import hashlib
import json
import os
import re
import uuid


def require(value, message):
    if not value:
        raise ValueError(message)


def sha(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def save(path, value):
    with Path(path).open("x", encoding="utf-8") as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write("\n")


def _pairs(items):
    result = {}
    for key, value in items:
        require(key not in result, "Duplicate JSON field: " + key)
        result[key] = value
    return result


def load(path):
    return json.loads(
        Path(path).read_text(encoding="utf-8-sig"),
        object_pairs_hook=_pairs,
        parse_constant=lambda value: (_ for _ in ()).throw(ValueError("Nonfinite JSON")),
    )


def atomic_new(path, value):
    """One publication, with a same-directory temporary file and no replacement."""
    path = Path(path)
    require(not path.exists(), "Preserve existing publication: " + str(path))
    temp = path.with_name(path.name + "." + uuid.uuid4().hex + ".tmp")
    with temp.open("x", encoding="utf-8") as stream:
        json.dump(value, stream, separators=(",", ":"), allow_nan=False)
        stream.flush()
        os.fsync(stream.fileno())
    # Request/session claims ensure a single writer for each final name.
    require(not path.exists(), "Publication appeared while writing")
    temp.rename(path)


def input_file(root, relative):
    root = Path(root)
    require(
        isinstance(relative, str)
        and re.fullmatch(r"[A-Za-z0-9_./-]+", relative)
        and not Path(relative).is_absolute()
        and ".." not in Path(relative).parts,
        "Invalid relative path",
    )
    path = root / relative
    require(
        path.is_file() and path.resolve().is_relative_to(root.resolve()),
        "Missing or escaped file: " + relative,
    )
    require(
        not root.is_symlink()
        and all(
            not p.is_symlink() for p in (path, *path.parents) if p == root or root in p.parents
        ),
        "Linked payload path",
    )
    return path


def hash_bytes(data):
    return hashlib.sha256(data).hexdigest()


def fingerprint(paths):
    result = {}
    for path in paths:
        path = Path(path).resolve()
        require(path.is_file(), "Evidence input is missing: " + str(path))
        result[str(path)] = sha(path)
    return dict(sorted(result.items()))


def seal(directory, inputs, output_names):
    directory = Path(directory)
    require(inputs == fingerprint(inputs), "Evidence inputs changed while computing the check")
    outputs = fingerprint(directory / name for name in output_names)
    save(directory / "binding.json", {"schema": 1, "inputs": inputs, "outputs": outputs})


def verify(directory, inputs, output_names):
    binding = json.loads((Path(directory) / "binding.json").read_text())
    require(
        binding.get("schema") == 1 and binding.get("inputs") == inputs,
        "Cached check is stale: inputs, implementation or configuration changed",
    )
    require(inputs == fingerprint(inputs), "Bound inputs changed during verification")
    expected_paths = {str((Path(directory) / name).resolve()) for name in output_names}
    require(
        expected_paths and set(binding.get("outputs", {})) == expected_paths,
        "Cached check has missing or unexpected sealed outputs",
    )
    require(binding["outputs"] == fingerprint(binding["outputs"]), "Cached output changed")
    return binding
