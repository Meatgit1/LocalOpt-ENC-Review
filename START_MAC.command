#!/bin/sh
set -eu
bundle_root=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
bundle_python=${INTRANCE_PYTHON:-"$bundle_root/.venv-mac/bin/python"}
if [ ! -x "$bundle_python" ]; then
  echo "Create .venv-mac with Python 3.12 first. The commands are in RUN.md." >&2
  exit 1
fi
if [ "$#" -eq 0 ]; then set -- run; fi
exec "$bundle_python" -B "$bundle_root/workflow/run.py" "$@"
