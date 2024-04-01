#!/usr/bin/env python3
"""Update .github/workflows/test.yml with list of notebooks to test."""

import fileinput
from pathlib import Path

nb_list = sorted(f.stem for f in Path("notebooks").glob("*.ipynb") if "nid_dib" not in f.stem)
for line in fileinput.input(".github/workflows/test.yml", inplace=True):
    if "notebook: [" in line:
        print(f"        notebook: {nb_list}", end="")
        print()
    else:
        print(line, end="")
