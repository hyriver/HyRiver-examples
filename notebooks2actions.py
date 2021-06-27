#!/usr/bin/env python3

import fileinput
from pathlib import Path

nb_list = [f.stem for f in Path("notebooks").glob("*.ipynb") if "nid_dib" not in f.stem]
for line in fileinput.input(".github/workflows/test.yml", inplace=True):
    if "notebook: [" in line:
        print(f"        notebook: {nb_list}", end="")
    else:
        print(line, end="")