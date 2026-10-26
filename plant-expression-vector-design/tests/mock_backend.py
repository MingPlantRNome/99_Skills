#!/usr/bin/env python3
"""PROTOCOL TEST DOUBLE ONLY. Does not load Evo2 or evaluate biology."""
import csv
import sys
from pathlib import Path

if __name__ == "__main__":
    input_file, output_file, weights = map(Path, sys.argv[1:])
    assert weights.is_file()
    names = [line[1:].split()[0] for line in input_file.read_text().splitlines() if line.startswith(">")]
    with output_file.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["sequence_name", "probability"])
        for name in names:
            writer.writerow([name, -1.0])
