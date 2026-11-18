"""Prepare a benchmark-only copy of the PGEN-enabled plugin source.

Usage: python3 instrument.py SOURCE_CHECKOUT NEW_DIRECTORY
Build the copy with the same Makefile, headers, compiler and flags as the release
binary. Only the profiled copies print timings or support the row-byte check.
"""
import shutil
import sys
from pathlib import Path

def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    source, dest = map(Path, sys.argv[1:])
    dest.mkdir()
    shutil.copytree(source / "src", dest / "src")
    shutil.copyfile(source / "Makefile", dest / "Makefile")
    shutil.copyfile(Path(__file__).with_name("pgen_probe.hpp"), dest / "src/pgen_probe.hpp")
    header = dest / "src/reference.hpp"
    text = header.read_text()
    probes = [
        ("#include <pgenlib_read.h>", '#include <pgenlib_read.h>\n#include "pgen_probe.hpp"'),
        ("    bool opened = false;", "    PgenReadProbe probe;\n    bool opened = false;"),
        ("        PglErr err = hardcalls ? PgrGet(", "        const auto begin = PgenReadProbe::Clock::now();\n        PglErr err = hardcalls ? PgrGet("),
        ('        require(err == kPglRetSuccess, "PGEN read failed:',
         '        const auto decoded = PgenReadProbe::Clock::now();\n        require(err == kPglRetSuccess, "PGEN read failed:'),
        ('        require(d == dosage_ct, "Inconsistent PGEN dosage count");',
         '        require(d == dosage_ct, "Inconsistent PGEN dosage count");\n'
         '        probe.observe(begin, decoded, PgenReadProbe::Clock::now(), row, out, dosage_ct);'),
    ]
    for marker, replacement in probes:
        if text.count(marker) != 1:
            sys.exit(f"Cannot locate unique probe site: {marker}")
        text = text.replace(marker, replacement)
    header.write_text(text)

if __name__ == "__main__":
    main()
