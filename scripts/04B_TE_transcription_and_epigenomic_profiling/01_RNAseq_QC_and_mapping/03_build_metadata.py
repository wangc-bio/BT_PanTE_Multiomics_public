#!/usr/bin/env python3
import os
# Public-release note: original executed source; replace local absolute paths with paths from your environment before running.
import re
from pathlib import Path

base = Path(os.path.join(os.environ.get("FB_ROOT", "PATH_TO_FB_ROOT"), "TE", "PanTE_v52"))
out = base / "05_TE_RNA/00_metadata/sample_metadata.tsv"

tissue_map = {
    "gan": "Liver",
    "jirou": "Muscle",
    "nao": "Brain",
    "sai": "Gill",
    "shen": "Kidney",
}

rows = []

for group in ["BSB", "TC", "BT"]:
    pathfile = base / f"{group}.bam.paths"

    with open(pathfile) as fh:
        for line in fh:
            bam = line.strip()
            if not bam:
                continue

            name = Path(bam).name.replace(
                ".Aligned.sortedByCoord.out.bam", ""
            )

            m = re.match(r"YPQ(\d+)(gan|jirou|nao|sai|shen)$", name)
            if not m:
                raise ValueError(f"Cannot parse sample: {name}")

            individual = "YPQ" + m.group(1)
            tissue_code = m.group(2)

            rows.append([
                name,
                group,
                individual,
                tissue_map[tissue_code],
                tissue_code,
                bam
            ])

with open(out, "w") as f:
    f.write(
        "Sample\tGroup\tIndividual\tTissue\tTissue_code\tBAM\n"
    )
    for r in rows:
        f.write("\t".join(r) + "\n")

print("Samples:", len(rows))
print("Output :", out)