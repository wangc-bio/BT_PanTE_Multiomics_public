#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import print_function

import os
import numpy as np
import pandas as pd


# ============================================================
# PATHS
# ============================================================

BASE = os.environ.get("KIMURA_BASE", "PATH_TO_RepeatMasker_unified_v5_2_Kimura_align")

KIM = os.path.join(
    BASE,
    "Kimura_divergence"
)

OUT = os.path.join(
    KIM,
    "01_parsed"
)

if not os.path.exists(OUT):
    os.makedirs(OUT)


FILES = {
    "BSB": os.path.join(
        KIM,
        "BSB",
        "BSB.divsum"
    ),

    "TC": os.path.join(
        KIM,
        "TC",
        "TC.divsum"
    ),

    "BTM": os.path.join(
        KIM,
        "BTM",
        "BTM.divsum"
    ),

    "BTP": os.path.join(
        KIM,
        "BTP",
        "BTP.divsum"
    )
}


for sample in FILES:

    if not os.path.isfile(
        FILES[sample]
    ):

        raise RuntimeError(
            "Missing divsum file: {}".format(
                FILES[sample]
            )
        )


# ============================================================
# TE CLASS MAPPING
# ============================================================

MAIN5 = [
    "DNA",
    "LTR",
    "LINE",
    "SINE",
    "RC"
]


def broad_class(cls):

    x = str(
        cls
    ).strip()

    xl = x.lower()

    if (
        x == "DNA"
        or
        x.startswith("DNA/")
    ):
        return "DNA"

    if (
        x == "LTR"
        or
        x.startswith("LTR/")
    ):
        return "LTR"

    if (
        x == "LINE"
        or
        x.startswith("LINE/")
    ):
        return "LINE"

    if (
        x == "SINE"
        or
        x.startswith("SINE/")
    ):
        return "SINE"

    if (
        x == "RC"
        or
        x.startswith("RC/")
        or
        xl.startswith("rolling-circle")
        or
        xl.startswith("rolling_circle")
    ):
        return "RC"

    return "Other"


# ============================================================
# PARSE ONE DIVSUM
# ============================================================

def parse_one(
    sample,
    path
):

    family_rows = []
    coverage_rows = []

    family_mode = False
    coverage_mode = False

    coverage_classes = None

    cpg_adjusted = False


    with open(
        path,
        "r"
    ) as fh:

        for raw in fh:

            line = raw.rstrip(
                "\n"
            )

            stripped = line.strip()


            # ----------------------------------------------
            # CpG correction flag
            # ----------------------------------------------

            if (
                "adjusted for CpG sites"
                in line
            ):

                cpg_adjusted = True


            # ----------------------------------------------
            # family-level section
            # ----------------------------------------------

            if line.startswith(
                "Class\tRepeat\tabsLen\twellCharLen\tKimura%"
            ):

                family_mode = True
                coverage_mode = False

                continue


            # ----------------------------------------------
            # Detect coverage matrix directly from "Div ..."
            # This is more robust than relying on section title.
            # ----------------------------------------------

            if stripped.startswith(
                "Div "
            ):

                fields = stripped.split()

                if len(fields) > 2:

                    coverage_classes = fields[
                        1:
                    ]

                    coverage_mode = True
                    family_mode = False

                    continue


            # ----------------------------------------------
            # parse family rows
            # ----------------------------------------------

            if family_mode:

                if not stripped:
                    continue

                if stripped.startswith(
                    "-----"
                ):
                    continue


                fields = line.split(
                    "\t"
                )


                if len(fields) != 5:
                    continue


                cls = fields[0]
                repeat = fields[1]


                try:

                    abs_len = int(
                        fields[2]
                    )

                    well_len = int(
                        fields[3]
                    )

                except ValueError:

                    continue


                try:

                    kimura = float(
                        fields[4]
                    )

                except ValueError:

                    kimura = np.nan


                family_rows.append({
                    "Sample": sample,
                    "Class": cls,
                    "Broad_class": broad_class(
                        cls
                    ),
                    "Repeat": repeat,
                    "absLen": abs_len,
                    "wellCharLen": well_len,
                    "Kimura_pct": kimura
                })


                continue


            # ----------------------------------------------
            # parse divergence coverage matrix
            # ----------------------------------------------

            if (
                coverage_mode
                and
                coverage_classes is not None
            ):

                if not stripped:
                    continue


                fields = stripped.split()


                if len(fields) < 2:
                    continue


                try:

                    divergence = int(
                        fields[0]
                    )

                except ValueError:

                    continue


                values = fields[
                    1:
                ]


                if len(values) != len(
                    coverage_classes
                ):

                    raise RuntimeError(
                        "{}: divergence bin {} has {} "
                        "values but header has {} classes".format(
                            sample,
                            divergence,
                            len(values),
                            len(
                                coverage_classes
                            )
                        )
                    )


                for cls, value in zip(
                    coverage_classes,
                    values
                ):

                    try:

                        bp = int(
                            value
                        )

                    except ValueError:

                        raise RuntimeError(
                            "Invalid coverage value: "
                            "{} Kimura={} Class={} Value={}".format(
                                sample,
                                divergence,
                                cls,
                                value
                            )
                        )


                    coverage_rows.append({
                        "Sample": sample,
                        "Kimura_bin": divergence,
                        "Class": cls,
                        "Broad_class": broad_class(
                            cls
                        ),
                        "Coverage_bp_wellChar": bp
                    })


    # ========================================================
    # convert
    # ========================================================

    family_df = pd.DataFrame(
        family_rows
    )

    coverage_df = pd.DataFrame(
        coverage_rows
    )


    if family_df.empty:

        raise RuntimeError(
            "No family-level rows parsed for {}".format(
                sample
            )
        )


    if coverage_df.empty:

        raise RuntimeError(
            "No Kimura coverage matrix parsed for {}".format(
                sample
            )
        )


    bins = sorted(
        coverage_df[
            "Kimura_bin"
        ].unique().tolist()
    )


    return (
        family_df,
        coverage_df,
        cpg_adjusted,
        bins,
        coverage_classes
    )


# ============================================================
# PARSE FOUR SAMPLES
# ============================================================

all_family = []
all_coverage = []

qc_rows = []


for sample in [
    "BSB",
    "TC",
    "BTM",
    "BTP"
]:

    print("")
    print(
        "Parsing {} ...".format(
            sample
        )
    )


    fam, cov, cpg, bins, classes = parse_one(
        sample,
        FILES[
            sample
        ]
    )


    all_family.append(
        fam
    )

    all_coverage.append(
        cov
    )


    main5_cov = cov[
        cov[
            "Broad_class"
        ].isin(
            MAIN5
        )
    ].copy()


    qc_rows.append({
        "Sample":
            sample,

        "CpG_adjusted":
            cpg,

        "N_repeat_families":
            int(
                len(
                    fam
                )
            ),

        "N_RM_classes":
            int(
                len(
                    classes
                )
            ),

        "N_Kimura_bins":
            int(
                len(
                    bins
                )
            ),

        "Kimura_min":
            int(
                min(
                    bins
                )
            ),

        "Kimura_max":
            int(
                max(
                    bins
                )
            ),

        "Coverage_all_classes_bp":
            int(
                cov[
                    "Coverage_bp_wellChar"
                ].sum()
            ),

        "Coverage_main5_bp":
            int(
                main5_cov[
                    "Coverage_bp_wellChar"
                ].sum()
            )
    })


    print(
        "  families        : {:,}".format(
            len(
                fam
            )
        )
    )

    print(
        "  RM classes      : {}".format(
            len(
                classes
            )
        )
    )

    print(
        "  Kimura bins     : {}-{} (N={})".format(
            min(
                bins
            ),
            max(
                bins
            ),
            len(
                bins
            )
        )
    )

    print(
        "  CpG adjusted    : {}".format(
            cpg
        )
    )


# ============================================================
# CONCATENATE
# ============================================================

family_all = pd.concat(
    all_family,
    ignore_index=True
)

coverage_all = pd.concat(
    all_coverage,
    ignore_index=True
)

qc = pd.DataFrame(
    qc_rows
)


# ============================================================
# CHECK SAME KIMURA BINS ACROSS FOUR SAMPLES
# ============================================================

bin_sets = {}


for sample in [
    "BSB",
    "TC",
    "BTM",
    "BTP"
]:

    bin_sets[
        sample
    ] = tuple(
        sorted(
            coverage_all.loc[
                coverage_all[
                    "Sample"
                ] == sample,
                "Kimura_bin"
            ].unique().tolist()
        )
    )


reference_bins = bin_sets[
    "BSB"
]


for sample in [
    "TC",
    "BTM",
    "BTP"
]:

    if (
        bin_sets[
            sample
        ]
        !=
        reference_bins
    ):

        print(
            "[WARNING] {} Kimura bins differ from BSB.".format(
                sample
            )
        )


# ============================================================
# AGGREGATE TO BROAD CLASS
# ============================================================

broad = (
    coverage_all
    .groupby(
        [
            "Sample",
            "Kimura_bin",
            "Broad_class"
        ],
        as_index=False
    )[
        "Coverage_bp_wellChar"
    ]
    .sum()
)


main5 = broad[
    broad[
        "Broad_class"
    ].isin(
        MAIN5
    )
].copy()


# ============================================================
# MAIN5 TOTALS
# ============================================================

sample_main5_totals = (
    main5
    .groupby(
        "Sample"
    )[
        "Coverage_bp_wellChar"
    ]
    .sum()
    .to_dict()
)


main5[
    "Within_main5_TE_pct"
] = main5.apply(
    lambda row:
        (
            100.0
            *
            row[
                "Coverage_bp_wellChar"
            ]
            /
            sample_main5_totals[
                row[
                    "Sample"
                ]
            ]
        ),
    axis=1
)


# ============================================================
# CLASS TOTALS
# ============================================================

class_totals = (
    main5
    .groupby(
        [
            "Sample",
            "Broad_class"
        ],
        as_index=False
    )[
        "Coverage_bp_wellChar"
    ]
    .sum()
)


class_totals[
    "Within_main5_TE_pct"
] = class_totals.apply(
    lambda row:
        (
            100.0
            *
            row[
                "Coverage_bp_wellChar"
            ]
            /
            sample_main5_totals[
                row[
                    "Sample"
                ]
            ]
        ),
    axis=1
)


# ============================================================
# RM CLASS LIST
# ============================================================

class_list = (
    coverage_all[
        [
            "Sample",
            "Class",
            "Broad_class"
        ]
    ]
    .drop_duplicates()
    .sort_values(
        [
            "Sample",
            "Broad_class",
            "Class"
        ]
    )
)


# ============================================================
# OUTPUT
# ============================================================

FAMILY_OUT = os.path.join(
    OUT,
    "Kimura_family_summary.all_samples.tsv"
)

COVERAGE_OUT = os.path.join(
    OUT,
    "Kimura_coverage_RMclass.long.tsv"
)

BROAD_OUT = os.path.join(
    OUT,
    "Kimura_coverage_broadclass.long.tsv"
)

MAIN5_OUT = os.path.join(
    OUT,
    "Kimura_coverage_main5.long.tsv"
)

TOTAL_OUT = os.path.join(
    OUT,
    "Kimura_main5_class_totals.tsv"
)

CLASS_OUT = os.path.join(
    OUT,
    "Kimura_RM_classes.tsv"
)

QC_OUT = os.path.join(
    OUT,
    "Kimura_parsing_QC.tsv"
)

LOG_OUT = os.path.join(
    OUT,
    "Kimura_parsing_QC.txt"
)


family_all.to_csv(
    FAMILY_OUT,
    sep="\t",
    index=False,
    na_rep="NA"
)

coverage_all.to_csv(
    COVERAGE_OUT,
    sep="\t",
    index=False
)

broad.to_csv(
    BROAD_OUT,
    sep="\t",
    index=False
)

main5.to_csv(
    MAIN5_OUT,
    sep="\t",
    index=False
)

class_totals.to_csv(
    TOTAL_OUT,
    sep="\t",
    index=False
)

class_list.to_csv(
    CLASS_OUT,
    sep="\t",
    index=False
)

qc.to_csv(
    QC_OUT,
    sep="\t",
    index=False
)


with open(
    LOG_OUT,
    "w"
) as out:

    out.write(
        "Kimura divsum parsing QC\n"
    )

    out.write(
        "=========================\n\n"
    )

    out.write(
        qc.to_string(
            index=False
        )
    )

    out.write(
        "\n\n"
    )

    out.write(
        "Main five TE class totals\n"
    )

    out.write(
        "=========================\n\n"
    )

    out.write(
        class_totals.to_string(
            index=False
        )
    )

    out.write(
        "\n"
    )


print("")
print("============================================")
print("Kimura parsing completed successfully.")
print("============================================")

print("")
print(qc.to_string(index=False))

print("")
print("Main five TE class totals:")
print(
    class_totals.to_string(
        index=False
    )
)

print("")
print("Output directory:")
print(
    OUT
)
