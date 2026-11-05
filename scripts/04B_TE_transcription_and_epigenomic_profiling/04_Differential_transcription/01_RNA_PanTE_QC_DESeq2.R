#!/usr/bin/env Rscript
# Public-release note: original executed source; replace local absolute paths with paths from your environment before running.

suppressPackageStartupMessages({
    library(DESeq2)
    library(ggplot2)
    library(pheatmap)
})

ROOT <- file.path(Sys.getenv("FB_ROOT", unset = "PATH_TO_FB_ROOT"), "TE", "PanTE_v52", "05_TE_RNA")

COUNT_FILE <- file.path(ROOT, "04_matrix/PanTE_raw_counts.tsv")
META_FILE  <- file.path(ROOT, "00_metadata/sample_metadata.tsv")
PRES_FILE  <- file.path(ROOT, "04_matrix/PanTE_genome_presence.tsv")

QC_DIR  <- file.path(ROOT, "05_DE/00_QC")
FLT_DIR <- file.path(ROOT, "05_DE/01_filter")
PCA_DIR <- file.path(ROOT, "05_DE/02_VST_PCA")
DE_DIR  <- file.path(ROOT, "05_DE/03_DESeq2")

dir.create(QC_DIR,  recursive=TRUE, showWarnings=FALSE)
dir.create(FLT_DIR, recursive=TRUE, showWarnings=FALSE)
dir.create(PCA_DIR, recursive=TRUE, showWarnings=FALSE)
dir.create(DE_DIR,  recursive=TRUE, showWarnings=FALSE)

# ============================================================
# 1. Input
# ============================================================

message("===== Reading input =====")

count_df <- read.delim(
    COUNT_FILE,
    header=TRUE,
    check.names=FALSE,
    stringsAsFactors=FALSE
)

rownames(count_df) <- count_df$PanTE_ID
count_df$PanTE_ID <- NULL

cts <- as.matrix(count_df)
storage.mode(cts) <- "integer"

meta <- read.delim(
    META_FILE,
    header=TRUE,
    stringsAsFactors=FALSE,
    check.names=FALSE
)

required_meta <- c("Sample", "Group", "Tissue")

if (!all(required_meta %in% colnames(meta))) {
    stop(
        paste(
            "Metadata must contain:",
            paste(required_meta, collapse=", ")
        )
    )
}

rownames(meta) <- meta$Sample

if (!all(meta$Sample %in% colnames(cts))) {
    stop("Some metadata samples are absent from count matrix.")
}

# exact metadata order
cts <- cts[, meta$Sample, drop=FALSE]

meta$Group <- factor(
    meta$Group,
    levels=c("BSB", "TC", "BT")
)

meta$Tissue <- factor(
    meta$Tissue,
    levels=c(
        "Brain",
        "Gill",
        "Kidney",
        "Liver",
        "Muscle"
    )
)

presence <- read.delim(
    PRES_FILE,
    header=TRUE,
    stringsAsFactors=FALSE
)

rownames(presence) <- presence$PanTE_ID

presence$Comparable_all3 <- ifelse(
    presence$BSB_present == 1 &
    presence$TC_present  == 1 &
    presence$BT_present  == 1,
    1, 0
)

message("Samples: ", ncol(cts))
message("PanTE:   ", nrow(cts))

stopifnot(ncol(cts) == 57)
stopifnot(nrow(cts) == 11950)

# ============================================================
# 2. Raw sample QC
# ============================================================

raw_qc <- data.frame(
    Sample=colnames(cts),
    Group=meta[colnames(cts), "Group"],
    Tissue=meta[colnames(cts), "Tissue"],
    Raw_PanTE_count=colSums(cts),
    Detected_PanTE=colSums(cts > 0),
    Count_ge_10=colSums(cts >= 10),
    stringsAsFactors=FALSE
)

write.table(
    raw_qc,
    file=file.path(QC_DIR, "RNA_PanTE_raw_sample_QC.tsv"),
    sep="\t",
    quote=FALSE,
    row.names=FALSE
)

# ============================================================
# Helper: PCA + correlation
# ============================================================

run_pca_corr <- function(vsd, sample_meta, prefix, outdir) {

    mat <- assay(vsd)

    # ------------------------
    # PCA
    # ------------------------
    pc <- prcomp(
        t(mat),
        center=TRUE,
        scale.=FALSE
    )

    pct <- (
        pc$sdev^2 /
        sum(pc$sdev^2)
    ) * 100

    pca_df <- data.frame(
        Sample=rownames(pc$x),
        Group=sample_meta[rownames(pc$x), "Group"],
        Tissue=sample_meta[rownames(pc$x), "Tissue"],
        PC1=pc$x[,1],
        PC2=pc$x[,2],
        stringsAsFactors=FALSE
    )

    write.table(
        pca_df,
        file=file.path(
            outdir,
            paste0(prefix, ".PCA_coordinates.tsv")
        ),
        sep="\t",
        quote=FALSE,
        row.names=FALSE
    )

    p <- ggplot(
        pca_df,
        aes(
            x=PC1,
            y=PC2,
            color=Group,
            shape=Tissue
        )
    ) +
        geom_point(size=3.2) +
        geom_text(
            aes(label=Sample),
            size=2.2,
            vjust=-0.7,
            check_overlap=TRUE
        ) +
        labs(
            x=sprintf("PC1 (%.1f%%)", pct[1]),
            y=sprintf("PC2 (%.1f%%)", pct[2]),
            title=paste0(prefix, " PanTE expression PCA")
        ) +
        theme_bw(base_size=11)

    ggsave(
        file.path(
            outdir,
            paste0(prefix, ".PCA.pdf")
        ),
        p,
        width=8,
        height=6
    )

    # ------------------------
    # Pearson correlation
    # ------------------------

    cor_mat <- cor(
        mat,
        method="pearson"
    )

    write.table(
        cor_mat,
        file=file.path(
            outdir,
            paste0(prefix, ".sample_correlation.tsv")
        ),
        sep="\t",
        quote=FALSE
    )

    ann <- data.frame(
        Group=sample_meta[colnames(mat), "Group"],
        Tissue=sample_meta[colnames(mat), "Tissue"]
    )

    rownames(ann) <- colnames(mat)

    pdf(
        file.path(
            outdir,
            paste0(prefix, ".sample_correlation.pdf")
        ),
        width=9,
        height=8
    )

    pheatmap(
        cor_mat,
        annotation_col=ann,
        annotation_row=ann,
        border_color=NA,
        fontsize=7
    )

    dev.off()
}

# ============================================================
# 3. Global QC
# ============================================================

message("===== Global VST/PCA =====")

# QC filter only:
# >=10 counts in at least 3 samples
global_keep <- rowSums(cts >= 10) >= 3

global_filter_summary <- data.frame(
    Total_PanTE=nrow(cts),
    Passed_PanTE=sum(global_keep),
    Removed_PanTE=sum(!global_keep)
)

write.table(
    global_filter_summary,
    file=file.path(
        FLT_DIR,
        "global_filter_summary.tsv"
    ),
    sep="\t",
    quote=FALSE,
    row.names=FALSE
)

dds_global <- DESeqDataSetFromMatrix(
    countData=cts[global_keep, , drop=FALSE],
    colData=meta,
    design=~1
)

dds_global <- estimateSizeFactors(dds_global)

global_sf <- data.frame(
    Sample=colnames(dds_global),
    Group=meta[colnames(dds_global), "Group"],
    Tissue=meta[colnames(dds_global), "Tissue"],
    SizeFactor=sizeFactors(dds_global),
    stringsAsFactors=FALSE
)

write.table(
    global_sf,
    file=file.path(
        QC_DIR,
        "global_DESeq2_size_factors.tsv"
    ),
    sep="\t",
    quote=FALSE,
    row.names=FALSE
)

vsd_global <- varianceStabilizingTransformation(
    dds_global,
    blind=TRUE
)

run_pca_corr(
    vsd_global,
    meta,
    "All57",
    PCA_DIR
)

# ============================================================
# 4. Tissue-wise DESeq2
# ============================================================

tissues <- c(
    "Brain",
    "Gill",
    "Kidney",
    "Liver",
    "Muscle"
)

de_summary <- list()

for (tissue in tissues) {

    message("")
    message("=====================================")
    message("Tissue: ", tissue)
    message("=====================================")

    samples <- rownames(meta)[meta$Tissue == tissue]

    sub_meta <- droplevels(
        meta[samples, , drop=FALSE]
    )

    # Restore 3-group level order
    sub_meta$Group <- factor(
        sub_meta$Group,
        levels=c("BSB", "TC", "BT")
    )

    sub_cts <- cts[, samples, drop=FALSE]

    # Main low-expression filter:
    # at least 3 biological samples have count >= 10
    keep <- rowSums(sub_cts >= 10) >= 3

    message(
        "Filter: ",
        sum(keep),
        " / ",
        nrow(sub_cts),
        " retained"
    )

    tissue_dir <- file.path(DE_DIR, tissue)
    dir.create(
        tissue_dir,
        recursive=TRUE,
        showWarnings=FALSE
    )

    filter_df <- data.frame(
        Tissue=tissue,
        Total_PanTE=nrow(sub_cts),
        Passed_PanTE=sum(keep),
        Removed_PanTE=sum(!keep),
        stringsAsFactors=FALSE
    )

    write.table(
        filter_df,
        file=file.path(
            tissue_dir,
            "filter_summary.tsv"
        ),
        sep="\t",
        quote=FALSE,
        row.names=FALSE
    )

    # Save retained IDs
    write.table(
        data.frame(
            PanTE_ID=rownames(sub_cts)[keep]
        ),
        file=file.path(
            tissue_dir,
            "retained_PanTE_ID.tsv"
        ),
        sep="\t",
        quote=FALSE,
        row.names=FALSE
    )

    # ------------------------
    # DESeq2
    # ------------------------

    dds <- DESeqDataSetFromMatrix(
        countData=sub_cts[keep, , drop=FALSE],
        colData=sub_meta,
        design=~Group
    )

    dds <- DESeq(
        dds,
        quiet=TRUE
    )

    # ------------------------
    # sample size factors
    # ------------------------

    sf_df <- data.frame(
        Sample=colnames(dds),
        Group=sub_meta[colnames(dds), "Group"],
        Tissue=tissue,
        SizeFactor=sizeFactors(dds),
        stringsAsFactors=FALSE
    )

    write.table(
        sf_df,
        file=file.path(
            tissue_dir,
            "DESeq2_size_factors.tsv"
        ),
        sep="\t",
        quote=FALSE,
        row.names=FALSE
    )

    # ------------------------
    # normalized counts
    # ------------------------

    norm <- counts(
        dds,
        normalized=TRUE
    )

    norm_out <- data.frame(
        PanTE_ID=rownames(norm),
        norm,
        check.names=FALSE
    )

    write.table(
        norm_out,
        file=file.path(
            tissue_dir,
            "PanTE_normalized_counts.tsv"
        ),
        sep="\t",
        quote=FALSE,
        row.names=FALSE
    )

    # ------------------------
    # normalized group means
    # ------------------------

    means <- data.frame(
        PanTE_ID=rownames(norm),
        BSB_mean=rowMeans(
            norm[, sub_meta$Group == "BSB", drop=FALSE]
        ),
        TC_mean=rowMeans(
            norm[, sub_meta$Group == "TC", drop=FALSE]
        ),
        BT_mean=rowMeans(
            norm[, sub_meta$Group == "BT", drop=FALSE]
        ),
        stringsAsFactors=FALSE
    )

    means$APE <- (
        means$BSB_mean +
        means$TC_mean
    ) / 2

    pres <- presence[
        match(
            means$PanTE_ID,
            presence$PanTE_ID
        ),
        c(
            "BSB_present",
            "TC_present",
            "BT_present",
            "Comparable_all3"
        ),
        drop=FALSE
    ]

    means <- cbind(
        means,
        pres
    )

    write.table(
        means,
        file=file.path(
            tissue_dir,
            "PanTE_normalized_group_means.tsv"
        ),
        sep="\t",
        quote=FALSE,
        row.names=FALSE
    )

    # ------------------------
    # VST PCA / correlation
    # ------------------------

    vsd <- varianceStabilizingTransformation(
        dds,
        blind=TRUE
    )

    run_pca_corr(
        vsd,
        sub_meta,
        tissue,
        tissue_dir
    )

    # ========================================================
    # Differential contrasts
    # ========================================================

    contrasts <- list(
        BSB_vs_TC=c("Group", "BSB", "TC"),
        BT_vs_BSB=c("Group", "BT", "BSB"),
        BT_vs_TC=c("Group", "BT", "TC")
    )

    for (nm in names(contrasts)) {

        rr <- results(
            dds,
            contrast=contrasts[[nm]],
            alpha=0.05
        )

        r <- as.data.frame(rr)

        r$PanTE_ID <- rownames(r)

        r <- r[
            ,
            c(
                "PanTE_ID",
                "baseMean",
                "log2FoldChange",
                "lfcSE",
                "stat",
                "pvalue",
                "padj"
            )
        ]

        # add group means
        m2 <- means[
            match(r$PanTE_ID, means$PanTE_ID),
            c(
                "BSB_mean",
                "TC_mean",
                "BT_mean",
                "APE",
                "BSB_present",
                "TC_present",
                "BT_present",
                "Comparable_all3"
            )
        ]

        r <- cbind(
            r,
            m2
        )

        # Main significance flag
        r$Significant <- ifelse(
            !is.na(r$padj) &
            r$padj < 0.05 &
            abs(r$log2FoldChange) >= 1,
            1,
            0
        )

        # sort
        r <- r[
            order(
                is.na(r$padj),
                r$padj
            ),
        ]

        outfile <- file.path(
            tissue_dir,
            paste0(nm, ".DESeq2.tsv")
        )

        write.table(
            r,
            file=outfile,
            sep="\t",
            quote=FALSE,
            row.names=FALSE
        )

        n_up <- sum(
            !is.na(r$padj) &
            r$padj < 0.05 &
            r$log2FoldChange >= 1
        )

        n_down <- sum(
            !is.na(r$padj) &
            r$padj < 0.05 &
            r$log2FoldChange <= -1
        )

        de_summary[[length(de_summary)+1]] <- data.frame(
            Tissue=tissue,
            Contrast=nm,
            Tested_PanTE=nrow(r),
            Significant_up=n_up,
            Significant_down=n_down,
            Significant_total=n_up+n_down,
            stringsAsFactors=FALSE
        )
    }
}

# ============================================================
# 5. Combined DE summary
# ============================================================

de_summary_df <- do.call(
    rbind,
    de_summary
)

write.table(
    de_summary_df,
    file=file.path(
        DE_DIR,
        "DESeq2_summary_all_tissues.tsv"
    ),
    sep="\t",
    quote=FALSE,
    row.names=FALSE
)

message("")
message("=====================================")
message("RNA PanTE QC + DESeq2 completed")
message("=====================================")
message("Global PCA: ", PCA_DIR)
message("DE results: ", DE_DIR)
message("")
print(de_summary_df)