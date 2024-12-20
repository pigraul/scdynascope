#!/bin/env python

import argparse
import gzip
import os
import sys

import pandas as pd
import pysam
import scipy.io
import scipy.sparse
import utils
from __init__ import ASSAY, BARCODE_FILE_NAME, FEATURE_FILE_NAME, MATRIX_FILE_NAME

DYNA_MATRIX_DIR_SUFFIX = ["labeled", "unlabeled"]


class Labeled:
    """
    Features
    - Quantify unlabeled and labeled RNA.

    Output
    - `{sample}_labeled_feature_bc_matrix` The labeled expression matrix of cell barcodes & all features in Matrix Market Exchange Formats. (will be deprecated in future versions)
    - `{sample}_unlabeled_feature_bc_matrix` The unlabeled expression matrix of cell barcodes & all features in Matrix Market Exchange Formats. (will be deprecated in future versions)
    - `{sample}_labeled_detail.txt`  tab-delimited  file:
        - Barcode: Cell barcode sequence
        - UMI: UMI sequence
        - geneID: gene ID
        - TC: TC site number in a read (backgroup snp removed)
    """

    def __init__(self, args):
        self.args = args

        # input files
        self.outdir = args.outdir
        self.sample = args.sample
        self.inbam = args.bam
        self.snp_file = args.bg_snp

        # set
        barcodes_file = os.path.join(args.filtered_matrix, BARCODE_FILE_NAME)
        features_file = os.path.join(args.filtered_matrix, FEATURE_FILE_NAME)
        self.features = pd.read_csv(features_file, sep="\t", header=None, index_col=0)
        self.barcodes = utils.read_one_col(barcodes_file)
        self.totaldf = pd.DataFrame()
        self.newdf, self.olddf = pd.DataFrame(), pd.DataFrame()
        self.bg = None
        self.report_dict = {}

        # output files
        # self.h5ad = f'{self.out_prefix}.labeled.h5ad'
        self.detail_txt = f"{self.outdir}/{self.sample}.labeled_detail.txt.gz"
        self.dir_labeled = f"{self.outdir}/{self.sample}.matrix/{DYNA_MATRIX_DIR_SUFFIX[0]}"
        self.dir_unlabeled = f"{self.outdir}/{self.sample}.matrix/{DYNA_MATRIX_DIR_SUFFIX[1]}"

    def run(self):
        # get backgroud snp
        self.bg = self.background_snp()
        # Labeled
        self.totaldf = self.modify_bam()
        self.totaldf.to_csv(self.detail_txt, sep="\t", index=False)
        self.run_quant()
        # report
        utils.write_multiqc(self.report_dict, self.sample, ASSAY, "labeled.stats")

    def run_quant(self):
        newdf = self.totaldf[self.totaldf["TC"] > 0]
        olddf = self.totaldf[self.totaldf["TC"] == 0]
        self.report_dict["Labeled rate"] = utils.get_frac(newdf.shape[0] / self.totaldf.shape[0])
        # matrix
        self.write_sparse_matrix(newdf.drop("TC", axis=1), self.dir_labeled)
        self.write_sparse_matrix(olddf.drop("TC", axis=1), self.dir_unlabeled)

    def modify_bam(self):
        save = pysam.set_verbosity(0)
        bamfile = pysam.AlignmentFile(self.inbam, "rb")
        pysam.set_verbosity(save)
        readdict = {}

        for read in bamfile.fetch(until_eof=True):
            try:
                cb = read.get_tag("CB")
                chro = read.reference_name
                ub = read.get_tag("UB")
                gene = read.get_tag("GX")
                tctag = 0
                true_tc = []

                if read.get_tag("ST") == "+":
                    stag = read.get_tag("TL")
                else:
                    stag = read.get_tag("AL")
                if len(stag) == 1 and stag[0] == 0:
                    tctag = 0
                    true_tc = stag
                else:
                    for si in range(0, len(stag)):
                        pos = chro + "_" + str(stag[si])
                        if pos not in self.bg:
                            true_tc.append(int(stag[si]))
                    tctag = len(true_tc)
                ## dedup: select the most TC read per UMI_gene
                readid = ":".join([cb, ub, gene])
                if readid not in readdict:
                    readdict[readid] = [tctag, read, true_tc]
                else:
                    if tctag > readdict[readid][0]:
                        readdict[readid] = [tctag, read, true_tc]

            except (ValueError, KeyError):
                continue
        bamfile.close()

        ## count df
        tc_df = pd.DataFrame.from_dict(readdict, orient="index", columns=["TC", "read", "loc"])
        tc_df.reset_index(inplace=True)
        ub_df = tc_df["index"].str.split(":", expand=True)
        ub_df.columns = ["Barcode", "UMI", "geneID"]
        outframe = pd.concat([ub_df, tc_df["TC"]], axis=1)

        return outframe

    def background_snp(self):
        outdict = {}
        bgs = []
        for bgargv in self.snp_file:
            if "," in bgargv:
                bgs += bgargv.strip().split(",")
            else:
                bgs.append(bgargv)

        for bgfile in bgs:
            if bgfile.endswith(".csv"):
                df = pd.read_csv(bgfile, dtype={"chrom": str})
                if "pos" in df.columns:
                    df["chrpos"] = df["chrom"] + "_" + df["pos"].astype(str)
                else:  # compatible with previous version
                    df["chrpos"] = df["chrom"] + "_" + df["pos2"].astype(str)
                df1 = df[["chrpos", "convs"]]
                df1.set_index("chrpos", inplace=True)
                for key1 in df1.index.to_list():
                    outdict[key1] = 1
                # outdict.update(df1.to_dict(orient='index'))
            elif bgfile.endswith(".vcf"):
                bcf_in = pysam.VariantFile(bgfile)
                for rec in bcf_in.fetch():
                    try:
                        chrom, pos = rec.chrom, rec.pos
                        chr_pos = chrom + "_" + str(pos - 1)
                        outdict[chr_pos] = 1
                    except (ValueError, KeyError):
                        continue
                bcf_in.close()
            else:
                raise ValueError(
                    "Background snp file format cannot be recognized! Only csv or vcf format."
                )
        return outdict

    def write_sparse_matrix(self, df, matrix_dir):
        count_matrix = self.dataframe_to_matrix(df, features=self.features.index, barcodes=self.barcodes)
        self.to_matrix_dir(count_matrix, matrix_dir)

    def dataframe_to_matrix(
        self, df, features, barcodes, barcode_column="Barcode", feature_column="geneID", value="UMI"
    ):
        if df.shape[0] > 0:
            series_grouped = df.groupby([barcode_column, feature_column], observed=True).size()
            series_grouped.name = value
            df_grouped = pd.DataFrame(series_grouped)
        else:
            empty_matrix = scipy.sparse.coo_matrix((len(features), len(barcodes)))
            return empty_matrix

        feature_index_dict = {}
        for index, gene_id in enumerate(features):
            feature_index_dict[gene_id] = index
        barcode_index_dict = {}
        for index, barcode in enumerate(barcodes):
            barcode_index_dict[barcode] = index

        # use all barcodes
        barcode_codes = [barcode_index_dict[barcode] for barcode in df_grouped.index.get_level_values(level=0)]
        # use all gene_id from features even if it is not in df
        gene_id_codes = [feature_index_dict[gene_id] for gene_id in df_grouped.index.get_level_values(level=1)]
        mtx = scipy.sparse.coo_matrix(
            (df_grouped[value], (gene_id_codes, barcode_codes)), shape=(len(features), len(barcodes))
        )

        return mtx

    def to_matrix_dir(self, count_matrix, matrix_dir):
        self.check_mkdir(dir_name=matrix_dir)
        self.features.to_csv(f"{matrix_dir}/{FEATURE_FILE_NAME}", sep="\t", header=False)
        pd.Series(self.barcodes).to_csv(f"{matrix_dir}/{BARCODE_FILE_NAME}", index=False, sep="\t", header=False)
        matrix_path = f"{matrix_dir}/{MATRIX_FILE_NAME}"
        with gzip.open(matrix_path, "wb") as f:
            scipy.io.mmwrite(f, count_matrix)

    def check_mkdir(self, dir_name):
        """if dir_name is not exist, make one"""
        if not os.path.exists(dir_name):
            os.system(f"mkdir -p {dir_name}")


def get_opts_labeled():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("--outdir", help="Output diretory.", required=True)
    parser.add_argument("--sample", help="Sample name.", required=True)
    parser.add_argument("--filtered_matrix", help="filtered matrix dir", required=True)
    parser.add_argument("--bam", help="tagged bam from conversion step", required=True)
    parser.add_argument("--bg_snp", nargs="+", required=False, help="background snp file, csv or vcf format")
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = get_opts_labeled()
    Labeled(args).run()
