#!/bin/env python

import argparse
import os
import re

import pandas as pd
import pysam
import utils
from __init__ import ASSAY


class Substitution:
    """
    ## Features
    - Computes the overall conversion rates in reads and plots a barplot.

    ## Output
    - `{sample}.substitution.txt` Tab-separated table of the overall conversion rates.
    """

    def __init__(self, args):
        self.args = args
        # input
        self.sample = args.sample
        self.outdir = args.outdir
        self.in_bam = args.bam
        self.snp_file = args.bg_snp
        self.all_type_plot = args.all_type_plot
        # set
        self.index_dict = {}
        self.report_dict = {"TC":{}, "All":{}}
        self.df = pd.DataFrame()
        # output files
        self.outstat = os.path.join(self.outdir, self.sample + ".substitution.tsv")
        self.tcstat = os.path.join(self.outdir, self.sample + ".TC_substitution.tsv")

    def run(self):
        bg_snp = self.background_snp()
        for_base, rev_base, is_forward, is_reverse, is_snp = self.get_sub_tag(
            self.in_bam, bg_snp
        )
        self.sub_stat(for_base, rev_base, is_forward, is_reverse, is_snp)
        utils.write_multiqc(self.report_dict, self.sample, ASSAY, "substitution")


    def get_sub_tag(self, bam, bg):
        save = pysam.set_verbosity(0)
        bamfile = pysam.AlignmentFile(bam, "rb", require_index=False)
        pysam.set_verbosity(save)
        is_reverse = {
            "cA": 0,
            "gA": 0,
            "tA": 0,
            "aC": 0,
            "gC": 0,
            "tC": 0,
            "aG": 0,
            "cG": 0,
            "tG": 0,
            "aT": 0,
            "cT": 0,
            "gT": 0,
        }
        is_forward = {
            "cA": 0,
            "gA": 0,
            "tA": 0,
            "aC": 0,
            "gC": 0,
            "tC": 0,
            "aG": 0,
            "cG": 0,
            "tG": 0,
            "aT": 0,
            "cT": 0,
            "gT": 0,
        }
        is_snp = {"reverse": 0, "forward": 0}
        for_base = {"a": 0, "c": 0, "g": 0, "t": 0}
        rev_base = {"a": 0, "c": 0, "g": 0, "t": 0}
        snp_tags = [
            "",
            "cA",
            "gA",
            "tA",
            "aC",
            "gC",
            "tC",
            "aG",
            "cG",
            "tG",
            "aT",
            "cT",
            "gT",
        ]
        ref_tags = ["", "a", "c", "g", "t"]

        for read in bamfile.fetch(until_eof=True):
            try:
                snpmatch = re.match(
                    r"cA(\d+);gA(\d+);tA(\d+);aC(\d+);gC(\d+);tC(\d+);aG(\d+);cG(\d+);tG(\d+);aT(\d+);cT(\d+);gT(\d+);",
                    read.get_tag("SC"),
                    re.M,
                )
                totmatch = re.match(
                    r"a(\d+);c(\d+);g(\d+);t(\d+)", read.get_tag("TC"), re.M
                )
                if snpmatch and totmatch:
                    if read.is_reverse:
                        for j in range(1, len(ref_tags)):
                            rev_base[ref_tags[j]] += int(totmatch.group(j))
                        for i in range(1, len(snp_tags)):
                            is_reverse[snp_tags[i]] += int(snpmatch.group(i))
                        stag = read.get_tag("AL")
                        rtag = "reverse"
                    else:
                        for j in range(1, len(ref_tags)):
                            for_base[ref_tags[j]] += int(totmatch.group(j))
                        for i in range(1, len(snp_tags)):
                            is_forward[snp_tags[i]] += int(snpmatch.group(i))
                        stag = read.get_tag("TL")
                        rtag = "forward"
                    if len(stag) == 1 and stag[0] == 0:
                        continue
                    else:
                        chro = read.reference_name
                        for si in range(0, len(stag)):
                            pos = chro + "_" + str(stag[si])
                            if pos in bg:
                                is_snp[rtag] += 1

            except (ValueError, KeyError):
                continue
        bamfile.close()

        return for_base, rev_base, is_forward, is_reverse, is_snp

    def sub_stat(self, for_base, rev_base, is_forward, is_reverse, is_snp):
        convertdict = {
            "a": ["aC", "aG", "aT"],
            "c": ["cA", "cG", "cT"],
            "g": ["gA", "gC", "gT"],
            "t": ["tA", "tC", "tG"],
        }
        subdict = {
            "a": "t",
            "t": "a",
            "c": "g",
            "g": "c",
            "aC": "tG",
            "aG": "tC",
            "aT": "tA",
            "cA": "gT",
            "cG": "gC",
            "cT": "gA",
            "gA": "cT",
            "gC": "cG",
            "gT": "cA",
            "tA": "aT",
            "tC": "aG",
            "tG": "aC",
        }
        outdict = {
            "aC": "A_to_C",
            "aG": "A_to_G",
            "aT": "A_to_T",
            "cA": "C_to_A",
            "cG": "C_to_G",
            "cT": "C_to_T",
            "gA": "G_to_A",
            "gC": "G_to_C",
            "gT": "G_to_T",
            "tA": "T_to_A",
            "tC": "T_to_C",
            "tG": "T_to_G",
        }

        if self.all_type_plot:
            outw = open(self.outstat, "w")
            for x in ["a", "c", "g", "t"]:
                fbase = for_base[x]
                rbase = rev_base[subdict[x]]
                for y in convertdict[x]:
                    fcov = is_forward[y] * 100 / float(fbase) if float(fbase) > 0 else 0
                    rcov = is_reverse[subdict[y]] * 100 / float(rbase) if float(rbase) > 0 else 0
                    self.report_dict["All"][outdict[y]] = fcov + rcov
                    outw.write(f"{outdict[y]}\t{fcov:.3f}\t{rcov:.3f}\n")
            outw.close()

        outw2 = open(self.tcstat, "w")
        fbase = for_base["t"]
        rbase = rev_base["a"]
        for_label = (is_forward["tC"] - is_snp["forward"]) * 100 / float(fbase) if float(fbase) > 0 else 0
        rev_label = (is_reverse["aG"] - is_snp["reverse"]) * 100 / float(rbase) if float(rbase) > 0 else 0
        for_snp = is_snp["forward"] * 100 / float(fbase) if float(fbase) > 0 else 0
        rev_snp = is_snp["reverse"] * 100 / float(rbase) if float(rbase) > 0 else 0
        outw2.write(f"Labeled\t{for_label:.3f}\t{rev_label:.3f}\n")
        outw2.write(f"Background\t{for_snp:.3f}\t{rev_snp:.3f}\n")
        outw2.close()
        self.report_dict["TC"]["Labeled"] = for_label + rev_label
        self.report_dict["TC"]["Background"] = for_snp + rev_snp


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


def get_opts_substitution():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("--outdir", help="Output diretory.", required=True)
    parser.add_argument("--sample", help="Sample name.", required=True)
    parser.add_argument("--bam", help="Required. bam file from step conversion.", required=True)
    parser.add_argument("--bg_snp", nargs="+", required=False, help="background snp file, csv or vcf format")
    parser.add_argument(
        "--all_type_plot",
        action="store_true",
        help="Plot subsititution rate for all conversion type",
    )
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = get_opts_substitution()
    Substitution(args).run()
