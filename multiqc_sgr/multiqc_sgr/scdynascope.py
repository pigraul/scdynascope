import json
import logging
from collections import defaultdict

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph, box, linegraph

# Initialise the logger
log = logging.getLogger("multiqc")


def get_int(x):
    return str(int(x))


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent object
        super().__init__(
            name="scdynascope",
            anchor="scdynascope",
            info="mapping, demultiplexing and quantification for single cell DynaSCOPE RNA-seq",
        )
        log.info(f"Running module: {self.name}")

        stat_data = self.parse_json(self.name, "stats")
        umi_count_data = self.parse_json(self.name, "umi_count")
        saturation_data = self.parse_json(self.name, "saturation")
        median_gene_data = self.parse_json(self.name, "median_gene")
        substitution_data = self.parse_json(self.name, "substitution")
        tor_data = self.parse_json(self.name, "tor")

        if all(
            len(x) == 0
            for x in [stat_data, umi_count_data, saturation_data, median_gene_data, substitution_data, tor_data]
        ):
            raise ModuleNoSamplesFound

        # Basic Stats Table
        self.general_stats_table(stat_data)

        # barcode rank plot
        self.add_section(
            name="Barcode Rank", anchor="scdynascope_barcode_rank", plot=self.barcode_rank_plot(umi_count_data)
        )

        # subsample
        if saturation_data:
            self.add_section(
                name="Saturation", anchor="scdynascope_subsample", plot=self.saturation_plot(saturation_data)
            )
        if median_gene_data:
            self.add_section(
                name="Median Gene", anchor="scdynascope_median_gene", plot=self.median_gene_plot(median_gene_data)
            )

        # substitution
        self.add_section(
            name="Substitution Rate", anchor="scdynascope_substitution", plot=self.substitution_plot(substitution_data)
        )

        if tor_data:
            self.add_section(name="Turn-over Rate", anchor="scdynascope_labeled_summary", plot=self.tor_plot(tor_data))

        # Superfluous function call to confirm that it is used in this module
        # Replace None with actual version if it is available
        self.add_software_version(None)

    def parse_json(self, assay, seg):
        data_dict = defaultdict(dict)
        n = 0
        for f in self.find_log_files(f"{assay}/{seg}"):
            log.debug(f"Found file: {f['fn']}")
            n += 1
            parsed_data = json.loads(f["f"])
            if parsed_data is not None:
                x = f["s_name"]
                s_name = x[: x.find(f".{assay}")]
                if s_name in data_dict:
                    log.info(f"Duplicate sample name found! Update: {s_name}")
                self.add_data_source(f, s_name=s_name, section=seg)
                data_dict[s_name].update(parsed_data)

        data_dict = self.ignore_samples(data_dict)

        log.info(f"Found {n} {assay} {seg} reports")
        # Write parsed report data to a file
        self.write_data_file(data_dict, f"multiqc_{assay}_{seg}")
        return data_dict

    def general_stats_table(self, summary_data):
        headers = {
            "Protocol": {
                "title": "Protocol",
                "description": "Predefined pattern of barcode and UMI",
                "scale": "purple",
                "hidden": True,
            },
            "Raw Reads": {
                "title": "Raw Reads",
                "description": "Number of reads in the input file",
                "scale": "blue",
                "format": "{:,.0f}",
                "hidden": True,
            },
            "Valid Reads": {
                "title": "Valid Reads",
                "description": "Percent of reads with valid barcodes",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "green",
            },
            "Corrected Barcodes": {
                "title": "Corrected Barcodes",
                "description": "Percent of corrected barcodes",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "green",
                "hidden": True,
            },
            "Reads Mapped To Unique Loci": {
                "title": "Unique Reads",
                "description": "Percent of valid reads mapped to unique loci on the genome",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "green",
                "hidden": True,
            },
            "Reads Mapped To Multiple Loci": {
                "title": "Multi Reads",
                "description": "Percent of valid reads mapped to multiple loci on the genome",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "green",
                "hidden": True,
            },
            "Reads Mapped Uniquely To Transcriptome": {
                "title": "Counted Unique Reads",
                "description": "Percent of valid reads mapped uniquely to transcriptome; These reads are used for UMI counting",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "green",
            },
            "Mapped Reads Assigned To Exonic Regions": {
                "title": "Exonic Reads",
                "description": "Percent of mapped reads assigned to exonic regions",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "green",
                "hidden": True,
            },
            "Mapped Reads Assigned To Intronic Regions": {
                "title": "Intronic Reads",
                "description": "Percent of mapped reads assigned to intronic regions",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "green",
                "hidden": True,
            },
            "Mapped Reads Assigned To Intergenic Regions": {
                "title": "Intergenic Reads",
                "description": "Percent of mapped reads assigned to intergenic regions",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "green",
                "hidden": True,
            },
            "Mapped Reads Assigned Antisense To Gene": {
                "title": "Antisense Reads",
                "description": "Percent of mapped reads assigned antisense to gene",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "green",
                "hidden": True,
            },
            "Estimated Number of Cells": {
                "title": "Number of Cells",
                "description": "Estimated number of cells",
                "scale": "blue",
                "format": "{:,.0f}",
            },
            "Fraction Reads in Cells": {
                "title": "Reads in Cells",
                "description": "Percent of unique reads in cells",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "green",
            },
            "Mean Used Reads per Cell": {
                "title": "Mean Used Reads",
                "description": "Mean number of reads used per cell",
                "scale": "blue",
                "format": "{:,.0f}",
            },
            "Median UMI per Cell": {
                "title": "Median UMI",
                "description": "Median number of UMIs per cell",
                "scale": "blue",
                "format": "{:,.0f}",
            },
            "Median Genes per Cell": {
                "title": "Median Genes",
                "description": "Median number of genes per cell",
                "scale": "blue",
                "format": "{:,.0f}",
            },
            "Total Genes": {
                "title": "Total Genes",
                "description": "Total number of genes detected",
                "scale": "blue",
                "format": "{:,.0f}",
            },
            "Saturation": {
                "title": "Saturation",
                "description": "Percent of reads originating from an already-observed UMI",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "green",
            },
            "Labeled rate": {
                "title": "Labeled rate",
                "description": "Percent of UMI labeled",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "orange",
            },
        }
        self.general_stats_addcols(summary_data, headers=headers)

    def barcode_rank_plot(self, umi_count_data):
        plot_data = {}
        colors = {}
        for sample in umi_count_data:
            for sub in umi_count_data[sample]:
                cur = umi_count_data[sample][sub]
                if not cur:
                    continue
                new = {}
                for k, v in cur.items():
                    new[int(k)] = v
                plot_data[sub] = new
                if "pure" in sub:
                    colors[sub] = "darkblue"
                elif "mix" in sub:
                    colors[sub] = "lightblue"
                elif "background" in sub:
                    colors[sub] = "lightgray"

        # Config for the plot
        pconfig = {
            "id": "scdynascope_barcode_rank_plot",
            "title": "scdynascope: Barcode Rank",
            "ylab": "UMI counts",
            "xlab": "Barcode Rank",
            "yLog": True,
            "xLog": True,
            "colors": colors,
            "ymin": 0,
            "height": 750,
        }

        return linegraph.plot(plot_data, pconfig)

    def saturation_plot(self, saturation_data):
        # Config for the plot
        pconfig = {
            "id": "scdynascope_saturation_plot",
            "title": "scdynascope: Saturation",
            "ylab": "Saturation",
            "xlab": "Percent of Reads",
            "height": 750,
        }

        return linegraph.plot(saturation_data, pconfig)

    def median_gene_plot(self, median_gene_data):
        # Config for the plot
        pconfig = {
            "id": "scdynascope_median_gene_plot",
            "title": "scdynascope: Median Gene",
            "ylab": "Median Gene",
            "xlab": "Percent of Reads",
            "height": 750,
        }
        return linegraph.plot(median_gene_data, pconfig)

    def substitution_plot(self, substitution_data):
        # Config for the plot
        pconfig = {
            "id": "scdynascope_substitution_plot",
            "title": "scdynascope: Substitution Rate",
            "ylab": "Substitution Rate %",
            "xlab": "Type",
            "stacking": "group",
            "cpswitch": False,
        }
        cats = [
            {
                "A_to_C": {"color": "#c2e699"},
                "A_to_G": {"color": "#78c679"},
                "A_to_T": {"color": "#238443"},
                "C_to_A": {"color": "#fbb4b9"},
                "C_to_G": {"color": "#f768a1"},
                "C_to_T": {"color": "#ae017e"},
                "G_to_A": {"color": "#bdd7e7"},
                "G_to_C": {"color": "#6baed6"},
                "G_to_T": {"color": "#2171b5"},
                "T_to_A": {"color": "#fcae91"},
                "T_to_C": {"color": "#fb6a4a"},
                "T_to_G": {"color": "#cb181d"},
            }
        ]
        return bargraph.plot(substitution_data, cats=cats, pconfig=pconfig)

    def tor_plot(self, tor_data):
        v_sample = defaultdict(dict)
        for sample in tor_data:
            for v in tor_data[sample]:
                v_sample[v].update({sample: tor_data[sample][v]})
        cats = sorted(list(v_sample.keys()))
        plot_data = [v_sample[v] for v in cats]

        pconfig = {
            "id": "scdynascope_tor_plot",
            "title": "scdynascope: Turn-over Rate",
            "ylab": "TOR",
            "data_labels": cats,
        }

        return box.plot(plot_data, pconfig)
