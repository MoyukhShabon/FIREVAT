#!/usr/bin/env Rscript

suppressPackageStartupMessages({
	library(optparse)
	library(FIREVAT)
})

# ---- Default Mutect2 config added inside the package ----
default.config <- system.file("config",
                              "GATK4_Mutect2_config_Tumor_Normal_SB.json",
                              package = "FIREVAT")

# ---- Argument parsing ----
option_list <- list(
	make_option(c("-v", "--vcf"), type = "character", default = NULL,
							help = "Path to input VCF file [required]", metavar = "FILE"),
	make_option(c("-c", "--config"), type = "character", default = default.config,
							help = "Path to config file [default: packaged GATK4_Mutect2_config_Tumor_Normal_SB.json]",
							metavar = "FILE"),
	make_option(c("-o", "--outdir"), type = "character", default = NULL,
							help = "Output directory [required]", metavar = "DIR"),
	make_option(c("-g", "--genome"), type = "character", default = "hg38",
							help = "Reference genome [default: %default]", metavar = "GENOME"),
	make_option(c("-n", "--num-cores"), type = "integer", default = 2,
							help = "Number of cores [default: %default]", metavar = "N")
)

opt <- parse_args(OptionParser(option_list = option_list))

# ---- Validate required arguments ----
required <- c("vcf", "outdir")
missing <- required[vapply(required, function(x) is.null(opt[[x]]), logical(1))]
if (length(missing) > 0) {
	stop("Missing required argument(s): ", paste0("--", missing, collapse = ", "),
			 call. = FALSE)
}

# config has a default, but make sure it actually resolved/exists
if (is.null(opt$config) || opt$config == "") {
	stop("No config supplied and packaged default could not be located. ",
			 "Pass one explicitly with --config.", call. = FALSE)
}

if (!file.exists(opt$vcf))    stop("VCF file not found: ", opt$vcf, call. = FALSE)
if (!file.exists(opt$config)) stop("Config file not found: ", opt$config, call. = FALSE)
if (!dir.exists(opt$outdir))  dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

# ---- Run FIREVAT ----
message("Running FIREVAT")
message("  VCF:     ", opt$vcf)
message("  Config:  ", opt$config)
message("  Outdir:  ", opt$outdir)
message("  Genome:  ", opt$genome)
message("  Cores:   ", opt$`num-cores`)

res <- RunFIREVAT(
	vcf.file = opt$vcf,
	vcf.file.genome = opt$genome,
	config.file = opt$config,
	df.ref.mut.sigs = GetPCAWGMutSigs(),
	target.mut.sigs = GetPCAWGMutSigsNames(),
	sequencing.artifact.mut.sigs = PCAWG.All.Sequencing.Artifact.Signatures,
	output.dir = opt$outdir,
	objective.fn = Default.Obj.Fn,
	num.cores = opt$`num-cores`,
	ga.pop.size = 100,
	ga.max.iter = 5,
	ga.run = 5,
	perform.strand.bias.analysis = TRUE,
	ref.forward.strand.var = "TumorDPRefForward",
	ref.reverse.strand.var = "TumorDPRefReverse",
	alt.forward.strand.var = "TumorDPAltForward",
	alt.reverse.strand.var = "TumorDPAltReverse",
	annotate = FALSE,
	write.vcf = FALSE,
	report = FALSE,
	save.rdata = FALSE,
	save.tsv = TRUE,
	save.scores = TRUE,
	report.format = "html",
	verbose = TRUE
)

message("FIREVAT complete. Results written to: ", opt$outdir)
