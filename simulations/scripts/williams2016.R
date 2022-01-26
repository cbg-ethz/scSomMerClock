#!/usr/bin/env Rscript

library(argparser)

p <- arg_parser('Run neutralitytestr from the Williams et al. 2016 paper.')
p <- add_argument(p, 'in_path', help = 'Input vcf file (str)')
p <- add_argument(p, 'out_path', help = 'output log file (str)')
p <- add_argument(p, '--depth', default = NULL,
    help = 'Sequencing depth (float)')
p <- add_argument(p, '--fmin', default = 0.01,
    help = 'Minimum VAF frequency (float)')
p <- add_argument(p, '--fmax', default = 0.25, help =
    'Maximum VAF frequency (float)')
p <- add_argument(p, '--cellularity', default = 1,
    help = 'Cellularity of sample (float)')
argv <- parse_args(p)

data <- read.csv(argv$in_path, sep = "\t")

# Run old Williams et al. 2016 frequentist test
library(neutralitytestr)
s <- neutralitytest(
    data$VAF,
    read_depth = as.numeric(argv$depth),
    fmin = argv$fmin,
    fmax = argv$fmax,
    cellularity = argv$cellularity,
)

# Run mobster test
library(mobster)

fit = mobster_fit(data)
evo <- function(fit){
    tryCatch(
        expr = {
            return(evolutionary_parameters(fit))
        },
        error = function(e){
            return(NA)
        }
    )
}

# Safe both outputs
sink(argv$out_path)
cat('MOBSTER Population Genetics statistics:\n')
print(data.frame(evo))
cat('\n\n')
print(summary(s))