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
p <- add_argument(p, '--plot', flag = TRUE, help = 'Generate plots')

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

if (argv$plot) {
    library(ggplot2)
    out.file1 = paste(as.character(argv$in_path), '_neutralitytestr.png', sep = '')
    gout1 <- plot_all(s)
    ggsave(out.file1,
        plot = gout1,
        devic = 'png',
        width = 12,
        height = 4,
        units = 'in',
        dpi = 300,
        bg = 'white'
    )
}

# Run mobster test
library(mobster)

fit = mobster_fit(
    data,
    K = 1
)
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

if (argv$plot) {
    out.file2 = paste(as.character(argv$in_path), '_mobster.png', sep = '')
    gout2 <- plot(fit$best)
    ggsave(out.file2,
        plot = gout2,
        devic = 'png',
        width = 12,
        height = 4,
        units = 'in',
        dpi = 300,
        bg = 'white'
    )
}



# Safe both outputs
sink(argv$out_path)
cat('MOBSTER Population Genetics statistics:\n')
print(data.frame(evo))
cat('\n\n')
print(summary(s))