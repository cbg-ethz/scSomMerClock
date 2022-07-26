#!/usr/bin/env Rscript

library(argparser)

p <- arg_parser('Run neutralitytestr from the Williams et al. 2016 paper.')
p <- add_argument(p, 'in_path', help = 'Input vcf file (str)')
p <- add_argument(p, 'out_path', help = 'output log file (str)')
p <- add_argument(p, '--fmin', default = 0.1, help = 'Min. VAF freq. (float)')
p <- add_argument(p, '--fmax', default = -1, help = 'Max. VAF freq. (float)')
p <- add_argument(p, '--depth', default = NULL, help = 'Seq. depth (float)')
p <- add_argument(p, '--cellularity', default = 1,
    help = 'Sample cellularity (float)')
p <- add_argument(p, '--ploidy', default = 2, help = 'Sample ploidy (float)')
p <- add_argument(p, '--K', default = 2, help = 'Cluster K (float)')
p <- add_argument(p, '--maxIter', default = 500, help = 'Mobster maxIter (float)')
p <- add_argument(p, '--bootstrap', flag = FALSE, help = 'Run bootstrapping')
p <- add_argument(p, '--plot', flag = FALSE, help = 'Generate plots')
p <- add_argument(p, '--stdout', flag = FALSE, help = 'Print summary to stdout')

argv <- parse_args(p)

data <- read.csv(argv$in_path, sep = "\t")

# Run old Williams et al. 2016 frequentist test
library(neutralitytestr)
if (argv$fmax > argv$fmin) {
    s <- neutralitytest(
        data$VAF,
        fmin = argv$fmin,
        fmax = as.numeric(argv$fmax),
    )
} else {
    s <- neutralitytest(
        data$VAF,
        read_depth = as.numeric(argv$depth),
        cellularity = argv$cellularity,
        ploidy = argv$ploidy,
    )
}

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

get.mobster.fit <- function(data, maxIter, maxK) {
    tryCatch(
        expr = {
            return(
                mobster_fit(
                    data, maxIter = maxIter, K = 1:maxK,
                )
            )
        },
        error = function(e) {
            return( FALSE )
        }
    )
}

fit <- get.mobster.fit(data, argv$maxIter, argv$K)

if (!is.logical(fit)) {
    if (argv$bootstrap) {
        bootstrap_results = mobster_bootstrap(
          fit$best,
          n.resamples = 20,
          auto_setup = 'FAST' # forwarded to mobster_fit
        )

        bootstrap_statistics = bootstrapped_statistics(
          fit$best,
          bootstrap_results = bootstrap_results
        )
    }

    mobster.evo <- function(fit){
        tryCatch(
            expr = {
                return( data.frame(evolutionary_parameters(fit)) )
            },
            error = function(e) {
                return( 'MOBSTER: no evolutionary parameters' )
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
}

# Safe both outputs
if (!argv$stdout) {
    sink(argv$out_path)
}

cat('MOBSTER Population Genetics statistics:\n')
if (!is.logical(fit)) {
    print(fit$best)
    print(mobster.evo(fit))
    if (argv$bootstrap) {
        cat('\nMOBSTER Model bootstrap:\n')
        print(data.frame(bootstrap_statistics$bootstrap_model))
    }
} else {
    cat('\nERROR\n')
}
cat('\n\n')
print(summary(s))

if (!argv$stdout) {
    sink()
}