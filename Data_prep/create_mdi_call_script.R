

#!/usr/bin/env Rscript

library(optparse)

# Load magrittr for the pipe %>%
library(magrittr)


input_arguments <- function() {
  option_list <- list(

    # Directory to read from
    optparse::make_option(c("-d", "--dir"),
      type = "character",
      default = ".",
      help = "directory to read files from [default= %default]",
      metavar = "character"
    ),

    # File extension of interest
    optparse::make_option(c("-e", "--ext"),
      type = "character",
      default = ".csv",
      help = "file extension to use for reading in files [default= %default]",
      metavar = "character"
    ),

    # Number of iterations for MDI
    optparse::make_option(c("-n", "--n_iter"),
      type = "integer",
      default = 1000,
      help = "number of iterations to run MDI for [default= %default]",
      metavar = "integer"
    ),

    # Thinning factor
    optparse::make_option(c("-t", "--thin"),
      type = "integer",
      default = 0,
      help = "number of iterations to thin MDI for [default= %default]",
      metavar = "integer"
    ),

    # Random seed for MDI
    optparse::make_option(c("-s", "--seed"),
      type = "integer",
      default = 1,
      help = "random seed for MDI [default= %default]",
      metavar = "integer"
    ),

    # Output file path
    optparse::make_option(c("-o", "--output"),
      type = "character",
      default = "output.csv",
      help = "file to feed MDI output into [default= %default]",
      metavar = "character"
    ),

    # File path to MDI
    optparse::make_option(c("--mdi_home"),
      type = "character",
      default = "mdipp-1.0.1",
      help = "path to call mdi from [default= %default]",
      metavar = "character"
    ),

    # Type of model to fit
    optparse::make_option(c("--type"),
      type = "character",
      default = "N",
      help = "type of model to fit (see MDI documenation) [default= %default]",
      metavar = "character"
    )
  )

  opt_parser <- optparse::OptionParser(option_list = option_list)
  opt <- optparse::parse_args(opt_parser)
}

# Input arguments
args <- input_arguments()

# Read in files
files_present <- list.files(path = args$dir, pattern = args$ext)

# Create parts of MDI code
shebang_code <- "#!/usr/bin/env bash"

# Call MDI itself
generic_call <- paste(args$mdi_home, "mdipp", sep = "/") # generic_call <-  "mdipp-1.0.1/mdipp"

# Arguments for MDI
tail <- paste(
  "-n", args$n_iter,
  "-t", args$thin,
  "-s", args$seed,
  ">", args$output
)

# Files to use
mdi_files <- files_present %>%
  paste(args$type, .) %>%
  paste(., collapse = " ")

# Put the main script together so we don't have a space at the beginning of the line
main_call <- paste(generic_call, mdi_files, tail)

# Output statement
cat(paste0(shebang_code, "\n", main_call))
