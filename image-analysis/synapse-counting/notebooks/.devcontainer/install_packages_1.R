# ---- Description
# This file installs the required R packages on the Docker image during the build time.
# The list of packages and their versions is set on the packages.json file
# Added support for Bioconductor packages

# ---- Dependencies
# To parse the json jq must be installed on the docker image.
# See: https://stedolan.github.io/jq/

# ---- Code starts here
# Set the working directory
setwd("./settings/")

# Set CRAN mirror explicitly
local({
  r <- getOption("repos")
  r["CRAN"] <- "https://cloud.r-project.org/"
  options(repos = r)
})

# Required packages
install.packages("remotes", repos = "https://cloud.r-project.org/")
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org/")
}

# Explicitly set Bioconductor repositories
BiocManager::repositories()

# Set the jq query
jq_command <- 'jq -r ".packages[] |  [.package, .version] | @tsv" packages.json'
# Debug mode
# jq_command <- 'jq -r ".package_dev[] |  [.package, .version] | @tsv" packages.json'

# Parse the json file with the list of package
raw <- system(command = jq_command, intern = TRUE)

package_list <- lapply(raw, function(i){
  x <- unlist(strsplit(x = i, split = "\t"))
  data.frame(package = x[1], version = x[2], stringsAsFactors = FALSE)
})

# Transform the list into a data.frame
packages_df <- as.data.frame(t(matrix(unlist(package_list), nrow = 2)),
                              stringsAsFactors = FALSE)
names(packages_df) <- c("package", "version")
packages_df$success <- FALSE

# Bioconductor package list
bioc_packages <- tryCatch({
  available_packages <- available.packages(repos = BiocManager::repositories())
  available_packages[,1]
}, error = function(e) {
  warning("Could not retrieve Bioconductor package list")
  character(0)
})

# ---- Install the packages
for(i in 1:nrow(packages_df)){
  # Check if package is from Bioconductor
  is_bioconductor <- packages_df$package[i] %in% bioc_packages
  
  if(is_bioconductor){
    cat("\033[0;92m", paste("Installing Bioconductor package", packages_df$package[i]), "\033[0m\n", sep = "")
    tryCatch({
      BiocManager::install(packages_df$package[i], ask = FALSE, update = FALSE)
    }, error = function(e) {
      warning(paste("Failed to install Bioconductor package", packages_df$package[i]))
    })
  } else {
    if(!packages_df$package[i] %in% rownames(installed.packages()) ||
       (packages_df$package[i] %in% rownames(installed.packages()) &&
        packageVersion(packages_df$package[i]) != packages_df$version[i])){
      cat("\033[0;92m", paste("Installing CRAN package", packages_df$package[i]), "\033[0m\n", sep = "")
      tryCatch({
        remotes::install_version(package = packages_df$package[i],
                                  version = packages_df$version[i],
                                  dependencies = c("Depends", "Imports"),
                                  upgrade = FALSE,
                                  verbose = FALSE,
                                  quiet = FALSE,
                                  repos = "https://cloud.r-project.org/")
      }, error = function(e) {
        install.packages(packages_df$package[i], 
                         repos = "https://cloud.r-project.org/")
      })
    }
  }
  
  # Check installation success
  if(!packages_df$package[i] %in% rownames(installed.packages()) ||
     (packages_df$package[i] %in% rownames(installed.packages()) &&
      packageVersion(packages_df$package[i]) != packages_df$version[i])){
    packages_df$success[i] <- FALSE
  } else {
    packages_df$success[i] <- TRUE
  }
}

# Report installation results
for(i in 1:nrow(packages_df)){
  if(packages_df$success[i]){
    cat("\033[0;92m", packages_df$package[i], "...OK", "\033[0m\n")
  } else {
    cat("\033[0;91m", packages_df$package[i], "...Fail", "\033[0m\n")
  }
}

# Final check
for(i in sort(c(packages_df$package))){
  if(!i %in% rownames(installed.packages())){
    cat(i, "...Failed\n")
  } else {
    cat(i, "...OK\n")
  }
}

if(!any(packages_df$success)){
  stop("One or more packages are missing...")
}