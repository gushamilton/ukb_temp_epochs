
# Load necessary libraries
if (!requireNamespace("ncdf4", quietly = TRUE)) {
  install.packages("ncdf4", repos = "https://cloud.r-project.org")
}
library(ncdf4)

# Define file paths
haduk_file <- "/Users/fh6520/R/temp_ukb/data/HadUK12_tmax_1981-2019_raw.nc"

# Open the NetCDF file
nc_data <- nc_open(haduk_file)

# Print the summary of the NetCDF file
print(nc_data)

# Close the NetCDF file
nc_close(nc_data)
