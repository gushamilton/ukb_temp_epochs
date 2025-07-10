pacman::p_load(tidyverse, ncdf4, lubridate, fuzzyjoin)

nc_data <- nc_open("data/HadUK12_tmax_1981-2019_raw.nc")
x <- ncvar_get(nc_data, "projection_x_coordinate") # easting (1D)
y <- ncvar_get(nc_data, "projection_y_coordinate") # northing (1D)
t <- ncvar_get(nc_data, "date") # date (1D)
tasmax.array <- ncvar_get(nc_data, "tasmax") # this is the main temperature data (3D)
y

latitude <- read_tsv("data/UKB_assessmentcentres_coords.txt") %>%
    mutate(easting = 100*gridx, northing = 100*gridy)

included_arrays<- latitude %>%
  rowwise() %>%
  mutate(easting_closest = which.min(abs(x - easting)),
         northing_closest =which.min(abs(y - northing))) %>%
  mutate(easting_closest_actual = x[which.min(abs(x - easting))],
         northing_closest_actual =y[which.min(abs(y - northing))]) %>%
  ungroup()

unique_combination <- included_arrays %>%
  distinct(easting_closest, northing_closest)


extract_temp <- function(x,y) {
tibble(temp = tasmax.array[x[1], y[1],],
       easting_closest = x,
       northing_closest = y,
       day = 1:14244)
}

#check it works

extract_temp(unique_combination$easting_closest[1], unique_combination$northing_closest[1]) 


temp_on_day <-map2_dfr(unique_combination$easting_closest, unique_combination$northing_closest, extract_temp)



final <- included_arrays %>%
  select(centre_id, centre_name, easting_closest, northing_closest) %>%
  right_join(temp_on_day) %>%
  mutate(date = as.Date("1980-01-01") + (day - 1)) %>%
  select(centre_name, centre_id, temp, day, date)  

final %>%
    write_tsv("data/ukb_temps_per_date.tsv.gz")
