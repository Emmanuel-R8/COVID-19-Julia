

library(rjson)
library(tidyverse)
library(lubridate)

######################################################################################################
##
## Import from the Neherlab case_counts.json file
##

result <- fromJSON(file = "case_counts.json")

res <- tibble(.rows = 0)

for (i in 1:length(result)) {
  locationName <- names(result[i])
  print(locationName)

  r <- result[[i]]
  if (length(result[[i]]) > 0) {
    for (j in 1:length(result[[i]])) {
      t <- r[[j]]$time
      c <- r[[j]]$cases
      d <- r[[j]]$deaths

      if (!is.null(t) & !is.null(c) & !is.null(d)) {
        df <- tibble(
          location = locationName,
          time = ymd(r[[j]]$time),
          cases = r[[j]]$cases,
          deaths = r[[j]]$deaths
        )
        res <- rbind(res, df)
      }
    }
  }
}


# Remove all rows with no deaths
deaths <- res %>% filter(deaths > 0)

# Only keep locations with more that 10 data points
deaths10 <- deaths %>%
  group_by(location) %>%
  mutate(n = n()) %>%
  filter(n >= 10) %>%
  view()


write_csv(deaths, "deaths.csv")
write_csv(deaths10, "deaths10.csv")


######################################################################################################
##
## Import from the ECDC website
##

ecdc <- read.csv("ecdc.csv")
ecdc <- as_tibble(ecdc)

ecdc_sum <- ecdc %>% group_by(countriesAndTerritories) %>%
  mutate(deaths = cumsum(deaths),
         cases = cumsum(deaths)) %>%
  ungroup() %>%
  filter(deaths > 0)

ecdc_sum %>% group_by(countriesAndTerritories) %>%
  summarise(n = n(), max = max(deaths)) %>% arrange(max) %>% view()
