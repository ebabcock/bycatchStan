#Function to calculate total bycatch from one stan model object
getBycatch <- function(mod1, logdat) {
  gg1 <- extract(mod1, pars = "StrataBycatch")$StrataBycatch
  gg1 <- as.data.frame(gg1) %>%
    rename(row = Var2)
  gg1$Year <- logdat$Year[gg1$row]
  modelyrSum1 <- gg1 %>% group_by(Year, iterations) %>%
    summarize(yearsum = sum(value)) %>%
    group_by(Year) %>%
    summarize(
      mean = mean(yearsum),
      se = sd(yearsum),
      lower = quantile(yearsum, 0.025),
      upper = quantile(yearsum, 0.975),
      median = quantile(yearsum, 0.5)
    )
  modelyrSum1
}

#Function to get mortality predictions 
getMortPred <- function(logdat,
                        mortPredDat,
                        sp,
                        mortMod,
                        bycatchMod,
                        codeName) {
  ggb <- extract(bycatchMod, pars = "StrataBycatch")$StrataBycatch %>%
    reshape2::melt() %>%
    as.data.frame() %>%
    rename(row = Var2, StrataBycatch = value) %>%
    mutate(Year = logdat$Year[row])
  ggm <- extract(mortMod, pars = "strataProb")$strataProb %>%
    reshape2::melt() %>%
    as.data.frame() %>%
    rename(row = Var2, strataProb = value) %>%
    mutate(Year = mortPredDat$Year[row], Species = mortPredDat$Species[row]) %>%
    filter(Species == sp) %>%
    mutate(row = row - min(row) + 1)
  gg1 <- left_join(ggb, ggm, by = c("Year", "row", "iterations"))
  modelyrSum1 <- gg1 %>% group_by(Year, iterations) %>%
    summarize(
      Bycatch = sum(StrataBycatch),
      Mortality = sum(StrataBycatch * strataProb)
    ) %>%
    ungroup() %>%
    pivot_longer(Bycatch:Mortality,
                 names_to = "Outcome",
                 values_to = "Total") %>%
    group_by(Year, Outcome) %>%
    summarize(
      mean = mean(Total),
      se = sd(Total),
      lower = quantile(Total, 0.025),
      upper = quantile(Total, 0.975),
      median = quantile(Total, 0.5)
    )
  modelyrSum1
}

