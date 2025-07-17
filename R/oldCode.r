#Simulate mortality data for fitting
mortData<-expand.grid(Year=1990:2018,
                      Species=c("SWO","BUM"),
                      Sample=1:8)
mortData<-mortData %>%
  mutate(logitS=rnorm(nrow(mortData),1+(Year-1990)/10+(Species=="BUM")*0.4))%>%
  mutate(S=exp(logitS)/(1+exp(logitS)))%>%
  mutate(alive=rbinom(nrow(mortData),1,S))
ggplot(mortData,aes(x=Year,y=alive,color=Species,fill=Species))+
  geom_jitter(width=0.1,height=0.05,alpha=0.5)+
  stat_smooth()
write_csv(mortData,"data/SimulatedMortality.csv")         

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
