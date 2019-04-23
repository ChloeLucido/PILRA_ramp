library(ggplot2)
library(readr)

pilra <- read_tsv("random_permutations_data.tsv")
pilra$Changed = pilra$Changed * 100
pilra$Created = pilra$Created * 100
pilra$Destroyed = pilra$Destroyed * 100
pilra$Location = as.integer(pilra$Location - 10)

pilra <- gather(pilra,key = "Type of Change", value = Percent, -c(Location,Mutation))
ggplot(pilra, aes(x = Location, y = Percent, colour = Mutation)) +
  geom_line() + 
  facet_wrap(~`Type of Change`) +
  ylab("Percent of Ramp Sequences Changed") +
  ylim(0,20) +
  scale_x_continuous(name="Position in Gene by Quantile",breaks=seq(0,100,10)) +
  ggtitle("Frequency of Ramp Sequences Changed By Random Mutations") +
  theme_bw()

 ggsave("random_permutations.pdf")

