# library(fishualize)
# library(tidyverse)
# library(cowplot)
# # select one species per family
# shapes<-fishapes() %>% group_by(family) %>% summarise(option=first(option))

# for(i in 1:length(shapes$option)){
# g1<-ggplot() + add_fishape(family = shapes$family[i], option = shapes$option[i]) +
#       theme_void()
# assign(paste('g1', i, sep='_'), g1)
# }

# pdf(file = 'fishualize_community.pdf', height=7, width=12)
  
# dev.off()