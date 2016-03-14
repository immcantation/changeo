# Imports
library(alakazam)
library(shazam)
library(ggplot2)

# Read in database file
db <- readChangeoDb('S43_db-pass_parse-select.tab')
# Calculate distance to nearest neighbor
db <- distToNearest(db, model='hs5f',symmetry='min')
# Plot resulting histogram with vertical threshold line
p1 <- ggplot() + theme_bw() + 
    ggtitle("Distance to nearest: hs5f") + xlab("distance") +
    geom_histogram(data=db, aes(x=DIST_NEAREST), binwidth=0.025, 
                   fill="steelblue", color="white") + 
    geom_vline(xintercept=0.165, linetype=2, color='red')
plot(p1)
