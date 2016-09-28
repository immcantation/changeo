# Imports
library(alakazam)
library(shazam)
library(ggplot2)

# Read in database file
db <- readChangeoDb("S43_db-pass_parse-select.tab")
# Calculate distance to nearest neighbor
db <- distToNearest(db, model="ham", symmetry="min")
# Find minimum between two modes of distribution
threshold <- findThreshold(db$DIST_NEAREST)
# Plot resulting histogram with vertical threshold line
p1 <- ggplot() + theme_bw() + 
    ggtitle("Distance to nearest: ham") + xlab("distance") +
    geom_histogram(data=db, aes(x=DIST_NEAREST), binwidth=0.025, 
                   fill="steelblue", color="white") + 
    geom_vline(xintercept=threshold, linetype=2, color="red")
plot(p1)
