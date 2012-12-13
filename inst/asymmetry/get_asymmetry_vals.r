setwd("/Users/ScottMac/Dropbox/CANPOLIN_treeshape_ms/data/empirical/rezendeetal2007nature/networks/matrices/binary")
path <- "/Users/ScottMac/Dropbox/CANPOLIN_treeshape_ms/data/empirical/rezendeetal2007nature/networks/matrices/binary/"
files<-dir()
files_path <- paste0(path, files)
datlist <- lapply(files_path, read.csv)

getasymm <- function(x) {
	plants <- names(x)[!names(x) %in% "X"]
	animals <- x[,1]
	c(p=length(plants), a=length(animals), 
		asymm = round((length(animals)-length(plants)) / (length(animals)+length(plants)), 2),
		ratio = round((length(animals)/length(plants)), 2) )
}

asymvalues <- ldply(datlist, getasymm)
summary(asymvalues)

makevec <- function(x, asymm){ c(start = x, plant = x, poll = x*asymm, tot = x + x*asymm) }
n <- c(10, 15, 20, 30, 40, 50, 60, 70, 80, 100, 110, 130, 150, 175, 200)

# These are the levels 
(things <- sapply(n, makevec, asymm=2.47))
things2 <- data.frame(t(things))
things2$asymm_calc <- things2$poll/things2$plant
things2