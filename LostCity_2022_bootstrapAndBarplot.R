RifR_LC <- read.csv("RifampicinResistanceAssay2022.csv")

library(boot)

# creates the RifR_LC.median function
RifR_LC.median <- function(RifR_LC.column, indices) {
  return(median(RifR_LC.column[indices], na.rm = "TRUE"))
} # creates the median function for bootstrapping 

# creates the RifR_LC.median function
RifR_LC.foldChange <- function(RifR_LC.frame, indices) {
  a <- median(RifR_LC.frame[indices,1], na.rm = "TRUE")
  b <- median(RifR_LC.frame[indices,2], na.rm = "TRUE")
  return(b/a)
}

# creates the RifR_LC.summary data frame with column names and column binding
RifR_LC.summary <- colnames(RifR_LC)
RifR_LC.summary <- cbind(RifR_LC.summary, 
                         data.frame("median" = numeric(length(colnames(RifR_LC))), 
                                    "ci.low" = numeric(length(colnames(RifR_LC))), 
                                    "ci.high" = numeric(length(colnames(RifR_LC))),
                                    "Fold.Change.from.EcMutY" = numeric(length(colnames(RifR_LC))),
                                    "ci.low.Fold.Change.from.EcMutY" = numeric(length(colnames(RifR_LC))),
                                    "ci.high.Fold.Change.from.EcMutY" = numeric(length(colnames(RifR_LC))),
                                    "n" = numeric(length(colnames(RifR_LC)))))

#  for-loop to walk through each set of RifR_LC frequencies for each of the different MutY variants
MutY.list <- colnames(RifR_LC)
index <- 1
for (MutY in MutY.list) {
  RifR_LC.summary[index,"median"] <- median(RifR_LC[,MutY], na.rm = "TRUE")
  RifR_LC.bootstrap = boot(RifR_LC[,MutY], RifR_LC.median, R=10000)
  RifR_LC.summary[index,"ci.low"] <- boot.ci(RifR_LC.bootstrap, type="perc")$percent[[4]]
  RifR_LC.summary[index,"ci.high"] <- boot.ci(RifR_LC.bootstrap, type="perc")$percent[[5]]
  index <- index + 1
}

# for-loop to walk through each set of RifR_LC frequencies for each LCHF MutY 
#variants calculate fold change
k <- 3
for (MutY in MutY.list[3:11]) {
  dataFrame <- data.frame("EcMutY" = RifR_LC[,2], RifR_LC[,MutY])
  RifR_LC.summary[k,"Fold.Change.from.EcMutY"] <- RifR_LC.foldChange(dataFrame)
  RifR_LC.bootstrap = boot(dataFrame, RifR_LC.foldChange, R=10000)
  RifR_LC.summary[k,"ci.low.Fold.Change.from.EcMutY"] <- boot.ci(RifR_LC.bootstrap, type="perc")$percent[[4]]
  RifR_LC.summary[k,"ci.high.Fold.Change.from.EcMutY"] <- boot.ci(RifR_LC.bootstrap, type="perc")$percent[[5]]
  i <- 0
  for (num in RifR_LC[,MutY]) {
    if (!is.na(num)) {
      i <- i + 1
    }
  }
  RifR_LC.summary[k, "n"] <- i
  k <- k + 1
}

# add fold change value for null MutY
dataFrame <- data.frame("EcMutY" = RifR_LC[,2], RifR_LC[,"MutY0"])
RifR_LC.summary[1,"Fold.Change.from.EcMutY"] <- RifR_LC.foldChange(dataFrame)
RifR_LC.bootstrap = boot(dataFrame, RifR_LC.foldChange, R=10000)
RifR_LC.summary[1,"ci.low.Fold.Change.from.EcMutY"] <- boot.ci(RifR_LC.bootstrap, type="perc")$percent[[4]]
RifR_LC.summary[1,"ci.high.Fold.Change.from.EcMutY"] <- boot.ci(RifR_LC.bootstrap, type="perc")$percent[[5]]
i <- 0
for (num in RifR_LC[,2]) {
  if (!is.na(num)) {
    i <- i + 1
  }
}
RifR_LC.summary[1, "n"] <- i

# add fold change value for EcMuty
dataFrame <- data.frame("EcMutY" = RifR_LC[,2], RifR_LC[,2])
RifR_LC.summary[2,"Fold.Change.from.EcMutY"] <- RifR_LC.foldChange(dataFrame)
RifR_LC.summary[2,"ci.low.Fold.Change.from.EcMutY"] <- as.numeric("NA")
RifR_LC.summary[2,"ci.high.Fold.Change.from.EcMutY"] <- as.numeric("NA")
i <- 0
for (num in RifR_LC[,1]) {
  if (!is.na(num)) {
    i <- i + 1
  }
}
RifR_LC.summary[2, "n"] <- i

# Creates function for plotting barplot
rotate_x <- function(data, column_to_plot, labels_vec, rot_angle) {
  plt <- barplot(data[[column_to_plot]],
                 main = "Rifampicin Resistance Assay",
                 #ylab = expression(Rif^R ~ Colony ~ Count),
                 ylim = c(0,200),
                 col = c("snow4", "snow2",
                         "coral","coral","coral",
                         "seagreen","seagreen","seagreen",
                         "skyblue","skyblue","skyblue"),
                 density = c(1000, 1000,
                             1000,50,50,
                             1000,50,50,
                             1000,50,50),
                 xaxt="n")
  title(ylab = expression(Rif^R ~ Colony ~ Count), line=2.5, cex.lab=1)
  text(plt, par("usr")[3], labels = labels_vec, 
       srt = rot_angle, adj = c(1,1.1), xpd = TRUE, cex=0.7) 
}

# Creates vector of labels for x-axis of barplot
v <- c("null", "EcMutY", "Ms MutY", "Ms MutY catalysis-",
       "Ms MutY recognition-", "Rb MutY", "Rb MutY catalysis-",
       "Rb MutY recognition-", "Tt MutY", "Tt MutY catalysis-",
       "Tt MutY recognition-")

# Creates barplot
rotate_x(RifR_LC.summary, 'median',v, 45)

# adds error bars to barplot based on confidence intervals.
i = 0.7
j = 1
for (ci in RifR_LC.summary$median) {
  segments(i, RifR_LC.summary[j,3], i, RifR_LC.summary[j,4])
  segments(i-.2, RifR_LC.summary[j,3], i+.2, RifR_LC.summary[j,3])
  segments(i-.2, RifR_LC.summary[j,4], i+.2, RifR_LC.summary[j,4])
  i <- i + 1.2
  j <- j + 1
}
text(x = c(1.3, 1.9, 3.6, 5.4), y = c(156, 170, 184, 198), "*")
segments(x0=0.7, y0 = 151 , x1 = 1.9, y1 = 151 )
segments(x0=0.7, y0 = 149 , x1 = 0.7, y1 = 153 )
segments(x0=1.9, y0 = 149 , x1 = 1.9, y1 = 153 )

segments(x0=0.7, y0 = 164 , x1 = 3.1, y1 = 164 )
segments(x0=0.7, y0 = 162 , x1 = 0.7, y1 = 166 )
segments(x0=3.1, y0 = 162 , x1 = 3.1, y1 = 166 )

segments(x0=0.7, y0 = 178 , x1 = 6.7, y1 = 178 )
segments(x0=0.7, y0 = 176 , x1 = 0.7, y1 = 180 )
segments(x0=6.7, y0 = 176 , x1 = 6.7, y1 = 180 )

segments(x0=0.7, y0 = 192 , x1 = 10.2, y1 = 192 )
segments(x0=0.7, y0 = 190 , x1 = 0.7, y1 = 194 )
segments(x0=10.2, y0 = 190 , x1 = 10.2, y1 = 194 )

# writes out the RifR_LC.summary table as a CSV 
write.csv(RifR_LC.summary, "AllExperiments2022data_5.csv")
