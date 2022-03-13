library(ggplot2)
library(glue)
library(WriteXLS)
library(reshape2)
library(xtable)

tools <- c("MOFA",
           "schema",
           "WNN",
           "scAI",
           "liger",
           "MOFARed",
           "DIABLORed",
           "symphonyRed",
           "MOJITOO")

source("../util/plot_util.R")


df.list <- lapply(tools, function(x) readRDS(glue("save/{x}_RAM_TIME.RDS")))
dff <- do.call(rbind, df.list)
dff$num <- as.numeric(dff$num)
head(dff)

pdf("viz/time_mem.pdf", width=5, height=4)

p1 <- ggplot(dff, aes(x=num, y=log(Elapsed_Time_sec), color=tool)) +
  geom_point() +
  geom_line() +
  ylab("log elapsed time seconds") +
  xlab("Cells") +
  theme_minimal() +
  scale_x_continuous(breaks=seq(3000, 30000, 3000))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_color_manual(values=color)

print(p1)
print(p1+theme(legend.position = "none"))

p2 <- ggplot(dff, aes(x=num, y=Peak_RAM_Used_MiB/1024.0, color=tool)) +
  geom_point() +
  geom_line() +
  ylab("Peak_Memory GB") +
  xlab("Cells") +
  theme_minimal() +
  scale_x_continuous(breaks=seq(3000, 30000, 3000))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_color_manual(values=color)

print(p2)
print(p2+theme(legend.position = "none"))

dev.off()

#system("cd viz/; rm *.svg; pdf2svgs time_mem.pdf time_mem.svg")


tdf <- dcast(dff[, c("Elapsed_Time_sec", "tool", "num")], num~tool, value.var="Elapsed_Time_sec")
rownames(tdf) <- tdf$num
tdf$num <- NULL
tdf <- round(tdf/60.0, 2)

mdf <- dcast(dff[, c("Peak_RAM_Used_MiB", "tool", "num")], num~tool, value.var="Peak_RAM_Used_MiB")
rownames(mdf) <- mdf$num
mdf$num <- NULL
mdf <- round(mdf/1024.0, 2)

df_list <- list(
  "time_minutes" = tdf,
  "PeakMem_Giga" = mdf
)

WriteXLS(df_list, "viz/time_mem.xlsx", SheetNames = names(df_list), row.names=T)


print(xtable(tdf))
print(xtable(mdf))

