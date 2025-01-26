library(ca)
library(devtools)
library(CAinterprTools)
library(cabootcrs)
library(readr)
library(ggplot2)
library(ggrepel)
library(paletteer)
library(extrafont)


#colour and text style
scale_colour_paletteer_d("MetBrewer::Greek")
scale_color_paletteer_d("MetBrewer::Greek")
scale_fill_paletteer_d("MetBrewer::Greek")
paletteer_d("MetBrewer::Greek")

font_import()
loadfonts(device = "win", quiet = TRUE)
choose_font(c("Consolas", "Cascadia-Code", "mono"), quiet = TRUE)

#read data as df
cadatawc <- read_delim("jang_data/CA Wachow cleaned+filled+noise removed.csv", 
                                                     delim = ";", escape_double = FALSE, trim_ws = TRUE)
#transform to matrix with first row as rownames
matrix.please <- function(cadatawc) {
  matrcadatawc <- as.matrix(cadatawc[,-1])
  rownames(matrcadatawc) <- cadatawc$Gräber
  matrcadatawc
}
matrcadatawc <- matrix.please(cadatawc)

#CA basic 
ca.fit <- ca(matrcadatawc)


#CA get variances
ca.sum <- summary(ca.fit)
dim.var.percs <- ca.sum$scree[,"values2"]
dim.var.percs

#CA plot data prep
ca.plot <- plot(ca.fit)

make.ca.plot.df <- function (ca.plot.obj,
                             row.lab = "Rows",
                             col.lab = "Columns") {
  df <- data.frame(Label = c(rownames(ca.plot.obj$rows),
                             rownames(ca.plot.obj$cols)),
                   Dim1 = c(ca.plot.obj$rows[,1], ca.plot.obj$cols[,1]),
                   Dim2 = c(ca.plot.obj$rows[,2], ca.plot.obj$cols[,2]),
                   Variable = c(rep(row.lab, nrow(ca.plot.obj$rows)),
                                rep(col.lab, nrow(ca.plot.obj$cols))))
  rownames(df) <- 1:nrow(df)
  df
}
ca.plot.df <- make.ca.plot.df(ca.plot,
                              row.lab = "Graves",
                              col.lab = "Types")
ca.plot.df$Size <- ifelse(ca.plot.df$Variable == "Graves", 2, 1)

#CA plot ggplot
p <- ggplot(ca.plot.df, aes(x = Dim1, y = Dim2,
                            col = Variable, shape = Variable,
                            label = Label, size = Size)) +
  geom_vline(xintercept = 0, lty = "dashed", alpha = .5) +
  geom_hline(yintercept = 0, lty = "dashed", alpha = .5) +
  geom_point()

#sizing and repel text
p <- p +
  scale_x_continuous(limits = range(ca.plot.df$Dim1) + c(diff(range(ca.plot.df$Dim1)) * -0.2,
                                                         diff(range(ca.plot.df$Dim1)) * 0.2)) +
  scale_y_continuous(limits = range(ca.plot.df$Dim2) + c(diff(range(ca.plot.df$Dim2)) * -0.2,
                                                         diff(range(ca.plot.df$Dim2)) * 0.2)) +
  scale_size(range = c(3, 3), guide = F) + 
  geom_text_repel(max.iter = 1000, max.overlaps = Inf, show.legend = F)

#graph labels, colours, and upside down flip plus theme
p <- p +
  labs(title ="Correspondence Analysis, Cemetery: Wachow, Kr. Havelland", caption = "based on Stetzuhn 2023 Fig. 10, modifications in style only" ,x = paste0("Dimension 1 (", signif(dim.var.percs[1], 3), "%)"),
       y = paste0("Dimension 2 (", signif(dim.var.percs[2], 3), "%)"),
       col = "", shape = "") +
  scale_y_reverse() +
  scale_x_reverse() +
  scale_colour_paletteer_d("MetBrewer::Greek") +
  theme_minimal()
plot(p)

#values
malinvaud(data=matrcadatawc)
rows.cntr(data=matrcadatawc)
cols.cntr(data=matrcadatawc)
aver.rule(matrcadatawc)
