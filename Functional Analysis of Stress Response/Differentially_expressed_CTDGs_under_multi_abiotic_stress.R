library(ggplot2)
library(ComplexUpset)
library(patchwork)
setwd("~/Desktop/project/CTDG_of_207_spe/draft_fig/")
colors <-c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a",
             "#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a",
             "#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a")

matrix_1 <- read.csv("~/R  project/207spe_ctdg/1.expression/16.Aty/expression_matrix_0.8.csv",
                                   sep = ',')
matrix_1$size <- factor(matrix_1$size, levels = c("member>10", "6<member<10", "3<member<5", "member=2"))
head(matrix_1,2)
treatment <- colnames(matrix_1)[3:7]

p <- upset(matrix_1,treatment,
      name = "treatment",
      width_ratio = 0.2,height_ratio = 0.5,
      base_annotations = list(
        "intersection size"=intersection_size(
          text = list(vjust=1,hjust=.5,angle=45),
          text_colors = c(on_background = "black",on_bar = "black"),
          mapping = aes(fill=size)
          )+
      scale_fill_manual(values=rev(c("#8dd3c7", "#ffffb3", "#bebada", "#fb8072")))+
      annotate(
            geom='text', x=Inf, y=Inf,
            label=paste('Total:', nrow(matrix_1)),
            vjust=1, hjust=1
          )
        ),
      min_degree=0);p

p1 <- upset(matrix_1,treatment,
            name = "treatment",
            width_ratio = 0.3,height_ratio = 0.5,
            base_annotations = list(
              "intersection size"=intersection_size(
                text = list(vjust=1,hjust=.5,angle=45),
                text_colors = c(on_background = "black",on_bar = "black"),
                mapping = aes(fill=size)
              )+
                scale_fill_manual(values=rev(c("#8dd3c7", "#ffffb3", "#bebada", "#fb8072")))+
                annotate(
                  geom='text', x=Inf, y=Inf,
                  label=paste('Total:', 3),
                  vjust=1, hjust=1
                )
            ),
            min_degree=4);p1

p2 <- upset(matrix_1,treatment,
            name = "treatment",
            width_ratio = 0.3,height_ratio = 0.5,
            base_annotations = list(
              "intersection size"=intersection_size(
                text = list(vjust=1,hjust=.5,angle=45),
                text_colors = c(on_background = "black",on_bar = "black"),
                mapping = aes(fill=size)
              )+
                scale_fill_manual(values=rev(c("#8dd3c7", "#ffffb3", "#bebada", "#fb8072")))+
                annotate(
                  geom='text', x=Inf, y=Inf,
                  label=paste('Total:', 41),
                  vjust=1, hjust=1
                )
            ),
            min_degree=3,max_degree=3);p2

p3 <- upset(matrix_1,treatment,
            name = "treatment",
            width_ratio = 0.3,height_ratio = 0.5,
            base_annotations = list(
              "intersection size"=intersection_size(
                text = list(vjust=1,hjust=.5,angle=45),
                text_colors = c(on_background = "black",on_bar = "black"),
                mapping = aes(fill=size)
              )+
                scale_fill_manual(values=rev(c("#8dd3c7", "#ffffb3", "#bebada", "#fb8072")))+
                annotate(
                  geom='text', x=Inf, y=Inf,
                  label=paste('Total:', 119),
                  vjust=1, hjust=1
                )
            ),
            min_degree=2,max_degree=2);p3

p4 <- upset(matrix_1,treatment,
            name = "treatment",
            width_ratio = 0.3,height_ratio = 0.5,
            base_annotations = list(
              "intersection size"=intersection_size(
                text = list(vjust=1,hjust=.5,angle=45),
                text_colors = c(on_background = "black",on_bar = "black"),
                mapping = aes(fill=size)
              )+
                scale_fill_manual(values=rev(c("#8dd3c7", "#ffffb3", "#bebada", "#fb8072")))+
                annotate(
                  geom='text', x=Inf, y=Inf,
                  label=paste('Total:', 354),
                  vjust=1, hjust=1
                )
            ),
            min_degree=1,max_degree=1);p4

T <- wrap_plots(p4, p3, p2, p1, ncol = 2);T
ggsave("DE-CTDGs_8.26*4.pdf", 
       plot = p, width = 8.26, height = 4, units = "in")

ggsave("sub_DE-CTDGs_8.26*11.pdf", 
       plot = T, width = 8.26, height = 11, units = "in")













