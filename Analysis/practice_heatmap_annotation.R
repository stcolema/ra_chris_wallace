library(pheatmap)
# Generate some data
test = matrix(rnorm(200), 20, 10)
test[1:10, seq(1, 10, 2)] = test[1:10, seq(1, 10, 2)] + 3
test[11:20, seq(2, 10, 2)] = test[11:20, seq(2, 10, 2)] + 2
test[15:20, seq(2, 10, 2)] = test[15:20, seq(2, 10, 2)] + 4
colnames(test) = paste("Test", 1:10, sep = "")
rownames(test) = paste("Gene", 1:20, sep = "")
# original figure
pheatmap(test)

# Add annotation as described above, and change the name of annotation
annotation <- data.frame(Var1 = factor(1:ncol(test) %% 2 == 0, labels = c("Exp1", "Exp2")))
annotation_col <- data.frame(Var2 = factor(1:nrow(test) %% 2 == 0, labels = c("Exp3", "Exp4")),
                             Var3 = factor(1:nrow(test) %% 3 == 0, labels = c("Exp3", "Exp4"))
                             )

rownames(annotation) <- colnames(test) # check out the row names of annotation
rownames(annotation_col) <- rownames(test) # check out the row names of annotation

pheatmap(test, annotation_col = annotation)
pheatmap(test, annotation_row = annotation_col)
pheatmap(test, annotation_col = annotation, annotation_row = annotation_col)

# change the color of annotation to what you want: (eg: "navy", "darkgreen")
Var1        <- c("navy", "darkgreen")
names(Var1) <- c("Exp1", "Exp2")
anno_colors <- list(Var1 = Var1)
pheatmap(test, annotation_col = annotation, annotation_colors = anno_colors)
