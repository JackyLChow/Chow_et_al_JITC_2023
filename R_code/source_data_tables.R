table(sce$class, sce$sample)
write.csv(table(sce$class, sce$sample), "tab_cellxsample.csv")
write.csv(round(prop.table(table(sce$class, sce$sample), 2), 3), "tab_cellxsample_prop.csv")

table(tnk$subclass_fine, tnk$sample)
write.csv(table(tnk$subclass_fine, tnk$sample), "tab_cellxsample_tnk.csv")
write.csv(round(prop.table(table(tnk$subclass_fine, tnk$sample), 2), 3), "tab_cellxsample_tnk_prop.csv")

table(mye$subclass_fine, mye$sample)
write.csv(table(mye$subclass_fine, mye$sample), "tab_cellxsample_mye.csv")
write.csv(round(prop.table(table(mye$subclass_fine, mye$sample), 2), 3), "tab_cellxsample_mye_prop.csv")
