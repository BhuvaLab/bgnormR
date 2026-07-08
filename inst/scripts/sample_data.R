devtools::load_all()
img = read_qptiff("inst/extdata/20240415_PA20P17131_HNC_Scan1.er.qptiff", level = 5)
markers = c("CD20", "CD3e", "CD8", "PanCK", "Vimentin")
write_qptiff(img[(1215-550+1):1215,51:850,markers], "inst/extdata/PA_HNC_sample.ome.tiff")
