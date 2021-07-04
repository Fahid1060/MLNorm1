#Convert Dataset into XPM
XPM <- function(x) {
    XPMCpp(x)
}
#Convert Dataset into XPM and transpose
tXPM <- function(x) {
    tXPMCpp(x)
}


BatchEffectLinnorm1 <- function(x,cleanData, minNonZeroPortion, BE_F_LC_Genes = 0.25,BE_F_HC_Genes = 0.05, BE_F_p = 0.5, BE_strength = 0.25, spikein = NULL) {
	#Save a copy of the raw matrix. x2 will not be filtered and used for normalization. x will be filtered and used as the model.
	x2 <- x

	CN <- colnames(x2)
	RN <- rownames(x2)

	x2 <- BatchEffect2(cleanData, x2, as.vector(sort(rowMeans(cleanData),decreasing = FALSE)), BE_strength)
	 colnames(x2) <- CN
	 rownames(x2) <- RN

	#Output
	return (x2)
}

BatchEffect2 <- function(x,y,z,z2) {
	BatchEffectCpp(x,y,z,z2)
}

