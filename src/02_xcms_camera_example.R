library (xcms)
library (CAMERA)
xset <- xcmsSet(method='centWave', ppm=20, snthresh=10, peakwidth=c(5,18))
xset <- group(xset, method=“density“, minfrac=0.5, minsamp=1, bw=10, mzwid=0.01, sleep=0.001)
xset2 <- retcor(xset, family= “s“, plottype= “m“, missing=1, extra=1, span=1)
xset2 <- retcor(xset, family= “s“, plottype= “m“, missing=1, extra=1, span=1)
xset3 <- group(xset2, method=“density“, mzwid=0.01, sleep=0.001, minfrac=0.5, minsamp=1, bw=5)
xset4 <- fillPeaks(xset3)
an <- xsAnnotate(xset4)
#Creation of an xsAnnotate object
anF <- groupFWHM(an, perfwhm = 0.6)
#Perfwhm = parameter defines the window width, which is used for matching
anI <- findIsotopes(anF, mzabs=0.01)
#Mzabs = the allowed m/z error
anIC <- groupCorr(anI, cor_eic_th=0.75)
anFA <- findAdducts(anIC, polarity=“positive“)
write.csv(getPeaklist(anFA), file=“Name.csv")