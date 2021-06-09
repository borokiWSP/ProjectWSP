#WCZYTANIE BIBLIOTEK
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("estrogen")
BiocManager::install("hgu95av2cdf")
BiocManager::install("affy")
BiocManager::install("Biobase")


library(affy)
library(Biobase)
library(hgu95av2cdf)
library(estrogen)

#WCZYTANIE SCIEZKI
sciezka = 'C:/Users/MagdalenaGierczyñska/Desktop/WSP/wsp'
setwd(sciezka)

#WCZYTANIE OPISU I PRZYNALEZNOSCI DO KLAS
pd = read.AnnotatedDataFrame("datasetB.txt", header = TRUE, sep = "\t", row.names = 1) #czytamy annotacje klas z wczesniej "ulozonego" pliku
pData(pd) #tworzymy nowy obiekt na podstawie wczytanych annotacji

#WCZYTANIE PLIKOW CEL Z OPISAMI I PRZYNALEZNOSCIA DO KLAS
celDat = ReadAffy(filenames = rownames(pData(pd)), phenoData = pd, verbose = TRUE, celfile.path = sciezka) #bierzemy wszystkie pliki z pliku pd, wczytujemy pliki z celfile.path, phenoData to informacje o annotacjach ; i na koniec dostajmy plik zmergowany czyli data + annotacje

allSets <- ls(hgu95av2cdf) #bierzemy probesety z brainarray
allSetDat <- mget(allSets, hgu95av2cdf) #tworzymy liste
hgu2 <- list2env(allSetDat) #tworzymy srodowisko (env) do prezentacji naszych informacji
celDat@cdfName <- "hgu2"

#NORMALIZACJA RMA
rma <- exprs(rma(celDat))

#ZAPISANIE PLIKU JAKO OBIEKT KLASY EXPRESSIONSET
save(rma, file="rma.Rdata") # zapisujemy nasz ExpressionSet jako najbardziej basic typ danych dla jezyka R czyli RData
