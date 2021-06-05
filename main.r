#WCZYTANIE BIBLIOTEK
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("hgu133plus2.db")
BiocManager::install("estrogen")
BiocManager::install("hgu95av2cdf")
BiocManager::install("affy")
BiocManager::install("Biobase")


library(affy)
library(Biobase)
#library(hgu95av2cdf)
library(estrogen)
library(seqinr)
library(affyPLM)

#WCZYTANIE SCIEZKI
sciezka = 'C:/Users/Aneta/Desktop/studia_mgr/semestr I/WSP/laboratorium'
setwd(sciezka)

#WCZYTANIE OPISU I PRZYNALEZNOSCI DO KLAS
pd = read.AnnotatedDataFrame("datasetA scans.txt", header = TRUE, sep = "\t", row.names = 1)
pData(pd)

#WCZYTANIE PLIKOW CEL Z OPISAMI I PRZYNALEZNOSCIA DO KLAS
celDat = ReadAffy(filenames = rownames(pData(pd)), phenoData = pd, verbose = TRUE, 
                  celfile.path = sciezka)
#NORMALIZACJA RMA
rma <- rma(celDat)

#ZAPISANIE PLIKU JAKO OBIEKT KLASY EXPRESSIONSET
save(rma, file="rma2.Rdata")

# ogólnie jakby ktoœ chcia³ to moze sie zastanowiæ dlaczego wczytuje tylo 195 rekordów, a nie wszystkie, czy to ju¿ jest jakaœ filtracja? xd
operacje_na_sondach(rma)