#WCZYTANIE BIBLIOTEK
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("hgu133plus2.db")
BiocManager::install("estrogen")
BiocManager::install("hgu95av2cdf")
BiocManager::install("affy")
BiocManager::install("Biobase")
install.packages("stringr")
install.packages("devtools")
devtools::install_github("tidyverse/stringr")


library(affy)
library(Biobase)
#library(hgu95av2cdf)
library(estrogen)
library(seqinr)
library(affyPLM)
library(stringr)

#WCZYTANIE SCIEZKI
sciezka = 'C:/Users/Aneta/Desktop/studia_mgr/semestr I/WSP/laboratorium'
setwd(sciezka)

#WCZYTANIE OPISU I PRZYNALEZNOSCI DO KLAS
pd = read.AnnotatedDataFrame("datasetA scans.txt", header = TRUE, sep = "\t", row.names = 1)
pData(pd)

#WCZYTANIE PLIKOW CEL Z OPISAMI I PRZYNALEZNOSCIA DO KLAS
celDat = ReadAffy(filenames = rownames(pData(pd)), phenoData = pd, verbose = TRUE, 
                  celfile.path = sciezka)
#NORMALIZACJA RMA - przemyślenia Anety: jak robiliśmy projekt to kazał nam robić z expresso i precyzowac te warunki, bo coś tam, myślicie, ze może zostać jednak funkcja RMA?
rma <- rma(celDat)

#ZAPISANIE PLIKU JAKO OBIEKT KLASY EXPRESSIONSET - kolejne przemyślnia: @Patrycja mi w sumie to wszystko dobrze działało na tym etapie więc to nie tak, ze olałąm Twój kod, ale działa. xd
save(rma, file="rma2.Rdata")

# ogólnie jakby ktoś chciał to moze sie zastanowić dlaczego wczytuje tylo 195 rekordów, a nie wszystkie, czy to już jest jakaś filtracja? xd
# sie dowiedziałam już, ze tak jakby ktos też się zastanawiał. xd
operacje_na_sondach(rma)
