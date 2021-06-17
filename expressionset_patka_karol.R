#WCZYTANIE BIBLIOTEK
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("estrogen")
BiocManager::install("hgu95av2cdf")
BiocManager::install("affy")
BiocManager::install("Biobase")
BiocManager::install("convert")

library(affy)
library(Biobase)
library(hgu95av2cdf)
library(estrogen)
library(convert)
library(seqinr)
library(affyPLM)
library(stringr)

#WCZYTANIE SCIEZKI

sciezka = 'C:/STUDIA/Stopieñ II/Semestr I/WSP'
setwd(sciezka)

#WCZYTANIE OPISU I PRZYNALEZNOSCI DO KLAS
pd = read.AnnotatedDataFrame("datasetB.txt", header = TRUE, sep = "\t", row.names = 1) #czytamy annotacje klas z wczesniej "ulozonego" pliku
pData(pd) #tworzymy nowy obiekt na podstawie wczytanych annotacji

sciezka = 'C:/STUDIA/Stopieñ II/Semestr I/WSP/wsp'
setwd(sciezka)

#EXPRESSIONSET
#ZALADOWANIE PLIKU
dataDirectory <- system.file("extdata", package = "Biobase")
exprsFile <- file.path(dataDirectory, "datasetB.txt")
#UTWORZENIE MACIERZY
exprs <- as.matrix(read.table(exprsFile))
colnames(exprs)[c(1,2,3,4)] <- c("simple_annotation","CLASS","Sample","scan")
#KLASA OBIEKTU EXPRS
class(exprs)
#WYMIARY MACIERZY EXPRS
dim(exprs)
#WYSWIETLENIE NAZW KOLUMN OBIEKTU EXPRS
colnames(exprs)
#WYSWIETLENIE KOLUMN 1-3
head(exprs[, 1:4])
#TWORZENIE MINIMALNEGO OBIEKTU EXPRESSIONSET
minimalSet <- new("ExpressionSet", exprs = exprs)
#TWORZENIE FENOTYPOWEGO OBIEKTU EXPRESSIONSET
pDataFile <- file.path(dataDirectory, "datasetB.txt")
pData <- read.table(pDataFile)
colnames(pData)[c(1,2,3,4)] <- c("simple_annotation","CLASS","Sample","scan")
dim(pData)
colnames(pData)
#PODSUMOWANIE FENOTYPU OBIEKTU EXPRESSIONSET
summary(pData)
#SPRAWDZA CZY NAZWY WIERSZY I KOLUMN S¥ IDENTYCZNIE UPORZ¥DKOWANE
all(rownames(pData) == colnames(exprs))

