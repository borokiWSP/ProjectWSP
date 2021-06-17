#WCZYTANIE BIBLIOTEK
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("hgu133plus2.db")
# BiocManager::install("estrogen")
# BiocManager::install("hgu95av2cdf")
# BiocManager::install("affy")
# BiocManager::install("Biobase")
# BiocManager::install("RColorBrewerr")
# BiocManager::install("gplots")
# 
# install.packages("RColorBrewer")
# install.packages("devtools")
# devtools::install_github("tidyverse/stringr")
# install.packages("stringr")
# install.packages("gplots")



#moje biblioteki
library(affy)
library(Biobase)
library(hgu95av2cdf)
library(estrogen)
library(seqinr)
library(affyPLM)
library(stringr)
library(RColorBrewer)
library(gplots)
library(devtools)
library(tidyverse)
library(stringr)
library(hgu133plus2.db)
library(matrixStats)
library(convert)

BiocManager::install("convert")

#WCZYTANIE SCIEZKI
sciezka = 'C:/Users/alicj/Desktop/WSP-aplikacja/wsp'
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
rma12 <- rma(celDat)

#ZAPISANIE PLIKU JAKO OBIEKT KLASY EXPRESSIONSET
save(rma12, file="rma12.Rdata") # zapisujemy nasz ExpressionSet jako najbardziej basic typ danych dla jezyka R czyli RData

# FUNKCJA ANETY + NATALII
operacje_na_sondach(rma12)

# Okreœlenie zale¿noœci pomiêdzy próbami
pca<-prcomp(rma12,scale=TRUE)
plot(pca$x[,1], pca$x[,2])

# Wyznaczenie odchylenia standardowego - wariancji
pca.var <- pca$sdev^2

# Procentowo
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)

# Wyœwietlenie histogramu reprezentuj¹cego wartoœci procentowe wariancji
barplot(pca.var.per, main="Histogram", xlab="Sk³adowe", ylab="Ró¿nica procentowa", col=rainbow(length(pca.var.per)))

# Wyznaczenie wartoœci potrzebnych do wyrysowania ggplot 
pca.data <- data.frame(Sample=rownames(pca$x),
                       X=pca$x[,1],
                       Y=pca$x[,2])
pca.data

# Wyrysowanie wykresu ggplot
gg_2<-ggplot(data=pca.data, aes(x=X, y=Y, label=Sample )) +
  geom_text() +
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep="")) +
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep="")) +
  theme_bw() +
  ggtitle("Analiza PCA")
gg_2 + stat_ellipse(aes(geom="polygon",level=0.99,alpha=0.2))

# Statystyka
wyniki<- pca$rotation[,1]
wyniki_abs <- abs(wyniki) #przyjêcie wartoœci bezwzglednej
ranking <- sort(wyniki_abs, decreasing=TRUE) #posortowanie malej¹co

# Wybranie 10 najbardziej istotnych 
top10 <- names(ranking[1:10])
top10

# Wyœwietlnie wartoœci 10 najbardziej istotnych
pca$rotation[top10,1]


######  Klasteryzacja hierarchiczna z map¹ ciep³a - Martyna 
klateryzacja_hierarchiczna(rma12)

# Patka i Karol
sciezka = 'C:/STUDIA/Stopieñ II/Semestr I/WSP'
setwd(sciezka)

#WCZYTANIE SCIEZKI
sciezka = 'C:/Users/alicj/Desktop/WSP-aplikacja/wsp'
setwd(sciezka)

#WCZYTANIE OPISU I PRZYNALEZNOSCI DO KLAS
pd = read.AnnotatedDataFrame("datasetB.txt", header = TRUE, sep = "\t", row.names = 1) #czytamy annotacje klas z wczesniej "ulozonego" pliku
pData(pd) #tworzymy nowy obiekt na podstawie wczytanych annotacji


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
#SPRAWDZA CZY NAZWY WIERSZY I KOLUMN S¼ IDENTYCZNIE UPORZ¼DKOWANE
all(rownames(pData) == colnames(exprs))

