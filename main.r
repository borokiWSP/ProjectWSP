#WCZYTANIE BIBLIOTEK
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("hgu133plus2.db")
BiocManager::install("estrogen")
BiocManager::install("hgu95av2cdf")
BiocManager::install("affy")
BiocManager::install("Biobase")
BiocManager::install("RColorBrewerr")
BiocManager::install("gplots")

install.packages("RColorBrewer")
install.packages("devtools")
devtools::install_github("tidyverse/stringr")
install.packages("stringr")
install.packages("gplots")

library(affy)
library(Biobase)
#library(hgu95av2cdf)
library(estrogen)
library(seqinr)
library(affyPLM)
library(stringr)
library(RColorBrewer)
library(gplots)

#WCZYTANIE SCIEZKI
sciezka = '/Users/madzia/Downloads/wsp'
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


# Określenie zależności pomiędzy próbami
pca<-prcomp(rma,scale=TRUE)

# Wyznaczenie odchylenia standardowego - wariancji
pca.var <- pca$sdev^2

# Procentowo
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)

# Wyświetlenie histogramu reprezentującego wartości procentowe wariancji
barplot(pca.var.per, main="Histogram", xlab="Składowe", ylab="Różnica procentowa", col=rainbow(length(pca.var.per)))

# Wyznaczenie wartości potrzebnych do wyrysowania ggplot 
pca.data <- data.frame(Sample=rownames(pca$x),
  X=pca$x[,1],
  Y=pca$x[,2])
pca.data
 
# Wyrysowanie wykresu ggplot
gg_2<-ggplot(data=pca.data, aes(x=X, y=Y, label=Sample, colour=red)) +
  geom_text() +
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep="")) +
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep="")) +
  theme_bw() +
  ggtitle("Analiza PCA")
#gg_2 + stat_ellipse(aes(fill = green), geom="polygon",level=0.99,alpha=0.2)

# Statystyka
wyniki<- pca$rotation[,1]
wyniki_abs <- abs(wyniki) #przyjęcie wartości bezwzglednej
ranking <- sort(wyniki_abs, decreasing=TRUE) #posortowanie malejąco

# Wybranie 10 najbardziej istotnych 
top10 <- names(ranking[1:10])
top10

# Wyświetlnie wartości 10 najbardziej istotnych
pca$rotation[top10,1]




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
#SPRAWDZA CZY NAZWY WIERSZY I KOLUMN SĄ IDENTYCZNIE UPORZĄDKOWANE
all(rownames(pData) == colnames(exprs))





