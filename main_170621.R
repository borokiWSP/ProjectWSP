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
BiocManager::install("genefilter") #dodałam nową bibliotekę więc jej jeszcze nie zakomendowałam 
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
library(genefilter)


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

# Określenie zależności pomiędzy próbami
pca<-prcomp(rma12,scale=TRUE)
plot(pca$x[,1], pca$x[,2])

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
gg_2<-ggplot(data=pca.data, aes(x=X, y=Y, label=Sample )) +
  geom_text() +
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep="")) +
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep="")) +
  theme_bw() +
  ggtitle("Analiza PCA")
gg_2 + stat_ellipse(aes(geom="polygon",level=0.99,alpha=0.2))

# Statystyka
wyniki<- pca$rotation[,1]
wyniki_abs <- abs(wyniki) #przyjęcie wartości bezwzglednej
ranking <- sort(wyniki_abs, decreasing=TRUE) #posortowanie malejąco

# Wybranie 10 najbardziej istotnych 
top10 <- names(ranking[1:10])
top10

# Wyświetlnie wartości 10 najbardziej istotnych
pca$rotation[top10,1]


### Klasteryzacja hierarchiczna z mapą ciepła - Martyna; DO WYŚWIETLENIA HEATMAPA

# Krok 1: Przygotowanie danych
#Według ustaleń biorę tylko klasę ADENO i NORMAL
ktory_adano_normal<- rma12@phenoData@data[["ADENO"]] %in% c("ADENO", "NORMAL") #sprawdza, który plik ma klasę anedo lub normal
adeno_normal <- rma12[, ktory_adano_normal] #tworze nową zmienną w której zapisane są tylko pliki z klasy adeno lub normal

#określam WARIANCJĘ dla wszystkich genów w próbkach i wybieram 
#200 genów o największej zmienności. (odrzucam te o najmniejsze wartości)
adeno_normal_wariancja <-apply(exprs(adeno_normal), 1, sd)
top200 <- names(sort(adeno_normal_wariancja , decreasing = TRUE))[1:200] #200 genów o największej zmienności
adeno_normal_200genow <- adeno_normal[top200, ] # nowa zmienna z wybranymi 200 genami 

# Krok 2: Metryka odległości

#funkcja odległosci 
dist_cor <- function(x) {
  as.dist(1 - cor(t(x), method = "pearson"))
}

# Krok 3: Metoda klasteryzacji 

# funkcja klasteryzacji 
clus_wd2 <- function(x) {
  hclust(x, method = "ward.D2")
}

# Krok 4: Mapa ceipła dla mikromacierzy

# zielony to geny 'down-regulated'
# czarny to geny neutralne (bez większych zmian w ekspresji)
# czerwony geny 'up-regulated' 
kolory <- colorRampPalette(c("green", "black", "red"))(n = 100) #wybieram kolory do heatmapy

#Definiowanie kolorów w pasku przynależności do klas- pierwszy "" to klasa dla której pasek będzie o kolorze 
#zdefiniowanym w drugim "", trzeci "" to kolor klasy niezdefiniowanej (normal)
klasy <- ifelse(adeno_normal_200genow@phenoData@data[["ADENO"]] == "ADENO","gray80" ,"gray40")

heatmap.2(exprs(adeno_normal_200genow), 
          #klasteryzacja
          distfun = dist_cor, 
          hclust = clus_wd2,
          
          # skalowanie sprawiające, że poszczególne geny będą w rzędach scaling 
          scale = "row",
          
          #kolory heatmapy
          col = kolory,
          
          # podpisy kolumn - bez określenia kolumn pojawia się nazwa pliku CEL
          # w zmienionej wersji może pojawić się: adeno_normal_200genow@phenoData@data[["AD10"]]
          #lub adeno_normal_200genow@phenoData@data[["ADA10T1_A389_7"]]. Ja wybrałam pierwszą opcję 
          #czyli "simple_annotation" bo wydawało mi się to bardziej przejrzyste
          labCol=(adeno_normal_200genow@phenoData@data[["AD10"]]),
          
          # pasek przynależności do klas
          ColSideColors = klasy, 
          
          # sprawia, że mamy ładą heatmapę a nie histogramy dla każdej z danych
          trace = "none",
          density.info = "none")
#--------------------

# Patka i Karol
sciezka = 'C:/STUDIA/Stopień II/Semestr I/WSP'
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
#SPRAWDZA CZY NAZWY WIERSZY I KOLUMN SĽ IDENTYCZNIE UPORZĽDKOWANE
all(rownames(pData) == colnames(exprs))

