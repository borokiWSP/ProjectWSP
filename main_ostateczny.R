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
# BiocManager::install("genefilter") #dodałam nową bibliotekę więc jej jeszcze nie zakomendowałam aby każdy widział, że ją też trzeba jeszcze wgrać 
# devtools::install_github("davidgohel/officer")
# install.packages("officer")
# BiocManager::install("convert")
# install.packages("RColorBrewer")
# install.packages("devtools")
# devtools::install_github("tidyverse/stringr")
# install.packages("stringr")
# install.packages("gplots")
install.packages("factoextra")
install.packages("doBy")

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
library(officer)
library(flextable)
library(factoextra)
library(doBy)

#WCZYTANIE SCIEZKI
sciezka = 'C:/Users/Aneta/Desktop/studia_mgr/semestr I/WSP/laboratorium'
setwd(sciezka)

#WCZYTANIE OPISU I PRZYNALEZNOSCI DO KLAS
pd = read.AnnotatedDataFrame("datasetB.txt", header = TRUE, sep = "\t", row.names = 1) #czytamy annotacje klas z wczesniej "ulozonego" pliku
pData(pd) #tworzymy nowy obiekt na podstawie wczytanych annotacji

sciezka = 'C:/Users/Aneta/Desktop/studia_mgr/semestr I/WSP/laboratorium/wsp'
setwd(sciezka)

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
wszystkie_klasy = operacje_na_sondach(rma12) # wyciągam wsyztkie 

adeno_min = data.frame(wszystkie_klasy[1]) # rozkładam je by się ładnie wyświetlały bo inaczej się niestety nie da, albo przynajmneij nie umiem
names_adeno_min = row.names(adeno_min)
adeno_min_df = data.frame(names_adeno_min, adeno_min)
names(adeno_min_df)[2] <- "Ekspresja" # zmieniam nazwy kolumn
names(adeno_min_df)[1] <- "Nazwa sondy w klasie ADENO"
# % sond o najmniejszej ekspresji w klasie adeno, 5% sond o największej ekspresji w klasie adeno,5% sond o najmniejszej ekspresji w klasie normal,5% sond o największej ekspresji w klasie normal 

adeno_max = data.frame(wszystkie_klasy[2])
names_adeno_max = row.names(adeno_max)
adeno_max_df = data.frame(names_adeno_max, adeno_max)
names(adeno_max_df)[1] <- "Nazwa sondy w klasie ADENO"
names(adeno_max_df)[2] <- "Ekspresja"

normal_min = data.frame(wszystkie_klasy[3])
names_normal_min = row.names(normal_min)
normal_min_df = data.frame(names_normal_min, normal_min)
names(normal_min_df)[1] <- "Nazwa sondy w klasie NORMAL"
names(normal_min_df)[2] <- "Ekspresja"

normal_max = data.frame(wszystkie_klasy[4])
names_normal_max = row.names(normal_max)
normal_max_df = data.frame(names_normal_max, normal_max)
names(normal_max_df)[1] <- "Nazwa sondy w klasie NORMAL"
names(normal_max_df)[2] <- "Ekspresja"

# jenda zbiorcza df:
all_sond = data.frame(adeno_min_df, adeno_max_df, normal_min_df, normal_max_df) # to do apki
all_adeno = data.frame(adeno_min_df, adeno_max_df) # to i to nizej do raportu, bo się nie mieści
all_normal = data.frame(normal_min_df, normal_max_df)


# ANALIZA PCA

#wyciagniecie indeksow z macierzy


ekspresja = exprs(rma12)
indeksy_do_grupy_adeno = which(rma12@phenoData@data[["ADENO"]] == "ADENO", arr.ind = TRUE)
indeksy_do_grupy_normal = which(rma12@phenoData@data[["ADENO"]] == "NORMAL", arr.ind = TRUE)
nowa_macierz=ekspresja[,c(indeksy_do_grupy_adeno,indeksy_do_grupy_normal)]


# Określenie zależnościci pomiędzy próbbami
rma12_exprs <- nowa_macierz
pca<-prcomp(t(rma12_exprs), scale=TRUE)
plot(pca$x[,1], pca$x[,2])

# Wyznaczenie odchylenia standardowego - wariancji
pca.var <- pca$sdev^2

# Procentowo
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)

# Przkładowo 5 najbardziej znaczczących
n_ogr <- 5

# Wyświetlenie histogramu reprezentującego wartości procentowe wariancji
fviz_eig(pca, ncp = n_ogr)

# Wyznaczenie wartości potrzebnych do wyrysowania ggplot
max_perc <- which.maxn(pca.var.per, n_ogr)
x_filt <- pca$x[max_perc,]
 
pca.data <- data.frame(Sample=rownames(pca$x),
                       X=pca$x[,1],
                       Y=pca$x[,2])
pca.data



length(pca.data)



rma12_exprs2 <- rma12_exprs[, 1:n_ogr]
Klasy <- colnames(rma12_exprs2)
nazwy <- unique(colnames(rma12_exprs2))
numeracja <- seq_len(ncol(rma12_exprs2))

for(j in seq_along(nazwy)){
	for(i in seq_len(ncol(rma12_exprs2))){
      if(Klasy[i]==nazwy[j]){
        numeracja[i] <- j}
    }
}

 # a=1:195
 # a[1:156]="ADENO"
 # a[157:170]="NORMAL"
 # a[171:195]="OTHERS"
 # P=matrix(c(1:195, a), nrow=195, ncol=2)
  Klasy=(rma12@phenoData@data[["ADENO"]])[c(indeksy_do_grupy_adeno,indeksy_do_grupy_normal)] 

# Wyrysowanie wykresu ggplot
gg_2<-ggplot(data=pca.data, aes(x=X, y=Y, label=".", color=Klasy)) +
  geom_text() +
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep="")) +
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep="")) +
  theme_bw() +
  ggtitle("Analiza PCA")
gg_2 + stat_ellipse(aes(fill = Klasy), geom="polygon",level=0.99,alpha=0.2)

# Statystyka
wyniki<- pca$rotation[,1]
wyniki_abs <- abs(wyniki) #przyjęcie wartości bezwzglednej
ranking <- sort(wyniki_abs, decreasing=TRUE) #posortowanie malejąco

# Wybranie 10 najbardziej istotnych 
top10 <- names(ranking[1:10])
top10

# Wyświetlnie wartości 10 najbardziej istotnych
pca$rotation[top10,1]

### TUTAJ ZACZYNA SIE MOJA CZĘŚĆ - NATALIA
# Wykorzystanie funkcji Anety
# operacje_na_sondach(rma12)
# Wyciągnięte z funkcji :
ekspresja = exprs(rma12)
# FILTROWANIE SOND KONTROLNYCH
nazwy_sond = row.names(ekspresja)
indeksy_sond_kontrolnych = which(str_detect(nazwy_sond, "^A"))
pierwsza_kontrolna = indeksy_sond_kontrolnych[[1]]

# ekspresja proby ADENO
indeksy_do_grupy_adeno = which(rma12@phenoData@data[["ADENO"]] == "ADENO", arr.ind = TRUE)
# kontrola = obiekt_expression_set@phenoData@data[["ADENO"]][which(obiekt_expression_set@phenoData@data[["ADENO"]] == "ADENO", arr.ind = TRUE)]
ekspresja_adeno = ekspresja[,indeksy_do_grupy_adeno[[1]]:indeksy_do_grupy_adeno[length(indeksy_do_grupy_adeno)]] # ADENO
srednia_ekspresja_adeno = rowSums(ekspresja_adeno)/length(indeksy_do_grupy_adeno)
srednia_ekspresja_adeno_po_filtracji_sond_kontrolnych = srednia_ekspresja_adeno[1:pierwsza_kontrolna-1]

# ekspresja proby NORMAL
indeksy_do_grupy_normal = which(rma12@phenoData@data[["ADENO"]] == "NORMAL", arr.ind = TRUE)
ekspresja_normal = ekspresja[,indeksy_do_grupy_normal[[1]]:indeksy_do_grupy_normal[length(indeksy_do_grupy_normal)]] # NORMAL
srednia_ekspresja_normal = rowSums(ekspresja_normal)/length(indeksy_do_grupy_normal)
srednia_ekspresja_normal_po_filtracji_sond_kontrolnych = srednia_ekspresja_normal[1:pierwsza_kontrolna-1]

# STATYSTYKA
# PORÓWNANIE GRUP NORMAL I ADENO 
# TRANSPOZYCJA MACIERZY
T_N = t(ekspresja_normal)
T_A = t(ekspresja_adeno)
# PĘTLA DO WYKONYWANIA TESTU WILCOXONA
wynik_testu=c() # pusty wektor na wyniki
for (i in seq(length(nazwy_sond))) {
  test=wilcox.test(T_N[,i], T_A[,i], paired = FALSE, alternative = "two.sided")
  wynik_testu<-c(wynik_testu,test$p.value)}
# PĘTLA DO WYCIĄGANIA NAZW SOND DLA KTÓRYCH P-WARTOŚĆ W TEŚCIE MA OKREŚLONĄ WARTOŚĆ
Nazwy = data.frame(nazwy_sond)
wynik=c() #pusty wektor na p-wartości
for(i in seq(wynik_testu)){
  if (wynik_testu[i]<= 0.0000000009){ # W TYM MIEJSCU DAJEMY WYBRANĄ P-WARTOŚĆ - WIEM ŻE NISKA ALE TO PRZEZ ILOŚĆ TESTÓW (POPRAWKA NA WIELOKROTNE TESTOWANIE)
    wynik<-c(wynik,nazwy_sond[i]) #wektor na nazwy sond
  }
}
# DO WYŚWIETLENIA - HISTOGRAM P-WARTOŚCI
# HISTOGRAM P- WARTOŚCI 
h =  hist(wynik_testu,
          main="Histogram p - wartości dla testu Wilcoxona", 
          xlab="p-wartości", 
          ylab = "Częstość występowania",
          border="black", 
          col="blue")
# DO WYŚWIETLENIA NAZWY SOND
Nazwy_sond_po_tescie = data.frame(wynik) # wyświetlenie nazw sond dla których p- wartość jest mniejsza niż zadana
# DO WYŚWIETLENIA OGÓLNA LICZBA WYNIKÓW 
l_wynikow = length(which(wynik_testu < 0.0000000009))# Liczba wyników które uzyskaliśmy 
# WYBRANE (po teście) GENY RÓŻNICUJĄCE
# BiocManager::install('hgu95av2')
library("hgu95av2")
# W SYMBOLACH SĄ TE GENY Z TESTU
symbole<- c("1001_at" , "1596_g_at" ,"1814_at", "268_at" ,"34708_at","35868_at", "36569_at", "37196_at", "37247_at" ,"37398_at","38044_at", "38177_at", "40841_at" )
Geny_roznicujace =mget(symbole, hgu95av2GENENAME,ifnotfound = NA)
# DO WYŚWIETLENIA NAZWY GENÓW RÓŻNICUJĄCYCH (GG)
GG1 = data.frame(Geny_roznicujace[1])
GG2 = data.frame(Geny_roznicujace[2])
GG3 = data.frame(Geny_roznicujace[3])
GG4 = data.frame(Geny_roznicujace[4])
GG5 = data.frame(Geny_roznicujace[5])
GG6 = data.frame(Geny_roznicujace[6])
GG7 = data.frame(Geny_roznicujace[7])
GG8 = data.frame(Geny_roznicujace[8])
GG9 = data.frame(Geny_roznicujace[9])
GG10 = data.frame(Geny_roznicujace[10])
Geny5 = data.frame(GG1, GG2, GG3, GG4, GG5)
Geny10 = data.frame(GG6, GG7, GG8, GG9, GG10)

# do pttx
stat = data.frame(pca$rotation[top10,1])
names_stat = row.names(stat)
stat_df = data.frame(names_stat, stat)
names(stat_df)[1] <- "Nazwy najistotniejszych" # zmieniam nazwy kolumn
names(stat_df)[2] <- "Wartości najistotniejszych "

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

hmap = heatmap.2(exprs(adeno_normal_200genow), 
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

# GENERACJA RAPORTU W PPTX
dok <- read_pptx()
loc_title <- ph_location_type(type = "title")
dok = add_slide(dok)
dok = ph_with(dok, value = "Raport z projektu z WSP - grupa II", location = loc_title)
dok <- add_slide(dok)
dok <- ph_with(dok, value = "5% sond o najmniejszej i największej ekspresji w klasie adeno", location = loc_title)
dok <- ph_with(dok, value = all_adeno, location = ph_location(width = 4, height = 4))
dok <- add_slide(dok)
dok <- ph_with(dok, value = "5% sond o najmniejszej i największej ekspresji w klasie normal", location = loc_title)
dok <- ph_with(dok, value = all_normal, location = ph_location_left())
dok = add_slide(dok)
dok = ph_with(dok, value = "Histogram z analizy PCA", location = loc_title)
dok = ph_with(dok, value = pca_hist, location = ph_location_left())
dok = add_slide(dok)
dok = ph_with(dok, value = "Wykres z analizy PCA", location = loc_title)
dok = ph_with(dok, value = pca_anal, location = ph_location_left())
dok = add_slide(dok)
dok = ph_with(dok, value = "Geny różnicujące cz1", location = loc_title)
dok = ph_with(dok, value = Geny5, location = ph_location_left())
dok = add_slide(dok)
dok = ph_with(dok, value = "Geny różnicujące cz2", location = loc_title)
dok = ph_with(dok, value = Geny10, location = ph_location_left())
dok = add_slide(dok)
dok = ph_with(dok, value = "Histogram Statystyka", location = loc_title)
#dok = ph_with(dok, value = h, location = ph_location_left())
dok = add_slide(dok)
dok = ph_with(dok, value = "Wyniki testu", location = loc_title)
dok = ph_with(dok, value = Nazwy_sond_po_tescie, location = ph_location_left())
dok = ph_with(dok, value = l_wynikow, location = ph_location_right())
dok = add_slide(dok)
dok = ph_with(dok, value = "Heatmapa", location = loc_title)
#dok = ph_with(dok, value = hmap, location = ph_location_left())
print(dok, target = "WSP_1.pptx")
