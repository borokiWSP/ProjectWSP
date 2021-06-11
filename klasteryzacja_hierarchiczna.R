# ----------------------------------------------------------
# Klasteryzacja hierarchiczna z map¹ ciep³a
# ----------------------------------------------------------

klateryzacja_hierarchiczna <- function (dane_po_normalizacji_rma) {
  
  
  BiocManager::install("gplots")
  BiocManager::install("RColorBrewer")
  
  library(gplots)
  library(RColorBrewer)
  
  
  #WCZYTANIE SCIEZKI
  sciezka = 'C:/Users/marty/Desktop/WSP/PROJEKT/wsp'
  setwd(sciezka)
  
  
  # ----------------------------------------------------------
  # Krok 1: Przygotowanie danych
  # ----------------------------------------------------------
  
  #Wed³ug ustaleñ biorê tylko klasê ADENO i NORMAL
  ktory_adano_normal<- dane_po_normalizacji_rma12@phenoData@data[["ADENO"]] %in% c("ADENO", "NORMAL") #sprawdza, który plik ma klasê anedo lub normal
  adeno_normal <- dane_po_normalizacji_rma[, ktory_adano_normal] #tworze now¹ zmienn¹ w której zapisane s¹ tylko pliki z klasy adeno lub normal
  
  #okreœlam WARIANCJÊ dla wszystkich genów w próbkach i wybieram 
  #200 genów o najwiêkszej zmiennoœci. (odrzucam te o najmniejsze wartoœci)
  adeno_normal_wariancja <- rowVars(exprs(adeno_normal)) #Wariancja 
  top200 <- names(sort(adeno_normal_wariancja, decreasing = TRUE))[1:200] #200 genów o najwiêkszej zmiennoœci
  adeno_normal_200genow <- adeno_normal[top200, ] # nowa zmienna z wybranymi 200 genami 
  
  # ----------------------------------------------------------
  # Krok 2: Metryka odleg³oœci
  # ----------------------------------------------------------
  
  #funkcja odleg³osci 
  dist_cor <- function(x) {
    as.dist(1 - cor(t(x), method = "pearson"))
  }
  
  # ----------------------------------------------------------
  # Krok 3: Metoda klasteryzacji 
  # ----------------------------------------------------------
  
  # funkcja klasteryzacji 
  clus_wd2 <- function(x) {
    hclust(x, method = "ward.D2")
  }
  
  # ----------------------------------------------------------
  # Krok 4: Mapa ceip³a dla mikromacierzy
  # ----------------------------------------------------------
  # zielony to geny 'down-regulated'
  # czarny to geny neutralne (bez wiêkszych zmian w ekspresji)
  # czerwony geny 'up-regulated' 
  kolory <- colorRampPalette(c("green", "black", "red"))(n = 100) #wybieram kolory do heatmapy
  
  #Definiowanie kolorów w pasku przynale¿noœci do klas- pierwszy "" to klasa dla której pasek bêdzie o kolorze 
  #zdefiniowanym w drugim "", trzeci "" to kolor klasy niezdefiniowanej (normal)
  klasy <- ifelse(adeno_normal_200genow@phenoData@data[["ADENO"]] == "ADENO","gray80" ,"gray40")
  
  heatmap.2(exprs(adeno_normal_200genow), 
            #klasteryzacja
            distfun = dist_cor, 
            hclust = clus_wd2,
            
            # skalowanie sprawiaj¹ce, ¿e poszczególne geny bêd¹ w rzêdach scaling 
            scale = "row",
            
            #kolory heatmapy
            col = kolory,
            
            # podpisy kolumn - bez okreœlenia kolumn pojawia siê nazwa pliku CEL
            # w zmienionej wersji mo¿e pojawiæ siê: adeno_normal_200genow@phenoData@data[["AD10"]]
            #lub adeno_normal_200genow@phenoData@data[["ADA10T1_A389_7"]]. Ja wybra³am pierwsz¹ opcjê 
            #czyli "simple_annotation" bo wydawa³o mi siê to bardziej przejrzyste
            labCol=(adeno_normal_200genow@phenoData@data[["AD10"]]),
            
            # pasek przynale¿noœci do klas
            ColSideColors = klasy, 
            
            # sprawia, ¿e mamy ³ad¹ heatmapê a nie histogramy dla ka¿dej z danych
            trace = "none",
            density.info = "none")
}



