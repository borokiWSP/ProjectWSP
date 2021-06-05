operacje_na_sondach <- function(obiekt_expression_set) {
  
  # EKSRTRAKCJE SOND
  ekspresja = exprs(obiekt_expression_set)
  
  # ekspresja kontroli (ADENO)
  indeksy_do_grupy_kontrolnej = which(obiekt_expression_set@phenoData@data[["ADENO"]] == "ADENO", arr.ind = TRUE)
  # kontrola = obiekt_expression_set@phenoData@data[["ADENO"]][which(obiekt_expression_set@phenoData@data[["ADENO"]] == "ADENO", arr.ind = TRUE)]
  ekspresja_kontroli = ekspresja[,indeksy_do_grupy_kontrolnej[[1]]:indeksy_do_grupy_kontrolnej[length(indeksy_do_grupy_kontrolnej)]] # ADENO
  srednia_ekspresja_konroli = rowSums(ekspresja_kontroli)/length(indeksy_do_grupy_kontrolnej)
  
  # ekspresja proby badawczej (CARCINOID)
  indeksy_do_grupy_carcinoid = which(obiekt_expression_set@phenoData@data[["ADENO"]] == "CARCINOID", arr.ind = TRUE)
  ekspresja_carcinoid = ekspresja[,indeksy_do_grupy_carcinoid[[1]]:indeksy_do_grupy_carcinoid[length(indeksy_do_grupy_carcinoid)]] # CARCINOID
  srednia_ekspresja_carcinoid = rowSums(ekspresja_carcinoid)/length(indeksy_do_grupy_carcinoid)
  
  # ekspresja proby badawczej (NORMAL)
  indeksy_do_grupy_normal = which(obiekt_expression_set@phenoData@data[["ADENO"]] == "NORMAL", arr.ind = TRUE)
  ekspresja_normal = ekspresja[,indeksy_do_grupy_normal[[1]]:indeksy_do_grupy_normal[length(indeksy_do_grupy_normal)]] # NORMAL
  srednia_ekspresja_normal = rowSums(ekspresja_normal)/length(indeksy_do_grupy_normal)
  
  # ekspresja proby badawczej (SMALL CELL)
  indeksy_do_grupy_smallcell = which(obiekt_expression_set@phenoData@data[["ADENO"]] == "SMALLCELL", arr.ind = TRUE)
  ekspresja_smallcell = ekspresja[,indeksy_do_grupy_smallcell[[1]]:indeksy_do_grupy_smallcell[length(indeksy_do_grupy_smallcell)]] # SMALL CELL
  srednia_ekspresa_smallcell = rowSums(ekspresja_smallcell)/length(indeksy_do_grupy_smallcell)
  
  # ekspresja proby badawczej (SQUAMOUS)
  indeksy_do_grupy_squamous = which(obiekt_expression_set@phenoData@data[["ADENO"]] == "SQUAMOUS", arr.ind = TRUE)
  ekspresja_squamous = ekspresja[,indeksy_do_grupy_squamous[[1]]:indeksy_do_grupy_squamous[length(indeksy_do_grupy_squamous)]] # SQUAMOUS
  srednia_ekspresja_squamous = rowSums(ekspresja_squamous)/length(indeksy_do_grupy_squamous)
  
  # FILTROWANIE SOND - wzrost o 30% wzglêdem kontroli (ADENO) ?? No chyba, ¿e chodzi o filtrowanie 
  roznica = srednia_ekspresja_konroli - srednia_ekspresja_carcinoid
  
  wzrost_nazwy_badanie = names(roznica[which(roznica >= log2(1.3))])
  wzrost_badanie = length(names(roznica[which(roznica >= log2(1.3))]))
  
  spadek_nazwy_badanie = (names(roznica[which(roznica <= -log2(1.3))]))
  spadek_badanie = length(names(roznica[which(roznica <= -log2(1.3))]))
  
  
  # 5% SOND O MAX EKSPRESJI
  ile_sond_to_5_procent = 12625*0.05/100
  ile_sond_to_5_procent = floor(ile_sond_to_5_procent)
  
  kolejno = sort(srednia_ekspresja_konroli)
  najwieksza_srednia_ekspresja_5procent = head(kolejno, ile_sond_to_5_procent)
  najmniejsza_srednia_ekspresja_5procent = tail(kolejno, ile_sond_to_5_procent)
  
  dane_wyjsciowe = list("5% sond o najmniejszej ekspresji" = najmniejsza_srednia_ekspresja_5procent,"5% sond o najwiêkszej ekspresji" = najwieksza_srednia_ekspresja_5procent)
  return(dane_wyjsciowe)
  
  #return(przefiltrowane sondy)
}
