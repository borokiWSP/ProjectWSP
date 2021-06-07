operacje_na_sondach <- function(obiekt_expression_set) {
  
  # EKSRTRAKCJE SOND 
  ekspresja = exprs(obiekt_expression_set)
  
  # FILTROWANIE SOND KONTROLNYCH
  nazwy_sond = row.names(ekspresja)
  indeksy_sond_kontrolnych = which(str_detect(nazwy_sond, "^A"))
  pierwsza_kontrolna = indeksy_sond_kontrolnych[[1]]
  #nowa_ekspresja = ekspresja[1:pierwsza_kontrolna-1,]  
  nowa_ekspresja = ekspresja[1:pierwsza_kontrolna-1] 
  #nowe_nazwy_sond = row.names(nowa_ekspresja)
  
  # ekspresja proby ADENO
  indeksy_do_grupy_adeno = which(obiekt_expression_set@phenoData@data[["ADENO"]] == "ADENO", arr.ind = TRUE)
  # kontrola = obiekt_expression_set@phenoData@data[["ADENO"]][which(obiekt_expression_set@phenoData@data[["ADENO"]] == "ADENO", arr.ind = TRUE)]
  ekspresja_adeno = ekspresja[,indeksy_do_grupy_adeno[[1]]:indeksy_do_grupy_adeno[length(indeksy_do_grupy_adeno)]] # ADENO
  srednia_ekspresja_adeno = rowSums(ekspresja_adeno)/length(indeksy_do_grupy_adeno)
  srednia_ekspresja_adeno_po_filtracji_sond_kontrolnych = srednia_ekspresja_adeno[1:pierwsza_kontrolna-1]
  
  # ekspresja proby CARCINOID
  indeksy_do_grupy_carcinoid = which(obiekt_expression_set@phenoData@data[["ADENO"]] == "CARCINOID", arr.ind = TRUE)
  ekspresja_carcinoid = ekspresja[,indeksy_do_grupy_carcinoid[[1]]:indeksy_do_grupy_carcinoid[length(indeksy_do_grupy_carcinoid)]] # CARCINOID
  srednia_ekspresja_carcinoid = rowSums(ekspresja_carcinoid)/length(indeksy_do_grupy_carcinoid)
  srednia_ekspresja_carcinoid_po_filtracji_sond_kontrolnych = srednia_ekspresja_carcinoid[1:pierwsza_kontrolna-1]
  
  # ekspresja proby NORMAL
  indeksy_do_grupy_normal = which(obiekt_expression_set@phenoData@data[["ADENO"]] == "NORMAL", arr.ind = TRUE)
  ekspresja_normal = ekspresja[,indeksy_do_grupy_normal[[1]]:indeksy_do_grupy_normal[length(indeksy_do_grupy_normal)]] # NORMAL
  srednia_ekspresja_normal = rowSums(ekspresja_normal)/length(indeksy_do_grupy_normal)
  srednia_ekspresja_normal_po_filtracji_sond_kontrolnych = srednia_ekspresja_normal[1:pierwsza_kontrolna-1]
  
  # ekspresja proby SMALL CELL
  indeksy_do_grupy_smallcell = which(obiekt_expression_set@phenoData@data[["ADENO"]] == "SMALLCELL", arr.ind = TRUE)
  ekspresja_smallcell = ekspresja[,indeksy_do_grupy_smallcell[[1]]:indeksy_do_grupy_smallcell[length(indeksy_do_grupy_smallcell)]] # SMALL CELL
  srednia_ekspresja_smallcell = rowSums(ekspresja_smallcell)/length(indeksy_do_grupy_smallcell)
  srednia_ekspresja_smallcell_po_filtracji_sond_kontrolnych = srednia_ekspresja_smallcell[1:pierwsza_kontrolna-1]
  
  # ekspresja proby SQUAMOUS
  indeksy_do_grupy_squamous = which(obiekt_expression_set@phenoData@data[["ADENO"]] == "SQUAMOUS", arr.ind = TRUE)
  ekspresja_squamous = ekspresja[,indeksy_do_grupy_squamous[[1]]:indeksy_do_grupy_squamous[length(indeksy_do_grupy_squamous)]] # SQUAMOUS
  srednia_ekspresja_squamous = rowSums(ekspresja_squamous)/length(indeksy_do_grupy_squamous)
  srednia_ekspresja_squamous_po_filtracji_sond_kontrolnych = srednia_ekspresja_squamous[1:pierwsza_kontrolna-1]


  # 5% SOND O MAX/MIN EKSPRESJI ADENO
  ile_sond_to_5_procent = length(nowa_ekspresja)*0.05/100
  ile_sond_to_5_procent = floor(ile_sond_to_5_procent)
  
  kolejno_adeno = sort(srednia_ekspresja_adeno_po_filtracji_sond_kontrolnych)
  adeno_najwieksza_srednia_ekspresja_5procent = head(kolejno_adeno, ile_sond_to_5_procent)
  adeno_najmniejsza_srednia_ekspresja_5procent = tail(kolejno_adeno, ile_sond_to_5_procent)
  
  # 5% SOND O MAX/MIN EKSPRESJI CARCINOID
  kolejno_carcinoid = sort(srednia_ekspresja_carcinoid_po_filtracji_sond_kontrolnych)
  carcinoid_najwieksza_srednia_ekspresja_5procent = head(kolejno_carcinoid, ile_sond_to_5_procent)
  carcinoid_najmniejsza_srednia_ekspresja_5procent = tail(kolejno_carcinoid, ile_sond_to_5_procent)
  
  # 5% SOND O MAX/MIN EKSPRESJI NORMAL
  kolejno_normal = sort(srednia_ekspresja_normal_po_filtracji_sond_kontrolnych)
  normal_najwieksza_srednia_ekspresja_5procent = head(kolejno_normal, ile_sond_to_5_procent)
  normal_najmniejsza_srednia_ekspresja_5procent = tail(kolejno_normal, ile_sond_to_5_procent)
  
  # 5% SOND O MAX/MIN EKSPRESJI SMALL CELL
  kolejno_smallcell = sort(srednia_ekspresja_smallcell_po_filtracji_sond_kontrolnych)
  smallcell_najwieksza_srednia_ekspresja_5procent = head(kolejno_smallcell, ile_sond_to_5_procent)
  smallcell_najmniejsza_srednia_ekspresja_5procent = tail(kolejno_smallcell, ile_sond_to_5_procent)
  
  # 5% SOND O MAX/MIN EKSPRESJI SQUAMOUS
  kolejno_squamous = sort(srednia_ekspresja_squamous_po_filtracji_sond_kontrolnych)
  squamous_najwieksza_srednia_ekspresja_5procent = head(kolejno_squamous, ile_sond_to_5_procent)
  squamous_najmniejsza_srednia_ekspresja_5procent = tail(kolejno_squamous, ile_sond_to_5_procent)
  
  
  dane_wyjsciowe = list("5% sond o najmniejszej ekspresji adeno" = adeno_najmniejsza_srednia_ekspresja_5procent,
                        "5% sond o najwiêkszej ekspresji adeno" = adeno_najwieksza_srednia_ekspresja_5procent,
                        "5% sond o najmniejszej ekspresji carcinoid" = carcinoid_najmniejsza_srednia_ekspresja_5procent,
                        "5% sond o najwiêkszej ekspresji carcinoid" = carcinoid_najwieksza_srednia_ekspresja_5procent,
                        "5% sond o najmniejszej ekspresji normal" = normal_najmniejsza_srednia_ekspresja_5procent,
                        "5% sond o najwiêkszej ekspresji normal" = normal_najwieksza_srednia_ekspresja_5procent,
                        "5% sond o najmniejszej ekspresji smallcell" = smallcell_najmniejsza_srednia_ekspresja_5procent,
                        "5% sond o najwiêkszej ekspresji smallcell" = smallcell_najwieksza_srednia_ekspresja_5procent,
                        "5% sond o najmniejszej ekspresji squamous" = squamous_najmniejsza_srednia_ekspresja_5procent,
                        "5% sond o najwiêkszej ekspresji squamous" = squamous_najwieksza_srednia_ekspresja_5procent)
  
  #dane_wyjsciowe = list("5% sond o najmniejszej ekspresji" = adeno_najmniejsza_srednia_ekspresja_5procent,"5% sond o najwiekszej ekspresji" = adeno_najwieksza_srednia_ekspresja_5procent)
  return(dane_wyjsciowe)
  
}
