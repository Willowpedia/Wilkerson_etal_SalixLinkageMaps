###### INTRO ###### 

# SCRIPT 1 OF 2: MARKER QC AND CODING
# WRITTEN BY DUSTIN G WILKERSON

# FOR WILKERSON ET AL 2021; Mapping the sex determination 
#         region in eight diverse Salix F1 hybrid families

# This script is a generalized version of the code used to 
# generate the .csv files required for R/qtl and ASMap 
# for each linkage map. 

# Please direct all questions to dgw65@cornell.edu

### REQUIRED PACKAGES ### 

library(stringr)

### FUNCTIONS ### 

parentconsensus <- function(Par) {
  
  P <- c()
  
  for (R in seq(1, nrow(Par))) { 
    if (length(unique(as.character(Par[R, ]))) == 1) { # if all runs are the same
      P[R] <- unique(as.character(Par[R, ]))
    } else {
      
      gfreq <- table(t(Par[R, ]))
      
      if (ncol(Par) == 2) {
        if (max(gfreq) == 2) {  
          P[R] <- names(gfreq[which.max(gfreq)])
        } else {
          if ("N" %in% names(gfreq)) {
            if (as.logical(gfreq[which(names(gfreq) == "N")] == 1)) {
              P[R] <- names(gfreq[which(names(gfreq) != "N")])
            } else {
              P[R] <- "N" 
            }
          } else {
            P[R] <- "N" 
          }
        }
      }
      
      if (ncol(Par) == 3) {
        if (max(gfreq) >= 2 & 
            names(gfreq[which.max(gfreq)]) != "N") {
          P[R] <- names(gfreq[which.max(gfreq)])
        } else {
          if ("N" %in% names(gfreq)) {
            if (as.logical(gfreq[which(names(gfreq) == "N")] == 2)) {
              P[R] <- names(gfreq[which(names(gfreq) != "N")])
            } else {
              P[R] <- "N" 
            }
          } else {
            P[R] <- "N" 
          }
        }
      }
      
      if (ncol(Par) == 4) {
        if (max(gfreq) >= 3 & 
            names(gfreq[which.max(gfreq)]) != "N") {
          P[R] <- names(gfreq[which.max(gfreq)])
        } else {
          if ("N" %in% names(gfreq)) {
            if (as.logical(gfreq[which(names(gfreq) == "N")] >= 2)) {
              if (length(gfreq) == 2) {
                P[R] <- names(gfreq[which(names(gfreq) != "N")])
              } else {
                P[R] <- "N" 
              }
            } else {
              P[R] <- "N" 
            }
          } else {
            P[R] <- "N"
          }
        }
      }
    }
  }
  return(P)
}

OPMSorter <- function(opm) { 
  # input file must only contain markers with one parent missing.
  Sorter <- opm$MarkType
  FEM <- opm$FEM
  MAL <- opm$MAL
  
  for (snp in 1:nrow(opm)) {
    
    freq <- as.data.frame(table(as.character(opm[ snp, -c(1:14)])))
    
    if (opm[ snp, "MAL"] == "N") {# male allele is missing
      
      if (opm[snp, "FEM"] %in% c("A", "T", "G", "C")) { # known female allele is homozygous
        
        if ((freq[which.max(freq$Freq), "Var1"] == opm[snp, "FEM"]) & 
            (freq[which.max(freq$Freq), "Freq"] >= ((ncol(opm) - 14) * 0.85))) { 
          # max allele matches known female allele and there's a lot of them
          Sorter[snp] <- "HOMOMATCH"
          MAL[snp] <- FEM[snp]
        } else {# max allele does NOT match known female allele OR doesn't have the prevalence
          
          if ((!(freq[which.max(freq$Freq), "Var1"] %in% c("A", "T", "G", "C"))) & 
              (freq[which.max(freq$Freq), "Freq"] >= ((ncol(opm) - 14) * 0.85))) { 
            # max freq is heterozygous AND there's a lot of them 
            Sorter[snp] <- "HOMODIFF"
            het <- as.character(freq[which.max(freq$Freq), "Var1"])
            known <- FEM[snp]
            
            if (het == "R" & known == "A") {
              MAL[snp] <- "G"
            }
            if (het == "R" & known == "G") {
              MAL[snp] <- "A"
            }
            if (het == "Y" & known == "C") {
              MAL[snp] <- "T"
            }
            if (het == "Y" & known == "T") {
              MAL[snp] <- "C"
            }
            if (het == "S" & known == "G") {
              MAL[snp] <- "C"
            }
            if (het == "S" & known == "C") {
              MAL[snp] <- "G"
            }
            if (het == "W" & known == "A") {
              MAL[snp] <- "T"
            }
            if (het == "W" & known == "T") {
              MAL[snp] <- "A"
            }
            if (het == "K" & known == "G") {
              MAL[snp] <- "T"
            }
            if (het == "K" & known == "T") {
              MAL[snp] <- "G"
            }
            if (het == "M" & known == "A") {
              MAL[snp] <- "C"
            }
            if (het == "M" & known == "C") {
              MAL[snp] <- "A"
            }
            
          } else {# missing parent is not homozygous
            
            if (opm[snp, "FEM"] %in% freq$Var1) {
              
              if (freq[which(freq$Var1 == opm[snp, "FEM"]), "Freq"] >= ((ncol(opm) - 14) * 0.3))  {
                
                freq <- droplevels(freq[which(freq$Var1 != opm[snp, "FEM"]), ])
                
                if ((!(freq[which.max(freq$Freq), "Var1"] %in% c("A", "T", "G", "C"))) & 
                    (freq[which.max(freq$Freq), "Freq"] >= ((ncol(opm) - 14) * 0.3))) {
                  
                  Sorter[snp] <- "MALEBC"
                  MAL[snp] <- as.character(freq[which.max(freq$Freq), "Var1"])
                  
                } else {
                  Sorter[snp] <- "ROGUE"
                }
              } else {
                Sorter[snp] <- "ROGUE"
              }
            } else {
              Sorter[snp] <- "ROGUE"
            }
          }
        }
      } else {# known female allele is heterozygous
        
        if ((!(freq[which.max(freq$Freq), "Var1"] %in% c("A", "T", "G", "C"))) & 
            freq[which.max(freq$Freq), "Var1"] == opm[snp, "FEM"] & 
            freq[which.max(freq$Freq), "Freq"] >= ((ncol(opm) - 14) * 0.3)) {
          ### this can still count for both ntrx and bc markers... :(
          
          het <- opm[snp, "FEM"]
          
          if (het == "R") {
            checks <- c("A", "G")
          } else {
            if (het == "Y") {
              checks <- c("T", "C")
            } else {
              if (het == "S") {
                checks <- c("G", "C")
              } else {
                if (het == "W") {
                  checks <- c("A", "T")
                } else {
                  if (het == "K") {
                    checks <- c("G", "T")
                  } else {
                    if (het == "M") {
                      checks <- c("A", "C")
                    }
                  }
                }
              }
            }
          }
          
          if (sum(checks %in% freq$Var1) == 2) {
            
            if (freq[which(freq$Var1 == checks[1]), "Freq"] >= ((ncol(opm) - 14) * 0.1) & 
                freq[which(freq$Var1 == checks[2]), "Freq"] >= ((ncol(opm) - 14) * 0.1)) {
              
              Sorter[snp] <- "NTRX"
              MAL[snp] <- FEM[snp]
            } 
            
            if (freq[which(freq$Var1 == checks[1]), "Freq"] <= ((ncol(opm) - 14) * 0.1) & 
                freq[which(freq$Var1 == checks[2]), "Freq"] >= ((ncol(opm) - 14) * 0.3)) {
              
              Sorter[snp] <- "FEMBC"
              MAL[snp] <- checks[2]
            }
            
            if (freq[which(freq$Var1 == checks[1]), "Freq"] >= ((ncol(opm) - 14) * 0.3) & 
                freq[which(freq$Var1 == checks[2]), "Freq"] <= ((ncol(opm) - 14) * 0.1)) {
              
              Sorter[snp] <- "FEMBC"
              MAL[snp] <- checks[1]
            }
            
            if (freq[which(freq$Var1 == checks[1]), "Freq"] <= ((ncol(opm) - 14) * 0.3) & 
                freq[which(freq$Var1 == checks[2]), "Freq"] <= ((ncol(opm) - 14) * 0.3)) {
              Sorter[snp] <- "ROGUE"
            }
          }
          if (sum(checks %in% freq$Var1) == 1) {
            
            if (freq[which(freq$Var1 %in% checks), "Freq"] >= ((ncol(opm) - 14) * 0.3)) {
              
              Sorter[snp] <- "FEMBC"
              MAL[snp] <- checks[which(checks %in% freq$Var1)]
            } else {
              Sorter[snp] <- "ROGUE"
            }
          }
          if (sum(checks %in% freq$Var1) == 0) {
            
            Sorter[snp] <- "ROGUE"
          }
          
        } else {
          
          if (opm[snp, "FEM"] %in% freq$Var1) {
            
            if (freq[which(freq$Var1 == opm[snp, "FEM"]), "Freq"] >= ((ncol(opm) - 14) * 0.3)) {
              
              freq <- droplevels(freq[which(freq$Var1 != opm[snp, "FEM"]), ])
              
              if ((freq[which.max(freq$Freq), "Var1"] %in% c("A", "T", "G", "C")) & 
                  (freq[which.max(freq$Freq), "Freq"] >= ((ncol(opm) - 14) * 0.3))) {
                
                Sorter[snp] <- "FEMBC"
                MAL[snp] <- as.character(freq[which.max(freq$Freq), "Var1"])
                
              } else {
                Sorter[snp] <- "ROGUE"
              }
            } else {
              Sorter[snp] <- "ROGUE"
            }
          } else {
            Sorter[snp] <- "ROGUE"
          }
          
        }
      }
      
    } else {# female allele is missing 
      
      if (opm[snp, "MAL"] %in% c("A", "T", "G", "C")) { # known male allele is homozygous
        
        if ((freq[which.max(freq$Freq), "Var1"] == opm[snp, "MAL"]) & 
            (freq[which.max(freq$Freq), "Freq"] >= ((ncol(opm) - 14) * 0.85))) { 
          # max allele matches known female allele and there's a lot of them
          Sorter[snp] <- "HOMOMATCH"
          FEM[snp] <- MAL[snp]
        } else {# max allele does NOT match known female allele OR doesn't have the prevalence
          
          if ((!(freq[which.max(freq$Freq), "Var1"] %in% c("A", "T", "G", "C"))) & 
              (freq[which.max(freq$Freq), "Freq"] >= ((ncol(opm) - 14) * 0.85))) { 
            # max freq is heterozygous AND there's a lot of them 
            Sorter[snp] <- "HOMODIFF"
            het <- as.character(freq[which.max(freq$Freq), "Var1"])
            known <- MAL[snp]
            
            if (het == "R" & known == "A") {
              FEM[snp] <- "G"
            }
            if (het == "R" & known == "G") {
              FEM[snp] <- "A"
            }
            if (het == "Y" & known == "C") {
              FEM[snp] <- "T"
            }
            if (het == "Y" & known == "T") {
              FEM[snp] <- "C"
            }
            if (het == "S" & known == "G") {
              FEM[snp] <- "C"
            }
            if (het == "S" & known == "C") {
              FEM[snp] <- "G"
            }
            if (het == "W" & known == "A") {
              FEM[snp] <- "T"
            }
            if (het == "W" & known == "T") {
              FEM[snp] <- "A"
            }
            if (het == "K" & known == "G") {
              FEM[snp] <- "T"
            }
            if (het == "K" & known == "T") {
              FEM[snp] <- "G"
            }
            if (het == "M" & known == "A") {
              FEM[snp] <- "C"
            }
            if (het == "M" & known == "C") {
              FEM[snp] <- "A"
            }
            
          } else {# missing parent is not homozygous
            
            if (opm[snp, "MAL"] %in% freq$Var1) {
              
              if (freq[which(freq$Var1 == opm[snp, "MAL"]), "Freq"] >= ((ncol(opm) - 14) * 0.3))  {
                
                freq <- droplevels(freq[which(freq$Var1 != opm[snp, "MAL"]), ])
                
                if ((!(freq[which.max(freq$Freq), "Var1"] %in% c("A", "T", "G", "C"))) & 
                    (freq[which.max(freq$Freq), "Freq"] >= ((ncol(opm) - 14) * 0.3))) {
                  
                  Sorter[snp] <- "FEMBC"
                  FEM[snp] <- as.character(freq[which.max(freq$Freq), "Var1"])
                  
                } else {
                  Sorter[snp] <- "ROGUE"
                }
              } else {
                Sorter[snp] <- "ROGUE"
              }
            } else {
              Sorter[snp] <- "ROGUE"
            }
          }
        }
      } else {# known male allele is heterozygous
        
        if ((!(freq[which.max(freq$Freq), "Var1"] %in% c("A", "T", "G", "C"))) & 
            freq[which.max(freq$Freq), "Var1"] == opm[snp, "MAL"] & 
            freq[which.max(freq$Freq), "Freq"] >= ((ncol(opm) - 14) * 0.3)) {
          ### this can still count for both ntrx and bc markers... :(
          
          het <- opm[snp, "MAL"]
          
          if (het == "R") {
            checks <- c("A", "G")
          } else {
            if (het == "Y") {
              checks <- c("T", "C")
            } else {
              if (het == "S") {
                checks <- c("G", "C")
              } else {
                if (het == "W") {
                  checks <- c("A", "T")
                } else {
                  if (het == "K") {
                    checks <- c("G", "T")
                  } else {
                    if (het == "M") {
                      checks <- c("A", "C")
                    }
                  }
                }
              }
            }
          }
          
          if (sum(checks %in% freq$Var1) == 2) {
            
            if (freq[which(freq$Var1 == checks[1]), "Freq"] >= ((ncol(opm) - 14) * 0.1) & 
                freq[which(freq$Var1 == checks[2]), "Freq"] >= ((ncol(opm) - 14) * 0.1)) {
              
              Sorter[snp] <- "NTRX"
              FEM[snp] <- MAL[snp]
            } 
            
            if (freq[which(freq$Var1 == checks[1]), "Freq"] <= ((ncol(opm) - 14) * 0.1) & 
                freq[which(freq$Var1 == checks[2]), "Freq"] >= ((ncol(opm) - 14) * 0.3)) {
              
              Sorter[snp] <- "MALEBC"
              FEM[snp] <- checks[2]
            }
            
            if (freq[which(freq$Var1 == checks[1]), "Freq"] >= ((ncol(opm) - 14) * 0.3) & 
                freq[which(freq$Var1 == checks[2]), "Freq"] <= ((ncol(opm) - 14) * 0.1)) {
              
              Sorter[snp] <- "MALEBC"
              FEM[snp] <- checks[1]
            }
            
            if (freq[which(freq$Var1 == checks[1]), "Freq"] <= ((ncol(opm) - 14) * 0.3) & 
                freq[which(freq$Var1 == checks[2]), "Freq"] <= ((ncol(opm) - 14) * 0.3)) {
              Sorter[snp] <- "ROGUE"
            }
            
            
          }
          if (sum(checks %in% freq$Var1) == 1) {
            
            if (freq[which(freq$Var1 %in% checks), "Freq"] >= ((ncol(opm) - 14) * 0.3)) {
              
              Sorter[snp] <- "MALEBC"
              FEM[snp] <- checks[which(checks %in% freq$Var1)]
            } else {
              Sorter[snp] <- "ROGUE"
            }
          }
          if (sum(checks %in% freq$Var1) == 0) {
            
            Sorter[snp] <- "ROGUE"
          }
        } else {
          
          if (opm[snp, "MAL"] %in% freq$Var1) {
            
            if (freq[which(freq$Var1 == opm[snp, "MAL"]), "Freq"] >= ((ncol(opm) - 14) * 0.3)) {
              
              freq <- droplevels(freq[which(freq$Var1 != opm[snp, "MAL"]), ])
              
              if ((freq[which.max(freq$Freq), "Var1"] %in% c("A", "T", "G", "C")) & 
                  (freq[which.max(freq$Freq), "Freq"] >= ((ncol(opm) - 14) * 0.3))) {
                
                Sorter[snp] <- "MALEBC"
                FEM[snp] <- as.character(freq[which.max(freq$Freq), "Var1"])
                
              } else {
                Sorter[snp] <- "ROGUE"
              }
            } else {
              Sorter[snp] <- "ROGUE"
            }
          } else {
            Sorter[snp] <- "ROGUE"
          }
          
        }
      }
      
    }
  }
  
  opm$MarkType <- Sorter
  opm$FEM <- FEM
  opm$MAL <- MAL
  
  return(opm)
}

MakeERRsGoAway <- function(df){
  
  codings <- lapply(1:nrow(df), function(x){
    
    alleles <- as.character(df[ x , c(15:ncol(df))])
    
    if (df[ x, "MarkType"] %in% c("FEMBC", "MALEBC")) {
      ifelse(alleles %in% c(df[ x , "MAL"], df[ x, "FEM"]), alleles, "N")
    } else {
      
      if (df[ x, "MarkType"] == "NTRX") {
        
        if (df[ x, "FEM"] == "W") {
          ifelse(alleles %in% c("W", "A", "T"), alleles, "N")
        } else {
          if (df[ x, "FEM"] == "R") {
            ifelse(alleles %in% c("R", "A", "G"), alleles, "N")
          } else {
            if (df[ x, "FEM"] == "Y") {
              ifelse(alleles %in% c("Y", "C", "T"), alleles, "N")
            } else {
              if (df[ x, "FEM"] == "S") {
                ifelse(alleles %in% c("S", "C", "G"), alleles, "N")
              } else {
                if (df[ x, "FEM"] == "K") {
                  ifelse(alleles %in% c("K", "G", "T"), alleles, "N")
                } else {
                  if (df[ x, "FEM"] == "M") {
                    ifelse(alleles %in% c("M", "A", "C"), alleles, "N")
                  }}}}}}
      } else {
        
        if (df[ x, "MarkType"] == "HOMOMATCH") {
          ifelse(alleles == df[ x , "FEM"], alleles, "N")
        } else {
          
          if (df[ x, "MarkType"] == "HOMODIFF") {
            if (df[ x, "FEM"] == "A") {
              if (df[ x, "MAL"] == "T") {
                ifelse(alleles == "W", alleles, "N")
              } else {
                if (df[ x, "MAL"] == "C") {
                  ifelse(alleles == "M", alleles, "N")
                } else {
                  if (df[ x, "MAL"] == "G") {
                    ifelse(alleles == "R", alleles, "N")
                  }}}
            } else {
              if (df[ x, "FEM"] == "T") {
                if (df[ x, "MAL"] == "A") {
                  ifelse(alleles == "W", alleles, "N")
                } else {
                  if (df[ x, "MAL"] == "C") {
                    ifelse(alleles == "Y", alleles, "N")
                  } else {
                    if (df[ x, "MAL"] == "G") {
                      ifelse(alleles == "K", alleles, "N")
                    }}}
              } else {
                if (df[ x, "FEM"] == "C") {
                  if (df[ x, "MAL"] == "A") {
                    ifelse(alleles == "M", alleles, "N")
                  } else {
                    if (df[ x, "MAL"] == "T") {
                      ifelse(alleles == "Y", alleles, "N")
                    } else {
                      if (df[ x, "MAL"] == "G") {
                        ifelse(alleles == "S", alleles, "N")
                      }}}
                } else {
                  if (df[ x, "FEM"] == "G") {
                    if (df[ x, "MAL"] == "A") {
                      ifelse(alleles == "R", alleles, "N")
                    } else {
                      if (df[ x, "MAL"] == "T") {
                        ifelse(alleles == "K", alleles, "N")
                      } else {
                        if (df[ x, "MAL"] == "C") {
                          ifelse(alleles == "S", alleles, "N")
                        }}}
                  }
                }
              }
            }
          }
        }
      }
    }
    
  })
  
  codings <- as.data.frame(do.call("rbind", codings), stringsAsFactors = F)
  colnames(codings) <- colnames(df)[c(15:ncol(df))]
  colnames(codings) <- toupper(str_replace_all(str_sub(colnames(codings), end = 12), "[[:punct:]]", "-"))
  codings <- cbind(df[ , c(1:14)], codings)
  
  return(codings)
}

InspectMarkers <- function(df, MarkType) {

    split <- data.frame("RS" = as.character(df$rs),
                        "x0" = rep(0, nrow(df)),
                        "x1" = rep(0, nrow(df)),
                        "x2" = rep(0, nrow(df)),
                        "x9" = rep(0, nrow(df)))
  
  for (x in df$rs) {
    
    alleles <- as.character(df[which(df$rs == x) , c(15:ncol(df))])
    
    if (MarkType == "FEMBC") {
      
      alleles <- ifelse(alleles != df[which(df$rs == x) , "MAL"] & alleles != df[ which(df$rs == x), "FEM"], "N", alleles)
      freq <- table(ifelse(alleles == df[which(df$rs == x) , "FEM"], 1, ifelse(alleles == df[which(df$rs == x) , "MAL"], 0, 9)))
      if ("0" %in% names(freq)) {
        split[which(split$RS == x) , "x0"] <- as.numeric(freq[which(names(freq) == "0")])
      }  
      if ("1" %in% names(freq)) {
        split[which(split$RS == x) , "x1"] <- as.numeric(freq[which(names(freq) == "1")])
      }
      if ("9" %in% names(freq)) {
        split[which(split$RS == x) , "x9"] <- as.numeric(freq[which(names(freq) == "9")])
      }
    } 
    
    if (MarkType == "MALEBC") {
      
      alleles <- ifelse(alleles != df[which(df$rs == x) , "MAL"] & alleles != df[ which(df$rs == x), "FEM"], "N", alleles)
      freq <- table(ifelse(alleles == df[which(df$rs == x) , "FEM"], 0, ifelse(alleles == df[which(df$rs == x), "MAL"], 1, 9)))
      if ("0" %in% names(freq)) {
        split[which(split$RS == x) , "x0"] <- as.numeric(freq[which(names(freq) == "0")])
      }
      if ("1" %in% names(freq)) {
        split[which(split$RS == x) , "x1"] <- as.numeric(freq[which(names(freq) == "1")])
      }
      if ("9" %in% names(freq)) {
        split[which(split$RS == x) , "x9"] <- as.numeric(freq[which(names(freq) == "9")])
      }
    } 
    
    if (MarkType == "NTRX") {
      
      het <- df[which(df$rs == x), "FEM"]
      
      if (het == "W") {
        alleles <- ifelse(alleles %in% c("W", "A", "T"), alleles, "N")
      }
      if (het == "R") {
        alleles <- ifelse(alleles %in% c("R", "A", "G"), alleles, "N")
      }
      if (het == "Y") {
        alleles <- ifelse(alleles %in% c("Y", "C", "T"), alleles, "N")
      }
      if (het == "S") {
        alleles <- ifelse(alleles %in% c("S", "C", "G"), alleles, "N")
      }
      if (het == "K") {
        alleles <- ifelse(alleles %in% c("K", "G", "T"), alleles, "N")
      }
      if (het == "M") {
        alleles <- ifelse(alleles %in% c("M", "A", "C"), alleles, "N")
      }
      
      freq <- table(alleles)[names(table(alleles)) %in% c("A", "T", "C", "G")]
      
      if (length(freq) <= 1) {
        split[which(split$RS == x), 2:5] <- as.numeric(404, 404, 404, 404)
      } else {
        
        if (freq[1] == freq[2]) {
          freq <- table(ifelse(alleles == names(freq[1]), 2, 
                               ifelse(alleles == names(freq[2]), 0,
                                      ifelse(alleles == het, 1, 9))))
        } else {
          freq <- table(ifelse(alleles == names(freq[which.min(freq)]), 2,
                               ifelse(alleles == names(freq[which.max(freq)]), 0,
                                      ifelse(alleles == het, 1, 9))))
        }
        
        if ("0" %in% names(freq)) {
          split[which(split$RS == x) , "x0"] <- as.numeric(freq[which(names(freq) == "0")])
        }
        if ("1" %in% names(freq)) {
          split[which(split$RS == x) , "x1"] <- as.numeric(freq[which(names(freq) == "1")])
        }
        if ("2" %in% names(freq)) {
          split[which(split$RS == x) , "x2"] <- as.numeric(freq[which(names(freq) == "2")])
        }
        if ("9" %in% names(freq)) {
          split[which(split$RS == x) , "x9"] <- as.numeric(freq[which(names(freq) == "9")])
        }
      }
    }
    
  }
  if (MarkType == "FEMBC" | MarkType == "MALEBC") {
    
    split$p <- .5
    split$q <- .5
    split$E0 <- .5*rowSums(split[ , c(2, 3, 5)])
    split$E1 <- .5*rowSums(split[ , c(2, 3, 5)])
    split$E2 <- 0
    split$Chi <- ((split$x0 - split$E0)^2/split$E0) + ((split$x1 - split$E1)^2/split$E1)
    split$pv <- 1 - pchisq(split$Chi, df = 1)
    
  } else {
    
    split$p <- ((split$x0*2) + split$x1)/(2*rowSums(split[ , c(2:5)]))
    split$q <- ((split$x2*2) + split$x1)/(2*rowSums(split[ , c(2:5)]))
    split$E0 <- (split$p*split$p)*rowSums(split[ , c(2:5)])
    split$E1 <- (2*split$p*split$q)*rowSums(split[ , c(2:5)])
    split$E2 <- (split$q*split$q)*rowSums(split[ , c(2:5)])
    split$Chi <- ((split$x0 - split$E0)^2/split$E0) + ((split$x1 - split$E1)^2/split$E1) + ((split$x2 - split$E2)^2/split$E2)
    split$pv <- 1 - pchisq(split$Chi, df = 1)
  }
  return(split)
}

AlleleCoding <- function(male_df, female_df, firstprog) {
  for (W in seq(1, nrow(female_df))) {
    P1 <- female_df[W, which(colnames(female_df) == "FEM")]
    P2 <- female_df[W, which(colnames(female_df) == "MAL")]
    for (C in seq(firstprog, ncol(female_df))) {
      female_df[W, C] <- ifelse(female_df[W, C] == P1, "AB", ifelse(female_df[W, C] == P2, "AA", NA))
    }
  }
  for (M in seq(1, nrow(male_df))) {
    P1 <- male_df[M, which(colnames(male_df) == "FEM")]
    P2 <- male_df[M, which(colnames(male_df) == "MAL")]
    for (C in seq(firstprog, ncol(male_df))) {
      male_df[M, C] <- ifelse(male_df[M, C] == P2, "AB", ifelse(male_df[M, C] == P1, "AA", NA))
    }  
  }
  return(list(female_df, male_df))
}

rQTLFormatter <- function(mal, fem) {
  
  insertrow <- function(dataframe, newrow, location) {
    
    dataframe[seq(location + 1, nrow(dataframe) + 1), ] <- dataframe[seq(location, nrow(dataframe)), ]
    dataframe[ location, ] <- newrow
    
    return(dataframe)
  }
  
  mal_rs <- mal[ , 1]
  fem_rs <- fem[ , 1]
  mal_cp <- mal[ , 3:4]
  fem_cp <- fem[ , 3:4]
  femGenotype <- colnames(fem)[-c(1:13)]
  femGenotype <- sub("X", "", toupper(gsub("[[:punct:]]", "-", femGenotype)))
  femGenotype <- c(NA, NA, femGenotype)
  malGenotype <- colnames(mal)[-c(1:13)]
  malGenotype <- sub("X", "", toupper(gsub("[[:punct:]]", "-", malGenotype)))
  malGenotype <- c(NA, NA, malGenotype)
  mal <- mal[ , 14:ncol(mal)]
  fem <- fem[ , 14:ncol(fem)]
  mal <- t(mal)
  fem <- t(fem)
  maldf <- as.data.frame(mal, row.names = F, stringsAsFactors = F)
  femdf <- as.data.frame(fem, row.names = F, stringsAsFactors = F)
  mal <- insertrow(insertrow(maldf, mal_cp$pos, 1), mal_cp$chrom, 1)
  fem <- insertrow(insertrow(femdf, fem_cp$pos, 1), fem_cp$chrom, 1)
  colnames(mal) <- mal_rs
  colnames(fem) <- fem_rs
  mal <- cbind("GENOTYPE" = malGenotype, mal)
  fem <- cbind("GENOTYPE" = femGenotype, fem)
  mal$GENOTYPE <- as.character(mal$GENOTYPE)
  fem$GENOTYPE <- as.character(fem$GENOTYPE)
  mal[1:2, 1] <- ""
  fem[1:2, 1] <- ""
  
  return(list("MAL" = mal, "FEM" = fem))
}

###### DERIVING PARENTAL CONSENSUS FOR GENOTYPE CODING ###### 

# read in a hapmap text file for a family with the F1 and multiple parents # 
F1 <- read.table("", header = T, sep = "\t", na.strings = "NA", stringsAsFactors = F) 

# pull the parents out into a separate file # 
parents <- c("94006", "94001", ".BN.", ".FF.", "MBG", "P336", "P295", "P294", "P63", "JORR")
parents <- droplevels(F1[ , grepl(paste(parents, collapse = "|"), colnames(F1))])

# remove the parents from the F1 file # 
F1noPar <- F1[ , (!(colnames(F1) %in% colnames(parents)))] 

# split the parent file into two # 
parents <- c("94006", "94001", ".BN.", ".FF.", "MBG", "P336", "P295", "P294", "P63", "JORR" )
FEM <- droplevels(F1[ , grepl(parents[], colnames(F1))])
MAL <- droplevels(F1[ , grepl(parents[], colnames(F1))])

# use the multiple runs of the parentts to derive their consensus genotypes # 
FEM <- parentconsensus(FEM)
MAL <- parentconsensus(MAL)

# reconstruct the family file, now with the consensus parent genotypes # 
F1 <- cbind(F1noPar[, 1:11], 
            "FEM" = FEM, "MAL" = MAL,
            F1noPar[, 12:ncol(F1noPar)],
            stringsAsFactors = F)

# Identify the marker types based on the parental consensus genotypes # 
F1$MarkType <- ifelse(((F1$FEM == F1$MAL) & 
                         (F1$FEM != "N" & F1$FEM != "A" & F1$FEM != "T" & F1$FEM != "C" & F1$FEM != "G")), "NTRX", 
                       ifelse(((F1$FEM != F1$MAL) & 
                                 (F1$FEM != "N" & F1$FEM != "A" & F1$FEM != "T" & F1$FEM != "C" & F1$FEM != "G") &
                                 (F1$MAL != "N") & 
                                 (F1$MAL == "A" | F1$MAL == "T" | F1$MAL == "C" | F1$MAL == "G")), "FEMBC",
                              ifelse(((F1$FEM != F1$MAL) &
                                        (F1$FEM != "N") & 
                                        (F1$FEM == "A" | F1$FEM == "T" | F1$FEM == "C" | F1$FEM == "G") &
                                        (F1$MAL != "N" & F1$MAL != "A" & F1$MAL != "T" & F1$MAL != "C" & F1$MAL != "G")), "MALEBC",
                                     "ROGUE"))) 

# Sort through the "ROGUE" type to pull out markers with one parent missing # 
F1$MarkType <- ifelse(F1$MarkType != "ROGUE", F1$MarkType,
                      ifelse(F1$MarkType == "ROGUE" &
                               (F1$FEM == "N" & F1$MAL != "N"), "FEMMISS", 
                             ifelse(F1$MarkType == "ROGUE" &
                                      (F1$FEM != "N" & F1$MAL == "N"), "MALMISS", "ROGUE")))

# Move MarkType column to before marker data # 
F1 <- F1[ , c(1:13, ncol(F1), 14:(ncol(F1) - 1))]

# Pull out the OneParentMissing markers from the others # 
OPM <- droplevels(F1[which(F1$MarkType %in% c("FEMMISS", "MALMISS")), ]) # 59901
F1 <- droplevels(F1[which(!(F1$MarkType %in% c("FEMMISS", "MALMISS"))), ]) # 44986

# Use the known parent and F1 genotypes to infer the missing parental genotype # 
OPM <- OPMSorter(OPM)

# Recombine the marker data # 
F1 <- rbind(F1, OPM) # 84913
F1 <- F1[ order(F1$chrom, F1$pos), ]

# Remove markers aligning to unmapped scaffolds in the reference # 
F1 <- droplevels(F1[which(F1$chrom <= 19), ])

# Keep only backcross markers # 
F1 <- droplevels(F1[which(F1$MarkType %in% c("MALEBC", "FEMBC")), ]) # 29403

# Use parental genotypes to set F1 genotypes to missing based on expected segregation
F1 <- MakeERRsGoAway(F1)

# Remove markers with more than 20% missing data # 
count <- lapply(seq.int(1, nrow(F1)), function(x) {
  sum(F1[x , -c(1:14)] == "N") })
sum(count <= ceiling(((ncol(F1) - 14)) * 0.2))
F1 <- droplevels(F1[which(count <= ceiling(((ncol(F1) - 14)) * 0.2)), ])

###### CHECKING UP ON SEGREGATION DISTORTION ###### 

# Split the markers into male and female informative maps 
FEM <- droplevels(F1[which(F1$MarkType == "FEMBC"), ]) 
MAL <- droplevels(F1[which(F1$MarkType == "MALEBC"), ]) 

# Use marker type to perform Chi-Square test for segregation distortion # 
FEM <- InspectMarkers(FEM, "FEMBC")
MAL <- InspectMarkers(MAL, "MALEBC")

# Remove markers significant for segregation distortion # 
FEM$MarkType <- "FEMBC"
MAL$MarkType <- "MALEBC"
MarkDat <- rbind(FEM, MAL)
MarkDat <- MarkDat[order(as.character(MarkDat$RS)), ] 
MarkDat <- droplevels(MarkDat[which(MarkDat$pv > 0), ])
F1 <- droplevels(F1[which(F1$rs %in% MarkDat$RS), ])

###### CODING and FORMATTING ###### 

# Split the male and female markers again # 
FEM <- droplevels(F1[which(F1$MarkType == "FEMBC"), ])
MAL <- droplevels(F1[which(F1$MarkType == "MALEBC"), ])

# Remove marker type column # 
FEM$MarkType <- NULL
MAL$MarkType <- NULL

# Code Alleles # 
catch <- AlleleCoding(MAL, FEM, 14)
FEM <- as.data.frame(catch[1])
MAL <- as.data.frame(catch[2])

# Format for R/qtl # 
catch <- rQTLFormatter(MAL, FEM)
FEM <- catch$FEM
MAL <- catch$MAL

# Save as .csv # 
write.table(FEM, ".csv", quote = F, sep = ",", row.names = F, col.names = T)
write.table(MAL, ".csv", quote = F, sep = ",", row.names = F, col.names = T)
