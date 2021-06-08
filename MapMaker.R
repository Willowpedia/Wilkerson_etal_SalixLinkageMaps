###### INTRO ###### 

# SCRIPT 2 OF 2: MAP MAKER 
# WRITTEN BY DUSTIN G WILKERSON

# FOR WILKERSON ET AL 2021; Mapping the sex determination 
#         region in eight diverse Salix F1 hybrid families

# This script is a generalized version of the code used to 
# create linkage maps using R/qtl and ASMap and the  
# functions defined below. Script is based on the creation of the 
# female parent map (FEM). Checkpoints were used to save the map at key
# stages in case of human error. Do not follow script blindly, each 
# population and therefore linkage maps will be different. Use the diagnostic 
# tools available in R/qtl and ASMap to make the best decision as you go. 

# Please direct all questions to dgw65@cornell.edu

### REQUIRED PACKAGES ###

library(qtl)
library(ASMap)
library(stringr)

### FUNCTIONS ### 

chrIDbymarks <- function(cross_obj, chrlist, Just19 = T) { 
# uses marker names to generate a list based on where the markers mapped to P294. 
# compares this list to the formed linkage groups to show if any chromosomes have
# grouped together
	
  lgnames <- c()
  snpnames <- c()
  
	for (chr in chrlist) {
		names <- markernames(cross_obj, chr = chr)
		names <- sub("S", "", word(sub("_", " ", names), 1))
		names <- as.numeric(names)
		if (Just19 == T) {
		  names <- names[which(names <= 19)]
		}
		
		snpnames <- append(snpnames, names)
		lgnames <- append(lgnames, as.numeric(rep(chr, length(names))))
	}
  
  
	return(table(snpnames, lgnames))
}

SwitchFinder <- function(df, lod) {
#locates 'switched' alleles via LOD scores above a set threshold "lod" by LG
#returns: G = list of LGs that show high linkage but low Recomb Freq
#		      M = list of LG locations of each marker within grouped LGs
#		      N = list of marker names within each grouped LGs
	names(df$geno) <- seq(1, length(df$geno))
	df <- est.rf(df)
	rf <- pull.rf(df, "lod")
	groups <- list()
	mlist <- list()
	namelist <- list()
	for (lgs in names(df$geno)) {
		rft <- subset(rf[rownames(rf) %in% markernames(df, lgs), ]) 
 		sigmarks <- c()
			for (C in seq(1, length(colnames(rft)))) {
				if (sum(rft[ , C] > lod, na.rm = T) > 0) {
					sigmarks <- append(sigmarks, C)
				}
			}
 		sigmarks <- colnames(rft[ , sigmarks])
		sigLGs <- c()
		for (M in sigmarks) {
      		for (G in seq(1, length(df$geno))) {
       			if (M %in% names(df$geno[[G]]$map)) {
					sigLGs <- append(sigLGs, G)
        			}
      		}
		}
 	groups <- append(groups, list(table(sigLGs)))
 	mlist <- append(mlist, list(sigLGs))
	namelist <- append(namelist, list(sigmarks))
  	}
return(list("G" = groups, "M" = mlist, "N" = namelist))
}

MessyMarks <- function(df_M, df_N) {
  # moves through LG groups from SwitchFinder and pulls out small, unlikely groups of 
  #  markers hiding with larger LGs,
  # need to find them as they likely carry extraneous error we don't need
  droppers <- c()
  for (G in seq(1, length(df_M))) {
    
    freq <- data.frame(table(df_M[[G]]))
    freq <- freq[order(-freq$Freq), ]
    ids <- df_N[[G]][which(df_M[[G]] %in% freq$Var1[-c(1:2)])]
    droppers <- append(droppers, ids)
    
  }
  return(droppers)
}

FindSmallerLGs <- function(df, df_M) {
  
  mins <- c()
  lgsizes <- unlist(lapply(1:length(df_M), function(x){length(df$geno[[x]]$map)}))
  
  for (G in seq(1, length(df_M))) {
    
    lkd <- as.data.frame(table(df_M[[G]]))
    lkd <- as.numeric(as.character(lkd[order(-lkd$Freq), ]$Var1[1:2]))
    
    min <- as.character(lkd[which.min(lgsizes[lkd])])
    
    mins <- append(mins, min)
  }
  return(mins)
}

MarkerOCD <- function(df) {
  
  output <- c()
  
  for (lg in names(df$geno)) {
    
    marks <- markernames(df, lg)
    chrs <- word(markernames(df, lg), start = 1, sep = "_")
    
    freq <- as.data.frame(table(chrs))
    
    toofew <- as.character(freq$chrs[which(freq$Freq < 5)])
    enough <- as.character(freq$chrs[which(freq$Freq >= 5)])
    
    droppers <- marks[which(chrs %in% toofew)]
    enough <- enough[which(enough != as.character(freq$chrs[which.max(freq$Freq)]))]
    
    
    if (length(enough) > 0) { 
      # if there were chrs with more than 5 markers w/n this LG
      
      for (grp in enough) { 
        # for each chr group
        pos <- which(chrs == grp) # creates a position list of all the markers from the chr grp
        ct <- 1 # starts a loop counter so I can reference other positions in the list
        streak <- c() # when positions are neighbors, they're stored here
        cat("STARTING:", lg, ":", grp, ";", length(pos), "POSITIONS", "\n", sep = "")
        
        for (loc in pos) {
          cat()
          if (!is.na(pos[ct + 1])) { 
            # checks whether or not this position is the last in the list
            
            if (loc + 1 != pos[ct + 1]) {
              # if this position + 1 is not the next position
              # means that either this snp is the last in a streak or is not
              # part of one
              
              if (loc == pos[1]) {
                # loc is the first position and
                # isn't the beginning of a streak
                cat("FIRST MARKER BAD", "\n")
                droppers <- append(droppers, marks[loc])
                
              } else {
                # loc is not the first position 
                # and can be checked if its the last
                # in a streak 
                
                if (loc - 1 == pos[ct - 1]) { 
                  # if this position - 1 equals the prior position 
                  # indicates this is the last SNP in a streak
                  # and appends this location to the streak list
                  cat("STREAK + 1", "\n")
                  streak <- append(streak, loc)
                  
                  if (length(streak) >= 5) {
                    # since this is the last loc in a streak
                    # checks the length of the streak
                    # if it is long enough, the streak list is reset 
                    streak <- c()
                    cat("STREAK ACCEPTED", "\n")
                    
                  } else {
                    # if its not long enough
                    # the streak locs are added to the drop list
                    # and the steak list is reset 
                    cat("STREAK TOO SHORT: ", length(streak), "\n")
                    
                    droppers <- append(droppers, marks[streak])
                    streak <- c()
                  }
                  
                } else {
                  # this loc is not part of a streak 
                  # added to the drop list
                  cat("REJECTED", "\n")
                  droppers <- append(droppers, marks[loc])
                }
              }
              
            } else {
              # this loc is part of a streak and is not the last one
              cat("STREAK + 1", "\n")
              streak <- append(streak, loc)
            } 
          } else {
            # loc is the last position
            
            if (loc - 1 == pos[ct - 1]) { 
              # if this position - 1 equals the prior position 
              # indicates this is the last SNP in a streak
              # and appends this location to the streak list
              cat("STREAK + 1", "\n")
              streak <- append(streak, loc)
              
              if (length(streak) >= 10) { 
                # since this is the last loc in a streak
                # checks the length of the streak
                # if it is long enough, the streak list is reset 
                cat("STREAK ACCEPTED", "\n")
                streak <- c()
                
              } else {
                # if its not long enough
                # the streak locs are added to the drop list
                # and the steak list is reset 
                cat("STREAK TOO SHORT: ", length(streak), "\n")
                droppers <- append(droppers, marks[streak])
                streak <- c()
              }
              
            } else {
              # this loc is not part of a streak 
              # added to the drop list
              cat("REJECTED", "\n")
              droppers <- append(droppers, marks[loc])
            }
          }
          ct <- ct + 1
        }
        cat("NUMBER DROPPED:", length(droppers), "\n")
        
      }
    }
    
    output <- append(output, droppers)
    cat("OUTPUT TOTAL: ", length(output), "\n")
    
    
  }
  
  return(output)
}

getcM <- function(cross_obj) {
	nexus <- list()
	counter <- 1	
	for (chr in chrnames(cross_obj)) {
		catch <- as.data.frame(cross_obj$geno[[chr]]$map)
		colnames(catch) <- c("cMdist")
		catch$physdist <- word(sub("_", " ", rownames(catch)), 2)
		catch$physdist <- as.numeric(sub("I", "", catch$physdist))
		catch$physdist <- catch$physdist/1000000 
		nexus[[counter]] <- catch
		counter <- counter + 1
	}
return(nexus)
}

XODXOCorrector <- function(df) {
  
  df$GENOTYPE[1:2] <- c("blank", "blank2")
  rownames(df) <- df$GENOTYPE
  
  corrected <- data.frame("GENOTYPE" = df$GENOTYPE, row.names = df$GENOTYPE)
  corrected$GENOTYPE[1:2] <- ""
  
  markorder <- data.frame("MARKNAME" = colnames(df)[-1], 
                          "LG" = as.numeric(df[1, ][-1]),
                          "CHR" = as.numeric(str_replace(word(colnames(df)[-1], 1, sep = "_"), "S", "")),
                          "POS" = as.numeric(word(colnames(df)[-1], 2, sep = "_")),
                          stringsAsFactors = F)
  
  if (sum(markorder$LG != markorder$CHR) != 0) {
    pullouts <- markorder[which(markorder$LG != markorder$CHR), ]
  } else {
    pullouts <- NULL
  }
  
  markorder <- markorder[order(markorder$LG, markorder$POS), ]
  df$GENOTYPE <- NULL
  
  pp <- df[ , match(markorder$MARKNAME, colnames(df))]
  
  for (lg in unique(as.character(pp[ 1, ]))) {
    grp <- droplevels(pp[which(pp[ 1, ] == lg)])
    nocor <- 0
    AvgLD <- c()
    
    for (geno in seq.int(3, nrow(grp))) {
      GenoLD <- c()
      LDLength <- 1
      
      for (pos in seq.int(2, (ncol(grp) - 1))) {
        
        #cat(paste("\n", pos, "\n", sep = ""))
        
        if (grp[geno, pos] == grp[geno, pos - 1] & 
            grp[geno, pos] == grp[geno, pos + 1]) { 
          # snp behind pos and snp in front of pos 
          # are the same as pos
          
          LDLength <- LDLength + 1
          #cat(paste("\n", "FLANK + 1", "\n", sep = ""))
          
          if (is.null(grp[geno, pos + 2])) { 
            LDLength <- LDLength + 1
            GenoLD <- append(GenoLD, LDLength)
            #cat(paste("\n", " LG END: ", LDLength, "\n", sep = ""))
          }
          
        } else {# flanking markers aren't the same 
          
          if (grp[geno, pos] != grp[geno, pos - 1] & 
              grp[geno, pos] != grp[geno, pos + 1]) { 
            # snp behind pos and snp in front of pos
            # are both different than pos; ie single SNP DXO
            
            grp[geno, pos] <- grp[geno, pos - 1]
            nocor <- nocor + 1
            LDLength <- LDLength + 1
            #cat(paste("\n", "DXO + 1", "\n", sep = ""))
            
            if (is.null(grp[geno, pos + 2])) { 
              LDLength <- LDLength + 1
              GenoLD <- append(GenoLD, LDLength)
              #cat(paste("\n", " LG END: ", LDLength, "\n", sep = ""))
            }
            
          } else {# not a single SNP DXO
            
            if ((!(is.null(grp[geno, pos + 3]))) & 
                (!(is.null(grp[geno, pos + 2])))) {
              # pos is not in relavent closeness to the end
              
              if (grp[geno, pos] == grp[geno, pos - 1] & 
                  grp[geno, pos] != grp[geno, pos + 1] & 
                  grp[geno, pos] != grp[geno, pos + 2] &
                  grp[geno, pos] != grp[geno, pos + 3]) {
                # snp behind is same
                # snp in front is different
                # snp 2 in front is different
                # end of LD block
                GenoLD <- append(GenoLD, LDLength)
                #cat(paste("\n", "LD Block End: ", LDLength, "\n", sep = ""))
                
              } 
              
              if (grp[geno, pos] != grp[geno, pos - 1] & 
                  grp[geno, pos] == grp[geno, pos + 1] & 
                  grp[geno, pos] == grp[geno, pos + 2]) {
                # snp behind is different
                # snp in front and snp 2 in front is same
                # start of new LD block
                
                LDLength <- 1
                #cat(paste("\n", "LD Block Begin: ", LDLength, "\n", sep = ""))
                
              }  
              
              if (grp[geno, pos] == grp[geno, pos - 1] & 
                  grp[geno, pos] != grp[geno, pos + 1] & 
                  grp[geno, pos] == grp[geno, pos + 2]) {
                # snp behind is same
                # snp in front is different
                # snp 2 in front is the same
                # right in front of a DXO
                LDLength <- LDLength + 1
                #cat(paste("\n", "Pre DXO + 1", "\n", sep = ""))
                
                
              }
              
              if (grp[geno, pos] == grp[geno, pos - 1] & 
                  grp[geno, pos] != grp[geno, pos + 1] & 
                  grp[geno, pos] != grp[geno, pos + 2] & 
                  grp[geno, pos] == grp[geno, pos + 3]) {
                # snp behind is same
                # snp in front is different
                # snp 2 in front is different
                # snp 3 in front in same
                # right in front two snp XO
                LDLength <- LDLength + 1
                #cat(paste("\n", "Pre Two SNP XO + 1", "\n", sep = ""))
                
              }
              
              if (grp[geno, pos] != grp[geno, pos - 1] & 
                  grp[geno, pos] == grp[geno, pos + 1] & 
                  grp[geno, pos] != grp[geno, pos + 2]) {
                # snp behind is different
                # snp in front is the same
                # snp 2 in front is different
                # xo of only 2 snps
                
                grp[geno, pos] <- grp[geno, pos - 1]
                nocor <- nocor + 1
                LDLength <- LDLength + 1
                
                #cat(paste("\n", "Two SNP XO + 1", "\n", sep = ""))
                
              }
            }
            
            if (is.null(grp[geno, pos + 3]) & 
                (!(is.null(grp[geno, pos + 2])))) {
              # pos is the next to last snp in the lg
              
              if (grp[geno, pos] != grp[geno, pos - 1] & 
                  grp[geno, pos] == grp[geno, pos + 1] & 
                  grp[geno, pos] == grp[geno, pos + 2]) {
                # snp behind is different
                # snp in front and snp 2 in front is same
                # start of new LD block
                
                LDLength <- 1
                #cat(paste("\n", "LD Block Begin: ", LDLength, "\n", sep = ""))
                
              }  
              
              if (grp[geno, pos] != grp[geno, pos - 1] & 
                  grp[geno, pos] == grp[geno, pos + 1] & 
                  grp[geno, pos] != grp[geno, pos + 2]) {
                # snp behind is different
                # snp in front is the same
                # snp 2 in front is different
                # xo of only 2 snps
                
                grp[geno, pos] <- grp[geno, pos - 1]
                nocor <- nocor + 1
                LDLength <- LDLength + 1
                
                #cat(paste("\n", "Two SNP XO + 1", "\n", sep = ""))
                
              }
              
              if (grp[geno, pos] == grp[geno, pos - 1] & 
                  grp[geno, pos] != grp[geno, pos + 1] & 
                  grp[geno, pos] == grp[geno, pos + 2]) {
                # snp behind is same
                # snp in front is different
                # snp 2 in front is the same
                # right in front of a DXO
                LDLength <- LDLength + 1
                #cat(paste("\n", "Pre DXO + 1", "\n", sep = ""))
                
                
              }
              
              if (grp[geno, pos] == grp[geno, pos - 1] & 
                  grp[geno, pos] != grp[geno, pos + 1] & 
                  grp[geno, pos] != grp[geno, pos + 2]) {
                # snp behind is same
                # snp in front is different
                # snp 2 in front is different
                # right in front of two snp XO
                LDLength <- LDLength + 1
                #cat(paste("\n", "Pre Two SNP XO + 1", "\n", sep = ""))
                
                
              }
              
            }
            
            if (is.null(grp[geno, pos + 3]) & 
                is.null(grp[geno, pos + 2])) {
              # pos is the next to last snp in the lg
              
              if (grp[geno, pos] == grp[geno, pos - 1] & 
                  grp[geno, pos] != grp[geno, pos + 1]) {
                
                grp[geno, pos + 1] <- grp[geno, pos]
                nocor <- nocor + 1
                LDLength <- LDLength + 2
                #cat(paste("\n", " LG END + 2: ", LDLength, "\n", sep = ""))
                GenoLD <- append(GenoLD, LDLength)
              }
              
              if (grp[geno, pos] != grp[geno, pos - 1] & 
                  grp[geno, pos] == grp[geno, pos + 1]) {
                
                grp[geno, pos] <- grp[geno, pos - 1]
                grp[geno, pos + 1] <- grp[geno, pos - 1]
                
                nocor <- nocor + 2
                LDLength <- LDLength + 2
                #cat(paste("\n", "LG END + 2: ",LDLength, "\n", sep = ""))
                GenoLD <- append(GenoLD, LDLength)
              }
            }
            
          }
        }
      }
      
      AvgLD <- append(AvgLD, mean(as.numeric(GenoLD)))
    }
    
    cat(paste("Finished: LG ", lg, "\n", 
              "Avg Corrected Per Ind.: ", ceiling(nocor/(nrow(grp) - 2)), "\n",
              "Avg LDLength Per Ind.: ", ceiling(mean(as.numeric(AvgLD))), "\n", sep = ""))
    
    if (lg %in% pullouts$LG) {
      marks <- pullouts$MARKNAME[which(pullouts$LG == lg)]
      marks <- droplevels(df[, which(colnames(df) %in% marks)])
      corrected <- cbind(corrected, grp, marks)
    } else {
      corrected <- cbind(corrected, grp)
    }
  }
  
  return(corrected) 
}

###### READING IN THE DATA ######

FEM <- read.cross("csv", "", ".csv", estimate.map = F, na.strings = c("NA", "-"),
			genotypes = c("AA", "AB"), alleles = c("A", "B")) 

###### DIAGNOSTICS ######

## Identify Duplicate Individuals ##

cgd <- genClones(FEM, tol = 0.95, id = "GENOTYPE")$cgd; cgd

## Removing Problematic Markers ## 

FEM <- pullCross(FEM, type = "co.located"); ncol(FEM$co.located$data)
FEM <- pullCross(FEM, type = "seg.distortion", 
                 pars = list(seg.thresh = (0.05/totmar(FEM)))); ncol(FEM$seg.distortion$data)

rm(cgd)

###### CHECKPOINT I ######

## OUT ##
write.cross(FEM, filestem = "_CK1", format = "csv")

## IN ## 
FEM <- read.cross("csv", "", "_CK1.csv", estimate.map = F, na.strings = c("NA", "-"),
			genotypes = c("AA", "AB"), alleles = c("A", "B")) 

###### INITIAL FORMATION OF LINKAGE GROUPS #######

# The p.value needs to be adjusted until the expected number of linkage groups is acheived 
FEM <- mstmap(FEM, p.value = 1e-12, objective.fun = "ML", dist.fun = "kosambi", id = "GENOTYPE", 
			suffix = "numeric", bychr = F, mvest.bc = F, detectBadData = F, trace = F, return.imputed = F)
length(FEM$geno)

# Removal of Singletons and Tiny LGs #
FEM <- subset(FEM, chr = names(FEM$geno)[nmar(FEM) > 10]) 
names(FEM$geno) <- c(1:length(FEM$geno))
length(FEM$geno)

# Splits marker ID to show a marker's physical alignment compared to how markers are linked 
chrIDbymarks(FEM, names(FEM$geno))

# check for switched alleles and remove high error marks #
catch <- SwitchFinder(FEM, 5)
drop <- MessyMarks(catch$M, catch$N); length(drop)
FEM <- drop.markers(FEM, unique(drop)); totmar(FEM)

# reform linkage groups  #
FEM <- mstmap(FEM, p.value = 1e-12, objective.fun = "ML", dist.fun = "kosambi", id = "GENOTYPE", 
			suffix = "numeric", bychr = F, mvest.bc = F, detectBadData = F, trace = F, return.imputed = F)
length(FEM$geno)
FEM <- subset(FEM, chr = names(FEM$geno)[nmar(FEM) >= 10])
names(FEM$geno) <- c(1:length(FEM$geno)); length(FEM$geno)

chrIDbymarks(FEM, names(FEM$geno)); totmar(FEM)

###### CHECKPOINT II ######

## OUT ##
write.cross(FEM, filestem = "_CK2", format = "csv")

## IN ##
FEM <- read.cross("csv", "", "_CK2.csv", estimate.map = F, na.strings = c("NA", "-"),
			genotypes = c("AA", "AB"), alleles = c("A", "B")) 

###### SWITCHING ALLELES FOR 19 LINKAGE GROUPS ######

# Resources frm R/qtl have great explanations of switched alleles #

# Identifies if separate linkage groups aligned to the same chromosome # 
lgnames <- c()
for (lg in seq(1:length(FEM$geno))) {
  c <- data.frame(table(as.numeric(str_replace(word(markernames(FEM, lg), 1, sep = "_"), "S", ""))))
  lgnames <- append(lgnames, as.character(c$Var1[which.max(c$Freq)]))
}
unique(lgnames[duplicated(lgnames)][duplicated(lgnames[duplicated(lgnames)])])

# Uses the assumption that the smaller of the linkage 
catch <- SwitchFinder(FEM, 5)
drop <- MessyMarks(catch$M, catch$N); length(drop)
FEM <- drop.markers(FEM, unique(drop)); totmar(FEM)

mins <- unique(FindSmallerLGs(FEM, catch$M))
toswitch <- markernames(FEM, chr = mins)
FEM <- switchAlleles(FEM, toswitch)

# As before, change the p.value until you have the expected number of linkage groups
FEM <- mstmap(FEM, p.value = 1e-12, objective.fun = "ML", dist.fun = "kosambi", id = "GENOTYPE", 
			suffix = "numeric", bychr = F, mvest.bc = F, detectBadData = F, trace = F, return.imputed = F)
length(FEM$geno)
FEM <- subset(FEM, chr = names(FEM$geno)[nmar(FEM) >= 10])
names(FEM$geno) <- c(1:length(FEM$geno)); length(FEM$geno)

chrIDbymarks(FEM, names(FEM$geno)); totmar(FEM)

###### CHECKPOINT III ###### 

# Just cleaning things up # 
rm(c, catch, drop, mins, toswitch, lg, lgnames)

## OUT ##
write.cross(FEM, filestem = "_CK3", format = "csv")

## IN ##
FEM <- read.cross("csv", "", "_CK3.csv", estimate.map = F, na.strings = c("NA", "-"),
                  genotypes = c("AA", "AB"), alleles = c("A", "B")) 

###### LOCKING DOWN LGs AND CLEANING UP SINGLETONS ###### 

FEM <- pullCross(FEM, type = "co.located"); ncol(FEM$co.located$data)
FEM <- drop.markers(FEM, MarkerOCD(FEM))
totmar(FEM)

lgnames <- c()
for (lg in seq(1:length(FEM$geno))) {
  c <- data.frame(table(as.numeric(str_replace(word(markernames(FEM, lg), 1, sep = "_"), "S", ""))))
  lgnames <- append(lgnames, as.character(c$Var1[which.max(c$Freq)]))
}

lgnames[duplicated(lgnames)]
which(lgnames == "")

# p.value greater than 1 will not break formed linkage groups # 
FEM <- mstmap(FEM, p.value = 10, objective.fun = "ML", dist.fun = "kosambi", id = "GENOTYPE", 
              suffix = "numeric", bychr = T, mvest.bc = F, detectBadData = F, trace = F, return.imputed = F)

chrIDbymarks(FEM, names(FEM$geno), F); totmar(FEM)

###### CHECKPOINT IV ###### 

rm(lgnames, c, lg)

## OUT ##
write.cross(FEM, filestem = "_CK4", format = "csv")

## IN ##
FEM <- read.cross("csv", "", "_CK4.csv", estimate.map = F, na.strings = c("NA", "-"),
                  genotypes = c("AA", "AB"), alleles = c("A", "B")) 

###### REDUCING DOUBLE CROSSOVERS & DROPPING ROGUE XO INDIVIDUALS ######

# Performs similar action as ABH genotypes, with simple error corrections based 
# on formed linkage groups and sequence alignments 
FEMX <- read.csv("_CK4.csv")
FEMX <- XODXOCorrector(FEMX)

write.table(FEMX, "_DXO.csv", quote = F, sep = ",", row.names = F, col.names = T)

FEMX <- read.cross("csv", "", "_DXO.csv", estimate.map = F, na.strings = c("NA", "-"),
                  genotypes = c("AA", "AB"), alleles = c("A", "B")) 

FEMX <- mstmap(FEMX, p.value = 10, objective.fun = "ML", dist.fun = "kosambi", id = "GENOTYPE", 
              suffix = "numeric", bychr = T, mvest.bc = F, detectBadData = F, trace = F, return.imputed = F)

profileMark(FEMX, stat.type = c("seg.dist", "dxo", "recomb"), layout = c(1, 3))

FEMX <- pullCross(FEMX, type = "co.located"); ncol(FEMX$co.located$data)
totmar(FEMX)

FEMX <- mstmap(FEMX, p.value = 10, objective.fun = "ML", dist.fun = "kosambi", id = "GENOTYPE", 
               suffix = "numeric", bychr = T, mvest.bc = F, detectBadData = F, trace = F, return.imputed = F)

profileGen(FEMX, stat.type = c("xo", "dxo"), bychr = F, id = "GENOTYPE", layout = c(1, 2))

# checks for genotypes that retained a large amount of double crossovers after error correction
# and removes them
Gen <- profileGen(FEMX, bychr = F, stat.type = c("xo", "dxo"), id = "GENOTYPE", 
                    layout = c(1, 2), lty = 2, cex = 0.7, xo.lambda = 100)
Gen <- Gen$stat$xo[names(Gen$stat$xo) %in% names(which(Gen$xo.lambda))]; names(Gen)

FEMX <- subsetCross(FEMX, ind = (!(FEMX$pheno$GENOTYPE %in% names(Gen)))); nind(FEMX)

FEMX <- mstmap(FEMX, p.value = 10, objective.fun = "ML", dist.fun = "kosambi", id = "GENOTYPE", 
               suffix = "numeric", bychr = T, mvest.bc = F, detectBadData = F, trace = F, return.imputed = F)

FEMX <- pullCross(FEMX, type = "co.located"); ncol(FEMX$co.located$data)
totmar(FEMX)

FEMX <- mstmap(FEMX, p.value = 10, objective.fun = "ML", dist.fun = "kosambi", id = "GENOTYPE", 
               suffix = "numeric", bychr = T, mvest.bc = F, detectBadData = F, trace = F, return.imputed = F)

###### CHECKPOINT V ###### 

## OUT ## 
write.cross(FEMX, filestem = "_CK5", format = "csv")

rm(FEMX, Gen)

## IN ##

FEM <- read.cross("csv", "", "_CK5.csv", estimate.map = F, na.strings = c("NA", "-"),
                  genotypes = c("AA", "AB"), alleles = c("A", "B")) 

###### MAP STATS ######

### CLIPPING TAILS ###

# Look at the linkage map to ID and remove long tails from linkage group
plotMap(FEM)

# Example # 
FEM$geno[[4]]$map 
FEM <- breakCross(FEM, split = list("4" = "S04_5434796"))

plotMap(FEM)

### ORIENTING LINKAGE GROUPS ###

# Look at the order of markers by cM and by physical distance 
# inverting the order if necessary
cm <- getcM(FEM)
par(mfrow = c(4, 5), mar = c(4, 4, 0.5, 0.5))

for (i in seq(1:19)) { # number of LGs
  plot(cm[[i]]$cMdist, cm[[i]]$physdist,
       pch = 19, xlab = "centiMorgan (cM)", 
       ylab = "Physical Distance (Mb)")
  legend("topleft", paste("Chr", as.character(i)), 
         bty = "n", cex = 1.5, inset = -0.05)
}

# Example # 
FEM <- flip.order(FEM, c(1))

cm <- getcM(FEM)
par(mfrow = c(4, 5), mar = c(4, 4, 0.5, 0.5))
for (i in seq(1:19)) {
  plot(cm[[i]]$cMdist, cm[[i]]$physdist,
       pch = 19, xlab = "centiMorgan (cM)", 
       ylab = "Physical Distance (Mb)")
  legend("topleft", paste("Chr", as.character(i)), 
         bty = "n", cex = 1.5, inset = -0.05)
}

### MAP STATS ### 

summaryMap(FEM)

###### CHECKPOINT FINAL ###### 


## OUT ## 
write.cross(FEM, filestem = "_FINAL", format = "csv")

## IN ##
FEM <- read.cross("csv", "", "_FINAL.csv", estimate.map = F, na.strings = c("NA", "-"),
                  genotypes = c("AA", "AB"), alleles = c("A", "B")) 
