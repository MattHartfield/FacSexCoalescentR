## AsexCoalescentSim.R
## Coalescent simulation with infrequent rates of sex
## Simulation originally written for use in Hartfield, Wright and Agrawal 2015 paper
## "Coalescent times and patterns of genetic diversity in species with facultative sex: effects of gene conversion, population structure and heterogeneity".
## Comments should be sent to Matthew Hartfield (matthew.hartfield@utoronto.ca)

## See related README.TXT file for instructions on execution

rm(list = ls()) ## Clearing space
require(compiler)	## For speeding up execution
invisible(enableJIT(3))

# Reading in arguments by command line
args <- commandArgs(trailingOnly = TRUE)

## Variables:
## THE FIRST EIGHT are: Nt; gen. conv.; theta; initial state; low->high sex prob; high->low sex prob; net migration rate; # demes.
Na <- as.integer(args[1])	# Population size
g <- as.double (args[2])	# Rate of gene conversion
theta <- as.double(args[3])	# Rate of mutation (4 NT mu)
pSTIN <- as.double(args[4])	# What to do (0=changing sex, 1=stepwise change, 2 = constant sex)
pLH <- as.double(args[5])	# Prob of low-sex to high-sex transition OR time of transition if stepwise change
pHL <- as.double(args[6])	# Prob of high-sex to low-sex transition
mig <- as.double(args[7])	# Migration rate (2 NT m)
mig <- mig/(2*Na)			# Rescaling migration for probability transitions
d <- as.integer(args[8])	# Number of demes in island model
if(d == 1){
	mig <- 0				# Set migration to zero if only one deme, as a precaution
}

# Initial Error checking
if(Na <= 0){
	stop("Total Population size N is zero or negative, not allowed.")
}
if(g < 0 || g > 1){
	stop("Rate of gene conversion has to lie between 0 and 1.")
}
if(g > (10/Na)){
	cat("WARNING: Analytical transitions assume gene conversion is weak.\n")
	cat("Current input is much larger than O(1/NT) - transitions may be inaccurate.\n")
}
if(theta < 0){
	stop("Mutation rate must be a positive (or zero) value.")
}
if(mig > (10/Na)){
	cat("WARNING: Analytical transitions assume migration is weak.\n")
	cat("Current input is much larger than O(1/NT) - transitions may be inaccurate.\n")
}
if(d <= 0){
	stop("Number of demes has to be a positive integer")
}
if(pSTIN != 1){
	if( any( c(pHL,pLH) < 0) || any( c(pHL,pLH) > 1)  ){
	stop("Sex transition probabilities have to lie between 0 and 1 (if there is no stepwise change in sex).")
}
}
if( any( c(pHL,pLH) == 0) && pSTIN == 0){
	stop("Sex transition probabilities have to lie between 0 and 1.")
}
if( !(pSTIN %in% c(0,1,2)) ){
	stop("pSTIN has to equal 0, 1, or 2.")
}

## Initial state of two samples;
## Iwith = No. of within host
## Ibet = No. of between host
## So total number of initial samples in a deme = 2*Iwith + Ibet
## THEN sexL[i], sexH[i] = low-sex, high-sex rate in each deme
sps <- args[9:(8+4*d)]
Iwith <- as.integer(sps[seq(1,4*d-3,by=4)])
Ibet <- as.integer(sps[seq(2,4*d-2,by=4)])
Itot <- (2*sum(Iwith) + sum(Ibet))
Iindv <- sum(Iwith) + sum(Ibet)
sexL <- as.double(sps[seq(3,4*d-1,by=4)])
sexH <- as.double(sps[seq(4,4*d,by=4)])

## Number of samples/reps to take
Nreps <- as.integer(args[4*d + 9])

Na <- Na/d	# Scaling NT to a demetic size, for more consistent use in calculations

# More error checking
if(any(sexL < 0) || any(sexL > 1) || any(sexH < 0) || any(sexH > 1)){
	stop("Rate of sexual reproduction has to lie between 0 and 1.")
}
if(any(sexH == 0) && pSTIN != 2){
	stop("Rate of sexual reproduction has to lie between 0 and 1.")
}
if(any(sexL > sexH) && pSTIN == 0){
	stop("All low-sex values must be less than or equal to high-sex values. Please re-check.")
}
if(any(Iwith < 0) || any(Ibet < 0)){
	stop("Number of within- and between-host samples has to be positive.")
}
if(Itot <= 1){
	stop("More than one sample must be initially present to execute the simulation.")
}
if(Nreps <= 0){
	stop("Must set positive number of repetitions.")
}

## Matrix of coalescent times for each rep
TCoal <- matrix(0,nrow=Nreps,ncol=(Itot - 1))
Twees <- vector(mode="character",length=Nreps) # Newick tree data per run
AllM <- vector("list", Nreps)	# List of mutation matrices

# Pre-defining triangle function for easing code alter on
# Also functions for each event
Trig <- function(x){ x*(x+1)/2 }
# Defining functions of each transition, for rapid per-deme calculations
P23 <- function(y,k,Na){ (k*y)/(Na) + choose(2*k,2)/(2*Na) }	# One of the paired samples is recreated: (x,y) -> (x-k+1, y + 2(k-1)). OR A unique samples coaleses with another (either pre-existing or new): (x,y) -> (x-k, y + 2k - 1)
P4 <- function(x,k,Na){((x-k)*2*k)/(Na)}	# A paired sample coaleses with new unique sample: (x,y) -> (x-k,y + 2k -1)
P56 <- function(y,Na){choose(y,2)/(2*Na)}	# Two pre-existing unique samples re-create a paired sample: (x,y) -> (x - k + 1, y + 2(k-1)). OR Two pre-existing paired samples coalesce: (x,y) -> (x - k, y + 2k-1)
P7 <- function(x,k,Na){(choose(x-k,2)/Na)}	# Two remaining paired samples doubly coalesce asexually: (x,y) -> (x-k-1,y+2k)
P8 <- function(x,y,k,Na){(x-k)*y/Na}		# One of the x - k remaining paired samples can coalesce with a unique sample: (x,y) -> (x-k,y+2k-1)
P9 <- function(x,k,g){g*(x-k)}				# Paired sample coaleses via gene conversion: (x,y) -> (x-k-1,y+2k+1)
P10 <- function(x,y,k,mig){mig*(x + y + k)}		# A sample migrates to another deme

# Updated function to calculate probability change vectors each time OVER EACH DEME
probset2 <- function(Na,g,mig,Nwith,Nbet,k,sw){
	pr <- matrix(data=0,nrow=9,ncol=d)		# Matrix of probabilities, if k splits per deme already determined
		
	pr[4,] <- mapply(P56,Nbet,Na)
	pr[5,] <- mapply(P56,Nbet,Na)
	pr[6,] <- mapply(P7,Nwith,k,Na)
	pr[7,] <- mapply(P8,Nwith,Nbet,k,Na)
	pr[8,] <- mapply(P9,Nwith,k,g)
	pr[9,] <- mapply(P10,Nwith,Nbet,k,mig)
	
	if(sw == 1){	#	Only activate the first three events if need to consider sex (fourth is 'split pairs remain split')
		if(sum(k) != 1){
			pr[1,] <- mapply(P23,Nbet,k,Na)
		}
		pr[2,] <- mapply(P23,Nbet,k,Na)
		pr[3,] <- mapply(P4,Nwith,k,Na)
	}

	return(pr)
}	# End of 'probset2' function

# Updated function to determine how to change state numbers following 'an event', taking into account events over all demes
stchange2 <- function(ev,deme,k,Nwith){
	
	oo1 <- -k
	oo2 <- 2*k
	
	## Now deciding extra events depending on deme and event
	oo3 <- switch(ev,
	c(0,0),
	c(1,-2),
	c(0,-1),
	c(0,-1),
	c(1,-2),
	c(0,-1),
	c(-1,0),
	c(0,-1),
	c(-1,1),
	c(0,0)
	)
	oo1[deme] <- oo1[deme] + oo3[1]
	oo2[deme] <- oo2[deme] + oo3[2]
	
	outlist <- list("WCH" = oo1, "BCH" = oo2)
	return(outlist)
}	# End of 'stchange' function

WHtab <- function(itab){	# Returning subtable of WH samples
	WH <- itab[itab[,3]==0,]
	if(is.null(dim(WH)) == 1){
		WH <- t(as.matrix(WH))
	}
	return(WH)
}	# End of 'WHtab' function

BHtab <- function(itab){	# Returning subtable of BH samples
	BH <- itab[itab[,3]==1,]
	if(is.null(dim(BH)) == 1){
		BH <- t(as.matrix(BH))
	}
	return(BH)
}	# End of 'BHtab' function

Ctab <- function(itab){	# Returning subtable of coalesced samples
	CT <- itab[itab[,3]==2,]
	if(is.null(dim(CT)) == 1){
		CT <- t(as.matrix(CT))
	}
	return(CT)
}	# End of 'Ctab' function

rntabs <- function(WH,BH,CT,Nwith,Nbet){	# Renumbering and merging tables after event

	WH <- rbind(WH,BH[BH[,3]==0,])	# Adding merged samples to WH table
	BH <- rbind(BH,WH[WH[,3]==1,])	# Adding split samples to BH table
	CT <- rbind(CT,rbind(WH[WH[,3]==2,],BH[BH[,3]==2,])) # Adding new coalesced samples to CT table
	
	WH <- WH[WH[,3]==0,]	# Removing split entries
	BH <- BH[BH[,3]==1,]	# Removing merged entries
	if(is.null(dim(WH)) == 1){
		WH <- t(as.matrix(WH))
	}
	if(is.null(dim(BH)) == 1){
		BH <- t(as.matrix(BH))
	}
	
	# Reordering then renumbering new subtables
	WH <- WH[order(WH[,6],WH[,2],WH[,1]),]
	WH[,2] <- rep(c(1:sum(Nwith)),each=2)

	BH <- BH[order(BH[,6],BH[,2]),]
	if(is.null(dim(BH)) != 1 && dim(BH)[1]!=0){
		BH[,2] <- c((sum(Nwith) + 1):(sum(Nwith) + sum(Nbet)))
	}else if(is.null(dim(BH)) == 1){
		BH[2] <- (sum(Nwith) + sum(Nbet))
	}
	
	itab <- rbind(WH,BH,CT)
	rownames(itab) <- NULL
	return(itab)
}	# End of 'rntabs' function

ssamp <- function(WH,BH,k,deme,state){	# Picking individuals for action based on action
	switch(state,
		{
			sps <- unique(WH[WH[,6]==deme,2])
			if(length(sps) > 1 | k == 0){
				sample(sps,k)
			}else if(length(sps) == 1 & k != 0){
				sps
			}
		},
		{
			sps <- unique(BH[BH[,6]==deme,1])
			if(length(sps) > 1){
				sample(sps,1)
			}else if(length(sps) == 1){
				sps
			}
		},
		sample(unique(BH[BH[,6]==deme,1]),2),
	)
}		# End of 'ssamp' function

# Function to determine which paired samples are split by sex
sexsamp <- function(WH,k){
	out <- numeric(0)
	for(j in 1:length(k)){
		sps <- unique(WH[WH[,6]==j,2])
		if(length(sps) > 1 | k[j] == 0){
			s2 <- sample(sps,k[j])
		}else if(length(sps) == 1 & k[j] != 0){
			s2 <- sps
		}
		out <- append(out,s2)
	}
	return(out)
}

cchange  <- function(tab,csamp,par,Tt){		# Renumbering tables after coalescent event
	
	for(i in 1:length(csamp)){
		tab[tab[,1]==csamp[i],2] <- 0			# 'Csamp' coalesces!
		tab[tab[,1]==csamp[i],3] <- 2			# 'Csamp' coalesces!
		tab[tab[,1]==csamp[i],4] <- Tt 			# Coalescent time
		tab[tab[,1]==csamp[i],5]  <- par[i]		# Noting parental sample
	}
	
	return(tab)
}	# End of 'cchange' function

# Function to change status of samples following event change
coalesce <- function(itab,Tt,Nwith,Nbet,din,kin,exin,dein,e2in){
	
	deme <- din
	k <- kin
	ex <- exin
	
	WH <- WHtab(itab)
	BH <- BHtab(itab)
	CT <- Ctab(itab)
	
	switch(ex,
	{	# Event 1: 2k new samples created from paired samples
		rsex <- sexsamp(WH,k)
		WH[WH[,2]%in%rsex,3] <- 1	# Setting samples as 'between-host' (due to split)
	},	
	{	# Event 2: One of the paired samples is recreated, no coalescence
		rsex <- sexsamp(WH,k)				# Paired samples that split...
		#rands2 <- ssamp(WH,BH,1,deme,1)	# ...except one that does not do so fully.
		rands2 <- sample(WH[(WH[,2]%in%rsex & WH[,6]==deme),2],1)			# ...except one that does not do so fully.
		nos <- (rbinom(1,1,0.5) + 1)		# Sub-sample that remained paired
		
		WH[WH[,2]%in%rsex,3] <- 1	# Setting samples as 'between-host' (due to split)
		WH[WH[,2]%in%rands2,3][nos] <- 0	# Except one that remained paired

		if(dim(rbind(WH[WH[,3]==1 & WH[,6]==deme,],BH[BH[,3]==1 & BH[,6]==deme,]))[1] != 1){
			bhc <- sample(rbind(WH[WH[,3]==1 & WH[,6]==deme,],BH[BH[,3]==1 & BH[,6]==deme,])[,1],1)	# Choosing unique sample that re-pairs
		}else if(dim(rbind(WH[WH[,3]==1 & WH[,6]==deme,],BH[BH[,3]==1 & BH[,6]==deme,]))[1] == 1){
			bhc <- rbind(WH[WH[,3]==1 & WH[,6]==deme,],BH[BH[,3]==1 & BH[,6]==deme,])[,1]
		}
		if(bhc %in% WH[,1]){
			WH[WH[,1]==bhc,3] <- 0	# Becomes within host
			WH[WH[,1]==bhc,2] <- rands2	# Ensuring same parents are made
		}else if(bhc %in% BH[,1]){
			BH[BH[,1]==bhc,3] <- 0	# Becomes within host
			BH[BH[,1]==bhc,2] <- rands2	# Ensuring same parents are made
		}
	},	
	{	# Event 3: One of the unique samples coaleses with another unique one (either pre-existing or new)
		rsex <- sexsamp(WH,k)
		WH[WH[,2]%in%rsex,3] <- 1	# Setting samples as 'between-host' (due to split)
						
		ord <- sample(WH[(WH[,2]%in%rsex & WH[,6]%in%deme),1],1)	# One sample involved in coalescence
		if(length(rbind(WH[WH[,3]==1 & WH[,6]==deme,],BH[BH[,3]==1 & BH[,6]==deme,])[,1]) != 2){
			ord <- c(ord,sample(setdiff(rbind(WH[WH[,3]==1 & WH[,6]==deme,],BH[BH[,3]==1 & BH[,6]==deme,])[,1],ord),1))		# Second sample involved		
		# }else if(length(rbind(WH[WH[,3]==1,] & WH[,6]==deme,BH[BH[,3]==1 & BH[,6]==deme,])[,1]) == 2){					
		}else if(length(rbind(WH[WH[,3]==1 & WH[,6]==deme,],BH[BH[,3]==1 & BH[,6]==deme,])[,1]) == 2){						# If only two samples exist, only they can pair, so deliberately choosing other sample
			ord <- c(ord,setdiff(rbind(WH[WH[,3]==1 & WH[,6]==deme,],BH[BH[,3]==1 & BH[,6]==deme,])[,1],ord))
		}
		jumb <- sample(ord)
		csamp <- jumb[1]	# The one that disappears
		par <- jumb[2]		# The one that does not
		
		itab2 <- cchange(rbind(WH,BH),csamp,par,Tt)	# Updating coalescent times
		WH <- WHtab(itab2)
		BH <- BHtab(itab2)
		CT <- rbind(CT,Ctab(itab2))
		# ol <- list("tw"=WH,"tb"=BH,"tc"=CT)
		# print(ol)
	},	
	{	# Event 4: A paired sample coaleses with new unique sample
		rsex <- sexsamp(WH,k)
		WH[WH[,2]%in%rsex,3] <- 1	# Setting samples as 'between-host' (due to split)
		
		ord <- sample(WH[WH[,3]==1 & WH[,6]==deme,1],1)	# Unique sample involved in coalescence
		ord <- c(ord,sample(WH[WH[,3]==0 & WH[,6]==deme,1],1))	# Paired sample involved in coalescence
		jumb <- sample(ord)
		csamp <- jumb[1]	# The one that disappears
		par <- jumb[2]	# The one that does not
		
		if(WH[WH[,1]==par,3] == 1){		# Correction if parential sample is unique sample
			WH[WH[,1]==par,3] <- 0			# Merged sample becomes WH
			WH[WH[,1]==par,2] <- WH[WH[,1]==csamp,2]			# Ensuring remaining sample enters same parent as coalesced sample
		}
		WH <- cchange(WH,csamp,par,Tt)	# Updating coalescent times
	},
	{	# Event 5: Two pre-existing unique samples re-create paired sample.
		rsex <- sexsamp(WH,k)
		WH[WH[,2]%in%rsex,3] <- 1	# Setting samples as 'between-host' (due to split)

		rands <- ssamp(WH,BH,k,deme,3)	# Two combining samples
		BH[BH[,1]%in%rands,3] <- 0	# Setting samples as 'within-host'
		BH[BH[,1]%in%rands,2][2] <- BH[BH[,1]%in%rands,2][1]	# Making sure they have same parents		
	},
	{	# Event 6: Two pre-existing unique samples coalesce
		rsex <- sexsamp(WH,k)
		WH[WH[,2]%in%rsex,3] <- 1	# Setting samples as 'between-host' (due to split)
		
		rands <- ssamp(WH,BH,k,deme,3)	# Two coalesced samples
		csamp <- sample(rands,1)		# The one that disappears
		par <- rands[rands!=csamp]		# The one that does not
		
		BH <- cchange(BH,csamp,par,Tt)	# Updating coalescent times
	},
	{	# Event 7: Two remaining paired samples doubly coalesce asexually.
		rsex <- sexsamp(WH,k)
		WH[WH[,2]%in%rsex,3] <- 1	# Setting samples as 'between-host' (due to split)
		
		rands <- ssamp(WH[!(WH[,2]%in%rsex),],BH,2,deme,1)			# 2 other within-host samples that coalesce
		lhs <- sample(WH[WH[,2]%in%rands,1][c(1,3)])		## Randomly resample the four samples; let them be parent/coalescent pairs respectively. These two are sample/parent for left 'Ceplitis' deme
		rhs <- sample(WH[WH[,2]%in%rands,1][c(2,4)])		## These two are sample/parent for right 'Ceplitis' deme
		csamp <- c(lhs[1],rhs[1])
		par <- c(lhs[2],rhs[2])
		
		WH <- cchange(WH,csamp,par,Tt)
		WH[WH[,1]==par[2],2] <- WH[WH[,1]==par[1],2]		# Making sure parental samples are in same individual
	},
	{	# Event 8: One of the x - k remaining paired samples coalesces with a unique sample
		rsex <- sexsamp(WH,k)
		WH[WH[,2]%in%rsex,3] <- 1	# Setting samples as 'between-host' (due to split)
		# if(k != 0){
			# WH[WH[,2]%in%rands[1:k],3] <- 1	# Setting samples as 'between-host' (due to split)
		# }
		
		rands <- ssamp(WH[!(WH[,2]%in%rsex),],BH,1,deme,1)	# A paired sample involved in coalescence
		ord <- sample(WH[WH[,2]==rands,1],1)	# A paired sample involved in coalescence
		ord <- c(ord,ssamp(WH,BH,0,deme,2))	# Existing Unique sample involved in coalescence
		jumb <- sample(ord)
		csamp <- jumb[1]	# The one that disappears
		par <- jumb[2]	# The one that does not
		
		itab2 <- cchange(rbind(WH,BH),csamp,par,Tt)	# Updating coalescent times
		WH <- WHtab(itab2)
		BH <- BHtab(itab2)
		CT <- rbind(CT,Ctab(itab2))
		
		if(par%in%BH[,1] == 1){		# Correction if parential sample is unique sample
			BH[BH[,1]==par,3] <- 0 	# Becomes within-host
			BH[BH[,1]==par,2] <- rands	# same parent as WH sample
		}
	},
	{	# Event 9: Paired sample coaleses via gene conversion: (x,y) -> (x-k-1,y+2k+1)
		rsex <- sexsamp(WH,k)
		WH[WH[,2]%in%rsex,3] <- 1	# Setting samples as 'between-host' (due to split)
		
	 	# if(k != 0){
			# WH[WH[,2]%in%rands[1:k],3] <- 1	# Setting samples as 'between-host' (due to split)
		# }
		
		rands <- ssamp(WH[!(WH[,2]%in%rsex),],BH,1,deme,1)	# A paired sample involved in coalescence
		ord <- sample(WH[WH[,2]==rands,1])	# Jumbing samples for coalescence
		csamp <- ord[1]				# The one that disappears
		par <- ord[2]					# The one that do not

		WH <- cchange(WH,csamp,par,Tt)	# Updating coalescent times
		WH[WH[,1]%in%par,3]	<- 1 # Making sure non-coalesced sample becomes BH
	},
	{	# Event 10: Migration of a sample
		rsex <- sexsamp(WH,k)
		WH[WH[,2]%in%rsex,3] <- 1	# Setting samples as 'between-host' (due to split)
		
		if(e2in == 1){
			bhc <- ssamp(WH[!(WH[,2]%in%rsex),],BH,1,deme,1)	# Choosing the WH sample that migrates
			WH[WH[,2]%in%bhc,6] <- dein							# Altering deme accordingly
		}else if(e2in == 2){
			if(length(rbind(WH[WH[,3]==1 & WH[,6]==deme,],BH[BH[,3]==1 & BH[,6]==deme,])[,1]) > 1){
				bhc <- sample(rbind(WH[WH[,3]==1 & WH[,6]==deme,],BH[BH[,3]==1 & BH[,6]==deme,])[,1],1)
			}else if(length(rbind(WH[WH[,3]==1 & WH[,6]==deme,],BH[BH[,3]==1 & BH[,6]==deme,])[,1]) == 1){
				bhc <- rbind(WH[WH[,3]==1 & WH[,6]==deme,],BH[BH[,3]==1 & BH[,6]==deme,])[,1]
			}
			if(bhc %in% WH[,1]){
				WH[WH[,1]%in%bhc,6] <- dein		# Altering deme accordingly
			}else if(bhc %in% BH[,1]){
				BH[BH[,1]%in%bhc,6] <- dein		# Altering deme accordingly
			}
		}
	},
	)
	
	itab <- rntabs(WH,BH,CT,Nwith,Nbet)	# Renumbering individuals
	return(itab)
	
}	# End of 'coalesce' function

# Function to reconstruct geneaology (based on and extending code from Jeremy Gray)
# Also to add mutation to samples (separate subfunction)
treemaker <- function(itab,theta){
	CT <- Ctab(itab)
	lct <- length(CT[,1]) + 1
	
	clades <- vector(mode="character",length=lct)	# Character vector of possible clades
	Cheight <- vector(mode="numeric",length=lct)	# Current 'height' (time) of each clade
	samps <- matrix(0,lct,lct)	# Table of samples present in each clade (row = each clade)
	nc <- 1	# {N}umber of {c}lades in current reconstruction
	MTab <- matrix(nrow=0,ncol=(lct+1))
	colnames(MTab) <- c("Position",c(1:(lct)))
	nmut <- 0
		
	for(i in 1:length(CT[,1])){
		birthtime <- CT[i,4]
	    child1 <- CT[i,1]
    	parent1 <- CT[i,5]
    	csum <- sum((child1 %in% samps) , (parent1 %in% samps))	# Testing how many of the pair have already been sampled, to decide tree reconstruction
    	ischild <- 0
    	
    	if(i == 1){
    		clades[1] <- paste("(",parent1,":",birthtime,",",child1,":",birthtime,")",sep="")
    		samps[1,1:2] <- c(parent1,child1)
    		Cheight[1] <- birthtime

			rmut <- rpois(2,lambda=(0.5*theta*birthtime))	# New mutations present in first and second sample respectively
			twos <- c(parent1,child1)
			for(a in 1:2){
				if(rmut[a] != 0){
					for(b in 1:rmut[a]){	# Adding mutations to table
			   			MTab <- rbind(MTab,rep(0,lct + 1))
    					MTab[nmut+b,1] <- runif(1)
	   					MTab[nmut+b,twos[a] + 1] <- 1		# Indicating location of mutants
 					}
 					nmut <- nmut + rmut[a]
				}
			}

    	}else if(i > 1){
    		# There can be three cases: Merge clades if child already listed; Add to clade if child new but parent already listed; Create new clade otherwise.
    		if(nc > 1 && csum==2){	# Merging existing clades if both samples already
    			cc <- which(samps==child1,arr.ind=T)[1]		# Choosing row (and therefore clade) containing existing child
    			pc <- which(samps==parent1,arr.ind=T)[1]		# Choosing row (and therefore clade) containing existing parent
    			ctext <- clades[cc]
    			ptext <- clades[pc]
    			tc <- paste("(",ptext,":",birthtime-Cheight[pc],",",ctext,":",birthtime-Cheight[cc],")",sep="")
    			
    			# Adding mutations
    			rmut <- rpois(2,lambda=(0.5*theta*c(birthtime-Cheight[pc],birthtime-Cheight[cc])))	# New mutations present in parential and child clade respectively
    			parsamps <- samps[pc,samps[pc,]!=0]
    			chisamps <- samps[cc,samps[cc,]!=0]
				for(a in 1:2){
					if(rmut[a] != 0){
						for(b in 1:rmut[a]){	# Adding mutations to table
			   				MTab <- rbind(MTab,rep(0,lct + 1))
    						MTab[nmut+b,1] <- runif(1)
    						if(a == 1){
    							MTab[nmut + b,parsamps + 1] <- 1
    						}else if(a == 2){
    							MTab[nmut + b,chisamps + 1] <- 1
    						}
	 					}
 						nmut <- nmut + rmut[a]
					}
				}
    			
    			# Now delete old clade data
    			j <- min(cc,pc)
    			k <- max(cc,pc)
    			clades[j] <- tc
    			clades[k] <- ""
    			Cheight[j] <- birthtime
    			Cheight[k] <- 0
    			ccs <- samps[cc,samps[cc,]!=0]
    			pcs <- samps[pc,samps[pc,]!=0]
    			ts <- c(ccs,pcs)
    			samps[j,] <- c(ts,rep(0,lct-length(ts)))
    			samps[k,] <- 0
    			
    			# Now to reorder (to prevent overflow/going out of array length)
    			# All I'm doing is extracting the non-zero entries and putting them together, so order should be maintained
    			# First the samples
    			s1 <- samps[apply(samps,1,function(x) !all(x==0)),]	# Non-zero rows
    			s2 <- samps[!(apply(samps,1,function(x) !all(x==0))),]	#zero rows
    			samps <- rbind(s1,s2)
    			rownames(samps) <- NULL
    			
    			# Then the clades and heights
    			tc2 <- clades[clades!=""]
    			ltc2 <- length(tc2)
    			tt <- Cheight[clades!=""]
    			clades <- c(tc2,rep(0,lct-ltc2))
    			Cheight <- c(tt,rep(0,lct-ltc2))
    			
    			nc <- nc - 1
    			
    		}else if(csum==1){	# Add to existing clade	
   				pc <- which(samps==parent1,arr.ind=T)[1]		# Choosing row (and therefore clade) containing existing parent (or child)
   				if(is.na(pc)==1){
   					ischild <- 1
   					pc <- which(samps==child1,arr.ind=T)[1]		
   				}    				
    			tc <- clades[pc]
    			if(ischild == 0){
    				tc <- paste("(",child1,":",birthtime,",",tc,":",birthtime-Cheight[pc],")",sep="")
	   				samps[pc,(lct - sum(samps[pc,]==0) + 1)] <- child1
					csamps <- samps[pc,(samps[pc,]!=0) & (samps[pc,]!=child1)]
	   				exsamps <- child1
    			}else if(ischild==1){
    				tc <- paste("(",parent1,":",birthtime,",",tc,":",birthtime-Cheight[pc],")",sep="")
	   				samps[pc,(lct - sum(samps[pc,]==0) + 1)] <- parent1
	   				csamps <- samps[pc,(samps[pc,]!=0) & (samps[pc,]!=parent1)]
	   				exsamps <- parent1
    			}
   				
   				# Adding mutations
    			rmut <- rpois(2,lambda=(0.5*theta*c(birthtime-Cheight[pc],birthtime)))	# New mutations present in new branch and existing clade respectively
				for(a in 1:2){
					if(rmut[a] != 0){
						for(b in 1:rmut[a]){	# Adding mutations to table
			   				MTab <- rbind(MTab,rep(0,lct + 1))
    						MTab[nmut+b,1] <- runif(1)
    						if(a == 1){
    							MTab[nmut+b,(csamps + 1)] <- 1
    						}else if(a == 2){
    							MTab[nmut+b,(exsamps + 1)] <- 1
    						}
	 					}
 						nmut <- nmut + rmut[a]
					}
				}
				
   				clades[pc] <- tc
   				Cheight[pc] <- birthtime
   				
    		}else if(csum==0){	# Create a new clade
    			nc <- nc + 1
    			clades[nc] <- paste("(",parent1,":",birthtime,",",child1,":",birthtime,")",sep="")
	    		samps[nc,1:2] <- c(parent1,child1)
	    		Cheight[nc] <- birthtime
	    		
				rmut <- rpois(2,lambda=(0.5*theta*birthtime))	# New mutations present in first and second sample respectively
				twos <- c(parent1,child1)
				for(a in 1:2){
					if(rmut[a] != 0){
						for(b in 1:rmut[a]){	# Adding mutations to table
			   				MTab <- rbind(MTab,rep(0,lct + 1))
    						MTab[nmut+b,1] <- runif(1)
	   						MTab[nmut+b,twos[a] + 1] <- 1		# Indicating location of mutants
	 					}
 						nmut <- nmut + rmut[a]
					}
				}
    		}    		
    	}
	}
	
	cout <- paste(clades[1],";",sep="")
	outlist <- list("tree" = cout, "muts" = MTab)
	return(outlist)
  
}	# End of 'treemaker' function

# Function to change rates of sex given a state change
rate_change <- function(pST,pLH,pHL,sexH,sexL,Na,d,switch1){
	# Setting up transition time (tts, or 'time to switch')
	tts <- NA
	if(pST == 0){
		sexCN <- sexL
		while(is.na(tts) == 1){
			tts <- rgeom(1, pLH/(1+pLH))/(2*Na*d)
			if(is.na(tts) == 1){
				cat("WARNING: Jump time was set to NA - transition probabilities may be too small.\n")
			}
		}
		npST <- 1
	}else if(pST == 1){
		sexCN <- sexH
		while(is.na(tts) == 1){
			tts <- rgeom(1, pHL/(1+pHL))/(2*Na*d)
			if(is.na(tts) == 1){
				cat("WARNING: Jump time was set to NA - transition probabilities may be too small.\n")
			}
		}
		npST <- 0
	}else if(pST == 2){	# If stepwise change, alter depending on whether already switched or not
		npST <- pST
		if(switch1 == 0){
			sexCN <- sexL
			tts <- pLH
		}else if(switch1 == 1){
			sexCN <- sexH
			tts <- Inf
		}
	}else if(pST == 3){	# If constant, no time to sex switch
		sexCN <- sexL
		tts <- Inf
		npST <- pST
	}
	outlist <- list("SCNew" = sexCN, "tswitch" = tts,"newpST"=as.integer(npST))
	return(outlist)
}

# Functions to save and restore seed (for debugging)
save_rng <- function(savefile=tempfile()) {
    if (exists(".Random.seed"))  {
        oldseed <- get(".Random.seed", .GlobalEnv)
    } else stop("don't know how to save before set.seed() or r*** call")
    oldRNGkind <- RNGkind()
    save("oldseed","oldRNGkind",file=savefile)
    invisible(savefile)
}

restore_rng <- function(savefile) {
    load(savefile)
    do.call("RNGkind",as.list(oldRNGkind))  ## must be first!
    assign(".Random.seed", oldseed, .GlobalEnv)
}

# The simulation, completed Nreps times
seed <- as.integer(runif(1,min=0,max=1000000))
print(seed)
write.table(as.matrix(seed),file="Seed.out",quote=F,row.names=F,col.names=F)
set.seed(seed)
for(i in 1:Nreps){
	# cat("Run ",i,":\n",sep='')
	if(pSTIN == 0){
		pST <- sample(c(0,1),size=1,prob=c(pHL/(pLH+pHL),pLH/(pLH+pHL)))	# Randomly assigning initial low-sex or high-sex state
	}else if(pSTIN == 1){
		pST <- 2
	}else if(pSTIN == 2){
		pST <- 3
	}
	Ttot <- 0			# Time in past, initiate at zero
	Nwith <- Iwith		# Resetting number of within-host samples
	Nbet <- Ibet		# Resetting number of between-host samples
	Ntot <- Itot		# Total number of samples at time
	Nindv <- Iindv		# Total number of individuals
	
	rcout <- rate_change(pST,pLH,pHL,sexH,sexL,Na,d,0)
	sexC <- as.double(rcout$SCNew)
	tts <- as.double(rcout$tswitch)
	pST <- as.integer(rcout$newpST)
	tls <- 0	# 'Time since Last Switch' or tls
	# print("Last switch, next switch")
	# print(c(tls,tts+tls))

	## Setting up table of individual samples (for drawing genealogy)
	indvs <- matrix(0,nrow=Itot,ncol=6)
	colnames(indvs) <- c("Sample","Individual","State","Coalescent time","Parent Sample","Deme")
	indvs[,1] <- c(1:Ntot)
	if(sum(Nwith) != 0){
		indvs[1:(2*sum(Nwith)),2] <- rep(c(1:sum(Nwith)),each=2)
		indvs[1:(2*sum(Nwith)),3] <- 0
		indvs[1:length(rep(c(1:d),(2*Nwith))),6] <- rep(c(1:d),(2*Nwith))
	}
	if(sum(Nbet) != 0){
		indvs[(2*sum(Nwith)+1):Ntot,2] <- c((sum(Nwith)+1):Nindv)
		indvs[(2*sum(Nwith)+1):Ntot,3] <- 1
		indvs[indvs[,6] == 0,6] <- rep(c(1:d),(Nbet))
	}
	
	while(Ntot > 1){	# Repeat until all samples coalesce to one MRCA
		
		# Setting up vector of state-change probabilities 
		probs <- probset2(Na,g,mig,Nwith,Nbet,rep(0,d),0)
		nosex <- prod((1-sexC)^Nwith)			# Probability of no sex splits, accounting for within-deme variation
		psum <- (1-nosex) + nosex*sum(probs)	# Sum of all event probabilites, for drawing random time
		# print(c(sexC,psum))
		if(psum > 1){
			stop("Summed probabilities exceed one, you need to double-check your algebra.")
		}
		if(psum <= 0 && all(sexC==0)!=1){
			stop("Summed probabilites are zero or negative, you need to double-check your algebra.")
		}
		
		# Drawing time to next event, SCALED TO 2NT GENERATIONS
		if(psum == 1){
			tjump <- 0
		}else if(psum == 0){
			tjump <- Inf
		}else{
			tjump <- rexp(1, rate = psum*(2*Na*d))
		}
		NextT <- (Ttot + tjump)
		
		# Outcomes depends on what's next: an event or change in rates of sex!
		if(NextT > (tls + tts)){ 	# If next event happens after a switch, change rates of sex
			tls <- (tls + tts)	# 'Time since Last Switch' or tls
			Ttot <- tls
			rcout <- rate_change(pST,pLH,pHL,sexH,sexL,Na,d,1)
			sexC <- as.double(rcout$SCNew)
			tts <- as.double(rcout$tswitch)
			pST <- as.integer(rcout$newpST)
			# print("Last switch, next switch")
			# print(c(tls,tts+tls))
		}else if (NextT <= (tls + tts)){	# If next event happens before a switch, draw an action
			
			Ttot <- NextT
						
			# Determines if sex occurs; if so, which samples are chosen
			# (deme-independent Binomial draws)
			esex <- 0
			evsex <- rep(0,d)
			csex <- match(1,rmultinom(1,1,c(nosex*sum(probs),(1-nosex))))-1
			if(csex == 1){				# Working out number of sex events IF it does occur
				while(esex == 0){
					evsex <- mapply(function(x,y){rbinom(1,x,y)},Nwith,sexC)
					esex <- sum(evsex)
				}
			}
	
			# Now redrawing probabilities with changed configuration
			if(esex >= 1){
				probs2 <- probset2(Na,g,mig,Nwith,Nbet,evsex,1)
				probs2 <- rbind(c(1-sum(probs2),rep(0,d-1)),probs2)
			}else if(esex == 0){
				probs2 <- probs
				probs2 <- rbind(rep(0,d),probs2)
			}
	
			# Given event happens, what is that event? Weighted average based on above probabilities.
			# Then drawing deme of event
			draw <- rmultinom(1,1,rowSums(probs2))[,1]
			event <- match(1,draw)
			draw2 <- rmultinom(1,1,probs2[event,])[,1]
			deme <- match(1,draw2)
			
			# Based on outcome, altering states accordingly
			ot <- stchange2(event,deme,evsex,Nwith)
			Nwith <- Nwith + ot$WCH
			Nbet <- Nbet + ot$BCH
			Ntot <- 2*sum(Nwith) + sum(Nbet)
			if(event == 10){	# Choosing demes to swap if there is a migration!
				if(d > 2){
					drec <- sample(setdiff(c(1:d),deme),1)
				}else if(d == 2){
					drec <- setdiff(c(1:d),deme)
				}
				draw3 <- rmultinom(1,1,c(Nwith[deme],Nbet[deme]))[,1]
				e2 <- match(1,draw3)
				if(e2 == 1){	# Paired sample migrates
					Nwith[deme] <- Nwith[deme] - 1
					Nwith[drec] <- Nwith[drec] + 1
				}else if(e2 == 2){	# Single sample migrates
					Nbet[deme] <- Nbet[deme] - 1
					Nbet[drec] <- Nbet[drec] + 1
				}
			}
			
			# Changing ancestry accordingly
			indvs <- coalesce(indvs,Ttot,Nwith,Nbet,deme,evsex,event,drec,e2)
		}
	}
	
	if(Ntot==1){	# Storing coalescent times and trees once all samples coalesce
		TCoal[i,] <- indvs[(2:Itot),4]
		treedat <- treemaker(indvs,theta)		# Makes Newick genealogy and mutation list based on output
		Twees[i] <- treedat$tree
		MTab <- treedat$muts
		MTab <- MTab[order(MTab[,1]),]	# Re-order by mutation position
		AllM[[i]] <- MTab
	}
}

## Printing output to screen and files
dir.create("Mutations", showWarnings = F)
cat("\n")
cat("Coalescent times:\n")
rownames(TCoal) <- c(1:Nreps)
colnames(TCoal) <- c(1:(Itot-1))
print(TCoal)
cat("\n")
cat("Genealogies:\n")
print(Twees)
cat("\n")
for(j in 1:Nreps){
	cat("\n")
	cat("Replicate ",j,":\n",sep='')
	cat("Mutation Positions:\n")
	MTab <- AllM[[j]]
	if(is.null(dim(MTab)) == 1){
		MTab <- t(as.matrix(MTab))
	}
	print(MTab[,1])
	cat("\n")
	cat("Individual States:\n")
	for(i in 1:Itot){
		cat(paste(MTab[,i + 1],collapse = ''))
		cat("\n")
	}
	write.table(MTab,file=paste0("Mutations/Muts_",j,".dat"),row.names=F,col.names=F,quote=F)
}

write.table(TCoal,file="CoalTimes.dat",quote=F)
write.table(Twees,file="Trees.dat",row.names=F,col.names=F)
		
# EOF