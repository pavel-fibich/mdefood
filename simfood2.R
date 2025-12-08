rm(list=ls())

if (FALSE) { #just for git commits/pulls
  library(gitcreds)
  system("cat ~/.Renviron")
  gitcreds_set() # paste token
}

# TO TRY -  pridat no spe a velikost matice, vary ranges
##
# static settings - may be varied later, check loop
outfilep<-"mdeFeas2_"
myseed<-12345; set.seed(myseed) # seed for random generator
noh <- 200 # no HI species - eaters
nol <- 100 # no LO species - trees
nos <- 100 # no gradient points
rangereds<-c(1,2,4)
rangered<-rangereds[1]

dropovers<-c("cut","border","regen","feasable")
dropover<-dropovers[4]
feasmore<-4 # generate 5 times more species to get enought feasable ones

# steps for binning
grstep<-10 # step for gradient binning (for how big interval is network statistic computed)
grstep2<-10 # step size
grstep3<-5 # second stepsize
grstep4<-20 # third stepsize

dospecs<-c("NOSPEC","DECREASING","PEAK","USHAPE")#"peak" # sort species according specialization (lower are generalist, upper specialists)
dospec<-dospecs[1]

# minimum and maximum of LO species that are HI specialized to (changed in loop)
speclo<-c(2,20,50)
spechi<-c(20,50,90)

noreps <- 500 # number of replicates -range 1:noreps
rep<-1 # replicate number - value
domodloop <- TRUE # if to run modelling loop

# funciton to make sort in peak - biggest values in the middle
makeTent = function(ints) {
  ints_o = ints[order(ints)]
  if((length(ints) %% 2) == 0) {
    # even number of observations
    ints_tent = c(ints_o[seq.int(from = 1, to = (length(ints) - 1), by = 2)],
                  rev(ints_o[seq.int(from = 2, to = length(ints), by = 2)]))
  } else {
    # odd number of observations
    ints_tent = c(ints_o[seq.int(from = 2, to = (length(ints) - 1), by = 2)],
                  rev(ints_o[seq.int(from = 1, to = length(ints), by = 2)]))
  }
  return(ints_tent)
}

###
library(bipartite)
library(foreach)
library(doParallel)
cl <- makeCluster(3) #not to overload your computer
registerDoParallel(cl)
##for(rangered in rangereds) 
rangered <- rangereds[2]
foreach( rangered=rangereds ) %dopar% {  
  for( dospec in dospecs[1:2]) {
    
    library(bipartite)
    # name used for storing the results  
    simsuf<-paste0("20_",dospec,"_",grstep,"grad_",rangered,"rangered",noreps,"reps")
    
    # uniformly generated ranges
    noh_r <- round( runif(noh,min=1,max=round(nos/rangered,0)),0 ) # for HI
    nol_r <- round( runif(nol,min=1,max=round(nos/rangered,0)),0 ) # for LO
    # mid domain centers
    noh_c <- round( runif(noh,min=1,max=nos),0 ) 
    nol_c <- round( runif(nol,min=1,max=nos),0 ) 
    #mid attractor for LO
    #nol_c <- rtrunc(nol,"norm",a = 1,b=nos,mean=mb_a,sd=mb_b)
    
    # species data frames - ifelse logic for the scenarios where species have 
    # the range overlaping border of gradient
    if ( dropover == "cut" ) { # drop species overlapping border
      speh <- data.frame(id=1:noh,range=noh_r, s=noh_c-noh_r/2, e=noh_c+noh_r/2 )
      spel <- data.frame(id=1:nol,range=nol_r, s=nol_c-nol_r/2, e=nol_c+nol_r/2 )
      
      speh$s <- round(speh$s,0); speh$e <- round(speh$e,0)
      spel$s <- round(spel$s,0); spel$e <- round(spel$e,0)
      
      speh$s <- ifelse( speh$s <1,1, speh$s); speh$e <- ifelse( speh$e >100,100, speh$e);  
      spel$s <- ifelse( spel$s <1,1, spel$s); spel$e <- ifelse( spel$e >100,100, spel$e);  
      
    } else if (dropover == "border") { # more overlaps inside gradient
      speh <- data.frame(id=1:noh,range=noh_r,
                         s=ifelse( (noh_c-noh_r/2)<1,1,
                                   ifelse((noh_c+noh_r/2)>nos, (noh_c-noh_r/2)-(noh_c+noh_r/2)+nos ,noh_c-noh_r/2 ) ),
                         e=ifelse((noh_c+noh_r/2)>nos,nos, 
                                  ifelse( (noh_c-noh_r/2)<1,(noh_r/2+noh_c)-(noh_c-noh_r/2)+1,noh_r/2+noh_c) ))
      spel <- data.frame(id=1:nol,range=nol_r,
                         s=ifelse( (nol_c-nol_r/2)<1,1,
                                   ifelse((nol_c+nol_r/2)>nos, (nol_c-nol_r/2)-(nol_c+nol_r/2)+nos ,nol_c-nol_r/2 ) ),
                         e=ifelse((nol_c+nol_r/2)>nos,nos, 
                                  ifelse( (nol_c-nol_r/2)<1,(nol_r/2+nol_c)-(nol_c-nol_r/2)+1,nol_r/2+nol_c) ))
    } else if ( dropover == "feasable" ) { # only feasable species goes in
      
      # uniformly generated ranges
      noh_r <- round( runif(noh*feasmore,min=1,max=round(nos/rangered,0)),0 ) # for HI
      nol_r <- round( runif(nol*feasmore,min=1,max=round(nos/rangered,0)),0 ) # for LO
      # mid domain centers
      noh_c <- round( runif(noh*feasmore,min=1,max=nos),0 ) 
      nol_c <- round( runif(nol*feasmore,min=1,max=nos),0 ) 
      
      speh <- data.frame(id=1:noh,range=noh_r, s=noh_c-noh_r/2, e=noh_c+noh_r/2 )
      spel <- data.frame(id=1:nol,range=nol_r, s=nol_c-nol_r/2, e=nol_c+nol_r/2 )
      
      good_h<- (speh$s>=1) & (speh$e<=100); #table(good_h);noh
      good_l<- (spel$s>=1) & (spel$e<=100); #table(good_l);nol
      
      speh<-speh[good_h,]; speh<-speh[1:min(nrow(speh),noh),]
      spel<-spel[good_l,]; spel<-spel[1:min(nrow(spel),nol),]
      #hist(speh$range);hist((speh$e+speh$s)/2)
      
    } else { # renerate overlapers but leave the range untouched
      speh <- data.frame(id=1:noh,range=noh_r, s=noh_c-noh_r/2, e=noh_c+noh_r/2 )
      spel <- data.frame(id=1:nol,range=nol_r, s=nol_c-nol_r/2, e=nol_c+nol_r/2 )
      
      speh$s <- round(speh$s,0); speh$e <- round(speh$e,0)
      spel$s <- round(spel$s,0); spel$e <- round(spel$e,0)
      
      spehToMove<-speh$s <1 | speh$e >100;table(spehToMove)
      speh[spehToMove,]
      for (si in speh[spehToMove,1] ) {
        si_range<-speh[speh[,1] == si,]$range
        si_dif<-nos-si_range
        if (si_dif == 0) {
          speh[speh[,1] == si,3:4]<-c(1,100)
        } else {
          si_newstart<-round( runif(1,min=1,max=si_dif),0 )
          speh[speh[,1] == si,3:4]<-c(si_newstart,si_newstart + si_range)
        }
      }
      #hist(speh$s +(speh$range/2) )
      spelToMove<-spel$s <1 | spel$e >100;table(spelToMove)
      spel[spelToMove,]
      for (si in spel[spelToMove,1] ) {
        si_range<-spel[spel[,1] == si,]$range
        si_dif<-nos-si_range
        if (si_dif == 0) {
          spel[spel[,1] == si,3:4]<-c(1,100)
        } else {
          si_newstart<-round( runif(1,min=1,max=si_dif),0 )
          spel[spel[,1] == si,3:4]<-c(si_newstart,si_newstart + si_range)
        }
      }
      
    }
    
    # just rounding to integer values
    speh$s <- round(speh$s,0); speh$e <- round(speh$e,0)
    spel$s <- round(spel$s,0); spel$e <- round(spel$e,0)
    
    # species number along gradient
    # colSums over gradient - species number
    #hi
    crh<-data.frame(gr=1:nos); for (i in 1:nrow(speh)) crh[[as.character(i)]] <- NA
    for (i in 1:nos) crh[i,2:ncol(crh)] <- (speh$s <= i) & (speh$e >=i)
    crah2 <- crah <- rowSums(crh[,2:ncol(crh)])
    #crah <- c(); for (i in 1:nos) crah<-c(crah, sum( ((speh$s <= i) & (speh$e >=i))))# + ((speh$s2 <=i) & (speh$e2 <=i)))  )
    #cral <- c(); for (i in 1:nos) cral<-c(cral, sum( ((spel$s <= i) & (spel$e >=i)))) #+ ((spel$s2 <=i) & (spel$e2 <=i)))  )
    #lo
    crl<-data.frame(gr=1:nos); for (i in 1:nrow(spel)) crl[[as.character(i)]] <- NA
    for (i in 1:nos) crl[i,2:ncol(crl)] <- (spel$s <= i) & (spel$e >=i)
    cral <- rowSums(crl[,2:ncol(crl)])
    
    #fast plot
    if (FALSE) {
      pdf (paste("gr_hilo_spe_",myseed,"_4oldborder.pdf",sep=""),6,14)
      par(mfcol=c(4,1))
      onx=1:(nos)
      # HI
      plot((seq(0,noh-2,by=2))~onx,type='n',main="", xlab="Environmental gradient", ylab="Higher level species")
      mtext(paste0("Range length ",nos/rangered),3,1,adj=0,cex=1.2)
      mtext(paste("n=",nrow(speh),sep=""),3,1.5,adj=1)
      for (i in 1:nrow(speh)) lines(x=c(speh[i,]$s,speh[i,]$e),y=c(i,i),col="gray" )
      lines(y=crah,x=1:nos,col="darkblue"); lines(y=crah2,x=1:nos,col="darkred")
      # LO
      plot((seq(0,nol-1,by=1))~onx,type='n',main="", xlab="Environmental gradient", ylab="Lower level species")
      mtext(paste("n=",nrow(spel),sep=""),3,1.5,adj=1)
      for (i in 1:nrow(spel)) lines(x=c(spel[i,]$s,spel[i,]$e),y=c(i,i),col="gray25" )
      lines(y=cral,x=1:nos,col="darkblue")
      #
      hist(speh$range,xlab="HI species gradient range",main="",col="gray")#;hist(spel$range,xlab="gradient range",main="",col="gray");
      hist(speh$s +(speh$range/2),xlab="HI species gradient center",main="",col="gray")
      dev.off()
    }
    
    # data frame for the results
    fw_res <- data.frame(fw=NA, rep=NA, minlo=NA, gr=NA, grstep=NA, allInts=NA, LOints=NA,HIints=NA,
                         connectance=NA, comps=NA, comps_div=NA, nestedness=NA, NODF=NA, spec_asym=NA,
                         H2=NA, generality=NA, vulnerability=NA)
    
    speci <- 1
    for(speci in 1:length(speclo) ) {
      minlo<-speclo[speci]; maxlo<-spechi[speci]
      
      # generate HI species for that are LO species specialized - just numbers
      spel_spec <- data.frame( id=1:noh )
      maxspec<- round(100*maxlo/nol,0); minspec<- round(100*minlo/nol,0)
      for (i in 1:maxspec) spel_spec[[as.character(i)]] <-NA
      # specialization decreasing according start of LO species (lo species lower had lower specialization)
      (nohs <- round( runif(noh,min=minspec,max=maxspec),0 )) # for LO
      
      if ( dospec == "DECREASING") { # sorted, generalists at the bottom, specialists at the top
        #nohs.d <- data.frame(x=(speh$s+speh$range/2) )
        nohs.d <- data.frame(x=speh$s )
        nohs.d$xo <- order( nohs.d$x ,decreasing=F )
        nohs.d$xx <- as.numeric(rownames(nohs.d))
        nohs.d <- nohs.d[nohs.d$xo,]
        nohs.d$nohs <- sort(nohs,decreasing = T)
        nohs.d <- nohs.d[order(nohs.d$xx,decreasing = F),]
        # plot(nohs.d$x~nohs.d$nohs);  nohs.d[1:10,]
        
        nohs<-nohs.d$nohs
      } else if (dospec =="PEAK"){
        nohs.d <- data.frame(x=speh$s )
        nohs.d$xo <- order( nohs.d$x ,decreasing=F )
        nohs.d$xx <- as.numeric(rownames(nohs.d))
        nohs.d <- nohs.d[nohs.d$xo,]
        nohs.d$nohs <- makeTent(nohs)#sort(nohs,decreasing = T)
        nohs.d <- nohs.d[order(nohs.d$xx,decreasing = F),]
        
        # plot(nohs.d$x~nohs.d$nohs);  nohs.d[1:10,]
        nohs<-nohs.d$nohs
      } else if (dospec =="USHAPE"){
        nohs.d <- data.frame(x=speh$s )
        nohs.d$xo <- order( nohs.d$x ,decreasing=F )
        nohs.d$xx <- as.numeric(rownames(nohs.d))
        nohs.d <- nohs.d[nohs.d$xo,]
        nohs.d$nohs <- makeTent(nohs)[c(100:200,1:99)]#makeTent(nohs)#sort(nohs,decreasing = T)
        nohs.d <- nohs.d[order(nohs.d$xx,decreasing = F),]
        
        # plot(nohs.d$x~nohs.d$nohs);  nohs.d[1:10,]
        nohs<-nohs.d$nohs
      }
      spel_spec <- spel_spec[,-c(1)]
      
      if (domodloop) for(rep in 1:noreps) {
        
        # generate HI species for that are LO species specialized - exactly choose HI species
        for (i in 1:noh) spel_spec[i,] <- c( sample(1:nol, nohs[i]), rep(NA,maxspec-nohs[i]) )
        
        if (FALSE){ # occasional plotting
          
          pdf (paste("gr_hilo_spe_",myseed,"_3.pdf",sep=""),12,6)
          par(mfcol=c(2,3))
          onx=1:(nos)
          # HI
          plot((seq(0,noh-2,by=2))~onx,type='n',main="", xlab="Environmental gradient", ylab="Higher level species")
          mtext(paste0("Range reduction by ",rangered),3,1,adj=0,cex=1.2)
          mtext(paste("n=",noh,sep=""),3,1.5,adj=1)
          #mtext(paste("HI mid att peak=",ma_a,", SD=",ma_b, sep=""),3,.5,adj=1)
          for (i in 1:nrow(speh)) lines(x=c(speh[i,]$s,speh[i,]$e),y=c(i,i),col="gray" )
          #colsums
          lines(y=crah,x=1:nos,col="darkblue"); lines(y=crah2,x=1:nos,col="darkred")
          
          # LO
          plot((seq(0,nol-1,by=1))~onx,type='n',main="", xlab="Environmental gradient", ylab="Lower level species")
          mtext(paste("n=",nol,sep=""),3,1.5,adj=1)
          #mtext(paste("LO mid att peak=",mb_a,", SD=",mb_b, sep=""),3,.5,adj=1)
          for (i in 1:nrow(spel)) lines(x=c(spel[i,]$s,spel[i,]$e),y=c(i,i),col="gray25" )
          #colsums
          lines(y=cral,x=1:nos,col="darkblue")
          #
          hist(speh$range,xlab="gradient range",main="",col="gray");hist(spel$range,xlab="gradient range",main="",col="gray");
          
          dev.off()
        }
        
        # overal res
        # does LO and HI overlap?
        #ra_hl <- sapply(speh$e, '>=', spel$s) & sapply(speh$s, '<=', spel$e)
        ## ion<-9;speh[ion,]; spel[ion,]; spel_spec[ion,]; ion %in% spel_spec[ion,]; ra_hl[ion,ion]
        ## how big is their overlap
        ##for ( i in 1:nol ) for( y in 1:noh ) ra_hl[i,y] <- if (i %in% spel_spec[y,]) sum(crl[,i+1]*crh[,y+1]) else 0
        ##hilo <- matrix(ra_hl,nrow=nol, ncol=noh, dimnames = list(c(1:nol),c(1:noh)))
        ##library(bipartite)
        ##( wnets <- networklevel(hilo,index=c("H2","generality","vulnerability","connectance","number of compartments","compartment diversity",
        #                                     "nestedness")) )
        #simcode=paste("alonggr","",sep="")
        #fw_res = rbind(fw_res, c(simcode,rep, minlo, "all", wnets) )
        
        ########
        ######## gradient loop
        ######## 
        gr<-grstep
        simcode=paste("alonggr","",sep="")
        ra_hl0 <- sapply(speh$e, '>=', spel$s) & sapply(speh$s, '<=', spel$e) # just to create matrix
        for ( gr in seq(grstep,nos,by=grstep2) ) {
          
          (gron<-(1:nos)[(gr-grstep+1):gr])
          ra_hl <- ra_hl0
          for ( i in 1:nol ) for( y in 1:noh ) 
            ra_hl[i,y] <- if (i %in% c(unlist(spel_spec[y,])) ) sum(crl[gron,i+1]*crh[gron,y+1]) else 0
          
          #nrow(ra_hl); colSums(ra_hl)
          hilo <- matrix(ra_hl,nrow=nol, ncol=noh, dimnames = list(c(1:nol),c(1:noh)))
          # check
          #hilo[1,]; table(hilo[1,]); 
          #xon=149;hilo[1,xon];speh[xon,];spel[1,]
          #xon=196;hilo[1,xon];speh[xon,];spel[1,]
          
          if ( any(hilo>0) ) {
            # get network measures
            ( wnets <- tryCatch(networklevel(hilo,index=
                                               c("H2","generality","vulnerability","connectance","number of compartments","compartment diversity",
                                                 "nestedness","NODF","SA")), error=function(e) return(e)) 
            )
            if ( length(wnets)>7 ) {
              fw_res = rbind(fw_res, c(simcode,rep,minlo, gr,grstep,sum(hilo!=0), sum(rowSums(hilo!=0)>0), sum(colSums(hilo!=0)>0), wnets) ) # save it the results data frame
            }
          }
        } # end of gradient loop
        
        # second step loop
        for ( gr in seq(grstep3,nos,by=grstep3) ) {
          (gron<-(1:nos)[(gr-grstep3+1):gr])
          ra_hl <- ra_hl0
          for ( i in 1:nol ) for( y in 1:noh ) 
            ra_hl[i,y] <- if (i %in% c(unlist(spel_spec[y,])) ) sum(crl[gron,i+1]*crh[gron,y+1]) else 0
          
          hilo <- matrix(ra_hl,nrow=nol, ncol=noh, dimnames = list(c(1:nol),c(1:noh)))
          # get network measures
          if ( any(hilo>0) ) {
            #( wnets <- networklevel(hilo,index=
            #                        c("H2","generality","vulnerability","connectance","number of compartments","compartment diversity",
            #                          "nestedness","NODF","SA")
            #) )
            
            wnets <- tryCatch(networklevel(hilo,index=
                                             c("H2","generality","vulnerability","connectance","number of compartments","compartment diversity",
                                               "nestedness","NODF","SA")), error=function(e) return(e))
            if ( length(wnets)>7 ) {
              fw_res = rbind(fw_res, c(simcode,rep,minlo, gr,grstep3,sum(hilo!=0), sum(rowSums(hilo!=0)>0), sum(colSums(hilo!=0)>0), wnets) ) # save it the results data frame
            }
          }
        } # end of gradient loop
        
        # third step loop
        for ( gr in seq(grstep4,nos,by=grstep4) ) {
          (gron<-(1:nos)[(gr-grstep4+1):gr])
          ra_hl <- ra_hl0
          for ( i in 1:nol ) for( y in 1:noh ) 
            ra_hl[i,y] <- if (i %in% c(unlist(spel_spec[y,])) ) sum(crl[gron,i+1]*crh[gron,y+1]) else 0
          
          hilo <- matrix(ra_hl,nrow=nol, ncol=noh, dimnames = list(c(1:nol),c(1:noh)))
          # get network measures
          if ( any(hilo>0) ) {
            ( wnets <- tryCatch(networklevel(hilo,index=
                                               c("H2","generality","vulnerability","connectance","number of compartments","compartment diversity",
                                                 "nestedness","NODF","SA")), error=function(e) return(e)))
            if ( length(wnets)>7 ) {
              fw_res = rbind(fw_res, c(simcode,rep,minlo, gr,grstep4,sum(hilo!=0), sum(rowSums(hilo!=0)>0), sum(colSums(hilo!=0)>0), wnets) ) # save it the results data frame
            }
          }
        } # end of gradient loop
        
      } # end of reps loop
      
    } # end of minlo
    fw_res <- fw_res[-c(1),]
    write.csv(fw_res,paste(outfilep,myseed,"_",simsuf,".csv",sep=""),row.names=F)
    
  } } # end of two for loops
save.image()
stopCluster(cl = cl)
