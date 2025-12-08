rm(list=ls())


myseed<-12345; set.seed(myseed) # seed for random generator
noh <- 200 # no HI species - eaters
nol <- 100 # no LO species - trees
nos <- 100 # no gradient points
rangereds<-c(1,2,4); rangered<-1
grstep<-10 # step for gradient binning (for how big interval is network statistic computed)
grstep2<-10 # step size

dropovers<-c("cut","border","regen")
dropover<-dropovers[3]

dospecs<-c("NOSPEC","DECREASING","PEAK","USHAPE")#"peak" # sort species according specialization (lower are generalist, upper specialists)
dospec<-dospecs[1]
maxlo <- 82 # maximum "
minlomi<-2; minloma<-maxlo # minimum of LO species that are HI specialized to -range
specstep<-20
minlo <- 18 # minimum of LO species that are HI specialized to -val (changed in loop)

noreps <- 20 # number of replicates -range 1:noreps
rep<-1 # replicate number - value
domodloop <- TRUE # if to run modelling loop

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


pdf (paste("gr_hilo_spe_",myseed,"_sche5.pdf",sep=""),width=13,height=7)
#par(mfcol=c(2,3))
#4 figs in one
#layout(matrix(c(1,1,2,2,3,3,3,4), 2, 4, byrow = TRUE),
#       widths=c(1,0.6,0.6,1), heights=c(1,1))
# 6 figs in one
par(mar=c(5,4.2,4,0.3))
layout(matrix(c(1,1,2,2,3,4, 5,5,5,6,6,7), 2, 6, byrow = TRUE),
       widths=c(1,1,1, 1,1,1), heights=c(1,1))

onx=1:(nos)

simsuf<-paste0("9_",dospec,"_",grstep,"grad_",rangered,"rangered",noreps,"reps")

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

# data frame for the results
fw_res <- data.frame(fw=NA, rep=NA, minlo=NA, gr=NA, allInts=NA, LOints=NA,HIints=NA,
                     connectance=NA, comps=NA, comps_div=NA, nestedness=NA, H2=NA, generality=NA, vulnerability=NA)

minlo
spel_spec <- data.frame( id=1:noh )
maxspec<- round(100*maxlo/nol,0); minspec<- round(100*minlo/nol,0)
for (i in 1:maxspec) spel_spec[[as.character(i)]] <-NA
# specialization decreasing according start of LO species (lo species lower had lower specialization)
(nohs <- round( runif(noh,min=minspec,max=maxspec),0 )) # for LO

spel_spec <- spel_spec[,-c(1)]
for (i in 1:noh) spel_spec[i,] <- c( sample(1:nol, nohs[i]), rep(NA,maxspec-nohs[i]) )

# HI
plot((seq(0,noh-2,by=2))~onx,type='n',main="", xlab="Artificial gradient", ylab="Higher level species")
axis(1,seq(10,90,by=20))
mtext("(a) Higher level species",3,1,adj=0)
#mtext(paste0("Range length ",round(100/as.numeric(rangered),0)) ,3,2,adj=0,cex=1.2)
mtext(paste("n=",nrow(speh),sep=""),3,0.5,adj=1)
#mtext(paste("HI mid att peak=",ma_a,", SD=",ma_b, sep=""),3,.5,adj=1)
abline(v=c(11,20),col="black",lwd=1.5)
#text(15,-15,"bin",xpd=T)
for (i in 1:nrow(speh)) lines(x=c(speh[i,]$s,speh[i,]$e),y=c(i,i),
                              col=ifelse( speh[i,]$e<=10 ,"gray",ifelse (speh[i,]$s > 20,"gray","red") ))
lines(y=crah,x=1:nos,col="darkblue"); lines(y=crah2,x=1:nos,col="darkred")


# LO
plot((seq(0,nol-1,by=1))~onx,type='n',main="", xlab="Artificial gradient", ylab="Lower level species")
axis(1,seq(10,90,by=20))
mtext("(b) Lower level species",3,1,adj=0)
mtext(paste("n=",nrow(spel),sep=""),3,0.5,adj=1)
#mtext(paste("LO mid att peak=",mb_a,", SD=",mb_b, sep=""),3,.5,adj=1)
abline(v=c(11,20),col="black",lwd=1.5)
#text(15,-15,"bin",xpd=T)
for (i in 1:nrow(spel)) lines(x=c(spel[i,]$s,spel[i,]$e),y=c(i,i),
                              col=ifelse( spel[i,]$e<=10 ,"gray",ifelse (spel[i,]$s > 20,"gray","blue") ))
lines(y=cral,x=1:nos,col="darkblue")

# histograms
hist(speh$range,xlab="Species gradient range",main="",col="gray")#;hist(spel$range,xlab="gradient range",main="",col="gray");
mtext("(c) Range",3,1,adj=0)
hist(speh$s +(speh$range/2),xlab="Species gradient center",main="",col="gray30",ylab="")
mtext("(d) Center",3,1,adj=0)


# HI & LO int + spec
gr<-20
ra_hl0 <- sapply(speh$e, '>=', spel$s) & sapply(speh$s, '<=', spel$e) # just to create matrix
(gron<-(1:nos)[(gr-grstep+1):gr])

#if ( all(ra_hl<2) ) {
ra_hl <- ra_hl0
ra_hli <- ra_hl0
# speciealized + length of overlap
for ( i in 1:nol ) for( y in 1:noh ) 
  ra_hl[i,y] <- if (i %in% c(unlist(spel_spec[y,])) ) sum(crl[gron,i+1]*crh[gron,y+1]) else 0
# only if they are specialized
for ( i in 1:nol ) for( y in 1:noh ) 
  ra_hli[i,y] <- if (i %in% c(unlist(spel_spec[y,])) ) 1 else 0
#}


dim(ra_hl)
#nrow(ra_hl); colSums(ra_hl)
hilo <- matrix(ra_hl,nrow=nol, ncol=noh, dimnames = list(c(1:nol),c(1:noh)))
#library(bipartite)
#networklevel(hilo>0,index="generality")

(grayc<-paste0("gray",rev( seq(9,95,by=8) ) ))
bluec<-colorRampPalette(c("lightgreen","darkgreen"))
(grayc<-c("lightblue",bluec(10)))
plot(y=seq(1,nol,length=noh),x=1:noh,ylab="Lower level species", xlab="Higher level species",
     type='n',axes=F)
axis(1);axis(2)
mtext("(e) Species specialization and overlap",3,1,adj=0)
xi=1;yi=1
for(xi in 1:nol) { 
  points(y=rep(xi,noh),x=1:noh, col=grayc[ 1+ra_hl[xi,] ], cex=0.5, pch=16)
  #points(y=rep(xi,noh),x=1:noh,col=ifelse(ra_hl[xi,]==0,"gray","green"),
  #                      cex=ifelse(ra_hl[xi,]==0,0.2,0.6),pch=16 )
  #text(y=rep(xi,noh),x=1:noh,labels=ifelse(ra_hl[xi,]==0,"", ra_hl[xi,]),cex=0.5)
}
for(xi in 1:nol) { 
  points(y=rep(xi,noh),x=1:noh, col=ifelse(ra_hli[xi,]==0,"gray95","NA"), cex=0.5, pch=16)
}

legend(x=5,y=112,c("not specialized, specialized with overlaps:",0,5,10),
       text.width=c(98,5,5,5),
       adj=0,pch=16,col=c("gray95",grayc[c(0,5,10)+1]),xpd=T,horiz = T,pt.cex=3,bty="n")

#lo 4 
plot( seq(1,8,length=30)~c(1:30),type='n',main="", xlab="Artificial gradient", ylab="", axes=F)
mtext("(f) Example of overlaps",3,1,adj=0)
axis(1)
abline(v=c(11,20),col="black",lwd=1.5)
#lo
onlo<-rowSums( (ra_hl>0) & (ra_hl<10) )
i=which(onlo>0)[1]; lines(x=c(spel[i,]$s,spel[i,]$e),y=c(1,1), col="blue",lwd=1.5)
text(x=0,y=1+0.4,labels=paste0("LO #",i),adj=0,col="blue")
#hi
onhi<-which( (ra_hl[i,]>0) & (ra_hl[i,]<10))
onhi[7]<-which( (ra_hl[i,]>0) & (ra_hl[i,]==10))[1]
points(x=rep(15,6),y=0.3+2:7,col=grayc[1+ (ra_hl[i,])[onhi[1:6]] ],pch=16,cex=3,xpd=T)
points(x=15,y=(0.3+2:8)[7],col="gray",pch=16,cex=3,xpd=T)
text(x=15,y=(0.3+2:8)[7],"0",xpd=T)

ii<-1; for (i in onhi[1:6]) {
  ii<-ii+1;
  lines(x=c(speh[i,]$s,speh[i,]$e),y=c(ii,ii),col="red")
  text(x=25,y=ii+0.4,labels=paste0("HI #",i),adj=0,col="red",xpd=T)
  text(x=15,y=ii+0.3, labels=ra_hl[which(onlo>0)[1], i],xpd=T )
}

ii<-ii+1;
lines(x=c(1,8),y=c(ii,ii),col="red")
text(x=25,y=ii+0.4,labels=paste0("HI #199"),adj=0,col="red",xpd=T)

library(bipartite)
par(mar=c(0,0,0,0) )
plotweb(hilo[1:7,onhi[1:5]])
mtext("(g) Network",3,1,adj=0)

dev.off()
