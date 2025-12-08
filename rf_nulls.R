# Function to run null models

###################################
# wmt<-data.frame(species_lo=NA,species_hi=NA,elevation=NA)
# for (i in grep("LOWER",wmtx$Level) ) for (y in 1:4) 
#   if ( (!is.na(wmtx[i,grep("D$",names(wmtx))[y]])) & (wmtx[i,grep("D$",names(wmtx))[y]] >0)  ) wmt<-rbind(wmt, c(wmtx[i,]$SpName,NA,y))
# for (i in grep("UPPER",wmtx$Level) ) for (y in 1:4) 
#   if ( (!is.na(wmtx[i,grep("D$",names(wmtx))[y]])) & (wmtx[i,grep("D$",names(wmtx))[y]] >0)  ) wmt<-rbind(wmt, c(NA,wmtx[i,]$SpName,y))
# wmt$elevation <- as.numeric(wmt$elevation)
# head(wmt)


{
  wnets0<-c(rep(NA,8))
names(wnets0)<-c("NOconnections","connectance","nestedness","NODF","SA","H2","generality.HL","vulnerability.LL")

if ( !ran ){
  pdf(paste0(onfw,"rangesx.pdf"),10,5)
  par(mfrow=c(1,2))
}


# if( ! pld ){ # loading data from original lifeweb
#   wmt<-wm[wm$dataset_id == onfw,]
#   wit<-wi[wi$dataset_id == onfw,]
#   plot(x=wit$longitude,y=wit$latitude)
#   names(wit)
#   #plot species ranges
#   wmt<-merge(wmt,wit[,c("dataset_id","network_id","elevation")],all.x=T)
#   wmt$species_lo <- gsub("_"," ",wmt$species_lo)
#   wmt$species_hi <- gsub("_"," ",wmt$species_hi)
# } 

his<-unique(wmt$species_hi); his<-his[!is.na(his)]; noh<-length(his)
los<-unique(wmt$species_lo); los<-los[!is.na(los)]; nol<-length(los)
onx=unique(sort(wmt$elevation))

# HI & LO ranges
hlra<-data.frame(hl=NA,species=NA,from=NA,to=NA)
hione<-his[1]; for(hione in his){
  hlra<-rbind(hlra, c("hi",hione,range(wmt[wmt$species_hi == hione,]$elevation,na.rm=T)))
}
hione<-los[1];for(hione in los){
  hlra<-rbind(hlra, c("lo",hione,range(wmt[wmt$species_lo == hione,]$elevation,na.rm=T)))
}
hlra<-hlra[ order(hlra$species),]
hlra$f<-NA;hlra$t<-NA
hlra<-hlra[!is.na(hlra$from),]
gr<-1
wnetsdall<-NULL
myran=1

# load(".RData")
for (myran in 0:500){ # random shifting
  print(paste(".................",myran))
  (maxra<-max(as.numeric(hlra$to),na.rm=T))
  (minra<-min(as.numeric(hlra$from),na.rm=T))
  if ( ran ){
    #shuffle ranges
    if(myran>0){ # do randomization  
      hlrai<-1; for (hlrai in 1:nrow(hlra) ){
        (onespera<-diff(as.numeric(hlra[hlrai,3:4])));hlra[hlrai,]
        
        if ( onfw == "Plowmanetalx" ) {
          #(ranran<-which( minra+100*(0:8)+onespera <=maxra) )
          (ranran<-which( minra+1*(0:8)+onespera <=maxra) )
        } else {
          (ranran<-which( minra+1*(0:3)+onespera <=maxra) )
        }
        
        (tomove<-sample(ranran-1,1))
        (tomovedir<-sample(c(-1,1),1))
        (tomovedirv<-if (tomovedir==1) { minra; } else {maxra;})
        
        if ( onfw == "Plowmanetalx" ) {
          (togo<-tomovedirv+tomovedir*100*tomove)
        }else {
          (togo<-tomovedirv+tomovedir*1*tomove)
        }
        
        hlra[hlrai,5:6] <- sort(c(togo,togo+tomovedir*onespera))
      }
      hlra$from <- hlra$f;hlra$to <- hlra$t
    }
  } else { # !ran, no randomization
    plot((seq(0,nrow(hlra),length=length(onx)))~onx,type='n',main="HI and LO species", xlab="gradient", ylab="species")
    mtext(paste("HI n=",noh,sep=""),3,1.5,col="green",adj=1)
    mtext(paste("LO n=",nol,sep=""),3,0,col="magenta",adj=1)
    mtext(onfw,3,1,adj=0,cex=1.2)
    #mtext(paste("HI mid att peak=",ma_a,", SD=",ma_b, sep=""),3,.5,adj=1)
    for (i in 1:nrow(hlra)) {
      lines(x=c(hlra[i,]$from,hlra[i,]$to),y=c(i,i),col=ifelse(hlra[i,]$hl =="hi","green","magenta"),lwd=2 )
      text(x=gr-0.1, y=i,labels=hlra$species[i],cex=0.3, col="black")#ifelse(hlra[i,]$hl =="hi","green","magenta"))
    }
    onx
    table (as.numeric(hlra$to) - as.numeric(hlra$from))
  }
  
  if ( pld ) {
    #hlra<-hlra[-c(nrow(hlra)),]; 
    hlra$from <- as.numeric(hlra$from); hlra$to <- as.numeric(hlra$to)
    
    (allgr<-sort(unique(c( hlra$from,hlra$to) )))
    for(alg in allgr) hlra[paste0("gr",alg)]<-FALSE
    for(algi in 1:nrow(hlra)) {
      hlra[algi,paste0("gr",allgr)] <- allgr %in% seq(hlra[algi,]$from,hlra[algi,]$to)
    }
    hlraA<-hlra
    
    #plow
    if ( onfw == "Plowmanetalx" ) {
      #gr<-700; grstep<-100
      #ongr<- c(700+grstep*c(0:7))
      gr<-1;grstep<-1
      ongr<-c(1+grstep*c(0:(mygr-1)))
    } else {
      #rob
      gr<-1;grstep<-1
      ongr<-c(1+grstep*c(0:(mygr-1)))
    }
    
    wnetsd<-as.data.frame(wnets0)
    for ( gr in ongr ) {
      #for (gr in c(1+grstep*c(0:(mygr-1)))) {
      
      print(gr)
      # are species in gradient bin?
      
      more2 <- FALSE
      if ( nrow(hlraA) > 2 ) {
        # 1-2, 2-3, 3-4
        #hlrai<-sapply(hlraA$to, '>=', gr) & sapply(hlraA$from, '<=', gr+grstep)
        # 1,2,3,4
        hlrai<-hlraA[,paste0("gr",gr)]
        
        # select species in gradient bin
        hlrah<-hlra[ hlrai & hlra$hl == "hi",]; hlral<-hlra[ hlrai & hlra$hl == "lo",];
        more2 <-TRUE
      }
      if( more2 & (nrow(hlrah) >1) & (nrow(hlral) >1) ){
        # have species overlaping ranges
        #ra_hl <- sapply(hlrah$to, '>=', hlral$from) & sapply(hlrah$from, '<=', hlral$to)
        ra_hl <- sapply(hlrah[,paste0("gr",gr)], '&', hlral[,paste0("gr",gr)]) 
        xl<-xh<-5; # xh=4;xl=1
        
        # how many times hi and lo species were together and in the gradient bin
        #plowman
        if ( onfw == "Plowmanetal" ) {
          xl=1;xh=1
          # how many times where observed together 
          for( xh in 1:nrow(hlrah) ) for( xl in 1:nrow(hlral) )
            ra_hl[xl,xh] <- ifelse ( any( (hlrah[xh,]$species == wmt$species_hi) & 
                                            (hlral[xl,]$species == wmt$species_lo) )  ,
                                     sum( (hlrah[xh,]$species == wmt$species_hi) & 
                                            (hlral[xl,]$species == wmt$species_lo)),# & 
                                     #( (wmt$elevation ==gr) | (wmt$elevation == (gr)))  ),
                                     0 )
        } else {
          #krkonose
          for(xh in 1:nrow(hlrah) ) for(xl in 1:nrow(hlral) )
            ra_hl[xl,xh] <- ifelse ( ra_hl[xl,xh] 
                                     & length(wmt[ (wmt$species_hi == hlrah[xh,]$species) 
                                                   & ((wmt$species_lo == hlral[xl,]$species)),]$value )>0,
                                     wmt[ (wmt$species_hi == hlrah[xh,]$species) 
                                          & ((wmt$species_lo == hlral[xl,]$species)),]$value,
                                     0)
        }
        
        #ra_hl[xl,xh] <- ifelse ( any( (hlrah[xh,]$species == wmt$species_hi) & 
        #                          (hlral[xl,]$species == wmt$species_lo) )  ,
        #                          sum( (hlrah[xh,]$species == wmt$species_hi) & 
        #                               (hlral[xl,]$species == wmt$species_lo) ),
        #                          ra_hl[xl,xh])
        
        nol<-nrow(hlral);noh<-nrow(hlrah)
        hilo <- matrix(ra_hl,nrow=nol, ncol=noh, dimnames = list(c(1:nol),c(1:noh)))
        
        library(bipartite)
        
        if ( (nrow(hilo) >= 2) & (ncol(hilo) >= 2) ) {
          ( wnets <- networklevel(hilo,index=
                                    c("H2","generality","vulnerability","connectance","nestedness",
                                      "nestedness","NODF","SA")
          ) )
          wnetsd<-cbind(wnetsd, c(sum(hilo>0,na.rm=T),wnets) )
        } else {
          wnetsd<-cbind(wnetsd, c(0,rep(NA,7)) )
        }
      } else {
        wnetsd<-cbind(wnetsd, c(0,rep(NA,7)) )
      }
      
    }
    #
  }
  
  wnetsd<-wnetsd[,-c(1)]
  
  if (is.null(wnetsdall)){
    wnetsdall<-as.data.frame.matrix(t(as.matrix(wnetsd)))
    wnetsdall$gr<-1:(mygr1)
    wnetsdall$null<-0
  } else {
    wnetsdallx<-as.data.frame.matrix(t(as.matrix(wnetsd)))
    wnetsdallx$gr<-1:(mygr1)
    wnetsdallx$null<-myran
    wnetsdall<-rbind(wnetsdall,wnetsdallx)
  }
}
write.csv(wnetsdall,paste0("realfood/netsgrad1000noi_",onfw,"x1.csv"),row.names=F)
#write.csv(wnetsdall,"netsgrad1000noi_krko.csv",row.names=F)
}
