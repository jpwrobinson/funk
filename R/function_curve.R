

## adapted from Fabio Pranovi November 2015

mycast<-function(biom,step,c0,c1,y0,y1){
  # matteo zucchetta
  # 15 maggio  2012
  # funzione che calcola le biomasse cumulate
  
  # biom : tabella biomasse (deve avere SP=specie; TL=livello trofico)
  # step : ampiezza delle classi
  # c0 : prima colonna biomasse
  # c1 : ultima colonna biomasse
  # y0 : primo anno
  # y1 : ultimo anno
  
  require(reshape)
  
  
  # definisco le classi
  my.lim<-seq(1-step,7+step,by=step) ### AH changed this form 5 to 7 to account for the much higher CV of trophic levels
  my.lab<-my.lim+step/2
  a<-cut( biom$TL,my.lim,labels=my.lab[-length(my.lim)])

  # creo un nuovo dataframe
  my.elab<-data.frame(sp=biom$SP,aveTL=biom$TL)
  my.elab$fTL<-cut(my.elab$aveTL,my.lim,labels=my.lab[-length(my.lim)])
  my.elab$fTLr<-as.numeric(levels(my.elab$fTL)[my.elab$fTL])
  
  B<-biom[,c(c0:c1)] 
  names(B)<-c(y0:y1)
  Btot<-cbind(my.elab,B)
  
  
  my.melt<-melt(Btot,id.vars=c("sp","aveTL","fTL","fTLr"))
  my.cast<-cast(my.melt,fTLr~variable,sum,add.missing=T,fill=0) 
  return(my.cast)
}


##################################

my.plot<-function(x,y){
  # matteo zucchetta
  # 15 maggio  2012
  # funzione che fitta le curve e ne calcola i parametri
  
  # x : TL
  # y : biomasse
  
  require(drc)
  
  
  if(sum(y)==0){ # se non ci sono dati non plotto la curva
    plot(x,y,ylim=c(0,1),xlim=c(2,6),xlab="",ylab="",xaxt="n",yaxt="n")
    text(2.4,0.8,paste(year))
    return(c(NA,NA,NA))
  }
  
  
  else if  # se i dati sono costanti negli anni non plotto la curva
  ( sum( diff(y))==0){
    plot(x,y,ylim=c(0,1),xlim=c(2,6),xlab="",ylab="",xaxt="n",yaxt="n")
    text(2.4,0.8,paste(year))
    return(c(NA,NA,NA))
  }
  
   #if(a < 1){
   #	bb<-c(2,4,6,8)
   	#return(bb)}
  
  else 
  {y<-y/max(y)
   ####B40F20, "#E2D200", "#46ACC8", "#E58601", "#B40F20"
   try(r <- drm(y ~ x,fct=baro5(fixed=c(NA,NA,NA,1,NA))))
   #plot(r, type="confidence")
   plot(x,y,type="n",ylim=c(0,1),xlim=c(2,6),xlab="",ylab="",xaxt="n",yaxt="n") #scatola
   abline(h=0,col="gray");abline(h=1,col="gray") # barre superiori ed inferiori
   plot(r,pch=20,lty=1,add=T,lwd=4, col="#B40F20") #curva
   
   xx<-NA;#length(pr)<-10000
   xx<-seq(min(x),max(x),length.out=10000)
   pr<-NA;length(pr)<-10000# preallocazione memoria
   pr<-predict(r,newdata=data.frame(x=xx))
   
   dpr<-diff(pr,1)/diff(xx)
      
   if(which.max(dpr)==1){
   	return(c(NA, NA, NA))
   }

   infl<-xx[which.max(dpr)-1]   
   
   sl<-max(dpr)
   int<- predict(r,newdata=data.frame(x=infl))-(sl * infl)
   bio<- predict(r, newdata=data.frame(x=infl)) ###RIGA AGGIUNTA###
   
   
   abline(v=infl,col="#E58601",lty=3,lwd=3)
   abline(a=int,b=sl,col="#B40F20",lty=3,lwd=3)
   abline(h=bio,col="#46ACC8",lty=3,lwd=3)
   abline(h=coef(r)[3],col="#E2D200",lty=3,lwd=3)
   
   
   text(2.4,0.8,paste(filelist[j]), cex=0.75)
   return(c(coef(r)[3],sl,infl,bio))  }
  
}

