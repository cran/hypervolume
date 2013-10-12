
plot.Hypervolume <- function(x, ...)
{
  templist = new("HypervolumeList")
  templist@HVList=list(x)
  plot.HypervolumeList(templist, ...)
}

plot.HypervolumeList <- function(x, npmax = 1000, colors=rainbow(length(x@HVList)), names=NULL, reshuffle=TRUE, showdensity=TRUE,showdata=TRUE,darkfactor=0.5,cex.random=0.5,cex.data=0.75,cex.axis=0.75,cex.names=1.0,cex.legend=0.75,legend=TRUE, varlims=NULL,pairplot=TRUE,whichaxes=NULL,...)
{
  alldims = sapply(x@HVList, function(z) { z@Dimensionality })
  allnames = sapply(x@HVList, function(z) { z@Name })
  stopifnot(all(alldims[1] == alldims))
  
  all <- NULL
  alldata <- NULL
  for (i in 1:length(x@HVList))
  {
    ivals = sample(nrow(x@HVList[[i]]@RandomUniformPointsThresholded), npmax, replace=T)
    subsampledpoints = data.frame(x@HVList[[i]]@RandomUniformPointsThresholded[ivals,])
    densityvals = x@HVList[[i]]@ProbabilityDensityAtRandomUniformPoints[ivals]
    
    if (nrow(subsampledpoints) > 0)
    {  
      subsampledpoints = cbind(subsampledpoints, ID=rep(i, nrow(subsampledpoints)), Density=densityvals/max(densityvals,na.rm=T))
    
      all <- rbind(all, subsampledpoints)
    }
    
    thisdata=x@HVList[[i]]@Data
    alldata <- rbind(alldata, cbind(thisdata, ID=rep(i,nrow(thisdata))))
  }  
  all <- unique(all)
  alldata <- as.data.frame(alldata)
  
  if (reshuffle==TRUE)
  {
    all <- all[sample(nrow(all)),] # reorder to shuffle colors
  }
  
  if (is.null(names))
  {
    names = names(all)[1:(ncol(all)-2)]
  }  
  
  if (!is.null(varlims) & !is.list(varlims))
  {
    varlimlist = vector('list',ncol(all)-2)
    for (i in 1:length(varlimlist))
    {
      varlimlist[[i]] <- varlims
    }
    varlims = varlimlist
  }
  
  colorlist <- colors[all$ID]
  alphavals <- all$Density
  if (showdensity)
  {
    colorlist <- rgb2rgba(colorlist, alphavals)
  }
  
  colorlistdata = colors[alldata$ID]
  colorlistdata <- rgb2rgbdark(colorlistdata, darkfactor)
  
  
  
  
  if (pairplot)
  {
    op = par(no.readonly = T)
    
    par(mfrow=c(ncol(all)-2, ncol(all)-2))
    par(mar=c(0,0,0,0))
    
    for (i in 1:(ncol(all)-2))
    {
      for (j in 1:(ncol(all)-2))  
      {
        if (j > i)
        {
          plot(all[,j], all[,i],type="n",axes=F,xlim=varlims[[j]], ylim=varlims[[i]])
          
          
          points(all[,j], all[,i], col=colorlist,cex=cex.random,pch=16)
          
          if (showdata)
          {
            points(alldata[,j], alldata[,i], col=colorlistdata,cex=cex.data,pch=16)
          }
          
          box()
        }
        else if (j == i)
        {
          plot(0,0,type="n",xlim=c(0,1),ylim=c(0,1),axes=F)
          text(0.5, 0.5, names[j],cex=cex.names)
        }
        else if (j==1 & i == (ncol(all) - 2))
        {
          plot(0,0,type="n",xlim=c(0,1),ylim=c(0,1),axes=F)
          
          if (legend == TRUE)
          {
            legend('topleft',legend=allnames,text.col=colors,bty='n',cex=cex.legend)
          }
        }
        else
        {
          plot(0,0,type="n",axes=F)    
        }
        
        if (j==i+1)
        {
          axis(side=1,cex.axis=cex.axis)
          axis(side=2,cex.axis=cex.axis)
        }
      }
    }  
    par(op)
  }
  else
  {
    if (is.null(whichaxes))
    {
      whichaxes=1:3  
    }
    if (is.null(names))
    {
      names <- names(data)
    }
    if(length(whichaxes)!=3) { stop('Must specify three axes') }
    
    if (all(is.numeric(whichaxes)))
    {
      axesnames <- names(data)[whichaxes]
    }
    
    plot3d(all[,whichaxes],col=colorlist,xlab=axesnames[1], ylab=axesnames[2], zlab=axesnames[3], xlim=varlims[[1]],ylim=varlims[[2]],zlim=varlims[[3]],size=cex.random,type='p',expand=1.05)
    
    if (legend==TRUE)
    {
      for (i in 1:length(allnames))
      {
        mtext3d(allnames[i],edge='x-+',line=1+i*cex.legend*1.25,color=colors[i],cex=cex.legend)  
      }
    }
    
    if (showdata)
    {
      if (!any(is.nan(as.matrix(alldata[,whichaxes]))))
      {
        points3d(x=alldata[,whichaxes[1]], y=alldata[,whichaxes[2]], z=alldata[,whichaxes[3]], col=colorlistdata,cex=cex.data,pch=16)
      }
    }
  }
}  
