# Adaption of the original function that can be found here: https://www.rdocumentation.org/packages/plotrix/versions/3.7-4/topics/axis.break

#' places a break marker at the position "breakpos" in user coordinates on the axis nominated
#'
#' @param axis numeric: which axis the break will be on
#' @param breakpos numeric: where the break will be
#' @param bgcol character string: the background color; white by default
#' @param breakcol character string: the color of the axis break; black by default
#' @param style list of strings: whether the break is a slash or a zigzag
#' @param brw numeric: the break width; 0.02 by default
#' @param Lwd numeric: line width for the plot; 1 by default
#'
#' @return The current plot with zigzag axis breaks at the desired place, to better fit data in a visualization
#' @import data.table
#' @import plyr
#' @import dplyr
#' @import tidyr
#' @import jsonlite
#' @import mongolite
#' @import stringr
#' @import tictoc
#' @import tibble
#'
#'
#' @export axis_break
#'
axis_break<-function(axis=1,breakpos,bgcol="white",breakcol="black",
                     style=c("slash"),brw=0.02, Lwd=1) {

  if(!missing(breakpos)) {
    # get the coordinates of the outside of the plot
    figxy<-par("usr")
    # flag if either axis is logarithmic
    xaxl<-par("xlog")
    yaxl<-par("ylog")
    # calculate the x and y offsets for the break
    xw<-(figxy[2]-figxy[1])*brw
    yw<-(figxy[4]-figxy[3])*brw
    if(xaxl && (axis == 1 || axis == 3)) breakpos<-log10(breakpos)
    if(yaxl && (axis == 2 || axis == 4)) breakpos<-log10(breakpos)
    # set up the "blank" rectangle (left, bottom, right, top)
    switch(axis,
           br<-c(breakpos-xw/2,figxy[3]-yw/2,breakpos+xw/2,figxy[3]+yw/2),
           br<-c(figxy[1]-xw/2,breakpos-yw/2,figxy[1]+xw/2,breakpos+yw/2),
           br<-c(breakpos-xw/2,figxy[4]-yw/2,breakpos+xw/2,figxy[4]+yw/2),
           br<-c(figxy[2]-xw/2,breakpos-yw/2,figxy[2]+xw/2,breakpos+yw/2),
           stop("Improper axis specification."))
    # get the current setting of xpd
    old.xpd<-par("xpd")
    # don't cut the break off at the edge of the plot
    par(xpd=T)
    # correct for logarithmic axes
    if(xaxl) br[c(1,3)]<-10^br[c(1,3)]
    if(yaxl) br[c(2,4)]<-10^br[c(2,4)]
    # draw the "blank" rectangle
    rect(br[1],br[2],br[3],br[4],col=bgcol,border=bgcol)
    if(style == "slash") {
      # calculate the slash ends
      if(axis == 1 || axis == 3) {
        xbegin<-c(breakpos-xw,breakpos)
        xend<-c(breakpos,breakpos+xw)
        ybegin<-c(br[2],br[2])
        yend<-c(br[4],br[4])
        if(xaxl) {
          xbegin<-10^xbegin
          xend<-10^xend
        }
      }
      else {
        xbegin<-c(br[1],br[1])
        xend<-c(br[3],br[3])
        ybegin<-c(breakpos-yw,breakpos)
        yend<-c(breakpos,breakpos+yw)
        if(yaxl) {
          ybegin<-10^ybegin
          yend<-10^yend
        }
      }
    }
    else {
      # calculate the zigzag ends
      if(axis == 1 || axis == 3) {
        xbegin<-c(breakpos-xw/2,breakpos-xw/4,breakpos+xw/4)
        xend<-c(breakpos-xw/4,breakpos+xw/4,breakpos+xw/2)
        ybegin<-c(ifelse(yaxl,10^figxy[3+(axis==3)],figxy[3+(axis==3)]),br[4],br[2])
        yend<-c(br[4],br[2],ifelse(yaxl,10^figxy[3+(axis==3)],figxy[3+(axis==3)]))
        if(xaxl) {
          xbegin<-10^xbegin
          xend<-10^xend
        }
      }
      else {
        xbegin<-c(ifelse(xaxl,10^figxy[1+(axis==4)],figxy[1+(axis==4)]),br[1],br[3])
        xend<-c(br[1],br[3],ifelse(xaxl,10^figxy[1+(axis==4)],figxy[1+(axis==4)]))
        ybegin<-c(breakpos-yw/2,breakpos-yw/4,breakpos+yw/4)
        yend<-c(breakpos-yw/4,breakpos+yw/4,breakpos+yw/2)
        if(yaxl) {
          ybegin<-10^ybegin
          yend<-10^yend
        }
      }
    }
    # draw the segments
    segments(xbegin,ybegin,xend,yend,col=breakcol,lty=1, lwd=Lwd)
    # restore xpd
    par(xpd=old.xpd)
  }
  else {
    cat("Usage: axis.break(axis=1,breakpos,bgcol=\"white\",breakcol=\"black\",\n")
    cat("\tstyle=c(\"slash\",\"zigzag\"),brw=0.02)\n")
  }
}

#' Replaces colons in chemical names with '-' so they will be read in correctly
#'
#' @param String character string: chem names to correct
#'
#' @return chem_names that fit the correct pattern
#'
#' @import data.table
#' @import plyr
#' @import dplyr
#' @import tidyr
#' @import jsonlite
#' @import mongolite
#' @import stringr
#' @import tictoc
#' @import tibble
#'
#' @export correctChemName
#'
#' @examples
#' correctChemName("test:chemical")
#'
correctChemName<-function(String){ #Converts chemical names to the proper format by replacing delimiters
  str_replace_all(String, pattern=":", replacement="-")
}

#' finds columns by name, if inv is true, only finds those columns that exist in the data
#'
#' @param Table A table of cell-painting data
#' @param inv A boolean, if it's TRUE, findDataCols only finds those columns that exist in the data
#' @param names A boolean, if it's TRUE, return only the column names, if it's FALSe, return the columns
#'
#' @return Either the column names, or the columns in Table, based on whether names is TRUE
#'
#' @import data.table
#' @import plyr
#' @import dplyr
#' @import tidyr
#' @import jsonlite
#' @import mongolite
#' @import stringr
#' @import tictoc
#' @import tibble
#'
#' @export findDataCols
#'
findDataCols<-function(Table, inv=F, names=F){
  iDataCols= c( grep("^Position", colnames(Table)),
                grep("^Shape", colnames(Table)),
                grep("^DNA_", colnames(Table)),
                grep("^RNA_", colnames(Table)),
                grep("^ER_", colnames(Table)),
                grep("^AGP_", colnames(Table)),
                grep("^Mito_", colnames(Table))  )
  if(inv){
    iDataCols=setdiff(1:dim(Table)[2], iDataCols)
  }
  if(names){
    return( colnames(Table)[iDataCols] )
  }else{
    return(iDataCols)
  }
}


#' Assigns a color to each column of data by name
#'
#' @param Data dataframe: the data being assigned color
#' @param dataform character string: whether the data is a "matrix" or a "vector"
#'
#' @return array of colors matching each of the columns in the data, to distinguish them in plots
#'
#' @import data.table
#' @import plyr
#' @import dplyr
#' @import tidyr
#' @import jsonlite
#' @import mongolite
#' @import stringr
#' @import tictoc
#' @import tibble
#'
#' @export ColumnColors
#'

ColumnColors<-function(Data, dataform="matrix"){ #
  if(dataform=="matrix"){
    ColumnColors=rep("gray80", times=dim(Data)[2])
    ColumnColors[grep("DNA", colnames(Data))]="blue"
    ColumnColors[grep("RNA", colnames(Data))]="green"
    ColumnColors[grep("^ER", colnames(Data))]="green4"
    ColumnColors[grep("AGP", colnames(Data))]="orange"
    ColumnColors[grep("Mito", colnames(Data))]="deeppink3"
  }else if (dataform=="vector"){
    ColumnColors=rep("gray80", times=length(Data))
    ColumnColors[grep("DNA", Data)]="blue"
    ColumnColors[grep("RNA", Data)]="green"
    ColumnColors[grep("^ER", Data)]="green4"
    ColumnColors[grep("AGP", Data)]="orange"
    ColumnColors[grep("Mito", Data)]="deeppink3"
  }else{
    stop("data must be in 'matrix' or 'vector' form")
  }

  return(ColumnColors)
}


#' Orders the columns to be a heatmap axis
#'
#' @param Data table to be reordered
#'
#' @import data.table
#' @import plyr
#' @import dplyr
#' @import tidyr
#' @import jsonlite
#' @import mongolite
#' @import stringr
#' @import tictoc
#' @import tibble
#'
#' @return newOrder
#'
#' @export Cluster.by.Channel
Cluster.by.Channel<-function(Data){ #
  newOrder<-NULL
  for(Parameter in c("Position",  "Shape", "DNA", "RNA", "^ER", "AGP", "Mito")){
    i=grep(Parameter, colnames(Data))
    if(length(i)>1){
      png("temp.png",  width=2000, height=1000)
      test=heatmap.2(Data[,i], scale="none", dendrogram="column", plot=F)
      dev.off()
      newOrder=c(newOrder, i[test$colInd])
      rm(i, test)
    }else{
      newOrder=c(newOrder, i)
      rm(i)
    }
  }
  file.remove("temp.png")
  return(newOrder)
}

#' gives one letter abbreviations to attributes
#'
#' @param Vector a vector of cell attributes
#'
#' @import data.table
#' @import plyr
#' @import dplyr
#' @import tidyr
#' @import jsonlite
#' @import mongolite
#' @import stringr
#' @import tictoc
#' @import tibble
#'
#' @export GOLetter
#'
#' @return vector of single-letter abbreviations of the attributes listed, to distinguish points on a plot
#'
#' @examples
#' GOLetter(c("test"))
#'
GOLetter<-function(Vector){ #
  Letter=rep("M", times=length(Vector) )
  Letter[grep("_Texture", Vector)]="T"
  Letter[grep("_Intensity", Vector)]="I"
  Letter[grep("Symmetry", Vector)]="S"
  Letter[grep("Compactness", Vector)]="C"
  Letter[grep("Axial", Vector)]="A"
  Letter[grep("Radial", Vector)]="R"
  Letter[grep("Profile", Vector)]="P"
  return(Letter)
}

#' PchShape
#'
#' @param Parameter.Name character string: the parameter names being assigned graph shapes
#'
#' @import data.table
#' @import plyr
#' @import dplyr
#' @import tidyr
#' @import jsonlite
#' @import mongolite
#' @import stringr
#' @import tictoc
#' @import tibble
#'
#' @export PchShape
#'
#' @return vector of 2-digit base r plots codes for single-letter abbreviations of the attributes listed, to distinguish points on a plot
#'
#' @examples
#' PchShape("AGP")
#'
PchShape<-function(Parameter.Name){ #assigns shapes to parameters for distinguishing them in plots
  PchShape=rep(19, times=length(Parameter.Name) ) #19=round
  PchShape[grep("_Nuclei",      Parameter.Name)]= 15 #square
  PchShape[grep("_Ring",        Parameter.Name)]= 18 #diamond
  PchShape[grep("_Membrane",    Parameter.Name)]= 17 #traingle
  PchShape[grep("Shape_Nuclei", Parameter.Name)]=15 #square #
  return(PchShape)
}

# PchLetter<-function(Vector){ #same as GOLetter but stores the letters as two-digit codes
#   PchLetter=rep(77, times=length(Vector) ) #M
#   PchLetter[grep("_Texture", Vector)]=84 #T
#   PchLetter[grep("_Intensity", Vector)]=73 #I
#   PchLetter[grep("_Symm", Vector)]=83 #S
#   PchLetter[grep("_Compact", Vector)]=67 #C
#   PchLetter[grep("_Axial", Vector)]=65 #A
#   PchLetter[grep("_Radial", Vector)]=82 #R
#   PchLetter[grep("_Profil", Vector)]=80 #P
#   return(PchLetter)
# }


#' Plots concentration-response data to visualize Benchmark Dose (BMD)
#'
#' @param Subset dataframe: the subset of the data to be plotted, for instance subset by a given chemical
#' @param xLim numeric: x limit of the plot; null by default.
#' @param TestedRange the dose range tested; null by default.
#' @param plotDatapoints boolean: whether the individual datapoints should be visualized on the plot
#' @param plotBMC boolean: whether the BMC should be highlighted on the plot; false by default
#' @param plotDoserange boolean: whether the doserange should be plotted; false by default
#' @param plotNoiseband  boolean: whether the noise band should be plotted; true by default
#' @param Lwd numeric: line width on the plot; 1 by default
#' @param cexAxis numeric: scale of the axis labels; 1 by default
#' @param yLim numeric: y limit of the plot; c(-5.5, 100) by default
#' @param yTicks numeric: the placement of y axis ticks; null by default.
#' @param yAxisSteps the frequency of y axis ticks, if yTicks is undefined
#'
#' @return hill plot of the data
#'
#' @import data.table
#' @import plyr
#' @import dplyr
#' @import tidyr
#' @import jsonlite
#' @import mongolite
#' @import stringr
#' @import tictoc
#' @import tibble
#'
#' @export plotCurves
#'
plotCurves <- function(Subset, xLim=NULL, TestedRange=NULL,
                       plotDatapoints=F, plotBMC=F, plotDoserange=F, plotNoiseband=T,
                       Lwd=1,
                       cexAxis=1,
                       yLim = c(-5.5, 100), yTicks=NULL, yAxisSteps=25){

  if(length(setdiff("a", colnames(Subset)))){ Subset = Subset %>% dplyr::mutate(a = NA)}
  if(length(setdiff("b", colnames(Subset)))){ Subset = Subset %>% dplyr::mutate(b = NA)}
  if(length(setdiff("p", colnames(Subset)))){ Subset =Subset %>%  dplyr::mutate(p = NA)}
  if(length(setdiff("q", colnames(Subset)))){ Subset = Subset %>% dplyr::mutate(q = NA)}
  if(length(setdiff("tp", colnames(Subset)))){ Subset = Subset %>% dplyr::mutate(tp = NA)}
  if(length(setdiff("ga", colnames(Subset)))){ Subset = Subset %>% dplyr::mutate(ga = NA)}
  if(length(setdiff("la", colnames(Subset)))){ Subset = Subset %>%  dplyr::mutate(la = NA)}

  #make concentrations log10
  Subset = Subset %>% dplyr::mutate_at(.vars=c("min_conc", "max_conc"), .funs=c("log10"))


  #define the plotting range
  if(is.null(xLim)){
    minlogX=floor  (min(Subset$min_conc))
    maxlogX=ceiling(max(Subset$max_conc))
    if(is.na(maxlogX)){maxlogX = minlogX+1}

    xLim = 10^c(minlogX, maxlogX) #in regular (non-log) units
  }else{
    minlogX = floor(log10(min(xLim)))
    maxlogX = floor(log10(max(xLim)))
  }

  if(is.null(yTicks)){
    yTicks = seq(0, max(yLim), by=yAxisSteps)
  }


  #empty plot to start with
  plot(c(0,0),type="n",xlab="",ylab="", xlim=xLim, ylim=yLim,
       log="x",
       xaxt="n", yaxt="n", bty="n",
       font.lab=2, las=1)

  if(plotDoserange){ #sets the plot range based on the data if no range is provided
    if(is.null(TestedRange)){
      rect(xleft=10^min(Subset$min_conc), xright=10^max(Subset$max_conc), ybottom=0, ytop=0.05*diff(yLim), col="gray80", border="transparent")
    }else{
      rect(xleft=TestedRange$min_conc, xright=TestedRange$max_conc, ybottom=0, ytop=0.05*diff(yLim), col="gray80", border="transparent")
    }
  }

  #Axis Ticks
  Ticks=NULL
  for (i in seq(minlogX, maxlogX-1)) {
    Ticks=c(Ticks, 10^i*c(1,2,3,4,5,6,7,8,9))
  }
  Ticks=c(Ticks, 10^maxlogX)
  Ticks=Ticks[Ticks>=xLim[1]]
  Ticks=Ticks[Ticks<=xLim[2]]

  TickLabelPos=10^seq(minlogX, maxlogX)

  ##plot x-axis ticks & labels
  axis(1, at=xLim, labels=F, lwd=Lwd, lwd.ticks=0, col.axis="black", cex.axis=cexAxis) #makes axis only if 2 values are entered
  axis(1, at=Ticks, labels=FALSE, lwd=Lwd, col.axis="black")
  axis(1, at=TickLabelPos, labels=FALSE, lwd=Lwd*1.5, col.axis="black")
  axis(1, at=c(TickLabelPos),               lwd=0,   labels=c(TickLabelPos),line=-0.5,  col.axis="black",  cex.axis=cexAxis)

  ##plot y-axis ticks & labels
  axis(2, at=c(yLim[1],yLim[2]*1.05),  labels=F, lwd=Lwd, lwd.ticks=0, col.axis="black", cex.axis=cexAxis)
  axis(2, at=yTicks, labels=F,  col.axis="black", las=1, lwd=Lwd, cex.axis=cexAxis)
  axis(2, at=yTicks, labels=TRUE, line=-0.25, col.axis="black", las=1, lwd=0, cex.axis=cexAxis)

  for(i in 1:dim(Subset)[1]){
    row = Subset[i,]
    print(row$endpoint)

    if(plotNoiseband){
      rect(xleft=10^min(row$min_conc), xright=10^max(row$max_conc), ybottom=row$ctr_median-row$ctr_nmad, ytop=row$ctr_median+row$ctr_nmad,
           col=adjustcolor(row$Color, alpha.f=0.2), border="transparent")

    }

    #reformat conc and resp as vectors
    ZF = tibble(conc =  as.numeric(str_split(row$conc,"\\|")[[1]]),
                resp = as.numeric(str_split(row$resp,"\\|")[[1]]) ) %>%
      filter(!is.na(resp)) %>%
      group_by(conc) %>%
      dplyr::summarise(Mean = mean(resp),
                       SD = sd(resp))

    conc_plot <- 10**seq(from=min(row$min_conc),to=max(row$max_conc),by=0.05)

    #get model parameters
    parnames = c("a", "tp", "b", "ga", "p", "la", "q")
    modpars = as.list(row[,parnames])
    modpars= modpars[!sapply(modpars, is.na)]

    fit_method = row$fit_method

    #calculate and plot model curves
    if(fit_method == "hill"){
      resp_plot <- do.call("hillfn",list(ps = unlist(modpars), x = conc_plot))
    } else if(!fit_method %in% c("cnst","none") ){
      resp_plot <- do.call(fit_method,list(ps = unlist(modpars), x = conc_plot))
    } else{
      resp_plot = rep(0, length.out=length(conc_plot))
    }

    resp_plot = resp_plot+row$ctr_median

    lines(resp_plot~conc_plot,col=row$Color, lwd=row$Lwd, lty=row$Lty)

    #plot data points and arrows
    if(plotDatapoints){
      ZF$Mean = ZF$Mean + row$ctr_median
      arrows(ZF$conc, ZF$Mean-ZF$SD,  ZF$conc, ZF$Mean+ZF$SD,
             code=3, angle=90, length=0.04, lwd=0.75*row$Lwd, col=row$Color)

      points(ZF$conc, ZF$Mean, col=row$Color, pch=row$Pch, cex=row$cexPch)
    }

    #plot a point for the BMC
    if(plotBMC){
      if(!is.na(row$BMC)){
        if(row$Assay=="cell count"){
          segments(x0 = row$BMC, x1=row$BMC, y0 = -3, y1=50, col=row$Color, lty="dashed")
          points(x=row$BMC, y=50,            col=row$Color,  pch=row$Pch, cex=row$PchCex)
        }else if(grepl(pattern="global", row$Assay)){
          segments(x0 = row$BMC, x1=row$BMC, y0 = -3, y1=row$ctr_median+row$bmr, col=row$Color, lty="dashed")
          points(x=row$BMC, y=row$ctr_median+row$bmr, col=row$Color,  pch=row$Pch, cex=row$PchCex)
        }else{
          segments(x0 = row$BMC, x1=row$BMC, y0 = -3, y1=row$ctr_median+3*row$ctr_nmad, col=row$Color, lty="dashed")
          points(x=row$BMC, y=row$ctr_median+3*row$ctr_nmad, col=row$Color,  pch=row$Pch, cex=row$PchCex)
        }
      }
    }#if plotBMC

  }#for each row

}#end of function

#' Concentration Response Plot
#'
#' @param row vector: Row containing response data to be plotted
#' @param ymin numeric: y axis minimum
#' @param ymax numeric: y axis maximum
#'
#' @import data.table
#' @import plyr
#' @import dplyr
#' @import tidyr
#' @import jsonlite
#' @import mongolite
#' @import stringr
#' @import tictoc
#' @import tibble
#'
#' @export concRespPlot_JN
#'
#' @return a log concentration-response plot
#'
#' @examples
#' conc <- list(.03, .1, .3, 1, 3, 10, 30, 100)
#' resp <- list(0, .2, .1, .4, .7, .9, .6, 1.2)
#' row <- list(conc = conc,
#'             resp = resp,
#'             bmed = 0,
#'             cutoff = 0.25,
#'             onesd = 0.125,
#'             name = "some chemical",
#'             assay = "some assay")
#' res <- tcplfit2::concRespCore(row, conthits = TRUE)
#' res <- dplyr::mutate(.data = res,
#'                      la = NA,
#'                      q = NA,
#'                      stype = "test sample",
#'                      endpoint = "test")
#' Subset <- dplyr::mutate(.data = res,
#'                         acc = ifelse(is.na(acc), bmd, acc))
#' Subset2 <- dplyr::filter(.data = Subset,
#'                          stype %in% c("test sample") & bmd < 100 & hitcall > 0.90)
#' concRespPlot_JN(Subset2,ymin=-0.5,ymax=2)
#'
concRespPlot_JN <- function(row, ymin=NULL, ymax=NULL) {

  dtxsid <- row[1,"dtxsid"]

  #every variable in PATHWAY_CR goes into the environment to make it easy
  #to update this function to use new PATHWAY_CR data.
  list2env(row,envir = environment())

  #reformat conc and resp as vectors
  conc <- as.numeric(str_split(row[1,"conc"],"\\|")[[1]])
  resp <- as.numeric(str_split(row[1,"resp"],"\\|")[[1]])

  #plotting points for curves (Johanna)
  plotrange_log = c(floor(min(log10(conc))),ceiling(max(log10(conc))))
  plotrange = 10^plotrange_log
  logc_plot <- seq(from=plotrange_log[1],to=plotrange_log[2],by=0.05)
  conc_plot <- 10**logc_plot

  #y range adapted (Johanna)
  if(is.null(ymax)){
    ymax <- max(1, cutoff, max(resp)*1.1)
  }
  if(is.null(ymin)){
    ymin <- min(-1, -cutoff, min(resp)*1.1)
  }

  col.list <- c("black","blue","red")

  #empty plot to start with
  plot(c(1,1),type="n",xlab="conc (uM)",ylab="Score",xlim=plotrange,ylim=c(ymin,ymax),
       log="x",main=paste(name,"\n",endpoint),cex.main=0.9)

  #cutoffs and gray rectangular noise region
  xmin = min(plotrange)
  xmax = max(plotrange)
  rect(xleft=xmin,ybottom=-bmr/1.349,xright=xmax,ytop=bmr/1.349,col="gray80", border=NULL)
  lines(c(xmin,xmax),c(cutoff,cutoff),lwd=1, col="gray50", lty="dashed")
  lines(c(xmin,xmax),c(-cutoff,-cutoff),lwd=1, col="gray50", lty="dashed")

  #thick line at 0 and bmrs
  lines(c(xmin,xmax),c(0,0),lwd=1, lty="dotted")

  #height for top labels
  yplot <- ymax*0.95
  #equally spaced labelling points pegged to plotrange
  xplot = 10^(seq(log10(plotrange[1]), log10(plotrange[2]), length.out = 8))[-8]

  #Top label headings
  text(xplot[1],yplot,"time",pos=4)
  text(xplot[2],yplot,"mthd",pos=4)
  text(xplot[3],yplot,"top_over_co",pos=4)
  text(xplot[4],yplot,"Top",pos=4)
  text(xplot[5],yplot,"BMD",pos=4)
  text(xplot[6],yplot,"ACC",pos=4)
  text(xplot[7],yplot,"Hitcall",pos=4)

  j = 1 #j = 1 is black
  conc <- conc[!is.na(resp)]
  resp <- resp[!is.na(resp)]
  points(resp~conc,pch=19,col=adjustcolor(col.list[j], alpha.f=0.5)) #plot actual data

  #get model parameters
  parnames = c("a", "tp", "b", "ga", "p", "la", "q")
  modpars = as.list(row[,parnames])
  modpars= modpars[!sapply(modpars, is.na)]

  #gcalculate and plot model curves
  if(fit_method == "hill"){
    resp_plot <- do.call("hillfn",list(ps = unlist(modpars), x = conc_plot))
    lines(resp_plot~conc_plot,col=col.list[2], lwd=2)
  } else if(!fit_method %in% c("cnst","none") ){
    resp_plot <- do.call(fit_method,list(ps = unlist(modpars), x = conc_plot))
    lines(resp_plot~conc_plot,col=col.list[2], lwd=2)
  }

  yplot <- yplot-(ymax-ymin)*0.05 #second row for top labels

  #Fill in top labels second row
  text(xplot[2],yplot,fit_method,pos=4)
  text(xplot[3],yplot,format(top_over_cutoff,digits=3),pos=4, col = "red")
  text(xplot[4],yplot,format(top,digits=3),pos=4)
  text(xplot[5],yplot,format(bmd,digits=2),pos=4, col = "green")
  text(xplot[6],yplot,format(acc,digits=2),pos=4, col = "blue")

  #Bottom left info
  text(xplot[1],0.85*ymin,paste("method:",fit_method,"\nCutoff=",format(cutoff,digits=2),sep=""),pos=4)

  #color hitcall based on whether it's a hit
  color <- "black"
  font <- 1
  if(hitcall==1) {
    color <- "red"
    font <- 2
  }
  text(xplot[7],yplot,format(hitcall, digits = 2),pos=4,col=color,cex=1,font=font)

  #plot green bmd with range
  if(hitcall>0.1 & acc<max(conc) & !is.na(bmd)) {
    lines(c(bmd,bmd),c(ymin/2,ymax/2),col="green",lwd=2, lty = isTRUE(bmd<min(conc)) + 1)
    if(is.na(bmdl)) xleft = plotrange[1]/10 else xleft = bmdl
    if(is.na(bmdu)) xright = plotrange[2]*10 else xright = bmdu

    rect(xleft=xleft,ybottom=ymin/2,xright=xright,ytop=ymax/2,col=rgb(0,1,0, alpha = .5), border = NA)
    lines(c(xleft,xleft),c(ymin/2,ymax/2),col="green",lwd=1)
    lines(c(xright,xright),c(ymin/2,ymax/2),col="green",lwd=1)
  }
}


# #' Plots heatmaps with finer control over color, breaks and level trace
# #'
# #' @param x a table to be converted to a matrix
# #' @param Rowv a boolean
# #' @param Colv
# #' @param distfun
# #' @param hclustfun
# #' @param dendrogram
# #' @param symm
# #' @param scale
# #' @param na.rm
# #' @param revC
# #' @param add.expr
# #' @param breaks
# #' @param symbreaks
# #' @param col
# #' @param colsep
# #' @param rowsep
# #' @param sepcolor
# #' @param sepwidth
# #' @param cellnote
# #' @param notecex
# #' @param notecol
# #' @param na.color
# #' @param trace
# #' @param tracecol
# #' @param hline
# #' @param vline
# #' @param linecol
# #' @param margins
# #' @param ColSideColors
# #' @param RowSideColors
# #' @param side.height.fraction
# #' @param cexRow
# #' @param cexCol
# #' @param labRow
# #' @param labCol
# #' @param key
# #' @param keysize
# #' @param density.info
# #' @param denscol
# #' @param symkey
# #' @param densadj
# #' @param main
# #' @param xlab
# #' @param ylab
# #' @param lmat
# #' @param lhei
# #' @param lwid
# #' @param ColSideColorsSize
# #' @param RowSideColorsSize
# #' @param KeyValueName
# #' @param ...
# #'
# #' @return heatmap of the input data
# #' @export
# #'
# #' @import data.table
# #' @import plyr
# #' @import dplyr
# #' @import tidyr
# #' @import jsonlite
# #' @import mongolite
# #' @import stringr
# #' @import tictoc
# #' @import tibble
# #'
# #' @examples
# heatmap.2plus <- function(x, #
#                           Rowv = TRUE,
#                           Colv = if (symm) "Rowv" else TRUE,
#                           distfun = dist,
#                           hclustfun = hclust,
#                           dendrogram = c("both","row", "column", "none"),
#                           symm = FALSE,
#                           scale = c("none","row", "column"),
#                           na.rm = TRUE,
#                           revC = identical(Colv,"Rowv"),
#                           add.expr,
#                           breaks,
#                           symbreaks = max(x < 0, na.rm = TRUE) || scale != "none",
#                           col = "heat.colors",
#                           colsep,
#                           rowsep,
#                           sepcolor = "white",
#                           sepwidth = c(0.05, 0.05),
#                           cellnote,
#                           notecex = 1,
#                           notecol = "cyan",
#                           na.color = par("bg"),
#                           trace = c("none", "column","row", "both"),
#                           tracecol = "cyan",
#                           hline = median(breaks),
#                           vline = median(breaks),
#                           linecol = tracecol,
#                           margins = c(5,5),
#                           ColSideColors,
#                           RowSideColors,
#                           side.height.fraction=0.3,
#                           cexRow = 0.2 + 1/log10(nr),
#                           cexCol = 0.2 + 1/log10(nc),
#                           labRow = NULL,
#                           labCol = NULL,
#                           key = TRUE,
#                           keysize = 1.5,
#                           density.info = c("none", "histogram", "density"),
#                           denscol = tracecol,
#                           symkey = max(x < 0, na.rm = TRUE) || symbreaks,
#                           densadj = 0.25,
#                           main = NULL,
#                           xlab = NULL,
#                           ylab = NULL,
#                           lmat = NULL,
#                           lhei = NULL,
#                           lwid = NULL,
#                           ColSideColorsSize = 1,
#                           RowSideColorsSize = 1,
#                           KeyValueName="Value",...){
#
#   invalid <- function (x) {
#     if (missing(x) || is.null(x) || length(x) == 0)
#       return(TRUE)
#     if (is.list(x))
#       return(all(sapply(x, invalid)))
#     else if (is.vector(x))
#       return(all(is.na(x)))
#     else return(FALSE)
#   }
#
#   x <- as.matrix(x)
#   scale01 <- function(x, low = min(x), high = max(x)) {
#     x <- (x - low)/(high - low)
#     x
#   }
#   retval <- list()
#   scale <- if (symm && missing(scale))
#     "none"
#   else match.arg(scale)
#   dendrogram <- match.arg(dendrogram)
#   trace <- match.arg(trace)
#   density.info <- match.arg(density.info)
#   if (length(col) == 1 && is.character(col))
#     col <- get(col, mode = "function")
#   if (!missing(breaks) && (scale != "none"))
#     warning("Using scale=\"row\" or scale=\"column\" when breaks are",
#             "specified can produce unpredictable results.", "Please consider using only one or the other.")
#   if (is.null(Rowv) || is.na(Rowv))
#     Rowv <- FALSE
#   if (is.null(Colv) || is.na(Colv))
#     Colv <- FALSE
#   else if (Colv == "Rowv" && !isTRUE(Rowv))
#     Colv <- FALSE
#   if (length(di <- dim(x)) != 2 || !is.numeric(x))
#     stop("`x' must be a numeric matrix")
#   nr <- di[1]
#   nc <- di[2]
#   if (nr <= 1 || nc <= 1)
#     stop("`x' must have at least 2 rows and 2 columns")
#   if (!is.numeric(margins) || length(margins) != 2)
#     stop("`margins' must be a numeric vector of length 2")
#   if (missing(cellnote))
#     cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
#   if (!inherits(Rowv, "dendrogram")) {
#     if (((!isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram %in%
#                                                  c("both", "row"))) {
#       if (is.logical(Colv) && (Colv))
#         dendrogram <- "column"
#       else dedrogram <- "none"
#       warning("Discrepancy: Rowv is FALSE, while dendrogram is `",
#               dendrogram, "'. Omitting row dendogram.")
#     }
#   }
#   if (!inherits(Colv, "dendrogram")) {
#     if (((!isTRUE(Colv)) || (is.null(Colv))) && (dendrogram %in%
#                                                  c("both", "column"))) {
#       if (is.logical(Rowv) && (Rowv))
#         dendrogram <- "row"
#       else dendrogram <- "none"
#       warning("Discrepancy: Colv is FALSE, while dendrogram is `",
#               dendrogram, "'. Omitting column dendogram.")
#     }
#   }
#   if (inherits(Rowv, "dendrogram")) {
#     ddr <- Rowv
#     rowInd <- order.dendrogram(ddr)
#   }
#   else if (is.integer(Rowv)) {
#     hcr <- hclustfun(distfun(x))
#     ddr <- as.dendrogram(hcr)
#     ddr <- reorder(ddr, Rowv)
#     rowInd <- order.dendrogram(ddr)
#     if (nr != length(rowInd))
#       stop("row dendrogram ordering gave index of wrong length")
#   }
#   else if (isTRUE(Rowv)) {
#     Rowv <- rowMeans(x, na.rm = na.rm)
#     hcr <- hclustfun(distfun(x))
#     ddr <- as.dendrogram(hcr)
#     ddr <- reorder(ddr, Rowv)
#     rowInd <- order.dendrogram(ddr)
#     if (nr != length(rowInd))
#       stop("row dendrogram ordering gave index of wrong length")
#   }
#   else {
#     rowInd <- nr:1
#   }
#   if (inherits(Colv, "dendrogram")) {
#     ddc <- Colv
#     colInd <- order.dendrogram(ddc)
#   }
#   else if (identical(Colv, "Rowv")) {
#     if (nr != nc)
#       stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
#     if (exists("ddr")) {
#       ddc <- ddr
#       colInd <- order.dendrogram(ddc)
#     }
#     else colInd <- rowInd
#   }
#   else if (is.integer(Colv)) {
#     hcc <- hclustfun(distfun(if (symm)
#       x
#       else t(x)))
#     ddc <- as.dendrogram(hcc)
#     ddc <- reorder(ddc, Colv)
#     colInd <- order.dendrogram(ddc)
#     if (nc != length(colInd))
#       stop("column dendrogram ordering gave index of wrong length")
#   }
#   else if (isTRUE(Colv)) {
#     Colv <- colMeans(x, na.rm = na.rm)
#     hcc <- hclustfun(distfun(if (symm)
#       x
#       else t(x)))
#     ddc <- as.dendrogram(hcc)
#     ddc <- reorder(ddc, Colv)
#     colInd <- order.dendrogram(ddc)
#     if (nc != length(colInd))
#       stop("column dendrogram ordering gave index of wrong length")
#   }
#   else {
#     colInd <- 1:nc
#   }
#   retval$rowInd <- rowInd
#   retval$colInd <- colInd
#   retval$call <- match.call()
#   x <- x[rowInd, colInd]
#   x.unscaled <- x
#   cellnote <- cellnote[rowInd, colInd]
#   if (is.null(labRow))
#     labRow <- if (is.null(rownames(x)))
#       (1:nr)[rowInd]
#   else rownames(x)
#   else labRow <- labRow[rowInd]
#   if (is.null(labCol))
#     labCol <- if (is.null(colnames(x)))
#       (1:nc)[colInd]
#   else colnames(x)
#   else labCol <- labCol[colInd]
#   if (scale == "row") {
#     retval$rowMeans <- rm <- rowMeans(x, na.rm = na.rm)
#     x <- sweep(x, 1, rm)
#     retval$rowSDs <- sx <- apply(x, 1, sd, na.rm = na.rm)
#     x <- sweep(x, 1, sx, "/")
#   }
#   else if (scale == "column") {
#     retval$colMeans <- rm <- colMeans(x, na.rm = na.rm)
#     x <- sweep(x, 2, rm)
#     retval$colSDs <- sx <- apply(x, 2, sd, na.rm = na.rm)
#     x <- sweep(x, 2, sx, "/")
#   }
#   if (missing(breaks) || is.null(breaks) || length(breaks) < 1) {
#     if (missing(col) || is.function(col))
#       breaks <- 16
#     else breaks <- length(col) + 1
#   }
#   if (length(breaks) == 1) {
#     if (!symbreaks)
#       breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm),
#                     length = breaks)
#     else {
#       extreme <- max(abs(x), na.rm = TRUE)
#       breaks <- seq(-extreme, extreme, length = breaks)
#     }
#   }
#   nbr <- length(breaks)
#   ncol <- length(breaks) - 1
#   if (class(col) == "function")
#     col <- col(ncol)
#   min.breaks <- min(breaks)
#   max.breaks <- max(breaks)
#   x[x < min.breaks] <- min.breaks
#   x[x > max.breaks] <- max.breaks
#   if (missing(lhei) || is.null(lhei))
#     lhei <- c(keysize, 4)
#   if (missing(lwid) || is.null(lwid))
#     lwid <- c(keysize, 4)
#   if (missing(lmat) || is.null(lmat)) {
#     lmat <- rbind(4:3, 2:1)
#
#     if (!missing(ColSideColors)) {
#       if (!is.character(ColSideColors) || nrow(ColSideColors) != nc)
#         stop("'ColSideColors' must be a matrix of nrow(x) rows")
#       lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
#       lhei=c(lhei[1], side.height.fraction*ColSideColorsSize/2, lhei[2])
#     }
#
#     if (!missing(RowSideColors)) {
#       if (!is.character(RowSideColors) || ncol(RowSideColors) != nr)
#         stop("'RowSideColors' must be a matrix of ncol(x) columns")
#       lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1), 1), lmat[,2] + 1)
#       lwid <- c(lwid[1], side.height.fraction*RowSideColorsSize/2, lwid[2])
#     }
#     lmat[is.na(lmat)] <- 0
#   }
#
#   if (length(lhei) != nrow(lmat))
#     stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
#   if (length(lwid) != ncol(lmat))
#     stop("lwid must have length = ncol(lmat) =", ncol(lmat))
#   op <- par(no.readonly = TRUE)
#   on.exit(par(op))
#
#   layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
#
#   if (!missing(RowSideColors)) {
#     if (!is.matrix(RowSideColors)){
#       par(mar = c(margins[1], 0, 0, 0.5))
#       image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
#     } else {
#       par(mar = c(margins[1], 0, 0, 0.5))
#       rsc = t(RowSideColors[,rowInd, drop=F])
#       rsc.colors = matrix()
#       rsc.names = names(table(rsc))
#       rsc.i = 1
#       for (rsc.name in rsc.names) {
#         rsc.colors[rsc.i] = rsc.name
#         rsc[rsc == rsc.name] = rsc.i
#         rsc.i = rsc.i + 1
#       }
#       rsc = matrix(as.numeric(rsc), nrow = dim(rsc)[1])
#       image(t(rsc), col = as.vector(rsc.colors), axes = FALSE)
#       if (length(rownames(RowSideColors)) > 0) {
#         axis(1, 0:(dim(rsc)[2] - 1)/max(1,(dim(rsc)[2] - 1)), rownames(RowSideColors), las = 2, tick = FALSE)
#       }
#     }
#   }
#
#   if (!missing(ColSideColors)) {
#
#     if (!is.matrix(ColSideColors)){
#       par(mar = c(0.5, 0, 0, margins[2]))
#       image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
#     } else {
#       par(mar = c(0.5, 0, 0, margins[2]))
#       csc = ColSideColors[colInd, , drop=F]
#       csc.colors = matrix()
#       csc.names = names(table(csc))
#       csc.i = 1
#       for (csc.name in csc.names) {
#         csc.colors[csc.i] = csc.name
#         csc[csc == csc.name] = csc.i
#         csc.i = csc.i + 1
#       }
#       csc = matrix(as.numeric(csc), nrow = dim(csc)[1])
#       image(csc, col = as.vector(csc.colors), axes = FALSE)
#       if (length(colnames(ColSideColors)) > 0) {
#         axis(2, 0:(dim(csc)[2] - 1)/max(1,(dim(csc)[2] - 1)), colnames(ColSideColors), las = 2, tick = FALSE)
#       }
#     }
#   }
#
#   par(mar = c(margins[1], 0, 0, margins[2]))
#   x <- t(x)
#   cellnote <- t(cellnote)
#   if (revC) {
#     iy <- nr:1
#     if (exists("ddr"))
#       ddr <- rev(ddr)
#     x <- x[, iy]
#     cellnote <- cellnote[, iy]
#   }
#   else iy <- 1:nr
#   image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, breaks = breaks, ...)
#   retval$carpet <- x
#   if (exists("ddr"))
#     retval$rowDendrogram <- ddr
#   if (exists("ddc"))
#     retval$colDendrogram <- ddc
#   retval$breaks <- breaks
#   retval$col <- col
#   if (!invalid(na.color) & any(is.na(x))) {
#     mmat <- ifelse(is.na(x), 1, NA)
#     image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "",
#           col = na.color, add = TRUE)
#   }
#   axis(1, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0,
#        cex.axis = cexCol)
#   if (!is.null(xlab))
#     mtext(xlab, side = 1, line = margins[1] - 1.25)
#   axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0,
#        cex.axis = cexRow)
#   if (!is.null(ylab))
#     mtext(ylab, side = 4, line = margins[2] - 1.25)
#   if (!missing(add.expr))
#     eval.parent(substitute(add.expr), n=1)
#   if (!missing(colsep))
#     for (csep in colsep) rect(xleft = csep + 0.5, ybottom = rep(0, length(csep)), xright = csep + 0.5 + sepwidth[1], ytop = rep(ncol(x) + 1, csep), lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
#   if (!missing(rowsep))
#     for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) + 1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) + 1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
#   min.scale <- min(breaks)
#   max.scale <- max(breaks)
#   x.scaled <- scale01(t(x), min.scale, max.scale)
#   if (trace %in% c("both", "column")) {
#     retval$vline <- vline
#     vline.vals <- scale01(vline, min.scale, max.scale)
#     for (i in colInd) {
#       if (!is.null(vline)) {
#         abline(v = i - 0.5 + vline.vals, col = linecol,
#                lty = 2)
#       }
#       xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
#       xv <- c(xv[1], xv)
#       yv <- 1:length(xv) - 0.5
#       lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
#     }
#   }
#   if (trace %in% c("both", "row")) {
#     retval$hline <- hline
#     hline.vals <- scale01(hline, min.scale, max.scale)
#     for (i in rowInd) {
#       if (!is.null(hline)) {
#         abline(h = i + hline, col = linecol, lty = 2)
#       }
#       yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
#       yv <- rev(c(yv[1], yv))
#       xv <- length(yv):1 - 0.5
#       lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
#     }
#   }
#   if (!missing(cellnote))
#     text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote),
#          col = notecol, cex = notecex)
#   par(mar = c(margins[1], 0, 0, 0))
#   if (dendrogram %in% c("both", "row")) {
#     plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
#   }
#   else plot.new()
#   par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
#   if (dendrogram %in% c("both", "column")) {
#     plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
#   }
#   else plot.new()
#   if (!is.null(main))
#     title(main, cex.main = 1.5 * op[["cex.main"]])
#   if (key) {
#     par(mar = c(5, 4, 2, 1), cex = 0.75)
#     tmpbreaks <- breaks
#     if (symkey) {
#       max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
#       min.raw <- -max.raw
#       tmpbreaks[1] <- -max(abs(x), na.rm = TRUE)
#       tmpbreaks[length(tmpbreaks)] <- max(abs(x), na.rm = TRUE)
#     }
#     else {
#       min.raw <- min(x, na.rm = TRUE)
#       max.raw <- max(x, na.rm = TRUE)
#     }
#
#     z <- seq(min.raw, max.raw, length = length(col))
#     image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks,
#           xaxt = "n", yaxt = "n")
#     par(usr = c(0, 1, 0, 1))
#     lv <- pretty(breaks)
#     xv <- scale01(as.numeric(lv), min.raw, max.raw)
#     axis(1, at = xv, labels = lv)
#     if (scale == "row")
#       mtext(side = 1, "Row Z-Score", line = 2)
#     else if (scale == "column")
#       mtext(side = 1, "Column Z-Score", line = 2)
#     else mtext(side = 1, KeyValueName, line = 2)
#     if (density.info == "density") {
#       dens <- density(x, adjust = densadj, na.rm = TRUE)
#       omit <- dens$x < min(breaks) | dens$x > max(breaks)
#       dens$x <- dens$x[-omit]
#       dens$y <- dens$y[-omit]
#       dens$x <- scale01(dens$x, min.raw, max.raw)
#       lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol,
#             lwd = 1)
#       axis(2, at = pretty(dens$y)/max(dens$y) * 0.95, pretty(dens$y))
#       title("Color Key\nand Density Plot")
#       par(cex = 0.5)
#       mtext(side = 2, "Density", line = 2)
#     }
#     else if (density.info == "histogram") {
#       h <- hist(x, plot = FALSE, breaks = breaks)
#       hx <- scale01(breaks, min.raw, max.raw)
#       hy <- c(h$counts, h$counts[length(h$counts)])
#       lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s",
#             col = denscol)
#       axis(2, at = pretty(hy)/max(hy) * 0.95, pretty(hy))
#       title("Color Key\nand Histogram")
#       par(cex = 0.5)
#       mtext(side = 2, "Count", line = 2)
#     }
#     else title("Color Key")
#   }
#   else plot.new()
#   retval$colorTable <- data.frame(low = retval$breaks[-length(retval$breaks)],
#                                   high = retval$breaks[-1], color = retval$col)
#   invisible(retval)
# }



