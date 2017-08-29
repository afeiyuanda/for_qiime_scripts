#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
my %opts;
my $VERSION = "2013-1-4";
GetOptions (\%opts,"i=s","n=i","m=i","g=s","ncol=i","bw=i","bh=i","pw=i","ph=i","hi=s","lcex=f","l1=f","lx=s","p=f","s=i","ltp=s","od=s","bar=s","pie=s","pcex=f","pxlim=f","pylim=f");

my $usage = <<"USAGE";
        Program : $0
        Version : $VERSION
        Contact : guo.yu\@majorbio.com
        Discription: plot bar and pie chart
        Usage:perl $0 [options] 
             base options:
                -i	    inputfile,a frame.
                -n	collect the max n in each sample
                -m	max m in all samples' sum
                -p	if num/all<=p push to other  defalt :0.01
                -s	if num<=s  push to other     defalt :0  
                -od	sort table with samples sum, [d/i/n],d:decreasing; i:increasing; n:not do anything; default:d
                -g	calculte sum of groups for samples eg. [   a1 A
                                              				a2 A
                                              				a3 A
                                              				b1 B
                                              				b2 B ...
                                           			    ]
                -bar  [T/F]plot bar or not . defult :T  
                -pie [T/F]plot pie for each sample or not . defult :T   
                
             bar options :                        			       
                -ncol  colums of legend defalt=3
                -lcex  legend cex defalt:0.4 #control legend's size
                -l1    make the width of legend colume1 large defalt:1.3
                -bw  width  defalt=8
                -bh  height    defalt=7
                -hi split of height 1:1.2 
                -lx [T/F] make legend looks good  defalt :F 

             pie options :
             	 -pw  	width      default=7
             	 -ph        height    default=5
             	 -pcex     label cex .defalt:0.8
                 -ltp          label type [t/l],t:text around pie ;l :legend on the pie right.defalt :t 
				-pxlim		xlim of pie;default:1
				-pylim	    ylim of pie:defualt:1
USAGE
die $usage if ( !$opts{i});

#define defalts
$opts{n}=defined $opts{n}?$opts{n}:-1;
$opts{m}=defined $opts{m}?$opts{m}:-1;
$opts{p}=defined $opts{p}?$opts{p}:0.01;
$opts{s}=defined $opts{s}?$opts{s}:-1;
$opts{ncol}=defined $opts{ncol}?$opts{ncol}:2;
$opts{lcex}=defined $opts{lcex}?$opts{lcex}:0.3;
$opts{pcex}=defined $opts{pcex}?$opts{pcex}:0.3;
$opts{l1}=defined $opts{l1}?$opts{l1}:1.3;
$opts{bw}=defined $opts{bw}?$opts{bw}:8;
$opts{bh}=defined $opts{bh}?$opts{bh}:7;
$opts{pw}=defined $opts{pw}?$opts{pw}:7;
$opts{ph}=defined $opts{ph}?$opts{ph}:4;
$opts{hi}=defined $opts{hi}?$opts{hi}:"1:1.2";
$opts{g}=defined $opts{g}?$opts{g}:"ALL";
$opts{lx}=defined $opts{lx}?$opts{lx}:"F";
$opts{od}=defined $opts{od}?$opts{od}:"d";
$opts{ltp}=$opts{ltp}?$opts{ltp}:"t";
$opts{bar}=$opts{bar}?$opts{bar}:"T";
$opts{pie}=$opts{pie}?$opts{pie}:"T";
$opts{pxlim}=defined $opts{pxlim}?$opts{pxlim}:1;
$opts{pylim}=defined $opts{pylim}?$opts{pylim}:1;

my $labtext;
my $legend;
if($opts{ltp}=~/^t$/){
          $labtext="TRUE";$legend="FALSE";
}elsif($opts{ltp}=~/^l$/){
          $labtext="FALSE";$legend="TRUE";
}else{
          print "opts -ltp must be 't' or 'l'!\n";exit;
}

$opts{g}=~/\/*([^\/]+)$/;
my $gs=$1;
my $del=0;
open RCMD, ">cmd.r";

print RCMD "
#library(RColorBrewer)
mycol <-c(119,132,147,454,89,404,123,529,463,104,552,28,54,84,256,100,558,43,652,31,610,477,588,99,81,503,562,76,96,495)
mycol <-colors()[rep(mycol,20)]
#mycol <- c( brewer.pal(8,\"Dark2\"),brewer.pal(12,\"Set3\"),brewer.pal(8,\"Set2\"),brewer.pal(8,\"Accent\"))

otu <-read.table(file=\"$opts{i}\",header=T,check.names=FALSE,comment.char=\"\",quote=\"\",sep=\"\t\")
rownames(otu) <- otu[,1]
rownames(otu) <-sapply(rownames(otu),function(x) gsub(\"_*{\.+}\",\" \",x,perl = TRUE)) 
otu <-otu[,-1]

#al <- which(rownames(otu) \%in% c(\"All\",\"No_Rank\",\"Trimed\"))
al <- which(rownames(otu) \%in% c(\"All\"))
if(length(al)) otu <-otu[-al,]

gs <-\"$opts{g}\"
if(gs!=\"ALL\"){
      group <- read.table(\"$opts{g}\")
      glst <- lapply(1:length(unique(group[,2])),function(x)group[which(group[,2] \%in\% unique(group[,2])[x]),1])
      names(glst) <-as.character(unique(group[,2]))
      tab <-sapply(1:length(glst),function(x) apply(otu[as.character(as.vector(glst[[x]]))],1,sum))
      otu <-tab[apply(tab,1,function(x)any(x>0)),]      
      colnames(otu) <-unique(group[,2])
}

n <-$opts{n}
m <-$opts{m}
p <-$opts{p}
s <-$opts{s}
od<-\"$opts{od}\"

rowsum <-sapply(1:nrow(otu),function(x) sum(otu[x,]))

if(od ==\"d\"){
	otu<-otu[order(rowsum,decreasing=TRUE),]
}else if(od ==\"i\"){
	otu<-otu[order(rowsum,decreasing=FALSE),]
}

if(n>0){
     #otu[order(otu[,1],decreasing=T)[1:10],1]
	
     max_names <- sapply(1:ncol(otu),function(x) rownames(otu)[order(otu[,x],decreasing=T)[1:n]])
     otu <-otu[which(rownames(otu) \%in% unique(as.vector(max_names))),]
     
     #unique(as.vector(max_names))
     #use_list <-which(rownames(otu) \%in% unique(as.vector(max_names)))
     #rownames(otu)[which(rownames(otu) \%in% unique(as.vector(max_names)))]
    
}else if(m>0){
        otu<-otu[order(rowsum,decreasing=TRUE),]
		otu<-otu[1:m,] 
}else if(p >0){
	 otu_pec <- otu
	 otu_pec <- sapply(1:ncol(otu),function(x) otu_pec[,x] <-otu[,x]/sum(otu[,x]))
	 minp <-sapply(1:nrow(otu_pec),function(y) all(otu_pec[y,]<=p)) 	 
	 otu_xp <-otu[minp,]
	 other <- sapply(1:ncol(otu_xp),function(y) sum(otu_xp[,y]))
	 otu <-rbind(otu[!minp,],other)
	 rownames(otu)[nrow(otu)] <-\"Others\"
}else if(s >0){
	 minp <-sapply(1:nrow(otu),function(y) all(otu[y,]<=s)) 	 
	 otu_xp <-otu[minp,]
	 other <- sapply(1:ncol(otu_xp),function(y) sum(otu_xp[,y]))
	 otu <-rbind(otu[!minp,],other)
	 rownames(otu)[nrow(otu)] <-\"Others\"	 
}

del <-unlist(sapply(1:ncol(otu),function(x) if(sum(otu[,x])==0) x))

if(!is.null(del)) {
  mydel <-colnames(otu)[c(del)]
  otu <-otu[,-c(del)]
}

bar <-\"$opts{bar}\"
pie <-\"$opts{pie}\"
################# plot  bar     ########################################
if(bar ==\"T\"){
	dat <-sapply(1:ncol(otu),function(x) otu[,x]/sum(otu[,x]))
	#print(dat),assignment sample name to the column
	#colnames(dat) <-colnames(otu)
	rownames(dat) <-rownames(otu)
	
	my_dat = dat
	colnames(my_dat) <-colnames(otu)
	write.table(my_dat,file=\"average.xls\",sep=\"\\t\",col.names=NA,row.names=TRUE,quote=FALSE)
	
	#plot.new()
	lab <-rownames(dat)
	#print(lab)
	str <-strwidth(lab, units = \"inches\")
	#print(str)
	#length(str)
	cnum <-ceiling(length(str)/$opts{ncol})
	cnum
	nstr <- lapply(1:$opts{ncol},function(x) if(cnum*x<=length(str)){str[(cnum*(x-1)+1):(cnum*x)]}else{str[(cnum*(x-1)+1):length(str)]})
	nlab <- lapply(1:$opts{ncol},function(x) if(cnum*x<=length(lab)){lab[(cnum*(x-1)+1):(cnum*x)]}else{lab[(cnum*(x-1)+1):length(lab)]})
	mw <-unlist(lapply(nstr,max))
	mw[1] <-mw[1]*$opts{l1}

	lx <-\"$opts{lx}\"
	#tiff(file=\"bar.$gs.$opts{i}.tiff\",width=$opts{bw},height=$opts{bh},pointsize=15)
	pdf(file=\"bar.$gs.$opts{i}.pdf\",width=$opts{bw},height=$opts{bh})
	if(lx==\"T\"){
	  la <-layout(matrix(c(rep(1,$opts{ncol}),2:($opts{ncol}+1)),2,$opts{ncol},byrow=TRUE),width=mw,heights=c($opts{hi}))
	}else{
	  layout(matrix(1:2,2,1),heights=c($opts{hi}))
	}
	par(mar=c(1,5,2,2))
	las <-1
	if(ncol(dat)>15) las <-3
	bar_plot <- barplot(dat*100,width=0.3,space=0.2,plot=T,las=las,yaxt=\"n\",col=mycol[1:nrow(dat)],cex.names=1,border=NA,offset=0)#cex.axis control the font of the axis's kedu.
	axis(side=2,lwd.ticks=0.5,tck=-0.02,lwd=0.5,labels=FALSE,line=-0.3)
	#line控制纵坐标轴与条形图的距离,正值是向左移动
	mtext(side=2,seq(0,100,20),at=seq(0,100,20),cex=0.3,line=-0.3)
	#tck指定轴上刻度长度的值,单位是百分比,以图形宽、高中最小一个作为基数; 如果tck=1则绘制grid
	#yaxt=\"n\"不显示y轴
	#add colunm name to the bar
	text(cex=0.2,x=bar_plot, y=-2, xpd=TRUE, lab=colnames(otu), srt=45, adj=0.5, font=2) # srt control the angle, y control the distance between column name and the x axis
	#1正常，2加粗，3斜体，4加粗，斜体，5符号
	#adj控制关于文字的对齐方式,0是左对齐,0.5是居中对齐,1是右对齐,值> 1时对齐位置在文本右边的地方,取负值时对齐位置在文本左边的地方
	#cex.axis控制坐标轴刻度数字大小,cex.lab控制坐标轴标签文字大小,cex.main控制标题文字大小,cex.sub控制副标题文字大小
	#las控制坐标轴刻度数字标记方向的整数(0: 平行于轴,1: 横排,2: 垂直于轴,3: 竖排)
	#side on which side of the plot (1=bottom, 2=left, 3=top, 4=right).
	mtext(1, text = \"Sample ID\", line = -0.1, cex=0.3) #cex control the font size, line control the distance between text and the x axis,negtive=up, positive=down
	mtext(\"Relative abundance (%)\",side=2,line=0.2,font=0.3,cex=0.3)  #add y lable
	
	#box()
	if(lx==\"T\"){
	 for(i in 1:$opts{ncol}){
	   if(i==1){par(mar=c(1,5,1,0), xpd=TRUE)}
	   else {par(mar=c(1,0,1,0), xpd=TRUE)}
	   #if(i==$opts{ncol}){par(mar=c(1,0,1,2))}
	   plot.new()
	   legend(\"topleft\",legend=nlab[[i]],fill=mycol[(cnum*(i-1)+1):(cnum*i)],cex=$opts{lcex},bty=\"n\",border=\"white\")
	 }
	}else{
	 par(mar=c(2,5,0,1), xpd=TRUE)
	 plot.new()
	 legend(\"topleft\",legend=rownames(dat),ncol=$opts{ncol},fill=mycol[1:nrow(dat)],cex=$opts{lcex},bty=\"n\",pt.cex=1,border=\"transparent\") #cex contrl the size of the legend
	 #legend(x=0,y=1.0,legend=rownames(dat),ncol=4,fill=mycol[1:nrow(dat)],cex=0.5,bty=\"n\"),border=\"white\"
	}
	dev.off()
	 
	if(!is.null(del)){
		prt <-paste(\"Warning: Samples:\",mydel,\"were none,disregard them in the plot.\")
		print(prt)
	}

}

######## bar finished ####################
if(pie ==\"T\"){

######## function spie ##########
spie <-function (x, labels = names(x), edges = 200, radius = 0.8, clockwise = FALSE, 
    init.angle = if (clockwise) 90 else 0, density = NULL, angle = 45, 
    col = NULL, border = NULL, lty = NULL, main = NULL, textlab=TRUE, legend = FALSE,...) 
{   
    if (!is.numeric(x) || any(is.na(x) | x < 0)) 
        stop(\"\'x\' values must be positive.\")
    if (is.null(labels)) 
        labels <- as.character(seq_along(x))
    else labels <- as.graphicsAnnot(labels)

    if(legend) { layout(matrix(1:2,1,2) ) 
                 par1 <-par(mar=c(1,1,1,0))
                 par2 <-par(mar=c(1,0,4,1))
    }

    x <- c(0, cumsum(x)/sum(x))
    dx <- diff(x)
    nx <- length(dx)
    plot.new()
    pin <- par(\"pin\")
    xlim <- ylim <- c(-1, 1)
    if (pin[1L] > pin[2L]) 
        xlim <- (pin[1L]/pin[2L]) * xlim
    else ylim <- (pin[2L]/pin[1L]) * ylim
    vx<-$opts{pxlim}+0.3
	vy<-$opts{pylim}+0.3
	ylim[2] <-ylim[2]*vy	
    if(legend) {v <- 0.8;par(par1)}
    plot.window(xlim*vx, ylim, \"\", asp = 1)
    if (is.null(col)) 
        col <- if (is.null(density)) 
            c(\"white\", \"lightblue\", \"mistyrose\", \"lightcyan\", 
                \"lavender\", \"cornsilk\")
        else par(\"fg\")
    col <- rep(col, length.out = nx)
    border <- rep(border, length.out = nx)
    lty <- rep(lty, length.out = nx)
    angle <- rep(angle, length.out = nx)
    density <- rep(density, length.out = nx)
    twopi <- if (clockwise) 
        -2 * pi
    else 2 * pi
    t2xy <- function(t) {
        t2p <- twopi * t + init.angle * pi/180
        list(x = radius * cos(t2p), y = radius * sin(t2p))
    }

    pl1 <- pl2 <- c(0,-radius*1.5)
    lb1 <- lb2 <-array()
    lb1[1] <- lb2[1] <- \"lab0\"
    li1 <- li2 <-1
    for (i in 1L:nx) {
        n <- max(2, floor(edges * dx[i]))
        P <- t2xy(seq.int(x[i], x[i + 1], length.out = n))
        polygon(c(P\$x, 0), c(P\$y, 0), density = density[i], angle = angle[i], 
            border = border[i], col = col[i], lty = lty[i])
        pm <- t2xy(mean(x[i + 0:1]))
        lab <- as.character(labels[i])
        if (!is.na(lab) && nzchar(lab)) {
            if(pm\$x >=0){ 
                 l1 <-c(pm\$x,pm\$y)
                 pl1 <-rbind(pl1,l1)
                 li1 <-li1+1
                 lb1[li1] <-lab
            }else { print (lab) 
                 print (pm\$x)
                 l2 <-c(pm\$x,pm\$y)
                 pl2 <-rbind(pl2,l2)
                 li2 <-li2+1
                 lb2[li2] <-lab
            }
        }
    }
    row.names(pl1) <-lb1
    row.names(pl2) <-lb2
    print(pl1)
    print(pl2)

############### labelxy ####
    labelxy <- function(pl) {
         pl <-pl[order(pl[,2]),]
         by <-pl[1,2]
         bx <-pl[1,1]
         d1 <-1.3
         for(j in 2:nrow(pl)) {
             ply <- pl[j,2]
             plx <- pl[j,1]
             pmy <- 1.3 * ply
             pmx <- 1.3 * plx
             d2 <-1.3
             if(ply<0 & abs(plx)<0.3*radius){ 
                  if(d1>1.1) {d1 <-d1-0.04}                  
                  while(abs(pmx-bx)<0.1 ){                   
                     d2 <-d2+0.01
                     pmx <- plx*d2}
             }else{ d1 <-1.2 }
             while(pmy-by<0.08){ pmy <- pmy+0.01 }
             d3 <-d2
             while(pmy*pmy+pmx*pmx <=radius*radius){
                  d3 <-d3+0.01
                  pmx <- plx*d3
             } 
             
             lines(c(1, d1) * plx,c(1, d1) * ply) #line1                                                  
             lines(c(d1* plx, pmx) , c(d1*ply, pmy) ) #line2
             lines(c(pmx, pmx*1.03) , c(pmy,pmy ))  #line3
             text(1.036*pmx, pmy, rownames(pl)[j], xpd = TRUE, adj = ifelse(plx < 0, 1, 0), ...)
             by <-pmy
             bx <-pmx
         }
    }
###############

  if(textlab) {       
            labelxy(pl1)
            labelxy(pl2)
  }
  title(main = main, ...)
  if(legend) { plot.new()
               par(par2)
               legend(\"topleft\",legend=paste(labels,dat,sep=\" \"),fill=col)

  }


  invisible(NULL)


}

#####################################################################

##### function ppie : plot a table of samples #########
ppie <-function(dat,label,col,smp){
	#tiff(paste(\"pie.\",smp,\".$opts{i}.tiff\",sep=\"\"),width=$opts{pw},height=$opts{ph},pointsize=15)
	pdf(paste(\"pie.\",smp,\".$opts{i}.pdf\",sep=\"\"),width=$opts{pw},height=$opts{ph})
	label <-sapply(1:nrow(dat),function(x) paste(label[x],\" \",round(dat[x,1]/sum(dat)*100,digits=2),\"%\",sep=\"\"))
	spie(dat,border=NULL,labels=label,col=col,main=smp,cex=$opts{pcex},textlab=$labtext,legend=$legend)
	dev.off()
}
###########################################

#########plot    pie             ##########
pcol=mycol[1:nrow(otu)]
dl <-lapply(1:ncol(otu),function(y) which(otu[,y]>0))
sapply(1:ncol(otu),function(y) ppie(dat=as.matrix(otu[dl[[y]],y]),label=rownames(otu)[dl[[y]]],col=pcol[dl[[y]]],smp=colnames(otu)[y]))
}


";


`R --restore --no-save < cmd.r`;
#system ('rm cmd.r');
















