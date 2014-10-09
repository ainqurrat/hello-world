############ Q_AIN (27-04-2014) ####################################
#####################################################################
################ HEATMAP for SIMILARITY MATRIX #####################
####################################################################
packages <- list('caret', 'kernlab','pROC','doMC','randomForest','elasticnet','gam','ipred','fastICA','pls','ggplot2')

# Install packages
lapply(packages, FUN = function(packages) {
  do.call("install.packages", list(packages)) 
})

# Load packages
lapply(packages, FUN = function(packages) {
  do.call("library", list(packages)) 
})
require("gplots")
###############################################################
###############################################################


# Tanimoto Distance between two vectors
tanimoto_bw <- function(a,b){
  
  if (length(a) != length(b)){print("vectors of unqual length");break}
  differentes = sum((a == b)*1) 
  comunes = length(a)
  tani =  differentes / comunes
  return(tani)
}

############################################################
###########################################################

################## Heatmap of Compounds #######################
df <- read.csv("input_file.csv") # input your fingerprint file
killset <- expression(c(Structure, Parent.ID, Target, Emax,Efficacy))  ######## Removes unnecessary columsn
x <- subset(df, select = -eval(killset))
partial_sim <- c()
m_sims <- matrix(0,156,156) # mention the dimension of your matrix (no. of rows and columns)

m_sims<-as.matrix(x)
for (i in 1:nrow(m_sims)){ 
   seq_now = x[i,]
   partial_sim <- c()
   for(j in 1:ncol(m_sims)){
     m_sims[i,j] <- tanimoto_bw(seq_now, x[j,]) #calculation of similarity matrix
     }
}



library(RColorBrewer)
pdf("~/Dropbox/Neusentis/11sep/functional_sim.pdf", pointsize=10, width=15, height=10)
result<-heatmap.2(x, Rowv=T,Colv=T,scale='none', dendrogram="none", symm=T, col=brewer.pal(9, "Reds"), tracecol=NULL, labRow=df$CMPD_CHEMBLID,labCol=df$CMPD_CHEMBLID, cexRow=0.5,key=F,main="")#m_sim is your similarity matrix
#labRow and labCol is your variable for labelling rows and columns

dev.off()
save.image('~/Desktop/output.RData')

#####################################################################
################ HIERARCHAL ClUSTERING #############################
####################################################################
require(graphics)
hc<-hclust(dist(x), "ward") ####################### Method could be changed
hc$labels<-df$Species ### You can set the labels of leaves of trees here and plot
plot(hc, hang=-1)

##################################################################
################ Principal Component Analysis#####################
####################################################################

######## Remove descriptors with zero-variance ##############
nzv.columns <- nearZeroVar(x, freqCut = 30/1)
nzv.names <- names(x)[nzv.columns]
x <- x[, -nzv.columns]
ncol(x)
hept_pca <- prcomp(x, retx=TRUE, center=TRUE,scale.=TRUE) 
plot(hept_pca)  #Scree plot , variance of each PC
scores <- hept_pca$x #Variance
loadings <- hept_pca$rotation #loading plot of descriptors
Targets<-df$subtype
Specie<-df$SPECIES

qplot(scores[,1], scores[,2],data=x, color=Targets,alpha=0.5,size=I(3),label=df$CMPD_CHEMBLID, xlab="PC1 (34%)", ylab="PC2 (26%)", main="")+theme_bw()+theme(legend.position="bottom")+guides(col = guide_legend(nrow = 1))+theme(legend.text=element_text(size=16))+theme(axis.text.x=element_text(size=14,face="bold"),axis.text.y=element_text(size=14,face="bold"),axis.title.x = element_text(colour="grey20",size=18,angle=0,hjust=.5,vjust=0,face="bold"),axis.title.y = element_text(colour="grey20",size=18,angle=90,hjust=.5,vjust=.5,face="bold"))+geom_text()
                                                                                                                                                                                  
biplot(hept_pca, col=c("red","black")) # this plot gives variance and loadings in one plot so might be useful to try


##################################################################
################ Multidimensional Scaling     #####################
####################################################################
d<-dist(x)
fit<-cmdscale(d, k = 5, eig = TRUE, add = FALSE, x.ret = FALSE)
#in the plot below, you can color accoridng to your vairables
qplot(fit$points[,2],fit$points[,3], colour=df$,size=I(6), xlab="Coordinate1", ylab="Coordiante2", main="")+scale_colour_gradient2(limit=c(-1,3))+opts(panel.background = theme_rect(fill='white', colour='white'))+theme(legend.position="right" )

#################################################################


library(ggplot2)
library(plyr)
library(reshape2)
df <- read.csv("input_file.csv")
melted <- melt(df, id.vars=c("Target"))
mean <- ddply(melted, c("Target", "variable"), summarise, mean=mean(value))
means.barplot <- qplot(x=Target, y=mean, fill=variable,
                       data=mean, geom="bar", stat="identity",
                       position="dodge", ylab="",xlab="")+theme_bw()+theme(axis.text.x=element_text(size=25,angle=45,face="bold"),axis.text.y=element_text(size=25,face="bold"))
means.sem <- ddply(melted, c("Target", "variable"), summarise,
                   mean=mean(value), sem=sd(value)/sqrt(length(value)))
#means.sem <- transform(means.sem, lower=mean-sem, upper=mean+sem)
means.sem <- read.csv("means.sem.csv") # I input means_sem separately (for Figure 3) as I have already calculated them in different file
means.barplot + geom_errorbar(aes(ymax=upper,ymin=lower),position=position_dodge(0.9),data=means.sem)
##################################################################
ggplot(df, aes(x=df$V1))+theme_bw() + ylab("Density") + xlab("Bioaffinity") + geom_density(alpha=.3)+ theme_bw()+theme(legend.position="right", legend.key=element_rect(size=5))+theme(axis.text.x=element_text(size=15,face="bold"),axis.text.y=element_text(size=15,face="bold"),axis.title.x = element_text(colour="grey20",size=25,angle=0,hjust=.5,vjust=0,face="plain"),
                                                                                                                                                                                                                       axis.title.y = element_text(colour="grey20",size=25,angle=90,hjust=.5,vjust=.5,face="plain"))
