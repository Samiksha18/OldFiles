####  K-mean Clustering in ggPlot
#     Nick Waters
#       20150121
#        Version 1
#Sub Source:http://www.r-bloggers.com/bot-botany-k-means-and-ggplot2/ 
#Main Source: http://rstudio-pubs-static.s3.amazonaws.com/2598_ccd642fc32854463844c6f9cc153983a.html



library(ggplot2)
library(plyr)
require(graphics)

#numeric variables of interest   TWEAK THIS LATER0
var1<-"carat"
var2<-"price"
var3<-"table"
#catagorical variables of interest
catvar<-"color"

#open file, but for now use diamonds dataset
file1<-diamonds
sub1<-file1[file1$cut=='Ideal', ]
sub1.5<-sub1[c(var1, var2, var3,catvar)]
sub2<-sub1[c(var1, var2, var3)]


#plot basic 
ggplot(sub2, aes(x=sub2[[var1]], y=sub2[[var2]])) + geom_point()+
  labs(title=(paste("Distribution of ",var1," vs ",var2)), x=var1, y=var2
)

#######################Graph to determine centers
dia.kmeans<-lapply(1:13, function(i){
  kmeans(sub2[,c(var1, var2)], centers = i)
})

lapply(dia.kmeans, function(x) x$withinss)

dia.within.ss <- sapply(dia.kmeans, function(x) sum(x$withinss))
ggplot(data.frame(cluster = 1:13, within.ss = dia.within.ss), aes(cluster, within.ss)) +
  geom_point() + geom_line() +
  scale_x_continuous(breaks = 0:13)
##################

clustRange<-1:3   ### Set this based on above plot's inflection point
dia.kmeans<-lapply(clustRange, function(i){
  kmeans(sub2[,c(var1, var2)], centers = i)})
cluster.colors <- lapply(dia.kmeans, function(x) x$cluster)


l_ply(cluster.colors,
      function(colors) {
        plot.dat <- cbind(sub2, cluster = factor(colors))
        gg.obj <- ggplot(plot.dat, aes(x=sub2[[var1]], y=sub2[[var2]], color = cluster)) +
          geom_point() + labs(title = paste(nlevels(factor(colors))),x=var1, y=var2)        
        print(gg.obj)
})

#####    Plotting centers for reference
sub3<-kmeans(sub2, 3)    ####SET RANGE####
sub4<-sub2
sub4[[catvar]]<-sub1.5[[catvar]]
sub4$cluster<-factor(sub3$cluster)
centers<-as.data.frame(sub3$centers)

#####Center Focus
ggplot(data=sub4, aes(x=sub4[[var1]], y=sub4[[var2]], color = sub4$cluster)) + 
  geom_point() + 
#  facet_wrap(~color)+
  geom_point(data=centers, aes(x=centers[[var1]], y=centers[[var2]], color='Center'))+
  geom_point(data=centers, aes(x=centers[[var1]], y=centers[[var2]], size=152, alpha=.3))


makeFacetClusters<- function(x){
  
}


ggplot(data=sub4, aes(x=sub4[[var1]], y=sub4[[var2]], color = sub4$cluster)) + 
  geom_point() + 
  facet_wrap(~color)

  