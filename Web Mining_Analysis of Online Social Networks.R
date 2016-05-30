install.packages("igraph");require("igraph");
install.packages("readr");require("readr");
dat=read.csv(file.choose("C://Assignments/out.youtube-links.csv"),header=FALSE,check.names=FALSE,sep=" ",skip=2)
dat<-dat[,1:2]
colnames(dat) <- c("From", "To")

g=graph.data.frame(dat[,1:2],directed=TRUE)
nodes<-V(g)
links<-E(g)

#plot(g, edge.arrow.size=.4,vertex.label=NA)
#Create Graph Object
g <- simplify(g, remove.multiple = F, remove.loops = T) 

###1.LINK SYMMETRY###http://kateto.net/networks-r-igraph
#The proportion of reciprocated ties (for a directed network).##
rec<-reciprocity(g,ignore.loops = TRUE)#Calculates the reciprocity of a directed graph.

#Alternate method

rec_calc<-2*dyad_census(g)$mut/ecount(g) # Calculating reciprocity. 
#dyad_census(g)  - Mutual, asymmetric, and null node pairs
#ecount(g) gives the number of edges in the graph.

#Output for Link Symmetry##
cat("Number of Users = ",vcount(g))
cat("Number of friend link = ",ecount(g))
cat("Fraction of Links symmetric = ",rec_calc*100,"%")

###2.Power Law Node degrees###

#Log-log plot of In degree of Complementary cumulative distribution functions
yin=degree.distribution(g,cumulative = TRUE,mode = c( "in"))
plot(yin,log="xy",col = "red",ylab = "P(Indegree>=x)", xlab = "Degree (log scale)", main = "Log-log plot of in-degree distribution")


#Log-log plot of Out degree of Complementary cumulative distribution functions
yout=degree.distribution(g,cumulative = TRUE,mode = c( "out"))
plot(yout,log="xy",col = "red",ylab = "P(outdegree>=x)", xlab = "Degree (log scale)", main = "Log-log plot of Out-degree distribution")

#Power-law coefficient estimates (??) and corresponding Kolmogorov-Smirnov goodness-of-fit metrics (D)
#have used plfit
#The 'plfit' implementation also uses the maximum likelihood principle to determine alpha for a given xmin; When xmin is not given in advance, the algorithm will attempt to find itsoptimal value for which the p-value of a Kolmogorov-Smirnov test between the fitted distribution and the original sample is the largest. The function uses the method of Clauset, Shalizi and Newman to calculate the parameters of the fitted distribution.

#In-degree
d <- degree(g, mode="in")
fit1 <- power.law.fit(d, 1)#1 - The lower bound for fitting the power-law
cat("Alpha = ")
fit1$alpha#the exponent of the fitted power-law distribution.
cat("test statistic of a Kolmogorov-Smirnov test = ")
fit1$KS.stat#the test statistic of a Kolmogorov-Smirnov test that compares the fitted distribution with the input vector. Smaller scores denote better fit.

#Out-degree
d <- degree(g, mode="out")
fit2 <- power.law.fit(d, 1)#1 - The lower bound for fitting the power-law
cat("Alpha = ")
fit2$alpha#the exponent of the fitted power-law distribution.
cat("test statistic of a Kolmogorov-Smirnov test = ")
fit2$KS.stat#the test statistic of a Kolmogorov-Smirnov test that compares the fitted distribution with the input vector. Smaller scores denote better fit.


###3. Corelation between indegree and outdegree
d.in <- degree(g, mode="in")
d.in<-sort(d.in,decreasing=TRUE)#sorting the degrees of nodes
d.out <- degree(g, mode="out")
d.out<-sort(d.out,decreasing=TRUE)#sorting the degrees of nodes

#calculation for x%
final<-NULL
fract_users<-seq(0.01,1.0,by=0.2)#x%
for(i in fract_users){
  d.in.top<-d.in[c(1:ceiling(i*length(d.in)))]
  d.out.top<-d.out[c(1:ceiling(i*length(d.out)))]
  count =0;
  for(n in d.in.top){
    if(is.element(n,d.out.top)==TRUE){
      count = count+1
    }
  }
  overlap1 = count/(length(d.in.top))
  final<-rbind(final,overlap1)
  #append(ratio_overlap,count/(length(d.in.top)))
  print(final)
  print(count)
  print(overlap1)
  
  
}
plot(x=fact,y=overlap,type="l",col = "green",ylab = "Overlap", xlab = "Fraction of Users(Ordered by In-degree)", main = "Plot of Overlap between top x% of Nodes(ranked by In-degree)")


###4. Average path length, radius and diamter
apl<-average.path.length(g, directed=TRUE)
apl
rad<-radius(g, mode = c("all"))
rad
diam<-diameter(g, directed = TRUE, unconnected = TRUE, weights = NULL)
diam

###5 Link Degree correlation
#5.1 Joint Degree Distribution

d.out.degree <- degree(g, mode="out")
d.out.degree<-sort(d.out.degree,decreasing=FALSE)#sorting the degrees of nodes

d.out.degree.sub<-head(d.out.degree,1000)
knn<-NULL
attr_vectr<-as.numeric(unlist(attributes(d.out.degree.sub)))
for(s in attr_vectr ){

  print(s)
  neigh_vector = neighbors(g, s, mode = "in")
  sum_1=0
  for(n in neigh_vector){
    sum_1 = sum_1 + degree(g,n, mode="in")
  }
  avg_indegree = sum_1/length(neigh_vector)
  knn<-rbind(knn,avg_indegree)
  
}

#incident(g, v, mode=c("all", "out", "in", "total"))

###Densly Connected Core
d.high.degree <- degree(g, mode="all")
d.high.degree.sort<-sort(d.out.degree,decreasing=TRUE)
high.degree.index<-as.numeric(unlist(attributes(d.high.degree.sort)))

#
total<-length(d.high.degree)
fract<-c(0.01,0.1,1,10)
final<-NULL
for (f in fract){
densly.connected.core<-high.degree.index[c(1:(f*total)/100)]
v<-c(densly.connected.core)
sub_g<-induced_subgraph(g, -v)
clusters_weak = no.clusters(sub_g,mode=c("weak"))
clusters_strong = no.clusters(sub_g,mode=c("strong"))
cluster.num<-cbind(f,clusters_weak,clusters_strong)
final<-rbind(final,cluster.num)
}

weak_cc_fract<-c(0.004,0.012,0.059,0.222)
strong_cc_fract<-c(0.996,0.988,0.941,0.778)
count<-cbind(weak_cc_fract,strong_cc_fract)
barplot(as.matrix(count))

strongclusters <- clusters(sub_g, mode="strong")$membership
clusters.with.size <- clusters(sub_g, mode="strong")$csize

#100 - 2298
#1000 - 6868
#10000 - 35966
#100000 - 175335

#gfam <- subgraph.edges(graph=g, eids=which(E(g)$label=="FAM"), delete.vertices = TRUE)