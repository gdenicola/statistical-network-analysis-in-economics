#loading necessary libraries
library(statnet)
library(dplyr)
library(amen)
library(abind)

#setting the path to current file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 

#loading the data
flandreau_jobst_internationalcurrencies_data <- read.delim("flandreau_jobst_internationalcurrencies_data.txt")
data <- flandreau_jobst_internationalcurrencies_data

#creating network object and adjacency matrix
df <- as.data.frame(cbind(flandreau_jobst_internationalcurrencies_data$country_A,flandreau_jobst_internationalcurrencies_data$country_B,flandreau_jobst_internationalcurrencies_data$quote1900))
el <- filter(df,V3==1)
el <- select(el,V1,V2)
net <- network(el, edgelist= T,lab=T)
Y2 <- as.matrix(net)

#create entries for self loops (with all NAs, used for creating covariate matrices)
countries <- unique(data$country_A)
for (i in countries) {
  data <- rbind.data.frame(data,c(i,i,rep(NA,23)))
}

data <- data %>% arrange(country_A,country_B)

#create edge covariate "colony"
colony_mat <- matrix(as.numeric(data$colony),nrow = 45,ncol=45)

#create edge covariate "distance"
dist_mat <- matrix(as.numeric(data$dist),nrow = 45,ncol=45)

#create edge covariate "bitrade"
bitrade_mat <- matrix(as.numeric(data$bitrade),nrow = 45,ncol=45)
bitrade_mat <- bitrade_mat 

#creating node covariates
data_temp <- data[-1,] # just for convenience
gold <- as.numeric(data_temp$gold[match(countries, data_temp$country_A)])
debtburden <- as.numeric(data_temp$debtburden[match(countries, data_temp$country_A)])
rlong <- as.numeric(data_temp$rlong[match(countries, data_temp$country_A)])
rshort <- as.numeric(data_temp$rshort1900[match(countries, data_temp$country_A)])
rgdp <- as.numeric(data_temp$rgdp[match(countries, data_temp$country_A)])
rgdpcap <- as.numeric(data_temp$rgdpcap[match(countries, data_temp$country_A)])
coverage <- as.numeric(data_temp$coverage[match(countries, data_temp$country_A)])
poldemo <- as.numeric(data_temp$poldemo[match(countries, data_temp$country_A)])

#replacing single NA in poldemo with 0: 
poldemo[34] <- 0


#putting edge covariates into an array (made of 2 matrices) to feed to AME
# we don't use "colony" as that is too sparse to estimate 
#(only 12 colonic relationships out of 1880)
Xd <- abind(dist_mat,log1p(bitrade_mat),along=3)
dimnames(Xd)[[3]]<-c("dist","bitrade")


#putting all node covariates into a matrix to feed to AME
#not using debtburden, rlong, rshort as too many NAs 
#  non-significant and creating singularities when imputing NAs to zeros
#(also just too many covariates)

#Xn <- cbind(gold,debtburden,rlong,rshort,rgdp,rgpcap,poldemo,coverage)
Xn <- cbind(gold,rgdpcap,poldemo,coverage)


#fitting the AME model
fitAME<-ame(Y2, Xdyad = Xd, Xcol = Xn, Xrow = Xn, R = 2, family = "bin",
            plot = T, print = F, intercept = T,seed = 123)

#fitting the Classical Probit model
probmod<-ame(Y2, Xdyad = Xd, Xcol = Xn, Xrow = Xn, R = 0, family = "bin",
             plot = T, print = F, intercept = T,seed = 123,rvar = F,cvar = F,
             dcor = F,nvar = F)




