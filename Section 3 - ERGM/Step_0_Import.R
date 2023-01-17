# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# ------                          Trade Of Arms                         ------ #
# ------                          (Import Data)                         ------ #
# ------                          Cornelius Fritz                       ------ #
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #


# This script preprocesses the data contained in the "Data" folder,
# producing the files (already uploaded) in the "RData" folder.
# The script was included for completeness, and it is not necessary for
# reproducing the ERGM modeling as in the paper (given that the
# preprocessed data is already contained in the "RData" folder).

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))  

library(reshape2)
library(sna)
library(foreign)
library(MASS)
library(ggplot2)
library(grid)
library(readstata13)
library(lubridate)
library(readxl)
library(states)
library(imputeTS)
library(data.table)
library(peacesciencer)
ccode_democracy

rm(list=ls(all=TRUE))
years = 2010:2020
exo_data = create_dyadyears(subset_years = years)
exo_data = add_gwcode_to_cow(exo_data)
exo_data = add_contiguity(exo_data)
exo_data = add_capital_distance(exo_data)
exo_data = add_atop_alliance(exo_data)
exo_data = add_igos(exo_data)
exo_data = add_minimum_distance(exo_data)
exo_data = add_nmc(exo_data)
exo_data = add_sdp_gdp(exo_data)
exo_data = add_cow_alliance(exo_data)
exo_data = add_democracy(exo_data)

# 9. Major Conventional Weapons by SIPRI -------

# Data Sources: 
# SIPRI Arms Transfers Database
# https://www.sipri.org/databases/armstransfers
dfl = fread("Data/DealsAndTIVs-2022-08-01-12 28 04.txt",nrows = 11043)
changeCols<- names(dfl)
dfl = dfl[dfl$Seller != ""]
# There seems to be Problem with deals whose description includes in the raw data a ; in the Description 
missing_tmp = which(suppressWarnings(is.na(as.numeric(as.character(dfl$Delivery.year)))))
dfl[missing_tmp,6:16] = dfl[missing_tmp,7:17]

dfl$Delivery.year = dfl$`Delivery year`
dfl$Order.date = dfl$`Order date`
dfl$Seller[dfl$Seller == "Cote d'Ivoire"] = "Cote dIvoire"
# Czechia = Czech Republic
dfl$Seller[dfl$Seller == "Czechia"] = "Czech Republic"
# Macedonia = Macedonia (FYROM)
dfl$Seller[dfl$Seller == "Macedonia"] = "Macedonia (FYROM)"
# East Germany (GDR) = German Democratic Republic
dfl$Seller[dfl$Seller == "East Germany (GDR)"] = "German Democratic Republic"
# Bosnia-Herzegovina = Bosnia and Herzegovina
dfl$Seller[dfl$Seller == "Bosnia-Herzegovina"] = "Bosnia and Herzegovina"
# Cote d'Ivoire = Cote dIvoire

dfl$Buyer[dfl$Buyer == "Cote d'Ivoire"] = "Cote dIvoire"
# Czechia = Czech Republic
dfl$Buyer[dfl$Buyer == "Czechia"] = "Czech Republic"
# Macedonia = Macedonia (FYROM)
dfl$Buyer[dfl$Buyer == "Macedonia"] = "Macedonia (FYROM)"
# East Germany (GDR) = German Democratic Republic
dfl$Buyer[dfl$Buyer == "East Germany (GDR)"] = "German Democratic Republic"
# Bosnia-Herzegovina = Bosnia and Herzegovina
dfl$Buyer[dfl$Buyer == "Bosnia-Herzegovina"] = "Bosnia and Herzegovina"

dfl$Seller[dfl$Seller == "Soviet Union"] = "Russia" # Change all transactions with the Soviet Union to Russia 
dfl$Buyer[dfl$Buyer == "Soviet Union"] = "Russia" 

dfl$TIV.delivery.values = as.numeric(dfl$`TIV delivery values`)
dfl = aggregate(dfl$TIV.delivery.values, by = list(paste(dfl$Seller,dfl$Buyer, dfl$Delivery.year,sep = "_")), sum) # Aggregate all TIV values by country combination and year 
tmp = unlist(strsplit(dfl$Group.1,split = "_")) # 1,4,7 ... are sender 2,5,8 ... are recevier 3,6,9 ... is the year  

dfl$from = tmp[seq(1,to = length(tmp),by = 3)]
dfl$to = tmp[seq(2,to = length(tmp),by = 3)]
dfl$year = as.numeric(tmp[seq(3,to = length(tmp),by = 3)])
country_list = fread("Data/country_list.csv")
dfl$from_id = match(dfl$from,country_list$SIPRI) # Match the country names with the ids in country_list
dfl$to_id = match(dfl$to,country_list$SIPRI)

# Try matching the unmatched names with the alternative name also given in country_list
dfl$from_id[is.na(dfl$from_id)] = match(dfl$from[is.na(dfl$from_id)],country_list$V1)
dfl$to_id[is.na(dfl$to_id)] = match(dfl$to[is.na(dfl$to_id)],country_list$V1)

unique(dfl$to[is.na(dfl$to_id)]) # Rest of the unmatched countries are not independent states and thus not part of the analysis 
unique(dfl$from[is.na(dfl$from_id)]) # Rest of the unmatched countries are not independent states and thus not part of the analysis 

dfl = dfl[!(is.na(dfl$from_id) | is.na(dfl$to_id)),]
setDT(dfl)
unique_actors = unique(c(dfl$from,dfl$to))

dfl$from_id = match(dfl$from,unique_actors)
dfl$to_id = match(dfl$to,unique_actors)
amk<- list()
# create a list of 69 adjacency matrices, one for each year from 1950:2017. amk[[1]] = 1950,...m amk[[67]]=2017

for (i in 1:length(years)){
  amk[[i]]<- matrix(0,length(unique_actors),length(unique_actors))
  colnames(amk[[i]])<-(unique_actors)
  rownames(amk[[i]])<-(unique_actors)
}

# Fill the matrix per year 
for(i in years){
  tmp_dfl = dfl[year == i,]
  amk[[i- min(years) +1]][cbind(tmp_dfl$from_id, tmp_dfl$to_id)] = 1
}


# todo: only study dyads where all exogenous info is available 
actors = unique(exo_data$gwcode1)
data(gwstates)
gwstates$start = ymd(gwstates$start)
gwstates$end = ymd(gwstates$end)
gwstates$start_year = year(gwstates$start)
gwstates$end_year = year(gwstates$end)
setDT(gwstates)
setDT(exo_data)
gwstates = gwstates[(start_year<= min(years)) & (end_year>= max(years))]

unique_actors_gw = unique_actors
unique_actors_gw = gwstates$gwcode[match(unique_actors,gwstates$country_name)]
unique_actors[is.na(match(unique_actors,gwstates$country_name))]
unique_actors[unique_actors == "Belarus"] = "Belarus (Byelorussia)"
unique_actors[unique_actors == "Bosnia and Herzegovina"] = "Bosnia-Herzegovina"
unique_actors[unique_actors == "Cote dIvoire"] = "Cote D'Ivoire"
unique_actors[unique_actors == "Germany"] = "German Federal Republic"
unique_actors[unique_actors == "Iran"] = "Iran (Persia)"
unique_actors[unique_actors == "Italy"] = "Italy/Sardinia"
unique_actors[unique_actors == "Kyrgyzstan"] = "Kyrgyz Republic"
unique_actors[unique_actors == "North Korea"] = "Korea, People's Republic of"
unique_actors[unique_actors == "Russia"] = "Russia (Soviet Union)"
unique_actors[unique_actors == "South Korea"] = "Korea, Republic of"
unique_actors[unique_actors == "Turkey"] = "Turkey (Ottoman Empire)"
unique_actors[unique_actors == "UAE"] = "United Arab Emirates"
unique_actors[unique_actors == "United States"] = "United States of America"
unique_actors[unique_actors == "Viet Nam"] = "Vietnam, Democratic Republic of"
unique_actors[unique_actors == "Burkina Faso"] = "Burkina Faso (Upper Volta)"
unique_actors[unique_actors == "Samoa"] = "Samoa/Western Samoa"
unique_actors[unique_actors == "Sri Lanka"] = "Sri Lanka (Ceylon)"
unique_actors[unique_actors == "Myanmar"] = "Myanmar (Burma)"
unique_actors[unique_actors == "Yemen"] = "Yemen (Arab Republic of Yemen)"
unique_actors[unique_actors == "Cambodia"] = "Cambodia (Kampuchea)"
unique_actors[unique_actors == "DR Congo"] = "Congo, Democratic Republic of (Zaire)"
unique_actors[unique_actors == "Tanzania"] = "Tanzania/Tanganyika"
unique_actors[unique_actors == "Timor-Leste"] = "East Timor"
unique_actors[unique_actors == "Suriname"] = "Surinam"
unique_actors_gw = gwstates$gwcode[match(unique_actors,gwstates$country_name)]
# We will only incorporate data from countries where the information is available for all years 
excluded = logical(length = length(unique_actors_gw))
excluded[is.na(unique_actors_gw)] = TRUE

exo_data = exo_data[(gwcode1 %in% unique_actors_gw) & (gwcode2 %in% unique_actors_gw)]
actors = actors[actors   %in%  exo_data$gwcode1]
gwstates = gwstates[gwcode %in% actors]
exo_data$id_1 = match(exo_data$gwcode1,unique_actors_gw)
exo_data$id_2 = match(exo_data$gwcode2,unique_actors_gw)
unique_actors_iso = gwstates$gwc[match(unique_actors_gw,gwstates$gwcode)]

# Save a list of matrices of the mcw in each year with ISO3 names instead of the full country names
iso3<- list()
# create a list of 69 adjacency matrices, one for each year from 1950:2017. amk[[1]] = 1950,...m amk[[67]]=2017

for (i in 1:length(years)){
  iso3[[i]]<- matrix(0,length(unique_actors),length(unique_actors))
  colnames(iso3[[i]])<-(unique_actors_iso)
  rownames(iso3[[i]])<-(unique_actors_iso)
}

# Fill the matrix per year 
for(i in years){
  tmp_dfl = dfl[year == i,]
  iso3[[i- min(years) +1]][cbind(tmp_dfl$from_id, tmp_dfl$to_id)] = 1
}




# Alliances 
alliance = list()
for (i in 1:length(years)){
  alliance[[i]]<- matrix(0,length(unique_actors),length(unique_actors))
  colnames(alliance[[i]])<-(unique_actors)
  rownames(alliance[[i]])<-(unique_actors)
}

# Fill the matrix per year 
exo_data$atop_any = exo_data$atop_consul + exo_data$atop_defense+ exo_data$atop_offense+ exo_data$atop_neutral+ exo_data$atop_nonagg
exo_data$atop_any = exo_data$atop_any>0
for(i in years){
  tmp_alliance = exo_data[(year == i) & (atop_defense),]
  alliance[[i- min(years) +1]][cbind(tmp_alliance$id_1, tmp_alliance$id_2)] = 1
}


# Min distance between countries 
mindist =  matrix(NA,length(unique_actors),length(unique_actors))
colnames(mindist)<-(unique_actors)
rownames(mindist)<-(unique_actors)
mindist[cbind(exo_data$id_1,exo_data$id_2)] = exo_data$mindist
mindist[(!is.na(unique_actors_gw)),(!is.na(unique_actors_gw))]

# Distance between capitals of countries 
capdist =  matrix(NA,length(unique_actors),length(unique_actors))
colnames(capdist)<-(unique_actors)
rownames(capdist)<-(unique_actors)
capdist[cbind(exo_data$id_1,exo_data$id_2)] = exo_data$capdist

# direct contiguity between countries 
contiguity =  matrix(NA,length(unique_actors),length(unique_actors))
colnames(contiguity)<-(unique_actors)
rownames(contiguity)<-(unique_actors)
contiguity[cbind(exo_data$id_1,exo_data$id_2)] = exo_data$conttype
contiguity[contiguity != 1] = 0


# Log GDP
gdp_data = matrix(ncol = length(years),nrow = length(unique_actors_gw))
for(tmp_year in years){
  data_tmp = exo_data[year == tmp_year] 
  log_gdp = data_tmp$wbgdp2011est1[match(unique_actors_gw,data_tmp$gwcode1)]
  gdp_data[,tmp_year-min(years)+1] = log_gdp
}

# Log GDPPC
gdppc_data = matrix(ncol = length(years),nrow = length(unique_actors_gw))
for(tmp_year in years){
  data_tmp = exo_data[year == tmp_year] 
  log_gdppc = data_tmp$wbgdppc2011est1[match(unique_actors_gw,data_tmp$gwcode1)]
  gdppc_data[,tmp_year-min(years)+1] = log_gdppc
}


# Polity
polity = matrix(ncol = length(years),nrow = length(unique_actors_gw))
for(tmp_year in years){
  data_tmp = exo_data[year == tmp_year] 
  polity_tmp = exo_data$polity21[match(unique_actors_gw,exo_data$gwcode1)]
  # polity_tmp = data_tmp$v2x_polyarchy1[match(unique_actors_gw,data_tmp$gwcode1)]
  polity[,tmp_year-min(years)+1] = polity_tmp
}

# Population
population = matrix(ncol = length(years),nrow = length(unique_actors_gw))
for(tmp_year in years){
  data_tmp = exo_data[year == tmp_year] 
  # polity_tmp = exo_data$polity21[match(unique_actors_gw,exo_data$gwcode1)]
  population_tmp = data_tmp$wbpopest1[match(unique_actors_gw,data_tmp$gwcode1)]
  population[,tmp_year-min(years)+1] = population_tmp
}

# Military Expenditure 
milexp = matrix(ncol = length(years),nrow = length(unique_actors_gw))
for(tmp_year in years){
  data_tmp = exo_data[year == tmp_year] 
  # polity_tmp = exo_data$polity21[match(unique_actors_gw,exo_data$gwcode1)]
  milexp_tmp = data_tmp$milex1[match(unique_actors_gw,data_tmp$gwcode1)]
  milexp[,tmp_year-min(years)+1] = log(milexp_tmp)
}
# excluded[rowSums(milexp,na.rm = T) == 0] = TRUE

# Cinc
cinc = matrix(ncol = length(years),nrow = length(unique_actors_gw))
for(tmp_year in years){
  data_tmp = exo_data[year == tmp_year] 
  # polity_tmp = exo_data$polity21[match(unique_actors_gw,exo_data$gwcode1)]
  cinc_tmp = data_tmp$cinc1[match(unique_actors_gw,data_tmp$gwcode1)]
  cinc[,tmp_year-min(years)+1] = log(cinc_tmp)
}

excluded[is.na(population[,6])] = TRUE
excluded[is.na(polity[,6])] = TRUE
excluded[is.na(gdppc_data[,6])] = TRUE
# excluded[rowSums(mindist,na.rm = T) == 0] = TRUE
excluded[rowSums(capdist,na.rm = T) == 0] = TRUE
excluded[is.na(contiguity[,6])] = TRUE
excluded[is.na(gdp_data[,6])] = TRUE
excluded[is.na(cinc[,6])] = TRUE

cinc = cinc[!excluded,]
milexp = milexp[!excluded,]
population = population[!excluded,]
polity = polity[!excluded,]
gdp_data = gdp_data[!excluded,]
gdppc_data = gdppc_data[!excluded,]
mindist = mindist[!excluded,!excluded]
capdist = capdist[!excluded,!excluded]
contiguity= contiguity[!excluded,!excluded]

mcw = amk
for(i in 1:length(years)){
  mcw[[i]] = amk[[i]][!excluded,!excluded]
}

mcw_iso3 = iso3
for(i in 1:length(years)){
  mcw_iso3[[i]] = iso3[[i]][!excluded,!excluded]
}

for(i in 1:length(years)){
  alliance[[i]] = alliance[[i]][!excluded,!excluded]
}
# 




# excluded[rowSums(population,na.rm = T) == 0] = TRUE
# excluded[rowSums(polity,na.rm = T) == 0] = TRUE
# excluded[rowSums(gdppc_data,na.rm = T) == 0] = TRUE
# excluded[rowSums(mindist,na.rm = T) == 0] = TRUE
# excluded[rowSums(capdist,na.rm = T) == 0] = TRUE
# excluded[rowSums(contiguity,na.rm = T) == 0] = TRUE
# excluded[rowSums(gdp_data,na.rm = T) == 0] = TRUE
# excluded[rowSums(cinc,na.rm = T) == 0] = TRUE


save(capdist,file = "RData/capdist.RData")
save(contiguity,file = "RData/contiguity.RData")
save(alliance,file = "RData/alliance.RData")
save(cinc,file = "RData/cinc.RData")
save(mcw,file = "RData/mcw.RData")
save(mcw_iso3,file = "RData/mcw_iso3.RData")
save(milexp,file = "RData/milexp.RData")
save(population,file = "RData/population.RData")
save(polity,file = "RData/polity.RData")
save(gdp_data,file = "RData/gdp_data.RData")
save(gdppc_data,file = "RData/gdppc_data.RData")
save(mindist,file = "RData/mindist.RData")

# rownames(gdp_data) = unique_actors_gw[!excluded]
# colnames(gdp_data) = years
# exo_data[(year == 2015) & (gwcode1 == 712),wbgdp2011est1][1]
# gdp_data


