library(statnet)
library(texreg)
library(PRROC)
library(ggpubr)

rm(list=ls())


setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 
load("RData/mcw.RData")
load("RData/cinc.RData")
load("RData/gdp_data.RData")
load("RData/gdppc_data.RData")
load("RData/polity.RData")
load("RData/population.RData")
load("RData/milexp.RData")
load("RData/mindist.RData")
load("RData/capdist.RData")
load("RData/alliance.RData")
load("RData/contiguity.RData")

# Response network ----
year = 8
cov_year = year-3
# 2 would
net<-network(mcw[[year]])

write.csv(x = mcw[[year]], file = "mcw.csv")

# Covariates----

# monadic attributes 
lgdp<-gdp_data[,cov_year]
tmp_cinc<-cinc[,cov_year]
tmp_population<-population[,cov_year]

# dyadic attributes 
poldiff<-outer(polity[,cov_year],polity[,cov_year],FUN=function(x,y) abs(x-y))
alliance_tmp<-alliance[[cov_year]]
capdist

# Lagged network
lag_net = mcw[[year-1]] + mcw[[year-2]] + mcw[[year-3]]
lag_net[lag_net>0] = 1
lnet<-network(lag_net)
lrecip<-network(t(lnet))

# set attributes
net <- set.vertex.attribute(net, "lgdp", lgdp)
net <- set.vertex.attribute(net, "cinc", tmp_cinc)
net <- set.vertex.attribute(net, "pop", tmp_population)
diag(mindist) = 0
diag(capdist) = 0

mindist = log1p(mindist)
capdist = log1p(capdist)

# Combine them in a list
nets<-network.list(list(lnet,net))

# Estimation ----
nsim=5000

model = ergm(net~ edges + 
               edgecov(lnet)+
               edgecov(capdist) +
               edgecov(poldiff) + 
               edgecov(alliance_tmp) +
               nodeocov("lgdp")+
               nodeicov("lgdp") +
               mutual + 
               gwidegree(log(2),fixed=T)+
               gwodegree(log(2),fixed=T) + 
               dgwesp(log(2),fixed=T, type = "OTP") + 
               dgwesp(log(2),fixed=T, type = "ISP") ,control = control.ergm(seed = 123))

gof_model = gof(model,control = control.gof.ergm(seed = 123))
pdf("Plots/gof_model.pdf")
plot(gof_model)
dev.off()
pdf("Plots/mcmc_diag_model.pdf")
mcmc.diagnostics(model)
dev.off()

# Model Comparison
model_logit = ergm(net~ edges + 
                     edgecov(lnet)+
                     edgecov(capdist) +
                     edgecov(poldiff) + 
                     edgecov(alliance_tmp) +
                     nodeocov("lgdp")+
                     nodeicov("lgdp") )

model_matrix_logit = ergmMPLE(net~ edges + 
                                edgecov(lnet)+
                                edgecov(capdist) +
                                edgecov(poldiff) + 
                                edgecov(alliance_tmp) +
                                nodeocov("lgdp")+
                                nodeicov("lgdp") )
logit=model_matrix_logit$predictor%*%coef(model_logit)
p_hat_logit=1/(1+exp(-logit))

tmp_simulations = simulate(model, nsim = 10000)
simulation_matrix = lapply(tmp_simulations, as.matrix)
p_hat_model = Reduce(simulation_matrix, f = "+")/10000
diag(p_hat_model) = NA
roc_logit<-roc.curve(scores.class0 = p_hat_logit,weights.class0 = model_matrix_logit$response,curve=T)
pr_logit<-pr.curve(scores.class0 = p_hat_logit,weights.class0 = model_matrix_logit$response,curve=T)

p_hat_model = as.vector(p_hat_model)
p_hat_model = p_hat_model[!is.na(p_hat_model)]
response =as.matrix(net) 
diag(response) = NA
response = as.vector(response)
response = response[!is.na(response)]
roc_ergm <-roc.curve(scores.class0 =p_hat_model,weights.class0 = response,curve=T)
pr_ergm <-pr.curve(scores.class0 =p_hat_model,weights.class0 = response,curve=T)

plot(pr_ergm)
plot(roc_ergm)

plot(pr_logit)
plot(roc_logit)

df_roc_ergm = data.frame("specificity" =1- roc_ergm$curve[,1],
                          "sensitivity" = roc_ergm$curve[,2], "Model" = "ERGM")
df_roc_logit = data.frame("specificity" =1- roc_logit$curve[,1],
                           "sensitivity" = roc_logit$curve[,2], "Model" = "Logit")
df_roc= rbind(df_roc_ergm, df_roc_logit)

pdf("Plots/roc_comparison.pdf",width = 9,height =9)
ggplot(data = df_roc ,mapping = aes(x = specificity,y = sensitivity, lty = Model)) +
  geom_step(alpha = 0.8, size=1) +
  scale_x_reverse() +
  geom_abline(slope = 1, intercept = 1, lty = 2, color = "grey40") +
  ylab("Sensitivity")+
  xlab("Specificity") +
  theme_pubr() +
  scale_linetype_manual(name = "", 
                        labels = c(paste0("ERGM (",round(roc_ergm$auc,digits = 3) ,")"),
                                   paste0("Logit (",round(roc_logit$auc,digits = 3) ,")")), 
                        values = c(1,4,2),guide = guide_legend(keywidth = 4)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.title.x = element_text(size = 30)) +
  theme(axis.title.y = element_text(size = 30)) +
  theme(axis.text.x = element_text(size = 20))  +
  theme(axis.text.y = element_text(size = 20)) +
  theme(plot.title = element_text(size = 30), 
        legend.text = element_text(size = 15) )

dev.off()


df_pr_ergm = data.frame("Precision" = pr_ergm$curve[,1],
                        "Recall" = pr_ergm$curve[,2], "Model" = "ERGM")
df_pr_logit = data.frame("Precision" =pr_logit$curve[,1],
                         "Recall" = pr_logit$curve[,2], "Model" = "Logit")
df_pr= rbind(df_pr_ergm,df_pr_logit)

pdf("Plots/pr_comparison.pdf",width = 9,height =9)
ggplot(data = df_pr ,mapping = aes(x = Precision,y = Recall, lty = Model)) +
  geom_step(alpha = 0.8, size=1) +
  ylab("Precision")+
  xlab("Recall") +
  theme_pubr() +
  scale_linetype_manual(name = "", 
                        labels = c(paste0("ERGM (",round(pr_ergm$auc.integral,digits = 3) ,")"),
                                   paste0("Logit (",round(pr_logit$auc.integral,digits = 3) ,")")), 
                        values = c(1,2),guide = guide_legend(keywidth = 4)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.title.x = element_text(size = 30)) +
  theme(axis.title.y = element_text(size = 30)) +
  theme(axis.text.x = element_text(size = 20))  +
  theme(axis.text.y = element_text(size = 20)) +
  theme(plot.title = element_text(size = 30), 
        legend.text = element_text(size = 15) )


dev.off()


texreg(list(model,model_logit),single.row = T,digits = 3)
