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


# 1. Analysis with four year lag ----

# Response network ----
year = 8
cov_year = year-4
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

# Estimation ----
nsim=5000

model_4 = ergm(net~ edges + 
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


gof_model_4 = gof(model_4,control = control.gof.ergm(seed = 123))
pdf("Plots/gof_model_4.pdf")
plot(gof_model_4)
dev.off()
pdf("Plots/mcmc_diag_model_4.pdf")
mcmc.diagnostics(model_4)
dev.off()

# Model Comparison
model_logit_4 = ergm(net~ edges + 
                     edgecov(lnet)+
                     edgecov(capdist) +
                     edgecov(poldiff) + 
                     edgecov(alliance_tmp) +
                     nodeocov("lgdp")+
                     nodeicov("lgdp") )

model_matrix_logit_4 = ergmMPLE(net~ edges + 
                                edgecov(lnet)+
                                edgecov(capdist) +
                                edgecov(poldiff) + 
                                edgecov(alliance_tmp) +
                                nodeocov("lgdp")+
                                nodeicov("lgdp") )
logit_4=model_matrix_logit_4$predictor%*%coef(model_logit_4)
p_hat_logit_4=1/(1+exp(-logit_4))

tmp_simulations_4 = simulate(model_4, nsim = 10000)
simulation_matrix_4 = lapply(tmp_simulations_4, as.matrix)
p_hat_model_4 = Reduce(simulation_matrix_4, f = "+")/10000
diag(p_hat_model_4) = NA
roc_logit_4 <-roc.curve(scores.class0 = p_hat_logit_4,weights.class0 = model_matrix_logit_4$response,curve=T)
pr_logit_4<-pr.curve(scores.class0 = p_hat_logit_4,weights.class0 = model_matrix_logit_4$response,curve=T)

p_hat_model_4 = as.vector(p_hat_model_4)
p_hat_model_4 = p_hat_model_4[!is.na(p_hat_model_4)]
response_4 =as.matrix(net) 
diag(response_4) = NA
response_4 = as.vector(response_4)
response_4 = response_4[!is.na(response_4)]
roc_ergm_4 <-roc.curve(scores.class0 =p_hat_model_4,weights.class0 = response_4,curve=T)
pr_ergm_4 <-pr.curve(scores.class0 =p_hat_model_4,weights.class0 = response_4,curve=T)

# You can also plot the pr and roc curves with the provided generic functions of the PRROC package
# plot(pr_ergm_4)
# plot(roc_ergm_4)
# 
# plot(pr_logit_4)
# plot(roc_logit_4)

df_roc_ergm_4 = data.frame("specificity" =1- roc_ergm_4$curve[,1],
                          "sensitivity" = roc_ergm_4$curve[,2], "Model" = "ERGM")
df_roc_logit_4 = data.frame("specificity" =1- roc_logit_4$curve[,1],
                           "sensitivity" = roc_logit_4$curve[,2], "Model" = "Logit")
df_roc_4= rbind(df_roc_ergm_4, df_roc_logit_4)

pdf("Plots/roc_comparison_4.pdf",width = 9,height =9)
ggplot(data = df_roc_4 ,mapping = aes(x = specificity,y = sensitivity, lty = Model)) +
  geom_step(alpha = 0.8, size=1) +
  scale_x_reverse() +
  geom_abline(slope = 1, intercept = 1, lty = 2, color = "grey40") +
  ylab("Sensitivity")+
  xlab("Specificity") +
  theme_pubr() +
  scale_linetype_manual(name = "", 
                        labels = c(paste0("ERGM (",round(roc_ergm_4$auc,digits = 3) ,")"),
                                   paste0("Logit (",round(roc_logit_4$auc,digits = 3) ,")")), 
                        values = c(1,4,2),guide = guide_legend(keywidth = 4)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.title.x = element_text(size = 30)) +
  theme(axis.title.y = element_text(size = 30)) +
  theme(axis.text.x = element_text(size = 20))  +
  theme(axis.text.y = element_text(size = 20)) +
  theme(plot.title = element_text(size = 30), 
        legend.text = element_text(size = 15) )

dev.off()


df_pr_ergm_4 = data.frame("Precision" = pr_ergm_4$curve[,1],
                        "Recall" = pr_ergm_4$curve[,2], "Model" = "ERGM")
df_pr_logit_4 = data.frame("Precision" =pr_logit_4$curve[,1],
                         "Recall" = pr_logit_4$curve[,2], "Model" = "Logit")
df_pr_4= rbind(df_pr_ergm_4,df_pr_logit_4)

pdf("Plots/pr_comparison_4.pdf",width = 9,height =9)
ggplot(data = df_pr_4 ,mapping = aes(x = Precision,y = Recall, lty = Model)) +
  geom_step(alpha = 0.8, size=1) +
  ylab("Precision")+
  xlab("Recall") +
  theme_pubr() +
  scale_linetype_manual(name = "", 
                        labels = c(paste0("ERGM (",round(pr_ergm_4$auc.integral,digits = 3) ,")"),
                                   paste0("Logit (",round(pr_logit_4$auc.integral,digits = 3) ,")")), 
                        values = c(1,2),guide = guide_legend(keywidth = 4)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.title.x = element_text(size = 30)) +
  theme(axis.title.y = element_text(size = 30)) +
  theme(axis.text.x = element_text(size = 20))  +
  theme(axis.text.y = element_text(size = 20)) +
  theme(plot.title = element_text(size = 30), 
        legend.text = element_text(size = 15) )

dev.off()
# Get the latex table 
texreg(list(model_4,model_logit_4),single.row = T,digits = 3)


# 2. Analysis with five year lag ----

# Response network ----
year = 8
cov_year = year-5
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

# Estimation ----
nsim=5000

model_5 = ergm(net~ edges + 
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


gof_model_5 = gof(model_5,control = control.gof.ergm(seed = 123))
pdf("Plots/gof_model_5.pdf")
plot(gof_model_5)
dev.off()
pdf("Plots/mcmc_diag_model_5.pdf")
mcmc.diagnostics(model_5)
dev.off()

# Model Comparison
model_logit_5 = ergm(net~ edges + 
                       edgecov(lnet)+
                       edgecov(capdist) +
                       edgecov(poldiff) + 
                       edgecov(alliance_tmp) +
                       nodeocov("lgdp")+
                       nodeicov("lgdp") )

model_matrix_logit_5 = ergmMPLE(net~ edges + 
                                  edgecov(lnet)+
                                  edgecov(capdist) +
                                  edgecov(poldiff) + 
                                  edgecov(alliance_tmp) +
                                  nodeocov("lgdp")+
                                  nodeicov("lgdp") )
logit_5=model_matrix_logit_5$predictor%*%coef(model_logit_5)
p_hat_logit_5=1/(1+exp(-logit_5))

tmp_simulations_5 = simulate(model_5, nsim = 10000)
simulation_matrix_5 = lapply(tmp_simulations_5, as.matrix)
p_hat_model_5 = Reduce(simulation_matrix_5, f = "+")/10000
diag(p_hat_model_5) = NA
roc_logit_5 <-roc.curve(scores.class0 = p_hat_logit_5,weights.class0 = model_matrix_logit_5$response,curve=T)
pr_logit_5<-pr.curve(scores.class0 = p_hat_logit_5,weights.class0 = model_matrix_logit_5$response,curve=T)

p_hat_model_5 = as.vector(p_hat_model_5)
p_hat_model_5 = p_hat_model_5[!is.na(p_hat_model_5)]
response_5 =as.matrix(net) 
diag(response_5) = NA
response_5 = as.vector(response_5)
response_5 = response_5[!is.na(response_5)]
roc_ergm_5 <-roc.curve(scores.class0 =p_hat_model_5,weights.class0 = response_5,curve=T)
pr_ergm_5 <-pr.curve(scores.class0 =p_hat_model_5,weights.class0 = response_5,curve=T)

# You can also plot the pr and roc curves with the provided generic functions of the PRROC package
# plot(pr_ergm_5)
# plot(roc_ergm_5)
# 
# plot(pr_logit_5)
# plot(roc_logit_5)

df_roc_ergm_5 = data.frame("specificity" =1- roc_ergm_5$curve[,1],
                           "sensitivity" = roc_ergm_5$curve[,2], "Model" = "ERGM")
df_roc_logit_5 = data.frame("specificity" =1- roc_logit_5$curve[,1],
                            "sensitivity" = roc_logit_5$curve[,2], "Model" = "Logit")
df_roc_5= rbind(df_roc_ergm_5, df_roc_logit_5)

pdf("Plots/roc_comparison_5.pdf",width = 9,height =9)
ggplot(data = df_roc_5 ,mapping = aes(x = specificity,y = sensitivity, lty = Model)) +
  geom_step(alpha = 0.8, size=1) +
  scale_x_reverse() +
  geom_abline(slope = 1, intercept = 1, lty = 2, color = "grey40") +
  ylab("Sensitivity")+
  xlab("Specificity") +
  theme_pubr() +
  scale_linetype_manual(name = "", 
                        labels = c(paste0("ERGM (",round(roc_ergm_5$auc,digits = 3) ,")"),
                                   paste0("Logit (",round(roc_logit_5$auc,digits = 3) ,")")), 
                        values = c(1,4,2),guide = guide_legend(keywidth = 4)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.title.x = element_text(size = 30)) +
  theme(axis.title.y = element_text(size = 30)) +
  theme(axis.text.x = element_text(size = 20))  +
  theme(axis.text.y = element_text(size = 20)) +
  theme(plot.title = element_text(size = 30), 
        legend.text = element_text(size = 15) )

dev.off()


df_pr_ergm_5 = data.frame("Precision" = pr_ergm_5$curve[,1],
                          "Recall" = pr_ergm_5$curve[,2], "Model" = "ERGM")
df_pr_logit_5 = data.frame("Precision" =pr_logit_5$curve[,1],
                           "Recall" = pr_logit_5$curve[,2], "Model" = "Logit")
df_pr_5= rbind(df_pr_ergm_5,df_pr_logit_5)

pdf("Plots/pr_comparison_5.pdf",width = 9,height =9)
ggplot(data = df_pr_5 ,mapping = aes(x = Precision,y = Recall, lty = Model)) +
  geom_step(alpha = 0.8, size=1) +
  ylab("Precision")+
  xlab("Recall") +
  theme_pubr() +
  scale_linetype_manual(name = "", 
                        labels = c(paste0("ERGM (",round(pr_ergm_5$auc.integral,digits = 3) ,")"),
                                   paste0("Logit (",round(pr_logit_5$auc.integral,digits = 3) ,")")), 
                        values = c(1,2),guide = guide_legend(keywidth = 4)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.title.x = element_text(size = 30)) +
  theme(axis.title.y = element_text(size = 30)) +
  theme(axis.text.x = element_text(size = 20))  +
  theme(axis.text.y = element_text(size = 20)) +
  theme(plot.title = element_text(size = 30), 
        legend.text = element_text(size = 15) )

dev.off()
# Get the latex table 
texreg(list(model_4,model_logit_4,model_5,model_logit_5),single.row = T,digits = 2,stars = 0,
       custom.coef.names = c("Edges", "Repetition", "Distance", "Abs. Diff. Polity",
                             "Alliance", "log-GDP (Sender)", "log-GDP (Receiver)", "Mutual", 
                             "GWIDEG", "GWODEG", "GWOTP", "GWISP"),
       custom.model.names = c("ERGM", "Logit","ERGM", "Logit"))
LR_statistic = -2*(logLik(model_4) -logLik(model_logit_4))
n_df = length(model_4$coefficients) - length(model_logit_4$coefficients)
qchisq(0.025,df = n_df)
