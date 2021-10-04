##Including humans, cattle, sheep and goats

library(deSolve)
library(MASS)

mat<-read.table("~/Desktop/Hybrid_schisto_pnas_code_datafiles/param_matrix_hybrid") #Parameter matrix
inits<-read.table("~/Desktop/Hybrid_schisto_pnas_code_datafiles/init_matrix_hybrid") #Inits matrix

tx_cover<-rep(0.5,1000)#vary these in simulations
tx_eff<-rep(0.95,1000)

mat<-as.matrix(mat)
inits<-as.matrix(inits)
mat<-cbind(mat,tx_cover,tx_eff) 


#### MATING FUNCTION #####

#Mating function as per May Woolhouse1993, allowing for unequal numbers of males and females. Probability that a female of species i will form a pair with a male of species j

mprob_iF_jM <- function(i,j,parameters,var,k_i){
  with(as.list(c(i,j,parameters,var)),{
    if(i>=0 & j>=0){
      k=k_i
      m=i+j
      q=ifelse(m==0,0,i/m)
      p=ifelse(m==0,0,j/m) ### Proportion male
      J=ifelse(q<0.5,1,(1-q)/q)
      al= m/(m+k)  ### Alpha as defined by equation 2 in May and Woolhouse
      ga = 2*((p*q)^0.5)       ### Gamma as defined by equation 11
      
      integrand <- function(theta)   ### Theta no definition here is just for solving the integral
      {(sin(theta)*sin(theta))/(((1+(al*ga*cos(theta)))^(1+k))*(1+(ga*cos(theta))))}
      intsol <- integrate(integrand, lower=0, upper=pi)
      
      int <- intsol$value
      
      I= (((1-al)^(1+k))/pi)*(int)
      
      prob = J- (2*(1-q)*I)
    }
    
    else(prob=0)
    
    return(prob)
    
  })
}

##### DENSITY-DEPENDENT FUNCTION ###

#Mean density-dependent reduction in fecundity per mated pair- exponential relationship. This is for human part of model
dd_exp <- function(parameters,var,Mean_F_H) 
{
  with(as.list(c(parameters,var,Mean_F_H)),{
    x <- 0:5000
    z <- exp(-b_H)
    M<-Mean_F_H
    p <- dnbinom(x = x, size = k_H, mu = M)
    dd <- z^((x/q_H)-1)
    f <- sum(p*dd*x)/M
    return(f)
  })
}

#Mean density-dependent reduction in fecundity per mated pair- exponential relationship, This is for livestock part of model
dd_exp_i <- function(parameters,var,Mean_F_i,b_i,k_i) ### 
{
  with(as.list(c(parameters,var,Mean_F_i,b_i,k_i)),{
    if(Mean_F_i>0){
      k<-k_i
      b<-b_i
      M<-Mean_F_i
      x <- 0:5000
      z <- exp(b)
      p <- dnbinom(x = x, size = k, mu = M)
      dd <- z^((x/q_Sb)-1)
      f <- sum(p*dd*x)/M
    }
    else(f=1)
    return(f)
  })
}


######################################## TRANSMISSION MODEL ###########################


MWB.dyn <- function(t, var, parameters){
  with(as.list(c(var,parameters)),{
    
    #1. ##SPECIFY VARIABLES TO BE TRACKED HERE 
    
    #HUMAN POPULATION#
    MW_H_Sh <- var[1] #Mean worm burden Sh per person etc.
    
    MW_H_Hyb<- var[2]
    
    MWH_H_F1a<-var[3]
    
    MW_H_F1b<- var[4] 
    
    MW_H_Sb<- var[5]
    
    #LIVESTOCK# 
    
    MW_C_Sb<-var[6] #Mean worm burden Sb per animal etc.
    
    MW_G_Sb<-var[7]
    
    MW_S_Sb<-var[8]
    
    #SNAILS/LARVAL POOL#  
    
    L_Sh <- var[9]   #Larval pool Sh etc
    
    L_Hyb<-var[10]
    
    L_Hyb2_F1<-var[11]
    
    L_Hyb3_F1 <-var[12]
    
    L_Sb <-var[13]
    
    MW_H_TOT<-MW_H_Sh + MW_H_Hyb +MW_H_F1a + MW_H_F1b +MW_H_Sb
    
    
    #2. ##Calculate mean females, males and pairs of each genotype######################### 
    
    #HUMANS 
    #Total Males and females of each species in humans
    
    F_Sh<-MW_H_Sh*q_H_sh
    M_Sh<-MW_H_Sh*(1-q_H_sh)
    
    F_Hyb<-MW_H_Hyb*q_H
    M_Hyb<-MW_H_Hyb*(1-q_H)
    
    F_F1a<-MW_H_F1a*q_H
    M_F1a<-MW_H_F1a*(1-q_H)
    
    F_F1b<-MW_H_F1b*(q_H)
    M_F1b<-MW_H_F1b*(1-q_H)
    
    F_Sb<-MW_H_Sb*q_H
    M_Sb<-MW_H_Sb*(1-q_H)
    
    Mean_F_H<- F_Sh+F_Hyb+F_F1a+F_F1b+F_Sb # This mean total females in human population used in the dd function
    
    #Call mating function sequentially: 7 possible pairing combinations in humans
    
    mprob_Hyb3_Back<-mprob_iF_jM(i=F_Sh, j=M_F1b, parameters, var,k_i=k_H)  #Female Sh, Male F1b, 
    mean_Hyb3_Back_pair<-F_Sh*mprob_Hyb3_Back                        #Mean pairs type 5
    F_Sh_nbackcross<-F_Sh-mean_Hyb3_Back_pair                        #From this calculate how many Female Sh available to mate with Sb
    
    mprob_Hyb3_F1<-mprob_iF_jM(i=F_Sh_nbackcross, j=M_Sb,parameters,var,k_i=k_H)           #Female Sh, Male Sb
    mean_Hyb3_F1_pair<-F_Sh_nbackcross*mprob_Hyb3_F1                #Mean pairs type 7
    F_Sh_spare<- F_Sh_nbackcross- mean_Hyb3_F1_pair                 #From this calculate how many Female Sh available to mate with Sh
    
    mprob_H_Sh_Sh <- mprob_iF_jM(i=F_Sh_spare, j=M_Sh, parameters, var,k_i=k_H)  ###Female Sh Male Sh
    mean_Sh_pair <- F_Sh_spare*mprob_H_Sh_Sh                              #Mean pairs type 1
    M_Sh_nSh<- M_Sh- mean_Sh_pair                                         ### From this calculate mean number of Sh M in humans not mated with Sh Females
    
    mprob_H_Hyb_Sh<-mprob_iF_jM(i=F_Hyb, j=M_Sh_nSh,parameters, var,k_i=k_H)                    #Female Hyb Male Sh
    mean_Hyb_Sh_pair<-F_Hyb*mprob_H_Hyb_Sh                              #Mean pairs type 2
    F_Hyb_nSh<- F_Hyb-mean_Hyb_Sh_pair                                  #From this calculate how many Female Hyb not mated with Sh Males
    M_Sh_nSh_nHyb <- M_Sh_nSh-mean_Hyb_Sh_pair                          #And how many Male Sh not mated with Sh or Hyb
    
    mprob_H_Hyb_Hyb<-mprob_iF_jM(i=F_Hyb_nSh, j=M_Hyb,parameters, var,k_i=k_H)                  #Female Hyb Male Hyb
    mean_Hyb_Hyb_pair<-F_Hyb_nSh*mprob_H_Hyb_Hyb                        #Mean pairs type 3
    
    mprob_Hyb2_Back<-mprob_iF_jM(i=F_F1a, j=M_Sh_nSh_nHyb,parameters, var,k_i=k_H)                  #Female F1a and Male Sh
    mean_Hyb2_Back_pair<-F_F1a*mprob_Hyb2_Back                          #Mean pairs type 4
    M_Sh_Spare<-M_Sh_nSh_nHyb-mean_Hyb2_Back_pair                       #Then mean number of Male Sh available for Sb pairing
    
    mprob_Hyb2_F1 <- mprob_iF_jM(i=F_Sb, j=M_Sh_Spare,parameters, var,k_i=k_H)                  #Female Sb and Male Sh
    mean_Hyb2_F1_pair<- F_Sb*mprob_Hyb2_F1                                   #Mean pairs type 6
    
    
    #LIVESTOCK 
    Mean_F_C_Sb<- MW_C_Sb* q_Sb
    Mean_M_C_Sb<- MW_C_Sb* (1-q_Sb)
    
    Mean_F_G_Sb<- MW_G_Sb* q_Sb
    Mean_M_G_Sb<- MW_G_Sb* (1-q_Sb)
    
    Mean_F_S_Sb<- MW_S_Sb* q_Sb
    Mean_M_S_Sb<- MW_S_Sb* (1-q_Sb)
    

    prob_C_Sb_Sb <- mprob_iF_jM(i=Mean_F_C_Sb, j= Mean_M_C_Sb, parameters, var,k_i=k_C)  ### Result from mating function, probability Female Sb in Cattle mated with male Sb
    prob_G_Sb_Sb <- mprob_iF_jM(i=Mean_F_G_Sb, j= Mean_M_G_Sb, parameters, var,k_i=k_G)  ### Result from mating function, probability Female Sb in Goats mated with male Sb
    prob_S_Sb_Sb <- mprob_iF_jM(i=Mean_F_S_Sb, j= Mean_M_S_Sb, parameters, var,k_i=k_S)  ### Result from mating function, probability Female Sb in sheep mated with male Sb
    
    
    Mean_couples_C <- MW_C_Sb*q_Sb*prob_C_Sb_Sb  ### Tracking number of pairs/mated females in cattle
    Mean_couples_G <- MW_G_Sb*q_Sb*prob_G_Sb_Sb  ### Tracking number of pairs/mated females in goats
    Mean_couples_S <- MW_S_Sb*q_Sb*prob_S_Sb_Sb  ### Tracking number of pairs/mated females in sheep
    
    
    #3. #### Calculating Shedding and FOI on snails/larval pool ####
    
    #Humans
    lambda_H<- dd_exp(parameters,var,Mean_F_H)*a_H*(u_d/u_s)
    
    #mira_k is the per pair type shedding rate per day for the population for each of the 7 pair types in human population
    
    mira_H_Sh <- lambda_H*N_H*mean_Sh_pair #Pair type 1
    mira_H_Hyb1<- lambda_H*N_H*(mean_Hyb_Sh_pair+mean_Hyb_Hyb_pair) #Pair type 2 & 3
    mira_H_Hyb2_Back<- lambda_H*N_H*mean_Hyb2_Back_pair #Pair type 4
    mira_H_Hyb3_Back <-lambda_H*N_H*mean_Hyb3_Back_pair #Pair type 5
    mira_H_Hyb2_F1<-lambda_H*N_H*mean_Hyb2_F1_pair #Pair type 6
    mira_H_Hyb3_F1 <- lambda_H*N_H*mean_Hyb3_F1_pair #Pair type 7
    
    mira_H_Hyb<- mira_H_Hyb1 + mira_H_Hyb2_Back + mira_H_Hyb3_Back # pooling these here for FOI on snails as going to pool at snail level
    
    #Livestock
    lambda_C<-a_C*gd_C
    lambda_G<- a_G*gd_G
    lambda_S<- a_S*gd_S
    
    mira_C_Sb<- lambda_C*N_C*Mean_couples_C *dd_exp_i(parameters,var,Mean_F_i = Mean_F_C_Sb,b_i=b_C,k_i=k_C) 
    mira_G_Sb<- lambda_G*N_G*Mean_couples_G*dd_exp_i(parameters,var,Mean_F_i = Mean_F_G_Sb,b_i=b_G,k_i=k_G)
    mira_S_Sb<- lambda_S*N_S*Mean_couples_S*dd_exp_i(parameters,var,Mean_F_i = Mean_F_S_Sb,b_i=b_S,k_i=k_S)
    
    
    mira_tot_Sb<-mira_C_Sb+mira_G_Sb+mira_S_Sb
    
    #4. Death rate worms
    
    del_tx<-(1/365)*tx_cover*tx_eff
    
    del_H<- del_nl + del_tx
    
    del_C<-del_nl+del_dcattle
    
    del_G<-del_nl+del_dg
    
    del_S<-del_nl+del_ds
    
    # 5. Dynamic equations for worms 
    
    #Humans
    
    d_MW_H_Sh <- (L_Sh*A_H_Sh)-(MW_H_Sh*del_H)
    
    #after 1 backcross, worms join general pool of hybrids and this is represented by pooling these snails for the Hyb FOI on humans 
    d_MW_H_Hyb <- (L_Hyb*A_H_Hyb)-(MW_H_Hyb*del_H) 
    
    d_MW_H_F1a <- (L_Hyb2_F1*A_H_F1a)-(MW_H_F1a*del_H)
    
    d_MW_H_F1b <- (L_Hyb3_F1*A_H_F1b)-(MW_H_F1b*del_H)
    
    d_MW_H_Sb <- (L_Sb*A_H_Sb)-(MW_H_Sb*del_H)
    
    
    #Livestock
    
    d_MW_C_Sb <- (L_Sb*A_C)-(MW_C_Sb*del_C) #rate of change in MWB Sb cattle
    
    d_MW_G_Sb<-(L_Sb*A_G)-(MW_G_Sb*del_G)
    
    d_MW_S_Sb <-(L_Sb*A_S)-(MW_S_Sb*del_S)
    
    
    # 6. dynamic equations for snails/larval pool 
    
    d_L_Sh <- (mira_H_Sh*B_H_tot)-(gamma_s_H*L_Sh) #rate of change Sh prevalence snails
    
    d_L_Hyb<- (mira_H_Hyb*B_H_tot)-(gamma_s_H*L_Hyb) #rate of change Hyb prevalence snails- this is pooling of Hyb1 and backcrosses
    
    d_L_Hyb2_F1 <- (mira_H_Hyb2_F1*B_H_tot)-(gamma_s_H*L_Hyb2_F1) #rate of change Hyb2 F1 prevalence snails
    
    d_L_Hyb3_F1 <- (mira_H_Hyb3_F1*B_H_tot)-(gamma_s_H*L_Hyb3_F1) #rate of change Hyb3 F1 prevalence snails
    
    d_L_Sb <- (mira_tot_Sb*Beta_livestock_Sb)-(gamma_s_Sb*L_Sb) #rate of change of Sb prevalence in snails set to zero pending adjustments in animal component
    
    list(c(d_MW_H_Sh, d_MW_H_Hyb, d_MW_H_F1a, d_MW_H_F1b, d_MW_H_Sb, d_MW_C_Sb, d_MW_G_Sb, d_MW_S_Sb, d_L_Sh, d_L_Hyb, d_L_Hyb2_F1, d_L_Hyb3_F1, d_L_Sb))})
  
}

####### Solving transmission dynamics model #############

parameters<-mat
MWB.init<-inits
MWB.par <- (parameters)
stepsize=7
times<-seq(0,365*15, stepsize) # YEARS

#Use this if just want to do one run to check anything, choose one set of params/inits
i=1
MWB.sol <- lsoda(MWB.init[i,], times, MWB.dyn, MWB.par[i,]) 
Output<-data.frame(MWB.sol)
plot(Output$time,Output$MW_H_Sh)
plot(Output$time,Output$MW_H_Hyb)
plot(Output$time,Output$MW_H_F1a)
plot(Output$time,Output$MW_H_F1b)
plot(Output$time,Output$MW_C_Sb)
plot(Output$time, Output$MW_G_Sb)
plot(Output$time, Output$MW_S_Sb)
plot(Output$time, Output$L_Sb)


#SIMULATIONS: Use this if want to simulate no zoonotic transmission by setting transmission of Sb to humans to zero
for(i in 1:1000){
  parameters[i,"A_H_Sb"]<-0
}


###################### Numerical solution: multiple runs #####################

data_store_all<-vector(length(1000),mode='list')

for(i in 1:1000){
  data_store_all[[i]]<-lsoda(MWB.init[i,], times, MWB.dyn, MWB.par[i,])
  
}

###Now have output from 1000 runs of model
## Extract outputs held in data_store_all
#create empty matrices and store median/confidence intervals for plotting
MW_H_Sh <- matrix(0,nrow=length(times), ncol=1000) #if want just worm tracking
MW_H_Hyb<-matrix(0,nrow=length(times), ncol=1000)
MW_H_F1a<-matrix(0,nrow=length(times),ncol=1000)
MW_H_F1b<-matrix(0,nrow = length(times),ncol=1000)
MW_H_Sb<-matrix(0,nrow=length(times),ncol=1000)
MW_TOT<-matrix(0,nrow=length(times), ncol=1000)
MW_H_Hyb_TOT<-matrix(0,nrow=length(times),ncol=1000)
MW_H_F1_TOT<-matrix(0,nrow=length(times),ncol=1000)
Prop_Hyb<-matrix(0,nrow=length(times),ncol=1000)

for(i in 1:1000){
  
  MW_H_Sh[,i]<-data_store_all[[i]][,"MW_H_Sh"]
  MW_H_Hyb[,i]<-data_store_all[[i]][,"MW_H_Hyb"]
  MW_H_F1a[,i]<-data_store_all[[i]][,"MW_H_F1a"]
  MW_H_F1b[,i]<-data_store_all[[i]][,"MW_H_F1b"]
  MW_H_Sb[,i]<-data_store_all[[i]][,"MW_H_Sb"]
  MW_TOT[,i]<-MW_H_Sh[,i]+MW_H_Hyb[,i]+MW_H_F1a[,i]+MW_H_F1b[,i]+MW_H_Sb[,i]
  MW_H_Hyb_TOT[,i]<-MW_H_Hyb[,i]+MW_H_F1a[,i]+MW_H_F1b[,i]
  MW_H_F1_TOT[,i]<-MW_H_F1a[,i]+MW_H_F1b[,i]
  Prop_Hyb[,i]<-MW_H_Hyb_TOT[,i]/MW_TOT[,i]
  
}


MW_H_Sh_Median<-apply(MW_H_Sh,1,median)
MW_H_Sh_lower<-apply(MW_H_Sh,1,quantile,probs=0.025)
MW_H_Sh_upper<-apply(MW_H_Sh,1,quantile,probs=0.975)
MW_H_Hyb_Median<-apply(MW_H_Hyb,1,median)
MW_H_Hyb_lower<-apply(MW_H_Hyb,1,quantile,probs=0.025)
MW_H_Hyb_upper<-apply(MW_H_Hyb,1,quantile,probs=0.975)

MW_H_F1_TOT_Median<-apply(MW_H_F1_TOT,1,median)
MW_H_F1_TOT_lower<-apply(MW_H_F1_TOT,1,quantile,probs=0.025)
MW_H_F1_TOT_upper<-apply(MW_H_F1_TOT,1,quantile,probs=0.975)
MW_H_Sb_Median<-apply(MW_H_Sb,1,median)
MW_H_Sb_lower<-apply(MW_H_Sb,1,quantile,probs=0.025)
MW_H_Sb_upper<-apply(MW_H_Sb,1,quantile,probs=0.975)


MW_TOT_median<-apply(MW_TOT, 1,median)
MW_TOT_upper<-apply(MW_TOT,1,quantile,probs=0.975)
MW_TOT_lower<-apply(MW_TOT,1,quantile,probs=0.025)
Prop_Hyb_median<-apply(Prop_Hyb,1,median)
Prop_Hyb_upper<-apply(Prop_Hyb,1,quantile,probs=0.975)
Prop_Hyb_lower<-apply(Prop_Hyb,1,quantile,probs=0.025)


#Then put into dataframe for ggplotting

dfMWB_All<-data.frame(times, MW_TOT_median, MW_TOT_lower, MW_TOT_upper,
                           Prop_Hyb_median, Prop_Hyb_lower, Prop_Hyb_upper,
                           MW_H_Sh_Median, MW_H_Sh_lower, MW_H_Sh_upper,
                           MW_H_F1_TOT_lower,MW_H_F1_TOT_Median,MW_H_F1_TOT_upper,
                           MW_H_Sb_Median,MW_H_Sb_lower,MW_H_Sb_upper,
                           MW_H_Hyb_Median,MW_H_Hyb_lower,MW_H_Hyb_upper)

library(ggplot2)

#PLOT SH +Hyb together 
Human_Sh_hyb_worms_v_time<- ggplot(data=dfMWB_All, aes(x=times/365))+
  geom_line(aes(y=MW_H_Sh_Median),colour="#FFCC66",size=0.8)+
  geom_line(aes(y=MW_H_Sh_upper), colour="#FFCC66", linetype='dashed')+
  geom_line(aes(y=MW_H_Sh_lower),colour="#FFCC66",linetype='dashed')+
  geom_line(aes(y=MW_H_Hyb_Median),colour="dark green",size=0.8)+
  geom_line(aes(y=MW_H_Hyb_upper), colour="dark green", linetype='dashed')+
  geom_line(aes(y=MW_H_Hyb_lower),colour="dark green",linetype='dashed')+
  xlab('Time (years)')+
  ylab('Mean Worm Burden')+
  ylim(0,45)+
  labs(title = "Multirun_Sh_Hyb")
Human_Sh_hyb_worms_v_time

Human_SB_F1_worms_v_time<- ggplot(data=dfMWB_All, aes(x=times/365))+
  geom_line(aes(y=MW_H_Sb_Median),colour="#336699",size=0.8)+
  geom_line(aes(y=MW_H_Sb_upper), colour="#336699", linetype='dashed')+
  geom_line(aes(y=MW_H_Sb_lower),colour="#336699",linetype='dashed')+
  geom_line(aes(y=MW_H_F1_TOT_Median),colour="#FFCC66",size=0.8)+
  geom_line(aes(y=MW_H_F1_TOT_Median),colour="#336699",size=0.8,linetype='dashed')+
  geom_line(aes(y=MW_H_F1_TOT_upper), colour="#FFCC66", linetype='dashed')+
  geom_line(aes(y=MW_H_F1_TOT_lower),colour="#FFCC66",linetype='dashed')+
  xlab('Time (years)')+
  ylab('Mean Worm Burden')+
  ylim(0,1)+
  labs(title = "Multirun_SB_F1")
Human_SB_F1_worms_v_time


Human_SB_v_time<- ggplot(data=dfMWB_All, aes(x=times/365))+
  geom_line(aes(y=MW_H_Sb_Median),colour="#336699",size=0.8)+
  geom_line(aes(y=MW_H_Sb_upper), colour="#336699", linetype='dashed')+
  geom_line(aes(y=MW_H_Sb_lower),colour="#336699",linetype='dashed')+
  xlab('Time (years)')+
  ylab('Mean Worm Burden')+
  ylim(0,0.25)+
  labs(title = "Multirun_SB")
Human_SB_v_time


All_worms_v_time<- ggplot(data=dfMWB_All, aes(x=times/365))+
  geom_line(aes(y=MW_TOT_median),colour="Black",size=0.8)+
  geom_line(aes(y=MW_TOT_upper), colour="Black", linetype='dashed')+
  geom_line(aes(y=MW_TOT_lower),colour="Black",linetype='dashed')+
  xlab('Time (years)')+
  ylab('Mean Worm Burden:All genotypes')+
  ylim(0,50)+
  labs(title = "Multirun_all_worms")
All_worms_v_time


Prop_hyb_v_time<- ggplot(data=dfMWB_All, aes(x=times/365))+
  geom_line(aes(y=Prop_Hyb_median),colour="Black",size=0.8)+
  geom_line(aes(y=Prop_Hyb_upper), colour="Black", linetype='dashed')+
  geom_line(aes(y=Prop_Hyb_lower),colour="Black",linetype='dashed')+
  xlab('Time (years)')+
  ylab('Proportion of worm burden hybrids')+
  ylim(0,1)+
  labs(title = "Multirun_prop_hybrids")
Prop_hyb_v_time
