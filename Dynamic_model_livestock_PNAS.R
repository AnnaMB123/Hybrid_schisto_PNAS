library(deSolve)
library(MASS)

### Code to run livestock Sb part of model ###

#Import matrices of parameter and initial values#
mat<-read.table("~/Desktop/Hybrid_schisto_pnas_code_datafiles/livestock_parameter_matrix")
inits<-read.table("~/Desktop/Hybrid_schisto_pnas_code_datafiles/livestock_init_matrix")

tx_cover<-rep(0,1000) #Specifiy treatment cover
tx_eff<-rep(0.95,1000) #Specify treatment efficacy

mat<-as.matrix(mat)
inits<-as.matrix(inits)
mat<-cbind(mat,tx_cover,tx_eff)

parameters<- mat
MWB.init<-inits

#### MATING FUNCTIONS #####

#Mating function as per May Woolhouse 1993, allowing for unequal numbers of males and females. Probability that a female of species i will form a pair with a male of species j

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


#Mean density-dependent reduction in fecundity per mated pair- exponential relationship
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
    
    #1.###SPECIFY VARIABLES TO BE TRACKED HERE Currently animal population static, worm& snail death is included but assumed snails and animals are immediately replaced
    
    #CATTLE
    MW_C_Sb <- var[1] ###Mean worm burden per cow
    
    MW_G_Sb <- var[2]  ###Mean worm burden per sheep
    
    MW_S_Sb <- var[3]  ###Mean worm burden per goat
    
    #LARVAL POOL
    
    L_Sb <- var[4]   #Baseline L_Sb fraction
    
    
    #2. ### Calculate mean females and males in each species ######
    
    ### CATTLE ###
    
    Mean_F_C_Sb<- MW_C_Sb* q_Sb
    Mean_M_C_Sb<- MW_C_Sb* (1-q_Sb)
    
    ### GOATS ###
    
    Mean_F_G_Sb<- MW_G_Sb* q_Sb
    Mean_M_G_Sb<- MW_G_Sb* (1-q_Sb)
    
    ### SHEEP ###
    
    Mean_F_S_Sb<- MW_S_Sb* q_Sb
    Mean_M_S_Sb<- MW_S_Sb* (1-q_Sb)
    
    
    ## 3. Call mating function for each species to calculate mating probability and then mean number of pairs
    
    ###############################SB MATING IN CATTLE #####################
    
    prob_C_Sb_Sb <- mprob_iF_jM(i=Mean_F_C_Sb, j= Mean_M_C_Sb, parameters, var,k_i=k_C)  ### Result from mating function, probability Female Sb in Cattle mated with male Sb
    
    Mean_couples_C <- MW_C_Sb*q_Sb*prob_C_Sb_Sb  ### Tracking number of mated females
    
    ###############################SB MATING IN GOATS #########################
    
    prob_G_Sb_Sb <- mprob_iF_jM(i=Mean_F_G_Sb, j= Mean_M_G_Sb, parameters, var,k_i=k_G)  ### Result from mating function, probability Female Sb in Goats mated with male Sb
    
    Mean_couples_G <- MW_G_Sb*q_Sb*prob_G_Sb_Sb  ### Tracking number of mated females
    
    
    ##############################SB MATING IN SHEEP ########################
    
    
    prob_S_Sb_Sb <- mprob_iF_jM(i=Mean_F_S_Sb, j= Mean_M_S_Sb, parameters, var,k_i=k_S)  ### Result from mating function, probability Female Sb in sheep mated with male Sb
    
    Mean_couples_S <- MW_S_Sb*q_Sb*prob_S_Sb_Sb  ### Tracking number of mated females
    
    
    #3. ### FOI ON SNAILS/LARVAL POOL/SHEDDING PARAMETERS  ####
    
    # lambda_i is per pair per day miracidia shedding rate for species i, in absence of dd effects
    
    lambda_C<-a_C*gd_C
    
    lambda_G<- a_G*gd_G
    
    lambda_S<- a_S*gd_S
    
    # mira here is the the per animal population daily shedding rate including mean dd effect
    
    mira_C<- lambda_C*N_C*Mean_couples_C *dd_exp_i(parameters,var,Mean_F_i = Mean_F_C_Sb,b_i=b_C,k_i=k_C)
    
    mira_G<- lambda_G*N_G*Mean_couples_G*dd_exp_i(parameters,var,Mean_F_i = Mean_F_G_Sb,b_i=b_G,k_i=k_G)
    
    mira_S<- lambda_S*N_S*Mean_couples_S*dd_exp_i(parameters,var,Mean_F_i = Mean_F_S_Sb,b_i=b_S,k_i=k_S)
    
    mira_tot<- mira_C+mira_G+mira_S
    
    #4.Death rate of worms overall including tx
    
    del_tx<-(1/365)*tx_cover*tx_eff
    
    del_C<-del_nl+del_dcattle+del_tx
    del_G<-del_nl+del_dg
    del_S<-del_nl+del_ds
    
    
    #5.## DERIVATIVES #####
    
    d_MW_C_Sb <- (L_Sb*A_C)-(MW_C_Sb*del_C) #rate of change in MWB Sb cattle
    
    d_MW_G_Sb <- (L_Sb*A_G)-(MW_G_Sb*del_G) #rate of change in MWB Sb goats
    
    d_MW_S_Sb <- (L_Sb*A_S)-(MW_S_Sb*del_S) #rate of change in MWB Sb sheep
    
    d_L_Sb <- (mira_tot*Beta_tot)-(gamma*L_Sb) #rate of change Sb prevalence snails
    
    
    
    list(c(d_MW_C_Sb, d_MW_G_Sb, d_MW_S_Sb,d_L_Sb))})
  
}


####### Solving transmission dynamics model #############

MWB.par <- parameters

stepsize=7
times<-seq(0,365*10, stepsize) #10 YEARS

###################### Numerical solution #####################


data_store_all_sim<-vector(length(1000),mode='list')


for(i in 1:1000){
  data_store_all_sim[[i]]<-lsoda(MWB.init[i,], times, MWB.dyn, MWB.par[i,])
  
}


#Use this if just want to do one run to check anything
i=1 #Choose any from parameter/init sets of 1000
MWB.sol <- lsoda(MWB.init[i,], times, MWB.dyn, MWB.par[i,]) 
Output<-data.frame(MWB.sol)

MW_C_Sb <- matrix(0,nrow=length(times), ncol=1000) 
MW_G_Sb<-matrix(0,nrow=length(times), ncol=1000)
MW_S_Sb<-matrix(0,nrow=length(times),ncol=1000)
L_Sb<-matrix(0,nrow=length(times),ncol=1000)
L_Sb_fraction<-(matrix(0,nrow=length(times),ncol=1000))

for(i in 1:1000){
  
  MW_C_Sb[,i]<-data_store_all_sim[[i]][,"MW_C_Sb"]
  MW_G_Sb[,i]<-data_store_all_sim[[i]][,"MW_G_Sb"]
  MW_S_Sb[,i]<-data_store_all_sim[[i]][,"MW_S_Sb"]
  L_Sb[,i]<-data_store_all_sim[[i]][,"L_Sb"]
  
}

MW_C_median<-apply(MW_C_Sb, 1,median)
MW_C_upper<-apply(MW_C_Sb,1,quantile,probs=0.975)
MW_C_lower<-apply(MW_C_Sb,1,quantile,probs=0.025)
MW_G_median<-apply(MW_G_Sb, 1,median)
MW_G_upper<-apply(MW_G_Sb,1,quantile,probs=0.975)
MW_G_lower<-apply(MW_G_Sb,1,quantile,probs=0.025)
MW_S_median<-apply(MW_S_Sb, 1,median)
MW_S_upper<-apply(MW_S_Sb,1,quantile,probs=0.975)
MW_S_lower<-apply(MW_S_Sb,1,quantile,probs=0.025)
L_Sb_median<-apply(L_Sb,1,median)
L_Sb_upper<-apply(L_Sb,1,quantile,probs=0.975)
L_Sb_lower<-apply(L_Sb,1,quantile,probs=0.025)
L_Fraction_median<-L_Sb_median/0.16*100
L_Fraction_upper<-L_Sb_upper/0.16*100
L_Fraction_lower<-L_Sb_lower/0.16*100

#Then put into dataframe for ggplotting

dfMWB_All_Sb<-data.frame(times, MW_C_median, MW_C_lower, 
                         MW_C_upper,MW_G_median, MW_G_lower,
                         MW_G_upper,MW_S_median, MW_S_lower, 
                         MW_S_upper,L_Fraction_median, L_Fraction_upper, L_Fraction_lower)

##PLOTTING##

library(ggplot2)
MW_C<- ggplot(data=dfMWB_All_Sb, aes(x=times/365))+
  geom_line(aes(y=MW_C_median),colour="#336699")+
  geom_line(aes(y=MW_C_lower),colour="#336699",linetype="dashed")+
  geom_line(aes(y=MW_C_upper),colour="#336699",linetype="dashed")+
   xlab('Time (years)')+
  ylab('Mean Worm burden')+
  labs(title = "Cattle")+
  ylim(0,100)
MW_C

MW_G<- ggplot(data=dfMWB_All_Sb, aes(x=times/365))+
  geom_line(aes(y=MW_G_median),colour="#336699")+
  geom_line(aes(y=MW_G_lower),colour="#336699",linetype="dashed")+
  geom_line(aes(y=MW_G_upper),colour="#336699",linetype="dashed")+
  xlab('Time (years)')+
  ylab('Mean Worm burden')+
  labs(title = "Goats")+
  ylim(0,3)
MW_G

MW_S<- ggplot(data=dfMWB_All_Sb, aes(x=times/365))+
  geom_line(aes(y=MW_S_median),colour="#336699")+
  geom_line(aes(y=MW_S_lower),colour="#336699",linetype="dashed")+
  geom_line(aes(y=MW_S_upper),colour="#336699",linetype="dashed")+
  xlab('Time (years)')+
  ylab('Mean Worm burden')+
  labs(title = "Sheep")+
  ylim(0,3)
MW_S

L_fraction<- ggplot(data=dfMWB_All_Sb, aes(x=times/365))+
  geom_line(aes(y=L_Fraction_median),colour="#336699")+
  geom_line(aes(y=L_Fraction_lower),colour="#336699",linetype="dashed")+
  geom_line(aes(y=L_Fraction_upper),colour="#336699",linetype="dashed")+
  xlab('Time (years)')+
  ylab('Percentage of Larval Pool from Baseline')+
  labs(title = "Larval fraction")+
 scale_y_continuous(limits=c(0,125),breaks=c(0,25,50,75,100,125))
  L_fraction


