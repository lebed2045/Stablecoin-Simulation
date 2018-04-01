set.seed(9)

maturity <- 5 #number of years the simulation runs for
dt <-  1/365 # the lengths in between simulations (ie period length)
simulation.length <- (maturity/dt) #number of days in the simulation
pB <- rep(1, times=simulation.length) #initialize the price vector with 1s
q <- b <- f <- pricing_mu <- pricing_sigma <- sigma_est <- mu_est <- M <- rep(0, times=simulation.length) #fill the q,b,f, and M vectors with zeroes
global_mu <- pricing_mu[1] <- -0.05 #annual drift in M
global_sigma <- pricing_sigma[1] <- 0.2 #annual volatility
M[1] <- M0 <- 1000 #intial market cap
b[1] <- b0 <- 300 #initial b of basebonds is 0 
q[1] <- q0 <- 1 #initial q of basecoin is 10 billion
bayesian <- 1  #set this to 1 to run the Bayesian model; set to 0 to run the perfect information model
is_weighted <- 0 #weight assigned to inial belief
mu_belief <- .05 #initial belief about drift in the Bayesian model
belief_weight <- 50 #weighting assigned to the intial belief about drift
Expiration <- 5/dt #expiration of bonds in days
ages <- list("0"=b0) #initialize ages
short_budget <- 500
auction_bundles <- 1 #sets number of bonds sold in each auction
CET <- 3654
#illiq_power <- log((short_budget+.08*(M0-short_budget))/(.08*(M0-short_budget)))/CET
illiq_power_v2 <- 1-(log( (.04*(M0-short_budget)) + (0*short_budget) )/5.43)
illiq_power_v2

#function that returns the key (number of days) of the oldest tranche of bonds
max_age_key <- function(ages){
  return(max(as.numeric(names(ages))))
}

#function that returns the value (number of bonds) of the oldest tranche of bonds
max_age_value <- function(ages){
  return(as.numeric(ages[toString(max(as.numeric(names(ages))))]))
}

#function that returns the total number of bonds in ages
sum_ages <- function(ages){
  return(sum(as.numeric(ages)))
}

#function that returns ages with the oldest tranche of bonds removed
remove_oldest <- function(ages){
  if(length(ages)>0) ages[toString(max(as.numeric(names(ages))))]<-NULL
  return(ages)
}

#this is the pricing function for bond auctions
price_B <- function(Mi, B, dt=1/365, Expiration=5/dt, pb_ages,mu,sigma){
  epv <- ddr <- rep(1,Expiration)
  Lwr <- Not_P <- rep(1,Expiration)
  P <- RP <- rep(0,Expiration)
  
  ages_mat <- cbind(as.numeric(names(pb_ages)), as.numeric(pb_ages))
  
  for(j in 1:Expiration){
    unexpired <- as.numeric(ages_mat[,1]<(1+Expiration-j)) #tests whether bonds in each row are unexpired
    outstanding_b <- sum(ages_mat[,2]*unexpired) #number of bonds outstanding in each period
    Mr <- outstanding_b+1+Mi #required increase in M for bond purchased to be paid back
    Lwr[j] <- (log((Mr/Mi)*exp(-(mu-sigma^2/2)*(j*dt)))*(1/sigma)*(j*dt)^(-1/2)) #lower bound of the normal PDF integral
    P[j] <-1-pnorm(Lwr[j])  #probability Z, a standard normal variable, is greater than than the lower bound
    RP[j] <- prod(Not_P)*P[j] #probability that the bond will pay back (payout probability given that it has not already been paid)
    Not_P[j] <- 1-P[j] #probability that the bond has not been redeemed in j
    ddr[j] <- ((1+0.0168)^(j*dt)) #daily discount rate
    epv[j] <- RP[j]/ddr[j] #expected present value of the bond for each period in j
  }
  #print(paste("epv",sum(epv),"Mi",Mi,"b",B,"Mr",Mr))
  return(sum(epv))
}

pB[1] <- price_B(M0,b0,dt,pb_ages=ages,mu=pricing_mu[i],sigma=pricing_sigma[i])

#for every day in the simulation
for(i in 2:simulation.length){
  
  if(i>1) { #if it's not the first day
    b[i] <- b[i-1] #set b[i] to be b[i-1]
  } else b[i] <- b0 #if it is the first day, set b[i] to be b0
  
  names(ages) <- as.numeric(names(ages))+1 #increase the age of all bonds by 1 day
  
  while(b[i]>0 && max_age_key(ages)>Expiration){ #while there are bonds older than expiration
    b[i] <- b[i] - max_age_value(ages) #reduce the number of bonds by the oldest bucket
    ages <- remove_oldest(ages) #remove the oldest bucket from ages
  }
  
  print(paste("i: ",i," ,b:", b[i]))
  
  f[i] <- f[i-1]+sqrt(dt)*rnorm(1) #f simulates brownian motion
  M[i] <- M[1]*exp((global_mu-(global_sigma^2)/2)*(i*dt)+global_sigma*f[i]) #market cap calc
  delta <- (M[i]-M[i-1]) #delta in market cap from last period
  
  sigma_est[i] <- sqrt(mean((diff(M[1:i]))^2/(M[2:i]^2*dt)))
  if(is_weighted==1){
    mu_est[i] <- ((.05*belief_weight)+(weighted.mean((log(M[2:i])-log(M[1:i-1])),seq(1,.1,-.9/(i-2)))/(dt))*i)/(i+belief_weight)
  }else {
    mu_est[i] <- ((.05*belief_weight)+(mean(log(M[2:i])-log(M[1:i-1]))/dt)*i)/(i+belief_weight)
  }
  
  if(bayesian==1){
    pricing_sigma[i] <- sigma_est[i]
    pricing_mu[i] <- mu_est[i]
  } else{
    pricing_sigma[i] <- global_sigma
    pricing_mu[i] <- global_mu
  }
  
  if(delta>0){ #if market cap (price) increased
    if(b[i]<delta){ #if the increase exceeds the number of bonds
      b[i]<-0 #should ideally issue coins to shareholders, but zeroes out the bonds
      ages <- list("0"=0) #resets ages
    } else { #otherwise if the increase does not exceed the number of bonds
      b[i] <- b[i] - delta #then decrease the number of bonds by the increase in price
      to_destroy <- delta #set the amount of bonds to destroy from ages equal to the increase in price
      while(sum_ages(ages) >=1 && to_destroy>0){ #while there are bonds to destroy from ages
        if(max_age_value(ages)>=to_destroy){ #if number of bonds in the oldest age bucket is >= the amount to destroy
          ages[toString(max_age_key(ages))] <- max_age_value(ages) - to_destroy #then reduce the oldest bucket by the amount to destroy
          to_destroy<-0 #and reset to_destroy to zero
        } else { #otherwise if number of bonds in the oldest bucket is < the amount to destroy
          to_destroy <- to_destroy - max_age_value(ages) #then reduce the amount to destroy by the number of bonds in the oldest bucket
          ages <-remove_oldest(ages) #and remove the oldest bucket
        } #end of number of bonds in the oldest bucket is < the amount to destroy case
      }#will re-enter the while loop if there are still bonds to destroy
    }
  } else{ #otherwise if market cap (price) decreased or stayed the same
    value_deficit <- -1*delta #set the amount of deficit to offset to be -1 times the delta
    while(value_deficit > 0){ #while the deficit is positive
      if(value_deficit < auction_bundles*price_B(M[i],b[i],dt,pb_ages=ages,mu=pricing_mu[i],sigma=pricing_sigma[i])) { #if the deficit is less than the price of auctioning a bond
        value_deficit<-0 #then reset the deficit and do nothing else (dont sell a bond)
      } else{ #otherwise if the deficit is >= the price of auctioning a bond
        value_deficit <- value_deficit - (auction_bundles*price_B(M[i],b[i],dt,pb_ages=ages,mu=pricing_mu[i],sigma=pricing_sigma[i])) #then reduce deficit by the price of a bond
        b[i] <- (b[i] + 1*auction_bundles) #and increase b by the amount being auctioned 
        if(ages["0"]=="NULL") ages["0"] <- 0 #if there are no new bundles in ages then set the amount equal to 0
        ages["0"] <- (as.numeric(ages["0"]) + auction_bundles) #increase the number of new bonds in ages by auction_bundles
      } #end of deficit is >= the price of auctioning a bond case
      short_power <- (short_budget*illiq_power_v2)
      if(short_budget>0 && price_B(M[i]-short_power,b[i],dt,pb_ages=ages,mu=pricing_mu[i],sigma=pricing_sigma[i])<=.1){
        print("he sells")
        M[i]<-M[i]-short_power
        value_deficit<-value_deficit+short_power
        short_budget <- 0
      }
    } #will re-enter the loop if there is still deficit to offset
  } #end of market cap (price) decreased case
  
  pB[i] <- price_B(M[i],b[i],dt,pb_ages=ages,mu=pricing_mu[i],sigma=pricing_sigma[i]) #save the price in the price vector
}

par(mfrow=c(1,3))
plot(M)
plot(b)
plot(pB)
summary(M)
summary(pB)