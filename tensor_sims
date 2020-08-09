## measles simulations and tensor decomp

iota = .1 
imports <- rbinom(Tmax,1,iota)
tsir <- function(I0=10,S0=1000,N=300000, beta=15,alpha=0.97,gamma=0.2, births=rep((30*300/26),26),
                 Tmax=2600, iota=.05, nsim = 10, growth = FALSE, shape = FALSE, inv_shape = FALSE, step = FALSE){
  Istore <- Sstore <- rep(NA,Tmax) 
  Istore[1] <- I0
  Sstore[1] <- S0
  seas.index = rep(1:26,Tmax/26)
  beta.seas <- beta*(1+gamma*cos(2*pi*(1:26)/26))
  imports <- imports
  
  dat <- data.frame(pop=rep(0,nsim*Tmax),
                    Istore=rep(0,nsim*Tmax), Sstore = rep(0,nsim*Tmax),
                    sim = rep(0, nsim*Tmax),time=rep(0,nsim*Tmax))
  
  if(growth == FALSE){
    for(n in 1:nsim){
      for (t in 2:Tmax){
        Istore[t] <- rpois(1, beta.seas[seas.index[t]]*Sstore[t-1]*((Istore[t-1]+imports[t])^alpha)/N)
        Istore[t] <- ifelse(Istore[t] > Sstore[t-1], Sstore[t-1], Istore[t])
        Sstore[t] <- Sstore[t-1]+(births[seas.index[t]])-Istore[t]}
      dat[((n-1)*Tmax+1):(n*Tmax), ] <- data.frame(pop=rep(N,Tmax),Istore=Istore,Sstore=Sstore,sim=rep(n,Tmax),time=1:Tmax)
    }
  }
  if(growth == TRUE & shape == TRUE){
    if(inv_shape == TRUE){births_growth = mean(births)*rep((1 - .75*sin(2*pi*1:260/260)),10)
    }else{births_growth = mean(births)*rep((1 + .75*sin(2*pi*1:260/260)),10)}
    
    for(n in 1:nsim){
      for (t in 2:Tmax){
        Istore[t] <- rpois(1, beta.seas[seas.index[t]]*Sstore[t-1]*((Istore[t-1]+imports[t])^alpha)/N)
        Istore[t] <- ifelse(Istore[t] > Sstore[t-1], Sstore[t-1], Istore[t])
        Sstore[t] <- Sstore[t-1]+(births_growth[t])-Istore[t]}
      dat[((n-1)*Tmax+1):(n*Tmax), ] <- data.frame(pop=rep(N,Tmax),Istore=Istore,Sstore=Sstore,sim=rep(n,Tmax),time=1:Tmax)
    } 
  }
  if(growth == TRUE & shape == FALSE){
    births_growth = mean(15*300/26) + c(rep(0, 1100), .15*(1:(Tmax- 500)))
    for(n in 1:nsim){
      for (t in 2:Tmax){
        Istore[t] <- rpois(1, beta.seas[seas.index[t]]*Sstore[t-1]*((Istore[t-1]+imports[t])^alpha)/N)
        Istore[t] <- ifelse(Istore[t] > Sstore[t-1], Sstore[t-1], Istore[t])
        Sstore[t] <- Sstore[t-1]+(births_growth[t])-Istore[t]}
      dat[((n-1)*Tmax+1):(n*Tmax), ] <- data.frame(pop=rep(N,Tmax),Istore=Istore,Sstore=Sstore,sim=rep(n,Tmax),time=1:Tmax)
    } 
  }
  if(growth == TRUE & step == TRUE){
    births_growth = mean(births)*(1 - c(rep(.6,2300), rep(-.2, 300)))
    for(n in 1:nsim){
      for (t in 2:Tmax){
        Istore[t] <- rpois(1, beta.seas[seas.index[t]]*Sstore[t-1]*((Istore[t-1]+imports[t])^alpha)/N)
        Istore[t] <- ifelse(Istore[t] > Sstore[t-1], Sstore[t-1], Istore[t])
        Sstore[t] <- Sstore[t-1]+(births_growth[t])-Istore[t]}
      dat[((n-1)*Tmax+1):(n*Tmax), ] <- data.frame(pop=rep(N,Tmax),Istore=Istore,Sstore=Sstore,sim=rep(n,Tmax),time=1:Tmax)
    } 
  }
  return(dat[dat$time > 1600,]) }


dat <- tsir(beta = 15, growth = TRUE, shape = TRUE, alpha = .97, iota = .1)
dat %>%
  filter(sim == 10) %>%
  ggplot(aes(y=Istore,x=time)) + geom_line()

dat <- tsir(beta = 15, growth = TRUE, shape = TRUE, inv_shape = TRUE, alpha = .97, iota = .1)
dat %>%
  filter(sim == 10) %>%
  ggplot(aes(y=Istore,x=time)) + geom_line()


dat <- tsir(beta = 18, growth = FALSE, alpha = .97, iota = .1)
dat %>%
  filter(sim == 10) %>%
  ggplot(aes(y=Istore,x=time)) + geom_line()


dat <- tsir(beta = 15, growth = TRUE, step = TRUE, alpha = .97, iota = .1)
dat %>%
  filter(sim == 4) %>%
  ggplot(aes(y=Istore,x=time)) + geom_line()

tmp <- tsir(beta = 15, iota = .1, growth = FALSE, births = rep((15*300/26),26))
tmp %>%
  filter(sim == 4) %>%
  ggplot(aes(y=Istore,x=time)) + geom_line()

case_data <- matrix(nrow = 520, ncol = 150)
for(i in 1:ncol(case_data)){
  if(i<=50){
    tmp <- tsir(beta = 15, iota = .1, growth = FALSE, nsim = 1, births = rep((15*300/26),26))
    case_data[,i] <- tmp$Istore[481:1000]
  }
  if(i>= 51 & i <= 100){
    tmp <- tsir(beta = 15, iota = .1, growth = TRUE, nsim = 1)
    case_data[,i] <- tmp$Istore[481:1000]
  }
  if(i>= 101 & i <= 150){
    tmp <- tsir(beta = 15, iota = .1, growth = TRUE, step = TRUE, nsim = 1)
    case_data[,i] <- tmp$Istore[481:1000]
  }
  if(i>= 151){
    tmp <- tsir(beta = 15, iota = .1, growth = TRUE, step = TRUE, nsim = 1)
    case_data[,i] <- tmp$Istore[480:1000]
  }
}

plot(case_data[,25], type = 'l')
lines(case_data[,5], col = 'blue')
lines(case_data[,105], col = 'red')


plot(case_data[,117], type = 'l')
lines(case_data[,10], col = 'blue')
lines(case_data[,75], col = 'red')
cases <- c(case_data[,10],case_data[,62],case_data[,112])
group <- c(rep(1, nrow(case_data)), rep(2, nrow(case_data)), rep(3,nrow(case_data)))
sample_data <- data.frame(cases = cases, group = group, time = time_vals )

cases <- c(ewMu4465[,442], ewMu4465[,'Manchester'], ewMu4465[,'Norwich'])
district <- c(rep('London', nrow(ewMu4465)), rep('Manchester', nrow(ewMu4465)), rep('Norwich', nrow(ewMu4465)))
sample_data_real <- data.frame(cases = cases, district = district, time = as.numeric(rownames(ewMu4465)))

cols <- c('black','dodgerblue','darkred')
sims<-ggplot(data = sample_data) + geom_line(aes(x=time, y = cases, group = group, color = as.factor(group)), size = 1, alpha = .85) + theme_minimal() + 
  scale_color_manual(values = cols, name = 'group') + ggtitle('Sample Simulations')
real <- ggplot(data = sample_data_real) + geom_line(aes(x=time, y = cases, group = district, color = as.factor(district)), size = 1, alpha = .85) + theme_minimal() + 
  scale_color_manual(values = cols, name = 'city') + ggtitle('Sample Data')

plot_grid(sims, real, labels=c('A','B'), nrow=2)

case_log <- log(case_data +1)
time_vals <- 1:nrow(case_log)/26
data <- case_log[,1]
wdat <- data.frame(time = time_vals, values = data)
test <- wt(wdat, pad = TRUE)
plot(test)

sims_array <- array(NA, dim = c(ncol(case_data), nrow(test$power), ncol(test$power)))
sims_array[1,,] <- test$power
for(i in 2:ncol(case_data)){
  data <- case_log[,i]
  wdat <-  data.frame(time = time_vals, values = data)
  twt <- wt(wdat, pad = TRUE)
  sims_array[i,,] <- twt$power
}

sim_ten <- as.tensor(sims_array)

sim_decom <- cp(sim_ten, num_components = 3, max_iter = 150)
sim_decom$conv

sim_decom$norm_percent


## can also attempt 4
sim_decom2 <- cp(sim_ten, num_components = 4, max_iter = 150)
sim_decom2$conv
sim_decom2$norm_percent


score_data <- data.frame(score1 = sim_decom$U[[1]][,1],
                         score2 = sim_decom$U[[1]][,2],
                         score3 = sim_decom$U[[1]][,3],
                         group = c(rep(1,50), rep(2,50), rep(3,50)))
obs <- kmeans(score_data[,1:3], centers = 3)
score_data$clus <- obs$cluster

box1 <- ggplot(data = score_data) + geom_boxplot(aes(x = group, y = score1, group = group),fill = 'grey70') + ggtitle('First Component') + ylab('score') +
  theme_minimal()
box2 <- ggplot(data = score_data) + geom_boxplot(aes(x = group, y = score2, group = group), fill = 'grey70') + ggtitle('Second Component') +
  ylab('score') + theme_minimal()
box3 <- ggplot(data = score_data) + geom_boxplot(aes(x = group, y = score4, group = group), fill = 'grey70')+ ggtitle('Third Component') +
  ylab('score') + theme_minimal()

sim_period_data <- data.frame(period1 = period1, period2= period2, period3 = period3, period = twt$period)
sim_time_data <- data.frame(time1=time1, time2=time2,time3=time3,time=time_vals)

p1 <- ggplot(sim_period_data) + geom_line(aes(x=period,y=period1)) + ylab('power') + theme_minimal() +
  geom_vline(xintercept = 1, lty = 3, col = 'red') + geom_vline(xintercept = 2, lty = 3, col = 'red')
p2 <- ggplot(sim_period_data) + geom_line(aes(x=period,y=period2)) + ylab('power') + theme_minimal() +
  geom_vline(xintercept = 1.5, lty = 3, col = 'red') + geom_vline(xintercept = 3, lty = 3, col = 'red')
p3 <- ggplot(sim_period_data) + geom_line(aes(x=period,y=period3)) + ylab('power') + theme_minimal() +
  geom_vline(xintercept = 1, lty = 3, col = 'red') + geom_vline(xintercept = 2, lty = 3, col = 'red') + geom_vline(xintercept = 3, lty = 3, col = 'red')

t1<- ggplot(data=sim_time_data) + geom_line(aes(x=time,y=time1)) + ylab('power') + theme_minimal()
t2<-ggplot(data=sim_time_data) + geom_line(aes(x=time,y=time2)) + ylab('power') + theme_minimal()
t3<-ggplot(data=sim_time_data) + geom_line(aes(x=time,y=time3)) + ylab('power') + theme_minimal() +
  geom_hline(yintercept = 0, col='dodgerblue', lty=3)

scores <- plot_grid(box1,box2,box3,nrow=1)
periods <- plot_grid(p1,p2,p3,nrow=1)
times <- plot_grid(t1,t2,t3,nrow=1)

plot_grid(scores, periods, times, labels = c('A','B','C'), ncol = 1)
