## tensor decomp for measles data

sub_meas <- ewMu4465
sub_meas <- log((sub_meas) + 1)
sub_meas <- apply(sub_meas, 2, function(c)(c - mean(c)))
sub_meas <- apply(sub_meas, 2, function(c)scale(c, center = 0, scale = 1))
wdat <- sub_meas[,1]
time <-  as.numeric(rownames(ewMu4465))
wdat <- data.frame(time = time, values = wdat)
test <- wt(wdat, pad = TRUE)

sub_array <- array(NA, dim = c(ncol(sub_meas), nrow(test$power.corr), ncol(test$power.corr)))
sub_array[1,,] <- test$power.corr

for(i in 2:ncol(sub_meas)){
  wdat <- sub_meas[,i]
  time <-  as.numeric(rownames(ewMu4465))
  wdat <- data.frame(time = time, values = wdat)
  wtran <- wt(wdat, pad = TRUE)
  #plot(wtran, main = paste(colnames(lml)[i]))
  sub_array[i, , ] <- (wtran$power.corr)
  #print(i)
}

wdat <- sub_meas[,442]
time <-  as.numeric(rownames(ewMu4465))
wdat <- data.frame(time = time, values = wdat)
test <- wt(wdat, pad = TRUE)

for(i in 1:ncol(sub_meas)){
  wdat <- sub_meas[,i]
  time <-  as.numeric(rownames(ewMu4465))
  wdat <- data.frame(time = time, values = wdat)
  wtran <- wt(wdat, pad = TRUE)
  #plot(wtran, main = paste(colnames(lml)[i]))
  mat <- (wtran$phase) - test$phase
  for(j in 1:nrow(mat)){
    for(k in 1:ncol(mat)){
      mat[j,k] <- ((((mat[j,k])*180/3.14159265359) + 540) %% 360) - 180
    }
  }
  sub_array[i, , ] <- mat
  print(i)
}

sub_X <- as.tensor(sub_array)

sub_decom <- cp(sub_X, num_components = 4, max_iter = 50)

sub_decom$conv
sub_decom$norm_percent

candidates <- 1:ncol(ewMu4465)
#centers = cbind(sub_decom$U[[1]][big_places,1],sub_decom$U[[1]][big_places,2],sub_decom$U[[1]][big_places,3])
ob <- kmeans(cbind(sub_decom$U[[1]][,1],sub_decom$U[[1]][,2],sub_decom$U[[1]][,3], sub_decom$U[[1]][,4]),centers=3)
earlybirths <- colMeans(ewBu4464[3:5,])/colMeans(ewPu4464[3:5,])
earlybirths <- colMeans(ewBu4464[3:5,]) - colMeans(ewBu4464[9:12,])
earlybirths <- apply(ewBu4464, 2, function(c)max(c)-min(c))
plot_data <- data.frame(x=as.numeric(ewXYu4465["Long",candidates]), y = as.numeric(ewXYu4465["Lat",candidates]),
                        #dist = mindist[candidates],
                        population= colMeans(ewPu4464)[candidates],births = colMeans(ewBu4464)[candidates],
                        score1 = sub_decom$U[[1]][,1], score2 = sub_decom$U[[1]][,2], 
                        score3 = sub_decom$U[[1]][,3], score4 = sub_decom$U[[1]][,4],
                        #cluster = ob$cluster,
                        name = colnames(ewMu4465), earlybirths = earlybirths)


### component tiles ###

comp1 <- outer(period1,time1)
comp2 <- outer(period2,time2)
comp3 <- outer(period3, time3)
comp4 <- outer(period4, time4)
period <- wtran$period
rownames(comp1) <-  rownames(comp2) <- rownames(comp4) <- rownames(comp3) <- wtran$period
colnames(comp1) <- colnames(comp2) <- colnames(comp3) <- colnames(comp4) <- (time)
#total <- as.data.frame(total)
comp1_long <- melt(comp1, id.vars = c(time,period))
comp3_long <- melt(comp3, id.vars = c(time,period))
comp2_long <- melt(comp2, id.vars = c(time,period))
comp4_long <- melt(comp4, id.vars = c(time,period))

ggplot() + geom_tile(data = comp1_long, aes(y=Var1,x=Var2,fill=value)) + stat_contour() + theme_classic() +
  scale_y_log10(breaks = c(.01,.5,1,2,5,10)) + scale_fill_viridis_c(name = 'power') +
  ylab('period') + xlab('time') + theme_void() + ggtitle('First Component')

ggplot() + geom_tile(data = comp2_long, aes(y=Var1,x=Var2,fill=value)) + stat_contour() + theme_classic() +
  scale_y_log10(breaks = c(.01,.5,1,2,5,10)) + scale_fill_viridis_c(name = 'power') +
  ylab('period') + xlab('time') + theme_void() + ggtitle('Second Component')

ggplot() + geom_tile(data = comp3_long, aes(y=Var1,x=Var2,fill=value)) + stat_contour() + theme_classic() +
  scale_y_log10(breaks = c(.01,.5,1,2,5,10)) + scale_fill_viridis_c(name = 'power') +
  ylab('period') + xlab('time') + theme_void() + ggtitle('Third Component')

ggplot() + geom_tile(data = comp4_long, aes(y=Var1,x=Var2,fill=value)) + stat_contour() + theme_classic() +
  scale_y_log10(breaks = c(.01,.5,1,2,5,10)) + scale_fill_viridis_c(name = 'power') +
  ylab('period') + xlab('time') + theme_void() + ggtitle('Fourth Component')



ggplot() + geom_tile(data = comp1_long, aes(y=Var1,x=Var2,fill=value), alpha = .75) + stat_contour() + theme_classic() +
  scale_y_log10(breaks = c(.01,.5,1,2,5,10)) + ggtitle('Component 1') + scale_fill_viridis_c(name = 'power') +
  ylab('period') + xlab('time')

ggplot() + geom_tile(data = comp1_long, aes(y=Var1,x=Var2,fill=value), alpha = .75) + stat_contour() + theme_classic() +
  scale_y_log10(breaks = c(.01,.5,1,2,5,10)) + ggtitle('Component 1') + scale_fill_viridis_c(name = 'power') +
  ylab('period') + xlab('time')


names <- c('London','Manchester','Liverpool','Norwich','Leeds')
plot_list <- cwt_list <- list(length = length(names))
i=1
for(name in names){
  num <- which(colnames(ewMu4465) == name)
  cid <- sub_decom$U[[1]][num,]
  total <- cid[1]*outer(period1,time1) + cid[2]*outer(period2,time2) + cid[3]*outer(period3,time3) + cid[4]*outer(period4,time4)
  #total <- cid[1]*outer(period1,time1) + cid[2]*outer(period2,time2) 
  
  period <- wtran$period
  rownames(total) <-  wtran$period
  colnames(total) <- (time)
  #total <- as.data.frame(total)
  total_long <- melt(total, id.vars = c(time,period))
  total_long$fill <- scale(total_long$value, center=TRUE, scale = TRUE)
  
  plot_list[[i]] <- ggplot() + geom_tile(data = total_long, aes(y=Var1,x=Var2,fill=value), alpha = .75) + stat_contour() + theme_classic() +
    scale_y_log10(breaks = c(.01,.5,1,2,5,10)) + ggtitle(paste('Reconstructed',name, sep = ': ')) + scale_fill_viridis_c(name = 'power') +
    ylab('period') + xlab('time')
  
  wdat <- sub_meas[,name]
  time <-  as.numeric(rownames(ewMu4465))
  wdat <- data.frame(time = time, values = wdat)
  wtran <- wt(wdat, pad = TRUE)
  #plot(wtran, main = paste(colnames(lml)[i]))
  power <- wtran$power.corr
  rownames(power) <-  wtran$period
  colnames(power) <- (time)
  total_long <- melt(power, id.vars = c(time,period))
  total_long$fill <- scale(total_long$value, center=TRUE, scale = TRUE)
  
  cwt_list[[i]] <- ggplot() + geom_tile(data = total_long, aes(y=Var1,x=Var2,fill=value), alpha = .75) + stat_contour() + theme_classic() +
    scale_y_log10(breaks = c(.01,.5,1,2,5,10)) + ggtitle(paste('Data',name, sep = ': ')) + scale_fill_viridis_c(name = 'power') +
    ylab('period') + xlab('time')
  
  i=i+1
}

place1 <- plot_list[[1]]
place1raw <- cwt_list[[1]]
place2 <- plot_list[[2]]
place2raw <- cwt_list[[2]]
place3 <- plot_list[[3]]
place3raw <- cwt_list[[3]]
place4 <- plot_list[[4]]
place4raw <- cwt_list[[4]]
place5 <- plot_list[[5]]
place5raw <- cwt_list[[5]]

lond_plots <- plot_grid(place1,place1raw,nrow=1)
man_plots <- plot_grid(place2,place2raw, nrow = 1)
nor_plots <- plot_grid(place4, place4raw, nrow = 1)
liv_plots <- plot_grid(place3,place3raw,nrow=1)
leeds_plots <- plot_grid(place5,place5raw,nrow=1)
plot_grid(lond_plots,man_plots, nor_plots, leeds_plots, nrow = 4, labels = c('A','B','C','D'))
plot_grid(place1,place1raw,place2,place2raw, nrow = 2)
plot_grid(place3,place3raw,place5,place5raw, nrow = 2)


time_data <- data.frame(time = time, power1 = (sub_decom$U[[3]][,1]), power2 = sub_decom$U[[3]][,2], power3 = sub_decom$U[[3]][,3],
                        power4 = (sub_decom$U[[3]][,4]))
period_data <- data.frame(period = wtran$period, power1 = (sub_decom$U[[2]][,1]), power2 = sub_decom$U[[2]][,2], power3 = sub_decom$U[[2]][,3],
                          power4 = (sub_decom$U[[2]][,4]))
time_data <- data.frame(time = time, power1 = (sub_decom$U[[3]][,1]), power2 = sub_decom$U[[3]][,2])

score1 <- ggplot(data = hist_data,aes(x = (score1), group = as.factor(cluster), fill = as.factor(cluster)))+
  geom_density(alpha = .6) +
  theme_minimal() + ggtitle('Score Density: Component 1') + xlab('score') + scale_fill_manual(values = viridis(4), name = 'cluster')

score2 <- ggplot(data = hist_data,aes(x = (score2), group = as.factor(cluster), fill = as.factor(cluster)))+
  geom_density(alpha = .6) +
  theme_minimal() + ggtitle('Score Density: Component 2') + xlab('score') + scale_fill_manual(values = viridis(4), name = 'cluster')

score3 <- ggplot(data = hist_data,aes(x = (score3), group = as.factor(cluster), fill = as.factor(cluster)))+
  geom_density(alpha = .6) +
  theme_minimal() + ggtitle('Score Density: Component 3') + xlab('score') + scale_fill_manual(values = viridis(4), name = 'cluster')

score4 <- ggplot(data = hist_data,aes(x = (score4), group = as.factor(cluster), fill = as.factor(cluster)))+
  geom_density(alpha = .6) +
  theme_minimal() + ggtitle('Score Density: Component 4') + xlab('score') + scale_fill_manual(values = viridis(4), name = 'cluster')

plot_data$birth_rate <- earlybirths
score4 <- ggplot(data = plot_data,
                 aes(x = birth_rate, group = as.factor(cluster), fill = as.factor(cluster)))+ geom_density(alpha = .6) +
  theme_minimal() + ggtitle('Population Density') + xlab('log pop') + scale_fill_manual(values = viridis(4), name = 'cluster')

plot_grid(score1,score2,score3,score4,labels=c('A','B','C',"D"), nrow = 2)

pop_plot <- ggplot(data = hist_data) + geom_point(aes(x=log10(pop), y = (score), col = as.factor(close))) +
  scale_color_manual(values = c('red','dodgerblue'), name = 'close')+ theme_minimal() + xlab('log population') + ylab('log score') +
  ggtitle('Score by Population')

ggplot(data = hist_data) + geom_point(aes(x=log10(pop), y = (score2), col = as.factor(close))) +
  scale_color_manual(values = c('red','dodgerblue'), name = 'close')+ theme_minimal() + xlab('log population') + ylab('log score') +
  ggtitle('Score by Population')


time_plot_one <- ggplot(data = time_data) + geom_line(aes(x=time, y = power1)) + theme_minimal() +
  ggtitle('First Component') + ylab('power')
time_plot_two <- ggplot(data = time_data) + geom_line(aes(x=time, y = power2)) + theme_minimal() +
  ggtitle('Second Component') + ylab('power')
time_plot_three <- ggplot(data = time_data) + geom_line(aes(x=time, y = power3)) + theme_minimal() +
  ggtitle('Third Component') + ylab('power')
time_plot_four <- ggplot(data = time_data) + geom_line(aes(x=time, y = power4)) + theme_minimal() +
  ggtitle('Fourth Component') + ylab('power')

period_plot_one <- ggplot(data = period_data) + geom_line(aes(x=period, y = power1)) + theme_minimal() +
  geom_vline(xintercept = 1, col = 'red', lty = 3) +
  geom_vline(xintercept = 2, col = 'red', lty = 3) +
  geom_hline(yintercept = 0, col = 'dodgerblue',lty=2) + ylab('power')
period_plot_two <- ggplot(data = period_data) + geom_line(aes(x=period, y = power2)) + theme_minimal() +
  geom_vline(xintercept = 1, col = 'red',lty =3) +
  geom_vline(xintercept = 2, col = 'red',lty=3) + 
  geom_hline(yintercept = 0, col = 'dodgerblue',lty=2) + ylab('power')
period_plot_three <- ggplot(data = period_data) + geom_line(aes(x=period, y = power3)) +
  theme_minimal() + geom_vline(xintercept = 1, col = 'red', lty=3) +
  geom_vline(xintercept = 2, col = 'red',lty=3) + 
  geom_hline(yintercept = 0, col = 'dodgerblue',lty=2) + ylab('power')
period_plot_four <- ggplot(data = period_data) + geom_line(aes(x=period, y = power4)) +
  theme_minimal() + geom_vline(xintercept = 1, col = 'red', lty=3) +
  geom_vline(xintercept = 2, col = 'red',lty=3) + 
  geom_hline(yintercept = 0, col = 'dodgerblue',lty=2) + ylab('power')

pop_plot_one <- ggplot(data=plot_data) + geom_point(aes(x=log10(population), y= (score1)), alpha = .75) + theme_minimal() +
  ylab('score') + ggtitle('')
pop_plot_two <- ggplot(data=plot_data) + geom_point(aes(x=log10(population), y= score2), alpha = .75) + theme_minimal() +
  ylab('score') + ggtitle('')
pop_plot_three <- ggplot(data=plot_data) + geom_point(aes(x=log10(population), y= score3), alpha = .5) + theme_minimal() +
  ylab('score') + ggtitle('')
pop_plot_four <- ggplot(data=plot_data) + geom_point(aes(x=log10(population), y= score4), alpha = .5) + theme_minimal() +
  ylab('score') + ggtitle('')


times <- plot_grid(time_plot_one, time_plot_two, time_plot_three, time_plot_four, nrow = 1)
periods <- plot_grid(period_plot_one, period_plot_two, period_plot_three, period_plot_four, nrow = 1)
pop_plots <-  plot_grid(pop_plot_one, pop_plot_two, pop_plot_three, pop_plot_four, nrow = 1)
time_period_sub <- plot_grid(times, periods, pop_plots, nrow = 3, labels = c('A','B','C'))

times <- plot_grid(time_plot_one, time_plot_two, nrow = 1)
periods <- plot_grid(period_plot_one, period_plot_two, nrow = 1)
time_period_sub <- plot_grid(times, periods, pop_plots, nrow = 3, labels = c('A','B','C'))

box_data <- data.frame(pos = ifelse(z>0,'pos','neg'), score = z, pop = colMeans(ewPu4464),
                       births = y)

box_data %>%
  filter(pop > 300000) %>%
  ggplot() + geom_boxplot(aes(x=pos,y=births)) + ggtitle('Baby Boom Birth Rates') + xlab('score') + ylab('crude birth rate') +
  scale_x_discrete(labels = c('negative','positive')) + theme_classic()

## other paper images
logmeas <- log(ewMu4465+1)
#logmeas <- log(unname(apply(ewMl4494,1,sum)))
wdat <- logmeas
wdat <- logmeas[,"London"]
time <-  as.numeric(rownames(logmeas))
wdat <- data.frame(time = time, values = wdat)
wave <- wt(wdat, pad = TRUE)

wave_res <- wave$signif
wave_power <- wave$power.corr
#wave_power <- wave_power/max(wave_power)
#wave_power <- (apply(wave_power, 2, function(c)c/(sum(c)*wave$scale)))
rownames(wave_res) <- rownames(wave_power) <- wave$period
colnames(wave_res) <- colnames(wave_power) <- time
wave_long <- melt(wave_res, id.vars = c(time,period))
wave_long_power <- melt(wave_power, id.vars = c(time,period))
wave_long$pval <- ifelse(wave_long$value >= 1, 1, 0)

wave_long$power <- wave_long_power$value
wave_long$fill <- ifelse(wave_long$pval == 1, NA, (wave_long$power))
wave_long$fill <- cut(log10(wave_long$fill), 9)
colfun <- colorRampPalette(c('dodgerblue','gold'))
colscale <- colfun(9)
colscale <- (c(brewer.pal(9, "Blues")))

ggplot() + geom_tile(data = wave_long, aes(y=Var1,x=Var2,fill=fill), alpha = .75) + stat_contour() + theme_classic() +
  scale_y_log10(breaks = c(.01,.5,1,2,5,10)) +
  geom_line(aes(x=time[which(log10(wave$coi)>-1.4)],y=(wave$coi)[which(log10(wave$coi)>-1.4)]),lwd=1)+
  scale_x_continuous(breaks = c(50,60,70,80,90), labels = 1900+c(50,60,70,80,90)) +
  scale_fill_manual("Power", values = (colscale), na.value = 'red',
                    labels = c('(-2.04,-0.991]','(-0.991,0.0509]','(0.0509,1.09]','(1.09,2.14]','(2.14,3.18]','(3.18,4.22]]',
                               '(4.22,5.26]','(5.26,6.31]','(6.31,7.36]','significant')) + ylab('Period') + xlab('Time') +
  theme(legend.position = "bottom") 

ggplot() + geom_tile(data = wave_long, aes(y=Var1,x=Var2,fill=fill), alpha = .75) + stat_contour() + theme_classic() +
  scale_y_log10(breaks = c(.01,.5,1,2,5,10)) +
  #geom_line(aes(x=time[which(log10(wave$coi)>-1.4)],y=(wave$coi)[which(log10(wave$coi)>-1.4)]),lwd=1)+
  scale_x_continuous(breaks = c(45,50,55,60,65), labels = 1900+c(45,50,55,60,65)) +
  scale_fill_manual("Power", values = (colscale), na.value = 'red',
                    labels = c('(-2.04,-0.991]','(-0.991,0.0509]','(0.0509,1.09]','(1.09,2.14]','(2.14,3.18]','(3.18,4.22]]',
                               '(4.22,5.26]','(5.26,6.31]','(6.31,7.36]','significant')) + ylab('Period') + xlab('Time') +
  theme_void() + theme(legend.position="none")


## London and Norwich

time1 <- sub_decom$U[[3]][,1]
period1 <- sub_decom$U[[2]][,1]
time2<- sub_decom$U[[3]][,2]
period2 <- sub_decom$U[[2]][,2]

names <- c('London','Norwich')
plot_list <- cwt_list <- list(length = length(names))
i=1
for(name in names){
  num <- which(names == name)
  cid <- sub_decom$U[[1]][num,]
  total <- cid[1]*outer(period1,time1) + cid[2]*outer(period2,time2)
  #total <- cid[1]*outer(period1,time1) + cid[2]*outer(period2,time2) 
  
  period <- wtran$period
  rownames(total) <-  wtran$period
  colnames(total) <- (time)
  #total <- as.data.frame(total)
  total_long <- melt(total, id.vars = c(time,period))
  total_long$fill <- scale(total_long$value, center=TRUE, scale = TRUE)
  
  plot_list[[i]] <- ggplot() + geom_tile(data = total_long, aes(y=Var1,x=Var2,fill=value), alpha = .75) + stat_contour() + theme_classic() +
    scale_y_log10(breaks = c(.01,.5,1,2,5,10)) + ggtitle(paste('Reconstructed',name, sep = ': ')) + scale_fill_viridis_c(name = 'power') +
    ylab('period') + xlab('time')
  
  wdat <- sub_meas[,name]
  time <-  as.numeric(rownames(ewMu4465))
  wdat <- data.frame(time = time, values = wdat)
  wtran <- wt(wdat, pad = TRUE)
  #plot(wtran, main = paste(colnames(lml)[i]))
  power <- wtran$power.corr
  rownames(power) <-  wtran$period
  colnames(power) <- (time)
  total_long <- melt(power, id.vars = c(time,period))
  total_long$fill <- scale(total_long$value, center=TRUE, scale = TRUE)
  
  cwt_list[[i]] <- ggplot() + geom_tile(data = total_long, aes(y=Var1,x=Var2,fill=value), alpha = .75) + stat_contour() + theme_classic() +
    scale_y_log10(breaks = c(.01,.5,1,2,5,10)) + ggtitle(paste('Data',name, sep = ': ')) + scale_fill_viridis_c(name = 'power') +
    ylab('period') + xlab('time')
  
  i=i+1
}

place1 <- plot_list[[1]]
place1raw <- cwt_list[[1]]
place2 <- plot_list[[2]]
place2raw <- cwt_list[[2]]

pdf(height = 8.5, width = 11, file = 'lon_nor.pdf')
plot_grid(place1, place1raw, place2, place2raw, nrow = 2, labels = c('A','B','C','D'))
dev.off()

