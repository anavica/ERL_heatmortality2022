
################################################################################
################################################################################
####### - "The footprint of human-induced climate change on ####################
####### heat-related deaths in the summer of 2022 in Switzerland"  #############
#######  By Vicedo-Cabrera et al. published in ERL #############################
################################################################################
################################################################################


# This code is to partially reproduce the analysis shown in the publication.
# This analysis aims at quantifying the heat-related mortality attributed to 
# climate change during the Summer of 2022 in Switzerland. We perform the analysis
# in 2 stages:    
# Step 1  - estimation of exposure-response (ER) function by Canton, sex & age (</>=65 years) 
#    using observed temperature-mortality (all cause) between 1990-2017.
# Step 2 - quantification of heat-related mortality by Canton in Summer 2022 (June-August) 
#   using observed temperature data (factual) and the counterfactual scenarios


# Note: Step (1) cannot be reproduced with this code since the original data to compute
# the ER functions cannot be distributed due to confidentiality issues, 
# in particular for small cantons with low counts per day/sex/age category.
# Thus, we provide the output from this step 1 separately that will be needed to 
# performed the two other steps ("output1step.Rdata).  

# Data needed (daily temperature-mortality in Summer 2022) can be found here: 
# https://doi.org/10.48620/315 

# LOAD PACKAGES
library(tidyverse); library(dlnm)
library(MASS)


################################################################################
# STEP 2
################################################################################

# LOAD TEMPERATURE-MORTALITY SUMMER 2022 DATA
data2022 <- read_csv("dataCHsummer2022.csv")
listcantonsname <- unique(data2022$canton_name)

# CREATE LIST DATA BY CANTON
datalist2022 <- lapply(listcantonsname, function(listcantonsname) 
  data2022[data2022$canton_name==listcantonsname,])
names(datalist2022) <- listcantonsname
rm(data2022)

# LOAD NECESSARY OBJECTS FROM STEP 1
load("output1step.Rdata")

# "output1step" contains:
# - "blupall.list" = list of lists of coef/vcov defining the ER by [[sex:age]][[canton]]
#     categories of sex:age --> 1: male_beloweq65 // 2: male_above65 
#                            // 3: female_beloweq65 // 4: female_above65
# - "datalist2022" = list of time series daily temperature-mortality (sex/age) 
#     by canton in Summer 2022. indsex: 1 male, 2 female // indage: 1 0-65, 2 above 65 years.
# - mintempcanton_list = list of lists of MMTT by [[sex:age]][[canton]]
#     categories of sex:age --> 1: male_beloweq65 // 2: male_above65 
#                            // 3: female_beloweq65 // 4: female_above65
# - predpoollist = list of coef/vcov defining the pooled ER functions by sex and age 
#     categories of sex:age --> 1: male_beloweq65 // 2: male_above65 
#                            // 3: female_beloweq65 // 4: female_above65
# - tmeandist_canton = distribution Tmean 1990-2017 by canton
# - tmeandist_all = overall Tmean 1990-2017 distribution 

                           
# SET NUMBER OF SIMULATIONS
nsim <- 1000

# DEFINE COMBINATIONS SEX/AGE
sexage <- rbind(c(1,1),c(1,2),c(2,1),c(2,2))

# CREATE EMPTY ARRAY TO STORE RESULTS
ancantonsim_attr <- array(0,dim=list(length(listcantonsname),4,5,2,nsim+1),
  dimnames=list(listcantonsname,c("male<=65y","male>65y","female<=65y","female>65y"),
                c("obs","mod1","mod2","mod3","mod4"),c("abs","dif"),
                c("est",paste0("sim",seq(nsim)))))

# CREATE EMPTY MATRIX TO STORE AN PER DAY AND SERIES
dailyan_attr <- array(NA,dim=list(length(datalist2022),4,5,2,
                             length(unique(as.character(datalist2022[[1]]$Date)))),
                 dimnames=list(listcantonsname,
                              c("male<=65y","male>65y","female<=65y","female>65y"),
                              c("obs","mod1","mod2","mod3","mod4"),c("abs","dif"),
                              unique(as.character(datalist2022[[1]]$Date))))
 
# FACTORS TO DEFINE SCENARIOS (FACTUAL + 4 COUNTERFACTUALS)
scen_fact <- c(0,2.16,2.75,1.19,2.27)

# RUN ACROSS CANTONS 
for (i in seq(length(datalist2022))){
  
  # SELECT DATA - CANTON
  data <- data.frame(datalist2022[[i]])

    # RUN ACROSS IND SUBGROUPS
    for (j in 1:4){
    
        # SELECT DATA
        datasel <-  data[data$indsex==sexage[j,1] & data$indage==sexage[j,2],] 
  
        # EXTRACT PARAMETERS
        coef <- blupall.list[[j]][[i]]$blup
        vcov <- blupall.list[[j]][[i]]$vcov
    
        # DEFINE ARGVAR
        argvar <- list(fun="ns",knots=tmeandist_canton[c("50.0%","90.0%"),i], 
                 Bound=tmeandist_canton[c("0.0%","100.0%"),i])

        for (g in seq(length(scen_fact))){

          # DEFINE THE TMEAN (SUBSTRACT FACTOR SCENARIO)
          tmeanmod <- datasel[,"tmean"] - scen_fact[g]
          
          # DERIVE THE CENTERED BASIS
          bvar <- do.call(onebasis,c(list(x=tmeanmod),argvar))
          cenvec <- do.call(onebasis,c(list(x=mintempcanton_list[[j]][[i]]),argvar))
          bvarcen <- scale(bvar,center=cenvec,scale=F)
      
          # INDICATOR FOR HEAT DAYS
          indheat <- tmeanmod>mintempcanton_list[[j]][[i]]
      
          # COMPUTE THE DAILY CONTRIBUTIONS OF ATTRIBUTABLE DEATHS
          an <- (1-exp(-bvarcen%*%coef))*datasel$deaths

          # STORE DAILY DEATHS
          dailyan_attr[i,j,g,1,] <- ifelse(indheat,an,0)

          # SUM 
          ancantonsim_attr[i,j,g,1,1] <- sum(an[indheat], na.rm=T)

          # SAMPLE COEF ASSUMING A MULTIVARIATE NORMAL DISTRIBUTION
          set.seed(13041975)
          coefsim <- mvrnorm(nsim,coef,vcov)
        
          # LOOP ACROSS ITERATIONS
          for(s in seq(nsim)) {
          
            # COMPUTE THE DAILY CONTRIBUTIONS OF ATTRIBUTABLE DEATHS
            ansim <- (1-exp(-bvarcen%*%coefsim[s,]))*datasel$deaths
                
            # COMPUTE THE RESULTS FOR EACH RANGE AND PERIOD AND SUM
            # NB: ACCOUNT FOR NO TEMPERATURE BELOW/ABOVE CEN FOR GIVEN PERIODS
            ancantonsim_attr[i,j,g,1,s+1] <- sum(ansim[indheat], na.rm=T)
          }
    }
  }
}

# COMPUTE DENOMINATORS
ndeaths_sum_total <- sapply(datalist2022, function(x) sum(x$deaths))
ndeaths_sum_malebelow65 <- sapply(datalist2022, function(x) sum(x$deaths[x$indsex==1 & x$indage==1]))
ndeaths_sum_femalebelow65 <- sapply(datalist2022, function(x) sum(x$deaths[x$indsex==2 & x$indage==1]))
ndeaths_sum_maleabove65 <- sapply(datalist2022, function(x) sum(x$deaths[x$indsex==1 & x$indage==2]))
ndeaths_sum_femaleabove65 <- sapply(datalist2022, function(x) sum(x$deaths[x$indsex==2 & x$indage==2]))


# ESTIMATE DIFFERENCE FACTUAL-COUNTERFACT
ancantonsim_attr[,,,2,] <- ancantonsim_attr[,,rep(1,length(scen_fact)),1,] - ancantonsim_attr[,,,1,]

# GET SUMMARIES
antot_attr <- aftot_attr <- array(NA, dim=list(2,2,3),
                                                 dimnames=list(c("abs","dif"),
                                                            c("fact","counterfact"),
                                                            c("est","ci.l","ci.u")))
ancanton_attr <- afcanton_attr <- array(NA, dim=list(length(listcantonsname),2,2,3),
                                                 dimnames=list(listcantonsname,c("abs","dif"),
                                                            c("fact","counterfact"),
                                                            c("est","ci.l","ci.u")))

ansexage_attr <- afsexage_attr <- 
  array(NA, dim=list(8,2,2,3),dimnames=list(c("male<=65y", "male>65y", "female<=65y",
            "female>65y","male","female", "<=65y", ">65y" ),c("abs","dif"),
            c("fact","counterfact"),c("est","ci.l","ci.u")))


# TOTAL AN / AF
antot_attr[,1,1] <- apply(ancantonsim_attr[,,1,,1],3,sum)
antot_attr[,1,2] <- apply(apply(ancantonsim_attr[,,1,,-1],c(3,4),sum),1,quantile,0.025)
antot_attr[,1,3] <- apply(apply(ancantonsim_attr[,,1,,-1],c(3,4),sum),1,quantile,0.975)

antot_attr[,2,1] <- apply(apply(ancantonsim_attr[,,-1,,1],c(3,4),sum),2,mean)
antot_attr[,2,2] <- apply(apply(ancantonsim_attr[,,-1,,-1],c(3,4,5),sum),2,quantile,0.025)
antot_attr[,2,3] <- apply(apply(ancantonsim_attr[,,-1,,-1],c(3,4,5),sum),2,quantile,0.975)

aftot_attr <- (antot_attr/sum(ndeaths_sum_total))*100
  
# AN / AF BY CANTON
ancanton_attr[,,1,1] <- apply(ancantonsim_attr[,,1,,1], c(1,3),sum)
ancanton_attr[,,1,2] <- apply(apply(ancantonsim_attr[,,1,,-1], c(1,3,4),sum),c(1,2) ,quantile,0.025)
ancanton_attr[,,1,3] <- apply(apply(ancantonsim_attr[,,1,,-1], c(1,3,4),sum),c(1,2) ,quantile,0.975)

ancanton_attr[,,2,1] <- apply(apply(ancantonsim_attr[,,-1,,1], c(1,3,4),sum), c(1,3),mean)
ancanton_attr[,,2,2] <- apply(apply(ancantonsim_attr[,,-1,,-1], c(1,3,4,5),sum), c(1,3),quantile,0.025)
ancanton_attr[,,2,3] <- apply(apply(ancantonsim_attr[,,-1,,-1], c(1,3,4,5),sum), c(1,3),quantile,0.975)

afcanton_attr <- (ancanton_attr/ndeaths_sum_total)*100

# AN / AF BY SEX/AGE
ansexage_attr[c(1:4),,1,1] <- apply(ancantonsim_attr[,,1,,1], c(2,3),sum)
ansexage_attr[c(1:4),,1,2] <- apply(apply(ancantonsim_attr[,,1,,-1], c(2,3,4),sum),c(1,2) ,quantile,0.025)
ansexage_attr[c(1:4),,1,3] <- apply(apply(ancantonsim_attr[,,1,,-1], c(2,3,4),sum),c(1,2) ,quantile,0.975)

ansexage_attr[c(1:4),,2,1] <- apply(apply(ancantonsim_attr[,,-1,,1], c(2,3,4),sum), c(1,3),mean)
ansexage_attr[c(1:4),,2,2] <- apply(apply(ancantonsim_attr[,,-1,,-1], c(2,3,4,5),sum),c(1,3) ,quantile,0.025)
ansexage_attr[c(1:4),,2,3] <- apply(apply(ancantonsim_attr[,,-1,,-1], c(2,3,4,5),sum),c(1,3) ,quantile,0.975)


ansexage_attr[5,,1,1] <- apply(ancantonsim_attr[,c(1,2),1,,1], c(3),sum)
ansexage_attr[5,,1,2] <- apply(apply(ancantonsim_attr[,c(1,2),1,,-1], c(3,4),sum), 1,quantile,0.025)
ansexage_attr[5,,1,3] <- apply(apply(ancantonsim_attr[,c(1,2),1,,-1], c(3,4),sum), 1,quantile,0.975)

ansexage_attr[6,,1,1] <- apply(ancantonsim_attr[,c(3,4),1,,1], c(3),sum)
ansexage_attr[6,,1,2] <- apply(apply(ancantonsim_attr[,c(3,4),1,,-1], c(3,4),sum), 1,quantile,0.025)
ansexage_attr[6,,1,3] <- apply(apply(ancantonsim_attr[,c(3,4),1,,-1], c(3,4),sum), 1,quantile,0.975)

ansexage_attr[7,,1,1] <- apply(ancantonsim_attr[,c(1,3),1,,1], c(3),sum)
ansexage_attr[7,,1,2] <- apply(apply(ancantonsim_attr[,c(1,3),1,,-1], c(3,4),sum), 1,quantile,0.025)
ansexage_attr[7,,1,3] <- apply(apply(ancantonsim_attr[,c(1,3),1,,-1], c(3,4),sum), 1,quantile,0.975)

ansexage_attr[8,,1,1] <- apply(ancantonsim_attr[,c(2,4),1,,1], c(3),sum)
ansexage_attr[8,,1,2] <- apply(apply(ancantonsim_attr[,c(2,4),1,,-1], c(3,4),sum), 1,quantile,0.025)
ansexage_attr[8,,1,3] <- apply(apply(ancantonsim_attr[,c(2,4),1,,-1], c(3,4),sum), 1,quantile,0.975)

ansexage_attr[5,,2,1] <- apply(apply(ancantonsim_attr[,c(1,2),-1,,1], c(3,4),sum),2,mean)
ansexage_attr[5,,2,2] <- apply(apply(ancantonsim_attr[,c(1,2),-1,,-1], c(3,4,5),sum), 2,quantile,0.025)
ansexage_attr[5,,2,3] <- apply(apply(ancantonsim_attr[,c(1,2),-1,,-1], c(3,4,5),sum), 2,quantile,0.975)

ansexage_attr[6,,2,1] <- apply(apply(ancantonsim_attr[,c(3,4),-1,,1], c(3,4),sum),2,mean)
ansexage_attr[6,,2,2] <- apply(apply(ancantonsim_attr[,c(3,4),-1,,-1], c(3,4,5),sum), 2,quantile,0.025)
ansexage_attr[6,,2,3] <- apply(apply(ancantonsim_attr[,c(3,4),-1,,-1], c(3,4,5),sum), 2,quantile,0.975)

ansexage_attr[7,,2,1] <- apply(apply(ancantonsim_attr[,c(1,3),-1,,1], c(3,4),sum),2,mean)
ansexage_attr[7,,2,2] <- apply(apply(ancantonsim_attr[,c(1,3),-1,,-1], c(3,4,5),sum), 2,quantile,0.025)
ansexage_attr[7,,2,3] <- apply(apply(ancantonsim_attr[,c(1,3),-1,,-1], c(3,4,5),sum), 2,quantile,0.975)

ansexage_attr[8,,2,1] <- apply(apply(ancantonsim_attr[,c(2,4),-1,,1], c(3,4),sum),2,mean)
ansexage_attr[8,,2,2] <- apply(apply(ancantonsim_attr[,c(2,4),-1,,-1], c(3,4,5),sum), 2,quantile,0.025)
ansexage_attr[8,,2,3] <- apply(apply(ancantonsim_attr[,c(2,4),-1,,-1], c(3,4,5),sum), 2,quantile,0.975)


ndeaths_sum_agesex <- datalist2022 %>% flatten_dfc() %>% 
  dplyr::select(starts_with("deaths")) %>% rowwise %>% 
  summarise(deaths=sum(c_across(everything())))

ndeaths_sum_agesex$indsex <- datalist2022[[1]]$indsex
ndeaths_sum_agesex$indage <- datalist2022[[1]]$indage
ndeaths_sum_agesex$date <- datalist2022[[1]]$Date

totaldeaths <- ndeaths_sum_agesex %>% group_by(indsex,indage) %>% summarize(sum=sum(deaths))

afsexage_attr[c(1,2),,,] <- (ansexage_attr[c(1,2),,,]/totaldeaths$sum[totaldeaths$indsex==1])*100
afsexage_attr[c(3,4),,,] <- (ansexage_attr[c(3,4),,,]/totaldeaths$sum[totaldeaths$indsex==2])*100
afsexage_attr[c(5),,,] <- (ansexage_attr[c(5),,,]/sum(totaldeaths$sum[totaldeaths$indsex==1]))*100
afsexage_attr[c(6),,,] <- (ansexage_attr[c(6),,,]/sum(totaldeaths$sum[totaldeaths$indsex==2]))*100
afsexage_attr[c(7),,,] <- (ansexage_attr[c(7),,,]/sum(totaldeaths$sum[totaldeaths$indage==1]))*100
afsexage_attr[c(8),,,] <- (ansexage_attr[c(8),,,]/sum(totaldeaths$sum[totaldeaths$indage==2]))*100


# PERCENTAGE
(ansexage_attr[,"dif","counterfact",1] / ansexage_attr[,"abs","fact",1])*100

# DAILY
dailyan_attr[,,,2,] <- dailyan_attr[,,rep(1,length(scen_fact)),1,] - dailyan_attr[,,,1,]

dailyan_attr_res <- array(NA,dim=list(length(dailyan_attr[1,1,1,2,]),2,2), 
                          dimnames=list(dimnames(dailyan_attr)[[5]],c("abs","dif"),
                                                            c("fact","counterfact")))
  
dailyan_attr_res[,,1] <- t(apply(dailyan_attr[,,1,,],c(3,4),sum))
dailyan_attr_res[,,2] <- t(apply(apply(dailyan_attr[,,-1,,],c(3,4,5),sum),c(2,3),mean))


################################################################################
#### FIGURES
################################################################################

# FIGURE 1 POOLED-BLUPS ER BY SEX/AGE
rangeplot <- c(6:114)
maintitle <- c("Male 0-65 years", "Male >65 years", "Female 0-65 years", "Female >65 years")
col <- c("blue","navy","red","maroon")

pdf("erc_poolblups.pdf",width=8,height=7)
layout(matrix(1:4,nrow=2),heights=c(rep(0.5,2)))

for (j in c(1:4)){
  
  predtoplot <- predpoollist[[j]]
  predblup <- blupall.list[[j]]

  par(mar=c(4,4,2,2),mgp=c(3,1,0))

  plot(predtoplot$predvar[rangeplot], predtoplot$allRRfit[rangeplot],
     ylim=c(0.8,1.8),xlab="", 
       ylab="", main=maintitle[j], cex.axis=0.8,cex.lab=0.8,
        col="white", xaxt="n")
     polygon(c(predtoplot$predvar[rangeplot], rev(predtoplot$predvar[rangeplot])), 
            c(predtoplot$allRRlow[rangeplot],rev(predtoplot$allRRhigh[rangeplot])), 
            density = 200, col =alpha(col[j],0.1))
     lines(predtoplot$predvar[rangeplot], predtoplot$allRRfit[rangeplot],
           lwd=2, col=col[j])
     axis(1, at=tmeandist_all[paste0(c(1,25,50,75,99),".0%")], 
        labels=c(1,25,50,75,99), cex.axis=0.8)
     mtext("Relative risk",side=2, line=2.2, cex=0.8)
     mtext("Temperature Percentile",side=1, line=2.2, cex=0.8)

    for(i in seq(length(listcantonsname))) {

        argvar <- list(x=tmeandist_all,fun="ns",
                 knots=tmeandist_all[paste0(c(50,90), ".0%")],
                 Bound=tmeandist_all[paste0(c(0,100), ".0%")])
        bvar <- do.call(onebasis,argvar)
        
        mmt <- tmeandist_all[(35:100)[which.min((bvar%*%predblup[[i]]$blup)[35:100])]]
  
        predplot <- crosspred(bvar,coef=predblup[[i]]$blup,vcov=predblup[[i]]$vcov,
                      model.link="log",at=tmeandist_all,cen=mmt)
        lines(predplot$predvar[rangeplot],predplot$allRRfit[rangeplot],
              col=alpha(col[j],0.5),lwd=1, lty="dotted")
        abline(h=1)
    }
}     
dev.off()


# FIGURE SUPPL (adapted) - ER CURVES BY SEX/AGE/CANTON 

pdf("erc_sexagecantons.pdf",width=9,height=13)
layout(matrix(seq(5*3),nrow=5,byrow=T))
par(mar=c(4,3.8,3,2.4),mgp=c(2.5,1,0),las=1)

for(i in seq(length(listcantonsname))) {
  predvar <- tmeandist_canton[,i]
  perc99 <- tmeandist_canton["99.0%",i]
    
  # REDEFINE THE FUNCTION USING ALL THE ARGUMENTS (BOUNDARY KNOTS INCLUDED)
  argvar <- list(x=predvar,fun="ns",
                 knots=tmeandist_canton[c("50.0%","90.0%"),i],
                 Bound=tmeandist_canton[c("0.0%","100.0%"),i])
  bvar <- do.call(onebasis,argvar)
  
  predblup11 <- crosspred(bvar,coef=blupall.list[[1]][[i]]$blup,vcov=blupall.list[[1]][[i]]$vcov,
                    model.link="log",by=0.1,cen=mintempcanton_list[[1]][i])
  predblup12 <- crosspred(bvar,coef=blupall.list[[2]][[i]]$blup,vcov=blupall.list[[2]][[i]]$vcov,
                    model.link="log",by=0.1,cen=mintempcanton_list[[2]][i])
  predblup21 <- crosspred(bvar,coef=blupall.list[[3]][[i]]$blup,vcov=blupall.list[[3]][[i]]$vcov,
                    model.link="log",by=0.1,cen=mintempcanton_list[[3]][i])
  predblup22 <- crosspred(bvar,coef=blupall.list[[4]][[i]]$blup,vcov=blupall.list[[4]][[i]]$vcov,
                    model.link="log",by=0.1,cen=mintempcanton_list[[4]][i])
  
  plot(predblup11,ylim=c(0,2.5),lab=c(6,5,7),xlab="Mean Temperature",
       ylab="Relative Risk", main=listcantonsname[i], col=col[1],
       ci.arg=list(col=alpha(col[1],0.1)),lwd=1.5)
  
  lines(predblup12$predvar,predblup12$allRRfit,col=col[2],lwd=1.5)
  polygon(c(rev(predblup12$predvar), predblup12$predvar), 
          c(rev(predblup12$allRRhigh), predblup12$allRRlow), 
          col = alpha(col[2],0.1), border = NA)
  
  lines(predblup21$predvar,predblup21$allRRfit,col=col[3],lwd=1.5)
  polygon(c(rev(predblup21$predvar), predblup21$predvar), 
          c(rev(predblup21$allRRhigh), predblup21$allRRlow), 
          col = alpha(col[3],0.1), border = NA)
  
  lines(predblup22$predvar,predblup22$allRRfit,col=col[4],lwd=1.5)
  polygon(c(rev(predblup22$predvar), predblup22$predvar), 
          c(rev(predblup22$allRRhigh), predblup22$allRRlow), 
          col = alpha(col[4],0.1), border = NA)

  abline(v=perc99, lty="dashed", col="grey")

}
dev.off()



# FIGURE 3 - TIME SERIES DAILY HEAT-MORTALITY & TMEAN

# BUILD DATASET FOR THE PLOT
date2022 <- seq(as.Date("2022/06/01"),as.Date("2022/08/31"), "days")
tmean_mean <- datalist2022 %>% flatten_dfc() %>% 
  dplyr::select(starts_with("tmean")) %>% rowwise %>% 
  summarise(aver_tmean=mean(c_across(everything())))

data_sum2022 <- data.frame(cbind(date2022,dailyan_attr_res[,,1],tmean_mean))

names(data_sum2022) <- c("date", "an", "tmean")

factscale <- 0.5

excess2022 <- ggplot(data=data_sum2022, aes(x=date, y=an)) +
  geom_segment(aes(x = date, xend = date,
            y = 0, yend = an),color = "lightgray", linewidth=1) +
  geom_point( aes(color=an),size =2, alpha=0.8) +
  geom_line(aes(y=tmean/factscale), col=alpha("black",0.7), lwd=0.6) +
  scale_color_gradient(low = "orange",
                            high = "red") +
  scale_x_date(date_labels = "%B") +
  scale_y_continuous(name = "Daily Heat-related Deaths",limits = c(0,50),
        sec.axis = sec_axis( trans=~.*factscale, name="Daily Mean Temperature (ºC)")) +
  geom_hline(yintercept = 0, color = "black", size=0.8)+
  theme_minimal() +
  theme(legend.position="none",
    axis.text.x = element_text(size=10),    
    axis.text.y = element_text(size=10),
    axis.title.y = element_text(size=12, vjust = 1.2),
    axis.title.x = element_blank(),
    axis.ticks = element_blank())

pdf("anexcess2022.pdf",width=10,height=6)
excess2022
dev.off()


# FIGURE 4 - BAR PLOT
 
res.plot <- matrix(NA, ncol=9, nrow=31)

res.plot[,c(1,4,7)] <- cbind(rbind(aftot_attr["abs",,"est"],
                        afsexage_attr[c(5:8),"abs",,"est"],
                        afcanton_attr[,"abs",,"est"]),
                        c(aftot_attr["dif","counterfact","est"],
                        afsexage_attr[c(5:8),"dif","counterfact","est"],
                        afcanton_attr[,"dif","counterfact","est"]))

res.plot[,c(2,5,8)] <- cbind(rbind(aftot_attr["abs",,"ci.l"],
                       afsexage_attr[c(5:8),"abs",,"ci.l"],
                        afcanton_attr[,"abs",,"ci.l"]),
                       c(aftot_attr["dif","counterfact","ci.l"],
                       afsexage_attr[c(5:8),"dif","counterfact","ci.l"],
                        afcanton_attr[,"dif","counterfact","ci.l"]))

res.plot[,c(3,6,9)] <- cbind(rbind(aftot_attr["abs",,"ci.u"],
                       afsexage_attr[c(5:8),"abs",,"ci.u"],
                        afcanton_attr[,"abs",,"ci.u"]),
                       c(aftot_attr["dif","counterfact","ci.u"],
                       afsexage_attr[c(5:8),"dif","counterfact","ci.u"],
                        afcanton_attr[,"dif","counterfact","ci.u"]))

colnames(res.plot) <- c("aaffact", "affact.l", "affact.u", 
                        "afcount", "afcount.l", "afcount.u",
                        "afdif", "afdif.l", "afdif.u")

res.plot <- rbind(res.plot[1,], rep(NA,9),
                  res.plot[c(2,3),],rep(NA,9),
                  res.plot[c(4,5),],rep(NA,9),
                  res.plot[c(6:31),])

res.plot <- data.frame(res.plot)
res.plot$idnum <- seq(1:nrow(res.plot))

res.plot$label <- c("Total", NA,"Male", "Female",NA, "0-65 years", ">65 years", NA,
                        listcantonsname)

res.plot.long1 <- res.plot[,c(1,4,7)] %>% pivot_longer(cols= aaffact:afdif, names_to = "scen", values_to = "af")
res.plot.long2 <- res.plot[,c(2,5,8)] %>% pivot_longer(cols= affact.l:afdif.l, names_to = "scen", values_to = "af.l")
res.plot.long3 <- res.plot[,c(3,6,9)] %>% pivot_longer(cols= affact.u:afdif.u, names_to = "scen", values_to = "af.u")

res.plot.long.comp <- cbind(res.plot.long1,res.plot.long2[,2],res.plot.long3[,2])
res.plot.long.comp$idnum <- c(1,2,3,4,4,4,
                             5,6,7,8,9,10,11,11,11,
                             12,13,14,15,16,17,18,18,18,
                             seq(19,96))

res.plot.long.comp$label <- NA
idnumlab <- which(res.plot.long.comp$scen=="afcount")

res.plot.long.comp$label[idnumlab] <- c("Total", NA,"Male", "Female",NA, "0-65 years", ">65 years", NA,
                        listcantonsname)

# PLOT VERSION 4
res.plot.main <- res.plot.long.comp[c(1:21),]

res.plot.main$idnum <- c(1,2,3,4,5,5,6,7,8,10,11,12,
                         13,14,14,15,16,17,
                         19,20,21)

barplot <- ggplot(res.plot.main, aes(x=idnum)) + 
      geom_crossbar(aes(x = idnum, y = 0, ymin = 0, ymax = af,fill=scen)
                    ,fatten = 0, colour = NA,width=1.5) +
      geom_errorbar(aes(ymin = af.l, ymax=af.u, x=idnum),
            width = 0.1, size=0.5) +
      geom_point(data=res.plot.main[res.plot.main$scen=="afdif",],
                 aes(y=af), fill="black",size = 4, shape = 21, 
             col = "white", show.legend = T) +
      scale_fill_manual(values = c("red", "orange", alpha("white",0)),
                        labels=c("Observed", "Counterfactual", "Attr. Climate change")) +
      guides(fill = guide_legend(override.aes = list(size=c(0,0,5)))) + 
      labs(x = "", y = "Heat-related Mortality (%)",
                         title = "") +
      geom_hline(yintercept=00, linetype="solid", color = "black") +
      scale_y_continuous("Heat-related Mortality (%)", breaks = seq(-1,8, 1),
                     limits = c(-1.7,8)) +
      geom_text(aes(x=idnum, y=-1.7, label=label),
                            alpha=0.6, size=4.5, col="black",
            inherit.aes = FALSE) +
      theme_minimal() +
      theme(legend.position="bottom",
            legend.title = element_blank(),
            legend.spacing.x= unit(0.3, 'cm'),
            legend.text = element_text(size=11),
                          plot.margin = unit(c(2, 1, 1, 0), units = "lines"),
                          #plot.title = element_textbox(size = 17, hjust = 0.5, face ="bold"),
                          axis.text.y = element_text(size=10), 
                          axis.text.x = element_blank(),
                          axis.title.y = element_text(size=12, vjust = 0.5),
                          panel.grid.major.y = element_line(linewidth=1),
                          panel.grid.major.x = element_blank(),
                          panel.grid.minor = element_blank())


pdf("barplot.pdf",width=10,height=6)
barplot
dev.off()

