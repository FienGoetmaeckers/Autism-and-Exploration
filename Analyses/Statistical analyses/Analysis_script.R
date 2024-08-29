library(ggplot2)
library(ggsignif)
library(factoextra)
library(tidyverse)
library(hrbrthemes)
library(viridis)
library(comprehenr)
library(afex)
library(brms)
library(effects)
library(sjPlot)


###
#step 1: read in all the data, make sure only relevant data is used
###

setwd("../Data")
data <- read.csv(file = "data.csv")
info <- read.csv(file = "info.csv")
quest <- read.csv(file = "quest.csv")

nr_participants <- nrow(info)
W = L = 11
nr_blocks = 10
nr_trials = 25

###
#step 2: add relevant measures to the dataframes
###

distance_cons <- c()
distances_top10 <- c()
new_clicks <- c() #boolean: true if this cell hasn't been opened before, 
prev_reward <- c()
counter <- 1
for (p in 1:nr_participants){
  for (r in 1:nr_blocks){
    
    o_cells <- c(data$initial_opened[counter]) #list with observed cells
    o_rewards <- c(2*data$average_reward[counter] - data$reward[counter]) #list with rewards of observed cells
    do <- data.frame(o_cells, o_rewards) #dataframe with all observations
    dh <- do[do$o_rewards > max(do$o_rewards)*0.9] #dataframe with only the top 10% highest reward observations
    dh["h_x"] <- dh["o_cells"]%%W
    dh["h_y"] <- (dh["o_cells"]-dh["h_x"])/W
    
    previous <- o_cells[1]
    
    for (t in 1:nr_trials){
      #save previous reward
      previous_reward <-  do[nrow(do),]$o_rewards
      
      now <- data$selected_choice[counter]
      n_x <- now%%W
      n_y <- (now - n_x)/W
      #calculate distances to all high value observations
      dh["distance"] <- abs(dh["h_x"] - n_x) + abs(dh["h_y"] - n_y)
      #find most nearby high value cell
      distances_top10 <- append(distances_top10, min(dh$distance))
      
      previous_x <- previous%%W
      previous_y <- (previous - previous_x)/W
      distance_cons <- append(distance_cons, abs(previous_x - n_x) + abs(previous_y - n_y))
      
      if (now %in% o_cells){
        new_clicks <- append(new_clicks, FALSE)
      }else{
        new_clicks <- append(new_clicks, TRUE)
      }
      
      #update
      prev_reward <- append(prev_reward, previous_reward)
      #add the current click to observations
      o_cells <- append(o_cells, now)
      new_row = c(o_cells = now, o_rewards=data$reward[counter])
      do = rbind(do,new_row)
      
      #recreate the dataframe with highet 10% observations
      dh <- do[do$o_rewards > max(do$o_rewards)*0.9,] #dataframe with only the top 10% highest reward observations
      dh["h_x"] <- dh["o_cells"]%%W
      dh["h_y"] <- (dh["o_cells"]-dh["h_x"])/W
      
      counter <- counter + 1
      
      previous <- now
    }
  }
}
data$prev_reward <- prev_reward
data$distance_prev <- distance_cons
data$distance <- distances_top10
data$new_click <- new_clicks
data$HV_click <- ifelse(data$distance == 0, TRUE, FALSE)

datalast <- data[data$trial_nr == 24,]

##add relevant behavioral participant measures to the info file
av_distance <- c()
av_distance_prev <- c()
reclicks <- c()
HVclicks <- c() #number of high-value clicks (== for this high-value distance measurement, a distance of 0)
Novclicks <- c() #number of novel clicks (hasn't been opened before)
performance <- c()
AQ <- c()
AQbin <- c()
CATI <- c()
CATI_rep <- c()
CATI_cam <- c()
SDS <- c()
SRS <- c()
PAQ <- c()
for (p in 1:nrow(info)){
  participant = info$subjectID[p]
  av_distance <- append(av_distance, 
                        mean(data[data$subjectID == participant,]$distance))
  av_distance_prev <- append(av_distance_prev, 
                             mean(data[data$subjectID == participant,]$distance_prev))
  reclicks <- append(reclicks,
                     nrow(data[data$subjectID == participant & data$distance == 0,]))
  HVclicks <- append(HVclicks, nrow(data[data$subjectID == participant & data$distance == 0,]))
  Novclicks <- append(Novclicks, nrow(data[data$subjectID == participant & data$new_click == TRUE,]))
  performance <- append(performance,
                        mean(datalast[datalast$subjectID == participant,]$average_reward))
  AQ <- append(AQ, quest[quest$subjectID == participant,]$aq_tot)
  AQbin <- append(AQbin, quest[quest$subjectID == participant,]$aq_tot_bin)
  CATI <- append(CATI, quest[quest$subjectID == participant,]$cati_tot)
  CATI_rep <- append(CATI_rep, quest[quest$subjectID == participant,]$cati_rep_beh)
  CATI_cam <- append(CATI_cam, quest[quest$subjectID == participant,]$cati_soc_cam)
  SDS <- append(SDS, quest[quest$subjectID == participant,]$sds_tot)
  SRS <- append(SRS, quest[quest$subjectID == participant,]$asrs_tot)
  PAQ <- append(PAQ, quest[quest$subjectID == participant,]$paq_tot)
}
info["av_distance"] <- av_distance
info["av_distance_prev"] <- av_distance_prev
info["reclicks"] <- reclicks
info["HVclicks"] <- HVclicks
info["Novclicks"] <- Novclicks
info["performance"] <- performance
info["AQ"] <- AQ
info["AQbin"] <- AQbin
info["CATI"] <- CATI
info["CATI_rep"] <- CATI_rep
info["CATI_cam"] <- CATI_cam
info["SDS"] <- SDS
info["SRS"] <- SRS
info["PAQ"] <- PAQ



#principal component analysis on AQ and ICAR
PCA <- prcomp(info[c("AQ", "CATI")], scale = TRUE)
summary(PCA)
fviz_eig(PCA)
fviz_pca_var(PCA,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)
info$PCA <- PCA$x[,1]
cor.test(info$PCA, info$CATI)
cor.test(info$PCA, info$CATI, method = "spearman")
cor.test(info$PCA, info$AQ)
cor.test(info$PCA, info$AQ, method = "spearman")

info["age"] = as.numeric(to_vec(for (i in 1:length(info["age"])) info["age"][i]))

#but first we need to add all AQ, age, gender data to the data dataframe
group <- c()
age <- c()
pc1 <- c()
cati_rep <- c()
cati_cam <- c()
gender <- c()
SDS <- c()
srs <- c()
IQ <- c()
PAQ <- c()
for (p in 1:length(unique(data$subjectID))){
  participant = unique(data$subjectID)[p]
  group <- append(group, to_vec(for (i in 1:(nr_trials*nr_blocks)) info[info$subjectID == participant,]$group))
  age <- append(age, to_vec(for (i in 1:(nr_trials*nr_blocks)) info[info$subjectID == participant,]$age))
  pc1 <- append(pc1, to_vec(for (i in 1:(nr_trials*nr_blocks)) info[info$subjectID == participant,]$PCA))
  cati_rep <- append(cati_rep, to_vec(for (i in 1:(nr_trials*nr_blocks)) info[info$subjectID == participant,]$CATI_rep))
  cati_cam <- append(cati_cam, to_vec(for (i in 1:(nr_trials*nr_blocks)) info[info$subjectID == participant,]$CATI_cam))
  gender <- append(gender, to_vec(for (i in 1:(nr_trials*nr_blocks)) info[info$subjectID == participant,]$gender))
  SDS <- append(SDS, to_vec(for (i in 1:(nr_trials*nr_blocks)) quest[quest$subjectID == participant,]$sds_tot))
  srs <- append(srs, to_vec(for (i in 1:(nr_trials*nr_blocks)) quest[quest$subjectID == participant,]$asrs_tot))
  IQ <- append(IQ, to_vec(for (i in 1:(nr_trials*nr_blocks)) info[info$subjectID == participant,]$ICAR_tot))
  PAQ <- append(PAQ, to_vec(for (i in 1:(nr_trials*nr_blocks)) info[info$subjectID == participant,]$PAQ))
}
data["group"] <- group
data["age"] <- age
data["PCA"] <- pc1
data["CATI_rep"] <- cati_rep
data["CATI_cam"] <- cati_cam
data["gender"] <- gender
data["SDS"] <- SDS
data["SRS"] <- srs
data["IQ"] <- IQ
data["PAQ"] <- PAQ
#for lindear trends, we use the logdistance and the logtrial
#add the log of trial_nr and distance to the dataframes
data["logtrial"] <- log(data["trial_nr"]+1.0,10) 
data["logdistance"] <- log(data["distance"]+1.0, 10)
data["logdistance_prev"] <- log(data["distance_prev"]+0.00001, 10)

options(contrasts = c("contr.sum", "contr.poly"))
info$group <- factor(info$group, levels = c("control", "autism"))
data$group <- factor(data$group, levels = c("control", "autism"))

#rescale factors for L(M)Models
info_Bfilter <- info[!(info$group == 'control' & info$AQbin > 32),] #delete all controls with too high AQbin
info_Bfilter$gender <- ifelse(info_Bfilter$gender == "Other" | info_Bfilter$gender == "Prefer not to say", "Female", info_Bfilter$gender)
info_Bfilter["gender"] <- as.factor(info_Bfilter$gender)
info_Bfilter["age"] <- scale(info_Bfilter$age)
info_Bfilter["SDS"] <- scale(info_Bfilter$SDS)
info_Bfilter["ICAR_total"] <- scale(info_Bfilter$ICAR_total)
info_Bfilter["SRS"] <- scale(info_Bfilter$SRS)
info_Bfilter["PAQ"] <- scale(info_Bfilter$PAQ)

info_Wfilter <- info[info$group == 'control',]
info_Wfilter$gender <- ifelse(info_Wfilter$gender == "Other" | info_Wfilter$gender == "Prefer not to say", "Female", info_Wfilter$gender)
info_Wfilter["PCA"] <- scale(info_Wfilter$PCA)
info_Wfilter["gender"] <- as.factor(info_Wfilter$gender)
info_Wfilter["age"] <- scale(info_Wfilter$age)
info_Wfilter["SDS"] <- scale(info_Wfilter$SDS)
info_Wfilter["ICAR_total"] <- scale(info_Wfilter$ICAR_total)
info_Wfilter["SRS"] <- scale(info_Wfilter$SRS)
info_Wfilter["PAQ"] <- scale(info_Wfilter$PAQ)

data_Bfilter <- data[data$subjectID %in% info_Bfilter$subjectID,]
data_Bfilter$gender <- ifelse(data_Bfilter$gender == "Other" | data_Bfilter$gender == "Prefer not to say", "Female", data_Bfilter$gender)
data_Bfilter["gender"] <- as.factor(data_Bfilter$gender)
data_Bfilter["age"] <- scale(data_Bfilter$age)
data_Bfilter["SDS"] <- scale(data_Bfilter$SDS)
data_Bfilter["IQ"] <- scale(data_Bfilter$IQ)
data_Bfilter["SRS"] <- scale(data_Bfilter$SRS)
data_Bfilter["PAQ"] <- scale(data_Bfilter$PAQ)
data_Bfilter["trial_nr"] <- scale(data_Bfilter$trial_nr)
data_Bfilter["logtrial"] <- scale(data_Bfilter$logtrial)
data_Bfilter["block_nr"] <- scale(data_Bfilter$block_nr)

data_Wfilter <- data[data$subjectID %in% info_Wfilter$subjectID,]
data_Wfilter$gender <- ifelse(data_Wfilter$gender == "Other" | data_Wfilter$gender == "Prefer not to say", "Female", data_Wfilter$gender)
data_Wfilter["PCA"] <- scale(data_Wfilter$PCA)
data_Wfilter["gender"] <- as.factor(data_Wfilter$gender)
data_Wfilter["age"] <- scale(data_Wfilter$age)
data_Wfilter["SDS"] <- scale(data_Wfilter$SDS)
data_Wfilter["IQ"] <- scale(data_Wfilter$IQ)
data_Wfilter["SRS"] <- scale(data_Wfilter$SRS)
data_Wfilter["PAQ"] <- scale(data_Wfilter$PAQ)
data_Wfilter["trial_nr"] <- scale(data_Wfilter$trial_nr)
data_Wfilter["logtrial"] <- scale(data_Wfilter$logtrial)
data_Wfilter["block_nr"] <- scale(data_Wfilter$block_nr)


####
#extra informative figures
###
#trend of behavioral measures over trials
#novel clicks
d <- aggregate(data$new_click, list(subjectID = data$subjectID, trial_nr = data$trial_nr), FUN=mean, digits=3)
d$x <- d$x * nr_blocks #to rescale from a percentage to a number
m <- d %>%
  group_by(trial_nr) %>% 
  summarise(mean_x = mean(x), lower_ci = mean_x - sd(x) / sqrt(nr_participants), upper_ci = mean_x + sd(x)/sqrt(nr_participants))
ggplot() +
  geom_point(data = d, aes(x = trial_nr, y = x), color = "grey", alpha = 0.08, position = position_jitter(width=1,height=.1)) +  # Individual data points
  geom_line(data = m, aes(x = trial_nr, y = mean_x), color = "darkolivegreen") +  # Mean line
  geom_ribbon(data = m, aes(x = trial_nr, ymin = lower_ci, ymax = upper_ci), fill = "darkolivegreen", alpha = 0.3) +  # Confidence interval ribbon
  theme_classic()

#high value clicks
d <- aggregate(data$HV_click, list(subjectID = data$subjectID, trial_nr = data$trial_nr), FUN=mean)
d$x <- d$x * nr_blocks #to rescale from a percentage to a number
m <- d %>%
  group_by(trial_nr) %>% 
  summarise(mean_x = mean(x), lower_ci = mean_x - sd(x) / sqrt(nr_participants), upper_ci = mean_x + sd(x)/sqrt(nr_participants))
ggplot() +
  geom_point(data = d, aes(x = trial_nr, y = x), color = "grey", alpha = 0.08, position = position_jitter(width=1,height=.1)) +  # Individual data points
  geom_line(data = m, aes(x = trial_nr, y = mean_x), color = "darkolivegreen") +  # Mean line
  geom_ribbon(data = m, aes(x = trial_nr, ymin = lower_ci, ymax = upper_ci), fill = "darkolivegreen", alpha = 0.3) +  # Confidence interval ribbon
  theme_classic()

#distance from previous click
d <- aggregate(data$distance_prev, list(data$subjectID, trial_nr = data$trial_nr), FUN=mean, digits = 4)
m <- d %>%
  group_by(trial_nr) %>% 
  summarise(mean_x = mean(x), lower_ci = mean_x - sd(x) / sqrt(nr_participants), upper_ci = mean_x + sd(x)/sqrt(nr_participants))
ggplot() +
  geom_point(data = d, aes(x = trial_nr, y = x), color = "grey", alpha = 0.08, position = position_jitter(width=1,height=.1)) +  # Individual data points
  geom_line(data = m, aes(x = trial_nr, y = mean_x), color = "darkolivegreen") +  # Mean line
  geom_ribbon(data = m, aes(x = trial_nr, ymin = lower_ci, ymax = upper_ci), fill = "darkolivegreen", alpha = 0.3) +  # Confidence interval ribbon
  theme_classic()

ggplot() +
  geom_point(data = d, aes(x = trial_nr, y = x), color = "grey", alpha = 0.2, position = position_jitter(width=1,height=.1)) +  # Individual data points
  geom_line(data = m, aes(x = trial_nr, y = mean_x), color = "darkolivegreen") +  # Mean line
  geom_ribbon(data = m, aes(x = trial_nr, ymin = lower_ci, ymax = upper_ci), fill = "darkolivegreen", alpha = 0.3) +  # Confidence interval ribbon
  theme_classic() + scale_x_log10() + scale_y_log10()

#distance from hv cell
d <- aggregate(data$distance, list(data$subjectID, trial_nr = data$trial_nr), FUN=mean, digits = 4)
m <- d %>%
  group_by(trial_nr) %>% 
  summarise(mean_x = mean(x), lower_ci = mean_x - sd(x) / sqrt(nr_participants), upper_ci = mean_x + sd(x)/sqrt(nr_participants))
ggplot() +
  geom_point(data = d, aes(x = trial_nr, y = x), color = "grey", alpha = 0.08, position = position_jitter(width=1,height=.1)) +  # Individual data points
  geom_line(data = m, aes(x = trial_nr, y = mean_x), color = "darkolivegreen") +  # Mean line
  geom_ribbon(data = m, aes(x = trial_nr, ymin = lower_ci, ymax = upper_ci), fill = "darkolivegreen", alpha = 0.3) +  # Confidence interval ribbon
  theme_classic()

ggplot() +
  geom_point(data = d, aes(x = trial_nr, y = x), color = "grey", alpha = 0.2, position = position_jitter(width=1,height=.1)) +  # Individual data points
  geom_line(data = m, aes(x = trial_nr, y = mean_x), color = "darkolivegreen") +  # Mean line
  geom_ribbon(data = m, aes(x = trial_nr, ymin = lower_ci, ymax = upper_ci), fill = "darkolivegreen", alpha = 0.3) +  # Confidence interval ribbon
  theme_classic() + scale_x_log10() + scale_y_log10()


###
###
#Section 1: Score
###
###

print(mean(info$performance))
print(sd(info$performance))
#A1: between group comparison
nrow(info_Bfilter)
print(mean(info_Bfilter[info_Bfilter$group == "control",]$performance))
print(sd(info_Bfilter[info_Bfilter$group == "control",]$performance))
print(mean(info_Bfilter[info_Bfilter$group == "autism",]$performance))
print(sd(info_Bfilter[info_Bfilter$group == "autism",]$performance))

t.test(info_Bfilter[info_Bfilter$group == 'control',]$performance, info_Bfilter[info_Bfilter$group == 'autism',]$performance, alternative = "two.sided", paired = FALSE)
ggplot(info_Bfilter, aes(x = factor(group), y = performance, fill= group)) + 
  ggdist::stat_halfeye(adjust = 0.5, justification = -.2, .width = 0, show.legend = FALSE) + 
  scale_fill_manual(values=c("grey", "orange")) + theme_classic() + 
  geom_boxplot(width = 0.15, outlier.color = NA, alpha = 0.5, show.legend = FALSE) +
  geom_signif(comparisons = list(c("control", "autism")), map_signif_level=TRUE) + 
  ggdist::stat_dots(side = "left", justification = 1.12, binwidth = 3, dotsize = 0.3, overflow = "compress", show.legend = FALSE) + ylim(35, 72)


#A2: within group comparison
nrow(info_Wfilter)
cor.test(info_Wfilter$performance,info_Wfilter$PCA , method="pearson")
cor.test(info_Wfilter$performance,info_Wfilter$PCA , method="spearman", exact = FALSE)
ggplot(info_Wfilter, aes(x = PCA, y = performance)) + geom_point() + geom_smooth(method = 'lm', color = "deepskyblue4", fill = "deepskyblue") +
  #annotate("text", x = max(info_Wfilter$PCA)-2, y = max(drop_na(info_Wfilter, performance)$performance - 1), label = paste("p = ", toString(round(c$p.value, digits = 3)))) +
  theme_classic() + ylim(35, 72)


#quick check: correlation with IQ?
cor.test(info$performance,info$ICAR_total , method="pearson")
cor.test(info$performance,info$ICAR_total , method="spearman", exact = FALSE)
ggplot(info, aes(x = ICAR_total, y = performance)) + geom_point() + geom_smooth(method = 'lm', color = "deepskyblue4", fill = "deepskyblue") +
  #annotate("text", x = max(info_Wfilter$ICAR_total)-2, y = max(drop_na(info_Wfilter, performance)$performance - 1), label = paste("p = ", toString(round(c$p.value, digits = 3))))
  theme_classic() + ylim(35, 72)

#correlation between autistic traits and ICAR
cor.test(info_Wfilter$PCA,info_Wfilter$ICAR_total , method="pearson", exact = FALSE)
cor.test(info_Wfilter$PCA,info_Wfilter$ICAR_total , method="spearman", exact = FALSE)
ggplot(info_Wfilter, aes(x = ICAR_total, y = PCA)) + geom_point() + geom_smooth(method = 'lm', color = "deepskyblue4", fill = "deepskyblue") +
  theme_classic()



t.test(info_Bfilter[info_Bfilter$group == 'control',]$ICAR_total, info_Bfilter[info_Bfilter$group == 'autism',]$ICAR_total, alternative = "two.sided", paired = FALSE)
ggplot(info_Bfilter, aes(x = group, y = ICAR_total, fill= group))+
  ggdist::stat_halfeye(adjust = 0.5, justification = -.2, .width = 0, show.legend = FALSE) + 
  scale_fill_manual(values=c("grey", "orange")) + theme_classic() + 
  geom_boxplot(width = 0.15, outlier.color = NA, alpha = 0.5, show.legend = FALSE) +
  geom_signif(comparisons = list(c("autism", "control")), map_signif_level=TRUE) + 
  ggdist::stat_dots(side = "left", justification = 1.12, dotsize = 0.4, binwidth = 0.5, overflow = "compress", show.legend = FALSE)


###
#Preregistered: trend of previous_reward vs distance varies with IQ?
###
#you can easily assess this by using a LMM to see interaction effect of previous_reward and IQ on distance_consecutive

data_temp <- data
data_temp$gender <- ifelse(data_temp$gender == "Other" | data_temp$gender == "Prefer not to say", "Female", data_temp$gender)
#scale needed factors
data_temp["age"] = scale(as.numeric(data_temp$age))
data_temp["gender"] = as.factor(data_temp$gender)
data_temp["prev_reward"] = scale(as.numeric(data_temp$prev_reward))
data_temp["IQ"] = scale(as.numeric(data_temp$IQ))

modelfull <- lmer(distance_prev ~ (prev_reward + IQ + gender + age)^2 + (1+prev_reward|subjectID), data = data_temp, control = lmerControl(optimizer = "bobyqa"))
summary(modelfull)
plot(effect("prev_reward:IQ", modelfull), ci.style="bands",
     multiline=T, key.args = list(x = 0.7, y = 0.62, corner = c(0.5, -0.4), cex = 0.5), 
     xaxt = "n",
     xaxp = c(0, 10, 20, 30, 40, 50, 60, 70, 80),
     yaxp = c(0, 2, 4, 6, 8),
     xlab="Previous reward", ylab="Distance from previous click", colors=c("grey", "darkcyan", "cyan", "chartreuse", "darkgreen"))


###
###
#Intermezzo: assessing effect of covariates
###
###
#for group
modelfull <- glm(group ~ (ICAR_total + SDS + SRS + PAQ)^2 + age + gender, data = info_Bfilter, family = binomial)
summary(modelfull)
#for traits
modelfull <- glm(PCA ~ (ICAR_total + SDS + SRS + PAQ)^2 + age + gender, data = info_Wfilter, family = gaussian)
summary(modelfull)


###
###
#Section 2: behavioral measurements for exploration
###
###
#for plotting purposes: create summary dataframe showing one trial per group/autistic traits bin
name_B <- data_Bfilter %>% group_by(group, trial_nr)
plot_dfB <- name_B %>% summarise(Novclicks = mean(new_click), HVclicks = mean(HV_click), D = mean(distance), Dprev = mean(distance_prev))

#bin the PCA data in same bins as model: -3, -1, 0 (0.1), 2, 3
bin_centers <- c(-3, -1, 0.1, 2, 3)
# find the closest bin center for each value in list
find_closest_bin <- function(value, bin_centers) {
  distances <- abs(bin_centers - value)
  closest_bin <- which.min(distances)
  return(closest_bin)
}
data_Wfilter$PCAbin <- v[sapply(data_Wfilter$PCA, find_closest_bin, bin_centers = bin_centers)]
name_W <- data_Wfilter %>% group_by(PCAbin, trial_nr)
plot_dfW <- name_W %>% summarise(Novclicks = mean(new_click), HVclicks = mean(HV_click), D = mean(distance), Dprev = mean(distance_prev))


###
#Novel clicks
###
#between group
modelfull <- glmer(new_click ~ (group + trial_nr + SDS + SRS + PAQ)^2 + gender + age + (1+trial_nr|subjectID), data = data_Bfilter, family = binomial, control = glmerControl(optimizer = "bobyqa"))
summary(modelfull)
z <- as.data.frame(effect("group:trial_nr", modelfull))
ggplot() +
  geom_line(data = z, aes(x = trial_nr, y = fit, group = group, color = group), size = 1, show.legend = FALSE) +
  geom_ribbon(data = z, aes(x = trial_nr, ymin = lower, ymax = upper, group = group, fill = group), alpha = 0.3, show.legend = FALSE) +
  labs(x = "Trial Number", y = "Fit", color = "Group", fill = "Group") +
  theme_classic() + scale_fill_manual(values=c("grey", "orange")) + scale_color_manual(values = c("grey", "orange"))


#within group
modelfull <- glmer(new_click ~ (PCA + trial_nr + SDS + SRS + PAQ)^2 + gender + age + (1+trial_nr|subjectID), data = data_Wfilter, family = binomial, control = glmerControl(optimizer = "bobyqa"))
summary(modelfull)
z <- as.data.frame(effect("PCA:trial_nr", modelfull))
ggplot() +
   geom_line(data = z, aes(x = trial_nr, y = fit, group = PCA, color = as.factor(PCA)), size = 1, show.legend = FALSE) +
   geom_ribbon(data = z, aes(x = trial_nr, ymin = lower, ymax = upper, group = PCA, fill = as.factor(PCA)), alpha = 0.18, show.legend = TRUE) +
   labs(x = "Trial Number", y = "Number of novel clicks", color = "PCA", fill = "PCA") +
   theme_classic() + scale_color_manual(values = c("deepskyblue", "grey", "yellow", "orange", "darkorange")) + scale_fill_manual(values=c("deepskyblue", "grey", "yellow", "orange", "darkorange")) 


###
#High value clicks
###
#between group
modelfull <- glmer(HV_click ~ (group + trial_nr + SDS + SRS + PAQ)^2 + gender + age + (1+trial_nr|subjectID), data = data_Bfilter, family = binomial, control = glmerControl(optimizer = "bobyqa"))
summary(modelfull)
z <- as.data.frame(effect("group:trial_nr", modelfull))
ggplot() +
  geom_line(data = z, aes(x = trial_nr, y = fit, group = group, color = group), size = 1, show.legend = FALSE) +
  geom_ribbon(data = z, aes(x = trial_nr, ymin = lower, ymax = upper, group = group, fill = group), alpha = 0.3, show.legend = FALSE) +
  labs(x = "Trial Number", y = "Number of high value clicks", color = "Group", fill = "Group") +
  theme_classic() + scale_fill_manual(values=c("grey", "orange")) + scale_color_manual(values = c("grey", "orange"))

#within group
modelfull <- glmer(HV_click ~ (PCA + trial_nr + SDS + SRS + PAQ)^2 + gender + age + (1+trial_nr|subjectID), data = data_Wfilter, family = binomial, control = glmerControl(optimizer = "bobyqa"))
summary(modelfull)
z <- as.data.frame(effect("PCA:trial_nr", modelfull))
ggplot() +
  geom_line(data = z, aes(x = trial_nr, y = fit, group = PCA, color = as.factor(PCA)), size = 1, show.legend = FALSE) +
  geom_ribbon(data = z, aes(x = trial_nr, ymin = lower, ymax = upper, group = PCA, fill = as.factor(PCA)), alpha = 0.18, show.legend = FALSE) +
  labs(x = "Trial Number", y = "Number of high value clicks", color = "PCA", fill = "PCA") +
  theme_classic() + scale_color_manual(values = c("deepskyblue", "grey", "yellow", "orange", "darkorange")) + scale_fill_manual(values=c("deepskyblue", "grey", "yellow", "orange", "darkorange")) 

###
###
#Section 3: Distance measures for exploration
###
###
###
#distance from most nearby high value cell
###

#between group
modelfull <- lmer(logdistance ~ (group + logtrial + block_nr + SDS + SRS + PAQ)^2 + gender + age + (1+logtrial*block_nr|subjectID), data = data_Bfilter, control = lmerControl(optimizer = "bobyqa"))
summary(modelfull)
z <- as.data.frame(effect("group:logtrial", modelfull))
ggplot() +
  geom_line(data = z, aes(x = logtrial, y = fit, group = group, color = group), size = 1, show.legend = FALSE) +
  geom_ribbon(data = z, aes(x = logtrial, ymin = lower, ymax = upper, group = group, fill = group), alpha = 0.3, show.legend = FALSE) +
  labs(x = "Log trial Number", y = "Distance from high value cell", color = "Group", fill = "Group") +
  theme_classic() + ylim(0.05,0.75) + scale_fill_manual(values=c("grey", "orange")) + scale_color_manual(values = c("grey", "orange"))

#within group
modelfull <- lmer(logdistance ~ (PCA + logtrial + block_nr + SDS + SRS + PAQ)^2 + gender + age + (1+logtrial*block_nr|subjectID), data = data_Wfilter, control = lmerControl(optimizer = "bobyqa"))
summary(modelfull)
z <- as.data.frame(effect("PCA:logtrial", modelfull))
ggplot() +
  geom_line(data = z, aes(x = logtrial, y = fit, group = PCA, color = as.factor(PCA)), size = 1, show.legend = FALSE) +
  geom_ribbon(data = z, aes(x = logtrial, ymin = lower, ymax = upper, group = PCA, fill = as.factor(PCA)), alpha = 0.3, show.legend = FALSE) +
  labs(x = "Log trial Number", y = "Distance from high value cell", color = "PCA", fill = "PCA") +
  theme_classic() + ylim(0.05,0.75) + scale_color_manual(values = c("deepskyblue", "grey", "yellow", "orange", "darkorange")) + scale_fill_manual(values=c("deepskyblue", "grey", "yellow", "orange", "darkorange")) 
###
#distance from previous click
###

#between group
modelfull <- lmer(logdistance_prev ~ (group + logtrial + block_nr + SDS + SRS + PAQ)^2 + gender + age + (1+logtrial*block_nr|subjectID), data = data_Bfilter, control = lmerControl(optimizer = "bobyqa"))
summary(modelfull)
z <- as.data.frame(effect("group:logtrial", modelfull))
ggplot() +
  geom_line(data = z, aes(x = logtrial, y = fit, group = group, color = group), size = 1, show.legend = FALSE) +
  geom_ribbon(data = z, aes(x = logtrial, ymin = lower, ymax = upper, group = group, fill = group), alpha = 0.3, show.legend = FALSE) +
  labs(x = "Log trial Number", y = "Distance from high value cell", color = "Group", fill = "Group") +
  theme_classic() + ylim(-3,2) + scale_fill_manual(values=c("grey", "orange")) + scale_color_manual(values = c("grey", "orange"))


#within group
modelfull <- lmer(logdistance_prev ~ (PCA + logtrial + block_nr + SDS + SRS + PAQ)^2 + gender + age + (1+logtrial*block_nr|subjectID), data = data_Wfilter, control = lmerControl(optimizer = "bobyqa"))
summary(modelfull)
z <- as.data.frame(effect("PCA:logtrial", modelfull))
ggplot() +
  geom_line(data = z, aes(x = logtrial, y = fit, group = PCA, color = as.factor(PCA)), size = 1, show.legend = FALSE) +
  geom_ribbon(data = z, aes(x = logtrial, ymin = lower, ymax = upper, group = PCA, fill = as.factor(PCA)), alpha = 0.3, show.legend = FALSE) +
  labs(x = "Log trial Number", y = "Distance from high value cell", color = "PCA", fill = "PCA") +
  theme_classic() + ylim(-3,2) + scale_color_manual(values = c("deepskyblue", "grey", "yellow", "orange", "darkorange")) + scale_fill_manual(values=c("deepskyblue", "grey", "yellow", "orange", "darkorange")) 


###
###
#Section 4: modelling results
###
###


#est <- read.csv(file = "1env_0405.csv") #estimates over all trials #classic model
est <- read.csv(file = "estimates.csv")                             #localized version
#only keep data for which we have an estimate
est <- est[est$Participant %in% data$subjectID,]
infoS <- info[info$subjectID %in% est$Participant,]

l <- c()
b <- c()
t <- c()
NLL <- c()

for (i in 1:length(unique(infoS$subjectID))){
  participant <- unique(infoS$subjectID)[i]
  l <- append(l, est[est$Participant == participant,]$l_fit)
  b <- append(b, est[est$Participant == participant,]$beta)
  t <- append(t, est[est$Participant == participant,]$tau)
  NLL <- append(NLL, est[est$Participant == participant,]$NLL)
}
infoS["l"] <- l
infoS["b"] <- b
infoS["t"] <- t
infoS["NLL"] <- NLL

info_B <- infoS[!(infoS$group == 'control' & infoS$AQbin > 32),] #delete all controls with too high AQbin
info_B$group <- factor(info_B$group, levels = c("control", "autism"))

#rescale factors for L(M)Models
info_B$gender <- ifelse(info_B$gender == "Other" | info_B$gender == "Prefer not to say", "Female", info_B$gender)
info_B["gender"] <- as.factor(info_B$gender)
info_B["age"] <- scale(info_B$age)
info_B["SDS"] <- scale(info_B$SDS)
info_B["ICAR_total"] <- scale(info_B$ICAR_total)
info_B["SRS"] <- scale(info_B$SRS)
info_B["PAQ"] <- scale(info_B$PAQ)
info_B["logl"] <- scale(log(info_B$l, 10))
info_B["logb"] <- scale(log(info_B$b, 10))
info_B["logt"] <- scale(log(info_B$t, 10))
info_B["NLL"] <- scale(info_B$NLL)

info_W <- infoS[infoS$subjectID %in% info_Wfilter$subjectID,]
info_W$gender <- ifelse(info_W$gender == "Other" | info_W$gender == "Prefer not to say", "Female", info_W$gender)
info_W["PCA"] <- scale(info_W$PCA)
info_W["gender"] <- as.factor(info_W$gender)
info_W["age"] <- scale(info_W$age)
info_W["SDS"] <- scale(info_W$SDS)
info_W["ICAR_total"] <- scale(info_W$ICAR_total)
info_W["SRS"] <- scale(info_W$SRS)
info_W["PAQ"] <- scale(info_W$PAQ)
info_W["logl"] <- scale(log(info_W$l, 10))
info_W["logb"] <- scale(log(info_W$b, 10))
info_W["logt"] <- scale(log(info_W$t, 10))
info_W["NLL"] <- scale(info_W$NLL)


#we type them here in the order that we discuss them:
#1: all in interaction, all covariates
m <- glm(formula = group ~ (gender + age + SDS + SRS + PAQ)^2 + logl + logb + logt, family = binomial, data = info_B)
summary(m)
plot(effect("logb", m),  ci.style="bands")


#only generalization
m <- glm(formula = group ~ (gender + age + ICAR_total + SDS + SRS + PAQ)^2 + logl, family = binomial, data = info_B)
summary(m)
m <- glm(formula = group ~ (gender + age)^2 + logl, family = binomial, data = info_B)
summary(m)
m <- glm(formula = group ~ logl, family = binomial, data = info_B)
summary(m)

#only Uncertainty guided exploration
m <- glm(formula = group ~ (gender + age + ICAR_total + SDS + SRS + PAQ)^2 + logb, family = binomial, data = info_B)
summary(m)
plot(effect("logb", m),  ci.style="bands")
m <- glm(formula = group ~ (gender + age)^2 + logb, family = binomial, data = info_B)
summary(m)
plot(effect("logb", m),  ci.style="bands")
m <- glm(formula = group ~ logb, family = binomial, data = info_B)
summary(m)
plot(effect("logb", m),  ci.style="bands")

#only random exploration
m <- glm(formula = group ~ (gender + age + ICAR_total + SDS + SRS + PAQ)^2 + logt, family = binomial, data = info_B)
summary(m)
plot(effect("logt", m),  ci.style="bands")
m <- glm(formula = group ~ (gender + age)^2 + logt, family = binomial, data = info_B)
summary(m)
plot(effect("logt", m),  ci.style="bands")
m <- glm(formula = group ~ logt, family = binomial, data = info_B)
summary(m)
plot(effect("logt", m),  ci.style="bands")



###
#for within group:
###
m <- glm(formula = PCA ~ (gender + age + SDS + SRS + PAQ)^2 + logl + logb + logt, family = gaussian, data = info_W)
summary(m)
plot(effect("logb", m),  ci.style="bands")

#generalization
m <- glm(formula = PCA ~ (gender + age + ICAR_total + SDS + SRS + PAQ)^2 + logl, family = gaussian, data = info_W)
summary(m)
m <- glm(formula = PCA ~ (gender + age)^2 + logl, family = gaussian, data = info_W)
summary(m)
m <- glm(formula = PCA ~ logl, family = gaussian, data = info_W)
summary(m)

#uncertainty guided exploration
m <- glm(formula = PCA ~ (gender + age + ICAR_total + SDS + SRS + PAQ)^2 + logb, family = gaussian, data = info_W)
summary(m)
m <- glm(formula = PCA ~ (gender + age )^2 + logb, family = gaussian, data = info_W)
summary(m)
m <- glm(formula = PCA ~  logb, family = gaussian, data = info_W)
summary(m)
plot(effect("logb", m), ci.style="bands")

#random exploration
m <- glm(formula = PCA ~ (gender + age + ICAR_total + SDS + SRS + PAQ)^2 + logt, family = gaussian, data = info_W)
summary(m)
m <- glm(formula = PCA ~ (gender + age)^2 + logt, family = gaussian, data = info_W)
summary(m)
m <- glm(formula = PCA ~ logt, family = gaussian, data = info_W)
summary(m)

###
###
#Section 5: matched subsets
###
###

#first, bin all the PCA (or round to 10, to allow for easier matching!)
scale <- max(info$PCA) - min(info$PCA)
range <- scale * 4/100
matched <- to_vec(for (i in 1:nrow(info)) 0)
for (p in 1:nrow(info)){
  print("next")
  if (info$group[p] == "autism"){ #only then continue
    participant <- info$subjectID[p]
    autism_traits <- info$PCA[p]
    print(autism_traits)
    #find a match by looping over all other pp, starting from where we are in the first loop
    p2 <- 1
    no_match <- TRUE
    while (no_match & p2 < nrow(info)){
      if (info$group[p2] == "control" & matched[p2] == 0){ #only check control participants that haven't been matched
        if (autism_traits - range < info$PCA[p2] & info$PCA[p2] < autism_traits + range){
          print(info$PCA[p2])
          matched[p] <- matched[p] + 1
          matched[p2] <- matched[p2] + 1
          no_match <- FALSE
        }
      }
      p2 <- p2 + 1
    }
  }
}
info["matched"] <- matched 
infoM <- info[info$matched == 1,]
#do they have the same PCA values? YES
t.test(infoM[infoM$group == 'autism',]$PCA, infoM[infoM$group == 'control',]$PCA, alternative = "two.sided", paired = FALSE)
step <- (max(infoM$PCA) - min(infoM$PCA))/25
ggplot(infoM, aes(x = factor(group), y = PCA, fill= group)) + 
  ggdist::stat_halfeye(adjust = 0.5, justification = -.2, .width = 0, show.legend = FALSE) + 
  scale_fill_manual(values=c("grey", "orange")) + theme_classic() + 
  geom_boxplot(width = 0.15, outlier.color = NA, alpha = 0.5, show.legend = FALSE) +
  geom_signif(comparisons = list(c("autism", "control")), map_signif_level=TRUE) + 
  ggdist::stat_dots(side = "left", justification = 1.12, binwidth = step, overflow = "compress", show.legend = FALSE)


#do they differ in behaviour?
t.test(infoM[infoM$group == 'autism',]$Novclicks, infoM[infoM$group == 'control',]$Novclicks, alternative = "two.sided", paired = FALSE)
step <- (max(infoM$Novclicks) - min(infoM$Novclicks))/25
ggplot(infoM, aes(x = factor(group), y = Novclicks, fill= group)) + 
  ggdist::stat_halfeye(adjust = 0.5, justification = -.2, .width = 0, show.legend = FALSE) + 
  scale_fill_manual(values=c("grey", "orange")) + theme_classic() + 
  geom_boxplot(width = 0.15, outlier.color = NA, alpha = 0.5, show.legend = FALSE) +
  geom_signif(comparisons = list(c("autism", "control")), map_signif_level=TRUE) + 
  ggdist::stat_dots(side = "left", justification = 1.12, binwidth = step, overflow = "compress", show.legend = FALSE)

t.test(infoM[infoM$group == 'autism',]$HVclicks, infoM[infoM$group == 'control',]$HVclicks, alternative = "two.sided", paired = FALSE)
step <- (max(infoM$HVclicks) - min(infoM$HVclicks))/25
ggplot(infoM, aes(x = factor(group), y = HVclicks, fill= group)) + 
  ggdist::stat_halfeye(adjust = 0.5, justification = -.2, .width = 0, show.legend = FALSE) + 
  scale_fill_manual(values=c("grey", "orange")) + theme_classic() + 
  geom_boxplot(width = 0.15, outlier.color = NA, alpha = 0.5, show.legend = FALSE) +
  geom_signif(comparisons = list(c("autism", "control")), map_signif_level=TRUE) + 
  ggdist::stat_dots(side = "left", justification = 1.12, binwidth = step, overflow = "compress", show.legend = FALSE)


#do they have different model parameters?
infoM <- info[info$subjectID %in% est$Participant,]
infoM["l"] <- l
infoM["b"] <- b
infoM["t"] <- t
infoM["logl"] = log(infoM$l, 10)
infoM["logb"] = log(infoM$b, 10)
infoM["logt"] = log(infoM$t, 10)
infoM <- infoM[infoM$matched == 1,]
#l
t.test(infoM[infoM$group == 'autism',]$logl, infoM[infoM$group == 'control',]$logl, alternative = "two.sided", paired = FALSE)
step <- (max(infoM$logl) - min(infoM$logl))/25
ggplot(infoM, aes(x = factor(group), y = logl, fill= group)) + 
  ggdist::stat_halfeye(adjust = 0.5, justification = -.2, .width = 0, show.legend = FALSE) + 
  scale_fill_manual(values=c("grey", "orange")) + theme_classic() + 
  geom_boxplot(width = 0.15, outlier.color = NA, alpha = 0.5, show.legend = FALSE) +
  geom_signif(comparisons = list(c("autism", "control")), map_signif_level=TRUE) + 
  ggdist::stat_dots(side = "left", justification = 1.12, binwidth = step, overflow = "compress", show.legend = FALSE)

#b
t.test(infoM[infoM$group == 'autism',]$logb, infoM[infoM$group == 'control',]$logb, alternative = "two.sided", paired = FALSE)
step <- (max(infoM$logb) - min(infoM$logb))/25
ggplot(infoM, aes(x = factor(group), y = logb, fill= group)) + 
  ggdist::stat_halfeye(adjust = 0.5, justification = -.2, .width = 0, show.legend = FALSE) + 
  scale_fill_manual(values=c("grey", "orange")) + theme_classic() + 
  geom_boxplot(width = 0.15, outlier.color = NA, alpha = 0.5, show.legend = FALSE) +
  geom_signif(comparisons = list(c("autism", "control")), map_signif_level=TRUE) + 
  ggdist::stat_dots(side = "left", justification = 1.12, binwidth = step, overflow = "compress", show.legend = FALSE)

#t
t.test(infoM[infoM$group == 'autism',]$logt, infoM[infoM$group == 'control',]$logt, alternative = "two.sided", paired = FALSE)
step <- (max(infoM$logt) - min(infoM$logt))/25
ggplot(infoM, aes(x = factor(group), y = logt, fill= group)) + 
  ggdist::stat_halfeye(adjust = 0.5, justification = -.2, .width = 0, show.legend = FALSE) + 
  scale_fill_manual(values=c("grey", "orange")) + theme_classic() + 
  geom_boxplot(width = 0.15, outlier.color = NA, alpha = 0.5, show.legend = FALSE) +
  geom_signif(comparisons = list(c("autism", "control")), map_signif_level=TRUE) + 
  ggdist::stat_dots(side = "left", justification = 1.12, binwidth = step, overflow = "compress", show.legend = FALSE) 

