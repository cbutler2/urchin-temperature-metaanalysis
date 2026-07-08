library(pacman)
p_load(tidyverse, metafor, splines)



# Preliminaries -----------------------------------------------------------

##  Basic plot formatting
theme_format <- 
  theme_bw()+
  theme(plot.title = element_text(size = 14, face = "bold"),
        plot.margin = margin(1, 1, 1, 1, "cm"),
        axis.title=element_text(size=14),
        axis.text = element_text(size = 12),
        axis.ticks = element_line(colour="black"),
        title = element_text(size = 12),
        legend.text=element_text(size=14),
        legend.title=element_text(size=14),
        legend.key.width=unit(1.2,"cm"),
        strip.text = element_text(size = 10, colour = "black"),
        strip.background = element_rect("white"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())





# Read in data ------------------------------------------------------------


d2 <- read_rds('data/UrchinReview_CompiledData_spatial.rds')

glimpse(d2)

# re categorise life stage
d2 <- d2 %>% 
  mutate(Life_Stage = factor(case_when(Development_category == 'gamete' ~ 'gamete',
                                       Development_category == 'embryo' ~ 'larvae',
                                       Development_category == 'embyro' ~ 'larvae',
                                       Development_category == 'larvae' ~ 'larvae',
                                       Development_category == 'embryo/larvae' ~ 'larvae',
                                       is.na(Development_category) ~ Urchin_Life_Stage,
                                       TRUE ~ 'oops')))

d2 <- d2 %>% 
  unnest(Effect_Size_Parameters) %>% 
  unnest(Variance_Measurements)



# Filter for species and covariate of interest ----------------------------
# the code below was used as a template to fit models to each combination or species and covariate.

dat <- d2 %>% 
  filter(Species == 'Heliocidaris erythrogramma' & Response_Type_Category == 'Development success') %>% 
  droplevels()


# Fit splines and select best-fit model ---------------------------------------------------------------


plot_list <- list()
mod_results <- list()

num_df <- c(1, 2, 3, 4)

for (i in num_df) {
  
  # fit model
  mod <- rma.mv(yi = effect_size, V = es_variance, mods = ~ ns(Temp_Anomaly, num_df[i]), random = ~1|Response_ID, data = dat)
  
  # create contrast matrix to predict over
  anomaly_preds <- seq(min(dat$Temp_Anomaly), max(dat$Temp_Anomaly), by = 1)
  
  # re-create the design matrix to use within the predict function
  # set-up the formula for the model
  form <- ~ ns(Temp_Anomaly, num_df[i])
  
  # create a new data frame that has input variables for our predictions
  newdat <- expand.grid(Temp_Anomaly = anomaly_preds)
  
  # create a new design matrix for the prediction matrix, rather than 
  # the raw data (minus the intercept)
  mod_design <- model.matrix(form, newdat)[,-1]
  
  
  # create the predictions
  preds <- 
    # depending on how many clusters there are, use the robust method if possible. If not use default methods
    predict(robust(mod, cluster = Response_ID), 
            # predict(mod, # can't use 'robust' - not enough clusters
            newmods = mod_design)
  
  mod_preds <- cbind(mod_design,
                     with(preds, data.frame(
                       pred = pred,
                       ci.lb = ci.lb,
                       ci.ub = ci.ub,
                       Temp_Anomaly = anomaly_preds)))
  # Plot
  plot <- ggplot(mod_preds) + 
    aes(x = Temp_Anomaly, y = pred) + 
    geom_line() + 
    geom_ribbon(aes(ymin = ci.lb, ymax = ci.ub), alpha =0.2, color = NA) +
    geom_point(data = dat, aes(x = Temp_Anomaly, y = effect_size)) +
    # if assessing larval development comment out the code for errorbars as these are too wide due to method of effect size calculations
    # geom_errorbar(data = dat, aes(x = Temp_Anomaly,
    #                              y = effect_size,
    #                              ymin = effect_size - es_variance,
    #                              ymax = effect_size + es_variance),
    #               width = 0.2) +
    theme_format + 
    geom_hline(yintercept = 0, linetype = 'dashed', colour = 'grey') +
    geom_vline(xintercept = 0, linetype = 'dashed', colour = 'grey') +
    scale_y_continuous(limits = c(-20, 10)) +
    labs(x = "Temperature Anomaly (\u00B0C)", 
         y = "Effect size",
         title = paste(num_df[i], ' DF', sep = ''))
  
  plot_list[[i]] <- plot
  
  plot
  
  
  print(paste(num_df[i], ' DF', sep = ''))
  
  # use appropriate print command depending on whether robust() method was used
  print(robust(mod, cluster = Response_ID))
  print(AIC(robust(mod, cluster = Response_ID)))
  # print(mod)
  # print(AIC(mod))
  
  mod_results[[i]] <- data.frame(fstat = mod$QM,
                                 df = mod$QMdf,
                                 pval = mod$QMp,
                                 aic = AIC(mod))
  
  
}

plot_list



# print out model information

mod_results <- bind_rows(lapply(mod_results, bind_rows))
mod_results

mod_results <- mod_results %>% 
  mutate(robust_method = 'yes',
         chosen_mod = 2, #change as appropriate
         species = 'heryth', #change as appropriate
         covariate = 'Metric', #change as appropriate
         level = 'Dev. success') %>% #change as appropriate
  select(species, covariate, level, fstat, df, pval, aic, robust_method, chosen_mod)


# Fit final model ---------------------------------------------------------

#change ns() DF as appropriate
mod <- rma.mv(yi = effect_size, V = es_variance, mods = ~ ns(Temp_Anomaly, 2), random = ~1|Response_ID, data = dat)

# use appropriate print command depending on whether robust() method was used 
robust(mod, cluster = Response_ID)
# mod

sink(file = "output/figures/model_outputs/heryth_devSuccess_output.txt") #change as appropriate
print('heryth_devSuccess') #change as appropriate
robust(mod, cluster = Response_ID)
sink(NULL)

# create contrast matrix to predict over
anomaly_preds <- seq(min(dat$Temp_Anomaly), max(dat$Temp_Anomaly), by = 1)

#Re-create the design matrix to use within the predict function
#Set-up the formula for the model
form <- ~ ns(Temp_Anomaly, 2)

#Create a new dataframe that has input variables for our predictions
newdat <- expand.grid(Temp_Anomaly = anomaly_preds)

#Create a new design matrix for the prediction matrix, rather than 
# the raw data (minus the intercept)
mod_design <- model.matrix(form, newdat)[,-1]


# Create the predictions
preds <- 
  predict(robust(mod, cluster = Response_ID),
          # predict(mod, 
          newmods = mod_design)

mod_preds <- cbind(mod_design,
                   with(preds, data.frame(
                     pred = pred,
                     ci.lb = ci.lb,
                     ci.ub = ci.ub,
                     Temp_Anomaly = anomaly_preds)))

mod_preds <- mod_preds %>% 
  mutate(species = 'H.eryth', #change as appropriate
         response_category = 'Development success') #change as appropriate

#saveRDS(mod_preds, file = 'data/urchin_responseCat_predictions/heryth_devSuccess.rds')


ggplot(mod_preds) + 
  aes(x = Temp_Anomaly, y = pred) + 
  geom_line() + 
  geom_ribbon(aes(ymin = ci.lb, ymax = ci.ub), alpha =0.2) +
  geom_point(data = dat, aes(x = Temp_Anomaly, y = effect_size)) +
  theme_format + 
  labs(x = "Temperature Anomaly (\u00B0C)", 
       y = "Effect size")
