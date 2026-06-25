library(tidyverse)
library(ggplot2)
library(cowplot)
library(GA)
site_sp_vec <- c("Corrigin_Eucalyptus_capillosa",
                 "ManyPeaksRange_Alphitonia_excelsa",
                 "ManyPeaksRange_Austromyrtus_bidwillii",
                 "ManyPeaksRange_Brachychiton_australis",
                 "ManyPeaksRange_Cochlospermum_gillivraei",
                 "Parque_Natural_Metropolitano_Calycophyllum_candidissimum",
                 "Puechabon_Quercus_ilex",
                 "Richmond_Eucalyptus_cladocalyx",
                 "Richmond_Eucalyptus_dunnii",
                 "Richmond_Eucalyptus_saligna",
                 "San_Lorenzo_Carapa_guianensis",
                 "San_Lorenzo_Tachigali_versicolor",
                 "San_Lorenzo_Tocoyena_pittieri",
                 "Sevilleta_Juniperus_monosperma",
                 "Sevilleta_Pinus_edulis",
                 "Vic_la_Gardiole_Quercus_ilex")
sp_vec <- c("Eucalyptus capillosa",
            "Alphitonia excelsa",
            "Austromyrtus bidwillii",
            "Brachychiton australis",
            "Cochlospermum gillivraei",
            "Calycophyllum candidissimum",
            "Quercus ilex",
            "Eucalyptus cladocalyx",
            "Eucalyptus dunnii",
            "Eucalyptus saligna",
            "Carapa guianensis",
            "Tachigali versicolor",
            "Tocoyena pittieri",
            "Juniperus monosperma",
            "Pinus edulis",
            "Quercus ilex")
site_vec <- c("Corrigin",
              "ManyPeaksRange",
              "ManyPeaksRange",
              "ManyPeaksRange",
              "ManyPeaksRange",
              "Parque_Natural",
              "Puechabon",
              "Richmond",
              "Richmond",
              "Richmond",
              "San_Lorenzo",
              "San_Lorenzo",
              "San_Lorenzo",
              "Sevilleta",
              "Sevilleta",
              "Vic_la_Gardiole")
valid = c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16)


# Relationship Gs_P50 vs Gs_slope -----------------------------------------
# From Martin-StPaul et al. (2017)
pgs90_Rev <- readxl::read_excel("data-raw/pgs90_Rev.xlsx") |>
  dplyr::select(Species, pgs90, pgs50, slope) |>
  dplyr::mutate(slope = 25/slope) |>
  dplyr::filter(slope < 300)
model_slope<-nls(slope ~ x1/pgs90^2 + x2/pgs90, start = list(x1= 1, x2 = 2), data = pgs90_Rev)
plot(pgs90_Rev$pgs90, pgs90_Rev$slope)
P90 <- seq(-6, -0.5, by = 0.1)
slope <- coefficients(model_slope)[1]/P90^2 + coefficients(model_slope)[2]/P90
lines(x = P90, y=slope)
plot(pgs90_Rev$pgs50, pgs90_Rev$pgs90)
model_pgs90<-lm(pgs90 ~ pgs50, data = pgs90_Rev)
coefficients(model_pgs90)
WFO_file <- "~/OneDrive/EMF_datasets/Taxonomy/WFO/WFO_Backbone/classification.csv"
DB_path <- "~/OneDrive/mcaceres_work/model_development/medfate_parameterization/traits_and_models/"
harmonized_trait_path <- paste0(DB_path,"data/harmonized_trait_sources")
db <- data.frame(originalName = sp_vec)
db_post <- traits4models::harmonize_taxonomy_WFO(db, WFO_backbone_file = WFO_file)
trait_obs  <- traits4models::get_taxon_trait_means(harmonized_trait_path, traits = c("Nleaf", "Vmax", "Ptlp", "Gs_P50", 
                                                                                     "VCstem_P50", "VCstem_slope", "LeafPI0", "LeafEPS", "Gswmin"))
## USE RELATIONSHIP BETWEEN P50 and slope
is_na <- is.na(trait_obs$VCstem_slope)
trait_obs$VCstem_slope[is_na] = 478 / (trait_obs$VCstem_P50[is_na]^2) - 149/trait_obs$VCstem_P50[is_na]

# ggplot(trait_obs)+
#   geom_point(aes(x=VCstem_P50, y=VCstem_slope))
# summary(trait_obs$VCstem_slope)
# P50 <- seq(-1, -8, by = -0.1)
# slope <- 478/(P50^2) - 149/P50
# plot(P50, slope, type="lines")

df_res <- data.frame(Site = site_vec,
                     Species_original = sp_vec,
                     Species_accepted = db_post$acceptedName) |>
  dplyr::left_join(trait_obs, by=c("Species_accepted"="acceptedName")) |>
  dplyr::mutate(n = NA,
                n_filt = NA,
                Pleaf_min = NA,
                Pleaf_max = NA,
                gs_max = NA,
                An_Cs_gs_max = NA,
                M1_g0 = NA,
                M2_g0 = NA,
                M3_g0 = NA,
                M4_g0 = NA,
                M1_g1 = NA,
                M2_g1 = NA,
                M3_g1 = NA,
                M4_g1 = NA,
                M1_gs_P50 = NA,
                M2_gs_P50 = NA,
                M3_gs_P50 = NA,
                M4_gs_P50 = NA,
                M1_gs_slope = NA,
                M2_gs_slope = NA,
                M3_gs_slope = NA,
                M4_gs_slope = NA,
                M1_mae = NA,
                M2_mae = NA,
                M3_mae = NA,
                M4_mae = NA,
                M1_rmse = NA,
                M2_rmse = NA,
                M3_rmse = NA,
                M4_rmse = NA,
                M1_r2 = NA,
                M2_r2 = NA,
                M3_r2 = NA,
                M4_r2 = NA)

## Sigmoid function
sigmoid <- function(y, P50, slope) {
  1 - 1/(1+exp((slope/25)*(y - P50)))
}
## Stomatal conductance function
gs <- function(x, y, g0, g1, gs_P50, gs_slope) {
  fl <- 1 - 1/(1+exp((gs_slope/25)*(y - gs_P50)))
  fl*pmax(g0, g0 + g1*x)
}
## Slope P50 relationship
f_slope <- function(P50) {
  P90 = -0.7497035 +  1.1525639*P50
  -40.20194/P90^2 + -208.65225/P90
}
## Optimization function (all parameters) - DISCARDED
f_opt_g0_g1_P50_slope <- function(x, xydata) {
  diff <- xydata$gs - gs(xydata$An_Cs, xydata$Pleaf, x[1], x[2], x[3], x[4])
  return(- mean(abs(diff)))
}
## Optimization function (all parameters except g0)
f_opt_g1_P50_slope <- function(x, xydata, g0) {
  diff <- xydata$gs - gs(xydata$An_Cs, xydata$Pleaf, g0, x[1], x[2], x[3])
  return(- mean(abs(diff)))
}
## Optimization function (all parameters except g0 and slope)
f_opt_g1_P50 <- function(x, xydata, g0) {
  diff <- xydata$gs - gs(xydata$An_Cs, xydata$Pleaf, g0, x[1], x[2], f_slope(x[2]))
  return(- mean(abs(diff)))
}
## Optimization function (only g1)
f_opt_g1 <- function(x, xydata, g0, gs_P50, gs_slope) {
  diff <- xydata$gs - gs(xydata$An_Cs, xydata$Pleaf, g0, x[1], gs_P50, gs_slope)
  return(- mean(abs(diff)))
}

 
for(index in valid) {
  print(site_sp_vec[index])
  xdata <- readr::read_csv(paste0("data-raw/ManonSabot-One_gs_model_to_rule_them_all-41cbd30/input/calibrations/obs_driven/", site_sp_vec[index], "_x.csv"), show_col_types = FALSE)
  ydata <- readr::read_csv(paste0("data-raw/ManonSabot-One_gs_model_to_rule_them_all-41cbd30/input/calibrations/obs_driven/", site_sp_vec[index], "_y.csv"), show_col_types = FALSE)

  xdata_day <- xdata[-1,]|>
    select(year, doy, hod, Ps, PPFD, Tair, Patm, CO2, u) |>
    mutate(year = as.integer(year),
           doy = as.integer(doy),
           hod = as.numeric(hod),
           Ps = as.numeric(Ps),
           PPFD = as.numeric(PPFD),
           Tair = as.numeric(Tair),
           Patm = as.numeric(Patm),
           CO2 = as.numeric(CO2),
           u = as.numeric(u),
           CO2_conc = CO2*1000/Patm)|>
    group_by(year, doy) |>
    summarise(Patm = mean(Patm, na.rm=TRUE), 
              CO2_conc = mean(CO2_conc, na.rm=TRUE), 
              .groups = "drop")
  
  ydata <- ydata[-1,] |>
    mutate(year = as.integer(year),
           doy = as.integer(doy),
           hod = as.numeric(hod),
           A = as.numeric(A),
           E = as.numeric(E),
           gs = as.numeric(gs),
           gb = as.numeric(gb),
           Tleaf = as.numeric(Tleaf),
           Pleaf = as.numeric(Pleaf))
  
  
  xydata <- ydata |>
    left_join(xdata_day, by=c("year", "doy"))|>
    filter(A >0) |>
    mutate(Cs = CO2_conc - A/gb,
           An_Cs = A/Cs) |>
    filter(!is.na(Pleaf) & !is.na(An_Cs))

  
  Q1 <- quantile(xydata$gs, 0.25)
  Q3 <- quantile(xydata$gs, 0.75)
  IQR_val <- IQR(xydata$gs)
  
  # Define bounds
  # lower_bound <- Q1 - 1.5 * IQR_val
  upper_bound <- Q3 + 1.5 * IQR_val
  
  df_res$n[index] <- nrow(xydata)
  xydata <- xydata |>
    filter(gs<=upper_bound )
    
  df_res$n_filt[index] <- nrow(xydata)
  df_res$gs_max[index] <- quantile(xydata$gs, 0.99)
  df_res$An_Cs_gs_max[index] <- mean(xydata$An_Cs[xydata$gs >= df_res$gs_max[index]])
  df_res$Pleaf_min[index] <- min(xydata$Pleaf)
  df_res$Pleaf_max[index] <- max(xydata$Pleaf)

  ## MODEL 1
  ga_M1 <- ga(type = "real-valued",
               fitness = f_opt_g1_P50_slope, xydata = xydata, g0 = df_res$gs_max[index]*0.05,
               popSize = 300,
               maxiter = 400,
               lower= c(1, df_res$Pleaf_min[index], 10),
               upper = c(20, df_res$Pleaf_max[index], 200),
               optim = TRUE, monitor = FALSE)
  
  df_res$M1_mae[index] <- -ga_M1@fitnessValue
  df_res$M1_g0[index] <- df_res$gs_max[index]*0.05 
  df_res$M1_g1[index] <-ga_M1@solution[1,1]
  df_res$M1_gs_P50[index] <-ga_M1@solution[1,2]
  df_res$M1_gs_slope[index] <- ga_M1@solution[1,3]
  xydata$M1_gs_pred <- gs(xydata$An_Cs,xydata$Pleaf, 
                          df_res$M1_g0[index], 
                          df_res$M1_g1[index], 
                          df_res$M1_gs_P50[index], 
                          df_res$M1_gs_slope[index])
  df_res$M1_r2[index] <- cor(xydata$gs, xydata$M1_gs_pred, use = "complete.obs")^2
  df_res$M1_rmse[index] <- sqrt(mean((xydata$gs-xydata$M1_gs_pred)^2))
  
  
  ## MODEL 2
  ga_M2 <- ga(type = "real-valued",
              fitness = f_opt_g1_P50, xydata = xydata, g0 = df_res$gs_max[index]*0.05,
              popSize = 300,
              maxiter = 400,
              lower= c(1, df_res$Pleaf_min[index]),
              upper = c(20, df_res$Pleaf_max[index]),
              optim = TRUE, monitor = FALSE)
  
  df_res$M2_mae[index] <- -ga_M2@fitnessValue
  df_res$M2_g0[index] <- df_res$gs_max[index]*0.05 
  df_res$M2_g1[index] <-ga_M2@solution[1,1]
  df_res$M2_gs_P50[index] <-ga_M2@solution[1,2]
  df_res$M2_gs_slope[index] <- f_slope(ga_M2@solution[1,2])
  xydata$M2_gs_pred <- gs(xydata$An_Cs,xydata$Pleaf, 
                          df_res$M2_g0[index], 
                          df_res$M2_g1[index], 
                          df_res$M2_gs_P50[index], 
                          df_res$M2_gs_slope[index])
  df_res$M2_r2[index] <- cor(xydata$gs, xydata$M2_gs_pred, use = "complete.obs")^2
  df_res$M2_rmse[index] <- sqrt(mean((xydata$gs-xydata$M2_gs_pred)^2))
  
  
  tlp <- df_res$Ptlp[index]
  P50 <- df_res$VCstem_P50[index]
  is_m3_m4 <- FALSE
  if(!is.na(tlp) || !is.na(P50)) {
    is_m3_m4 <- TRUE
    if(is.na(tlp)) {
      P_gs88 <- max(-3.5, P50*0.23 - 1.5)
      P_gs12 <- P_gs88 + 1
    } else {
      P_gs88 <- tlp
      P_gs12 <- tlp/2
    }
    slope <- (88.0 - 12.0)/(abs(P_gs88) - abs(P_gs12))
    P_gs50 <- (P_gs88+P_gs12)/2.0
  }
  if(is_m3_m4) {
    ## MODEL 3
    ga_M3 <- ga(type = "real-valued",
                fitness = f_opt_g1, xydata = xydata, g0 = df_res$gs_max[index]*0.05, gs_P50 = P_gs50, gs_slope = slope,
                popSize = 300,
                maxiter = 400,
                lower= c(1, df_res$Pleaf_min[index]),
                upper = c(20, df_res$Pleaf_max[index]),
                optim = TRUE, monitor = FALSE)
    df_res$M3_g0[index] <- df_res$gs_max[index]*0.05 
    df_res$M3_g1[index] <-ga_M3@solution[1,1]
    df_res$M3_gs_P50[index] <- P_gs50
    df_res$M3_gs_slope[index] <- slope
    xydata$M3_gs_pred <- gs(xydata$An_Cs,xydata$Pleaf, 
                            df_res$M3_g0[index], 
                            df_res$M3_g1[index], 
                            df_res$M3_gs_P50[index], 
                            df_res$M3_gs_slope[index])
    df_res$M3_mae[index] <- mean(abs(xydata$gs-xydata$M3_gs_pred))
    df_res$M3_r2[index] <- cor(xydata$gs, xydata$M3_gs_pred, use = "complete.obs")^2
    df_res$M3_rmse[index] <- sqrt(mean((xydata$gs-xydata$M3_gs_pred)^2))
    
    ## MODEL 4
    df_res$M4_g0[index] <- df_res$gs_max[index]*0.05 
    df_res$M4_g1[index] <- (df_res$gs_max[index] - df_res$M4_g0[index])/df_res$An_Cs_gs_max[index]
    df_res$M4_gs_P50[index] <- P_gs50
    df_res$M4_gs_slope[index] <- slope
    xydata$M4_gs_pred <- gs(xydata$An_Cs,xydata$Pleaf, 
                            df_res$M4_g0[index], 
                            df_res$M4_g1[index], 
                            df_res$M4_gs_P50[index], 
                            df_res$M4_gs_slope[index])
    df_res$M4_mae[index] <- mean(abs(xydata$gs-xydata$M4_gs_pred))
    df_res$M4_r2[index] <- cor(xydata$gs, xydata$M4_gs_pred, use = "complete.obs")^2
    df_res$M4_rmse[index] <- sqrt(mean((xydata$gs-xydata$M4_gs_pred)^2))
  }
  print(df_res[index,])
  
  g1 <- ggplot(xydata)+
    geom_point(aes(x = gs, y = M1_gs_pred), col = "black", alpha= 0.5)+
    geom_point(aes(x = gs, y = M2_gs_pred), col = "red", alpha= 0.5)+
    geom_text(x =max(xydata$gs)*0.05, y = max(xydata$gs)*0.95,
              label = paste0("M1 - R2 = ", round(100*df_res$M1_r2[index],1), 
                             "% MAE = ", round(df_res$M1_mae[index],3), 
                             " RMSE = ", round(df_res$M1_rmse[index],3)),
              size = 3, col = "black",
              hjust = "left")+
    geom_text(x =max(xydata$gs)*0.05, y = max(xydata$gs)*0.90,
              label = paste0("M2 - R2 = ", round(100*df_res$M2_r2[index],1), 
                             "% MAE = ", round(df_res$M2_mae[index],3), 
                             " RMSE = ", round(df_res$M2_rmse[index],3)),
              size = 3, col = "red",
              hjust = "left")
    
  if(is_m3_m4) {
    g1 <- g1+
      geom_point(aes(x = gs, y = M3_gs_pred), col = "blue", alpha= 0.5)+
      geom_point(aes(x = gs, y = M4_gs_pred), col = "darkgreen", alpha= 0.5)+
      geom_abline(intercept = 0, slope = 1, col = "gray")+
      geom_text(x =max(xydata$gs)*0.05, y = max(xydata$gs)*0.85,
                label = paste0("M3 - R2 = ", round(100*df_res$M3_r2[index],1), 
                               "% MAE = ", round(df_res$M3_mae[index],3), 
                               " RMSE = ", round(df_res$M3_rmse[index],3)),
                size = 3, col = "blue",
                hjust = "left")+
      geom_text(x =max(xydata$gs)*0.05, y = max(xydata$gs)*0.80,
                label = paste0("M4 - R2 = ", round(100*df_res$M4_r2[index],1), 
                               "% MAE = ", round(df_res$M4_mae[index],3), 
                               " RMSE = ", round(df_res$M4_rmse[index],3)),
                size = 3, col = "darkgreen",
                hjust = "left")
  }
  g1 <- g1+
    xlab("Observed stomatal conductance")+
    ylab("Predicted stomatal conductance")+
    scale_x_continuous(limits = c(0,max(xydata$gs)))+
    scale_y_continuous(limits = c(0,max(xydata$gs)))+
    theme_bw() +
    labs(subtitle = " ")
  
  mdata1 <- data.frame(Pleaf = seq(-9, 0, by=0.01),
                       An_Cs = max(xydata$An_Cs)) |>
    dplyr::mutate(M1_gs_pred = gs(An_Cs, Pleaf, 
                               df_res$M1_g0[index], df_res$M1_g1[index], 
                               df_res$M1_gs_P50[index], df_res$M1_gs_slope[index])) |>
    dplyr::mutate(M2_gs_pred = gs(An_Cs, Pleaf, 
                               df_res$M2_g0[index], df_res$M2_g1[index], 
                               df_res$M2_gs_P50[index], df_res$M2_gs_slope[index]))
  if(is_m3_m4) {
    mdata1 <- mdata1 |>
    dplyr::mutate(M3_gs_pred = gs(An_Cs, Pleaf, 
                                 df_res$M3_g0[index], df_res$M3_g1[index], 
                                 df_res$M3_gs_P50[index], df_res$M3_gs_slope[index]))|>
    dplyr::mutate(M4_gs_pred = gs(An_Cs, Pleaf, 
                                  df_res$M4_g0[index], df_res$M4_g1[index], 
                                  df_res$M4_gs_P50[index], df_res$M4_gs_slope[index]))
  }
  
  g2 <- ggplot(xydata)+
    geom_point(aes(x = Pleaf, y = gs, col = An_Cs))+
    geom_line(aes(x=Pleaf, y=M1_gs_pred), data = mdata1, col = "black")+
    geom_vline(xintercept = df_res$M1_gs_P50[index], col = "black", size = 2, alpha = 0.5)+
    geom_line(aes(x=Pleaf, y=M2_gs_pred), data = mdata1, col = "red")+
    geom_vline(xintercept = df_res$M2_gs_P50[index], col = "red", size = 2, alpha = 0.5)
  if(is_m3_m4) {
    g2 <- g2+
      geom_line(aes(x=Pleaf, y=M3_gs_pred), data = mdata1, col = "blue")+
      geom_vline(xintercept = df_res$M3_gs_P50[index], col = "blue", size = 2, alpha = 0.5)+
      geom_line(aes(x=Pleaf, y=M4_gs_pred), data = mdata1, col = "darkgreen")+
      geom_vline(xintercept = df_res$M4_gs_P50[index], col = "darkgreen", size = 2, alpha = 0.5)
  }
  g2 <- g2+
    scale_x_continuous(limits = c(df_res$Pleaf_min[index]-0.5,0))+
    xlab("Leaf water potential (MPa)")+
    ylab("Observed stomatal conductance")
  if(!is.na(df_res$VCstem_P50[index]) && !is.na(df_res$VCstem_slope[index])) {
    mdata1$PLC <- sigmoid(mdata1$Pleaf, df_res$VCstem_P50[index], df_res$VCstem_slope[index])
    fact <- 100/max(mdata1$M1_gs_pred)
    g2 <- g2 +
      geom_line(aes(x=Pleaf, y=PLC*max(M1_gs_pred)), data = mdata1, col ="gray")+
      scale_x_continuous(limits = c(df_res$VCstem_P50[index]-1.0,0))+
      scale_y_continuous(sec.axis = sec_axis(~ . * fact, name = "PLC (%)")) +
      geom_vline(xintercept = df_res$VCstem_P50[index], col = "gray", size = 2, alpha = 0.5)
  }
  g2 <- g2 +
    theme_bw()+
    theme(legend.position = c(0.9, 0.2)) +
    labs(subtitle = paste(site_vec[index], "/", df_res$Species_accepted[index]))

  mdata2 <- data.frame(Pleaf = 0.0,
                       An_Cs = xydata$An_Cs) |>
    dplyr::mutate(M1_gs_pred = gs(An_Cs, Pleaf, 
                                  df_res$M1_g0[index], df_res$M1_g1[index], 
                                  df_res$M1_gs_P50[index], df_res$M1_gs_slope[index])) |>
    dplyr::mutate(M2_gs_pred = gs(An_Cs, Pleaf, 
                                  df_res$M2_g0[index], df_res$M2_g1[index], 
                                  df_res$M2_gs_P50[index], df_res$M2_gs_slope[index]))
  if(is_m3_m4) {
    mdata2 <- mdata2 |>
      dplyr::mutate(M3_gs_pred = gs(An_Cs, Pleaf, 
                                    df_res$M3_g0[index], df_res$M3_g1[index], 
                                    df_res$M3_gs_P50[index], df_res$M3_gs_slope[index])) |>
      dplyr::mutate(M4_gs_pred = gs(An_Cs, Pleaf, 
                                    df_res$M4_g0[index], df_res$M4_g1[index], 
                                    df_res$M4_gs_P50[index], df_res$M4_gs_slope[index]))
  }
  g3<- ggplot(xydata)+
    geom_point(aes(x = An_Cs, y = gs, col = Pleaf))+
    geom_line(aes(x=An_Cs, y=M1_gs_pred), data = mdata2, col="black") +
    geom_line(aes(x=An_Cs, y=M2_gs_pred), data = mdata2, col="red")
  if(is_m3_m4) {
    g3 <- g3 +
     geom_line(aes(x=An_Cs, y=M3_gs_pred), data = mdata2, col="blue")
  }
  g3 <- g3+
    xlab("An/Cs")+
    ylab("Observed stomatal conductance")+
    theme_bw()+
    theme(legend.position = c(0.9, 0.2)) +
    labs(subtitle = " ")

  g <- plot_grid(g2, g3, g1, nrow = 1)
  ggsave2(paste0("plots/", site_sp_vec[index], ".png"), g,
          width = 15, height  = 5)
}

df_res_out <- df_res  |>
  dplyr::filter(!is.na(n))

write.csv2(df_res_out, file = "data/gs_model_fit.csv", row.names = FALSE)
