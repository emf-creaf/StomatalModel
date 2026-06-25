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


WFO_file <- "~/OneDrive/EMF_datasets/Taxonomy/WFO/WFO_Backbone/classification.csv"
DB_path <- "~/OneDrive/mcaceres_work/model_development/medfate_parameterization/traits_and_models/"
harmonized_trait_path <- paste0(DB_path,"data/harmonized_trait_sources")
db <- data.frame(originalName = sp_vec)
db_post <- traits4models::harmonize_taxonomy_WFO(db, WFO_backbone_file = WFO_file)
trait_obs  <- traits4models::get_taxon_trait_means(harmonized_trait_path, traits = c("Nleaf", "Vmax", "Ptlp", "Gs_P50", "VCstem_P50", "LeafPI0", "LeafEPS", "Gswmin"))

df_res <- data.frame(Site = site_vec,
                     Species_original = sp_vec,
                     Species_accepted = db_post$acceptedName,
                     n = NA,
                     Pleaf_min = NA,
                     Pleaf_max = NA,
                     mae = NA,
                     valid = NA,
                     gs_max = NA,
                     g0 = NA,
                     g1 = NA,
                     gs_P50 = NA,
                     gs_slope = NA)

## Stomatal conductance function
gs <- function(x, y, g0, g1, gs_P50, gs_slope) {
  fl <- 1 - 1/(1+exp((gs_slope/25)*(y - gs_P50)))
  fl*pmax(g0, g0 + g1*x)
}

## Optimization function (all parameters)
f_opt_g0_g1_P50_slope <- function(x, xydata) {
  diff <- xydata$gs - gs(xydata$An_Cs, xydata$Pleaf, x[1], x[2], x[3], x[4])
  return(- mean(abs(diff)))
}
f_opt_g1_P50_slope <- function(x, xydata, g0) {
  diff <- xydata$gs - gs(xydata$An_Cs, xydata$Pleaf, g0, x[1], x[2], x[3])
  return(- mean(abs(diff)))
}



for(index in valid) {
  print(site_sp_vec[index])
  xdata <- readr::read_csv(paste0("raw-data/ManonSabot-One_gs_model_to_rule_them_all-41cbd30/input/calibrations/obs_driven/", site_sp_vec[index], "_x.csv"), show_col_types = FALSE)
  ydata <- readr::read_csv(paste0("raw-data/ManonSabot-One_gs_model_to_rule_them_all-41cbd30/input/calibrations/obs_driven/", site_sp_vec[index], "_y.csv"), show_col_types = FALSE)

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

  df_res$n[index] <- nrow(xydata)
  df_res$gs_max[index] <- quantile(xydata$gs, 0.95)
  df_res$Pleaf_min[index] <- min(xydata$Pleaf)
  df_res$Pleaf_max[index] <- max(xydata$Pleaf)

  ga_res <- ga(type = "real-valued",
               fitness = f_opt_g1_P50_slope, xydata = xydata, g0 = df_res$gs_max[index]*0.05,
               popSize = 300,
               maxiter = 400,
               lower= c(1, df_res$Pleaf_min[index], 25), 
               upper = c(20, df_res$Pleaf_max[index], 90),
               optim = TRUE, monitor = FALSE)
  

  df_res$mae[index] <- -ga_res@fitnessValue
  df_res$g0[index] <-df_res$gs_max[index]*0.05 #ga_res@solution[1]
  df_res$g1[index] <-ga_res@solution[1,1]
  df_res$gs_P50[index] <-ga_res@solution[1,2]
  df_res$gs_slope[index] <- ga_res@solution[1,3]
  
  df_res$valid[index] <- (df_res$gs_P50[index]> df_res$Pleaf_min[index]) & (df_res$gs_slope[index]>5) 
  print(df_res[index,])

  
  
  xydata$gs_pred <- gs(xydata$An_Cs,xydata$Pleaf, 
                       df_res$g0[index], df_res$g1[index], df_res$gs_P50[index], df_res$gs_slope[index])
  
  g1 <- ggplot(xydata)+
    geom_point(aes(x = gs, y = gs_pred))+
    geom_abline(intercept = 0, slope = 1, col = "gray")+
    scale_x_continuous(limits = c(0,max(xydata$gs)))+
    scale_y_continuous(limits = c(0,max(xydata$gs)))+
    theme_bw()

  
  mdata1 <- data.frame(Pleaf = seq(-9, 0, by=0.01)) |>
    dplyr::mutate(sigmoid = (1 - 1/(1+exp((df_res$gs_slope[index]/25)*(Pleaf -  df_res$gs_P50[index])))),
                  gs_pred = max(xydata$gs)*sigmoid)
  
  g2 <- ggplot(xydata)+
    geom_point(aes(x = Pleaf, y = gs, col = An_Cs))+
    geom_line(aes(x=Pleaf, y=gs_pred), data = mdata1)+
    scale_x_continuous(limits = c(-6,0))+
    theme_bw()+
    theme(legend.position = c(0.9, 0.2))

  mdata2 <- data.frame(An_Cs = seq(0, max(xydata$An_Cs), by=0.001)) |>
    dplyr::mutate(gs_pred = pmax(df_res$g0[index], df_res$g0[index] + df_res$g1[index]*An_Cs))
  
  g3<- ggplot(xydata)+
    geom_point(aes(x = An_Cs, y = gs, col = Pleaf))+
    geom_line(aes(x=An_Cs, y=gs_pred), data = mdata2) +
    theme_bw()+
    theme(legend.position = c(0.9, 0.2))
  
  g <- plot_grid(g1,g2, g3, nrow = 1)
  ggsave2(paste0("plots/", site_sp_vec[index], ".png"), g,
          width = 12, height  = 6)
}



write.csv2(df_res  |>
             dplyr::left_join(trait_obs, by=c("Species_accepted"="acceptedName")) |>
             dplyr::filter(!is.na(n)), file = "data/gs_model_fit.csv", row.names = FALSE)
