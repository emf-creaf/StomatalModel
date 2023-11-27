library(tidyverse)
library(ggplot2)

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
valid = c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16)
index = 9
print(site_sp_vec[index])
xdata <- readr::read_csv(paste0("Data/ManonSabot-One_gs_model_to_rule_them_all-41cbd30/input/calibrations/obs_driven/", site_sp_vec[index], "_x.csv"))
ydata <- readr::read_csv(paste0("Data/ManonSabot-One_gs_model_to_rule_them_all-41cbd30/input/calibrations/obs_driven/", site_sp_vec[index], "_y.csv"))
xunits <- xdata[1,]
yunits <- ydata[1,]

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
         An_Cs = A/Cs)
  
ggplot(xydata)+
  geom_point(aes(x = An_Cs, y = gs, col = Pleaf)) 
ggplot(xydata)+
  geom_point(aes(x = Pleaf, y = gs, col = An_Cs)) 
summary(lm(gs ~  An_Cs, data = xydata))

m <- function(x, y, g0, g1, gs_P50, gs_slope) {
  fl <- 1/(1+exp((gs_slope/25)*(y - gs_P50)))
  fl*max(g0, g0 + g1*x)
}

m(0.02, -1)