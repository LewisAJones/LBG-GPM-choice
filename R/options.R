# Header ----------------------------------------------------------------
# Project: LBG-GPM-choice
# File name: options.R
# Last updated: 2023-03-06
# Repository: https://github.com/LewisAJones/LBG-GPM-choice
# Parameters ------------------------------------------------------------
params <- list(
  # Which Global Plate Models should be used?
  models = c("PALEOMAP", "GOLONKA", "MERDITH2021"),
  # Naming convention for longitude
  lng = "lng",
  # Naming convention for latitude
  lat = "lat",
  # Naming convention for palaeolongitude
  p_lng = "p_lng",
  # Naming convention for palaeolatitude
  p_lat = "p_lat",
  # Number of latitudinal bins desired
  lat_bin = 12,
  # Naming convention for temporal bin midpoint
  age = "bin_midpoint",
  # The geological rank for conducting analyses
  rank = "stage",
  # The Geological Time Scale to be used
  GTS = "GTS2020",
  # How should occurrences be temporally binned?
  method = "majority",
  # Threshold for majority binning rule
  threshold = 90,
  # PBDB API link for occurrences
  pbdb_api = paste0("https://paleobiodb.org/data1.2/occs/list.csv?base_name=",
                    "Bivalvia,Brachiopoda,Cephalopoda,Gastropoda,Trilobita",
                    "&taxon_reso=genus",
                    "&ident=latest",
                    "&taxon_status=valid",
                    "&idqual=genus_certain",
                    "&pres=regular",
                    "&interval=Fortunian,Holocene",
                    "&envtype=marine",
                    "&show=genus,pres,strat,coll,coords,paleoloc,loc,class"),
  # Should a new version of PBDB data be downloaded?
  download = FALSE,
  # Should results be plotted?
  plot = FALSE,
  # Notify when script finished running?
  notify = TRUE
)
