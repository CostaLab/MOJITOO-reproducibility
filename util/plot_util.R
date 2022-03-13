reductions <- c("MOFARed",
                "schema",
                "WNN",
                "scAI",
                "liger",
                "MOJITOO",
                "DIABLORed",
                "symphonyRed"
)

color <- c(
  "MOFARed"             = "#FED439FF",
  "schema"              = "#709AE1FF",
  "WNN"                 = "#8A9197FF",
  "scAI"                = "#D2AF81FF",
  "liger"               = "#FD7446FF",
  "MOJITOO"             = "#D5E4A2FF",
  "DIABLORed"           = "#C80813FF",
  "symphonyRed"         = "#1A9993FF",
  "MOJITOOAll"          = "#075149FF",
  "MOFA"                = "#46732EFF",
  "MOJITOORaw"          = "#91331FFF"

)


pbmc_color <- c("CD4 Naive"                = "#456883",
                "CD4 Memory"               = "#f0e685",
                "CD8 Naive"                = "#5cb1dd",
                "CD16+ Monocytes"          = "#739b57",
                "NK cell"                  = "#d595a7",
                "CD14+ Monocytes"          = "#ce3c31",
                "B cell progenitor"        = "#4f4ffe",
                "Dendritic cell"           = "#7f2167",
                "CD8 effector"             = "#ba6237",
                "Double negative T cell"   = "#6ad76a",
                "pre-B cell"               = "#c75026",
                "pDC"                      = "#924721",
                "Platelets"                = "#837a8d"
)


asterisk <- function(x){
  if(x<0.01){
    return("***")
  }else if(x < 0.05){
    return("**")
  }else if(x < 0.1){
    return("*")
  }else{
    return("")
  }
  return("")
}

