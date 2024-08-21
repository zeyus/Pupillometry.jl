# check if pacman is installed
if (!require("pacman")) install.packages("pacman")
library("pacman")

pacman::p_load(
  "tidyverse",
  "devtools",
  "readxl",
  "interp"
)
if (!require("pupillometry")) devtools::install_github("dr-JT/pupillometry")

library("pupillometry")


# Load the data
data_path <- "D:/data/pupillometry/"

# load the data from
#  - data_path/[SUBJ]/[SUBJ]_dx.xlsx (left eye)
#  - data_path/[SUBJ]/[SUBJ]_sx.xlsx (right eye)
# subjects 1M/F to 5M/F are LL lighting, 6M/F to 10M/F are HL lighting

# excel sheet "baseline" contains the baseline data
# excel sheet "audio" contains the task data
# load the data
data <- list()

rename_cols <- function(df, rename_map) {
  for (old_name in names(rename_map)) {
    new_name <- rename_map[[old_name]]
    if (old_name %in% names(df)) {
      names(df)[names(df) == old_name] <- new_name
    }
  }
  return(df)
}

possibl_col_renames <- list(
  "diameter (pixel)" = "diameter_px",
  "blink (1 yes/0 no)" = "blink",
  "blink (1 sì/0 no)" = "blink",
  "audio (1 yes/0 no)" = "audio",
  "audio (1 sì/0 no)" = "audio",
  "artifact (1 yes/0 no)" = "artifact",
  "artefatto (1 sì/0 no)" = "artifact"
)

for (subj in 1:10) {
  for (s in c("M", "F")) {
    subj_str <- sprintf("%d%s", subj, s)
    dx_file <- sprintf("%s/%s/%s_dx.xlsx", data_path, subj_str, subj_str)
    sx_file <- sprintf("%s/%s/%s_sx.xlsx", data_path, subj_str, subj_str)
    lighting <- ifelse(subj <= 5, "LL", "HL")
    dxa <- readxl::read_excel(dx_file, sheet = "audio")
    dxa <- rename_cols(dxa, possibl_col_renames)
    # add lighting column
    dxa$lighting <- lighting
    if (s == "M") {
      subj <- subj + 10
    }
    # add subject column
    dxa$subj <- subj
    # rename diameter_px to L_Pupil_Diameter.px
    dxa <- dxa %>% rename(
      L_Pupil_Diameter.px = diameter_px,
      L_Blink = blink
    )
    # add condition column
    dxa$condition <- "audio"

    sxa <- readxl::read_excel(sx_file, sheet = "audio")
    sxa <- rename_cols(sxa, possibl_col_renames)
    sxa <- sxa %>% rename(
      R_Pupil_Diameter.px = diameter_px,
      R_Blink = blink
    )

    # get dxa milliseconds not in sxa
    missing_l <- setdiff(dxa$milliseconds, sxa$milliseconds)

    linterp <- interp::interp(
      x = sxa$milliseconds,
      y = sxa$R_Pupil_Diameter.px,
      xo = missing_l,
      output = "point",
    )

    # get sxa milliseconds not in dxa
    missing_r <- setdiff(sxa$milliseconds, dxa$milliseconds)
    rinterp <- interp::interp(
      x = dxa$milliseconds,
      y = dxa$L_Pupil_Diameter.px,
      xo = missing_r,
      output = "point",
    )
    sxa$L_Pupil_Diameter.px <- rinterp$y



    


    dxb <- readxl::read_excel(dx_file, sheet = "baseline")
    dxb <- rename_cols(dxb, possibl_col_renames)
    dxb$lighting <- lighting
    dxb$subj <- subj
    dxb$condition <- "baseline"
    dxb <- dxb %>% rename(
      L_Pupil_Diameter.px = diameter_px,
      L_Blink = blink
    )

    sxb <- readxl::read_excel(sx_file, sheet = "baseline")
    sxb <- rename_cols(sxb, possibl_col_renames)
    sxb <- sxb %>% rename(
      R_Pupil_Diameter.px = diameter_px,
      R_Blink = blink
    )

    # rename the columns if they exist
    data[[sprintf("%s_audio", subj_str)]] <- dxa
    data[[sprintf("%s_baseline", subj_str)]] <- dxb
  }
}


# combine the data
all_data <- bind_rows(data)

# update to 
# left eye diameter_px -> L_Pupil_Diameter.px
# left eye blink -> L_Blink
# right eye diameter_px -> R_Pupil_Diameter.px
# right eye blink -> R_Blink
# milliseconds/1000 -> Time

all_data <- all_data %>%
  mutate(
    Time = milliseconds / 1000
  ) %>%
  select(-milliseconds) %>%
  rename(
    L_Pupil_Diameter.px = diameter_px,
    R_Pupil_Diameter.px = diameter_px,
    L_Blink = blink,
    R_Blink = blink
  )







