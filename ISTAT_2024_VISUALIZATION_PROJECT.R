
RUN_MODELLING <- TRUE 
RUN_CLUSTERING <- TRUE  
SHOW_FINAL_PLOTS <- TRUE  
EXPORT_FINAL_PLOTS <- TRUE 

# 00 | SETUP -----------------------------------------------------------------

required_packages <- c("readr", "dplyr", "tidyr", "ggplot2", "stringr", "scales",
                       "sf", "patchwork", "cluster", "ggrepel", "fmsb")
missing_packages <- setdiff(required_packages, rownames(installed.packages()))
if (length(missing_packages)) {
  stop("Missing R packages: ", paste(missing_packages, collapse = ", "),
       ". Install them with install.packages().")
}
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(ggplot2)
  library(stringr); library(scales); library(sf); library(patchwork)
})

DATA_FILE <- "data_raw/INCSTRAD_Microdati_2024.txt"
GEO_ZIP <- "data_raw/Limiti01012024_g.zip"
GEO_DIR <- "data_raw/istat_2024_geography"
GEO_URL <- paste0("https://www.istat.it/storage/cartografia/",
                  "confini_amministrativi/generalizzati/2024/Limiti01012024_g.zip")
MIN_N_PLOT <- 100L       # Plotting threshold

label_istat <- function(x, dictionary, name) {
  values <- unname(dictionary[as.character(x)])
  unknown <- unique(as.character(x[!is.na(x) & is.na(values)]))
  if (length(unknown)) stop(name, " has unmapped codes: ", paste(unknown, collapse = ", "))
  values
}

# vehicle grouping
vehicle_group_from_code <- function(x) {
  dplyr::case_when(
    x %in% sprintf("%02d", 1:4) ~ "Cars",
    x %in% c("05", "06") ~ "Buses",
    x %in% c("08", "09", "10", "12", "18") ~ "Goods/commercial vehicles",
    x %in% c("14", "23") ~ "Bicycles/e-bikes",
    x %in% c("15", "16", "17") ~ "Motorcycles/mopeds",
    x == "22" ~ "E-scooters",
    TRUE ~ "Other/unknown"
  )
}

vehicle_names <- c(
  "01"="Private car", "02"="Private car with trailer",
  "03"="Public car/taxi", "04"="Emergency/police car",
  "05"="Urban bus/trolleybus", "06"="Extra-urban bus",
  "07"="Tram", "08"="Truck", "09"="Truck with trailer",
  "10"="Articulated goods vehicle", "11"="Special vehicle",
  "12"="Road tractor", "13"="Agricultural tractor",
  "14"="Bicycle", "15"="Moped", "16"="Motorcycle (solo)",
  "17"="Motorcycle with passenger", "18"="Motor-truck/van",
  "19"="Animal-/hand-drawn vehicle", "20"="Unidentified hit-and-run vehicle",
  "21"="Quadricycle", "22"="Electric kick scooter", "23"="E-bike"
)

codebooks <- list(
  localizzazione_incidente = c(
    "0"="Regional road within built-up area", "1"="Urban road",
    "2"="Provincial road within built-up area",
    "3"="State road within built-up area", "4"="Extra-urban municipal road",
    "5"="Provincial road", "6"="State road", "7"="Motorway",
    "8"="Other road", "9"="Regional road"),
  tipo_di_strada = c(
    "1"="One-way single carriageway", "2"="Two-way single carriageway",
    "3"="Two carriageways", "4"="More than two carriageways"),
  intersezione_o_non_interse3 = c(
    "01"="Crossroads", "02"="Roundabout", "03"="Signed intersection",
    "04"="Traffic-light/police-controlled intersection",
    "05"="Unsigned intersection", "06"="Level crossing", "07"="Straight road",
    "08"="Bend", "09"="Crest/narrowing", "10"="Slope",
    "11"="Lit tunnel", "12"="Unlit tunnel"),
  pavimentazione = c("1"="Paved road", "2"="Damaged paved road", "3"="Unpaved road"),
  fondo_stradale = c("1"="Dry", "2"="Wet", "3"="Slippery", "4"="Icy",
                    "5"="Snow-covered"),
  segnaletica = c("1"="Absent", "2"="Vertical signs", "3"="Road markings",
                 "4"="Both signs and road markings", "5"="Temporary roadwork signs"),
  condizioni_meteorologiche = c(
    "1"="Clear", "2"="Fog", "3"="Rain", "4"="Hail", "5"="Snow",
    "6"="Strong wind", "7"="Other"),
  natura_incidente = c(
    "01"="Head-on collision", "02"="Front-to-side collision",
    "03"="Side collision", "04"="Rear-end collision",
    "05"="Pedestrian collision", "06"="Collision with stopping vehicle",
    "07"="Collision with parked vehicle", "08"="Collision with obstacle",
    "09"="Collision with train", "10"="Run-off-road/skidding",
    "11"="Injury from sudden braking", "12"="Injury from falling from vehicle")
)

weekday_levels <- c("Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday", "Sunday")
quarter_levels <- c("Jan–Mar", "Apr–Jun", "Jul–Sep", "Oct–Dec")
time_band_levels <- c("00–05", "06–11", "12–17", "18–23", "Not reported")
age_levels <- c("0-5", "6-9", "10-14", "15-17", "18-29", "30-44",
                "45-54", "55-64", "65+", "Unknown")
selected_roads <- c("Urban road", "Provincial road", "State road", "Motorway")
selected_natures <- c("Front-to-side collision", "Rear-end collision",
                      "Side collision", "Head-on collision", "Pedestrian collision",
                      "Run-off-road/skidding", "Collision with obstacle")
selected_vehicle_groups <- c("Cars", "Motorcycles/mopeds",
                             "Goods/commercial vehicles", "Bicycles/e-bikes")
# 01 | DATA PREPARATION ------------------------------------------------------


if (!file.exists(DATA_FILE)) {
  stop("Cannot find ISTAT microdata at '", DATA_FILE,
       "'. Set the RStudio working directory to the project root.")
}
accidents_raw <- readr::read_tsv(DATA_FILE,
  col_types = cols(.default = col_character()),
  na = "", trim_ws = TRUE, show_col_types = FALSE)
stopifnot(ncol(accidents_raw) == 120L)

variable_overview <- tibble(
  variable = names(accidents_raw),
  observed = as.integer(colSums(!is.na(accidents_raw))),
  distinct = vapply(accidents_raw, function(x) length(unique(stats::na.omit(x))), integer(1)),
  missing_pct = round(100 * colMeans(is.na(accidents_raw)), 2)
)
empty_cols <- variable_overview$variable[variable_overview$observed == 0L]

# 01.1 | Incident-level table ------------------------------------------------
accidents_main <- accidents_raw %>%
  mutate(
    incident_id = row_number(),
    provincia, comune, comune_istat = paste0(provincia, comune),
    giorno, trimestre,
    hour = as.integer(na_if(Ora, "25")),  # Missing/imprecise, NOT an extra hour.
    weekend = giorno %in% c("6", "7"),
    localizzazione_incidente, tipo_di_strada,
    intersezione_o_non_interse3, pavimentazione,
    fondo_stradale, segnaletica, condizioni_meteorologiche, natura_incidente,
    tipo_veicolo_a, tipo_veicoli__b_, tipo_veicolo__c_,
    other_vehicles = as.integer(altri_veicoli_coinvolti),
    deaths_24h = as.integer(morti_entro_24_ore),
    deaths_2_30d = as.integer(morti_entro_30_giorni),
    total_injured = as.integer(feriti)
  ) %>%
  mutate(
    total_deaths = deaths_24h + deaths_2_30d,
    fatal_accident = total_deaths > 0,
    n_vehicles_abc = as.integer(!is.na(tipo_veicolo_a)) +
      as.integer(!is.na(tipo_veicoli__b_)) +
      as.integer(!is.na(tipo_veicolo__c_)),
    n_vehicles = if_else(is.na(other_vehicles), NA_integer_,
                         n_vehicles_abc + other_vehicles),
    day_label = factor(as.integer(giorno), levels = 1:7,
                       labels = weekday_levels, ordered = TRUE),
    quarter_label = factor(as.integer(trimestre), levels = 1:4,
                           labels = quarter_levels, ordered = TRUE),
    hour_label = if_else(is.na(hour), "Not reported", sprintf("%02d", hour)),
    # Final corrected bands, defined ONCE, not redefined in the EDA.
    time_band = case_when(
      is.na(hour) ~ "Not reported",
      hour <= 5 ~ "00–05", hour <= 11 ~ "06–11",
      hour <= 17 ~ "12–17", hour <= 23 ~ "18–23",
      TRUE ~ NA_character_
    ),
    time_band = factor(time_band, levels = time_band_levels, ordered = TRUE)
  )
stopifnot(all(is.na(accidents_main$hour) |
              accidents_main$hour %in% 0:23))
for (v in names(codebooks)) {
  accidents_main[[paste0(v, "_label")]] <-
    label_istat(accidents_main[[v]], codebooks[[v]], v)
}
accidents_main <- accidents_main %>% mutate(
  vehicle_A_label = label_istat(tipo_veicolo_a, vehicle_names, "vehicle A"),
  vehicle_B_label = label_istat(tipo_veicoli__b_, vehicle_names, "vehicle B"),
  vehicle_C_label = label_istat(tipo_veicolo__c_, vehicle_names, "vehicle C")
)

# 01.2 | Vehicle-level table -------------------------------------------------
vehicles_long <- accidents_main %>%
  select(incident_id, A = tipo_veicolo_a, B = tipo_veicoli__b_, C = tipo_veicolo__c_) %>%
  pivot_longer(c(A, B, C), names_to = "vehicle_slot",
               values_to = "vehicle_code", values_drop_na = TRUE) %>%
  mutate(
    vehicle_group = vehicle_group_from_code(vehicle_code),
    vehicle_label = label_istat(vehicle_code, vehicle_names, "vehicle code")
  )
vehicle_flags <- vehicles_long %>%
  group_by(incident_id) %>%
  summarise(
    has_car = any(vehicle_group == "Cars"),
    has_motorcycle = any(vehicle_group == "Motorcycles/mopeds"),
    has_bicycle = any(vehicle_group == "Bicycles/e-bikes"),
    has_goods_vehicle = any(vehicle_group == "Goods/commercial vehicles"),
    .groups = "drop"
  )
accidents_main <- accidents_main %>% left_join(vehicle_flags, by = "incident_id")

# 01.3 | Driver-level table --------------------------------------------------
driver_fields <- list(
  A = c(vehicle="tipo_veicolo_a", age="veicolo__a___et__conducente",
        sex="veicolo__a___sesso_conducente", outcome="veicolo__a___esito_conducente",
        license="veicolo__a___anno_rilascio30"),
  B = c(vehicle="tipo_veicoli__b_", age="veicolo__b___et__conducente",
        sex="veicolo__b___sesso_conducente", outcome="veicolo__b___esito_conducente",
        license="veicolo__b___anno_rilascio52"),
  C = c(vehicle="tipo_veicolo__c_", age="veicolo__c___et__conducente",
        sex="veicolo__c___sesso_conducente", outcome="veicolo__c___esito_conducente",
        license="veicolo__c___anno_rilascio74")
)
drivers_long <- bind_rows(lapply(names(driver_fields), function(slot) {
  f <- driver_fields[[slot]]
  accidents_raw %>% transmute(
    incident_id = row_number(), vehicle_slot = slot,
    vehicle_code = .data[[f[["vehicle"]]]],
    age_group = .data[[f[["age"]]]], sex_code = .data[[f[["sex"]]]],
    outcome_code = .data[[f[["outcome"]]]],
    license_issue_group = .data[[f[["license"]]]]
  )
})) %>%
  filter(!is.na(age_group) | !is.na(sex_code) | !is.na(outcome_code)) %>%
  mutate(
    sex = case_when(sex_code == "1" ~ "Male", sex_code == "2" ~ "Female",
                    TRUE ~ NA_character_),
    outcome = case_when(outcome_code == "1" ~ "Uninjured",
                        outcome_code == "2" ~ "Injured",
                        outcome_code %in% c("3", "4") ~ "Killed",
                        TRUE ~ NA_character_),
    age_group_clean = na_if(age_group, "n.i."),
    vehicle_group = vehicle_group_from_code(vehicle_code)
  ) %>%
  left_join(vehicles_long %>% select(incident_id, vehicle_slot, vehicle_label),
            by = c("incident_id", "vehicle_slot"))

# Preserve the original licence class, but DO NOT analyse its duration as actual driving experience.

# 01.4 | Pedestrian casualties -----------------------------------------------
pedestrians_long <- bind_rows(lapply(c("morto", "ferito"), function(status) {
  bind_rows(lapply(1:4, function(slot) {
    prefix <- paste0("pedone_", status, "_", slot, "__")
    tibble(
      incident_id = seq_len(nrow(accidents_raw)), pedestrian_slot = slot,
      outcome = if (status == "morto") "Killed" else "Injured",
      sex_code = accidents_raw[[paste0(prefix, "sesso")]],
      age_group = accidents_raw[[paste0(prefix, "et_")]]
    )
  }))
})) %>%
  filter(!is.na(sex_code) | !is.na(age_group)) %>%
  mutate(
    sex = case_when(
      outcome == "Killed" & sex_code == "1" ~ "Male",
      outcome == "Killed" & sex_code == "2" ~ "Female",
      outcome == "Injured" & sex_code == "3" ~ "Male",
      outcome == "Injured" & sex_code == "4" ~ "Female",
      TRUE ~ NA_character_),
    age_group_clean = na_if(age_group, "n.i.")
  )

# 01.5 | Passenger and other-occupant totals ---------------------------------
# Explicit variable names replace brittle raw positions [31,34,...].
passenger_outcome_slots <- c(
  "veicolo__a___passeggeri_an35", "veicolo__a___esito_passegg38",
  "veicolo__a___esito_passegg41", "veicolo__a___esito_passegg44",
  "veicolo__b___passeggeri_an57", "veicolo__b___esito_passegg60",
  "veicolo__b___esito_passegg63", "veicolo__b___esito_passegg66",
  "veicolo__c___passeggeri_an79", "veicolo__c___esito_passegg82",
  "veicolo__c___esito_passegg85", "veicolo__c___esito_passegg88"
)
extra_passenger_deaths <- c(
  "veicolo__a___altri_passegg47", "veicolo__a___altri_passegg48",
  "veicolo__b___altri_passegg69", "veicolo__b___altri_passegg70",
  "veicolo__c___altri_passegg91", "veicolo__c___altri_passegg92"
)
extra_passenger_injuries <- c(
  "veicolo__a___altri_passegg49", "veicolo__a___altri_passegg50",
  "veicolo__b___altri_passegg71", "veicolo__b___altri_passegg72",
  "veicolo__c___altri_passegg93", "veicolo__c___altri_passegg94"
)
other_deaths <- c("morti_maschi_coinvolti_su_95", "morti_femmine_coinvolti_su96")
other_injuries <- c("feriti_maschi_coinvolti_su97", "feriti_femmine_coinvolti_s98")
expected_fields <- c(passenger_outcome_slots, extra_passenger_deaths,
                     extra_passenger_injuries, other_deaths, other_injuries)
stopifnot(all(expected_fields %in% names(accidents_raw)))
count_code <- function(cols, code) {
  rowSums(as.matrix(accidents_raw[, cols, drop = FALSE]) == code, na.rm = TRUE)
}
sum_counts <- function(cols) {
  x <- sapply(accidents_raw[, cols, drop = FALSE], as.integer)
  if (is.null(dim(x))) x <- matrix(x, ncol = 1L)
  rowSums(x, na.rm = TRUE)
}
passenger_counts <- tibble(
  incident_id = seq_len(nrow(accidents_raw)),
  passengers_killed = count_code(passenger_outcome_slots, "1") +
    sum_counts(extra_passenger_deaths),
  passengers_injured = count_code(passenger_outcome_slots, "2") +
    sum_counts(extra_passenger_injuries)
)
other_person_counts <- tibble(
  incident_id = seq_len(nrow(accidents_raw)),
  other_vehicle_people_killed = sum_counts(other_deaths),
  other_vehicle_people_injured = sum_counts(other_injuries)
)

driver_counts <- drivers_long %>% group_by(incident_id) %>% summarise(
  drivers_killed = sum(outcome == "Killed", na.rm = TRUE),
  drivers_injured = sum(outcome == "Injured", na.rm = TRUE), .groups = "drop")
pedestrian_counts <- pedestrians_long %>% group_by(incident_id) %>% summarise(
  pedestrians_killed = sum(outcome == "Killed"),
  pedestrians_injured = sum(outcome == "Injured"), .groups = "drop")

accidents_main <- accidents_main %>%
  left_join(driver_counts, by = "incident_id") %>%
  left_join(pedestrian_counts, by = "incident_id") %>%
  left_join(passenger_counts, by = "incident_id") %>%
  left_join(other_person_counts, by = "incident_id") %>%
  mutate(across(c(drivers_killed, drivers_injured,
                  pedestrians_killed, pedestrians_injured), ~ coalesce(.x, 0L)),
         reconstructed_deaths = drivers_killed + passengers_killed +
           pedestrians_killed + other_vehicle_people_killed,
         reconstructed_injured = drivers_injured + passengers_injured +
           pedestrians_injured + other_vehicle_people_injured)

# 01.6 | ISTAT geographical labels ------------------------------------------
# Statistics preserve 2024 municipality/region codes; three municipalities
# created after 01/01/2024 are not in the early-2024 map geometries.
if (!file.exists(GEO_ZIP)) {
  dir.create(dirname(GEO_ZIP), recursive = TRUE, showWarnings = FALSE)
  tryCatch(utils::download.file(GEO_URL, GEO_ZIP, mode = "wb", quiet = TRUE),
           error = function(e) message("Map download failed: ", conditionMessage(e)))
}
if (!dir.exists(GEO_DIR)) dir.create(GEO_DIR, recursive = TRUE)
if (file.exists(GEO_ZIP) && !length(list.files(GEO_DIR, "\\.shp$", recursive = TRUE))) {
  utils::unzip(GEO_ZIP, exdir = GEO_DIR)
}
shape_files <- list.files(GEO_DIR, pattern = "\\.shp$", full.names = TRUE,
                          recursive = TRUE)
find_shape <- function(prefix) {
  match <- shape_files[grepl(prefix, basename(shape_files))]
  if (!length(match)) stop("Missing ISTAT 2024 shapefile: ", prefix,
                           ". Check GEO_ZIP / GEO_DIR.")
  match[1]
}
municipal_file <- find_shape("^Com01012024")
province_file <- find_shape("^Prov(CM)?01012024")
region_file <- find_shape("^Reg01012024")
municipalities_sf <- st_read(municipal_file, quiet = TRUE)
provinces_sf <- st_read(province_file, quiet = TRUE)
regions_sf <- st_read(region_file, quiet = TRUE)

municipality_code_column <- intersect(c("PRO_COM_T", "PRO_COM"),
                                      names(municipalities_sf))[1]
province_code_column <- intersect(c("COD_PROV", "COD_UTS"),
                                  names(provinces_sf))[1]
province_name_column <- intersect(c("DEN_UTS", "DEN_PROV", "PROVINCIA"),
                                  names(provinces_sf))[1]
stopifnot(!is.na(municipality_code_column), !is.na(province_code_column),
          !is.na(province_name_column))

municipality_lookup <- municipalities_sf %>% st_drop_geometry() %>%
  transmute(
    comune_istat = str_pad(as.character(.data[[municipality_code_column]]), 6, pad = "0"),
    municipality_name = as.character(COMUNE),
    region_code = sprintf("%02d", as.integer(COD_REG)),
    macro_code = as.character(COD_RIP)
  ) %>% distinct(comune_istat, .keep_all = TRUE)
new_municipalities <- tibble(
  comune_istat = c("025075", "024128", "028108"),
  municipality_name = c("Setteville", "Sovizzo", "Santa Caterina d'Este"),
  region_code = c("05", "05", "05"), macro_code = c("2", "2", "2")
)
municipality_lookup <- municipality_lookup %>%
  filter(!comune_istat %in% new_municipalities$comune_istat) %>%
  bind_rows(new_municipalities)
province_lookup <- provinces_sf %>% st_drop_geometry() %>%
  transmute(
    provincia = sprintf("%03d", as.integer(.data[[province_code_column]])),
    province_name = as.character(.data[[province_name_column]])
  ) %>% distinct(provincia, .keep_all = TRUE)
region_names <- c(
  "01"="Piemonte", "02"="Valle d'Aosta", "03"="Lombardia",
  "04"="Trentino-Alto Adige", "05"="Veneto", "06"="Friuli-Venezia Giulia",
  "07"="Liguria", "08"="Emilia-Romagna", "09"="Toscana", "10"="Umbria",
  "11"="Marche", "12"="Lazio", "13"="Abruzzo", "14"="Molise",
  "15"="Campania", "16"="Puglia", "17"="Basilicata", "18"="Calabria",
  "19"="Sicilia", "20"="Sardegna"
)
macro_names <- c("1"="North-West", "2"="North-East", "3"="Centre",
                 "4"="South", "5"="Islands")
accidents_main <- accidents_main %>%
  left_join(municipality_lookup, by = "comune_istat") %>%
  left_join(province_lookup, by = "provincia") %>%
  mutate(region_name = unname(region_names[region_code]),
         macro_area = unname(macro_names[macro_code]))

# 02 | INTEGRITY TESTS AND DATA QUALITY ------------------------------------
quality_summary <- accidents_main %>% summarise(
  incidents = n(), deaths = sum(total_deaths), injured = sum(total_injured),
  fatal_accidents = sum(fatal_accident),
  zero_casualty_incidents = sum(total_deaths == 0 & total_injured == 0),
  death_count_mismatches = sum(reconstructed_deaths != total_deaths),
  injury_count_mismatches = sum(reconstructed_injured != total_injured),
  missing_total_vehicle_count = sum(is.na(n_vehicles))
)
stopifnot(
  nrow(accidents_main) == nrow(accidents_raw),
  nrow(accidents_main) == 173364L,
  quality_summary$deaths == 3030L,
  quality_summary$injured == 233853L,
  quality_summary$death_count_mismatches == 0L,
  quality_summary$injury_count_mismatches == 0L,
  quality_summary$missing_total_vehicle_count == 0L,
  quality_summary$zero_casualty_incidents == 0L,
  !anyDuplicated(accidents_main$incident_id),
  !anyDuplicated(drivers_long[c("incident_id", "vehicle_slot")]),
  !anyNA(accidents_main$municipality_name),
  !anyNA(accidents_main$province_name),
  !anyNA(accidents_main$region_name),
  !any(is.na(pedestrians_long$sex))
)
accidents_main <- accidents_main %>%
  select(-reconstructed_deaths, -reconstructed_injured)

# Identify mismatches
vehicles_without_drivers <- vehicles_long %>%
  anti_join(drivers_long %>% select(incident_id, vehicle_slot),
            by = c("incident_id", "vehicle_slot"))
drivers_without_vehicles <- drivers_long %>%
  anti_join(vehicles_long %>% select(incident_id, vehicle_slot),
            by = c("incident_id", "vehicle_slot"))
preparation_audit <- tibble(
  check = c("Duplicate incident IDs", "Duplicate driver ID/slot pairs",
            "Vehicles without driver details", "Drivers without a matching vehicle",
            "Unidentified hit-and-run vehicles without a driver",
            "Drivers with unknown age", "Drivers with missing vehicle code",
            "Pedestrians with unknown age", "Pedestrians with unknown sex",
            "Accidents without municipality name"),
  value = c(
    sum(duplicated(accidents_main$incident_id)),
    sum(duplicated(drivers_long[c("incident_id", "vehicle_slot")])),
    nrow(vehicles_without_drivers), nrow(drivers_without_vehicles),
    sum(vehicles_without_drivers$vehicle_code == "20", na.rm = TRUE),
    sum(is.na(drivers_long$age_group_clean)),
    sum(is.na(drivers_long$vehicle_code)),
    sum(is.na(pedestrians_long$age_group_clean)),
    sum(is.na(pedestrians_long$sex)),
    sum(is.na(accidents_main$municipality_name))
  )
)
# Expected from the supplied outputs: 2265 without driver, of which 2264
# are code 20; one driver in accident 80490 has an absent vehicle code.
stopifnot(nrow(vehicles_without_drivers) == 2265L,
          sum(vehicles_without_drivers$vehicle_code == "20") == 2264L,
          nrow(drivers_without_vehicles) == 1L)

# Retain useful diagnosis; 'n.i.' counts as unknown, not a valid known age.
driver_completeness <- drivers_long %>% summarise(
  drivers = n(), age_unknown = sum(is.na(age_group_clean)),
  sex_missing = sum(is.na(sex)), outcome_missing = sum(is.na(outcome_code)),
  license_group_missing = sum(is.na(license_issue_group)))
pedestrian_sex_summary <- pedestrians_long %>% count(outcome, sex, name = "people")
license_age_crosscheck <- drivers_long %>%
  filter(!is.na(age_group_clean)) %>%
  count(age_group_clean, license_issue_group, name = "drivers") %>%
  group_by(age_group_clean) %>%
  mutate(percentage = round(100 * drivers / sum(drivers), 2)) %>%
  ungroup()

cat("\nPREPARATION COMPLETE\n")
cat("Incidents:", nrow(accidents_main), "| Vehicles:", nrow(vehicles_long),
    "| Drivers:", nrow(drivers_long), "| Pedestrian casualties:",
    nrow(pedestrians_long), "\n")
print(quality_summary, width = Inf)
print(preparation_audit, n = Inf)

# 03 | DESCRIPTIVE EDA TABLES ----------------------------------------------

# 03.1 | Incident outcomes ---------------------------------------------------
death_distribution <- accidents_main %>% count(total_deaths, name = "incidents") %>%
  mutate(percentage = round(100 * incidents / sum(incidents), 2))
injury_distribution <- accidents_main %>%
  mutate(injury_group = if_else(total_injured >= 5L, "5+", as.character(total_injured)),
         injury_group = factor(injury_group, levels = c("0", "1", "2", "3", "4", "5+"))) %>%
  count(injury_group, name = "incidents") %>%
  mutate(percentage = round(100 * incidents / sum(incidents), 2))
vehicle_count_stats <- accidents_main %>% group_by(n_vehicles) %>% summarise(
  incidents = n(), fatal_accidents = sum(fatal_accident),
  fatal_rate = round(100 * mean(fatal_accident), 2), .groups = "drop")

# 03.2 | Time ---------------------------------------------------------------
day_stats <- accidents_main %>% group_by(day_label) %>% summarise(
  incidents = n(), fatal_accidents = sum(fatal_accident),
  total_deaths = sum(total_deaths), fatal_rate = round(100 * mean(fatal_accident), 2),
  .groups = "drop") %>% mutate(percentage = round(100 * incidents / sum(incidents), 2))
hour_stats <- accidents_main %>% filter(!is.na(hour)) %>% group_by(hour) %>% summarise(
  incidents = n(), fatal_accidents = sum(fatal_accident),
  fatal_rate = round(100 * mean(fatal_accident), 2), .groups = "drop")
quarter_stats <- accidents_main %>% group_by(quarter_label) %>% summarise(
  incidents = n(), fatal_accidents = sum(fatal_accident),
  total_deaths = sum(total_deaths), fatal_rate = round(100 * mean(fatal_accident), 2),
  .groups = "drop") %>% mutate(percentage = round(100 * incidents / sum(incidents), 2))
time_band_stats <- accidents_main %>% group_by(time_band) %>% summarise(
  incidents = n(), fatal_accidents = sum(fatal_accident),
  fatal_rate = round(100 * mean(fatal_accident), 2), .groups = "drop")
day_time_stats <- accidents_main %>% filter(time_band != "Not reported") %>%
  group_by(day_label, time_band) %>% summarise(
    incidents = n(), fatal_accidents = sum(fatal_accident),
    fatal_rate = round(100 * mean(fatal_accident), 2), .groups = "drop")

# Calendar-day normalization matters: 2024 has 53 Mondays and Tuesdays.
weekday_calendar <- tibble(
  day_label = factor(weekday_levels, levels = weekday_levels, ordered = TRUE),
  days_in_2024 = c(53, 53, rep(52, 5)))
weekday_normalized <- day_stats %>%
  left_join(weekday_calendar, by = "day_label") %>%
  mutate(accidents_per_day = round(incidents / days_in_2024, 1))
quarter_calendar <- tibble(
  quarter_label = factor(quarter_levels, levels = quarter_levels, ordered = TRUE),
  days_in_quarter = c(91, 91, 92, 92))
quarter_normalized <- quarter_stats %>%
  left_join(quarter_calendar, by = "quarter_label") %>%
  mutate(accidents_per_day = round(incidents / days_in_quarter, 1))
# Missing hour is retained as 'Not reported'; analyse reporting bias separately.
hour_missing_check <- accidents_main %>%
  mutate(hour_status = if_else(is.na(hour), "Not reported", "Reported")) %>%
  group_by(hour_status) %>% summarise(
    incidents = n(), fatal_accidents = sum(fatal_accident),
    fatal_rate = round(100 * mean(fatal_accident), 2), .groups = "drop")

# 03.3 | Environment, roads, accident nature --------------------------------
environment_variables <- c(
  "Weather"="condizioni_meteorologiche_label",
  "Road surface"="fondo_stradale_label", "Pavement"="pavimentazione_label",
  "Road signs"="segnaletica_label", "Road location"="localizzazione_incidente_label",
  "Road type"="tipo_di_strada_label",
  "Intersection / road layout"="intersezione_o_non_interse3_label",
  "Accident type"="natura_incidente_label")
environment_stats <- lapply(environment_variables, function(v) {
  accidents_main %>%
    mutate(category = coalesce(as.character(.data[[v]]), "Not reported")) %>%
    group_by(category) %>% summarise(
      incidents = n(), fatal_accidents = sum(fatal_accident),
      total_deaths = sum(total_deaths),
      fatal_rate = round(100 * mean(fatal_accident), 2), .groups = "drop") %>%
    mutate(percentage = round(100 * incidents / sum(incidents), 2)) %>%
    arrange(desc(incidents))
})
weather_surface_stats <- accidents_main %>%
  mutate(weather = coalesce(condizioni_meteorologiche_label, "Not reported"),
         road_surface = coalesce(fondo_stradale_label, "Not reported")) %>%
  group_by(weather, road_surface) %>% summarise(
    incidents = n(), fatal_accidents = sum(fatal_accident),
    fatal_rate = round(100 * mean(fatal_accident), 2), .groups = "drop") %>%
  arrange(desc(incidents))

# 03.4 | Vehicle counts, incident presence and exactly-two-vehicle pairs ----
vehicle_frequency <- vehicles_long %>% count(vehicle_group, name = "vehicles") %>%
  mutate(percentage = round(100 * vehicles / sum(vehicles), 2)) %>%
  arrange(desc(vehicles))
vehicle_severity <- vehicles_long %>% distinct(incident_id, vehicle_group) %>%
  left_join(accidents_main %>% select(incident_id, fatal_accident), by = "incident_id") %>%
  group_by(vehicle_group) %>% summarise(
    incidents = n(), fatal_accidents = sum(fatal_accident),
    fatal_rate = round(100 * mean(fatal_accident), 2), .groups = "drop") %>%
  arrange(desc(fatal_rate))
# Group incidence is overlapping: a collision can feature multiple groups.
vehicle_pairs <- vehicles_long %>%
  inner_join(accidents_main %>% filter(n_vehicles == 2) %>%
               select(incident_id, fatal_accident), by = "incident_id") %>%
  group_by(incident_id) %>% filter(n() == 2L) %>%
  summarise(vehicle_pair = paste(sort(vehicle_group), collapse = " + "),
            fatal_accident = first(fatal_accident), .groups = "drop")
vehicle_pair_stats <- vehicle_pairs %>% group_by(vehicle_pair) %>% summarise(
  incidents = n(), fatal_accidents = sum(fatal_accident),
  fatal_rate = round(100 * mean(fatal_accident), 2), .groups = "drop") %>%
  arrange(desc(incidents))

# 03.5 | People  ---------------------
driver_eda <- drivers_long %>%
  filter(outcome %in% c("Uninjured", "Injured", "Killed")) %>%
  mutate(age_class = factor(coalesce(age_group_clean, "Unknown"),
                            levels = age_levels),
         gender = coalesce(sex, "Unknown"))
driver_outcome_stats <- driver_eda %>% count(outcome, name = "drivers") %>%
  mutate(percentage = round(100 * drivers / sum(drivers), 2))
driver_age_stats <- driver_eda %>% group_by(age_class) %>% summarise(
  drivers = n(), uninjured = sum(outcome == "Uninjured"),
  injured = sum(outcome == "Injured"), killed = sum(outcome == "Killed"),
  death_pct = round(100 * killed / drivers, 2), .groups = "drop")
driver_gender_stats <- driver_eda %>% group_by(gender) %>% summarise(
  drivers = n(), injured = sum(outcome == "Injured"),
  killed = sum(outcome == "Killed"),
  death_pct = round(100 * killed / drivers, 2), .groups = "drop")

pedestrian_eda <- pedestrians_long %>%
  filter(outcome %in% c("Injured", "Killed")) %>%
  mutate(age_class = factor(coalesce(age_group_clean, "Unknown"),
                            levels = age_levels),
         gender = coalesce(sex, "Unknown"))
pedestrian_outcome_stats <- pedestrian_eda %>% count(outcome, name = "pedestrians") %>%
  mutate(percentage = round(100 * pedestrians / sum(pedestrians), 2))
pedestrian_age_stats <- pedestrian_eda %>% group_by(age_class) %>% summarise(
  casualties = n(), injured = sum(outcome == "Injured"),
  killed = sum(outcome == "Killed"),
  fatal_share = round(100 * killed / casualties, 2), .groups = "drop")
pedestrian_gender_stats <- pedestrian_eda %>% group_by(gender) %>% summarise(
  casualties = n(), injured = sum(outcome == "Injured"),
  killed = sum(outcome == "Killed"),
  fatal_share = round(100 * killed / casualties, 2), .groups = "drop")

casualty_role_stats <- tibble(
  role = c("Drivers", "Passengers", "Pedestrians", "Other vehicle occupants"),
  injured = c(sum(accidents_main$drivers_injured), sum(accidents_main$passengers_injured),
              sum(accidents_main$pedestrians_injured),
              sum(accidents_main$other_vehicle_people_injured)),
  killed = c(sum(accidents_main$drivers_killed), sum(accidents_main$passengers_killed),
             sum(accidents_main$pedestrians_killed),
             sum(accidents_main$other_vehicle_people_killed))
) %>% mutate(casualties = injured + killed,
             fatal_share = round(100 * killed / casualties, 2))
stopifnot(sum(casualty_role_stats$killed) == sum(accidents_main$total_deaths),
          sum(casualty_role_stats$injured) == sum(accidents_main$total_injured))

# 04 | GEOGRAPHICAL EDA -----------------------------------------------------
# All regional/provincial/municipal summaries retain every incident.
region_stats <- accidents_main %>% group_by(region_code, region_name) %>%
  summarise(incidents = n(), fatal_accidents = sum(fatal_accident),
            total_deaths = sum(total_deaths), fatal_rate = 100 * mean(fatal_accident),
            .groups = "drop") %>% arrange(desc(incidents))
province_stats <- accidents_main %>% group_by(provincia, province_name, region_name) %>%
  summarise(incidents = n(), fatal_accidents = sum(fatal_accident),
            total_deaths = sum(total_deaths), fatal_rate = 100 * mean(fatal_accident),
            .groups = "drop") %>% arrange(desc(incidents))
municipality_stats <- accidents_main %>%
  group_by(comune_istat, municipality_name, provincia, region_name) %>%
  summarise(incidents = n(), fatal_accidents = sum(fatal_accident),
            total_deaths = sum(total_deaths), fatal_rate = 100 * mean(fatal_accident),
            .groups = "drop") %>% arrange(desc(incidents))
stopifnot(nrow(region_stats) == 20L,
          sum(region_stats$incidents) == nrow(accidents_main),
          sum(province_stats$incidents) == nrow(accidents_main),
          sum(municipality_stats$incidents) == nrow(accidents_main))

regions_map <- regions_sf %>%
  mutate(region_code = sprintf("%02d", as.integer(COD_REG))) %>%
  left_join(region_stats, by = "region_code")
provinces_map <- provinces_sf %>%
  mutate(provincia = sprintf("%03d", as.integer(.data[[province_code_column]]))) %>%
  left_join(province_stats %>%
              select(provincia, incidents, fatal_accidents, fatal_rate),
            by = "provincia")
municipalities_map <- municipalities_sf %>%
  mutate(comune_istat = str_pad(as.character(.data[[municipality_code_column]]),
                                width = 6, pad = "0")) %>%
  left_join(municipality_stats %>%
              select(comune_istat, incidents, fatal_accidents, fatal_rate),
            by = "comune_istat")
unmapped_municipalities <- municipality_stats %>%
  filter(!comune_istat %in% municipalities_map$comune_istat)
milan_municipalities <- municipalities_map %>%
  filter(substr(comune_istat, 1, 3) == "015") %>%
  mutate(incidents = replace_na(incidents, 0L))

# 05 | CONDITIONAL + MULTIDIMENSIONAL DESCRIPTIVE EDA ----------------------

# 05.1 | Road location x time band (non-reported hour remains in main data) ----
road_time_stats <- accidents_main %>%
  filter(localizzazione_incidente_label %in% selected_roads,
         !is.na(time_band), time_band != "Not reported") %>%
  group_by(localizzazione_incidente_label, time_band) %>%
  summarise(incidents = n(), fatal_accidents = sum(fatal_accident),
            fatal_rate = round(100 * mean(fatal_accident), 2), .groups = "drop")

# 05.2 | Accident nature x number of vehicles -------------------------------
top_accident_types <- environment_stats[["Accident type"]] %>%
  arrange(desc(incidents)) %>% slice_head(n = 7) %>% pull(category)
accident_vehicle_stats <- accidents_main %>%
  mutate(vehicle_count_group = case_when(
    n_vehicles == 1 ~ "1", n_vehicles == 2 ~ "2", n_vehicles == 3 ~ "3",
    n_vehicles >= 4 ~ "4+"),
    vehicle_count_group = factor(vehicle_count_group,
                                 levels = c("1", "2", "3", "4+"))) %>%
  filter(natura_incidente_label %in% top_accident_types) %>%
  group_by(natura_incidente_label, vehicle_count_group) %>%
  summarise(incidents = n(), fatal_accidents = sum(fatal_accident),
            fatal_rate = round(100 * mean(fatal_accident), 2), .groups = "drop") %>%
  mutate(accident_type = factor(natura_incidente_label,
                                levels = rev(top_accident_types)))

# 05.3 | Driver age x vehicle (person-level, all outcomes documented) -------
driver_age_vehicle_stats <- driver_eda %>%
  filter(vehicle_group %in% selected_vehicle_groups,
         as.character(age_class) %in% c("15-17", "18-29", "30-44", "45-54",
                                      "55-64", "65+")) %>%
  group_by(vehicle_group, age_class) %>%
  summarise(drivers = n(), killed = sum(outcome == "Killed"),
            fatal_pct = round(100 * killed / drivers, 2), .groups = "drop")
driver_age_vehicle_plot <- driver_age_vehicle_stats %>%
  filter(as.character(age_class) %in% c("18-29", "30-44", "45-54", "55-64", "65+"),
         drivers >= 200) %>%
  mutate(age_class = factor(as.character(age_class),
                            levels = c("18-29", "30-44", "45-54", "55-64", "65+")))

# 05.4 | Road x accident nature ---------------------------------------------
road_nature_stats <- accidents_main %>%
  filter(localizzazione_incidente_label %in% selected_roads) %>%
  mutate(road = localizzazione_incidente_label,
         nature = if_else(natura_incidente_label %in% selected_natures,
                          natura_incidente_label, "Other/unspecified")) %>%
  group_by(road, nature) %>%
  summarise(incidents = n(), fatal_accidents = sum(fatal_accident),
            fatal_pct = round(100 * mean(fatal_accident), 2), .groups = "drop") %>%
  group_by(road) %>%
  mutate(share_within_road = round(100 * incidents / sum(incidents), 2)) %>%
  ungroup() %>% arrange(road, desc(incidents))
road_nature_plot <- road_nature_stats %>%
  filter(nature %in% c("Head-on collision", "Front-to-side collision",
                       "Run-off-road/skidding", "Pedestrian collision"),
         incidents >= 100) %>%
  mutate(road = factor(road, levels = selected_roads))

# 05.5 | Driver age x sex ----------------------------------------------------
selected_ages <- c("15-17", "18-29", "30-44", "45-54", "55-64", "65+")
driver_age_gender_stats <- drivers_long %>%
  filter(age_group_clean %in% selected_ages, sex %in% c("Male", "Female"),
         outcome %in% c("Uninjured", "Injured", "Killed")) %>%
  mutate(age_class = factor(age_group_clean, levels = selected_ages, ordered = TRUE)) %>%
  group_by(age_class, sex) %>%
  summarise(drivers = n(), killed = sum(outcome == "Killed"),
            fatal_pct = round(100 * killed / drivers, 2), .groups = "drop") %>%
  arrange(age_class, sex)

# 05.6 | Road x accident nature x night/other hours x outcome --------------
multivariate_stats <- accidents_main %>%
  filter(!is.na(hour), localizzazione_incidente_label %in% selected_roads,
         natura_incidente_label %in% c(
           "Head-on collision", "Front-to-side collision", "Rear-end collision",
           "Run-off-road/skidding", "Pedestrian collision")) %>%
  mutate(period = if_else(hour < 6, "night", "other_hours")) %>%
  group_by(localizzazione_incidente_label, natura_incidente_label, period) %>%
  summarise(incidents = n(), fatal_accidents = sum(fatal_accident),
            fatal_pct = 100 * mean(fatal_accident), .groups = "drop")
multivariate_comparison <- multivariate_stats %>%
  pivot_wider(names_from = period,
              values_from = c(incidents, fatal_accidents, fatal_pct)) %>%
  filter(incidents_night >= 100, incidents_other_hours >= 100) %>%
  mutate(difference_pp = round(fatal_pct_night - fatal_pct_other_hours, 2),
         fatal_pct_night = round(fatal_pct_night, 2),
         fatal_pct_other_hours = round(fatal_pct_other_hours, 2)) %>%
  arrange(desc(difference_pp))
multivariate_plot <- multivariate_comparison %>%
  arrange(difference_pp) %>%
  mutate(combination = paste(localizzazione_incidente_label,
                              natura_incidente_label, sep = " | "),
         combination = factor(combination, levels = combination))

# 06 | PRESUMED CIRCUMSTANCES -- EXPLORATORY, NOT FULLY VALIDATED ----------
circumstances_A_codes <- accidents_raw %>%
  count(code = veicolo__a____circostanze_10, sort = TRUE, name = "incidents")
circumstances_B_codes <- accidents_raw %>%
  count(code = veicolo__b____circostanze_13, sort = TRUE, name = "incidents")
circumstance_availability <- accidents_raw %>% summarise(
  A_reported = sum(!is.na(veicolo__a____circostanze_10)),
  B_reported = sum(!is.na(veicolo__b____circostanze_13)),
  B_vehicles = sum(!is.na(tipo_veicoli__b_)),
  B_without_vehicle_B = sum(is.na(tipo_veicoli__b_) &
                              !is.na(veicolo__b____circostanze_13)))
circumstances_B_without_vehicle <- accidents_raw %>%
  filter(is.na(tipo_veicoli__b_)) %>%
  count(code = veicolo__b____circostanze_13, sort = TRUE, name = "incidents")

circumstances_verified <- accidents_main %>%
  select(incident_id, natura_incidente, fatal_accident) %>%
  left_join(accidents_raw %>% transmute(
    incident_id = row_number(), code_A = veicolo__a____circostanze_10,
    code_B = veicolo__b____circostanze_13), by = "incident_id") %>%
  mutate(
    moving_vehicles = natura_incidente %in% c("01", "02", "03", "04"),
    speed_related = code_A %in% c("11", "12", "23", "24", "41", "42", "64", "65", "72") |
      (moving_vehicles & code_B %in% c("11", "12", "23", "24")),
    distracted_driving = code_A %in% c("02", "21", "61", "71") |
      (moving_vehicles & code_B %in% c("02", "21")),
    right_of_way_or_signals = code_A %in% c("04", "05", "06", "08", "46", "49", "50") |
      (moving_vehicles & code_B %in% c("04", "05", "06", "08")),
    following_too_closely = code_A %in% c("03", "22", "62") |
      (moving_vehicles & code_B %in% c("03", "22"))
  )
stopifnot(nrow(circumstances_verified) == nrow(accidents_main),
          !anyDuplicated(circumstances_verified$incident_id))
circumstances_long <- circumstances_verified %>%
  select(incident_id, speed_related, distracted_driving,
         right_of_way_or_signals, following_too_closely) %>%
  pivot_longer(-incident_id, names_to = "circumstance", values_to = "present") %>%
  filter(present) %>% distinct(incident_id, circumstance)
circumstances_corrected_stats <- circumstances_long %>%
  left_join(accidents_main %>% select(incident_id, fatal_accident),
            by = "incident_id") %>%
  group_by(circumstance) %>% summarise(
    incidents = n(), fatal_accidents = sum(fatal_accident),
    fatal_rate = round(100 * mean(fatal_accident), 2), .groups = "drop") %>%
  arrange(desc(incidents))
# Keep the original plotting object name for compatibility.
circumstance_stats <- circumstances_corrected_stats

cat("\nEDA TABLES READY | Geographic regions:", nrow(region_stats),
    "| Road x nature cells:", nrow(road_nature_stats),
    "| Multidimensional comparisons:", nrow(multivariate_comparison), "\n")


# 09 | FEATURE ENGINEERING & VALIDATED DESCRIPTIVE EXTENSIONS
# ============================================================================
# 09.1 | Multi-injury(MC) target and descriptive target-encoded rates -------------
accidents_main <- accidents_main %>%
  mutate(multi_casualty = total_injured >= 2L)

te_smooth <- function(df, group_col, target_col, smooth = 30) {
  global_mean <- mean(df[[target_col]], na.rm = TRUE)
  df %>%
    group_by(.data[[group_col]]) %>%
    summarise(te_sum = sum(.data[[target_col]], na.rm = TRUE),
              te_count = sum(!is.na(.data[[target_col]])),
              te_raw = mean(.data[[target_col]], na.rm = TRUE),
              .groups = "drop") %>%
    mutate(te_smoothed = (te_sum + global_mean * smooth) /
             (te_count + smooth), global_mean = global_mean)
}

te_natura <- te_smooth(accidents_main, "natura_incidente", "multi_casualty")
te_location <- te_smooth(accidents_main, "localizzazione_incidente", "multi_casualty")
te_road <- te_smooth(accidents_main, "tipo_di_strada", "multi_casualty")
te_vehicle <- te_smooth(accidents_main, "tipo_veicolo_a", "multi_casualty")

# EDA tables
te_natura_table <- te_natura %>%
  left_join(accidents_main %>% distinct(natura_incidente, natura_incidente_label),
            by = "natura_incidente") %>%
  arrange(desc(te_smoothed)) %>%
  transmute(code = natura_incidente, label = natura_incidente_label,
            n_accidents = te_count, raw_rate_pct = 100 * te_raw,
            encoded_value = te_smoothed)
te_vehicle_table <- te_vehicle %>%
  mutate(vehicle_label = label_istat(tipo_veicolo_a, vehicle_names, "vehicle A")) %>%
  arrange(desc(te_smoothed)) %>%
  transmute(code = tipo_veicolo_a, label = vehicle_label,
            n_accidents = te_count, raw_rate_pct = 100 * te_raw,
            encoded_value = te_smoothed) %>%
  filter(n_accidents >= 50)

national_mc_rate <- 100 * mean(accidents_main$multi_casualty, na.rm = TRUE)
national_fatal_rate <- 100 * mean(accidents_main$fatal_accident, na.rm = TRUE)

# 10 | TEMPORAL FEATURES -----------------------------------------------------
stopifnot(all(is.na(accidents_main$hour) | accidents_main$hour %in% 0:23))
accidents_main <- accidents_main %>%
  mutate(is_night = case_when(
           is.na(hour) ~ NA,
           TRUE ~ hour >= 22L | hour <= 5L),
         is_weekend = weekend)


night_table <- accidents_main %>%
  filter(!is.na(is_night)) %>%
  mutate(group = if_else(is_night,
                         "Night (22-05)",
                         "Other hours (06-21)")) %>%
  group_by(group) %>%
  summarise(
    accidents = n(),
    multi_injury_accidents = sum(multi_casualty, na.rm = TRUE),
    mc_rate = 100 * mean(multi_casualty, na.rm = TRUE),
    fatal_rate = 100 * mean(fatal_accident, na.rm = TRUE),
    .groups = "drop"
  )

weekend_table <- accidents_main %>%
  mutate(group = if_else(is_weekend,
                         "Weekend (Sat-Sun)",
                         "Weekday (Mon-Fri)")) %>%
  group_by(group) %>%
  summarise(
    accidents = n(),
    multi_injury_accidents = sum(multi_casualty, na.rm = TRUE),
    mc_rate = 100 * mean(multi_casualty, na.rm = TRUE),
    fatal_rate = 100 * mean(fatal_accident, na.rm = TRUE),
    .groups = "drop"
  )

unknown_hour_table <- accidents_main %>%
  group_by(hour_reported = !is.na(hour)) %>%
  summarise(accidents = n(), mc_rate = 100 * mean(multi_casualty),
            fatal_rate = 100 * mean(fatal_accident), .groups = "drop")

hourly_mc <- accidents_main %>%
  filter(!is.na(hour)) %>% group_by(hour) %>%
  summarise(accidents = n(), mc_rate = 100 * mean(multi_casualty),
            .groups = "drop") %>%
  mutate(is_night_hour = hour >= 22L | hour <= 5L)
temporal_comparison_data <- bind_rows(
    night_table %>% mutate(comparison = "Night vs other hours"),
    weekend_table %>% mutate(comparison = "Weekend vs weekday")) %>%
  pivot_longer(c(mc_rate, fatal_rate), names_to = "metric", values_to = "rate") %>%
  mutate(metric = recode(metric, "mc_rate" = "Multi-injury accident rate",
                         "fatal_rate" = "Fatal accident rate"),
         highlighted = grepl("^Night|^Weekend", group))
stopifnot(all(is.na(temporal_comparison_data$rate) |
              (temporal_comparison_data$rate >= 0 & temporal_comparison_data$rate <= 100)))
# 11 | VALIDATED BINARY FLAGS ------------------------------------------------
# Preserve unknown values as NA instead of declaring them FALSE.
flag_from_code <- function(x, codes) {
  dplyr::if_else(is.na(x), NA, x %in% codes)
}

accidents_main <- accidents_main %>%
  mutate(
    is_motorway = flag_from_code(localizzazione_incidente, "7"),
    is_open = flag_from_code(intersezione_o_non_interse3, "07"),
    at_junc = flag_from_code(intersezione_o_non_interse3,
                            c("01", "02", "03", "04", "05", "06")),
    bad_weather = flag_from_code(condizioni_meteorologiche, c("3", "4", "5", "6")),
    bad_surface = flag_from_code(fondo_stradale, c("2", "3", "4", "5")),
    is_headon = flag_from_code(natura_incidente, "01"),
    is_rearend = flag_from_code(natura_incidente, "04"),
    is_ped_hit = flag_from_code(natura_incidente, "05"),
    is_runoff = flag_from_code(natura_incidente, "10"),
    is_obstacle_hit = flag_from_code(natura_incidente, "08"),
    is_moto_a = flag_from_code(tipo_veicolo_a, c("15", "16", "17")),
    is_vru_a = flag_from_code(tipo_veicolo_a, c("14", "15", "16", "17", "22", "23")),
    cil_a_high = case_when(
      is.na(Cilindrata_A_Classi) ~ NA,
      suppressWarnings(as.integer(Cilindrata_A_Classi)) %in% 3:17 ~ TRUE,
      suppressWarnings(as.integer(Cilindrata_A_Classi)) %in% 1:2 ~ FALSE,
      TRUE ~ NA),
    has_veh_b = !is.na(tipo_veicoli__b_),
    is_single = if_else(is.na(n_vehicles), NA, n_vehicles == 1L),
    vb_is_heavy = flag_from_code(tipo_veicoli__b_,
                                 c("08", "09", "10", "12", "18")),
    male_a = flag_from_code(veicolo__a___sesso_conducente, "1"),
    young_a = flag_from_code(veicolo__a___et__conducente, c("15-17", "18-29")),
    senior_a = flag_from_code(veicolo__a___et__conducente, "65+"),
    no_licence_a = flag_from_code(veicolo__a___patente_condu29, "9")
  )

step3_flags <- c("is_motorway", "is_open", "at_junc", "bad_weather",
                 "bad_surface", "is_headon", "is_rearend", "is_ped_hit",
                 "is_runoff", "is_obstacle_hit", "is_moto_a", "is_vru_a",
                 "cil_a_high", "has_veh_b", "is_single", "vb_is_heavy",
                 "male_a", "young_a", "senior_a", "no_licence_a")
step3_labels <- c("Motorway", "Straight road", "Junction / level crossing", "Adverse weather",
                  "Non-dry road surface", "Head-on", "Rear-end", "Pedestrian collision",
                  "Run-off/skidding", "Collision with obstacle", "Motorcycle (A)",
                  "Vulnerable road user (A)", "Engine >250cc (A)",
                  "Vehicle B recorded", "Exactly one vehicle", "Heavy vehicle B",
                  "Male driver A", "Driver A aged 15-29", "Driver A aged 65+",
                  "Driver A without licence")
step3_groups <- c(rep("Road", 3), rep("Conditions", 2),
                  rep("Accident nature", 5), rep("Vehicle A", 3),
                  rep("Vehicle configuration", 3), rep("Driver A", 4))
stopifnot(length(step3_flags) == length(step3_labels),
          length(step3_flags) == length(step3_groups))
step3_summary <- bind_rows(lapply(seq_along(step3_flags), function(i) {
  x <- accidents_main[[step3_flags[i]]]
  tibble(feature = step3_flags[i], label = step3_labels[i],
         group = step3_groups[i],
         n_known = sum(!is.na(x)),
         n_yes = sum(x %in% TRUE), n_no = sum(x %in% FALSE),
         pct_yes = 100 * mean(x, na.rm = TRUE),
         mc_yes = 100 * mean(accidents_main$multi_casualty[x %in% TRUE]),
         mc_no = 100 * mean(accidents_main$multi_casualty[x %in% FALSE]),
         fatal_yes = 100 * mean(accidents_main$fatal_accident[x %in% TRUE]),
         fatal_no = 100 * mean(accidents_main$fatal_accident[x %in% FALSE]))
})) %>%
  mutate(delta_mc = mc_yes - mc_no,
         delta_fatal = fatal_yes - fatal_no)

nature_flags_long <- step3_summary %>%
  filter(group == "Accident nature") %>%
  select(label, delta_mc, delta_fatal) %>%
  pivot_longer(c(delta_mc, delta_fatal), names_to = "metric", values_to = "delta") %>%
  mutate(metric = recode(metric,
    "delta_mc" = "Multi-injury difference (pp)",
    "delta_fatal" = "Fatal accident difference (pp)"))
road_flags_data <- accidents_main %>%
  mutate(road_condition = case_when(
    is_motorway %in% TRUE ~ "Motorway",
    is_open %in% TRUE ~ "Straight road (not motorway)",
    at_junc %in% TRUE ~ "Junction (neither above)",
    TRUE ~ "Other/unknown")) %>%
  group_by(road_condition) %>%
  summarise(incidents = n(), mc_rate = 100 * mean(multi_casualty),
            fatal_rate = 100 * mean(fatal_accident), .groups = "drop")
# The descriptive grouping above is PRIORITIZED and mutually exclusive;
# it does not represent mutually exclusive original ISTAT road variables.
# 12 | VEHICLE-EXPOSURE PROXIES VS OBSERVED CASUALTIES -------------------------
# ISTAT records injured/deceased passengers, NOT total occupants.
accidents_main <- accidents_main %>%
  mutate(vehicle_count_group = case_when(
    is.na(n_vehicles) ~ NA_character_,
    n_vehicles >= 4L ~ "4+", TRUE ~ as.character(n_vehicles)),
    vehicle_count_group = factor(vehicle_count_group, levels = c("1","2","3","4+")))

vehicle_exposure_table <- accidents_main %>%
  filter(!is.na(vehicle_count_group)) %>%
  group_by(vehicle_count_group) %>%
  summarise(accidents = n(), multi_injury_accidents = sum(multi_casualty),
            mc_rate = 100 * mean(multi_casualty),
            fatal_rate = 100 * mean(fatal_accident), .groups = "drop")
# Outcome variables used only for describing the composition of casualties.
casualty_profile_table <- casualty_role_stats %>%
  mutate(role = factor(role, levels = role))
# 13 | COMPUTED FEATURE OVERVIEW (NO HARD-CODED 'IMPORTANCE') -----------------
feature_overview <- step3_summary %>%
  filter(is.finite(delta_mc), n_yes >= MIN_N_PLOT, n_no >= MIN_N_PLOT) %>%
  arrange(delta_mc) %>%
  mutate(label = factor(label, levels = unique(label)))
# Safe baseline for future modelling: no outcome-derived passenger/person fields;
# no globally target-encoded variables; no fatalities or injuries as predictors.
accidents_ml_safe <- accidents_main %>%
  select(incident_id, multi_casualty,
         giorno, trimestre, hour, localizzazione_incidente, tipo_di_strada,
         intersezione_o_non_interse3, pavimentazione, fondo_stradale,
         segnaletica, condizioni_meteorologiche, natura_incidente,
         tipo_veicolo_a, tipo_veicoli__b_, tipo_veicolo__c_,
         n_vehicles, is_night, is_weekend,
         is_motorway, is_open, at_junc, bad_weather, bad_surface,
         is_headon, is_rearend, is_ped_hit, is_runoff, is_obstacle_hit,
         is_moto_a, is_vru_a, cil_a_high,
         has_veh_b, is_single, vb_is_heavy,
         male_a, young_a, senior_a, no_licence_a,
         has_car, has_motorcycle, has_bicycle, has_goods_vehicle)
stopifnot(nrow(accidents_ml_safe) == nrow(accidents_main),
          !any(c("passengers_injured", "passengers_killed", "total_injured",
                 "fatal_accident", "natura_te", "vehicle_a_te",
                 "pedestrians_injured", "total_potential_victims") %in%
               names(accidents_ml_safe)))


# 14 | ADDITIONAL MULTIDIMENSIONAL PLOTS --------------------------------------
# All summaries are recomputed from existing incident/person-level tables.
# Missing or low-count combinations are excluded explicitly; they do NOT imply 0.

day_time_mc <- accidents_main %>%
  filter(!is.na(hour), !is.na(time_band), time_band != "Not reported") %>%
  group_by(day_label, time_band) %>%
  summarise(incidents = n(), fatal_rate = 100 * mean(fatal_accident),
            mc_rate = 100 * mean(multi_casualty), .groups = "drop")

quarter_day_stats <- accidents_main %>%
  group_by(quarter_label, day_label) %>%
  summarise(incidents = n(), fatal_rate = 100 * mean(fatal_accident),
            mc_rate = 100 * mean(multi_casualty), .groups = "drop")

road_nature_hm <- accidents_main %>%
  filter(localizzazione_incidente_label %in% selected_roads,
         natura_incidente_label %in% selected_natures) %>%
  group_by(road = localizzazione_incidente_label,
           nature = natura_incidente_label) %>%
  summarise(incidents = n(), fatal_rate = 100 * mean(fatal_accident),
            mc_rate = 100 * mean(multi_casualty), .groups = "drop") %>%
  filter(incidents >= MIN_N_PLOT)

weather_surface_mc <- accidents_main %>%
  filter(!is.na(condizioni_meteorologiche_label),
         !is.na(fondo_stradale_label)) %>%
  group_by(weather = condizioni_meteorologiche_label,
           surface = fondo_stradale_label) %>%
  summarise(incidents = n(), fatal_rate = 100 * mean(fatal_accident),
            mc_rate = 100 * mean(multi_casualty), .groups = "drop") %>%
  filter(incidents >= MIN_N_PLOT)

# Use already-aggregated driver-age-by-vehicle statistics from primary EDA.
driver_age_vehicle_hm <- driver_age_vehicle_stats %>%
  filter(vehicle_group %in% selected_vehicle_groups,
         as.character(age_class) %in% c("18-29","30-44","45-54","55-64","65+"),
         drivers >= MIN_N_PLOT)

# Exactly-two-vehicle accidents only; no comparison against all-vehicle totals.
vehicle_pair_mc <- vehicle_pairs %>%
  left_join(accidents_main %>%
              select(incident_id, multi_casualty), by = "incident_id") %>%
  group_by(vehicle_pair) %>%
  summarise(incidents = n(), fatal_rate = 100 * mean(fatal_accident),
            mc_rate = 100 * mean(multi_casualty), .groups = "drop") %>%
  filter(incidents >= 150)

sm_data <- driver_age_vehicle_hm %>%
  mutate(age_class = factor(as.character(age_class),
    levels = c("18-29", "30-44", "45-54", "55-64", "65+")))

# Geography: avoid repeated 'incidents' columns when joining to provinces_map.
province_mc <- accidents_main %>%
  group_by(provincia) %>%
  summarise(mc_rate = 100 * mean(multi_casualty), .groups = "drop")
provinces_map_full <- provinces_map %>%
  left_join(province_mc, by = "provincia")

# Retain comparison table even though the discarded dumbbell plot is omitted.
if (nrow(multivariate_comparison) > 0) {
  dumbbell_data <- multivariate_comparison %>%
    arrange(difference_pp) %>%
    mutate(combination = paste(localizzazione_incidente_label,
                               natura_incidente_label, sep = " | "),
           combination = factor(combination, levels = combination))
} else {
  dumbbell_data <- multivariate_comparison
}


# ==========================================================================
# 14.1 | OPTIONAL DESCRIPTIVE TABLE CATALOGUE (NO ADDITIONAL PLOTS)
# ==========================================================================
# Use e.g. print(master_eda_tables$road_nature_stats).
# For geographic sf tables (regions_map, provinces_map, municipalities_map),
master_table_names <- c(
  "variable_overview", "quality_summary", "preparation_audit",
  "driver_completeness", "pedestrian_sex_summary", "license_age_crosscheck",
  "death_distribution", "injury_distribution", "vehicle_count_stats",
  "day_stats", "hour_stats", "quarter_stats", "time_band_stats",
  "day_time_stats", "weekday_normalized", "quarter_normalized",
  "hour_missing_check", "weather_surface_stats", "vehicle_frequency",
  "vehicle_severity", "vehicle_pair_stats", "driver_outcome_stats",
  "driver_age_stats", "driver_gender_stats", "pedestrian_outcome_stats",
  "pedestrian_age_stats", "pedestrian_gender_stats", "casualty_role_stats",
  "region_stats", "province_stats", "municipality_stats",
  "unmapped_municipalities", "road_time_stats", "accident_vehicle_stats",
  "driver_age_vehicle_stats", "road_nature_stats", "driver_age_gender_stats",
  "multivariate_stats", "multivariate_comparison", "circumstances_A_codes",
  "circumstances_B_codes", "circumstance_availability",
  "circumstances_B_without_vehicle", "circumstances_corrected_stats",
  "te_natura_table", "te_vehicle_table", "night_table", "weekend_table",
  "unknown_hour_table", "hourly_mc", "temporal_comparison_data",
  "step3_summary", "nature_flags_long", "road_flags_data",
  "vehicle_exposure_table", "casualty_profile_table", "feature_overview",
  "day_time_mc", "quarter_day_stats", "road_nature_hm",
  "weather_surface_mc", "driver_age_vehicle_hm", "vehicle_pair_mc",
  "dumbbell_data", "sm_data", "province_mc"
)
stopifnot(all(vapply(master_table_names, exists, logical(1))))
master_eda_tables <- mget(master_table_names, inherits = TRUE)
cat("\nMaster EDA retained tabular objects:", length(master_eda_tables), "\n")

# ============================================================================
# 15 | PRESENTATION EDA: FINAL SELECTED CHARTS — ITALY 2024
# ============================================================================
required_objects <- c("accidents_main", "weekday_normalized", "regions_sf",
                      "drivers_long", "weekday_levels")
missing_objects <- required_objects[!vapply(required_objects, exists, logical(1))]
if (length(missing_objects)) {
  stop("Dataset preparation did not create: ",
       paste(missing_objects, collapse = ", "))
}
stopifnot("multi_casualty" %in% names(accidents_main))

ppt_eda <- list()

# Presentation palette
BLUE   <- "#2865A3"
DARK   <- "#173E72"
GREY   <- "#929292"
LIGHT  <- "#CBD5E1"
ORANGE <- "#ED8619"
RED    <- "#C83D3D"
GREEN  <- "#34863C"

MC_AVG <- 100 * mean(accidents_main$multi_casualty, na.rm = TRUE)

ppt_theme <- theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 11, colour = "grey35"),
    plot.caption = element_text(size = 9, colour = "grey45"),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "bottom"
  )

# ================================================================
# PLOT 1 — DISTRIBUTION OF INJURED PEOPLE
# ================================================================

injury_order <- c(
  "1 injured",
  "2 injured",
  "3 injured",
  "4+ injured",
  "Deaths only"
)

injury_colours <- c(
  "1 injured" = GREY,
  "2 injured" = ORANGE,
  "3 injured" = RED,
  "4+ injured" = "#B72C42",
  "Deaths only" = "#B8BEC7"
)

injury_ppt_data <- accidents_main %>%
  mutate(
    outcome_group = case_when(
      total_injured == 0 & total_deaths > 0 ~ "Deaths only",
      total_injured == 1 ~ "1 injured",
      total_injured == 2 ~ "2 injured",
      total_injured == 3 ~ "3 injured",
      total_injured >= 4 ~ "4+ injured",
      TRUE ~ NA_character_
    )
  )

stopifnot(!anyNA(injury_ppt_data$outcome_group))

injury_ppt_data <- injury_ppt_data %>%
  count(outcome_group, name = "accidents") %>%
  mutate(
    pct = 100 * accidents / sum(accidents),
    outcome_group = factor(
      outcome_group,
      levels = rev(injury_order)
    )
  )

ppt_eda$injury_distribution <- ggplot(
  injury_ppt_data,
  aes(
    x = outcome_group,
    y = pct,
    fill = outcome_group
  )
) +
  geom_col(width = 0.68, show.legend = FALSE) +
  geom_text(
    aes(label = sprintf("%.1f%%", pct)),
    hjust = -0.12,
    size = 4.3,
    fontface = "bold"
  ) +
  coord_flip(clip = "off") +
  scale_fill_manual(values = injury_colours) +
  scale_y_continuous(
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0, 0.18))
  ) +
  labs(
    title = "Accidents by Number of Injured People",
    subtitle = "Distribution of recorded accidents | Italy 2024",
    x = NULL,
    y = "Share of accidents (%)",
    caption = paste(
      "Deaths only = zero injured and at least one death.",
      "Fatal accidents with injured people remain in their injury category."
    )
  ) +
  ppt_theme

# ================================================================
# PLOT 2 — WEEKDAY: FREQUENCY + MULTI-INJURY RATE
# ================================================================

weekday_mc_data <- accidents_main %>%
  group_by(day_label) %>%
  summarise(
    mc_rate = 100 * mean(multi_casualty, na.rm = TRUE),
    .groups = "drop"
  )

weekday_ppt_data <- weekday_normalized %>%
  select(day_label, accidents_per_day) %>%
  left_join(weekday_mc_data, by = "day_label")

weekday_factor <- max(
  weekday_ppt_data$accidents_per_day
) / max(weekday_ppt_data$mc_rate)

ppt_eda$weekday <- ggplot(
  weekday_ppt_data,
  aes(x = day_label)
) +
  geom_col(
    aes(y = accidents_per_day),
    fill = LIGHT,
    width = 0.70
  ) +
  geom_line(
    aes(y = mc_rate * weekday_factor, group = 1),
    colour = BLUE,
    linewidth = 1.25
  ) +
  geom_point(
    aes(y = mc_rate * weekday_factor),
    colour = BLUE,
    size = 3
  ) +
  geom_text(
    aes(
      y = mc_rate * weekday_factor,
      label = sprintf("%.1f%%", mc_rate)
    ),
    vjust = -0.85,
    colour = BLUE,
    fontface = "bold",
    size = 3.2
  ) +
  scale_y_continuous(
    name = "Average accidents per weekday",
    sec.axis = sec_axis(
      ~ . / weekday_factor,
      name = "Multi-injury accident rate (%)"
    ),
    expand = expansion(mult = c(0, 0.18))
  ) +
  labs(
    title = "Road Accidents by Weekday",
    subtitle = "Frequency and multi-injury rate | Italy 2024",
    x = NULL,
    caption = paste(
      "Bars: mean accidents per calendar weekday.",
      "Blue line: share of recorded accidents with 2+ injured.",
      "Left and right axes have different units."
    )
  ) +
  ppt_theme +
  theme(panel.grid.major.x = element_blank())

# ================================================================
# PLOT 3 — HEATMAP: WEEKDAY x TIME BAND
# Four time bands, NOT 24 separate hours
# ================================================================

time_band_order <- c(
  "00–05", "06–11", "12–17", "18–23"
)

heatmap_ppt_data <- accidents_main %>%
  filter(
    !is.na(time_band),
    as.character(time_band) %in% time_band_order
  ) %>%
  group_by(day_label, time_band) %>%
  summarise(
    accidents = n(),
    mc_rate = 100 * mean(multi_casualty),
    .groups = "drop"
  ) %>%
  mutate(
    weekday_plot = factor(
      as.character(day_label),
      levels = rev(weekday_levels)
    ),
    band_plot = factor(
      as.character(time_band),
      levels = time_band_order
    )
  )

heatmap_mid <- mean(range(
  heatmap_ppt_data$mc_rate,
  na.rm = TRUE
))

ppt_eda$time_heatmap <- ggplot(
  heatmap_ppt_data,
  aes(
    x = band_plot,
    y = weekday_plot,
    fill = mc_rate
  )
) +
  geom_tile(
    colour = "white",
    linewidth = 1.2
  ) +
  geom_text(
    aes(
      label = sprintf("%.1f%%", mc_rate),
      colour = mc_rate > heatmap_mid
    ),
    size = 4,
    fontface = "bold"
  ) +
  scale_colour_manual(
    values = c(
      "FALSE" = "#1E293B",
      "TRUE" = "white"
    ),
    guide = "none"
  ) +
  scale_fill_gradient(
    low = "#E7F0FA",
    high = DARK,
    name = "2+ injured (%)"
  ) +
  labs(
    title = "Multi-Injury Rate by Weekday and Time Band",
    subtitle = "Four six-hour periods | Italy 2024",
    x = "Time of day",
    y = NULL,
    caption = "Rates conditional on recorded accidents; unknown hours excluded."
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    panel.grid = element_blank(),
    legend.position = "right"
  )
# ================================================================
# PLOT 4 — ACCIDENT NATURE (FINAL: RATE-ORDERED WITH BENCHMARK)
# ================================================================

nature_ppt_data <- accidents_main %>%
  mutate(
    category = case_when(
      natura_incidente_label == "Head-on collision" ~ "Head-on",
      natura_incidente_label == "Front-to-side collision" ~ "Front-lateral",
      natura_incidente_label == "Rear-end collision" ~ "Rear-end",
      natura_incidente_label == "Pedestrian collision" ~ "Pedestrian hit",
      natura_incidente_label == "Run-off-road/skidding" ~ "Run-off road",
      natura_incidente_label == "Collision with obstacle" ~ "Hit obstacle",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(category)) %>%
  group_by(category) %>%
  summarise(
    accidents = n(),
    mc_rate = 100 * mean(multi_casualty),
    .groups = "drop"
  ) %>%
  mutate(
    highlight = category %in% c("Head-on", "Front-lateral", "Rear-end"),
    category = reorder(category, mc_rate)
  )

ppt_eda$accident_nature <- ggplot(
  nature_ppt_data,
  aes(
    x = category,
    y = mc_rate
  )
) +
  geom_col(
    aes(fill = highlight),
    width = 0.67,
    show.legend = FALSE
  ) +
  geom_hline(
    yintercept = MC_AVG,
    colour = "#D94841",
    linetype = "dashed",
    linewidth = 1
  ) +
  geom_text(
    aes(label = sprintf("%.1f%%", mc_rate)),
    hjust = -0.12,
    size = 4,
    fontface = "bold"
  ) +
  coord_flip(clip = "off") +
  scale_fill_manual(
    values = c(
      "TRUE" = BLUE,
      "FALSE" = GREY
    )
  ) +
  scale_y_continuous(
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0, 0.20))
  ) +
  labs(
    title = "Multi-Injury Rate by Accident Nature",
    subtitle = "Collision types ordered by multi-injury rate",
    x = NULL,
    y = "Accidents with 2+ injured (%)",
    caption = paste0(
      "Red dashed line: national average (",
      round(MC_AVG, 1), "%)."
    )
  ) +
  ppt_theme +
  theme(panel.grid.major.y = element_blank())


# ================================================================
# PLOT 5 — DRIVER A AGE
# ================================================================

age_order <- c(
  "15-17", "18-29", "30-44",
  "45-54", "55-64", "65+"
)

age_colours <- c(
  "15-17" = GREY,
  "18-29" = RED,
  "30-44" = BLUE,
  "45-54" = BLUE,
  "55-64" = GREEN,
  "65+" = GREEN
)

age_ppt_data <- drivers_long %>%
  filter(
    vehicle_slot == "A",
    age_group_clean %in% age_order
  ) %>%
  distinct(incident_id, .keep_all = TRUE) %>%
  select(incident_id, age_group_clean) %>%
  inner_join(
    accidents_main %>%
      select(incident_id, multi_casualty),
    by = "incident_id"
  ) %>%
  group_by(age_group_clean) %>%
  summarise(
    accidents = n(),
    mc_rate = 100 * mean(multi_casualty),
    .groups = "drop"
  ) %>%
  mutate(
    age_group_clean = factor(
      age_group_clean,
      levels = age_order
    )
  )

ppt_eda$driver_age <- ggplot(
  age_ppt_data,
  aes(
    x = age_group_clean,
    y = mc_rate,
    fill = age_group_clean
  )
) +
  geom_col(width = 0.70, show.legend = FALSE) +
  geom_hline(
    yintercept = MC_AVG,
    colour = ORANGE,
    linetype = "dashed",
    linewidth = 0.9
  ) +
  geom_text(
    aes(label = sprintf("%.1f%%", mc_rate)),
    vjust = -0.55,
    size = 4,
    fontface = "bold"
  ) +
  scale_fill_manual(values = age_colours) +
  scale_y_continuous(
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0, 0.20))
  ) +
  labs(
    title = "Multi-Injury Rate by Driver Age",
    subtitle = "Age of driver A | Italy 2024",
    x = "Driver age group",
    y = "Accidents with 2+ injured (%)",
    caption = paste(
      "Dashed line: overall accident rate.",
      "Accidents without a recognised driver A age excluded."
    )
  ) +
  ppt_theme +
  theme(panel.grid.major.x = element_blank())

# ================================================================
# PLOT 6 — ITALY REGIONAL CHOROPLETH
# 20 administrative regions
# ================================================================

region_mc_ppt <- accidents_main %>%
  group_by(region_code, region_name) %>%
  summarise(
    accidents = n(),
    multi_injury = sum(multi_casualty, na.rm = TRUE),
    mc_rate = 100 * mean(multi_casualty, na.rm = TRUE),
    .groups = "drop"
  )

regions_ppt_map <- regions_sf %>%
  mutate(
    region_code = sprintf("%02d", as.integer(COD_REG))
  ) %>%
  left_join(
    region_mc_ppt,
    by = "region_code"
  )

stopifnot(
  nrow(regions_ppt_map) == 20L,
  !anyNA(regions_ppt_map$mc_rate)
)

ppt_eda$region_map <- ggplot(regions_ppt_map) +
  geom_sf(
    aes(fill = mc_rate),
    colour = "white",
    linewidth = 0.35
  ) +
  scale_fill_gradient(
    low = "#E3EFF9",
    high = DARK,
    name = "2+ injured (%)"
  ) +
  coord_sf(datum = NA) +
  labs(
    title = "Multi-Injury Accident Rate by Region",
    subtitle = "Italy 2024 | 20 regions",
    caption = paste(
      "Rate = accidents with 2+ injured / recorded accidents",
      "within the region. Not population-adjusted."
    )
  ) +
  theme_void(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 11),
    legend.position = "right"
  )

# Companion chart: top five regions
region_top5_ppt <- region_mc_ppt %>%
  slice_max(
    order_by = mc_rate,
    n = 5,
    with_ties = FALSE
  )

ppt_eda$region_top5 <- ggplot(
  region_top5_ppt,
  aes(
    x = reorder(region_name, mc_rate),
    y = mc_rate
  )
) +
  geom_col(fill = BLUE, width = 0.68) +
  geom_text(
    aes(label = sprintf("%.1f%%", mc_rate)),
    hjust = -0.12,
    size = 3.6,
    fontface = "bold"
  ) +
  coord_flip(clip = "off") +
  scale_y_continuous(
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0, 0.25))
  ) +
  labs(
    title = "Top 5 Regions",
    subtitle = "Highest multi-injury accident rate",
    x = NULL,
    y = "Multi-injury rate (%)"
  ) +
  ppt_theme

# Map and ranking together, for one PowerPoint slide
ppt_eda$geography_slide <-
  ppt_eda$region_map +
  ppt_eda$region_top5 +
  patchwork::plot_layout(widths = c(1.6, 1))

## ================================================================
# PLOT 7 — ROAD CONTEXT
# Blue = above national average | Grey = below national average
# ================================================================

MC_AVG <- 100 * mean(
  accidents_main$multi_casualty,
  na.rm = TRUE
)

road_context_data <- accidents_main %>%
  mutate(
    category = case_when(
      localizzazione_incidente %in% c("0", "1", "2", "3") ~
        "Built-up roads",
      localizzazione_incidente == "7" ~
        "Motorway",
      localizzazione_incidente %in% c("5", "6", "9") ~
        "Other extra-urban main roads",
      localizzazione_incidente == "4" ~
        "Extra-urban municipal roads",
      TRUE ~ "Other / unknown"
    )
  ) %>%
  group_by(category) %>%
  summarise(
    accidents = n(),
    mc_rate = 100 * mean(multi_casualty, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    category = reorder(category, mc_rate),
    above_avg = mc_rate >= MC_AVG
  )

ppt_eda$road_context <- ggplot(
  road_context_data,
  aes(
    x = category,
    y = mc_rate,
    fill = above_avg
  )
) +
  geom_col(width = 0.68) +
  geom_hline(
    yintercept = MC_AVG,
    colour = "#D94841",
    linetype = "dashed",
    linewidth = 1
  ) +
  geom_text(
    aes(label = sprintf("%.1f%%", mc_rate)),
    hjust = -0.15,
    size = 4,
    fontface = "bold"
  ) +
  coord_flip(clip = "off") +
  scale_fill_manual(
    values = c(
      "TRUE" = BLUE,
      "FALSE" = GREY
    ),
    guide = "none"
  ) +
  scale_y_continuous(
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0, 0.18))
  ) +
  labs(
    title = "Multi-Injury Rate by Road Context",
    subtitle = "Location-based road groups | Italy 2024",
    x = NULL,
    y = "Accidents with 2+ injured (%)",
    caption = paste0(
      "Red dashed line = national multi-injury rate (",
      round(MC_AVG, 1), "%)."
    )
  ) +
  ppt_theme +
  theme(
    panel.grid.major.y = element_blank()
  )
# ================================================================
# DESCRIPTIVE TABLES: KEEP DATA, OMIT UNUSED PRESENTATION PLOTS
# ================================================================

weather_ppt_data <- accidents_main %>%
  group_by(condizioni_meteorologiche_label) %>%
  summarise(
    accidents = n(),
    multi_injury_accidents = sum(multi_casualty, na.rm = TRUE),
    mc_rate = 100 * mean(multi_casualty, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(accidents))

surface_ppt_data <- accidents_main %>%
  mutate(
    category = case_when(
      fondo_stradale_label == "Dry" ~ "Dry",
      fondo_stradale_label == "Wet" ~ "Wet",
      fondo_stradale_label %in% c(
        "Slippery", "Icy", "Snow-covered"
      ) ~ "Slippery / icy / snowy",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(category)) %>%
  group_by(category) %>%
  summarise(
    accidents = n(),
    mc_rate = 100 * mean(multi_casualty),
    .groups = "drop"
  )


# Vehicle-count table
vehicles_ppt_data <- accidents_main %>%
  filter(!is.na(vehicle_count_group)) %>%
  group_by(vehicle_count_group) %>%
  summarise(
    accidents = n(),
    multi_injury_accidents = sum(multi_casualty, na.rm = TRUE),
    mc_rate = 100 * mean(multi_casualty, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    accident_share = 100 * accidents / sum(accidents)
  )



# Fatality summary table
fatal_ppt_table <- accidents_main %>%
  summarise(
    recorded_accidents = n(),
    fatal_accidents = sum(fatal_accident, na.rm = TRUE),
    fatal_accident_rate = 100 * mean(
      fatal_accident, na.rm = TRUE
    ),
    people_killed = sum(total_deaths, na.rm = TRUE),
    deaths_only_accidents = sum(
      total_injured == 0 & total_deaths > 0,
      na.rm = TRUE
    ),
    fatal_and_multi_injury = sum(
      fatal_accident & multi_casualty,
      na.rm = TRUE
    )
  )

eda_tables <- list(
  injury_distribution = injury_ppt_data,
  weekday = weekday_ppt_data,
  time_heatmap = heatmap_ppt_data,
  accident_nature = nature_ppt_data,
  driver_age = age_ppt_data,
  regional_rates = region_mc_ppt,
  road_context = road_context_data,
  weather = weather_ppt_data,
  road_surface = surface_ppt_data,
  vehicle_count = vehicles_ppt_data,
  fatality_summary = fatal_ppt_table
)

cat("\n--- PRESENTATION EDA TABLES ---\n")
for (nm in names(eda_tables)) {
  cat("\n[", nm, "]\n", sep = "")
  print(eda_tables[[nm]], width = Inf)
}

cat("\nVehicle-count coverage: ", sum(vehicles_ppt_data$accidents),
    " / ", nrow(accidents_main), " incidents\n", sep = "")
cat("Vehicle-count shares: ",
    round(sum(vehicles_ppt_data$accident_share), 2), "%\n", sep = "")
cat("National multi-injury rate: ", round(MC_AVG, 2), "%\n", sep = "")

# =================
# OUTPUT
# =================

selected_plot_names <- c(
  "injury_distribution",
  "weekday",
  "time_heatmap",
  "accident_nature",
  "driver_age",
  "region_map",
  "region_top5",
  "road_context"
)
stopifnot(all(selected_plot_names %in% names(ppt_eda)))

cat("\n--- SELECTED PRESENTATION CHARTS ---\n")
print(selected_plot_names)

# Print selected charts 
print(ppt_eda$injury_distribution)
print(ppt_eda$weekday)
print(ppt_eda$time_heatmap)
print(ppt_eda$accident_nature)
print(ppt_eda$driver_age)
print(ppt_eda$geography_slide)
print(ppt_eda$road_context)

out_dir <- "outputs/presentation_eda2024"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Export
plots_to_export <- c(selected_plot_names, "geography_slide")
for (nm in plots_to_export) {
  is_wide <- nm %in% c("time_heatmap", "geography_slide")
  ggplot2::ggsave(
    filename = file.path(out_dir, paste0(nm, ".png")),
    plot = ppt_eda[[nm]],
    width = if (is_wide) 12 else 10,
    height = if (is_wide) 6.5 else 5.5,
    dpi = 300,
    bg = "white"
  )
}

# Keep all informative tables
for (nm in names(eda_tables)) {
  utils::write.csv(
    as.data.frame(eda_tables[[nm]]),
    file = file.path(out_dir, paste0(nm, ".csv")),
    row.names = FALSE
  )
}
cat("\nPresentation EDA complete: ", length(selected_plot_names),
    " final charts and ", length(eda_tables), " supporting tables.\n", sep = "")

# ==========================================================================
# 16 | SUPERVISED MODELLING 
# ==========================================================================

if (RUN_MODELLING) {
# ── Package check ──────────────────────────────────────────────────────────
ml_packages <- c(
  "caret",        # train/test split + cross-validation framework
  "glmnet",       # logistic regression with regularisation
  "ranger",       # fast Random Forest
  "xgboost",      # gradient boosting
  "pROC",         # AUC-ROC curve
  "PRROC",        # Precision-Recall curve (better for imbalanced targets)
  "yardstick",    # tidy metric functions
  "cluster",      # silhouette score
  "factoextra",   # cluster visualisation helpers
  "fmsb",         # radar/spider charts
  "ggrepel",      # non-overlapping text labels
  "patchwork"     # combine plots
)
missing_ml <- setdiff(ml_packages, rownames(installed.packages()))
if (length(missing_ml) > 0) {
  message("Installing missing ML packages: ", paste(missing_ml, collapse = ", "))
  install.packages(missing_ml, quiet = TRUE)
}
suppressPackageStartupMessages({
  library(caret);    library(glmnet);  library(ranger)
  library(xgboost);  library(pROC);    library(PRROC)
  library(yardstick);library(cluster); library(factoextra)
  library(fmsb);     library(ggrepel); library(patchwork)
})

cat("\n============================================================\n")
cat("SECTION 16: SUPERVISED + UNSUPERVISED LEARNING\n")
cat("============================================================\n")

# ── Colour palette (consistent with EDA script) ────────────────────────────
ML_COLS <- c(
  LR   = "#2563A6",
  RF   = "#1a7f37",
  XGB  = "#bc4c00",
  base = "#636c76"
)
CLUSTER_COLS <- c("#bc4c00", "#0969da", "#636c76", "#8250df", "#1a7f37")
CLUSTER_NAMES <- c(
  "C0: Motorcycle/VRU",
  "C1: Car multi-vehicle",
  "C2: Pedestrian/single",
  "C3: Night run-off",
  "C4: Senior driver"
)

# ============================================================================
# 16.1 | FEATURE PREPARATION
# ============================================================================

cat("\n--- 16.1  Feature preparation ---\n")

# Binary flags 
BINARY_FEATS <- c(
  "is_night", "is_weekend",
  "is_motorway", "is_open", "at_junc",
  "bad_weather", "bad_surface",
  "is_headon", "is_rearend", "is_ped_hit", "is_runoff", "is_obstacle_hit",
  "is_moto_a", "is_vru_a", "cil_a_high",
  "has_veh_b", "is_single", "vb_is_heavy",
  "male_a", "young_a", "senior_a", "no_licence_a",
  "has_car", "has_motorcycle", "has_bicycle", "has_goods_vehicle"
)

# Numeric / ordinal
NUMERIC_FEATS <- c("n_vehicles")

# Raw categorical codes
CAT_FEATS <- c(
  "natura_incidente",
  "localizzazione_incidente",
  "tipo_di_strada",
  "tipo_veicolo_a"
)

# Hour: preserve unknown times as NA
accidents_ml_safe <- accidents_ml_safe %>%
  mutate(
    hour_num = suppressWarnings(as.integer(as.character(hour))),
    hour_num = ifelse(is.na(hour_num) | hour_num < 0 | hour_num > 23,
                      NA_integer_, hour_num),
    giorno_num = suppressWarnings(as.integer(as.character(giorno))),
    trimestre_num = suppressWarnings(as.integer(as.character(trimestre)))
  )

NUMERIC_FEATS <- c(NUMERIC_FEATS, "hour_num", "giorno_num", "trimestre_num")

# Critical: as.numeric(as.character(TRUE)) returns NA, not 1.
# Convert logical TRUE/FALSE directly to 1/0. Text 'TRUE'/'FALSE' is also handled.
to_binary_01 <- function(x) {
  if (is.logical(x) || is.numeric(x)) return(as.numeric(x))
  z <- tolower(trimws(as.character(x)))
  dplyr::case_when(
    z %in% c("true", "t", "1") ~ 1,
    z %in% c("false", "f", "0") ~ 0,
    TRUE ~ NA_real_
  )
}

ml_data <- accidents_ml_safe %>%
  mutate(across(all_of(BINARY_FEATS), to_binary_01)) %>%
  mutate(across(all_of(NUMERIC_FEATS),
                ~ suppressWarnings(as.numeric(as.character(.))))) %>%
  mutate(across(all_of(CAT_FEATS), as.character)) %>%
  mutate(y = factor(ifelse(multi_casualty, "MultiCas", "SingleCas")))

stopifnot(all(vapply(ml_data[BINARY_FEATS],
                     function(x) all(is.na(x) | x %in% c(0, 1)), logical(1))))

cat("  Feature matrix rows:", nrow(ml_data), "\n")
cat("  Binary features:   ", length(BINARY_FEATS), "\n")
cat("  Numeric features:  ", length(NUMERIC_FEATS), "\n")
cat("  Cat (TE):          ", length(CAT_FEATS), "\n")
cat("  Positive rate:     ", round(mean(ml_data$multi_casualty) * 100, 1), "%\n")

# ── Target encoding: category codes stay CHARACTER (preserve leading zeros). ─
SMOOTH <- 30

te_encode <- function(train_df, test_df, col, target_col, smooth = SMOOTH) {
  gm <- mean(as.numeric(train_df[[target_col]]), na.rm = TRUE)
  map_df <- train_df %>%
    group_by(cat = as.character(.data[[col]])) %>%
    summarise(s = sum(as.numeric(.data[[target_col]]), na.rm = TRUE),
              n = sum(!is.na(.data[[target_col]])), .groups = "drop") %>%
    mutate(te = (s + gm * smooth) / (n + smooth))
  lookup <- match(as.character(test_df[[col]]), map_df$cat)
  out <- map_df$te[lookup]
  out[is.na(out)] <- gm  # Unknown categories fall back to train prevalence.
  as.numeric(out)
}

# Cross-fitted TE for the outer-training rows: no row sees its own outcome.
# This safeguards the external CV predictions, though a fully nested
# hyperparameter-search pipeline would additionally re-fit TE in each inner fold.
te_encode_train_oof <- function(train_df, col, target_col, smooth = SMOOTH,
                                k = 5L) {
  fold <- caret::createFolds(as.factor(train_df[[target_col]]),
                            k = k, list = FALSE, returnTrain = FALSE)
  out <- numeric(nrow(train_df))
  for (j in seq_len(k)) {
    i <- which(fold == j)
    out[i] <- te_encode(train_df[-i, , drop = FALSE],
                        train_df[i, , drop = FALSE],
                        col, target_col, smooth)
  }
  out
}

# Replace missing/Inf values using the OUTER TRAINING fold's statistics.
# Completely unobserved columns become constants (0), never artificial signals.
impute_train_test <- function(X_tr, X_te) {
  for (j in seq_len(ncol(X_tr))) {
    observed <- X_tr[is.finite(X_tr[, j]), j]
    fill <- if (length(observed)) stats::median(observed) else 0
    X_tr[!is.finite(X_tr[, j]), j] <- fill
    X_te[!is.finite(X_te[, j]), j] <- fill
  }
  stopifnot(all(is.finite(X_tr)), all(is.finite(X_te)))
  list(train = X_tr, test = X_te)
}

# Constant columns have sd=0; scale(..., scale=TRUE) would create NaN.
scale_train_test <- function(X_tr, X_te) {
  center <- colMeans(X_tr)
  scale_factor <- apply(X_tr, 2, stats::sd)
  scale_factor[!is.finite(scale_factor) | scale_factor == 0] <- 1
  tr <- sweep(sweep(X_tr, 2, center, "-"), 2, scale_factor, "/")
  te <- sweep(sweep(X_te, 2, center, "-"), 2, scale_factor, "/")
  stopifnot(all(is.finite(tr)), all(is.finite(te)))
  list(train = tr, test = te)
}

# ============================================================================
# 16.2 | 5-FOLD STRATIFIED CROSS-VALIDATION SETUP (2024 only)
# ============================================================================
# Because we have only one year, there is no temporal split.
# We use 5-fold stratified CV: each fold uses 80% as train and 20% as test.
# Target encoding is recomputed inside every fold so there is no leakage.
# Out-of-fold predictions are collected and used for all final evaluation.
# ============================================================================

cat("\n--- 16.2  5-fold stratified CV (2024, n=", nrow(ml_data), ") ---\n")

N_FOLDS <- 5
set.seed(42)
fold_ids <- createFolds(ml_data$y, k = N_FOLDS, list = FALSE, returnTrain = FALSE)
# fold_ids[i] = which fold observation i belongs to (1..5)

cat("  Folds created. Distribution:\n")
print(table(fold_ids))

# Helper: build imputed numeric matrices, with all preprocessing fit on train.
build_X_fold <- function(tr, te) {
  tr_bin <- tr[, BINARY_FEATS, drop = FALSE]
  te_bin <- te[, BINARY_FEATS, drop = FALSE]
  tr_num <- tr[, NUMERIC_FEATS, drop = FALSE]
  te_num <- te[, NUMERIC_FEATS, drop = FALSE]

  for (col in CAT_FEATS) {
    te_col <- paste0(col, "_te")
    tr_bin[[te_col]] <- te_encode_train_oof(tr, col, "multi_casualty")
    te_bin[[te_col]] <- te_encode(tr, te, col, "multi_casualty")
  }

  X_tr <- as.matrix(cbind(tr_bin, tr_num))
  X_te <- as.matrix(cbind(te_bin, te_num))
  storage.mode(X_tr) <- "double"
  storage.mode(X_te) <- "double"
  imputed <- impute_train_test(X_tr, X_te)

  list(
    X_tr = imputed$train,
    X_te = imputed$test,
    # ranger's class.weights expect levels '0'/'1', not 'FALSE'/'TRUE'.
    y_tr = as.integer(tr$multi_casualty),
    y_te = as.integer(te$multi_casualty)
  )
}

# Determine feature count from fold 1
f1_tr  <- ml_data[fold_ids != 1, ]
f1_te  <- ml_data[fold_ids == 1, ]
f1_mat <- build_X_fold(f1_tr, f1_te)
FEAT_NAMES <- colnames(f1_mat$X_tr)
N_FEATS    <- length(FEAT_NAMES)
cat("  Features per fold:", N_FEATS, "\n")

# Global pos_wt based on full dataset (used for model fitting inside folds)
pos_wt <- sum(ml_data$multi_casualty == 0) / sum(ml_data$multi_casualty == 1)
cat("  scale_pos_weight:", round(pos_wt, 2), "\n")
cat("  Overall MC rate: ", round(mean(ml_data$multi_casualty) * 100, 1), "%\n\n")

# ── Run 5-fold CV for all three models ───────────────────────────────────────
# For each fold we:
#   1. Build train/test matrices with fold-specific TE
#   2. Fit LR, RF, XGB on train
#   3. Collect out-of-fold (OOF) probabilities for test rows
# After all folds: evaluate on the full OOF probability vector

oof_lr  <- numeric(nrow(ml_data))
oof_rf  <- numeric(nrow(ml_data))
oof_xgb <- numeric(nrow(ml_data))

# Store per-fold AUC for stability reporting
fold_auc_lr  <- numeric(N_FOLDS)
fold_auc_rf  <- numeric(N_FOLDS)
fold_auc_xgb <- numeric(N_FOLDS)

# Also keep the last fold's matrices for Platt calibration (section 16.6)
last_fold_tr <- NULL
last_fold_te <- NULL
last_fold_mat <- NULL

cat("  Running 5-fold CV...\n")
cat(sprintf("  %-6s %-8s %-8s %-8s\n", "Fold", "LR AUC", "RF AUC", "XGB AUC"))
cat("  ", rep("-", 36), "\n", sep="")

for (fold in seq_len(N_FOLDS)) {
  tr_df <- ml_data[fold_ids != fold, ]
  te_df <- ml_data[fold_ids == fold, ]
  mf    <- build_X_fold(tr_df, te_df)

  X_tr <- mf$X_tr;  y_tr <- mf$y_tr
  X_te <- mf$X_te;  y_te <- mf$y_te
  pw   <- sum(y_tr == 0) / sum(y_tr == 1)
  te_idx <- which(fold_ids == fold)

  # ── LR (glmnet ridge) ──────────────────────────────────────────────────────
  scaled_fold <- scale_train_test(X_tr, X_te)
  X_tr_sc <- scaled_fold$train
  X_te_sc <- scaled_fold$test
  w_lr <- ifelse(y_tr == 1, pw, 1)
  cv_lr_fold <- cv.glmnet(X_tr_sc, y_tr, family = "binomial",
                           alpha = 0, nfolds = 3,
                           weights = w_lr, type.measure = "auc")
  oof_lr[te_idx] <- as.numeric(predict(cv_lr_fold, newx = X_te_sc,
                                        s = "lambda.1se", type = "response"))

  # ── RF (ranger) ────────────────────────────────────────────────────────────
  rf_fold <- ranger(x = X_tr, y = as.factor(y_tr),
                    num.trees = 300, max.depth = 12, min.node.size = 50,
                    class.weights = c("0"=1,"1"=pw),
                    importance = if(fold==1)"impurity" else "none",
                    probability = TRUE, seed = fold)
  oof_rf[te_idx] <- predict(rf_fold, data = X_te)$predictions[, "1"]
  if (fold == 1) rf_model_fold1 <- rf_fold   # save for feature importance

  # ── XGBoost ────────────────────────────────────────────────────────────────
  dtr <- xgb.DMatrix(X_tr, label = y_tr)
  dte <- xgb.DMatrix(X_te, label = y_te)
  xgb_fold <- xgb.train(
    params    = list(booster="gbtree", objective="binary:logistic",
                     eval_metric="logloss", max_depth=6, eta=0.05,
                     subsample=0.80, colsample_bytree=0.80,
                     scale_pos_weight=pw, seed=fold),
    data      = dtr,
    nrounds   = 500,
    watchlist = list(train=dtr, test=dte),
    verbose   = 0
  )
  oof_xgb[te_idx] <- predict(xgb_fold, dte)
  if (fold == N_FOLDS) {
    xgb_model_last <- xgb_fold
    last_fold_tr   <- tr_df
    last_fold_te   <- te_df
    last_fold_mat  <- mf
    xgb_fi <- xgb.importance(model = xgb_fold) %>%
      as.data.frame() %>% arrange(desc(Gain)) %>%
      mutate(Gain_pct = Gain / sum(Gain) * 100)
  }

  # Per-fold AUC
  fold_auc_lr[fold]  <- as.numeric(pROC::auc(
    pROC::roc(y_te, oof_lr[te_idx],  quiet=TRUE)))
  fold_auc_rf[fold]  <- as.numeric(pROC::auc(
    pROC::roc(y_te, oof_rf[te_idx],  quiet=TRUE)))
  fold_auc_xgb[fold] <- as.numeric(pROC::auc(
    pROC::roc(y_te, oof_xgb[te_idx], quiet=TRUE)))

  cat(sprintf("  Fold %d  %.3f    %.3f    %.3f\n",
              fold, fold_auc_lr[fold], fold_auc_rf[fold], fold_auc_xgb[fold]))
}

cat(sprintf("  Mean   %.3f±%.3f %.3f±%.3f %.3f±%.3f\n",
            mean(fold_auc_lr),  sd(fold_auc_lr),
            mean(fold_auc_rf),  sd(fold_auc_rf),
            mean(fold_auc_xgb), sd(fold_auc_xgb)))

# OOF vectors are now complete — use them as y_pred for all evaluation
# y_test = full ml_data labels; oof_* = out-of-fold predicted probabilities
y_test  <- as.integer(ml_data$multi_casualty)
lr_probs  <- oof_lr
rf_probs  <- oof_rf
xgb_probs <- oof_xgb

# Keep train/test from last fold for Platt calibration (section 16.6)
X_train <- last_fold_mat$X_tr;  X_test <- last_fold_mat$X_te
y_train <- last_fold_mat$y_tr

# Scaled X for LR coefficient extraction (using last fold)
scaled_last <- scale_train_test(X_train, X_test)
X_train_sc <- scaled_last$train
X_test_sc  <- scaled_last$test

# LR coefficients from last fold
w_last <- ifelse(y_train == 1, sum(y_train==0)/sum(y_train==1), 1)
cv_lr_last <- cv.glmnet(X_train_sc, y_train, family="binomial",
                         alpha=0, nfolds=3, weights=w_last,
                         type.measure="auc")
lr_coef <- coef(cv_lr_last, s="lambda.1se")
lr_coef_df <- data.frame(
  feature = rownames(lr_coef)[-1],
  coef    = as.numeric(lr_coef)[-1]
) %>% filter(coef != 0) %>% arrange(desc(abs(coef)))

cat("\n  LR non-zero coefficients:", nrow(lr_coef_df), "\n")
cat("  RF OOB error (fold 1):   ", round(rf_model_fold1$prediction.error*100,2), "%\n")

# ============================================================================
# 16.3 | LOGISTIC REGRESSION — summary from CV
# ============================================================================
cat("\n--- 16.3  Logistic Regression summary (from CV) ---\n")
cat("  OOF predictions collected across 5 folds.\n")
cat("  Ridge (L2) regularisation, lambda by 3-fold inner CV per outer fold.\n")
cat("  Top 5 coefficients (last fold, absolute value):\n")
print(head(lr_coef_df, 5))

# ============================================================================
# 16.4 | RANDOM FOREST — summary from CV
# ============================================================================
cat("\n--- 16.4  Random Forest summary (from CV) ---\n")
cat("  300 trees per fold, max_depth=12, min_node_size=50.\n")
cat("  OOB error (fold 1):", round(rf_model_fold1$prediction.error*100,2), "%\n")

# RF feature importance from fold 1
rf_fi <- data.frame(
  feature    = FEAT_NAMES,
  importance = rf_model_fold1$variable.importance
) %>% arrange(desc(importance)) %>%
  mutate(importance_pct = importance / sum(importance) * 100)

cat("  Top 5 RF features:\n")
print(head(rf_fi[, c("feature","importance_pct")], 5))

# ============================================================================
# 16.5 | XGBOOST — summary from CV
# ============================================================================
cat("\n--- 16.5  XGBoost summary (from CV) ---\n")
cat("  500 trees, eta=0.05, max_depth=6, subsample=0.80.\n")
cat("  Top 5 XGBoost features (last fold Gain):\n")
print(head(xgb_fi[, c("Feature","Gain_pct")], 5))

# Store XGB params for reference
xgb_params <- list(
  booster="gbtree", objective="binary:logistic", eval_metric="logloss",
  max_depth=6, eta=0.05, subsample=0.80, colsample_bytree=0.80,
  scale_pos_weight=pos_wt, seed=42
)

# ============================================================================
# 16.6 | PLATT SCALING CALIBRATION
# ============================================================================
# With 5-fold CV and 2024-only data the OOF XGBoost predictions are already
# produced on unseen data, so we can fit Platt directly on the full OOF vector.
# The Platt logistic regression maps raw OOF scores → calibrated probabilities.
# IMPORTANT: evaluating calibration, Brier and chosen threshold on the SAME
# OOF rows used to fit Platt can be optimistic. For a final publication use
# nested/cross-fitted calibration or a separate untouched holdout.
# ============================================================================

cat("\n--- 16.6  Platt scaling calibration (XGBoost OOF) ---\n")

# Fit Platt on the full out-of-fold XGBoost predictions
platt_fit     <- glm(y_test ~ xgb_probs, family = binomial())
xgb_cal_probs <- predict(platt_fit,
                          newdata = data.frame(xgb_probs = xgb_probs),
                          type    = "response")

cat("  Mean predicted before:", round(mean(xgb_probs) * 100, 1), "%\n")
cat("  Mean predicted after: ", round(mean(xgb_cal_probs) * 100, 1), "%\n")
cat("  True positive rate:   ", round(mean(y_test) * 100, 1), "%\n")

brier_before <- mean((xgb_probs - y_test)^2)
brier_after  <- mean((xgb_cal_probs - y_test)^2)
cat("  Brier score before:", round(brier_before, 4), "\n")
cat("  Brier score after: ", round(brier_after,  4), "\n")

# ============================================================================
# 16.7 | EVALUATION METRICS
# ============================================================================

cat("\n--- 16.7  Evaluation metrics ---\n")

THRESHOLD <- 0.25   # Working decision threshold; evaluate independently

eval_model <- function(probs, y_true, name, thr = THRESHOLD) {
  roc_obj  <- pROC::roc(y_true, probs, quiet = TRUE)
  auc_val  <- as.numeric(pROC::auc(roc_obj))
  pr_obj   <- PRROC::pr.curve(scores.class0 = probs,
                               weights.class0 = y_true, curve = FALSE)
  ap_val   <- pr_obj$auc.integral
  base     <- mean(y_true)
  yhat     <- as.integer(probs >= thr)
  cm       <- table(Actual = y_true, Predicted = yhat)
  tn <- cm["0","0"]; fp <- cm["0","1"]
  fn <- cm["1","0"]; tp <- cm["1","1"]
  prec  <- ifelse((tp+fp)>0, tp/(tp+fp), 0)
  rec   <- ifelse((tp+fn)>0, tp/(tp+fn), 0)
  spec  <- ifelse((tn+fp)>0, tn/(tn+fp), 0)
  f1    <- ifelse((prec+rec)>0, 2*prec*rec/(prec+rec), 0)
  f2    <- ifelse((4*prec+rec)>0, 5*prec*rec/(4*prec+rec), 0)
  ba    <- (rec + spec) / 2
  brier <- mean((probs - y_true)^2)
  mcc_num <- tp*tn - fp*fn
  mcc_den <- sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
  mcc   <- ifelse(mcc_den > 0, mcc_num / mcc_den, 0)
  list(name=name, auc=auc_val, ap=ap_val, ap_lift=ap_val/base,
       prec=prec, rec=rec, spec=spec, f1=f1, f2=f2, ba=ba,
       brier=brier, mcc=mcc,
       tp=tp, fp=fp, fn=fn, tn=tn, thr=thr,
       roc_obj=roc_obj, probs=probs)
}

metrics_lr  <- eval_model(lr_probs,       y_test, "Logistic Regression")
metrics_rf  <- eval_model(rf_probs,       y_test, "Random Forest")
metrics_xgb <- eval_model(xgb_probs,     y_test, "XGBoost (raw)")
metrics_cal <- eval_model(xgb_cal_probs, y_test, "XGBoost + Platt")

# Summary table
metrics_all <- list(metrics_lr, metrics_rf, metrics_xgb, metrics_cal)
metrics_tbl <- do.call(rbind, lapply(metrics_all, function(m) {
  data.frame(Model     = m$name,
             AUC       = round(m$auc,  3),
             AP        = round(m$ap,   3),
             AP_lift   = round(m$ap_lift, 2),
             F1        = round(m$f1,   3),
             F2        = round(m$f2,   3),
             Bal_Acc   = round(m$ba,   3),
             MCC       = round(m$mcc,  3),
             Precision = paste0(round(m$prec*100, 1), "%"),
             Recall    = paste0(round(m$rec*100,  1), "%"),
             Brier     = round(m$brier, 4))
}))
cat("\n  Full metric table (threshold =", THRESHOLD, "):\n")
print(metrics_tbl, row.names = FALSE)

# Detailed confusion matrix — calibrated XGBoost
cat("\n  Confusion matrix — XGBoost + Platt (thr =", THRESHOLD, "):\n")
cat("  ┌──────────────────────────────────────────┐\n")
cat("  │                Pred NO      Pred YES     │\n")
cat(sprintf("  │  Actual NO   %9d    %9d     │\n", metrics_cal$tn, metrics_cal$fp))
cat(sprintf("  │  Actual YES  %9d    %9d     │\n", metrics_cal$fn, metrics_cal$tp))
cat("  └──────────────────────────────────────────┘\n")

# Threshold analysis
cat("\n  Threshold analysis — XGBoost + Platt:\n")
cat(sprintf("  %-10s %-10s %-8s %-8s %-10s\n",
            "Threshold","Precision","Recall","F1","Flagged"))
for (thr in c(0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50)) {
  yhat  <- as.integer(xgb_cal_probs >= thr)
  cm_t  <- table(y_test, yhat)
  if (all(c("0","1") %in% rownames(cm_t)) && "1" %in% colnames(cm_t)) {
    tp_t <- cm_t["1","1"]; fp_t <- cm_t["0","1"]
    fn_t <- cm_t["1","0"]
    p <- tp_t/(tp_t+fp_t); r <- tp_t/(tp_t+fn_t)
    f <- 2*p*r/(p+r)
    cat(sprintf("  %-10.2f %-10s %-8s %-8s %-10d\n",
                thr, paste0(round(p*100,1),"%"), paste0(round(r*100,1),"%"),
                round(f,3), sum(yhat)))
  }
}

# ============================================================================
# 16.8 | SUPERVISED PLOTS
# ============================================================================

cat("\n--- 16.8  Supervised plots ---\n")

plots_ml <- list()

# ── Plot S1: Model comparison bar chart ──────────────────────────────────────
comp_long <- metrics_tbl %>%
  select(Model, AUC, AP, Bal_Acc, MCC) %>%
  tidyr::pivot_longer(-Model, names_to = "Metric", values_to = "Value") %>%
  mutate(col = case_when(
    grepl("Logistic",  Model) ~ ML_COLS["LR"],
    grepl("Random",    Model) ~ ML_COLS["RF"],
    grepl("raw",       Model) ~ ML_COLS["XGB"],
    TRUE                      ~ "#C44949"
  ))

plots_ml$model_comparison <- ggplot(comp_long,
  aes(x = Model, y = Value, fill = col)) +
  geom_col(width = 0.65, show.legend = FALSE) +
  geom_text(aes(label = round(Value, 3)),
            vjust = -0.4, size = 3.2) +
  scale_fill_identity() +
  scale_y_continuous(limits = c(0, 1.08)) +
  facet_wrap(~Metric, nrow = 1) +
  labs(title    = "Supervised Model Comparison — Key Metrics",
       subtitle = paste("5-fold out-of-fold CV (2024) · threshold =", THRESHOLD,
                        "· target = feriti ≥ 2 (multi_casualty)"),
       x = NULL, y = "Score") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x  = element_text(angle = 25, hjust = 1, size = 8),
        panel.grid.major.x = element_blank(),
        strip.text   = element_text(face = "bold"))

# ── Plot S2: ROC curves (all 3 models) ───────────────────────────────────────
roc_data <- bind_rows(lapply(list(metrics_lr, metrics_rf, metrics_cal), function(m) {
  roc_df <- data.frame(
    fpr  = 1 - m$roc_obj$specificities,
    tpr  = m$roc_obj$sensitivities,
    model= paste0(m$name, " (AUC=", round(m$auc, 3), ")")
  )
  roc_df[order(roc_df$fpr), ]
}))

plots_ml$roc_curves <- ggplot(roc_data, aes(x = fpr, y = tpr, colour = model)) +
  geom_line(linewidth = 1.1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              colour = "grey60", linewidth = 0.8) +
  scale_colour_manual(values = c(ML_COLS["LR"], ML_COLS["RF"], "#C44949")) +
  scale_x_continuous(labels = scales::percent_format()) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(title    = "ROC Curves — All Models",
       subtitle = "Dashed diagonal = random classifier (AUC = 0.50)",
       x        = "False Positive Rate (1 − Specificity)",
       y        = "True Positive Rate (Sensitivity)",
       colour   = NULL) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

# ── Plot S3: Precision-Recall curve (calibrated XGBoost) ─────────────────────
pr_curve_data <- data.frame(
  threshold  = seq(0.05, 0.95, by = 0.01)
) %>%
  rowwise() %>%
  mutate(
    yhat  = list(as.integer(xgb_cal_probs >= threshold)),
    tp    = sum(unlist(yhat) == 1 & y_test == 1),
    fp    = sum(unlist(yhat) == 1 & y_test == 0),
    fn    = sum(unlist(yhat) == 0 & y_test == 1),
    prec  = ifelse((tp+fp)>0, tp/(tp+fp), NA_real_),
    rec   = ifelse((tp+fn)>0, tp/(tp+fn), NA_real_),
    f1    = ifelse(!is.na(prec) & !is.na(rec) & (prec+rec)>0,
                   2*prec*rec/(prec+rec), NA_real_)
  ) %>%
  ungroup() %>%
  filter(!is.na(prec), !is.na(rec))

best_f1_row <- pr_curve_data %>% slice_max(f1, n = 1)

plots_ml$pr_curve <- ggplot(pr_curve_data,
  aes(x = rec, y = prec)) +
  geom_line(colour = ML_COLS["XGB"], linewidth = 1.3) +
  geom_hline(yintercept = mean(y_test),
             linetype = "dashed", colour = "grey50", linewidth = 0.8) +
  geom_point(data = best_f1_row,
             aes(x = rec, y = prec),
             colour = "#C44949", size = 4) +
  geom_label(data = best_f1_row,
             aes(label = paste0("Best F1=", round(f1, 3),
                                "\nthr=", round(threshold, 2))),
             nudge_x = 0.04, size = 3.2, colour = "#C44949") +
  annotate("text", x = 0.9, y = mean(y_test) + 0.02,
           label = paste0("Random baseline = ", round(mean(y_test)*100, 1), "%"),
           colour = "grey50", size = 3.2, hjust = 1) +
  scale_x_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
  scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
  labs(title    = "Precision-Recall Curve — XGBoost + Platt Calibration",
       subtitle = paste0("AP = ", round(metrics_cal$ap, 3),
                         " vs random baseline = ", round(mean(y_test), 3),
                         " · AP lift = ", round(metrics_cal$ap_lift, 2), "×"),
       x = "Recall", y = "Precision") +
  theme_minimal(base_size = 12)

# ── Plot S4: Calibration curves before / after Platt ─────────────────────────
cal_bins <- 10
make_cal_data <- function(probs, y, name) {
  cuts  <- cut(probs, breaks = seq(0, 1, length.out = cal_bins + 1),
               include.lowest = TRUE)
  data.frame(probs = probs, y = y, bin = cuts) %>%
    group_by(bin) %>%
    summarise(mean_pred = mean(probs),
              frac_pos  = mean(y),
              n         = n(), .groups = "drop") %>%
    mutate(model = name)
}

cal_data <- bind_rows(
  make_cal_data(xgb_probs,     y_test, paste0("XGBoost raw\n(mean=",
                round(mean(xgb_probs)*100,1),"%)")),
  make_cal_data(xgb_cal_probs, y_test, paste0("XGBoost + Platt\n(mean=",
                round(mean(xgb_cal_probs)*100,1),"%)"))
)

plots_ml$calibration <- ggplot(cal_data,
  aes(x = mean_pred, y = frac_pos, colour = model, size = n)) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", colour = "grey60", linewidth = 0.8) +
  geom_line(linewidth = 1, show.legend = TRUE) +
  geom_point(show.legend = FALSE) +
  scale_colour_manual(values = c(ML_COLS["XGB"], "#C44949")) +
  scale_size_continuous(range = c(2, 6), guide = "none") +
  scale_x_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
  scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
  labs(title    = "Calibration Curves — Before and After Platt Scaling",
       subtitle = paste0("Brier: before = ", round(brier_before, 4),
                         " → after = ", round(brier_after, 4),
                         " (", round((1 - brier_after/brier_before)*100, 0),
                         "% improvement)"),
       x = "Mean predicted probability", y = "Observed fraction positive",
       colour = NULL) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

# ── Plot S5: Feature importance (XGBoost) ─────────────────────────────────────
plots_ml$feature_importance <- ggplot(
  head(xgb_fi, 20),
  aes(x = reorder(Feature, Gain_pct), y = Gain_pct)) +
  geom_col(fill = ML_COLS["XGB"], width = 0.7) +
  geom_text(aes(label = paste0(round(Gain_pct, 1), "%")),
            hjust = -0.15, size = 3.2) +
  coord_flip() +
  scale_y_continuous(limits = c(0, max(head(xgb_fi, 20)$Gain_pct) * 1.2)) +
  labs(title    = "XGBoost Feature Importance — Top 20 (Gain)",
       subtitle = "★ = target-encoded categorical · higher = more important",
       x = NULL, y = "Gain (%)") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.major.y = element_blank())

# ── Plot S6: Logistic Regression coefficients ─────────────────────────────────
plots_ml$lr_coefficients <- ggplot(
  head(lr_coef_df, 15),
  aes(x = reorder(feature, coef), y = coef,
      fill = coef > 0)) +
  geom_col(width = 0.7, show.legend = FALSE) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  geom_text(aes(label = round(coef, 3),
                hjust = ifelse(coef >= 0, -0.15, 1.15)),
            size = 3.2) +
  coord_flip() +
  scale_fill_manual(values = c("TRUE" = ML_COLS["LR"], "FALSE" = "#C44949")) +
  labs(title    = "Logistic Regression — Standardised Coefficients (top 15)",
       subtitle = "Positive = increases multi-casualty probability · Ridge (L2) regularisation",
       x = NULL, y = "Coefficient") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.major.y = element_blank())

# ── Plot S7: Threshold analysis ───────────────────────────────────────────────
thr_seq <- seq(0.10, 0.70, by = 0.05)
thr_rows <- lapply(thr_seq, function(thr) {
  yhat <- as.integer(xgb_cal_probs >= thr)
  cm_t <- table(factor(y_test, levels=c(0,1)),
                factor(yhat,   levels=c(0,1)))
  tp_t <- cm_t["1","1"]; fp_t <- cm_t["0","1"]; fn_t <- cm_t["1","0"]
  p <- ifelse((tp_t+fp_t)>0, tp_t/(tp_t+fp_t), 0)
  r <- ifelse((tp_t+fn_t)>0, tp_t/(tp_t+fn_t), 0)
  f1<- ifelse((p+r)>0, 2*p*r/(p+r), 0)
  sp<- 1 - fp_t/(fp_t + cm_t["0","0"])
  data.frame(threshold=thr, Precision=p, Recall=r, F1=f1, Specificity=sp)
})
thr_df <- bind_rows(thr_rows) %>%
  tidyr::pivot_longer(-threshold, names_to="Metric", values_to="Value")

plots_ml$threshold_analysis <- ggplot(thr_df,
  aes(x = threshold, y = Value, colour = Metric)) +
  geom_line(linewidth = 1.1) +
  geom_vline(xintercept = THRESHOLD, linetype = "dashed",
             colour = "grey50", linewidth = 0.8) +
  annotate("text", x = THRESHOLD + 0.02, y = 0.1,
           label = paste0("default\nthr=", THRESHOLD),
           size = 3, colour = "grey50", hjust = 0) +
  scale_colour_manual(values = c(
    Precision = ML_COLS["LR"], Recall = ML_COLS["RF"],
    F1 = ML_COLS["XGB"], Specificity = "#8250df")) +
  scale_x_continuous(breaks = seq(0.1, 0.7, 0.1)) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(title    = "Threshold Analysis — XGBoost + Platt Calibration",
       subtitle = "Raising the threshold increases precision but reduces recall",
       x = "Classification threshold", y = "Metric value",
       colour   = NULL) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

# ── Plot S8: Combined supervised overview ─────────────────────────────────────
plots_ml$supervised_overview <-
  (plots_ml$roc_curves | plots_ml$pr_curve) /
  (plots_ml$calibration | plots_ml$threshold_analysis) +
  patchwork::plot_annotation(
    title    = "Supervised Learning — Model Evaluation Overview",
    subtitle = "Top: ranking quality | Bottom: calibration and threshold choice",
    theme    = theme(plot.title    = element_text(size = 13, face = "bold"),
                     plot.subtitle = element_text(size = 10, colour = "grey40"))
  )

# ============================================================================
# 16.8B | SUPERVISED VISUALISATIONS
# ============================================================================

# Precision–Recall curve for each original OOF probability vector.
pr_points_ooof <- function(probs, truth, model_label, max_points = 3500L) {
  stopifnot(length(probs) == length(truth),
            all(is.finite(probs)), all(truth %in% c(0L, 1L)))
  o <- order(probs, decreasing = TRUE)
  p <- probs[o]
  y <- truth[o]
  cum_tp <- cumsum(y == 1L)
  last_in_tie <- c(which(diff(p) != 0), length(p))
  keep <- unique(as.integer(round(seq(1, length(last_in_tie),
                                      length.out = min(max_points, length(last_in_tie))))))
  at <- last_in_tie[keep]
  rbind(
    data.frame(recall = 0, precision = 1, model = model_label),
    data.frame(recall = cum_tp[at] / sum(y == 1L),
               precision = cum_tp[at] / at,
               model = model_label)
  )
}

pr_names <- c(
  LR  = sprintf("Logistic regression (AP %.3f)", metrics_lr$ap),
  RF  = sprintf("Random forest (AP %.3f)", metrics_rf$ap),
  XGB = sprintf("XGBoost raw (AP %.3f)", metrics_xgb$ap)
)
pr_all_models_data <- dplyr::bind_rows(
  pr_points_ooof(lr_probs,  y_test, pr_names[["LR"]]),
  pr_points_ooof(rf_probs,  y_test, pr_names[["RF"]]),
  pr_points_ooof(xgb_probs, y_test, pr_names[["XGB"]])
)
plots_ml$pr_models <- ggplot(pr_all_models_data,
                             aes(x = recall, y = precision, colour = model)) +
  geom_line(linewidth = 1.0) +
  geom_hline(yintercept = mean(y_test), colour = "#7C8792",
             linewidth = 0.8, linetype = "dashed") +
  annotate("label", x = 0.79, y = mean(y_test) + 0.035,
           label = sprintf("Baseline %.1f%%", 100 * mean(y_test)),
           size = 3.1, colour = "#596575") +
  scale_colour_manual(values = stats::setNames(
    unname(ML_COLS[c("LR", "RF", "XGB")]),
    unname(pr_names[c("LR", "RF", "XGB")])
  )) +
  scale_x_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
  scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
  labs(title = "Precision–Recall Curves — Three Classifiers",
       subtitle = "2024, 5-fold out-of-fold predictions | dashed = observed positive rate",
       x = "Recall", y = "Precision", colour = NULL) +
  theme_minimal(base_size = 12) + theme(legend.position = "bottom")

# Confusion matrix from the EXISTING reported calibrated XGBoost decision rule.
# Note: calibration is evaluated in-sample on OOF scores in the original script.
cm_actual <- factor(y_test, levels = c(0, 1),
                    labels = c("0–1 injured", "2+ injured"))
cm_pred <- factor(as.integer(xgb_cal_probs >= THRESHOLD),
                  levels = c(0, 1), labels = c("Predict 0–1", "Predict 2+"))
confusion_plot_data <- as.data.frame(table(Actual = cm_actual, Predicted = cm_pred)) %>%
  group_by(Actual) %>%
  mutate(row_pct = if (sum(Freq) > 0) 100 * Freq / sum(Freq) else 0,
         tile_label = sprintf("%s\n%.1f%% of actual class",
                              scales::comma(Freq), row_pct)) %>%
  ungroup()

plots_ml$confusion_matrix <- ggplot(confusion_plot_data,
                                    aes(x = Predicted, y = Actual, fill = row_pct)) +
  geom_tile(colour = "white", linewidth = 1.8) +
  geom_text(aes(label = tile_label), size = 4, lineheight = 1.1,
            fontface = "bold", colour = "#172B42") +
  scale_fill_gradient(low = "#EAF2FA", high = "#83B4E0", name = "Row share (%)") +
  scale_y_discrete(limits = rev(levels(cm_actual))) +
  coord_equal() +
  labs(title = "Confusion Matrix — XGBoost + Platt",
       subtitle = sprintf("Existing OOF-based results | working threshold = %.2f", THRESHOLD),
       caption = "Percentages are calculated within each actual class. Calibration fit uses the same OOF rows (not an independent calibration test).",
       x = "Predicted class", y = "Actual class") +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(face = "bold"),
        legend.position = "right")

# The original Gain feature importance is computed on the LAST CV fold only.
# State this prominently: it is not a stable across-fold importance estimate.
xgb_top10_data <- xgb_fi %>%
  arrange(desc(Gain_pct)) %>%
  slice_head(n = 10) %>%
  mutate(display_feature = gsub("_", " ", Feature, fixed = TRUE),
         display_feature = ifelse(grepl("_te$", Feature),
                                  paste0(display_feature, " (TE)"), display_feature),
         display_feature = factor(display_feature,
                                  levels = rev(display_feature)))
plots_ml$xgb_importance_top10 <- ggplot(xgb_top10_data,
      aes(x = display_feature, y = Gain_pct)) +
  geom_col(fill = unname(ML_COLS[["XGB"]]), width = 0.68) +
  geom_text(aes(label = sprintf("%.1f%%", Gain_pct)),
            hjust = -0.15, fontface = "bold", size = 3.6) +
  coord_flip(clip = "off") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
  labs(title = "XGBoost — Top 10 Predictors by Gain",
       subtitle = "Importance from the last cross-validation fold only | TE = target encoding",
       x = NULL, y = "Share of total feature Gain (%)") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.major.y = element_blank(),
        plot.title = element_text(face = "bold"))

if (SHOW_FINAL_PLOTS) {
  print(plots_ml$model_comparison)
  print(plots_ml$pr_models)
  print(plots_ml$confusion_matrix)
  print(plots_ml$xgb_importance_top10)
}
if (EXPORT_FINAL_PLOTS) {
  supervised_output_dir <- "outputs/presentation_ml_2024"
  dir.create(supervised_output_dir, recursive = TRUE, showWarnings = FALSE)
  supervised_to_save <- list(
    model_comparison = plots_ml$model_comparison,
    precision_recall_three_models = plots_ml$pr_models,
    confusion_matrix_xgb = plots_ml$confusion_matrix,
    xgb_feature_importance_top10 = plots_ml$xgb_importance_top10
  )
  for (nm in names(supervised_to_save)) {
    ggplot2::ggsave(file.path(supervised_output_dir, paste0(nm, ".png")),
                    plot = supervised_to_save[[nm]],
                    width = if (nm == "model_comparison") 13 else 10,
                    height = 6.3, dpi = 300, bg = "white")
  }
  utils::write.csv(metrics_tbl, file.path(supervised_output_dir, "model_metrics.csv"),
                   row.names = FALSE)
  utils::write.csv(as.data.frame(confusion_plot_data),
                   file.path(supervised_output_dir, "confusion_matrix.csv"),
                   row.names = FALSE)
  utils::write.csv(xgb_top10_data,
                   file.path(supervised_output_dir, "xgboost_gain_top10.csv"),
                   row.names = FALSE)
  cat("\nSaved supervised Canva graphics in ", supervised_output_dir, "\n", sep = "")
}


}  

if (RUN_CLUSTERING) {
# ============================================================================
# 16.9 | FINAL K=6 CLUSTERING — 16 BINARY FEATURES
# ============================================================================.

cat("============================================================\n")
cat("  ISTAT 2024 — K-MEANS CLUSTERING (k = 6, 16 binary features)\n")
cat("============================================================\n\n")

# ── 1. Build the 16 binary features 
cat("--- 1. Building 16 binary features ---\n")

# The validated main data contain raw ISTAT fields that accidents_ml_safe omits.
# accidents_main for clustering; the 16 predictors below do NOT include injuries.
stopifnot(exists("accidents_main"), exists("accidents_ml_safe"))
df <- accidents_main

num <- function(x) suppressWarnings(as.numeric(as.character(x)))

F <- data.frame(
  Motorcycle        = as.integer(num(df$tipo_veicolo_a) %in% c(15,16,17)),
  Cyclist_escooter  = as.integer(num(df$tipo_veicolo_a) %in% c(14,22)),
  Heavy_vehicle     = as.integer(num(df$tipo_veicolo_a) %in% c(8,9,10,12)),
  X2nd_vehicle      = as.integer(!is.na(df$tipo_veicoli__b_) &
                                   as.character(df$tipo_veicoli__b_) != ""),
  Pedestrian_hit    = as.integer(num(df$natura_incidente) == 5),
  Run_off           = as.integer(num(df$natura_incidente) == 10),
  Rear_side         = as.integer(num(df$natura_incidente) %in% c(3,4)),
  Head_on           = as.integer(num(df$natura_incidente) %in% c(1,2)),
  Night             = as.integer({
                        h <- num(df$Ora)
                        !is.na(h) & (h >= 22 | h <= 5)
                      }),
  Weekend           = as.integer(num(df$giorno) %in% c(6,7)),
  Young_18_29       = as.integer(trimws(as.character(
                        df$veicolo__a___et__conducente)) == "18-29"),
  Senior_65         = as.integer(trimws(as.character(
                        df$veicolo__a___et__conducente)) == "65+"),
  Urban             = as.integer(num(df$localizzazione_incidente) %in% c(0,1,2,3)),
  Motorway          = as.integer(num(df$localizzazione_incidente) == 7),
  Intersection      = as.integer(num(df$intersezione_o_non_interse3) %in% c(1,2,3,4,5)),
  Wet_road          = as.integer(num(df$fondo_stradale) != 1)
)
F[is.na(F)] <- 0L

# Outcomes 
df$mc    <- as.integer(df$multi_casualty)
df$fatal <- as.integer(df$fatal_accident)

cat("  Feature matrix:", nrow(F), "rows ×", ncol(F), "binary columns\n")
cat("  Outcome: mc =", sum(df$mc, na.rm=TRUE), "(",
    round(mean(df$mc, na.rm=TRUE)*100,1), "%),",
    "fatal =", sum(df$fatal, na.rm=TRUE), "\n\n")

# ── 2. Standardise---------------------------------------------------
cat("--- 2. StandardScaler (mean=0, sd=1) ---\n")
feature_sd <- vapply(F, stats::sd, numeric(1))
if (any(!is.finite(feature_sd) | feature_sd == 0)) {
  stop("Constant/invalid clustering features: ",
       paste(names(feature_sd)[!is.finite(feature_sd) | feature_sd == 0], collapse=", "))
}
X <- scale(F)
stopifnot(nrow(X) == nrow(df), ncol(X) == 16L, all(is.finite(X)))

# ── 3. Validity metrics —---------------------------------------------------
cat("--- 3. Cluster validity sweep k=2..10 ---\n")
set.seed(0)
smp_idx <- sample(nrow(X), 12000)
Xs      <- X[smp_idx, ]
# 12,000-point distance vector uses approximately 576 MB; compute only once.
dist_Xs <- stats::dist(Xs)

sweep_res <- data.frame(k=2:10, inertia=NA, silhouette=NA, davies_bouldin=NA)
for (i in seq_len(nrow(sweep_res))) {
  k_val <- sweep_res$k[i]
  km    <- kmeans(Xs, centers=k_val, nstart=4, iter.max=300)
  sil   <- cluster::silhouette(km$cluster, dist_Xs)
  db    <- tryCatch(
    clusterSim::index.DB(Xs, km$cluster)$DB,
    error = function(e) {
      # manual Davies-Bouldin if clusterSim unavailable
      k    <- k_val
      labs <- km$cluster
      cents<- km$centers
      s    <- sapply(1:k, function(ci) {
        pts <- Xs[labs==ci, , drop=FALSE]
        if (nrow(pts)<2) return(0)
        mean(sqrt(rowSums((pts - matrix(cents[ci,],nrow(pts),ncol(pts),byrow=TRUE))^2)))
      })
      rij <- outer(1:k,1:k, Vectorize(function(a,b) {
        if (a==b) return(0)
        (s[a]+s[b])/sqrt(sum((cents[a,]-cents[b,])^2))
      }))
      mean(apply(rij,1,max))
    })
  sweep_res$inertia[i]      <- round(km$tot.withinss)
  sweep_res$silhouette[i]   <- round(mean(sil[,3]),4)
  sweep_res$davies_bouldin[i]<- round(db,4)
  cat("  k=", k_val,"  inertia=", round(km$tot.withinss),
      "  sil=", round(mean(sil[,3]),3),
      "  DB=", round(db,3), "\n")
}

# ── 4. Stability (ARI)---------------------------------------------------
cat("\n--- 4. Stability ARI ---\n")
if (requireNamespace("mclust", quietly=TRUE)) {
  for (k_val in c(5,6,7,8)) {
    labs_list <- lapply(0:3, function(s) {
      set.seed(s); kmeans(X, centers=k_val, nstart=1, iter.max=200)$cluster
    })
    pairs <- combn(4,2)
    ari   <- mean(apply(pairs, 2, function(p)
      mclust::adjustedRandIndex(labs_list[[p[1]]], labs_list[[p[2]]])))
    cat("  k=", k_val, "  ARI=", round(ari,3), "\n")
  }
} else {
  cat("  (mclust not installed — skipping ARI; install.packages('mclust') to enable)\n")
}

# ── 5. Final model: k=6---------------------------------------------------
cat("\n--- 5. Final model k=6 (nstart=10, seed=42) ---\n")
set.seed(42)
km6 <- kmeans(X, centers=6, nstart=10, iter.max=500)
df$cluster <- as.factor(km6$cluster)
F$cluster  <- df$cluster
accidents_main$cluster <- df$cluster
accidents_ml_safe$cluster <- df$cluster
if (RUN_MODELLING && exists("ml_data")) ml_data$cluster <- df$cluster

sil6 <- cluster::silhouette(km6$cluster[smp_idx], dist_Xs)
cat("  k=6 silhouette:", round(mean(sil6[,3]),4), "\n\n")

# ── 6. Cluster profiles ---------------------------------------------------
cat("--- 6. Cluster profiles ---\n")
feat_cols <- names(F)[names(F) != "cluster"]
profiles  <- df %>%
  dplyr::bind_cols(F[, feat_cols]) %>%
  dplyr::group_by(cluster) %>%
  dplyr::summarise(
    n          = n(),
    pct        = round(n()/nrow(df)*100,1),
    mc_rate    = round(mean(mc,  na.rm=TRUE)*100,1),
    fatal_rate = round(mean(fatal,na.rm=TRUE)*100,2),
    across(all_of(feat_cols), ~round(mean(.x, na.rm=TRUE)*100,1)),
    .groups="drop"
  ) %>%
  dplyr::arrange(dplyr::desc(mc_rate))

# ── Verified cluster labels---------------------------------------------------
# Profiles audited using within-cluster prevalence and national-rate deviations.
CLUSTER_LABELS <- c(
  "1" = "Heavy-Vehicle Collisions",
  "2" = "Urban Pedestrian Crashes",
  "3" = "Rear/Side Multi-Vehicle Collisions",
  "4" = "Cyclist / E-Scooter Crashes",
  "5" = "Run-Off Road Crashes",
  "6" = "Frontal / Front-Lateral Collisions"
)

profile_idx <- setNames(seq_len(nrow(profiles)), as.character(profiles$cluster))
stopifnot(setequal(names(profile_idx), names(CLUSTER_LABELS)))
profile_signatures_match <- c(
  "1" = profiles$Heavy_vehicle[profile_idx[["1"]]] >= 95,
  "2" = profiles$Pedestrian_hit[profile_idx[["2"]]] >= 95,
  "3" = profiles$Rear_side[profile_idx[["3"]]] >= 60 &&
        profiles$X2nd_vehicle[profile_idx[["3"]]] >= 70,
  "4" = profiles$Cyclist_escooter[profile_idx[["4"]]] >= 95,
  "5" = profiles$Run_off[profile_idx[["5"]]] >= 95,
  "6" = profiles$Head_on[profile_idx[["6"]]] >= 95 &&
        profiles$X2nd_vehicle[profile_idx[["6"]]] >= 95
)
if (!all(profile_signatures_match)) {
  stop("Cluster label audit mismatch for IDs: ",
       paste(names(profile_signatures_match)[!profile_signatures_match], collapse = ", "),
       ". Re-audit the six profiles before using these presentation labels.")
}
profiles$name <- unname(CLUSTER_LABELS[as.character(profiles$cluster)])
stopifnot(!anyNA(profiles$name), !anyDuplicated(profiles$name))

# Keep the six existing, visually distinct cluster colours (C1 pink, C4 teal).
CL_COLS <- c(
  "Heavy-Vehicle Collisions"          = "#E98EA5",
  "Urban Pedestrian Crashes"          = "#636c76",
  "Rear/Side Multi-Vehicle Collisions"= "#2563A6",
  "Cyclist / E-Scooter Crashes"       = "#00BDB5",
  "Run-Off Road Crashes"              = "#8250DF",
  "Frontal / Front-Lateral Collisions"= "#C83C3C"
)
profiles$colour <- unname(CL_COLS[profiles$name])

cluster_slide_table <- profiles %>%
  dplyr::mutate(
    defining_characteristics = dplyr::case_when(
      as.character(cluster) == "1" ~ sprintf(
        "Heavy vehicle A %.0f%%; second vehicle %.1f%%; motorway %.1f%%",
        Heavy_vehicle, X2nd_vehicle, Motorway),
      as.character(cluster) == "2" ~ sprintf(
        "Pedestrian collision %.0f%%; urban %.1f%%; driver 65+ %.1f%%",
        Pedestrian_hit, Urban, Senior_65),
      as.character(cluster) == "3" ~ sprintf(
        "Rear/side collision %.1f%%; second vehicle %.1f%%; urban %.1f%%",
        Rear_side, X2nd_vehicle, Urban),
      as.character(cluster) == "4" ~ sprintf(
        "Cyclist/e-scooter vehicle A %.0f%%; urban %.1f%%; second vehicle %.1f%%",
        Cyclist_escooter, Urban, X2nd_vehicle),
      as.character(cluster) == "5" ~ sprintf(
        "Run-off road %.0f%%; wet road %.1f%%; night %.1f%%",
        Run_off, Wet_road, Night),
      as.character(cluster) == "6" ~ sprintf(
        "Frontal/front-lateral %.0f%%; second vehicle %.0f%%; intersection %.1f%%",
        Head_on, X2nd_vehicle, Intersection)
    )
  ) %>%
  dplyr::transmute(
    cluster_id = as.character(cluster),
    profile = name,
    accidents = n,
    share_of_accidents_pct = pct,
    multi_injury_rate_pct = mc_rate,
    fatal_accident_rate_pct = fatal_rate,
    defining_characteristics
  ) %>%
  dplyr::arrange(dplyr::desc(multi_injury_rate_pct))

cat("\n--- Verified six-cluster profile table (slide-ready) ---\n")
print(cluster_slide_table, n = Inf, width = Inf)


# ── 7. Plots ---------------------------------------------------
cat("\n--- 7. Generating plots ---\n")
plots <- list()

THEME <- theme_minimal(base_size=13) +
  theme(panel.grid.minor=element_blank(),
        plot.title=element_text(face="bold",size=14),
        plot.subtitle=element_text(size=11,colour="#636c76"),
        legend.position="bottom")

## 7a. Silhouette sweep
plots$silhouette_sweep <- ggplot(sweep_res, aes(x=k, y=silhouette)) +
  geom_line(colour="#2563A6", linewidth=1.2) +
  geom_point(aes(colour=ifelse(k==6,"chosen","other")), size=4, show.legend=FALSE) +
  scale_colour_manual(values=c(chosen="#E8821A", other="#2563A6")) +
  geom_text(aes(label=round(silhouette,3)), vjust=-0.8, size=3.2, colour="#444") +
  annotate("point", x=6, y=sweep_res$silhouette[sweep_res$k==6],
           colour="#E8821A", size=6, shape=8) +
  scale_x_continuous(breaks=2:10) +
  labs(title="Cluster validity — silhouette score (full 0–1 scale)",
       subtitle=sprintf("Scores: %.3f–%.3f | k = 6 selected", min(sweep_res$silhouette), max(sweep_res$silhouette)),
       x="Number of clusters (k)", y="Mean silhouette score (0–1 scale)") +
  coord_cartesian(ylim=c(0,1)) +
  THEME

## 7b. Elbow
plots$elbow <- ggplot(sweep_res, aes(x=k, y=inertia/1000)) +
  geom_line(colour="#1baf7a", linewidth=1.2) +
  geom_point(aes(colour=ifelse(k==6,"chosen","other")), size=4, show.legend=FALSE) +
  scale_colour_manual(values=c(chosen="#E8821A", other="#1baf7a")) +
  annotate("point", x=6, y=sweep_res$inertia[sweep_res$k==6]/1000,
           colour="#E8821A", size=6, shape=8) +
  scale_x_continuous(breaks=2:10) +
  labs(title="Elbow method — within-cluster inertia",
       subtitle="From 0 to show true scale · k=6 ★ chosen",
       x="Number of clusters (k)", y="Total within-cluster SS (thousands)") +
  coord_cartesian(ylim=c(0,NA)) +
  THEME

## 7c. MC rate by cluster
prof_plot <- profiles %>%
  dplyr::mutate(name=factor(name, levels=rev(profiles$name)))

plots$mc_by_cluster <- ggplot(prof_plot,
    aes(x=name, y=mc_rate, fill=name)) +
  geom_col(width=0.7, show.legend=FALSE) +
  geom_text(aes(label=paste0(mc_rate,"%")), hjust=-0.1, size=4, fontface="bold") +
  scale_fill_manual(values=CL_COLS) +
  coord_flip(ylim=c(0,45)) +
  labs(title="Multi-casualty rate by cluster",
       subtitle="Outcomes computed AFTER clustering — never used as features",
       x=NULL, y="Multi-casualty rate (%)") +
  THEME

## 7d. Fatal rate by cluster
plots$fatal_by_cluster <- ggplot(prof_plot,
    aes(x=name, y=fatal_rate, fill=name)) +
  geom_col(width=0.7, show.legend=FALSE) +
  geom_text(aes(label=paste0(fatal_rate,"%")), hjust=-0.1, size=4, fontface="bold") +
  scale_fill_manual(values=CL_COLS) +
  coord_flip(ylim=c(0,5)) +
  labs(title="Fatal accident rate by cluster",
       subtitle="Fatal outcomes reported separately from clustering features",
       x=NULL, y="Fatal rate (%)") +
  THEME

## 7e. Feature heatmap
feat_long <- profiles %>%
  dplyr::select(name, all_of(feat_cols)) %>%
  tidyr::pivot_longer(-name, names_to="feature", values_to="pct") %>%
  dplyr::mutate(
    feature = gsub("_"," ", gsub("X2nd","2nd", feature)),
    name    = factor(name, levels=profiles$name)
  )

plots$feature_heatmap <- ggplot(feat_long, aes(x=feature, y=name, fill=pct)) +
  geom_tile(colour="white", linewidth=0.5) +
  geom_text(aes(label=ifelse(pct>=20, paste0(round(pct),"%"), "")),
            size=3, colour="white", fontface="bold") +
  scale_fill_gradient(low="#EBF5FB", high="#1A5276",
                      name="% of cluster", limits=c(0,100)) +
  scale_x_discrete(guide=guide_axis(angle=40)) +
  labs(title="Feature profile heatmap — k=6 clusters",
       subtitle="Percentage of accidents in each cluster with this feature = 1",
       x=NULL, y=NULL) +
  theme_minimal(base_size=11) +
  theme(plot.title=element_text(face="bold"), legend.position="right")

## 7f. PCA scatter
pca_res <- prcomp(X, center=FALSE, scale.=FALSE)
ev      <- round(pca_res$sdev^2 / sum(pca_res$sdev^2) * 100, 1)
pca_df  <- data.frame(
  PC1     = pca_res$x[,1],
  PC2     = pca_res$x[,2],
  cluster = df$cluster
)
# Map cluster IDs to names
cl_name_map <- setNames(profiles$name, profiles$cluster)
pca_df$cluster_name <- cl_name_map[as.character(pca_df$cluster)]

# Sample 200 per cluster for plotting
set.seed(42)
pca_samp <- pca_df %>%
  dplyr::group_by(cluster_name) %>%
  dplyr::slice_sample(n=200) %>%
  dplyr::ungroup()

# Centroids
cents <- pca_df %>%
  dplyr::group_by(cluster_name) %>%
  dplyr::summarise(PC1=mean(PC1), PC2=mean(PC2), .groups="drop")

plots$pca_scatter <- ggplot(pca_samp, aes(x=PC1, y=PC2, colour=cluster_name)) +
  geom_point(alpha=0.35, size=1.8, show.legend=FALSE) +
  geom_point(data=cents, aes(x=PC1, y=PC2, fill=cluster_name),
             shape=21, size=7, colour="white", stroke=2, show.legend=TRUE) +
  ggrepel::geom_text_repel(data=cents, aes(x=PC1, y=PC2, label=cluster_name),
                            size=3.5, fontface="bold", show.legend=FALSE,
                            box.padding=0.5, max.overlaps=20) +
  scale_colour_manual(values=CL_COLS) +
  scale_fill_manual(values=CL_COLS, name="Cluster") +
  labs(
    title="PCA scatter — k=6 clusters (2024 data, 16 binary features)",
    subtitle=paste0("PC1 (", ev[1], "%) × PC2 (", ev[2], "%) — ",
                    round(ev[1]+ev[2],1), "% of variance shown · ",
                    "2D projection; interpret separation together with cluster profiles"),
    x=paste0("PC1 — ", ev[1], "% variance"),
    y=paste0("PC2 — ", ev[2], "% variance")
  ) +
  THEME + theme(legend.position="right")

# ── PCA INTERPRETATION: ORIGINAL 16 STANDARDISED BINARY FEATURES ---------------------------------------
loadings_pc12 <- data.frame(
  feature = rownames(pca_res$rotation),
  PC1 = pca_res$rotation[, 1], PC2 = pca_res$rotation[, 2],
  row.names = NULL
)
pca_feature_contributions <- loadings_pc12 %>%
  tidyr::pivot_longer(cols = c(PC1, PC2),
                      names_to = "component", values_to = "loading") %>%
  group_by(component) %>%
  mutate(contribution_pct = 100 * loading^2 / sum(loading^2),
         loading_direction = ifelse(loading >= 0, "Positive", "Negative")) %>%
  ungroup()
stopifnot(all(abs(pca_feature_contributions %>% group_by(component) %>%
                    summarise(total = sum(contribution_pct), .groups = "drop") %>%
                    pull(total) - 100) < 1e-8))

pc_contribution_panel <- function(pc_name, variance_pct) {
  dat <- pca_feature_contributions %>%
    filter(component == pc_name) %>%
    arrange(desc(contribution_pct)) %>%
    slice_head(n = 8) %>%
    mutate(feature_label = gsub("_", " ", feature, fixed = TRUE),
           feature_label = factor(feature_label, levels = rev(feature_label)))
  ggplot(dat, aes(x = feature_label, y = contribution_pct,
                  fill = loading_direction)) +
    geom_col(width = 0.68, show.legend = FALSE) +
    geom_text(aes(label = sprintf("%.1f%%", contribution_pct)),
              hjust = -0.12, size = 3.6, fontface = "bold") +
    coord_flip(clip = "off") +
    scale_fill_manual(values = c(Positive = "#2865A3", Negative = "#ED8619")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.24))) +
    labs(title = paste0(pc_name, " — ", round(variance_pct, 1), "% variance explained"),
         x = NULL, y = "Feature contribution to component (%)") +
    theme_minimal(base_size = 11) +
    theme(panel.grid.major.y = element_blank(),
          plot.title = element_text(face = "bold"))
}
plots$pca_feature_contributions <-
  (pc_contribution_panel("PC1", ev[1]) |
   pc_contribution_panel("PC2", ev[2])) +
  patchwork::plot_annotation(
    title = "Which Features Define the First Two Principal Components?",
    subtitle = "Top 8 features by squared loading (share %) | Blue = positive, orange = negative loading",
    caption = "Contribution percentages sum to 100% across all 16 features for each PC, not across only those shown.")

# Optional companion: position/direction of ALL features in PC1–PC2 loading space.
plots$pca_feature_loading_map <- ggplot(loadings_pc12, aes(x = PC1, y = PC2)) +
  geom_hline(yintercept = 0, colour = "grey75", linetype = "dashed") +
  geom_vline(xintercept = 0, colour = "grey75", linetype = "dashed") +
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2),
               linewidth = 0.55, colour = "#2865A3",
               arrow = grid::arrow(length = grid::unit(0.11, "inches"))) +
  ggrepel::geom_text_repel(aes(label = gsub("_", " ", feature, fixed = TRUE)),
                           size = 3.25, max.overlaps = Inf) +
  coord_equal() +
  labs(title = "Feature Loadings on PC1 and PC2",
       subtitle = "Arrow direction and length represent each standardised feature's loadings",
       x = sprintf("PC1 (%.1f%% explained variance)", ev[1]),
       y = sprintf("PC2 (%.1f%% explained variance)", ev[2])) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

## 7g. Radar chart (fmsb) — one per cluster
make_radar <- function(prof_row, col_hex) {
  feats_to_show <- c("X2nd_vehicle","Head_on","Rear_side","Pedestrian_hit",
                     "Run_off","Night","Senior_65","Motorcycle","Motorway","Heavy_vehicle")
  feat_labels   <- c("2nd veh","Head-on","Rear/side","Ped","Run-off",
                     "Night","Senior","Moto","Motorway","Heavy")
  vals <- as.numeric(prof_row[, feats_to_show])
  radar_df <- rbind(rep(100,length(vals)), rep(0,length(vals)), vals)
  colnames(radar_df) <- feat_labels
  radar_df <- as.data.frame(radar_df)
  col_fill  <- paste0(col_hex, "40")
  fmsb::radarchart(radar_df, axistype=1, pcol=col_hex, pfcol=col_fill,
                   plwd=2.5, cglcol="#DDDDDD", cglty=1, axislabcol="#888",
                   vlcex=0.85, calcex=0.7, vlabels=feat_labels,
                   title=paste0(prof_row$name,
                                "\nMC=",prof_row$mc_rate,"% · fatal=",prof_row$fatal_rate,
                                "% · n=",format(prof_row$n,big.mark=",")))
}

## 7h. Cluster size bubble
plots$cluster_size <- ggplot(profiles,
    aes(x=reorder(name,-n), y=mc_rate, size=n, fill=name)) +
  geom_point(shape=21, colour="white", stroke=1.5, alpha=0.9) +
  geom_text(aes(label=paste0(pct,"%")), vjust=-1.6, size=3.5, fontface="bold") +
  scale_size_continuous(range=c(6,22), guide="none") +
  scale_fill_manual(values=CL_COLS, guide="none") +
  scale_x_discrete(guide=guide_axis(angle=20)) +
  labs(title="Cluster size (bubble area) vs multi-casualty rate",
       subtitle="Label = share of all 2024 accidents",
       x=NULL, y="Multi-casualty rate (%)") +
  THEME

# ── 8. Output all plots ---------------------------------------------------
cat("\n--- 8. Printing plots ---\n")

if (SHOW_FINAL_PLOTS) {
  print(plots$silhouette_sweep)
  print(plots$mc_by_cluster)
  print(plots$pca_scatter)
  print(plots$pca_feature_contributions)
  print(plots$pca_feature_loading_map)
}

# Radar charts are exported together and displayed LAST below, so they are
# not hidden by the summary plot on the RStudio Plots tab.

# ── 9. Patchwork summary ─────────────────────────────────────────────────────
cat("\n--- 9. Summary patchwork ---\n")
summary_plot <- (plots$silhouette_sweep | plots$elbow) /
                (plots$mc_by_cluster   | plots$fatal_by_cluster) /
                (plots$pca_scatter     | plots$cluster_size) +
  patchwork::plot_annotation(
    title    = "ISTAT 2024 — K-Means Clustering (k=6)",
    subtitle = "16 binary features · no injury data used during clustering · outcomes shown for validation only",
    theme    = theme(plot.title    = element_text(face="bold", size=16),
                     plot.subtitle = element_text(size=12, colour="#636c76"))
  )
# Optional: print(summary_plot) for a consolidated diagnostics dashboard.

# ── 10. Console summary ───────────────────────────────────────────────────────
cat("\n============================================================\n")
cat("  CLUSTERING SUMMARY\n")
cat("============================================================\n")
cat("  Algorithm:   K-Means\n")
cat("  k chosen:    6\n")
cat("  Features:    16 binary (no outcome data)\n")
cat("  Data:        2024 ISTAT microdata (n =", nrow(df), ")\n")
cat("  Silhouette: ", round(mean(sil6[,3]), 3), "\n\n", sep="")
cat("  Cluster profiles:\n")
cat(sprintf("  %-30s %6s %7s %7s\n","Name","Share","MC%","Fatal%"))
cat("  ", strrep("-",55), "\n", sep="")
for (i in seq_len(nrow(profiles))) {
  cat(sprintf("  %-30s %5.1f%% %6.1f%% %6.2f%%\n",
      profiles$name[i], profiles$pct[i],
      profiles$mc_rate[i], profiles$fatal_rate[i]))
}
cat("\n  PCA is a two-dimensional projection; check profiles for distinctions\n")
cat("============================================================\n")


# Compatibility alias for project workflows that previously used plots_cl.
plots_cl <- plots
cat("\nK=6 plots are available as plots$... or plots_cl$...\n")

# Draw ALL SIX radar profiles on one device/page, with a deliberate layout.
# Base graphics cannot be stored in a ggplot object; use this function or the PNG.
draw_six_radars <- function() {
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(mfrow = c(2, 3), mar = c(2.7, 2.2, 4.1, 2.2),
                oma = c(0.5, 0.5, 2.0, 0.5))
  for (i in seq_len(nrow(profiles))) {
    make_radar(profiles[i, , drop = FALSE], profiles$colour[i])
  }
  graphics::mtext("ISTAT 2024 — Six K-means Cluster Profiles",
                  side = 3, outer = TRUE, line = 0.7, font = 2, cex = 1.15)
  invisible(NULL)
}

if (EXPORT_FINAL_PLOTS) {
  cluster_output_dir <- "outputs/presentation_clusters_2024"
  dir.create(cluster_output_dir, recursive = TRUE, showWarnings = FALSE)
  cluster_to_save <- list(
    silhouette_k2_k10 = plots$silhouette_sweep,
    multi_injury_by_cluster = plots$mc_by_cluster,
    pca_cluster_scatter = plots$pca_scatter,
    pca_pc1_pc2_feature_contributions = plots$pca_feature_contributions,
    pca_pc1_pc2_feature_loading_map = plots$pca_feature_loading_map
  )
  for (nm in names(cluster_to_save)) {
    ggplot2::ggsave(
      filename = file.path(cluster_output_dir, paste0(nm, ".png")),
      plot = cluster_to_save[[nm]],
      width = if (nm == "pca_pc1_pc2_feature_contributions") 14 else 11,
      height = if (nm == "pca_pc1_pc2_feature_contributions") 7 else 7,
      dpi = 300, bg = "white"
    )
  }
  grDevices::png(filename = file.path(cluster_output_dir, "six_radar_profiles.png"),
                 width = 4200, height = 2750, res = 300, bg = "white")
  tryCatch(draw_six_radars(), finally = grDevices::dev.off())
  utils::write.csv(as.data.frame(profiles),
                   file.path(cluster_output_dir, "cluster_profiles.csv"),
                   row.names = FALSE)
  utils::write.csv(as.data.frame(cluster_slide_table),
                   file.path(cluster_output_dir, "cluster_profiles_for_slides.csv"),
                   row.names = FALSE)
  utils::write.csv(as.data.frame(sweep_res),
                   file.path(cluster_output_dir, "cluster_validity.csv"),
                   row.names = FALSE)
  utils::write.csv(as.data.frame(pca_feature_contributions),
                   file.path(cluster_output_dir, "pca_feature_contributions.csv"),
                   row.names = FALSE)
  cat("\nSaved six radar plots and PCA figures in ", cluster_output_dir, "\n", sep = "")
}

# Show the radar composite LAST so it remains visible in RStudio's Plots tab.
# When using Rscript (non-interactive), the PNG above is the reliable output.
if (SHOW_FINAL_PLOTS && interactive()) {
  draw_six_radars()
}
cat("\nTo revisit the plots: print(plots$pca_scatter); ",
    "print(plots$pca_feature_contributions); draw_six_radars()\n", sep = "")

}
