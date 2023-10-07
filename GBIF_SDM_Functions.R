# ---- Download Bioclimatic data from Worldclim ----
bioclim_download = function(resolution = 10) {
  # Create the "data" directory if it doesn't exist
  data_dir = here::here("data2")
  if (!dir.exists(data_dir)) {
    dir.create(data_dir, recursive = TRUE)
  }
  
  worldclim_global(var = "bio",
                   res = resolution,
                   path = data_dir)
  
  cat("Bioclimatic data downloaded and saved in:", data_dir, "\n")
}

# ---- Automated Batch GBIF Download ----
gbif_occ_data = function(scientific_name, countries) {
  # Create a semicolon-separated string of country codes
  country_codes = paste(unique(unlist(strsplit(countries, ";"))), collapse = ";")
  
  sp = occ_data(scientificName = scientific_name, 
                 country = country_codes,
                 limit = 100000) # Need to set limit or it defaults to 500
  
  sp_df = data.frame(sp$data)
  
  # Make sure the data frame has the required columns
  if (!all(c("decimalLatitude", "decimalLongitude", "countryCode") %in% names(sp_df))) {
    return(data.frame(Latitude = NA, Longitude = NA, Country = NA))
  }
  
  # Filter out NA points & Rename Columns
  occ = na.omit(sp_df[, c("decimalLatitude", "decimalLongitude", "countryCode")]) 
  names(occ) = c("Latitude", "Longitude", "Country")
  return(occ)
}

# ---- Run MaxEnt Batch ----
run_maxent = function(occ_sf, env) {
  # Get unique Scientific Names
  print(class(occ_sf$Scientific.Name))
  species = unique(occ_sf$Scientific.Name)
  
  # Initialize lists to store results
  enm_results = list()
  maxent_models = list()
  
  # Loop through each species
  for(sp in species) {
    cat(paste("Processing", sp, "\n"))
    
    # Subset occurrence data for the species
    occ_sp = occ_sf[occ_sf$Scientific.Name == sp, c("Longitude", "Latitude")]
    
    # Run ENMevaluate
    enmeval_results = ENMevaluate(occ_sp, env, 
                                   bg = NULL, 
                                   tune.args = list(fc = c("L","LQ","H", "LQH", "LQHP", "LQHPT"), rm = 1:5), 
                                   partitions = "randomkfold", partition.settings = list(kfolds = 2), 
                                   algorithm = "maxnet", 
                                   taxon.name = sp)  # specify the taxon name here
    
    enmeval_df = enmeval_results@results
    
    # Store ENMevaluate results
    enm_results[[sp]] = enmeval_df
    
    # Subset the ENMeval results to get the best model
    enmeval_bestm = subset(enmeval_df, delta.AICc == 0)
    
    # Decode the features
    maxent_feats = as.character(enmeval_bestm$fc)
    maxent_rm = as.character(enmeval_bestm$rm)
    
    # Print out the results
    cat(paste("Best features for", sp, ":", maxent_feats, "\n"))
    cat(paste("Best regularization multiplier for", sp, ":", maxent_rm, "\n"))
    
    # Run MaXent SDM
    sp_maxent_model = dismo::maxent(env, as.matrix(occ_sp), features = maxent_feats, betamultiplier = maxent_rm)
    
    # Store maxent model
    maxent_models[[sp]] = sp_maxent_model
  }
  
  return(list(enm = enm_results, models = maxent_models))
}

# ---- Replace General Regions w/ List of Countries ----
region_to_country = function(location_string, mapping) {
  for(region in names(mapping)) {
    if(grepl(region, location_string, ignore.case = TRUE)) {
      location_string = gsub(region, paste(mapping[[region]], collapse = "; "), location_string, ignore.case = TRUE)
    }
  }
  return(location_string)
}

# ---- Convert Source Location to countryCode ----
flexible_country_code = function(location_string) {
  # List of common two-word countries
  two_word_countries = c("Antigua and Barbuda", "Bosnia and Herzegovina", "Central African Republic", 
                          "Costa Rica", "Czech Republic", "Dominican Republic", "Timor-Leste", "El Salvador",
                          "Equatorial Guinea", "Marshall Islands", "New Zealand", "North Korea", 
                          "Papua New Guinea", "San Marino", "Sao Tome and Principe", "Sierra Leone", 
                          "Saudi Arabia", "Solomon Islands", "South Korea", "Sri Lanka", 
                          "Trinidad and Tobago", "United Arab Emirates", "United Kingdom", "United States")
  
  # Match two-word countries and store their codes
  matched_two_word_codes = sapply(two_word_countries, function(country) {
    if(str_detect(location_string, regex(country, ignore_case = TRUE))) {
      countrycode(country, "country.name", "iso2c", nomatch = NA)
    } else {
      NA
    }
  })
  
  # Remove matched two-word countries from the location string
  for(country in two_word_countries) {
    location_string = str_remove(location_string, regex(country, ignore_case = TRUE))
  }
  
  # Split remaining string by non-word sequences
  potential_countries = unlist(strsplit(location_string, split = "\\W+"))
  
  # Convert each potential country name to its country code
  country_codes = sapply(potential_countries, function(country) {
    country_trimmed = str_trim(country)
    countrycode(country_trimmed, "country.name", "iso2c", nomatch = NA)
  })
  
  # Combine codes and filter out NA values
  all_codes = c(matched_two_word_codes, country_codes)
  valid_codes = paste(all_codes[!is.na(all_codes)], collapse = ";")
  
  return(valid_codes)
}

# ---- Region Countries Repository ----
# Define a mapping between region names and the countries they contain
region_mapping = list(
  
  # Asia
  
  'Asia' = c('Afghanistan', 'Armenia', 'Azerbaijan', 'Bahrain', 'Bangladesh', 'Bhutan', 'Brunei', 'Cambodia', 'China', 'Cyprus', 'Georgia', 'India', 'Indonesia', 'Iran', 'Iraq', 'Israel', 'Japan', 'Jordan', 'Kazakhstan', 'Kuwait', 'Kyrgyzstan', 'Laos', 'Lebanon', 'Malaysia', 'Maldives', 'Mongolia', 'Myanmar (Burma)', 'Nepal', 'North Korea', 'Oman', 'Pakistan', 'Palestine', 'Philippines', 'Qatar', 'Russia', 'Saudi Arabia', 'Singapore', 'South Korea', 'Sri Lanka', 'Syria', 'Taiwan', 'Tajikistan', 'Thailand', 'Timor-Leste', 'Turkey', 'Turkmenistan', 'United Arab Emirates', 'Uzbekistan', 'Vietnam', 'Yemen'),
  
  'East_Asia' = c('China', 'Japan', 'Mongolia', 'North Korea', 'South Korea', 'Taiwan'),
  
  'West_Asia (Middle East)' = c('Bahrain', 'Iraq', 'Iran', 'Israel', 'Jordan', 'Kuwait', 
                                'Lebanon', 'Oman', 'Qatar', 'Saudi Arabia', 'Syria', 'United Arab Emirates', 
                                'Yemen'),
  
  'SE_Asia' = c('Indonesia', 'Philippines', 'Vietnam', 'Thailand', 'Myanmar (Burma)', 
                'Malaysia', 'Cambodia', 'Laos', 'Singapore', 'Timor-Leste', 
                'Brunei Darussalam'),
  'SW_Asia' = c('Armenia', 'Azerbaijan', 'Bahrain', 'Cyprus', 'Georgia', 'Iran', 'Iraq', 'Israel', 'Jordan', 'Kuwait', 'Lebanon', 'Oman', 'Palestine', 'Qatar', 'Saudi Arabia', 'Syria', 'Turkey', 'United Arab Emirates', 'Yemen'),
  
  'NE_Asia' = c('China', 'Japan', 'Mongolia', 'North Korea', 'South Korea'),
  
  'NW_Asia' = c('Russia', 'Kazakhstan'),
  
  'Central_Asia' = c('Kazakhstan', 'Kyrgyzstan', 'Tajikistan', 'Turkmenistan', 'Uzbekistan'),
  
  'South_Asia' = c('India', 'Pakistan', 'Bangladesh', 'Sri Lanka', 'Nepal', 'Bhutan', 'Maldives'),
  
  'Tropical_Asia' = c(
    "Bangladesh", "Brunei", "Cambodia", "Timor-Leste", "India",
    "Indonesia", "Laos",
    "Malaysia", "Maldives", "Myanmar (Burma)", "Papua New Guinea (part of it is in Southeast Asia)",
    "Philippines", "Singapore", "Sri Lanka", "Thailand", "Vietnam"),
  
  'Eurasia' = c('Albania', 'Andorra', 'Armenia', 'Austria', 'Azerbaijan', 'Belarus', 'Belgium', 
                'Bosnia and Herzegovina', 'Bulgaria', 'Croatia', 'Cyprus', 'Czech Republic', 
                'Denmark', 'Estonia', 'Finland', 'France', 'Georgia', 'Germany', 'Greece', 
                'Hungary', 'Iceland', 'Ireland', 'Italy', 'Kazakhstan', 'Kosovo', 'Latvia', 
                'Liechtenstein', 'Lithuania', 'Luxembourg', 'Macedonia', 
                'Malta', 'Moldova', 'Monaco', 'Montenegro', 'Netherlands', 'Norway', 'Poland', 
                'Portugal', 'Romania', 'Russia', 'San Marino', 'Serbia', 'Slovakia', 'Slovenia', 
                'Spain', 'Sweden', 'Switzerland', 'Turkey', 'Ukraine', 'United Kingdom', 'Vatican City', 
                'Kazakhstan', 'Kyrgyzstan', 'Tajikistan', 'Turkmenistan', 'Uzbekistan',
                'China', 'Japan', 'Mongolia', 'North Korea', 'South Korea', 'Taiwan',
                'Afghanistan', 'Bangladesh', 'Bhutan', 'India', 'Maldives', 'Nepal', 'Pakistan', 'Sri Lanka',
                'Bahrain', 'Iraq', 'Iran', 'Israel', 'Jordan', 'Kuwait', 
                'Lebanon', 'Oman', 'Qatar', 'Saudi Arabia', 'Syria', 'United Arab Emirates', 'Yemen'),
  
  # Europe
  
  'Europe' = c('Albania', 'Andorra', 'Armenia', 'Austria', 'Azerbaijan', 'Belarus', 'Belgium', 'Bosnia and Herzegovina', 'Bulgaria', 'Croatia', 'Cyprus', 'Czech Republic', 'Denmark', 'Estonia', 'Finland', 'France', 'Georgia', 'Germany', 'Greece', 'Hungary', 'Iceland', 'Ireland', 'Italy', 'Kazakhstan', 'Latvia', 'Liechtenstein', 'Lithuania', 'Luxembourg', 'Malta', 'Moldova', 'Monaco', 'Montenegro', 'Netherlands', 'North Macedonia', 'Norway', 'Poland', 'Portugal', 'Romania', 'Russia', 'San Marino', 'Serbia', 'Slovakia', 'Slovenia', 'Spain', 'Sweden', 'Switzerland', 'Turkey', 'Ukraine', 'United Kingdom', 'Vatican City'),
  
  'Western_Europe' = c('Austria', 'Belgium', 'France', 'Germany', 'Liechtenstein', 'Luxembourg', 'Monaco', 'Netherlands', 'Switzerland'),
  
  'Eastern_Europe' = c('Bulgaria', 'Czech Republic', 'Hungary', 'Poland', 'Romania', 'Russia', 'Slovakia'),
  
  'Northern_Europe' = c('Denmark', 'Estonia', 'Finland', 'Iceland', 'Ireland', 'Latvia', 'Lithuania', 'Norway', 'Sweden', 'United Kingdom'),
  
  'Southern_Europe' = c('Albania', 'Andorra', 'Bosnia and Herzegovina', 'Croatia', 'Greece', 'Italy', 'Malta', 'Montenegro', 'North Macedonia', 'Portugal', 'San Marino', 'Serbia', 'Slovenia', 'Spain', 'Vatican City'),
  
  'NW_Europe' = c(
    "Belgium", "Channel Islands (Crown Dependencies of the United Kingdom)", "France",
    "Ireland", "Isle of Man (a Crown Dependency of the United Kingdom)", "Luxembourg", "Monaco",
    "Netherlands", "United Kingdom"
  ),
  
  'SW_Europe' = c("Andorra", "Gibraltar (a British Overseas Territory)", "Portugal", "Spain"),
  
  'SE_Europe' = c(
    "Albania", "Bosnia and Herzegovina", "Bulgaria", "Croatia", "Greece", "Kosovo", "Montenegro",
    "North Macedonia", "Romania", "Serbia", "Slovenia", "Turkey"
  ),
  
  'NE_Europe' = c(
    "Belarus", "Estonia", "Latvia", "Lithuania", "Russia"
  ),
  
  # Middle East
  
  'Middle_East' = c('Bahrain', 'Cyprus', 'Egypt', 'Iran', 'Iraq', 'Israel', 'Jordan', 'Kuwait', 'Lebanon', 'Oman', 'Palestine', 'Qatar', 'Saudi Arabia', 'Syria', 'Turkey', 'United Arab Emirates', 'Yemen'),
  
  'Mediterranean' = c(
    "Albania", "Algeria", "Bosnia and Herzegovina", "Croatia", "Cyprus",
    "Egypt", "France", "Greece", "Israel", "Italy", "Lebanon", "Libya",
    "Malta", "Monaco", "Montenegro", "Morocco", "North Macedonia", "Palestine", 
    "Slovenia", "Spain", "Syria", "Tunisia", "Turkey"),
  
  
  # Americas
  
  'Americas' = c('Antigua and Barbuda', 'Argentina', 'Bahamas', 'Barbados', 'Belize', 'Bolivia', 'Brazil', 'Canada', 'Chile', 'Colombia', 'Costa Rica', 'Cuba', 'Dominica', 'Dominican Republic', 'Ecuador', 'El Salvador', 'Grenada', 'Guatemala', 'Guyana', 'Haiti', 'Honduras', 'Jamaica', 'Mexico', 'Nicaragua', 'Panama', 'Paraguay', 'Peru', 'Saint Kitts and Nevis', 'Saint Lucia', 'Saint Vincent and the Grenadines', 'Suriname', 'Trinidad and Tobago', 'United States', 'Uruguay', 'Venezuela'),
  
  'North_America' = c('Canada', 'United States', 'Mexico'),
  
  'Central_America' = c('Belize', 'Costa Rica', 'El Salvador', 'Guatemala', 'Honduras', 'Nicaragua', 'Panama'),
  
  'South_America' = c('Argentina', 'Bolivia', 'Brazil', 'Chile', 'Colombia', 'Ecuador', 'Guyana', 'Paraguay', 'Peru', 'Suriname', 'Uruguay', 'Venezuela'),
  
  'Caribbean' = c(
    "Antigua and Barbuda", "Bahamas", "Barbados", "Cayman Islands ",
    "Cuba", "Dominica", "Dominican Republic", "Grenada", "Guadeloupe",
    "Haiti", "Jamaica", "Martinique", "Montserrat",
    "Puerto Rico", "Saint Kitts and Nevis", "Saint Lucia",
    "Saint Vincent", "Trinidad", "Tobago", "Turks", "Caicos Islands",
    "Virgin Islands"),
  
  'Caroline_Islands' = c(
    "Federated States of Micronesia",
    "Palau",
    "Marshall Islands",
    "Papua New Guinea"),
  
  # Africa
  
  'Africa' = c(
    "Algeria", "Angola", "Benin", "Botswana", "Burkina Faso", "Burundi",
    "Cabo Verde", "Cameroon", "Central African Republic", "Chad", "Comoros",
    "Congo (Brazzaville)", "Congo (Kinshasa)", "Cote d'Ivoire", "Djibouti",
    "Egypt", "Equatorial Guinea", "Eritrea", "Eswatini", "Ethiopia", "Gabon",
    "Gambia", "Ghana", "Guinea", "Guinea-Bissau", "Kenya", "Lesotho", "Liberia",
    "Libya", "Madagascar", "Malawi", "Mali", "Mauritania", "Mauritius", "Morocco",
    "Mozambique", "Namibia", "Niger", "Nigeria", "Rwanda", "Sao Tome and Principe",
    "Senegal", "Seychelles", "Sierra Leone", "Somalia", "Sudan",
    "Sudan", "Tanzania", "Togo", "Tunisia", "Uganda", "Zambia", "Zimbabwe"),
  
  
  'Sub-Saharan_Africa' = c('Angola', 'Benin', 'Botswana', 'Burkina Faso', 'Burundi', 'Cape Verde', 'Cameroon', 'Central African Republic', 'Chad', 'Comoros', 'Congo (DRC)', 'Republic of Congo', 'Djibouti', 'Equatorial Guinea', 'Eritrea', 'Eswatini', 'Ethiopia', 'Gabon', 'Gambia', 'Ghana', 'Guinea', 'Guinea-Bissau', 'Ivory Coast', 'Kenya', 'Lesotho', 'Liberia', 'Madagascar', 'Malawi', 'Mali', 'Mauritius', 'Mozambique', 'Namibia', 'Niger', 'Nigeria', 'Rwanda', 'São Tomé and Príncipe', 'Senegal', 'Seychelles', 'Sierra Leone', 'Somalia', 'Sudan', 'Sudan', 'Tanzania', 'Togo', 'Uganda', 'Zambia', 'Zimbabwe'),
  
  'North_Africa' = c(
    "Algeria", "Egypt", "Libya", "Mauritania", "Morocco", "Sudan", "Tunisia", "Western Sahara"
  ),
  'South_Africa' = c(
    "Angola", "Botswana", "Eswatini", "Lesotho", "Malawi", "Mozambique", 
    "Namibia", "Zambia", "Zimbabwe"
  ),
  'East_Africa' = c(
    "Burundi", "Comoros", "Djibouti", "Eritrea", "Ethiopia", "Kenya", "Madagascar",
    "Mauritius", "Rwanda", "Seychelles", "Somalia", "Sudan", "Tanzania", "Uganda"
  ),
  'West_Africa' = c(
    "Benin", "Burkina Faso", "Cote d'Ivoire", "Gambia", "Ghana", "Guinea", "Guinea-Bissau",
    "Liberia", "Mali", "Niger", "Nigeria", "Senegal", "Sierra Leone", "Togo"),
  
  'Central_Africa' = c(
    "Angola", "Burundi", "Cameroon", "Central African Republic", "Chad", "Congo (Brazzaville)",
    "Congo (Kinshasa)", "Equatorial Guinea", "Gabon", "Rwanda", "Sao Tome and Principe"),
  
  'NW_Africa' = c(
    "Algeria", "Egypt", "Libya", "Mauritania", "Morocco", "Sudan", "Tunisia", "Western Sahara"),
  
  'SW_Africa' = c(
    "Angola", "Namibia", "Botswana", "Eswatini", "Lesotho"),
  
  'SE_Africa' = c(
    "Kenya", "Uganda", "Tanzania", "Rwanda", "Burundi", "Comoros", "Djibouti", 
    "Eritrea", "Ethiopia", "Madagascar", "Malawi", "Mauritius", "Mozambique", "Seychelles", 
    "Somalia", "South Sudan", "Zambia", "Zimbabwe"),
  
  'NE_Africa' = c(
    "Benin", "Burkina Faso", "Cote d'Ivoire", "Gambia", "Ghana", "Guinea",
    "Guinea-Bissau", "Liberia", "Mali", "Niger", "Nigeria", "Senegal", "Sierra Leone", "Togo", 
    "Cameroon", "Central African Republic", "Chad", "Congo (Brazzaville)", "Congo (Kinshasa)", 
    "Equatorial Guinea", "Gabon", "Sao Tome and Principe"),
  
  # Biogeographic Realms
  
  'Nearctic' = c("Canada", "United States", "Greenland", "Mexico"),
  
  'Neotropical' = c(
    "Argentina", "Belize", "Bolivia", "Brazil", "Chile", "Colombia", "Costa Rica",
    "Cuba", "Dominican Republic", "Ecuador", "El Salvador", "French Guiana",
    "Guatemala", "Guyana", "Honduras", "Jamaica", "Mexico", "Nicaragua", "Panama", "Paraguay", 
    "Peru", "Suriname", "Trinidad and Tobago", "Uruguay", "Venezuela"),
  
  'Palearctic' = c(
    "Afghanistan", "Albania", "Algeria", "Andorra", "Armenia", "Austria", 
    "Azerbaijan", "Bahrain", "Belarus", "Belgium", "Bhutan", 
    "Bosnia and Herzegovina", "Bulgaria", "China",
    "Croatia", "Cyprus", "Czech Republic", "Denmark", "Egypt",
    "Estonia", "Faroe Islands ", "Finland", "France", "Georgia", "Germany",
    "Greece", "Greenland", "Hungary", "Iceland", "India",
    "Iran", "Iraq", "Ireland", "Israel", "Italy", "Jordan", "Kazakhstan",
    "Kuwait", "Kyrgyzstan", "Latvia", "Lebanon", "Libya",
    "Liechtenstein", "Lithuania", "Luxembourg", "Macedonia",
    "Malta", "Moldova", "Monaco", "Mongolia", "Montenegro",
    "Morocco", "Nepal", "Netherlands",
    "Norway", "Oman (northern and central regions)", "Pakistan",
    "Poland", "Portugal", "Qatar", "Romania", "Russia",
    "Saudi Arabia", "Serbia", "Slovakia", "Slovenia",
    "Spain", "Sri Lanka", "Sweden", "Switzerland", "Syria",
    "Tajikistan", "Tunisia", "Turkey", "Turkmenistan",
    "Ukraine", "United Arab Emirates", "United Kingdom",
    "Uzbekistan", "Vatican City"),
  
  'Ethiopian' = c(
    "Angola", "Benin", "Botswana", "Burkina Faso", "Burundi", "Cameroon", "Central African Republic",
    "Chad", "Comoros", "Congo (Brazzaville)", "Congo (Kinshasa)", "Djibouti", "Equatorial Guinea",
    "Eritrea", "Eswatini", "Ethiopia", "Gabon", "Gambia", "Ghana", "Guinea", "Guinea-Bissau", "Ivory Coast",
    "Kenya", "Lesotho", "Liberia", "Madagascar", "Malawi", "Mali", "Mauritania", "Mauritius", "Mozambique",
    "Namibia", "Niger", "Nigeria", "Rwanda", "Sao Tome and Principe", "Senegal", "Seychelles", "Sierra Leone",
    "Somalia", "Sudan", "Tanzania", "Togo", "Uganda", "Zambia", "Zimbabwe"),
  
  'Australian' = c(
    "Australia", "Fiji", "Indonesia", "Kiribati", "Marshall Islands", "Micronesia", 
    "Nauru", "New Zealand", "Palau", "Papua New Guinea", "Samoa", "Solomon Islands", 
    "Tonga", "Tuvalu", "Vanuatu"),
  
  'Oceanian' = c(
    "American Samoa", "Cook Islands (a self-governing territory in free association with New Zealand)",
    "Fiji", "French Polynesia", "Guam", "Kiribati",
    "Marshall Islands", "Micronesia", "Nauru", "New Caledonia", "New Zealand", "Niue",
    "Norfolk Island", "Northern Mariana Islands",
    "Palau", "Papua New Guinea", "Pitcairn Islands", "Samoa", "Solomon Islands", "Tonga", "Tuvalu", "Vanuatu", "Wallis and Futuna"),
  
  'Antarctic' = c(
    "Antarctica", "Bouvet Island", "French Southern and Antarctic Lands", "Heard Island and McDonald Islands", "South Georgia and the South Sandwich Islands")
)

# ---- Clip Prediction Raster by Species ----
spp_clip_raster <- function(spp, worldbound, env_rs) {
  
  # Ensure the worldbound is in the same CRS as env_rs
  worldbound <- st_transform(worldbound, crs = crs(env_rs))
  
  # Convert the string of country codes in spp dataframe to a list
  spp <- spp %>%
    mutate(country_list = str_split(countryCode, ";"))
  
  # Initialize an empty list to store the clipped rasters
  clipped_rasters <- list()
  
  # Iterate over each unique Scientific.Name
  for (sci_name in unique(spp$Scientific.Name)) {
    # Extract the country codes associated with the current Scientific.Name
    country_codes <- spp %>%
      filter(Scientific.Name == sci_name) %>%
      pull(country_list) %>%
      unlist()
    
    # Subset the worldbound based on the country codes
    subset_polygons <- worldbound %>%
      filter(iso_3166_1_ %in% country_codes)
    
    # Union the subsetted polygons to get a single polygon in an sf object
    combined_polygon <- st_union(subset_polygons) %>%
      st_sf()
    
    # Clip the env_rs raster stack using the combined polygon
    clipped_raster <- mask(env_rs, combined_polygon)
    
    # Store the clipped raster in the list
    clipped_rasters[[sci_name]] <- clipped_raster
  }
  
  return(clipped_rasters)
}
# ---- Predict MaxEnt Models ----
maxent_predict <- function() {
  for (i in seq_along(maxent_results$models)) {
    # If the models have names, use those; otherwise, use the index
    model_name <- names(maxent_results$models)[i]
    if (is.null(model_name) || model_name == "") {
      model_name <- paste0("model_", i)
    }
    
    # Check if the raster clip for this model exists
    if (!model_name %in% names(clipped_rasters_list)) {
      cat("No raster clip found for model:", model_name, "\n")
      next
    }
    
    # Use the corresponding clipped raster for prediction
    prediction_raster <- clipped_rasters_list[[model_name]]
    
    # Predict the model
    prediction <- predict(prediction_raster, maxent_results$models[[i]], progress = 'text')
    
    # Define the filename for the TIFF using the model name
    model_name <- gsub(" ", "_", model_name)
    tif_file <- paste0("MaxEnt_Target_Predictions/", model_name, "_MaxEnt_Pred", ".tif")
    
    # Save the prediction as a TIFF
    writeRaster(prediction, filename = tif_file, format = "GTiff", overwrite = TRUE)
    
    cat("Saved prediction for", model_name, "as", tif_file, "\n")
  }
}