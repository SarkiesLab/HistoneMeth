#start with empty workspace

rm(list = ls(all = TRUE))

# turn off scientific notation for plots

options(scipen=10000)

#### set working directory ####

# here create new folder and set working directory within it

dir.create("~/NNMT_manuscript")
setwd("~/NNMT_manuscript")

# create subfolders for input, output and graphics

dir.create("input")

# into input folder, add input files 

dir.create("output")
dir.create("graphics")
dir.create("plot_data")

# make directories to deposit data underlying plots in main figures
sapply(1:4,
       function(i){
         
         dir.create(paste0("plot_data/Fig ", i))
         
         NULL
         
       })

sapply(2:13,
       function(i){
         
         dir.create(paste0("plot_data/Fig S", i))
         
         NULL
         
       })

#### DEFINE FUNCTIONS ####

apply.lm.to.counts <- function(countdata, lookuptable = LUT) {
  
  # create data frame with expression data and explanatory variables (age, sex, death) for linear regression
  # data will be transposed to have genes in columns and sample in rows
  # 3 new columns will be added: age, sex, cause of death for each sample
  temp_df <- as.data.frame(matrix(nrow = ncol(countdata), ncol = nrow(countdata) + 5))
  
  colnames(temp_df) <- c(row.names(countdata), "age_bracket", "sex", "death", "ischemic_time", "batch")
  row.names(temp_df) <- colnames(countdata)
  
  # insert tranposed data into new data frame
  temp_df[, 1:nrow(countdata)] <- t(countdata)
  
  # insert age bracket for each sample
  # use lookup table created above (lookuptable)
  temp_df[, "age_bracket"] <- sapply(str_split(lookuptable[row.names(temp_df), "age_bracket"], 
                                               pattern = "-"), function(x){
                                                 
                                                 mean(as.numeric(x))
                                                 
                                               })
  
  # from lookup table, add Sex and Death variables
  temp_df[, "sex"] <- lookuptable[row.names(temp_df), "Sex"]
  temp_df[, "death"] <- lookuptable[row.names(temp_df), "DeathScale"]
  
  # replace death = NA (unknown cause of death) with a 5th factor level
  temp_df[is.na(temp_df$death), "death"] <- 5
  
  temp_df[, "ischemic_time"] <- lookuptable[row.names(temp_df), "IschemicTime"]
  
  temp_df[, "batch"] <- lookuptable[row.names(temp_df), "Batch"]
  
  # check for any missing values and remove entry if so
  # here additional clause to work with GTEX v6, which has 4 brain tissues with no values for ischemic time
  # we aim to check if the tissue has NO values for ischemic time. in which case we will leave it
  # if it has some values, then will remove any entries with missing values
  
  if(sum(!is.na(temp_df$ischemic_time)) > 0){
    
    if(sum(is.na(temp_df$ischemic_time)) > 0){
      
      temp_df <- temp_df[!(is.na(temp_df$ischemic_time)), ]
      
    }
    
  }
  
  # do linear regression across all samples (all tissues combined)
  # perform linear regression with age (bracket midpoint), sex and cause of death
  
  # create dataframe to deposit residuals
  residuals_df <- as.data.frame(matrix(nrow = nrow(temp_df), ncol = ncol(temp_df) - 5))
  row.names(residuals_df) <- row.names(temp_df)
  colnames(residuals_df) <- colnames(temp_df)[1:(ncol(temp_df) - 5)]
  
  # loop over columns (genes); for each gene perform linear regression and deposit residuals in new data frame
  # sex is considered a factor automatically as a character string; as.factor() must be specified for cause of death
  
  if(sum(!is.na(temp_df$ischemic_time)) > 0){
    
    ischemicvariable <- "ischemic_time"
    
  } else {
    
    ischemicvariable <- NULL
    
  }
  
  if(length(unique(temp_df$sex)) == 1) {
    
    sexvariable <- NULL
    
  } else {
    
    sexvariable <- "as.factor(sex)"
    
  }
  
  if(length(unique(temp_df$batch)) == 1) {
    
    batchvariable <- NULL
    
  } else {
    
    batchvariable <- "as.factor(batch)"
    
  }
  
  modelvariables <- c("as.factor(death)", "age_bracket", ischemicvariable, batchvariable, sexvariable)
  
  for(j in 1:(ncol(residuals_df))){
    
    outcome <- paste0("temp_df[, ", j, "]")
    
    f <- as.formula(
      paste(outcome,
            paste(modelvariables, collapse = " + "),
            sep = " ~ "))
    
    residuals_df[, colnames(temp_df)[j]] <- lm(formula = f, data = temp_df)$residuals
    
  } # end of residuals loop
  
  return(t(residuals_df))
  
}

TCGA.apply.lm.to.counts <- function(countdata, lookuptable = TCGA_LUT) {
  
  # create data frame with expression data and explanatory variables (age, sex, death) for linear regression
  # data will be transposed to have genes in columns and sample in rows
  # 5 new columns will be added for variables of interest
  
  # restrict to samples present in lookuptable 
  tempcountdata <- countdata[, colnames(countdata) %in% lookuptable$SAMPID]
  
  temp_df <- as.data.frame(matrix(nrow = ncol(tempcountdata), ncol = nrow(tempcountdata) + 5))
  
  colnames(temp_df) <- c(row.names(tempcountdata), "race", "gender", "tumour_stage", "sequencing_centre", "days_to_birth")
  row.names(temp_df) <- colnames(tempcountdata)
  
  # insert tranposed data into new data frame
  temp_df[, 1:nrow(tempcountdata)] <- t(tempcountdata)
  
  # from lookup table, add variables
  temp_df[, "days_to_birth"] <- as.numeric(lookuptable[row.names(temp_df), "days_to_birth"])
  temp_df[, "gender"] <- lookuptable[row.names(temp_df), "gender"]
  temp_df[, "race"] <- lookuptable[row.names(temp_df), "race"]
  temp_df[, "tumour_stage"] <- lookuptable[row.names(temp_df), "tumour_stage"]
  temp_df[, "sequencing_centre"] <- lookuptable[row.names(temp_df), "sequencing_centre"]
  
  # check for any missing values and remove entry if so
  temp_df <- temp_df[!(is.na(temp_df$days_to_birth)), ]
  temp_df <- temp_df[!(is.na(temp_df$gender)), ]
  temp_df <- temp_df[!(is.na(temp_df$race)), ]
  temp_df <- temp_df[!(is.na(temp_df$tumour_stage)), ]
  temp_df <- temp_df[!(is.na(temp_df$sequencing_centre)), ]
  
  # do linear regression across all samples (all tissues combined)
  # perform linear regression with age (bracket midpoint), sex and cause of death
  
  # create dataframe to deposit residuals
  residuals_df <- as.data.frame(matrix(nrow = nrow(temp_df), ncol = ncol(temp_df) - 5))
  row.names(residuals_df) <- row.names(temp_df)
  colnames(residuals_df) <- colnames(temp_df)[1:(ncol(temp_df) - 5)]
  
  # loop over columns (genes); for each gene perform linear regression and deposit residuals in new data frame
  # sex is considered a factor automatically as a character string; sequencing centre is already considered as factor
  # here we note that in some projects, all sequencing centres are equal. Likewise gender for some cancer tpyes (e.g. ovarian)
  # need conditional to ignore seq centre, gender or tumour_stage if only one level exists among considered data in order to avoid an error
  
  if(length(unique(temp_df$race)) == 1){
    
    racevariable <- NULL
    
  } else {
    
    racevariable <- "race"
    
  }
  
  if(length(unique(temp_df$tumour_stage)) == 1){
    
    tumourstagevariable <- NULL
    
  } else {
    
    tumourstagevariable <- "tumour_stage"
    
  }
  
  if(length(unique(temp_df$gender)) == 1) {
    
    gendervariable <- NULL
    
  } else {
    
    gendervariable <- "gender"
  }
  
  if(length(unique(temp_df$sequencing_centre)) == 1) {
    
    seqcentrevariable <- NULL
    
  } else {
    
    seqcentrevariable <- "sequencing_centre"
    
  }
  
  modelvariables <- c("days_to_birth", gendervariable, seqcentrevariable, tumourstagevariable, racevariable)
  
  for(j in 1:(ncol(residuals_df))){
    
    outcome <- paste0("temp_df[, ", j, "]")
    
    f <- as.formula(
      paste(outcome,
            paste(modelvariables, collapse = " + "),
            sep = " ~ "))
    
    residuals_df[, colnames(temp_df)[j]] <- lm(formula = f, data = temp_df)$residuals
    
  } # end of residuals loop
  
  
  return(t(residuals_df))
  
}

GTEX.apply.lm.to.combined.counts <- function(countdata, lookuptable = LUT) {
  
  # create data frame with expression data and explanatory variables (age, sex, death) for linear regression
  # data will be transposed to have genes in columns and sample in rows
  # 3 new columns will be added: age, sex, cause of death for each sample
  temp_df <- as.data.frame(matrix(nrow = ncol(countdata), ncol = nrow(countdata) + 7))
  
  colnames(temp_df) <- c(row.names(countdata), "age_bracket", "sex", "death", "ischemic_time", "batch", "tissue", "donorID")
  row.names(temp_df) <- colnames(countdata)
  
  # insert tranposed data into new data frame
  temp_df[, 1:nrow(countdata)] <- t(countdata)
  
  # insert age bracket for each sample
  # use lookup table created above (lookuptable)
  temp_df[, "age_bracket"] <- sapply(str_split(lookuptable[row.names(temp_df), "age_bracket"], 
                                               pattern = "-"), function(x){
                                                 
                                                 mean(as.numeric(x))
                                                 
                                               })
  
  # from lookup table, add Sex and Death variables
  temp_df[, "sex"] <- lookuptable[row.names(temp_df), "Sex"]
  temp_df[, "death"] <- lookuptable[row.names(temp_df), "DeathScale"]
  
  # replace death = NA (unknown cause of death) with a 5th factor level
  temp_df[is.na(temp_df$death), "death"] <- 5
  
  # ischemic time rescaled so lmer doesn't complain about scale issues
  temp_df[, "ischemic_time"] <- lookuptable[row.names(temp_df), "IschemicTime"]/60 
  
  temp_df[, "batch"] <- lookuptable[row.names(temp_df), "Batch"]
  
  temp_df[, "tissue"] <- lookuptable[row.names(temp_df), "Tissue"]
  
  temp_df[, "donorID"] <- lookuptable[row.names(temp_df), "SUBJID"]
  
  # check for any missing values and remove entry if so
  # here additional clause to work with GTEX v6, which has 4 brain tissues with no values for ischemic time
  # we aim to check if the tissue has NO values for ischemic time. in which case we will leave it
  # if it has some values, then will remove any entries with missing values
  
  if(sum(!is.na(temp_df$ischemic_time)) > 0){
    
    if(sum(is.na(temp_df$ischemic_time)) > 0){
      
      temp_df <- temp_df[!(is.na(temp_df$ischemic_time)), ]
      
    }
    
  }
  
  # do linear regression across all samples (all tissues combined)
  # perform linear regression with age (bracket midpoint), sex and cause of death
  
  # create dataframe to deposit residuals
  residuals_df <- as.data.frame(matrix(nrow = nrow(temp_df), ncol = ncol(temp_df) - 7))
  row.names(residuals_df) <- row.names(temp_df)
  colnames(residuals_df) <- colnames(temp_df)[1:(ncol(temp_df) - 7)]
  
  # loop over columns (genes); for each gene perform linear regression and deposit residuals in new data frame
  # sex is considered a factor automatically as a character string; as.factor() must be specified for cause of death
  
  for(j in 1:(ncol(residuals_df))){
    
    residuals_df[, colnames(temp_df)[j]] <- resid(lmer(temp_df[, j] ~ (sex + as.factor(death) + age_bracket + ischemic_time)*tissue + (1|batch) + (1|donorID), data = temp_df))
    
  } # end of residuals loop
  
  return(t(residuals_df))
  
}

CCLE.correlate.metabolites.with.gene <- function(hgnc_gene = NULL,
                                                 ensembl_gene = NULL){

  if(is.null(ensembl_gene) & is.null(hgnc_gene)){
    
    stop("You must provide a gene ID")
    
  }
  
  if(!is.null(ensembl_gene) & !is.null(hgnc_gene)){
    
    stop("You must provide a single gene ID only")
    
  }
  
  if(!is.null(hgnc_gene)){
    
    message("Retrieving gene ID")
    
    ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    
    convertedensembl <- getBM(attributes = "ensembl_gene_id",
                              filters = "hgnc_symbol",
                              values = hgnc_gene,
                              mart = ensembl)
    
    if(length(convertedensembl) > 1){
      
      stop("HGNC symbol provided corresponds to multiple ENSEMBL gene IDs. Please provide an ENSEMBL ID instead")
      
    }
    
    ensembl_gene <- unlist(convertedensembl)
    
  }
  
  favgene_MOR_list <- lapply(CCLE_MOR_list, function(x){
    
    x[ensembl_gene, ]
    
  })
  
  favgene_z_list <- lapply(favgene_MOR_list, function(x){
    
    mean <- mean(log10(x + 1), na.rm = TRUE)
    sd <- sd(log10(x + 1), na.rm = TRUE)
    
    (log10(x + 1) - mean) / sd
    
  })       
  
  favgene_z_vec <- do.call(c, favgene_z_list)
  
  names(favgene_z_vec) <- str_extract(names(favgene_z_vec), pattern = "ACH.*")
  
  favgene_z_vec <- favgene_z_vec[match(row.names(metab_CCLE_z), names(favgene_z_vec))]
  
  matching <- names(favgene_z_vec)[!is.na(favgene_z_vec)]
  
  metabolomics_z <- metab_CCLE_z[match(matching, row.names(metab_CCLE_z)), ]
  favgene_z_vec <- favgene_z_vec[match(matching, names(favgene_z_vec))]
  
  results <- apply(metabolomics_z, 2, function(x){
    
    corr <- cor.test(unlist(x), favgene_z_vec, method = "pearson")
    
    c(corr$estimate, corr$p.value)
    
  })
  
  results <- t(results)
  colnames(results) <- c("corr", "p_value")
  
  padjust_vec <- p.adjust(results[, "p_value"], method = "BH")
  
  results <- cbind(results, padjust_vec)
  
  return(results)
  
}

CCLE.correlate.mRNA.with.protein <- function(hgnc_gene = NULL,
                                                 ensembl_gene = NULL){

  if(is.null(ensembl_gene) & is.null(hgnc_gene)){
    
    stop("You must provide a gene ID")
    
  }
  
  if(!is.null(ensembl_gene) & !is.null(hgnc_gene)){
    
    stop("You must provide a single gene ID only")
    
  }
  
  if(!is.null(hgnc_gene)){
    
    message("Retrieving gene ID")
    
    ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    
    convertedensembl <- getBM(attributes = "ensembl_gene_id",
                              filters = "hgnc_symbol",
                              values = hgnc_gene,
                              mart = ensembl)
    
    if(length(convertedensembl) > 1){
      
      stop("HGNC symbol provided corresponds to multiple ENSEMBL gene IDs. Please provide an ENSEMBL ID instead")
      
    }
    
    ensembl_gene <- unlist(convertedensembl)
    
  }
  
  favgene_MOR_list <- lapply(CCLE_MOR_list, function(x){
    
    x[ensembl_gene, ]
    
  })
  
  favgene_z_list <- lapply(favgene_MOR_list, function(x){
    
    mean <- mean(log10(x + 1), na.rm = TRUE)
    sd <- sd(log10(x + 1), na.rm = TRUE)
    
    (log10(x + 1) - mean) / sd
    
  })       
  
  favgene_z_vec <- do.call(c, favgene_z_list)
  
  names(favgene_z_vec) <- str_extract(names(favgene_z_vec), pattern = "ACH.*")
  
  matching <- row.names(proteo_CCLE_z)[match(names(favgene_z_vec), row.names(proteo_CCLE_z))]
  matching <- matching[!is.na(matching)]
  
  favgene_z_vec <- favgene_z_vec[match(matching, names(favgene_z_vec))]
  proteo_z <- proteo_CCLE_z[match(matching, row.names(proteo_CCLE_z)), ]
  
  favprot_z_vec <- proteo_CCLE_z[match(matching, row.names(proteo_CCLE_z)), ensembl_gene]
  
  results <- c(cor.test(favprot_z_vec, favgene_z_vec, method = "pearson")$estimate, cor.test(favprot_z_vec, favgene_z_vec, method = "pearson")$p.value)

  names(results) <- c("corr", "p_value")
  
  return(results)
  
}


CCLE.correlate.metabolites.with.protein <- function(hgnc_gene = NULL,
                                                    ensembl_gene = NULL){

  if(is.null(ensembl_gene) & is.null(hgnc_gene)){
    
    stop("You must provide a gene ID")
    
  }
  
  if(!is.null(ensembl_gene) & !is.null(hgnc_gene)){
    
    stop("You must provide a single gene ID only")
    
  }
  
  if(!is.null(hgnc_gene)){
    
    message("Retrieving gene ID")
    
    ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    
    convertedensembl <- getBM(attributes = "ensembl_gene_id",
                              filters = "hgnc_symbol",
                              values = hgnc_gene,
                              mart = ensembl)
    
    if(length(convertedensembl) > 1){
      
      stop("HGNC symbol provided corresponds to multiple ENSEMBL gene IDs. Please provide an ENSEMBL ID instead")
      
    }
    
    ensembl_gene <- unlist(convertedensembl)
    
  }
  
  
  matching <- row.names(proteo_CCLE_z)[match(row.names(metab_CCLE_z), row.names(proteo_CCLE_z))]
  matching <- matching[!is.na(matching)]
  
  metabolomics_z <- metab_CCLE_z[match(matching, row.names(metab_CCLE_z)), ]
  proteo_z <- proteo_CCLE_z[match(matching, row.names(proteo_CCLE_z)), ]
  
  favprot_z_vec <- proteo_CCLE_z[match(matching, row.names(proteo_CCLE_z)), ensembl_gene]
  
  results <- apply(metabolomics_z, 2, function(x){
    
    corr <- cor.test(unlist(x), favprot_z_vec, method = "pearson")
    
    c(corr$estimate, corr$p.value)
    
  })
  
  results <- t(results)
  colnames(results) <- c("corr", "p_value")
  
  padjust_vec <- p.adjust(results[, "p_value"], method = "BH")
  
  results <- cbind(results, padjust_vec)
  
  return(results)
  
}

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

convert.metabolite.labels.for.plot <- function(labels_vec){
  plotlabels <- str_replace_all(labels_vec, pattern = "\\.", replacement = "-")
  plotlabels <- str_remove(plotlabels, "^X")
  plotlabels[str_detect(plotlabels, pattern = "^C")] <- str_replace(plotlabels[str_detect(plotlabels, pattern = "^C")], pattern = "-", replacement = ":")
  plotlabels[str_detect(plotlabels, pattern = "^C")] <- str_replace(plotlabels[str_detect(plotlabels, pattern = "^C")], pattern = "-", replacement = " ")
  plotlabels[str_detect(plotlabels, pattern = "HILIC")] <- str_replace_all(plotlabels[str_detect(plotlabels, pattern = "HILIC")], pattern = "-", replacement = "")
  plotlabels[str_detect(plotlabels, pattern = "HILIC")] <- str_replace_all(plotlabels[str_detect(plotlabels, pattern = "HILIC")], pattern = "HILICneg", replacement = " (HILIC neg)")
  plotlabels[str_detect(plotlabels, pattern = "HILIC")] <- str_replace_all(plotlabels[str_detect(plotlabels, pattern = "HILIC")], pattern = "HILICpos", replacement = " (HILIC pos)")
  return(plotlabels)
}

convert.GTEX.tissue.labels.for.plot <- function(GTEX_labels,
                                                supershort = FALSE){
  
  # Here replace tissue names merely to shorten long names for plotting
  tissue_select_vec <- c("Adipose \\- Subcutaneous",
                         "Adipose \\- Visceral \\(Omentum\\)",
                         "Adrenal Gland",
                         "Artery \\- Aorta",
                         "Artery \\- Coronary",
                         "Artery  \\- Tibial",
                         "Brain \\- Amygdala",
                         "Brain \\- Anterior cingulate cortex \\(BA24\\)",
                         "Brain \\- Caudate \\(basal ganglia\\)",
                         "Brain \\- Cerebellar Hemisphere",
                         "Brain \\- Cerebellum",
                         "Brain \\- Cortex",
                         "Brain \\- Frontal Cortex \\(BA9\\)",
                         "Brain \\- Hippocampus",
                         "Brain \\- Hypothalamus",
                         "Brain \\- Nucleus accumbens \\(basal ganglia\\)",
                         "Brain \\- Putamen \\(basal ganglia\\)",
                         "Brain \\- Spinal cord \\(cervical c\\-1\\)",
                         "Brain \\- Substantia nigra",
                         "Breast \\- Mammary Tissue",
                         "Cells \\- Cultured fibroblasts",
                         "Cells \\- EBV-transformed lymphocytes",
                         "Colon \\- Sigmoid",
                         "Colon \\- Transverse",
                         "Esophagus \\- Gastroesophageal Junction",
                         "Esophagus \\- Mucosa",
                         "Esophagus \\- Muscularis",
                         "Heart \\- Atrial Appendage",
                         "Heart \\- Left Ventricle",
                         "Liver",
                         "Lung",
                         "Minor Salivary Gland",
                         "Muscle \\- Skeletal",
                         "Nerve \\- Tibial",
                         "Ovary",
                         "Pancreas",
                         "Pituitary",
                         "Prostate",
                         "Skin \\- Not Sun Exposed \\(Suprapubic\\)",
                         "Skin \\- Sun Exposed \\(Lower leg\\)",
                         "Small Intestine \\- Terminal Ileum",
                         "Spleen",
                         "Stomach",
                         "Testis",
                         "Thyroid",
                         "Uterus",
                         "Vagina",
                         "Whole Blood")
  
  if(supershort == FALSE){
  tissue_replace_vec <- c("Adipose (Subcut.)",
                          "Adipose (Visc.)",
                          "Adrenal Gland",
                          "Artery (Aorta)",
                          "Artery (Coronary)",
                          "Artery (Tibial)",
                          "Brain (Amygdala)",
                          "Brain (Ant.cing. cortex)",
                          "Brain (Caudate)",
                          "Brain (Cereb. Hemsph.)",
                          "Brain (Cerebellum)",
                          "Brain (Cortex)",
                          "Brain (Frontal Cortex)",
                          "Brain (Hippocampus)",
                          "Brain (Hypothalamus)",
                          "Brain (Nucl. acc.)",
                          "Brain (Putamen)",
                          "Brain (Spinal cord)",
                          "Brain (Subst. nigra)",
                          "Breast (Mammary)",
                          "Cultured fibroblasts",
                          "EBV-transf. lymphocytes",
                          "Colon (Sigmoid)",
                          "Colon (Transverse)",
                          "Esophagus (Gastr. Junc.)",
                          "Esophagus (Mucosa)",
                          "Esophagus (Muscularis)",
                          "Heart (Atrial Appendage)",
                          "Heart (Left Ventricle)",
                          "Liver",
                          "Lung",
                          "Min. Saliv. Gland",
                          "Muscle (Skeletal)",
                          "Nerve (Tibial)",
                          "Ovary",
                          "Pancreas",
                          "Pituitary",
                          "Prostate",
                          "Skin (Suprapubic)",
                          "Skin (Lower leg)",
                          "Small Intestine",
                          "Spleen",
                          "Stomach",
                          "Testis",
                          "Thyroid",
                          "Uterus",
                          "Vagina",
                          "Whole Blood")
  } else {
    tissue_replace_vec <- c("Adip.(Subc)",
                            "Adip. (Visc)",
                            "Ad. Gl.",
                            "Aorta",
                            "Cor.Art.",
                            "Tib.Art",
                            "Brain(Amyg.)",
                            "Brain(ACC)",
                            "Brain(Caud.)",
                            "Brain(Cer.Hem.)",
                            "Brain(Cereb.)",
                            "Brain(Cortex)",
                            "Brain(Fr.Crtx)",
                            "Brain(Hipp.)",
                            "Brain(Hyp.)",
                            "Brain(Nucl.acc.)",
                            "Brain(Put.)",
                            "Spinal cord",
                            "Brain(S.nigra)",
                            "Breast",
                            "Fibrob.",
                            "Lymph.",
                            "Colon(Sig.)",
                            "Colon(Trans.)",
                            "Esoph(GJ)",
                            "Esoph(Muc.)",
                            "Esoph(Musc.)",
                            "Heart(Atr.)",
                            "Heart(Vent.)",
                            "Liver",
                            "Lung",
                            "Saliv.Gland",
                            "Muscle",
                            "Nerve",
                            "Ovary",
                            "Pancreas",
                            "Pituitary",
                            "Prostate",
                            "Skin (hidden)",
                            "Skin (exposed)",
                            "Intestine",
                            "Spleen",
                            "Stomach",
                            "Testis",
                            "Thyroid",
                            "Uterus",
                            "Vagina",
                            "Blood")
  }
  
  for (i in 1:length(tissue_select_vec)){
    GTEX_labels <- str_replace_all(GTEX_labels, 
                                             pattern = tissue_select_vec[i], 
                                             replacement = tissue_replace_vec[i])
  }
  
  return(GTEX_labels)
  
}

metabolite.to.allgenes.correlation <- function(metabolite,
                                               extrageneset,
                                               genesetname){
  # extrageneset = SET_HMTs$ensembl_gene_id
  # genesetname = "allHMTs"
  
  # as we compute Z scores restrict to nuclear genes expressed set
  
  nuclear_genes_expressed <- nuclear_genes_expressed[nuclear_genes_expressed %in% row.names(CCLE_MOR_list[[1]])]
  
  nuclear_genes_expressed_and_geneset <- unique(c(nuclear_genes_expressed, extrageneset))
  
  geneset_MOR_list <- lapply(CCLE_MOR_list, function(x){
    
    x[nuclear_genes_expressed_and_geneset, ]
    
  })
  
  geneset_MOR_list <- lapply(geneset_MOR_list, function(x){
    
    x <- data.frame(x)
    x[genesetname, ] <- colSums(x[extrageneset, ])
    return(x)
    
  })
  
  geneset_z_list <- lapply(geneset_MOR_list, function(x){
    
    mean_vec <- apply(x, 1, function(thisrow){
      
      mean(log10(thisrow + 1), na.rm = TRUE)
      
    })
    
    sd_vec <- apply(x, 1, function(thisrow){
      
      sd(log10(thisrow + 1), na.rm = TRUE)
      
    })
    
    Zscores <- apply(x, 2, function(thiscolumn){
      
      (log10(thiscolumn + 1) - mean_vec) / sd_vec  
      
    })
    
  }) 
  
  geneset_z_df <- data.frame(do.call(cbind, geneset_z_list))
  colnames(geneset_z_df) <- str_replace(colnames(geneset_z_df), pattern = "\\.", replacement = "-")
  geneset_z_df <- geneset_z_df[, colnames(geneset_z_df) %in% row.names(metab_CCLE_z)]
  
  metab_CCLE_z_cens <- metab_CCLE_z[row.names(metab_CCLE_z) %in% colnames(geneset_z_df), ]
  metabolite_CCLE_z_cens <- metab_CCLE_z_cens[, metabolite]
  
  geneset_z_df <- geneset_z_df[, match(row.names(metab_CCLE_z_cens), colnames(geneset_z_df))]
  
  correlations <- apply(geneset_z_df, 1, function(x){
    
    gene_vec <- unlist(x)
    
    tryCatch({tempcor <- cor.test(gene_vec, metabolite_CCLE_z_cens)
    
    return(tempcor)}, error = function(e){return(list(estimate = NA, p.value = NA))})
    
  })
  
  generesults <- sapply(correlations, function(x){x$estimate})
  names(generesults) <- str_remove(names(generesults), pattern = "\\.cor")
  
  genepvals <- sapply(correlations, function(x){x$p.value})
  names(genepvals) <- str_remove(names(generesults), pattern = "\\.cor")
  
  resultstogether <- cbind(generesults, genepvals)
  
  return(resultstogether)
  
}

geneset.with.everything <- function(database = c("GTEX", "TCGA", "CCLE", "TARGET"),
                                    supplied_datalist = NULL,
                                    geneset_ensembl,
                                    basename = "mygeneset",
                                    save = TRUE,
                                    iterations = 100,
                                    chosensamples = NULL,
                                    filename = NULL){
  
  if(length(geneset_ensembl) < 2){
    
    stop("You must provide multiple genes. If you would like to use a single gene, use the anything.with.everything function!")
    
  }
  
  match.arg(database)
  
  # load in MRN/MOR-corrected pseudocounts generated elsewhere
  
  message("Reading in data")
  
  if(!is.null(supplied_datalist)){
    
    message("Using supplied datalist")
    
    if(database == "TCGA"){
      
    TCGA_LUT <- readRDS("output/TCGA_LUT.rds")
    datalist <- lapply(supplied_datalist, function(x){x[, colnames(x) %in% TCGA_LUT$SAMPID]})
      
    } else {
    
    datalist <- supplied_datalist
    
    }
    
  } else {
    
    if(database == "TCGA") {
      
      datalist <- readRDS("output/TCGA_MOR_list.rds")
      TCGA_LUT <- readRDS("output/TCGA_LUT.rds")
      datalist <- lapply(datalist, function(x){x[, colnames(x) %in% TCGA_LUT$SAMPID]})
      
    }
    
    if(database == "GTEX") {
      
      datalist <- readRDS("output/GTEX_MOR_list.rds")
      
    }
    
    if(database == "CCLE") {
      
      datalist <- readRDS("output/CCLE_MOR_list.rds")
      
    }
    
    if(database == "TARGET") {
      
      datalist <- readRDS("output/TARGET_MOR_list.rds")
      TARGETmetadata <- readRDS("output/TARGET_metadata_list.rds")
      
      datalist <- lapply(datalist, function(x){
        
        x <- x[, colnames(x) %in% TARGETmetadata[str_detect(TARGETmetadata$sample_type, "rimary"), "cases" ]]
        
      })
      
      datalist <- datalist[sapply(datalist, function(x){ncol(x) > 50})]
      
    }
    
  }
  
  if(any(!(geneset_ensembl %in% row.names(datalist[[1]])))){
    
    stop(paste0("The following gene is not present in chosen database:", geneset_ensembl[!(geneset_ensembl %in% row.names(datalist[[1]]))]))
    
  }
  
  datalist <- lapply(datalist, function(x){
    
    geneset_sum <- colSums(x[geneset_ensembl, ])
    df <- rbind(x, geneset_sum)
    row.names(df)[nrow(df)] <- basename
    return(df)
    
  })
  
  types <- names(datalist)
  
  message("Starting all ranks")
  
  allelse_ranks <- lapply(datalist, function(x){
    
    # add crumb - small random noise to break up ties. ranks will then be randomised among samples with 0 expression.
    crumb_ranks <- apply(x, 1, FUN = function(y){y + rnorm(y, 0, 0.01)})
    
    allelse_ranks <- apply(crumb_ranks, 2, rank)
    
    allelse_rankspercent <- allelse_ranks / (nrow(allelse_ranks) + 1)
    
    t(allelse_rankspercent)
    
  })
  
  allelse_allrankspercent <- data.frame(do.call(cbind, allelse_ranks))
  colnames(allelse_allrankspercent) <- str_replace_all(colnames(allelse_allrankspercent), pattern= "\\.", replacement = "-")
  
  if(!is.null(chosensamples)) {
    
    chosensamples_list <- chosensamples
    
  } else {
    
    message(paste0("Object '", database, "_randomsamples' not defined. Attempting to load from output..."))
    
    if(file.exists(paste0("output/chosensamples_", database, "_", iterations, "iterations.rds"))) {
      
      message(paste0("Now loading output/chosensamples_", database, "_", iterations, "iterations.rds"))
      
      chosensamples_list <- readRDS(paste0("output/chosensamples_", database, "_", iterations, "iterations.rds"))
      
    } else {
      
      message(paste0("No matching file found in output. Generating new random samples, which will be attached to output object"))
      
      chosensamples_list <- define.random.samples.for.iterations(database = database,
                                                                 supplied_datalist = supplied_datalist,
                                                                 iterations = iterations,
                                                                 save = FALSE)
      
    }
    
  }
  
  alliterations <- lapply(1:iterations, function(thisround) {
    
    chosensamples_vec <- chosensamples_list[[thisround]]
    
    tempranks <- allelse_allrankspercent[, chosensamples_vec]
    
    allcorrs <- vector(length = nrow(tempranks))
    names(allcorrs) <- row.names(tempranks)
    
    allpvals <- vector(length = nrow(tempranks))
    names(allpvals) <- row.names(tempranks)
    
    b4 <- Sys.time()
    message("Starting correlations")
    
    favranks <- unlist(tempranks[paste(basename), ])
    
    for(i in 1:nrow(tempranks)){
      
      ranks_vec <- unlist(tempranks[i, ])
      
      tempcor <- cor.test(favranks,
                          ranks_vec,
                          method = "spearman")
      
      allcorrs[i] <- tempcor$estimate
      allpvals[i] <- tempcor$p.value
      
      # give progress
      
      if (i %% 1000 == 0) {
        
        message(paste0("Iteration: ", thisround))
        
        print(paste("Combination", i, sep = " "))
        
        time_elapsed <- Sys.time() - b4
        print(paste("Time elapsed:", round(as.numeric(time_elapsed), digits = 2), attr(time_elapsed, "units"), sep = " "))
        
        timeperiteration <- time_elapsed / i
        
        remaining <- nrow(allelse_allrankspercent) - i
        print(paste("Estimated time remaining:", round((as.numeric(timeperiteration, units = "hours") * remaining), digits = 2), "hours", sep = " "))
        
      }
      
    }
    
    temp_list <- list(allcorrs = allcorrs, allpvals = allpvals)
    
    return(temp_list)
    
  }) # end of sapply for iterations
  
  allcorrs_df <- do.call(cbind, lapply(alliterations, function(x){unlist(x["allcorrs"])}))
  row.names(allcorrs_df) <- str_remove(row.names(allcorrs_df), "allcorrs\\.")
  
  allpvals_df <- do.call(cbind, lapply(alliterations, function(x){unlist(x["allpvals"])}))
  row.names(allpvals_df) <- str_remove(row.names(allpvals_df), "allpvals\\.")
  
  output_corrs_vec <- apply(allcorrs_df, 1, median)
  
  output_p_vec_raw <- apply(allpvals_df, 1, gm_mean)
  output_p_vec_FDR <- p.adjust(output_p_vec_raw, method = "BH")
  
  output_list <- list(medians = output_corrs_vec,
                      FDR = output_p_vec_FDR,
                      fullcorrs = allcorrs_df,
                      fullpvalues = allpvals_df,
                      randomsamples = chosensamples_list)
  
  if(save == TRUE){
    
    if(!is.null(filename)){
      
      saveRDS(output_list, paste0("output/", filename,".rds"))
      
    } else {
    
    saveRDS(output_list, paste0("output/", database, "_", basename, "_correlate_with_everything.rds"))
      
    }
  }
  
  return(output_list)
  
}

add.geneset.to.everything <- function(database = c("GTEX", "TCGA", "CCLE"),
                                      supplied_datalist = NULL,
                                      geneset_ensembl,
                                      genesetname = "mygeneset",
                                      everything_object){
  
  match.arg(database)
  
  everything_gene <- names(everything_object$medians[everything_object$medians == 1])
  
  random_samples <- everything_object$randomsamples
  
  if(!is.null(supplied_datalist)){
    
    message("Using supplied datalist")
    
    if(database == "TCGA"){
      
      TCGA_LUT <- readRDS("output/TCGA_LUT.rds")
      datalist <- lapply(supplied_datalist, function(x){x[, colnames(x) %in% TCGA_LUT$SAMPID]})
      
    } else {
    
    datalist <- supplied_datalist
    
    }
    
  } else {
  
  if(database == "TCGA") {
    
    datalist <- readRDS("output/TCGA_MOR_list.rds")
    TCGA_LUT <- readRDS("output/TCGA_LUT.rds")
    datalist <- lapply(datalist, function(x){x[, colnames(x) %in% TCGA_LUT$SAMPID]})
    
  }
  
  if(database == "GTEX") {
    
    datalist <- readRDS("output/GTEX_MOR_list.rds")
    
  }
  
  if(database == "CCLE") {
    
    datalist <- readRDS("output/CCLE_MOR_list.rds")
    
  }
  
  }
    
  datalist <- lapply(datalist, as.data.frame)
  
  for(i in 1:length(datalist)){

    datalist[[i]][paste0(genesetname), ] <- colSums(datalist[[i]][geneset_ensembl, ])
    
  }
  
  everything_and_geneset_ranks <- lapply(datalist, function(x){
    
    # add crumb - small random noise to break up ties. ranks will then be randomised among samples with 0 expression.
    everything_gene_crumb <- x[everything_gene, ] + rnorm(x[everything_gene, ], 0, 0.01)
    
    everything_gene_ranks <- rank(everything_gene_crumb)
    
    everything_gene_rankspercent <- everything_gene_ranks / (length(everything_gene_ranks) + 1)
    
    geneset_crumb <- x[genesetname, ] + rnorm(x[genesetname, ], 0, 0.01)
    
    geneset_ranks <- rank(geneset_crumb)
    
    geneset_rankspercent <- geneset_ranks / (length(geneset_ranks) + 1)
    
    return(data.frame(rbind(everything_gene_rankspercent, geneset_rankspercent)))
    
  })
  
  everything_and_geneset_ranks_df <- data.frame(do.call(cbind, everything_and_geneset_ranks))
  colnames(everything_and_geneset_ranks_df) <- str_replace_all(do.call(c, lapply(everything_and_geneset_ranks, colnames)), pattern= "\\.", replacement = "-")
  
  tempcorrs <- sapply(1:length(random_samples), function(thisround){

    thisround_samples <- random_samples[[thisround]]

    everything_gene_vec <- unlist(everything_and_geneset_ranks_df["everything_gene_rankspercent", thisround_samples])
    geneset_vec <- unlist(everything_and_geneset_ranks_df["geneset_rankspercent", thisround_samples])
    
    cor.test(everything_gene_vec, geneset_vec, method = "spearman")$estimate
    
  })
  
  temppvals <- sapply(1:length(random_samples), function(thisround){
    
    thisround_samples <- random_samples[[thisround]]
    
    everything_gene_vec <- unlist(everything_and_geneset_ranks_df["everything_gene_rankspercent", thisround_samples])
    geneset_vec <- unlist(everything_and_geneset_ranks_df["geneset_rankspercent", thisround_samples])
    
    cor.test(everything_gene_vec, geneset_vec, method = "spearman")$p.value
    
  })
  
  fullcorrs_added <- rbind(everything_object$fullcorrs, tempcorrs)
  row.names(fullcorrs_added)[nrow(fullcorrs_added)] <- genesetname
  
  medians_added <- apply(fullcorrs_added, 1, median)
  
  allpvals_added <- rbind(everything_object$fullpvalues, temppvals)
  row.names(allpvals_added)[nrow(allpvals_added)] <- genesetname
  
  output_p_vec_raw_added <- apply(allpvals_added, 1, gm_mean)
  output_p_vec_FDR_added <- p.adjust(output_p_vec_raw_added, method = "BH")
  
  output_list <- list(medians = medians_added,
                      FDR = output_p_vec_FDR_added,
                      fullcorrs = fullcorrs_added,
                      fullpvalues = allpvals_added,
                      randomsamples = everything_object$randomsamples)
  
  return(output_list)
  
}

add.geneset.to.everything.subtype <- function(database = c("GTEX", "TCGA", "CCLE"),
                                      geneset_ensembl,
                                      genesetname = "mygeneset",
                                      everything_object_list){
  
  match.arg(database)
  
  everything_gene <- row.names(everything_object_list[[1]][everything_object_list[[1]][, "allcorrs"] > 0.9999, ])[1]
  
  if(database == "TCGA") {
    
    datalist <- readRDS("output/TCGA_MOR_list.rds")
    TCGA_LUT <- readRDS("output/TCGA_LUT.rds")
    datalist <- lapply(datalist, function(x){x[, colnames(x) %in% TCGA_LUT$SAMPID]})
    
  }
  
  if(database == "GTEX") {
    
    datalist <- readRDS("output/GTEX_MOR_list.rds")
    
  }
  
  if(database == "CCLE") {
    
    datalist <- readRDS("output/CCLE_MOR_list.rds")
    
  }
  
  datalist <- lapply(datalist, as.data.frame)
  
  for(i in 1:length(datalist)){

    datalist[[i]][paste0(genesetname), ] <- colSums(datalist[[i]][row.names(datalist[[i]]) %in% geneset_ensembl, ])
    
  }
  
  everything_object_list_added <- lapply(1:length(datalist), function(i){

  thistype <- names(datalist)[[i]]
     
  everything_gene_vec <- unlist(datalist[[thistype]][everything_gene, ]) + rnorm(ncol(datalist[[thistype]]), 0, 0.01)
  
  geneset_vec <- unlist(datalist[[thistype]][genesetname, ]) + rnorm(ncol(datalist[[thistype]]), 0, 0.01)
  
  thiscorr <- c(allcorrs = cor.test(everything_gene_vec, geneset_vec, method = "spearman")$estimate, 
                allpvals = cor.test(everything_gene_vec, geneset_vec, method = "spearman")$p.value)
  
  output <- rbind(everything_object_list[[thistype]], thiscorr)
  row.names(output)[nrow(output)] <- genesetname
  
  return(output)
  
  })
  
  names(everything_object_list_added) <- names(datalist)
  
  return(everything_object_list_added)
  
}

anything.with.everything <- function(database = c("GTEX", "TCGA", "CCLE"),
                                     supplied_datalist = NULL,
                                     ensembl_gene = NULL,
                                     hgnc_gene = NULL,
                                     basename = "mygene",
                                     iterations = 100,
                                     save = TRUE,
                                     extrageneset = NULL,
                                     genesetname = "mygeneset",
                                     filename = NULL,
                                     chosensamples = NULL){
  
  if(is.null(ensembl_gene) & is.null(hgnc_gene)){
    
    stop("You must provide a gene ID")
    
  }
  
  if(!is.null(ensembl_gene) & !is.null(hgnc_gene)){
    
    stop("You must provide a single gene ID only")
    
  }
  
  match.arg(database)
  
  if(!is.null(hgnc_gene)){
    
    require(biomaRt)
    
    message("Retrieving gene ID")
    
    ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    
    convertedensembl <- getBM(attributes = "ensembl_gene_id",
                              filters = "hgnc_symbol",
                              values = hgnc_gene,
                              mart = ensembl)
    
    if(length(convertedensembl) > 1){
      
      stop("HGNC symbol provided corresponds to multiple ENSEMBL gene IDs. Please provide an ENSEMBL ID instead")
      
    }
    
    ensembl_gene <- convertedensembl
    
  }
  
  # load in MRN/MOR-corrected pseudocounts generated elsewhere
  
  message("Reading in data")
  
  if(!is.null(supplied_datalist)){
    
    message("Using supplied datalist...")
    
    if(database == "TCGA") {
      
      TCGA_LUT <- readRDS("output/TCGA_LUT.rds")
      datalist <- lapply(supplied_datalist, function(x){x[, colnames(x) %in% TCGA_LUT$SAMPID]})
      
    } else {
    
    datalist <- supplied_datalist
    
    }
    
  } else {
    
    if(database == "TCGA") {
      
      datalist <- readRDS("output/TCGA_MOR_list.rds")
      TCGA_LUT <- readRDS("output/TCGA_LUT.rds")
      datalist <- lapply(datalist, function(x){x[, colnames(x) %in% TCGA_LUT$SAMPID]})
      
    }
    
    if(database == "GTEX") {
      
      datalist <- readRDS("output/GTEX_MOR_list.rds")
      
    }
    
    if(database == "CCLE") {
      
      datalist <- readRDS("output/CCLE_MOR_list.rds")
      
    }
    
  }
  
  if(!(ensembl_gene %in% row.names(datalist[[1]]))){
    
    stop("Gene is not present in chosen database")
    
  }
  
  if(any(!(extrageneset %in% row.names(datalist[[1]])))){
    
    stop(paste0("The following gene is not present in chosen database:", geneset_ensembl[!(geneset_ensembl %in% row.names(datalist[[1]]))]))
    
  }
  
  types <- names(datalist)
  
  datalist <- lapply(datalist, as.data.frame)
  
  if(!is.null(extrageneset)){
    
    for(i in 1:length(datalist)){
      
      datalist[[i]][paste0(genesetname), ] <- colSums(datalist[[i]][extrageneset, ])
      
    }
    
  }
  
  message("Starting all ranks")
  
  allelse_ranks <- lapply(datalist, function(x){
    
    # add crumb - small random noise to break up ties. ranks will then be randomised among samples with 0 expression.
    crumb_ranks <- apply(x, 1, FUN = function(y){y + rnorm(y, 0, 0.01)})
    
    allelse_ranks <- apply(crumb_ranks, 2, rank)
    
    allelse_rankspercent <- allelse_ranks / (nrow(allelse_ranks) + 1)
    
    t(allelse_rankspercent)
    
  })
  
  allelse_allrankspercent <- data.frame(do.call(cbind, allelse_ranks))
  colnames(allelse_allrankspercent) <- str_replace_all(colnames(allelse_allrankspercent), pattern= "\\.", replacement = "-")
  
  if(!is.null(chosensamples)) {
    
    chosensamples_list <- chosensamples
    
  } else {
    
    message(paste0("Object '", database, "_randomsamples' not defined. Attempting to load from output..."))
    
    if(file.exists(paste0("output/chosensamples_", database, "_", iterations, "iterations.rds"))) {
      
      message(paste0("Now loading output/chosensamples_", database, "_", iterations, "iterations.rds"))
      
      chosensamples_list <- readRDS(paste0("output/chosensamples_", database, "_", iterations, "iterations.rds"))
      
    } else {
      
      message(paste0("No matching file found in output. Generating new random samples, which will be attached to output object"))
      
      chosensamples_list <- define.random.samples.for.iterations(database = database,
                                                                 supplied_datalist = supplied_datalist,
                                                                 iterations = iterations,
                                                                 save = FALSE)
      
    }
    
  }
  
  alliterations <- lapply(1:iterations, function(thisround) {
    
    chosensamples_vec <- chosensamples_list[[thisround]]
    
    tempranks <- allelse_allrankspercent[, chosensamples_vec]
    
    allcorrs <- vector(length = nrow(tempranks))
    
    names(allcorrs) <- row.names(tempranks) 
    
    allpvals <- vector(length = nrow(tempranks))
    
    names(allpvals) <- row.names(tempranks) 
    
    b4 <- Sys.time()
    message("Starting correlations")
    
    favranks <- unlist(tempranks[paste(ensembl_gene), ])
    
    for(i in 1:nrow(tempranks)){
      
      ranks_vec <- unlist(tempranks[i, ])
      
      tempcor <- cor.test(favranks,
                          ranks_vec,
                          method = "spearman")
      
      allcorrs[i] <- tempcor$estimate
      allpvals[i] <- tempcor$p.value
      
      # give progress
      
      if (i %% 1000 == 0) {
        
        message(paste0("Iteration: ", thisround))
        
        print(paste("Combination", i, sep = " "))
        
        time_elapsed <- Sys.time() - b4
        print(paste("Time elapsed:", round(as.numeric(time_elapsed), digits = 2), attr(time_elapsed, "units"), sep = " "))
        
        timeperiteration <- time_elapsed / i
        
        remaining <- nrow(allelse_allrankspercent) - i
        print(paste("Estimated time remaining:", round((as.numeric(timeperiteration, units = "hours") * remaining), digits = 2), "hours", sep = " "))
        
      }
      
    } 
    
    temp_list <- list(allcorrs = allcorrs, allpvals = allpvals)
    
    return(temp_list)
    
  }) # end of sapply for iterations
  
  allcorrs_df <- do.call(cbind, lapply(alliterations, function(x){x$allcorrs}))
  allpvals_df <- do.call(cbind, lapply(alliterations, function(x){x$allpvals}))
  
  output_corrs_vec <- apply(allcorrs_df, 1, median)
  
  output_p_vec_raw <- apply(allpvals_df, 1, gm_mean)
  output_p_vec_FDR <- p.adjust(output_p_vec_raw, method = "BH")
  
  output_list <- list(medians = output_corrs_vec,
                      FDR = output_p_vec_FDR,
                      fullcorrs = allcorrs_df,
                      fullpvalues = allpvals_df,
                      randomsamples = chosensamples_list)
  
  if(save == TRUE){
    
    if(!is.null(filename)){
      
      saveRDS(output_list, paste0("output/", filename,".rds"))
      
    } else {
      
      saveRDS(output_list, paste0("output/", database, "_", basename, "_correlate_with_everything.rds"))
      
    }
  }
  
  return(output_list)
  
}

analyse.by.subtype <- function(database = c("GTEX", "TCGA", "CCLE", "TARGET"),
                                       geneset_ensembl = NULL,
                                       extrageneset_ensembl = NULL,
                                       singlegene = NULL,
                                       genebasename,
                               extragenesetbasename = NULL) {
  
  match.arg(database)
  
  if(database == "TCGA") {
    
    datalist <- readRDS("output/TCGA_MOR_list.rds")
    TCGA_LUT <- readRDS("output/TCGA_LUT.rds")
    
    datalist <- lapply(datalist, function(x){x[, colnames(x) %in% TCGA_LUT$SAMPID]})  
    
  }
  
  if(database == "GTEX") {
    
    datalist <- readRDS("output/GTEX_MOR_list.rds")
    
  }
  
  if(database == "CCLE") {
    
    datalist <- readRDS("output/CCLE_MOR_list.rds")
    
  }
  
  
  if(!is.null(geneset_ensembl)){
    
    if(sum(row.names(datalist[[1]]) %in% geneset_ensembl) < length(geneset_ensembl)){
      
      message("The following genes are missing from the data:")
      message(paste0(geneset_ensembl[!(geneset_ensembl %in% row.names(datalist[[1]]))], collapse = ","))
      
    }
    
    datalist <- lapply(datalist, function(x){
      
      geneset_sum <- colSums(x[row.names(x) %in% geneset_ensembl, ])
      df <- rbind(x, geneset_sum)
      row.names(df)[nrow(df)] <- genebasename
      return(df)
      
    })
    
  }
  
  if(!is.null(extrageneset_ensembl)){
    
    datalist <- lapply(datalist, function(x){
      
      geneset_sum <- colSums(x[extrageneset_ensembl, ])
      df <- rbind(x, geneset_sum)
      row.names(df)[nrow(df)] <- extragenesetbasename
      return(df)
      
    })
    
  }
  
  if(!is.null(singlegene)){
    
    datalist <- lapply(datalist, function(x){
      
      geneset_sum <- x[singlegene, ]
      df <- rbind(x, geneset_sum)
      row.names(df)[nrow(df)] <- genebasename
      return(df)
      
    })
    
  }
  
  correlations <- lapply(1:length(datalist), function(i){

    #  message(paste0("Starting ", names(datalist)[i]))
    
    allcorrs <- vector(length = nrow(datalist[[i]]))
    names(allcorrs) <- row.names(datalist[[i]])
    
    allpvals <- vector(length = nrow(datalist[[i]]))
    names(allpvals) <- row.names(datalist[[i]])
    
    fav_vec <- unlist(datalist[[i]][paste(genebasename), ])
    
    for(j in 1:nrow(datalist[[i]])){

      temp_vec <- unlist(datalist[[i]][j, ])
      
      crumb_ranks <- temp_vec + rnorm(temp_vec, 0, 0.01)
      
      tempcor <- cor.test(fav_vec,
                          crumb_ranks,
                          method = "spearman")
      
      allcorrs[j] <- tempcor$estimate
      allpvals[j] <- tempcor$p.value
      
    }
    
    output_df <- data.frame(cbind(allcorrs, allpvals))
    
    return(output_df)
    
  })
  
  names(correlations) <- names(datalist)
  
  return(correlations)
  
}

calculate.relative.reciprocal.scores <- function(everything_object1,
                                                 everything_object2,
                                                 gene1_name = NULL,
                                                 gene2_name = NULL){

  focusgene1 <- names(everything_object1[everything_object1 > 0.9999])
  focusgene2 <- names(everything_object2[everything_object2 > 0.9999])[1]
  
  if(is.null(gene1_name)){
    
    gene1_name = focusgene1
    
  }
  
  
  if(is.null(gene2_name)){
    
    gene2_name = focusgene2
    
  }
  
  if(!(focusgene1 %in% names(everything_object2))){
    
    stop(paste0(gene1_name, " is not present in everything_object2"))
    
  }
  
  if(!(focusgene2 %in% names(everything_object1))){
    
    stop(paste0(gene2_name, "is not present in everything_object1"))
    
  }
  
  focusgene1_rank_in_object2 <- which(names(everything_object2[order(everything_object2)]) == focusgene1)
  focusgene2_rank_in_object1 <- which(names(everything_object1[order(everything_object1)]) == focusgene2)
  
  reciprocal_score <- (focusgene1_rank_in_object2)^2 + (focusgene2_rank_in_object1)^2
  
  outputvec <- c(focusgene1_rank_in_object2,
                     focusgene2_rank_in_object1,
                     reciprocal_score)
  
  names(outputvec) <- c(paste0(gene1_name, "_rank_in_", gene2_name, "_everything"),
                         paste0(gene2_name, "_rank_in_", gene1_name, "_everything"),
                         "reciprocal_score")
  
  return(outputvec)
  
}

define.random.samples.for.iterations <- function(database = c("GTEX", "TCGA", "CCLE", "TARGET"),
                                                 supplied_datalist = NULL, 
                                                 iterations = 100,
                                                 save = TRUE,
                                                 filename = NULL){
  
  if(!is.null(supplied_datalist)){
    
    if(database == "TCGA"){
      
      datalist <- lapply(supplied_datalist, function(x){x[, colnames(x) %in% TCGA_LUT$SAMPID]})  
      
    } else {
    
    datalist <- supplied_datalist
    
    }
    
  } else {
    
    if(database == "TCGA") {
      
      datalist <- readRDS("output/TCGA_MOR_list.rds")
      TCGA_LUT <- readRDS("output/TCGA_LUT.rds")
      
      datalist <- lapply(datalist, function(x){x[, colnames(x) %in% TCGA_LUT$SAMPID]})  
      
    }
    
    if(database == "GTEX") {
      
      datalist <- readRDS("output/GTEX_MOR_list.rds")
      
    }
    
    if(database == "CCLE") {
      
      datalist <- readRDS("output/CCLE_MOR_list.rds")
      
    }
    
    if(database == "TARGET") {
      
      datalist <- readRDS("output/TARGET_MOR_list.rds")
      TARGETmetadata <- readRDS("output/TARGET_metadata_list.rds")
      
      datalist <- lapply(datalist, function(x){
        
        x <- x[, colnames(x) %in% TARGETmetadata[str_detect(TARGETmetadata$sample_type, "rimary"), "cases" ]]
        
      })
      
      datalist <- datalist[sapply(datalist, function(x){ncol(x) > 50})]
      
    } 
    
  }
  
  chosensamples_list <- lapply(1:iterations, function(thisround){
    
    chosensamples <- lapply(datalist, function(x){
      
      if(database == "GTEX"){
        n = 100
      }
      
      if(database == "TCGA"){
        n = 36
      }
      
      if(database == "CCLE"){
        n = 20
      }
      
      if(database == "TARGET"){
        n = 60
      }
      
      sample(colnames(x), n)
      
    })  
    
    chosensamples_vec <- do.call(c, chosensamples)
    chosensamples_vec <- str_replace_all(chosensamples_vec, pattern= "\\.", replacement = "-")
    
    return(chosensamples_vec)
    
  })
  
  if(save == TRUE){
    
    if(!is.null(filename)){
      
      saveRDS(chosensamples_list, paste0("output/", filename, ".rds"))
      
    } else {
      
      saveRDS(chosensamples_list, paste0("output/chosensamples_", database, "_", iterations, "iterations.rds"))
      
    }
    
  }
  
  return(chosensamples_list)
  
}

whole.geneset.metabolite.correlation <- function(geneset){
  
  geneset_MOR_list <- lapply(CCLE_MOR_list, function(x){
    
    colSums(x[geneset, ])
    
  })
  
  geneset_z_list <- lapply(geneset_MOR_list, function(x){
    
    mean <- mean(log10(x + 1), na.rm = TRUE)
    sd <- sd(log10(x + 1), na.rm = TRUE)
    
    (log10(x + 1) - mean) / sd
    
  }) 
  
  geneset_z_vec <- do.call(c, geneset_z_list)
  
  names(geneset_z_vec) <- str_extract(names(geneset_z_vec), pattern = "ACH.*")
  
  geneset_z_vec <- geneset_z_vec[match(row.names(metab_CCLE_z), names(geneset_z_vec))]
  
  matching <- names(geneset_z_vec)[!is.na(geneset_z_vec)]
  
  metabolomics_z <- metab_CCLE_z[match(matching, row.names(metab_CCLE_z)), ]
  geneset_z_vec <- geneset_z_vec[match(matching, names(geneset_z_vec))]
  
  genesetresults <- apply(metabolomics_z, 2, function(x){
    
    cor.test(geneset_z_vec, x, method = "pearson")$estimate
    
  })
  
  genesetpvalues <- apply(metabolomics_z, 2, function(x){
    
    cor.test(geneset_z_vec, x, method = "pearson")$p.value
    
  })
  
  resultstogether <- cbind(genesetresults, genesetpvalues)
  
  return(resultstogether)
  
}

map2color <- function(x, 
                      pal,
                      limits = NULL){
  
  if(is.null(limits)){
    limits = range(x)
  }
  
  pal[findInterval(x, seq(limits[1], limits[2], length.out = length(pal) + 1), all.inside = TRUE)]
  
}

plot.selection.against.everything <- function(everything_object,
                                              selection_ensembl,
                                              title = "Everything vs selection"){
  
  selvec <- names(everything_object) %in% selection_ensembl
  rhovec <- data.frame(everything_object[selvec])
  colnames(rhovec) <- "rho"
  
  everything_df <- data.frame(everything_object)
  colnames(everything_df) <- "rho"

unique(row.names(everything_df))

  densityforplot <- with(density(everything_df$rho, na.rm = TRUE), data.frame(x, y))
  low_extreme <- densityforplot[densityforplot$x < quantile(everything_df$rho, 0.025, na.rm = TRUE), ]
  hi_extreme <- densityforplot[densityforplot$x > quantile(everything_df$rho, 0.975, na.rm = TRUE), ]
  
  print(ggplot(data = everything_df, aes(x = rho)) + 
          geom_density() + 
          geom_area(data = low_extreme, aes(x = x, y = y), alpha = 0.8, fill = "dodgerblue") +
          geom_area(data = hi_extreme, aes(x = x, y = y), alpha = 0.8, fill = "dodgerblue") +
          # ggtitle(title) +
          geom_density(data = data.frame(rhovec), (aes(x = rho)), color = "red") + 
          coord_cartesian(xlim = c(-0.6, 0.6)) +
          geom_point(data = rhovec,
                     aes(x = rho, y = -0.25),
                     color = "red",
                     size = 2,
                     alpha = 0.7) + 
          xlab(substitute("Spearman's"~rho)) +
          ylab("Density")+
          theme_classic()) + 
    theme(axis.text.x = element_text(colour = "black", 
                                     size = 10),
          axis.text.y = element_text(colour = "black"),
          axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10))

  
}

plot.selection.against.everything.facet <- function(everything_object,
                                                    selection_ensembl,
                                                    title = "Everything vs selection",
                                                    pointsize = 0.1){
  
  everything_object_ids_df <- lapply(everything_object, function(x){
    
    x[, "row_names"] <- row.names(x)
    x <- x[, !(colnames(x) == "p_val")]
    return(x)
    
  })
  
  everything_melted <- melt(everything_object_ids_df)
  colnames(everything_melted) <- c("ensembl_gene_id", "variable", "rho", "subtype")
  
  sel_df <- everything_melted[everything_melted$ensembl_gene_id %in% selection_ensembl, c("rho", "subtype")]
  
  low_df <- lapply(everything_object, function(subobject){
    
    my_df <- subobject[!(subobject$rho == 1), ]
    my_df <- my_df[!is.na(my_df$rho), ]
    densityforplot <- with(density(my_df$rho, na.rm = TRUE), data.frame(x, y))
    low_extreme <- densityforplot[densityforplot$x < quantile(my_df$rho, 0.025, na.rm = TRUE), ]
    
  })
  
  low_docalled <- do.call(rbind, low_df)
  low_docalled[, "subtype"] <- str_remove(row.names(low_docalled), pattern = "\\..*$")
  
  high_df <- lapply(everything_object, function(subobject){
    
    my_df <- subobject[!(subobject$rho == 1), ]
    my_df <- my_df[!is.na(my_df$rho), ]
    densityforplot <- with(density(my_df$rho, na.rm = TRUE), data.frame(x, y))
    low_extreme <- densityforplot[densityforplot$x > quantile(my_df$rho, 0.975, na.rm = TRUE), ]
    
  })
  
  high_docalled <- do.call(rbind, high_df)
  high_docalled[, "subtype"] <- str_remove(row.names(high_docalled), pattern = "\\..*$")
  
  print(ggplot(data = everything_melted, aes(x = rho)) + 
          geom_density(stat = "density", col = "red") +
          geom_density(data = sel_df, aes(x = rho), col = "blue") +
          geom_area(data = low_docalled, aes(x = x, y = y), alpha = 0.8, fill = "grey") +
          geom_area(data = high_docalled, aes(x = x, y = y), alpha = 0.8, fill = "grey") +
          geom_jitter(data = sel_df,
                      aes(x = rho, y = 0),
                      height = 0.5,
                      width = 0,
                      color = "blue",
                      size = pointsize) +
          facet_wrap(~ subtype) +
          ggtitle(title) +
          theme_classic2()+
          theme(strip.text.x = element_text(size = 8, margin = margin(0, 0.2, 0, 0.2, "cm"))) +
          ylab("Density")+
          xlab(substitute(rho)))
  
}

plot.selection.against.everything.volcano <- function(everything_object,
                                                      selection_ensembl,
                                                      selection_colours,
                                                      filename,
                                                      ylim = NULL,
                                                      xlim = NULL,
                                                      plotheight = 2.5,
                                                      plotwidth = 2.5,
                                                      xaxislabel = substitute("Pearson's"~italic(r)),
                                                      yaxislabel = substitute(-log[10]~"(false discovery rate)",),
                                                      yaxisbreaks = NULL,
                                                      xaxisbreaks = NULL,
                                                      fileoutput = "pdf"
                                                      ){

  if(!is.list(selection_ensembl)){
    selection_ensembl <- list(selection_ensembl)
  }
  
  everything_df <- data.frame(everything_object)
  colnames(everything_df) <- c("rho", "FDR")
  
  everything_df[, "gene"] <- row.names(everything_df)

  everything_df$gene <- factor(everything_df$gene, levels = everything_df[order(everything_df$rho, decreasing = FALSE), "gene"])
  everything_df[, "barcolour"] <- "grey"
  everything_df[everything_df$rho < quantile(everything_df$rho, 0.025, na.rm = TRUE), "barcolour"] <- "gray55"
  everything_df[everything_df$rho > quantile(everything_df$rho, 0.975, na.rm = TRUE), "barcolour"] <- "gray55"
  
  selection_df <- everything_df[unlist(selection_ensembl), ]
  selection_df[, "hgnc_symbol"] <- SET_HMTs[match(selection_df$gene, SET_HMTs$ensembl_gene_id), "hgnc_symbol"]
  
  for(i in 1:length(selection_ensembl)){
    
    selection_df[selection_ensembl[[i]], "barcolour"] <- selection_colours[i]
    
  }
  
  if(fileoutput == "pdf"){
  pdf(paste0("graphics/", filename, ".pdf"),
      height = plotheight,
      width = plotwidth)} else {
        png(paste0("graphics/", filename, ".png"),
            height = plotheight,
            width = plotwidth,
            units = "in",
            res = 1200)
      }
  
  print(ggplot(aes(y = -log10(FDR), x = rho), data = everything_df) + 
          geom_jitter(colour = everything_df$barcolour, width = 0.02, height = 0.02, size = 0.05, alpha = 0.6) +
          theme_classic() +
          theme(axis.text.y = element_text(colour = "black",
                                           size = 10),
                axis.text.x = element_text(colour = "black",
                                           size = 10),
                axis.title.x = element_text(size = 10),
                axis.title.y = element_text(size = 10),
                axis.line.y = element_line(colour = "black"),
                plot.margin = margin(t = 5, r = 10, b = 5, l = 5)) +
          geom_jitter(data = selection_df, width = 0.01, height = 0.06, colour = selection_df$barcolour, size = 0.7) + 
          # geom_label_repel(aes(label = hgnc_symbol), data = selection_df, size = 0.5, max.overlaps = 40) +
          coord_cartesian(xlim = xlim,
                          ylim = ylim) + 
          scale_y_continuous(breaks = yaxisbreaks,
                             limits = ylim,
                             expand = c(0,0)) +
          scale_x_continuous(breaks = xaxisbreaks,
                             limits = xlim,
                             expand = c(0,0),
                             labels = label_number(drop0trailing = TRUE)) +
          geom_vline(xintercept = 0, linetype = "dashed", colour = "grey") + 
          ylab(yaxislabel) +
          xlab(xaxislabel))
  
  dev.off()
  
}

reverselog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  trans_new(paste0("reverselog-", format(base)), trans, inv, 
            log_breaks(base = base), 
            domain = c(1e-100, Inf))
}

return.sites.over.peak.threshold <- function(chiptarget,
                                             sites_gr = genes,
                                             chosenthresh = 0.75) {
  
  counts_tissue_list <- lapply(ENCODEtissuetypes, function(thistissue){
    
    tissue_file_list <- list.files(paste0("input/", thistissue))
    tissue_files <- tissue_file_list[str_detect(tissue_file_list, "bed$")]
    
    tissue_chip_metadata <- read.table(paste0("input/", thistissue, "/", thistissue, "_chip_metadata.tsv"), 
                                       sep = "\t", 
                                       header = TRUE,
                                       quote = "",
                                       fill = TRUE)  
    
    tissue_CHIPtarget_metadata <- tissue_chip_metadata[str_detect(tissue_chip_metadata$Experiment.target, chiptarget), ]
    tissue_CHIPtarget_files <- tissue_files[str_remove(tissue_files, "\\.bed$") %in% tissue_CHIPtarget_metadata$File.accession]
    
    tissue_chip_data <- lapply(tissue_CHIPtarget_files, function(x){
      
      tempdata <- read.table(paste0("input/", thistissue, "/", x),
                             header = FALSE,
                             sep = "\t",
                             stringsAsFactors = FALSE,
                             quote = "")
      
      colnames(tempdata)[1:3] <- c("chr", "start", "end")
      
      return(tempdata)
      
    })
    
    patients <- str_remove(str_remove(tissue_chip_metadata[match(str_remove(tissue_CHIPtarget_files, pattern = "\\.bed$"), tissue_chip_metadata$File.accession), "Donor.s."], pattern = "\\/human-donors\\/"), pattern = "\\/$")
    names(tissue_chip_data) <- paste0(patients, "-", str_remove(tissue_chip_metadata[match(str_remove(tissue_CHIPtarget_files, pattern = "\\.bed$"), tissue_chip_metadata$File.accession), "Experiment.target"], pattern = "-human$"))
    
    tissue_chiptarget_gr_list <- lapply(tissue_chip_data, function(x){
      
      makeGRangesFromDataFrame(x)
      
    })
    
    counts <- sapply(tissue_chiptarget_gr_list, function(x){
      # this finds overlapping pairs of peaks and promoters
      # duplications appear on both sides i.e. peaks that overlap multiple promoters, promoters that overlap multiple peaks.
      # from here on we will take a promoter view
      overlaps <- countOverlaps(sites_gr, x)
      
    })
    
    return(counts)
    
  })
  
  counts_tissue_df <- do.call(cbind, counts_tissue_list)
  
  # restrict to present in RNA
  counts_tissue_df <- counts_tissue_df[, str_remove(colnames(counts_tissue_df), "-H.*") %in% str_remove(colnames(normalised_counts), "-.*")]
  
  number_of_samples <- ncol(counts_tissue_df)
  cutoff <- ceiling(chosenthresh * number_of_samples)
  
  peak_logical <- apply(counts_tissue_df, 1, function(x){x > 0})
  peak_totals <- colSums(peak_logical)
  
  sites_exceeding_threshold <- names(peak_totals)[peak_totals >= cutoff]
  
  return(sites_exceeding_threshold)
  
}

fit.models.to.chip <- function(computematrix_outputdirectory,
                               selectedtag,
                               removetag,
                               selectedsites,
                               genename,
                               gene_ensembl,
                               save,
                               filename){
  
  computematrixfiles <- list.files(paste0("input/", computematrix_outputdirectory))
  selectedfiles <- computematrixfiles[str_detect(computematrixfiles, selectedtag)]
  
  readmatrix <- lapply(selectedfiles, function(x){
    
    tempfile <- read.table(paste0("input/", computematrix_outputdirectory, "/", x),  skip = 1)
    
  })
  
  names(readmatrix) <- str_remove(selectedfiles, pattern = paste0(removetag, ".*\\.gz$"))
  
  readmatrix_df <- sapply(readmatrix, function(x){
    
    x$V7
    
  })
  
  row.names(readmatrix_df) <- readmatrix[[1]]$V4
  
  # restrict to samples in the metadata file
  readmatrix_df <- readmatrix_df[, sapply(colnames(readmatrix_df), function(x){
    str_detect(string = paste(chipmetadata_df$Files, sep = '', collapse = ''), pattern = x)}
  )]
  
  readmatrixt_df <- data.frame(t(readmatrix_df))
  row.names(readmatrixt_df) <- str_replace_all(paste0(sapply(1:nrow(readmatrixt_df), function(i){chipmetadata_df[str_detect(chipmetadata_df$Files, colnames(readmatrix_df)[i]), "Donor"]}),
                                                      "-",
                                                      sapply(1:nrow(readmatrixt_df), function(i){chipmetadata_df[str_detect(chipmetadata_df$Files, colnames(readmatrix_df)[i]), "Biosample.term.name"]})), pattern = " ", replacement = "_")
  
  matrix_for_lm <- readmatrixt_df
  
  matrix_for_lm_peaks <- matrix_for_lm[, str_remove(colnames(matrix_for_lm), "X") %in% selectedsites]
  
  matrix_for_lm_peaks[ , "tissue"] <- str_replace_all(sapply(1:nrow(readmatrixt_df), function(i){chipmetadata_df[str_detect(chipmetadata_df$Files, colnames(readmatrix_df)[i]), "Biosample.term.name"]}), pattern = " ", replacement = "_")
  matrix_for_lm_peaks[, "donor"] <- sapply(1:nrow(readmatrixt_df), function(i){chipmetadata_df[str_detect(chipmetadata_df$Files, colnames(readmatrix_df)[i]), "Donor"]})
  matrix_for_lm_peaks[, genename] <- normalised_counts[gene_ensembl, match(row.names(matrix_for_lm), colnames(normalised_counts) )]
  
  # exclude tissue types with fewer than three instances as cannot reasonably calculate residuals for those samples
  matrix_for_lm_peaks <- matrix_for_lm_peaks[!(matrix_for_lm_peaks$tissue %in% names(table(matrix_for_lm_peaks[, "tissue"]))[table(matrix_for_lm_peaks[, "tissue"]) < 3]), ]
  
  gene_lm <- lmer(formula = as.formula(paste0("log10(", genename, ")", " ~  + tissue + (1|donor)")), data = matrix_for_lm_peaks)
  matrix_for_lm_peaks[, paste0(genename, "resid")] <- residuals(gene_lm)
  
  temp_signal_model <- lapply(colnames(matrix_for_lm_peaks)[1:(ncol(matrix_for_lm_peaks) - 4)], function(x){
    
    tryCatch(expr = {
      templm <- lmer(formula = as.formula(paste0(x, " ~ ", genename, "resid + tissue + (1|donor)")), data = matrix_for_lm_peaks)
      templmsum <- summary(templm)
      
      t <- templmsum$coefficients[paste0(genename, "resid"), "t value"]
      
      output_vec <- c(t)
      names(output_vec) <- c("t_value")
      
      return(output_vec)
      
    }, error = function(e){
      
      output_vec <- c(NA)
      names(output_vec) <- c("t_value")
      
      return(output_vec)
      
    }
    )
    
  })
  
  names(temp_signal_model) <- colnames(matrix_for_lm_peaks)[1:(ncol(matrix_for_lm_peaks)-4)]
  
  if (save == TRUE) {
    
    saveRDS(temp_signal_model, paste0("output/", filename, ".rds"))
    
  }
  
  return(temp_signal_model)
  
}

spit.out.signal.matrix <- function(computematrix_outputdirectory,
                                   selectedtag,
                                   removetag,
                                   selectedsites = NULL
){
  
  computematrixfiles <- list.files(paste0("input/", computematrix_outputdirectory))
  selectedfiles <- computematrixfiles[str_detect(computematrixfiles, selectedtag)]
  
  readmatrix <- lapply(selectedfiles, function(x){
    
    tempfile <- read.table(paste0("input/", computematrix_outputdirectory, "/", x),  skip = 1)
    
  })
  
  names(readmatrix) <- str_remove(selectedfiles, pattern = paste0(removetag, ".*\\.gz$"))
  
  readmatrix_df <- sapply(readmatrix, function(x){
    
    x$V7
    
  })
  
  row.names(readmatrix_df) <- readmatrix[[1]]$V4
  
  # restrict to samples in the metadata file
  readmatrix_df <- readmatrix_df[, sapply(colnames(readmatrix_df), function(x){
    str_detect(string = paste(chipmetadata_df$Files, sep = '', collapse = ''), pattern = x)}
  )]
  
  readmatrixt_df <- data.frame(t(readmatrix_df))

  row.names(readmatrixt_df) <- str_replace_all(paste0(sapply(1:nrow(readmatrixt_df), function(i){chipmetadata_df[str_detect(chipmetadata_df$Files, colnames(readmatrix_df)[i]), "Donor"]}),
                                                      "-",
                                                      sapply(1:nrow(readmatrixt_df), function(i){chipmetadata_df[str_detect(chipmetadata_df$Files, colnames(readmatrix_df)[i]), "Biosample.term.name"]})), pattern = " ", replacement = "_")
  
  matrix_for_lm <- readmatrixt_df
  
  if(!is.null(selectedsites)){
    
  matrix_for_lm_peaks <- matrix_for_lm[, str_remove(colnames(matrix_for_lm), "X") %in% selectedsites]
  
  } else {
    
  matrix_for_lm_peaks <- matrix_for_lm  
    
  }
  
  return(matrix_for_lm_peaks)
  
}

chip.models.random.samples <- function(computematrix_outputdirectory,
                                       selectedtag,
                                       removetag,
                                       selectedsites,
                                       generandomsample,
                                       save,
                                       filename){
  
  computematrixfiles <- list.files(paste0("input/", computematrix_outputdirectory))
  selectedfiles <- computematrixfiles[str_detect(computematrixfiles, selectedtag)]
  
  readmatrix <- lapply(selectedfiles, function(x){
    
    tempfile <- read.table(paste0("input/", computematrix_outputdirectory, "/", x),  skip = 1)
    
  })
  
  names(readmatrix) <- str_remove(selectedfiles, pattern = paste0(removetag, ".*\\.gz$"))
  
  readmatrix_df <- sapply(readmatrix, function(x){
    
    x$V7
    
  })
  
  row.names(readmatrix_df) <- readmatrix[[1]]$V4
  
  # restrict to samples in the metadata file
  readmatrix_df <- readmatrix_df[, sapply(colnames(readmatrix_df), function(x){
    str_detect(string = paste(chipmetadata_df$Files, sep = '', collapse = ''), pattern = x)}
  )]
  
  readmatrixt_df <- data.frame(t(readmatrix_df))
  row.names(readmatrixt_df) <- str_replace_all(paste0(sapply(1:nrow(readmatrixt_df), function(i){chipmetadata_df[str_detect(chipmetadata_df$Files, colnames(readmatrix_df)[i]), "Donor"]}),
                                                      "-",
                                                      sapply(1:nrow(readmatrixt_df), function(i){chipmetadata_df[str_detect(chipmetadata_df$Files, colnames(readmatrix_df)[i]), "Biosample.term.name"]})), pattern = " ", replacement = "_")
  
  matrix_for_lm <- readmatrixt_df
  
  matrix_for_lm_peaks <- matrix_for_lm[, str_remove(colnames(matrix_for_lm), "X") %in% selectedsites]
  
  matrix_for_lm_peaks[ , "tissue"] <- str_replace_all(sapply(1:nrow(readmatrixt_df), function(i){chipmetadata_df[str_detect(chipmetadata_df$Files, colnames(readmatrix_df)[i]), "Biosample.term.name"]}), pattern = " ", replacement = "_")
  matrix_for_lm_peaks[, "donor"] <- sapply(1:nrow(readmatrixt_df), function(i){chipmetadata_df[str_detect(chipmetadata_df$Files, colnames(readmatrix_df)[i]), "Donor"]})
  
  # exclude tissue types with fewer than three instances as cannot reasonably calculate residuals
  matrix_for_lm_peaks <- matrix_for_lm_peaks[!(matrix_for_lm_peaks$tissue %in% names(table(matrix_for_lm_peaks[, "tissue"]))[table(matrix_for_lm_peaks[, "tissue"]) < 3]), ]
  
  for(i in 1:length(generandomsample)){
    
    matrix_for_lm_peaks[, paste(generandomsample[i])] <- normalised_counts[generandomsample[i], match(str_remove(row.names(matrix_for_lm_peaks), pattern = "-H.*"), colnames(normalised_counts))]
    
    temp_lm <- lmer(formula = paste0("log10(", generandomsample[i], " + 0.001) ~ tissue + (1|donor)"), data = matrix_for_lm_peaks)
    temp_resid <- residuals(temp_lm)
    
    matrix_for_lm_peaks[, paste0(generandomsample[i], "_resid")] <- temp_resid
    
  }
  
  sample_t_values <- sapply(generandomsample, function(thisgene){
    
    sapply(colnames(matrix_for_lm_peaks)[1:(ncol(matrix_for_lm_peaks) - 2 - (2 * length(generandomsample)))], function(x){
      
      tryCatch(expr = {
        
        glm <- lmer(formula = as.formula(paste0(x, " ~ ", thisgene, "_resid + tissue + (1|donor)")), data = matrix_for_lm_peaks)
        glmsum <- summary(glm)
        
        t <- glmsum$coefficients[paste0(thisgene, "_resid"), "t value"]
        
        output_vec <- c(t)
        names(output_vec) <- c("t_value")
        
        return(output_vec)
        
      }, error = function(e){
        
        output_vec <- c(NA)
        names(output_vec) <- c("t_value")
        
        return(output_vec)
        
      }
      )
      
    })
    
  })
  
  colnames(sample_t_values) <- generandomsample
  
  if(save == TRUE){
    
    saveRDS(sample_t_values, paste0("output/", filename, ".rds"))
    
  }
  
  return(sample_t_values)
  
}


promoter.peak.widths <- function(thistissue,
                                 chiptarget) {
  
  tissue_file_list <- list.files(paste0("input/", thistissue))
  tissue_files <- tissue_file_list[str_detect(tissue_file_list, "bed$")]
  
  tissue_chip_metadata <- read.table(paste0("input/", thistissue, "/", thistissue, "_chip_metadata.tsv"), 
                                     sep = "\t", 
                                     header = TRUE,
                                     quote = "",
                                     fill = TRUE)  
  
  tissue_CHIPtarget_metadata <- tissue_chip_metadata[str_detect(tissue_chip_metadata$Experiment.target, chiptarget), ]
  tissue_CHIPtarget_files <- tissue_files[str_remove(tissue_files, "\\.bed$") %in% tissue_CHIPtarget_metadata$File.accession]
  
  tissue_chip_data <- lapply(tissue_CHIPtarget_files, function(x){
    
    tempdata <- read.table(paste0("input/", thistissue, "/", x),
                           header = FALSE,
                           sep = "\t",
                           stringsAsFactors = FALSE,
                           quote = "")
    
    colnames(tempdata)[1:3] <- c("chr", "start", "end")
    
    return(tempdata)
    
  })
  
  patients <- str_remove(str_remove(tissue_chip_metadata[match(str_remove(tissue_CHIPtarget_files, pattern = "\\.bed$"), tissue_chip_metadata$File.accession), "Donor.s."], pattern = "\\/human-donors\\/"), pattern = "\\/$")
  names(tissue_chip_data) <- paste0(patients, "-", str_remove(tissue_chip_metadata[match(str_remove(tissue_CHIPtarget_files, pattern = "\\.bed$"), tissue_chip_metadata$File.accession), "Experiment.target"], pattern = "-human$"))
  
  tissue_chiptarget_gr_list <- lapply(tissue_chip_data, function(x){
    
    makeGRangesFromDataFrame(x)
    
  })
  
  tempprom <- promoters(genes(txdb_hg38))
  
  promoter_widths <- lapply(tissue_chiptarget_gr_list, function(x){

    # this finds overlapping pairs of peaks and promoters
    # duplications appear on both sides i.e. peaks that overlap multiple promoters, promoters that overlap multiple peaks.
    # from here on we will take a promoter view
    temppairs <- IRanges::findOverlapPairs(x, tempprom)
    
    # calculate width of peaks
    peak_widths <- width(S4Vectors::first(temppairs))
    
    # take out entrez gene ids for corresponding overlapping promoters
    promoter_entrez <- second(temppairs)$gene_id
    
    tempdf <- data.frame(peak_width = peak_widths, gene_promoter_entrez_id = promoter_entrez)
    
    # I aggregate it by summing the widths of all peaks overlapping the promoter.
    # this will give us one number per promoter
    agg_df <- aggregate(peak_width ~ gene_promoter_entrez_id, data = tempdf, FUN = sum)
    
    # here I have clear duplications in genes with consecutive identifiers and identical peak width
    # delete these duplicates (keep one) to not generate an artificial enrichment later on
    
    cens_agg_df <- data.frame(matrix(nrow = 1, ncol = 2))
    cens_agg_df[1, ] <- agg_df[1, ]
    
    for(i in 2:nrow(agg_df)){
      
      if(agg_df[i, 2] != agg_df[i-1, 2]){
        
        cens_agg_df[nrow(cens_agg_df) + 1, ] <- agg_df[i, ]
        
      }
      
    }
    
    colnames(cens_agg_df) <- colnames(agg_df)
    
    return(cens_agg_df)
    
  })
  
  # we add the tissue in the colname for later
  names(promoter_widths) <- paste0(names(promoter_widths), "-", thistissue)
  
  return(promoter_widths)
  
}


region.peak.counts <- function(NCI60_list, region = c("genes", "promoters", "repeats")){
  
  regioninput <- match.arg(region)
  
  if(regioninput == "genes"){
    chosenregion = genes
  }
  
  if(regioninput == "promoters"){
    chosenregion = gene_promoters
  }
  
  if(regioninput == "repeats"){
    chosenregion = rmask_gr
  }
  
  if(regioninput == "repeats"){
    
    chosenregion_peaks <- lapply(NCI60_list, function(x){
      
      temp_gr <- data.frame(seqnames = seqnames(chosenregion),
                            starts = start(chosenregion),
                            ends = end(chosenregion),
                            gene_id = chosenregion$myID)
      
      temp_gr$overlaps <- countOverlaps(chosenregion, x)
      
      return(temp_gr)
      
    })} else {
      
      chosenregion_peaks <- lapply(NCI60_list, function(x){
        
        temp_gr <- data.frame(seqnames = seqnames(chosenregion),
                              starts = start(chosenregion),
                              ends = end(chosenregion),
                              gene_id = chosenregion$gene_id)
        
        temp_gr$overlaps <- countOverlaps(chosenregion, x)
        
        return(temp_gr)
        
      })
    }
  
  if(regioninput == "repeats"){
    
    dfbase <- data.frame(seqnames = seqnames(chosenregion),
                         starts = start(chosenregion),
                         ends = end(chosenregion),
                         gene_id = chosenregion$myID)
    
    chosenregion_peaks_df <- cbind(dfbase, do.call(cbind, lapply(chosenregion_peaks, function(x){x$overlaps})))
    
    return(chosenregion_peaks_df)
    
  } else {
    
    dfbase <- data.frame(seqnames = seqnames(chosenregion),
                         starts = start(chosenregion),
                         ends = end(chosenregion),
                         gene_id = chosenregion$gene_id)
    
    chosenregion_peaks_df <- cbind(dfbase, do.call(cbind, lapply(chosenregion_peaks, function(x){x$overlaps})))
    
    return(chosenregion_peaks_df)
    
  }}

write.bed.file.for.computematrix <- function(sampleGR){
  
  temp_predbed_df <- data.frame(seqnames = seqnames(sampleGR),
                                starts = start(sampleGR)-1,
                                ends = end(sampleGR),
                                gene_id = names(sampleGR),
                                scores = c(rep(".", length(sampleGR))),
                                strands = strand(sampleGR))
  
  options(scipen=999)
  
  write.table(temp_predbed_df, file = paste0("output/", deparse(substitute(sampleGR)), ".bed"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  
  
}

plot.selection.against.everything.volcano.facet <- function(everything_object,
                                                            selection_ensembl, #a vector or list of vectors
                                                            selection_colours,  # a vector
                                                            filename = "myfacetvolcano",
                                                            xlim = NULL,
                                                            ylim = NULL,
                                                            plotheight = 11,
                                                            plotwidth = 8,
                                                            xaxislabel = substitute("Pearson's"~italic(r)),
                                                            yaxislabel = substitute(-log[10]~"(false discovery rate)"),
                                                            fileoutput = "pdf"){
  
  if(!is.list(selection_ensembl)){
    selection_ensembl <- list(selection_ensembl)
  }
  
  everything_object_ids_df <- lapply(everything_object, function(x){
    
    x[, "row_names"] <- row.names(x)
    x[, "FDR"] <- p.adjust(x$allpvals, method = "BH")
    x <- x[, !(colnames(x) == "allpvals")]
    return(x)
    
  })
  
  everything_melted <- melt(everything_object_ids_df)
  everything_melted <- everything_melted[!is.na(everything_melted$value), ]
  
  FDR_melt <- everything_melted[everything_melted$variable == "FDR", ]
  rho_melt <- everything_melted[everything_melted$variable == "allcorrs", ]
  
  meltyplot <- rho_melt
  meltyplot[, "FDR"] <- FDR_melt$value
  meltyplot <- meltyplot[, c(1,3,4,5)]
  colnames(meltyplot) <- c("ensembl_gene_id", "rho", "subtype", "FDR")
  
  meltyplot <- meltyplot[!(meltyplot$rho > 0.999), ]
  
  for(i in 1:length(unique(meltyplot$subtype))){
    
    FDRs_pos <- meltyplot[meltyplot$subtype == unique(meltyplot$subtype)[i] & meltyplot$rho > 0, "FDR"]
    FDRs_neg <- meltyplot[meltyplot$subtype == unique(meltyplot$subtype)[i] & meltyplot$rho < 0, "FDR"]
    
    FDRs_pos[FDRs_pos == 0] <- min(FDRs_pos[!(FDRs_pos == 0)])
    FDRs_neg[FDRs_neg == 0] <- min(FDRs_neg[!(FDRs_neg == 0)])
    
    meltyplot[meltyplot$subtype == unique(meltyplot$subtype)[i] & meltyplot$rho > 0, "FDR"] <- FDRs_pos
    meltyplot[meltyplot$subtype == unique(meltyplot$subtype)[i] & meltyplot$rho < 0, "FDR"] <- FDRs_neg
    
  }
  
  # not sure what i need this one for?? ah it wa a legacy from bar chart
  # meltyplot$ensembl_gene_id <- factor(meltyplot$ensembl_gene_id, levels = meltyplot[order(meltyplot$rho, decreasing = FALSE), "ensembl_gene_id"])
  meltyplot[, "barcolour"] <- "grey"
  meltyplot[meltyplot$rho < quantile(meltyplot$rho, 0.025, na.rm = TRUE), "barcolour"] <- "gray55"
  meltyplot[meltyplot$rho > quantile(meltyplot$rho, 0.975, na.rm = TRUE), "barcolour"] <- "gray55"
  
  selection_df <- meltyplot[meltyplot$ensembl_gene_id %in% unlist(selection_ensembl), ]
  
  for(i in 1:length(selection_ensembl)){
    
    selection_df[selection_df$ensembl_gene_id %in% selection_ensembl[[i]], "barcolour"] <- selection_colours[i]
    
  }
  
  # now make the plot below into volcano
  
  
  if(fileoutput == "pdf"){
    pdf(paste0("graphics/", filename, ".pdf"),
        height = plotheight,
        width = plotwidth)} else {
          png(paste0("graphics/", filename, ".png"),
              height = plotheight,
              width = plotwidth,
              units = "in",
              res = 300)
        }
  
  print(ggplot(aes(y = -log10(FDR), x = rho), data = meltyplot) + 
          geom_jitter(colour = meltyplot$barcolour, width = 0.02, height = 0.02, size = 0.05, alpha = 0.6) +
          theme_classic() +
          theme(axis.text.y = element_text(colour = "black"),
                axis.text.x = element_text(colour = "black"),
                axis.title.y = element_text(size = 10),
                axis.line.y = element_line(colour = "black")) +
          geom_jitter(data = selection_df, width = 0.01, height = 0.06, colour = selection_df$barcolour, size = 0.7) + 
          # geom_label_repel(aes(label = hgnc_symbol), data = selection_df, size = 0.5, max.overlaps = 40) +
          coord_cartesian(xlim = xlim,
                          ylim = ylim) +
          geom_vline(xintercept = 0, linetype = "dashed", colour = "grey") + 
          facet_wrap(~ subtype) +
          ylab(yaxislabel) +
          xlab(xaxislabel))
  
  dev.off()
  
}

plot.selection.against.everything.volcano.full.plots <- function(everything_object,
                                                            selection_ensembl, #a vector or list of vectors
                                                            selection_colours,  # a vector
                                                            filename_prefix = "myvolcano",
                                                            plotheight = 11,
                                                            plotwidth = 8,
                                                            xaxislabel = substitute("Pearson's"~italic(r)),
                                                            yaxislabel = substitute(-log[10]~"(false discovery rate)"),
                                                            fileoutput = "pdf"){
  
  if(!is.list(selection_ensembl)){
    selection_ensembl <- list(selection_ensembl)
  }
  
  everything_object_ids_df <- lapply(everything_object, function(x){
    
    x[, "row_names"] <- row.names(x)
    x[, "FDR"] <- p.adjust(x$allpvals, method = "BH")
    x <- x[, !(colnames(x) == "allpvals")]
    return(x)
    
  })
  
  everything_melted <- melt(everything_object_ids_df)
  everything_melted <- everything_melted[!is.na(everything_melted$value), ]
  
  FDR_melt <- everything_melted[everything_melted$variable == "FDR", ]
  rho_melt <- everything_melted[everything_melted$variable == "allcorrs", ]
  
  meltyplot <- rho_melt
  meltyplot[, "FDR"] <- FDR_melt$value
  meltyplot <- meltyplot[, c(1,3,4,5)]
  colnames(meltyplot) <- c("ensembl_gene_id", "rho", "subtype", "FDR")
  
  meltyplot <- meltyplot[!(meltyplot$rho > 0.999), ]
  
  for(i in 1:length(unique(meltyplot$subtype))){
    
    FDRs_pos <- meltyplot[meltyplot$subtype == unique(meltyplot$subtype)[i] & meltyplot$rho > 0, "FDR"]
    FDRs_neg <- meltyplot[meltyplot$subtype == unique(meltyplot$subtype)[i] & meltyplot$rho < 0, "FDR"]
    
    FDRs_pos[FDRs_pos == 0] <- min(FDRs_pos[!(FDRs_pos == 0)])
    FDRs_neg[FDRs_neg == 0] <- min(FDRs_neg[!(FDRs_neg == 0)])
    
    meltyplot[meltyplot$subtype == unique(meltyplot$subtype)[i] & meltyplot$rho > 0, "FDR"] <- FDRs_pos
    meltyplot[meltyplot$subtype == unique(meltyplot$subtype)[i] & meltyplot$rho < 0, "FDR"] <- FDRs_neg
    
  }
  
  # not sure what i need this one for?? ah it wa a legacy from bar chart
  # meltyplot$ensembl_gene_id <- factor(meltyplot$ensembl_gene_id, levels = meltyplot[order(meltyplot$rho, decreasing = FALSE), "ensembl_gene_id"])
  
  selection_df <- meltyplot[meltyplot$ensembl_gene_id %in% unlist(selection_ensembl), ]
  
  for(i in 1:length(selection_ensembl)){
    
    selection_df[selection_df$ensembl_gene_id %in% selection_ensembl[[i]], "barcolour"] <- selection_colours[i]
    
  }
  
  # now make the plot below into volcano

  for(i in 1:length(unique(meltyplot$subtype))){

    meltyplot_temp <- meltyplot[meltyplot$subtype == unique(meltyplot$subtype)[i], ]
    meltyplot_temp[, "barcolour"] <- "grey"
    meltyplot_temp[meltyplot_temp$rho < quantile(meltyplot_temp$rho, 0.025, na.rm = TRUE), "barcolour"] <- "gray55"
    meltyplot_temp[meltyplot_temp$rho > quantile(meltyplot_temp$rho, 0.975, na.rm = TRUE), "barcolour"] <- "gray55"
    
    selection_df_temp <- selection_df[selection_df$subtype == unique(meltyplot$subtype)[i], ]
    
    selection_df_temp[, "hgnc_symbol"] <- SET_HMTs_annotated[match(selection_df_temp$ensembl_gene_id, SET_HMTs_annotated$ensembl_gene_id), "hgnc_symbol"]
    selection_df_temp[(!selection_df_temp$rho < quantile(meltyplot_temp$rho, 0.025, na.rm = TRUE)) & (!selection_df_temp$rho > quantile(meltyplot_temp$rho, 0.975, na.rm = TRUE)), "hgnc_symbol"] <- ""
    
    selection_df_temp[selection_df_temp$ensembl_gene_id == "allHMTs", "hgnc_symbol"] <- "all HMTs"
    
    meltyplot_temp$title <- unique(meltyplot$subtype)[i]
    
  if(fileoutput == "pdf"){
    pdf(paste0("graphics/", filename_prefix, "_", unique(meltyplot$subtype)[i],".pdf"),
        height = plotheight,
        width = plotwidth)} else {
          png(paste0("graphics/", filename_prefix, "_", unique(meltyplot$subtype)[i], ".png"),
              height = plotheight,
              width = plotwidth,
              units = "in",
              res = 300)
        }
    
  print(ggplot(aes(y = -log10(FDR), x = rho), data = meltyplot_temp) + 
          geom_jitter(colour = meltyplot_temp$barcolour, width = 0.02, height = 0.02, size = 0.05, alpha = 0.6) +
          theme_classic() +
          theme(axis.text.y = element_text(colour = "black"),
                axis.text.x = element_text(colour = "black"),
                axis.title.y = element_text(size = 10),
                axis.line.y = element_line(colour = "black"),
                strip.text = element_text(size = 12, color = "white", face = "bold"),
                strip.background = element_rect(fill = "black")) +
          geom_jitter(data = selection_df_temp, width = 0.01, height = 0.06, colour = selection_df_temp$barcolour, size = 1) + 
          geom_label_repel(aes(label = hgnc_symbol), data = selection_df_temp, max.overlaps = 60) +
          geom_vline(xintercept = 0, linetype = "dashed", colour = "grey") +
          xlab(substitute("Spearman's"~italic(rho))) +
          ylab(substitute(-log[10]~"(FDR)")) +
    facet_grid(. ~ title))

  dev.off()
  
  }
  
}

#### LOAD PACKAGES & FUNCTIONS ####

## First specify the packages of interest
packages <- c("stringr", 
             "matrixStats", 
             "biomaRt", 
             "ggplot2", 
             "reshape2",
             "gtools",
             "ggrepel", # ggrepel for bubble plots
             "ggdendro", # ggdendro for metabolite clustering dendrograms
             "dendextend", # for dendrograms too
             "pals", # for metabolomics pathways colours
             "openxlsx",
             "dichromat", # for plot colours
             "ppcor", # ppcor for partial correlations
             "ggpubr", #ggpubr for theme_classic2 function,
             "recount", # for TCGA metadata
             "scales", # for logarithmic scale labelling without trailing 0s
             "enrichR",
             "lme4",
             "ProliferativeIndex",
             "readxl",
             "gplots",
             "corrplot", # for HMT correlation visualisation
             "Hmisc", # for p values from correlation matrices
             "car", # for LIN35 ANOVA
             "rgexf", # for preparing network plots for import into Gephi
             "EnvStats", # for co-regulation modelling
             "MASS",
             "igraph",
             "plyr",
             "CePa" # for read.gct
             )

## Now load or install&load all
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)

# Here define packages which need to be loaded through biocmanager

biocmanager_packages <- c("DESeq2",
                          "KEGGREST",
                          "GenomicRanges", # for intersecting
                          "GenomicAlignments",
                          "TxDb.Hsapiens.UCSC.hg38.knownGene",
                          "TxDb.Celegans.UCSC.ce11.refGene",
                          "TCGAbiolinks",
                          "recount",
                          "depmap",
                          "limma",
                          "GEOquery",
                          # "org.Hs.es.db",
                          "decoupleR",
                          "dorothea",
                          "memes",
                          "BSgenome.Celegans.UCSC.ce11",
                          "OmnipathR")

bioc_package.check <- lapply(
  biocmanager_packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      
      if (!requireNamespace("BiocManager", quietly = TRUE)){
        install.packages("BiocManager")
      }
      
      BiocManager::install(x, dependencies = TRUE)
      
      library(x, character.only = TRUE)
      
    }
  }
)

source("source/plot_ggdendro.R")

# here load RAPToR from Github

github_packages <- c("LBMC/RAPToR",
                     "LBMC/wormRef")

if (!requireNamespace("devtools", quietly = TRUE)){
  install.packages("devtools")}

gh_package.check <- lapply(
  github_packages,
  FUN = function(x) {
    if (!require(str_remove(x, ".*\\/"), character.only = TRUE)) {
      
      if (!requireNamespace("devtools", quietly = TRUE)){
        install.packages("devtools")}
      
      devtools::install_github(x, build_vignettes = TRUE)
      
      library(str_remove(x, ".*\\/"), character.only = TRUE)
      
    }
  }
)

#### INPUT DATA ####

# provide short names for genome objects
txdb_hg38 <- TxDb.Hsapiens.UCSC.hg38.knownGene
genes <- genes(txdb_hg38)

# repetitive elements file for hg38
rmask <- read.table("input/rmsk.txt", sep = "\t", header = FALSE, quote = "")
colnames(rmask)[c(6,7,8)] <- c("chr", "start", "end")
rmask[, "myID"] <- paste0(rmask$chr, rmask$start, rmask$end)
row.names(rmask) <- rmask$myID

# set ensembl MART
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# read in HMTs identities and annotations
SET_HMTs <- read.table("input/HMTs_list.txt", header = TRUE)
SET_HMTs_annotated <- read.table("input/HMTs_list_annotated.txt", header = TRUE, quote = "", sep = "\t")

act_rep_list <- list(activating_HMTs = SET_HMTs_annotated[SET_HMTs_annotated$Transcription == "activating", "ensembl_gene_id"],
                     repressing_HMTs = SET_HMTs_annotated[SET_HMTs_annotated$Transcription == "repressing", "ensembl_gene_id"],
                     unclear_HMTs = SET_HMTs_annotated[SET_HMTs_annotated$Transcription == "unclear", "ensembl_gene_id"],
                     allHMTs = "allHMTs")

act_rep_list_hgnc <- list(activating_HMTs = SET_HMTs_annotated[SET_HMTs_annotated$Transcription == "activating", "hgnc_symbol"],
                     repressing_HMTs = SET_HMTs_annotated[SET_HMTs_annotated$Transcription == "repressing", "hgnc_symbol"],
                     unclear_HMTs = SET_HMTs_annotated[SET_HMTs_annotated$Transcription == "unclear", "hgnc_symbol"],
                     allHMTs = "allHMTs")

# read in demethylase identitites and annotations
deMTs <- read.table("input/HdeMTs_list_annotated.txt", header = TRUE)

# save NNMT and PEMT common IDs for usage later
NNMT_ensembl <- "ENSG00000166741"
NNMT_entrezgene_id <- "4837"
PEMT_ensembl <- "ENSG00000133027"
PEMT_entrezgene_id <- "10400"

# read in CCLE metabolomics
metabolomics <- read.table("input/CCLE_metabolomics_20190502.csv",
                           sep = ",",
                           header = TRUE)

# transform the metabolomics data frame
# first eliminate one entry without a depmap id
metabolomics <- metabolomics[!(is.na(metabolomics$DepMap_ID)), ]

# then set depmap id as row names
row.names(metabolomics) <- metabolomics$DepMap_ID

# then eliminate first two columns (cell line name and ID)
metabolomics <- metabolomics[, 3:ncol(metabolomics)]

CCLE_proteomics <- read.table("input/Nusinov2020_Table_S2_Protein_Quant_Normalized.txt",
                              sep = "\t",
                              quote = "",
                              header = TRUE)

# C. elegans CeNDR data (supplementary files to GEO accession number GSE186719)

CeNDR_raw_counts <- read.table(file = "input/GSE186719_Celegans_208strains_609samples_rawCounts.tsv",
                         sep = "\t",
                         fill = TRUE,
                         quote = "",
                         header = TRUE)

CeNDR_TPMs <- read.table(file = "input/GSE186719_Celegans_208strains_609samples_rawTPM.tsv",
                   sep = "\t",
                   fill = TRUE,
                   quote = "",
                   header = TRUE)

#### generate nuclear genes expressed gold standard gene set ####

# read in expression data from GTEX v8 downloaded directly as transcript per million (TPM) values
# file accessed through GTEX data portal: [https://gtexportal.org/home/datasets]
# file can be downloaded here: [https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz]
TPMdata <- read.gct("input/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct")
row.names(TPMdata) <- str_remove(row.names(TPMdata), pattern = "\\..*$")

nuclear_genes_expressed <- row.names(TPMdata[matrixStats::rowMedians(TPMdata) > 5, ])

# retrieve ensembl gene id and HGNC gene symbols for genes on mitochondrial DNA from biomart
mitogenes <- getBM(mart = ensembl,
                   attributes = c("ensembl_gene_id", "hgnc_symbol"),
                   filters = "chromosome_name",
                   values = "MT")

nuclear_genes_expressed <- nuclear_genes_expressed[!(nuclear_genes_expressed %in% mitogenes$ensembl_gene_id)]

saveRDS(nuclear_genes_expressed, "output/GTEX_nuclear_genes_expressed.rds")
# nuclear_genes_expressed <- readRDS("output/GTEX_nuclear_genes_expressed.rds")

#### DOWNLOAD AND PROCESS TCGA AND GTEX RNASEQ DATA ####

#### QUERY DATA WITH TCGA BIOLINKS ####
# n.b. after recent update to TCGA biolinks, it is possible that this code as written would need modification

#project names downloaded from GDC Data Portal: [https://portal.gdc.cancer.gov/projects]

# change file name to match
proj_table <- read.table("input/projects-table.2020-12-03.tsv", sep="\t", header=TRUE)
proj_names <- proj_table$Project
proj_names <- proj_names[1:67]

# restrict to only TCGA projects
proj_names <- proj_names[str_detect(proj_names, pattern = "^TCGA-")]

dir.create("output/fulldata-counts")

# Query platform 
TCGA_counts_list <- list()

for (i in 1:length(proj_names)){
  
  favfile <- paste0("output/fulldata-counts/", proj_names[i], "-counts.rds")
  
  if (file.exists(favfile)) {
    
    message(paste0("We got this one (", proj_names[i], ") already"))
    next
  } else {
    
    tryCatch(
      expr = {
        query <- GDCquery(project = paste(proj_names[i]),
                          data.category = "Transcriptome Profiling",
                          data.type = "Gene Expression Quantification",
                          workflow.type = "HTSeq - Counts",
                          experimental.strategy = "RNA-Seq",
                          legacy = FALSE)
        
        GDCdownload(query, method = "api", files.per.chunk = 10)
        
        tempdata <- GDCprepare(query, summarizedExperiment = FALSE)
        
        tempdata_mat <- as.matrix(tempdata[,2:ncol(tempdata)])
        row.names(tempdata_mat) = sub("\\.[0-9]+$", "", unlist(tempdata[,1]))
        
        saveRDS(tempdata_mat, file = paste0("output/fulldata-counts/", proj_names[i] ,"-counts.rds"))
        
        TCGA_counts_list[[i]] <- tempdata_mat
        
        rm(query)
        rm(tempdata)
        rm(tempmatrix)
        rm(tempfav)
        
        gc()
        
      }, error = function(e) {
        message(paste("An error occurred for project ", proj_names[i], sep=""))
        print(e)
        
        rm(query)
        rm(tempdata)
        
        gc()
        
      }
    )
  }
}


# if count files already downloaded and saved as RDS #

# TCGA_counts_list <- list()

# TCGAcountsfilelist <- list.files(path = "output/fulldata-counts/")
#for(i in 1:length(TCGAcountsfilelist)){

#TCGA_counts_list[[i]] <- readRDS(paste0("output/fulldata-counts/", TCGAcountsfilelist[i]))
#}

names(TCGA_counts_list) <- str_remove(TCGAcountsfilelist, pattern = "-counts.rds")

saveRDS(TCGA_counts_list, "output/TCGA_counts_list.rds")
#  TCGA_counts_list <-readRDS("output/TCGA_counts_list.rds")

#### TCGA LOOKUP TABLES ####

# download TCGA metadata using 'recount' package
TCGA_metadata <- all_metadata(subset = "tcga", verbose = TRUE) 

# compile fields of interest into data frame
metadata_for_lm <- data.frame(cbind(TCGA_metadata$gdc_cases.samples.portions.analytes.aliquots.submitter_id,
                                    TCGA_metadata$gdc_cases.project.project_id,
                                    TCGA_metadata$gdc_cases.samples.sample_type,
                                    TCGA_metadata$gdc_cases.demographic.race,
                                    TCGA_metadata$gdc_cases.demographic.gender,
                                    TCGA_metadata$gdc_cases.diagnoses.tumor_stage,
                                    TCGA_metadata$gdc_cases.diagnoses.days_to_birth
))

# name columns
colnames(metadata_for_lm) <- c("SAMPID",
                               "cancer",
                               "sample_type",
                               "race",
                               "gender",
                               "tumour_stage",
                               "days_to_birth"
)

# Add field from barcode for sequencing centre
metadata_for_lm[, "sequencing_centre"] <- as.factor(str_extract(string = metadata_for_lm$SAMPID, pattern = "[:digit:]{2}$"))

# simplify tumour stage classification into stage 0, 1, 2, 3 or 4 (or not reported)
metadata_for_lm$tumour_stage <- str_replace_all(metadata_for_lm$tumour_stage, pattern = "stage i$", replacement = "stage_1")
metadata_for_lm$tumour_stage <- str_replace_all(metadata_for_lm$tumour_stage, pattern = "stage ia$", replacement = "stage_1")
metadata_for_lm$tumour_stage <- str_replace_all(metadata_for_lm$tumour_stage, pattern = "stage ib$", replacement = "stage_1")

metadata_for_lm$tumour_stage <- str_replace_all(metadata_for_lm$tumour_stage, pattern = "stage ii$", replacement = "stage_2")
metadata_for_lm$tumour_stage <- str_replace_all(metadata_for_lm$tumour_stage, pattern = "stage iia$", replacement = "stage_2")
metadata_for_lm$tumour_stage <- str_replace_all(metadata_for_lm$tumour_stage, pattern = "stage iib$", replacement = "stage_2")
metadata_for_lm$tumour_stage <- str_replace_all(metadata_for_lm$tumour_stage, pattern = "stage iic$", replacement = "stage_2")

metadata_for_lm$tumour_stage <- str_replace_all(metadata_for_lm$tumour_stage, pattern = "stage iii$", replacement = "stage_3")
metadata_for_lm$tumour_stage <- str_replace_all(metadata_for_lm$tumour_stage, pattern = "stage iiia$", replacement = "stage_3")
metadata_for_lm$tumour_stage <- str_replace_all(metadata_for_lm$tumour_stage, pattern = "stage iiib$", replacement = "stage_3")
metadata_for_lm$tumour_stage <- str_replace_all(metadata_for_lm$tumour_stage, pattern = "stage iiic$", replacement = "stage_3")

metadata_for_lm$tumour_stage <- str_replace_all(metadata_for_lm$tumour_stage, pattern = "stage iv$", replacement = "stage_4")
metadata_for_lm$tumour_stage <- str_replace_all(metadata_for_lm$tumour_stage, pattern = "stage iva$", replacement = "stage_4")
metadata_for_lm$tumour_stage <- str_replace_all(metadata_for_lm$tumour_stage, pattern = "stage ivb$", replacement = "stage_4")
metadata_for_lm$tumour_stage <- str_replace_all(metadata_for_lm$tumour_stage, pattern = "stage ivc$", replacement = "stage_4")

#Keep only primary tumour samples; exclude stages i/ii nos, is, stage 0 and stage x; latter excludes only 84 samples across all TCGA
TCGA_LUT <- metadata_for_lm[metadata_for_lm$sample_type %in% c("Primary Tumor","Primary Blood Derived Cancer - Peripheral Blood"), ]
TCGA_LUT <- TCGA_LUT[!(TCGA_LUT$tumour_stage %in% c("i/ii nos", "is", "stage x")), ]

# remove duplicate entries (duplicated sample ID)
TCGA_LUT <- TCGA_LUT[!(duplicated(TCGA_LUT$SAMPID)),]

# remove entries where multiple samples are derived from one donor. keep only donors that provide one cancer sample.
TGA_LUT <- TCGA_LUT[(!duplicated(str_extract(string = TCGA_LUT$SAMPID, pattern = "^TCGA-[:alnum:]{2}-[:alnum:]{4}")) & !duplicated(str_extract(string = TCGA_LUT$SAMPID, pattern = "^TCGA-[:alnum:]{2}-[:alnum:]{4}"), fromLast = TRUE)), ]

# set SAMPID as row names
row.names(TCGA_LUT) <- TCGA_LUT$SAMPID

saveRDS(TCGA_LUT, "output/TCGA_LUT.rds")
# TCGA_LUT <- readRDS("output/TCGA_LUT.rds")

cancertypes <- str_sort(unique(TCGA_LUT$cancer))

#### TCGA RNA-SEQ MRN NORMALISATION ####

TCGA_MOR_list <- list()

for (i in 1:length(cancertypes)){
  
  # retrieve sample IDs for this tissue from the lookup table. use metadata_for_lm to not restrict to cancer samples; include normal samples.
  sampleIDs <- metadata_for_lm[metadata_for_lm$cancer == cancertypes[i], "SAMPID"]
  
  # subset the dataframe of raw counts by sample ids corresponding to tissue i
  cancercounts <- TCGA_counts_list[[paste(cancertypes[i])]][, colnames(TCGA_counts_list[[paste(cancertypes[i])]]) %in% sampleIDs]
  
  # DESeq2 wants a colData object. Not actually used for the normalisation. Here we can use the sample IDs with tissue.
  col_cancer <- metadata_for_lm[colnames(cancercounts), "cancer"]
  col_cancer <- as.matrix(col_cancer)
  
  rownames(col_cancer) <- colnames(cancercounts)
  colnames(col_cancer) <- c("cancer")
  
  # need to create a DESeq2 object. Design set to ~1 allows for use of estimateSizeFactors
  tempdds <- DESeqDataSetFromMatrix(countData = cancercounts, colData = col_cancer, design = ~ 1)
  
  # this function estimates the scaling factors from the samples from the median of ratios wrt to the geometric mean for each gene across samples
  tempdds <- estimateSizeFactors(tempdds)
  
  # put the counts normalised by the scaling factors in a new object
  cancernormalised_counts <- counts(tempdds, normalized = TRUE)
  
  # put object in list
  TCGA_MOR_list[[i]] <- cancernormalised_counts
  
}

# set the names in the list to be the tissue names
names(TCGA_MOR_list) <- cancertypes

# save the list in output folder
saveRDS(TCGA_MOR_list, "output/TCGA_MOR_list.rds")
# TCGA_MOR_list <- readRDS("output/TCGA_MOR_list.rds")

# perform MOR normalisation across all cancers (for later correlation with e.g. Proliferative Index)

# DESeq2 wants a colData object. Not actually used for the normalisation. Here we can use the sample IDs with tissue.
col_data <- TCGA_LUT[unlist(sapply(TCGA_counts_list, colnames)), "cancer"]

col_data <- as.matrix(col_data)

row.names(col_data) <- unlist(sapply(TCGA_counts_list, colnames))
colnames(col_data) <- c("cancer")

TCGA_counts_df <- do.call(cbind, TCGA_counts_list)

col_data <- col_data[colnames(TCGA_counts_df), ]

col_data <- col_data[!is.na(col_data)]
col_data <- data.frame(col_data)
TCGA_counts_df <- TCGA_counts_df[, colnames(TCGA_counts_df) %in% row.names(col_data)]

# Normalise reads by scaling factors using DESeq2

# need to create a DESeq2 object. Design set to ~1 allows for use of estimateSizeFactors
dds <- DESeqDataSetFromMatrix(countData = TCGA_counts_df, colData = col_data, design = ~ 1)

# this function estimates the scaling factors from the samples from the median of ratios wrt to the geometric mean for each gene across samples
dds <- estimateSizeFactors(dds)

# save the counts normalised by the scaling factors in a new object
TCGA_MOR_across_cancers <- counts(dds, normalized = TRUE)

# save it. this is for all genes. note this normalisation across all samples from all tissues
saveRDS(TCGA_MOR_across_cancers, "output/TCGA-MOR-normalisation-across-cancers.rds")
# TCGA_MOR_across_cancers <- readRDS("~/manuscript/output/TCGA-MOR-normalisation-across-cancers.rds")

#### GTEX IMPORT RNA-SEQ DATA ####

# read in expression data from GTEX v8 downloaded as read counts
# file accessed through GTEX data portal: [https://gtexportal.org/home/datasets]
# file can be downloaded here: [https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz]
countsdata <- read.gct("input/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct")
row.names(countsdata) <- str_remove(row.names(countsdata), pattern = "\\..*$")

# save as RDS for rapid loading later if code is run in separate sessions
saveRDS(countsdata, "output/countsdata.rds")
# countsdata <- readRDS("output/countsdata.rds")

#### GTEX LOOKUP TABLE ####

# create lookup table for GTEX v8
# load sample annotation and phenotype files and make lookup table for samples and tissues for later use

# file accessed through GTEX data portal: [https://gtexportal.org/home/datasets]
# file can be downloaded here: [https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt]
sampleattributes <- read.table("input/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt", sep = "\t", header = TRUE, quote = "")
row.names(sampleattributes) <- sampleattributes$SAMPID

# file accessed through GTEX data portal: [https://gtexportal.org/home/datasets]
# file can be downloaded here: [https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt]
subjectphenotypes <- read.table("input/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt", sep = "\t", header = TRUE, quote = "")
row.names(subjectphenotypes) <- subjectphenotypes$SUBJID

# Create new data frame in which to build our lookup table for the samples present in the expression data table
GTEX_LUT <- as.data.frame(matrix(nrow = ncol(TPMdata), ncol = 2))
colnames(GTEX_LUT) <- c("SUBJID", "SAMPID")

# for the purposes of constructing the lookup table, will use column names from the TPM data
# this code recreates the SUBJID from the column/sample names in the RPKMtable
# first splitting by the separator (in this case, ".")
sampid_split_list <- str_split(colnames(TPMdata), pattern = "\\.")

# then pasting back the first two fields with a "-" separator to get the subject ID, as SUBJID appears in sample annotation file.
SUBJID_vec <- c()

for(i in 1:length(sampid_split_list)){
  SUBJID_vec[i] <- paste(sampid_split_list[[i]][1:2], collapse = "-")  
}

# fill SUBJID column with SUBJIDs recreated from RPKM file
GTEX_LUT$SUBJID <- SUBJID_vec

# fill SAMPID column with sample IDs from RPKM file
GTEX_LUT$SAMPID <- colnames(TPMdata)

# add age bracket as column to look up table
GTEX_LUT[, "age_bracket"] <- subjectphenotypes[GTEX_LUT$SUBJID, "AGE"]

# set row names as sample IDs to allow for easy subsetting later
rownames(GTEX_LUT) <- GTEX_LUT$SAMPID

# add tissue as column to look up table (in sample attributes SMTSD column)
# note str_replace command is needed here to match the sampleattributes file as SAMPID from TPM table has "." separator, not "-"
GTEX_LUT[, "Tissue"] <- sampleattributes[str_replace_all(GTEX_LUT$SAMPID, pattern = "\\.", replacement = "-"), "SMTSD"]

# add gender as column to look up table 
GTEX_LUT[, "Sex"] <- subjectphenotypes[GTEX_LUT$SUBJID, "SEX"]

# replace gender codes (0,1) with appropriate strings
GTEX_LUT[, "Sex"] <- str_replace_all(as.character(GTEX_LUT$Sex), pattern = "1", replacement = "MALE")
GTEX_LUT[, "Sex"] <- str_replace_all(as.character(GTEX_LUT$Sex), pattern = "2", replacement = "FEMALE")

# add Hardy deathscale ('cause of death') to lookup table
GTEX_LUT[, "DeathScale"] <- subjectphenotypes[GTEX_LUT$SUBJID, "DTHHRDY" ]

# add sample ischemic time to lookup table
GTEX_LUT[, "IschemicTime"] <- sampleattributes[str_replace_all(GTEX_LUT$SAMPID, pattern = "\\.", replacement = "-"), "SMTSISCH"]

# add sample sequencing batch to lookup table
GTEX_LUT[, "Batch"] <- sampleattributes[str_replace_all(GTEX_LUT$SAMPID, pattern = "\\.", replacement = "-"), "SMGEBTCH"]

# save the lookup table in the output folder
saveRDS(GTEX_LUT, file = "output/GTEX-version8-sampleID-LUT.rds")
# GTEX_LUT <- readRDS("output/GTEX-version8-sampleID-LUT.rds")

# Make vector of tissue names to use downstream
# limited to tissue types with >100 samples

tissues <- names(table(LUT$Tissue))[table(LUT$Tissue) > 100]

#### GTEX RNA-SEQ MRN NORMALISATION ####

# perform MRN normalisation by tissue
# NB in this script the MRN normalisation is referred to as 'MOR' (median of ratios) instead of MRN (median ratio normalisation) as in the accompanying manuscript
# MRN/MOR normalisation uses 'DESeq2' package

GTEX_MOR_list <- list()

for (i in 1:length(tissues)){
  
  # retrieve sample IDs for this tissue from the lookup table
  sampleIDs <- LUT[LUT$Tissue == tissues[i], "SAMPID"]
  
  # subset the dataframe of raw counts by sample ids corresponding to tissue i
  tissuecounts <- countsdata[, colnames(countsdata) %in% sampleIDs]
  
  # DESeq2 wants a colData object. Not actually used for the normalisation. Here we can use the sample IDs with tissue.
  col_tissue <- LUT[colnames(tissuecounts), "Tissue"]
  col_tissue <- as.matrix(col_tissue)
  
  rownames(col_tissue) <- colnames(tissuecounts)
  colnames(col_tissue) <- c("Tissue")
  
  # need to create a DESeq2 object. Design set to ~1 allows for use of estimateSizeFactors
  tempdds <- DESeqDataSetFromMatrix(countData = tissuecounts, colData = col_tissue, design = ~ 1)
  
  # this function estimates the scaling factors from the samples from the median of ratios wrt to the geometric mean for each gene across samples
  tempdds <- estimateSizeFactors(tempdds)
  
  # put the counts normalised by the scaling factors in a new object
  tissuenormalised_counts <- counts(tempdds, normalized = TRUE)
  
  # put object in list
  GTEX_MOR_list[[i]] <- tissuenormalised_counts
  
}

# set the names in the list to be the tissue names
names(GTEX_MOR_list) <- tissues

# save the tissue MOR-normalised values in the output folder
saveRDS(GTEX_MOR_list, "output/GTEX_MOR_list.rds")
# GTEX_MOR_list <- readRDS("output/GTEX_MOR_list.rds")

# perform MOR normalisation across all tissues combined, rather than individually (we will use this for combined analysis and for correlating to NFKB/immune cell fraction/Proliferative Index)

# DESeq2 wants a colData object. Not actually used for the normalisation. Here we can use the sample IDs with tissue.
col_data <- LUT[colnames(countsdata), "Tissue"]
col_data <- as.matrix(col_data)

rownames(col_data) <- colnames(countsdata)
colnames(col_data) <- c("Tissue")

# Normalise reads by scaling factors using DESeq2

# need to create a DESeq2 object. Design set to ~1 allows for use of estimateSizeFactors
dds <- DESeqDataSetFromMatrix(countData = countsdata, colData = col_data, design = ~ 1)

# this function estimates the scaling factors from the samples from the median of ratios wrt to the geometric mean for each gene across samples
dds <- estimateSizeFactors(dds)

# save the counts normalised by the scaling factors in a new object
GTEX_MOR_across_tissues <- counts(dds, normalized = TRUE)

saveRDS(GTEX_MOR_across_tissues, "output/GTEX_MOR_across_tissues.rds")
GTEX_MOR_across_tissues <- readRDS("output/GTEX_MOR_across_tissues.rds")

#### NCI60 RNA-seq read counts ####


#### use Genomic Alignments to count RNA-seq reads ####

exbygene <- exonsBy(TxDb.Hsapiens.UCSC.hg38.knownGene, by = "gene")

# point to directory with aligned sequencing bam files
thebams <- list.files("insert_bam_directory", full = TRUE)

olap <- summarizeOverlaps(exbygene, thebams)
deseq <- DESeqDataSet(olap, design= ~ 1)

deseq <- estimateSizeFactors(deseq)
NCI60_normalised_counts <- counts(deseq, normalized = TRUE)

NCI60RNA_metadata <- read.table("input/NCI60_RNAseq_SRA.csv", sep = ",", header = TRUE)
colnames(NCI60_normalised_counts) <- NCI60RNA_metadata[match(str_remove(colnames(NCI60_normalised_counts), pattern = "\\.bam"), NCI60RNA_metadata$Run), "LibraryName"]

names_changed_to_match_format <- str_remove(str_remove(str_replace(str_replace(str_remove(colnames(NCI60_normalised_counts), "^[A-Z]{2,3}:"), "/", "_"), " ", "_"), "\\(TB\\)"), "_ATCC")
names_changed_to_match_format <- str_replace(names_changed_to_match_format, pattern = "MDA-N", replacement = "MDA-MB-468")

colnames(NCI60_normalised_counts) <- names_changed_to_match_format

saveRDS(NCI60_normalised_counts, "output/NCI60_RNAseq_normalised_counts.rds")
# NCI60_normalised_counts <- readRDS("output/NCI60_RNAseq_normalised_counts_old.rds")

#### ANNOTATE METABOLITES ####

# annotate metabolites

metabolites_for_KEGG <- colnames(metabolomics)

metabolites_for_KEGG <- str_remove(metabolites_for_KEGG, pattern = "^X[0-9]{1}\\.")
metabolites_for_KEGG <- str_replace_all(metabolites_for_KEGG, pattern = "\\.", replacement = " ")
metabolites_for_KEGG <- str_replace_all(metabolites_for_KEGG, pattern = "_", replacement = " ")

metabolites_for_KEGG <- str_replace(metabolites_for_KEGG, pattern = "^C[0-9]{2} [0-9]{1} ", replacement = "")
metabolites_for_KEGG <- str_replace_all(metabolites_for_KEGG, pattern = "TAG", replacement = "triacylglycerol")
metabolites_for_KEGG <- str_replace_all(metabolites_for_KEGG, pattern = "LPC", replacement = "lysophosphatidylcholine")
metabolites_for_KEGG <- str_replace_all(metabolites_for_KEGG, pattern = "PC", replacement = "phosphatidylcholine")
metabolites_for_KEGG <- str_replace_all(metabolites_for_KEGG, pattern = "DAG", replacement = "diacylglycerol")
metabolites_for_KEGG <- str_replace_all(metabolites_for_KEGG, pattern = "CE", replacement = "cholesterol ester")
metabolites_for_KEGG <- str_replace_all(metabolites_for_KEGG, pattern = "SM", replacement = "sphingomyelin")
metabolites_for_KEGG <- str_replace_all(metabolites_for_KEGG, pattern = "LPE", replacement = "lysophosphatidylethanolamine")

KEGGfind_list <- lapply(metabolites_for_KEGG, function(x){
  keggFind("compound", x)
})

names(KEGGfind_list) <- metabolites_for_KEGG

KEGG_metabolites_df <- data.frame(cbind(colnames(metabolomics), metabolites_for_KEGG))
KEGG_metabolites_df[, "KEGG_result"] <- sapply(KEGGfind_list, function(x){
  
  if(length(x) == 1){
    return(x)
  } else {
    return("")
  }
})

KEGG_metabolites_df[, "KEGG_ID"] <- sapply(KEGGfind_list, function(x){
  
  if(length(x) == 1){
    return(names(x))
  } else {
    return("")
  }
})

KEGGchooselist <- list(c("aminoadipate", 1),
                       c("phosphoglycerate", 2),
                       c("pyridoxate", 1),
                       c("aconitate", 1),
                       c("adenine", 6),
                       c("adipate", 14),
                       c("AMP", 1),
                       c("citrate", 1),
                       c("isocitrate", 1),
                       c("CMP", 1),
                       c("cystathionine", 1),
                       c("cytidine", 8),
                       c("dCMP", 1),
                       c("erythrose 4 phosphate", 1),
                       c("glucuronate", 2),
                       c("GMP", 1),
                       c("guanosine", 8),
                       c("hypoxanthine", 1),
                       c("inosine", 4),
                       c("kynurenine", 2),
                       c("lactate", 6),
                       c("lactose", 3),
                       c("malate", 3),
                       c("NAD", 1),
                       c("NADP", 2),
                       c("oxalate", 2),
                       c("pantothenate", 1),
                       c("sorbitol", 1),
                       c("sucrose", 1),
                       c("thymine", 1),
                       c("UMP", 1),
                       c("uracil", 1),
                       c("urate", 1), 
                       c("uridine", 5),
                       c("xanthine", 2),
                       c("hydroxyglutarate", 1),
                       c("inositol", 1),
                       c("glycine", 1),
                       c("alanine", 18),
                       c("serine", 3),
                       c("threonine", 1),
                       c("methionine", 4),
                       c("aspartate", 22),
                       c("glutamate", 3),
                       c("asparagine", 9),
                       c("glutamine", 2),
                       c("histidine", 4),
                       c("arginine", 10),
                       c("lysine", 71),
                       c("valine", 8),
                       c("leucine", 17),
                       c("isoleucine", 3),
                       c("phenylalanine", 3),
                       c("tyrosine", 9),
                       c("tryptophan", 4),
                       c("proline", 23),
                       c("ornithine", 4),
                       c("citrulline", 1),
                       c("taurine", 1),
                       c("serotonin", 1),
                       c("GABA", 1),
                       c("dimethylglycine", 1),
                       c("homocysteine", 2),
                       c("allantoin", 1),
                       c("anthranilic acid", 1),
                       c("kynurenic acid", 1),
                       c("carnosine", 1),
                       c("N carbamoyl beta alanine", 1),
                       c("thiamine", 2),
                       c("betaine", 4),
                       c("choline", 1),
                       c("acetylcholine", 1),
                       c("creatine", 1),
                       c("thyroxine", 1),
                       c("adenosine", 9),
                       c("thymidine", 1),
                       c("xanthosine", 4),
                       c("deoxyadenosine", 4),
                       c("deoxycytidine", 4),
                       c("cAMP", 1),
                       c("pipecolic acid", 1),
                       c("pyroglutamic acid", 1),
                       c("butyrobetaine", 3),
                       c("putrescine", 1),
                       c("carnitine", 1),
                       c("sarcosine", 1),
                       c("beta alanine", 2),
                       c("phosphatidylcholine", 1),
                       c("phosphatidylcholine A", 1),
                       c("diacylglycerol", 1),
                       c("cholesterol ester", 4),
                       c("lysophosphatidylcholine", 1),
                       c("lysophosphatidylethanolamine", 1))

for(i in 1:length(KEGGchooselist)){
  
  tempone <- KEGGchooselist[[i]]
  
  KEGG_metabolites_df[KEGG_metabolites_df$metabolites_for_KEGG == tempone[1], "KEGG_ID"] <- names(KEGGfind_list[[tempone[1]]][as.numeric(tempone[2])])
  KEGG_metabolites_df[KEGG_metabolites_df$metabolites_for_KEGG == tempone[1], "KEGG_result"] <- KEGGfind_list[[tempone[1]]][as.numeric(tempone[2])]
  
  
}

carnitine <- keggFind("compound", "carnitine")
KEGG_metabolites_df[str_detect(KEGG_metabolites_df$metabolites_for_KEGG, "valerylcarnitine"), "KEGG_ID"] <- names(carnitine[12])
KEGG_metabolites_df[str_detect(KEGG_metabolites_df$metabolites_for_KEGG, "valerylcarnitine"), "KEGG_result"] <- carnitine[12]

alpha_ketoglutarate <- keggFind("compound", "ketoglutaric")
KEGG_metabolites_df[str_detect(KEGG_metabolites_df$metabolites_for_KEGG, "alpha.ketoglutarate"), "KEGG_ID"] <- names(alpha_ketoglutarate[1])
KEGG_metabolites_df[str_detect(KEGG_metabolites_df$metabolites_for_KEGG, "alpha.ketoglutarate"), "KEGG_result"] <- alpha_ketoglutarate[1]

alpha_glycerophosphate <- keggFind("compound", "glycerophosphoric")
KEGG_metabolites_df[str_detect(KEGG_metabolites_df$metabolites_for_KEGG, "alpha.glycerophosphate"), "KEGG_ID"] <- names(alpha_glycerophosphate [1])
KEGG_metabolites_df[str_detect(KEGG_metabolites_df$metabolites_for_KEGG, "alpha.glycerophosphate"), "KEGG_result"] <- alpha_glycerophosphate [1]

glyceraldeyde_3P <- keggFind("compound", "glycerophosphoric")
KEGG_metabolites_df[str_detect(KEGG_metabolites_df$metabolites_for_KEGG, "DHAP.glyceraldehyde"), "KEGG_ID"] <- names(glyceraldeyde_3P [1])
KEGG_metabolites_df[str_detect(KEGG_metabolites_df$metabolites_for_KEGG, "DHAP.glyceraldehyde"), "KEGG_result"] <- glyceraldeyde_3P [1]

F1P <- keggFind("compound", "fructose phosphate")
KEGG_metabolites_df[str_detect(KEGG_metabolites_df$metabolites_for_KEGG, "F1P"), "KEGG_ID"] <- names(F1P[4])
KEGG_metabolites_df[str_detect(KEGG_metabolites_df$metabolites_for_KEGG, "F1P"), "KEGG_result"] <- F1P[4]

hexose <- keggFind("compound", "hexose")
KEGG_metabolites_df[str_detect(KEGG_metabolites_df$metabolites_for_KEGG, "hexose"), "KEGG_ID"] <- names(hexose[2])
KEGG_metabolites_df[str_detect(KEGG_metabolites_df$metabolites_for_KEGG, "hexose"), "KEGG_result"] <- hexose[2]

fumarate <- keggFind("compound", "fumarate")
KEGG_metabolites_df[str_detect(KEGG_metabolites_df$metabolites_for_KEGG, "fumarate"), "KEGG_ID"] <- names(fumarate[1])
KEGG_metabolites_df[str_detect(KEGG_metabolites_df$metabolites_for_KEGG, "fumarate"), "KEGG_result"] <- fumarate[1]

PEP <- keggFind("compound", "phosphoenolpyruvic")
KEGG_metabolites_df[str_detect(KEGG_metabolites_df$metabolites_for_KEGG, "PEP"), "KEGG_ID"] <- names(PEP[1])
KEGG_metabolites_df[str_detect(KEGG_metabolites_df$metabolites_for_KEGG, "PEP"), "KEGG_result"] <- PEP[1]

ribuloseP <- keggFind("compound", "ribulose phosphate")
KEGG_metabolites_df[str_detect(KEGG_metabolites_df$metabolites_for_KEGG, "ribulose"), "KEGG_ID"] <- names(ribuloseP[1])
KEGG_metabolites_df[str_detect(KEGG_metabolites_df$metabolites_for_KEGG, "ribulose"), "KEGG_result"] <- ribuloseP[1]

succinate <- keggFind("compound", "succinate")
KEGG_metabolites_df[str_detect(KEGG_metabolites_df$metabolites_for_KEGG, "succinate"), "KEGG_ID"] <- names(succinate[1])
KEGG_metabolites_df[str_detect(KEGG_metabolites_df$metabolites_for_KEGG, "succinate"), "KEGG_result"] <- succinate[1]

UDP_galactose <- keggFind("compound", "UDP galactose")
KEGG_metabolites_df[str_detect(KEGG_metabolites_df$metabolites_for_KEGG, "UDP.galactose"), "KEGG_ID"] <- names(UDP_galactose[1])
KEGG_metabolites_df[str_detect(KEGG_metabolites_df$metabolites_for_KEGG, "UDP.galactose"), "KEGG_result"] <- UDP_galactose[1]

glycodeoxycholate <- keggFind("compound", "glycodeoxycholate")
KEGG_metabolites_df[str_detect(KEGG_metabolites_df$metabolites_for_KEGG, "glycodeoxycholate"), "KEGG_ID"] <- names(glycodeoxycholate[1])
KEGG_metabolites_df[str_detect(KEGG_metabolites_df$metabolites_for_KEGG, "glycodeoxycholate"), "KEGG_result"] <- glycodeoxycholate[1]

taurodeoxycholate <- keggFind("compound", "taurodeoxycholate")
KEGG_metabolites_df[str_detect(KEGG_metabolites_df$metabolites_for_KEGG, "taurodeoxycholate"), "KEGG_ID"] <- names(taurodeoxycholate[1])
KEGG_metabolites_df[str_detect(KEGG_metabolites_df$metabolites_for_KEGG, "taurodeoxycholate"), "KEGG_result"] <- taurodeoxycholate[1]

pimelate <- keggFind("compound", "pimelate")
KEGG_metabolites_df[str_detect(KEGG_metabolites_df$metabolites_for_KEGG, "pimelate"), "KEGG_ID"] <- names(pimelate[3])
KEGG_metabolites_df[str_detect(KEGG_metabolites_df$metabolites_for_KEGG, "pimelate"), "KEGG_result"] <- pimelate[3]

hydroxybutyrate <- keggFind("compound", "hydroxybutyrate")
KEGG_metabolites_df[str_detect(KEGG_metabolites_df$metabolites_for_KEGG, "hydroxybutyrate"), "KEGG_ID"] <- names(hydroxybutyrate[1])
KEGG_metabolites_df[str_detect(KEGG_metabolites_df$metabolites_for_KEGG, "hydroxybutyrate"), "KEGG_result"] <- hydroxybutyrate[1]

hydroxybutyrate <- keggFind("compound", "hydroxybutyrate")
KEGG_metabolites_df[str_detect(KEGG_metabolites_df$metabolites_for_KEGG, "hydroxybutyrate"), "KEGG_ID"] <- names(hydroxybutyrate[1])
KEGG_metabolites_df[str_detect(KEGG_metabolites_df$metabolites_for_KEGG, "hydroxybutyrate"), "KEGG_result"] <- hydroxybutyrate[1]

hydroxyproline <- keggFind("compound", "hydroxyproline")
KEGG_metabolites_df[str_detect(KEGG_metabolites_df$metabolites_for_KEGG, "hydroxyproline"), "KEGG_ID"] <- names(hydroxyproline[2])
KEGG_metabolites_df[str_detect(KEGG_metabolites_df$metabolites_for_KEGG, "hydroxyproline"), "KEGG_result"] <- hydroxyproline[2]

hydroxyindoleacetate <- keggFind("compound", "hydroxyindoleacetate")
KEGG_metabolites_df[str_detect(KEGG_metabolites_df$metabolites_for_KEGG, "HIAA"), "KEGG_ID"] <- names(hydroxyindoleacetate[1])
KEGG_metabolites_df[str_detect(KEGG_metabolites_df$metabolites_for_KEGG, "HIAA"), "KEGG_result"] <- hydroxyindoleacetate[1]

ADMA <- keggFind("compound", "ADMA")
KEGG_metabolites_df[str_detect(KEGG_metabolites_df$metabolites_for_KEGG, "ADMA"), "KEGG_ID"] <- names(ADMA[1])
KEGG_metabolites_df[str_detect(KEGG_metabolites_df$metabolites_for_KEGG, "ADMA"), "KEGG_result"] <- ADMA[1]

NMMA <- keggFind("compound", "methylmalonic")
KEGG_metabolites_df[str_detect(KEGG_metabolites_df$metabolites_for_KEGG, "NMMA"), "KEGG_ID"] <- names(NMMA[1])
KEGG_metabolites_df[str_detect(KEGG_metabolites_df$metabolites_for_KEGG, "NMMA"), "KEGG_result"] <- NMMA[1]

methionine_sulfoxide <- keggFind("compound", "methionine S oxide")
KEGG_metabolites_df[str_detect(KEGG_metabolites_df$metabolites_for_KEGG, "methionine.sulfoxide"), "KEGG_ID"] <- names(methionine_sulfoxide[1])
KEGG_metabolites_df[str_detect(KEGG_metabolites_df$metabolites_for_KEGG, "methionine.sulfoxide"), "KEGG_result"] <- methionine_sulfoxide[1]

# lastly, phosphatidylcholine B
KEGG_metabolites_df[str_detect(KEGG_metabolites_df$metabolites_for_KEGG, "phosphatidylcholine B"), "KEGG_ID"] <- KEGG_metabolites_df[str_detect(KEGG_metabolites_df$metabolites_for_KEGG, "^phosphatidylcholine"), "KEGG_ID"][1]
KEGG_metabolites_df[str_detect(KEGG_metabolites_df$metabolites_for_KEGG, "phosphatidylcholine B"), "KEGG_result"] <- KEGG_metabolites_df[str_detect(KEGG_metabolites_df$metabolites_for_KEGG, "^phosphatidylcholine"), "KEGG_result"][1]

# now to retrieve pathways

KEGG_pathways <- lapply(KEGG_metabolites_df$KEGG_ID, function(x){
  print(x)
  tryCatch({keggGet(x)}, error = function(e){return(NULL)})
})

pathways_list <- sapply(KEGG_pathways, function(biglist){
  sapply(biglist, function(littlelist){
    littlelist$PATHWAY
  })
})

names(pathways_list) <- KEGG_metabolites_df$metabolites_for_KEGG

## COMMENTED OUT CODE USED TO IDENTIFY SUITABLE SINGLE PATHWAYS TO COVER AS MANY METABOLITES AS POSSIBLE. NECESSARILY MANUAL AND HEURISTIC.
# sofar <- c(KEGG_metabolites_df$metabolites_for_KEGG[sapply(pathways_list, function(X){
#   any(str_detect(X, "Tryptophan"))
# })],
# 
# KEGG_metabolites_df$metabolites_for_KEGG[sapply(pathways_list, function(X){
#   any(str_detect(X, "Glutathione"))
# })],
# 
# KEGG_metabolites_df$metabolites_for_KEGG[sapply(pathways_list, function(X){
#   any(str_detect(X, "proline metabolism"))
# })],
# 
# KEGG_metabolites_df$metabolites_for_KEGG[sapply(pathways_list, function(X){
#   any(str_detect(X, "Vitamin digestion"))
# })],
# 
# KEGG_metabolites_df$metabolites_for_KEGG[sapply(pathways_list, function(X){
#   any(str_detect(X, "beta-Alanine metabolism"))
# })],
# 
# KEGG_metabolites_df$metabolites_for_KEGG[sapply(pathways_list, function(X){
#   any(str_detect(X, "cofactors"))
# })],
# 
# KEGG_metabolites_df$metabolites_for_KEGG[sapply(pathways_list, function(X){
#   any(str_detect(X, "Glycine, serine and threonine metabolism"))
# })],
# KEGG_metabolites_df$metabolites_for_KEGG[sapply(pathways_list, function(X){
#   any(str_detect(X, "Biosynthesis of amino acids"))
# })],
# 
# KEGG_metabolites_df$metabolites_for_KEGG[sapply(pathways_list, function(X){
#   any(str_detect(X, "Sphingolipid metabolism"))
# })],
# 
# KEGG_metabolites_df$metabolites_for_KEGG[sapply(pathways_list, function(X){
#   any(str_detect(X, "Cholesterol metabolism"))
# })],
# 
#   KEGG_metabolites_df$metabolites_for_KEGG[sapply(pathways_list, function(X){
#   any(str_detect(X, "TCA cycle"))
# })],
# 
# KEGG_metabolites_df$metabolites_for_KEGG[sapply(pathways_list, function(X){
#   any(str_detect(X, "Purine metabolism"))
# })],
# 
# KEGG_metabolites_df$metabolites_for_KEGG[sapply(pathways_list, function(X){
#   any(str_detect(X, "Glycerolipid metabolism"))
# })],
# 
# KEGG_metabolites_df$metabolites_for_KEGG[sapply(pathways_list, function(X){
#   any(str_detect(X, "Pyrimidine metabolism"))
# })],
# 
# KEGG_metabolites_df$metabolites_for_KEGG[sapply(pathways_list, function(X){
#   any(str_detect(X, "Pentose phosphate pathway"))
# })],
# 
# KEGG_metabolites_df$metabolites_for_KEGG[sapply(pathways_list, function(X){
#   any(str_detect(X, "Glycerophospholipid metabolism"))
# })],
# 
# KEGG_metabolites_df$metabolites_for_KEGG[sapply(pathways_list, function(X){
#   any(str_detect(X, "nicotinamide"))
# })],
# 
# KEGG_metabolites_df$metabolites_for_KEGG[sapply(pathways_list, function(X){
#   any(str_detect(X, "Galactose"))
# })]
# 
# )
# 
# unique(sofar)
# 
# KEGG_metabolites_df$metabolites_for_KEGG[sapply(pathways_list, function(X){
#   any(str_detect(X, "Amino acid"))
# })] %in% sofar
# sum(KEGG_metabolites_df$metabolites_for_KEGG == "lysophosphatidylethanolamine")
# 
# 
# tabme <- table(unlist(pathways_list))
# tabme[order(tabme, decreasing = FALSE)]

# now to implement it

pathways_vec <- c("Tryptophan metabolism",
                  "Glutathione metabolism",
                  "Arginine and proline metabolism",
                  "beta-Alanine metabolism",
                  "Glycine, serine and threonine metabolism",
                  "Biosynthesis of amino acids",
                  "Sphingolipid metabolism",
                  "Cholesterol metabolism",
                  "TCA cycle",
                  "Purine metabolism",
                  "Glycerolipid metabolism",
                  "Pyrimidine metabolism",
                  "Pentose phosphate pathway",
                  "Glycerophospholipid metabolism",
                  "Nicotinate and nicotinamide metabolism",
                  "Galactose metabolism")

for(i in 1:length(pathways_vec)){
  
  KEGG_metabolites_df[sapply(pathways_list, function(X){
    any(str_detect(X, pathways_vec[i]))}), "pathway"] <- pathways_vec[i]
}

# and some manual annotation

KEGG_metabolites_df[str_detect(KEGG_metabolites_df$metabolites_for_KEGG, "arnitine"), "pathway"] <- "Carnitine metabolism"
KEGG_metabolites_df[str_detect(KEGG_metabolites_df$metabolites_for_KEGG, "diacylglycerol"), "pathway"] <- "Diacylglycerols"

amino_acids <- c("tryptophan",
                 "phenylalanine",
                 "glycine",
                 "alanine",
                 "valine",
                 "isoleucine",
                 "leucine",
                 "methionine",
                 "proline",
                 "tyrosine",
                 "serine",
                 "threonine",
                 "asparagine",
                 "glutamine",
                 "cysteine",
                 "lysine",
                 "arginine",
                 "histidine",
                 "aspartate",
                 "glutamate",
                 "ornithine",
                 "taurine",
                 "pipecolic acid")

KEGG_metabolites_df[KEGG_metabolites_df$metabolites_for_KEGG %in% amino_acids, "pathway"] <- "Amino acids"

KEGG_metabolites_df[KEGG_metabolites_df$metabolites_for_KEGG %in% c("phosphoglycerate", "F1P F6P G1P G6P", "hexoses  HILIC neg ", "hexoses  HILIC pos "), "pathway"] <- "Glycolysis"

KEGG_metabolites_df[str_detect(KEGG_metabolites_df$metabolites_for_KEGG, "malate"), "pathway"] <- "TCA cycle"
KEGG_metabolites_df[str_detect(KEGG_metabolites_df$metabolites_for_KEGG, "deoxycholate"), "pathway"] <- "Bile acids"
KEGG_metabolites_df[str_detect(KEGG_metabolites_df$metabolites_for_KEGG, "cystathionine"), "pathway"] <- "Glutathione metabolism"
KEGG_metabolites_df[str_detect(KEGG_metabolites_df$metabolites_for_KEGG, "hydroxybutyrate"), "pathway"] <- "Glutathione metabolism"
KEGG_metabolites_df[str_detect(KEGG_metabolites_df$metabolites_for_KEGG, "hydroxyglutarate"), "pathway"] <- "TCA cycle"
KEGG_metabolites_df[str_detect(KEGG_metabolites_df$metabolites_for_KEGG, "acetylglycine"), "pathway"] <- "Glycine, serine and threonine metabolism"
KEGG_metabolites_df[str_detect(KEGG_metabolites_df$metabolites_for_KEGG, "allantoin"), "pathway"] <- "Purine metabolism"
KEGG_metabolites_df[str_detect(KEGG_metabolites_df$metabolites_for_KEGG, "butyrobetaine"), "pathway"] <- "Carnitine metabolism"
KEGG_metabolites_df[str_detect(KEGG_metabolites_df$metabolites_for_KEGG, "kynurenine"), "pathway"] <- "Tryptophan metabolism"
KEGG_metabolites_df[str_detect(KEGG_metabolites_df$metabolites_for_KEGG, "lactate"), "pathway"] <- "Glycolysis"

# consolidate purine and pyrimidine metabolism into one class

KEGG_metabolites_df[(KEGG_metabolites_df$pathway == "Pyrimidine metabolism" | KEGG_metabolites_df$pathway == "Purine metabolism"), "pathway"] <- "Purine/pyrimidine metabolism"

KEGG_metabolites_df[is.na(KEGG_metabolites_df$pathway), "pathway"] <- "undetermined"

write.table(KEGG_metabolites_df, "output/KEGG_metabolites_annotated.txt")
saveRDS(KEGG_metabolites_df, "output/KEGG_metabolites_annotated.rds")

# KEGG_metabolites_df <- readRDS("output/KEGG_metabolites_annotated.rds")

#### RNA-seq CCLE ####

# reads final obtained from DepMap data portal [https://depmap.org/portal/download/all/]
CCLEexpression_reads <- read.table("input/CCLE_RNAseq_reads.csv",
                                   sep = ",",
                                   header = TRUE)

row.names(CCLEexpression_reads) <- CCLEexpression_reads[, "X"]
CCLEexpression_reads <- CCLEexpression_reads[, 2:ncol(CCLEexpression_reads)]

colnames(CCLEexpression_reads) <- str_extract(colnames(CCLEexpression_reads), pattern = "ENSG\\d+")

saveRDS(CCLEexpression_reads, "output/CCLEexpression_reads.rds")
# CCLEexpression_reads <- readRDS("output/CCLEexpression_reads.rds")

CCLEdiseasetypes <-  unique(metadata$primary_disease)

# this gets metadata for all CCLE cell lines
metadata <- depmap_metadata()
useful_metadata <- metadata[(metadata$depmap_id %in% row.names(metabolomics))&(metadata$depmap_id %in% row.names(CCLEexpression_reads)), ]

# we have some primary diseases wth too few samples to work with.
usablecounts <- sapply(CCLEdiseasetypes, function(x){
  
  templines <- metadata[metadata$primary_disease == x, "depmap_id"]
  return(sum(templines$depmap_id %in% row.names(CCLEexpression_reads)))
  
})

# we set an arbitrary cutoff at 20 cell lines
sufficientcounts <- usablecounts[!(usablecounts < 20)]

usableIds <- sapply(names(sufficientcounts), function(x){
  
  templines <- metadata[metadata$primary_disease == x, "depmap_id"]
  return(templines[templines$depmap_id %in% row.names(CCLEexpression_reads), "depmap_id"])
  
})

readcountlist <- list()

for(i in 1:length(usableIds)){
  
  readcountlist[[i]] <- t(CCLEexpression_reads[usableIds[[i]],])
  names(readcountlist)[i] <- str_remove(names(usableIds)[i], ".depmap_id$")
  
}

CCLE_MOR_list <- lapply(readcountlist, function(x){

  intreads <- apply(x, 2, as.integer)

  # DESeq2 wants a colData object. Not actually used for the normalisation. Here we can use the sample IDs with tissue.
  col_tissue <- metadata[match(colnames(intreads), metadata$depmap_id), "primary_disease"]
  col_tissue <- as.matrix(col_tissue)

  rownames(col_tissue) <- colnames(intreads)
  colnames(col_tissue) <- c("primary_disease")

  # need to create a DESeq2 object. Design set to ~1 allows for use of estimateSizeFactors
  tempdds <- DESeqDataSetFromMatrix(countData = intreads, colData = col_tissue, design = ~ 1)

  # this function estimates the scaling factors from the samples from the median of ratios wrt to the geometric mean for each gene across samples
  tempdds <- estimateSizeFactors(tempdds)

  # put the counts normalised by the scaling factors in a new object
  normcounts <- counts(tempdds, normalized = TRUE)

  row.names(normcounts) <- rownames(x)

  return(normcounts)

})

# set the names in the list to be the tissue names
saveRDS(CCLE_MOR_list, "output/CCLE_MOR_list.rds")
# CCLE_MOR_list <- readRDS("output/CCLE_MOR_list.rds")

## CCLE across cancer types

CCLE_readintegers <- apply(CCLEexpression_reads, 1, as.integer)

# DESeq2 wants a colData object. Not actually used for the normalisation. Here we can use the sample IDs with tissue.
CCLEacross_col_tissue <- metadata[match(row.names(CCLEexpression_reads), metadata$depmap_id), "primary_disease"]
CCLEacross_col_tissue <- as.matrix(CCLEacross_col_tissue)

rownames(CCLEacross_col_tissue) <- row.names(CCLEexpression_reads)
colnames(CCLEacross_col_tissue) <- c("primary_disease")

# need to create a DESeq2 object. Design set to ~1 allows for use of estimateSizeFactors
CCLEacross_tempdds <- DESeqDataSetFromMatrix(countData = CCLE_readintegers, colData = CCLEacross_col_tissue, design = ~ 1)

# this function estimates the scaling factors from the samples from the median of ratios wrt to the geometric mean for each gene across samples
CCLEacross_tempdds <- estimateSizeFactors(CCLEacross_tempdds)

# put the counts normalised by the scaling factors in a new object
CCLE_MOR_across_cancers <- counts(CCLEacross_tempdds, normalized = TRUE)

row.names(CCLE_MOR_across_cancers) <- colnames(CCLEexpression_reads)

saveRDS(CCLE_MOR_across_cancers, "output/CCLE_MOR_across_cancers.rds")
# CCLE_MOR_across_cancers <- readRDS("output/CCLE_MOR_across_cancers.rds")

# select chosen samples for 100 iterations of pan-cancer or pan-tissue correlations

GTEX_randomsamples <- define.random.samples.for.iterations(database = "GTEX",
                                     iterations = 100,
                                     save = TRUE)

TCGA_randomsamples <- define.random.samples.for.iterations(database = "TCGA",
                                     iterations = 100,
                                     save = TRUE)

CCLE_randomsamples <- define.random.samples.for.iterations(database = "CCLE",
                                     iterations = 100,
                                     save = TRUE)

#### METABOLOMICS FOR FIGURE 1 ####

#### hierarchical clustering ####

# cluster metabolites to see what generally goes together

metab_dist <- dist(t(metabolomics))
metab_clust <- hclust(metab_dist)

# make full dendrogram for Suppplementary Figure 1 / Fig S1

metab_clust_copy2 <- metab_clust

metab_clust_copy2$labels <- convert.metabolite.labels.for.plot(metab_clust$labels)

metabdatafull <- dendro_data_k(metab_clust_copy2, k = 1)

pdf("graphics/metabolite_clustering_full.pdf", height = 11, width = 8)

plot_ggdendro(hcdata = metabdatafull, scale.color = c("#000000", "#FF0000"),
              branch.size = 0.2,
              expand.y = 0.25, label.size = 1.5) + 
  theme_classic() + 
  theme(axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())

dev.off()

# make sidebar

dend_raw <- as.dendrogram(metab_clust)

pathway_df <- KEGG_metabolites_df[, c("V1", "pathway")]

mycolourpallette <- pals::alphabet(n = length(unique(pathway_df$pathway)))
names(mycolourpallette) <- unique(pathway_df$pathway)

pathway_df[, "mycolours"] <- mycolourpallette[pathway_df$pathway]

labels_colors(dend_raw) <- "grey"
colourlabels <- labels_colors(dend_raw)
labels_colors(dend_raw) <- pathway_df[match(names(colourlabels), pathway_df$V1), "mycolours"]
labels(dend_raw) <- convert.metabolite.labels.for.plot(labels(dend_raw))

forsidebar_df <- data.frame(cbind(labels = colnames(metabolomics)[order.hclust(metab_clust)], mycolours = labels_colors(dend_raw)[convert.metabolite.labels.for.plot(colnames(metabolomics)[order.hclust(metab_clust)])]))
forsidebar_df[, "pseudox"] <- "pseudox"
forsidebar_df[, "order"] <- seq.int(nrow(forsidebar_df))
forsidebar_df[, "rev_order"] <- seq(from = 225, to = 1, length.out = 225)
forsidebar_df[, "pathway"] <- KEGG_metabolites_df[match(forsidebar_df$labels, KEGG_metabolites_df$V1), "pathway"]
forsidebar_df[forsidebar_df$pathway == "undetermined", "mycolours"] <- "white"

pdf("graphics/metabolomics_raw_sidebar.pdf",
    width = 1, 
    height = 9)

ggplot(aes(x= pseudox, y = order), data = forsidebar_df) + 
  geom_tile(fill = forsidebar_df$mycolours) + 
  theme_classic() + 
  theme(axis.text = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank()) + 
  scale_size_manual(values = 0)

dev.off()

#### Principal component analysis for Fig 1C ####

prcomp_metabolomics_raw <- prcomp(t(metabolomics))

# call summary for variance associated with first two PCs for plot
summary(prcomp_metabolomics_raw)

prcomp_raw_for_ggplot <- data.frame(PC1 = prcomp_metabolomics_raw$x[, 1], PC2 = prcomp_metabolomics_raw$x[, 2])
prcomp_raw_for_ggplot[, "pathway"] <- KEGG_metabolites_df[match(row.names(prcomp_raw_for_ggplot), KEGG_metabolites_df$V1), "pathway"]
prcomp_raw_for_ggplot[, "mycolours"] <- mycolourpallette[prcomp_raw_for_ggplot$pathway]

prcomp_NNMTalone <- prcomp_raw_for_ggplot[str_detect(row.names(prcomp_raw_for_ggplot), "methylnicotinamide"), ]

write.table(prcomp_raw_for_ggplot,
            file = "plot_data/Fig 1/Fig_1C_plotdata.txt",
            row.names = TRUE,
            col.names = TRUE,
            sep = "\t")

pdf("graphics/PC_plot.pdf",
    height = 2.5,
    width = 2.5)

ggplot(aes(x = PC1, y = PC2, colour = pathway), data = prcomp_raw_for_ggplot) + 
  geom_point(color = prcomp_raw_for_ggplot$mycolours, size = 0.8) + 
  geom_point(data = prcomp_NNMTalone, colour = "red", size = 5, shape = 1, stroke = 2) +
  theme(legend.position = "none") +
  theme_classic() +
  xlab("PC1 (21.96 % of variance)") + 
  ylab("PC2 (8.25 % of variance)") +
  scale_y_continuous(breaks = seq(-10, 10, 5),
                     expand = c(0, 0),
                     limits = c(-10, 10)) +
  theme(axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10))

dev.off()

#### turn metabolite concentrations into primary disease adusted Z scores ####

metab_for_lm <- cbind(data.frame(metabolomics[match(useful_metadata$depmap_id, row.names(metabolomics)), ]), useful_metadata[, "primary_disease"])

metab_CCLE_z <- data.frame(matrix(nrow = nrow(metab_for_lm), ncol = ncol(metab_for_lm) - 1))
row.names(metab_CCLE_z) <- row.names(metab_for_lm)
colnames(metab_CCLE_z) <- colnames(metab_for_lm)[1:(ncol(metabolomics))]

for (i in 1:length(unique(metab_for_lm$primary_disease))){

  diseasetype <- unique(metab_for_lm$primary_disease)[i]
  disease_metabolites <- metab_for_lm[metab_for_lm$primary_disease == diseasetype, colnames(metab_for_lm) != "primary_disease"]
  mean_metab <- apply(disease_metabolites, 2, mean)
  sd_metab <- apply(disease_metabolites, 2, sd)

  disease_metab_z <- apply(disease_metabolites, 1, function(x){
    (x - mean_metab) / sd_metab
  })

  metab_CCLE_z[match(colnames(disease_metab_z), row.names(metab_CCLE_z)), ] <- t(disease_metab_z)

}

# exclude PRDM16 from individual analysis because it has 0 expression in many samples and produces spurious results
# include NNMT in this analysis so that supplementary table includes NNMT correlations

HMT_metab_corrlist <- lapply(c(SET_HMTs$ensembl_gene_id, NNMT_ensembl), function(gene_id){

  tryCatch({CCLE.correlate.metabolites.with.gene(ensembl_gene = gene_id)},
           error = function(e){})
  
})

names(HMT_metab_corrlist) <- c(SET_HMTs$ensembl_gene_id, NNMT_ensembl)

HMT_metab_corrlist <- HMT_metab_corrlist[!sapply(HMT_metab_corrlist, is.null)]

HMT_metab_corr_df <- data.frame(do.call(cbind, lapply(HMT_metab_corrlist, function(x){x[, "corr"]})))
HMT_metab_FDR_df <- data.frame(do.call(cbind, lapply(HMT_metab_corrlist, function(x){x[, "padjust_vec"]})))
HMT_metab_p_df <- data.frame(do.call(cbind, lapply(HMT_metab_corrlist, function(x){x[, "p_value"]})))

# how many genes are significant at FDR < 0.05? 18 for methylnicotinamide. 
totalHMTssig <- apply(HMT_metab_p_df[, 1:34], 1, function(x){ 
  sum(x < 0.05)
})

totalHMTssig[order(totalHMTssig, decreasing = TRUE)]

HMT_metab_p_df[, "geometric_mean_across_HMTs"] <- apply(HMT_metab_p_df[, 1:(ncol(HMT_metab_p_df) - 1)], 1, function(x){
  
  gm_mean(x)
  
})

# find range across HMTs (cited in manuscript text)

min(HMT_metab_corr_df["X1.methylnicotinamide", !(colnames(HMT_metab_p_df) %in% c("geometric_mean_across_HMTs", NNMT_ensembl))])
max(HMT_metab_corr_df["X1.methylnicotinamide", !(colnames(HMT_metab_p_df) %in% c("geometric_mean_across_HMTs", NNMT_ensembl))])
mean(unlist(HMT_metab_corr_df["X1.methylnicotinamide", !(colnames(HMT_metab_p_df) %in% c("geometric_mean_across_HMTs", NNMT_ensembl))]))

#### now with pooling the HMT read counts ####

HMTs_metab_together <- data.frame(whole.geneset.metabolite.correlation(SET_HMTs$ensembl_gene_id))
HMTs_metab_together[, "FDR"] <- p.adjust(HMTs_metab_together$genesetpvalues, method = "BH")

HMTs_metab_together[, "label"] <- convert.metabolite.labels.for.plot(row.names(HMTs_metab_together))

colours_for_plot <- map2color(-log10(HMTs_metab_together$genesetpvalues), pal = colorRampPalette(c("grey", "grey", "red", "red", "red"))(100))

write.table(HMTs_metab_together,
            file = "plot_data/Fig 1/Fig_1A_plotdata.txt",
            row.names = TRUE,
            col.names = TRUE,
            sep = "\t")

#### figure 1A (metabolite correlation bubble plot) ####

pdf("~/NNMT_manuscript/graphics/allHMT_metabcorr.pdf",
    height = 2.5,
    width = 2.5)

ggplot(aes(x = genesetresults, y = -log10(FDR), size = -log10(FDR), label = label), data = HMTs_metab_together) + 
  geom_point(col = colours_for_plot, alpha = 0.7) + 
  theme_classic() + 
  geom_label_repel(size = 2, max.overlaps = 40) +
  geom_vline(xintercept = 0, linetype = "dashed", col = "grey") + 
  theme(legend.position = "none",
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10)) + 
  ylab(substitute("log"[10]~"(false discovery rate)")) + 
  xlab(substitute("Metabolite-HMTs Pearson's"~italic(r))) + 
  coord_cartesian(xlim = c(-0.3, 0.3),
                  ylim = c(0, 16))

dev.off()

#### Supplementary Table S2 (metabolomics correlations excel table) ####

output_table_corr <- data.frame(matrix(nrow = nrow(HMT_metab_corr_df), ncol = ncol(HMT_metab_corr_df) + 3))
output_table_corr[ , 1] <- convert.metabolite.labels.for.plot(row.names(HMT_metab_corr_df))
output_table_corr[, 2:(ncol(output_table_corr)-2)] <- HMT_metab_corr_df
colnames(output_table_corr) <- c("Metabolite", SET_HMTs[match(colnames(HMT_metab_corr_df)[1:(ncol(HMT_metab_corr_df) - 1)], SET_HMTs$ensembl_gene_id), "hgnc_symbol"], "NNMT", "Average correlation across HMTs", "Correlation to all HMTs pooled")
output_table_corr[, "Average correlation across HMTs"] <- rowMeans(output_table_corr[, 2:(ncol(output_table_corr) - 3)])
output_table_corr[, "Correlation to all HMTs pooled"] <- HMTs_metab_together[match(output_table_corr$Metabolite, HMTs_metab_together$label), "genesetresults"]

# sort alphabetically for easy perusing by eye
output_table_corr[ , 2:(ncol(output_table_corr)-2)] <- output_table_corr[, match(sort(colnames(output_table_corr)[2:(ncol(output_table_corr)-2)]), colnames(output_table_corr))]
colnames(output_table_corr)[2:(ncol(output_table_corr)-2)] <- sort(colnames(output_table_corr)[2:(ncol(output_table_corr)-2)])

# put NNMT last to avoid confusion
output_table_corr <- output_table_corr %>% relocate(NNMT, .after = last_col())

output_table_FDR <- data.frame(matrix(nrow = nrow(HMT_metab_FDR_df), ncol = ncol(HMT_metab_FDR_df) + 3))
output_table_FDR[ , 1] <- convert.metabolite.labels.for.plot(row.names(HMT_metab_FDR_df))
output_table_FDR[, 2:(ncol(output_table_FDR)-2)] <- HMT_metab_FDR_df
colnames(output_table_FDR) <- c("Metabolite", SET_HMTs[match(colnames(HMT_metab_FDR_df)[!str_detect(colnames(HMT_metab_FDR_df), pattern = "mean") & !str_detect(colnames(HMT_metab_FDR_df), NNMT_ensembl)], SET_HMTs$ensembl_gene_id), "hgnc_symbol"], "NNMT", "Geometric mean FDR across HMTs", "FDR correlation to all HMTs pooled")
output_table_FDR[, "Geometric mean FDR across HMTs"] <- apply(output_table_FDR[, 2:(ncol(output_table_FDR) - 3)], 1, gm_mean)
output_table_FDR[, "FDR correlation to all HMTs pooled"] <- HMTs_metab_together[match(output_table_FDR$Metabolite, HMTs_metab_together$label), "FDR"]

# sort alphabetically for easy perusing by eye
output_table_FDR[ , 2:(ncol(output_table_FDR)-2)] <- output_table_FDR[, match(sort(colnames(output_table_FDR)[2:(ncol(output_table_FDR)-2)]), colnames(output_table_FDR))]
colnames(output_table_FDR)[2:(ncol(output_table_FDR)-2)] <- sort(colnames(output_table_FDR)[2:(ncol(output_table_FDR)-2)])

# put NNMT last to avoid confusion
output_table_FDR <- output_table_FDR %>% relocate(NNMT, .after = last_col())

output_table_pval <- data.frame(matrix(nrow = nrow(HMT_metab_p_df), ncol = ncol(HMT_metab_p_df) + 3))
output_table_pval[ , 1] <- convert.metabolite.labels.for.plot(row.names(HMT_metab_p_df))
output_table_pval[, 2:(ncol(output_table_pval)-2)] <- HMT_metab_p_df
colnames(output_table_pval) <- c("Metabolite", SET_HMTs[match(colnames(HMT_metab_p_df)[!str_detect(colnames(HMT_metab_p_df), pattern = "mean") & !str_detect(colnames(HMT_metab_p_df), NNMT_ensembl)], SET_HMTs$ensembl_gene_id), "hgnc_symbol"], "NNMT", "Geometric mean p-value across HMTs", "p-value correlation to all HMTs pooled")
output_table_pval[, "Geometric mean p-value across HMTs"] <- apply(output_table_FDR[, 2:(ncol(output_table_FDR) - 3)], 1, gm_mean)
output_table_pval[, "p-value correlation to all HMTs pooled"] <- HMTs_metab_together[match(output_table_pval$Metabolite, HMTs_metab_together$label), "FDR"]

# sort alphabetically for easy perusing by eye
output_table_pval[ , 2:(ncol(output_table_pval)-2)] <- output_table_pval[, match(sort(colnames(output_table_pval)[2:(ncol(output_table_pval)-2)]), colnames(output_table_pval))]
colnames(output_table_pval)[2:(ncol(output_table_pval)-2)] <- sort(colnames(output_table_pval)[2:(ncol(output_table_pval)-2)])

# put NNMT last to avoid confusion
output_table_pval <- output_table_pval %>% relocate(NNMT, .after = last_col())

write.xlsx(list("Pearson's correlations" = output_table_corr,
                "raw p-values" = output_table_pval,
                "FDR-adjusted p-values" = output_table_FDR),
           file = "output/CCLE_HMT_metabolite_correlations.xlsx",
           rowNames = FALSE, colnames = TRUE)

#### correlating NNMT expression with MNA in CCLE for Fig S2A / Supplementary Figure 2A ####

CCLE_NNMT_list <- lapply(CCLE_MOR_list, function(x){x[NNMT_ensembl, ]})

CCLE_NNMT_z_list <-  lapply(CCLE_NNMT_list, function(x){
  
  mean <- mean(log10(x + 1), na.rm = TRUE)
  sd <- sd(log10(x + 1), na.rm = TRUE)
  
  (log10(x) - mean) / sd
  
})       

CCLE_NNMT_z_vec <- do.call(c, CCLE_NNMT_z_list)

names(CCLE_NNMT_z_vec) <- str_extract(names(CCLE_NNMT_z_vec), pattern = "ACH.*")

# remove infinite values
CCLE_NNMT_z_vec <- CCLE_NNMT_z_vec[is.finite(CCLE_NNMT_z_vec)]

CCLE_NNMT_z_vec <- CCLE_NNMT_z_vec[match(row.names(metab_CCLE_z), names(CCLE_NNMT_z_vec))]
CCLE_NNMT_z_vec <- CCLE_NNMT_z_vec[!is.na(CCLE_NNMT_z_vec)]

CCLE_MNA_z_vec <- metab_CCLE_z[match(names(CCLE_NNMT_z_vec), row.names(metab_CCLE_z)), "X1.methylnicotinamide"]
names(CCLE_MNA_z_vec) <- names(CCLE_NNMT_z_vec)

CCLE_NNMT_MNA_forplot <- data.frame(cbind(CCLE_NNMT_z_vec, CCLE_MNA_z_vec))

CCLE_NNMT_MNA_cortest <- cor.test(CCLE_NNMT_z_vec, CCLE_MNA_z_vec)       
CCLE_NNMT_MNA_cortest$p.value

CCLE_NNMT_MNA_linefit <- lm(CCLE_MNA_z_vec ~ CCLE_NNMT_z_vec)
sum_CCLE_NNMT_MNA_linefit <- summary(CCLE_NNMT_MNA_linefit)

write.table(CCLE_NNMT_MNA_forplot,
            file = "plot_data/Fig S2/Fig_S2A_plotdata.txt",
            row.names = TRUE,
            col.names = TRUE,
            sep = "\t")

pdf("~/NNMT_manuscript/graphics/CCLE_NNMT_vs_MNA_scatterplot.pdf", 
    height = 2.5,
    width = 2.5)

ggplot(aes(x = CCLE_NNMT_z_vec, 
           y = CCLE_MNA_z_vec), 
       data = CCLE_NNMT_MNA_forplot) + 
  geom_point(size = 0.4, 
             alpha = 0.3) + 
  theme_classic() + 
  ylab("Cell line 1MNA level (z-score)") + 
  xlab(substitute(italic(NNMT)~"expression (RNA-seq z)")) + 
  theme(axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10)) + 
  scale_y_continuous(breaks = seq(-4, 4, 2),
                     limits = c(-4, 4),
                     expand = c(0,0)) + 
  scale_x_continuous(breaks = seq(-4, 4, 2),
                     limits = c(-4, 4),
                     expand = c(0,0)) + 
  annotate(geom = "text", 
           x = -1.2, 
           y = 3.7, 
           label = substitute(italic(r)~"= 0.729")) +
  annotate(geom = "text", 
           x = -1.2, 
           y = 2.9, 
           label = substitute(italic(p)~"= 3.19 x"~10^-148),
           size = 4) +
  # coord_cartesian(xlim = c(-4.1, 4.1), 
  #                 ylim = c(-4.1, 4.1)) +
  geom_smooth(method = "lm", 
              se = FALSE,
              colour = "red")

dev.off()

#### correlating 1MNA to total HMT in CCLE ####

CCLE_allHMT_list <- lapply(CCLE_MOR_list, function(x){colSums(x[SET_HMTs$ensembl_gene_id, ])})

CCLE_allHMT_z_list <-  lapply(CCLE_allHMT_list, function(x){
  
  mean <- mean(log10(x + 1), na.rm = TRUE)
  sd <- sd(log10(x + 1), na.rm = TRUE)
  
  (log10(x + 1) - mean) / sd
  
})       

CCLE_allHMT_z_vec <- do.call(c, CCLE_allHMT_z_list)

names(CCLE_allHMT_z_vec) <- str_extract(names(CCLE_allHMT_z_vec), pattern = "ACH.*")

# remove infinite values

CCLE_MNA_z_vec_forHMTs <- metab_CCLE_z[row.names(metab_CCLE_z) %in% names(CCLE_allHMT_z_vec), "X1.methylnicotinamide"]
names(CCLE_MNA_z_vec_forHMTs) <- row.names(metab_CCLE_z)[row.names(metab_CCLE_z) %in% names(CCLE_allHMT_z_vec)]

MNA_allHMTs_matching_samples <- intersect(names(CCLE_MNA_z_vec_forHMTs), names(CCLE_allHMT_z_vec))

CCLE_allHMT_z_vec <- CCLE_allHMT_z_vec[MNA_allHMTs_matching_samples]
CCLE_MNA_z_vec_forHMTs <- CCLE_MNA_z_vec_forHMTs[MNA_allHMTs_matching_samples]

CCLE_allHMT_MNA_forplot <- data.frame(cbind(CCLE_allHMT_z_vec, CCLE_MNA_z_vec_forHMTs))

CCLE_allHMT_MNA_cortest <- cor.test(CCLE_allHMT_z_vec, CCLE_MNA_z_vec_forHMTs)       
CCLE_allHMT_MNA_cortest$p.value

CCLE_allHMT_MNA_linefit <- lm(CCLE_MNA_z_vec_forHMTs ~ CCLE_allHMT_z_vec)
sum_CCLE_allHMT_MNA_linefit <- summary(CCLE_allHMT_MNA_linefit)

pdf("~/NNMT_manuscript/graphics/CCLE_allHMT_vs_MNA_scatterplot.pdf", 
    height = 2.5,
    width = 2.5)

ggplot(aes(x = CCLE_allHMT_z_vec, 
           y = CCLE_MNA_z_vec_forHMTs), 
       data = CCLE_allHMT_MNA_forplot) + 
  geom_point(size = 0.4, 
             alpha = 0.3) + 
  theme_classic() + 
  ylab("Cell line 1MNA level (z-score)") + 
  xlab("total HMT expression\n(RNA-seq z-score)") + 
  theme(axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10)) + 
  annotate(geom = "text", 
           x = -2.2, 
           y = 3.7, 
           label = substitute(italic(r)~"= 0.284"),
           size = 3) +
  annotate(geom = "text", 
           x = -2.2, 
           y = 2.9, 
           label = substitute(italic(p)~"= 5.50 x"~10^-18),
           size = 3) +
  annotate(geom = "text",
           x = 3,
           y = -3,
           label = "CCLE",
           colour = "dodgerblue4",
           size = 6,
           fontface = 2) +
  coord_cartesian(xlim = c(-4.1, 4.1), 
                  ylim = c(-4.1, 4.1)) +
  geom_smooth(method = "lm", 
              se = FALSE,
              colour = "red")

dev.off()

#### correlating NNMT expression with protein in NCI60 for Fig S2D / Supplementary Figure 2D ####

NCI60_RNAseq <- readRDS("output/NCI60_RNAseq_normalised_counts.rds")

# NNMT proteomics values downloaded from Cell Miner Database [https://discover.nci.nih.gov/rsconnect/cellminercdb/]
NNMT_proteomics <- read.table("input/NNMT_proteomics.txt",
                              quote = "",
                              sep = "\t",
                              header = TRUE)

names_changed_to_match_format <- str_remove(str_remove(str_replace(str_replace(str_remove(NNMT_proteomics$Cell.Line, "^[A-Z]{2,3}:"), "/", "_"), " ", "_"), "\\(TB\\)"), "_ATCC")
NNMT_proteomics$Cell.Line <- names_changed_to_match_format

NNMTproteomics_vec <- NNMT_proteomics$swaNNMT_nci60
NNMT_RNAseq_vec <- NCI60_RNAseq["4837", match(NNMT_proteomics$Cell.Line, colnames(NCI60_RNAseq))]

# NNMT_prot_RNAseq_df <- data.frame(cbind(NNMTproteomics_vec, NNMT_RNAseq_vec))
# NNMT_prot_RNAseq_df[, "tissues"] <- NNMT_proteomics[match(row.names(NNMT_prot_RNAseq_df), NNMT_proteomics$Cell.Line), "tissues"]
# 
# NNMT_prot_RNAseq_df[, "protresid"] <- lm(NNMT_prot_RNAseq_df$NNMTproteomics_vec ~ NNMT_prot_RNAseq_df$tissues)$residuals
# NNMT_prot_RNAseq_df[, "RNAresid"] <- lm(log10(NNMT_prot_RNAseq_df$NNMT_RNAseq_vec + 1) ~ NNMT_prot_RNAseq_df$tissues)$residuals
# 
# plot(NNMT_prot_RNAseq_df$protresid, log10(NNMT_prot_RNAseq_df$NNMT_RNAseq_vec + 1))
# cor.test(NNMT_prot_RNAseq_df$RNAresid, NNMT_prot_RNAseq_df$protresid)

cor.test(NNMT_proteomics$swaNNMT_nci60, NNMT_proteomics$expNNMT_nci60)

write.table(NNMT_proteomics,
            file = "plot_data/Fig S2/Fig_S2D_plotdata.txt",
            row.names = FALSE,
            col.names = TRUE,
            sep = "\t")

pdf("graphics/NCI60_NNMTmRNA_vs_protein_scatterplot.pdf",
    height = 2.5,
    width = 2.5)

ggplot(aes(x = expNNMT_nci60, y = swaNNMT_nci60), data = NNMT_proteomics) + 
  geom_point(size = 0.4, 
             alpha = 0.3) +
  theme_classic() + 
  theme(legend.position = "none",
        axis.text.x = element_text(colour = "black",
                                   size = 10),
        axis.text.y = element_text(colour = "black",
                                   size = 10),
        axis.title.x = element_text(size = 8.5),
        axis.title.y = element_text(size = 8.5)) + 
  ylab("NNMT protein (SWATH-MS signal)") + 
  xlab(substitute(italic(NNMT)~"mRNA level (microarray z)")) + 
  scale_y_continuous(labels = label_number(drop0trailing = TRUE)) +
  annotate(geom = "text",
           x = -0.4,
           y = 5.5,
           label = substitute(italic(r)~"= 0.635")) +
  annotate(geom = "text",
           x = -0.3,
           y = 5.3,
           label = substitute(italic(p)~"= 6.46 x"~10^-8)) +
  geom_smooth(method = "lm",
              se = FALSE,
              colour = "red")

dev.off()

#### CCLE proteomics ####

# manually find missing HMTs since I can see some of them are there. they are:
# KMT5A  aka SETD8
# KMT5B  aka SUV420H1
# KMT5C   aka SUV420H2
# NSD2    aka WHSC1
# NSD3   aka WHSC1L1
# PRDM13  not found 
# PRDM7   not found
# PRDM9   not found
# SETD4   not found

HMTsnames_for_proteomics <- SET_HMTs[, c("hgnc_symbol", "ensembl_gene_id") ]
HMTsnames_for_proteomics[, "proteomics_name"] <- HMTsnames_for_proteomics$hgnc_symbol
pattern_Vec <- SET_HMTs_annotated$hgnc_symbol[!(SET_HMTs_annotated$hgnc_symbol %in% CCLE_proteomics$Gene_Symbol)][1:5]
replace_Vec <- c("SETD8", 
                 "SUV420H1", 
                 "SUV420H2",
                 "WHSC1",
                 "WHSC1L1")

for(i in 1:length(pattern_Vec)){
  
  HMTsnames_for_proteomics$proteomics_name <- str_replace(HMTsnames_for_proteomics$proteomics_name,
                                                          pattern = pattern_Vec[i],
                                                          replacement = replace_Vec[i])
  
}

# how many entries for each gene? mostly 1. exceptions are WHSC1 / NSD2 and NSD1
sapply(HMTsnames_for_proteomics$proteomics_name, function(x){sum(str_detect(CCLE_proteomics$Gene_Symbol, pattern = x))})

CCLE_proteomics_of_interest <- CCLE_proteomics[match(c(HMTsnames_for_proteomics$proteomics_name, "NNMT"), CCLE_proteomics$Gene_Symbol), 49:426]
row.names(CCLE_proteomics_of_interest) <- c(HMTsnames_for_proteomics$ensembl_gene_id, NNMT_ensembl)

CCLE_proteomics_of_interest <- CCLE_proteomics_of_interest[, str_remove(colnames(CCLE_proteomics_of_interest), pattern = "_Ten.*") %in% useful_metadata$cell_line]
colnames(CCLE_proteomics_of_interest) <- unlist(useful_metadata[match(str_remove(colnames(CCLE_proteomics_of_interest), pattern = "_Ten.*"), useful_metadata$cell_line), "depmap_id"])

CCLE_proteomics_of_interest <- CCLE_proteomics_of_interest[apply(CCLE_proteomics_of_interest,1, function(x){any(!is.na(x))}), ]

# see how HMTs ad NNMT protein levels correlate with RNA levels (uncorrected)
expression_protein_correlations <- sapply(row.names(CCLE_proteomics_of_interest), function(thisgene){
  
  cor.test(unlist(CCLE_proteomics_of_interest[thisgene, ]), unlist(CCLE_MOR_across_cancers[thisgene, match(colnames(CCLE_proteomics_of_interest), colnames(CCLE_MOR_across_cancers))]), method = "spearman", na.rm = T)$estimate
  
})  

expression_protein_correlations[order(expression_protein_correlations, decreasing = T)]

individualcorrs <- apply(CCLE_proteomics_of_interest, 1, function(x){
  
  cor.test(x + rnorm(ncol(CCLE_proteomics_of_interest), 0, 0.01), unlist(CCLE_proteomics_of_interest[NNMT_ensembl, ]) + rnorm(ncol(CCLE_proteomics_of_interest), 0, 0.01))$estimate
  
})

names(individualcorrs) <- SET_HMTs[match(names(individualcorrs), SET_HMTs$ensembl_gene_id), "hgnc_symbol"]

individualcorrs[order(individualcorrs)]


# converting proteomics data to cancer type Z score

proteomics_for_lm <- cbind(data.frame(useful_metadata[match(colnames(CCLE_proteomics_of_interest), useful_metadata$depmap_id), "primary_disease"]), t(CCLE_proteomics_of_interest))

# restrict to cancer types with at least 12 samples for calculating Z score
proteomics_for_lm_cut1 <- proteomics_for_lm[proteomics_for_lm$primary_disease %in% names(table(proteomics_for_lm$primary_disease)[table(proteomics_for_lm$primary_disease) > 12]), ]

# restrict to genes with missing data in < 30 samples (~10%)
proteomics_for_lm_cut2 <- proteomics_for_lm_cut1[, !(colnames(proteomics_for_lm_cut1) %in% names(apply(proteomics_for_lm_cut1, 2, function(x){sum(is.na(x))})[apply(proteomics_for_lm_cut1, 2, function(x){sum(is.na(x))}) > 30]))]

# leaves 298 samples for 20 HMTs + NNMT

proteo_CCLE_z <- data.frame(matrix(nrow = nrow(proteomics_for_lm_cut2), ncol = ncol(proteomics_for_lm_cut2) - 1))
row.names(proteo_CCLE_z) <- row.names(proteomics_for_lm_cut2)
colnames(proteo_CCLE_z) <- colnames(proteomics_for_lm_cut2)[2:(ncol(proteomics_for_lm_cut2))]

for (i in 1:length(unique(proteomics_for_lm_cut2$primary_disease))){
  
  diseasetype <- unique(proteomics_for_lm_cut2$primary_disease)[i]
  disease_proteoolites <- proteomics_for_lm_cut2[proteomics_for_lm_cut2$primary_disease == diseasetype, colnames(proteomics_for_lm_cut2) != "primary_disease"]
  mean_proteo <- apply(disease_proteoolites, 2, function(x){mean(x, na.rm = TRUE)})
  sd_proteo <- apply(disease_proteoolites, 2, function(x){sd(x, na.rm = TRUE)})
  
  disease_proteo_z <- apply(disease_proteoolites, 1, function(x){
    (x - mean_proteo) / sd_proteo
  })
  
  proteo_CCLE_z[match(colnames(disease_proteo_z), row.names(proteo_CCLE_z)), ] <- t(disease_proteo_z)
  
}

row.names(proteo_CCLE_z) <- str_replace(row.names(proteo_CCLE_z), pattern = "\\.", replacement = "-")

Zscoreindividualcorrs <- apply(proteo_CCLE_z, 2, function(x){
  
  cor.test(x + rnorm(nrow(proteo_CCLE_z), 0, 0.01), unlist(proteo_CCLE_z[, NNMT_ensembl]) + rnorm(nrow(proteo_CCLE_z), 0, 0.01))$estimate
  
})

names(Zscoreindividualcorrs) <- SET_HMTs[match(names(Zscoreindividualcorrs), SET_HMTs$ensembl_gene_id), "hgnc_symbol"]

Zscoreindividualcorrs[order(Zscoreindividualcorrs)]

# do proteomics vs RNA-seq correlations for HMTs and NNMT with Z scores instead of absolute

protein_with_mRNA_CCLE <- lapply(colnames(proteo_CCLE_z), function(x){
  CCLE.correlate.mRNA.with.protein(ensembl_gene = x)
})

names(protein_with_mRNA_CCLE) <- colnames(proteo_CCLE_z)

# from the plot we can see that the result is very similar regardless of whether we correct to Z score or not
plot(sapply(protein_with_mRNA_CCLE, function(x){x[1]}), expression_protein_correlations[str_remove(names(expression_protein_correlations), "\\.rho") %in% names(protein_with_mRNA_CCLE)])

#### correlate HMT NNMT protein levels with metabolites for Supplementary Table S3 / Table S3 ####

protein_with_metabolites <- lapply(colnames(proteo_CCLE_z), function(x){
  CCLE.correlate.metabolites.with.protein(ensembl_gene = x)
})

names(protein_with_metabolites) <- colnames(proteo_CCLE_z)

CCLEprotein_metab_corr_df <- data.frame(do.call(cbind, lapply(protein_with_metabolites, function(x){x[, "corr"]})))
CCLEprotein_metab_p_df <- data.frame(do.call(cbind, lapply(protein_with_metabolites, function(x){x[, "padjust_vec"]})))
CCLEprotein_metab_p_unadjusted_df <- data.frame(do.call(cbind, lapply(protein_with_metabolites, function(x){x[, "p_value"]})))

# very few significant correlations for proteins and metabolites.
totalHMTssigprotmet <- apply(CCLEprotein_metab_p_df[, 1:20], 1, function(x){ 
  sum(x < 0.05)
})

min(CCLEprotein_metab_corr_df["X1.methylnicotinamide", 1:20])
max(CCLEprotein_metab_corr_df["X1.methylnicotinamide", 1:20])
mean(unlist(CCLEprotein_metab_corr_df["X1.methylnicotinamide", 1:20]))

avg_CCLE_metab_protein_corrs <- apply(CCLEprotein_metab_corr_df[, 1:20], 1, mean)

# let's try metabolites to overall sample HMT

mean_HMT_protZ <- rowMeans(proteo_CCLE_z[, 1:20], na.rm = TRUE)

overallsampleHMTprot_metab <-  apply(metab_CCLE_z, 2, function(x){
  
  cor.test(x[match(names(mean_HMT_protZ), row.names(metab_CCLE_z))], mean_HMT_protZ, na.rm = TRUE)$estimate
  
})

overallsampleHMTprot_metab_p <-  apply(metab_CCLE_z, 2, function(x){
  
  cor.test(x[match(names(mean_HMT_protZ), row.names(metab_CCLE_z))], mean_HMT_protZ, na.rm = TRUE)$p.value
  
})

overallsampleHMTprot_metab[order(overallsampleHMTprot_metab, decreasing = F)]

# output to Excel table
output_table_corr_prot <- data.frame(matrix(nrow = nrow(CCLEprotein_metab_corr_df), ncol = ncol(CCLEprotein_metab_corr_df) + 3))
output_table_corr_prot[ , 1] <- convert.metabolite.labels.for.plot(row.names(CCLEprotein_metab_corr_df))
output_table_corr_prot[, 2:(ncol(output_table_corr_prot)-2)] <- CCLEprotein_metab_corr_df
colnames(output_table_corr_prot) <- c("Metabolite", SET_HMTs[match(colnames(CCLEprotein_metab_corr_df[1:20]), SET_HMTs$ensembl_gene_id), "hgnc_symbol"], "NNMT", "Average correlation across HMTs", "Correlation to all HMTs pooled")
output_table_corr_prot[, "Average correlation across HMTs"] <- rowMeans(output_table_corr_prot[, 2:(ncol(output_table_corr_prot) - 3)])
output_table_corr_prot[, "Correlation to all HMTs pooled"] <- overallsampleHMTprot_metab[match(output_table_corr_prot$Metabolite, convert.metabolite.labels.for.plot(names(overallsampleHMTprot_metab)))]

# sort alphabetically for easy perusing by eye
output_table_corr_prot[ , 2:(ncol(output_table_corr_prot)-2)] <- output_table_corr_prot[, match(sort(colnames(output_table_corr_prot)[2:(ncol(output_table_corr_prot)-2)]), colnames(output_table_corr_prot))]
colnames(output_table_corr_prot)[2:(ncol(output_table_corr_prot)-2)] <- sort(colnames(output_table_corr_prot)[2:(ncol(output_table_corr_prot)-2)])

# but stick NNMT last to avoid confusion
output_table_corr_prot <- output_table_corr_prot %>% relocate(NNMT, .after = last_col())

output_table_pval_prot <- data.frame(matrix(nrow = nrow(CCLEprotein_metab_p_unadjusted_df), ncol = ncol(CCLEprotein_metab_p_unadjusted_df) + 3))
output_table_pval_prot[ , 1] <- convert.metabolite.labels.for.plot(row.names(CCLEprotein_metab_p_unadjusted_df))
output_table_pval_prot[, 2:(ncol(output_table_pval_prot)-2)] <- CCLEprotein_metab_p_unadjusted_df
colnames(output_table_pval_prot) <- c("Metabolite", SET_HMTs[match(colnames(CCLEprotein_metab_p_unadjusted_df)[!str_detect(colnames(CCLEprotein_metab_p_df), pattern = "mean") & !str_detect(colnames(CCLEprotein_metab_p_unadjusted_df), NNMT_ensembl)], SET_HMTs$ensembl_gene_id), "hgnc_symbol"], "NNMT", "Geometric mean p-value across HMTs", "p-value correlation to all HMTs mean Z-score")
output_table_pval_prot[, "Geometric mean p-value across HMTs"] <- apply(output_table_pval_prot[, 2:(ncol(output_table_pval_prot)-3)], 1, gm_mean)
output_table_pval_prot[, "p-value correlation to all HMTs mean Z-score"] <- overallsampleHMTprot_metab_p[match(output_table_pval_prot$Metabolite, convert.metabolite.labels.for.plot(names(overallsampleHMTprot_metab_p)))]

# sort alphabetically for easy perusing by eye
output_table_pval_prot[ , 2:(ncol(output_table_pval_prot)-2)] <- output_table_pval_prot[, match(sort(colnames(output_table_pval_prot)[2:(ncol(output_table_pval_prot)-2)]), colnames(output_table_pval_prot))]
colnames(output_table_pval_prot)[2:(ncol(output_table_pval_prot)-2)] <- sort(colnames(output_table_pval_prot)[2:(ncol(output_table_pval_prot)-2)])

# but stick NNMT last to avoid confusion
output_table_pval_prot <- output_table_pval_prot %>% relocate(NNMT, .after = last_col())

output_table_FDR_prot <- data.frame(matrix(nrow = nrow(CCLEprotein_metab_p_df), ncol = ncol(CCLEprotein_metab_p_df) + 3))
output_table_FDR_prot[ , 1] <- convert.metabolite.labels.for.plot(row.names(CCLEprotein_metab_p_df))
output_table_FDR_prot[, 2:(ncol(output_table_FDR_prot)-2)] <- CCLEprotein_metab_p_df
colnames(output_table_FDR_prot) <- c("Metabolite", SET_HMTs[match(colnames(CCLEprotein_metab_p_df)[!str_detect(colnames(CCLEprotein_metab_p_df), pattern = "mean") & !str_detect(colnames(CCLEprotein_metab_p_df), NNMT_ensembl)], SET_HMTs$ensembl_gene_id), "hgnc_symbol"], "NNMT", "Geometric mean FDR across HMTs", "FDR correlation to all HMTs mean Z-score")
output_table_FDR_prot[, "Geometric mean FDR across HMTs"] <- apply(output_table_FDR_prot[, 2:(ncol(output_table_FDR_prot)-3)], 1, gm_mean)
output_table_FDR_prot[, "FDR correlation to all HMTs mean Z-score"] <- p.adjust(overallsampleHMTprot_metab_p[match(output_table_FDR_prot$Metabolite, convert.metabolite.labels.for.plot(names(overallsampleHMTprot_metab_p)))], method = "BH")

# sort alphabetically for easy perusing by eye
output_table_FDR_prot[ , 2:(ncol(output_table_FDR_prot)-2)] <- output_table_FDR_prot[, match(sort(colnames(output_table_FDR_prot)[2:(ncol(output_table_FDR_prot)-2)]), colnames(output_table_FDR_prot))]
colnames(output_table_FDR_prot)[2:(ncol(output_table_FDR_prot)-2)] <- sort(colnames(output_table_FDR_prot)[2:(ncol(output_table_FDR_prot)-2)])

# but stick NNMT last to avoid confusion
output_table_FDR_prot <- output_table_FDR_prot %>% relocate(NNMT, .after = last_col())

write.xlsx(list("Pearson's correlations" = output_table_corr_prot,
                "raw p-values" = output_table_pval_prot,
                "FDR-adjusted p-values" = output_table_FDR_prot),
           file = "output/CCLE_HMTproteomics_metabolite_correlations.xlsx",
           rowNames = FALSE, colnames = TRUE)

#### NNMT protein to HMT protein for Supplementary Figure 2E / Fig S2E  ####

CCLE_NNMTtoHMTprotein_forplot <- data.frame(HMT = mean_HMT_protZ, NNMT = proteo_CCLE_z[names(mean_HMT_protZ), NNMT_ensembl])

cor.test(CCLE_NNMTtoHMTprotein_forplot$HMT, CCLE_NNMTtoHMTprotein_forplot$NNMT)

write.table(CCLE_NNMTtoHMTprotein_forplot,
            file = "plot_data/Fig S2/Fig_S2E_plotdata.txt",
            row.names = TRUE,
            col.names = TRUE,
            sep = "\t")

pdf("graphics/CCLE_NNMTprotein_vs_HMTprotein_scatterplot.pdf",
    width = 2.5,
    height = 2.5)

ggplot(aes(y = NNMT, x = HMT), data = CCLE_NNMTtoHMTprotein_forplot) +
  geom_point(size = 0.4,
             alpha = 0.3) +
  theme_classic() + 
  theme(legend.position = "none",
        axis.text.x = element_text(colour = "black",
                                   size = 10),
        axis.text.y = element_text(colour = "black",
                                   size = 10),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        plot.margin = margin(t = 5, 
                             l = 5,
                             b = 5,
                             r = 8)) + 
  ylab("NNMT protein level (z-score)") + 
  xlab("Mean HMT protein (z-score)") + 
  scale_x_continuous(breaks = seq(-1.5, 1.5, 0.75),
                     limits = c(-1.5, 1.5),
                     expand = c(0,0),
                     labels = label_number(drop0trailing = TRUE)) + 
  scale_y_continuous(breaks = seq(-2, 3, 1),
                     limits = c(-2, 3),
                     expand = c(0,0)) +
  annotate(geom = "text",
           x = -0.6,
           y = 2.7,
           label = substitute(italic(r)~"= -0.238")) +
  annotate(geom = "text",
           x = -0.6,
           y = 2.3,
           label = substitute(italic(p)~"= 3.32 x"~10^-5)) +
  geom_smooth(method = "lm",
              se = FALSE,
              colour = "red")

dev.off()


CCLE_NNMTtoHMTprotein_raw <- data.frame(HMT = mean_HMT_protZ, NNMT = proteo_CCLE_z[, NNMT_ensembl])

cor.test(CCLE_NNMTtoHMTprotein_forplot$HMT, CCLE_NNMTtoHMTprotein_forplot$NNMT)

pdf("graphics/CCLE_NNMTprotein_vs_HMTprotein_scatterplot.pdf",
    width = 2.5,
    height = 2.5)

ggplot(aes(y = NNMT, x = HMT), data = CCLE_NNMTtoHMTprotein_forplot) +
  geom_point(size = 0.4,
             alpha = 0.3) +
  theme_classic() + 
  theme(legend.position = "none",
        axis.text.x = element_text(colour = "black",
                                   size = 10),
        axis.text.y = element_text(colour = "black",
                                   size = 10),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        plot.margin = margin(t = 5, 
                             l = 5,
                             b = 5,
                             r = 8)) + 
  ylab("NNMT protein level (z-score)") + 
  xlab("Mean HMT protein (z-score)") + 
  scale_x_continuous(breaks = seq(-1.5, 1.5, 0.75),
                     limits = c(-1.5, 1.5),
                     expand = c(0,0),
                     labels = label_number(drop0trailing = TRUE)) + 
  scale_y_continuous(breaks = seq(-2, 3, 1),
                     limits = c(-2, 3),
                     expand = c(0,0)) +
  annotate(geom = "text",
           x = -0.6,
           y = 2.7,
           label = substitute(italic(r)~"= -0.238")) +
  annotate(geom = "text",
           x = -0.6,
           y = 2.3,
           label = substitute(italic(p)~"= 3.32 x"~10^-5)) +
  geom_smooth(method = "lm",
              se = FALSE,
              colour = "red")

dev.off()

#### NNMT protein to 1MNA Supplementary Figure 2B / Fig S2B ####

NNMTprotein_vec <- proteo_CCLE_z[, NNMT_ensembl]
names(NNMTprotein_vec) <- row.names(proteo_CCLE_z)

MNA_vec <- metab_CCLE_z[, "X1.methylnicotinamide"]
names(MNA_vec) <- row.names(metab_CCLE_z)

NNMT_MNA_forplot <- data.frame(NNMT = NNMTprotein_vec, MNA = MNA_vec[names(NNMTprotein_vec)])
NNMT_MNA_forplot <- NNMT_MNA_forplot[!apply(NNMT_MNA_forplot, 1, function(x){any(is.na(x))}), ]

cor.test(NNMT_MNA_forplot$NNMT, NNMT_MNA_forplot$MNA)$p.value

write.table(NNMT_MNA_forplot,
            file = "plot_data/Fig S2/Fig_S2B_plotdata.txt",
            row.names = TRUE,
            col.names = TRUE,
            sep = "\t")

pdf("graphics/CCLE_NNMTprotein_vs_MNA_scatterplot.pdf",
    width = 2.5,
    height = 2.5)

ggplot(aes(x = NNMT, y = MNA), data = NNMT_MNA_forplot) +
  geom_point(size = 0.4,
             alpha = 0.3) +
  theme_classic() + 
  theme(legend.position = "none",
        axis.text.x = element_text(colour = "black",
                                   size = 10),
        axis.text.y = element_text(colour = "black",
                                   size = 10),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10)) + 
  xlab("NNMT protein level (z-score)") + 
  ylab("Cell line 1MNA level (z-score)") + 
  scale_y_continuous(breaks = seq(-4, 4, 2),
                     limits = c(-4, 4),
                     expand = c(0,0)) + 
  scale_x_continuous(breaks = seq(-4, 4, 2),
                     limits = c(-4, 4),
                     expand = c(0,0)) +
  annotate(geom = "text",
           x = -1.2,
           y = 3.7,
           label = substitute(italic(r)~"= 0.551")) +
  annotate(geom = "text",
           x = -1.2,
           y = 2.9,
           label = substitute(italic(p)~"= 8.29 x"~10^-25)) +
  geom_smooth(method = "lm",
              se = FALSE,
              colour = "red")

dev.off()

#### CCLE NNMT protein to gene expression Supplementary Figure 2C / Fig S2C ####

CCLE_NNMTprotein_to_expression_forplot <- data.frame(NNMTprotein = NNMTprotein_vec, NNMTexpression = CCLE_NNMT_z_vec[names(NNMTprotein_vec)])
CCLE_NNMTprotein_to_expression_forplot <- CCLE_NNMTprotein_to_expression_forplot[!apply(CCLE_NNMTprotein_to_expression_forplot, 1, function(x){any(is.na(x))}), ]

cor.test(CCLE_NNMTprotein_to_expression_forplot$NNMTprotein, CCLE_NNMTprotein_to_expression_forplot$NNMTexpression)$p.value

write.table(CCLE_NNMTprotein_to_expression_forplot,
            file = "plot_data/Fig S2/Fig_S2C_plotdata.txt",
            row.names = TRUE,
            col.names = TRUE,
            sep = "\t")

pdf("graphics/CCLE_NNMTprotein_vs_NNMTexpression_scatterplot.pdf",
    width = 2.5,
    height = 2.5)

ggplot(aes(x = NNMTexpression, y = NNMTprotein), data = CCLE_NNMTprotein_to_expression_forplot) +
  geom_point(size = 0.4,
             alpha = 0.3) +
  theme_classic() + 
  theme(legend.position = "none",
        axis.text.x = element_text(colour = "black",
                                   size = 10),
        axis.text.y = element_text(colour = "black",
                                   size = 10),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10)) + 
  ylab("NNMT protein level (z-score)") + 
  xlab(substitute(italic(NNMT)~"expression (RNA-seq z)")) + 
  scale_y_continuous(breaks = seq(-3, 3, 1.5),
                     limits = c(-3, 3),
                     expand = c(0,0),
                     labels = label_number(drop0trailing = TRUE)) + 
  scale_x_continuous(breaks = seq(-3, 3, 1.5),
                     limits = c(-3, 3),
                     expand = c(0,0),
                     labels = label_number(drop0trailing = TRUE)) +
  annotate(geom = "text",
           x = -0.9,
           y = 2.6,
           label = substitute(italic(r)~"= 0.694")) +
  annotate(geom = "text",
           x = -0.9,
           y = 2,
           label = substitute(italic(p)~"= 1.33 x"~10^-43)) +
  geom_smooth(method = "lm",
              se = FALSE,
              colour = "red")

dev.off()

#### HMT expression to NNMT in CCLE ####

HMT_MOR_list <- lapply(CCLE_MOR_list, function(x){
  
  colSums(x[SET_HMTs$ensembl_gene_id, ])
  
})

HMT_z_list <- lapply(HMT_MOR_list, function(x){
  
  mean <- mean(log10(x + 1), na.rm = TRUE)
  sd <- sd(log10(x + 1), na.rm = TRUE)
  
  (log10(x + 1) - mean) / sd
  
}) 

HMT_z_vec <- do.call(c, HMT_z_list)
names(HMT_z_vec) <- str_extract(names(HMT_z_vec), pattern = "ACH.*")

NNMT_z_list <- lapply(lapply(CCLE_MOR_list, function(x){x[NNMT_ensembl, ]}), function(x){
  
  mean <- mean(log10(x + 1), na.rm = TRUE)
  sd <- sd(log10(x + 1), na.rm = TRUE)
  
  (log10(x + 1) - mean) / sd
  
}) 

NNMT_z_vec <- do.call(c, NNMT_z_list)
names(NNMT_z_vec) <- str_extract(names(NNMT_z_vec), pattern = "ACH.*")

CCLE_NNMT_HMT_df <- data.frame(cbind(NNMT_z_vec, HMT_z_vec))

cor.test(HMT_z_vec, NNMT_z_vec)$p.value

pdf("graphics/CCLE_NNMT_vs_HMT_scatterplot.pdf",
    width = 2.5,
    height = 2.5)

ggplot(aes(x = NNMT_z_vec, y = HMT_z_vec), data = CCLE_NNMT_HMT_df) +
  geom_point(size = 0.4,
             alpha = 0.3) +
  theme_classic() + 
  theme(legend.position = "none",
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10)) + 
  xlab(substitute(italic(NNMT)~"RNA-seq z-score")) + 
  ylab("HMTs (pooled) RNA-seq z-score") + 
  annotate(geom = "text",
           x = -2.3,
           y = 4,
           label = substitute(italic(r)~"= -0.253")) +
  annotate(geom = "text",
           x = -2.3,
           y = 3.4,
           label = substitute(italic(p)~"= 3.25 x"~10^-20)) + 
  annotate(geom = "text",
           x = 3.4,
           y = -2.8,
           label = "CCLE",
           colour = "dodgerblue4",
           size = 6,
           fontface = 2) +
  geom_smooth(method = "lm",
              se = FALSE,
              colour = "red")

dev.off()

#### MNA to all genes in CCLE Fig 1B ####

MNA_to_everything <- metabolite.to.allgenes.correlation(metabolite = "X1.methylnicotinamide",
                                   extrageneset = SET_HMTs$ensembl_gene_id,
                                   genesetname = "allHMTs")

MNA_to_everything <- data.frame(MNA_to_everything)
colnames(MNA_to_everything)[1] <- "rho"

MNA_to_everything <- MNA_to_everything[!is.na(MNA_to_everything$rho), ]

MNA_to_everything[, "FDR"] <- p.adjust(MNA_to_everything$genepvals, method = "BH")

MNA_to_HMTs <- MNA_to_everything[c(SET_HMTs$ensembl_gene_id, "allHMTs"), ]

selection_colours <- c("green", "magenta", "dodgerblue", "black")

write.table(MNA_to_everything[, c(1, 3)],
            file = "plot_data/Fig 1/Fig_1B_plotdata.txt",
            row.names = TRUE,
            col.names = TRUE,
            sep = "\t")

plot.selection.against.everything.volcano(everything_object = MNA_to_everything[, c(1,3)],
                                          selection_ensembl = act_rep_list,
                                          selection_colours = selection_colours,
                                          filename = "MNA_vs_HMTs",
                                          xlim = c(-0.5, 0.5),
                                          ylim = c(0, 40),
                                          yaxisbreaks = seq(0, 40, 10),
                                          xaxislabel = substitute("Pearson's"~italic(r)~" to 1MNA"),
                                          fileoutput = "png")

ggplot(aes(x = generesults, y = -log10(FDR)), data = MNA_to_everything) + 
  geom_jitter(width = 0.02, height = 0.02, size = 0.1, alpha = 0.6) + 
  theme_classic() +
  theme(axis.text.y = element_text(colour = "black"),
        axis.title.y = element_text(size = 10),
        axis.line.y = element_line(colour = "black")) +
  geom_jitter(data = MNA_to_HMTs, width = 0.01, height = 0.01, colour = "red", size = 2) + 
  coord_cartesian(xlim = c(-0.5, 0.5),
                  ylim = c(0, 25)) + 
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey") + 
  ylab(substitute(-log[10]~"(p value)")) +
  xlab (substitute("Spearman's"~rho))

ggplot(aes(y = abs(rho), x = rho), data = everything_df) + 
  geom_jitter(colour = everything_df$barcolour, width = 0.02, height = 0.02, size = 0.1, alpha = 0.6) + 
  theme_classic() +
  theme(axis.text.y = element_text(colour = "black"),
        axis.title.y = element_text(size = 10),
        axis.line.y = element_line(colour = "black")) +
  geom_jitter(data = selection_df, width = 0.01, height = 0.01, colour = selection_df$barcolour, size = 2) + 
  coord_cartesian(xlim = c(-0.6, 0.6),
                  ylim = c(0, 0.6)) + 
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey") + 
  ylab(substitute(-log[10]~"(p value)")) +
  xlab (substitute("Spearman's"~rho))

#### partial correlations MNA NNMT HMTs Supplementary Figure 2F / Fig S2F ####

MNAM <- metab_CCLE_z$X1.methylnicotinamide
names(MNAM) <- row.names(metab_CCLE_z)

SET_HMTs_list <- lapply(CCLE_MOR_list, function(x){
  
  sums <- colSums(x[SET_HMTs$ensembl_gene_id, ])
  mean <- mean(log10(sums), na.rm = TRUE)
  sd <- sd(log10(sums), na.rm = TRUE)
  
  (log10(sums) - mean) / sd
  
})

SET_HMTs_z_vec <- do.call(c, SET_HMTs_list)
names(SET_HMTs_z_vec) <- str_extract(names(SET_HMTs_z_vec), pattern = "ACH.*")

NNMT_HMTs_list <- lapply(CCLE_MOR_list, function(x){
  
  nnmt <- x[NNMT_ensembl, ]
  mean <- mean(log10(nnmt[nnmt != 0]), na.rm = TRUE)
  sd <- sd(log10(nnmt[nnmt != 0]), na.rm = TRUE)
  
  (log10(nnmt) - mean) / sd
  
})

NNMT_HMTs_z_vec <- do.call(c, NNMT_HMTs_list)
names(NNMT_HMTs_z_vec) <- str_extract(names(NNMT_HMTs_z_vec), pattern = "ACH.*")

matching_forpcor <- intersect(intersect(names(NNMT_HMTs_z_vec), names(SET_HMTs_z_vec)), names(MNAM))

NNMT_HMTs_z_vec <- NNMT_HMTs_z_vec[match(matching_forpcor, names(NNMT_HMTs_z_vec))]
SET_HMTs_z_vec <- SET_HMTs_z_vec[match(matching_forpcor, names(SET_HMTs_z_vec))]
MNAM <- MNAM[match(matching_forpcor, names(MNAM))]

alltogether_forpcor <- data.frame(cbind(NNMT_HMTs_z_vec, MNAM, SET_HMTs_z_vec))
alltogether_forpcor <- alltogether_forpcor[!is.infinite(alltogether_forpcor$NNMT_HMTs_z_vec), ]

pcor(alltogether_forpcor)$estimate
cor.test(alltogether_forpcor$NNMT_HMTs_z_vec, alltogether_forpcor$SET_HMTs_z_vec)
cor.test(alltogether_forpcor$MNAM, alltogether_forpcor$SET_HMTs_z_vec)
cor.test(alltogether_forpcor$MNAM, alltogether_forpcor$NNMT_HMTs_z_vec)

#### NNMT expression across cancers Supplementary Figure 2H / Fig S2H ####

TCGA_NNMTexpression_acrosscancers <- TCGA_MOR_across_cancers[NNMT_ensembl, ]

TCGA_NNMTexpression_acrosscancers_df <- data.frame(NNMTpseudocounts = TCGA_NNMTexpression_acrosscancers)
# its already limited to primary cancers

TCGA_NNMTexpression_acrosscancers_df[, "disease"] <- TCGA_LUT[match(row.names(TCGA_NNMTexpression_acrosscancers_df), TCGA_LUT$SAMPID), "cancer"]

TCGA_NNMTexpression_acrosscancers_geomean <- TCGA_NNMTexpression_acrosscancers_df %>% 
  group_by(disease) %>%
  summarise(geomean = gm_mean(NNMTpseudocounts))

TCGA_NNMTexpression_acrosscancers_df[, "disease"] <- factor(str_remove(TCGA_NNMTexpression_acrosscancers_df$disease, "^TCGA-"), levels = str_remove(unlist(TCGA_NNMTexpression_acrosscancers_geomean[order(TCGA_NNMTexpression_acrosscancers_geomean$geomean), "disease"]), "^TCGA-"))

write.table(TCGA_NNMTexpression_acrosscancers_df,
            file = "plot_data/Fig S2/Fig_S2H_plotdata.txt",
            row.names = TRUE,
            col.names = TRUE,
            sep = "\t")

pdf("graphics/TCGA_NNMTexpression_acrosscancers.pdf",
    width = 5,
    height = 2.5)

ggplot(aes(y = NNMTpseudocounts, x = disease), data = TCGA_NNMTexpression_acrosscancers_df) + 
  geom_boxplot(aes(fill = disease),
               outlier.size = 0.1) + 
  scale_y_log10(breaks = c(1, 10, 100, 1000, 10000, 100000, 1000000),
                expand = c(0, 0),
                limits = c(1, 1000000),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90,
                                   colour = "black",
                                   vjust = 0.5,
                                   hjust = 1),
        axis.text.y = element_text(colour = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10,
                                    hjust = 0.8),
        legend.position = "none",
        plot.margin = margin(t = 6, 
                              r = 3,
                              l = 3,
                              b = 3)) + 
  ylab(substitute(italic(NNMT)~"expression (pseudocounts)"))

dev.off()

#### HMT-NNNMT CORRELATIONS Fig 1E and Fig 1F ####

CCLE_allHMTs_with_everything <- geneset.with.everything(database = "CCLE",
                                                         geneset_ensembl = SET_HMTs$ensembl_gene_id,
                                                         basename = "allHMTs")

saveRDS(CCLE_allHMTs_with_everything, "output/CCLE_allHMTs_correlate_with_everything.rds")
CCLE_allHMTs_with_everything <- readRDS("output/CCLE_allHMTs_correlate_with_everything.rds")

pdf("output/CCLEtest.pdf",
    width = 2.5,
    height = 2.5)

plot.selection.against.everything(everything_object = CCLE_allHMTs_with_everything$medians,
                                  selection_ensembl = NNMT_ensembl,
                                  title = "")

dev.off()

# NNMT_with_everything_CCLE <-anything.with.everything(database = "CCLE",
#                                                      ensembl_gene = NNMT_ensembl,
#                                                      basename = "NNMT",
#                                                      iterations = 100,
#                                                      save = TRUE,
#                                                      extrageneset = SET_HMTs$ensembl_gene_id,
#                                                      genesetname = "allHMTs")

NNMT_with_everything_CCLEfull <- readRDS("output/CCLE_NNMT_correlate_with_everything.rds")
NNMT_with_everything_CCLE <- rowMedians(NNMT_with_everything_CCLEfull)
names(NNMT_with_everything_CCLE) <- row.names(NNMT_with_everything_CCLEfull)

# NNMT_with_everything_TCGA <-anything.with.everything(database = "TCGA",
#                                                      ensembl_gene = NNMT_ensembl,
#                                                      basename = "NNMT",
#                                                      iterations = 100,
#                                                      save = TRUE,
#                                                      extrageneset = SET_HMTs$ensembl_gene_id,
#                                                      genesetname = "allHMTs")

NNMT_with_everything_TCGAfull <- readRDS("output/TCGA_NNMT_correlate_with_everything.rds")
NNMT_with_everything_TCGA <- rowMedians(NNMT_with_everything_TCGAfull$medians)
names(NNMT_with_everything_TCGA) <- row.names(NNMT_with_everything_TCGAfull)

NNMT_with_everything_GTEX <-anything.with.everything(database = "GTEX",
                                                     ensembl_gene = NNMT_ensembl,
                                                     basename = "NNMT",
                                                     iterations = 100,
                                                     save = TRUE,
                                                     extrageneset = SET_HMTs$ensembl_gene_id,
                                                     genesetname = "allHMTs")

PEMT_with_everything_CCLE <-anything.with.everything(database = "CCLE",
                                                     ensembl_gene = PEMT_ensembl,
                                                     basename = "PEMT",
                                                     iterations = 100,
                                                     save = TRUE,
                                                     extrageneset = SET_HMTs$ensembl_gene_id,
                                                     genesetname = "allHMTs")

PEMT_with_everything_TCGA <-anything.with.everything(database = "TCGA",
                                                     ensembl_gene = PEMT_ensembl,
                                                     basename = "PEMT",
                                                     iterations = 100,
                                                     save = TRUE,
                                                     extrageneset = SET_HMTs$ensembl_gene_id,
                                                     genesetname = "allHMTs")

PEMT_with_everything_GTEX <-anything.with.everything(database = "GTEX",
                                                     ensembl_gene = PEMT_ensembl,
                                                     basename = "PEMT",
                                                     iterations = 100,
                                                     save = TRUE,
                                                     extrageneset = SET_HMTs$ensembl_gene_id,
                                                     genesetname = "allHMTs")

HMT_with_everything_CCLE <- geneset.with.everything(database = "CCLE",
                                                    geneset_ensembl = SET_HMTs$ensembl_gene_id,
                                                    basename = "allHMTs",
                                                    iterations = 100,
                                                    save = TRUE)

HMT_with_everything_TCGA <- geneset.with.everything(database = "TCGA",
                                                    geneset_ensembl = SET_HMTs$ensembl_gene_id,
                                                    basename = "allHMTs",
                                                    iterations = 100,
                                                    save = TRUE)

HMT_with_everything_GTEX <- geneset.with.everything(database = "GTEX",
                                                    geneset_ensembl = SET_HMTs$ensembl_gene_id,
                                                    basename = "allHMTs",
                                                    iterations = 100,
                                                    save = TRUE)

pdf("graphics/NNMT_vs_HMT_everything_CCLE.pdf",
    width = 1.8,
    height = 1.8)

plot.selection.against.everything(everything_object = HMT_with_everything_CCLE,
                                  selection_ensembl = NNMT_ensembl,
                                  title = "")

dev.off()

pdf("graphics/HMT_vs NNMT_everything_CCLE.pdf",
    width = 1.8,
    height = 1.8)

plot.selection.against.everything(everything_object = NNMT_with_everything_CCLE,
                                  selection_ensembl = "allHMTs",
                                  title = "")

dev.off()

pdf("graphics/NNMT_vs_HMT_everything_TCGA.pdf",
    width = 1.8,
    height = 1.8)

plot.selection.against.everything(everything_object = HMT_with_everything_TCGA,
                                  selection_ensembl = NNMT_ensembl,
                                  title = "")

dev.off()

pdf("graphics/HMT_vs_NNMT_everything_TCGA.pdf",
    width = 1.8,
    height = 1.8)

plot.selection.against.everything(everything_object = NNMT_with_everything_TCGA,
                                  selection_ensembl = "allHMTs",
                                  title = "")

dev.off()

#### HMTs NNMT scatterplot TCGA ####

TCGA_across_cancers_HMT_NNMT <- cbind(colSums(TCGA_MOR_across_cancers[SET_HMTs$ensembl_gene_id, ]), TCGA_MOR_across_cancers[NNMT_ensembl, ])
colnames(TCGA_across_cancers_HMT_NNMT) <- c("allHMTs", "NNMT")

TCGA_across_cancers_HMT_NNMT <- TCGA_across_cancers_HMT_NNMT[row.names(TCGA_across_cancers_HMT_NNMT) %in% TCGA_LUT$SAMPID, ]

TCGA_across_cancers_HMT_NNMT <- data.frame(cbind(TCGA_across_cancers_HMT_NNMT, TCGA_LUT[match(row.names(TCGA_across_cancers_HMT_NNMT), TCGA_LUT$SAMPID), c("race", "gender", "tumour_stage", "days_to_birth", "sequencing_centre", "cancer")]))

TCGA_across_cancers_HMT_NNMT_residuals <- data.frame(NNMT_residuals = lm(log10(NNMT + 1) ~ cancer + as.factor(race) + as.numeric(days_to_birth) + as.factor(gender) + as.factor(tumour_stage) + as.factor(sequencing_centre), data = TCGA_across_cancers_HMT_NNMT)$residuals,
                                                           HMTresiduals = lm(log10(allHMTs + 1) ~ cancer + as.factor(race) + as.numeric(days_to_birth) + as.factor(gender) + as.factor(tumour_stage) + as.factor(sequencing_centre), data = TCGA_across_cancers_HMT_NNMT)$residuals)

cor.test(TCGA_across_cancers_HMT_NNMT_residuals$NNMT_residuals, TCGA_across_cancers_HMT_NNMT_residuals$HMTresiduals)$p.value

pdf("graphics/TCGA_allcancers_HMT_NNMT_residual_scatterplot.pdf",
    height = 2.5,
    width = 2.5)

ggplot(aes(x = HMTresiduals, y = NNMT_residuals), data = TCGA_across_cancers_HMT_NNMT_residuals) + 
  geom_point(size = 0.1, alpha = 0.1) + 
  theme_classic() + 
  theme(axis.text.x = element_text(colour = "black",
                                   size = 10),
        axis.text.y = element_text(colour = "black",
                                   size = 10)) +
  geom_smooth(method = "lm",
              se = FALSE,
              colour = "red") +
  xlab("total HMTs expression (resid.)") + 
  ylab(italic(NNMT)~"expression (residuals)") +
  annotate(geom = "text",
           x = -0.25,
           y = -2,
           label = substitute(italic(r)~"= -0.383")) + 
  annotate(geom = "text",
           x = -0.25,
           y = -2.5,
           label = substitute(italic(p)~"= 0"))

dev.off()

#### individual cancer/tissue types for Supplementary File 1 and Supplementary File 2 / File S1 and File S2 ####

TCGA_SET_HMTs_bycancer <- analyse.by.subtype(database = "TCGA",
                                             geneset_ensembl = SET_HMTs$ensembl_gene_id,
                                             genebasename = "allHMTs")

saveRDS(TCGA_SET_HMTs_bycancer, "output/TCGA_SET_HMTs_with_everything_bycancer.rds")
# TCGA_SET_HMTs_bycancer <- readRDS("output/TCGA_SET_HMTs_with_everything_bycancer.rds")

TCGA_NNMT_bycancer <- analyse.by.subtype(database = "TCGA",
                                             singlegene = NNMT_ensembl,
                                             genebasename = NNMT_ensembl,
                                             extrageneset_ensembl = SET_HMTs$ensembl_gene_id,
                                             extragenesetbasename = "allHMTs")

saveRDS(TCGA_NNMT_bycancer, "output/TCGA_NNMT_with_everything_bycancer.rds")
# TCGA_NNMT_bycancer <- readRDS("output/TCGA_NNMT_with_everything_bycancer.rds")

CCLE_allHMTs_bycancer <- analyse.by.subtype(database = "CCLE",
                                             geneset_ensembl = SET_HMTs$ensembl_gene_id,
                                             genebasename = "allHMTs")

saveRDS(CCLE_allHMTs_bycancer, "output/CCLE_allHMTs_with_everything_bycancer.rds")
# CCLE_allHMTs_bycancer <- readRDS("output/CCLE_allHMTs_with_everything_bycancer.rds")

CCLE_NNMT_bycancer <- analyse.by.subtype(database = "CCLE",
                                         singlegene = NNMT_ensembl,
                                         genebasename = NNMT_ensembl,
                                         extrageneset_ensembl = SET_HMTs$ensembl_gene_id,
                                         extragenesetbasename = "allHMTs")

saveRDS(CCLE_NNMT_bycancer, "output/CCLE_NNMT_with_everything_bycancer.rds")
# CCLE_NNMT_bycancer <- readRDS("output/CCLE_NNMT_with_everything_bycancer.rds")

#### compute reciprocal scores for TCGA cancer types for Fig 1I #### 

TCGA_cancer_reciprocal_scores <- sapply(names(TCGA_SET_HMTs_bycancer), function(tempcancertype){

  NNMT_in_HMTs <- which(row.names(TCGA_SET_HMTs_bycancer[[tempcancertype]][order(TCGA_SET_HMTs_bycancer[[tempcancertype]]$allcorrs), ]) == NNMT_ensembl)
  HMTs_in_NNMT <- which(row.names(TCGA_NNMT_bycancer[[tempcancertype]][order(TCGA_NNMT_bycancer[[tempcancertype]]$allcorrs), ]) == "allHMTs")
  
  reciprocal_score <- (NNMT_in_HMTs^2 + HMTs_in_NNMT^2)
  
  return(c(NNMT_in_HMTs, HMTs_in_NNMT, reciprocal_score))
  
})

TCGA_cancer_reciprocal_scores_df <- data.frame(t(TCGA_cancer_reciprocal_scores))
colnames(TCGA_cancer_reciprocal_scores_df) <- c("NNMT_in_HMTs", "HMTs_in_NNMT", "reciprocal_score")

TCGA_cancer_reciprocal_scores_df[, "cancer"] <- str_remove(row.names(TCGA_cancer_reciprocal_scores_df), "^TCGA-")

TCGA_cancer_reciprocal_scores_df[, "HMTs_in_NNMT_percent"] <- (TCGA_cancer_reciprocal_scores_df$HMTs_in_NNMT * 100)/60489
TCGA_cancer_reciprocal_scores_df[, "NNMT_in_HMTs_percent"] <- (TCGA_cancer_reciprocal_scores_df$NNMT_in_HMTs * 100)/60489

write.table(TCGA_cancer_reciprocal_scores_df,
            file = "plot_data/Fig 1/Fig_1I_plotdata.txt",
            row.names = TRUE,
            col.names = TRUE,
            sep = "\t")

# plot out

pdf("graphics/TCGA-cancers-reciprocal_scores.pdf",
    width = 2,
    height = 2)

ggplot(aes(x = NNMT_in_HMTs_percent, y = HMTs_in_NNMT_percent, label = cancer), data = TCGA_cancer_reciprocal_scores_df) + 
  geom_rect(aes(xmax = 2.5,
                xmin = 0,
                ymax = 2.5,
                ymin = 0),
                fill = "white",
            linetype = "dashed",
            color = "grey") +
  geom_point(aes(size = 10/log10(reciprocal_score)),
             colour = "red",
             alpha = 0.7) + 
  scale_x_log10(breaks = c(0.01, 0.1, 1, 10, 100),
                expand = c(0, 0),
                limits = c(0.001, 100),
                labels = label_number(drop0trailing = TRUE)) +
  scale_y_log10(breaks = c(0.01, 0.1, 1, 10, 100),
                expand = c(0, 0),
                limits = c(0.001, 100),
                labels = label_number(drop0trailing = TRUE)) +
  geom_label_repel(size = 1.5,
                   max.overlaps = 40,
                   label.size = 0.1,
                   box.padding = 0.1,
                   label.padding = 0.1) + 
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(colour = "black",
                                   size = 8),
        axis.text.y = element_text(colour = "black",
                                   size = 8),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        plot.margin = margin(t = 3,
                             r = 8,
                             b = 3,
                             l = 3)) + 
  xlab(substitute(italic(NNMT)~"%"^"ile"~"in HMTs correlations")) +
  ylab(substitute("HMTs"~"%"^"ile"~"in"~italic(NNMT)~"correlations"))

dev.off()

#### plot reciprocal scores vs NNMT expressionc for Supplementary Figure 2I / Fig S2I #### 

NNMT_expression_vs_score <- data.frame(reciprocal_score = TCGA_cancer_reciprocal_scores_df$reciprocal_score, NNMTexpression = unlist(TCGA_NNMTexpression_acrosscancers_geomean[match(TCGA_cancer_reciprocal_scores_df$cancer, str_remove(TCGA_NNMTexpression_acrosscancers_geomean$disease, pattern = "^TCGA-")), "geomean"]))
NNMT_expression_vs_score[, "cancer"] <- TCGA_cancer_reciprocal_scores_df$cancer

cor.test(log10(NNMT_expression_vs_score$reciprocal_score), log10(NNMT_expression_vs_score$NNMTexpression))

pdf("graphics/TCGA_cancer_NNMTexpression_vs_reciprocalscore_scatterplot.pdf",
    width = 2.5,
    height = 2.5)

ggplot(aes(x = NNMTexpression, y = reciprocal_score), data = NNMT_expression_vs_score) + 
  geom_point() + 
  scale_x_log10(breaks = c(10, 100, 1000, 10000, 10000, 100000),
                expand = c(0,0),
                limits = c(10, 100000),
                labels = trans_format("log10", math_format(10^.x))) + 
  scale_y_continuous(trans = reverselog_trans(10), 
                     breaks = c(10, 1000, 100000, 10000000, 1000000000),
                     expand = c(0,0),
                     limits = c(1000000000, 10),
                     labels = trans_format("log10", math_format(10^.x))) +
  geom_smooth(method = "lm",
              se = FALSE) +
  theme_classic() +
  theme(axis.text.x = element_text(colour = "black",
                                   size = 10),
        axis.text.y = element_text(colour = "black",
                                   size = 10),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10)) + 
  ylab(substitute("HMTs-"~italic(NNMT)~"reciprocal score")) + 
  xlab(substitute(italic(NNMT)~"expression (pseudocounts)")) + 
  annotate(geom = "text",
           x = 100,
           y = 100,
           label = substitute(italic(r)~"= -0.405")) +
  annotate(geom = "text",
           x = 100,
           y = 1000,
           label = substitute(italic(p)~"= 0.0193"))
  
dev.off()

TCGA_NNMTexpression_acrosscancers_df[, "reciprocal_score"] <- NNMT_expression_vs_score[match(TCGA_NNMTexpression_acrosscancers_df$disease, NNMT_expression_vs_score$cancer), "reciprocal_score"]

cor.test(log10(TCGA_NNMTexpression_acrosscancers_df$NNMTpseudocounts + 1), log10(TCGA_NNMTexpression_acrosscancers_df$reciprocal_score))$estimate^2
cor.test(log10(TCGA_NNMTexpression_acrosscancers_df$NNMTpseudocounts + 1), log10(TCGA_NNMTexpression_acrosscancers_df$reciprocal_score))$p.value

write.table(TCGA_NNMTexpression_acrosscancers_df,
            file = "plot_data/Fig S2/Fig_S2I_plotdata.txt",
            row.names = TRUE,
            col.names = TRUE,
            sep = "\t")

pdf("graphics/TCGA_cancer_NNMTexpression_vs_reciprocalscore_violin.pdf",
    height = 2.5,
    width = 3)

ggplot(data = TCGA_NNMTexpression_acrosscancers_df, aes(y = NNMTpseudocounts, x = reciprocal_score)) +
  geom_violin(aes(group = disease, fill = disease), 
              show.legend = FALSE, 
              scale = 1, 
              width = 0.3, 
              position = "identity", 
              lwd = 0.2) +
  scale_x_continuous(trans = reverselog_trans(10), 
                breaks = c(10, 1000, 100000, 10000000, 1000000000),
                expand = c(0,0),
                limits = c(1000000000, 10),
                labels = trans_format("log10", math_format(10^.x))) + 
  scale_y_log10(breaks = c(1, 10, 100, 1000, 1000, 10000,  100000, 1000000),
                expand = c(0,0),
                limits = c(1, 1000000),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_classic() +
  theme(axis.text.x = element_text(colour = "black",
                                   size = 10),
        axis.text.y = element_text(colour = "black",
                                   size = 10),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10,
                                    hjust = 0.8)) + 
  xlab(substitute("HMTs-"~italic(NNMT)~"reciprocal score")) + 
  ylab(substitute(italic(NNMT)~"expression (pseudocounts)")) + 
  geom_smooth(method = "lm",
              se = FALSE) + 
  annotate(geom = "text",
           x = 10000000,
           y = 600000,
           label = substitute(R^2~"= 0.0631")) +
  annotate(geom = "text",
           x = 10000000,
           y = 200000,
           label = substitute(italic(p)~"= 3.38 x"~10^-141))

dev.off()
  
#### compute reciprocal scores for TCGA cancer types for PEMT for Fig 2G ####

TCGA_PEMT_cancer_reciprocal_scores <- sapply(names(TCGA_SET_HMTs_bycancer), function(tempcancertype){

  PEMT_in_HMTs <- which(row.names(TCGA_SET_HMTs_bycancer[[tempcancertype]][order(TCGA_SET_HMTs_bycancer[[tempcancertype]]$allcorrs), ]) == PEMT_ensembl)
  HMTs_in_PEMT <- which(row.names(TCGA_PEMT_bycancer[[tempcancertype]][order(TCGA_PEMT_bycancer[[tempcancertype]]$allcorrs), ]) == "allHMTs")
  TCGA_SET_HMTs_bycancer[[tempcancertype]][PEMT_ensembl, ]
  TCGA_PEMT_bycancer[[tempcancertype]]["allHMTs", ]
  reciprocal_score <- (PEMT_in_HMTs^2 + HMTs_in_PEMT^2)
  
  return(c(PEMT_in_HMTs, HMTs_in_PEMT, reciprocal_score))
  
})

TCGA_PEMT_cancer_reciprocal_scores_df <- data.frame(t(TCGA_PEMT_cancer_reciprocal_scores))
colnames(TCGA_PEMT_cancer_reciprocal_scores_df) <- c("PEMT_in_HMTs", "HMTs_in_PEMT", "reciprocal_score")

TCGA_PEMT_cancer_reciprocal_scores_df[, "cancer"] <- str_remove(row.names(TCGA_PEMT_cancer_reciprocal_scores_df), "^TCGA-")

TCGA_PEMT_cancer_reciprocal_scores_df[, "HMTs_in_PEMT_percent"] <- (TCGA_PEMT_cancer_reciprocal_scores_df$HMTs_in_PEMT * 100)/60489
TCGA_PEMT_cancer_reciprocal_scores_df[, "PEMT_in_HMTs_percent"] <- (TCGA_PEMT_cancer_reciprocal_scores_df$PEMT_in_HMTs * 100)/60489

saveRDS(TCGA_PEMT_cancer_reciprocal_scores_df, "output/TCGA_PEMT_cancer_reciprocal_scores_df.rds")
# TCGA_PEMT_cancer_reciprocal_scores_df <- readRDS("output/TCGA_PEMT_cancer_reciprocal_scores_df.rds")

TCGA_PEMT_cancer_reciprocal_scores_df[TCGA_PEMT_cancer_reciprocal_scores_df$HMTs_in_PEMT_percent > 2.5 | TCGA_PEMT_cancer_reciprocal_scores_df$PEMT_in_HMTs_percent > 2.5, "cancer"] <- ""

write.table(TCGA_PEMT_cancer_reciprocal_scores_df,
            file = "plot_data/Fig 2/Fig_2G_plotdata.txt",
            row.names = TRUE,
            col.names = TRUE,
            sep = "\t")

# plot out

pdf("graphics/TCGA-PEMT_cancers-reciprocal_scores.pdf",
    width = 2,
    height = 2)

ggplot(aes(x = PEMT_in_HMTs_percent, y = HMTs_in_PEMT_percent, label = cancer), data = TCGA_PEMT_cancer_reciprocal_scores_df) + 
  geom_rect(aes(xmax = 2.5,
                xmin = 0,
                ymax = 2.5,
                ymin = 0),
            fill = "white",
            linetype = "dashed",
            color = "grey") +
  geom_point(aes(size = 10/log10(reciprocal_score)),
             colour = "red",
             alpha = 0.7) + 
  scale_x_log10(breaks = c(0.01, 0.1, 1, 10, 100),
                expand = c(0, 0),
                limits = c(0.001, 100),
                labels = label_number(drop0trailing = TRUE)) +
  scale_y_log10(breaks = c(0.01, 0.1, 1, 10, 100),
                expand = c(0, 0),
                limits = c(0.001, 100),
                labels = label_number(drop0trailing = TRUE)) +
  geom_label_repel(size = 1.5,
                   max.overlaps = 40,
                   label.size = 0.1,
                   box.padding = 0.1,
                   label.padding = 0.1) + 
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(colour = "black",
                                   size = 8),
        axis.text.y = element_text(colour = "black",
                                   size = 8),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        plot.margin = margin(t = 3,
                             r = 8,
                             b = 3,
                             l = 3)) + 
  xlab(substitute(italic(PEMT)~"%"^"ile"~"in HMTs correlations")) +
  ylab(substitute("HMTs"~"%"^"ile"~"in"~italic(PEMT)~"correlations"))

dev.off()

# compute reciprocal scores for CCLE cancer types

CCLE_cancer_reciprocal_scores <- sapply(names(CCLE_allHMTs_bycancer), function(tempcancertype){
  
  NNMT_in_HMTs <- which(row.names(CCLE_allHMTs_bycancer[[tempcancertype]][order(CCLE_allHMTs_bycancer[[tempcancertype]]$allcorrs), ]) == NNMT_ensembl)
  HMTs_in_NNMT <- which(row.names(CCLE_NNMT_bycancer[[tempcancertype]][order(CCLE_NNMT_bycancer[[tempcancertype]]$allcorrs), ]) == "allHMTs")
  
  reciprocal_score <- (NNMT_in_HMTs^2 + HMTs_in_NNMT^2)
  
  return(c(NNMT_in_HMTs, HMTs_in_NNMT, reciprocal_score))
  
})

CCLE_cancer_reciprocal_scores_df <- data.frame(t(CCLE_cancer_reciprocal_scores))
colnames(CCLE_cancer_reciprocal_scores_df) <- c("NNMT_in_HMTs", "HMTs_in_NNMT", "reciprocal_score")

CCLE_cancer_reciprocal_scores_df[, "cancer"] <- str_remove(row.names(CCLE_cancer_reciprocal_scores_df), "^CCLE-")

CCLE_cancer_reciprocal_scores_df[, "HMTs_in_NNMT_percent"] <- (CCLE_cancer_reciprocal_scores_df$HMTs_in_NNMT * 100)/60489
CCLE_cancer_reciprocal_scores_df[, "NNMT_in_HMTs_percent"] <- (CCLE_cancer_reciprocal_scores_df$NNMT_in_HMTs * 100)/60489

# expand slightly to label sarcoma, which is marginal
CCLE_cancer_reciprocal_scores_df[CCLE_cancer_reciprocal_scores_df$HMTs_in_NNMT_percent > 3 | CCLE_cancer_reciprocal_scores_df$NNMT_in_HMTs_percent > 3, "cancer"] <- ""

write.table(CCLE_cancer_reciprocal_scores_df,
            file = "plot_data/Fig S2/Fig_S2G_plotdata.txt",
            row.names = FALSE,
            col.names = TRUE,
            sep = "\t")

# plot out

pdf("graphics/CCLE-cancers-reciprocal_scores.pdf",
    width = 2.5,
    height = 2.5)

ggplot(aes(x = NNMT_in_HMTs_percent, y = HMTs_in_NNMT_percent, label = cancer), data = CCLE_cancer_reciprocal_scores_df) + 
  geom_rect(aes(xmax = 2.5,
                xmin = 0,
                ymax = 2.5,
                ymin = 0),
            fill = "white",
            linetype = "dashed",
            color = "grey") +
  geom_point(aes(size = 10/log10(reciprocal_score)),
             colour = "red",
             alpha = 0.7) + 
  scale_x_log10(breaks = c(0.01, 0.1, 1, 10, 100),
                expand = c(0, 0),
                limits = c(0.001, 100),
                labels = label_number(drop0trailing = TRUE)) +
  scale_y_log10(breaks = c(0.01, 0.1, 1, 10, 100),
                expand = c(0, 0),
                limits = c(0.001, 100),
                labels = label_number(drop0trailing = TRUE)) +
  geom_label_repel(size = 1.5,
                   max.overlaps = 40,
                   label.size = 0.1,
                   box.padding = 0.1,
                   label.padding = 0.1) + 
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(colour = "black",
                                   size = 10),
        axis.text.y = element_text(colour = "black",
                                   size = 10),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        plot.margin = margin(t = 3,
                             r = 8,
                             b = 3,
                             l = 3)) + 
  xlab(substitute(italic(NNMT)~"%"^"ile"~"in HMTs correlations")) +
  ylab(substitute("HMTs"~"%"^"ile"~"in"~italic(NNMT)~"correlations"))

dev.off()

TCGA_PEMT_bycancer <- analyse.by.subtype(database = "TCGA",
                                             singlegene = PEMT_ensembl,
                                             genebasename = PEMT_ensembl,
                                             extrageneset_ensembl = SET_HMTs$ensembl_gene_id,
                                             extragenesetbasename = "allHMTs")

saveRDS(TCGA_PEMT_bycancer, "output/TCGA_PEMT_with_everything_bycancer.rds")
# TCGA_PEMT_bycancer <- readRDS("output/TCGA_PEMT_with_everything_bycancer.rds")

GTEX_SET_HMTs_bytissue <- analyse.by.subtype(database = "GTEX",
                                             geneset_ensembl = SET_HMTs$ensembl_gene_id,
                                             genebasename = "allHMTs")

saveRDS(GTEX_SET_HMTs_bytissue, "output/GTEX_SET_HMTs_with_everything_bytissue.rds")
# GTEX_SET_HMTs_bytissue <- readRDS("output/GTEX_SET_HMTs_with_everything_bytissue.rds")

GTEX_NNMT_bytissue <- analyse.by.subtype(database = "GTEX",
                                             singlegene = NNMT_ensembl,
                                             genebasename = NNMT_ensembl,
                                             extrageneset_ensembl = SET_HMTs$ensembl_gene_id,
                                             extragenesetbasename = "allHMTs")

saveRDS(GTEX_NNMT_bytissue, "output/GTEX_NNMT_with_everything_bytissue.rds")
# GTEX_NNMT_bytissue <- readRDS("output/GTEX_NNMT_with_everything_bytissue.rds")

GTEX_PEMT_bytissue <- analyse.by.subtype(database = "GTEX",
                                         singlegene = PEMT_ensembl,
                                         genebasename = PEMT_ensembl,
                                         extrageneset_ensembl = SET_HMTs$ensembl_gene_id,
                                         extragenesetbasename = "allHMTs")

saveRDS(GTEX_PEMT_bytissue, "output/GTEX_PEMT_with_everything_bytissue.rds")
# GTEX_PEMT_bytissue <- readRDS("output/GTEX_PEMT_with_everything_bytissue.rds")

#### compute reciprocal scores for GTEX tissues Fig 2C ####

GTEX_tissue_HMT_PEMT_reciprocal_scores <- sapply(names(GTEX_SET_HMTs_bytissue), function(temptissue){
  
  PEMT_in_HMTs <- which(row.names(GTEX_SET_HMTs_bytissue[[temptissue]][order(GTEX_SET_HMTs_bytissue[[temptissue]]$allcorrs), ]) == PEMT_ensembl)
  HMTs_in_PEMT <- which(row.names(GTEX_PEMT_bytissue[[temptissue]][order(GTEX_PEMT_bytissue[[temptissue]]$allcorrs), ]) == "allHMTs")
  
  reciprocal_score <- (PEMT_in_HMTs^2 + HMTs_in_PEMT^2)
  
  return(c(PEMT_in_HMTs, HMTs_in_PEMT, reciprocal_score))
  
})

GTEX_tissue_HMT_PEMT_reciprocal_scores_df <- data.frame(t(GTEX_tissue_HMT_PEMT_reciprocal_scores))
colnames(GTEX_tissue_HMT_PEMT_reciprocal_scores_df) <- c("PEMT_in_HMTs", "HMTs_in_PEMT", "reciprocal_score")

GTEX_tissue_HMT_PEMT_reciprocal_scores_df[, "tissue"] <- convert.GTEX.tissue.labels.for.plot(str_remove(row.names(GTEX_tissue_HMT_PEMT_reciprocal_scores_df), "^GTEX-"), supershort = TRUE)

GTEX_tissue_HMT_PEMT_reciprocal_scores_df[, "HMTs_in_PEMT_percent"] <- (GTEX_tissue_HMT_PEMT_reciprocal_scores_df$HMTs_in_PEMT * 100)/56201
GTEX_tissue_HMT_PEMT_reciprocal_scores_df[, "PEMT_in_HMTs_percent"] <- (GTEX_tissue_HMT_PEMT_reciprocal_scores_df$PEMT_in_HMTs * 100)/56201

saveRDS(GTEX_tissue_HMT_PEMT_reciprocal_scores_df, "output/GTEX_tissue_HMT_PEMT_reciprocal_scores_df.rds")
# GTEX_tissue_HMT_PEMT_reciprocal_scores_df <- readRDS("output/GTEX_tissue_HMT_PEMT_reciprocal_scores_df.rds")

write.table(GTEX_tissue_HMT_PEMT_reciprocal_scores_df,
            file = "plot_data/Fig 2/Fig_2C_plotdata.txt",
            row.names = TRUE,
            col.names = TRUE,
            sep = "\t")

GTEX_tissue_HMT_PEMT_reciprocal_scores_df[GTEX_tissue_HMT_PEMT_reciprocal_scores_df$HMTs_in_PEMT_percent > 2.5 | GTEX_tissue_HMT_PEMT_reciprocal_scores_df$PEMT_in_HMTs_percent > 2.5, "tissue"] <- ""

GTEX_PEMT_HMT_strongtissues <- row.names(GTEX_tissue_HMT_PEMT_reciprocal_scores_df[GTEX_tissue_HMT_PEMT_reciprocal_scores_df$HMTs_in_PEMT_percent < 2.5 & GTEX_tissue_HMT_PEMT_reciprocal_scores_df$PEMT_in_HMTs_percent < 2.5, ])

# plot out as bubble plot

pdf("graphics/GTEX-tissues-reciprocal_scores.pdf",
    width = 2,
    height = 2)

ggplot(aes(x = PEMT_in_HMTs_percent, y = HMTs_in_PEMT_percent, label = tissue), data = GTEX_tissue_HMT_PEMT_reciprocal_scores_df) + 
  geom_rect(aes(xmax = 2.5,
                xmin = 0,
                ymax = 2.5,
                ymin = 0),
            fill = "white",
            linetype = "dashed",
            color = "grey") +
  geom_point(aes(size = 10/log10(reciprocal_score)),
             colour = "red",
             alpha = 0.7) + 
  scale_x_log10(breaks = c(0.01, 0.1, 1, 10, 100),
                expand = c(0, 0),
                limits = c(0.001, 100),
                labels = label_number(drop0trailing = TRUE)) +
  scale_y_log10(breaks = c(0.01, 0.1, 1, 10, 100),
                expand = c(0, 0),
                limits = c(0.001, 100),
                labels = label_number(drop0trailing = TRUE)) +
  geom_label_repel(size = 1.5,
                   max.overlaps = 40,
                   label.size = 0.1,
                   box.padding = 0.1,
                   label.padding = 0.1) + 
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(colour = "black",
                                   size = 8),
        axis.text.y = element_text(colour = "black",
                                   size = 8),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        plot.margin = margin(t = 5,
                             r = 8,
                             b = 3,
                             l = 3)) + 
  xlab(substitute(italic(PEMT)~"%"^"ile"~"in HMTs correlations")) +
  ylab(substitute("HMTs"~"%"^"ile"~"in"~italic(PEMT)~"correlations"))

dev.off()

#### compute reciprocal scores for NNMT for GTEX tissues Fig 2A and Fig 2B ####

GTEX_tissue_HMT_NNMT_reciprocal_scores <- sapply(names(GTEX_SET_HMTs_bytissue), function(temptissue){
  
  NNMT_in_HMTs <- which(row.names(GTEX_SET_HMTs_bytissue[[temptissue]][order(GTEX_SET_HMTs_bytissue[[temptissue]]$allcorrs), ]) == NNMT_ensembl)
  HMTs_in_NNMT <- which(row.names(GTEX_NNMT_bytissue[[temptissue]][order(GTEX_NNMT_bytissue[[temptissue]]$allcorrs), ]) == "allHMTs")
  
  reciprocal_score <- (NNMT_in_HMTs^2 + HMTs_in_NNMT^2)
  
  return(c(NNMT_in_HMTs, HMTs_in_NNMT, reciprocal_score))
  
})

GTEX_tissue_HMT_NNMT_reciprocal_scores_df <- data.frame(t(GTEX_tissue_HMT_NNMT_reciprocal_scores))
colnames(GTEX_tissue_HMT_NNMT_reciprocal_scores_df) <- c("NNMT_in_HMTs", "HMTs_in_NNMT", "reciprocal_score")

GTEX_tissue_HMT_NNMT_reciprocal_scores_df[, "tissue"] <- convert.GTEX.tissue.labels.for.plot(str_remove(row.names(GTEX_tissue_HMT_NNMT_reciprocal_scores_df), "^GTEX-"), supershort = TRUE)

GTEX_tissue_HMT_NNMT_reciprocal_scores_df[, "HMTs_in_NNMT_percent"] <- (GTEX_tissue_HMT_NNMT_reciprocal_scores_df$HMTs_in_NNMT * 100)/56201
GTEX_tissue_HMT_NNMT_reciprocal_scores_df[, "NNMT_in_HMTs_percent"] <- (GTEX_tissue_HMT_NNMT_reciprocal_scores_df$NNMT_in_HMTs * 100)/56201

dir.create("plot_data/Fig 2")

write.table(GTEX_tissue_HMT_NNMT_reciprocal_scores_df,
            file = "plot_data/Fig 2/Fig_2A_plotdata.txt",
            row.names = TRUE,
            col.names = TRUE,
            sep = "\t")

# will restrict labelling to 6 with lowest reciprocal scores

GTEX_tissue_HMT_NNMT_reciprocal_scores_df[GTEX_tissue_HMT_NNMT_reciprocal_scores_df$HMTs_in_NNMT_percent > 2.5 | GTEX_tissue_HMT_NNMT_reciprocal_scores_df$NNMT_in_HMTs_percent > 2.5, "tissue"] <- ""

# plot out as bubble plot

pdf("graphics/GTEX-tissues-NNMT-reciprocal_scores.pdf",
    width = 2,
    height = 2)

ggplot(aes(x = NNMT_in_HMTs_percent, y = HMTs_in_NNMT_percent, label = tissue), data = GTEX_tissue_HMT_NNMT_reciprocal_scores_df) + 
  geom_rect(aes(xmax = 2.5,
                xmin = 0,
                ymax = 2.5,
                ymin = 0),
            fill = "white",
            linetype = "dashed",
            color = "grey") +
  geom_point(aes(size = 10/log10(reciprocal_score)),
             colour = "red",
             alpha = 0.7) + 
  scale_x_log10(breaks = c(0.01, 0.1, 1, 10, 100),
                expand = c(0, 0),
                limits = c(0.001, 100),
                labels = label_number(drop0trailing = TRUE)) +
  scale_y_log10(breaks = c(0.01, 0.1, 1, 10, 100),
                expand = c(0, 0),
                limits = c(0.001, 100),
                labels = label_number(drop0trailing = TRUE)) +
  geom_label_repel(size = 1.5,
                   max.overlaps = 40,
                   label.size = 0.1,
                   box.padding = 0.1,
                   label.padding = 0.1) +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(colour = "black",
                                   size = 8),
        axis.text.y = element_text(colour = "black",
                                   size = 8),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        plot.margin = margin(t = 5,
                             r = 8,
                             b = 3,
                             l = 3)) + 
  xlab(substitute(italic(NNMT)~"%"^"ile"~"in HMTs correlations")) +
  ylab(substitute("HMTs"~"%"^"ile"~"in"~italic(NNMT)~"correlations"))

dev.off()

GTEX_NNMT_with_everything_bytissue <- readRDS("output/GTEX_NNMT_with_everything_bytissue.rds")

names(GTEX_NNMT_with_everything_bytissue) <- convert.GTEX.tissue.labels.for.plot(names(GTEX_NNMT_with_everything_bytissue), supershort = TRUE)

plot.selection.against.everything.volcano.facet(everything_object = GTEX_NNMT_with_everything_bytissue,
                                        selection_ensembl = act_rep_list,
                                        selection_colours = c("green", "magenta", "dodgerblue", "black"),
                                        filename = "GTEX_NNMT_HMTs_bytissue",
                                        fileoutput = "png")

plot.selection.against.everything.volcano.facet(everything_object = GTEX_PEMT_bytissue,
                                                selection_ensembl = act_rep_list,
                                                selection_colours = c("green", "magenta", "dodgerblue", "black"),
                                                filename = "GTEX_PEMT_HMTs_bytissue",
                                                fileoutput = "png")

plot.selection.against.everything.volcano.full.plots(everything_object = GTEX_PEMT_bytissue,
                                                selection_ensembl = act_rep_list,
                                                selection_colours = c("green", "magenta", "dodgerblue", "black"),
                                                filename_prefix = "GTEX_PEMT_HMTs",
                                                fileoutput = "png")

plot.selection.against.everything.volcano.full.plots(everything_object = GTEX_NNMT_bytissue,
                                                     selection_ensembl = act_rep_list,
                                                     selection_colours = c("green", "magenta", "dodgerblue", "black"),
                                                     filename_prefix = "GTEX_NNMT_HMTs",
                                                     fileoutput = "png")

plot.selection.against.everything.volcano.facet(everything_object = GTEX_SET_HMTs_bytissue,
                                                selection_ensembl = PEMT_ensembl,
                                                selection_colours = "red",
                                                filename = "GTEX_HMTs_PEMT_bytissue")

plot.selection.against.everything.volcano.facet(everything_object = GTEX_SET_HMTs_bytissue,
                                                selection_ensembl = NNMT_ensembl,
                                                selection_colours = "red",
                                                filename = "GTEX_HMTs_NNMT_bytissue")

plot.selection.against.everything.volcano.facet(everything_object = TCGA_SET_HMTs_bycancer,
                                                selection_ensembl = NNMT_ensembl,
                                                selection_colours = "red",
                                                filename = "TCGA_HMTs_NNMT_bytissue"
                                                )

TCGA_NNMT_bycancer_plot <- TCGA_NNMT_bycancer
names(TCGA_NNMT_bycancer_plot) <- str_remove(names(TCGA_NNMT_bycancer_plot), pattern = "TCGA-")

plot.selection.against.everything.volcano.facet(everything_object = TCGA_NNMT_bycancer_plot,
                                                selection_ensembl = act_rep_list,
                                                selection_colours = c("green", "magenta", "dodgerblue", "black"),
                                                filename = "TCGA_HMTs_NNMT_bytissue",
                                                fileoutput = "png",
                                                xaxislabel = substitute("Spearman's"~rho~"to"~italic(NNMT)))

plot.selection.against.everything.volcano.full.plots(everything_object = TCGA_NNMT_bycancer_plot,
                                                     selection_ensembl = act_rep_list,
                                                     selection_colours = c("green", "magenta", "dodgerblue", "black"),
                                                     filename_prefix = "TCGA_NNMT_HMTs",
                                                     fileoutput = "png")

TCGA_PEMT_bycancer_plot <- TCGA_PEMT_bycancer
names(TCGA_PEMT_bycancer_plot) <- str_remove(names(TCGA_PEMT_bycancer_plot), pattern = "TCGA-")

plot.selection.against.everything.volcano.facet(everything_object = TCGA_PEMT_bycancer_plot,
                                                selection_ensembl = act_rep_list,
                                                selection_colours = c("green", "magenta", "dodgerblue", "black"),
                                                filename = "TCGA_HMTs_PEMT_bytissue",
                                                fileoutput = "png",
                                                xaxislabel = substitute("Spearman's"~rho~"to"~italic(PEMT)))

plot.selection.against.everything.volcano.facet(everything_object = CCLE_NNMT_bycancer,
                                                selection_ensembl = act_rep_list,
                                                selection_colours = c("green", "magenta", "dodgerblue", "black"),
                                                filename = "CCLE_HMTs_NNMT_bytissue",
                                                fileoutput = "png",
                                                xaxislabel = substitute("Spearman's"~rho~"to"~italic(NNMT)),
                                                ylim = c(0, 12))

#### Supplementary File 1 / File S1 ####

CCLE_NNMT_bycancer_plot <- CCLE_NNMT_bycancer
names(CCLE_NNMT_bycancer_plot) <- str_replace_all(names(CCLE_NNMT_bycancer), "/", " & ")

plot.selection.against.everything.volcano.full.plots(everything_object = CCLE_NNMT_bycancer_plot,
                                                     selection_ensembl = act_rep_list,
                                                     selection_colours = c("green", "magenta", "dodgerblue", "black"),
                                                     filename_prefix = "CCLE_NNMT_HMTs",
                                                     fileoutput = "png")

#### Cancer type violins for example Fig 1G ####

make.cancer.type.violins <- function(cancertypedata,
                                     cancername){

thisone_NNMT <- TCGA_NNMT_bycancer[[cancername]]

thisone_NNMT[, "group"] <- "Other genes" 
thisone_NNMT[row.names(thisone_NNMT) %in% SET_HMTs$ensembl_gene_id, "group"] <- "HMTs"

thisone_NNMT[, "group"] <- factor(thisone_NNMT$group, levels = c("Other genes", "HMTs"))

thisone_NNMT_allHMTs <- thisone_NNMT["allHMTs", ]
thisone_NNMT_allHMTs[, "group"] <- "HMTs"

thisone_NNMT <- thisone_NNMT[row.names(thisone_NNMT) != "allHMTs", ]
thisone_NNMT <- thisone_NNMT[!str_detect(row.names(thisone_NNMT), NNMT_ensembl), ]

thisone_NNMT_HMTsubset <- thisone_NNMT[row.names(thisone_NNMT) %in% SET_HMTs$ensembl_gene_id, ]

for(i in 1:length(act_rep_list)){
  
  thisone_NNMT_HMTsubset[act_rep_list[[i]], "colour"] <- selection_colours[i]
  
}

thisone_NNMT_HMTsubset <- thisone_NNMT_HMTsubset[SET_HMTs$ensembl_gene_id, ]

pdf(paste0("graphics/", cancername,"_violin.pdf"),
    height = 2.5,
    width = 2.5)

print(ggplot(aes(x = group, y = allcorrs, fill = group), data = thisone_NNMT) + 
  geom_violin(trim = TRUE) + 
  scale_fill_manual(values = c("grey", "black")) +
  geom_point(data = thisone_NNMT_allHMTs, colour = "black", size = 2) + 
  geom_jitter(data = thisone_NNMT_HMTsubset, width = 0.05, colour = thisone_NNMT_HMTsubset$colour, size = 0.5) +
  geom_hline(yintercept = 0, 
             linetype = "dashed",
             colour = "grey55") +
  theme_classic() + 
  theme(axis.text.x = element_text(colour = "black", 
                                   size = 10),
        axis.text.y = element_text(colour = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10),
        legend.position = "none") +
  ylab("Spearman's"~rho~"to"~italic(NNMT)))

dev.off()

}

lapply(names(TCGA_MOR_list), function(x){
  
  make.cancer.type.violins(cancertypedata = TCGA_MOR_list[[x]],
                           cancername = x)
  
})

# how many PAAD samples?
sum(colnames(TCGA_MOR_list[["TCGA-PAAD"]]) %in% TCGA_LUT$SAMPID)

# respective ranks
which(row.names(TCGA_NNMT_bycancer[["TCGA-PAAD"]][order(TCGA_NNMT_bycancer$`TCGA-PAAD`$allcorrs), ]) == "allHMTs")
which(row.names(TCGA_SET_HMTs_bycancer[["TCGA-PAAD"]][order(TCGA_SET_HMTs_bycancer$`TCGA-PAAD`$allcorrs), ]) == NNMT_ensembl)

# PAAD plot data

PAAD_NNMT <- TCGA_NNMT_bycancer[["TCGA-PAAD"]]

PAAD_NNMT[, "group"] <- "Other genes" 
PAAD_NNMT[row.names(PAAD_NNMT) %in% SET_HMTs$ensembl_gene_id, "group"] <- "HMTs"

PAAD_NNMT[, "group"] <- factor(PAAD_NNMT$group, levels = c("Other genes", "HMTs"))

PAAD_NNMT_allHMTs <- PAAD_NNMT["allHMTs", ]
PAAD_NNMT_allHMTs[, "group"] <- "HMTs"

PAAD_NNMT <- PAAD_NNMT[row.names(PAAD_NNMT) != "allHMTs", ]
PAAD_NNMT <- PAAD_NNMT[!str_detect(row.names(PAAD_NNMT), NNMT_ensembl), ]

write.table(PAAD_NNMT,
            file = "plot_data/Fig 1/Fig_1G_plotdata.txt",
            row.names = TRUE,
            col.names = TRUE,
            sep = "\t")

#### GTEX violins for example - Fig 2E ####

# GTEX_MOR_list <- readRDS("output/GTEX_MOR_list.rds")

make.GTEX.type.violins <- function(tissuename,
                                   selection_colours){
  
  thisone_PEMT <- GTEX_PEMT_bytissue[[tissuename]]
  
  thisone_PEMT[, "group"] <- "Other genes" 
  thisone_PEMT[row.names(thisone_PEMT) %in% SET_HMTs$ensembl_gene_id, "group"] <- "HMTs"
  
  thisone_PEMT[, "group"] <- factor(thisone_PEMT$group, levels = c("Other genes", "HMTs"))
  
  thisone_PEMT_allHMTs <- thisone_PEMT["allHMTs", ]
  thisone_PEMT_allHMTs[, "group"] <- "HMTs"
  
  thisone_PEMT <- thisone_PEMT[row.names(thisone_PEMT) != "allHMTs", ]
  thisone_PEMT <- thisone_PEMT[!str_detect(row.names(thisone_PEMT), PEMT_ensembl), ]
  
  thisone_PEMT_HMTsubset <- thisone_PEMT[row.names(thisone_PEMT) %in% SET_HMTs$ensembl_gene_id, ]
  
  for(i in 1:length(act_rep_list)){
    
    thisone_PEMT_HMTsubset[act_rep_list[[i]], "colour"] <- selection_colours[i]
    
  }
  
  thisone_PEMT_HMTsubset <- thisone_PEMT_HMTsubset[SET_HMTs$ensembl_gene_id, ]
  
  pdf(paste0("graphics/", tissuename,"_violin.pdf"),
      height = 2.5,
      width = 2.5)
  
  print(ggplot(aes(x = group, y = allcorrs, fill = group), data = thisone_PEMT) + 
          geom_violin(trim = TRUE) + 
          scale_fill_manual(values = c("grey", "black")) +
          geom_point(data = thisone_PEMT_allHMTs, colour = "black", size = 2) + 
          geom_jitter(data = thisone_PEMT_HMTsubset, width = 0.05, colour = thisone_PEMT_HMTsubset$colour, size = 0.5) +
          geom_hline(yintercept = 0, 
                     linetype = "dashed",
                     colour = "grey55") +
          theme_classic() + 
          theme(axis.text.x = element_text(colour = "black", 
                                           size = 10),
                axis.text.y = element_text(colour = "black"),
                axis.title.x = element_blank(),
                axis.title.y = element_text(size = 10),
                legend.position = "none") +
          ylab("Spearman's"~rho~"to"~italic(PEMT)))
  
  dev.off()
  
}

lapply(names(GTEX_MOR_list), function(x){
  
  make.GTEX.type.violins(tissuename = x)
  
})

make.GTEX.type.violins(tissuename = "Esophagus - Gastroesophageal Junction",
                       selection_colours = c("green", "magenta", "dodgerblue", "black"))

# plot data to save

EsGasJ_PEMT <- GTEX_PEMT_bytissue[["Esophagus - Gastroesophageal Junction"]]

EsGasJ_PEMT[, "group"] <- "Other genes" 
EsGasJ_PEMT[row.names(EsGasJ_PEMT) %in% SET_HMTs$ensembl_gene_id, "group"] <- "HMTs"

EsGasJ_PEMT[, "group"] <- factor(EsGasJ_PEMT$group, levels = c("Other genes", "HMTs"))

write.table(EsGasJ_PEMT,
            file = "plot_data/Fig 2/Fig_2E_plotdata.txt",
            row.names = TRUE,
            col.names = TRUE,
            sep = "\t")

#### CCLE and TCGA pancancer volcanos - Fig 1D, Fig 1F ####

CCLE_NNMT_with_everything <- readRDS("output/CCLE_NNMT_correlate_with_everything.rds")
CCLE_NNMT_everything_forplot <- data.frame(cbind(CCLE_NNMT_with_everything$medians, CCLE_NNMT_with_everything$FDR))
colnames(CCLE_NNMT_everything_forplot) <- c("rho", "FDR")

write.table(CCLE_NNMT_everything_forplot,
            "plot_data/Fig 1/Fig_1E_plotdata.txt",
            row.names = TRUE,
            col.names = TRUE,
            sep = "\t")

plot.selection.against.everything.volcano(everything_object = CCLE_NNMT_everything_forplot,
                                          selection_ensembl = act_rep_list,
                                          selection_colours = selection_colours,
                                          filename = "CCLE_NNMT_vs_HMTs",
                                          xlim = c(-0.5, 0.5),
                                          ylim = c(0, 30),
                                          yaxisbreaks = seq(0, 30, 10),
                                          xaxislabel = substitute("Spearman's"~italic(rho)~"to"~italic(NNMT)),
                                          yaxislabel = substitute(-log[10]~"(false discovery rate)"),
                                          plotheight = 2.5, 
                                          plotwidth = 2.5,
                                          fileoutput = "png")

TCGA_NNMT_with_everything <- readRDS("output/TCGA_NNMT_correlate_with_everything.rds")
TCGA_NNMT_everything_forplot <- data.frame(cbind(TCGA_NNMT_with_everything$medians, TCGA_NNMT_with_everything$FDR))
colnames(TCGA_NNMT_everything_forplot) <- c("rho", "FDR")

write.table(TCGA_NNMT_everything_forplot,
            "plot_data/Fig 1/Fig_1F_plotdata.txt",
            row.names = TRUE,
            col.names = TRUE,
            sep = "\t")

plot.selection.against.everything.volcano(everything_object = TCGA_NNMT_everything_forplot,
                                          selection_ensembl = act_rep_list,
                                          selection_colours = c("green", "magenta", "dodgerblue", "black"),
                                          filename = "TCGA_NNMT_vs_HMTs",
                                          xlim = c(-0.5, 0.5),
                                          ylim = c(0, 75),
                                          yaxisbreaks = seq(0, 75, 15),
                                          xaxislabel = substitute("Spearman's"~italic(rho)~"to"~italic(NNMT)),
                                          yaxislabel = substitute(-log[10]~"(false discovery rate)"),
                                          plotheight = 2.5, 
                                          plotwidth = 2.5,
                                          fileoutput = "png")

TCGA_HMTs_with_everything <- readRDS("output/TCGA_allHMTs_correlate_with_everything.rds")
TCGA_HMTs_everything_forplot <- data.frame(cbind(TCGA_HMTs_with_everything$medians, TCGA_HMTs_with_everything$FDR))
colnames(TCGA_HMTs_everything_forplot) <- c("rho", "FDR")

plot.selection.against.everything.volcano(everything_object = TCGA_HMTs_everything_forplot,
                                          selection_ensembl = NNMT_ensembl,
                                          selection_colours = "black",
                                          filename = "TCGA_NNMT_vs_allHMTs",
                                          xlim = c(-0.5, 0.5),
                                          ylim = c(0, 75),
                                          xaxislabel = substitute("Spearman's"~italic(rho)~"to"~italic(HMTs)),
                                          yaxislabel = substitute(-log[10]~"(FDR)"),
                                          plotheight = 2, 
                                          plotwidth = 2)

TCGA_PEMT_with_everything <- readRDS("output/TCGA_PEMT_correlate_with_everything.rds")
TCGA_PEMT_everything_forplot <- data.frame(cbind(TCGA_PEMT_with_everything$medians, TCGA_PEMT_with_everything$FDR))
colnames(TCGA_PEMT_everything_forplot) <- c("rho", "FDR")

plot.selection.against.everything.volcano(everything_object = TCGA_PEMT_everything_forplot,
                                          selection_ensembl = act_rep_list,
                                          selection_colours = c("green", "magenta", "dodgerblue", "red"),
                                          filename = "TCGA_PEMT_vs_HMTs",
                                          xlim = c(-0.5, 0.5),
                                          ylim = c(0, 50),
                                          xaxislabel = substitute("Spearman's"~italic(rho)),
                                          yaxislabel = substitute(-log[10]~"(FDR)"),
                                          plotheight = 2, 
                                          plotwidth = 2)

GTEX_NNMT_with_everything <- readRDS("output/GTEX_NNMT_correlate_with_everything.rds")
GTEX_NNMT_everything_forplot <- data.frame(cbind(GTEX_NNMT_with_everything$medians, GTEX_NNMT_with_everything$FDR))
colnames(GTEX_NNMT_everything_forplot) <- c("rho", "FDR")

write.table(GTEX_NNMT_everything_forplot,
            file = "plot_data/Fig 2/Fig_2B_plotdata.txt",
            row.names = TRUE,
            col.names = TRUE,
            sep = "\t")

GTEX_individualHMT_NNMTranks <- match(SET_HMTs$ensembl_gene_id, names(GTEX_NNMT_with_everything$medians[order(GTEX_NNMT_with_everything$medians)]))
names(GTEX_individualHMT_NNMTranks) <- SET_HMTs$hgnc_symbol

plot.selection.against.everything.volcano(everything_object = GTEX_NNMT_everything_forplot,
                                          selection_ensembl = act_rep_list,
                                          selection_colours = c("green", "magenta", "dodgerblue", "black"),
                                          filename = "GTEX_NNMT_vs_HMTs",
                                          xlim = c(-0.6, 0.6),
                                          ylim = c(0, 400),
                                          xaxislabel = substitute("Spearman's"~italic(rho)),
                                          yaxislabel = substitute(-log[10]~"(FDR)"),
                                          plotheight = 2, 
                                          plotwidth = 2)

GTEX_PEMT_with_everything <- readRDS("output/GTEX_PEMT_correlate_with_everything.rds")
GTEX_PEMT_everything_forplot <- data.frame(cbind(GTEX_PEMT_with_everything$medians, GTEX_PEMT_with_everything$FDR))
colnames(GTEX_PEMT_everything_forplot) <- c("rho", "FDR")

plot.selection.against.everything.volcano(everything_object = GTEX_PEMT_everything_forplot,
                                          selection_ensembl = act_rep_list,
                                          selection_colours = c("green", "magenta", "dodgerblue", "black"),
                                          filename = "GTEX_PEMT_vs_HMTs",
                                          xlim = c(-0.6, 0.6),
                                          ylim = c(-5, 350),
                                          xaxislabel = substitute("Spearman's"~italic(rho)~"to"~italic(PEMT)),
                                          yaxislabel = substitute(-log[10]~"(false discovery rate)"),
                                          yaxisbreaks = seq(0, 350, 50), 
                                          plotheight = 2.5, 
                                          plotwidth = 2.5,
                                          fileoutput = "png")

GTEX_PEMT_HMTranks <- match(SET_HMTs$ensembl_gene_id, row.names(GTEX_PEMT_everything_forplot[order(GTEX_PEMT_everything_forplot$rho), ]))
names(GTEX_PEMT_HMTranks) <- SET_HMTs$hgnc_symbol

GTEX_NNMT_with_everything <- readRDS("output/GTEX_NNMT_correlate_with_everything.rds")
GTEX_NNMT_everything_forplot <- data.frame(cbind(GTEX_NNMT_with_everything$medians, GTEX_NNMT_with_everything$FDR))
colnames(GTEX_NNMT_everything_forplot) <- c("rho", "FDR")

plot.selection.against.everything.volcano(everything_object = GTEX_NNMT_everything_forplot,
                                          selection_ensembl = act_rep_list,
                                          selection_colours = c("green", "magenta", "dodgerblue", "black"),
                                          filename = "GTEX_NNMT_vs_HMTs",
                                          xlim = c(-0.6, 0.6),
                                          ylim = c(-5, 350),
                                          xaxislabel = substitute("Spearman's"~italic(rho)~"to"~italic(NNMT)),
                                          yaxislabel = substitute(-log[10]~"(false discovery rate)"),
                                          yaxisbreaks = seq(0, 350, 50),
                                          plotheight = 2.5, 
                                          plotwidth = 2.5,
                                          fileoutput = "png")

GTEX_individualHMT_PEMTranks <- match(SET_HMTs$ensembl_gene_id, names(GTEX_PEMT_with_everything$medians[order(GTEX_PEMT_with_everything$medians)]))
names(GTEX_individualHMT_PEMTranks) <- SET_HMTs$hgnc_symbol

#### SPECIFIC LYSINE RESIDUE SETS ####

# define random samples for GTEX for strong tissues
define.random.samples.for.iterations(database = "GTEX",
                                     supplied_datalist = GTEX_MOR_list[GTEX_PEMT_HMT_strongtissues],
                                     iterations = 100,
                                     save = TRUE,
                                     filename = "chosensamples_GTEX_strongtissues_100iterations")

allHMTs <- SET_HMTs_annotated$ensembl_gene_id

H3K4mts <- c(SET_HMTs_annotated[str_detect(SET_HMTs_annotated$Mark, "H3K4"), "ensembl_gene_id"])

H3K9mts <- c(SET_HMTs_annotated[str_detect(SET_HMTs_annotated$Mark, "H3K9"), "ensembl_gene_id"])

H3K27mts <- c(SET_HMTs_annotated[str_detect(SET_HMTs_annotated$Mark, "H3K27"), "ensembl_gene_id"])

H3K36mts <- c(SET_HMTs_annotated[str_detect(SET_HMTs_annotated$Mark, "H3K36"), "ensembl_gene_id"])

H4K20mts <- c(SET_HMTs_annotated[str_detect(SET_HMTs_annotated$Mark, "H4K20"), "ensembl_gene_id"])

geneset.with.everything(database = c("GTEX"),
                        geneset_ensembl = SET_HMTs_annotated[str_detect(SET_HMTs_annotated$Mark, "H3K4"), "ensembl_gene_id"],
                        basename = "H3K4mts",
                        save = TRUE,
                        iterations = 100)

geneset.with.everything(database = c("TCGA"),
                        geneset_ensembl = SET_HMTs_annotated[str_detect(SET_HMTs_annotated$Mark, "H3K4"), "ensembl_gene_id"],
                        basename = "H3K4mts",
                        save = TRUE,
                        iterations = 100)

geneset.with.everything(database = c("GTEX"),
                        geneset_ensembl = SET_HMTs_annotated[str_detect(SET_HMTs_annotated$Mark, "H3K9"), "ensembl_gene_id"],
                        basename = "H3K9mts",
                        save = TRUE,
                        iterations = 100)

geneset.with.everything(database = c("TCGA"),
                        geneset_ensembl = SET_HMTs_annotated[str_detect(SET_HMTs_annotated$Mark, "H3K9"), "ensembl_gene_id"],
                        basename = "H3K9mts",
                        save = TRUE,
                        iterations = 100)

geneset.with.everything(database = c("GTEX"),
                        geneset_ensembl = SET_HMTs_annotated[str_detect(SET_HMTs_annotated$Mark, "H3K27"), "ensembl_gene_id"],
                        basename = "H3K27mts",
                        save = TRUE,
                        iterations = 100)

geneset.with.everything(database = c("TCGA"),
                        geneset_ensembl = SET_HMTs_annotated[str_detect(SET_HMTs_annotated$Mark, "H3K27"), "ensembl_gene_id"],
                        basename = "H3K27mts",
                        save = TRUE,
                        iterations = 100)

geneset.with.everything(database = c("GTEX"),
                        geneset_ensembl = SET_HMTs_annotated[str_detect(SET_HMTs_annotated$Mark, "H3K36"), "ensembl_gene_id"],
                        basename = "H3K36mts",
                        save = TRUE,
                        iterations = 100)

geneset.with.everything(database = c("TCGA"),
                        geneset_ensembl = SET_HMTs_annotated[str_detect(SET_HMTs_annotated$Mark, "H3K36"), "ensembl_gene_id"],
                        basename = "H3K36mts",
                        save = TRUE,
                        iterations = 100)

## TCGA ##

geneset.with.everything(database = c("GTEX"),
                        supplied_datalist = GTEX_MOR_list[GTEX_PEMT_HMT_strongtissues],
                        geneset_ensembl = SET_HMTs_annotated[str_detect(SET_HMTs_annotated$Mark, "H3K4"), "ensembl_gene_id"],
                        basename = "H3K4mts_strongtissues",
                        save = TRUE,
                        iterations = 100,
                        chosensamples = chosensamples_GTEX_strongtissues_100iterations)

geneset.with.everything(database = c("GTEX"),
                        supplied_datalist = GTEX_MOR_list[GTEX_PEMT_HMT_strongtissues],
                        geneset_ensembl = SET_HMTs_annotated[str_detect(SET_HMTs_annotated$Mark, "H3K9"), "ensembl_gene_id"],
                        basename = "H3K9mts_strongtissues",
                        save = TRUE,
                        iterations = 100,
                        chosensamples = chosensamples_GTEX_strongtissues_100iterations)

geneset.with.everything(database = c("GTEX"),
                        supplied_datalist = GTEX_MOR_list[GTEX_PEMT_HMT_strongtissues],
                        geneset_ensembl = SET_HMTs_annotated[str_detect(SET_HMTs_annotated$Mark, "H3K27"), "ensembl_gene_id"],
                        basename = "H3K27mts_strongtissues",
                        save = TRUE,
                        iterations = 100,
                        chosensamples = chosensamples_GTEX_strongtissues_100iterations)

geneset.with.everything(database = c("GTEX"),
                        supplied_datalist = GTEX_MOR_list[GTEX_PEMT_HMT_strongtissues],
                        geneset_ensembl = SET_HMTs_annotated[str_detect(SET_HMTs_annotated$Mark, "H3K36"), "ensembl_gene_id"],
                        basename = "H3K36mts_strongtissues",
                        save = TRUE,
                        iterations = 100,
                        chosensamples = chosensamples_GTEX_strongtissues_100iterations)

geneset.with.everything(database = c("GTEX"),
                        supplied_datalist = GTEX_MOR_list[GTEX_PEMT_HMT_strongtissues],
                        geneset_ensembl = SET_HMTs_annotated[str_detect(SET_HMTs_annotated$Mark, "H4K20"), "ensembl_gene_id"],
                        basename = "H4K20mts_strongtissues",
                        save = TRUE,
                        iterations = 100,
                        chosensamples = chosensamples_GTEX_strongtissues_100iterations)

# let's do with TCGA_PEMT_strongtissues

TCGA_PEMT_HMT_strongcancers <- row.names(TCGA_PEMT_cancer_reciprocal_scores_df[TCGA_PEMT_cancer_reciprocal_scores_df$HMTs_in_PEMT_percent < 2.5 & TCGA_PEMT_cancer_reciprocal_scores_df$PEMT_in_HMTs_percent < 2.5, ])

# define random samples for TCGA for strong cancers
define.random.samples.for.iterations(database = "TCGA",
                                     supplied_datalist = TCGA_MOR_list[TCGA_PEMT_HMT_strongcancers],
                                     iterations = 100,
                                     save = TRUE,
                                     filename = "chosensamples_TCGA_strongcancers_100iterations")

# just MTs

geneset.with.everything(database = c("TCGA"),
                        supplied_datalist = TCGA_MOR_list[TCGA_PEMT_HMT_strongcancers],
                        geneset_ensembl = SET_HMTs_annotated[str_detect(SET_HMTs_annotated$Mark, "H3K4"), "ensembl_gene_id"],
                        basename = "H3K4mts_strongcancers",
                        save = TRUE,
                        iterations = 100,
                        chosensamples = chosensamples_TCGA_strongcancers_100iterations)

geneset.with.everything(database = c("TCGA"),
                        supplied_datalist = TCGA_MOR_list[TCGA_PEMT_HMT_strongcancers],
                        geneset_ensembl = SET_HMTs_annotated[str_detect(SET_HMTs_annotated$Mark, "H3K9"), "ensembl_gene_id"],
                        basename = "H3K9mts_strongcancers",
                        save = TRUE,
                        iterations = 100,
                        chosensamples = chosensamples_TCGA_strongcancers_100iterations)

geneset.with.everything(database = c("TCGA"),
                        supplied_datalist = TCGA_MOR_list[TCGA_PEMT_HMT_strongcancers],
                        geneset_ensembl = SET_HMTs_annotated[str_detect(SET_HMTs_annotated$Mark, "H3K27"), "ensembl_gene_id"],
                        basename = "H3K27mts_strongcancers",
                        save = TRUE,
                        iterations = 100,
                        chosensamples = chosensamples_TCGA_strongcancers_100iterations)

geneset.with.everything(database = c("TCGA"),
                        supplied_datalist = TCGA_MOR_list[TCGA_PEMT_HMT_strongcancers],
                        geneset_ensembl = SET_HMTs_annotated[str_detect(SET_HMTs_annotated$Mark, "H3K36"), "ensembl_gene_id"],
                        basename = "H3K36mts_strongcancers",
                        save = TRUE,
                        iterations = 100,
                        chosensamples = chosensamples_TCGA_strongcancers_100iterations)

geneset.with.everything(database = c("TCGA"),
                        supplied_datalist = TCGA_MOR_list[TCGA_PEMT_HMT_strongcancers],
                        geneset_ensembl = SET_HMTs_annotated[str_detect(SET_HMTs_annotated$Mark, "H4K20"), "ensembl_gene_id"],
                        basename = "H4K20mts_strongcancers",
                        save = TRUE,
                        iterations = 100,
                        chosensamples = chosensamples_TCGA_strongcancers_100iterations)

geneset.with.everything(database = c("TCGA"),
                        supplied_datalist = TCGA_MOR_list[TCGA_PEMT_HMT_strongcancers],
                        geneset_ensembl = SET_HMTs_annotated$ensembl_gene_id,
                        basename = "allHMTs_strongcancers",
                        save = TRUE,
                        iterations = 100,
                        chosensamples = chosensamples_TCGA_strongcancers_100iterations)

#### RELATIVE RECIPROCAL SCORES FOR LYSINE RESIDUES - for Fig 2I, Fig 2J, Fig 2K ####

TCGA_NNMT_all_mark_sets_list <- lapply(c("H3K4mts",
                                         "H3K9mts",
                                         "H3K27mts",
                                         "H3K36mts",
                                         "H4K20mts",
                                         "allHMTs"), function(x){
                        
  temp_everything <- readRDS(paste0("output/TCGA_", x, "_correlate_with_everything.rds"))
  
  if(x != "allHMTs"){
  
  NNMT_added <- add.geneset.to.everything(database = "TCGA",
                                          geneset_ensembl = get(x),
                                          genesetname = x,
                                          everything_object = TCGA_NNMT_with_everything)

  } else {
    
    NNMT_added <- TCGA_NNMT_with_everything
    
  }

  tempscore <- calculate.relative.reciprocal.scores(everything_object1 = NNMT_added$medians,
                                                    everything_object2 = temp_everything$medians,
                                                    gene1_name = "NNMT")
  
  output_df <- as.data.frame(tempscore)
  row.names(output_df) <- c("NNMT_rank_in_geneset_everything", "Geneset_rank_in_NNMT_everything", "reciprocal_score")
  colnames(output_df) <- x
  
  return(output_df)
  
})

TCGA_NNMT_all_mark_sets_df <- do.call(cbind, TCGA_NNMT_all_mark_sets_list)
TCGA_NNMT_all_mark_sets_df <- data.frame(t(TCGA_NNMT_all_mark_sets_df))

saveRDS(TCGA_NNMT_all_mark_sets_df, "output/TCGA_NNMT_all_mark_sets_df.rds")
TCGA_NNMT_all_mark_sets_df <- readRDS("output/TCGA_NNMT_all_mark_sets_df.rds")

TCGA_NNMT_all_mark_sets_df[, "set"] <- row.names(TCGA_NNMT_all_mark_sets_df)

TCGA_NNMT_all_mark_sets_df[, "set_in_NNMT_percent"] <- (TCGA_NNMT_all_mark_sets_df$Geneset_rank_in_NNMT_everything * 100)/60489
TCGA_NNMT_all_mark_sets_df[, "NNMT_in_set_percent"] <- (TCGA_NNMT_all_mark_sets_df$NNMT_rank_in_geneset_everything * 100)/60489

TCGA_NNMT_all_mark_sets_df[, "set"] <- str_replace(TCGA_NNMT_all_mark_sets_df$set, "mts", " MTs")
TCGA_NNMT_all_mark_sets_df[, "set"] <- str_replace(TCGA_NNMT_all_mark_sets_df$set, "set", " MTs + deMTs")
TCGA_NNMT_all_mark_sets_df[, "set"] <- str_replace(TCGA_NNMT_all_mark_sets_df$set, "allHMTs", "all HMTs")

TCGA_NNMT_all_mark_MTs_df <- TCGA_NNMT_all_mark_sets_df[!str_detect(TCGA_NNMT_all_mark_sets_df$set, "deMTs"), ]

write.table(TCGA_NNMT_all_mark_MTs_df,
            file = "plot_data/Fig 2/Fig_2I_plotdata.txt",
            row.names = TRUE,
            col.names = TRUE,
            sep = "\t")

pdf("graphics/TCGA-MTsets-reciprocal_scores.pdf",
    width = 2,
    height = 2)

ggplot(aes(x = NNMT_in_set_percent, y = set_in_NNMT_percent, label = set), data = TCGA_NNMT_all_mark_sets_df) + 
  geom_rect(aes(xmax = 2.5,
                xmin = 0,
                ymax = 2.5,
                ymin = 0),
            fill = "white",
            linetype = "dashed",
            color = "grey") +
  geom_point(aes(size = 10/log10(reciprocal_score)),
             colour = "red",
             alpha = 0.7) + 
  scale_x_log10(breaks = c(0.01, 0.1, 1, 10, 100),
                expand = c(0, 0),
                limits = c(0.001, 100),
                labels = label_number(drop0trailing = TRUE)) +
  scale_y_log10(breaks = c(0.01, 0.1, 1, 10, 100),
                expand = c(0, 0),
                limits = c(0.001, 100),
                labels = label_number(drop0trailing = TRUE)) +
  geom_label_repel(size = 1.5,
                   max.overlaps = 40,
                   label.size = 0.1,
                   box.padding = 0.1,
                   label.padding = 0.1) + 
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(colour = "black",
                                   size = 8),
        axis.text.y = element_text(colour = "black",
                                   size = 8),
        axis.title.x = element_text(size = 8,
                                    hjust = 1),
        axis.title.y = element_text(size = 8,
                                    hjust = 1),
        plot.margin = margin(t = 5,
                             r = 8,
                             b = 3,
                             l = 3)) + 
  xlab(substitute(italic(NNMT)~"%"^"ile"~"in geneset correlations")) +
  ylab(substitute("Geneset"~"%"^"ile"~"in"~italic(NNMT)~"correlations"))

dev.off()

pdf("graphics/TCGA-MTsalone-reciprocal_scores.pdf",
    width = 2.2,
    height = 2.2)

ggplot(aes(x = NNMT_in_set_percent, y = set_in_NNMT_percent, label = set), data = TCGA_NNMT_all_mark_MTs_df) + 
  geom_rect(aes(xmax = 2.5,
                xmin = 0,
                ymax = 2.5,
                ymin = 0),
            fill = "white",
            linetype = "dashed",
            color = "grey") +
  geom_point(aes(size = 10/log10(reciprocal_score)),
             colour = "red",
             alpha = 0.7) + 
  scale_x_log10(breaks = c(0.01, 0.1, 1, 10, 100),
                expand = c(0, 0),
                limits = c(0.001, 100),
                labels = label_number(drop0trailing = TRUE)) +
  scale_y_log10(breaks = c(0.01, 0.1, 1, 10, 100),
                expand = c(0, 0),
                limits = c(0.001, 100),
                labels = label_number(drop0trailing = TRUE)) +
  geom_label_repel(size = 1.5,
                   max.overlaps = 40,
                   label.size = 0.1,
                   box.padding = 0.1,
                   label.padding = 0.1) + 
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(colour = "black",
                                   size = 8),
        axis.text.y = element_text(colour = "black",
                                   size = 8),
        axis.title.x = element_text(size = 8,
                                    hjust = 1),
        axis.title.y = element_text(size = 8,
                                    hjust = 1),
        plot.margin = margin(t = 5,
                             r = 8,
                             b = 3,
                             l = 3)) + 
  xlab(substitute(italic(NNMT)~"%"^"ile"~"in geneset correlations")) +
  ylab(substitute("Geneset"~"%"^"ile"~"in"~italic(NNMT)~"correlations"))

dev.off()

# now GTEX with PEMT - Fig 2C

GTEX_strongtissues_PEMT_with_everything <- geneset.with.everything(database = "GTEX",
                                                                   supplied_datalist = GTEX_MOR_list[GTEX_PEMT_HMT_strongtissues],
                                                                   geneset_ensembl = PEMT_ensembl,
                                                                   basename = "PEMT",
                                                                   save = TRUE,
                                                                   iterations = 100,
                                                                   chosensamples = chosensamples_GTEX_strongtissues_100iterations,
                                                                   filename = "GTEXstrongtissues_PEMT_with_everything")

GTEX_strongtissues_PEMT_with_everything <- readRDS("output/GTEX_strongtissues_PEMT_with_everything.rds")

# plot out volcano for strong tissues

GTEX_strongtissues_PEMT_forvolcano <- add.geneset.to.everything(database = "GTEX",
                                                                geneset_ensembl = SET_HMTs_annotated$ensembl_gene_id,
                                                                genesetname = "allHMTs",
                                                                everything_object = GTEX_strongtissues_PEMT_with_everything)

GTEX_strongtissues_PEMT_with_everything_forplot <- data.frame(cbind(GTEX_strongtissues_PEMT_forvolcano$medians, GTEX_strongtissues_PEMT_forvolcano$FDR))
colnames(GTEX_strongtissues_PEMT_with_everything_forplot) <- c("rho", "FDR")
GTEX_strongtissues_PEMT_with_everything_forplot[order(GTEX_strongtissues_PEMT_with_everything_forplot$rho), ][1:10, "rho"]

write.table(GTEX_strongtissues_PEMT_with_everything_forplot,
            file = "plot_data/Fig 2/Fig_2D_plotdata.txt",
            row.names = TRUE,
            col.names = TRUE,
            sep = "\t")

plot.selection.against.everything.volcano(everything_object = GTEX_strongtissues_PEMT_with_everything_forplot,
                                          selection_ensembl = act_rep_list,
                                          selection_colours = c("green", "magenta", "dodgerblue", "black"),
                                          filename = "GTEX_PEMT_vs_HMTs_strongtissues",
                                          xlim = c(-0.75, 0.75),
                                          ylim = c(0, 200),
                                          yaxisbreaks = seq(0, 200, 50),
                                          xaxisbreaks = seq(-0.75, 0.75, 0.25),
                                          xaxislabel = substitute("Spearman's"~italic(rho)~"to"~italic(PEMT)),
                                          yaxislabel = substitute(-log[10]~"(false discovery rate)"),
                                          plotheight = 2.5, 
                                          plotwidth = 2.5,
                                          fileoutput = "png")

GTEX_strongtissues_PEMT_all_mark_sets_list <- lapply(c("H3K4set_strongtissues",
                                                       "H3K4mts_strongtissues",
                                         "H3K9set_strongtissues",
                                         "H3K9mts_strongtissues",
                                         "H3K27set_strongtissues",
                                         "H3K27mts_strongtissues",
                                         "H3K36set_strongtissues",
                                         "H3K36mts_strongtissues",
                                         "H4K20set_strongtissues",
                                         "H4K20mts_strongtissues",
                                         "allHMTs_strongtissues"
                                         ), function(x){

                                           temp_everything <- readRDS(paste0("output/GTEX_", x, "_correlate_with_everything.rds"))
                                             
                                             PEMT_added <- add.geneset.to.everything(database = "GTEX",
                                                                                     geneset_ensembl = get(str_remove(x, "_strongtissues")),
                                                                                     genesetname = x,
                                                                                     everything_object = GTEX_strongtissues_PEMT_with_everything)
                                           
                                           tempscore <- calculate.relative.reciprocal.scores(everything_object1 = PEMT_added$medians,
                                                                                             everything_object2 = temp_everything$medians,
                                                                                             gene1_name = "PEMT")
                                           
                                           output_df <- as.data.frame(tempscore)
                                           row.names(output_df) <- c("PEMT_rank_in_geneset_everything", "Geneset_rank_in_PEMT_everything", "reciprocal_score")
                                           colnames(output_df) <- str_remove(x, "_strongtissues")
                                           
                                           return(output_df)
                                           
                                         })

GTEX_strongtissues_PEMT_all_mark_sets_df <- do.call(cbind, GTEX_strongtissues_PEMT_all_mark_sets_list)
GTEX_strongtissues_PEMT_all_mark_sets_df <- data.frame(t(GTEX_strongtissues_PEMT_all_mark_sets_df))

saveRDS(GTEX_strongtissues_PEMT_all_mark_sets_df, "output/GTEX_strongtissues_PEMT_all_mark_sets_df.rds")
# GTEX_strongtissues_PEMT_all_mark_sets_df <- readRDS("output/GTEX_strongtissues_PEMT_all_mark_sets_df.rds")

GTEX_strongtissues_PEMT_all_mark_sets_df[, "set"] <- row.names(GTEX_strongtissues_PEMT_all_mark_sets_df)
GTEX_strongtissues_PEMT_all_mark_sets_df[, "set"] <- str_replace(GTEX_strongtissues_PEMT_all_mark_sets_df$set, "mts", " MTs")
GTEX_strongtissues_PEMT_all_mark_sets_df[, "set"] <- str_replace(GTEX_strongtissues_PEMT_all_mark_sets_df$set, "set", " MTs + deMTs")
GTEX_strongtissues_PEMT_all_mark_sets_df[, "set"] <- str_replace(GTEX_strongtissues_PEMT_all_mark_sets_df$set, "allHMTs", "all HMTs")

GTEX_strongtissues_PEMT_all_mark_sets_df[, "set_in_PEMT_percent"] <- (GTEX_strongtissues_PEMT_all_mark_sets_df$Geneset_rank_in_PEMT_everything * 100)/60489
GTEX_strongtissues_PEMT_all_mark_sets_df[, "PEMT_in_set_percent"] <- (GTEX_strongtissues_PEMT_all_mark_sets_df$PEMT_rank_in_geneset_everything * 100)/60489

GTEX_strongtissues_PEMT_all_mark_MTs_df <- GTEX_strongtissues_PEMT_all_mark_sets_df[!str_detect(GTEX_strongtissues_PEMT_all_mark_sets_df$set, "deMTs"), ]

write.table(GTEX_strongtissues_PEMT_all_mark_MTs_df,
            file = "plot_data/Fig 2/Fig_2J_plotdata.txt",
            row.names = TRUE,
            col.names = TRUE,
            sep = "\t")

pdf("graphics/GTEX-strongtissues-MTsets-reciprocal_scores.pdf",
    width = 2,
    height = 2)

ggplot(aes(x = PEMT_in_set_percent, y = set_in_PEMT_percent, label = set), data = GTEX_strongtissues_PEMT_all_mark_sets_df) + 
  geom_rect(aes(xmax = 2.5,
                xmin = 0,
                ymax = 2.5,
                ymin = 0),
            fill = "white",
            linetype = "dashed",
            color = "grey") +
  geom_point(aes(size = 10/log10(reciprocal_score)),
             colour = "red",
             alpha = 0.7) + 
  scale_x_log10(breaks = c(0.01, 0.1, 1, 10, 100),
                expand = c(0, 0),
                limits = c(0.001, 100),
                labels = label_number(drop0trailing = TRUE)) +
  scale_y_log10(breaks = c(0.01, 0.1, 1, 10, 100),
                expand = c(0, 0),
                limits = c(0.001, 100),
                labels = label_number(drop0trailing = TRUE)) +
  geom_label_repel(size = 1.5,
                   max.overlaps = 40,
                   label.size = 0.1,
                   box.padding = 0.1,
                   label.padding = 0.1) + 
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(colour = "black",
                                   size = 8),
        axis.text.y = element_text(colour = "black",
                                   size = 8),
        axis.title.x = element_text(size = 8,
                                    hjust = 1),
        axis.title.y = element_text(size = 8,
                                    hjust = 1),
        plot.margin = margin(t = 5,
                             r = 8,
                             b = 3,
                             l = 3)) + 
  xlab(substitute(italic(PEMT)~"%"^"ile"~"in geneset correlations")) +
  ylab(substitute("Geneset"~"%"^"ile"~"in"~italic(PEMT)~"correlations"))

dev.off()

pdf("graphics/GTEX-strongtissues-MTsalone-reciprocal_scores.pdf",
    width = 2.2,
    height = 2.2)

ggplot(aes(x = PEMT_in_set_percent, y = set_in_PEMT_percent, label = set), data = GTEX_strongtissues_PEMT_all_mark_MTs_df) + 
  geom_rect(aes(xmax = 2.5,
                xmin = 0,
                ymax = 2.5,
                ymin = 0),
            fill = "white",
            linetype = "dashed",
            color = "grey") +
  geom_point(aes(size = 10/log10(reciprocal_score)),
             colour = "red",
             alpha = 0.7) + 
  scale_x_log10(breaks = c(0.01, 0.1, 1, 10, 100),
                expand = c(0, 0),
                limits = c(0.001, 100),
                labels = label_number(drop0trailing = TRUE)) +
  scale_y_log10(breaks = c(0.01, 0.1, 1, 10, 100),
                expand = c(0, 0),
                limits = c(0.001, 100),
                labels = label_number(drop0trailing = TRUE)) +
  geom_label_repel(size = 1.5,
                   max.overlaps = 40,
                   label.size = 0.1,
                   box.padding = 0.1,
                   label.padding = 0.1) + 
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(colour = "black",
                                   size = 8),
        axis.text.y = element_text(colour = "black",
                                   size = 8),
        axis.title.x = element_text(size = 8,
                                    hjust = 1),
        axis.title.y = element_text(size = 8,
                                    hjust = 1),
        plot.margin = margin(t = 5,
                             r = 8,
                             b = 3,
                             l = 3)) + 
  xlab(substitute(italic(PEMT)~"%"^"ile"~"in geneset correlations")) +
  ylab(substitute("Geneset"~"%"^"ile"~"in"~italic(PEMT)~"correlations"))

dev.off()

# TCGA with strong cancers for PEMT 

TCGA_strongtissues_PEMT_with_everything <- anything.with.everything(database = "TCGA",
                                                                    supplied_datalist = TCGA_MOR_list[TCGA_PEMT_HMT_strongcancers],
                                                                    ensembl_gene = PEMT_ensembl,
                                                                    basename = "PEMT",
                                                                    save = TRUE,
                                                                    iterations = 100,
                                                                    chosensamples = chosensamples_TCGA_strongcancers_100iterations,
                                                                    filename = "TCGA_strongcancers_PEMT_with_everything")

# TCGA_strongcancers_PEMT_with_everything <- readRDS("output/TCGA_strongcancers_PEMT_with_everything.rds")

TCGA_strongcancers_PEMT_all_mark_sets_list <- lapply(c("H3K4mts_strongcancers",
                                                       "H3K9mts_strongcancers",
                                                       "H3K27mts_strongcancers",
                                                       "H3K36mts_strongcancers",
                                                       "H4K20mts_strongcancers",
                                                       "allHMTs_strongcancers"), function(x){

                                                                                                    
  temp_everything <- readRDS(paste0("output/TCGA_", x, "_correlate_with_everything.rds"))

  PEMT_added <- add.geneset.to.everything(database = "TCGA",
                                          geneset_ensembl = get(str_remove(x, "_strongcancers")),
                                          genesetname = x,
                                          everything_object = TCGA_strongcancers_PEMT_with_everything)

    tempscore <- calculate.relative.reciprocal.scores(everything_object1 = PEMT_added$medians,
                                                    everything_object2 = temp_everything$medians,
                                                    gene1_name = "PEMT")
  
  output_df <- as.data.frame(tempscore)
  row.names(output_df) <- c("PEMT_rank_in_geneset_everything", "Geneset_rank_in_PEMT_everything", "reciprocal_score")
  colnames(output_df) <- str_remove(x, "_strongcancers")
  
  return(output_df)
  
})

TCGA_strongcancers_PEMT_all_mark_sets_df <- do.call(cbind, TCGA_strongcancers_PEMT_all_mark_sets_list)
TCGA_strongcancers_PEMT_all_mark_sets_df <- data.frame(t(TCGA_strongcancers_PEMT_all_mark_sets_df))

saveRDS(TCGA_strongcancers_PEMT_all_mark_sets_df, "output/TCGA_strongcancers_PEMT_all_mark_sets_df.rds")
TCGA_strongcancers_PEMT_all_mark_sets_df <- readRDS("output/TCGA_strongcancers_PEMT_all_mark_sets_df.rds")

TCGA_strongcancers_PEMT_all_mark_sets_df[, "set"] <- row.names(TCGA_strongcancers_PEMT_all_mark_sets_df)
TCGA_strongcancers_PEMT_all_mark_sets_df[, "set"] <- str_replace(TCGA_strongcancers_PEMT_all_mark_sets_df$set, "mts", " MTs")
TCGA_strongcancers_PEMT_all_mark_sets_df[, "set"] <- str_replace(TCGA_strongcancers_PEMT_all_mark_sets_df$set, "allHMTs", "all HMTs")

TCGA_strongcancers_PEMT_all_mark_sets_df[, "set_in_PEMT_percent"] <- (TCGA_strongcancers_PEMT_all_mark_sets_df$Geneset_rank_in_PEMT_everything * 100)/60489
TCGA_strongcancers_PEMT_all_mark_sets_df[, "PEMT_in_set_percent"] <- (TCGA_strongcancers_PEMT_all_mark_sets_df$PEMT_rank_in_geneset_everything * 100)/60489

write.table(TCGA_strongcancers_PEMT_all_mark_sets_df,
            file = "plot_data/Fig 2/Fig_2K_plotdata.txt",
            row.names = TRUE,
            col.names = TRUE,
            sep = "\t")

pdf("graphics/TCGA-strongcancers-MTsets-reciprocal_scores.pdf",
    width = 2.2,
    height = 2.2)

ggplot(aes(x = PEMT_in_set_percent, y = set_in_PEMT_percent, label = set), data = TCGA_strongcancers_PEMT_all_mark_sets_df) + 
  geom_rect(aes(xmax = 2.5,
                xmin = 0,
                ymax = 2.5,
                ymin = 0),
            fill = "white",
            linetype = "dashed",
            color = "grey") +
  geom_point(aes(size = 10/log10(reciprocal_score)),
             colour = "red",
             alpha = 0.7) + 
  scale_x_log10(breaks = c(0.01, 0.1, 1, 10, 100),
                expand = c(0, 0),
                limits = c(0.001, 100),
                labels = label_number(drop0trailing = TRUE)) +
  scale_y_log10(breaks = c(0.01, 0.1, 1, 10, 100),
                expand = c(0, 0),
                limits = c(0.001, 100),
                labels = label_number(drop0trailing = TRUE)) +
  geom_label_repel(size = 1.5,
                   max.overlaps = 40,
                   label.size = 0.1,
                   box.padding = 0.1,
                   label.padding = 0.1) + 
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(colour = "black",
                                   size = 8),
        axis.text.y = element_text(colour = "black",
                                   size = 8),
        axis.title.x = element_text(size = 8,
                                    hjust = 1),
        axis.title.y = element_text(size = 8,
                                    hjust = 1),
        plot.margin = margin(t = 5,
                             r = 8,
                             b = 3,
                             l = 3)) + 
  xlab(substitute(italic(PEMT)~"%"^"ile"~"in geneset correlations")) +
  ylab(substitute("Geneset"~"%"^"ile"~"in"~italic(PEMT)~"correlations"))

dev.off()

#### PEMT expression across tissues ####

GTEX_PEMTexpression_acrosstissues <- GTEX_MOR_across_tissues[PEMT_ensembl, ]

GTEX_PEMTexpression_acrosstissues_df <- data.frame(PEMTpseudocounts = GTEX_PEMTexpression_acrosstissues)

GTEX_PEMTexpression_acrosstissues_df[, "Tissue"] <- GTEX_LUT[match(row.names(GTEX_PEMTexpression_acrosstissues_df), GTEX_LUT$SAMPID), "Tissue"]
GTEX_PEMTexpression_acrosstissues_df[, "Tissue"] <- convert.GTEX.tissue.labels.for.plot(GTEX_PEMTexpression_acrosstissues_df$Tissue, supershort = TRUE)

GTEX_PEMTexpression_acrosstissues_geomean <- GTEX_PEMTexpression_acrosstissues_df %>% 
  group_by(Tissue) %>%
  summarise(geomean = gm_mean(PEMTpseudocounts))

GTEX_PEMTexpression_acrosstissues_df[, "Tissue"] <- factor(str_remove(GTEX_PEMTexpression_acrosstissues_df$Tissue, "^GTEX-"), levels = str_remove(unlist(GTEX_PEMTexpression_acrosstissues_geomean[order(GTEX_PEMTexpression_acrosstissues_geomean$geomean), "Tissue"]), "^GTEX-"))

pdf("graphics/GTEX_PEMTexpression_acrosstissues.pdf",
    width = 5,
    height = 2.5)

ggplot(aes(y = PEMTpseudocounts, x = Tissue), data = GTEX_PEMTexpression_acrosstissues_df) + 
  geom_boxplot(aes(fill = Tissue),
               outlier.size = 0.1) + 
  scale_y_log10(breaks = c(1, 10, 100, 1000, 10000, 100000, 1000000),
                expand = c(0, 0),
                limits = c(1, 1000000),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90,
                                   colour = "black",
                                   vjust = 0.5,
                                   hjust = 1),
        axis.text.y = element_text(colour = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10,
                                    hjust = 0.8),
        legend.position = "none",
        plot.margin = margin(t = 6, 
                             r = 3,
                             l = 3,
                             b = 3)) + 
  ylab(substitute(italic(PEMT)~"expression (pseudocounts)"))

dev.off()

#### plot PEMT expression vs reciprocal score #### 

PEMT_expression_vs_score <- data.frame(reciprocal_score = GTEX_tissue_HMT_PEMT_reciprocal_scores_df$reciprocal_score, PEMTexpression = unlist(GTEX_PEMTexpression_acrosstissues_geomean[match(GTEX_tissue_HMT_PEMT_reciprocal_scores_df$tissue, GTEX_PEMTexpression_acrosstissues_geomean$Tissue), "geomean"]))
PEMT_expression_vs_score[, "tissue"] <- GTEX_tissue_HMT_PEMT_reciprocal_scores_df$tissue

cor.test(log10(PEMT_expression_vs_score$reciprocal_score), log10(PEMT_expression_vs_score$PEMTexpression), method = "spearman")

pdf("graphics/GTEX_tissue_PEMTexpression_vs_reciprocalscore_scatterplot.pdf",
    width = 2.5,
    height = 2.5)

ggplot(aes(x = PEMTexpression, y = reciprocal_score), data = PEMT_expression_vs_score) + 
  geom_point() + 
  scale_x_log10(breaks = c(1, 10, 100, 1000, 10000, 100000, 1000000),
                expand = c(0, 0),
                limits = c(1, 1000000),
                labels = trans_format("log10", math_format(10^.x))) + 
  scale_y_log10(breaks = c(100, 1000, 10000, 100000),
                expand = c(0, 0),
                limits = c(50, 100000),
                labels = trans_format("log10", math_format(10^.x))) +
  geom_smooth(method = "lm",
              se = FALSE) +
  theme_classic() +
  theme(axis.text.x = element_text(colour = "black",
                                   size = 10),
        axis.text.y = element_text(colour = "black",
                                   size = 10),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10)) + 
  ylab(substitute("PEMTs-"~italic(PEMT)~"reciprocal score")) + 
  xlab(substitute(italic(PEMT)~"expression (pseudocounts)"))

dev.off()

GTEX_PEMTexpression_acrosstissues_df[, "reciprocal_score"] <- PEMT_expression_vs_score[match(GTEX_PEMTexpression_acrosstissues_df$Tissue, PEMT_expression_vs_score$tissue), "reciprocal_score"]

cor.test(log10(GTEX_PEMTexpression_acrosstissues_df$PEMTpseudocounts + 1), log10(GTEX_PEMTexpression_acrosstissues_df$reciprocal_score))$estimate^2
cor.test(log10(GTEX_PEMTexpression_acrosstissues_df$PEMTpseudocounts + 1), log10(GTEX_PEMTexpression_acrosstissues_df$reciprocal_score))$p.value

pdf("graphics/GTEX_tissue_PEMTexpression_vs_reciprocalscore_violin.pdf",
    height = 2.5,
    width = 5)

ggplot(data = GTEX_PEMTexpression_acrosstissues_df, aes(y = PEMTpseudocounts, x = reciprocal_score)) +
  geom_violin(aes(group = Tissue, fill = Tissue), 
              show.legend = FALSE, 
              scale = 1, 
              width = 0.3, 
              position = "identity", 
              lwd = 0.2) +
  scale_x_continuous(trans = reverselog_trans(10), 
                     breaks = c(1000, 100000, 10000000, 1000000000),
                     expand = c(0,0),
                     limits = c(1000000000, 1000),
                     labels = trans_format("log10", math_format(10^.x))) + 
  scale_y_log10(breaks = c(100, 1000, 10000, 100000),
                expand = c(0, 0),
                limits = c(50, 100000),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_classic() +
  theme(axis.text.x = element_text(colour = "black",
                                   size = 10),
        axis.text.y = element_text(colour = "black",
                                   size = 10),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10,
                                    hjust = 0.8)) + 
  xlab(substitute("PEMTs-"~italic(PEMT)~"reciprocal score")) + 
  ylab(substitute(italic(PEMT)~"expression (pseudocounts)")) + 
  geom_smooth(method = "lm",
              se = FALSE) + 
  annotate(geom = "text",
           x = 100000,
           y = 5000,
           label = substitute(italic(r)~"= 0.236")) +
  annotate(geom = "text",
           x = 100000,
           y = 10000,
           label = substitute(italic(p)~"= 0"))

dev.off()

#### HMT expression across tissues ####

GTEX_HMTexpression_acrosstissues <- colSums(GTEX_MOR_across_tissues[SET_HMTs_annotated$ensembl_gene_id, ])

GTEX_HMTexpression_acrosstissues_df <- data.frame(HMTpseudocounts = GTEX_HMTexpression_acrosstissues)

GTEX_HMTexpression_acrosstissues_df[, "Tissue"] <- GTEX_LUT[match(row.names(GTEX_HMTexpression_acrosstissues_df), GTEX_LUT$SAMPID), "Tissue"]
GTEX_HMTexpression_acrosstissues_df[, "Tissue"] <- convert.GTEX.tissue.labels.for.plot(GTEX_HMTexpression_acrosstissues_df$Tissue, supershort = TRUE)

GTEX_HMTexpression_acrosstissues_geomean <- GTEX_HMTexpression_acrosstissues_df %>% 
  group_by(Tissue) %>%
  summarise(geomean = gm_mean(HMTpseudocounts))

GTEX_HMTexpression_acrosstissues_df[, "Tissue"] <- factor(str_remove(GTEX_HMTexpression_acrosstissues_df$Tissue, "^GTEX-"), levels = str_remove(unlist(GTEX_HMTexpression_acrosstissues_geomean[order(GTEX_HMTexpression_acrosstissues_geomean$geomean), "Tissue"]), "^GTEX-"))

pdf("graphics/GTEX_HMTexpression_acrosstissues.pdf",
    width = 5,
    height = 2.5)

ggplot(aes(y = HMTpseudocounts, x = Tissue), data = GTEX_HMTexpression_acrosstissues_df) + 
  geom_boxplot(aes(fill = Tissue),
               outlier.size = 0.1) + 
  scale_y_log10() +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90,
                                   colour = "black",
                                   vjust = 0.5,
                                   hjust = 1),
        axis.text.y = element_text(colour = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10,
                                    hjust = 0.8),
        legend.position = "none",
        plot.margin = margin(t = 6, 
                             r = 3,
                             l = 3,
                             b = 3)) + 
  ylab(substitute(italic(HMT)~"expression (pseudocounts)"))

dev.off()

# plot HMT expression vs reciprocal score

HMT_expression_vs_score <- data.frame(reciprocal_score = GTEX_tissue_HMT_PEMT_reciprocal_scores_df$reciprocal_score, HMTexpression = unlist(GTEX_HMTexpression_acrosstissues_geomean[match(GTEX_tissue_HMT_PEMT_reciprocal_scores_df$tissue, GTEX_HMTexpression_acrosstissues_geomean$Tissue), "geomean"]))
HMT_expression_vs_score[, "tissue"] <- GTEX_tissue_HMT_PEMT_reciprocal_scores_df$tissue

cor.test(log10(HMT_expression_vs_score$reciprocal_score), log10(HMT_expression_vs_score$HMTexpression), method = "spearman")

pdf("graphics/GTEX_tissue_HMTexpression_vs_reciprocalscore_scatterplot.pdf",
    width = 2.5,
    height = 2.5)

ggplot(aes(x = HMTexpression, y = reciprocal_score), data = HMT_expression_vs_score) + 
  geom_point() + 
  scale_x_log10() + 
  scale_y_log10() +
  geom_smooth(method = "lm",
              se = FALSE) +
  theme_classic() +
  theme(axis.text.x = element_text(colour = "black",
                                   size = 10),
        axis.text.y = element_text(colour = "black",
                                   size = 10),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10)) + 
  ylab(substitute("HMTs-"~italic(HMT)~"reciprocal score")) + 
  xlab(substitute(italic(HMT)~"expression (pseudocounts)"))
dev.off()

GTEX_HMTexpression_acrosstissues_df[, "reciprocal_score"] <- HMT_expression_vs_score[match(GTEX_HMTexpression_acrosstissues_df$Tissue, HMT_expression_vs_score$tissue), "reciprocal_score"]

cor.test(log10(GTEX_HMTexpression_acrosstissues_df$HMTpseudocounts + 1), log10(GTEX_HMTexpression_acrosstissues_df$reciprocal_score))$estimate^2
cor.test(log10(GTEX_HMTexpression_acrosstissues_df$HMTpseudocounts + 1), log10(GTEX_HMTexpression_acrosstissues_df$reciprocal_score))$p.value

pdf("graphics/GTEX_tissue_HMTexpression_vs_reciprocalscore_violin.pdf",
    height = 2.5,
    width = 5)

ggplot(data = GTEX_HMTexpression_acrosstissues_df, aes(y = HMTpseudocounts, x = reciprocal_score)) +
  geom_violin(aes(group = Tissue, fill = Tissue), 
              show.legend = FALSE, 
              scale = 1, 
              width = 0.3, 
              position = "identity", 
              lwd = 0.2) +
  scale_x_log10() + 
  scale_y_log10() +
  theme_classic() +
  theme(axis.text.x = element_text(colour = "black",
                                   size = 10),
        axis.text.y = element_text(colour = "black",
                                   size = 10),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10,
                                    hjust = 0.8)) + 
  xlab(substitute("HMTs-"~italic(HMT)~"reciprocal score")) + 
  ylab(substitute(italic(HMT)~"expression (pseudocounts)")) + 
  geom_smooth(method = "lm",
              se = FALSE) + 
  annotate(geom = "text",
           x = 80000000,
           y = 500000,
           label = substitute(italic(r)~"= -0.251")) +
  annotate(geom = "text",
           x = 80000000,
           y = 100000,
           label = substitute(italic(p)~"= 3.38 x"~10^-141))

dev.off()

#### CANCER VS NORMAL for Fig S6A ####

# download TCGA metadata using 'recount' package
TCGA_metadata <- all_metadata(subset = "tcga", verbose = TRUE) 

# compile fields of interest into data frame
metadata_for_lm <- data.frame(cbind(TCGA_metadata$gdc_cases.samples.portions.analytes.aliquots.submitter_id,
                                    TCGA_metadata$gdc_cases.project.project_id,
                                    TCGA_metadata$gdc_cases.samples.sample_type,
                                    TCGA_metadata$gdc_cases.demographic.race,
                                    TCGA_metadata$gdc_cases.demographic.gender,
                                    TCGA_metadata$gdc_cases.diagnoses.tumor_stage,
                                    TCGA_metadata$gdc_cases.diagnoses.days_to_birth
))

# name columns
colnames(metadata_for_lm) <- c("SAMPID",
                               "cancer",
                               "sample_type",
                               "race",
                               "gender",
                               "tumour_stage",
                               "days_to_birth"
)

normsamples <- metadata_for_lm[metadata_for_lm$sample_type == "Solid Tissue Normal", "SAMPID"]

norm_by_cancer <- sapply(TCGA_MOR_list, function(x){
  sum(colnames(x) %in% normsamples)
})

normal_cancers <- names(norm_by_cancer[norm_by_cancer > 30])

normal_and_cancer_NNMT_HMT_corrs_and_scores <- lapply(normal_cancers, function(tempcancertype){

  normal_samples <- metadata_for_lm[metadata_for_lm$sample_type == "Solid Tissue Normal" &
                                      metadata_for_lm$cancer == tempcancertype , "SAMPID"]
  
  normal_donors <- str_extract(string = normal_samples, pattern = "^TCGA-[:alnum:]{2}-[:alnum:]{4}")
  
  matched_cancer_samples <- metadata_for_lm[metadata_for_lm$sample_type == "Primary Tumor" 
                                            & str_extract(string = metadata_for_lm$SAMPID, pattern = "^TCGA-[:alnum:]{2}-[:alnum:]{4}") %in% normal_donors, "SAMPID"]
  
  # remove duplicated donors (choose first one)
  matched_cancer_samples <- matched_cancer_samples[!duplicated(str_extract(matched_cancer_samples, pattern = "^TCGA-[:alnum:]{2}-[:alnum:]{4}"))]
  
  # ensure all samples are present in data
  matched_cancer_samples <- matched_cancer_samples[matched_cancer_samples %in% colnames(TCGA_MOR_list[[tempcancertype]])]
  normal_samples <- normal_samples[normal_samples %in% colnames(TCGA_MOR_list[[tempcancertype]])]
  
  finaldonors <- intersect(str_extract(string = matched_cancer_samples, pattern = "^TCGA-[:alnum:]{2}-[:alnum:]{4}"), 
            str_extract(string = normal_samples, pattern = "^TCGA-[:alnum:]{2}-[:alnum:]{4}"))
  
  matched_cancer_samples <- matched_cancer_samples[str_extract(string = matched_cancer_samples, pattern = "^TCGA-[:alnum:]{2}-[:alnum:]{4}") %in% finaldonors]
  normal_samples <- normal_samples[str_extract(string = normal_samples, pattern = "^TCGA-[:alnum:]{2}-[:alnum:]{4}") %in% finaldonors]
  
  normal_samples_MOR <- data.frame(TCGA_MOR_list[[tempcancertype]][, normal_samples])
  cancer_samples_MOR <- data.frame(TCGA_MOR_list[[tempcancertype]][, matched_cancer_samples])
  
  normal_samples_MOR["allHMTs", ] <- colSums(normal_samples_MOR[SET_HMTs$ensembl_gene_id, ])
  cancer_samples_MOR["allHMTs", ] <- colSums(cancer_samples_MOR[SET_HMTs$ensembl_gene_id, ])
    
  # now HMTs for normal samples
  
    all_normal_HMT_corrs <- vector(length = nrow(normal_samples_MOR))
    names(all_normal_HMT_corrs) <- row.names(normal_samples_MOR)
    
    all_normal_HMT_pvals <- vector(length = nrow(normal_samples_MOR))
    names(all_normal_HMT_pvals) <- row.names(normal_samples_MOR)
    
    normal_HMT_vec <- unlist(normal_samples_MOR["allHMTs", ])
    
    for(j in 1:nrow(normal_samples_MOR)){

      temp_vec <- unlist(normal_samples_MOR[j, ])
      
      crumb_ranks <- temp_vec + rnorm(temp_vec, 0, 0.01)
      
      tempcor <- cor.test(normal_HMT_vec,
                          crumb_ranks,
                          method = "spearman")
      
      all_normal_HMT_corrs[j] <- tempcor$estimate
      all_normal_HMT_pvals[j] <- tempcor$p.value
      
    }
    
    # now HMTs for cancer
    
    all_cancer_HMT_corrs <- vector(length = nrow(cancer_samples_MOR))
    names(all_cancer_HMT_corrs) <- row.names(cancer_samples_MOR)
    
    all_cancer_HMT_pvals <- vector(length = nrow(cancer_samples_MOR))
    names(all_cancer_HMT_pvals) <- row.names(cancer_samples_MOR)
    
    cancer_HMT_vec <- unlist(cancer_samples_MOR["allHMTs", ])
    
    for(j in 1:nrow(cancer_samples_MOR)){
      
      temp_vec <- unlist(cancer_samples_MOR[j, ])
      
      crumb_ranks <- temp_vec + rnorm(temp_vec, 0, 0.01)
      
      tempcor <- cor.test(cancer_HMT_vec,
                          crumb_ranks,
                          method = "spearman")
      
      all_cancer_HMT_corrs[j] <- tempcor$estimate
      all_cancer_HMT_pvals[j] <- tempcor$p.value
      
    }
    
    # now NNMT for normal samples
    
    all_normal_NNMT_corrs <- vector(length = nrow(normal_samples_MOR))
    names(all_normal_NNMT_corrs) <- row.names(normal_samples_MOR)
    
    all_normal_NNMT_pvals <- vector(length = nrow(normal_samples_MOR))
    names(all_normal_NNMT_pvals) <- row.names(normal_samples_MOR)
    
    normal_NNMT_vec <- unlist(normal_samples_MOR[NNMT_ensembl, ])
    
    for(j in 1:nrow(normal_samples_MOR)){
      
      temp_vec <- unlist(normal_samples_MOR[j, ])
      
      crumb_ranks <- temp_vec + rnorm(temp_vec, 0, 0.01)
      
      tempcor <- cor.test(normal_NNMT_vec,
                          crumb_ranks,
                          method = "spearman")
      
      all_normal_NNMT_corrs[j] <- tempcor$estimate
      all_normal_NNMT_pvals[j] <- tempcor$p.value
      
    }
    
    # now NNMT for cancer
    
    all_cancer_NNMT_corrs <- vector(length = nrow(cancer_samples_MOR))
    names(all_cancer_NNMT_corrs) <- row.names(cancer_samples_MOR)
    
    all_cancer_NNMT_pvals <- vector(length = nrow(cancer_samples_MOR))
    names(all_cancer_NNMT_pvals) <- row.names(cancer_samples_MOR)
    
    cancer_NNMT_vec <- unlist(cancer_samples_MOR[NNMT_ensembl, ])
    
    for(j in 1:nrow(cancer_samples_MOR)){
      
      temp_vec <- unlist(cancer_samples_MOR[j, ])
      
      crumb_ranks <- temp_vec + rnorm(temp_vec, 0, 0.01)
      
      tempcor <- cor.test(cancer_NNMT_vec,
                          crumb_ranks,
                          method = "spearman")
      
      all_cancer_NNMT_corrs[j] <- tempcor$estimate
      all_cancer_NNMT_pvals[j] <- tempcor$p.value
      
    }
    
    cancertypecorrelations <- list(normal_HMT_correlations = data.frame(cbind(all_normal_HMT_corrs, all_normal_HMT_pvals)),
                                  normal_NNMT_correlations = data.frame(cbind(all_normal_NNMT_corrs, all_normal_NNMT_pvals)),
                                  cancer_HMT_correlations = data.frame(cbind(all_cancer_HMT_corrs, all_cancer_HMT_pvals)),
                                  cancer_NNMT_correlations = data.frame(cbind(all_cancer_NNMT_corrs, all_cancer_NNMT_pvals)))

  normal_NNMTrank_in_HMTs <- which(names(all_normal_HMT_corrs[order(all_normal_HMT_corrs)]) == NNMT_ensembl)
  normal_HMTrank_in_NNMT <- which(names(all_normal_NNMT_corrs[order(all_normal_NNMT_corrs)]) == "allHMTs")
  
  normal_reciprocal_score = normal_NNMTrank_in_HMTs^2 + normal_HMTrank_in_NNMT^2
  
  cancer_NNMTrank_in_HMTs <- which(names(all_cancer_HMT_corrs[order(all_cancer_HMT_corrs)]) == NNMT_ensembl)
  cancer_HMTrank_in_NNMT <- which(names(all_cancer_NNMT_corrs[order(all_cancer_NNMT_corrs)]) == "allHMTs")
  
  cancer_reciprocal_score = cancer_NNMTrank_in_HMTs^2 + cancer_HMTrank_in_NNMT^2
  
  output_list = list(correlations = cancertypecorrelations,
                     normal_reciprocal_score = normal_reciprocal_score,
                     cancer_reciprocal_score = cancer_reciprocal_score)
  
})

names(normal_and_cancer_NNMT_HMT_corrs_and_scores) <- normal_cancers
saveRDS(normal_and_cancer_NNMT_HMT_corrs_and_scores, "output/normal_and_cancer_NNMT_HMT_corrs_and_scores.rds")
# normal_and_cancer_NNMT_HMT_corrs_and_scores <- readRDS("output/normal_and_cancer_NNMT_HMT_corrs_and_scores.rds")

normal_cancer_reciprocal_scores <- sapply(normal_and_cancer_NNMT_HMT_corrs_and_scores, function(x){
  c(x["normal_reciprocal_score"], x["cancer_reciprocal_score"])
})

# paired wilcox test for p value quoted in text
wilcox.test(log10(unlist(normal_cancer_reciprocal_scores[1, ])), log10(unlist(normal_cancer_reciprocal_scores[2, ])), paired = TRUE)

normal_cancer_reciprocal_scores_for_plot <- data.frame(t(normal_cancer_reciprocal_scores))
normal_cancer_reciprocal_scores_for_plot_df <- data.frame(apply(normal_cancer_reciprocal_scores_for_plot, 2, unlist))
normal_cancer_reciprocal_scores_for_plot_df[, "cancer"] <- str_remove(row.names(normal_cancer_reciprocal_scores_for_plot_df), "^TCGA-")

write.table(normal_cancer_reciprocal_scores_for_plot_df,
            file = "plot_data/Fig S6/Fig_S6A_plodata.txt",
            row.names = FALSE,
            col.names = TRUE,
            sep = "\t")

pdf("graphics/TCGA_cancer_normal_scatterplot.pdf",
    height = 2.5,
    width = 2.5)

ggplot(data = normal_cancer_reciprocal_scores_for_plot_df, aes(x = normal_reciprocal_score, y = cancer_reciprocal_score, label = cancer)) + 
  geom_point(colour = "red", size = 3, alpha = 0.7) +
  scale_x_continuous(trans = reverselog_trans(10),
                breaks = c(10000, 100000, 1000000, 10000000, 100000000, 1000000000, 10000000000),
                expand = c(0, 0),
                limits = c(10000000000, 10000),
                labels = trans_format("log10", math_format(10^.x))) + 
  scale_y_continuous(trans = reverselog_trans(10),
                breaks = c(10000, 100000, 1000000, 10000000, 100000000, 1000000000, 10000000000),
                expand = c(0, 0),
                limits = c(10000000000, 10000),
                labels = trans_format("log10", math_format(10^.x))) + 
  # coord_cartesian(xlim = c(10000, 1000000000),
  #                 ylim = c(10000, 1000000000)) +
  geom_label_repel(size = 2) +
  geom_abline(slope = 1,
              color = "grey",
              linetype = "dashed") + 
  theme_classic()+
  theme(axis.text.x = element_text(colour = "black",
                                           size = 10),
                axis.text.y = element_text(colour = "black",
                                           size = 10),
                axis.title.x = element_text(size = 10),
                axis.title.y = element_text(size = 10)) + 
  xlab("Normal relative reciprocal score") + 
  ylab("Cancer relative reciprocal score")

dev.off()

#### Cancer vs normal panels Fig S6 ####

normal_and_cancer_HMT_expression <- sapply(normal_cancers, function(tempcancertype){
  
  normal_samples <- metadata_for_lm[metadata_for_lm$sample_type == "Solid Tissue Normal" &
                                      metadata_for_lm$cancer == tempcancertype , "SAMPID"]
  
  normal_donors <- str_extract(string = normal_samples, pattern = "^TCGA-[:alnum:]{2}-[:alnum:]{4}")
  
  matched_cancer_samples <- metadata_for_lm[metadata_for_lm$sample_type == "Primary Tumor" 
                                            & str_extract(string = metadata_for_lm$SAMPID, pattern = "^TCGA-[:alnum:]{2}-[:alnum:]{4}") %in% normal_donors, "SAMPID"]
  
  # remove duplicated donors (choose first one)
  matched_cancer_samples <- matched_cancer_samples[!duplicated(str_extract(matched_cancer_samples, pattern = "^TCGA-[:alnum:]{2}-[:alnum:]{4}"))]
  
  # ensure all samples are present in data
  matched_cancer_samples <- matched_cancer_samples[matched_cancer_samples %in% colnames(TCGA_MOR_list[[tempcancertype]])]
  normal_samples <- normal_samples[normal_samples %in% colnames(TCGA_MOR_list[[tempcancertype]])]
  
  finaldonors <- intersect(str_extract(string = matched_cancer_samples, pattern = "^TCGA-[:alnum:]{2}-[:alnum:]{4}"), 
                           str_extract(string = normal_samples, pattern = "^TCGA-[:alnum:]{2}-[:alnum:]{4}"))
  
  matched_cancer_samples <- matched_cancer_samples[str_extract(string = matched_cancer_samples, pattern = "^TCGA-[:alnum:]{2}-[:alnum:]{4}") %in% finaldonors]
  normal_samples <- normal_samples[str_extract(string = normal_samples, pattern = "^TCGA-[:alnum:]{2}-[:alnum:]{4}") %in% finaldonors]
  
  normal_samples_MOR <- data.frame(TCGA_MOR_list[[tempcancertype]][, normal_samples])
  cancer_samples_MOR <- data.frame(TCGA_MOR_list[[tempcancertype]][, matched_cancer_samples])
  
  normal_samples_MOR["allHMTs", ] <- colSums(normal_samples_MOR[SET_HMTs$ensembl_gene_id, ])
  cancer_samples_MOR["allHMTs", ] <- colSums(cancer_samples_MOR[SET_HMTs$ensembl_gene_id, ])
  
  paired_test <- t.test(log10(unlist(normal_samples_MOR["allHMTs", ])), log10(unlist(cancer_samples_MOR["allHMTs", ])), paired = TRUE)
  
  boxplot(unlist(normal_samples_MOR["allHMTs", ]), unlist(cancer_samples_MOR["allHMTs", ]))
  
  output_vec <- c()
  output_vec["paired_t"] <- paired_test$statistic
  output_vec["p_val"] <- paired_test$p.value
  
  return(output_vec)
  
})

normal_and_cancer_HMT_expression_forplot <- lapply(normal_cancers, function(tempcancertype){
  
  normal_samples <- metadata_for_lm[metadata_for_lm$sample_type == "Solid Tissue Normal" &
                                      metadata_for_lm$cancer == tempcancertype , "SAMPID"]
  
  normal_donors <- str_extract(string = normal_samples, pattern = "^TCGA-[:alnum:]{2}-[:alnum:]{4}")
  
  matched_cancer_samples <- metadata_for_lm[metadata_for_lm$sample_type == "Primary Tumor" 
                                            & str_extract(string = metadata_for_lm$SAMPID, pattern = "^TCGA-[:alnum:]{2}-[:alnum:]{4}") %in% normal_donors, "SAMPID"]
  
  # remove duplicated donors (choose first one)
  matched_cancer_samples <- matched_cancer_samples[!duplicated(str_extract(matched_cancer_samples, pattern = "^TCGA-[:alnum:]{2}-[:alnum:]{4}"))]
  
  # ensure all samples are present in data
  matched_cancer_samples <- matched_cancer_samples[matched_cancer_samples %in% colnames(TCGA_MOR_list[[tempcancertype]])]
  normal_samples <- normal_samples[normal_samples %in% colnames(TCGA_MOR_list[[tempcancertype]])]
  
  finaldonors <- intersect(str_extract(string = matched_cancer_samples, pattern = "^TCGA-[:alnum:]{2}-[:alnum:]{4}"), 
                           str_extract(string = normal_samples, pattern = "^TCGA-[:alnum:]{2}-[:alnum:]{4}"))
  
  matched_cancer_samples <- matched_cancer_samples[str_extract(string = matched_cancer_samples, pattern = "^TCGA-[:alnum:]{2}-[:alnum:]{4}") %in% finaldonors]
  normal_samples <- normal_samples[str_extract(string = normal_samples, pattern = "^TCGA-[:alnum:]{2}-[:alnum:]{4}") %in% finaldonors]
  
  normal_samples_MOR <- data.frame(TCGA_MOR_list[[tempcancertype]][, normal_samples])
  cancer_samples_MOR <- data.frame(TCGA_MOR_list[[tempcancertype]][, matched_cancer_samples])
  
  output_list <- list()
  
  output_list[["normal"]] <- colSums(normal_samples_MOR[SET_HMTs$ensembl_gene_id, ])
  output_list[["cancer"]] <- colSums(cancer_samples_MOR[SET_HMTs$ensembl_gene_id, ])
  
  return(output_list)
  
})

names(normal_and_cancer_HMT_expression_forplot) <- normal_cancers

normal_and_cancer_HMT_expression_forplot_melt <- reshape2::melt(normal_and_cancer_HMT_expression_forplot)
colnames(normal_and_cancer_HMT_expression_forplot_melt) <- c("HMT",
                                                             "state",
                                                             "cancer_type")

normal_and_cancer_HMT_expression_forplot_melt$state <- factor(normal_and_cancer_HMT_expression_forplot_melt$state, 
                                                              level = c("normal", "cancer"))

normal_and_cancer_HMT_expression_forplot_melt$cancer_type <- str_remove(normal_and_cancer_HMT_expression_forplot_melt$cancer_type,
                                                                        "^TCGA-")

saveRDS(normal_and_cancer_HMT_expression_forplot_melt,
        "plot_data/normal_and_cancer_HMT_expression_forplot_melt.rds")

normal_and_cancer_HMT_expression_forplot_melt <- readRDS("plot_data/normal_and_cancer_HMT_expression_forplot_melt.rds")

write.table(normal_and_cancer_HMT_expression_forplot_melt,
            file = "plot_data/Fig S6/Fig_S6C_plotdata.txt",
            row.names = FALSE,
            col.names = TRUE,
            sep = "\t")

scientific_10 = function(x) {
  ifelse(
    x==0, "0",
    parse(text = sub("e[+]?", " %*% 10^", scientific_format()(x)))
  )
}

# Fig S6C

pdf("graphics/normal_cancer_HMT_expression.pdf",
    width = 3.75,
    height = 3)

ggplot(data = normal_and_cancer_HMT_expression_forplot_melt, aes(x = state, y = HMT)) + 
  geom_boxplot(fill = rep(c("grey", "red"), times = 12),
               lwd = 0.2,
               outlier.size = 0.2) + 
  facet_wrap(~ cancer_type) + 
  scale_y_log10(breaks = c(30000, 100000, 300000),
                expand = c(0, 0),
                limits = c(20000, 350000),
                label = scientific_10) +
  theme_classic() + 
  theme(strip.background = element_rect(fill = "black"),
        strip.text = element_text(colour = "white",
                                  face = "bold"),
        axis.text = element_text(colour = "black",
                                 size = 7)) +
  xlab("") + 
  ylab("Total HMT expression (pseudocounts)")

dev.off()

normal_and_cancer_HMT_expression <- data.frame(t(normal_and_cancer_HMT_expression))
normal_and_cancer_HMT_expression[, "FDR"] <- p.adjust(normal_and_cancer_HMT_expression$p_val, 
                                                      method = "BH")

normal_and_cancer_HMT_expression[, "cancer_type"] <- row.names(normal_and_cancer_HMT_expression)

# difference significant in 7 of 12 cancers
# of those 7, 5 significantly upregulated in cancer, 2 significantly downregulated 

normal_and_cancer_NNMT_expression <- sapply(normal_cancers, function(tempcancertype){
  
  normal_samples <- metadata_for_lm[metadata_for_lm$sample_type == "Solid Tissue Normal" &
                                      metadata_for_lm$cancer == tempcancertype , "SAMPID"]
  
  normal_donors <- str_extract(string = normal_samples, pattern = "^TCGA-[:alnum:]{2}-[:alnum:]{4}")
  
  matched_cancer_samples <- metadata_for_lm[metadata_for_lm$sample_type == "Primary Tumor" 
                                            & str_extract(string = metadata_for_lm$SAMPID, pattern = "^TCGA-[:alnum:]{2}-[:alnum:]{4}") %in% normal_donors, "SAMPID"]
  
  # remove duplicated donors (choose first one)
  matched_cancer_samples <- matched_cancer_samples[!duplicated(str_extract(matched_cancer_samples, pattern = "^TCGA-[:alnum:]{2}-[:alnum:]{4}"))]
  
  # ensure all samples are present in data
  matched_cancer_samples <- matched_cancer_samples[matched_cancer_samples %in% colnames(TCGA_MOR_list[[tempcancertype]])]
  normal_samples <- normal_samples[normal_samples %in% colnames(TCGA_MOR_list[[tempcancertype]])]
  
  finaldonors <- intersect(str_extract(string = matched_cancer_samples, pattern = "^TCGA-[:alnum:]{2}-[:alnum:]{4}"), 
                           str_extract(string = normal_samples, pattern = "^TCGA-[:alnum:]{2}-[:alnum:]{4}"))
  
  matched_cancer_samples <- matched_cancer_samples[str_extract(string = matched_cancer_samples, pattern = "^TCGA-[:alnum:]{2}-[:alnum:]{4}") %in% finaldonors]
  normal_samples <- normal_samples[str_extract(string = normal_samples, pattern = "^TCGA-[:alnum:]{2}-[:alnum:]{4}") %in% finaldonors]
  
  normal_samples_MOR <- data.frame(TCGA_MOR_list[[tempcancertype]][, normal_samples])
  cancer_samples_MOR <- data.frame(TCGA_MOR_list[[tempcancertype]][, matched_cancer_samples])
  
  normal_samples_MOR["NNMT", ] <- normal_samples_MOR[NNMT_ensembl, ]
  cancer_samples_MOR["NNMT", ] <- cancer_samples_MOR[NNMT_ensembl, ]
  
  paired_test <- t.test(log10(unlist(normal_samples_MOR["NNMT", ])), log10(unlist(cancer_samples_MOR["NNMT", ])), paired = TRUE)
  
  boxplot(unlist(normal_samples_MOR["NNMT", ]), unlist(cancer_samples_MOR["NNMT", ]))
  
  output_vec <- c()
  output_vec["paired_t"] <- paired_test$statistic
  output_vec["p_val"] <- paired_test$p.value
  
  return(output_vec)
  
})

normal_and_cancer_NNMT_expression_forplot <- lapply(normal_cancers, function(tempcancertype){
  
  normal_samples <- metadata_for_lm[metadata_for_lm$sample_type == "Solid Tissue Normal" &
                                      metadata_for_lm$cancer == tempcancertype , "SAMPID"]
  
  normal_donors <- str_extract(string = normal_samples, pattern = "^TCGA-[:alnum:]{2}-[:alnum:]{4}")
  
  matched_cancer_samples <- metadata_for_lm[metadata_for_lm$sample_type == "Primary Tumor" 
                                            & str_extract(string = metadata_for_lm$SAMPID, pattern = "^TCGA-[:alnum:]{2}-[:alnum:]{4}") %in% normal_donors, "SAMPID"]
  
  # remove duplicated donors (choose first one)
  matched_cancer_samples <- matched_cancer_samples[!duplicated(str_extract(matched_cancer_samples, pattern = "^TCGA-[:alnum:]{2}-[:alnum:]{4}"))]
  
  # ensure all samples are present in data
  matched_cancer_samples <- matched_cancer_samples[matched_cancer_samples %in% colnames(TCGA_MOR_list[[tempcancertype]])]
  normal_samples <- normal_samples[normal_samples %in% colnames(TCGA_MOR_list[[tempcancertype]])]
  
  finaldonors <- intersect(str_extract(string = matched_cancer_samples, pattern = "^TCGA-[:alnum:]{2}-[:alnum:]{4}"), 
                           str_extract(string = normal_samples, pattern = "^TCGA-[:alnum:]{2}-[:alnum:]{4}"))
  
  matched_cancer_samples <- matched_cancer_samples[str_extract(string = matched_cancer_samples, pattern = "^TCGA-[:alnum:]{2}-[:alnum:]{4}") %in% finaldonors]
  normal_samples <- normal_samples[str_extract(string = normal_samples, pattern = "^TCGA-[:alnum:]{2}-[:alnum:]{4}") %in% finaldonors]
  
  normal_samples_MOR <- data.frame(TCGA_MOR_list[[tempcancertype]][, normal_samples])
  cancer_samples_MOR <- data.frame(TCGA_MOR_list[[tempcancertype]][, matched_cancer_samples])
  
  output_list <- list()
  
  output_list[["normal"]] <- unlist(normal_samples_MOR[NNMT_ensembl, ])
  output_list[["cancer"]] <- unlist(cancer_samples_MOR[NNMT_ensembl, ])
  
  return(output_list)
  
})

names(normal_and_cancer_NNMT_expression_forplot) <- normal_cancers

normal_and_cancer_NNMT_expression_forplot_melt <- reshape2::melt(normal_and_cancer_NNMT_expression_forplot)
colnames(normal_and_cancer_NNMT_expression_forplot_melt) <- c("NNMT",
                                                              "state",
                                                              "cancer_type")

normal_and_cancer_NNMT_expression_forplot_melt$state <- factor(normal_and_cancer_NNMT_expression_forplot_melt$state, 
                                                               level = c("normal", "cancer"))

normal_and_cancer_NNMT_expression_forplot_melt$cancer_type <- str_remove(normal_and_cancer_NNMT_expression_forplot_melt$cancer_type, 
                                                                         "^TCGA-")

saveRDS(normal_and_cancer_NNMT_expression_forplot_melt,
        "plot_data/normal_and_cancer_NNMT_expression_forplot_melt.rds")

write.table(normal_and_cancer_NNMT_expression_forplot_melt,
            file = "plot_data/Fig S6/Fig_S6B_plotdata.txt",
            row.names = FALSE,
            col.names = TRUE,
            sep = "\t")

# Fig S6B

pdf("graphics/normal_cancer_NNMT_expression.pdf",
    width = 3.75,
    height = 3)

ggplot(data = normal_and_cancer_NNMT_expression_forplot_melt, aes(x = state, y = NNMT)) + 
  geom_boxplot(fill = rep(c("grey", "red"), times = 12),
               outlier.size = 0.2,
               lwd = 0.2) + 
  facet_wrap(~ cancer_type) + 
  scale_y_log10(breaks = c(10, 100, 1000, 10000, 100000, 1000000),
                expand = c(0, 0),
                limits = c(10, 1000000),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_classic() + 
  theme(strip.background = element_rect(fill = "black"),
        strip.text = element_text(colour = "white",
                                  face = "bold"),
        axis.text = element_text(colour = "black",
                                 size = 7)) +
  xlab("") + 
  ylab(italic(NNMT)~"expression (pseudocounts)")

dev.off()

normal_and_cancer_NNMT_expression <- data.frame(t(normal_and_cancer_NNMT_expression))
normal_and_cancer_NNMT_expression[, "FDR"] <- p.adjust(normal_and_cancer_NNMT_expression$p_val, 
                                                       method = "BH")

normal_cancer_HMT_NNMT_t_cor <- cor.test(normal_and_cancer_HMT_expression$paired_t, normal_and_cancer_NNMT_expression$paired_t)

normal_cancer_HMT_NNMT_t_foroplot <- data.frame(HMT = normal_and_cancer_HMT_expression$paired_t, 
                                                NNMT = normal_and_cancer_NNMT_expression$paired_t)

write.table(normal_cancer_HMT_NNMT_t_foroplot,
            file = "plot_data/Fig S6/Fig_S6D_plotdata.txt",
            row.names = TRUE,
            col.names = TRUE,
            sep = "\t")

# Fig S6D
pdf("graphics/HMTvsNNMT_normalvscancer_t.pdf",
    width = 3,
    height = 3)

ggplot(data = normal_cancer_HMT_NNMT_t_foroplot, aes(x = HMT, y = NNMT)) + 
  geom_point() + 
  geom_smooth(method = "lm",
              se = FALSE) + 
  ylab("NNMT normal vs cancer"~italic(t)) + 
  xlab("HMT normal vs cancer"~italic(t)) +
  theme_classic() + 
  coord_cartesian(ylim = c(-20, 20)) +
  geom_hline(yintercept = 0,
             colour = "grey",
             linetype = "dashed") +
  geom_vline(xintercept = 0,
             colour = "grey",
             linetype = "dashed") +
  annotate(geom = 'text',
           x = -3,
           y = -15,
           label = paste0("Pearson's r = ", round(normal_cancer_HMT_NNMT_t_cor$estimate, 3))) +
  annotate(geom = 'text',
           x = -3,
           y = -18,
           label = paste0("p = ", round(normal_cancer_HMT_NNMT_t_cor$p.value, 3)))

dev.off()

#### STEPWISE HMTs CORRELATIONS Fig 1H Fig 2F ####

# one_thousand_samples <- sapply(1:1000, function(x){sample(SET_HMTs$ensembl_gene_id, size = nrow(SET_HMTs), replace = FALSE)})
# saveRDS(one_thousand_samples, "output/one_thousand_HMTsamples.rds")
one_thousand_samples <- readRDS("output/one_thousand_HMTsamples.rds")

# samples are in columns

stepwise.samples <- function(database = c("GTEX", "TCGA"),
                             supplied_datalist = NULL,
                             subtype = NULL,
                             ensembl_gene = NULL,
                             stepwise_set = NULL){
  
  if(is.null(ensembl_gene)){
    
    stop("You must provide a gene ID")
    
  }
  
  if(length(stepwise_set) < 2){
    
    stop("You must provide a set of multiple genes for stepwise analysis")
    
  }
  
  match.arg(database)
  
  # load in MRN/MOR-corrected pseudocounts generated elsewhere
  
  message("Reading in data")
  
  if(!is.null(supplied_datalist)){
    
    datalist <- supplied_datalist
    
  } else {
  
  if(database == "TCGA") {
    
    datalist <- readRDS("output/TCGA_MOR_list.rds")
    TCGA_LUT <- readRDS("output/TCGA_LUT.rds")
    datalist <- lapply(datalist, function(x){x[, colnames(x) %in% TCGA_LUT$SAMPID]})
    
  }
  
  if(database == "GTEX") {
    
    datalist <- readRDS("output/GTEX_MOR_list.rds")

  }
    
  }
  
  if(!is.null(subtype)){
    
    datalist <- datalist[subtype]
    
  }
  
  if(!(ensembl_gene %in% row.names(datalist[[1]]))){
    
    stop("Gene is not present in chosen database")
    
  }
  
  if(any(!(stepwise_set %in% row.names(datalist[[1]])))){
    
    stop("At least one stepwise set member is not present in chosen database")
    
  }
  
  types <- names(datalist)
  
  favgene_ranks <- lapply(datalist, function(x){
    
    favgene_vec <- x[ensembl_gene, ]
    
    # add crumb - small random noise to break up ties. ranks will then be randomised among samples with 0 expression.
    favgene_crumb <- favgene_vec + rnorm(favgene_vec, 0, 0.01)
    
    favgene_ranks <- rank(favgene_crumb)
    
    favgene_rankspercent <- favgene_ranks / (length(favgene_ranks) + 1)
    
  })
  
  stepwise_list <- lapply(datalist, function(x){
    
    output_df <- data.frame(matrix(nrow = length(stepwise_set), ncol = ncol(x)))
    rownames(output_df) <- paste0("stepwise_", 1:length(stepwise_set))
    colnames(output_df) <- colnames(x)
    
    for (i in 1:length(stepwise_set)){
      
      if(i == 1){stepwise_vec <- x[stepwise_set[1:i], ]} else {
        stepwise_vec <- colSums(x[stepwise_set[1:i], ])}
      
      stepwise_vec_crumb <- stepwise_vec + rnorm(stepwise_vec, 0, 0.01)
      stepwise_vec_ranks <- rank(stepwise_vec_crumb)
      
      stepwise_vec_rankspercent <- stepwise_vec_ranks / (length(stepwise_vec_ranks) + 1)
      
      output_df[i, ] <- stepwise_vec_rankspercent
      
    }
    
    return(output_df)
    
  })

  if(!is.null(subtype)){
    
    favranks <- unlist(favgene_ranks[[1]])
    
    stepwise_ranks <- do.call(cbind, stepwise_list)
    colnames(stepwise_ranks) <- str_remove(colnames(stepwise_ranks), pattern = paste0("^", subtype, "."))
    
    stepwise_correlation <- apply(stepwise_ranks, 1, function(x){
      cor.test(x, favranks[match(colnames(stepwise_ranks), names(favranks))], method = "spearman")$estimate
    })
    
    output <- stepwise_correlation
    
  } else {
  
  iteratedlist <- sapply(1:100, function(i){
    
    chosensamples <- lapply(datalist, function(x){
      
      if(database == "GTEX"){
        n = 100
      }
      
      if(database == "TCGA"){
        n = 36
      }
      
      sample(colnames(x), n)
      
    })  
    
    for(i in 1:length(favgene_ranks)){
      
      favgene_ranks[[i]] <- favgene_ranks[[i]][chosensamples[[i]]]
      
    }
    
    favranks <- unlist(favgene_ranks)
    
    for(i in 1:length(stepwise_list)){
      
      stepwise_list[[i]] <- stepwise_list[[i]][chosensamples[[i]]]
      
    }
    
    stepwise_ranks <- do.call(cbind, stepwise_list)
    
    stepwise_correlation <- apply(stepwise_ranks, 1, function(x){
      cor.test(x, favranks[match(colnames(stepwise_ranks), names(favranks))], method = "spearman")$estimate
    })
    
    return(stepwise_correlation)
    
  }) 
  
  output <- rowMedians(iteratedlist)
  
  }
  
  return(output)
  
}

stepwise_NNMT_HMT_onethousand <- apply(one_thousand_samples, 2, FUN = function(x){
  
  stepwise.samples(database = "TCGA",
                   ensembl_gene = NNMT_ensembl,
                   stepwise_set = x)
  
})

saveRDS(stepwise_NNMT_HMT_onethousand, "output/stepwise_NNMT_HMT_onethousand.rds")
# stepwise_NNMT_HMT_onethousand  <- readRDS("output/stepwise_NNMT_HMT_onethousand.rds")

stepwisemelt <- melt(stepwise_NNMT_HMT_onethousand)

stepwisemelt_plot <- stepwisemelt
colnames(stepwisemelt_plot) <- c("pooledHMTs", "simulation", "rho")

write.table(stepwisemelt_plot,
            file = "plot_data/Fig 1/Fig_1H_plotdata.txt",
            row.names = TRUE,
            col.names = TRUE,
            sep = "\t")

pdf("graphics/TCGA_stepwise_1000.pdf", 
    width = 2.5, 
    height = 2.5)

ggplot(aes(x = Var1, y = value), data = stepwisemelt) + 
  geom_line(aes(group = Var2), alpha = 0.03) + 
  geom_smooth(col = "red", se = FALSE) + 
  theme_classic() + 
  theme(axis.text.x = element_text(colour = "black",
                                   size = 10),
        axis.text.y = element_text(colour = "black",
                                   size = 10),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10)) + 
  xlab("Number of pooled HMTs") + 
  ylab(substitute("Spearman's"~rho~"of HMT to"~paste(italic("NNMT")))) +
  geom_hline(yintercept = 0, linetype = "dashed", col = "grey")

dev.off()

GTEXstrongtissues_stepwise_PEMT_HMT_onethousand <- apply(one_thousand_samples, 2, FUN = function(x){
  print(x)
  stepwise.samples(database = "GTEX",
                   supplied_datalist = GTEX_MOR_list[GTEX_PEMT_HMT_strongtissues],
                   ensembl_gene = PEMT_ensembl,
                   stepwise_set = x)
  
})

saveRDS(GTEXstrongtissues_stepwise_PEMT_HMT_onethousand, "output/GTEXstrongtissues_stepwise_PEMT_HMT_onethousand.rds")
# GTEXstrongtissues_stepwise_PEMT_HMT_onethousand  <- readRDS("output/GTEXstrongtissues_stepwise_PEMT_HMT_onethousand.rds")

GTEX_strongtissues_stepwisemelt <- melt(GTEXstrongtissues_stepwise_PEMT_HMT_onethousand)

GTEX_strongtissues_stepwisemelt_plot <- GTEX_strongtissues_stepwisemelt
colnames(GTEX_strongtissues_stepwisemelt_plot) <- c("pooledHMTs", "simulation", "rho")

write.table(GTEX_strongtissues_stepwisemelt_plot,
            file = "plot_data/Fig 2/Fig_2F_plotdata.txt",
            row.names = FALSE,
            col.names = TRUE,
            sep = "\t")

pdf("graphics/GTEX_strongtissues_stepwise_1000.pdf", 
    width = 2.5, 
    height = 2.5)

ggplot(aes(x = Var1, y = value), data = GTEX_strongtissues_stepwisemelt) + 
  geom_line(aes(group = Var2), alpha = 0.03) + 
  geom_smooth(col = "red", se = FALSE) + 
  theme_classic() + 
  theme(axis.text.x = element_text(colour = "black",
                                   size = 10),
        axis.text.y = element_text(colour = "black",
                                   size = 10),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10)) + 
  xlab("Number of pooled HMTs") + 
  ylab(substitute("Spearman's"~rho~"of HMT to"~paste(italic("PEMT")))) +
  geom_hline(yintercept = 0, linetype = "dashed", col = "grey")

dev.off()

#### SETDB1 in melanoma Fig S4 ####

## violin

SKCM_NNMT <- TCGA_NNMT_bycancer[["TCGA-SKCM"]]

SKCM_NNMT[, "group"] <- "Other genes" 
SKCM_NNMT[row.names(SKCM_NNMT) %in% SET_HMTs$ensembl_gene_id, "group"] <- "HMTs"

SKCM_NNMT[, "group"] <- factor(SKCM_NNMT$group, levels = c("Other genes", "HMTs"))

# SKCM_NNMT_allHMTs <- SKCM_NNMT["allHMTs", ]
# SKCM_NNMT_allHMTs[, "group"] <- "HMTs"

# SKCM_NNMT <- SKCM_NNMT[row.names(SKCM_NNMT) != "allHMTs", ]
SKCM_NNMT[row.names(SKCM_NNMT) == "allHMTs", "group"] <- "HMTs"

SKCM_NNMT <- SKCM_NNMT[!str_detect(row.names(SKCM_NNMT), NNMT_ensembl), ]

SKCM_NNMT_HMTsubset <- SKCM_NNMT[row.names(SKCM_NNMT) %in% SET_HMTs$ensembl_gene_id, ]

for(i in 1:length(act_rep_list)){
  
  SKCM_NNMT_HMTsubset[act_rep_list[[i]], "colour"] <- selection_colours[i]
  
}

SKCM_NNMT_HMTsubset <- SKCM_NNMT_HMTsubset[SET_HMTs$ensembl_gene_id, ]

write.table(SKCM_NNMT,
            file = "plot_data/Fig S4/Fig_S4A_plotdata.txt",
            row.names = TRUE,
            col.names = TRUE,
            sep = "\t")

## scatterplots

SKCM_NNMT <- TCGA_NNMT_bycancer[["TCGA-SKCM"]]
SKCM_NNMT <- SKCM_NNMT[order(SKCM_NNMT$allcorrs), ]
SKCM_NNMT_ranks <- match(SET_HMTs$ensembl_gene_id, row.names(SKCM_NNMT))
names(SKCM_NNMT_ranks) <- SET_HMTs$hgnc_symbol

TCGA_SETDB1_bycancer <- analyse.by.subtype(database = "TCGA",
                                           singlegene = "ENSG00000143379",
                                           genebasename = "ENSG00000143379",
                                           extrageneset_ensembl = SET_HMTs$ensembl_gene_id,
                                           extragenesetbasename = "allHMTs")

saveRDS(TCGA_SETDB1_bycancer, "output/TCGA_SETDB1_with_everything_bycancer.rds")
# TCGA_SETDB1_bycancer <- readRDS("output/TCGA_SETDB1_with_everything_bycancer.rds")

# TCGA_MOR_list <- readRDS("output/TCGA_MOR_list.rds")
SKCM_MOR <- TCGA_MOR_list[["TCGA-SKCM"]]

SKCM_MOR_tumours <- SKCM_MOR[, colnames(SKCM_MOR) %in% TCGA_LUT$SAMPID]
ncol(SKCM_MOR_tumours)

SKCM_tumours_for_lm <- cbind(SKCM_MOR_tumours[NNMT_ensembl, ], SKCM_MOR_tumours["ENSG00000143379", ])

colnames(SKCM_tumours_for_lm) <- c("NNMT", "SETDB1")

SKCM_tumours_for_lm <- data.frame(cbind(SKCM_tumours_for_lm, TCGA_LUT[match(row.names(SKCM_tumours_for_lm), TCGA_LUT$SAMPID), c("race", "gender", "tumour_stage", "days_to_birth", "sequencing_centre")]))

SKCM_tumours_NNMTcorrected_vec <- lm(log10(NNMT) ~ race + gender + tumour_stage + as.numeric(days_to_birth), data = SKCM_tumours_for_lm)$residuals
SKCM_tumours_SETDB1corrected_vec <- lm(log10(SETDB1) ~ race + gender + tumour_stage + as.numeric(days_to_birth), data = SKCM_tumours_for_lm)$residuals

SKCM_NNMT_SETDB1_forplot <- data.frame(cbind(SKCM_tumours_NNMTcorrected_vec, SKCM_tumours_SETDB1corrected_vec))
colnames(SKCM_NNMT_SETDB1_forplot) <- c("NNMT_residuals", "SETDB1_residuals")

write.table(SKCM_NNMT_SETDB1_forplot,
            file = "plot_data/Fig S4/Fig_S4B_plotdata.txt",
            row.names = TRUE,
            col.names = TRUE,
            sep = "\t")

pdf("graphics/SETDB1_NNMT_melanomaprimary_scatter.pdf",
    height = 2.5,
    width = 2.5)

ggplot(aes(x = NNMT_residuals, y = SETDB1_residuals), data = SKCM_NNMT_SETDB1_forplot) + 
  geom_point(size = 0.5) + 
  theme_classic() +
  theme(axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10)) + 
  ylab(substitute(italic(SETDB1)~"expression residuals")) + 
  xlab(substitute(italic(NNMT)~"expression residuals")) + 
  annotate(geom = "text",
           x = -0.5,
           y = -0.325,
           label = substitute(italic(r)~"= -0.591")) +
  annotate(geom = "text",
           x = -0.5,
           y = -0.4,
           label = substitute(italic(p)~"= 6.28 x"~10^-11)) + 
  geom_smooth(method = "lm",
              se = FALSE,
              colour = "red")

dev.off()


# what about in the abundant metastasis samples, rather than only primary tumours? 

SKCM_MOR_metastasis <- SKCM_MOR[, str_detect(colnames(SKCM_MOR), "TCGA-[0-9A-Z]{2}-[0-9A-Z]{4}-06")]

# download TCGA metadata using 'recount' package
TCGA_metadata <- all_metadata(subset = "tcga", verbose = TRUE) 

# compile fields of interest into data frame
metadata_for_lm <- data.frame(cbind(TCGA_metadata$gdc_cases.samples.portions.analytes.aliquots.submitter_id,
                                    TCGA_metadata$gdc_cases.project.project_id,
                                    TCGA_metadata$gdc_cases.samples.sample_type,
                                    TCGA_metadata$gdc_cases.demographic.race,
                                    TCGA_metadata$gdc_cases.demographic.gender,
                                    TCGA_metadata$gdc_cases.diagnoses.tumor_stage,
                                    TCGA_metadata$gdc_cases.diagnoses.days_to_birth
))

# name columns
colnames(metadata_for_lm) <- c("SAMPID",
                               "cancer",
                               "sample_type",
                               "race",
                               "gender",
                               "tumour_stage",
                               "days_to_birth"
)

# Add field from barcode for sequencing centre
metadata_for_lm[, "sequencing_centre"] <- as.factor(str_extract(string = metadata_for_lm$SAMPID, pattern = "[:digit:]{2}$"))

# simplify tumour stage classification into stage 0, 1, 2, 3 or 4 (or not reported)
metadata_for_lm$tumour_stage <- str_replace_all(metadata_for_lm$tumour_stage, pattern = "stage i$", replacement = "stage_1")
metadata_for_lm$tumour_stage <- str_replace_all(metadata_for_lm$tumour_stage, pattern = "stage ia$", replacement = "stage_1")
metadata_for_lm$tumour_stage <- str_replace_all(metadata_for_lm$tumour_stage, pattern = "stage ib$", replacement = "stage_1")

metadata_for_lm$tumour_stage <- str_replace_all(metadata_for_lm$tumour_stage, pattern = "stage ii$", replacement = "stage_2")
metadata_for_lm$tumour_stage <- str_replace_all(metadata_for_lm$tumour_stage, pattern = "stage iia$", replacement = "stage_2")
metadata_for_lm$tumour_stage <- str_replace_all(metadata_for_lm$tumour_stage, pattern = "stage iib$", replacement = "stage_2")
metadata_for_lm$tumour_stage <- str_replace_all(metadata_for_lm$tumour_stage, pattern = "stage iic$", replacement = "stage_2")

metadata_for_lm$tumour_stage <- str_replace_all(metadata_for_lm$tumour_stage, pattern = "stage iii$", replacement = "stage_3")
metadata_for_lm$tumour_stage <- str_replace_all(metadata_for_lm$tumour_stage, pattern = "stage iiia$", replacement = "stage_3")
metadata_for_lm$tumour_stage <- str_replace_all(metadata_for_lm$tumour_stage, pattern = "stage iiib$", replacement = "stage_3")
metadata_for_lm$tumour_stage <- str_replace_all(metadata_for_lm$tumour_stage, pattern = "stage iiic$", replacement = "stage_3")

metadata_for_lm$tumour_stage <- str_replace_all(metadata_for_lm$tumour_stage, pattern = "stage iv$", replacement = "stage_4")
metadata_for_lm$tumour_stage <- str_replace_all(metadata_for_lm$tumour_stage, pattern = "stage iva$", replacement = "stage_4")
metadata_for_lm$tumour_stage <- str_replace_all(metadata_for_lm$tumour_stage, pattern = "stage ivb$", replacement = "stage_4")
metadata_for_lm$tumour_stage <- str_replace_all(metadata_for_lm$tumour_stage, pattern = "stage ivc$", replacement = "stage_4")

# remove duplicate entries (duplicated sample ID)
metadata_for_lm <- metadata_for_lm[!(duplicated(metadata_for_lm$SAMPID)),]

# remove entries where multiple samples are derived from one donor. keep only donors that provide one cancer sample.
metadata_for_lm <- metadata_for_lm[(!duplicated(str_extract(string = metadata_for_lm$SAMPID, pattern = "^TCGA-[:alnum:]{2}-[:alnum:]{4}-[:alnum:]{2}")) & !duplicated(str_extract(string = metadata_for_lm$SAMPID, pattern = "^TCGA-[:alnum:]{2}-[:alnum:]{4}-[:alnum:]{2}"), fromLast = TRUE)), ]

# set SAMPID as row names
row.names(metadata_for_lm) <- metadata_for_lm$SAMPID

SKCM_metastasis_for_lm <- data.frame(cbind(SKCM_MOR_metastasis[NNMT_ensembl, ], SKCM_MOR_metastasis["ENSG00000143379", ]))
colnames(SKCM_metastasis_for_lm) <- c("NNMT", "SETDB1")

SKCM_metastasis_for_lm <- data.frame(cbind(SKCM_metastasis_for_lm, metadata_for_lm[match(row.names(SKCM_metastasis_for_lm), metadata_for_lm$SAMPID), c("race", "gender", "tumour_stage", "days_to_birth", "sequencing_centre")]))

SKCM_metastasis_for_lm$tumour_stage <- str_replace(SKCM_metastasis_for_lm$tumour_stage, pattern = "i/ii nos", replacement = "stage_1.5")
SKCM_metastasis_for_lm <- SKCM_metastasis_for_lm[SKCM_metastasis_for_lm$tumour_stage !="not reported", ]
SKCM_metastasis_for_lm <- SKCM_metastasis_for_lm[!is.na(SKCM_metastasis_for_lm$days_to_birth), ]


SKCM_metastasis_NNMTcorrected_vec <- lm(log10(NNMT) ~ race + gender + tumour_stage + as.numeric(days_to_birth), data = SKCM_metastasis_for_lm)$residuals
SKCM_metastasis_SETDB1corrected_vec <- lm(log10(SETDB1) ~ race + gender + tumour_stage + as.numeric(days_to_birth), data = SKCM_metastasis_for_lm)$residuals

SKCMmetastasis_NNMT_SETDB1_forplot <- data.frame(cbind(SKCM_metastasis_NNMTcorrected_vec, SKCM_metastasis_SETDB1corrected_vec))
colnames(SKCMmetastasis_NNMT_SETDB1_forplot) <- c("NNMT_residuals", "SETDB1_residuals")

cor.test(SKCM_metastasis_NNMTcorrected_vec, SKCM_metastasis_SETDB1corrected_vec)$p.value

write.table(SKCMmetastasis_NNMT_SETDB1_forplot,
            file = "plot_data/Fig S4/Fig_S4C_plotdata.txt",
            row.names = TRUE,
            col.names = TRUE,
            sep = "\t")

pdf("graphics/SETDB1_NNMT_melanomametastasis_scatter.pdf",
    height = 2.5,
    width = 2.5)

ggplot(aes(x = NNMT_residuals, y = SETDB1_residuals), data = SKCMmetastasis_NNMT_SETDB1_forplot) + 
  geom_point(size = 0.5) + 
  theme_classic() +
  theme(axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10)) + 
  ylab(substitute(italic(SETDB1)~"expression residuals")) + 
  xlab(substitute(italic(NNMT)~"expression residuals")) + 
  annotate(geom = "text",
           x = -0.8,
           y = -0.325,
           label = substitute(italic(r)~"= -0.472")) +
  annotate(geom = "text",
           x = -0.8,
           y = -0.43,
           label = substitute(italic(p)~"= 1.43 x"~10^-19)) + 
  geom_smooth(method = "lm",
              se = FALSE,
              colour = "red")

dev.off()

# what about the expression level of SETDB1 in SKCM - vs other HMTs and vs other cancers

SKCM_MOR_HMTs <- SKCM_MOR[SET_HMTs$ensembl_gene_id, colnames(SKCM_MOR) %in% TCGA_LUT$SAMPID]
SKCM_MOR_HMTs_melt <- melt(SKCM_MOR_HMTs)

SKCM_HMTs_mean_expression <- SKCM_MOR_HMTs_melt %>%
  group_by(Var1) %>%
  summarise(geomean = gm_mean(value))

SKCM_HMTs_mean_expression[, "hgnc_symbol"] <- SET_HMTs[match(SKCM_HMTs_mean_expression$Var1, SET_HMTs$ensembl_gene_id), "hgnc_symbol"]

SKCM_MOR_HMTs_melt$Var1 <- SET_HMTs[match(SKCM_MOR_HMTs_melt$Var1, SET_HMTs$ensembl_gene_id), "hgnc_symbol"]

SKCM_MOR_HMTs_melt$Var1 <- factor(SKCM_MOR_HMTs_melt$Var1, levels = unlist(SKCM_HMTs_mean_expression[order(SKCM_HMTs_mean_expression$geomean), "hgnc_symbol"]))

SKCM_MOR_HMTs_melt_copy <- SKCM_MOR_HMTs_melt
colnames(SKCM_MOR_HMTs_melt_copy) <- c("HMT", "SAMPID", "pseudocounts")

write.table(SKCM_MOR_HMTs_melt_copy,
            file = "plot_data/Fig S4/Fig_S4D_plotdata.txt",
            row.names = FALSE,
            col.names = TRUE,
            sep = "\t")

pdf("graphics/HMTs_expression_in_primary_melanoma.pdf",
    height = 2.5,
    width = 5)

ggplot(aes(y = value, x = Var1), data = SKCM_MOR_HMTs_melt) + 
  geom_boxplot(aes(fill = Var1),
               outlier.size = 0.1,
               lwd = 0.3) + 
  scale_y_log10(breaks = c(1, 10, 100, 1000, 10000, 100000),
                expand = c(0, 0),
                limits = c(1, 100000),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90,
                                   colour = "black",
                                   face = "italic",
                                   size = 10,
                                   vjust = 0.5,
                                   hjust = 1
  ),
  axis.text.y = element_text(colour = "black",
                             size = 10),
  axis.title.x = element_blank(),
  axis.title.y = element_text(size = 10,
                              hjust = 0.8),
  legend.position = "none") + 
  ylab("Gene expression (pseudocounts)")

dev.off()

# now for all cancers

# TCGA_MOR_across_cancers <- readRDS("~/manuscript/output/TCGA-MOR-normalisation-across-cancers.rds")
SETDB1_across_cancers <- TCGA_MOR_across_cancers["ENSG00000143379", ]

SETDB1_across_cancers <- data.frame(SETDB1_across_cancers)
SETDB1_across_cancers[, "cancer"] <- TCGA_LUT[match(row.names(SETDB1_across_cancers), TCGA_LUT$SAMPID), "cancer"]

SETDB1_across_cancers_geomean <- SETDB1_across_cancers %>% 
  group_by(cancer) %>%
  summarise(geomean = gm_mean(SETDB1_across_cancers))

SETDB1_across_cancers$cancer <- factor(str_remove(SETDB1_across_cancers$cancer, "^TCGA-"), levels = str_remove(unlist(SETDB1_across_cancers_geomean[order(SETDB1_across_cancers_geomean$geomean), "cancer"]), "^TCGA-"))

write.table(SETDB1_across_cancers,
            file = "plot_data/Fig S4/Fig_S4E_plotdata.txt",
            row.names = TRUE,
            col.names = TRUE,
            sep = "\t")

pdf("graphics/SETDB1_expression_across_cancers.pdf",
    height = 2.5,
    width = 5)

ggplot(aes(y = SETDB1_across_cancers, x = cancer), data = SETDB1_across_cancers) + 
  geom_boxplot(aes(fill = cancer),
               outlier.size = 0.1,
               lwd = 0.3) + 
  scale_y_log10(breaks = c(100, 1000, 10000),
                expand = c(0, 0),
                limits = c(100, 10000),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90,
                                   colour = "black",
                                   size = 10,
                                   vjust = 0.5,
                                   hjust = 1
  ),
  axis.text.y = element_text(colour = "black",
                             size = 10),
  axis.title.x = element_blank(),
  axis.title.y = element_text(size = 10,
                              hjust = 0.8),
  legend.position = "none") +
  ylab(substitute(italic(SETDB1)~"expression (pseudocounts)"))

dev.off()

#### OTHER METHYLTRANSFERASE SETS for Fig S5 Fig S7 ####

# RNA methyltransferases from here:
# https://pubmed.ncbi.nlm.nih.gov/26566070/
mRNA_MTs <- c("RNMT",
              "METTL14",
              "METTL3",
              "CMTR1",
              "CMTR2")

smRNA_MTs <- c("HENMT1",
               "MEPCE",
               "TGS1",
               "BCDIN3D",
               "METTL4")

rRNA_MTs <- c("DIMT1",
              "EMG1",
              "FBL",
              "FBLL1",
              "FTSJ2",
              "FTSJ3",
              "MRM1",
              "NSUN4",
              "NSUN5",
              "RNMTL1",
              "TFB1M",
              "TFB2M",
              "TRMT112",
              "WBSCR22")

tRNA_MTs <- c("ALKBH8",
              "CDK5RAP1",
              "CDKAL1",
              "ELP3",
              "FTSJ1",
              "KIAA1456",
              "METTL1",
              "NSUN2",
              "NSUN6",
              "TRDMT1",
              "TRMT1",
              "TRMT10A",
              "TRMT10B",
              "TRMT10C",
              "TRMT11",
              "TRMT12",
              "TRMT2B",
              "TRMT44",
              "TRMT5",
              "TRMT61A",
              "TRMT61B",
              "TYW3")

unknownRNA_MTs <- c("SPOUT1",
                    "NSUN3",
                    "NSUN5P1",
                    "NSUN5P2",
                    "NSUN7",
                    "RSAD1",
                    "TARBP1",
                    "THUMPD2",
                    "THUMPD3",
                    "TRMT1L",
                    "TRMT2A")

DNA_MTs <- c("DNMT1",
             "DNMT3A",
             "DNMT3B")

PRMTs <- c(paste0("PRMT", 1:11))

RNA_MT_list <- list(mRNA_MTs = mRNA_MTs,
                    smRNA_MTs = smRNA_MTs,
                    rRNA_MTs = rRNA_MTs,
                    tRNA_MTs = tRNA_MTs,
                    unknownRNA_MTs = unknownRNA_MTs,
                    DNA_MTs = DNA_MTs,
                    PRMTs = PRMTs
                    )

RNA_MT_ensembl_list <- lapply(RNA_MT_list, function(x){
  getBM(attributes = c("ensembl_gene_id", "description", "hgnc_symbol"),
        filters = "hgnc_symbol",
        values = x,
        mart = ensembl)
})

saveRDS(RNA_MT_ensembl_list, "output/RNA_MT_ensembl_list.rds")
RNA_MT_ensembl_list <- readRDS("output/RNA_MT_ensembl_list.rds")

RNA_MT_ensembl_list_melt <- melt(RNA_MT_ensembl_list)
colnames(RNA_MT_ensembl_list_melt)[ncol(RNA_MT_ensembl_list_melt)] <- "category"

write.xlsx(RNA_MT_ensembl_list_melt,
           "input/RNA_MT.xlsx")

RNA_MT_NNMT_reciprocal_scores <- lapply(names(RNA_MT_ensembl_list), function(x){

  if(!file.exists(paste0("output/TCGA_", x, "_everything_object.rds"))){

  TCGA_thisone_everything_object <- analyse.by.subtype(database = "TCGA",
                                                             geneset_ensembl = RNA_MT_ensembl_list[[x]]$ensembl_gene_id,
                                                             genebasename = x)
  
  saveRDS(TCGA_thisone_everything_object, paste0("output/TCGA_", x, "_everything_object.rds"))
  
  } else {
    
    TCGA_thisone_everything_object <- readRDS(paste0("output/TCGA_", x, "_everything_object.rds"))
    
  }
  
  if(!file.exists(paste0("output/TCGA_NNMT_everything_object_", x, "added.rds"))){
  
  thisoneadded <- add.geneset.to.everything.subtype(database = "TCGA",
                                          geneset_ensembl = RNA_MT_ensembl_list[[x]]$ensembl_gene_id,
                                          genesetname = x,
                                          everything_object_list = TCGA_NNMT_bycancer)
  
  saveRDS(thisoneadded, paste0("output/TCGA_NNMT_everything_object_", x, "added.rds"))
  
  } else {
    
    thisoneadded <- readRDS(paste0("output/TCGA_NNMT_everything_object_", x, "added.rds"))
    
  }
  
  tempscorelist <- lapply(1:length(TCGA_thisone_everything_object), function(i){

    TCGA_thisone_vec <- TCGA_thisone_everything_object[[i]]$allcorrs
    names(TCGA_thisone_vec) <- row.names(TCGA_thisone_everything_object[[i]])
    
    thisoneadded_vec <- thisoneadded[[i]]$allcorrs
    names(thisoneadded_vec) <- row.names(thisoneadded[[i]])

    calculate.relative.reciprocal.scores(everything_object1 = TCGA_thisone_vec,
                                         everything_object2 = thisoneadded_vec)
    
  })
  
  tempscore_df <- do.call(rbind, tempscorelist)
  row.names(tempscore_df) <- names(TCGA_thisone_everything_object)
  
  return(tempscore_df)
  
})

names(RNA_MT_NNMT_reciprocal_scores) <- names(RNA_MT_ensembl_list)
saveRDS(RNA_MT_NNMT_reciprocal_scores, "output/RNA_MT_NNMT_reciprocal_scores.rds")
RNA_MT_NNMT_reciprocal_scores <- readRDS("output/RNA_MT_NNMT_reciprocal_scores.rds")

RNA_MT_ensembl_list <- lapply(RNA_MT_ensembl_list, function(x){
  
  x[x$ensembl_gene_id %in% row.names(TCGA_MOR_list[[1]]), ]
  
})

RNA_MT_NNMT_reciprocal_scores <- lapply(names(RNA_MT_NNMT_reciprocal_scores), function(thisone){
  
  x <- data.frame(RNA_MT_NNMT_reciprocal_scores[[thisone]])
  x[, "cancer"] <- str_remove(row.names(x), "^TCGA-")
  
  colnames(x) <- str_replace(colnames(x), thisone, "MTset")
  
  x[, "MTs_in_NNMT_percent"] <- (unlist(x[, 1]) * 100)/60489
  x[, "NNMT_in_MTs_percent"] <- (unlist(x[, 2]) * 100)/60489
  
  x[x$NNMT_in_MTs_percent > 2.5 | x$MTs_in_NNMT_percent > 2.5, "cancer"] <- ""
  
  return(x)
  
})

names(RNA_MT_NNMT_reciprocal_scores) <- names(RNA_MT_ensembl_list)

RNA_MT_NNMT_reciprocal_scores_for_facet <- lapply(names(RNA_MT_NNMT_reciprocal_scores), function(thisone){
  
  x <- RNA_MT_NNMT_reciprocal_scores[[thisone]]
  
  x[, "MTset"] <- thisone
  
  return(x)
  
})

RNA_MT_NNMT_reciprocal_scores_for_facet <- do.call(rbind, RNA_MT_NNMT_reciprocal_scores_for_facet)

write.table(RNA_MT_NNMT_reciprocal_scores_for_facet,
            file = "plot_data/Fig S5/Fig_S5_bubble_plotdata.txt",
            row.names = TRUE,
            col.names = TRUE,
            sep = "\t")

pdf("graphics/TCGA_NNMT_otherMTs_facet.pdf",
    height = 8,
    width = 5)

ggplot(aes(x = NNMT_in_MTs_percent, y = MTs_in_NNMT_percent, label = cancer), data = RNA_MT_NNMT_reciprocal_scores_for_facet) + 
  geom_rect(aes(xmax = 2.5,
                xmin = 0,
                ymax = 2.5,
                ymin = 0),
            fill = "white",
            linetype = "dashed",
            color = "grey") +
  geom_point(aes(size = 10/log10(reciprocal_score)),
             colour = "red",
             alpha = 0.7) + 
  scale_x_log10(breaks = c(0.01, 0.1, 1, 10, 100),
                expand = c(0, 0),
                limits = c(0.001, 100),
                labels = label_number(drop0trailing = TRUE)) +
  scale_y_log10(breaks = c(0.01, 0.1, 1, 10, 100),
                expand = c(0, 0),
                limits = c(0.001, 100),
                labels = label_number(drop0trailing = TRUE)) +
  geom_label_repel(size = 1.5,
                   max.overlaps = 40,
                   label.size = 0.1,
                   box.padding = 0.1,
                   label.padding = 0.1) + 
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(colour = "black",
                                   size = 10),
        axis.text.y = element_text(colour = "black",
                                   size = 10),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.line = element_line(),
        plot.margin = margin(t = 8,
                             r = 8,
                             b = 8,
                             l = 8)) + 
  xlab(substitute(italic(NNMT)~"%"^"ile"~"in"~thisonelabel~"correlations")) +
  ylab(substitute(thisonelabel~"%"^"ile"~"in"~italic(NNMT)~"correlations")) + 
  facet_wrap(~ MTset, scales = "free")

dev.off()

pdf("graphics/TCGA_NNMT_otherMTs_facet_nostrip.pdf",
    height = 4,
    width = 4)

ggplot(aes(x = NNMT_in_MTs_percent, y = MTs_in_NNMT_percent, label = cancer), data = RNA_MT_NNMT_reciprocal_scores_for_facet) + 
  geom_rect(aes(xmax = 2.5,
                xmin = 0,
                ymax = 2.5,
                ymin = 0),
            fill = "white",
            linetype = "dashed",
            color = "grey") +
  geom_point(aes(size = 10/log10(reciprocal_score)),
             colour = "red",
             alpha = 0.7) + 
  scale_x_log10(breaks = c(0.01, 0.1, 1, 10, 100),
                expand = c(0, 0),
                limits = c(0.001, 100),
                labels = label_number(drop0trailing = TRUE)) +
  scale_y_log10(breaks = c(0.01, 0.1, 1, 10, 100),
                expand = c(0, 0),
                limits = c(0.001, 100),
                labels = label_number(drop0trailing = TRUE)) +
  geom_label_repel(size = 1.5,
                   max.overlaps = 40,
                   label.size = 0.1,
                   box.padding = 0.1,
                   label.padding = 0.1) + 
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(colour = "black",
                                   size = 8,
                                   angle = 30),
        axis.text.y = element_text(colour = "black",
                                   size = 8),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.line = element_line(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.margin = margin(t = 8,
                             r = 8,
                             b = 8,
                             l = 8)) + 
  xlab(substitute(italic(NNMT)~"%"^"ile"~"in"~thisonelabel~"correlations")) +
  ylab(substitute(thisonelabel~"%"^"ile"~"in"~italic(NNMT)~"correlations")) + 
  facet_wrap(~ MTset, scales = "free")

dev.off()

lapply(names(RNA_MT_NNMT_reciprocal_scores), function(thisone){

  thisonelabel <- str_replace(thisone, "_", " ")
  
  thisdata <- RNA_MT_NNMT_reciprocal_scores[[thisone]]
  thisdata[(thisdata$MTs_in_NNMT_percent > 2.5) | (thisdata$NNMT_in_MTs_percent > 2.5), "cancer"] <- ""
  
  pdf(paste0("graphics/TCGA-", thisone, "-NNMT-cancers-reciprocal_scores.pdf"),
    width = 1.75,
    height = 1.75)

print(
ggplot(aes(x = NNMT_in_MTs_percent, y = MTs_in_NNMT_percent, label = cancer), data = thisdata) + 
  geom_rect(aes(xmax = 2.5,
                xmin = 0,
                ymax = 2.5,
                ymin = 0),
            fill = "white",
            linetype = "dashed",
            color = "grey") +
  geom_point(aes(size = 10/log10(reciprocal_score)),
             colour = "red",
             alpha = 0.7) + 
  scale_x_log10(breaks = c(0.01, 0.1, 1, 10, 100),
                expand = c(0, 0),
                limits = c(0.001, 100),
                labels = label_number(drop0trailing = TRUE)) +
  scale_y_log10(breaks = c(0.01, 0.1, 1, 10, 100),
                expand = c(0, 0),
                limits = c(0.001, 100),
                labels = label_number(drop0trailing = TRUE)) +
  geom_label_repel(size = 1.5,
                   max.overlaps = 40,
                   label.size = 0.1,
                   box.padding = 0.1,
                   label.padding = 0.1) + 
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(colour = "black",
                                   size = 10),
        axis.text.y = element_text(colour = "black",
                                   size = 10),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.margin = margin(t = 8,
                             r = 8,
                             b = 8,
                             l = 8)) + 
  xlab(substitute(italic(NNMT)~"%"^"ile"~"in"~thisonelabel~"correlations")) +
  ylab(substitute(thisonelabel~"%"^"ile"~"in"~italic(NNMT)~"correlations"))

)

dev.off()

})

# do stepwise hockeysticks for these

RNA_MT_samples_list <- lapply(RNA_MT_ensembl_list, function(x){

  number_of_genes <- nrow(x)

  if(number_of_genes <= 5){
  
  permutations <- permutations(number_of_genes, number_of_genes)
  
  samples <- apply(permutations, 1, function(y){

    x$ensembl_gene_id[y]
    
  })
  
  } else {
    
    samples <- sapply(1:120, function(i){sample(x$ensembl_gene_id, size = nrow(x), replace = FALSE)})
    
  }
  
  return(samples)
  
})

saveRDS(RNA_MT_samples_list, "output/otherMT_stepwise_samples_list.rds")
RNA_MT_samples_list <- readRDS("output/otherMT_stepwise_samples_list.rds")

RNA_hockeystick <- lapply(names(RNA_MT_samples_list), function(thisname){

stepwise_NNMT_MT <- apply(RNA_MT_samples_list[[thisname]], 2, FUN = function(genesamples){

 stepwise.samples(database = "TCGA",
                   ensembl_gene = NNMT_ensembl,
                   stepwise_set = genesamples)
  
})

saveRDS(stepwise_NNMT_MT, paste0("output/", thisname, "_stepwise_samples_list.rds"))

})

outputfiles <- list.files("output")
hockeystickfiles <- c(outputfiles[str_detect(outputfiles, "_MTs_stepwise_")], "PRMTs_stepwise_samples_list.rds")
hockeystickfiles <- hockeystickfiles[!str_detect(hockeystickfiles, "GTEX")]

for(i in 1:length(hockeystickfiles)){

  temphock <- readRDS(paste0("output/", hockeystickfiles[i]))
  
  colnames(temphock) <- paste0("permutation_", 1:ncol(temphock))
  
  write.table(temphock,
              file = paste0("plot_data/Fig S5/Fig_S5_", str_remove(hockeystickfiles[i], "_stepwise_samples_list.rds"), "_hockeystick_plotdata.txt"),
              row.names = FALSE,
              col.names = TRUE,
              sep = "\t")
  
}
  
# plot them out

lapply(hockeystickfiles, function(x){

  tempdata <- readRDS(paste0("output/",x))
  tempname <- str_remove(x, "_stepwise_samples_list.rds")
  
  tempmelt <- melt(tempdata)
  
  pdf(paste0("graphics/TCGA_stepwise_", tempname, ".pdf"), 
      width = 2.5, 
      height = 2.5)
  
  print(
  ggplot(aes(x = Var1, y = value), data = tempmelt) + 
    geom_line(aes(group = Var2), alpha = 0.03) + 
    geom_smooth(col = "red", se = FALSE) + 
    theme_classic() + 
    theme(axis.text.x = element_text(colour = "black",
                                     size = 10),
          axis.text.y = element_text(colour = "black",
                                     size = 10),
          axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10)) + 
    xlab("Number of pooled MTs") + 
    ylab(substitute("Spearman's"~rho~"of MTs to"~paste(italic("NNMT")))) +
    geom_hline(yintercept = 0, linetype = "dashed", col = "grey")
  )
  
  dev.off()
  
})

# checking out unknwon RNA MTs

TCGA_LUSC <- TCGA_MOR_list[["TCGA-LUSC"]][, colnames(TCGA_MOR_list[["TCGA-LUSC"]]) %in% TCGA_LUT$SAMPID]
TCGA_STAD <- TCGA_MOR_list[["TCGA-STAD"]][, colnames(TCGA_MOR_list[["TCGA-STAD"]]) %in% TCGA_LUT$SAMPID]
TCGA_COAD <- TCGA_MOR_list[["TCGA-COAD"]][, colnames(TCGA_MOR_list[["TCGA-COAD"]]) %in% TCGA_LUT$SAMPID]

cbind(c(sapply(RNA_MT_ensembl_list$unknownRNA_MTs$ensembl_gene_id, function(x){
  cor.test(TCGA_LUSC[x, ], TCGA_LUSC[NNMT_ensembl, ])$estimate
}), total = cor.test(colSums(TCGA_LUSC[RNA_MT_ensembl_list$unknownRNA_MTs$ensembl_gene_id, ]), TCGA_LUSC[NNMT_ensembl, ])$estimate),
c(sapply(RNA_MT_ensembl_list$unknownRNA_MTs$ensembl_gene_id, function(x){
  cor.test(TCGA_STAD[x, ], TCGA_STAD[NNMT_ensembl, ])$estimate
}), total = cor.test(colSums(TCGA_STAD[RNA_MT_ensembl_list$unknownRNA_MTs$ensembl_gene_id, ]), TCGA_STAD[NNMT_ensembl, ])$estimate),
c(sapply(RNA_MT_ensembl_list$unknownRNA_MTs$ensembl_gene_id, function(x){
  cor.test(TCGA_COAD[x, ], TCGA_COAD[NNMT_ensembl, ])$estimate
}), total = cor.test(colSums(TCGA_COAD[RNA_MT_ensembl_list$unknownRNA_MTs$ensembl_gene_id, ]), TCGA_COAD[NNMT_ensembl, ])$estimate)

)

# now for GTEX - Fig S7

GTEX_RNA_MT_PEMT_reciprocal_scores <- lapply(names(RNA_MT_ensembl_list), function(x){
  
  if(!file.exists(paste0("output/GTEX_", x, "_everything_object.rds"))){
    
    GTEX_thisone_everything_object <- analyse.by.subtype(database = "GTEX",
                                                         geneset_ensembl = RNA_MT_ensembl_list[[x]]$ensembl_gene_id,
                                                         genebasename = x)
    
    saveRDS(GTEX_thisone_everything_object, paste0("output/GTEX_", x, "_everything_object.rds"))
    
  } else {
    
    GTEX_thisone_everything_object <- readRDS(paste0("output/GTEX_", x, "_everything_object.rds"))
    
  }
  
  if(!file.exists(paste0("output/GTEX_PEMT_everything_object_", x, "added.rds"))){
    
    thisoneadded <- add.geneset.to.everything.subtype(database = "GTEX",
                                                      geneset_ensembl = RNA_MT_ensembl_list[[x]]$ensembl_gene_id,
                                                      genesetname = x,
                                                      everything_object_list = GTEX_PEMT_bytissue)
    
    saveRDS(thisoneadded, paste0("output/GTEX_PEMT_everything_object_", x, "added.rds"))
    
  } else {
    
    thisoneadded <- readRDS(paste0("output/GTEX_PEMT_everything_object_", x, "added.rds"))
    
  }
  
  tempscorelist <- lapply(1:length(GTEX_thisone_everything_object), function(i){
    
    GTEX_thisone_vec <- GTEX_thisone_everything_object[[i]]$allcorrs
    names(GTEX_thisone_vec) <- row.names(GTEX_thisone_everything_object[[i]])
    
    thisoneadded_vec <- thisoneadded[[i]]$allcorrs
    names(thisoneadded_vec) <- row.names(thisoneadded[[i]])
    
    calculate.relative.reciprocal.scores(everything_object1 = GTEX_thisone_vec,
                                         everything_object2 = thisoneadded_vec)
    
  })
  
  tempscore_df <- do.call(rbind, tempscorelist)
  row.names(tempscore_df) <- names(GTEX_thisone_everything_object)
  
  return(tempscore_df)
  
})

names(GTEX_RNA_MT_PEMT_reciprocal_scores) <- names(RNA_MT_ensembl_list)
saveRDS(GTEX_RNA_MT_PEMT_reciprocal_scores, "output/GTEX_RNA_MT_PEMT_reciprocal_scores.rds")
GTEX_RNA_MT_PEMT_reciprocal_scores <- readRDS("output/GTEX_RNA_MT_PEMT_reciprocal_scores.rds")

GTEX_RNA_MT_PEMT_reciprocal_scores <- lapply(names(GTEX_RNA_MT_PEMT_reciprocal_scores), function(thisone){

  x <- data.frame(GTEX_RNA_MT_PEMT_reciprocal_scores[[thisone]])
  
  x[, "MTs_in_PEMT_percent"] <- (unlist(x[, 1]) * 100)/60489
  x[, "PEMT_in_MTs_percent"] <- (unlist(x[, 2]) * 100)/60489
  
  colnames(x) <- str_replace(colnames(x), thisone, "MTset")
  
  x[x$PEMT_in_MTs_percent < 2.5 & x$MTs_in_PEMT_percent < 2.5, "tissue"] <- convert.GTEX.tissue.labels.for.plot(row.names(x[x$PEMT_in_MTs_percent < 2.5 & x$MTs_in_PEMT_percent < 2.5, ]), supershort = TRUE)
  x[x$PEMT_in_MTs_percent > 2.5 | x$MTs_in_PEMT_percent > 2.5, "tissue"] <- ""
  
  return(x)
  
})

names(GTEX_RNA_MT_PEMT_reciprocal_scores) <- names(RNA_MT_ensembl_list)


lapply(names(GTEX_RNA_MT_PEMT_reciprocal_scores), function(thisone){

  thisonelabel <- str_replace(thisone, "_", " ")
  
  thisdata <- GTEX_RNA_MT_PEMT_reciprocal_scores[[thisone]]
  thisdata[(thisdata$MTs_in_NNMT_percent > 2.5) | (thisdata$NNMT_in_MTs_percent > 2.5), "tissue"] <- ""
  
  pdf(paste0("graphics/GTEX-", thisone, "-PEMT-tissues-reciprocal_scores.pdf"),
      width = 1.75,
      height = 1.75)
  
  print(
    ggplot(aes(x = PEMT_in_MTs_percent, y = MTs_in_PEMT_percent, label = cancer), data = thisdata) + 
      geom_rect(aes(xmax = 2.5,
                    xmin = 0,
                    ymax = 2.5,
                    ymin = 0),
                fill = "white",
                linetype = "dashed",
                color = "grey") +
      geom_point(aes(size = 10/log10(reciprocal_score)),
                 colour = "red",
                 alpha = 0.7) + 
      scale_x_log10(breaks = c(0.01, 0.1, 1, 10, 100),
                    expand = c(0, 0),
                    limits = c(0.001, 100),
                    labels = label_number(drop0trailing = TRUE)) +
      scale_y_log10(breaks = c(0.01, 0.1, 1, 10, 100),
                    expand = c(0, 0),
                    limits = c(0.001, 100),
                    labels = label_number(drop0trailing = TRUE)) +
      geom_label_repel(size = 1.5,
                       max.overlaps = 40,
                       label.size = 0.1,
                       box.padding = 0.1,
                       label.padding = 0.1) + 
      theme_classic() +
      theme(legend.position = "none",
            axis.text.x = element_text(colour = "black",
                                       size = 10),
            axis.text.y = element_text(colour = "black",
                                       size = 10),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            plot.margin = margin(t = 8,
                                 r = 8,
                                 b = 8,
                                 l = 8)) + 
      xlab(substitute(italic(NNMT)~"%"^"ile"~"in"~thisonelabel~"correlations")) +
      ylab(substitute(thisonelabel~"%"^"ile"~"in"~italic(NNMT)~"correlations"))
    
  )
  
  dev.off()
  
})

GTEX_RNA_MT_PEMT_reciprocal_scores_for_facet <- lapply(names(GTEX_RNA_MT_PEMT_reciprocal_scores), function(thisone){

  x <- GTEX_RNA_MT_PEMT_reciprocal_scores[[thisone]]
  
  x[, "MTset"] <- thisone
  
  return(x)
  
})

GTEX_RNA_MT_PEMT_reciprocal_scores_for_facet <- do.call(rbind, GTEX_RNA_MT_PEMT_reciprocal_scores_for_facet)

write.table(GTEX_RNA_MT_PEMT_reciprocal_scores_for_facet,
            file = "plot_data/Fig S7/Fig_S7_bubble_plotdata.txt",
            row.names = TRUE,
            col.names = TRUE,
            sep = "\t")

pdf("graphics/GTEX_PEMT_otherMTs_facet.pdf")

ggplot(aes(x = PEMT_in_MTs_percent, y = MTs_in_PEMT_percent, label = tissue), data = GTEX_RNA_MT_PEMT_reciprocal_scores_for_facet) + 
  geom_rect(aes(xmax = 2.5,
                xmin = 0,
                ymax = 2.5,
                ymin = 0),
            fill = "white",
            linetype = "dashed",
            color = "grey") +
  geom_point(aes(size = 10/log10(reciprocal_score)),
             colour = "red",
             alpha = 0.7) + 
  scale_x_log10(breaks = c(0.01, 0.1, 1, 10, 100),
                expand = c(0, 0),
                limits = c(0.001, 100),
                labels = label_number(drop0trailing = TRUE)) +
  scale_y_log10(breaks = c(0.01, 0.1, 1, 10, 100),
                expand = c(0, 0),
                limits = c(0.001, 100),
                labels = label_number(drop0trailing = TRUE)) +
  geom_label_repel(size = 1.5,
                   max.overlaps = 40,
                   label.size = 0.1,
                   box.padding = 0.1,
                   label.padding = 0.1) + 
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(colour = "black",
                                   size = 10),
        axis.text.y = element_text(colour = "black",
                                   size = 10),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.line = element_line(),
        plot.margin = margin(t = 8,
                             r = 8,
                             b = 8,
                             l = 8)) + 
  xlab(substitute(italic(NNMT)~"%"^"ile"~"in"~thisonelabel~"correlations")) +
  ylab(substitute(thisonelabel~"%"^"ile"~"in"~italic(NNMT)~"correlations")) + 
  facet_wrap(~ MTset, scales = "free")

dev.off()

pdf("graphics/GTEX_PEMT_otherMTs_facet_nostrip.pdf",
    height = 4,
    width = 4)

ggplot(aes(x = PEMT_in_MTs_percent, y = MTs_in_PEMT_percent, label = tissue), data = GTEX_RNA_MT_PEMT_reciprocal_scores_for_facet) + 
  geom_rect(aes(xmax = 2.5,
                xmin = 0,
                ymax = 2.5,
                ymin = 0),
            fill = "white",
            linetype = "dashed",
            color = "grey") +
  geom_point(aes(size = 10/log10(reciprocal_score)),
             colour = "red",
             alpha = 0.7) + 
  scale_x_log10(breaks = c(0.01, 0.1, 1, 10, 100),
                expand = c(0, 0),
                limits = c(0.001, 100),
                labels = label_number(drop0trailing = TRUE)) +
  scale_y_log10(breaks = c(0.01, 0.1, 1, 10, 100),
                expand = c(0, 0),
                limits = c(0.001, 100),
                labels = label_number(drop0trailing = TRUE)) +
  geom_label_repel(size = 1.5,
                   max.overlaps = 40,
                   label.size = 0.1,
                   box.padding = 0.1,
                   label.padding = 0.1) + 
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(colour = "black",
                                   size = 8,
                                   angle = 30),
        axis.text.y = element_text(colour = "black",
                                   size = 8),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.line = element_line(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.margin = margin(t = 8,
                             r = 8,
                             b = 8,
                             l = 8)) + 
  xlab(substitute(italic(NNMT)~"%"^"ile"~"in"~thisonelabel~"correlations")) +
  ylab(substitute(thisonelabel~"%"^"ile"~"in"~italic(NNMT)~"correlations")) + 
  facet_wrap(~ MTset, scales = "free")

dev.off()

GTEX_RNA_hockeystick <- lapply(names(RNA_MT_samples_list), function(thisname){
  
  stepwise_PEMT_MT <- apply(RNA_MT_samples_list[[thisname]], 2, FUN = function(genesamples){
    
    stepwise.samples(database = "GTEX",
                     ensembl_gene = PEMT_ensembl,
                     stepwise_set = genesamples)
    
  })
  
  saveRDS(stepwise_PEMT_MT, paste0("output/GTEX_", thisname, "_stepwise_samples_list.rds"))
  
})

outputfiles <- list.files("output")
GTExhockeystickfiles <- c(outputfiles[str_detect(outputfiles, "_MTs_stepwise_")][str_detect(outputfiles[str_detect(outputfiles, "_MTs_stepwise_")], "GTEX")], "GTEX_PRMTs_stepwise_samples_list.rds")

for(i in 1:length(GTExhockeystickfiles)){
  
  temphock <- readRDS(paste0("output/", GTExhockeystickfiles[i]))
  
  colnames(temphock) <- paste0("permutation_", 1:ncol(temphock))
  
  write.table(temphock,
              file = paste0("plot_data/Fig S7/Fig_S7_", str_remove(hockeystickfiles[i], "_stepwise_samples_list.rds"), "_hockeystick_plotdata.txt"),
              row.names = FALSE,
              col.names = TRUE,
              sep = "\t")
  
}

# plot them out

lapply(GTExhockeystickfiles, function(x){
  
  tempdata <- readRDS(paste0("output/",x))
  tempname <- str_remove(x, "_stepwise_samples_list.rds")
  
  tempmelt <- melt(tempdata)
  
  pdf(paste0("graphics/GTEX_stepwise_", tempname, ".pdf"), 
      width = 2.5, 
      height = 2.5)
  
  print(
    ggplot(aes(x = Var1, y = value), data = tempmelt) + 
      geom_line(aes(group = Var2), alpha = 0.03) + 
      geom_smooth(col = "red", se = FALSE) + 
      theme_classic() + 
      theme(axis.text.x = element_text(colour = "black",
                                       size = 10),
            axis.text.y = element_text(colour = "black",
                                       size = 10),
            axis.title.x = element_text(size = 10),
            axis.title.y = element_text(size = 10)) + 
      xlab("Number of pooled MTs") + 
      ylab(substitute("Spearman's"~rho~"of MTs to"~paste(italic("NNMT")))) +
      geom_hline(yintercept = 0, linetype = "dashed", col = "grey")
  )
  
  dev.off()
  
})

#### PROLIFERATIVE INDEX Fig 2L Fig 2M ####

# PROLIFERATIVE INDEX FOR GTEX

# DESeq2 wants a colData object. Not actually used for the normalisation. Here we can use the sample IDs with tissue.
col_tissue <- LUT[colnames(countsdata), "Tissue"]
col_tissue <- as.matrix(col_tissue)

rownames(col_tissue) <- colnames(countsdata)
colnames(col_tissue) <- c("Tissue")

# need to create a DESeq2 object. Design set to ~1 allows for use of estimateSizeFactors
tempdds <- DESeqDataSetFromMatrix(countData = countsdata, colData = col_tissue, design = ~ 1)

tempvsd <- varianceStabilizingTransformation(tempdds, blind = TRUE)

# retrieve ensembl gene id and HGNC gene symbols for genes on mitochondrial DNA from biomart
replace_gene_symbols <- getBM(mart = ensembl,
                              attributes = c("hgnc_symbol", "ensembl_gene_id"),
                              filters = "ensembl_gene_id",
                              values = row.names(tempvsd))

replace_gene_symbols <- replace_gene_symbols[!(duplicated(replace_gene_symbols$ensembl_gene_id)), ]

row.names(replace_gene_symbols) <- replace_gene_symbols$ensembl_gene_id

MOR_vst_forPI <- data.frame(assay(tempvsd[row.names(tempvsd) %in% replace_gene_symbols$ensembl_gene_id, ]))
names <- replace_gene_symbols[row.names(MOR_vst_forPI), "hgnc_symbol"]
names[is.na(names)] <- paste("na", 1:sum(is.na(names)))
names[names == ""] <- paste(1:sum(names == ""))
names[duplicated(names)] <- paste0(names[duplicated(names)], 1:length(names[duplicated(names)]))
row.names(MOR_vst_forPI) <- names

GTEX_PIobject <- readDataForPI(MOR_vst_forPI, modelIDs = c("SCYL3"))
GTEX_PIs <- calculatePI(GTEX_PIobject)

GTEX_PIs_long <- data.frame(matrix(ncol = 2, nrow = 0))
colnames(GTEX_PIs_long) <- c("PI", "tissue")

for(i in 1:length(tissues)){
  
  tempsampIDs <- LUT[LUT$Tissue == tissues[i], "SAMPID"]
  
  tissuelabelvec <- vector(length = length(tempsampIDs))
  tissuelabelvec[] <- tissues[i]
  
  temp_df <- cbind(GTEX_PIs[tempsampIDs], tissuelabelvec)
  colnames(temp_df) <- c("PI", "tissue")
  
  GTEX_PIs_long <- rbind(GTEX_PIs_long, temp_df)
  
}

GTEX_PIs_long$PI <- as.numeric(GTEX_PIs_long$PI)
saveRDS(GTEX_PIs_long, "output/GTEX_PIs_long.rds")
# GTEX_PIs_long <- readRDS("output/GTEX_PIs_long.rds")

# TCGA PROLIFERATIVE INDEX

# TCGA_counts_list <- readRDS("output/TCGA_counts_list.rds")
# TCGA_LUT <- readRDS("output/TCGA_LUT.rds")

TCGA_counts_df <- do.call(cbind, TCGA_counts_list)

TCGA_counts_df <- TCGA_counts_df[, colnames(TCGA_counts_df) %in% TCGA_LUT$SAMPID]
col_data <- TCGA_LUT[colnames(TCGA_counts_df), "cancer"]
names(col_data) <- TCGA_LUT[colnames(TCGA_counts_df), "SAMPID"]
col_data <- data.frame(col_data)

# need to create a DESeq2 object. Design set to ~1 allows for use of estimateSizeFactors
TCGAdds <- DESeqDataSetFromMatrix(countData = TCGA_counts_df, colData = col_data, design = ~ 1)

TCGAvsd <- varianceStabilizingTransformation(TCGAdds, blind = TRUE)

TCGA_vst_df <- data.frame(assay(TCGAvsd))
# TCGA_vst_df <- readRDS("output/TCGA_vst_df.rds")

# retrieve ensembl gene id and HGNC gene symbols for genes on mitochondrial DNA from biomart
TCGA_replace_gene_symbols <- getBM(mart = ensembl,
                                   attributes = c("hgnc_symbol", "ensembl_gene_id"),
                                   filters = "ensembl_gene_id",
                                   values = row.names(TCGA_vst_df))

TCGA_replace_gene_symbols <- TCGA_replace_gene_symbols[!(duplicated(TCGA_replace_gene_symbols$ensembl_gene_id)), ]
row.names(TCGA_replace_gene_symbols) <- TCGA_replace_gene_symbols$ensembl_gene_id

TCGA_vst_df <- TCGA_vst_df[row.names(TCGA_vst_df) %in% TCGA_replace_gene_symbols$ensembl_gene_id, ]
names <- TCGA_replace_gene_symbols[row.names(TCGA_vst_df), "hgnc_symbol"]
names[is.na(names)] <- paste("na.", 1:sum(is.na(names)))
names[names == ""] <- paste(1:sum(names == ""))
names[duplicated(names)] <- paste0(names[duplicated(names)], 1:length(names[duplicated(names)]))
row.names(TCGA_vst_df) <- names

TCGA_PIobject <- readDataForPI(TCGA_vst_df, modelIDs = c("SCYL3"))
TCGA_PI <- calculatePI(TCGA_PIobject)

TCGAforbind <- TCGA_MOR_mtOXPHOS_nuOXPHOS$cancer_mean_corr_df$median_mtnuOXPHOS_spearman
names(TCGAforbind) <- TCGA_MOR_mtOXPHOS_nuOXPHOS$cancer_mean_corr_df$cancertype

TCGA_PIs_long <- data.frame(matrix(ncol = 2, nrow = 0))
colnames(TCGA_PIs_long) <- c("PI", "cancer")

for(i in 1:length(TCGA_MOR_mtOXPHOS_nuOXPHOS$cancer_mean_corr_df$cancertype)){
  
  tempsampIDs <- TCGA_LUT[TCGA_LUT$cancer == TCGA_MOR_mtOXPHOS_nuOXPHOS$cancer_mean_corr_df$cancertype[i], "SAMPID"]
  presentsampIDs <- tempsampIDs[tempsampIDs %in% str_replace_all(names(TCGA_PI), pattern = "\\.", replacement = "-")]
  
  cancerlabelvec <- vector(length = length(presentsampIDs))
  cancerlabelvec[] <- TCGA_MOR_mtOXPHOS_nuOXPHOS$cancer_mean_corr_df$cancertype[i]
  
  temp_df <- cbind(TCGA_PI[str_replace_all(presentsampIDs, pattern = "-", replacement = "\\.")], cancerlabelvec)
  colnames(temp_df) <- c("PI", "cancer")
  
  TCGA_PIs_long <- rbind(TCGA_PIs_long, temp_df)
  
}

TCGA_PIs_long$PI <- as.numeric(TCGA_PIs_long$PI)
saveRDS(TCGA_PIs_long, "output/TCGA_PIs_long.rds")
# TCGA_PIs_long <- readRDS("output/TCGA_PIs_long.rds")

GTEX_PIs_long[, "reciprocal_score"] <- GTEX_tissue_HMT_PEMT_reciprocal_scores_df[match(GTEX_PIs_long$tissue, row.names(GTEX_tissue_HMT_PEMT_reciprocal_scores_df)), "reciprocal_score"] 
TCGA_PIs_long[, "reciprocal_score"] <- TCGA_cancer_reciprocal_scores_df[match(TCGA_PIs_long$cancer, row.names(TCGA_cancer_reciprocal_scores_df)), "reciprocal_score"] 

cor.test(GTEX_PIs_long$PI, log10(GTEX_PIs_long$reciprocal_score))
cor.test(as.numeric(TCGA_PIs_long$PI), log10(TCGA_PIs_long$reciprocal_score))$p.value

write.table(GTEX_PIs_long[, c(1:2, 4)],
            file = "plot_data/Fig 2/Fig_2L_plotdata.txt",
            row.names = TRUE,
            col.names = TRUE,
            sep = "\t")

write.table(TCGA_PIs_long[, c(1:2, 4)],
            file = "plot_data/Fig 2/Fig_2M_plotdata.txt",
            row.names = TRUE,
            col.names = TRUE,
            sep = "\t")

pdf("graphics/GTEX_tissue_PI_vs_reciprocalscore_violin.pdf",
    height = 3,
    width = 4)

ggplot(data = GTEX_PIs_long, aes(y = PI, x = reciprocal_score)) +
  geom_violin(aes(group = tissue, fill = tissue), 
              show.legend = FALSE, 
              scale = 1, 
              width = 0.3, 
              position = "identity", 
              lwd = 0.2) +
  scale_x_continuous(trans = reverselog_trans(10), 
                     breaks = c(10, 1000, 100000, 10000000, 1000000000),
                     expand = c(0,0),
                     limits = c(1000000000, 10),
                     labels = trans_format("log10", math_format(10^.x))) + 
  scale_y_continuous() +
  theme_classic() +
  theme(axis.text.x = element_text(colour = "black",
                                   size = 10),
        axis.text.y = element_text(colour = "black",
                                   size = 10),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10,
                                    hjust = 0.4)) + 
  ylab("Proliferative Index (AU)") + 
  xlab(substitute(italic(PEMT)~"- HMT reciprocal score")) + 
  geom_smooth(method = "lm",
              se = FALSE) +
  annotate(geom = "text",
           x = 1000,
           y = 11,
           label = substitute(R^2~"= 0.116")) +
  annotate(geom = "text",
           x = 1000,
           y = 10,
           label = substitute(italic(p)~"= 0"))

dev.off()

pdf("graphics/TCGA_cancer_PI_vs_reciprocalscore_violin.pdf",
    height = 3,
    width = 4)

ggplot(data = TCGA_PIs_long, aes(y = PI, x = reciprocal_score)) +
  geom_violin(aes(group = cancer, fill = cancer), 
              show.legend = FALSE, 
              scale = 1, 
              width = 0.3, 
              position = "identity", 
              lwd = 0.2) +
  scale_x_continuous(trans = reverselog_trans(10), 
                     breaks = c(10, 1000, 100000, 10000000, 1000000000),
                     expand = c(0,0),
                     limits = c(1000000000, 10),
                     labels = trans_format("log10", math_format(10^.x))) + 
  scale_y_continuous() +
  theme_classic() +
  theme(axis.text.x = element_text(colour = "black",
                                   size = 10),
        axis.text.y = element_text(colour = "black",
                                   size = 10),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10,
                                    hjust = 0.4)) + 
  ylab("Proliferative Index (AU)") + 
  xlab(substitute(italic(NNMT)~"-HMT reciprocal score")) + 
  geom_smooth(method = "lm",
              se = FALSE) +
  annotate(geom = "text",
           x = 1000,
           y = 7.5,
           label = substitute(R^2~"= 0.038")) +
  annotate(geom = "text",
           x = 1000,
           y = 7,
           label = substitute(italic(p)~"= 1.05 x"~10^-84))

dev.off()

#### EPIC TIMER DECONVOLUTIONS Fig S3 ####

# this file downloaded from here [http://timer.cistrome.org/]
tcga_infiltration <- read.table("input/infiltration_estimation_for_tcga.csv", sep = ",", header = TRUE, quote = "")

# 14 is LAML
TIMERinfiltration_NNMT_cor <- sapply(TCGA_MOR_list[-14], function(x){
  
  samples <- colnames(x)[colnames(x) %in% TCGA_LUT$SAMPID]
  samples <- samples[str_extract(samples, pattern = "TCGA-[0-9A-Z]{2}-[0-9A-Z]{4}-[0-9A-Z]{2}") %in% tcga_infiltration$cell_type]
  
  if(length(samples) == 0){
    return()
  } else {
    
    NNMT_vec <- x["ENSG00000166741", samples]
    
    infiltration_matrix <- tcga_infiltration[match(str_extract(samples, pattern = "TCGA-[0-9A-Z]{2}-[0-9A-Z]{4}-[0-9A-Z]{2}"), tcga_infiltration$cell_type), str_detect(colnames(tcga_infiltration), pattern = "TIMER$")]
    
    if(any(is.na(infiltration_matrix))){
      return()
    }else{
      
      
      apply(infiltration_matrix, 2, function(y){
        cor.test(NNMT_vec, y, method = "spearman")$estimate
      })
      
    }
  }
  
})

TIMERinfiltration_HMT_cor <- sapply(TCGA_MOR_list[-14], function(x){
  
  samples <- colnames(x)[colnames(x) %in% TCGA_LUT$SAMPID]
  samples <- samples[str_extract(samples, pattern = "TCGA-[0-9A-Z]{2}-[0-9A-Z]{4}-[0-9A-Z]{2}") %in% tcga_infiltration$cell_type]
  
  if(length(samples) == 0){
    return()
  } else {
    
    HMT_vec <- colSums(x[SET_HMTs$ensembl_gene_id, samples])
    
    infiltration_matrix <- tcga_infiltration[match(str_extract(samples, pattern = "TCGA-[0-9A-Z]{2}-[0-9A-Z]{4}-[0-9A-Z]{2}"), tcga_infiltration$cell_type), str_detect(colnames(tcga_infiltration), pattern = "TIMER$")]
    
    if(any(is.na(infiltration_matrix))){
      return()
    }else{
      
      
      apply(infiltration_matrix, 2, function(y){
        cor.test(HMT_vec, y, method = "spearman")$estimate
      })
      
    }
  }
  
})

TIMERinfiltration_NNMTcor_melt <- melt(TIMERinfiltration_NNMT_cor)
TIMERinfiltration_NNMTcor_melt$Var1 <- str_remove(TIMERinfiltration_NNMTcor_melt$Var1, "_TIMER$")
TIMERinfiltration_NNMTcor_melt$Var1 <- str_replace(TIMERinfiltration_NNMTcor_melt$Var1, "\\.$", "+")
TIMERinfiltration_NNMTcor_melt$Var1 <- str_replace_all(TIMERinfiltration_NNMTcor_melt$Var1, "\\.", " ")

TIMERinfiltration_NNMTcor_celltypes <- TIMERinfiltration_NNMTcor_melt %>% group_by(Var1) %>% summarize(mean = mean(value))

TIMERinfiltration_NNMTcor_melt$Var1 <- factor(TIMERinfiltration_NNMTcor_melt$Var1, levels = unlist(TIMERinfiltration_NNMTcor_celltypes[order(TIMERinfiltration_NNMTcor_celltypes$mean), "Var1"]))

TIMERinfiltration_HMTcor_melt <- melt(TIMERinfiltration_HMT_cor)
TIMERinfiltration_HMTcor_melt$Var1 <- str_remove(TIMERinfiltration_HMTcor_melt$Var1, "_TIMER$")
TIMERinfiltration_HMTcor_melt$Var1 <- str_replace(TIMERinfiltration_HMTcor_melt$Var1, "\\.$", "+")
TIMERinfiltration_HMTcor_melt$Var1 <- str_replace_all(TIMERinfiltration_HMTcor_melt$Var1, "\\.", " ")

TIMERinfiltration_HMTcor_celltypes <- TIMERinfiltration_HMTcor_melt %>% group_by(Var1) %>% summarize(mean = mean(value))

TIMERinfiltration_HMTcor_melt$Var1 <- factor(str_remove(TIMERinfiltration_HMTcor_melt$Var1, "_TIMER$"), levels = unlist(TIMERinfiltration_NNMTcor_celltypes[order(TIMERinfiltration_NNMTcor_celltypes$mean), "Var1"]))

TIMERinfiltration_NNMTcor_melt_copy <- TIMERinfiltration_NNMTcor_melt
colnames(TIMERinfiltration_NNMTcor_melt_copy) <- c("Cell_type", "cancer", "correlation")

TIMERinfiltration_HMTcor_melt_copy <- TIMERinfiltration_HMTcor_melt
colnames(TIMERinfiltration_HMTcor_melt_copy) <- c("Cell_type", "cancer", "correlation")

write.table(TIMERinfiltration_NNMTcor_melt_copy,
            file = "plot_data/Fig S3/Fig_S3A_plotdata.txt",
            row.names = FALSE,
            col.names = TRUE,
            sep = "\t")

write.table(TIMERinfiltration_HMTcor_melt_copy,
            file = "plot_data/Fig S3/Fig_S3C_plotdata.txt",
            row.names = FALSE,
            col.names = TRUE,
            sep = "\t")

pdf("graphics/TIMER_NNMT_correlations.pdf",
    height = 2.5,
    width = 3.5)

ggplot(data = TIMERinfiltration_NNMTcor_melt, aes(x = Var1, y = value, fill = Var1)) + 
  geom_boxplot(outlier.size = 0.5,
               fill = c("salmon", "goldenrod2", "darkolivegreen3", "dark turquoise", "dodger blue", "magenta")) + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45,
                                   size = 8, 
                                   color = "black",
                                   hjust = 1),
        axis.text.y = element_text(size = 10,
                                   color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10,
                                    hjust = 1),
        legend.position = "none") + 
  scale_y_continuous(breaks = seq(-0.6, 0.6, 0.3),
                     expand = c(0,0),
                     limits = c(-0.6, 0.6)) +
  geom_hline(yintercept = 0,
             linetype = "dashed",
             color = "grey") +
  ylab(substitute("Cancer type infiltration -"~italic(NNMT)~rho))

dev.off()

pdf("graphics/TIMER_HMT_correlations.pdf",
    height = 2.5,
    width = 3.5)

ggplot(data = TIMERinfiltration_HMTcor_melt, aes(x = Var1, y = value, fill = Var1)) + 
  geom_boxplot(outlier.size = 0.5,
               fill = c("salmon", "goldenrod2", "darkolivegreen3", "dark turquoise", "dodger blue", "magenta")) + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45,
                                   size = 8, 
                                   color = "black",
                                   hjust = 1),
        axis.text.y = element_text(size = 10,
                                   color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10,
                                    hjust = 1),
        legend.position = "none") + 
  scale_y_continuous(breaks = seq(-0.6, 0.6, 0.3),
                     expand = c(0,0),
                     limits = c(-0.6, 0.6)) +
  geom_hline(yintercept = 0,
             linetype = "dashed",
             color = "grey") +
  ylab(substitute("Cancer type infiltration - HMT"~rho))

dev.off()

rowMeans(infiltration_HMT_cor)
rowMeans(infiltration_NNMT_cor)

EPICinfiltration_NNMT_cor <- sapply(TCGA_MOR_list[-14], function(x){
  
  samples <- colnames(x)[colnames(x) %in% TCGA_LUT$SAMPID]
  samples <- samples[str_extract(samples, pattern = "TCGA-[0-9A-Z]{2}-[0-9A-Z]{4}-[0-9A-Z]{2}") %in% tcga_infiltration$cell_type]
  
  if(length(samples) == 0){
    return()
  } else {
    
    NNMT_vec <- x[NNMT_ensembl, samples]
    
    infiltration_matrix <- tcga_infiltration[match(str_extract(samples, pattern = "TCGA-[0-9A-Z]{2}-[0-9A-Z]{4}-[0-9A-Z]{2}"), tcga_infiltration$cell_type), str_detect(colnames(tcga_infiltration), pattern = "EPIC$")]
    infiltration_matrix <- infiltration_matrix[, !str_detect(colnames(infiltration_matrix), "uncharacterized")]
    
    if(any(is.na(infiltration_matrix))){
      return()
    }else{
      
      
      apply(infiltration_matrix, 2, function(y){
        cor.test(NNMT_vec, y, method = "spearman")$estimate
      })
      
    }
  }
  
})

EPICinfiltration_HMT_cor <- sapply(TCGA_MOR_list[-14], function(x){

  samples <- colnames(x)[colnames(x) %in% TCGA_LUT$SAMPID]
  samples <- samples[str_extract(samples, pattern = "TCGA-[0-9A-Z]{2}-[0-9A-Z]{4}-[0-9A-Z]{2}") %in% tcga_infiltration$cell_type]
  
  if(length(samples) == 0){
    return()
  } else {
    
    HMT_vec <- colSums(x[SET_HMTs$ensembl_gene_id, samples])
    
    infiltration_matrix <- tcga_infiltration[match(str_extract(samples, pattern = "TCGA-[0-9A-Z]{2}-[0-9A-Z]{4}-[0-9A-Z]{2}"), tcga_infiltration$cell_type), str_detect(colnames(tcga_infiltration), pattern = "EPIC$")]
    
    infiltration_matrix <- infiltration_matrix[, !str_detect(colnames(infiltration_matrix), "uncharacterized")]
    
    if(any(is.na(infiltration_matrix))){
      return()
    }else{
      
      
      apply(infiltration_matrix, 2, function(y){
        cor.test(HMT_vec, y, method = "spearman")$estimate
      })
      
    }
  }
  
})

EPICinfiltration_NNMTcor_melt <- melt(EPICinfiltration_NNMT_cor)
EPICinfiltration_NNMTcor_melt$Var1 <- str_remove(EPICinfiltration_NNMTcor_melt$Var1, "_EPIC$")
EPICinfiltration_NNMTcor_melt$Var1 <- str_replace(EPICinfiltration_NNMTcor_melt$Var1, "\\.$", "+")
EPICinfiltration_NNMTcor_melt$Var1 <- str_replace_all(EPICinfiltration_NNMTcor_melt$Var1, "\\.", " ")
EPICinfiltration_NNMTcor_melt$Var1 <- str_replace_all(EPICinfiltration_NNMTcor_melt$Var1, "Cancer associated fibroblast", "CAFs")

EPICinfiltration_NNMTcor_celltypes <- EPICinfiltration_NNMTcor_melt %>% group_by(Var1) %>% summarize(mean = mean(value))

EPICinfiltration_NNMTcor_melt$Var1 <- factor(EPICinfiltration_NNMTcor_melt$Var1, levels = unlist(EPICinfiltration_NNMTcor_celltypes[order(EPICinfiltration_NNMTcor_celltypes$mean), "Var1"]))

EPICinfiltration_HMTcor_melt <- melt(EPICinfiltration_HMT_cor)
EPICinfiltration_HMTcor_melt$Var1 <- str_remove(EPICinfiltration_HMTcor_melt$Var1, "_EPIC$")
EPICinfiltration_HMTcor_melt$Var1 <- str_replace(EPICinfiltration_HMTcor_melt$Var1, "\\.$", "+")
EPICinfiltration_HMTcor_melt$Var1 <- str_replace_all(EPICinfiltration_HMTcor_melt$Var1, "\\.", " ")
EPICinfiltration_HMTcor_melt$Var1 <- str_replace_all(EPICinfiltration_HMTcor_melt$Var1, "Cancer associated fibroblast", "CAFs")

EPICinfiltration_HMTcor_celltypes <- EPICinfiltration_HMTcor_melt %>% group_by(Var1) %>% summarize(mean = mean(value))

EPICinfiltration_HMTcor_melt$Var1 <- factor(EPICinfiltration_HMTcor_melt$Var1, levels = unlist(EPICinfiltration_NNMTcor_celltypes[order(EPICinfiltration_NNMTcor_celltypes$mean), "Var1"]))

EPICinfiltration_NNMTcor_melt_copy <- EPICinfiltration_NNMTcor_melt
colnames(EPICinfiltration_NNMTcor_melt_copy) <- c("Cell_type", "cancer", "correlation")

EPICinfiltration_HMTcor_melt_copy <- EPICinfiltration_HMTcor_melt
colnames(EPICinfiltration_HMTcor_melt_copy) <- c("Cell_type", "cancer", "correlation")

write.table(EPICinfiltration_NNMTcor_melt_copy,
            file = "plot_data/Fig S3/Fig_S3B_plotdata.txt",
            row.names = FALSE,
            col.names = TRUE,
            sep = "\t")

write.table(EPICinfiltration_HMTcor_melt_copy,
            file = "plot_data/Fig S3/Fig_S3D_plotdata.txt",
            row.names = FALSE,
            col.names = TRUE,
            sep = "\t")

pdf("graphics/EPIC_NNMT_correlations.pdf",
    height = 2.5,
    width = 3.5)

ggplot(data = EPICinfiltration_NNMTcor_melt, aes(x = Var1, y = value, fill = Var1)) + 
  geom_boxplot(outlier.size = 0.5,
               fill = c("goldenrod2", "darkolivegreen3", "salmon", "purple", "red", "dark turquoise", "grey")) + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45,
                                   size = 8, 
                                   color = "black",
                                   hjust = 1),
        axis.text.y = element_text(size = 10,
                                   color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10,
                                    hjust = 1),
        legend.position = "none") + 
  scale_y_continuous(breaks = seq(-0.6, 0.6, 0.3),
                     expand = c(0,0),
                     limits = c(-0.6, 0.6)) +
  geom_hline(yintercept = 0,
             linetype = "dashed",
             color = "grey") +
  ylab(substitute("Cancer type infiltration -"~italic(NNMT)~rho))

dev.off()

pdf("graphics/EPIC_HMT_correlations.pdf",
    height = 2.5,
    width = 3.5)

ggplot(data = EPICinfiltration_HMTcor_melt, aes(x = Var1, y = value, fill = Var1)) + 
  geom_boxplot(outlier.size = 0.5,
               fill = c("goldenrod2", "darkolivegreen3", "salmon", "purple", "red", "dark turquoise", "grey")) + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45,
                                   size = 8, 
                                   color = "black",
                                   hjust = 1),
        axis.text.y = element_text(size = 10,
                                   color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10,
                                    hjust = 1),
        legend.position = "none") + 
  scale_y_continuous(breaks = seq(-0.6, 0.6, 0.3),
                     expand = c(0,0),
                     limits = c(-0.6, 0.6)) +
  geom_hline(yintercept = 0,
             linetype = "dashed",
             color = "grey") +
  ylab(substitute("Cancer type infiltration - HMT"~rho))

dev.off()

NNMT_SET_immune <- sapply(TCGA_MOR_list[-14], function(x){
  
  samples <- colnames(x)[colnames(x) %in% TCGA_LUT$SAMPID]
  samples <- samples[str_extract(samples, pattern = "TCGA-[0-9A-Z]{2}-[0-9A-Z]{4}-[0-9A-Z]{2}") %in% tcga_infiltration$cell_type]
  
  if(length(samples) == 0){
    return()
  } else {
    
    tempdata <- data.frame(rbind(NNMT_expr = as.numeric(x["ENSG00000166741", samples]),
                                 HMT = colSums(x[SET_HMTs$ensembl_gene_id, samples]),
                                 t(tcga_infiltration[match(str_extract(samples, pattern = "TCGA-[0-9A-Z]{2}-[0-9A-Z]{4}-[0-9A-Z]{2}"), tcga_infiltration$cell_type), str_detect(colnames(tcga_infiltration), pattern = "TIMER$")])))
    
    NNMT_residuals <- lm(formula = unlist(tempdata["NNMT_expr", ]) ~ unlist(tempdata["B.cell_TIMER",]) + unlist(tempdata["T.cell.CD4._TIMER",]) + unlist(tempdata["T.cell.CD8._TIMER",]) + unlist(tempdata["Neutrophil_TIMER",]) + unlist(tempdata["Macrophage_TIMER",]) + unlist(tempdata["Myeloid.dendritic.cell_TIMER",]))$residuals
    
    HMT_residuals <- lm(formula = unlist(tempdata["HMT", ]) ~ unlist(tempdata["B.cell_TIMER",]) + unlist(tempdata["T.cell.CD4._TIMER",]) + unlist(tempdata["T.cell.CD8._TIMER",]) + unlist(tempdata["Neutrophil_TIMER",]) + unlist(tempdata["Macrophage_TIMER",]) + unlist(tempdata["Myeloid.dendritic.cell_TIMER",]))$residuals
    
    cor.test(NNMT_residuals, HMT_residuals, method = "spearman")$estimate
    
  }
  
})

NNMT_SET_noimmune <- sapply(TCGA_MOR_list[-14], function(x){
  
  samples <- colnames(x)[colnames(x) %in% TCGA_LUT$SAMPID]
  samples <- samples[str_extract(samples, pattern = "TCGA-[0-9A-Z]{2}-[0-9A-Z]{4}-[0-9A-Z]{2}") %in% tcga_infiltration$cell_type]
  
  if(length(samples) == 0){
    return()
  } else {
    
    tempdata <- data.frame(rbind(NNMT_expr = as.numeric(x["ENSG00000166741", samples]),
                                 HMT = colSums(x[SET_HMTs$ensembl_gene_id, samples]),
                                 t(tcga_infiltration[match(str_extract(samples, pattern = "TCGA-[0-9A-Z]{2}-[0-9A-Z]{4}-[0-9A-Z]{2}"), tcga_infiltration$cell_type), str_detect(colnames(tcga_infiltration), pattern = "TIMER$")])))
    
    cor.test(unlist(tempdata["NNMT_expr",]), unlist(tempdata["HMT",]), method = "spearman")$estimate
    
  }
  
})

TIMER_correction_forplot <- data.frame(cbind(immune = NNMT_SET_immune, noimmune = NNMT_SET_noimmune))
TIMER_correction_forplot[, "cancer"] <- str_remove(str_remove(row.names(TIMER_correction_forplot), "\\.rho"), "TCGA-")

write.table(TIMER_correction_forplot,
            file = "plot_data/Fig S3/Fig_S3E_plotdata.txt",
            row.names = FALSE,
            col.names = TRUE,
            sep = "\t")

TIMER_correction_forplot[!(TIMER_correction_forplot$cancer %in% c("KICH", "ACC", "LGG")), "cancer"] <- ""

pdf("graphics/TIMER_scatterplot.pdf",
    height = 2.5,
    width = 2.5)

ggplot(data = TIMER_correction_forplot, aes(x = noimmune, y = immune, label = cancer)) + 
  geom_point() + 
  theme_classic() + 
  theme(axis.text.x = element_text(size = 10,
                                   colour = "black"),
        axis.text.y = element_text(size = 10,
                                   colour = "black"),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        plot.margin = margin(t = 6,
                             r = 10,
                             b = 6,
                             l = 6)) +
  geom_smooth(method = "lm", 
              se = FALSE) +
  geom_label_repel(size = 2) +
  xlab(substitute("HMT -"~italic(NNMT)~rho~"no correction")) + 
  ylab(substitute("HMT -"~italic(NNMT)~rho~"+ infiltration corr.")) +
  coord_cartesian(ylim = c(-1, 0.25),
                  xlim= c(-1, 0.25)) +
  geom_abline(slope = 1, linetype = "dashed") + 
  scale_y_continuous(breaks = seq(-1, 0.25, 0.25),
                     expand = c(0,0),
                     limits = c(-1, 0.25)) +
    scale_x_continuous(breaks = seq(-1, 0.25, 0.25),
                       expand = c(0,0),
                       limits = c(-1, 0.25)) +
  geom_hline(yintercept = 0,
             colour = "grey",
             linetype = "dashed") + 
  geom_vline(xintercept = 0,
             colour = "grey",
             linetype = "dashed") 

dev.off()

NNMT_SET_immune_EPIC <- sapply(TCGA_MOR_list[-14], function(x){

  samples <- colnames(x)[colnames(x) %in% TCGA_LUT$SAMPID]
  samples <- samples[str_extract(samples, pattern = "TCGA-[0-9A-Z]{2}-[0-9A-Z]{4}-[0-9A-Z]{2}") %in% tcga_infiltration$cell_type]
  
  if(length(samples) == 0){
    return()
  } else {
    
    EPIC_cells <- colnames(tcga_infiltration)[str_detect(colnames(tcga_infiltration), pattern = "EPIC$")]
    EPIC_cells <- EPIC_cells[!str_detect(EPIC_cells, "uncharacterized")]
    
    tempdata <- data.frame(cbind(NNMT_expr = as.numeric(x["ENSG00000166741", samples]),
                                 HMT = colSums(x[SET_HMTs$ensembl_gene_id, samples]),
                                 tcga_infiltration[match(str_extract(samples, pattern = "TCGA-[0-9A-Z]{2}-[0-9A-Z]{4}-[0-9A-Z]{2}"), tcga_infiltration$cell_type), str_detect(colnames(tcga_infiltration), pattern = "EPIC$")]))
    
    NNMT_residuals <- lm(formula = as.formula(paste("NNMT_expr", paste(EPIC_cells, collapse = " + "), sep = "~")), data = tempdata)$residuals
    
    HMT_residuals <- lm(formula = as.formula(paste("HMT", paste(EPIC_cells, collapse = " + "), sep = "~")), data = tempdata)$residuals
    
    cor.test(NNMT_residuals, HMT_residuals, method = "spearman")$estimate
    
  }
  
})

# ACC not robust result in EPIC, but LGG and KICH are

EPIC_correction_forplot <- data.frame(cbind(immune = NNMT_SET_immune_EPIC, noimmune = NNMT_SET_noimmune))
EPIC_correction_forplot[, "cancer"] <- str_remove(str_remove(row.names(EPIC_correction_forplot), "\\.rho"), "TCGA-")

write.table(EPIC_correction_forplot,
            file = "plot_data/Fig S3/Fig_3F_plotdata.txt",
            row.names = FALSE,
            col.names = TRUE,
            sep = "\t")

EPIC_correction_forplot[!(EPIC_correction_forplot$cancer %in% c("KICH", "LGG")), "cancer"] <- ""

pdf("graphics/EPIC_scatterplot.pdf",
    height = 2.5,
    width = 2.5)

ggplot(data = EPIC_correction_forplot, aes(x = noimmune, y = immune, label = cancer)) + 
  geom_point() + 
  theme_classic() + 
  theme(axis.text.x = element_text(size = 10,
                                   colour = "black"),
        axis.text.y = element_text(size = 10,
                                   colour = "black"),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        plot.margin = margin(t = 6,
                             r = 10,
                             b = 6,
                             l = 6)) +
  geom_smooth(method = "lm", 
              se = FALSE) +
  geom_label_repel(size = 2,
                   max.overlaps = 40) +
  xlab(substitute("HMT -"~italic(NNMT)~rho~"no correction")) + 
  ylab(substitute("HMT -"~italic(NNMT)~rho~"+ infiltration corr.")) +  
  coord_cartesian(ylim = c(-1, 0.25),
                  xlim= c(-1, 0.25)) +
  geom_abline(slope = 1, linetype = "dashed") + 
  scale_y_continuous(breaks = seq(-1, 0.25, 0.25),
                     expand = c(0,0),
                     limits = c(-1, 0.25)) +
  scale_x_continuous(breaks = seq(-1, 0.25, 0.25),
                     expand = c(0,0),
                     limits = c(-1, 0.25)) +
  geom_hline(yintercept = 0,
             colour = "grey",
             linetype = "dashed") + 
  geom_vline(xintercept = 0,
             colour = "grey",
             linetype = "dashed") 

dev.off()

#### ENCODE CHIPSEQ DATA Fig 3A ####

ENCODEtissuetypes <- c("esophagus_muscularis_mucosa",
                       "esophagus_squamous_epithelium",
                       "gastroesophageal_sphincter",
                       "sigmoid_colon",
                       "spleen")

# get file metadata

filesabove <- list.files("input")
metadatafile <- filesabove[str_detect(filesabove, "metadata")]
chipmetadatafile <- metadatafile[str_detect(metadatafile, "chip")]

chipmetadata_list <- lapply(chipmetadatafile, function(x){
  
  read.table(paste0("input/", x), quote = "", sep = "\t", header = TRUE, skip = 1)
  
})

chipmetadata_df <- do.call(rbind, chipmetadata_list)

RNAmetadatafiles <- lapply(ENCODEtissuetypes, function(thistissue){
  
  read.table(paste0("input/", thistissue, "/", thistissue, "_RNA_seq_metadata.tsv"), sep = "\t", header = TRUE, fill = TRUE)
  
})

RNAmetadata_df <- do.call(rbind, RNAmetadatafiles)

altRNAmetadatafiles <- metadatafile[str_detect(metadatafile, "RNA")]
altRNAmetadatafiles <- altRNAmetadatafiles[!str_detect(altRNAmetadatafiles, "NCI60")]

altRNAmetadata_list <- lapply(altRNAmetadatafiles, function(x){
  
  read.table(paste0("input/", x), quote = "", sep = "\t", header = TRUE, skip = 1)
  
})

altRNAmetadata_df <- do.call(rbind, altRNAmetadata_list)

RNAmetadata_df[, "biosample_summary"] <- altRNAmetadata_df[match(RNAmetadata_df$Experiment.accession, altRNAmetadata_df$Accession), "Biosample.summary"]

chipmetadata_df[, "Donor"] <- str_remove(str_remove(RNAmetadata_df[match(chipmetadata_df$Biosample.summary, RNAmetadata_df$biosample_summary), "Donor.s."], "/human-donors/"), "/$")

# get RNAseq data

process.RNA.data <- function(thistissue){
  
  tissue_file_list <- list.files(paste0("input/", thistissue))
  
  tissue_RNA_files <- tissue_file_list[str_detect(tissue_file_list, "tsv$")]
  tissue_RNA_files <- tissue_RNA_files[!str_detect(tissue_RNA_files, pattern = "metadata")]
  
  tissue_RNAseq_metadata <- read.table(paste0("input/", thistissue, "/", thistissue, "_RNA_seq_metadata.tsv"),
                                       sep = "\t",
                                       header = TRUE,
                                       quote = "",
                                       fill = TRUE)
  
  tissue_RNA_data <- lapply(tissue_RNA_files, function(x){
    read.table(paste0("input/", thistissue, "/", x),
               header = TRUE,
               sep = "\t",
               stringsAsFactors = FALSE,
               quote = "")
  })
  
  names(tissue_RNA_data) <- str_remove(str_remove(tissue_RNAseq_metadata[match(str_remove(tissue_RNA_files, pattern = "\\.tsv$"), tissue_RNAseq_metadata$File.accession), "Donor.s."], pattern = "\\/human-donors\\/"), pattern = "\\/$")
  
  countstable <- sapply(tissue_RNA_data, function(x){
    x$expected_count
  })
  
  row.names(countstable) <- str_remove(tissue_RNA_data[[1]]$gene_id, "\\..*$")
  mode(countstable) <- "integer"
  
  colnames(countstable) <- paste0(colnames(countstable), "-", thistissue)
  
  return(countstable)
  
}

RNAcounts <- lapply(ENCODEtissuetypes, process.RNA.data)
RNAcounts_df <- do.call(cbind, RNAcounts)

coldata <- data.frame(matrix(nrow = ncol(RNAcounts_df), ncol = 1))
coldata[] <- str_remove(str_extract(colnames(RNAcounts_df), pattern = "-.*$"), pattern = "^-")

tempdds <- DESeqDataSetFromMatrix(countData = RNAcounts_df, colData = coldata, design = ~ 1)

# this function estimates the scaling factors from the samples from the median of ratios wrt to the geometric mean for each gene across samples
tempdds <- estimateSizeFactors(tempdds)

# put the counts normalised by the scaling factors in a new object
normalised_counts <- counts(tempdds, normalized = TRUE)

saveRDS(normalised_counts, "output/ENCODE_normalised_counts.rds")
normalised_counts <- readRDS("output/ENCODE_normalised_counts.rds")

#### make beds for computematrix ####

genes <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene)

gene_bodies_df <- data.frame(seqnames = seqnames(genes),
                             starts = start(genes)-1,
                             ends = end(genes),
                             gene_id = genes$gene_id,
                             scores = c(rep(".", length(genes))),
                             strands = strand(genes))

gene_bodies_200bpflank_df <- gene_bodies_df
gene_bodies_200bpflank_df$starts <- gene_bodies_df$starts - 200
gene_bodies_200bpflank_df$ends <- gene_bodies_df$ends + 200

write.table(gene_bodies_df, file = "output/gene_bodies.bed", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(gene_bodies_200bpflank_df, file = "output/gene_bodies_200bpflanks.bed", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

gene_promoters <- promoters(genes)

promoters_df <- data.frame(seqnames = seqnames(gene_promoters),
                           starts = start(gene_promoters)-1,
                           ends = end(gene_promoters),
                           gene_id = gene_promoters$gene_id,
                           scores = c(rep(".", length(gene_promoters))),
                           strands = strand(gene_promoters))

write.table(promoters_df, file = "output/promoters.bed", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

# ENCODE BEDS FOR REPETITIVE ELEMENTS (FOR H3K9me3)

# limit repetitive elements to clearly annotated categories and convert to genomic ranges object
# repetitive elements file for hg38
rmask <- read.table("input/rmsk.txt", sep = "\t", header = FALSE, quote = "")
colnames(rmask)[c(6,7,8)] <- c("chr", "start", "end")
rmask[, "myID"] <- paste0(rmask$chr, rmask$start, rmask$end)
row.names(rmask) <- rmask$myID

rmask_gr <- makeGRangesFromDataFrame(rmask[rmask$V12 %in% c("Satellite", "LINE", "SINE", "LTR", "DNA", "Low_complexity", "snRNA", "tRNA", "Retroposon"), c("chr", "start", "end", "myID", "V11", "V12", "V13")], keep.extra.columns = TRUE)

ENCODE_H3K9_repeats_sites <- return.sites.over.peak.threshold(chiptarget = "H3K9me3",
                                                              sites_gr = rmask_gr,
                                                              chosenthresh = 0.75)

rmask_H3K9_select <- rmask_gr[ENCODE_H3K9_repeats_sites, ]

saveRDS(rmask_H3K9_select, "output/rmask_gr_ENCODE_H3K9_18of24peaks.rds")
# rmask_select <- readRDS("output/rmask_gr_18of24peaks.rds")

rmask_H3K9_select_df <- data.frame(seqnames = seqnames(rmask_H3K9_select),
                                   starts = start(rmask_H3K9_select)-1,
                                   ends = end(rmask_H3K9_select),
                                   gene_id = rmask_H3K9_select$myID,
                                   scores = c(rep(".", length(rmask_H3K9_select))),
                                   strands = strand(rmask_H3K9_select))

options(scipen=999)

write.table(rmask_H3K9_select_df, file = "output/rmask_ENCODE_H3K9_18of24peaks.bed", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

# get computematrix output

ENCODE_H3K27me3_promoter_sites <- return.sites.over.peak.threshold(chiptarget = "H3K27me3",
                                                                   sites_gr = promoters(genes),
                                                                   chosenthresh = 0.75)

ENCODE_H3K9me3_promoter_sites <- return.sites.over.peak.threshold(chiptarget = "H3K9me3",
                                                                  sites_gr = promoters(genes),
                                                                  chosenthresh = 0.75)

ENCODE_H3K36me3_promoter_sites <- return.sites.over.peak.threshold(chiptarget = "H3K36me3",
                                                                   sites_gr = promoters(genes),
                                                                   chosenthresh = 0.75)

ENCODE_H3K4me3_promoter_sites <- return.sites.over.peak.threshold(chiptarget = "H3K4me3",
                                                                  sites_gr = promoters(genes),
                                                                  chosenthresh = 1)

ENCODE_H3K27me3_genebody_sites <- return.sites.over.peak.threshold(chiptarget = "H3K27me3",
                                                                   sites_gr = genes,
                                                                   chosenthresh = 0.75)

ENCODE_H3K9me3_genebody_sites <- return.sites.over.peak.threshold(chiptarget = "H3K9me3",
                                                                  sites_gr = genes,
                                                                  chosenthresh = 0.75)

ENCODE_H3K36me3_genebody_sites <- return.sites.over.peak.threshold(chiptarget = "H3K36me3",
                                                                   sites_gr = genes,
                                                                   chosenthresh = 0.75)

ENCODE_H3K4me3_genebody_sites <- return.sites.over.peak.threshold(chiptarget = "H3K4me3",
                                                                  sites_gr = genes,
                                                                  chosenthresh = 1)

ENCODE_PEMT_H3K27_promoter_model <- fit.models.to.chip(computematrix_outputdirectory = "ENCODE_H3K27_computematrix",
                                                       selectedtag = "promoter",
                                                       removetag = "_log2.bw",
                                                       selectedsites = ENCODE_H3K27me3_promoter_sites,
                                                       genename = "PEMT",
                                                       gene_ensembl = PEMT_ensembl,
                                                       save = TRUE,
                                                       filename = "ENCODE_PEMT_H3K27_promoter_model")

ENCODE_PEMT_H3K4_promoter_model <- fit.models.to.chip(computematrix_outputdirectory = "ENCODE_H3K4_computematrix",
                                                      selectedtag = "promoter",
                                                      removetag = "_log2.bw",
                                                      selectedsites = ENCODE_H3K4me3_promoter_sites,
                                                      genename = "PEMT",
                                                      gene_ensembl = PEMT_ensembl,
                                                      save = TRUE,
                                                      filename = "ENCODE_PEMT_H3K4_promoter_model")

ENCODE_PEMT_H3K9_promoter_model <- fit.models.to.chip(computematrix_outputdirectory = "ENCODE_H3K9_computematrix",
                                                      selectedtag = "promoter",
                                                      removetag = "_log2.bw",
                                                      selectedsites = ENCODE_H3K9me3_promoter_sites,
                                                      genename = "PEMT",
                                                      gene_ensembl = PEMT_ensembl,
                                                      save = TRUE,
                                                      filename = "ENCODE_PEMT_H3K9_promoter_model")

ENCODE_PEMT_H3K36_promoter_model <- fit.models.to.chip(computematrix_outputdirectory = "ENCODE_H3K36_computematrix",
                                                       selectedtag = "promoter",
                                                       removetag = "_log2.bw",
                                                       selectedsites = ENCODE_H3K36me3_promoter_sites,
                                                       genename = "PEMT",
                                                       gene_ensembl = PEMT_ensembl,
                                                       save = TRUE,
                                                       filename = "ENCODE_PEMT_H3K36_promoter_model")

ENCODE_PEMT_H3K27_genebody_model <- fit.models.to.chip(computematrix_outputdirectory = "ENCODE_H3K27_computematrix",
                                                       selectedtag = "genebodies200bpflanks",
                                                       removetag = "_log2.bw",
                                                       selectedsites = ENCODE_H3K27me3_genebody_sites,
                                                       genename = "PEMT",
                                                       gene_ensembl = PEMT_ensembl,
                                                       save = TRUE,
                                                       filename = "ENCODE_PEMT_H3K27_genebody_model")

ENCODE_PEMT_H3K4_genebody_model <- fit.models.to.chip(computematrix_outputdirectory = "ENCODE_H3K4_computematrix",
                                                      selectedtag = "genebodies200bpflanks",
                                                      removetag = "_log2.bw",
                                                      selectedsites = ENCODE_H3K4me3_genebody_sites,
                                                      genename = "PEMT",
                                                      gene_ensembl = PEMT_ensembl,
                                                      save = TRUE,
                                                      filename = "ENCODE_PEMT_H3K4_genebody_model")

ENCODE_PEMT_H3K9_genebody_model <- fit.models.to.chip(computematrix_outputdirectory = "ENCODE_H3K9_computematrix",
                                                      selectedtag = "genebodies200bpflanks",
                                                      removetag = "_log2.bw",
                                                      selectedsites = ENCODE_H3K9me3_genebody_sites,
                                                      genename = "PEMT",
                                                      gene_ensembl = PEMT_ensembl,
                                                      save = TRUE,
                                                      filename = "ENCODE_PEMT_H3K9_genebody_model")

ENCODE_PEMT_H3K36_genebody_model <- fit.models.to.chip(computematrix_outputdirectory = "ENCODE_H3K36_computematrix",
                                                       selectedtag = "genebodies200bpflanks",
                                                       removetag = "_log2.bw",
                                                       selectedsites = ENCODE_H3K36me3_genebody_sites,
                                                       genename = "PEMT",
                                                       gene_ensembl = PEMT_ensembl,
                                                       save = TRUE,
                                                       filename = "ENCODE_PEMT_H3K36_genebody_model")

ENCODE_PEMT_H3K9_repeats_model <- fit.models.to.chip(computematrix_outputdirectory = "ENCODE_H3K9_computematrix",
                                                     selectedtag = "repeats",
                                                     removetag = "_log2.bw",
                                                     selectedsites = ENCODE_H3K9_repeats_sites,
                                                     genename = "PEMT",
                                                     gene_ensembl = PEMT_ensembl,
                                                     save = TRUE,
                                                     filename = "ENCODE_PEMT_H3K9_repeats_model")

# now for NNMT

ENCODE_NNMT_H3K27_promoter_model <- fit.models.to.chip(computematrix_outputdirectory = "ENCODE_H3K27_computematrix",
                                                       selectedtag = "promoter",
                                                       removetag = "_log2.bw",
                                                       selectedsites = ENCODE_H3K27me3_promoter_sites,
                                                       genename = "NNMT",
                                                       gene_ensembl = NNMT_ensembl,
                                                       save = TRUE,
                                                       filename = "ENCODE_NNMT_H3K27_promoter_model")

ENCODE_NNMT_H3K9_promoter_model <- fit.models.to.chip(computematrix_outputdirectory = "ENCODE_H3K9_computematrix",
                                                      selectedtag = "promoter",
                                                      removetag = "_log2.bw",
                                                      selectedsites = ENCODE_H3K9me3_promoter_sites,
                                                      genename = "NNMT",
                                                      gene_ensembl = NNMT_ensembl,
                                                      save = TRUE,
                                                      filename = "ENCODE_NNMT_H3K9_promoter_model")

ENCODE_NNMT_H3K36_promoter_model <- fit.models.to.chip(computematrix_outputdirectory = "ENCODE_H3K36_computematrix",
                                                       selectedtag = "promoter",
                                                       removetag = "_log2.bw",
                                                       selectedsites = ENCODE_H3K36me3_promoter_sites,
                                                       genename = "NNMT",
                                                       gene_ensembl = NNMT_ensembl,
                                                       save = TRUE,
                                                       filename = "ENCODE_NNMT_H3K36_promoter_model")

ENCODE_NNMT_H3K27_genebody_model <- fit.models.to.chip(computematrix_outputdirectory = "ENCODE_H3K27_computematrix",
                                                       selectedtag = "genebodies200bpflanks",
                                                       removetag = "_log2.bw",
                                                       selectedsites = ENCODE_H3K27me3_genebody_sites,
                                                       genename = "NNMT",
                                                       gene_ensembl = NNMT_ensembl,
                                                       save = TRUE,
                                                       filename = "ENCODE_NNMT_H3K27_genebody_model")

ENCODE_NNMT_H3K9_genebody_model <- fit.models.to.chip(computematrix_outputdirectory = "ENCODE_H3K9_computematrix",
                                                      selectedtag = "genebodies200bpflanks",
                                                      removetag = "_log2.bw",
                                                      selectedsites = ENCODE_H3K9me3_genebody_sites,
                                                      genename = "NNMT",
                                                      gene_ensembl = NNMT_ensembl,
                                                      save = TRUE,
                                                      filename = "ENCODE_NNMT_H3K9_genebody_model")

ENCODE_NNMT_H3K36_genebody_model <- fit.models.to.chip(computematrix_outputdirectory = "ENCODE_H3K36_computematrix",
                                                       selectedtag = "genebodies200bpflanks",
                                                       removetag = "_log2.bw",
                                                       selectedsites = ENCODE_H3K36me3_genebody_sites,
                                                       genename = "NNMT",
                                                       gene_ensembl = NNMT_ensembl,
                                                       save = TRUE,
                                                       filename = "ENCODE_NNMT_H3K36_genebody_model")

ENCODE_NNMT_H3K9_repeats_model <- fit.models.to.chip(computematrix_outputdirectory = "ENCODE_H3K9_computematrix",
                                                     selectedtag = "repeats",
                                                     removetag = "_log2.bw",
                                                     selectedsites = ENCODE_H3K9_repeats_sites,
                                                     genename = "NNMT",
                                                     gene_ensembl = NNMT_ensembl,
                                                     save = TRUE,
                                                     filename = "ENCODE_NNMT_H3K9_repeats_model")

ENCODE_NNMT_H3K4_promoter_model <- fit.models.to.chip(computematrix_outputdirectory = "ENCODE_H3K4_computematrix",
                                                      selectedtag = "promoter",
                                                      removetag = "_log2.bw",
                                                      selectedsites = ENCODE_H3K4me3_promoter_sites,
                                                      genename = "NNMT",
                                                      gene_ensembl = NNMT_ensembl,
                                                      save = TRUE,
                                                      filename = "ENCODE_NNMT_H3K4_promoter_model")

ENCODE_NNMT_H3K4_genebody_model <- fit.models.to.chip(computematrix_outputdirectory = "ENCODE_H3K4_computematrix",
                                                      selectedtag = "genebodies200bpflanks",
                                                      removetag = "_log2.bw",
                                                      selectedsites = ENCODE_H3K4me3_genebody_sites,
                                                      genename = "NNMT",
                                                      gene_ensembl = NNMT_ensembl,
                                                      save = TRUE,
                                                      filename = "ENCODE_NNMT_H3K4_genebody_model")

# make random sample of 1000 expressed genes (from genes with TPM > 5 across GTEX) to do random distribution

# first remove the HMTs and deMTs from the list
# nuclear_genes_exclude_epi <- nuclear_genes_expressed[!(nuclear_genes_expressed %in% c(SET_HMTs$ensembl_gene_id, deMTs$ensembl_gene_id))]
# 
# sample1000 <- sample(nuclear_genes_exclude_epi, size = 1000, replace = FALSE)
# 
# write.table(sample1000, "output/sample_1000randomgenes.txt", row.names = FALSE)
sample1000 <- read.table("output/sample_1000randomgenes.txt", header = TRUE)

ENCODE_H3K4_promoters_randommodel_1st250 <- chip.models.random.samples(computematrix_outputdirectory = "ENCODE_H3K4_computematrix",
                                                                       selectedtag = "promoter",
                                                                       removetag = "_log2.bw",
                                                                       selectedsites = ENCODE_H3K4me3_promoter_sites,
                                                                       generandomsample = unlist(sample1000)[1:250],
                                                                       save = TRUE,
                                                                       filename = "ENCODE_H3K4_promoters_randommodel_1st250")

ENCODE_H3K4_promoters_randommodel_2nd250 <- chip.models.random.samples(computematrix_outputdirectory = "ENCODE_H3K4_computematrix",
                                                                       selectedtag = "promoter",
                                                                       removetag = "_log2.bw",
                                                                       selectedsites = ENCODE_H3K4me3_promoter_sites,
                                                                       generandomsample = unlist(sample1000)[251:500],
                                                                       save = TRUE,
                                                                       filename = "ENCODE_H3K4_promoters_randommodel_2nd250")

ENCODE_H3K4_promoters_randommodel_3rd250 <- chip.models.random.samples(computematrix_outputdirectory = "ENCODE_H3K4_computematrix",
                                                                       selectedtag = "promoter",
                                                                       removetag = "_log2.bw",
                                                                       selectedsites = ENCODE_H3K4me3_promoter_sites,
                                                                       generandomsample = unlist(sample1000)[501:750],
                                                                       save = TRUE,
                                                                       filename = "ENCODE_H3K4_promoters_randommodel_3rd250")

ENCODE_H3K4_promoters_randommodel_4th250 <- chip.models.random.samples(computematrix_outputdirectory = "ENCODE_H3K4_computematrix",
                                                                       selectedtag = "promoter",
                                                                       removetag = "_log2.bw",
                                                                       selectedsites = ENCODE_H3K4me3_promoter_sites,
                                                                       generandomsample = unlist(sample1000)[751:1000],
                                                                       save = TRUE,
                                                                       filename = "ENCODE_H3K4_promoters_randommodel_4th250")

ENCODE_H3K9_promoters_randommodel <- chip.models.random.samples(computematrix_outputdirectory = "ENCODE_H3K9_computematrix",
                                                                selectedtag = "promoter",
                                                                removetag = "_log2.bw",
                                                                selectedsites = ENCODE_H3K9me3_promoter_sites,
                                                                generandomsample = unlist(sample1000),
                                                                save = TRUE,
                                                                filename = "ENCODE_H3K9_promoters_randommodel")

ENCODE_H3K27_promoters_randommodel <- chip.models.random.samples(computematrix_outputdirectory = "ENCODE_H3K27_computematrix",
                                                                 selectedtag = "promoter",
                                                                 removetag = "_log2.bw",
                                                                 selectedsites = ENCODE_H3K27me3_promoter_sites,
                                                                 generandomsample = unlist(sample1000),
                                                                 save = TRUE,
                                                                 filename = "ENCODE_H3K27_promoters_randommodel")

ENCODE_H3K36_promoters_randommodel <- chip.models.random.samples(computematrix_outputdirectory = "ENCODE_H3K36_computematrix",
                                                                 selectedtag = "promoter",
                                                                 removetag = "_log2.bw",
                                                                 selectedsites = ENCODE_H3K36me3_promoter_sites,
                                                                 generandomsample = unlist(sample1000),
                                                                 save = TRUE,
                                                                 filename = "ENCODE_H3K36_promoters_randommodel")

ENCODE_H3K4_genebody_randommodel_1st250 <- chip.models.random.samples(computematrix_outputdirectory = "ENCODE_H3K4_computematrix",
                                                                      selectedtag = "genebodies",
                                                                      removetag = "_log2.bw",
                                                                      selectedsites = ENCODE_H3K4me3_genebody_sites,
                                                                      generandomsample = unlist(sample1000)[1:250],
                                                                      save = TRUE,
                                                                      filename = "ENCODE_H3K4_genebody_randommodel_1st250")

ENCODE_H3K4_genebody_randommodel_2nd250 <- chip.models.random.samples(computematrix_outputdirectory = "ENCODE_H3K4_computematrix",
                                                                      selectedtag = "genebodies",
                                                                      removetag = "_log2.bw",
                                                                      selectedsites = ENCODE_H3K4me3_genebody_sites,
                                                                      generandomsample = unlist(sample1000)[251:500],
                                                                      save = TRUE,
                                                                      filename = "ENCODE_H3K4_genebody_randommodel_2nd250")

ENCODE_H3K4_genebody_randommodel_3rd250 <- chip.models.random.samples(computematrix_outputdirectory = "ENCODE_H3K4_computematrix",
                                                                      selectedtag = "genebodies",
                                                                      removetag = "_log2.bw",
                                                                      selectedsites = ENCODE_H3K4me3_genebody_sites,
                                                                      generandomsample = unlist(sample1000)[501:750],
                                                                      save = TRUE,
                                                                      filename = "ENCODE_H3K4_genebody_randommodel_3rd250")

ENCODE_H3K4_genebody_randommodel_4th250 <- chip.models.random.samples(computematrix_outputdirectory = "ENCODE_H3K4_computematrix",
                                                                      selectedtag = "genebodies",
                                                                      removetag = "_log2.bw",
                                                                      selectedsites = ENCODE_H3K4me3_genebody_sites,
                                                                      generandomsample = unlist(sample1000)[751:1000],
                                                                      save = TRUE,
                                                                      filename = "ENCODE_H3K4_genebody_randommodel_4th250")

ENCODE_H3K9_genebody_randommodel <- chip.models.random.samples(computematrix_outputdirectory = "ENCODE_H3K9_computematrix",
                                                               selectedtag = "genebodies",
                                                               removetag = "_log2.bw",
                                                               selectedsites = ENCODE_H3K9me3_genebody_sites,
                                                               generandomsample = unlist(sample1000),
                                                               save = TRUE,
                                                               filename = "ENCODE_H3K9_genebody_randommodel")

ENCODE_H3K27_genebody_randommodel <- chip.models.random.samples(computematrix_outputdirectory = "ENCODE_H3K27_computematrix",
                                                                selectedtag = "genebodies",
                                                                removetag = "_log2.bw",
                                                                selectedsites = ENCODE_H3K27me3_genebody_sites,
                                                                generandomsample = unlist(sample1000),
                                                                save = TRUE,
                                                                filename = "ENCODE_H3K27_genebody_randommodel")

ENCODE_H3K36_genebody_randommodel_1st500 <- chip.models.random.samples(computematrix_outputdirectory = "ENCODE_H3K36_computematrix",
                                                                selectedtag = "genebodies",
                                                                removetag = "_log2.bw",
                                                                selectedsites = ENCODE_H3K36me3_genebody_sites,
                                                                generandomsample = unlist(sample1000)[1:500],
                                                                save = TRUE,
                                                                filename = "ENCODE_H3K36_genebody_randommodel_1st500")

ENCODE_H3K36_genebody_randommodel_2nd500 <- chip.models.random.samples(computematrix_outputdirectory = "ENCODE_H3K36_computematrix",
                                                                       selectedtag = "genebodies",
                                                                       removetag = "_log2.bw",
                                                                       selectedsites = ENCODE_H3K36me3_genebody_sites,
                                                                       generandomsample = unlist(sample1000)[501:1000],
                                                                       save = TRUE,
                                                                       filename = "ENCODE_H3K36_genebody_randommodel_2nd500")

ENCODE_H3K9_repeats_randommodel <- chip.models.random.samples(computematrix_outputdirectory = "ENCODE_H3K9_computematrix",
                                                              selectedtag = "repeats",
                                                              removetag = "_log2.bw",
                                                              selectedsites = ENCODE_H3K9_repeats_sites,
                                                              generandomsample = unlist(sample1000),
                                                              save = TRUE,
                                                              filename = "ENCODE_H3K9_repeats_randommodel")

ENCODE_H3K4_promoter_signal <- spit.out.signal.matrix(computematrix_outputdirectory = "ENCODE_H3K4_computematrix",
                                                      selectedtag = "promoter",
                                                      removetag = "_log2.bw",
                                                      selectedsites = ENCODE_H3K4me3_promoter_sites)

ENCODE_H3K9_promoter_signal <- spit.out.signal.matrix(computematrix_outputdirectory = "ENCODE_H3K9_computematrix",
                                                      selectedtag = "promoter",
                                                      removetag = "_log2.bw",
                                                      selectedsites = ENCODE_H3K9me3_promoter_sites)

ENCODE_H3K27_promoter_signal <- spit.out.signal.matrix(computematrix_outputdirectory = "ENCODE_H3K27_computematrix",
                                                       selectedtag = "promoter",
                                                       removetag = "_log2.bw",
                                                       selectedsites = ENCODE_H3K27me3_promoter_sites)

ENCODE_H3K36_promoter_signal <- spit.out.signal.matrix(computematrix_outputdirectory = "ENCODE_H3K36_computematrix",
                                                       selectedtag = "promoter",
                                                       removetag = "_log2.bw",
                                                       selectedsites = ENCODE_H3K36me3_promoter_sites)

ENCODE_H3K4_genebody_signal <- spit.out.signal.matrix(computematrix_outputdirectory = "ENCODE_H3K4_computematrix",
                                                      selectedtag = "genebodies200bpflanks",
                                                      removetag = "_log2.bw",
                                                      selectedsites = ENCODE_H3K4me3_genebody_sites)

ENCODE_H3K9_genebody_signal <- spit.out.signal.matrix(computematrix_outputdirectory = "ENCODE_H3K9_computematrix",
                                                      selectedtag = "genebodies200bpflanks",
                                                      removetag = "_log2.bw",
                                                      selectedsites = ENCODE_H3K9me3_genebody_sites)

ENCODE_H3K27_genebody_signal <- spit.out.signal.matrix(computematrix_outputdirectory = "ENCODE_H3K27_computematrix",
                                                       selectedtag = "genebodies200bpflanks",
                                                       removetag = "_log2.bw",
                                                       selectedsites = ENCODE_H3K27me3_genebody_sites)

ENCODE_H3K36_genebody_signal <- spit.out.signal.matrix(computematrix_outputdirectory = "ENCODE_H3K36_computematrix",
                                                       selectedtag = "genebodies200bpflanks",
                                                       removetag = "_log2.bw",
                                                       selectedsites = ENCODE_H3K36me3_genebody_sites)

ENCODE_H3K9_repeats_signal <- spit.out.signal.matrix(computematrix_outputdirectory = "ENCODE_H3K9_computematrix",
                                                     selectedtag = "repeats",
                                                     removetag = "_log2.bw",
                                                     selectedsites = ENCODE_H3K9_repeats_sites)

ENCODE_H3K9_repeats_signal <- ENCODE_H3K9_repeats_signal[, apply(ENCODE_H3K9_repeats_signal, 2, function(x){any(!is.na(x))})]

ENCODE_H3K4_promoter_signal_sites <- names(ENCODE_H3K4_promoter_signal[colMeans(ENCODE_H3K4_promoter_signal) > 0])
ENCODE_H3K9_promoter_signal_sites <- names(ENCODE_H3K9_promoter_signal[colMeans(ENCODE_H3K9_promoter_signal) > 0])
ENCODE_H3K27_promoter_signal_sites <- names(ENCODE_H3K27_promoter_signal[colMeans(ENCODE_H3K27_promoter_signal) > 0])
ENCODE_H3K36_promoter_signal_sites <- names(ENCODE_H3K36_promoter_signal[colMeans(ENCODE_H3K36_promoter_signal) > 0])

ENCODE_H3K4_genebody_signal_sites <- names(ENCODE_H3K4_genebody_signal[colMeans(ENCODE_H3K4_genebody_signal) > 0])
ENCODE_H3K9_genebody_signal_sites <- names(ENCODE_H3K9_genebody_signal[colMeans(ENCODE_H3K9_genebody_signal) > 0])
ENCODE_H3K27_genebody_signal_sites <- names(ENCODE_H3K27_genebody_signal[colMeans(ENCODE_H3K27_genebody_signal) > 0])
ENCODE_H3K36_genebody_signal_sites <- names(ENCODE_H3K36_genebody_signal[colMeans(ENCODE_H3K36_genebody_signal) > 0])

ENCODE_H3K9_repeats_signal_sites <- names(ENCODE_H3K9_repeats_signal[colMeans(ENCODE_H3K9_repeats_signal) > 0])

ENCODE_H3K4random_promoters_signal_model <- readRDS("output/ENCODE_noTC_H3K4_promoters_randommodel.rds")
ENCODE_H3K9random_promoters_signal_model <- readRDS("output/ENCODE_noTC_H3K9_promoters_randommodel.rds")
ENCODE_H3K27random_promoters_signal_model <- readRDS("output/ENCODE_noTC_H3K27_promoters_randommodel.rds")
ENCODE_H3K36random_promoters_signal_model <- readRDS("output/ENCODE_noTC_H3K36_promoters_randommodel.rds")

ENCODE_H3K4random_genebody_signal_model <- readRDS("output/ENCODE_noTC_H3K4_genebody_randommodel.rds")
ENCODE_H3K9random_genebody_signal_model <- readRDS("output/ENCODE_noTC_H3K9_genebody_randommodel.rds")
ENCODE_H3K27random_genebody_signal_model <- readRDS("output/ENCODE_noTC_H3K27_genebody_randommodel.rds")
ENCODE_H3K36random_genebody_signal_model <- readRDS("output/ENCODE_noTC_H3K36_genebody_randommodel.rds")

ENCODE_H3K9random_repeats_signal_model <- readRDS("output/ENCODE_noTC_H3K9_repeats_randommodel.rds")

ENCODE_PEMT_H3K4_genebody_model <- readRDS("output/ENCODE_PEMT_H3K4_genebody_model.rds")
ENCODE_PEMT_H3K9_genebody_model <- readRDS("output/ENCODE_PEMT_H3K9_genebody_model.rds")
ENCODE_PEMT_H3K27_genebody_model <- readRDS("output/ENCODE_PEMT_H3K27_genebody_model.rds")
ENCODE_PEMT_H3K36_genebody_model <- readRDS("output/ENCODE_PEMT_H3K36_genebody_model.rds")

ENCODE_PEMT_H3K4_promoter_model <- readRDS("output/ENCODE_PEMT_H3K4_promoter_model.rds")
ENCODE_PEMT_H3K9_promoter_model <- readRDS("output/ENCODE_PEMT_H3K9_promoter_model.rds")
ENCODE_PEMT_H3K27_promoter_model <- readRDS("output/ENCODE_PEMT_H3K27_promoter_model.rds")
ENCODE_PEMT_H3K36_promoter_model <- readRDS("output/ENCODE_PEMT_H3K36_promoter_model.rds")

ENCODE_PEMT_H3K9_repeats_model <- readRDS("output/ENCODE_PEMT_H3K9_repeats_model.rds")

ENCODE_RNAseq_plot_list <- list(
  repeats = list(H3K9_PEMT = ENCODE_PEMT_H3K9_repeats_model[ENCODE_H3K9_repeats_signal_sites],
                 H3K9_random = rowMeans(ENCODE_H3K9random_repeats_signal_model[str_remove(row.names(ENCODE_H3K9random_repeats_signal_model), "\\.t_value") %in%  ENCODE_H3K9_repeats_signal_sites, ], na.rm = TRUE)),
  genes = list(H3K4_PEMT = ENCODE_PEMT_H3K4_genebody_model[names(ENCODE_PEMT_H3K4_genebody_model) %in% ENCODE_H3K4_genebody_signal_sites],
               H3K4_random = rowMeans(ENCODE_H3K4random_genebody_signal_model[str_remove(row.names(ENCODE_H3K4random_genebody_signal_model), "\\.t_value") %in%  ENCODE_H3K4_genebody_signal_sites, ], na.rm = TRUE),
               H3K9_PEMT = ENCODE_PEMT_H3K9_genebody_model[names(ENCODE_PEMT_H3K9_genebody_model) %in% ENCODE_H3K9_genebody_signal_sites], 
               H3K9_random = rowMeans(ENCODE_H3K9random_genebody_signal_model[str_remove(row.names(ENCODE_H3K9random_genebody_signal_model), "\\.t_value") %in% ENCODE_H3K9_genebody_signal_sites, ], na.rm = TRUE),
               H3K27_PEMT = ENCODE_PEMT_H3K27_genebody_model[names(ENCODE_PEMT_H3K27_genebody_model) %in% ENCODE_H3K27_genebody_signal_sites],
               H3K27_random = rowMeans(ENCODE_H3K27random_genebody_signal_model[str_remove(row.names(ENCODE_H3K27random_genebody_signal_model), "\\.t_value") %in% ENCODE_H3K27_genebody_signal_sites, ], na.rm = TRUE),
               H3K36_PEMT = ENCODE_PEMT_H3K36_genebody_model[names(ENCODE_PEMT_H3K36_genebody_model) %in% ENCODE_H3K36_genebody_signal_sites],
               H3K36_random = rowMeans(ENCODE_H3K36random_genebody_signal_model[str_remove(row.names(ENCODE_H3K36random_genebody_signal_model), "\\.t_value") %in% ENCODE_H3K36_genebody_signal_sites, ], na.rm = TRUE)),
  promoters = list(H3K4_PEMT = ENCODE_PEMT_H3K4_promoter_model[names(ENCODE_PEMT_H3K4_promoter_model) %in% ENCODE_H3K4_promoter_signal_sites],
                   H3K4_random = rowMeans(ENCODE_H3K4random_promoters_signal_model[str_remove(row.names(ENCODE_H3K4random_promoters_signal_model), "\\.t_value") %in%  ENCODE_H3K4_promoter_signal_sites, ], na.rm = TRUE),
                   H3K9_PEMT = ENCODE_PEMT_H3K9_promoter_model[names(ENCODE_PEMT_H3K9_promoter_model) %in% ENCODE_H3K9_promoter_signal_sites],
                   H3K9_random = rowMeans(ENCODE_H3K9random_promoters_signal_model[str_remove(row.names(ENCODE_H3K9random_promoters_signal_model), "\\.t_value") %in%  ENCODE_H3K9_promoter_signal_sites, ], na.rm = TRUE),
                   H3K27_PEMT = ENCODE_PEMT_H3K27_promoter_model[names(ENCODE_PEMT_H3K27_promoter_model) %in% ENCODE_H3K27_promoter_signal_sites],
                   H3K27_random = rowMeans(ENCODE_H3K27random_promoters_signal_model[str_remove(row.names(ENCODE_H3K27random_promoters_signal_model), "\\.t_value") %in%  ENCODE_H3K27_promoter_signal_sites, ], na.rm = TRUE),
                   H3K36_PEMT = ENCODE_PEMT_H3K36_promoter_model[names(ENCODE_PEMT_H3K36_promoter_model) %in% ENCODE_H3K36_promoter_signal_sites],
                   H3K36_random = rowMeans(ENCODE_H3K36random_promoters_signal_model[str_remove(row.names(ENCODE_H3K36random_promoters_signal_model), "\\.t_value") %in%  ENCODE_H3K36_promoter_signal_sites, ], na.rm = TRUE))
)         

ENCODE_RNAseq_plot_melt <- reshape2::melt(ENCODE_RNAseq_plot_list)
colnames(ENCODE_RNAseq_plot_melt) <- c("signal", "gene", "group", "location")



write.table(ENCODE_RNAseq_plot_melt,
            file = "plot_data/Fig 3/Fig_3A_plotdata.txt",
            col.names = TRUE,
            row.names = FALSE,
            sep = "\t")

factet_labs <- c("Gene bodies",
                 "Promoters",
                 "Repetitive elements")
names(factet_labs) <- c("genes",
                        "promoters",
                        "repeats")

options(scipen = 3)

ENCODE_samplesizes <- table(paste0(ENCODE_RNAseq_plot_melt$group, ENCODE_RNAseq_plot_melt$location))
ENCODE_samplesizes <- as.data.frame(ENCODE_samplesizes)
ENCODE_samplesizes[, "group"] <- str_remove(str_remove(str_remove(ENCODE_samplesizes$Var1, pattern = "promoters"), pattern = "genes"), pattern = "repeats")

for(i in 1:nrow(ENCODE_samplesizes)){
  ENCODE_samplesizes[i, "location"] <- str_remove(ENCODE_samplesizes[i, "Var1"], pattern = ENCODE_samplesizes[i, "group"])
}

ENCODE_samplesizes$Freq = paste0("n = ", ENCODE_samplesizes$Freq)

ENCODE_samplesizes[!str_detect(ENCODE_samplesizes$Var1, pattern = "random"), "Freq"] <- ""

ENCODE_samplesizes[ENCODE_samplesizes$Var1 == "H3K4_PEMTgenes", "p_value"] <- signif(wilcox.test(ENCODE_RNAseq_plot_list$genes$H3K4_random, unlist(ENCODE_RNAseq_plot_list$genes$H3K4_PEMT)[names(ENCODE_RNAseq_plot_list$genes$H3K4_random)], paired = TRUE)$p.value, digits = 2)
ENCODE_samplesizes[ENCODE_samplesizes$Var1 == "H3K4_randomgenes", "p_value_parse"] <- paste0(str_remove(ENCODE_samplesizes[ENCODE_samplesizes$Var1 == "H3K4_PEMTgenes", "p_value"], pattern = "e.*"), " * 10^", str_remove(ENCODE_samplesizes[ENCODE_samplesizes$Var1 == "H3K4_PEMTgenes", "p_value"], pattern = "p = .*e"))

ENCODE_samplesizes[ENCODE_samplesizes$Var1 == "H3K4_PEMTpromoters", "p_value"] <- signif(wilcox.test(ENCODE_RNAseq_plot_list$promoters$H3K4_random, unlist(ENCODE_RNAseq_plot_list$promoters$H3K4_PEMT)[names(ENCODE_RNAseq_plot_list$promoters$H3K4_random)], paired = TRUE)$p.value, digits = 2)
ENCODE_samplesizes[ENCODE_samplesizes$Var1 == "H3K4_randompromoters", "p_value_parse"] <- paste0(str_remove(ENCODE_samplesizes[ENCODE_samplesizes$Var1 == "H3K4_PEMTpromoters", "p_value"], pattern = "e.*"), " * 10^", str_remove(ENCODE_samplesizes[ENCODE_samplesizes$Var1 == "H3K4_PEMTpromoters", "p_value"], pattern = "p = .*e"))

ENCODE_samplesizes[ENCODE_samplesizes$Var1 == "H3K9_PEMTgenes", "p_value"] <- signif(wilcox.test(ENCODE_RNAseq_plot_list$genes$H3K9_random, unlist(ENCODE_RNAseq_plot_list$genes$H3K9_PEMT)[names(ENCODE_RNAseq_plot_list$genes$H3K9_random)], paired = TRUE)$p.value, digits = 2)
ENCODE_samplesizes[ENCODE_samplesizes$Var1 == "H3K9_randomgenes", "p_value_parse"] <- paste0(str_remove(ENCODE_samplesizes[ENCODE_samplesizes$Var1 == "H3K9_PEMTgenes", "p_value"], pattern = "e.*"), " * 10^", str_remove(ENCODE_samplesizes[ENCODE_samplesizes$Var1 == "H3K9_PEMTgenes", "p_value"], pattern = "p = .*e"))

ENCODE_samplesizes[ENCODE_samplesizes$Var1 == "H3K9_PEMTpromoters", "p_value"] <- signif(wilcox.test(ENCODE_RNAseq_plot_list$promoters$H3K9_random, unlist(ENCODE_RNAseq_plot_list$promoters$H3K9_PEMT)[names(ENCODE_RNAseq_plot_list$promoters$H3K9_random)], paired = TRUE)$p.value, digits = 2)
ENCODE_samplesizes[ENCODE_samplesizes$Var1 == "H3K9_randompromoters", "p_value_parse"] <- paste0(str_remove(ENCODE_samplesizes[ENCODE_samplesizes$Var1 == "H3K9_PEMTpromoters", "p_value"], pattern = "e.*"), " * 10^", str_remove(ENCODE_samplesizes[ENCODE_samplesizes$Var1 == "H3K9_PEMTpromoters", "p_value"], pattern = "p = .*e"))

ENCODE_samplesizes[ENCODE_samplesizes$Var1 == "H3K9_PEMTrepeats", "p_value"] <- signif(wilcox.test(ENCODE_RNAseq_plot_list$repeats$H3K9_random, unlist(ENCODE_RNAseq_plot_list$repeats$H3K9_PEMT)[names(ENCODE_RNAseq_plot_list$repeats$H3K9_random)], paired = TRUE)$p.value, digits = 2)
ENCODE_samplesizes[ENCODE_samplesizes$Var1 == "H3K9_randomrepeats", "p_value_parse"] <- paste0(str_remove(ENCODE_samplesizes[ENCODE_samplesizes$Var1 == "H3K9_PEMTrepeats", "p_value"], pattern = "e.*"), " * 10^", str_remove(ENCODE_samplesizes[ENCODE_samplesizes$Var1 == "H3K9_PEMTrepeats", "p_value"], pattern = "p = .*e"))

ENCODE_samplesizes[ENCODE_samplesizes$Var1 == "H3K27_PEMTgenes", "p_value"] <- signif(wilcox.test(ENCODE_RNAseq_plot_list$genes$H3K27_random, unlist(ENCODE_RNAseq_plot_list$genes$H3K27_PEMT)[names(ENCODE_RNAseq_plot_list$genes$H3K27_random)], paired = TRUE)$p.value, digits = 2)
ENCODE_samplesizes[ENCODE_samplesizes$Var1 == "H3K27_randomgenes", "p_value_parse"] <- paste0(str_remove(ENCODE_samplesizes[ENCODE_samplesizes$Var1 == "H3K27_PEMTgenes", "p_value"], pattern = "e.*"), " * 10^", str_remove(ENCODE_samplesizes[ENCODE_samplesizes$Var1 == "H3K27_PEMTgenes", "p_value"], pattern = "p = .*e"))

ENCODE_samplesizes[ENCODE_samplesizes$Var1 == "H3K27_PEMTpromoters", "p_value"] <- signif(wilcox.test(ENCODE_RNAseq_plot_list$promoters$H3K27_random, unlist(ENCODE_RNAseq_plot_list$promoters$H3K27_PEMT)[names(ENCODE_RNAseq_plot_list$promoters$H3K27_random)], paired = TRUE)$p.value, digits = 2)
ENCODE_samplesizes[ENCODE_samplesizes$Var1 == "H3K27_randompromoters", "p_value_parse"] <- paste0(str_remove(ENCODE_samplesizes[ENCODE_samplesizes$Var1 == "H3K27_PEMTpromoters", "p_value"], pattern = "e.*"), " * 10^", str_remove(ENCODE_samplesizes[ENCODE_samplesizes$Var1 == "H3K27_PEMTpromoters", "p_value"], pattern = "p = .*e"))

ENCODE_samplesizes[ENCODE_samplesizes$Var1 == "H3K36_PEMTgenes", "p_value"] <- signif(wilcox.test(ENCODE_RNAseq_plot_list$genes$H3K36_random, unlist(ENCODE_RNAseq_plot_list$genes$H3K36_PEMT)[names(ENCODE_RNAseq_plot_list$genes$H3K36_random)], paired = TRUE)$p.value, digits = 2)
ENCODE_samplesizes[ENCODE_samplesizes$Var1 == "H3K36_randomgenes", "p_value_parse"] <- paste0(str_remove(ENCODE_samplesizes[ENCODE_samplesizes$Var1 == "H3K36_PEMTgenes", "p_value"], pattern = "e.*"), " * 10^", str_remove(ENCODE_samplesizes[ENCODE_samplesizes$Var1 == "H3K36_PEMTgenes", "p_value"], pattern = "p = .*e"))

ENCODE_samplesizes[ENCODE_samplesizes$Var1 == "H3K36_PEMTpromoters", "p_value"] <- signif(wilcox.test(ENCODE_RNAseq_plot_list$promoters$H3K36_random, unlist(ENCODE_RNAseq_plot_list$promoters$H3K36_PEMT)[names(ENCODE_RNAseq_plot_list$promoter$H3K36_random)], paired = TRUE)$p.value, digits = 2)
ENCODE_samplesizes[ENCODE_samplesizes$Var1 == "H3K36_randompromoters", "p_value_parse"] <- paste0(str_remove(ENCODE_samplesizes[ENCODE_samplesizes$Var1 == "H3K36_PEMTpromoters", "p_value"], pattern = "e.*"), " * 10^", str_remove(ENCODE_samplesizes[ENCODE_samplesizes$Var1 == "H3K36_PEMTpromoters", "p_value"], pattern = "p = .*e"))

ENCODE_dat_text <- data.frame(
  label = ENCODE_samplesizes[!(is.na(ENCODE_samplesizes$p_value)) & ENCODE_samplesizes$p_value > -0.001, "p_value"],
  location   = ENCODE_samplesizes[!(is.na(ENCODE_samplesizes$p_value)) & ENCODE_samplesizes$p_value > -0.001, "location"],
  group = ENCODE_samplesizes[!(is.na(ENCODE_samplesizes$p_value)) & ENCODE_samplesizes$p_value > -0.001, "group"]
)

ENCODE_dat_text[, "p_value_parse"] <- paste0(str_remove(ENCODE_dat_text$label, pattern = "e.*"), " %*% 10^", str_remove(ENCODE_dat_text$label, pattern = ".*e"))
ENCODE_dat_text[!is.na(ENCODE_dat_text$label) & ENCODE_dat_text$label == "0", "p_value_parse"] <- "0"

ENCODE_RNAseq_plot_melt$group <- factor(ENCODE_RNAseq_plot_melt$group, levels = c("H3K36_random", "H3K36_PEMT", "H3K27_random", "H3K27_PEMT", "H3K9_random", "H3K9_PEMT", "H3K4_random", "H3K4_PEMT"))

pdf("graphics/ENCODE_allmarks.pdf")

ggplot(aes(x = signal, y = group, group = location), data = ENCODE_RNAseq_plot_melt ) + 
  geom_boxplot(aes(group = group), fill = c("grey",
                                            "darkseagreen1",
                                            "grey",
                                            "darkorchid1",
                                            "grey",
                                            "magenta",
                                            "grey",
                                            "darkolivegreen3",
                                            "grey",
                                            "darkseagreen1",
                                            "grey",
                                            "darkorchid1",
                                            "grey",
                                            "magenta",
                                            "grey",
                                            "darkolivegreen3",
                                            "grey",
                                            "magenta"), notch = TRUE, outlier.shape = NA) +
  theme_classic() +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey") +
  xlab("Linear model t-value") + 
  ylab(NULL) +
  scale_y_discrete(label = c("H3K36me3 (random genes)",
                             "H3K36me3 (PEMT RNA-seq)",
                             "H3K27me3 (random genes)",
                             "H3K27me3 (PEMT RNA-seq)",
                             "H3K9me3 (random genes)",
                             "H3K9me3 (PEMT RNA-seq)",
                             "H3K4me3 (random genes)",
                             "H3K4me3 (PEMT RNA-seq)")) +
  facet_wrap(~ location, ncol = 1, labeller = labeller(location = factet_labs), scales = "free") + 
  coord_cartesian(xlim = c(-4, 5.5)) + 
  theme(strip.text = element_text(size = 12, color = "white", face = "bold"),
        strip.background = element_rect(fill = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.title.x = element_text(size = 10, color = "black", face = "bold")) +
  geom_text(aes(x = 4.75, group = location, label = p_value_parse), data = ENCODE_dat_text, parse = TRUE) +
  geom_text(aes(x = 4.75, group = location, label = Freq), data = ENCODE_samplesizes)

dev.off()

pdf("graphics/ENCODE_allmarks.pdf",
    width = 3.75,
    height = 5)

ggplot(aes(x = signal, y = group, group = location), data = ENCODE_RNAseq_plot_melt ) + 
  geom_boxplot(aes(group = group), fill = c("grey",
                                            "darkseagreen1",
                                            "grey",
                                            "darkorchid1",
                                            "grey",
                                            "magenta",
                                            "grey",
                                            "darkolivegreen3",
                                            "grey",
                                            "darkseagreen1",
                                            "grey",
                                            "darkorchid1",
                                            "grey",
                                            "magenta",
                                            "grey",
                                            "darkolivegreen3",
                                            "grey",
                                            "magenta"), notch = TRUE, outlier.shape = NA) +
  theme_classic() +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey") +
  xlab("Linear model t-value") + 
  ylab(NULL) +
  facet_wrap(~ location, ncol = 1, labeller = labeller(location = factet_labs), scales = "free") + 
  coord_cartesian(xlim = c(-4, 5.5)) + 
  theme(strip.text = element_text(size = 12, color = "white", face = "bold"),
        strip.background = element_rect(fill = "black"),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.title.x = element_text(size = 10, color = "black", face = "bold")) +
  geom_text(aes(x = 4.75, group = location, label = p_value_parse), data = ENCODE_dat_text, parse = TRUE, size =  3) +
  geom_text(aes(x = 4.75, group = location, label = Freq), data = ENCODE_samplesizes, size = 3)

dev.off()

pdf("graphics/ENCODE_allmarks_nopval.pdf",
    width = 3.75,
    height = 5)

ggplot(aes(x = signal, y = group, group = location), data = ENCODE_RNAseq_plot_melt ) + 
  geom_boxplot(aes(group = group), fill = c("grey",
                                            "darkseagreen1",
                                            "grey",
                                            "darkorchid1",
                                            "grey",
                                            "magenta",
                                            "grey",
                                            "darkolivegreen3",
                                            "grey",
                                            "darkseagreen1",
                                            "grey",
                                            "darkorchid1",
                                            "grey",
                                            "magenta",
                                            "grey",
                                            "darkolivegreen3",
                                            "grey",
                                            "magenta"), 
               notch = TRUE, 
               outlier.shape = NA,
               lwd = 0.3,
               fatten = 0.5) +
  theme_classic() +
  geom_vline(xintercept = 0,
             linetype = "dashed",
             colour = "grey") +
  xlab("Linear model t-value") + 
  ylab(NULL) +
  facet_wrap(~ location, ncol = 1, labeller = labeller(location = factet_labs), scales = "free") + 
  coord_cartesian(xlim = c(-3, 4.5)) + 
  theme(strip.text = element_text(size = 12, color = "white", face = "bold"),
        strip.background = element_rect(fill = "black"),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.title.x = element_text(size = 10, color = "black", face = "bold")) +
  geom_text(aes(x = 3, group = location, label = Freq), data = ENCODE_samplesizes, size = 3, hjust = 0)

dev.off()

#### ENCODE NNMT Fig S8C ####

ENCODE_NNMT_H3K4_genebody_model <- readRDS("ENCODE_redo_noTC/output/ENCODE_NNMT_H3K4_genebody_model.rds")
ENCODE_NNMT_H3K9_genebody_model <- readRDS("ENCODE_redo_noTC/output/ENCODE_NNMT_H3K9_genebody_model.rds")
ENCODE_NNMT_H3K27_genebody_model <- readRDS("ENCODE_redo_noTC/output/ENCODE_NNMT_H3K27_genebody_model.rds")
ENCODE_NNMT_H3K36_genebody_model <- readRDS("ENCODE_redo_noTC/output/ENCODE_NNMT_H3K36_genebody_model.rds")

ENCODE_NNMT_H3K4_promoter_model <- readRDS("ENCODE_redo_noTC/output/ENCODE_NNMT_H3K4_promoter_model.rds")
ENCODE_NNMT_H3K9_promoter_model <- readRDS("ENCODE_redo_noTC/output/ENCODE_NNMT_H3K9_promoter_model.rds")
ENCODE_NNMT_H3K27_promoter_model <- readRDS("ENCODE_redo_noTC/output/ENCODE_NNMT_H3K27_promoter_model.rds")
ENCODE_NNMT_H3K36_promoter_model <- readRDS("ENCODE_redo_noTC/output/ENCODE_NNMT_H3K36_promoter_model.rds")

ENCODE_NNMT_H3K9_repeats_model <- readRDS("ENCODE_redo_noTC/output/ENCODE_NNMT_H3K9_repeats_model.rds")

ENCODE_NNMT_RNAseq_plot_list <- list(
  repeats = list(H3K9_NNMT = ENCODE_NNMT_H3K9_repeats_model[ENCODE_H3K9_repeats_signal_sites],
                 H3K9_random = rowMeans(ENCODE_H3K9random_repeats_signal_model[str_remove(row.names(ENCODE_H3K9random_repeats_signal_model), "\\.t_value") %in%  ENCODE_H3K9_repeats_signal_sites, ], na.rm = TRUE)),
  genes = list(H3K4_NNMT = ENCODE_NNMT_H3K4_genebody_model[names(ENCODE_NNMT_H3K4_genebody_model) %in% ENCODE_H3K4_genebody_signal_sites],
               H3K4_random = rowMeans(ENCODE_H3K4random_genebody_signal_model[str_remove(row.names(ENCODE_H3K4random_genebody_signal_model), "\\.t_value") %in%  ENCODE_H3K4_genebody_signal_sites, ], na.rm = TRUE),
               H3K9_NNMT = ENCODE_NNMT_H3K9_genebody_model[names(ENCODE_NNMT_H3K9_genebody_model) %in% ENCODE_H3K9_genebody_signal_sites], 
               H3K9_random = rowMeans(ENCODE_H3K9random_genebody_signal_model[str_remove(row.names(ENCODE_H3K9random_genebody_signal_model), "\\.t_value") %in% ENCODE_H3K9_genebody_signal_sites, ], na.rm = TRUE),
               H3K27_NNMT = ENCODE_NNMT_H3K27_genebody_model[names(ENCODE_NNMT_H3K27_genebody_model) %in% ENCODE_H3K27_genebody_signal_sites],
               H3K27_random = rowMeans(ENCODE_H3K27random_genebody_signal_model[str_remove(row.names(ENCODE_H3K27random_genebody_signal_model), "\\.t_value") %in% ENCODE_H3K27_genebody_signal_sites, ], na.rm = TRUE),
               H3K36_NNMT = ENCODE_NNMT_H3K36_genebody_model[names(ENCODE_NNMT_H3K36_genebody_model) %in% ENCODE_H3K36_genebody_signal_sites],
               H3K36_random = rowMeans(ENCODE_H3K36random_genebody_signal_model[str_remove(row.names(ENCODE_H3K36random_genebody_signal_model), "\\.t_value") %in% ENCODE_H3K36_genebody_signal_sites, ], na.rm = TRUE)),
  promoters = list(H3K4_NNMT = ENCODE_NNMT_H3K4_promoter_model[names(ENCODE_NNMT_H3K4_promoter_model) %in% ENCODE_H3K4_promoter_signal_sites],
                   H3K4_random = rowMeans(ENCODE_H3K4random_promoters_signal_model[str_remove(row.names(ENCODE_H3K4random_promoters_signal_model), "\\.t_value") %in%  ENCODE_H3K4_promoter_signal_sites, ], na.rm = TRUE),
                   H3K9_NNMT = ENCODE_NNMT_H3K9_promoter_model[names(ENCODE_NNMT_H3K9_promoter_model) %in% ENCODE_H3K9_promoter_signal_sites],
                   H3K9_random = rowMeans(ENCODE_H3K9random_promoters_signal_model[str_remove(row.names(ENCODE_H3K9random_promoters_signal_model), "\\.t_value") %in%  ENCODE_H3K9_promoter_signal_sites, ], na.rm = TRUE),
                   H3K27_NNMT = ENCODE_NNMT_H3K27_promoter_model[names(ENCODE_NNMT_H3K27_promoter_model) %in% ENCODE_H3K27_promoter_signal_sites],
                   H3K27_random = rowMeans(ENCODE_H3K27random_promoters_signal_model[str_remove(row.names(ENCODE_H3K27random_promoters_signal_model), "\\.t_value") %in%  ENCODE_H3K27_promoter_signal_sites, ], na.rm = TRUE),
                   H3K36_NNMT = ENCODE_NNMT_H3K36_promoter_model[names(ENCODE_NNMT_H3K36_promoter_model) %in% ENCODE_H3K36_promoter_signal_sites],
                   H3K36_random = rowMeans(ENCODE_H3K36random_promoters_signal_model[str_remove(row.names(ENCODE_H3K36random_promoters_signal_model), "\\.t_value") %in%  ENCODE_H3K36_promoter_signal_sites, ], na.rm = TRUE))
)         

ENCODE_NNMT_RNAseq_plot_melt <- reshape2::melt(ENCODE_NNMT_RNAseq_plot_list)
colnames(ENCODE_NNMT_RNAseq_plot_melt) <- c("signal", "gene", "group", "location")

write.table(ENCODE_NNMT_RNAseq_plot_melt[, c(1:2, 4)],
            file = "plot_data/Fig S8/Fig_S8C_plotdata.txt",
            row.names = FALSE,
            col.names = TRUE,
            sep = "\t")

factet_labs <- c("Gene bodies",
                 "Promoters",
                 "Repetitive elements")
names(factet_labs) <- c("genes",
                        "promoters",
                        "repeats")

options(scipen = 3)

ENCODE_NNMT_samplesizes <- table(paste0(ENCODE_NNMT_RNAseq_plot_melt$group, ENCODE_NNMT_RNAseq_plot_melt$location))
ENCODE_NNMT_samplesizes <- as.data.frame(ENCODE_NNMT_samplesizes)
ENCODE_NNMT_samplesizes[, "group"] <- str_remove(str_remove(str_remove(ENCODE_NNMT_samplesizes$Var1, pattern = "promoters"), pattern = "genes"), pattern = "repeats")

for(i in 1:nrow(ENCODE_NNMT_samplesizes)){
  ENCODE_NNMT_samplesizes[i, "location"] <- str_remove(ENCODE_NNMT_samplesizes[i, "Var1"], pattern = ENCODE_NNMT_samplesizes[i, "group"])
}

ENCODE_NNMT_samplesizes$Freq = paste0("n = ", ENCODE_NNMT_samplesizes$Freq)

ENCODE_NNMT_samplesizes[!str_detect(ENCODE_NNMT_samplesizes$Var1, pattern = "random"), "Freq"] <- ""

ENCODE_NNMT_samplesizes[ENCODE_NNMT_samplesizes$Var1 == "H3K4_NNMTgenes", "p_value"] <- signif(wilcox.test(ENCODE_NNMT_RNAseq_plot_list$genes$H3K4_random, unlist(ENCODE_NNMT_RNAseq_plot_list$genes$H3K4_NNMT)[names(ENCODE_NNMT_RNAseq_plot_list$genes$H3K4_random)], paired = TRUE)$p.value, digits = 2)
ENCODE_NNMT_samplesizes[ENCODE_NNMT_samplesizes$Var1 == "H3K4_randomgenes", "p_value_parse"] <- paste0(str_remove(ENCODE_NNMT_samplesizes[ENCODE_NNMT_samplesizes$Var1 == "H3K4_NNMTgenes", "p_value"], pattern = "e.*"), " * 10^", str_remove(ENCODE_NNMT_samplesizes[ENCODE_NNMT_samplesizes$Var1 == "H3K4_NNMTgenes", "p_value"], pattern = "p = .*e"))

ENCODE_NNMT_samplesizes[ENCODE_NNMT_samplesizes$Var1 == "H3K4_NNMTpromoters", "p_value"] <- signif(wilcox.test(ENCODE_NNMT_RNAseq_plot_list$promoters$H3K4_random, unlist(ENCODE_NNMT_RNAseq_plot_list$promoters$H3K4_NNMT)[names(ENCODE_NNMT_RNAseq_plot_list$promoters$H3K4_random)], paired = TRUE)$p.value, digits = 2)
ENCODE_NNMT_samplesizes[ENCODE_NNMT_samplesizes$Var1 == "H3K4_randompromoters", "p_value_parse"] <- paste0(str_remove(ENCODE_NNMT_samplesizes[ENCODE_NNMT_samplesizes$Var1 == "H3K4_NNMTpromoters", "p_value"], pattern = "e.*"), " * 10^", str_remove(ENCODE_NNMT_samplesizes[ENCODE_NNMT_samplesizes$Var1 == "H3K4_NNMTpromoters", "p_value"], pattern = "p = .*e"))

ENCODE_NNMT_samplesizes[ENCODE_NNMT_samplesizes$Var1 == "H3K9_NNMTgenes", "p_value"] <- signif(wilcox.test(ENCODE_NNMT_RNAseq_plot_list$genes$H3K9_random, unlist(ENCODE_NNMT_RNAseq_plot_list$genes$H3K9_NNMT)[names(ENCODE_NNMT_RNAseq_plot_list$genes$H3K9_random)], paired = TRUE)$p.value, digits = 2)
ENCODE_NNMT_samplesizes[ENCODE_NNMT_samplesizes$Var1 == "H3K9_randomgenes", "p_value_parse"] <- paste0(str_remove(ENCODE_NNMT_samplesizes[ENCODE_NNMT_samplesizes$Var1 == "H3K9_NNMTgenes", "p_value"], pattern = "e.*"), " * 10^", str_remove(ENCODE_NNMT_samplesizes[ENCODE_NNMT_samplesizes$Var1 == "H3K9_NNMTgenes", "p_value"], pattern = "p = .*e"))

ENCODE_NNMT_samplesizes[ENCODE_NNMT_samplesizes$Var1 == "H3K9_NNMTpromoters", "p_value"] <- signif(wilcox.test(ENCODE_NNMT_RNAseq_plot_list$promoters$H3K9_random, unlist(ENCODE_NNMT_RNAseq_plot_list$promoters$H3K9_NNMT)[names(ENCODE_NNMT_RNAseq_plot_list$promoters$H3K9_random)], paired = TRUE)$p.value, digits = 2)
ENCODE_NNMT_samplesizes[ENCODE_NNMT_samplesizes$Var1 == "H3K9_randompromoters", "p_value_parse"] <- paste0(str_remove(ENCODE_NNMT_samplesizes[ENCODE_NNMT_samplesizes$Var1 == "H3K9_NNMTpromoters", "p_value"], pattern = "e.*"), " * 10^", str_remove(ENCODE_NNMT_samplesizes[ENCODE_NNMT_samplesizes$Var1 == "H3K9_NNMTpromoters", "p_value"], pattern = "p = .*e"))

ENCODE_NNMT_samplesizes[ENCODE_NNMT_samplesizes$Var1 == "H3K9_NNMTrepeats", "p_value"] <- signif(wilcox.test(ENCODE_NNMT_RNAseq_plot_list$repeats$H3K9_random, unlist(ENCODE_NNMT_RNAseq_plot_list$repeats$H3K9_NNMT)[names(ENCODE_NNMT_RNAseq_plot_list$repeats$H3K9_random)], paired = TRUE)$p.value, digits = 2)
ENCODE_NNMT_samplesizes[ENCODE_NNMT_samplesizes$Var1 == "H3K9_randomrepeats", "p_value_parse"] <- paste0(str_remove(ENCODE_NNMT_samplesizes[ENCODE_NNMT_samplesizes$Var1 == "H3K9_NNMTrepeats", "p_value"], pattern = "e.*"), " * 10^", str_remove(ENCODE_NNMT_samplesizes[ENCODE_NNMT_samplesizes$Var1 == "H3K9_NNMTrepeats", "p_value"], pattern = "p = .*e"))

ENCODE_NNMT_samplesizes[ENCODE_NNMT_samplesizes$Var1 == "H3K27_NNMTgenes", "p_value"] <- signif(wilcox.test(ENCODE_NNMT_RNAseq_plot_list$genes$H3K27_random, unlist(ENCODE_NNMT_RNAseq_plot_list$genes$H3K27_NNMT)[names(ENCODE_NNMT_RNAseq_plot_list$genes$H3K27_random)], paired = TRUE)$p.value, digits = 2)
ENCODE_NNMT_samplesizes[ENCODE_NNMT_samplesizes$Var1 == "H3K27_randomgenes", "p_value_parse"] <- paste0(str_remove(ENCODE_NNMT_samplesizes[ENCODE_NNMT_samplesizes$Var1 == "H3K27_NNMTgenes", "p_value"], pattern = "e.*"), " * 10^", str_remove(ENCODE_NNMT_samplesizes[ENCODE_NNMT_samplesizes$Var1 == "H3K27_NNMTgenes", "p_value"], pattern = "p = .*e"))

ENCODE_NNMT_samplesizes[ENCODE_NNMT_samplesizes$Var1 == "H3K27_NNMTpromoters", "p_value"] <- signif(wilcox.test(ENCODE_NNMT_RNAseq_plot_list$promoters$H3K27_random, unlist(ENCODE_NNMT_RNAseq_plot_list$promoters$H3K27_NNMT)[names(ENCODE_NNMT_RNAseq_plot_list$promoters$H3K27_random)], paired = TRUE)$p.value, digits = 2)
ENCODE_NNMT_samplesizes[ENCODE_NNMT_samplesizes$Var1 == "H3K27_randompromoters", "p_value_parse"] <- paste0(str_remove(ENCODE_NNMT_samplesizes[ENCODE_NNMT_samplesizes$Var1 == "H3K27_NNMTpromoters", "p_value"], pattern = "e.*"), " * 10^", str_remove(ENCODE_NNMT_samplesizes[ENCODE_NNMT_samplesizes$Var1 == "H3K27_NNMTpromoters", "p_value"], pattern = "p = .*e"))

ENCODE_NNMT_samplesizes[ENCODE_NNMT_samplesizes$Var1 == "H3K36_NNMTgenes", "p_value"] <- signif(wilcox.test(ENCODE_NNMT_RNAseq_plot_list$genes$H3K36_random, unlist(ENCODE_NNMT_RNAseq_plot_list$genes$H3K36_NNMT)[names(ENCODE_NNMT_RNAseq_plot_list$genes$H3K36_random)], paired = TRUE)$p.value, digits = 2)
ENCODE_NNMT_samplesizes[ENCODE_NNMT_samplesizes$Var1 == "H3K36_randomgenes", "p_value_parse"] <- paste0(str_remove(ENCODE_NNMT_samplesizes[ENCODE_NNMT_samplesizes$Var1 == "H3K36_NNMTgenes", "p_value"], pattern = "e.*"), " * 10^", str_remove(ENCODE_NNMT_samplesizes[ENCODE_NNMT_samplesizes$Var1 == "H3K36_NNMTgenes", "p_value"], pattern = "p = .*e"))

ENCODE_NNMT_samplesizes[ENCODE_NNMT_samplesizes$Var1 == "H3K36_NNMTpromoters", "p_value"] <- signif(wilcox.test(ENCODE_NNMT_RNAseq_plot_list$promoters$H3K36_random, unlist(ENCODE_NNMT_RNAseq_plot_list$promoters$H3K36_NNMT)[names(ENCODE_NNMT_RNAseq_plot_list$promoter$H3K36_random)], paired = TRUE)$p.value, digits = 2)
ENCODE_NNMT_samplesizes[ENCODE_NNMT_samplesizes$Var1 == "H3K36_randompromoters", "p_value_parse"] <- paste0(str_remove(ENCODE_NNMT_samplesizes[ENCODE_NNMT_samplesizes$Var1 == "H3K36_NNMTpromoters", "p_value"], pattern = "e.*"), " * 10^", str_remove(ENCODE_NNMT_samplesizes[ENCODE_NNMT_samplesizes$Var1 == "H3K36_NNMTpromoters", "p_value"], pattern = "p = .*e"))

ENCODE_NNMT_dat_text <- data.frame(
  label = ENCODE_NNMT_samplesizes[!(is.na(ENCODE_NNMT_samplesizes$p_value)) & ENCODE_NNMT_samplesizes$p_value > -0.001, "p_value"],
  location   = ENCODE_NNMT_samplesizes[!(is.na(ENCODE_NNMT_samplesizes$p_value)) & ENCODE_NNMT_samplesizes$p_value > -0.001, "location"],
  group = ENCODE_NNMT_samplesizes[!(is.na(ENCODE_NNMT_samplesizes$p_value)) & ENCODE_NNMT_samplesizes$p_value > -0.001, "group"]
)

ENCODE_NNMT_dat_text[, "p_value_parse"] <- paste0(str_remove(ENCODE_NNMT_dat_text$label, pattern = "e.*"), " %*% 10^", str_remove(ENCODE_NNMT_dat_text$label, pattern = ".*e"))
ENCODE_NNMT_dat_text[!is.na(ENCODE_NNMT_dat_text$label) & ENCODE_NNMT_dat_text$label == "0", "p_value_parse"] <- "0"

ENCODE_NNMT_RNAseq_plot_melt$group <- factor(ENCODE_NNMT_RNAseq_plot_melt$group, levels = c("H3K36_random", "H3K36_NNMT", "H3K27_random", "H3K27_NNMT", "H3K9_random", "H3K9_NNMT", "H3K4_random", "H3K4_NNMT"))

pdf("graphics/ENCODE_NNMT_allmarks_nopval.pdf",
    width = 3.75,
    height = 5)

ggplot(aes(x = signal, y = group, group = location), data = ENCODE_NNMT_RNAseq_plot_melt ) + 
  geom_boxplot(aes(group = group), fill = c("grey",
                                            "darkseagreen1",
                                            "grey",
                                            "darkorchid1",
                                            "grey",
                                            "magenta",
                                            "grey",
                                            "darkolivegreen3",
                                            "grey",
                                            "darkseagreen1",
                                            "grey",
                                            "darkorchid1",
                                            "grey",
                                            "magenta",
                                            "grey",
                                            "darkolivegreen3",
                                            "grey",
                                            "magenta"), 
               notch = TRUE, 
               outlier.shape = NA,
               lwd = 0.3,
               fatten = 0.7) +
  theme_classic() +
  geom_vline(xintercept = 0,
             linetype = "dashed",
             colour = "grey") +
  xlab("Linear model t-value") + 
  ylab(NULL) +
  facet_wrap(~ location, ncol = 1, labeller = labeller(location = factet_labs), scales = "free") + 
  coord_cartesian(xlim = c(-3, 4.5)) + 
  theme(strip.text = element_text(size = 12, color = "white", face = "bold"),
        strip.background = element_rect(fill = "black"),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.title.x = element_text(size = 10, color = "black", face = "bold")) +
  geom_text(aes(x = 3, group = location, label = Freq), data = ENCODE_NNMT_samplesizes, size = 3, hjust = 0)

dev.off()

pdf("graphics/ENCODE_NNMT_allmarks_pval.pdf",
    width = 3.75,
    height = 5)

ggplot(aes(x = signal, y = group, group = location), data = ENCODE_NNMT_RNAseq_plot_melt ) + 
  geom_boxplot(aes(group = group), fill = c("grey",
                                            "darkseagreen1",
                                            "grey",
                                            "darkorchid1",
                                            "grey",
                                            "magenta",
                                            "grey",
                                            "darkolivegreen3",
                                            "grey",
                                            "darkseagreen1",
                                            "grey",
                                            "darkorchid1",
                                            "grey",
                                            "magenta",
                                            "grey",
                                            "darkolivegreen3",
                                            "grey",
                                            "magenta"), 
               notch = TRUE, 
               outlier.shape = NA,
               lwd = 0.3,
               fatten = 0.7) +
  theme_classic() +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey") +
  xlab("Linear model t-value") + 
  ylab(NULL) +
  facet_wrap(~ location, ncol = 1, labeller = labeller(location = factet_labs), scales = "free") + 
  coord_cartesian(xlim = c(-3, 4.5)) + 
  theme(strip.text = element_text(size = 12, color = "white", face = "bold"),
        strip.background = element_rect(fill = "black"),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.title.x = element_text(size = 10, color = "black", face = "bold")) +
  geom_text(aes(x = 3, group = location, label = Freq), data = ENCODE_NNMT_samplesizes, size = 3, hjust = 0)+
  geom_text(aes(x = 3, group = location, label = p_value_parse), data = ENCODE_NNMT_dat_text, parse = TRUE)


dev.off()

#### ENCODE H3K4me3 WIDTH Fig S8A Fig S8B ####

alltissue_peak_widths <- lapply(ENCODEtissuetypes, function(tissue){
 
  promoter.peak.widths(thistissue = tissue,
                       chiptarget = "H3K4me3")
  
})

names(alltissue_peak_widths) <- ENCODEtissuetypes

# for constructing table, will be limited to promoters with peaks for all? no. but could limit to promoters with peaks in at least half the samples?

promoters_present <-  lapply(alltissue_peak_widths, function(smalllist){
  
  lapply(smalllist, function(x){x[, "gene_promoter_entrez_id"]})
  
})

unique_promoters <- unique(unlist(promoters_present))

allsample_promoter_width_list <- lapply(alltissue_peak_widths, function(smalllist){
  
  tempdata <- sapply(unique_promoters, function(x){
    
    sapply(smalllist, function(y){
      
      y[y$gene_promoter_entrez_id == x, "peak_width"]
      
    })
    
  })
  
  tempdatadf <- data.frame(tempdata)
  
  newdatadf <- apply(tempdatadf, 2, function(x){
    
    sapply(x, function(y){
      
      if(length(y) == 0){
        y <- NA
      } else {
        y <- y
      }
      
    })
    
  })
  
  return(newdatadf)
  
})

allsample_promoter_width_df <- do.call(rbind, allsample_promoter_width_list)

# filter to samples present in RNA data
allsample_promoter_width_df <- data.frame(allsample_promoter_width_df[str_remove(row.names(allsample_promoter_width_df), pattern = "-H3K4me3") %in% colnames(normalised_counts), ])

promoters_contain_NA <- apply(allsample_promoter_width_df, 2, function(x){
  
  tempvec <- sapply(x, function(y){
    is.na(y)
  })
  
  any(sapply(tempvec, isTRUE))
  
})

allsample_completepromoter_width_df <- allsample_promoter_width_df[, !promoters_contain_NA]

allsample_completepromoter_width_df[, "PEMT"] <- normalised_counts["ENSG00000133027", match(str_remove(row.names(allsample_completepromoter_width_df), pattern = "-H3K4me3"), colnames(normalised_counts))]

allsample_completepromoter_width_df[, "tissue"] <- str_remove(row.names(allsample_completepromoter_width_df), pattern = "[A-Z0-9]{11}-H3K4me3-")

allsample_completepromoter_width_df[, "donor"] <- str_remove(row.names(allsample_completepromoter_width_df), pattern = "-H3K4me3-.*")

allsample_completepromoter_width_df[, "PEMT_resid"] <- residuals(lmer(formula = log10(PEMT) ~ tissue + (1|donor), data = allsample_completepromoter_width_df))

ENCODE_H3K4_width_signal_model <- lapply(colnames(allsample_completepromoter_width_df)[1:(ncol(allsample_completepromoter_width_df) - 4)], function(x){
  
  tryCatch(expr = {
    templm <- lmer(formula = as.formula(paste0(x, " ~ PEMT_resid + tissue + (1|donor)")), data = allsample_completepromoter_width_df)
    templmsum <- summary(templm)
    
    t <- templmsum$coefficients["PEMT_resid", "t value"]
    
    output_vec <- c(t)
    names(output_vec) <- c("t_value")
    
    return(output_vec)
    
  }, error = function(e){
    
    output_vec <- c(NA)
    names(output_vec) <- c("t_value")
    
    return(output_vec)
    
  }
  )
  
})

names(ENCODE_H3K4_width_signal_model) <- colnames(allsample_completepromoter_width_df)[1:(ncol(allsample_completepromoter_width_df)-4)]

saveRDS(ENCODE_H3K4_width_signal_model, "output/ENCODE_H3K4_width_signal_model.rds")
# ENCODE_H3K4_width_signal_model <- readRDS("output/ENCODE_H3K4_width_signal_model.rds")

# what percentage is negatively correlated? 99.3%
sum(unlist(ENCODE_H3K4_width_signal_model) < 0) / length(ENCODE_H3K4_width_signal_model)

#### H3K4 promoter width random distribution ####

matrix_for_width_distribution <- allsample_completepromoter_width_df

for(i in 1:length(unlist(sample1000))){

  matrix_for_width_distribution[, unlist(sample1000)[i]] <- normalised_counts[unlist(sample1000)[i], match(str_remove(row.names(matrix_for_width_distribution), pattern = "H3.*-"), colnames(normalised_counts))]
  
  temp_lm <- lmer(formula = paste0("log10(", unlist(sample1000)[i], "+ 0.001) ~ tissue + (1|donor)"), data = matrix_for_width_distribution)
  temp_resid <- residuals(temp_lm)
  
  matrix_for_width_distribution[, paste0(unlist(sample1000)[i], "_resid")] <- temp_resid
  
}

ENCODE_H3K4promoter_width_randomsample_t_values <- sapply(unlist(sample1000), function(thisgene){
  
  sapply(colnames(matrix_for_width_distribution)[1:(ncol(matrix_for_width_distribution) - 2 - (2 * length(unlist(sample1000))))], function(x){
    
    tryCatch(expr = {
      
      glm <- lmer(formula = as.formula(paste0(x, " ~ ", thisgene, "_resid + tissue + (1|donor)")), data = matrix_for_width_distribution)
      glmsum <- summary(glm)
      
      t <- glmsum$coefficients[paste0(thisgene, "_resid"), "t value"]
      
      output_vec <- c(t)
      names(output_vec) <- c("t_value")
      
      return(output_vec)
      
    }, error = function(e){
      
      output_vec <- c(NA)
      names(output_vec) <- c("t_value")
      
      return(output_vec)
      
    }
    )
    
  })
  
})

colnames(ENCODE_H3K4promoter_width_randomsample_t_values) <- unlist(sample1000)

saveRDS(ENCODE_H3K4promoter_width_randomsample_t_values, "output/ENCODE_noTC_H3K4promoter_width_randomsample_t_values.rds")
ENCODE_H3K4promoter_width_randomsample_t_values <- readRDS("output/ENCODE_noTC_H3K4promoter_width_randomsample_t_values.rds")

width_match_vec <- intersect(names(ENCODE_H3K4_width_signal_model),
                             str_remove(row.names(ENCODE_H3K4promoter_width_randomsample_t_values), "\\.t_value"))

width_plot_list <- list(width = unlist(ENCODE_H3K4_width_signal_model[width_match_vec]),
                        random = rowMeans(ENCODE_H3K4promoter_width_randomsample_t_values[paste0(width_match_vec, ".t_value"), ]))

wilcox.test(width_plot_list[[1]], width_plot_list[[2]], paired = TRUE)$p.value

width_plot_df <- melt(width_plot_list)

width_plot_copy <- width_plot_df
width_plot_copy[width_plot_copy$L1 == "width", "L1"] <- "PEMT"
colnames(width_plot_copy) <- c("t_value", "group")

write.table(width_plot_copy,
            file = "plot_data/Fig S8/Fig_S8B_plotdata.txt",
            row.names = FALSE,
            col.names = TRUE,
            sep = "\t")

pdf("graphics/ENCODE_H3K4me3_width_model.pdf",
    width = 2,
    height = 2.2)

ggplot(data = width_plot_df, aes(y = L1, x = value)) + 
  geom_boxplot(fill = c("grey", "darkolivegreen1"),
               outlier.shape = NA,
               notch = TRUE,
               lwd = 0.3,
               fatten = 0.7) + 
  theme_classic2() + 
  theme(axis.title.y = element_blank(),
        axis.title.x = element_text(size = 10),
        axis.text.y = element_blank(),
        axis.text.x = element_text(colour = "black",
                                   size = 8)) +
  geom_vline(xintercept = 0,
             linetype = "dashed",
             colour = "grey") + 
  coord_cartesian(xlim = c(-3, 1)) + 
  xlab("Linear model t-value")

dev.off()

#### correlate width measure to signal Fig S8A ####

# generate matrix of H3K4 signal

H3K4promoters_computematrixfiles <- list.files("input/ENCODE_H3K4_computematrix/")[str_detect(list.files("input/ENCODE_H3K4_computematrix/"), "promoter")]

H3K4promoters_signalmatrix <- lapply(H3K4promoters_computematrixfiles, function(x){
  
  tempfile <- read.table(paste0("input/ENCODE_H3K4_computematrix/", x),  skip = 1)
  
})

names(H3K4promoters_signalmatrix) <- str_remove(H3K4promoters_computematrixfiles, pattern = paste0("_log2.bw.*\\.gz$"))

H3K4promoters_signalmatrix_df <- sapply(H3K4promoters_signalmatrix, function(x){
  
  x$V7
  
})

row.names(H3K4promoters_signalmatrix_df) <- H3K4promoters_signalmatrix[[1]]$V4

# restrict to samples in the metadata file
H3K4promoters_signalmatrix_df <- H3K4promoters_signalmatrix_df[, sapply(colnames(H3K4promoters_signalmatrix_df), function(x){
  str_detect(string = paste(chipmetadata_df$Files, sep = '', collapse = ''), pattern = x)}
)]

H3K4promoters_signalmatrix_tdf <- data.frame(t(H3K4promoters_signalmatrix_df))
row.names(H3K4promoters_signalmatrix_tdf) <- str_replace_all(paste0(sapply(1:nrow(H3K4promoters_signalmatrix_tdf), function(i){chipmetadata_df[str_detect(chipmetadata_df$Files, row.names(H3K4promoters_signalmatrix_tdf)[i]), "Donor"]}),
                                                                    "-",
                                                                    sapply(1:nrow(H3K4promoters_signalmatrix_tdf), function(i){chipmetadata_df[str_detect(chipmetadata_df$Files, row.names(H3K4promoters_signalmatrix_tdf)[i]), "Biosample.term.name"]})), pattern = " ", replacement = "_")

matching_vec <- colnames(allsample_completepromoter_width_df[1, 1:(ncol(allsample_completepromoter_width_df) - 4)])[match(colnames(H3K4promoters_signalmatrix_tdf[1, ]), colnames(allsample_completepromoter_width_df[1, 1:(ncol(allsample_completepromoter_width_df) - 4)]))]
matching_vec <- matching_vec[!is.na(matching_vec)]

allsample_width_match <- allsample_completepromoter_width_df[, matching_vec]
H3K4promoters_signal_match <- H3K4promoters_signalmatrix_tdf[, matching_vec]

row.names(allsample_width_match) <- str_remove(row.names(allsample_width_match), "H3K4me3-")
allsample_width_match <- allsample_width_match[match(row.names(H3K4promoters_signal_match), row.names(allsample_width_match)), ]

allsites_cor <- sapply(1:ncol(allsample_width_match), function(i){
  
  width_vec <- allsample_width_match[, i]
  signal_vec <- H3K4promoters_signal_match[, i]
  
  cor.test(width_vec, signal_vec)$estimate
  
})

width_signal_example_plot <- data.frame(width = allsample_width_match[, 1], 
                                        signal = H3K4promoters_signal_match[, 1])

write.table(width_signal_example_plot,
            file = "plot_data/Fig S8/Fig_S8A_plotdata.txt",
            row.names = FALSE,
            col.names = TRUE,
            sep = "\t")

cor.test(width_signal_example_plot$width, width_signal_example_plot$signal)

pdf("graphics/H3K4me3_width_vs_signal_example.pdf",
    width = 2.2,
    height = 2.2)

ggplot(data = width_signal_example_plot, aes(x = width, y = signal)) + 
  geom_point() + 
  geom_smooth(se = FALSE, 
              method = "lm") + 
  theme_classic2() + 
  theme(axis.text = element_text(colour = "black",
                                 size = 10),
        axis.title = element_text(size = 8),
        plot.margin = margin(t = 5, 
                             l = 5,
                             b = 5,
                             r = 10),
        axis.title.y = element_text(hjust = 0.8)) + 
  xlab(substitute("Promoter H3K4me3 peak width (bp)")) + 
  ylab(substitute("H3K4me3 signal ("~log[2]~"F.C. over control)")) +
  annotate(geom = "text",
           x = 2000, 
           y = 0.5, 
           label = substitute(italic(r)~"= 0.788"))

dev.off()

median(allsites_cor)
max(allsites_cor)
min(allsites_cor)
quantile(allsites_cor, 0.25)
quantile(allsites_cor, 0.75)
quantile(allsites_cor, 0.50)

# range of correlation in individual samples is 0.508 - 0.608, mean of 0.544

#### NNMT / ChIP RELATIONSHIP TO EXPRESSION Fig S8D Fig S8E Fig S8F ####

ENCODE_H3K4_promoter_signal_for_lm <- ENCODE_H3K4_promoter_signal

ENCODE_H3K4_promoter_signal_for_lm[, "tissue"] <- str_remove(row.names(ENCODE_H3K4_promoter_signal), "^ENCD.*-")
ENCODE_H3K4_promoter_signal_for_lm[, "donor"] <- str_remove(row.names(ENCODE_H3K4_promoter_signal), "-.*$")

ENCODE_H3K4_promoter_signal_for_lm[, "PEMT"] <- normalised_counts[PEMT_ensembl, match(row.names(ENCODE_H3K4_promoter_signal), colnames(normalised_counts))]

# exclude tissue types with fewer than three instances as cannot reasonably calculate residuals for those samples
ENCODE_H3K4_promoter_signal_for_lm <- ENCODE_H3K4_promoter_signal_for_lm[!(ENCODE_H3K4_promoter_signal_for_lm$tissue %in% names(table(ENCODE_H3K4_promoter_signal_for_lm[, "tissue"]))[table(ENCODE_H3K4_promoter_signal_for_lm[, "tissue"]) < 3]), ]

PEMT_lm <- lmer(formula = as.formula(paste0("log10(PEMT)", " ~  + tissue + (1|donor)")), data = ENCODE_H3K4_promoter_signal_for_lm)
ENCODE_H3K4_promoter_signal_for_lm[, "PEMT_resid"] <- residuals(PEMT_lm)

ENCODE_H3K4_promoter_signal_for_lm <- ENCODE_H3K4_promoter_signal_for_lm[, colnames(ENCODE_H3K4_promoter_signal_for_lm) %in% c(ENCODE_H3K4_promoter_signal_sites, "tissue", "donor", "PEMT", "PEMT_resid")]

ENCODE_PEMT_H3K4_promoter_signal_model <- lapply(colnames(ENCODE_H3K4_promoter_signal_for_lm)[1:(ncol(ENCODE_H3K4_promoter_signal_for_lm) - 4)], function(x){
  
  tryCatch(expr = {
    templm <- lmer(formula = as.formula(paste0(x, " ~ PEMT_resid + tissue + (1|donor)")), data = ENCODE_H3K4_promoter_signal_for_lm)
    templmsum <- summary(templm)
    
    t <- templmsum$coefficients["PEMT_resid", "t value"]
    
    output_vec <- c(t)
    names(output_vec) <- c("t_value")
    
    return(output_vec)
    
  }, error = function(e){
    
    output_vec <- c(NA)
    names(output_vec) <- c("t_value")
    
    return(output_vec)
    
  }
  )
  
})

names(ENCODE_PEMT_H3K4_promoter_signal_model) <- colnames(ENCODE_H3K4_promoter_signal_for_lm)[1:(ncol(ENCODE_H3K4_promoter_signal_for_lm)-4)]

saveRDS(ENCODE_PEMT_H3K4_promoter_signal_model, "output/ENCODE_PEMT_H3K4_promoter_signal_model.rds")
# ENCODE_PEMT_H3K4_promoter_signal_model <- readRDS("output/ENCODE_PEMT_H3K4_promoter_signal_model.rds")

ensembl_to_entrez_table <- getBM(attributes = c("ensembl_gene_id", "entrezgene_id"),
                                 values = row.names(normalised_counts),
                                 filters = "ensembl_gene_id",
                                 mart = ensembl ,
)

normalised_counts_entrez <- normalised_counts
row.names(normalised_counts_entrez) <- ensembl_to_entrez_table[match(row.names(normalised_counts), ensembl_to_entrez_table$ensembl_gene_id), "entrezgene_id"]

normalised_counts_entrez <- normalised_counts_entrez[!is.na(row.names(normalised_counts_entrez)), ]

ENCODE_H3K4_expression_matrix_for_lm <- normalised_counts_entrez[match(ENCODE_H3K4_promoter_signal_sites , paste0("X", row.names(normalised_counts_entrez))), ]
ENCODE_H3K4_expression_matrix_for_lm <- data.frame(t(ENCODE_H3K4_expression_matrix_for_lm[!is.na(row.names(ENCODE_H3K4_expression_matrix_for_lm)), colnames(ENCODE_H3K4_expression_matrix_for_lm) %in% row.names(ENCODE_H3K4_promoter_signal_for_lm)]))

ENCODE_H3K4_expression_matrix_for_lm <- ENCODE_H3K4_expression_matrix_for_lm[, colnames(ENCODE_H3K4_expression_matrix_for_lm) %in% ENCODE_H3K4_promoter_signal_sites]

# exclude any genes which don't have expression in some samples
ENCODE_H3K4_expression_matrix_for_lm <- ENCODE_H3K4_expression_matrix_for_lm[, !apply(ENCODE_H3K4_expression_matrix_for_lm, 2, function(x){any(x < 0.001)})]

ENCODE_H3K4_expression_matrix_for_lm[, "tissue"] <- str_remove(row.names(ENCODE_H3K4_expression_matrix_for_lm ), "^ENCD.*-")
ENCODE_H3K4_expression_matrix_for_lm[, "donor"] <- str_remove(row.names(ENCODE_H3K4_expression_matrix_for_lm ), "-.*$")

ENCODE_H3K4_expression_matrix_for_lm[, "PEMT"] <- normalised_counts[PEMT_ensembl, match(row.names(ENCODE_H3K4_expression_matrix_for_lm), colnames(normalised_counts))]

# exclude tissue types with fewer than three instances as cannot reasonably calculate residuals for those samples
ENCODE_H3K4_expression_matrix_for_lm <- ENCODE_H3K4_expression_matrix_for_lm[!(ENCODE_H3K4_expression_matrix_for_lm$tissue %in% names(table(ENCODE_H3K4_expression_matrix_for_lm[, "tissue"]))[table(ENCODE_H3K4_expression_matrix_for_lm[, "tissue"]) < 3]), ]

PEMT_lm <- lmer(formula = as.formula(paste0("log10(PEMT)", " ~  + tissue + (1|donor)")), data = ENCODE_H3K4_expression_matrix_for_lm)
ENCODE_H3K4_expression_matrix_for_lm[, "PEMT_resid"] <- residuals(PEMT_lm)

ENCODE_PEMT_promoter_expression_model <- lapply(colnames(ENCODE_H3K4_expression_matrix_for_lm)[1:(ncol(ENCODE_H3K4_expression_matrix_for_lm) - 4)], function(x){
  
  tryCatch(expr = {
    templm <- lmer(formula = as.formula(paste0(x, " ~ PEMT_resid + tissue + (1|donor)")), data = ENCODE_H3K4_expression_matrix_for_lm)
    templmsum <- summary(templm)
    
    t <- templmsum$coefficients["PEMT_resid", "t value"]
    
    output_vec <- c(t)
    names(output_vec) <- c("t_value")
    
    return(output_vec)
    
  }, error = function(e){
    
    output_vec <- c(NA)
    names(output_vec) <- c("t_value")
    
    return(output_vec)
    
  }
  )
  
})

names(ENCODE_PEMT_promoter_expression_model) <- colnames(ENCODE_H3K4_expression_matrix_for_lm)[1:(ncol(ENCODE_H3K4_expression_matrix_for_lm)-4)]

saveRDS(ENCODE_PEMT_promoter_expression_model, "output/ENCODE_PEMT_promoter_expression_model.rds")
# ENCODE_PEMT_promoter_expression_model <- readRDS("output/ENCODE_PEMT_promoter_expression_model.rds")


matching_sites <- names(ENCODE_PEMT_H3K4_promoter_model)[match(names(ENCODE_PEMT_promoter_expression_model), names(ENCODE_PEMT_H3K4_promoter_model))]
matching_sites <- matching_sites[!is.na(matching_sites)]

ENCODE_H3K4_promoter_expression_vs_signal_models <- data.frame(expression_t_value = unlist(ENCODE_PEMT_promoter_expression_model[matching_sites]),
                                                               signal_t_value = unlist(ENCODE_PEMT_H3K4_promoter_model[matching_sites]))

melted_expression_vs_signal <- melt(ENCODE_H3K4_promoter_expression_vs_signal_models )

unique(melted_expression_vs_signal$variable)
melted_expression_vs_signal$variable <- factor(melted_expression_vs_signal$variable, levels = c("signal_t_value", "expression_t_value"))

write.table(melted_expression_vs_signal,
            file = "plot_data/Fig S8/Fig_S8D_promoterpanel_plotdata.txt",
            row.names = FALSE,
            col.names = TRUE,
            sep = "\t")

pdf("graphics/ENCODE_H3K4promoter_PEMTmodels_expression_signal.pdf",
    height = 2.5,
    width = 2.5)

ggplot(melted_expression_vs_signal, aes(x = variable, y = value)) + 
  geom_boxplot(fill = c("darkolivegreen3", "orange"),
               outlier.shape = NA) + 
  theme_classic() + 
  theme(axis.title.x = element_blank(),
        axis.text = element_text(size = 10,
                                 colour = "black"),
        axis.title.y = element_text(size = 10,
                                    colour = "black")) + 
  geom_hline(yintercept = 0 , linetype = "dashed", colour = "grey") + 
  coord_cartesian(ylim = c(-5, 5)) + 
  ylab(substitute(italic(PEMT)~"expression t-value")) +
  xlab("") + 
  scale_x_discrete(labels = c("signal_t_value"  = "Promoter\nH3K4 signal", "expression_t_value" = "Gene\nexpression"))

dev.off()

#### H3K9me3 at genebodies to expression ####

ENCODE_H3K9_genebody_signal_for_lm <- ENCODE_H3K9_genebody_signal

ENCODE_H3K9_genebody_signal_for_lm[, "tissue"] <- str_remove(row.names(ENCODE_H3K9_genebody_signal), "^ENCD.*-")
ENCODE_H3K9_genebody_signal_for_lm[, "donor"] <- str_remove(row.names(ENCODE_H3K9_genebody_signal), "-.*$")

ENCODE_H3K9_genebody_signal_for_lm[, "PEMT"] <- normalised_counts[PEMT_ensembl, match(row.names(ENCODE_H3K9_genebody_signal), colnames(normalised_counts))]

# exclude tissue types with fewer than three instances as cannot reasonably calculate residuals for those samples
ENCODE_H3K9_genebody_signal_for_lm <- ENCODE_H3K9_genebody_signal_for_lm[!(ENCODE_H3K9_genebody_signal_for_lm$tissue %in% names(table(ENCODE_H3K9_genebody_signal_for_lm[, "tissue"]))[table(ENCODE_H3K9_genebody_signal_for_lm[, "tissue"]) < 3]), ]

PEMT_lm <- lmer(formula = as.formula(paste0("log10(PEMT)", " ~  + tissue + (1|donor)")), data = ENCODE_H3K9_genebody_signal_for_lm)
ENCODE_H3K9_genebody_signal_for_lm[, "PEMT_resid"] <- residuals(PEMT_lm)

ENCODE_H3K9_genebody_signal_for_lm <- ENCODE_H3K9_genebody_signal_for_lm[, colnames(ENCODE_H3K9_genebody_signal_for_lm) %in% c(ENCODE_H3K9_genebody_signal_sites, "tissue", "donor", "PEMT", "PEMT_resid")]

ENCODE_PEMT_H3K9_genebody_signal_model <- lapply(colnames(ENCODE_H3K9_genebody_signal_for_lm)[1:(ncol(ENCODE_H3K9_genebody_signal_for_lm) - 4)], function(x){
  
  tryCatch(expr = {
    
    templm <- lmer(formula = as.formula(paste0(x, " ~ PEMT_resid + tissue + (1|donor)")), data = ENCODE_H3K9_genebody_signal_for_lm)
    templmsum <- summary(templm)
    
    t <- templmsum$coefficients["PEMT_resid", "t value"]
    
    output_vec <- c(t)
    names(output_vec) <- c("t_value")
    
    return(output_vec)
    
  }, error = function(e){
    
    output_vec <- c(NA)
    names(output_vec) <- c("t_value")
    
    return(output_vec)
    
  }
  )
  
})

names(ENCODE_PEMT_H3K9_genebody_signal_model) <- colnames(ENCODE_H3K9_genebody_signal_for_lm)[1:(ncol(ENCODE_H3K9_genebody_signal_for_lm)-4)]
# 
# ensembl_to_entrez_table <- getBM(attributes = c("ensembl_gene_id", "entrezgene_id"),
#                                  values = row.names(normalised_counts),
#                                  filters = "ensembl_gene_id",
#                                  mart = ensembl ,
# )
# 
# normalised_counts_entrez <- normalised_counts
# row.names(normalised_counts_entrez) <- ensembl_to_entrez_table[match(row.names(normalised_counts), ensembl_to_entrez_table$ensembl_gene_id), "entrezgene_id"]

ENCODE_H3K9_expression_matrix_for_lm <- normalised_counts_entrez[match(ENCODE_H3K9_genebody_signal_sites , paste0("X", row.names(normalised_counts_entrez))), ]
ENCODE_H3K9_expression_matrix_for_lm <- data.frame(t(ENCODE_H3K9_expression_matrix_for_lm[!is.na(row.names(ENCODE_H3K9_expression_matrix_for_lm)), colnames(ENCODE_H3K9_expression_matrix_for_lm) %in% row.names(ENCODE_H3K9_genebody_signal_for_lm)]))

ENCODE_H3K9_expression_matrix_for_lm <- ENCODE_H3K9_expression_matrix_for_lm[, colnames(ENCODE_H3K9_expression_matrix_for_lm) %in% ENCODE_H3K9_genebody_signal_sites]

# exclude any genes which don't have expression in some samples
ENCODE_H3K9_expression_matrix_for_lm <- ENCODE_H3K9_expression_matrix_for_lm[, !apply(ENCODE_H3K9_expression_matrix_for_lm, 2, function(x){any(x < 0.001)})]

ENCODE_H3K9_expression_matrix_for_lm[, "tissue"] <- str_remove(row.names(ENCODE_H3K9_expression_matrix_for_lm ), "^ENCD.*-")
ENCODE_H3K9_expression_matrix_for_lm[, "donor"] <- str_remove(row.names(ENCODE_H3K9_expression_matrix_for_lm ), "-.*$")

ENCODE_H3K9_expression_matrix_for_lm[, "PEMT"] <- normalised_counts[PEMT_ensembl, match(row.names(ENCODE_H3K9_expression_matrix_for_lm), colnames(normalised_counts))]

# exclude tissue types with fewer than three instances as cannot reasonably calculate residuals for those samples
ENCODE_H3K9_expression_matrix_for_lm <- ENCODE_H3K9_expression_matrix_for_lm[!(ENCODE_H3K9_expression_matrix_for_lm$tissue %in% names(table(ENCODE_H3K9_expression_matrix_for_lm[, "tissue"]))[table(ENCODE_H3K9_expression_matrix_for_lm[, "tissue"]) < 3]), ]

PEMT_lm <- lmer(formula = as.formula(paste0("log10(PEMT)", " ~  + tissue + (1|donor)")), data = ENCODE_H3K9_expression_matrix_for_lm)
ENCODE_H3K9_expression_matrix_for_lm[, "PEMT_resid"] <- residuals(PEMT_lm)

ENCODE_PEMT_genebody_expression_model <- lapply(colnames(ENCODE_H3K9_expression_matrix_for_lm)[1:(ncol(ENCODE_H3K9_expression_matrix_for_lm) - 4)], function(x){
  
  tryCatch(expr = {
    templm <- lmer(formula = as.formula(paste0(x, " ~ PEMT_resid + tissue + (1|donor)")), data = ENCODE_H3K9_expression_matrix_for_lm)
    templmsum <- summary(templm)
    
    t <- templmsum$coefficients["PEMT_resid", "t value"]
    
    output_vec <- c(t)
    names(output_vec) <- c("t_value")
    
    return(output_vec)
    
  }, error = function(e){
    
    output_vec <- c(NA)
    names(output_vec) <- c("t_value")
    
    return(output_vec)
    
  }
  )
  
})

names(ENCODE_PEMT_genebody_expression_model) <- colnames(ENCODE_H3K9_expression_matrix_for_lm)[1:(ncol(ENCODE_H3K9_expression_matrix_for_lm)-4)]

matching_sites <- names(ENCODE_PEMT_H3K9_genebody_signal_model)[match(names(ENCODE_PEMT_genebody_expression_model), names(ENCODE_PEMT_H3K9_genebody_signal_model))]
matching_sites <- matching_sites[!is.na(matching_sites)]

ENCODE_H3K9_genebody_expression_vs_signal_models <- data.frame(expression_t_value = unlist(ENCODE_PEMT_genebody_expression_model[matching_sites]),
                                                               signal_t_value = unlist(ENCODE_PEMT_H3K9_genebody_signal_model[matching_sites]))

melted_expression_vs_signal <- melt(ENCODE_H3K9_genebody_expression_vs_signal_models )

unique(melted_expression_vs_signal$variable)
melted_expression_vs_signal$variable <- factor(melted_expression_vs_signal$variable, levels = c("signal_t_value", "expression_t_value"))

write.table(melted_expression_vs_signal,
            file = "plot_data/Fig S8/Fig_S8D_genebodypanel_plotdata.txt",
            row.names = FALSE,
            col.names = TRUE,
            sep = "\t")

pdf("graphics/ENCODE_H3K9genebody_PEMTmodels_expression_signal.pdf",
    height = 2.5,
    width = 2.5)

ggplot(melted_expression_vs_signal, aes(x = variable, y = value)) + 
  geom_boxplot(fill = c("darkolivegreen3", "orange"),
               outlier.shape = NA) + 
  theme_classic() + 
  theme(axis.title.x = element_blank(),
        axis.text = element_text(size = 10,
                                 colour = "black"),
        axis.title.y = element_text(size = 10,
                                    colour = "black")) + 
  geom_hline(yintercept = 0 , linetype = "dashed", colour = "grey") + 
  coord_cartesian(ylim = c(-5, 5)) + 
  ylab(substitute(italic(PEMT)~"expression t-value")) +
  xlab("") + 
  scale_x_discrete(labels = c("signal_t_value"  = "genebody\nH3K9 signal", "expression_t_value" = "Gene\nexpression"))

dev.off()

#### correlate H3K4me3 promoter signal to expression ####

H3K4_promoter_signal_expression_tvals <- sapply(1:(ncol(ENCODE_H3K4_expression_matrix_for_lm) - 4), function(i){

  tempgene <- colnames(ENCODE_H3K4_expression_matrix_for_lm)[i]

  tempmatrix_for_lm <- data.frame(signal_resid = residuals(lmer(formula = as.formula(paste0(tempgene, " ~ tissue + (1|donor)")), data = ENCODE_H3K4_promoter_signal_for_lm)))

  tempmatrix_for_lm[, "expression"] <- ENCODE_H3K4_expression_matrix_for_lm[row.names(tempmatrix_for_lm), tempgene]
  tempmatrix_for_lm[, "tissue"] <- ENCODE_H3K4_promoter_signal_for_lm[row.names(tempmatrix_for_lm), "tissue"]
  tempmatrix_for_lm[, "donor"] <- ENCODE_H3K4_promoter_signal_for_lm[row.names(tempmatrix_for_lm), "donor"]
  
  tempmodelfit <- lmer(log10(expression + 1) ~ signal_resid + tissue + (1|donor), data = tempmatrix_for_lm)
  
  tempsumm <- summary(tempmodelfit)
  
  return(tempsumm$coefficients["signal_resid", "t value"])
  
})

saveRDS(H3K4_promoter_signal_expression_tvals, "output/H3K4_promoter_signal_expression_tvals.rds")
H3K4_promoter_signal_expression_tvals <- readRDS("output/H3K4_promoter_signal_expression_tvals.rds")

#### correlate H3K4me3 promoter peak width to expression ####

ENCODE_H3K4_expression_for_width_cor <- ENCODE_H3K4_expression_matrix_for_lm[, colnames(ENCODE_H3K4_expression_matrix_for_lm) %in% colnames(matrix_for_width_distribution)]

H3K4_promoter_width_vs_expresion <- sapply(1:(ncol(ENCODE_H3K4_expression_for_width_cor) - 4), function(i){

  tempgene <- colnames(ENCODE_H3K4_expression_for_width_cor)[i]
  
  ENCODE_H3K4_expression_for_width_cor[, paste0(tempgene, "_widthresid")] <- residuals(lmer(formula = as.formula(paste0(tempgene, " ~ tissue + (1|donor)")), data = matrix_for_width_distribution))
  
  expression_model <- (lmer(formula = as.formula(paste0("log10(", tempgene, " +1) ~ ", tempgene,"_widthresid + tissue + (1|donor)")), data = ENCODE_H3K4_expression_for_width_cor))
  
  width_tvalue <- summary(expression_model)$coefficients[paste0(tempgene, "_widthresid"), "t value"]
  
  return(width_tvalue)
  
})

H3K4_promoter_width_vs_expresion_df <- list(promoters = data.frame(cor = H3K4_promoter_width_vs_expresion,
                                                  group = "H3K4me3"))

H3K4_promoter_width_vs_expresion_meltdf <- reshape2::melt( H3K4_promoter_width_vs_expresion_df )

write.table(H3K4_promoter_width_vs_expresion_meltdf[, c(1, 3:4)],
            file = "plot_data/Fig S8/Fig_S8F_plotdata.txt",
            row.names = FALSE,
            col.names = TRUE,
            sep = "\t")

H3K4widthfacet_labels <- c("Promoters")
names(H3K4widthfacet_labels) <- c("promoters")

pdf("graphics/ENCODE_H3K4promoter_width_vs_expression_corr.pdf",
    height = 3,
    width = 2)

ggplot(H3K4_promoter_width_vs_expresion_meltdf, aes(x = group, y = value)) + 
  geom_boxplot(fill = "darkolivegreen3",
               outlier.shape = NA,
               notch = TRUE)+ 
  facet_wrap(~ L1, labeller = labeller(L1 = H3K4widthfacet_labels)) +
  theme_classic() +
  theme(strip.text = element_text(size = 12, color = "white", face = "bold"),
        strip.background = element_rect(fill = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.title.x = element_blank(),
        panel.border = element_rect(colour = "black",
                                    fill = "transparent"),
        axis.title.y = element_text(size = 9, color = "black")) +
  geom_hline(yintercept = 0,
             linetype = "dashed",
             colour = "grey") + 
  geom_text(aes(y = 4.75, x = 1), 
            label = "n = 9512",
            size = 3) +
  coord_cartesian(ylim = c(-5, 5)) + 
  ylab("Expression ~ peak width linear model t-value")

dev.off()

#### correlate H3K9me3 promoter signal to expression ####

H3K9_genebody_signal_expression_tvals <- sapply(1:(ncol(ENCODE_H3K9_expression_matrix_for_lm) - 4), function(i){
  
  tempgene <- colnames(ENCODE_H3K9_expression_matrix_for_lm)[i]
  
  tempmatrix_for_lm <- data.frame(signal_resid = residuals(lmer(formula = as.formula(paste0(tempgene, " ~ tissue + (1|donor)")), data = ENCODE_H3K9_genebody_signal_for_lm)))
  
  tempmatrix_for_lm[, "expression"] <- ENCODE_H3K9_expression_matrix_for_lm[row.names(tempmatrix_for_lm), tempgene]
  tempmatrix_for_lm[, "tissue"] <- ENCODE_H3K9_genebody_signal_for_lm[row.names(tempmatrix_for_lm), "tissue"]
  tempmatrix_for_lm[, "donor"] <- ENCODE_H3K9_genebody_signal_for_lm[row.names(tempmatrix_for_lm), "donor"]
  
  tempmodelfit <- lmer(log10(expression + 1) ~ signal_resid + tissue + (1|donor), data = tempmatrix_for_lm)
  
  tempsumm <- summary(tempmodelfit)
  
  return(tempsumm$coefficients["signal_resid", "t value"])
  
})

saveRDS(H3K9_genebody_signal_expression_tvals, "output/H3K9_genebody_signal_expression_tvals.rds")
# H3K9_genebody_signal_expression_tvals <- readRDS("output/H3K9_genebody_signal_expression_tvals.rds")

#### plots for ENCODE related to expression ####

ENCODE_H3K4promoter_expression_vs_signal_plot <- list(promoters = list("H3K4me3" = H3K4_promoter_signal_expression_tvals))
ENCODE_H3K4promoter_expression_vs_signal_plot <- reshape2::melt(ENCODE_H3K4promoter_expression_vs_signal_plot)

write.table(ENCODE_H3K4promoter_expression_vs_signal_plot, 
            file = "plot_data/Fig S8/Fig_S8E_H3K4promoterpanel_plotdata.txt",
            row.names = FALSE,
            col.names = TRUE,
            sep = "\t")

H3K4facet_labels <- c("Promoters")
names(H3K4facet_labels) <- c("promoters")

pdf("graphics/ENCODE_H3K4_promoter_signal_vs_expression.pdf",
        width = 2,
    height = 3)

ggplot(data = ENCODE_H3K4promoter_expression_vs_signal_plot, aes(x = L2, y = value)) +
  geom_boxplot(fill = "darkolivegreen3",
               notch = TRUE,
               outlier.shape = NA) +
  facet_wrap(~ L1, labeller = labeller(L1 = H3K4facet_labels)) +
  theme_classic() +
  theme(strip.text = element_text(size = 12, color = "white", face = "bold"),
        strip.background = element_rect(fill = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.text.x = element_text(size = 10, color = "black"),
        panel.border = element_rect(colour = "black",
                                    fill = "transparent"),
        axis.title.y = element_text(size = 9, color = "black"),
        axis.title.x = element_text(size = 9, color = "black", face = "bold")) +
  geom_hline(yintercept = 0,
             linetype = "dashed",
             colour = "grey") +
  xlab(NULL) +
  geom_text(aes(y = 4.75, x = 1), 
            label = "n = 9873",
            size = 3) +
  coord_cartesian(ylim = c(-5, 5)) + 
  ylab("Expression ~ ChIP signal linear model t-value")

dev.off()

ENCODE_H3K9genebody_expression_vs_signal_plot <- list(genebodys = list("H3K9me3" = H3K9_genebody_signal_expression_tvals))
ENCODE_H3K9genebody_expression_vs_signal_plot <- reshape2::melt(ENCODE_H3K9genebody_expression_vs_signal_plot)

write.table(ENCODE_H3K9genebody_expression_vs_signal_plot, 
            file = "plot_data/Fig S8/Fig_S8E_H3K9genebodypanel_plotdata.txt",
            row.names = FALSE,
            col.names = TRUE,
            sep = "\t")

H3K9facet_labels <- c("Gene bodies")
names(H3K9facet_labels) <- c("genebodys")

pdf("graphics/ENCODE_H3K9_genebody_signal_vs_expression.pdf",
    width = 1.55,
    height = 3)

ggplot(data = ENCODE_H3K9genebody_expression_vs_signal_plot, aes(x = L2, y = value)) +
  geom_boxplot(fill = "magenta",
               notch = TRUE,
               outlier.shape = NA) +
  facet_wrap(~ L1, labeller = labeller(L1 = H3K9facet_labels)) +
  theme_classic() +
  theme(strip.text = element_text(size = 12, color = "white", face = "bold"),
        strip.background = element_rect(fill = "black"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = 10, color = "black"),
        panel.border = element_rect(colour = "black",
                                    fill = "transparent"),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 9, color = "black", face = "bold")) +
  geom_hline(yintercept = 0,
             linetype = "dashed",
             colour = "grey") +
  xlab(NULL) +
  geom_text(aes(y = 4.75, x = 1), 
            label = "n = 280",
            size = 3) +
  coord_cartesian(ylim = c(-5, 5))

dev.off()

# plot for PEMT signal

H3K9genebodymatching_sites <- names(ENCODE_PEMT_H3K9_genebody_signal_model)[match(names(ENCODE_PEMT_genebody_expression_model), names(ENCODE_PEMT_H3K9_genebody_signal_model))]
H3K9genebodymatching_sites <- H3K9genebodymatching_sites[!is.na(H3K9genebodymatching_sites)]

ENCODE_H3K9_genebody_expression_vs_signal_forplot <- list(gene_bodies = list("expression" = unlist(ENCODE_PEMT_genebody_expression_model[H3K9genebodymatching_sites]),
                                                               "signal" = unlist(ENCODE_PEMT_H3K9_genebody_signal_model[H3K9genebodymatching_sites])))

ENCODE_H3K9_genebody_expression_vs_signal_forplotmelt <- reshape2::melt(ENCODE_H3K9_genebody_expression_vs_signal_forplot)
colnames(ENCODE_H3K9_genebody_expression_vs_signal_forplotmelt) <- c("value", "group", "location")

ENCODE_H3K9_genebody_expression_vs_signal_forplotmelt[ENCODE_H3K9_genebody_expression_vs_signal_forplotmelt$group == "signal", "group"] <- "H3K9me3\nsignal"
ENCODE_H3K9_genebody_expression_vs_signal_forplotmelt[ENCODE_H3K9_genebody_expression_vs_signal_forplotmelt$group == "expression", "group"] <- "Gene\nexpression"

H3K9facet_labels <- c("Gene bodies")
names(H3K9facet_labels) <- c("gene_bodies")

ENCODE_H3K9_genebody_expression_vs_signal_forplotmelt$group <- factor(ENCODE_H3K9_genebody_expression_vs_signal_forplotmelt$group, levels = c("H3K9me3\nsignal", "Gene\nexpression"))

table(ENCODE_H3K9_genebody_expression_vs_signal_forplotmelt$group)
samplesizes <- table(paste0(ENCODE_H3K9_genebody_expression_vs_signal_forplotmelt$group, ENCODE_H3K9_genebody_expression_vs_signal_forplotmelt$location))

samplesizes_edit <- c("gene_bodies" = max(samplesizes[str_detect(names(samplesizes), "gene_bodies")]))

samplesizes_edit <- data.frame(cbind(samplesizes_edit,
                                     names(samplesizes_edit)))

names(samplesizes_edit) <- c("Freq", "location")

samplesizes_edit$Freq = paste0("n = ", samplesizes_edit$Freq)

pdf("graphics/ENCODE_H3K9me3_expression_by_PEMT_boxplot.pdf",
    width = 1.6, height = 3)

ggplot(aes(y = value, x = group, group = location), data = ENCODE_H3K9_genebody_expression_vs_signal_forplotmelt) + 
  geom_boxplot(aes(group = group),
               notch = TRUE,
               outlier.shape = NA,
               fill = c("magenta",
                        "red")) +
  facet_wrap(~ location, labeller = labeller(location = factet_labs)) +
  theme_classic() +
  theme(strip.text = element_text(size = 12, color = "white", face = "bold"),
        strip.background = element_rect(fill = "black"),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 10, color = "black"),
        panel.border = element_rect(colour = "black",
                                    fill = "transparent"),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_blank()) +
  geom_hline(yintercept = 0,
             linetype = "dashed",
             colour = "grey") +
  xlab("Response variable") + 
  geom_text(aes(x = 1.5, y = 3, group = location, label = Freq), data = samplesizes_edit) + 
  coord_cartesian(ylim = c(-4, 3.5)) + 
  scale_y_continuous(breaks = c(-3, -1.5, 0, 1.5, 3))

dev.off()

H3K4promotermatching_sites <- names(ENCODE_PEMT_H3K4_promoter_model)[match(names(ENCODE_PEMT_promoter_expression_model), names(ENCODE_PEMT_H3K4_promoter_model))]
H3K4promotermatching_sites <- H3K4promotermatching_sites[!is.na(H3K4promotermatching_sites)]

ENCODE_H3K4_promoter_expression_vs_signal_forplot <- list(promoters = list("expression" = unlist(ENCODE_PEMT_promoter_expression_model[H3K4promotermatching_sites]),
                                                                             "signal" = unlist(ENCODE_PEMT_H3K4_promoter_signal_model[H3K4promotermatching_sites])))

ENCODE_H3K4_promoter_expression_vs_signal_forplotmelt <- reshape2::melt(ENCODE_H3K4_promoter_expression_vs_signal_forplot)
colnames(ENCODE_H3K4_promoter_expression_vs_signal_forplotmelt) <- c("value", "group", "location")

ENCODE_H3K4_promoter_expression_vs_signal_forplotmelt[ENCODE_H3K4_promoter_expression_vs_signal_forplotmelt$group == "signal", "group"] <- "H3K4me3\nsignal"
ENCODE_H3K4_promoter_expression_vs_signal_forplotmelt[ENCODE_H3K4_promoter_expression_vs_signal_forplotmelt$group == "expression", "group"] <- "Gene\nexpression"

H3K4facet_labels <- c("Promoters")
names(H3K4facet_labels) <- c("promoters")

ENCODE_H3K4_promoter_expression_vs_signal_forplotmelt$group <- factor(ENCODE_H3K4_promoter_expression_vs_signal_forplotmelt$group, levels = c("H3K4me3\nsignal", "Gene\nexpression"))

table(ENCODE_H3K4_promoter_expression_vs_signal_forplotmelt$group)
samplesizes <- table(paste0(ENCODE_H3K4_promoter_expression_vs_signal_forplotmelt$group, ENCODE_H3K4_promoter_expression_vs_signal_forplotmelt$location))

samplesizes_edit <- c("promoters" = max(samplesizes[str_detect(names(samplesizes), "promoters")]))

samplesizes_edit <- data.frame(cbind(samplesizes_edit,
                                     names(samplesizes_edit)))

names(samplesizes_edit) <- c("Freq", "location")

samplesizes_edit$Freq = paste0("n = ", samplesizes_edit$Freq)

pdf("graphics/ENCODE_H3K4me3_expression_by_PEMT_boxplot.pdf",
    width = 2, height = 3)

ggplot(aes(y = value, x = group, group = location), data = ENCODE_H3K4_promoter_expression_vs_signal_forplotmelt) + 
  geom_boxplot(aes(group = group),
               notch = TRUE,
               outlier.shape = NA,
               fill = c("darkolivegreen3",
                        "red")) +
  facet_wrap(~ location, labeller = labeller(location = H3K4facet_labels)) +
  theme_classic() +
  theme(strip.text = element_text(size = 12, color = "white", face = "bold"),
        strip.background = element_rect(fill = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.text.x = element_text(size = 10, color = "black"),
        panel.border = element_rect(colour = "black",
                                    fill = "transparent"),
        axis.title.y = element_text(size = 9, color = "black"),
        axis.title.x = element_blank()) +
  geom_hline(yintercept = 0,
             linetype = "dashed",
             colour = "grey") +
  ylab(italic(PEMT)~" linear model t-value") + 
  xlab("Response variable") + 
  geom_text(aes(x = 1.5, y = 3, group = location, label = Freq), data = samplesizes_edit) + 
  coord_cartesian(ylim = c(-4, 3.5)) + 
  scale_y_continuous(breaks = c(-3, -1.5, 0, 1.5, 3))

dev.off()

#### NCI60 ChIP-seq for Fig 3C Fig 3D Fig S7 ####

#### H3K4me3 ####

#### combine narrow and broad peak options####

# point to directory with genrich output
my_directory <- "input/genrich_peak_calls/"

genrich_file_list <- list.files(my_directory)

# choose the narrowPeak files
genrich_file_list_narrowPeak <- genrich_file_list[str_detect(genrich_file_list, pattern = "\\.narrowPeak")]

# choose the H3K4me3
H3K4me3_file_list_narrowPeak <- genrich_file_list_narrowPeak[str_detect(genrich_file_list_narrowPeak, pattern = "H3K4me3")]

# loop through unique experiments to pool the two replicates

H3K4me3_ranges_list <- lapply(unique(str_remove(H3K4me3_file_list_narrowPeak, pattern = "_genrich.*")), function(x){
  
  uniqueexperiment <- str_remove(x, pattern = "_genrich_.*")
  
  matchingfiles <- H3K4me3_file_list_narrowPeak[str_detect(H3K4me3_file_list_narrowPeak, paste0(uniqueexperiment, "_"))]
  
  if(length(matchingfiles) != 2){
    message(paste0("Warning! ", uniqueexperiment, " has ", length(matchingfiles),"  matching files!"))
  }
  
  peaks_1 <- read.table(paste0(my_directory, matchingfiles[1]), header = FALSE)
  colnames(peaks_1)[1:3] <- c("chr", "start", "end")
  
  peaks_2 <- read.table(paste0(my_directory, matchingfiles[2]), header = FALSE)
  colnames(peaks_2)[1:3] <- c("chr", "start", "end")
  
  peaks_1_GR <- makeGRangesFromDataFrame(peaks_1)
  peaks_2_GR <- makeGRangesFromDataFrame(peaks_2)
  
  all_peaks_GR <- c(peaks_1_GR, peaks_2_GR)
  
  all_peaks_reduced_GR <- reduce(all_peaks_GR)
  
  return(all_peaks_reduced_GR)
  
})

names(H3K4me3_ranges_list) <- unique(str_remove(H3K4me3_file_list_narrowPeak, pattern = "_genrich.*"))

saveRDS(H3K4me3_ranges_list, "output/NCI60_H3K4me3_combinedpeaks_list.rds")
# H3K4me3_ranges_list <- readRDS("output/NCI60_H3K4me3_combinedpeaks_list.rds")

#### SELECT H3K4me3 REGIONS WITH CONSISTENT PEAKS ####

H3K4me3_gene_peak_counts <- region.peak.counts(NCI60_list = H3K4me3_ranges_list,
                                               region = "genes")

H3K4genes_in_40of60_number <- apply(H3K4me3_gene_peak_counts[, 5:ncol(H3K4me3_gene_peak_counts)], 1, function(x){
  
  sum(x > 0)
  
})

H3K4genes_in_40of60_selec <- H3K4genes_in_40of60_number > 39
H3K4genes_in_40of60 <- H3K4me3_gene_peak_counts[H3K4genes_in_40of60_selec, "gene_id"]

H3K4genes_in_all <- H3K4me3_gene_peak_counts[H3K4genes_in_40of60_number == 60, "gene_id"]

# total of 15231 genes with H3K4me3 peaks in most samples
# total of 10555 genes with H3K4me3 peaks in most samples

H3K4me3_promoter_peak_counts <- region.peak.counts(NCI60_list = H3K4me3_ranges_list,
                                                   region = "promoters")

H3K4promoters_in_40of60_number <- apply(H3K4me3_promoter_peak_counts[, 5:ncol(H3K4me3_promoter_peak_counts)], 1, function(x){
  
  sum(x > 0)
  
})

H3K4promoters_in_40of60_selec <- H3K4promoters_in_40of60_number > 39
H3K4promoters_in_40of60 <- H3K4me3_promoter_peak_counts[H3K4promoters_in_40of60_selec, "gene_id"]

H3K4promoters_in_all <- H3K4me3_promoter_peak_counts[H3K4promoters_in_40of60_number == 60, "gene_id"]

# total of 12872 promoters with H3K4me3 peaks in most samples
# total of 8986 promoters with H3K4me3 peaks in every sample

#### MAKE H3K4me3 MATRIX FOR LM (GENEBODIES) ####

computematrix_files <- list.files("input/computematrix_out/")
H3K4me3_computematrix_files <- computematrix_files[str_detect(computematrix_files, "H3K4me3")] 
H3K4me3_genes_computematrix_files <- H3K4me3_computematrix_files[str_detect(H3K4me3_computematrix_files, "genebodies")]

uniqueexperiments <- unique(str_remove(str_remove(H3K4me3_genes_computematrix_files, "_bamcompare.bigWig_genebodies200bpflanks_computematrix.gz"), "_rep2"))

H3K4me3_genes_experimentavg_list <- lapply(uniqueexperiments, function(x){
  
  replicatexperiments <- H3K4me3_genes_computematrix_files[str_detect(H3K4me3_genes_computematrix_files, paste0(x, "_"))]
  
  if(length(replicatexperiments) > 2){message(paste0("Warning! More than 2 files detected for this experiment (", x, ")"))}
  
  rep1 <- read.table(paste0("input/computematrix_out/", replicatexperiments[1]), skip = 1)
  rep2 <- read.table(paste0("input/computematrix_out/", replicatexperiments[2]), skip = 1)
  
  combined <- rep1
  combined$V7 <- rowMeans(cbind(rep1$V7, rep2$V7))
  
  return(combined)
  
})

names(H3K4me3_genes_experimentavg_list) <- uniqueexperiments

H3K4genesreadmatrix_df <- sapply(H3K4me3_genes_experimentavg_list, function(x){
  
  x$V7
  
})

row.names(H3K4genesreadmatrix_df) <- H3K4me3_genes_experimentavg_list[[1]]$V4

H3K4genesmatrix_for_lm <- data.frame(t(H3K4genesreadmatrix_df))

H3K4genesmatrix_for_lm <- H3K4genesmatrix_for_lm[, str_remove(colnames(H3K4genesmatrix_for_lm), "^X") %in% H3K4genes_in_all]

H3K4genesmatrix_for_lm[ , "tissue"] <- str_extract(str_remove(row.names(H3K4genesmatrix_for_lm), "H3K4me3_"), "[A-Z]{2,3}")

NNMT_proteomics <- read.table("input/NCI60_NNMT_proteomics.txt", sep = "\t", header = TRUE)
NNMT_mRNA <- read.table("input/NCI60_NNMT_mRNA_zscores.txt", sep = "\t", header = TRUE)

# here I am changing the names from the NNMT proteomics and gene expression file to match the format
# replace slashes with underscores
# replace spaces with underscores
# remove "(TB)" from HL-60
# remove "_ATCC" from A549
names_changed_to_match_format_protein <- str_remove(str_remove(str_replace(str_replace(str_remove(NNMT_proteomics$Cell.Line, "^[A-Z]{2,3}:"), "/", "_"), " ", "_"), "\\(TB\\)"), "_ATCC")
names_changed_to_match_format_RNA <- str_remove(str_remove(str_replace(str_replace(str_remove(NNMT_mRNA$Cell.Line, "^[A-Z]{2,3}:"), "/", "_"), " ", "_"), "\\(TB\\)"), "_ATCC")

# why 59?
# only 59 qvailable data points for proteomics
# 60 for expression so will download expression data separately from CellMinerDB

# the 60th in the RNA dataset is "MDA-N". Will assume it is the last one in the panel which is MDA-MB-468
names_changed_to_match_format_RNA[!sapply(names_changed_to_match_format_RNA, function(x){
  any(str_detect(row.names(H3K4genesmatrix_for_lm), x))
})]

names_changed_to_match_format_RNA <- str_replace(names_changed_to_match_format_RNA, pattern = "MDA-N", replacement = "MDA-MB-468")

NNMT_proteomics$Cell.Line <- names_changed_to_match_format_protein
NNMT_mRNA$Cell.Line <- names_changed_to_match_format_RNA

H3K4genesmatrix_for_lm[, "NNMTprotein"] <- NNMT_proteomics[match(str_remove(str_remove(row.names(H3K4genesmatrix_for_lm), pattern = "H3K4me3_"), pattern = "^[A-Z]{2,3}_"), NNMT_proteomics$Cell.Line), "swaNNMT_nci60"]

H3K4genesmatrix_for_lm[, "NNMTmRNA_z"] <- NNMT_mRNA[match(str_remove(str_remove(row.names(H3K4genesmatrix_for_lm), pattern = "H3K4me3_"), pattern = "^[A-Z]{2,3}_"), NNMT_mRNA$Cell.Line), "expNNMT_nci60"]

H3K4genesmatrix_for_lm[, "NNMT_RNAseq"] <- NCI60_normalised_counts[NNMT_entrezgene_id, match(str_remove(str_remove(row.names(H3K4genesmatrix_for_lm), pattern = "H3K4me3_"), pattern = "^[A-Z]{2,3}_"), colnames(NCI60_normalised_counts))]
H3K4genesmatrix_for_lm[, "PEMT_RNAseq"] <- NCI60_normalised_counts[PEMT_entrezgene_id, match(str_remove(str_remove(row.names(H3K4genesmatrix_for_lm), pattern = "H3K4me3_"), pattern = "^[A-Z]{2,3}_"), colnames(NCI60_normalised_counts))]

H3K4genesnumber_of_sites <- ncol(H3K4genesmatrix_for_lm) - 5

saveRDS(H3K4genesmatrix_for_lm, "output/NCI60_H3K4me3_genes_signal_matrix_for_lm.rds")
# H3K4genesmatrix_for_lm <- readRDS("output/NCI60_H3K4me3_genes_signal_matrix_for_lm.rds")

#### MAKE H3K4me3 MATRIX FOR LM (PROMOTERS) ####

# point to deeptools computematrix output directory
computematrix_files <- list.files("input/computematrix_out/")
H3K4me3_computematrix_files <- computematrix_files[str_detect(computematrix_files, "H3K4me3")] 
H3K4me3_promoters_computematrix_files <- H3K4me3_computematrix_files[str_detect(H3K4me3_computematrix_files, "promoter")]

uniqueexperiments <- unique(str_remove(str_remove(H3K4me3_promoters_computematrix_files, "_bamcompare.bigWig_promoters_computematrix.gz"), "_rep2"))

H3K4me3_promoters_experimentavg_list <- lapply(uniqueexperiments, function(x){
  
  replicatexperiments <- H3K4me3_promoters_computematrix_files[str_detect(H3K4me3_promoters_computematrix_files, paste0(x, "_"))]
  
  if(length(replicatexperiments) > 2){message(paste0("Warning! More than 2 files detected for this experiment (", x, ")"))}
  
  rep1 <- read.table(paste0("input/computematrix_out/", replicatexperiments[1]), skip = 1)
  rep2 <- read.table(paste0("input/computematrix_out/", replicatexperiments[2]), skip = 1)
  
  combined <- rep1
  combined$V7 <- rowMeans(cbind(rep1$V7, rep2$V7))
  
  return(combined)
  
})

names(H3K4me3_promoters_experimentavg_list) <- uniqueexperiments

H3K4promotersreadmatrix_df <- sapply(H3K4me3_promoters_experimentavg_list, function(x){
  
  x$V7
  
})

row.names(H3K4promotersreadmatrix_df) <- H3K4me3_promoters_experimentavg_list[[1]]$V4

H3K4promotersmatrix_for_lm <- data.frame(t(H3K4promotersreadmatrix_df))

H3K4promotersmatrix_for_lm <- H3K4promotersmatrix_for_lm[, str_remove(colnames(H3K4promotersmatrix_for_lm), "^X") %in% H3K4promoters_in_all]

H3K4promotersmatrix_for_lm[ , "tissue"] <- str_extract(str_remove(row.names(H3K4promotersmatrix_for_lm), "H3K4me3_"), "[A-Z]{2,3}")

H3K4promotersmatrix_for_lm[, "NNMTprotein"] <- NNMT_proteomics[match(str_remove(str_remove(row.names(H3K4promotersmatrix_for_lm), pattern = "H3K4me3_"), pattern = "^[A-Z]{2,3}_"), NNMT_proteomics$Cell.Line), "swaNNMT_nci60"]

H3K4promotersmatrix_for_lm[, "NNMTmRNA_z"] <- NNMT_mRNA[match(str_remove(str_remove(row.names(H3K4promotersmatrix_for_lm), pattern = "H3K4me3_"), pattern = "^[A-Z]{2,3}_"), NNMT_mRNA$Cell.Line), "expNNMT_nci60"]

H3K4promotersmatrix_for_lm[, "NNMT_RNAseq"] <- NCI60_normalised_counts[NNMT_entrezgene_id, match(str_remove(str_remove(row.names(H3K4promotersmatrix_for_lm), pattern = "H3K4me3_"), pattern = "^[A-Z]{2,3}_"), colnames(NCI60_normalised_counts))]
H3K4promotersmatrix_for_lm[, "PEMT_RNAseq"] <- NCI60_normalised_counts[PEMT_entrezgene_id, match(str_remove(str_remove(row.names(H3K4promotersmatrix_for_lm), pattern = "H3K4me3_"), pattern = "^[A-Z]{2,3}_"), colnames(NCI60_normalised_counts))]

H3K4promotersnumber_of_sites <- ncol(H3K4promotersmatrix_for_lm) - 5

saveRDS(H3K4promotersmatrix_for_lm, "output/NCI60_H3K4me3_promoters_signal_matrix_for_lm.rds")
# H3K4promotersmatrix_for_lm <- readRDS("output/NCI60_H3K4me3_promoters_signal_matrix_for_lm.rds")

#### FIT H3K4 NNMT mRNA PROMOTERS MODEL ####

fit.model.to.NCI60.chip <- function(histone_mark,
                                    siteofchoice = c("promoters", "genes", "repeats"),
                                    gene_enquiry = c("PEMT", "NNMT"),
                                    expressiontype = c("mRNA_z", "protein", "_RNAseq")
){
  
  datamatrix <- get(paste0(histone_mark, siteofchoice, "matrix_for_lm"))
  datasitenumber <- get(paste0(histone_mark, siteofchoice, "number_of_sites"))
  modelterm <- paste0(gene_enquiry, expressiontype)
  
  output <- lapply(colnames(datamatrix)[1:datasitenumber], function(x){

    tryCatch(expr = {

      templm <- glm(formula = as.formula(paste0(x, " ~ ", modelterm, " + tissue + ", modelterm, " * tissue")), data = datamatrix)
      templmsum <- summary(templm)
      
      t <- templmsum$coefficients[paste0(modelterm), "t value"]
      
      output_vec <- c(t)
      names(output_vec) <- c("t_value")
      
      return(output_vec)
      
    }, error = function(e){
      
      output_vec <- c(NA)
      names(output_vec) <- c("t_value")
      
      return(output_vec)
      
    }
    )
    
  })
  
  names(output) <- colnames(datamatrix)[1:datasitenumber]
  
  return(output)
  
}

H3K4NNMTmRNA_promoters_signal_model <- fit.model.to.NCI60.chip(histone_mark = "H3K4",
                                                               siteofchoice = "promoters",
                                                               gene_enquiry = "NNMT",
                                                               expressiontype = "mRNA_z")

saveRDS(H3K4NNMTmRNA_promoters_signal_model, "output/NCI60_H3K4me3_promoters_signal_model_NNMTmRNA.rds")

H3K4NNMTmRNA_genes_signal_model <- fit.model.to.NCI60.chip(histone_mark = "H3K4",
                                                           siteofchoice = "genes",
                                                           gene_enquiry = "NNMT",
                                                           expressiontype = "mRNA_z")

saveRDS(H3K4NNMTmRNA_genes_signal_model, "output/NCI60_H3K4me3_genes_signal_model_NNMTmRNA.rds")

H3K4NNMTmRNA_promoters_signal_model <- fit.model.to.NCI60.chip(histone_mark = "H3K4",
                                                               siteofchoice = "promoters",
                                                               gene_enquiry = "NNMT",
                                                               expressiontype = "protein")

saveRDS(H3K4NNMTmRNA_promoters_signal_model, "output/NCI60_H3K4me3_promoters_signal_model_NNMTprotein.rds")

H3K4NNMTprotein_genes_signal_model <- fit.model.to.NCI60.chip(histone_mark = "H3K4",
                                                              siteofchoice = "genes",
                                                              gene_enquiry = "NNMT",
                                                              expressiontype = "protein")

saveRDS(H3K4NNMTprotein_genes_signal_model, "output/NCI60_H3K4me3_genes_signal_model_NNMTprotein.rds")

H3K4NNMT_RNAseq_promoters_signal_model <- fit.model.to.NCI60.chip(histone_mark = "H3K4",
                                                                  siteofchoice = "promoters",
                                                                  gene_enquiry = "NNMT",
                                                                  expressiontype = "_RNAseq")

saveRDS(H3K4NNMT_RNAseq_promoters_signal_model, "output/NCI60_H3K4me3_promoters_signal_model_NNMT_RNAseq.rds")

H3K4NNMT_RNAseq_genes_signal_model <- fit.model.to.NCI60.chip(histone_mark = "H3K4",
                                                              siteofchoice = "genes",
                                                              gene_enquiry = "NNMT",
                                                              expressiontype = "_RNAseq")

saveRDS(H3K4NNMT_RNAseq_genes_signal_model, "output/NCI60_H3K4me3_genes_signal_model_NNMT_RNAseq.rds")

H3K4PEMT_RNAseq_promoters_signal_model <- fit.model.to.NCI60.chip(histone_mark = "H3K4",
                                                                  siteofchoice = "promoters",
                                                                  gene_enquiry = "PEMT",
                                                                  expressiontype = "_RNAseq")

saveRDS(H3K4PEMT_RNAseq_promoters_signal_model, "output/NCI60_H3K4me3_promoters_signal_model_PEMT_RNAseq.rds")

H3K4PEMT_RNAseq_genes_signal_model <- fit.model.to.NCI60.chip(histone_mark = "H3K4",
                                                              siteofchoice = "genes",
                                                              gene_enquiry = "PEMT",
                                                              expressiontype = "_RNAseq")

saveRDS(H3K4PEMT_RNAseq_genes_signal_model, "output/NCI60_H3K4me3_genes_signal_model_PEMT_RNAseq.rds")

#### H3K4 random distributions ####

#### RANDOM DISTRIBUTIONS H3K4me3 RNAseq GENE BODIES ####

nuclear_genes_expressed <- readRDS("output/GTEX_nuclear_genes_expressed.rds")
nuclear_genes_expressed_entrez <- getBM(attributes = c("ensembl_gene_id", "entrezgene_id"),
                                        values = nuclear_genes_expressed,
                                        filter = "ensembl_gene_id",
                                        mart = ensembl)

saveRDS(nuclear_genes_expressed_entrez, "output/nuclear_genes_expressed_entrez.rds")
# nuclear_genes_expressed_entrez <- readRDS("output/nuclear_genes_expressed_entrez.rds")

RNAseq_data_expressed <- NCI60_normalised_counts[row.names(NCI60_normalised_counts) %in% nuclear_genes_expressed_entrez$entrezgene_id, ]
row.names(RNAseq_data_expressed) <- nuclear_genes_expressed_entrez[match(row.names(RNAseq_data_expressed), nuclear_genes_expressed_entrez$entrezgene_id), "ensembl_gene_id"]

# will use same random genes for mRNA and for RNA-seq
random_expressed_mRNA_sample <- sample(colnames(mRNA_expressed_datat), 1000, replace = FALSE)
saveRDS(random_expressed_mRNA_sample, "output/random_expressed_mRNA_sample.rds")
# random_expressed_mRNA_sample <- readRDS("output/random_expressed_mRNA_sample.rds")

for(i in 1:length(random_expressed_mRNA_sample)){
  
  H3K4genesmatrix_for_lm[, paste0(random_expressed_mRNA_sample[i], "_RNAseq")] <- as.numeric(RNAseq_data_expressed[random_expressed_mRNA_sample[i], match(str_remove(row.names(H3K4genesmatrix_for_lm), "H3K4me3_[A-Z]{2,3}_"), colnames(RNAseq_data_expressed))])
  
}

RNAseqrandomsample_H3K4genes_t_values <- sapply(random_expressed_mRNA_sample, function(thisgene){
  
  sapply(colnames(H3K4genesmatrix_for_lm)[1:H3K4genesnumber_of_sites], function(x){
    
    tryCatch(expr = {
      glmtemp <- glm(formula = as.formula(paste0(x, " ~ log10(", as.character(thisgene), "_RNAseq + 1) + tissue + log10(", thisgene, "_RNAseq + 1) * tissue")), data = H3K4genesmatrix_for_lm)
      glmsum <- summary(glmtemp)
      
      t <- glmsum$coefficients[paste0("log10(", thisgene, "_RNAseq + 1)"), "t value"]
      
      output_vec <- c(t)
      names(output_vec) <- c("t_value")
      
      return(output_vec)
      
    }, error = function(e){
      
      output_vec <- c(NA)
      names(output_vec) <- c("t_value")
      
      return(output_vec)
      
    }
    )
    
  })
  
})

saveRDS(RNAseqrandomsample_H3K4genes_t_values, "output/RNAseqrandomsample_H3K4me3genes_t_values.rds")

#### RANDOM DISTRIBUTIONS H3K4me3 RNAseq PROMOTERS ####

RNAseq_data_expressed <- NCI60_normalised_counts[row.names(NCI60_normalised_counts) %in% nuclear_promoters_expressed_entrez$entrezgene_id, ]
row.names(RNAseq_data_expressed) <- nuclear_promoters_expressed_entrez[match(row.names(RNAseq_data_expressed), nuclear_promoters_expressed_entrez$entrezgene_id), "ensembl_gene_id"]

# random_expressed_mRNA_sample <- sample(colnames(mRNA_expressed_datat), 1000, replace = FALSE)
# saveRDS(random_expressed_mRNA_sample, "output/random_expressed_mRNA_sample.rds")
random_expressed_mRNA_sample <- readRDS("output/random_expressed_mRNA_sample.rds")

for(i in 1:length(random_expressed_mRNA_sample)){
  
  H3K4promotersmatrix_for_lm[, paste0(random_expressed_mRNA_sample[i], "_RNAseq")] <- as.numeric(RNAseq_data_expressed[random_expressed_mRNA_sample[i], match(str_remove(row.names(H3K4promotersmatrix_for_lm), "H3K4me3_[A-Z]{2,3}_"), colnames(RNAseq_data_expressed))])
  
}

RNAseqrandomsample_H3K4promoters_t_values <- sapply(random_expressed_mRNA_sample, function(thisgene){
  
  sapply(colnames(H3K4promotersmatrix_for_lm)[1:H3K4promotersnumber_of_sites], function(x){
    
    tryCatch(expr = {
      glmtemp <- glm(formula = as.formula(paste0(x, " ~ log10(", as.character(thisgene), "_RNAseq + 1) + tissue + log10(", thisgene, "_RNAseq + 1) * tissue")), data = H3K4promotersmatrix_for_lm)
      glmsum <- summary(glmtemp)
      
      t <- glmsum$coefficients[paste0("log10(", thisgene, "_RNAseq + 1)"), "t value"]
      
      output_vec <- c(t)
      names(output_vec) <- c("t_value")
      
      return(output_vec)
      
    }, error = function(e){
      
      output_vec <- c(NA)
      names(output_vec) <- c("t_value")
      
      return(output_vec)
      
    }
    )
    
  })
  
})

saveRDS(RNAseqrandomsample_H3K4promoters_t_values, "output/RNAseqrandomsample_H3K4me3promoters_t_values.rds")

#### RANDOM DISTRIBUTION OF H3K4 PROMOTERS MRNA ####

# downloaded from [https://discover.nci.nih.gov/cellminer/loadDownload.do]
all_mRNA_data <- read_excel("input/RNA__5_Platform_Gene_Transcript_Average_z_scores.xls", skip = 10)
names_changed_to_match_format_RNAsheet <- str_remove(str_remove(str_replace(str_replace(str_remove(colnames(all_mRNA_data)[7:ncol(all_mRNA_data)], "^[A-Z]{2,3}:"), "/", "_"), " ", "_"), "\\(TB\\)"), "_ATCC")
names_changed_to_match_format_RNAsheet <- str_replace(names_changed_to_match_format_RNAsheet, pattern = "MDA-N", replacement = "MDA-MB-468")
colnames(all_mRNA_data)[7:ncol(all_mRNA_data)] <- names_changed_to_match_format_RNAsheet

mRNAdatat <- t(all_mRNA_data)
colnames(mRNAdatat) <- mRNAdatat[1,]

all_mRNA_data_expressed <- all_mRNA_data[all_mRNA_data$`Entrez gene id e` %in% nuclear_genes_expressed_entrez$entrezgene_id, ]
all_mRNA_data_expressed[, "ensembl_gene_id"] <- nuclear_genes_expressed_entrez[match(all_mRNA_data_expressed$`Entrez gene id e`, nuclear_genes_expressed_entrez$entrezgene_id), "ensembl_gene_id"]

mRNA_expressed_datat <- t(all_mRNA_data_expressed)
colnames(mRNA_expressed_datat) <- mRNA_expressed_datat["ensembl_gene_id",]

# random_expressed_mRNA_sample <- sample(colnames(mRNA_expressed_datat), 1000, replace = FALSE)
# saveRDS(random_expressed_mRNA_sample, "output/random_expressed_mRNA_sample.rds")
random_expressed_mRNA_sample <- readRDS("output/random_expressed_mRNA_sample.rds")

for(i in 1:length(random_expressed_mRNA_sample)){
  
  H3K4promotersmatrix_for_lm[, paste0(random_expressed_mRNA_sample[i], "_mRNA")] <- as.numeric(mRNA_expressed_datat[match(str_remove(row.names(H3K4promotersmatrix_for_lm), "H3K4me3_[A-Z]{2,3}_"), row.names(mRNA_expressed_datat)), random_expressed_mRNA_sample[i]])
  
}

mRNArandomsample_H3K4me3promoters_t_values <- sapply(random_expressed_mRNA_sample, function(thisgene){
  
  sapply(colnames(H3K4promotersmatrix_for_lm)[1:H3K4promotersnumber_of_sites], function(x){
    
    tryCatch(expr = {
      glmtemp <- glm(formula = as.formula(paste0(x, " ~ ", as.character(thisgene), "_mRNA + tissue + ", thisgene, "_mRNA * tissue")), data = H3K4promotersmatrix_for_lm)
      glmsum <- summary(glmtemp)
      
      t <- glmsum$coefficients[paste0(thisgene, "_mRNA"), "t value"]
      
      output_vec <- c(t)
      names(output_vec) <- c("t_value")
      
      return(output_vec)
      
    }, error = function(e){
      
      output_vec <- c(NA)
      names(output_vec) <- c("t_value")
      
      return(output_vec)
      
    }
    )
    
  })
  
})

saveRDS(mRNArandomsample_H3K4me3promoters_t_values, "output/mRNArandomsample_H3K4me3promoters_t_values.rds")
# mRNArandomsample_H3K4me3promoters_t_values <- readRDS("output/mRNArandomsample_H3K4me3promoters_t_values.rds")

#### RANDOM DISTRIBUTION OF H3K4 GENES MRNA ####

# random_expressed_mRNA_sample <- sample(colnames(mRNA_expressed_datat), 1000, replace = FALSE)
# saveRDS(random_expressed_mRNA_sample, "output/random_expressed_mRNA_sample.rds")
random_expressed_mRNA_sample <- readRDS("output/random_expressed_mRNA_sample.rds")

for(i in 1:length(random_expressed_mRNA_sample)){
  
  H3K4genesmatrix_for_lm[, paste0(random_expressed_mRNA_sample[i], "_mRNA")] <- as.numeric(mRNA_expressed_datat[match(str_remove(row.names(H3K4genesmatrix_for_lm), "H3K4me3_[A-Z]{2,3}_"), row.names(mRNA_expressed_datat)), random_expressed_mRNA_sample[i]])
  
}

mRNArandomsample_H3K4me3genes_t_values <- sapply(random_expressed_mRNA_sample, function(thisgene){
  
  sapply(colnames(H3K4genesmatrix_for_lm)[1:H3K4genesnumber_of_sites], function(x){
    
    tryCatch(expr = {
      glmtemp <- glm(formula = as.formula(paste0(x, " ~ ", as.character(thisgene), "_mRNA + tissue + ", thisgene, "_mRNA * tissue")), data = H3K4genesmatrix_for_lm)
      glmsum <- summary(glmtemp)
      
      t <- glmsum$coefficients[paste0(thisgene, "_mRNA"), "t value"]
      
      output_vec <- c(t)
      names(output_vec) <- c("t_value")
      
      return(output_vec)
      
    }, error = function(e){
      
      output_vec <- c(NA)
      names(output_vec) <- c("t_value")
      
      return(output_vec)
      
    }
    )
    
  })
  
})

saveRDS(mRNArandomsample_H3K4me3genes_t_values, "output/mRNArandomsample_H3K4me3genes_t_values.rds")
mRNArandomsample_H3K4me3genes_t_values <- readRDS("output/mRNArandomsample_H3K4me3genes_t_values.rds")

#### RANDOM DISTRIBUTION OF H3K4 PROMOTERS PROTEIN ####

# downloaded from [https://discover.nci.nih.gov/cellminer/loadDownload.do]
all_protein_data <- read_excel("input/Protein__SWATH_(Mass_spectrometry)_Protein.xls", skip = 10)

names_changed_to_match_format_protsheet <- str_remove(str_remove(str_replace(str_replace(str_remove(colnames(all_protein_data)[10:ncol(all_protein_data)], "^[A-Z]{2,3}:"), "/", "_"), " ", "_"), "\\(TB\\)"), "_ATCC")
colnames(all_protein_data)[10:ncol(all_protein_data)] <- names_changed_to_match_format_protsheet

protdatat <- t(all_protein_data)
colnames(protdatat) <- protdatat[2, ]

random_detected_protein_sample <- sample(all_protein_data$`Gene name d`, 1000, replace = FALSE)
saveRDS(random_detected_protein_sample, "output/random_detected_protein_sample.rds" )
# random_detected_protein_sample <- readRDS("output/random_detected_protein_sample.rds" )

for(i in 1:length(random_detected_protein_sample)){
  
  H3K4promotersmatrix_for_lm[, paste0(random_detected_protein_sample[i], "_protein")] <- as.numeric(protdatat[match(str_remove(row.names(H3K4promotersmatrix_for_lm), "H3K4me3_[A-Z]{2,3}_"), row.names(protdatat)), random_detected_protein_sample[i]])
  
}

proteinrandomsample_H3K4me3promoters_t_values <- sapply(random_detected_protein_sample, function(thisgene){
  
  sapply(colnames(H3K4promotersmatrix_for_lm)[1:H3K4promotersnumber_of_sites], function(x){
    
    tryCatch(expr = {
      glmtemp <- glm(formula = as.formula(paste0(x, " ~ ", as.character(thisgene), "_protein + tissue + ", thisgene, "_protein * tissue")), data = H3K4promotersmatrix_for_lm)
      glmsum <- summary(glmtemp)
      
      t <- glmsum$coefficients[paste0(thisgene, "_protein"), "t value"]
      
      output_vec <- c(t)
      names(output_vec) <- c("t_value")
      
      return(output_vec)
      
    }, error = function(e){
      
      output_vec <- c(NA)
      names(output_vec) <- c("t_value")
      
      return(output_vec)
      
    }
    )
    
  })
  
})

saveRDS(proteinrandomsample_H3K4me3promoters_t_values, "output/proteinrandomsample_H3K4me3promoters_t_values.rds")
# proteinrandomsample_H3K4me3promoters_t_values <- readRDS("output/proteinrandomsample_H3K4me3promoters_t_values.rds")

#### RANDOM DISTRIBUTION OF H3K4 GENE BODIES PROTEIN ####

# random_detected_protein_sample <- sample(all_protein_data$`Gene name d`, 1000, replace = FALSE)
# saveRDS(random_detected_protein_sample, "output/random_detected_protein_sample.rds" )
random_detected_protein_sample <- readRDS("output/random_detected_protein_sample.rds" )

for(i in 1:length(random_detected_protein_sample)){
  
  H3K4genesmatrix_for_lm[, paste0(random_detected_protein_sample[i], "_protein")] <- as.numeric(protdatat[match(str_remove(row.names(H3K4genesmatrix_for_lm), "H3K4me3_[A-Z]{2,3}_"), row.names(protdatat)), random_detected_protein_sample[i]])
  
}

proteinrandomsample_H3K4me3genes_t_values <- sapply(random_detected_protein_sample, function(thisgene){
  
  sapply(colnames(H3K4genesmatrix_for_lm)[1:H3K4genesnumber_of_sites], function(x){
    
    tryCatch(expr = {
      glmtemp <- glm(formula = as.formula(paste0(x, " ~ ", as.character(thisgene), "_protein + tissue + ", thisgene, "_protein * tissue")), data = H3K4genesmatrix_for_lm)
      glmsum <- summary(glmtemp)
      
      t <- glmsum$coefficients[paste0(thisgene, "_protein"), "t value"]
      
      output_vec <- c(t)
      names(output_vec) <- c("t_value")
      
      return(output_vec)
      
    }, error = function(e){
      
      output_vec <- c(NA)
      names(output_vec) <- c("t_value")
      
      return(output_vec)
      
    }
    )
    
  })
  
})

saveRDS(proteinrandomsample_H3K4me3genes_t_values, "output/proteinrandomsample_H3K4me3genes_t_values.rds")
# proteinrandomsample_H3K4me3genes_t_values <- readRDS("output/proteinrandomsample_H3K4me3genes_t_values.rds")

#### H3K9me3 ####

#### combine narrow and broad peak options####

# point to directory with genrich output
my_directory <- "input/NCI60_genrich_peak_calls/"

genrich_file_list <- list.files(my_directory)

# choose the narrowPeak files
genrich_file_list_narrowPeak <- genrich_file_list[str_detect(genrich_file_list, pattern = "\\.narrowPeak")]

# choose the H3K9me3 files
H3K9me3_file_list_narrowPeak <- genrich_file_list_narrowPeak[str_detect(genrich_file_list_narrowPeak, pattern = "H3K9me3")]

# loop through unique experiments to pool the two replicates

H3K9me3_ranges_list <- lapply(H3K9me3_file_list_narrowPeak, function(x){
  
  peaks_1 <- read.table(paste0(my_directory, x), header = FALSE)
  colnames(peaks_1)[1:3] <- c("chr", "start", "end")
  
  peaks_1_GR <- makeGRangesFromDataFrame(peaks_1)
  
  
})

names(H3K9me3_ranges_list) <- unique(str_remove(H3K9me3_file_list_narrowPeak, pattern = "_genrich.*"))

saveRDS(H3K9me3_ranges_list, "output/NCI60_H3K9me3_combinedpeaks_list.rds")
# H3K9me3_ranges_list <- readRDS("output/NCI60_H3K9me3_combinedpeaks_list.rds")

#### SELECT H3K9me3 REGIONS WITH CONSISTENT PEAKS ####

H3K9me3_gene_peak_counts <- region.peak.counts(NCI60_list = H3K9me3_ranges_list,
                                               region = "genes")

H3K9genes_in_40of60_number <- apply(H3K9me3_gene_peak_counts[, 5:ncol(H3K9me3_gene_peak_counts)], 1, function(x){
  
  sum(x > 0)
  
})

H3K9genes_in_40of60_selec <- H3K9genes_in_40of60_number > 39

H3K9genes_in_40of60 <- H3K9me3_gene_peak_counts[H3K9genes_in_40of60_selec, "gene_id"]

H3K9genes_in_all <- H3K9me3_gene_peak_counts[H3K9genes_in_40of60_number == 60, "gene_id"]

# total of 307 genes with H3K9me3 peaks in most samples
# total of 11 genes with H3K9me3 peaks in most samples

H3K9me3_promoter_peak_counts <- region.peak.counts(NCI60_list = H3K9me3_ranges_list,
                                                   region = "promoters")

H3K9promoters_in_40of60_number <- apply(H3K9me3_promoter_peak_counts[, 5:ncol(H3K9me3_promoter_peak_counts)], 1, function(x){
  
  sum(x > 0)
  
})

H3K9promoters_in_40of60_selec <- H3K9promoters_in_40of60_number > 39
H3K9promoters_in_40of60 <- H3K9me3_promoter_peak_counts[H3K9promoters_in_40of60_selec, "gene_id"]

H3K9promoters_in_all <- H3K9me3_promoter_peak_counts[H3K9promoters_in_40of60_number == 60, "gene_id"]

# total of 25 promoters with H3K9me3 peaks in most samples
# total of 1 promoters with H3K9me3 peaks in every sample

#### MAKE H3K9me3 MATRIX FOR LM (GENEBODIES) ####

computematrix_files <- list.files("input/NCI60_H3K9_computematrix/")
H3K9me3_computematrix_files <- computematrix_files[str_detect(computematrix_files, "H3K9me3")] 
H3K9me3_genes_computematrix_files <- H3K9me3_computematrix_files[str_detect(H3K9me3_computematrix_files, "genebodies")]

uniqueexperiments <- unique(str_remove(str_remove(H3K9me3_genes_computematrix_files, "_bamcompare.bigWig_genebodies200bpflanks_computematrix.gz"), "_rep2"))

H3K9me3_genes_experimentavg_list <- lapply(uniqueexperiments, function(x){
  
  replicatexperiments <- H3K9me3_genes_computematrix_files[str_detect(H3K9me3_genes_computematrix_files, paste0(x, "_"))]
  
  if(length(replicatexperiments) > 2){message(paste0("Warning! More than 2 files detected for this experiment (", x, ")"))}
  
  rep1 <- read.table(paste0("input/NCI60_H3K9_computematrix/", replicatexperiments[1]), skip = 1)
  rep2 <- read.table(paste0("input/NCI60_H3K9_computematrix/", replicatexperiments[2]), skip = 1)
  
  combined <- rep1
  combined$V7 <- rowMeans(cbind(rep1$V7, rep2$V7))
  
  return(combined)
  
})

names(H3K9me3_genes_experimentavg_list) <- uniqueexperiments

H3K9genesreadmatrix_df <- sapply(H3K9me3_genes_experimentavg_list, function(x){
  
  x$V7
  
})

row.names(H3K9genesreadmatrix_df) <- H3K9me3_genes_experimentavg_list[[1]]$V4

H3K9genesmatrix_for_lm <- data.frame(t(H3K9genesreadmatrix_df))

H3K9genesmatrix_for_lm <- H3K9genesmatrix_for_lm[, str_remove(colnames(H3K9genesmatrix_for_lm), "^X") %in% H3K9genes_in_40of60]

H3K9genesmatrix_for_lm[ , "tissue"] <- str_extract(str_remove(row.names(H3K9genesmatrix_for_lm), "H3K9me3_"), "[A-Z]{2,3}")

NNMT_proteomics <- read.table("input/NCI60_NNMT_proteomics.txt", sep = "\t", header = TRUE)
NNMT_mRNA <- read.table("input/NCI60_NNMT_mRNA_zscores.txt", sep = "\t", header = TRUE)

# here I am changing the names from the NNMT proteomics and gene expression file to match the format
# replace slashes with underscores
# replace spaces with underscores
# remove "(TB)" from HL-60
# remove "_ATCC" from A549
names_changed_to_match_format_protein <- str_remove(str_remove(str_replace(str_replace(str_remove(NNMT_proteomics$Cell.Line, "^[A-Z]{2,3}:"), "/", "_"), " ", "_"), "\\(TB\\)"), "_ATCC")
names_changed_to_match_format_RNA <- str_remove(str_remove(str_replace(str_replace(str_remove(NNMT_mRNA$Cell.Line, "^[A-Z]{2,3}:"), "/", "_"), " ", "_"), "\\(TB\\)"), "_ATCC")

# why 59?
# only 59 qvailable data points for proteomics
# 60 for expression so will download expression data separately from CellMinerDB

# the 60th in the RNA dataset is "MDA-N". Will assume it is the last one in the panel which is MDA-MB-468
names_changed_to_match_format_RNA[!sapply(names_changed_to_match_format_RNA, function(x){
  any(str_detect(row.names(H3K9genesmatrix_for_lm), x))
})]

names_changed_to_match_format_RNA <- str_replace(names_changed_to_match_format_RNA, pattern = "MDA-N", replacement = "MDA-MB-468")

NNMT_proteomics$Cell.Line <- names_changed_to_match_format_protein
NNMT_mRNA$Cell.Line <- names_changed_to_match_format_RNA

H3K9genesmatrix_for_lm[, "NNMTprotein"] <- NNMT_proteomics[match(str_remove(str_remove(row.names(H3K9genesmatrix_for_lm), pattern = "H3K9me3_"), pattern = "^[A-Z]{2,3}_"), NNMT_proteomics$Cell.Line), "swaNNMT_nci60"]

H3K9genesmatrix_for_lm[, "NNMTmRNA_z"] <- NNMT_mRNA[match(str_remove(str_remove(row.names(H3K9genesmatrix_for_lm), pattern = "H3K9me3_"), pattern = "^[A-Z]{2,3}_"), NNMT_mRNA$Cell.Line), "expNNMT_nci60"]

H3K9genesmatrix_for_lm[, "NNMT_RNAseq"] <- NCI60_normalised_counts[NNMT_entrezgene_id, match(str_remove(str_remove(row.names(H3K9genesmatrix_for_lm), pattern = "H3K9me3_"), pattern = "^[A-Z]{2,3}_"), colnames(NCI60_normalised_counts))]
H3K9genesmatrix_for_lm[, "PEMT_RNAseq"] <- NCI60_normalised_counts[PEMT_entrezgene_id, match(str_remove(str_remove(row.names(H3K9genesmatrix_for_lm), pattern = "H3K9me3_"), pattern = "^[A-Z]{2,3}_"), colnames(NCI60_normalised_counts))]

H3K9genesnumber_of_sites <- ncol(H3K9genesmatrix_for_lm) - 5

saveRDS(H3K9genesmatrix_for_lm, "output/NCI60_H3K9me3_genes_signal_matrix_for_lm.rds")
# H3K9genesmatrix_for_lm <- readRDS("output/NCI60_H3K9me3_genes_signal_matrix_for_lm.rds")

#### MAKE H3K9me3 MATRIX FOR LM (PROMOTERS) ####

# point to deeptools computematrix output directory
computematrix_files <- list.files("input/computematrix_out/")
H3K9me3_computematrix_files <- computematrix_files[str_detect(computematrix_files, "H3K9me3")] 
H3K9me3_promoters_computematrix_files <- H3K9me3_computematrix_files[str_detect(H3K9me3_computematrix_files, "promoter")]

uniqueexperiments <- unique(str_remove(str_remove(H3K9me3_promoters_computematrix_files, "_bamcompare.bigWig_promoters_computematrix.gz"), "_rep2"))

H3K9me3_promoters_experimentavg_list <- lapply(uniqueexperiments, function(x){
  
  replicatexperiments <- H3K9me3_promoters_computematrix_files[str_detect(H3K9me3_promoters_computematrix_files, paste0(x, "_"))]
  
  if(length(replicatexperiments) > 2){message(paste0("Warning! More than 2 files detected for this experiment (", x, ")"))}
  
  rep1 <- read.table(paste0("input/computematrix_out/", replicatexperiments[1]), skip = 1)
  rep2 <- read.table(paste0("input/computematrix_out/", replicatexperiments[2]), skip = 1)
  
  combined <- rep1
  combined$V7 <- rowMeans(cbind(rep1$V7, rep2$V7))
  
  return(combined)
  
})

names(H3K9me3_promoters_experimentavg_list) <- uniqueexperiments

H3K9promotersreadmatrix_df <- sapply(H3K9me3_promoters_experimentavg_list, function(x){
  
  x$V7
  
})

row.names(H3K9promotersreadmatrix_df) <- H3K9me3_promoters_experimentavg_list[[1]]$V4

H3K9promotersmatrix_for_lm <- data.frame(t(H3K9promotersreadmatrix_df))

H3K9promotersmatrix_for_lm <- H3K9promotersmatrix_for_lm[, str_remove(colnames(H3K9promotersmatrix_for_lm), "^X") %in% H3K9promoters_in_40of60]

H3K9promotersmatrix_for_lm[ , "tissue"] <- str_extract(str_remove(row.names(H3K9promotersmatrix_for_lm), "H3K9me3_"), "[A-Z]{2,3}")

H3K9promotersmatrix_for_lm[, "NNMTprotein"] <- NNMT_proteomics[match(str_remove(str_remove(row.names(H3K9promotersmatrix_for_lm), pattern = "H3K9me3_"), pattern = "^[A-Z]{2,3}_"), NNMT_proteomics$Cell.Line), "swaNNMT_nci60"]

H3K9promotersmatrix_for_lm[, "NNMTmRNA_z"] <- NNMT_mRNA[match(str_remove(str_remove(row.names(H3K9promotersmatrix_for_lm), pattern = "H3K9me3_"), pattern = "^[A-Z]{2,3}_"), NNMT_mRNA$Cell.Line), "expNNMT_nci60"]

H3K9promotersmatrix_for_lm[, "NNMT_RNAseq"] <- NCI60_normalised_counts[NNMT_entrezgene_id, match(str_remove(str_remove(row.names(H3K9promotersmatrix_for_lm), pattern = "H3K9me3_"), pattern = "^[A-Z]{2,3}_"), colnames(NCI60_normalised_counts))]
H3K9promotersmatrix_for_lm[, "PEMT_RNAseq"] <- NCI60_normalised_counts[PEMT_entrezgene_id, match(str_remove(str_remove(row.names(H3K9promotersmatrix_for_lm), pattern = "H3K9me3_"), pattern = "^[A-Z]{2,3}_"), colnames(NCI60_normalised_counts))]

H3K9promotersnumber_of_sites <- ncol(H3K9promotersmatrix_for_lm) - 5

saveRDS(H3K9promotersmatrix_for_lm, "output/NCI60_H3K9me3_promoters_signal_matrix_for_lm.rds")
# H3K9promotersmatrix_for_lm <- readRDS("output/NCI60_H3K9me3_promoters_signal_matrix_for_lm.rds")

#### FIT H3K9 NNMT MODELS ####

H3K9NNMTmRNA_promoters_signal_model <- fit.model.to.NCI60.chip(histone_mark = "H3K9",
                                                               siteofchoice = "promoters",
                                                               gene_enquiry = "NNMT",
                                                               expressiontype = "mRNA_z")

saveRDS(H3K9NNMTmRNA_promoters_signal_model, "output/NCI60_H3K9me3_promoters_signal_model_NNMTmRNA.rds")

H3K9NNMTmRNA_genes_signal_model <- fit.model.to.NCI60.chip(histone_mark = "H3K9",
                                                           siteofchoice = "genes",
                                                           gene_enquiry = "NNMT",
                                                           expressiontype = "mRNA_z")

saveRDS(H3K9NNMTmRNA_genes_signal_model, "output/NCI60_H3K9me3_genes_signal_model_NNMTmRNA.rds")

H3K9NNMTmRNA_repeats_signal_model <- fit.model.to.NCI60.chip(histone_mark = "H3K9",
                                                             siteofchoice = "repeats",
                                                             gene_enquiry = "NNMT",
                                                             expressiontype = "mRNA_z")

saveRDS(H3K9NNMTmRNA_repeats_signal_model, "output/NCI60_H3K9me3_repeats_signal_model_NNMT_mRNA.rds")

H3K9NNMTprotein_promoters_signal_model <- fit.model.to.NCI60.chip(histone_mark = "H3K9",
                                                                  siteofchoice = "promoters",
                                                                  gene_enquiry = "NNMT",
                                                                  expressiontype = "protein")

saveRDS(H3K9NNMTprotein_promoters_signal_model, "output/NCI60_H3K9me3_promoters_signal_model_NNMTprotein.rds")

H3K9NNMTprotein_genes_signal_model <- fit.model.to.NCI60.chip(histone_mark = "H3K9",
                                                              siteofchoice = "genes",
                                                              gene_enquiry = "NNMT",
                                                              expressiontype = "protein")

saveRDS(H3K9NNMTprotein_genes_signal_model, "output/NCI60_H3K9me3_genes_signal_model_NNMTprotein.rds")


H3K9NNMTprotein_repeats_signal_model <- fit.model.to.NCI60.chip(histone_mark = "H3K9",
                                                                siteofchoice = "repeats",
                                                                gene_enquiry = "NNMT",
                                                                expressiontype = "protein")

saveRDS(H3K9NNMTprotein_repeats_signal_model, "output/NCI60_H3K9me3_repeats_signal_model_NNMT_protein.rds")

H3K9NNMT_RNAseq_promoters_signal_model <- fit.model.to.NCI60.chip(histone_mark = "H3K9",
                                                                  siteofchoice = "promoters",
                                                                  gene_enquiry = "NNMT",
                                                                  expressiontype = "_RNAseq")

saveRDS(H3K9NNMT_RNAseq_promoters_signal_model, "output/NCI60_H3K9me3_promoters_signal_model_NNMT_RNAseq.rds")

H3K9NNMT_RNAseq_genes_signal_model <- fit.model.to.NCI60.chip(histone_mark = "H3K9",
                                                              siteofchoice = "genes",
                                                              gene_enquiry = "NNMT",
                                                              expressiontype = "_RNAseq")

saveRDS(H3K9NNMT_RNAseq_genes_signal_model, "output/NCI60_H3K9me3_genes_signal_model_NNMT_RNAseq.rds")
H3K9NNMT_RNAseq_genes_signal_model <- readRDS("output/NCI60_H3K9me3_genes_signal_model_NNMT_RNAseq.rds")

H3K9NNMT_RNAseq_repeats_signal_model <- fit.model.to.NCI60.chip(histone_mark = "H3K9",
                                                                siteofchoice = "repeats",
                                                                gene_enquiry = "NNMT",
                                                                expressiontype = "_RNAseq")

saveRDS(H3K9NNMT_RNAseq_repeats_signal_model, "output/NCI60_H3K9me3_repeats_signal_model_NNMT_RNAseq.rds")

H3K9PEMT_RNAseq_promoters_signal_model <- fit.model.to.NCI60.chip(histone_mark = "H3K9",
                                                                  siteofchoice = "promoters",
                                                                  gene_enquiry = "PEMT",
                                                                  expressiontype = "_RNAseq")

saveRDS(H3K9PEMT_RNAseq_promoters_signal_model, "output/NCI60_H3K9me3_promoters_signal_model_PEMT_RNAseq.rds")

H3K9PEMT_RNAseq_genes_signal_model <- fit.model.to.NCI60.chip(histone_mark = "H3K9",
                                                              siteofchoice = "genes",
                                                              gene_enquiry = "PEMT",
                                                              expressiontype = "_RNAseq")

saveRDS(H3K9PEMT_RNAseq_genes_signal_model, "output/NCI60_H3K9me3_genes_signal_model_PEMT_RNAseq.rds")

H3K9PEMT_RNAseq_repeats_signal_model <- fit.model.to.NCI60.chip(histone_mark = "H3K9",
                                                                siteofchoice = "repeats",
                                                                gene_enquiry = "PEMT",
                                                                expressiontype = "_RNAseq")

saveRDS(H3K9PEMT_RNAseq_repeats_signal_model, "output/NCI60_H3K9me3_repeats_signal_model_PEMT_RNAseq.rds")

#### H3K9 random distributions ####

#### RANDOM DISTRIBUTIONS H3K9me3 RNAseq GENE BODIES ####

nuclear_genes_expressed_entrez <- readRDS("output/nuclear_genes_expressed_entrez.rds")

RNAseq_data_expressed <- NCI60_normalised_counts[row.names(NCI60_normalised_counts) %in% nuclear_genes_expressed_entrez$entrezgene_id, ]
row.names(RNAseq_data_expressed) <- nuclear_genes_expressed_entrez[match(row.names(RNAseq_data_expressed), nuclear_genes_expressed_entrez$entrezgene_id), "ensembl_gene_id"]

# random_expressed_mRNA_sample <- sample(colnames(mRNA_expressed_datat), 1000, replace = FALSE)
# saveRDS(random_expressed_mRNA_sample, "output/random_expressed_mRNA_sample.rds")
random_expressed_mRNA_sample <- readRDS("output/random_expressed_mRNA_sample.rds")

for(i in 1:length(random_expressed_mRNA_sample)){
  
  H3K9genesmatrix_for_lm[, paste0(random_expressed_mRNA_sample[i], "_RNAseq")] <- as.numeric(RNAseq_data_expressed[random_expressed_mRNA_sample[i], match(str_remove(row.names(H3K9genesmatrix_for_lm), "H3K9me3_[A-Z]{2,3}_"), colnames(RNAseq_data_expressed))])
  
}

RNAseqrandomsample_H3K9genes_t_values <- sapply(random_expressed_mRNA_sample, function(thisgene){
  
  sapply(colnames(H3K9genesmatrix_for_lm)[1:H3K9genesnumber_of_sites], function(x){
    
    tryCatch(expr = {
      glmtemp <- glm(formula = as.formula(paste0(x, " ~ log10(", as.character(thisgene), "_RNAseq + 1) + tissue + log10(", thisgene, "_RNAseq + 1) * tissue")), data = H3K9genesmatrix_for_lm)
      glmsum <- summary(glmtemp)
      
      t <- glmsum$coefficients[paste0("log10(", thisgene, "_RNAseq + 1)"), "t value"]
      
      output_vec <- c(t)
      names(output_vec) <- c("t_value")
      
      return(output_vec)
      
    }, error = function(e){
      
      output_vec <- c(NA)
      names(output_vec) <- c("t_value")
      
      return(output_vec)
      
    }
    )
    
  })
  
})

saveRDS(RNAseqrandomsample_H3K9genes_t_values, "output/RNAseqrandomsample_H3K9me3genes_t_values.rds")

#### RANDOM DISTRIBUTIONS H3K9me3 RNAseq PROMOTERS ####

RNAseq_data_expressed <- NCI60_normalised_counts[row.names(NCI60_normalised_counts) %in% nuclear_promoters_expressed_entrez$entrezgene_id, ]
row.names(RNAseq_data_expressed) <- nuclear_promoters_expressed_entrez[match(row.names(RNAseq_data_expressed), nuclear_promoters_expressed_entrez$entrezgene_id), "ensembl_gene_id"]

# random_expressed_mRNA_sample <- sample(colnames(mRNA_expressed_datat), 1000, replace = FALSE)
# saveRDS(random_expressed_mRNA_sample, "output/random_expressed_mRNA_sample.rds")
random_expressed_mRNA_sample <- readRDS("output/random_expressed_mRNA_sample.rds")

for(i in 1:length(random_expressed_mRNA_sample)){
  
  H3K9promotersmatrix_for_lm[, paste0(random_expressed_mRNA_sample[i], "_RNAseq")] <- as.numeric(RNAseq_data_expressed[random_expressed_mRNA_sample[i], match(str_remove(row.names(H3K9promotersmatrix_for_lm), "H3K9me3_[A-Z]{2,3}_"), colnames(RNAseq_data_expressed))])
  
}

RNAseqrandomsample_H3K9promoters_t_values <- sapply(random_expressed_mRNA_sample, function(thisgene){
  
  sapply(colnames(H3K9promotersmatrix_for_lm)[1:H3K9promotersnumber_of_sites], function(x){
    
    tryCatch(expr = {
      glmtemp <- glm(formula = as.formula(paste0(x, " ~ log10(", as.character(thisgene), "_RNAseq + 1) + tissue + log10(", thisgene, "_RNAseq + 1) * tissue")), data = H3K9promotersmatrix_for_lm)
      glmsum <- summary(glmtemp)
      
      t <- glmsum$coefficients[paste0("log10(", thisgene, "_RNAseq + 1)"), "t value"]
      
      output_vec <- c(t)
      names(output_vec) <- c("t_value")
      
      return(output_vec)
      
    }, error = function(e){
      
      output_vec <- c(NA)
      names(output_vec) <- c("t_value")
      
      return(output_vec)
      
    }
    )
    
  })
  
})

saveRDS(RNAseqrandomsample_H3K9promoters_t_values, "output/RNAseqrandomsample_H3K9me3promoters_t_values.rds")

#### RANDOM DISTRIBUTION OF H3K9 PROMOTERS MRNA ####

# downloaded from [https://discover.nci.nih.gov/cellminer/loadDownload.do]
all_mRNA_data <- read_excel("input/RNA__5_Platform_Gene_Transcript_Average_z_scores.xls", skip = 10)
names_changed_to_match_format_RNAsheet <- str_remove(str_remove(str_replace(str_replace(str_remove(colnames(all_mRNA_data)[7:ncol(all_mRNA_data)], "^[A-Z]{2,3}:"), "/", "_"), " ", "_"), "\\(TB\\)"), "_ATCC")
names_changed_to_match_format_RNAsheet <- str_replace(names_changed_to_match_format_RNAsheet, pattern = "MDA-N", replacement = "MDA-MB-468")
colnames(all_mRNA_data)[7:ncol(all_mRNA_data)] <- names_changed_to_match_format_RNAsheet

mRNAdatat <- t(all_mRNA_data)
colnames(mRNAdatat) <- mRNAdatat[1,]

all_mRNA_data_expressed <- all_mRNA_data[all_mRNA_data$`Entrez gene id e` %in% nuclear_genes_expressed_entrez$entrezgene_id, ]
all_mRNA_data_expressed[, "ensembl_gene_id"] <- nuclear_genes_expressed_entrez[match(all_mRNA_data_expressed$`Entrez gene id e`, nuclear_genes_expressed_entrez$entrezgene_id), "ensembl_gene_id"]

mRNA_expressed_datat <- t(all_mRNA_data_expressed)
colnames(mRNA_expressed_datat) <- mRNA_expressed_datat["ensembl_gene_id",]

# random_expressed_mRNA_sample <- sample(colnames(mRNA_expressed_datat), 1000, replace = FALSE)
# saveRDS(random_expressed_mRNA_sample, "output/random_expressed_mRNA_sample.rds")
random_expressed_mRNA_sample <- readRDS("output/random_expressed_mRNA_sample.rds")

for(i in 1:length(random_expressed_mRNA_sample)){
  
  H3K9promotersmatrix_for_lm[, paste0(random_expressed_mRNA_sample[i], "_mRNA")] <- as.numeric(mRNA_expressed_datat[match(str_remove(row.names(H3K9promotersmatrix_for_lm), "H3K9me3_[A-Z]{2,3}_"), row.names(mRNA_expressed_datat)), random_expressed_mRNA_sample[i]])
  
}

mRNArandomsample_H3K9me3promoters_t_values <- sapply(random_expressed_mRNA_sample, function(thisgene){
  
  sapply(colnames(H3K9promotersmatrix_for_lm)[1:H3K9promotersnumber_of_sites], function(x){
    
    tryCatch(expr = {
      glmtemp <- glm(formula = as.formula(paste0(x, " ~ ", as.character(thisgene), "_mRNA + tissue + ", thisgene, "_mRNA * tissue")), data = H3K9promotersmatrix_for_lm)
      glmsum <- summary(glmtemp)
      
      t <- glmsum$coefficients[paste0(thisgene, "_mRNA"), "t value"]
      
      output_vec <- c(t)
      names(output_vec) <- c("t_value")
      
      return(output_vec)
      
    }, error = function(e){
      
      output_vec <- c(NA)
      names(output_vec) <- c("t_value")
      
      return(output_vec)
      
    }
    )
    
  })
  
})

saveRDS(mRNArandomsample_H3K9me3promoters_t_values, "output/mRNArandomsample_H3K9me3promoters_t_values.rds")
# mRNArandomsample_H3K9me3promoters_t_values <- readRDS("output/mRNArandomsample_H3K9me3promoters_t_values.rds")

#### RANDOM DISTRIBUTION OF H3K9 GENES MRNA ####

# random_expressed_mRNA_sample <- sample(colnames(mRNA_expressed_datat), 1000, replace = FALSE)
# saveRDS(random_expressed_mRNA_sample, "output/random_expressed_mRNA_sample.rds")
random_expressed_mRNA_sample <- readRDS("output/random_expressed_mRNA_sample.rds")

for(i in 1:length(random_expressed_mRNA_sample)){
  
  H3K9genesmatrix_for_lm[, paste0(random_expressed_mRNA_sample[i], "_mRNA")] <- as.numeric(mRNA_expressed_datat[match(str_remove(row.names(H3K9genesmatrix_for_lm), "H3K9me3_[A-Z]{2,3}_"), row.names(mRNA_expressed_datat)), random_expressed_mRNA_sample[i]])
  
}

mRNArandomsample_H3K9me3genes_t_values <- sapply(random_expressed_mRNA_sample, function(thisgene){
  
  sapply(colnames(H3K9genesmatrix_for_lm)[1:H3K9genesnumber_of_sites], function(x){
    
    tryCatch(expr = {
      glmtemp <- glm(formula = as.formula(paste0(x, " ~ ", as.character(thisgene), "_mRNA + tissue + ", thisgene, "_mRNA * tissue")), data = H3K9genesmatrix_for_lm)
      glmsum <- summary(glmtemp)
      
      t <- glmsum$coefficients[paste0(thisgene, "_mRNA"), "t value"]
      
      output_vec <- c(t)
      names(output_vec) <- c("t_value")
      
      return(output_vec)
      
    }, error = function(e){
      
      output_vec <- c(NA)
      names(output_vec) <- c("t_value")
      
      return(output_vec)
      
    }
    )
    
  })
  
})

saveRDS(mRNArandomsample_H3K9me3genes_t_values, "output/mRNArandomsample_H3K9me3genes_t_values.rds")
mRNArandomsample_H3K9me3genes_t_values <- readRDS("output/mRNArandomsample_H3K9me3genes_t_values.rds")

#### RANDOM DISTRIBUTION OF H3K9 PROMOTERS PROTEIN ####

# downloaded from [https://discover.nci.nih.gov/cellminer/loadDownload.do]
all_protein_data <- read_excel("input/Protein__SWATH_(Mass_spectrometry)_Protein.xls", skip = 10)

names_changed_to_match_format_protsheet <- str_remove(str_remove(str_replace(str_replace(str_remove(colnames(all_protein_data)[10:ncol(all_protein_data)], "^[A-Z]{2,3}:"), "/", "_"), " ", "_"), "\\(TB\\)"), "_ATCC")
colnames(all_protein_data)[10:ncol(all_protein_data)] <- names_changed_to_match_format_protsheet

protdatat <- t(all_protein_data)
colnames(protdatat) <- protdatat[2, ]

# random_detected_protein_sample <- sample(all_protein_data$`Gene name d`, 1000, replace = FALSE)
# saveRDS(random_detected_protein_sample, "output/random_detected_protein_sample.rds" )
random_detected_protein_sample <- readRDS("output/random_detected_protein_sample.rds" )

for(i in 1:length(random_detected_protein_sample)){
  
  H3K9promotersmatrix_for_lm[, paste0(random_detected_protein_sample[i], "_protein")] <- as.numeric(protdatat[match(str_remove(row.names(H3K9promotersmatrix_for_lm), "H3K9me3_[A-Z]{2,3}_"), row.names(protdatat)), random_detected_protein_sample[i]])
  
}

proteinrandomsample_H3K9me3promoters_t_values <- sapply(random_detected_protein_sample, function(thisgene){
  
  sapply(colnames(H3K9promotersmatrix_for_lm)[1:H3K9promotersnumber_of_sites], function(x){
    
    tryCatch(expr = {
      glmtemp <- glm(formula = as.formula(paste0(x, " ~ ", as.character(thisgene), "_protein + tissue + ", thisgene, "_protein * tissue")), data = H3K9promotersmatrix_for_lm)
      glmsum <- summary(glmtemp)
      
      t <- glmsum$coefficients[paste0(thisgene, "_protein"), "t value"]
      
      output_vec <- c(t)
      names(output_vec) <- c("t_value")
      
      return(output_vec)
      
    }, error = function(e){
      
      output_vec <- c(NA)
      names(output_vec) <- c("t_value")
      
      return(output_vec)
      
    }
    )
    
  })
  
})

saveRDS(proteinrandomsample_H3K9me3promoters_t_values, "output/proteinrandomsample_H3K9me3promoters_t_values.rds")
# proteinrandomsample_H3K9me3promoters_t_values <- readRDS("output/proteinrandomsample_H3K9me3promoters_t_values.rds")

#### RANDOM DISTRIBUTION OF H3K9 GENE BODIES PROTEIN ####

random_detected_protein_sample <- readRDS("output/random_detected_protein_sample.rds" )

for(i in 1:length(random_detected_protein_sample)){
  
  H3K9genesmatrix_for_lm[, paste0(random_detected_protein_sample[i], "_protein")] <- as.numeric(protdatat[match(str_remove(row.names(H3K9genesmatrix_for_lm), "H3K9me3_[A-Z]{2,3}_"), row.names(protdatat)), random_detected_protein_sample[i]])
  
}

proteinrandomsample_H3K9me3genes_t_values <- sapply(random_detected_protein_sample, function(thisgene){
  
  sapply(colnames(H3K9genesmatrix_for_lm)[1:H3K9genesnumber_of_sites], function(x){
    
    tryCatch(expr = {
      glmtemp <- glm(formula = as.formula(paste0(x, " ~ ", as.character(thisgene), "_protein + tissue + ", thisgene, "_protein * tissue")), data = H3K9genesmatrix_for_lm)
      glmsum <- summary(glmtemp)
      
      t <- glmsum$coefficients[paste0(thisgene, "_protein"), "t value"]
      
      output_vec <- c(t)
      names(output_vec) <- c("t_value")
      
      return(output_vec)
      
    }, error = function(e){
      
      output_vec <- c(NA)
      names(output_vec) <- c("t_value")
      
      return(output_vec)
      
    }
    )
    
  })
  
})

saveRDS(proteinrandomsample_H3K9me3genes_t_values, "output/proteinrandomsample_H3K9me3genes_t_values.rds")
# proteinrandomsample_H3K9me3genes_t_values <- readRDS("output/proteinrandomsample_H3K9me3genes_t_values.rds")

#### H4K20me3 ####

#### combine narrow and broad peak options####

# point to directory with genrich output
my_directory <- "input/NCI60_genrich_peak_calls/"

genrich_file_list <- list.files(my_directory)

# choose the narrowPeak files
genrich_file_list_narrowPeak <- genrich_file_list[str_detect(genrich_file_list, pattern = "\\.narrowPeak")]

# choose the H4K20me3 files
H4K20me3_file_list_narrowPeak <- genrich_file_list_narrowPeak[str_detect(genrich_file_list_narrowPeak, pattern = "H4K20me3")]

# loop through unique experiments to pool the two replicates

H4K20me3_ranges_list <- lapply(H4K20me3_file_list_narrowPeak, function(x){
  
  peaks_1 <- read.table(paste0(my_directory, x), header = FALSE)
  colnames(peaks_1)[1:3] <- c("chr", "start", "end")
  
  peaks_1_GR <- makeGRangesFromDataFrame(peaks_1)
  
  
})

names(H4K20me3_ranges_list) <- unique(str_remove(H4K20me3_file_list_narrowPeak, pattern = "_genrich.*"))

saveRDS(H4K20me3_ranges_list, "output/NCI60_H4K20me3_combinedpeaks_list.rds")
# H4K20me3_ranges_list <- readRDS("output/NCI60_H4K20me3_combinedpeaks_list.rds")

#### SELECT H4K20me3 REGIONS WITH CONSISTENT PEAKS ####

H4K20me3_gene_peak_counts <- region.peak.counts(NCI60_list = H4K20me3_ranges_list,
                                                region = "genes")

H4K20genes_in_40of60_number <- apply(H4K20me3_gene_peak_counts[, 5:ncol(H4K20me3_gene_peak_counts)], 1, function(x){
  
  sum(x > 0)
  
})

H4K20genes_in_40of60_selec <- H4K20genes_in_40of60_number > 39
H4K20genes_in_40of60 <- H4K20me3_gene_peak_counts[H4K20genes_in_40of60_selec, "gene_id"]

H4K20genes_in_all <- H4K20me3_gene_peak_counts[H4K20genes_in_40of60_number == 60, "gene_id"]

# total of 92 genes with H4K20me3 peaks in most samples
# total of 11 genes with H4K20me3 peaks in most samples

H4K20me3_promoter_peak_counts <- region.peak.counts(NCI60_list = H4K20me3_ranges_list,
                                                    region = "promoters")

H4K20promoters_in_40of60_number <- apply(H4K20me3_promoter_peak_counts[, 5:ncol(H4K20me3_promoter_peak_counts)], 1, function(x){
  
  sum(x > 0)
  
})

H4K20promoters_in_40of60_selec <- H4K20promoters_in_40of60_number > 39
H4K20promoters_in_40of60 <- H4K20me3_promoter_peak_counts[H4K20promoters_in_40of60_selec, "gene_id"]

H4K20promoters_in_all <- H4K20me3_promoter_peak_counts[H4K20promoters_in_40of60_number == 60, "gene_id"]

# total of 40 promoters with H4K20me3 peaks in most samples
# total of 1 promoters with H4K20me3 peaks in every sample

#### MAKE H4K20me3 MATRIX FOR LM (GENEBODIES) ####

H4K20me3_computematrix_files <- list.files("input/NCI60_H4K20_computematrix/")
H4K20me3_genes_computematrix_files <- H4K20me3_computematrix_files[str_detect(H4K20me3_computematrix_files, "genebodies")]

uniqueexperiments <- unique(str_remove(str_remove(H4K20me3_genes_computematrix_files, "_bamcompare.bigWig_genebodies200bpflanks_computematrix.gz"), "_rep2"))

H4K20me3_genes_experimentavg_list <- lapply(uniqueexperiments, function(x){
  
  replicatexperiments <- H4K20me3_genes_computematrix_files[str_detect(H4K20me3_genes_computematrix_files, paste0(x, "_"))]
  
  if(length(replicatexperiments) > 2){message(paste0("Warning! More than 2 files detected for this experiment (", x, ")"))}
  
  rep1 <- read.table(paste0("input/NCI60_H4K20_computematrix/", replicatexperiments[1]), skip = 1)
  rep2 <- read.table(paste0("input/NCI60_H4K20_computematrix/", replicatexperiments[2]), skip = 1)
  
  combined <- rep1
  combined$V7 <- rowMeans(cbind(rep1$V7, rep2$V7))
  
  return(combined)
  
})

names(H4K20me3_genes_experimentavg_list) <- uniqueexperiments

H4K20genesreadmatrix_df <- sapply(H4K20me3_genes_experimentavg_list, function(x){
  
  x$V7
  
})

row.names(H4K20genesreadmatrix_df) <- H4K20me3_genes_experimentavg_list[[1]]$V4

H4K20genesmatrix_for_lm <- data.frame(t(H4K20genesreadmatrix_df))

H4K20genesmatrix_for_lm <- H4K20genesmatrix_for_lm[, str_remove(colnames(H4K20genesmatrix_for_lm), "^X") %in% H4K20genes_in_40of60]

H4K20genesmatrix_for_lm[ , "tissue"] <- str_extract(str_remove(row.names(H4K20genesmatrix_for_lm), "H4K20me3_"), "[A-Z]{2,3}")

NNMT_proteomics <- read.table("input/NCI60_NNMT_proteomics.txt", sep = "\t", header = TRUE)
NNMT_mRNA <- read.table("input/NCI60_NNMT_mRNA_zscores.txt", sep = "\t", header = TRUE)

# here I am changing the names from the NNMT proteomics and gene expression file to match the format
# replace slashes with underscores
# replace spaces with underscores
# remove "(TB)" from HL-60
# remove "_ATCC" from A549
names_changed_to_match_format_protein <- str_remove(str_remove(str_replace(str_replace(str_remove(NNMT_proteomics$Cell.Line, "^[A-Z]{2,3}:"), "/", "_"), " ", "_"), "\\(TB\\)"), "_ATCC")
names_changed_to_match_format_RNA <- str_remove(str_remove(str_replace(str_replace(str_remove(NNMT_mRNA$Cell.Line, "^[A-Z]{2,3}:"), "/", "_"), " ", "_"), "\\(TB\\)"), "_ATCC")

# why 59?
# only 59 qvailable data points for proteomics
# 60 for expression so will download expression data separately from CellMinerDB

# the 60th in the RNA dataset is "MDA-N". Will assume it is the last one in the panel which is MDA-MB-468
names_changed_to_match_format_RNA[!sapply(names_changed_to_match_format_RNA, function(x){
  any(str_detect(row.names(H4K20genesmatrix_for_lm), x))
})]

names_changed_to_match_format_RNA <- str_replace(names_changed_to_match_format_RNA, pattern = "MDA-N", replacement = "MDA-MB-468")

NNMT_proteomics$Cell.Line <- names_changed_to_match_format_protein
NNMT_mRNA$Cell.Line <- names_changed_to_match_format_RNA

H4K20genesmatrix_for_lm[, "NNMTprotein"] <- NNMT_proteomics[match(str_remove(str_remove(row.names(H4K20genesmatrix_for_lm), pattern = "H4K20me3_"), pattern = "^[A-Z]{2,3}_"), NNMT_proteomics$Cell.Line), "swaNNMT_nci60"]

H4K20genesmatrix_for_lm[, "NNMTmRNA_z"] <- NNMT_mRNA[match(str_remove(str_remove(row.names(H4K20genesmatrix_for_lm), pattern = "H4K20me3_"), pattern = "^[A-Z]{2,3}_"), NNMT_mRNA$Cell.Line), "expNNMT_nci60"]

H4K20genesmatrix_for_lm[, "NNMT_RNAseq"] <- NCI60_normalised_counts[NNMT_entrezgene_id, match(str_remove(str_remove(row.names(H4K20genesmatrix_for_lm), pattern = "H4K20me3_"), pattern = "^[A-Z]{2,3}_"), colnames(NCI60_normalised_counts))]
H4K20genesmatrix_for_lm[, "PEMT_RNAseq"] <- NCI60_normalised_counts[PEMT_entrezgene_id, match(str_remove(str_remove(row.names(H4K20genesmatrix_for_lm), pattern = "H4K20me3_"), pattern = "^[A-Z]{2,3}_"), colnames(NCI60_normalised_counts))]

H4K20genesnumber_of_sites <- ncol(H4K20genesmatrix_for_lm) - 5

saveRDS(H4K20genesmatrix_for_lm, "output/NCI60_H4K20me3_genes_signal_matrix_for_lm.rds")
# H4K20genesmatrix_for_lm <- readRDS("output/NCI60_H4K20me3_genes_signal_matrix_for_lm.rds")

#### MAKE H4K20me3 MATRIX FOR LM (PROMOTERS) ####

# point to deeptools computematrix output directory
computematrix_files <- list.files("input/computematrix_out/")
H4K20me3_computematrix_files <- computematrix_files[str_detect(computematrix_files, "H4K20me3")] 
H4K20me3_promoters_computematrix_files <- H4K20me3_computematrix_files[str_detect(H4K20me3_computematrix_files, "promoter")]

uniqueexperiments <- unique(str_remove(str_remove(H4K20me3_promoters_computematrix_files, "_bamcompare.bigWig_promoters_computematrix.gz"), "_rep2"))

H4K20me3_promoters_experimentavg_list <- lapply(uniqueexperiments, function(x){
  
  replicatexperiments <- H4K20me3_promoters_computematrix_files[str_detect(H4K20me3_promoters_computematrix_files, paste0(x, "_"))]
  
  if(length(replicatexperiments) > 2){message(paste0("Warning! More than 2 files detected for this experiment (", x, ")"))}
  
  rep1 <- read.table(paste0("input/computematrix_out/", replicatexperiments[1]), skip = 1)
  rep2 <- read.table(paste0("input/computematrix_out/", replicatexperiments[2]), skip = 1)
  
  combined <- rep1
  combined$V7 <- rowMeans(cbind(rep1$V7, rep2$V7))
  
  return(combined)
  
})

names(H4K20me3_promoters_experimentavg_list) <- uniqueexperiments

H4K20promotersreadmatrix_df <- sapply(H4K20me3_promoters_experimentavg_list, function(x){
  
  x$V7
  
})

row.names(H4K20promotersreadmatrix_df) <- H4K20me3_promoters_experimentavg_list[[1]]$V4

H4K20promotersmatrix_for_lm <- data.frame(t(H4K20promotersreadmatrix_df))

H4K20promotersmatrix_for_lm <- H4K20promotersmatrix_for_lm[, str_remove(colnames(H4K20promotersmatrix_for_lm), "^X") %in% H4K20promoters_in_40of60]

H4K20promotersmatrix_for_lm[ , "tissue"] <- str_extract(str_remove(row.names(H4K20promotersmatrix_for_lm), "H4K20me3_"), "[A-Z]{2,3}")

H4K20promotersmatrix_for_lm[, "NNMTprotein"] <- NNMT_proteomics[match(str_remove(str_remove(row.names(H4K20promotersmatrix_for_lm), pattern = "H4K20me3_"), pattern = "^[A-Z]{2,3}_"), NNMT_proteomics$Cell.Line), "swaNNMT_nci60"]

H4K20promotersmatrix_for_lm[, "NNMTmRNA_z"] <- NNMT_mRNA[match(str_remove(str_remove(row.names(H4K20promotersmatrix_for_lm), pattern = "H4K20me3_"), pattern = "^[A-Z]{2,3}_"), NNMT_mRNA$Cell.Line), "expNNMT_nci60"]

H4K20promotersmatrix_for_lm[, "NNMT_RNAseq"] <- NCI60_normalised_counts[NNMT_entrezgene_id, match(str_remove(str_remove(row.names(H4K20promotersmatrix_for_lm), pattern = "H4K20me3_"), pattern = "^[A-Z]{2,3}_"), colnames(NCI60_normalised_counts))]
H4K20promotersmatrix_for_lm[, "PEMT_RNAseq"] <- NCI60_normalised_counts[PEMT_entrezgene_id, match(str_remove(str_remove(row.names(H4K20promotersmatrix_for_lm), pattern = "H4K20me3_"), pattern = "^[A-Z]{2,3}_"), colnames(NCI60_normalised_counts))]

H4K20promotersnumber_of_sites <- ncol(H4K20promotersmatrix_for_lm) - 5

saveRDS(H4K20promotersmatrix_for_lm, "output/NCI60_H4K20me3_promoters_signal_matrix_for_lm.rds")
# H4K20promotersmatrix_for_lm <- readRDS("output/NCI60_H4K20me3_promoters_signal_matrix_for_lm.rds")

#### NCI60 BEDS FOR REPETITIVE ELEMENTS (FOR H4K20me3) ####

# limit repetitive elements to clearly annotated categories and convert to genomic ranges object
# repetitive elements file for hg38
rmask <- read.table("input/rmsk.txt", sep = "\t", header = FALSE, quote = "")
colnames(rmask)[c(6,7,8)] <- c("chr", "start", "end")
rmask[, "myID"] <- paste0(rmask$chr, rmask$start, rmask$end)
row.names(rmask) <- rmask$myID

rmask_gr <- makeGRangesFromDataFrame(rmask[rmask$V12 %in% c("Satellite", "LINE", "SINE", "LTR", "DNA", "Low_complexity", "snRNA", "tRNA", "Retroposon"), c("chr", "start", "end", "myID", "V11", "V12", "V13")], keep.extra.columns = TRUE)

NCI60_H4K20_repeats_sites <- region.peak.counts(NCI60_list = H4K20me3_ranges_list,
                                                region = "repeats")

rmask_NCI60_H4K20_select <- rmask_gr[NCI60_H4K20_repeats_sites , ]

saveRDS(rmask_NCI60_H4K20_select, "output/rmask_40ofNCI60peaks_H4K20me3.rds")
# rmask_NCI60_H4K20_select <- readRDS("output/rmask_40ofNCI60peaks_H4K20me3.rds")

rmask_NCI60_H4K20_select_df <- data.frame(seqnames = seqnames(rmask_NCI60_H4K20_select),
                                          starts = start(rmask_NCI60_H4K20_select)-1,
                                          ends = end(rmask_NCI60_H4K20_select),
                                          gene_id = rmask_NCI60_H4K20_select$myID,
                                          scores = c(rep(".", length(rmask_NCI60_H4K20_select))),
                                          strands = strand(rmask_NCI60_H4K20_select))

options(scipen=999)

write.table(rmask_NCI60_H4K20_select_df, file = "output/rmask_40ofNCI60peaks_H4K20me3.bed", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

#### FIT H4K20 NNMT MODELS ####

H4K20NNMTmRNA_promoters_signal_model <- fit.model.to.NCI60.chip(histone_mark = "H4K20",
                                                                siteofchoice = "promoters",
                                                                gene_enquiry = "NNMT",
                                                                expressiontype = "mRNA_z")

saveRDS(H4K20NNMTmRNA_promoters_signal_model, "output/NCI60_H4K20me3_promoters_signal_model_NNMTmRNA.rds")

H4K20NNMTmRNA_genes_signal_model <- fit.model.to.NCI60.chip(histone_mark = "H4K20",
                                                            siteofchoice = "genes",
                                                            gene_enquiry = "NNMT",
                                                            expressiontype = "mRNA_z")

saveRDS(H4K20NNMTmRNA_genes_signal_model, "output/NCI60_H4K20me3_genes_signal_model_NNMTmRNA.rds")

H4K20NNMTmRNA_repeats_signal_model <- fit.model.to.NCI60.chip(histone_mark = "H4K20",
                                                              siteofchoice = "repeats",
                                                              gene_enquiry = "NNMT",
                                                              expressiontype = "mRNA_z")

saveRDS(H4K20NNMTmRNA_repeats_signal_model, "output/NCI60_H4K20me3_repeats_signal_model_NNMT_mRNA.rds")

H4K20NNMTprotein_promoters_signal_model <- fit.model.to.NCI60.chip(histone_mark = "H4K20",
                                                                   siteofchoice = "promoters",
                                                                   gene_enquiry = "NNMT",
                                                                   expressiontype = "protein")

saveRDS(H4K20NNMTprotein_promoters_signal_model, "output/NCI60_H4K20me3_promoters_signal_model_NNMTprotein.rds")

H4K20NNMTprotein_genes_signal_model <- fit.model.to.NCI60.chip(histone_mark = "H4K20",
                                                               siteofchoice = "genes",
                                                               gene_enquiry = "NNMT",
                                                               expressiontype = "protein")

saveRDS(H4K20NNMTprotein_genes_signal_model, "output/NCI60_H4K20me3_genes_signal_model_NNMTprotein.rds")


H4K20NNMTprotein_repeats_signal_model <- fit.model.to.NCI60.chip(histone_mark = "H4K20",
                                                                 siteofchoice = "repeats",
                                                                 gene_enquiry = "NNMT",
                                                                 expressiontype = "protein")

saveRDS(H4K20NNMTprotein_repeats_signal_model, "output/NCI60_H4K20me3_repeats_signal_model_NNMT_protein.rds")

H4K20NNMT_RNAseq_promoters_signal_model <- fit.model.to.NCI60.chip(histone_mark = "H4K20",
                                                                   siteofchoice = "promoters",
                                                                   gene_enquiry = "NNMT",
                                                                   expressiontype = "_RNAseq")

saveRDS(H4K20NNMT_RNAseq_promoters_signal_model, "output/NCI60_H4K20me3_promoters_signal_model_NNMT_RNAseq.rds")

H4K20NNMT_RNAseq_genes_signal_model <- fit.model.to.NCI60.chip(histone_mark = "H4K20",
                                                               siteofchoice = "genes",
                                                               gene_enquiry = "NNMT",
                                                               expressiontype = "_RNAseq")

saveRDS(H4K20NNMT_RNAseq_genes_signal_model, "output/NCI60_H4K20me3_genes_signal_model_NNMT_RNAseq.rds")
H4K20NNMT_RNAseq_genes_signal_model <- readRDS("output/NCI60_H4K20me3_genes_signal_model_NNMT_RNAseq.rds")

H4K20NNMT_RNAseq_repeats_signal_model <- fit.model.to.NCI60.chip(histone_mark = "H4K20",
                                                                 siteofchoice = "repeats",
                                                                 gene_enquiry = "NNMT",
                                                                 expressiontype = "_RNAseq")

saveRDS(H4K20NNMT_RNAseq_repeats_signal_model, "output/NCI60_H4K20me3_repeats_signal_model_NNMT_RNAseq.rds")

H4K20PEMT_RNAseq_promoters_signal_model <- fit.model.to.NCI60.chip(histone_mark = "H4K20",
                                                                   siteofchoice = "promoters",
                                                                   gene_enquiry = "PEMT",
                                                                   expressiontype = "_RNAseq")

saveRDS(H4K20PEMT_RNAseq_promoters_signal_model, "output/NCI60_H4K20me3_promoters_signal_model_PEMT_RNAseq.rds")

H4K20PEMT_RNAseq_genes_signal_model <- fit.model.to.NCI60.chip(histone_mark = "H4K20",
                                                               siteofchoice = "genes",
                                                               gene_enquiry = "PEMT",
                                                               expressiontype = "_RNAseq")

saveRDS(H4K20PEMT_RNAseq_genes_signal_model, "output/NCI60_H4K20me3_genes_signal_model_PEMT_RNAseq.rds")

H4K20PEMT_RNAseq_repeats_signal_model <- fit.model.to.NCI60.chip(histone_mark = "H4K20",
                                                                 siteofchoice = "repeats",
                                                                 gene_enquiry = "PEMT",
                                                                 expressiontype = "_RNAseq")

saveRDS(H4K20PEMT_RNAseq_repeats_signal_model, "output/NCI60_H4K20me3_repeats_signal_model_PEMT_RNAseq.rds")

#### H4K20 random distributions ####

#### RANDOM DISTRIBUTIONS H4K20me3 RNAseq GENE BODIES ####

nuclear_genes_expressed_entrez <- readRDS("output/nuclear_genes_expressed_entrez.rds")

RNAseq_data_expressed <- NCI60_normalised_counts[row.names(NCI60_normalised_counts) %in% nuclear_genes_expressed_entrez$entrezgene_id, ]
row.names(RNAseq_data_expressed) <- nuclear_genes_expressed_entrez[match(row.names(RNAseq_data_expressed), nuclear_genes_expressed_entrez$entrezgene_id), "ensembl_gene_id"]

# random_expressed_mRNA_sample <- sample(colnames(mRNA_expressed_datat), 1000, replace = FALSE)
# saveRDS(random_expressed_mRNA_sample, "output/random_expressed_mRNA_sample.rds")
random_expressed_mRNA_sample <- readRDS("output/random_expressed_mRNA_sample.rds")

for(i in 1:length(random_expressed_mRNA_sample)){
  
  H4K20genesmatrix_for_lm[, paste0(random_expressed_mRNA_sample[i], "_RNAseq")] <- as.numeric(RNAseq_data_expressed[random_expressed_mRNA_sample[i], match(str_remove(row.names(H4K20genesmatrix_for_lm), "H4K20me3_[A-Z]{2,3}_"), colnames(RNAseq_data_expressed))])
  
}

RNAseqrandomsample_H4K20genes_t_values <- sapply(random_expressed_mRNA_sample, function(thisgene){
  
  sapply(colnames(H4K20genesmatrix_for_lm)[1:H4K20genesnumber_of_sites], function(x){
    
    tryCatch(expr = {
      glmtemp <- glm(formula = as.formula(paste0(x, " ~ log10(", as.character(thisgene), "_RNAseq + 1) + tissue + log10(", thisgene, "_RNAseq + 1) * tissue")), data = H4K20genesmatrix_for_lm)
      glmsum <- summary(glmtemp)
      
      t <- glmsum$coefficients[paste0("log10(", thisgene, "_RNAseq + 1)"), "t value"]
      
      output_vec <- c(t)
      names(output_vec) <- c("t_value")
      
      return(output_vec)
      
    }, error = function(e){
      
      output_vec <- c(NA)
      names(output_vec) <- c("t_value")
      
      return(output_vec)
      
    }
    )
    
  })
  
})

saveRDS(RNAseqrandomsample_H4K20genes_t_values, "output/RNAseqrandomsample_H4K20me3genes_t_values.rds")

#### RANDOM DISTRIBUTIONS H4K20me3 RNAseq PROMOTERS ####

nuclear_genes_expressed_entrez <- getBM(attributes = c("ensembl_gene_id", "entrezgene_id"),
                                        values = nuclear_promoters_expressed,
                                        filter = "ensembl_gene_id",
                                        mart = ensembl)

RNAseq_data_expressed <- NCI60_normalised_counts[row.names(NCI60_normalised_counts) %in% nuclear_promoters_expressed_entrez$entrezgene_id, ]
row.names(RNAseq_data_expressed) <- nuclear_promoters_expressed_entrez[match(row.names(RNAseq_data_expressed), nuclear_promoters_expressed_entrez$entrezgene_id), "ensembl_gene_id"]

random_expressed_mRNA_sample <- sample(colnames(mRNA_expressed_datat), 1000, replace = FALSE)
saveRDS(random_expressed_mRNA_sample, "output/random_expressed_mRNA_sample.rds")
# random_expressed_mRNA_sample <- readRDS("output/random_expressed_mRNA_sample.rds")

for(i in 1:length(random_expressed_mRNA_sample)){
  
  H4K20promotersmatrix_for_lm[, paste0(random_expressed_mRNA_sample[i], "_RNAseq")] <- as.numeric(RNAseq_data_expressed[random_expressed_mRNA_sample[i], match(str_remove(row.names(H4K20promotersmatrix_for_lm), "H4K20me3_[A-Z]{2,3}_"), colnames(RNAseq_data_expressed))])
  
}

RNAseqrandomsample_H4K20promoters_t_values <- sapply(random_expressed_mRNA_sample, function(thisgene){
  
  sapply(colnames(H4K20promotersmatrix_for_lm)[1:H4K20promotersnumber_of_sites], function(x){
    
    tryCatch(expr = {
      glmtemp <- glm(formula = as.formula(paste0(x, " ~ log10(", as.character(thisgene), "_RNAseq + 1) + tissue + log10(", thisgene, "_RNAseq + 1) * tissue")), data = H4K20promotersmatrix_for_lm)
      glmsum <- summary(glmtemp)
      
      t <- glmsum$coefficients[paste0("log10(", thisgene, "_RNAseq + 1)"), "t value"]
      
      output_vec <- c(t)
      names(output_vec) <- c("t_value")
      
      return(output_vec)
      
    }, error = function(e){
      
      output_vec <- c(NA)
      names(output_vec) <- c("t_value")
      
      return(output_vec)
      
    }
    )
    
  })
  
})

saveRDS(RNAseqrandomsample_H4K20promoters_t_values, "output/RNAseqrandomsample_H4K20me3promoters_t_values.rds")

#### RANDOM DISTRIBUTION OF H4K20 PROMOTERS MRNA ####

# downloaded from [https://discover.nci.nih.gov/cellminer/loadDownload.do]
all_mRNA_data <- read_excel("input/RNA__5_Platform_Gene_Transcript_Average_z_scores.xls", skip = 10)
names_changed_to_match_format_RNAsheet <- str_remove(str_remove(str_replace(str_replace(str_remove(colnames(all_mRNA_data)[7:ncol(all_mRNA_data)], "^[A-Z]{2,3}:"), "/", "_"), " ", "_"), "\\(TB\\)"), "_ATCC")
names_changed_to_match_format_RNAsheet <- str_replace(names_changed_to_match_format_RNAsheet, pattern = "MDA-N", replacement = "MDA-MB-468")
colnames(all_mRNA_data)[7:ncol(all_mRNA_data)] <- names_changed_to_match_format_RNAsheet

mRNAdatat <- t(all_mRNA_data)
colnames(mRNAdatat) <- mRNAdatat[1,]

all_mRNA_data_expressed <- all_mRNA_data[all_mRNA_data$`Entrez gene id e` %in% nuclear_genes_expressed_entrez$entrezgene_id, ]
all_mRNA_data_expressed[, "ensembl_gene_id"] <- nuclear_genes_expressed_entrez[match(all_mRNA_data_expressed$`Entrez gene id e`, nuclear_genes_expressed_entrez$entrezgene_id), "ensembl_gene_id"]

mRNA_expressed_datat <- t(all_mRNA_data_expressed)
colnames(mRNA_expressed_datat) <- mRNA_expressed_datat["ensembl_gene_id",]

# random_expressed_mRNA_sample <- sample(colnames(mRNA_expressed_datat), 1000, replace = FALSE)
# saveRDS(random_expressed_mRNA_sample, "output/random_expressed_mRNA_sample.rds")
random_expressed_mRNA_sample <- readRDS("output/random_expressed_mRNA_sample.rds")

for(i in 1:length(random_expressed_mRNA_sample)){
  
  H4K20promotersmatrix_for_lm[, paste0(random_expressed_mRNA_sample[i], "_mRNA")] <- as.numeric(mRNA_expressed_datat[match(str_remove(row.names(H4K20promotersmatrix_for_lm), "H4K20me3_[A-Z]{2,3}_"), row.names(mRNA_expressed_datat)), random_expressed_mRNA_sample[i]])
  
}

mRNArandomsample_H4K20me3promoters_t_values <- sapply(random_expressed_mRNA_sample, function(thisgene){
  
  sapply(colnames(H4K20promotersmatrix_for_lm)[1:H4K20promotersnumber_of_sites], function(x){
    
    tryCatch(expr = {
      glmtemp <- glm(formula = as.formula(paste0(x, " ~ ", as.character(thisgene), "_mRNA + tissue + ", thisgene, "_mRNA * tissue")), data = H4K20promotersmatrix_for_lm)
      glmsum <- summary(glmtemp)
      
      t <- glmsum$coefficients[paste0(thisgene, "_mRNA"), "t value"]
      
      output_vec <- c(t)
      names(output_vec) <- c("t_value")
      
      return(output_vec)
      
    }, error = function(e){
      
      output_vec <- c(NA)
      names(output_vec) <- c("t_value")
      
      return(output_vec)
      
    }
    )
    
  })
  
})

saveRDS(mRNArandomsample_H4K20me3promoters_t_values, "output/mRNArandomsample_H4K20me3promoters_t_values.rds")
# mRNArandomsample_H4K20me3promoters_t_values <- readRDS("output/mRNArandomsample_H4K20me3promoters_t_values.rds")

#### RANDOM DISTRIBUTION OF H4K20 GENES MRNA ####

random_expressed_mRNA_sample <- readRDS("output/random_expressed_mRNA_sample.rds")

for(i in 1:length(random_expressed_mRNA_sample)){
  
  H4K20genesmatrix_for_lm[, paste0(random_expressed_mRNA_sample[i], "_mRNA")] <- as.numeric(mRNA_expressed_datat[match(str_remove(row.names(H4K20genesmatrix_for_lm), "H4K20me3_[A-Z]{2,3}_"), row.names(mRNA_expressed_datat)), random_expressed_mRNA_sample[i]])
  
}

mRNArandomsample_H4K20me3genes_t_values <- sapply(random_expressed_mRNA_sample, function(thisgene){
  
  sapply(colnames(H4K20genesmatrix_for_lm)[1:H4K20genesnumber_of_sites], function(x){
    
    tryCatch(expr = {
      glmtemp <- glm(formula = as.formula(paste0(x, " ~ ", as.character(thisgene), "_mRNA + tissue + ", thisgene, "_mRNA * tissue")), data = H4K20genesmatrix_for_lm)
      glmsum <- summary(glmtemp)
      
      t <- glmsum$coefficients[paste0(thisgene, "_mRNA"), "t value"]
      
      output_vec <- c(t)
      names(output_vec) <- c("t_value")
      
      return(output_vec)
      
    }, error = function(e){
      
      output_vec <- c(NA)
      names(output_vec) <- c("t_value")
      
      return(output_vec)
      
    }
    )
    
  })
  
})

saveRDS(mRNArandomsample_H4K20me3genes_t_values, "output/mRNArandomsample_H4K20me3genes_t_values.rds")
mRNArandomsample_H4K20me3genes_t_values <- readRDS("output/mRNArandomsample_H4K20me3genes_t_values.rds")

#### RANDOM DISTRIBUTION OF H4K20 PROMOTERS PROTEIN ####

# downloaded from [https://discover.nci.nih.gov/cellminer/loadDownload.do]
all_protein_data <- read_excel("input/Protein__SWATH_(Mass_spectrometry)_Protein.xls", skip = 10)

names_changed_to_match_format_protsheet <- str_remove(str_remove(str_replace(str_replace(str_remove(colnames(all_protein_data)[10:ncol(all_protein_data)], "^[A-Z]{2,3}:"), "/", "_"), " ", "_"), "\\(TB\\)"), "_ATCC")
colnames(all_protein_data)[10:ncol(all_protein_data)] <- names_changed_to_match_format_protsheet

protdatat <- t(all_protein_data)
colnames(protdatat) <- protdatat[2, ]

# random_detected_protein_sample <- sample(all_protein_data$`Gene name d`, 1000, replace = FALSE)
# saveRDS(random_detected_protein_sample, "output/random_detected_protein_sample.rds" )
random_detected_protein_sample <- readRDS("output/random_detected_protein_sample.rds" )

for(i in 1:length(random_detected_protein_sample)){
  
  H4K20promotersmatrix_for_lm[, paste0(random_detected_protein_sample[i], "_protein")] <- as.numeric(protdatat[match(str_remove(row.names(H4K20promotersmatrix_for_lm), "H4K20me3_[A-Z]{2,3}_"), row.names(protdatat)), random_detected_protein_sample[i]])
  
}

proteinrandomsample_H4K20me3promoters_t_values <- sapply(random_detected_protein_sample, function(thisgene){
  
  sapply(colnames(H4K20promotersmatrix_for_lm)[1:H4K20promotersnumber_of_sites], function(x){
    
    tryCatch(expr = {
      glmtemp <- glm(formula = as.formula(paste0(x, " ~ ", as.character(thisgene), "_protein + tissue + ", thisgene, "_protein * tissue")), data = H4K20promotersmatrix_for_lm)
      glmsum <- summary(glmtemp)
      
      t <- glmsum$coefficients[paste0(thisgene, "_protein"), "t value"]
      
      output_vec <- c(t)
      names(output_vec) <- c("t_value")
      
      return(output_vec)
      
    }, error = function(e){
      
      output_vec <- c(NA)
      names(output_vec) <- c("t_value")
      
      return(output_vec)
      
    }
    )
    
  })
  
})

saveRDS(proteinrandomsample_H4K20me3promoters_t_values, "output/proteinrandomsample_H4K20me3promoters_t_values.rds")
# proteinrandomsample_H4K20me3promoters_t_values <- readRDS("output/proteinrandomsample_H4K20me3promoters_t_values.rds")

#### RANDOM DISTRIBUTION OF H4K20 GENE BODIES PROTEIN ####

# random_detected_protein_sample <- sample(all_protein_data$`Gene name d`, 1000, replace = FALSE)
# saveRDS(random_detected_protein_sample, "output/random_detected_protein_sample.rds" )
random_detected_protein_sample <- readRDS("output/random_detected_protein_sample.rds" )

for(i in 1:length(random_detected_protein_sample)){
  
  H4K20genesmatrix_for_lm[, paste0(random_detected_protein_sample[i], "_protein")] <- as.numeric(protdatat[match(str_remove(row.names(H4K20genesmatrix_for_lm), "H4K20me3_[A-Z]{2,3}_"), row.names(protdatat)), random_detected_protein_sample[i]])
  
}

proteinrandomsample_H4K20me3genes_t_values <- sapply(random_detected_protein_sample, function(thisgene){
  
  sapply(colnames(H4K20genesmatrix_for_lm)[1:H4K20genesnumber_of_sites], function(x){
    
    tryCatch(expr = {
      glmtemp <- glm(formula = as.formula(paste0(x, " ~ ", as.character(thisgene), "_protein + tissue + ", thisgene, "_protein * tissue")), data = H4K20genesmatrix_for_lm)
      glmsum <- summary(glmtemp)
      
      t <- glmsum$coefficients[paste0(thisgene, "_protein"), "t value"]
      
      output_vec <- c(t)
      names(output_vec) <- c("t_value")
      
      return(output_vec)
      
    }, error = function(e){
      
      output_vec <- c(NA)
      names(output_vec) <- c("t_value")
      
      return(output_vec)
      
    }
    )
    
  })
  
})

saveRDS(proteinrandomsample_H4K20me3genes_t_values, "output/proteinrandomsample_H4K20me3genes_t_values.rds")
# proteinrandomsample_H4K20me3genes_t_values <- readRDS("output/proteinrandomsample_H4K20me3genes_t_values.rds")

#### MAKE BOXPLOTS ####

#### check signal output in sites ####

NCI60_H3K4_promoters_signal_sites <- names(H3K4promotersmatrix_for_lm[colMeans(H3K4promotersmatrix_for_lm[, 1:H3K4promotersnumber_of_sites]) > 0])
NCI60_H3K4_promoters_signal_sites <- NCI60_H3K4_promoter_signal_sites[!(NCI60_H3K4_promoters_signal_sites %in% c("tissue")) & !str_detect(NCI60_H3K4_promoters_signal_sites, "MT")]

NCI60_H3K4_genes_signal_sites <- names(H3K4genesmatrix_for_lm[colMeans(H3K4genesmatrix_for_lm[, 1:H3K4genesnumber_of_sites]) > 0])
NCI60_H3K4_genes_signal_sites <- NCI60_H3K4_genes_signal_sites[!(NCI60_H3K4_genes_signal_sites %in% c("tissue")) & !str_detect(NCI60_H3K4_genes_signal_sites, "MT")]


NCI60_H3K9_promoters_signal_sites <- names(H3K9promotersmatrix_for_lm[colMeans(H3K9promotersmatrix_for_lm[, 1:H3K9promotersnumber_of_sites]) > 0])
NCI60_H3K9_promoters_signal_sites <- NCI60_H3K9_promoters_signal_sites[!(NCI60_H3K9_promoters_signal_sites %in% c("tissue")) & !str_detect(NCI60_H3K9_promoters_signal_sites, "MT")]

NCI60_H3K9_genes_signal_sites <- names(H3K9genesmatrix_for_lm[colMeans(H3K9genesmatrix_for_lm[, 1:H3K9genesnumber_of_sites]) > 0])
NCI60_H3K9_genes_signal_sites <- NCI60_H3K9_genes_signal_sites[!(NCI60_H3K9_genes_signal_sites %in% c("tissue")) & !str_detect(NCI60_H3K9_genes_signal_sites, "MT")]

NCI60_H4K20_promoters_signal_sites <- names(H4K20promotersmatrix_for_lm[colMeans(H4K20promotersmatrix_for_lm[, 1:H4K20promotersnumber_of_sites]) > 0])
NCI60_H4K20_promoters_signal_sites <- NCI60_H4K20_promoters_signal_sites[!(NCI60_H4K20_promoters_signal_sites %in% c("tissue")) & !str_detect(NCI60_H4K20_promoters_signal_sites, "MT")]

NCI60_H4K20_genes_signal_sites <- names(H4K20genesmatrix_for_lm[colMeans(H4K20genesmatrix_for_lm[, 1:H4K20genesnumber_of_sites]) > 0])
NCI60_H4K20_genes_signal_sites <- NCI60_H4K20_genes_signal_sites[!(NCI60_H4K20_genes_signal_sites %in% c("tissue")) & !str_detect(NCI60_H4K20_genes_signal_sites, "MT")]

#### RNA-SEQ BOXPLOT ####

# load RNA-seq models

H4K20NNMTRNAseq_genes_signal_model <- readRDS("output/NCI60_H4K20me3_genes_signal_model_NNMT_RNAseq.rds")
H3K9NNMTRNAseq_genes_signal_model <- readRDS("output/NCI60_H3K9me3_genes_signal_model_NNMT_RNAseq.rds")
H3K4NNMTRNAseq_genes_signal_model <- readRDS("output/NCI60_H3K4me3_genes_signal_model_NNMT_RNAseq.rds")

H4K20NNMTRNAseq_promoters_signal_model <- readRDS("output/NCI60_H4K20me3_promoters_signal_model_NNMT_RNAseq.rds")
H3K9NNMTRNAseq_promoters_signal_model <- readRDS("output/NCI60_H3K9me3_promoters_signal_model_NNMT_RNAseq.rds")
H3K4NNMTRNAseq_promoters_signal_model <- readRDS("output/NCI60_H3K4me3_promoters_signal_model_NNMT_RNAseq.rds")

H4K20random_genes_signal_model <- readRDS("output/RNAseqrandomsample_H4K20me3genes_t_values.rds")
H3K9random_genes_signal_model <- readRDS("output/RNAseqrandomsample_H3K9me3genes_t_values.rds")
H3K4random_genes_signal_model <- readRDS("output/RNAseqrandomsample_H3K4me3genes_t_values.rds")

H4K20random_promoters_signal_model <- readRDS("output/RNAseqrandomsample_H4K20me3promoters_t_values.rds")
H3K9random_promoters_signal_model <- readRDS("output/RNAseqrandomsample_H3K9me3promoters_t_values.rds")
H3K4random_promoters_signal_model <- readRDS("output/RNAseqrandomsample_H3K4me3promoters_t_values.rds")

NCI60_RNAseq_plot_list <- list(genes = list(H3K4_NNMT = unlist(H3K4NNMTRNAseq_genes_signal_model[NCI60_H3K4_genes_signal_sites]),
                                            H3K4_random = rowMeans(H3K4random_genes_signal_model[str_remove(row.names(H3K4random_genes_signal_model), "\\.t_value") %in% NCI60_H3K4_genes_signal_sites, ], na.rm = TRUE),
                                            H4K20_NNMT = unlist(H4K20NNMTRNAseq_genes_signal_model[NCI60_H4K20_genes_signal_sites]),
                                            H4K20_random = rowMeans(H4K20random_genes_signal_model[str_remove(row.names(H4K20random_genes_signal_model), "\\.t_value") %in% NCI60_H4K20_genes_signal_sites, ], na.rm = TRUE),
                                            H3K9_NNMT = unlist(H3K9NNMTRNAseq_genes_signal_model[NCI60_H3K9_genes_signal_sites]),
                                            H3K9_random = rowMeans(H3K9random_genes_signal_model[str_remove(row.names(H3K9random_genes_signal_model), "\\.t_value") %in% NCI60_H3K9_genes_signal_sites, ], na.rm = TRUE)),
                               promoters = list(H3K4_NNMT = unlist(H3K4NNMTRNAseq_promoters_signal_model[NCI60_H3K4_promoters_signal_sites]),
                                                H3K4_random = rowMeans(H3K4random_promoters_signal_model[str_remove(row.names(H3K4random_promoters_signal_model), "\\.t_value") %in% NCI60_H3K4_promoters_signal_sites, ], na.rm = TRUE),
                                                H4K20_NNMT = unlist(H4K20NNMTRNAseq_promoters_signal_model[NCI60_H4K20_promoters_signal_sites]),
                                                H4K20_random = rowMeans(H4K20random_promoters_signal_model[str_remove(row.names(H4K20random_promoters_signal_model), "\\.t_value") %in% NCI60_H4K20_promoters_signal_sites, ], na.rm = TRUE),
                                                H3K9_NNMT = unlist(H3K9NNMTRNAseq_promoters_signal_model[NCI60_H3K9_promoters_signal_sites]),
                                                H3K9_random = rowMeans(H3K9random_promoters_signal_model[str_remove(row.names(H3K9random_promoters_signal_model), "\\.t_value") %in% NCI60_H3K9_promoters_signal_sites, ], na.rm = TRUE)))           

RNAseq_plot_melt <- melt(NCI60_RNAseq_plot_list)
colnames(RNAseq_plot_melt) <- c("signal", "group", "location")

write.table(RNAseq_plot_melt,
            file = "plot_data/Fig 3/Fig_3C_plotdata.txt",
            col.names = TRUE,
            row.names = FALSE,
            sep = "\t")

factet_labs <- c("Gene bodies", "Promoters")
names(factet_labs) <- c("genes", "promoters")

samplesizes <- table(paste0(RNAseq_plot_melt$group, RNAseq_plot_melt$location))
samplesizes <- as.data.frame(samplesizes)
samplesizes[, "group"] <- str_remove(str_remove(samplesizes$Var1, pattern = "promoters"), pattern = "genes")

for(i in 1:nrow(samplesizes)){
  samplesizes[i, "location"] <- str_remove(samplesizes[i, "Var1"], pattern = samplesizes[i, "group"])
}

samplesizes$Freq = paste0("n = ", samplesizes$Freq)

samplesizes[!str_detect(samplesizes$Var1, pattern = "random"), "Freq"] <- ""

samplesizes[samplesizes$Var1 == "H4K20_NNMTgenes", "p_value"] <- signif(wilcox.test(NCI60_RNAseq_plot_list$genes$H4K20_random, unlist(NCI60_RNAseq_plot_list$genes$H4K20_NNMT), paired = TRUE)$p.value, digits = 2)
samplesizes[samplesizes$Var1 == "H4K20_randomgenes", "p_value_parse"] <- paste0(str_remove(samplesizes[samplesizes$Var1 == "H4K20_NNMTgenes", "p_value"], pattern = "e.*"), " * 10^", str_remove(samplesizes[samplesizes$Var1 == "H4K20_NNMTgenes", "p_value"], pattern = "p = .*e"))

samplesizes[samplesizes$Var1 == "H4K20_NNMTpromoters", "p_value"] <- signif(wilcox.test(NCI60_RNAseq_plot_list$promoters$H4K20_random, unlist(NCI60_RNAseq_plot_list$promoters$H4K20_NNMT), paired = TRUE)$p.value, digits = 2)
samplesizes[samplesizes$Var1 == "H4K20_randompromoters", "p_value_parse"] <- paste0(str_remove(samplesizes[samplesizes$Var1 == "H4K20_NNMTpromoters", "p_value"], pattern = "e.*"), " * 10^", str_remove(samplesizes[samplesizes$Var1 == "H4K20_NNMTpromoters", "p_value"], pattern = "p = .*e"))

samplesizes[samplesizes$Var1 == "H3K9_NNMTgenes", "p_value"] <- signif(wilcox.test(NCI60_RNAseq_plot_list$genes$H3K9_random, unlist(NCI60_RNAseq_plot_list$genes$H3K9_NNMT), paired = TRUE)$p.value, digits = 2)
samplesizes[samplesizes$Var1 == "H3K9_randomgenes", "p_value_parse"] <- paste0(str_remove(samplesizes[samplesizes$Var1 == "H3K9_NNMTgenes", "p_value"], pattern = "e.*"), " * 10^", str_remove(samplesizes[samplesizes$Var1 == "H3K9_NNMTgenes", "p_value"], pattern = "p = .*e"))

samplesizes[samplesizes$Var1 == "H3K9_NNMTpromoters", "p_value"] <- signif(wilcox.test(NCI60_RNAseq_plot_list$promoters$H3K9_random, unlist(NCI60_RNAseq_plot_list$promoters$H3K9_NNMT), paired = TRUE)$p.value, digits = 2)
samplesizes[samplesizes$Var1 == "H3K9_randompromoters", "p_value_parse"] <- paste0(str_remove(samplesizes[samplesizes$Var1 == "H3K9_NNMTpromoters", "p_value"], pattern = "e.*"), " * 10^", str_remove(samplesizes[samplesizes$Var1 == "H3K9_NNMTpromoters", "p_value"], pattern = "p = .*e"))

samplesizes[samplesizes$Var1 == "H3K4_NNMTgenes", "p_value"] <- signif(wilcox.test(NCI60_RNAseq_plot_list$genes$H3K4_random, unlist(NCI60_RNAseq_plot_list$genes$H3K4_NNMT), paired = TRUE)$p.value, digits = 2)
samplesizes[samplesizes$Var1 == "H3K4_randomgenes", "p_value_parse"] <- paste0(str_remove(samplesizes[samplesizes$Var1 == "H3K4_NNMTgenes", "p_value"], pattern = "e.*"), " * 10^", str_remove(samplesizes[samplesizes$Var1 == "H3K4_NNMTgenes", "p_value"], pattern = "p = .*e"))

samplesizes[samplesizes$Var1 == "H3K4_NNMTpromoters", "p_value"] <- signif(wilcox.test(NCI60_RNAseq_plot_list$promoters$H3K4_random, unlist(NCI60_RNAseq_plot_list$promoters$H3K4_NNMT), paired = TRUE)$p.value, digits = 2)
samplesizes[samplesizes$Var1 == "H3K4_randompromoters", "p_value_parse"] <- paste0(str_remove(samplesizes[samplesizes$Var1 == "H3K4_NNMTpromoters", "p_value"], pattern = "e.*"), " * 10^", str_remove(samplesizes[samplesizes$Var1 == "H3K4_NNMTpromoters", "p_value"], pattern = "p = .*e"))

samplesizes[is.na(samplesizes$p_value), "p_value"] <- ""

dat_text <- data.frame(
  label = samplesizes[samplesizes$p_value > -0.001, "p_value"],
  location   = samplesizes[samplesizes$p_value > -0.001, "location"],
  group = samplesizes[samplesizes$p_value > -0.001, "group"]
)

RNAseq_plot_melt$group <- factor(RNAseq_plot_melt$group, levels = c("H4K20_random", "H4K20_NNMT", "H3K9_random", "H3K9_NNMT", "H3K4_random", "H3K4_NNMT"))

pdf("graphics/NCI60_NNMT_RNAseq_nopval.pdf",
    width = 3.75,
    height = 5)

ggplot(aes(x = signal, y = group, group = location), data = RNAseq_plot_melt ) + 
  geom_boxplot(aes(group = group), 
               fill = c("grey",
                        "darkolivegreen3",
                        "grey",
                        "purple",
                        "grey",
                        "magenta",
                        "grey",
                        "darkolivegreen3", 
                        "grey",
                        "purple",
                        "grey",
                        "magenta"), 
               notch = TRUE,
               outlier.shape = NA,
               lwd = 0.3,
               fatten = 0.7) +
  theme_classic() +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey") +
  xlab("Linear model t-value") + 
  ylab(NULL) +
  # scale_y_discrete(label = c("H4K20me3 (random genes)",
  #                            "H4K20me3 (NNMT RNA-seq)",
  #                            "H3K9me3 (random genes)",
  #                            "H3K9me3 (NNMT RNA-seq)",
  #                            "H3K4me3 (random genes)",
  #                            "H3K4me3 (NNMT RNA-seq)")) +
  facet_wrap(~ location, ncol = 1, labeller = labeller(location = factet_labs), scales = "free") + 
  coord_cartesian(xlim = c(-4, 5.5)) + 
  theme(strip.text = element_text(size = 12, color = "white", face = "bold"),
        strip.background = element_rect(fill = "black"),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.title.x = element_text(size = 10, color = "black", face = "bold")) +
  geom_text(aes(x = 3.8, group = location, label = Freq), data = samplesizes, size = 3, hjust = 0)

dev.off()


pdf("graphics/NCI60_NNMT_RNAseq_pval.pdf",
    width = 3.75,
    height = 5)

ggplot(aes(x = signal, y = group, group = location), data = RNAseq_plot_melt ) + 
  geom_boxplot(aes(group = group), 
               fill = c("grey",
                        "darkolivegreen3",
                        "grey",
                        "purple",
                        "grey",
                        "magenta",
                        "grey",
                        "darkolivegreen3", 
                        "grey",
                        "purple",
                        "grey",
                        "magenta"), 
               notch = TRUE,
               outlier.shape = NA,
               lwd = 0.3,
               fatten = 0.7) +
  theme_classic() +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey") +
  xlab("Linear model t-value") + 
  ylab(NULL) +
  scale_y_discrete(label = c("H4K20me3 (random genes)",
                             "H4K20me3 (NNMT RNA-seq)",
                             "H3K9me3 (random genes)",
                             "H3K9me3 (NNMT RNA-seq)",
                             "H3K4me3 (random genes)",
                             "H3K4me3 (NNMT RNA-seq)")) +
  facet_wrap(~ location, ncol = 1, labeller = labeller(location = factet_labs), scales = "free") +
  coord_cartesian(xlim = c(-4, 5.5)) + 
  theme(strip.text = element_text(size = 12, color = "white", face = "bold"),
        strip.background = element_rect(fill = "black"),
        # axis.text.y = element_blank(),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.title.x = element_text(size = 10, color = "black", face = "bold")) +
  geom_text(aes(x = 3.8, group = location, label = Freq), data = samplesizes, size = 3, hjust = 0) +
  geom_text(aes(x = 3.8, group = location, label = label), data = dat_text, parse = TRUE, size =  3)

dev.off()

#### mRNA BOXPLOT ####

# load mRNA models

H4K20NNMTmRNA_genes_signal_model <- readRDS("output/NCI60_H4K20me3_genes_signal_model_NNMTmRNA.rds")
H3K9NNMTmRNA_genes_signal_model <- readRDS("output/NCI60_H3K9me3_genes_signal_model_NNMTmRNA.rds")
H3K4NNMTmRNA_genes_signal_model <- readRDS("output/NCI60_H3K4me3_genes_signal_model_NNMTmRNA.rds")

H4K20NNMTmRNA_promoters_signal_model <- readRDS("output/NCI60_H4K20me3_promoters_signal_model_NNMTmRNA.rds")
H3K9NNMTmRNA_promoters_signal_model <- readRDS("output/NCI60_H3K9me3_promoters_signal_model_NNMTmRNA.rds")
H3K4NNMTmRNA_promoters_signal_model <- readRDS("output/NCI60_H3K4me3_promoters_signal_model_NNMTmRNA.rds")

H4K20random_genes_signal_model <- readRDS("output/mRNArandomsample_H4K20me3genes_t_values.rds")
H3K9random_genes_signal_model <- readRDS("output/mRNArandomsample_H3K9me3genes_t_values.rds")
H3K4random_genes_signal_model <- readRDS("output/mRNArandomsample_H3K4me3genes_t_values.rds")

H4K20random_promoters_signal_model <- readRDS("output/mRNArandomsample_H4K20promoters_t_values.rds")
H3K9random_promoters_signal_model <- readRDS("output/mRNArandomsample_H3K9me3promoters_t_values.rds")
H3K4random_promoters_signal_model <- readRDS("output/mRNArandomsample_H3K4me3promoters_t_values.rds")

NCI60mRNA_plot_list <- list(genes = list(H3K4_NNMT = unlist(H3K4NNMTmRNA_genes_signal_model[NCI60_H3K4_genes_signal_sites]),
                                         H3K4_random = rowMeans(H3K4random_genes_signal_model[str_remove(row.names(H3K4random_genes_signal_model), "\\.t_value") %in% NCI60_H3K4_genes_signal_sites, ], na.rm = TRUE),
                                         H4K20_NNMT = unlist(H4K20NNMTmRNA_genes_signal_model[NCI60_H4K20_genes_signal_sites]),
                                         H4K20_random = rowMeans(H4K20random_genes_signal_model[str_remove(row.names(H4K20random_genes_signal_model), "\\.t_value") %in% NCI60_H4K20_genes_signal_sites, ], na.rm = TRUE),
                                         H3K9_NNMT = unlist(H3K9NNMTmRNA_genes_signal_model[NCI60_H3K9_genes_signal_sites]),
                                         H3K9_random = rowMeans(H3K9random_genes_signal_model[str_remove(row.names(H3K9random_genes_signal_model), "\\.t_value") %in% NCI60_H3K9_genes_signal_sites, ], na.rm = TRUE)),
                            promoters = list(H3K4_NNMT = unlist(H3K4NNMTmRNA_promoters_signal_model[NCI60_H3K4_promoters_signal_sites]),
                                             H3K4_random = rowMeans(H3K4random_promoters_signal_model[str_remove(row.names(H3K4random_promoters_signal_model), "\\.t_value") %in% NCI60_H3K4_promoters_signal_sites, ], na.rm = TRUE),
                                             H4K20_NNMT = unlist(H4K20NNMTmRNA_promoters_signal_model[NCI60_H4K20_promoters_signal_sites]),
                                             H4K20_random = rowMeans(H4K20random_promoters_signal_model[str_remove(row.names(H4K20random_promoters_signal_model), "\\.t_value") %in% NCI60_H4K20_promoters_signal_sites, ], na.rm = TRUE),
                                             H3K9_NNMT = unlist(H3K9NNMTmRNA_promoters_signal_model[NCI60_H3K9_promoters_signal_sites]),
                                             H3K9_random = rowMeans(H3K9random_promoters_signal_model[str_remove(row.names(H3K9random_promoters_signal_model), "\\.t_value") %in% NCI60_H3K9_promoters_signal_sites, ], na.rm = TRUE)))           


mRNA_plot_melt <- melt(NCI60mRNA_plot_list)
colnames(mRNA_plot_melt) <- c("signal", "group", "location")

write.table(mRNA_plot_melt,
            file = "plot_data/Fig S9/Fig_S9A_plotdata.txt",
            col.names = TRUE,
            row.names = FALSE,
            sep = "\t")

factet_labs <- c("Gene bodies", "Promoters")
names(factet_labs) <- c("genes", "promoters")

samplesizes <- table(paste0(mRNA_plot_melt$group, mRNA_plot_melt$location))
samplesizes <- as.data.frame(samplesizes)
samplesizes[, "group"] <- str_remove(str_remove(samplesizes$Var1, pattern = "promoters"), pattern = "genes")

for(i in 1:nrow(samplesizes)){
  samplesizes[i, "location"] <- str_remove(samplesizes[i, "Var1"], pattern = samplesizes[i, "group"])
}

samplesizes$Freq = paste0("n = ", samplesizes$Freq)

samplesizes[!str_detect(samplesizes$Var1, pattern = "random"), "Freq"] <- ""

samplesizes[samplesizes$Var1 == "H4K20_NNMTgenes", "p_value"] <- signif(wilcox.test(NCI60mRNA_plot_list$genes$H4K20_random, unlist(NCI60mRNA_plot_list$genes$H4K20_NNMT), paired = TRUE)$p.value, digits = 2)
samplesizes[samplesizes$Var1 == "H4K20_randomgenes", "p_value_parse"] <- paste0(str_remove(samplesizes[samplesizes$Var1 == "H4K20_NNMTgenes", "p_value"], pattern = "e.*"), " * 10^", str_remove(samplesizes[samplesizes$Var1 == "H4K20_NNMTgenes", "p_value"], pattern = "p = .*e"))

samplesizes[samplesizes$Var1 == "H4K20_NNMTpromoters", "p_value"] <- signif(wilcox.test(NCI60mRNA_plot_list$promoters$H4K20_random, unlist(NCI60mRNA_plot_list$promoters$H4K20_NNMT), paired = TRUE)$p.value, digits = 2)
samplesizes[samplesizes$Var1 == "H4K20_randompromoters", "p_value_parse"] <- paste0(str_remove(samplesizes[samplesizes$Var1 == "H4K20_NNMTpromoters", "p_value"], pattern = "e.*"), " * 10^", str_remove(samplesizes[samplesizes$Var1 == "H4K20_NNMTpromoters", "p_value"], pattern = "p = .*e"))

samplesizes[samplesizes$Var1 == "H3K9_NNMTgenes", "p_value"] <- signif(wilcox.test(NCI60mRNA_plot_list$genes$H3K9_random, unlist(NCI60mRNA_plot_list$genes$H3K9_NNMT), paired = TRUE)$p.value, digits = 2)
samplesizes[samplesizes$Var1 == "H3K9_randomgenes", "p_value_parse"] <- paste0(str_remove(samplesizes[samplesizes$Var1 == "H3K9_NNMTgenes", "p_value"], pattern = "e.*"), " * 10^", str_remove(samplesizes[samplesizes$Var1 == "H3K9_NNMTgenes", "p_value"], pattern = "p = .*e"))

samplesizes[samplesizes$Var1 == "H3K9_NNMTpromoters", "p_value"] <- signif(wilcox.test(NCI60mRNA_plot_list$promoters$H3K9_random, unlist(NCI60mRNA_plot_list$promoters$H3K9_NNMT), paired = TRUE)$p.value, digits = 2)
samplesizes[samplesizes$Var1 == "H3K9_randompromoters", "p_value_parse"] <- paste0(str_remove(samplesizes[samplesizes$Var1 == "H3K9_NNMTpromoters", "p_value"], pattern = "e.*"), " * 10^", str_remove(samplesizes[samplesizes$Var1 == "H3K9_NNMTpromoters", "p_value"], pattern = "p = .*e"))

samplesizes[samplesizes$Var1 == "H3K4_NNMTgenes", "p_value"] <- signif(wilcox.test(NCI60mRNA_plot_list$genes$H3K4_random, unlist(NCI60mRNA_plot_list$genes$H3K4_NNMT), paired = TRUE)$p.value, digits = 2)
samplesizes[samplesizes$Var1 == "H3K4_randomgenes", "p_value_parse"] <- paste0(str_remove(samplesizes[samplesizes$Var1 == "H3K4_NNMTgenes", "p_value"], pattern = "e.*"), " * 10^", str_remove(samplesizes[samplesizes$Var1 == "H3K4_NNMTgenes", "p_value"], pattern = "p = .*e"))

samplesizes[samplesizes$Var1 == "H3K4_NNMTpromoters", "p_value"] <- signif(wilcox.test(NCI60mRNA_plot_list$promoters$H3K4_random, unlist(NCI60mRNA_plot_list$promoters$H3K4_NNMT), paired = TRUE)$p.value, digits = 2)
samplesizes[samplesizes$Var1 == "H3K4_randompromoters", "p_value_parse"] <- paste0(str_remove(samplesizes[samplesizes$Var1 == "H3K4_NNMTpromoters", "p_value"], pattern = "e.*"), " * 10^", str_remove(samplesizes[samplesizes$Var1 == "H3K4_NNMTpromoters", "p_value"], pattern = "p = .*e"))

samplesizes[is.na(samplesizes$p_value), "p_value"] <- ""

dat_text <- data.frame(
  label = samplesizes[samplesizes$p_value > -0.001, "p_value"],
  location   = samplesizes[samplesizes$p_value > -0.001, "location"],
  group = samplesizes[samplesizes$p_value > -0.001, "group"]
)


mRNA_plot_melt$group <- factor(mRNA_plot_melt$group, levels = c("H4K20_random", "H4K20_NNMT", "H3K9_random", "H3K9_NNMT", "H3K4_random", "H3K4_NNMT"))

pdf("graphics/NCI60_NNMTmRNA_nopval.pdf",
    width = 2.25,
    height = 5)

ggplot(aes(x = signal, y = group, group = location), data = mRNA_plot_melt ) + 
  geom_boxplot(aes(group = group), 
               fill = c("grey",
                        "darkolivegreen3",
                        "grey",
                        "purple",
                        "grey",
                        "magenta",
                        "grey",
                        "darkolivegreen3", 
                        "grey",
                        "purple",
                        "grey",
                        "magenta"), 
               notch = TRUE,
               outlier.shape = NA,
               lwd = 0.3,
               fatten = 0.7) +
  theme_classic() +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey") +
  xlab("Linear model t-value") + 
  ylab(NULL) +
  # scale_y_discrete(label = c("H4K20me3 (random genes)",
  #                            "H4K20me3 (NNMT RNA-seq)",
  #                            "H3K9me3 (random genes)",
  #                            "H3K9me3 (NNMT RNA-seq)",
  #                            "H3K4me3 (random genes)",
  #                            "H3K4me3 (NNMT RNA-seq)")) +
  facet_wrap(~ location, ncol = 1, labeller = labeller(location = factet_labs), scales = "free") + 
  coord_cartesian(xlim = c(-4, 4)) + 
  theme(strip.text = element_text(size = 12, color = "white", face = "bold"),
        strip.background = element_rect(fill = "black"),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.title.x = element_text(size = 10, color = "black", face = "bold")) +
  geom_text(aes(x = 2, group = location, label = Freq), data = samplesizes, size = 3, hjust = 0)

dev.off()

pdf("graphics/NCI60_NNMTmRNA_pval.pdf",
    width = 3.75,
    height = 5)

ggplot(aes(x = signal, y = group, group = location), data = mRNA_plot_melt ) + 
  geom_boxplot(aes(group = group), 
               fill = c("grey",
                        "darkolivegreen3",
                        "grey",
                        "purple",
                        "grey",
                        "magenta",
                        "grey",
                        "darkolivegreen3", 
                        "grey",
                        "purple",
                        "grey",
                        "magenta"), 
               notch = TRUE,
               outlier.shape = NA,
               lwd = 0.3,
               fatten = 0.7) +
  theme_classic() +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey") +
  xlab("Linear model t-value") + 
  ylab(NULL) +
  scale_y_discrete(label = c("H4K20me3 (random genes)",
                             "H4K20me3 (NNMT RNA-seq)",
                             "H3K9me3 (random genes)",
                             "H3K9me3 (NNMT RNA-seq)",
                             "H3K4me3 (random genes)",
                             "H3K4me3 (NNMT RNA-seq)")) +
  facet_wrap(~ location, ncol = 1, labeller = labeller(location = factet_labs), scales = "free") + 
  coord_cartesian(xlim = c(-4, 5.5)) + 
  theme(strip.text = element_text(size = 12, color = "white", face = "bold"),
        strip.background = element_rect(fill = "black"),
        # axis.text.y = element_blank(),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.title.x = element_text(size = 10, color = "black", face = "bold")) +
  geom_text(aes(x = 3.8, group = location, label = Freq), data = samplesizes, size = 3, hjust = 0) +
  geom_text(aes(x = 3.8, group = location, label = label), data = dat_text, parse = TRUE, size =  3)

dev.off()

#### PROTEIN BOXPLOT ####

# load protein models

H4K20NNMTprotein_genes_signal_model <- readRDS("output/NCI60_H4K20me3_genes_signal_model_NNMTprotein.rds")
H3K9NNMTprotein_genes_signal_model <- readRDS("output/NCI60_H3K9me3_genes_signal_model_NNMTprotein.rds")
H3K4NNMTprotein_genes_signal_model <- readRDS("output/NCI60_H3K4me3_genes_signal_model_NNMTprotein.rds")

H4K20NNMTprotein_promoters_signal_model <- readRDS("output/NCI60_H4K20me3_promoters_signal_model_NNMTprotein.rds")
H3K9NNMTprotein_promoters_signal_model <- readRDS("output/NCI60_H3K9me3_promoters_signal_model_NNMTprotein.rds")
H3K4NNMTprotein_promoters_signal_model <- readRDS("output/NCI60_H3K4me3_promoters_signal_model_NNMTprotein.rds")

H4K20random_genes_signal_model <- readRDS("output/proteinrandomsample_H4K20me3genes_t_values.rds")
H3K9random_genes_signal_model <- readRDS("output/proteinrandomsample_H3K9me3genes_t_values.rds")
H3K4random_genes_signal_model <- readRDS("output/proteinrandomsample_H3K4me3genes_t_values.rds")

H4K20random_promoters_signal_model <- readRDS("output/proteinrandomsample_H4K20me3promoters_t_values.rds")
H3K9random_promoters_signal_model <- readRDS("output/proteinrandomsample_H3K9me3promoters_t_values.rds")
H3K4random_promoters_signal_model <- readRDS("output/proteinrandomsample_H3K4me3promoters_t_values.rds")

NCI60protein_plot_list <- list(genes = list(H3K4_NNMT = unlist(H3K4NNMTprotein_genes_signal_model[NCI60_H3K4_genes_signal_sites]),
                                            H3K4_random = rowMeans(H3K4random_genes_signal_model[str_remove(row.names(H3K4random_genes_signal_model), "\\.t_value") %in% NCI60_H3K4_genes_signal_sites, ], na.rm = TRUE),
                                            H4K20_NNMT = unlist(H4K20NNMTprotein_genes_signal_model[NCI60_H4K20_genes_signal_sites]),
                                            H4K20_random = rowMeans(H4K20random_genes_signal_model[str_remove(row.names(H4K20random_genes_signal_model), "\\.t_value") %in% NCI60_H4K20_genes_signal_sites, ], na.rm = TRUE),
                                            H3K9_NNMT = unlist(H3K9NNMTprotein_genes_signal_model[NCI60_H3K9_genes_signal_sites]),
                                            H3K9_random = rowMeans(H3K9random_genes_signal_model[str_remove(row.names(H3K9random_genes_signal_model), "\\.t_value") %in% NCI60_H3K9_genes_signal_sites, ], na.rm = TRUE)),
                               promoters = list(H3K4_NNMT = unlist(H3K4NNMTprotein_promoters_signal_model[NCI60_H3K4_promoters_signal_sites]),
                                                H3K4_random = rowMeans(H3K4random_promoters_signal_model[str_remove(row.names(H3K4random_promoters_signal_model), "\\.t_value") %in% NCI60_H3K4_promoters_signal_sites, ], na.rm = TRUE),
                                                H4K20_NNMT = unlist(H4K20NNMTprotein_promoters_signal_model[NCI60_H4K20_promoters_signal_sites]),
                                                H4K20_random = rowMeans(H4K20random_promoters_signal_model[str_remove(row.names(H4K20random_promoters_signal_model), "\\.t_value") %in% NCI60_H4K20_promoters_signal_sites, ], na.rm = TRUE),
                                                H3K9_NNMT = unlist(H3K9NNMTprotein_promoters_signal_model[NCI60_H3K9_promoters_signal_sites]),
                                                H3K9_random = rowMeans(H3K9random_promoters_signal_model[str_remove(row.names(H3K9random_promoters_signal_model), "\\.t_value") %in% NCI60_H3K9_promoters_signal_sites, ], na.rm = TRUE)))           

protein_plot_melt <- melt(NCI60protein_plot_list)
colnames(protein_plot_melt) <- c("signal", "group", "location")

write.table(protein_plot_melt,
            file = "plot_data/Fig S9/Fig_S9B_plotdata.txt",
            row.names = FALSE,
            col.names = TRUE,
            sep = "\t")

factet_labs <- c("Gene bodies", "Promoters")
names(factet_labs) <- c("genes", "promoters")

samplesizes <- table(paste0(protein_plot_melt$group, protein_plot_melt$location))
samplesizes <- as.data.frame(samplesizes)
samplesizes[, "group"] <- str_remove(str_remove(samplesizes$Var1, pattern = "promoters"), pattern = "genes")

for(i in 1:nrow(samplesizes)){
  samplesizes[i, "location"] <- str_remove(samplesizes[i, "Var1"], pattern = samplesizes[i, "group"])
}

samplesizes$Freq = paste0("n = ", samplesizes$Freq)

samplesizes[!str_detect(samplesizes$Var1, pattern = "random"), "Freq"] <- ""

samplesizes[samplesizes$Var1 == "H4K20_NNMTgenes", "p_value"] <- signif(wilcox.test(NCI60protein_plot_list$genes$H4K20_random, unlist(NCI60protein_plot_list$genes$H4K20_NNMT), paired = TRUE)$p.value, digits = 2)
samplesizes[samplesizes$Var1 == "H4K20_randomgenes", "p_value_parse"] <- paste0(str_remove(samplesizes[samplesizes$Var1 == "H4K20_NNMTgenes", "p_value"], pattern = "e.*"), " * 10^", str_remove(samplesizes[samplesizes$Var1 == "H4K20_NNMTgenes", "p_value"], pattern = "p = .*e"))

samplesizes[samplesizes$Var1 == "H4K20_NNMTpromoters", "p_value"] <- signif(wilcox.test(NCI60protein_plot_list$promoters$H4K20_random, unlist(NCI60protein_plot_list$promoters$H4K20_NNMT), paired = TRUE)$p.value, digits = 2)
samplesizes[samplesizes$Var1 == "H4K20_randompromoters", "p_value_parse"] <- paste0(str_remove(samplesizes[samplesizes$Var1 == "H4K20_NNMTpromoters", "p_value"], pattern = "e.*"), " * 10^", str_remove(samplesizes[samplesizes$Var1 == "H4K20_NNMTpromoters", "p_value"], pattern = "p = .*e"))

samplesizes[samplesizes$Var1 == "H3K9_NNMTgenes", "p_value"] <- signif(wilcox.test(NCI60protein_plot_list$genes$H3K9_random, unlist(NCI60protein_plot_list$genes$H3K9_NNMT), paired = TRUE)$p.value, digits = 2)
samplesizes[samplesizes$Var1 == "H3K9_randomgenes", "p_value_parse"] <- paste0(str_remove(samplesizes[samplesizes$Var1 == "H3K9_NNMTgenes", "p_value"], pattern = "e.*"), " * 10^", str_remove(samplesizes[samplesizes$Var1 == "H3K9_NNMTgenes", "p_value"], pattern = "p = .*e"))

samplesizes[samplesizes$Var1 == "H3K9_NNMTpromoters", "p_value"] <- signif(wilcox.test(NCI60protein_plot_list$promoters$H3K9_random, unlist(NCI60protein_plot_list$promoters$H3K9_NNMT), paired = TRUE)$p.value, digits = 2)
samplesizes[samplesizes$Var1 == "H3K9_randompromoters", "p_value_parse"] <- paste0(str_remove(samplesizes[samplesizes$Var1 == "H3K9_NNMTpromoters", "p_value"], pattern = "e.*"), " * 10^", str_remove(samplesizes[samplesizes$Var1 == "H3K9_NNMTpromoters", "p_value"], pattern = "p = .*e"))

samplesizes[samplesizes$Var1 == "H3K4_NNMTgenes", "p_value"] <- signif(wilcox.test(NCI60protein_plot_list$genes$H3K4_random, unlist(NCI60protein_plot_list$genes$H3K4_NNMT), paired = TRUE)$p.value, digits = 2)
samplesizes[samplesizes$Var1 == "H3K4_randomgenes", "p_value_parse"] <- paste0(str_remove(samplesizes[samplesizes$Var1 == "H3K4_NNMTgenes", "p_value"], pattern = "e.*"), " * 10^", str_remove(samplesizes[samplesizes$Var1 == "H3K4_NNMTgenes", "p_value"], pattern = "p = .*e"))

samplesizes[samplesizes$Var1 == "H3K4_NNMTpromoters", "p_value"] <- signif(wilcox.test(NCI60protein_plot_list$promoters$H3K4_random, unlist(NCI60protein_plot_list$promoters$H3K4_NNMT), paired = TRUE)$p.value, digits = 2)
samplesizes[samplesizes$Var1 == "H3K4_randompromoters", "p_value_parse"] <- paste0(str_remove(samplesizes[samplesizes$Var1 == "H3K4_NNMTpromoters", "p_value"], pattern = "e.*"), " * 10^", str_remove(samplesizes[samplesizes$Var1 == "H3K4_NNMTpromoters", "p_value"], pattern = "p = .*e"))

samplesizes[is.na(samplesizes$p_value), "p_value"] <- ""

dat_text <- data.frame(
  label = samplesizes[samplesizes$p_value > -0.001, "p_value"],
  location   = samplesizes[samplesizes$p_value > -0.001, "location"],
  group = samplesizes[samplesizes$p_value > -0.001, "group"]
)


protein_plot_melt$group <- factor(protein_plot_melt$group, levels = c("H4K20_random", "H4K20_NNMT", "H3K9_random", "H3K9_NNMT", "H3K4_random", "H3K4_NNMT"))

pdf("graphics/NCI60_NNMTprotein_nopval.pdf",
    width = 2.25,
    height = 5)

ggplot(aes(x = signal, y = group, group = location), data = protein_plot_melt ) + 
  geom_boxplot(aes(group = group), 
               fill = c("grey",
                        "darkolivegreen3",
                        "grey",
                        "purple",
                        "grey",
                        "magenta",
                        "grey",
                        "darkolivegreen3", 
                        "grey",
                        "purple",
                        "grey",
                        "magenta"), 
               notch = TRUE,
               outlier.shape = NA,
               lwd = 0.3,
               fatten = 0.7) +
  theme_classic() +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey") +
  xlab("Linear model t-value") + 
  ylab(NULL) +
  # scale_y_discrete(label = c("H4K20me3 (random genes)",
  #                            "H4K20me3 (NNMT RNA-seq)",
  #                            "H3K9me3 (random genes)",
  #                            "H3K9me3 (NNMT RNA-seq)",
  #                            "H3K4me3 (random genes)",
  #                            "H3K4me3 (NNMT RNA-seq)")) +
  facet_wrap(~ location, ncol = 1, labeller = labeller(location = factet_labs), scales = "free") + 
  coord_cartesian(xlim = c(-4, 4)) + 
  theme(strip.text = element_text(size = 12, color = "white", face = "bold"),
        strip.background = element_rect(fill = "black"),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.title.x = element_text(size = 10, color = "black", face = "bold")) +
  geom_text(aes(x = 2, group = location, label = Freq), data = samplesizes, size = 3, hjust = 0)

dev.off()

pdf("graphics/NCI60_NNMTprotein_pval.pdf",
    width = 3.75,
    height = 5)

ggplot(aes(x = signal, y = group, group = location), data = protein_plot_melt ) + 
  geom_boxplot(aes(group = group), 
               fill = c("grey",
                        "darkolivegreen3",
                        "grey",
                        "purple",
                        "grey",
                        "magenta",
                        "grey",
                        "darkolivegreen3", 
                        "grey",
                        "purple",
                        "grey",
                        "magenta"), 
               notch = TRUE,
               outlier.shape = NA,
               lwd = 0.3,
               fatten = 0.7) +
  theme_classic() +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey") +
  xlab("Linear model t-value") + 
  ylab(NULL) +
  scale_y_discrete(label = c("H4K20me3 (random genes)",
                             "H4K20me3 (NNMT RNA-seq)",
                             "H3K9me3 (random genes)",
                             "H3K9me3 (NNMT RNA-seq)",
                             "H3K4me3 (random genes)",
                             "H3K4me3 (NNMT RNA-seq)")) +
  facet_wrap(~ location, ncol = 1, labeller = labeller(location = factet_labs), scales = "free") + 
  coord_cartesian(xlim = c(-4, 5.5)) + 
  theme(strip.text = element_text(size = 12, color = "white", face = "bold"),
        strip.background = element_rect(fill = "black"),
        # axis.text.y = element_blank(),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.title.x = element_text(size = 10, color = "black", face = "bold")) +
  geom_text(aes(x = 3.8, group = location, label = Freq), data = samplesizes, size = 3, hjust = 0) +
  geom_text(aes(x = 3.8, group = location, label = label), data = dat_text, parse = TRUE, size =  3)

dev.off()

#### PEMT RNA-seq BOXPLOT Fig S8C ####

# load RNA-seq models

H4K20PEMTRNAseq_genes_signal_model <- readRDS("output/NCI60_H4K20me3_genes_signal_model_PEMT_RNAseq.rds")
H3K9PEMTRNAseq_genes_signal_model <- readRDS("output/NCI60_H3K9me3_genes_signal_model_PEMT_RNAseq.rds")
H3K4PEMTRNAseq_genes_signal_model <- readRDS("output/NCI60_H3K4me3_genes_signal_model_PEMT_RNAseq.rds")

H4K20PEMTRNAseq_promoters_signal_model <- readRDS("output/NCI60_H4K20me3_promoters_signal_model_PEMT_RNAseq.rds")
H3K9PEMTRNAseq_promoters_signal_model <- readRDS("output/NCI60_H3K9me3_promoters_signal_model_PEMT_RNAseq.rds")
H3K4PEMTRNAseq_promoters_signal_model <- readRDS("output/NCI60_H3K4me3_promoters_signal_model_PEMT_RNAseq.rds")

H4K20random_genes_signal_model <- readRDS("output/RNAseqrandomsample_H4K20me3genes_t_values.rds")
H3K9random_genes_signal_model <- readRDS("output/RNAseqrandomsample_H3K9me3genes_t_values.rds")
H3K4random_genes_signal_model <- readRDS("output/RNAseqrandomsample_H3K4me3genes_t_values.rds")

H4K20random_promoters_signal_model <- readRDS("output/RNAseqrandomsample_H4K20me3promoters_t_values.rds")
H3K9random_promoters_signal_model <- readRDS("output/RNAseqrandomsample_H3K9me3promoters_t_values.rds")
H3K4random_promoters_signal_model <- readRDS("output/RNAseqrandomsample_H3K4me3promoters_t_values.rds")

NCI60_PEMT_RNAseq_plot_list <- list(genes = list(H3K4_PEMT = unlist(H3K4PEMTRNAseq_genes_signal_model[NCI60_H3K4_genes_signal_sites]),
                                                 H3K4_random = rowMeans(H3K4random_genes_signal_model[str_remove(row.names(H3K4random_genes_signal_model), "\\.t_value") %in% NCI60_H3K4_genes_signal_sites, ], na.rm = TRUE),
                                                 H4K20_PEMT = unlist(H4K20PEMTRNAseq_genes_signal_model[NCI60_H4K20_genes_signal_sites]),
                                                 H4K20_random = rowMeans(H4K20random_genes_signal_model[str_remove(row.names(H4K20random_genes_signal_model), "\\.t_value") %in% NCI60_H4K20_genes_signal_sites, ], na.rm = TRUE),
                                                 H3K9_PEMT = unlist(H3K9PEMTRNAseq_genes_signal_model[NCI60_H3K9_genes_signal_sites]),
                                                 H3K9_random = rowMeans(H3K9random_genes_signal_model[str_remove(row.names(H3K9random_genes_signal_model), "\\.t_value") %in% NCI60_H3K9_genes_signal_sites, ], na.rm = TRUE)),
                                    promoters = list(H3K4_PEMT = unlist(H3K4PEMTRNAseq_promoters_signal_model[NCI60_H3K4_promoters_signal_sites]),
                                                     H3K4_random = rowMeans(H3K4random_promoters_signal_model[str_remove(row.names(H3K4random_promoters_signal_model), "\\.t_value") %in% NCI60_H3K4_promoters_signal_sites, ], na.rm = TRUE),
                                                     H4K20_PEMT = unlist(H4K20PEMTRNAseq_promoters_signal_model[NCI60_H4K20_promoters_signal_sites]),
                                                     H4K20_random = rowMeans(H4K20random_promoters_signal_model[str_remove(row.names(H4K20random_promoters_signal_model), "\\.t_value") %in% NCI60_H4K20_promoters_signal_sites, ], na.rm = TRUE),
                                                     H3K9_PEMT = unlist(H3K9PEMTRNAseq_promoters_signal_model[NCI60_H3K9_promoters_signal_sites]),
                                                     H3K9_random = rowMeans(H3K9random_promoters_signal_model[str_remove(row.names(H3K9random_promoters_signal_model), "\\.t_value") %in% NCI60_H3K9_promoters_signal_sites, ], na.rm = TRUE)))           


RNAseq_plot_melt <- melt(NCI60_PEMT_RNAseq_plot_list)
colnames(RNAseq_plot_melt) <- c("signal", "group", "location")

write.table(RNAseq_plot_melt,
            file = "plot_data/Fig S9/Fig_S9C_plotdata.txt",
            col.names = TRUE,
            row.names = FALSE,
            sep = "\t")

factet_labs <- c("Gene bodies", "Promoters")
names(factet_labs) <- c("genes", "promoters")

samplesizes <- table(paste0(RNAseq_plot_melt$group, RNAseq_plot_melt$location))
samplesizes <- as.data.frame(samplesizes)
samplesizes[, "group"] <- str_remove(str_remove(str_remove(samplesizes$Var1, pattern = "promoters"), pattern = "genes"), pattern = "repeats")

for(i in 1:nrow(samplesizes)){
  samplesizes[i, "location"] <- str_remove(samplesizes[i, "Var1"], pattern = samplesizes[i, "group"])
}

samplesizes$Freq = paste0("n = ", samplesizes$Freq)

samplesizes[!str_detect(samplesizes$Var1, pattern = "random"), "Freq"] <- ""

samplesizes[samplesizes$Var1 == "H4K20_PEMTgenes", "p_value"] <- signif(wilcox.test(NCI60_PEMT_RNAseq_plot_list$genes$H4K20_random, unlist(NCI60_PEMT_RNAseq_plot_list$genes$H4K20_PEMT), paired = TRUE)$p.value, digits = 2)
samplesizes[samplesizes$Var1 == "H4K20_randomgenes", "p_value_parse"] <- paste0(str_remove(samplesizes[samplesizes$Var1 == "H4K20_PEMTgenes", "p_value"], pattern = "e.*"), " * 10^", str_remove(samplesizes[samplesizes$Var1 == "H4K20_PEMTgenes", "p_value"], pattern = "p = .*e"))

samplesizes[samplesizes$Var1 == "H4K20_PEMTpromoters", "p_value"] <- signif(wilcox.test(NCI60_PEMT_RNAseq_plot_list$promoters$H4K20_random, unlist(NCI60_PEMT_RNAseq_plot_list$promoters$H4K20_PEMT), paired = TRUE)$p.value, digits = 2)
samplesizes[samplesizes$Var1 == "H4K20_randompromoters", "p_value_parse"] <- paste0(str_remove(samplesizes[samplesizes$Var1 == "H4K20_PEMTpromoters", "p_value"], pattern = "e.*"), " * 10^", str_remove(samplesizes[samplesizes$Var1 == "H4K20_PEMTpromoters", "p_value"], pattern = "p = .*e"))

samplesizes[samplesizes$Var1 == "H3K9_PEMTgenes", "p_value"] <- signif(wilcox.test(NCI60_PEMT_RNAseq_plot_list$genes$H3K9_random, unlist(NCI60_PEMT_RNAseq_plot_list$genes$H3K9_PEMT), paired = TRUE)$p.value, digits = 2)
samplesizes[samplesizes$Var1 == "H3K9_randomgenes", "p_value_parse"] <- paste0(str_remove(samplesizes[samplesizes$Var1 == "H3K9_PEMTgenes", "p_value"], pattern = "e.*"), " * 10^", str_remove(samplesizes[samplesizes$Var1 == "H3K9_PEMTgenes", "p_value"], pattern = "p = .*e"))

samplesizes[samplesizes$Var1 == "H3K9_PEMTpromoters", "p_value"] <- signif(wilcox.test(NCI60_PEMT_RNAseq_plot_list$promoters$H3K9_random, unlist(NCI60_PEMT_RNAseq_plot_list$promoters$H3K9_PEMT), paired = TRUE)$p.value, digits = 2)
samplesizes[samplesizes$Var1 == "H3K9_randompromoters", "p_value_parse"] <- paste0(str_remove(samplesizes[samplesizes$Var1 == "H3K9_PEMTpromoters", "p_value"], pattern = "e.*"), " * 10^", str_remove(samplesizes[samplesizes$Var1 == "H3K9_PEMTpromoters", "p_value"], pattern = "p = .*e"))

samplesizes[samplesizes$Var1 == "H3K4_PEMTgenes", "p_value"] <- signif(wilcox.test(NCI60_PEMT_RNAseq_plot_list$genes$H3K4_random, unlist(NCI60_PEMT_RNAseq_plot_list$genes$H3K4_PEMT), paired = TRUE)$p.value, digits = 2)
samplesizes[samplesizes$Var1 == "H3K4_randomgenes", "p_value_parse"] <- paste0(str_remove(samplesizes[samplesizes$Var1 == "H3K4_PEMTgenes", "p_value"], pattern = "e.*"), " * 10^", str_remove(samplesizes[samplesizes$Var1 == "H3K4_PEMTgenes", "p_value"], pattern = "p = .*e"))

samplesizes[samplesizes$Var1 == "H3K4_PEMTpromoters", "p_value"] <- signif(wilcox.test(NCI60_PEMT_RNAseq_plot_list$promoters$H3K4_random, unlist(NCI60_PEMT_RNAseq_plot_list$promoters$H3K4_PEMT), paired = TRUE)$p.value, digits = 2)
samplesizes[samplesizes$Var1 == "H3K4_randompromoters", "p_value_parse"] <- paste0(str_remove(samplesizes[samplesizes$Var1 == "H3K4_PEMTpromoters", "p_value"], pattern = "e.*"), " * 10^", str_remove(samplesizes[samplesizes$Var1 == "H3K4_PEMTpromoters", "p_value"], pattern = "p = .*e"))

samplesizes[is.na(samplesizes$p_value), "p_value"] <- ""

dat_text <- data.frame(
  label = samplesizes[samplesizes$p_value > -0.001, "p_value"],
  location   = samplesizes[samplesizes$p_value > -0.001, "location"],
  group = samplesizes[samplesizes$p_value > -0.001, "group"]
)

dat_text[, "p_value_parse"] <- paste0(str_remove(dat_text$label, pattern = "e.*"), " %*% 10^", str_remove(dat_text$label, pattern = ".*e"))
dat_text[!is.na(dat_text$label) & dat_text$label == "0", "p_value_parse"] <- "0"

RNAseq_plot_melt$group <- factor(RNAseq_plot_melt$group, levels = c("H4K20_random", "H4K20_PEMT", "H3K9_random", "H3K9_PEMT", "H3K4_random", "H3K4_PEMT"))

pdf("graphics/NCI60_PEMT_RNAseq_nopval.pdf",
    width = 2.25,
    height = 5)

ggplot(aes(x = signal, y = group, group = location), data = RNAseq_plot_melt ) + 
  geom_boxplot(aes(group = group), 
               fill = c("grey",
                        "darkolivegreen3",
                        "grey",
                        "purple",
                        "grey",
                        "magenta",
                        "grey",
                        "darkolivegreen3", 
                        "grey",
                        "purple",
                        "grey",
                        "magenta"), 
               notch = TRUE,
               outlier.shape = NA,
               lwd = 0.3,
               fatten = 0.7) +
  theme_classic() +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey") +
  xlab("Linear model t-value") + 
  ylab(NULL) +
  # scale_y_discrete(label = c("H4K20me3 (random genes)",
  #                            "H4K20me3 (PEMT RNA-seq)",
  #                            "H3K9me3 (random genes)",
  #                            "H3K9me3 (PEMT RNA-seq)",
  #                            "H3K4me3 (random genes)",
  #                            "H3K4me3 (PEMT RNA-seq)")) +
  facet_wrap(~ location, ncol = 1, labeller = labeller(location = factet_labs), scales = "free") + 
  coord_cartesian(xlim = c(-2, 4)) + 
  theme(strip.text = element_text(size = 12, color = "white", face = "bold"),
        strip.background = element_rect(fill = "black"),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.title.x = element_text(size = 10, color = "black", face = "bold")) +
  geom_text(aes(x = 2.5, group = location, label = Freq), data = samplesizes, size = 3, hjust = 0)

dev.off()

pdf("graphics/NCI60_PEMT_RNAseq_pval.pdf",
    width = 3.75,
    height = 5)

ggplot(aes(x = signal, y = group, group = location), data = RNAseq_plot_melt ) + 
  geom_boxplot(aes(group = group), 
               fill = c("grey",
                        "darkolivegreen3",
                        "grey",
                        "purple",
                        "grey",
                        "magenta",
                        "grey",
                        "darkolivegreen3", 
                        "grey",
                        "purple",
                        "grey",
                        "magenta"), 
               notch = TRUE,
               outlier.shape = NA,
               lwd = 0.3,
               fatten = 0.7) +
  theme_classic() +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey") +
  xlab("Linear model t-value") + 
  ylab(NULL) +
  scale_y_discrete(label = c("H4K20me3 (random genes)",
                             "H4K20me3 (PEMT RNA-seq)",
                             "H3K9me3 (random genes)",
                             "H3K9me3 (PEMT RNA-seq)",
                             "H3K4me3 (random genes)",
                             "H3K4me3 (PEMT RNA-seq)")) +
  facet_wrap(~ location, ncol = 1, labeller = labeller(location = factet_labs), scales = "free") + 
  coord_cartesian(xlim = c(-4, 5.5)) + 
  theme(strip.text = element_text(size = 12, color = "white", face = "bold"),
        strip.background = element_rect(fill = "black"),
        # axis.text.y = element_blank(),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.title.x = element_text(size = 10, color = "black", face = "bold")) +
  geom_text(aes(x = 3.8, group = location, label = Freq), data = samplesizes, size = 3, hjust = 0) +
  geom_text(aes(x = 3.8, group = location, label = label), data = dat_text, parse = TRUE, size =  3)

dev.off()

#### REPETITIVE ELEMENTS IN THE NCI60 WITH BAYESIAN REMAPPING Fig 3D Fig S9D ####

TEfolders <- list.dirs("input/TELESCOPE/")
TEfolders <- TEfolders[2:length(TEfolders)]

TEreports <- lapply(TEfolders, function(thisfolder){

  read.table(paste0(thisfolder, "/telescope-telescope_report.tsv"),
             header = TRUE)

})

names(TEreports) <- str_extract(TEfolders, "SRR.*$")

unique_mapped_elements <- lapply(TEreports, function(x){

  unique(x[, 1])

})

unique_mapped_elements <- unique(unlist(unique_mapped_elements))

allelements_counts <- sapply(unique_mapped_elements, function(thiselement){

  sapply(TEreports, function(thissample){

    if(thiselement %in% thissample$transcript){
      return(thissample[thissample$transcript == thiselement, "final_count"])
    } else {
      return(0)
    }

  })

})

allelements_counts_df <- data.frame(t(allelements_counts))

filetable <- read.table("input/NCI60file.txt", header = TRUE, fill = TRUE)

colnames(allelements_counts_df) <- filetable[match(colnames(allelements_counts_df), filetable$Run), "LibrarySelection"]

# change names to match RNA-seq expression file
names_changed_to_match_format <- str_remove(str_remove(str_replace(str_replace(str_remove(colnames(allelements_counts_df), "^[A-Z]{2,3}:"), "/", "_"), " ", "_"), "\\(TB\\)"), "-ATCC")

names_changed_to_match_format[names_changed_to_match_format == "NCI-ADR-RES"] <- "NCI_ADR-RES"
names_changed_to_match_format[names_changed_to_match_format == "K562"] <- "K-562"
names_changed_to_match_format[names_changed_to_match_format == "IGR-OV1"] <- "IGROV1"
names_changed_to_match_format[names_changed_to_match_format == "Hs_578T"] <- "HS_578T"
names_changed_to_match_format[names_changed_to_match_format == "MDA-N"] <- "MDA-MB-468"
names_changed_to_match_format[names_changed_to_match_format == "HCC2998"] <- "HCC-2998"
names_changed_to_match_format[names_changed_to_match_format == "COLO-205"] <- "COLO_205"
names_changed_to_match_format[names_changed_to_match_format == "RXF-393"] <- "RXF_393"
names_changed_to_match_format[names_changed_to_match_format == "DU145"] <- "DU-145"

colnames(allelements_counts_df) <- names_changed_to_match_format

saveRDS(allelements_counts_df, file = "output/telescope_counts.rds")
# allelements_counts_df <- readRDS(file = "output/telescope_counts.rds")

#### match with repetitive elements for ChIP-seq ####

# point to directory with genrich output
my_directory <- "input/NCI60_SmartMap_bdgsubtract_MACS_peak_calls//"

SMpeaks_file_list <- list.files(my_directory)
SMpeaks_file_list <- SMpeaks_file_list[!str_detect(SMpeaks_file_list, "\\.gz")]

# data files for mark 

NCI60_ChIP_metadata <- read.table("input/NCI60_SRA_fordownload_single.txt",
                             sep = ",",
                             fill = TRUE,
                             quote = "")

NCI60_ChIP_metadata_colnames <- NCI60_ChIP_metadata[1, ]

NCI60_ChIP_metadata2 <- NCI60_ChIP_metadata[2:nrow(NCI60_ChIP_metadata), c(1:14, 20:ncol(NCI60_ChIP_metadata_colnames))]
colnames(NCI60_ChIP_metadata2) <- NCI60_ChIP_metadata_colnames[1:28]

#### do bed files for H3K9-marked peaks ####

H3K9_metadata <- NCI60_ChIP_metadata2[str_detect(paste0(NCI60_ChIP_metadata2$Antibody, NCI60_ChIP_metadata2$chip_antibody), "H3K9"), ]

# choose the H3K9me3 files
H3K9me3_file_list <- SMpeaks_file_list[str_detect(SMpeaks_file_list, "H3K9me3")]

# # loop through unique experiments to read in files

H3K9me3_ranges_list <- lapply(H3K9me3_file_list, function(x){

  tryCatch({

    peaks <- read.table(paste0(my_directory, x), header = FALSE, sep = "\t", skip = 1)
    colnames(peaks)[1:3] <- c("chr", "start", "end")

    peaks_GR <- makeGRangesFromDataFrame(peaks)

  }, error = function(e){

    return(NULL)

  })

})

names(H3K9me3_ranges_list) <- str_remove(str_remove(H3K9me3_file_list, "^H3K9me3_"), "_MACS2broad")

saveRDS(H3K9me3_ranges_list, "output/NCI60_SmartMap_bdgsubtract_MACS2_H3K9me3_peaks_list.rds")
# H3K9me3_ranges_list <- readRDS("output/NCI60_SmartMap_bdgsubtract_MACS2_H3K9me3_peaks_list.rds")

H3K9me3_repeats_peak_counts <- region.peak.counts(NCI60_list = H3K9me3_ranges_list,
                                                  region = "repeats")

H3K9repeats_in_40of60_number <- apply(H3K9me3_repeats_peak_counts[, 5:ncol(H3K9me3_repeats_peak_counts)], 1, function(x){
  
  sum(x > 0)
  
})

H3K9repeats_in_40of60_selec <- H3K9repeats_in_40of60_number > 39
H3K9repeats_in_40of60 <- H3K9me3_repeats_peak_counts[H3K9repeats_in_40of60_selec, "gene_id"]

H3K9repeats_in40of60_gr <- rmask_gr[match(H3K9repeats_in_40of60, mcols(rmask_gr)$myID), ]

H3K9repeats_in40of60_gr_noalt <- H3K9repeats_in40of60_gr[seqnames(H3K9repeats_in40of60_gr) %in% c(paste0("chr", 1:23), "chrX", "chrY"), ]

H3K9repeats_in40of60_gr_noalt_LINEsample <- H3K9repeats_in40of60_gr_noalt[H3K9repeats_in40of60_gr_noalt$V12 == "LINE", ][sample(names(H3K9repeats_in40of60_gr_noalt[H3K9repeats_in40of60_gr_noalt$V12 == "LINE", ]), 10000), ]
write.bed.file.for.computematrix(sampleGR = H3K9repeats_in40of60_gr_noalt_LINEsample)

H3K9repeats_in40of60_gr_noalt_SINEsample <- H3K9repeats_in40of60_gr_noalt[H3K9repeats_in40of60_gr_noalt$V12 == "SINE", ][sample(names(H3K9repeats_in40of60_gr_noalt[H3K9repeats_in40of60_gr_noalt$V12 == "SINE", ]), 10000), ]
write.bed.file.for.computematrix(sampleGR = H3K9repeats_in40of60_gr_noalt_SINEsample)

H3K9repeats_in40of60_gr_noalt_LTRsample <- H3K9repeats_in40of60_gr_noalt[H3K9repeats_in40of60_gr_noalt$V12 == "LTR", ][sample(names(H3K9repeats_in40of60_gr_noalt[H3K9repeats_in40of60_gr_noalt$V12 == "LTR", ]), 10000), ]
write.bed.file.for.computematrix(sampleGR = H3K9repeats_in40of60_gr_noalt_LTRsample)

H3K9repeats_in40of60_gr_noalt_Centromeresample <- H3K9repeats_in40of60_gr_noalt[H3K9repeats_in40of60_gr_noalt$V13 == "centr", ][sample(names(H3K9repeats_in40of60_gr_noalt[H3K9repeats_in40of60_gr_noalt$V13 == "centr", ]), 1000), ]
write.bed.file.for.computematrix(sampleGR = H3K9repeats_in40of60_gr_noalt_Centromeresample)

#### now for H4K20 ####

H4K20_metadata <- NCI60_ChIP_metadata2[str_detect(paste0(NCI60_ChIP_metadata2$Antibody, NCI60_ChIP_metadata2$chip_antibody), "H4K20"), ]

# choose the H4K20me3 files
H4K20me3_file_list <- SMpeaks_file_list[str_detect(SMpeaks_file_list, "H4K20me3")]

# loop through unique experiments to read in files

H4K20me3_ranges_list <- lapply(H4K20me3_file_list, function(x){

  tryCatch({

    peaks <- read.table(paste0(my_directory, x), header = FALSE, sep = "\t", skip = 1)
    colnames(peaks)[1:3] <- c("chr", "start", "end")

    peaks_GR <- makeGRangesFromDataFrame(peaks)

  }, error = function(e){

    return(NULL)

  })

})

names(H4K20me3_ranges_list) <- str_remove(str_remove(H4K20me3_file_list, "^H4K20me3_"), "_MACS2broad")

saveRDS(H4K20me3_ranges_list, "output/NCI60_SmartMap_bdgsubtract_MACS2_H4K20me3_peaks_list.rds")
# H4K20me3_ranges_list <- readRDS("output/NCI60_SmartMap_bdgsubtract_MACS2_H4K20me3_peaks_list.rds")

H4K20me3_repeats_peak_counts <- region.peak.counts(NCI60_list = H4K20me3_ranges_list,
                                                  region = "repeats")

H4K20repeats_in_40of60_number <- apply(H4K20me3_repeats_peak_counts[, 5:ncol(H4K20me3_repeats_peak_counts)], 1, function(x){
  
  sum(x > 0)
  
})

H4K20repeats_in_40of60_selec <- H4K20repeats_in_40of60_number > 39
H4K20repeats_in_40of60 <- H4K20me3_repeats_peak_counts[H4K20repeats_in_40of60_selec, "gene_id"]

H4K20repeats_in40of60_gr <- rmask_gr[match(H4K20repeats_in_40of60, mcols(rmask_gr)$myID), ]

H4K20repeats_in40of60_gr_noalt <- H4K20repeats_in40of60_gr[seqnames(H4K20repeats_in40of60_gr) %in% c(paste0("chr", 1:23), "chrX", "chrY"), ]

H4K20repeats_in40of60_gr_noalt_LINEsample <- H4K20repeats_in40of60_gr_noalt[H4K20repeats_in40of60_gr_noalt$V12 == "LINE", ][sample(names(H4K20repeats_in40of60_gr_noalt[H4K20repeats_in40of60_gr_noalt$V12 == "LINE", ]), 10000), ]
write.bed.file.for.computematrix(sampleGR = H4K20repeats_in40of60_gr_noalt_LINEsample)

H4K20repeats_in40of60_gr_noalt_SINEsample <- H4K20repeats_in40of60_gr_noalt[H4K20repeats_in40of60_gr_noalt$V12 == "SINE", ][sample(names(H4K20repeats_in40of60_gr_noalt[H4K20repeats_in40of60_gr_noalt$V12 == "SINE", ]), 10000), ]
write.bed.file.for.computematrix(sampleGR = H4K20repeats_in40of60_gr_noalt_SINEsample)

H4K20repeats_in40of60_gr_noalt_LTRsample <- H4K20repeats_in40of60_gr_noalt[H4K20repeats_in40of60_gr_noalt$V12 == "LTR", ][sample(names(H4K20repeats_in40of60_gr_noalt[H4K20repeats_in40of60_gr_noalt$V12 == "LTR", ]), 10000), ]
write.bed.file.for.computematrix(sampleGR = H4K20repeats_in40of60_gr_noalt_LTRsample)

H4K20repeats_in40of60_gr_noalt_Centromeresample <- H4K20repeats_in40of60_gr_noalt[H4K20repeats_in40of60_gr_noalt$V13 == "centr", ][sample(names(H4K20repeats_in40of60_gr_noalt[H4K20repeats_in40of60_gr_noalt$V13 == "centr", ]), 1000), ]
write.bed.file.for.computematrix(sampleGR = H4K20repeats_in40of60_gr_noalt_Centromeresample)

#### read in TELESCOPE HERV annotation ####

# telescope annotation obtained from [https://github.com/mlbendall/telescope_annotation_db]

readgtf <- read.table("input/retrohg38.v1.gtf", fill = TRUE)
colnames(readgtf)[c(1, 4:5, 7, 10)] <- c("Chr", "start", "end", "strand", "locusID")

readgtf_mini <- readgtf[, !str_detect(colnames(readgtf), "^V")]
readgtf_mini <- readgtf_mini[apply(readgtf_mini, 1, function(x){!any(is.na(x))}), ]

retrohg_GR <- makeGRangesFromDataFrame(readgtf_mini, keep.extra.columns = TRUE)

# H3K9 #

# find overlaps
H3K9overlaps <- data.frame(first = first(findOverlapPairs(H3K9repeats_in40of60_gr, retrohg_GR)), second = second(findOverlapPairs(H3K9repeats_in40of60_gr, retrohg_GR)))
H3K9overlaps$second.start <- H3K9overlaps$second.start-1

H3K9overlaps_match <- H3K9overlaps[apply(H3K9overlaps, 1, function(x){x["first.start"] == x["second.start"]}), ]

H3K9countsforonesIvegot <- allelements_counts_df[H3K9overlaps_match$second.locusID, ]

H3K9overlaps_names <- str_extract(row.names(H3K9countsforonesIvegot), "[0-9A-Z]+_[0-9]{1,2}[p,q]{1}[0-9]{1,2}")

H3K9finaloverlappingcounts <- H3K9countsforonesIvegot[!duplicated(H3K9overlaps_names), ]

# 1 NA out of the ~5000. save subset vector and remove to use (NA throws error otherwise)
retrohg_H3K9_selecvec <- match(row.names(H3K9finaloverlappingcounts), retrohg_GR$locusID)
retrohg_H3K9_selecvec  <- retrohg_H3K9_selecvec[!is.na(retrohg_H3K9_selecvec)]

retrohg_H3K9_select_GR <- retrohg_GR[retrohg_H3K9_selecvec, ]

retrohg_H3K9_select_df <- data.frame(seqnames = seqnames(retrohg_H3K9_select_GR),
                                     starts = start(retrohg_H3K9_select_GR)-1,
                                     ends = end(retrohg_H3K9_select_GR),
                                     gene_id = retrohg_H3K9_select_GR$locusID,
                                     scores = c(rep(".", length(retrohg_H3K9_select_GR))),
                                     strands = strand(retrohg_H3K9_select_GR))

options(scipen=999)
write.table(retrohg_H3K9_select_df, file = "output/retrohg_NCI60_H3K9_select.bed", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

# H4K20 #

H4K20overlaps <- data.frame(first = S4Vectors::first(findOverlapPairs(H4K20repeats_in40of60_gr, retrohg_GR)), second = second(findOverlapPairs(H4K20repeats_in40of60_gr, retrohg_GR)))
H4K20overlaps$second.start <- H4K20overlaps$second.start-1

H4K20overlaps_match <- H4K20overlaps[apply(H4K20overlaps, 1, function(x){x["first.start"] == x["second.start"]}), ]

H4K20countsforonesIvegot <- allelements_counts_df[H4K20overlaps_match$second.locusID, ]

H4K20overlaps_names <- str_extract(row.names(H4K20countsforonesIvegot), "[0-9A-Z]+_[0-9]{1,2}[p,q]{1}[0-9]{1,2}")

H4K20finaloverlappingcounts <- H4K20countsforonesIvegot[!duplicated(H4K20overlaps_names), ]

retrohg_H4K20_selecvec <- match(row.names(H4K20finaloverlappingcounts), retrohg_GR$locusID)
retrohg_H4K20_selecvec  <- retrohg_H4K20_selecvec[!is.na(retrohg_H4K20_selecvec)]

retrohg_H4K20_select_GR <- retrohg_GR[retrohg_H4K20_selecvec, ]

retrohg_H4K20_select_df <- data.frame(seqnames = seqnames(retrohg_H4K20_select_GR),
                                     starts = start(retrohg_H4K20_select_GR)-1,
                                     ends = end(retrohg_H4K20_select_GR),
                                     gene_id = retrohg_H4K20_select_GR$locusID,
                                     scores = c(rep(".", length(retrohg_H4K20_select_GR))),
                                     strands = strand(retrohg_H4K20_select_GR))

options(scipen=999)
write.table(retrohg_H4K20_select_df, file = "output/retrohg_NCI60_H4K20_select.bed", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

#### now have computematrix output ####

make.lm.matrix.from.files <- function(signalfiles_vec,
                                      histonemark,
                                      output_name,
                                      fileindicator,
                                      saveoutput = FALSE){

temp_signal_list <- lapply(signalfiles_vec, function(x){
  
  rep1 <- read.table(paste0("input/bigwigcompare_out/computematrix_out/", x), skip = 1)
  
})

names(temp_signal_list) <- str_remove(signalfiles_vec, paste0("_bigwigcompare.bigWig_", fileindicator,"_computematrix.gz"))

readmatrix_df <- sapply(temp_signal_list, function(x){
  
  x$V7
  
})

row.names(readmatrix_df) <- temp_signal_list[[1]]$V4

temp_matrix_for_lm <- data.frame(t(readmatrix_df))

temp_matrix_for_lm[ , "tissue"] <- str_extract(str_remove(row.names(temp_matrix_for_lm), paste0(histonemark, "_")), "[A-Z]{2,3}")

if(!exists("NNMT_proteomics")){
  NNMT_proteomics <- read.table("input/NNMT_proteomics.txt", sep = "\t", header = TRUE)
}

if(!exists("NNMT_mRNA")){
  NNMT_mRNA <- read.table("input/NNMT_mRNA_zscores.txt", sep = "\t", header = TRUE)
}

# here I am changing the names from the NNMT proteomics and gene expression file to match the format
# replace slashes with underscores
# replace spaces with underscores
# remove "(TB)" from HL-60
# remove "_ATCC" from A549
names_changed_to_match_format <- str_remove(str_remove(str_replace(str_replace(str_remove(NNMT_proteomics$Cell.Line, "^[A-Z]{2,3}:"), "/", "_"), " ", "_"), "\\(TB\\)"), "_ATCC")
names_changed_to_match_format_RNA <- str_remove(str_remove(str_replace(str_replace(str_remove(NNMT_mRNA$Cell.Line, "^[A-Z]{2,3}:"), "/", "_"), " ", "_"), "\\(TB\\)"), "_ATCC")

# why 59?
# only 59 qvailable data points for proteomics
# 60 for expression so will download expression data separately from CellMinerDB

# the 60th in the RNA dataset is "MDA-N". Will assume it is the last one in the panel which is MDA-MB-468
names_changed_to_match_format_RNA[!sapply(names_changed_to_match_format_RNA, function(x){
  any(str_detect(row.names(temp_matrix_for_lm), x))
})]

names_changed_to_match_format_RNA <- str_replace(names_changed_to_match_format_RNA, pattern = "MDA-N", replacement = "MDA-MB-468")

NNMT_proteomics$Cell.Line <- names_changed_to_match_format
NNMT_mRNA$Cell.Line <- names_changed_to_match_format_RNA

temp_matrix_for_lm[, "NNMTprotein"] <- NNMT_proteomics[match(str_remove(str_remove(row.names(temp_matrix_for_lm), pattern = paste0(histonemark, "_")), pattern = "^[A-Z]{2,3}_"), NNMT_proteomics$Cell.Line), "swaNNMT_nci60"]

temp_matrix_for_lm[, "NNMTmRNA_z"] <- NNMT_mRNA[match(str_remove(str_remove(row.names(temp_matrix_for_lm), pattern = paste0(histonemark, "_")), pattern = "^[A-Z]{2,3}_"), NNMT_mRNA$Cell.Line), "expNNMT_nci60"]

temp_matrix_for_lm[, "NNMT_RNAseq"] <- NCI60_normalised_counts[NNMT_entrezgene_id, match(str_remove(str_remove(row.names(temp_matrix_for_lm), pattern = paste0(histonemark, "_")), pattern = "^[A-Z]{2,3}_"), colnames(NCI60_normalised_counts))]
temp_matrix_for_lm[, "PEMT_RNAseq"] <- NCI60_normalised_counts[PEMT_entrezgene_id, match(str_remove(str_remove(row.names(temp_matrix_for_lm), pattern = paste0(histonemark, "_")), pattern = "^[A-Z]{2,3}_"), colnames(NCI60_normalised_counts))]

# remove columns with only NAs
temp_matrix_for_lm <- temp_matrix_for_lm[, apply(temp_matrix_for_lm, 2, function(x){any(!is.na(x))})]

temp_number_of_sites <- ncol(temp_matrix_for_lm) - 5

output_list <- list(matrix_for_lm = temp_matrix_for_lm,
                    number_of_sites = temp_number_of_sites)

if(saveoutput == TRUE){
saveRDS(output_list, paste0("output/", output_name,".rds"))
}

return(output_list)
  
}

computematrix_files <- list.files("input/bigwigcompare_out/computematrix_out/")

myNCI60histonemark_vec <- c("H3K9",
                            "H4K20")

myrepelements_vec <- c("HERVs",
                       "LINEs",
                       "SINEs",
                       "LTRs",
                       "Centromeres")

NCI60marks_NNMT_RNAseq_matrix_for_lm_list <- lapply(myNCI60histonemark_vec, function(thismark){
  
  temp_list <- lapply(myrepelements_vec, function(thiselement){
    
    thismarkfiles <- computematrix_files[str_detect(computematrix_files, thismark)]
   temp_files <- thismarkfiles[str_detect(thismarkfiles, str_remove(thiselement, "s$"))]
    
   temp_matrix_for_lm_out <- make.lm.matrix.from.files(signalfiles_vec = temp_files,
                                                              histonemark = paste0(thismark, "me3"),
                                                              output_name = paste0("NCI60_", thismark, "_", thiselement,"_matrix_for_lm_out"),
                                                              fileindicator = paste0(str_remove(thiselement, "s$"), "SAMPLE"),
                                                              saveoutput = TRUE)
    
   names(temp_matrix_for_lm_out) <- c(paste0(thismark, thiselement, "_matrix_for_lm"),
                                      paste0(thismark, thiselement, "number_of_sites"))
   
   saveRDS(temp_matrix_for_lm_out, paste0("output/", thismark, thiselement, "_matrix_for_lm_out.rds"))
   
   return(temp_matrix_for_lm_out)
    
  })
  
  names(temp_list) <- paste0(thismark, myrepelements_vec, "_matrix_for_lm_out")
  
  return(temp_list)
  
})

# # if already calculated load here
# 
# NCI60marks_NNMT_RNAseq_matrix_for_lm_list <- lapply(myNCI60histonemark_vec, function(thismark){
# 
#   lapply(myrepelements_vec, function(thiselement){
# 
#     readRDS(paste0("output/", thismark, thiselement, "_matrix_for_lm_out.rds"))
# 
#   })
# 
# })

names(NCI60marks_NNMT_RNAseq_matrix_for_lm_list) <- myNCI60histonemark_vec

for(i in 1:length(myNCI60histonemark_vec)){
  
  for(j in 1:length(myrepelements_vec)){

    assign(names(NCI60marks_NNMT_RNAseq_matrix_for_lm_list[[i]][[j]][1]), NCI60marks_NNMT_RNAseq_matrix_for_lm_list[[i]][[j]][[1]])
    assign(names(NCI60marks_NNMT_RNAseq_matrix_for_lm_list[[i]][[j]][2]), NCI60marks_NNMT_RNAseq_matrix_for_lm_list[[i]][[j]][[2]])
    
  }
  
}

#### fit models ####

myNCI60histonemark_vec <- c("H3K9",
                            "H4K20")

myrepelements_vec <- c("HERVs",
                       "LINEs",
                       "SINEs",
                       "LTRs",
                       "Centromeres")

NCI60marks_NNMT_RNAseq_elements_signal_model_list <- lapply(myNCI60histonemark_vec, function(thismark){
  
  temp_list <- lapply(myrepelements_vec, function(thiselement){
    
    temp_model <- fit.model.to.NCI60.chip(histone_mark = thismark,
                            siteofchoice = thiselement,
                            gene_enquiry = "NNMT",
                            expressiontype = "_RNAseq")
    
    saveRDS(temp_model, 
            paste0("output/", thismark, "NNMT_RNAseq_", thiselement, "_signal_model.rds"))
    
  })
  
  names(temp_list) <- paste0(thismark, gene, "_RNAseq_", myrepelements_vec, "_signal_model")
  
  return(temp_list)
  
})

# # if already calculated load here
# 
# NCI60marks_NNMT_RNAseq_elements_signal_model_list <- lapply(myNCI60histonemark_vec, function(thismark){
# 
#   lapply(myrepelements_vec, function(thiselement){
# 
#     readRDS(paste0("output/", thismark, "NNMT_RNAseq_", thiselement, "_signal_model.rds"))
# 
#   })
# 
# })
# 

names(NCI60marks_NNMT_RNAseq_elements_signal_model_list) <- myNCI60histonemark_vec

for(i in 1:length(names(NCI60marks_NNMT_RNAseq_elements_signal_model_list))){
  
  names(NCI60marks_NNMT_RNAseq_elements_signal_model_list[[i]]) <- paste0(names(NCI60marks_NNMT_RNAseq_elements_signal_model_list)[[i]], "NNMT_RNAseq_", myrepelements_vec, "_signal_model")
  
}

for(i in 1:length(myNCI60histonemark_vec)){
  
  for(j in 1:length(myrepelements_vec)){

    assign(names(NCI60marks_NNMT_RNAseq_elements_signal_model_list[[i]][j]), NCI60marks_NNMT_RNAseq_elements_signal_model_list[[i]][[j]])
    
  }
  
}

#### signal vs t dfs ####

do.signal.vs.t.df <- function(histonemark,
                              elementname,
                              expressiontype,
                              gene = "NNMT",
                              no_of_bins = 5,
                              save = TRUE){
  
  if(!exists(paste0(histonemark, elementname, "matrix_for_lm"))){
    
    dataobject <- readRDS(paste0("output/NCI60_", histonemark, "_", elementname, "_matrix_for_lm_out.rds"))
    datamatrix <- dataobject$matrix_for_lm
    assign(paste0(histonemark, elementname, "matrix_for_lm"), datamatrix)
    
  }
  
  n = get(paste0(histonemark, elementname, "number_of_sites"))
  
  tempsample_expression_vs_t <- data.frame(t = unlist(get(paste0(histonemark, gene, "_", expressiontype, "_", elementname, "_signal_model"))),
                                       signal = colMeans(get(paste0(histonemark, elementname, "matrix_for_lm"))[, 1:n])[str_remove(names(unlist(get(paste0(histonemark, gene, "_", expressiontype, "_", elementname, "_signal_model")))), "\\.t_value")])

  tempsample_expression_vs_t_bins <- tempsample_expression_vs_t %>% mutate(points_bin = ntile(signal,  n = no_of_bins))
  
  if(save == TRUE){

        saveRDS(tempsample_expression_vs_t_bins, 
            paste0("output/", histonemark, gene, "_", expressiontype, "_", elementname, "_signal_vs_t_df.rds"))
    
  }
  
  return(tempsample_expression_vs_t_bins)
  
}

# H3K9

myNCI60histonemark_vec <- c("H3K9",
                            "H4K20")

myrepelements_vec <- c("HERVs",
                       "LINEs",
                       "SINEs",
                       "LTRs",
                       "Centromeres")

NCI60marks_NNMT_RNAseq_elements_signal_vs_t_list <- lapply(myNCI60histonemark_vec, function(thismark){
  
  temp_list <- lapply(myrepelements_vec, function(thiselement){
    
    do.signal.vs.t.df(histonemark = thismark,
                      gene = "NNMT",
                      expressiontype = "RNAseq",
                      elementname = thiselement,
                      save = TRUE)
    
  })
  
  names(temp_list) <- paste0(thismark, "NNMT_RNAseq_", myrepelements_vec, "_signal_vs_t_df")
  
  return(temp_list)
  
})

names(NCI60marks_NNMT_RNAseq_elements_signal_vs_t_list) <- myNCI60histonemark_vec

for(i in 1:length(myNCI60histonemark_vec)){
  
  for(j in 1:length(myrepelements_vec)){

    assign(names(NCI60marks_NNMT_RNAseq_elements_signal_vs_t_list[[i]][j]), NCI60marks_NNMT_RNAseq_elements_signal_vs_t_list[[i]][[j]])
    
  }
  
}

# NCI60marks_NNMT_RNAseq_elements_signal_vs_t_all_df <- reshape2::melt(NCI60marks_NNMT_RNAseq_elements_signal_vs_t_list)
NCI60marks_NNMT_RNAseq_elements_signal_vs_t_H3K9_df <- do.call(rbind, NCI60marks_NNMT_RNAseq_elements_signal_vs_t_list[["H3K9"]])
NCI60marks_NNMT_RNAseq_elements_signal_vs_t_H4K20_df <- do.call(rbind, NCI60marks_NNMT_RNAseq_elements_signal_vs_t_list[["H4K20"]])

NCI60marks_NNMT_RNAseq_elements_signal_vs_t_all_df <- rbind(NCI60marks_NNMT_RNAseq_elements_signal_vs_t_H3K9_df ,
                                                            NCI60marks_NNMT_RNAseq_elements_signal_vs_t_H4K20_df)

NCI60marks_NNMT_RNAseq_elements_signal_vs_t_all_df[, "location"] <- str_remove(str_remove(row.names(NCI60marks_NNMT_RNAseq_elements_signal_vs_t_all_df), "^[A-Z0-9]{8,9}_RNAseq_"), "_signal_vs_t.*$")
NCI60marks_NNMT_RNAseq_elements_signal_vs_t_all_df[, "mark"] <- paste0(str_extract(row.names(NCI60marks_NNMT_RNAseq_elements_signal_vs_t_all_df), "^H[3-4]K[0-9]{1,2}"), "me3")

breaks <- c(-3, -1.5, 0, 1.5, 3)
labels <- as.character(breaks)
labels[!(breaks %% 3 == 0)] <- ''

saveRDS(NCI60marks_NNMT_RNAseq_elements_signal_vs_t_all_df, "plot_data/NCI60marks_NNMT_RNAseq_elements_signal_vs_t_all_df.rds")

pdf("graphics/NCI60_repelements_HXKXsignal_by_bin2.pdf",
    width = 4, 
    height = 6)

ggplot(data = NCI60marks_NNMT_RNAseq_elements_signal_vs_t_all_df, aes(x = points_bin, group = points_bin, y = t)) + 
  geom_boxplot(fill = rep(c(rep("magenta", times = 5),
                        rep("dark orchid", times = 5)), times = 5),
               outlier.size = 0.1) +
  theme_classic() + 
  xlab('Histone methylation mark signal bin') + 
  ylab(italic(NNMT)~"linear model t-value") + 
  theme(strip.text = element_text(size = 12, color = "white", face = "bold"),
        strip.background = element_rect(fill = "black"),
        axis.text.x = element_blank(),
        axis.text.y = element_text(colour = "black",
                                   size = 10),
        panel.border = element_rect(colour = "black",
                                    fill = "transparent"),
        panel.spacing.y = unit(0, "line")) +
  geom_hline(yintercept = 0,
             linetype = "dashed",
             colour = "dark grey") +
  coord_cartesian(ylim = c(-4, 4)) +
  scale_y_continuous(breaks = c(-3, -1.5, 0, 1.5, 3),
                     labels = labels) +
  facet_grid(location ~ mark) 

dev.off()

#### comparison of relationship by signal Fig S9D Fig S9E ####

NCI60marks_NNMT_RNAseq_elements_signal_vs_t_list_copy <- NCI60marks_NNMT_RNAseq_elements_signal_vs_t_list

NCI60marks_NNMT_RNAseq_elements_signal_vs_t_list_copy <- lapply(names(NCI60marks_NNMT_RNAseq_elements_signal_vs_t_list_copy), function(x){

  templist <- lapply(NCI60marks_NNMT_RNAseq_elements_signal_vs_t_list_copy[[x]], function(y){
      
      y[, "mark"] <- paste0(x, "me3")
      
      return(y)
    
  })
  
  return(templist)
  
})

NCI60marks_NNMT_RNAseq_elements_signal_vs_t_list_copy <- lapply(NCI60marks_NNMT_RNAseq_elements_signal_vs_t_list_copy, function(thismark){

  lapply(names(thismark), function(thiselement){
    
    thismark[[thiselement]][, "location"] <- str_remove(str_remove(thiselement, "H.*_RNAseq_"), "_signal_vs_t_df")
    
    return(thismark[[thiselement]])
    
  })
  
})

all_sample_expression_vs_t_melt_df <- do.call(c, NCI60marks_NNMT_RNAseq_elements_signal_vs_t_list_copy)
all_sample_expression_vs_t_melt_df <- do.call(rbind, all_sample_expression_vs_t_melt_df)
all_sample_expression_vs_t_melt_df[, "element"] <- str_remove(row.names(all_sample_expression_vs_t_melt_df), "\\.chr.*")

saveRDS(all_sample_expression_vs_t_melt_df, "plot_data/all_sample_expression_vs_t_melt_df.rds")
# all_sample_expression_vs_t_melt_df <- readRDS("plot_data/all_sample_expression_vs_t_melt_df.rds")

write.table(all_sample_expression_vs_t_melt_df,
            file = "plot_data/Fig S9/Fig_S9E_plotdata.txt",
            row.names = TRUE,
            col.names = TRUE,
            sep = "\t")

pdf("graphics/NCI60_repelements_t_by_signal.pdf",
    width = 3.5,
    height = 2)

ggplot(data = all_sample_expression_vs_t_melt_df, aes(x = signal, y = t, group = location, colour = location)) + 
  geom_smooth(method =  lm) + 
  theme_classic() + 
  theme(strip.text = element_text(size = 12, color = "white", face = "bold"),
        strip.background = element_rect(fill = "black"), 
        axis.text = element_text(size = 10, 
                                  colour = "black"),
        legend.position = "none") +
  xlab("Signal ("~log[2]~"fold change over input )") + 
  ylab("Linear model t value") + 
  facet_wrap(~ mark) + 
  geom_vline(xintercept = 0,
             linetype = "dashed",
             colour = "grey") +
  geom_hline(yintercept = 0,
             linetype = "dashed",
             colour = "grey")

dev.off()

#### RANDOM DISTRIBUTIONS ####

nuclear_genes_expressed <- readRDS("output/GTEX_nuclear_genes_expressed.rds")
nuclear_genes_expressed_entrez <- getBM(attributes = c("ensembl_gene_id", "entrezgene_id"),
                                        values = nuclear_genes_expressed,
                                        filter = "ensembl_gene_id",
                                        mart = ensembl)

saveRDS(nuclear_genes_expressed_entrez, "output/nuclear_genes_expressed_entrez.rds")
# nuclear_genes_expressed_entrez <- readRDS("output/nuclear_genes_expressed_entrez.rds")

RNAseq_data_expressed <- NCI60_normalised_counts[row.names(NCI60_normalised_counts) %in% nuclear_genes_expressed_entrez$entrezgene_id, ]
row.names(RNAseq_data_expressed) <- nuclear_genes_expressed_entrez[match(row.names(RNAseq_data_expressed), nuclear_genes_expressed_entrez$entrezgene_id), "ensembl_gene_id"]

# random_expressed_RNAseq_sample <- sample(row.names(RNAseq_data_expressed), 1000, replace = FALSE)
# saveRDS(random_expressed_RNAseq_sample, "output/random_expressed_RNAseq_sample.rds")
random_expressed_RNAseq_sample <- readRDS("output/random_expressed_RNAseq_sample.rds")

sapply(myNCI60histonemark_vec, function(thismark){
  
  sapply(myrepelements_vec)
  
})

random_modified_matrices_list <- lapply(myNCI60histonemark_vec, function(thismark){
  
  temp_list <- lapply(myrepelements_vec, function(thiselement){
    
    tempobject <- get(paste0(thismark, thiselement, "matrix_for_lm"))
    
    for(i in 1:length(random_expressed_RNAseq_sample)){
      
      tempobject[, paste0(random_expressed_RNAseq_sample[i], "_RNAseq")] <- as.numeric(RNAseq_data_expressed[random_expressed_RNAseq_sample[i], match(str_remove(row.names(tempobject), paste0(thismark, "me3_[A-Z]{2,3}_")), colnames(RNAseq_data_expressed))])
      
    }
    
    return(tempobject)
    
  })
  
  names(temp_list) <- paste0(thismark, "_", myrepelements_vec, "matrix_for_lm_plusrandom")

  return(temp_list)
  
})

for(i in 1:length(myNCI60histonemark_vec)){
  
  for(j in 1:length(myrepelements_vec)){
    
    assign(names(random_modified_matrices_list[[i]][j]), random_modified_matrices_list[[i]][[j]])
    
  }
  
}

NCI60_repelements_random_t_list <- lapply(myNCI60histonemark_vec, function(thismark){
  
  temp_list <- lapply(myrepelements_vec, function(thiselement){
    
    temp_t_values <- sapply(random_expressed_RNAseq_sample[1:400], function(thisgene){
      
      sapply(str_remove(row.names(get(paste0(thismark, "NNMT_RNAseq_", thiselement, "_signal_vs_t_df"))[get(paste0(thismark, "NNMT_RNAseq_", thiselement, "_signal_vs_t_df"))$points_bin == 5, ]), "\\.t_value$"), function(x){
        
        tryCatch(expr = {
          
          glmtemp <- glm(formula = as.formula(paste0(x, " ~ log10(", as.character(thisgene), "_RNAseq + 1) + tissue + log10(", thisgene, "_RNAseq + 1) * tissue")), data = get(paste0(thismark, "_", thiselement, "matrix_for_lm_plusrandom")))
          glmsum <- summary(glmtemp)
          
          t <- glmsum$coefficients[paste0("log10(", thisgene, "_RNAseq + 1)"), "t value"]
          
          output_vec <- c(t)
          names(output_vec) <- c("t_value")
          
          return(output_vec)
          
        }, error = function(e){
          
          output_vec <- c(NA)
          names(output_vec) <- c("t_value")
          
          return(output_vec)
          
        }
        )
        
      })
      
    }) # end of internal loop
    
  })
  
  names(temp_list) <- paste0(thismark, "NNMT","_", myrepelements_vec, "_random_genes_signal_model")

  sapply(myrepelements_vec, function(thiselement){

    saveRDS(temp_list[[paste0(thismark, "NNMT","_", thiselement, "_random_genes_signal_model")]], paste0("output/", thismark, "NNMT","_", thiselement, "_random_genes_signal_model.rds"))
    
  })
  
  return(temp_list)
  
})

# ## Use comented code below to reconstitute this list run previously
# 
# NCI60_repelements_random_t_list <- lapply(myNCI60histonemark_vec, function(thismark){
# 
#   templist <- lapply(myrepelements_vec, function(thiselement){
# 
#     readRDS(paste0("output/", thismark, "NNMT","_", thiselement, "_random_genes_signal_model.rds"))
# 
#   })
# 
#   names(templist) <-  paste0(thismark, "NNMT","_", myrepelements_vec, "_random_genes_signal_model")
# 
#   templist
# 
# })
# 
# names(NCI60_repelements_random_t_list) <- myNCI60histonemark_vec

NCI60_repelements_random_mean_t_list <- lapply(NCI60_repelements_random_t_list, function(x){
  
  lapply(x, rowMeans, na.rm = TRUE)
  
})

for(i in 1:length(myNCI60histonemark_vec)){

  for(j in 1:length(myrepelements_vec)){
    
    assign(paste0(names(NCI60_repelements_random_mean_t_list[[i]][j]), "_avg"), NCI60_repelements_random_mean_t_list[[i]][[j]])
    
  }
  
}

#### make plot for NCI60 repetitive elements Fig 3E ####

H4K20NNMT_RNAseq_Centromeres_signal_vs_t_df <- readRDS("output/H4K20NNMT_RNAseq_Centromeres_signal_vs_t_df.rds")
H3K9NNMT_RNAseq_Centromeres_signal_vs_t_df <- readRDS("output/H3K9NNMT_RNAseq_Centromeres_signal_vs_t_df.rds")

H4K20NNMT_RNAseq_HERVs_signal_vs_t_df <- readRDS("output/H4K20NNMT_RNAseq_HERVs_signal_vs_t_df.rds")
H3K9NNMT_RNAseq_HERVs_signal_vs_t_df <- readRDS("output/H3K9NNMT_RNAseq_HERVs_signal_vs_t_df.rds")

H4K20NNMT_RNAseq_LINEs_signal_vs_t_df <- readRDS("output/H4K20NNMT_RNAseq_LINEs_signal_vs_t_df.rds")
H3K9NNMT_RNAseq_LINEs_signal_vs_t_df <- readRDS("output/H3K9NNMT_RNAseq_LINEs_signal_vs_t_df.rds")

H4K20NNMT_RNAseq_LTRs_signal_vs_t_df <- readRDS("output/H4K20NNMT_RNAseq_LTRs_signal_vs_t_df.rds")
H3K9NNMT_RNAseq_LTRs_signal_vs_t_df <- readRDS("output/H3K9NNMT_RNAseq_LTRs_signal_vs_t_df.rds")

H4K20NNMT_RNAseq_SINEs_signal_vs_t_df <- readRDS("output/H4K20NNMT_RNAseq_SINEs_signal_vs_t_df.rds")
H3K9NNMT_RNAseq_SINEs_signal_vs_t_df <- readRDS("output/H3K9NNMT_RNAseq_SINEs_signal_vs_t_df.rds")

H4K20NNMT_Centromeres_random_genes_signal_model <- readRDS("output/H4K20NNMT_Centromeres_random_genes_signal_model.rds")

NCI60_RNAseq_repeats_plot_list <- list(Centromeres = list(H4K20_NNMT = H4K20NNMT_RNAseq_Centromeres_signal_vs_t_df[H4K20NNMT_RNAseq_Centromeres_signal_vs_t_df$points_bin == 5, "t"],
                                                      H4K20_random = H4K20NNMT_Centromeres_random_genes_signal_model_avg,
                                                      H3K9_NNMT = H3K9NNMT_RNAseq_Centromeres_signal_vs_t_df[H3K9NNMT_RNAseq_Centromeres_signal_vs_t_df$points_bin == 5, "t"],
                                                      H3K9_random = H3K9NNMT_Centromeres_random_genes_signal_model_avg),
                                       HERVs = list(H4K20_NNMT = H4K20NNMT_RNAseq_HERVs_signal_vs_t_df[H4K20NNMT_RNAseq_HERVs_signal_vs_t_df$points_bin == 5, "t"],
                                                    H4K20_random = H4K20NNMT_HERVs_random_genes_signal_model_avg,
                                                    H3K9_NNMT = H3K9NNMT_RNAseq_HERVs_signal_vs_t_df[H3K9NNMT_RNAseq_HERVs_signal_vs_t_df$points_bin == 5, "t"],
                                                    H3K9_random = H3K9NNMT_HERVs_random_genes_signal_model_avg),
                                       LINEs = list(H4K20_NNMT = H4K20NNMT_RNAseq_LINEs_signal_vs_t_df[H4K20NNMT_RNAseq_LINEs_signal_vs_t_df$points_bin == 5, "t"],
                                                    H4K20_random = H4K20NNMT_LINEs_random_genes_signal_model_avg,
                                                    H3K9_NNMT = H3K9NNMT_RNAseq_LINEs_signal_vs_t_df[H3K9NNMT_RNAseq_LINEs_signal_vs_t_df$points_bin == 5, "t"],
                                                    H3K9_random = H3K9NNMT_LINEs_random_genes_signal_model_avg),
                                       LTRs = list(H4K20_NNMT = H4K20NNMT_RNAseq_LTRs_signal_vs_t_df[H4K20NNMT_RNAseq_LTRs_signal_vs_t_df$points_bin == 5, "t"],
                                                    H4K20_random = H4K20NNMT_LTRs_random_genes_signal_model_avg,
                                                    H3K9_NNMT = H3K9NNMT_RNAseq_LTRs_signal_vs_t_df[H3K9NNMT_RNAseq_LTRs_signal_vs_t_df$points_bin == 5, "t"],
                                                    H3K9_random = H3K9NNMT_LTRs_random_genes_signal_model_avg),
                                       SINEs = list(H4K20_NNMT = H4K20NNMT_RNAseq_SINEs_signal_vs_t_df[H4K20NNMT_RNAseq_SINEs_signal_vs_t_df$points_bin == 5, "t"],
                                                    H4K20_random = H4K20NNMT_SINEs_random_genes_signal_model_avg,
                                                    H3K9_NNMT = H3K9NNMT_RNAseq_SINEs_signal_vs_t_df[H3K9NNMT_RNAseq_SINEs_signal_vs_t_df$points_bin == 5, "t"],
                                                    H3K9_random = H3K9NNMT_SINEs_random_genes_signal_model_avg))           


RNAseq_plot_melt <- melt(NCI60_RNAseq_repeats_plot_list)
colnames(RNAseq_plot_melt) <- c("signal", "group", "location")

write.table(RNAseq_plot_melt,
            file = "plot_data/Fig 3/Fig_3D_plotdata.txt",
            col.names = TRUE,
            row.names = FALSE,
            sep = "\t")

# factet_labs <- c("Centromeres", "HERVs", "LINEs", "LTRs", "SINEs")
# names(factet_labs) <- c("genes", "promoters", "repeats")

samplesizes <- table(paste0(RNAseq_plot_melt$group, RNAseq_plot_melt$location))
samplesizes <- as.data.frame(samplesizes)
samplesizes[, "group"] <- str_remove(str_remove(str_remove(str_remove(str_remove(samplesizes$Var1, pattern = "Centromeres"), pattern = "HERVs"), pattern = "LINEs"), pattern = "LTRs"), pattern = "SINEs")

for(i in 1:nrow(samplesizes)){
  samplesizes[i, "location"] <- str_remove(samplesizes[i, "Var1"], pattern = samplesizes[i, "group"])
}

samplesizes$Freq = paste0("n = ", samplesizes$Freq)

samplesizes[!str_detect(samplesizes$Var1, pattern = "random"), "Freq"] <- ""

for(i in 1:length(myNCI60histonemark_vec)){
  
  for(j in 1:length(myrepelements_vec)){

samplesizes[samplesizes$Var1 == paste0(myNCI60histonemark_vec[i], "_NNMT", myrepelements_vec[j]), "p_value"] <- signif(wilcox.test(unlist(NCI60_RNAseq_repeats_plot_list[[myrepelements_vec[j]]][[paste0(myNCI60histonemark_vec[i], "_random")]]),
                                                                                                                                   unlist(NCI60_RNAseq_repeats_plot_list[[myrepelements_vec[j]]][[paste0(myNCI60histonemark_vec[i], "_NNMT")]]),
                                                                                                                                   paired = TRUE)$p.value,
                                                                                                                       digits = 2)

samplesizes[samplesizes$Var1 == paste0(myNCI60histonemark_vec[i], "_random", myrepelements_vec[j]), "p_value_parse"] <- paste0(as.character(round(as.numeric(str_remove(samplesizes[samplesizes$Var1 == paste0(myNCI60histonemark_vec[i], "_NNMT", myrepelements_vec[j]), "p_value"], pattern = "e.*")), digits = 2)),
                                                                                                                               " * 10^",
                                                                                                                               str_remove(samplesizes[samplesizes$Var1 == paste0(myNCI60histonemark_vec[i], "_NNMT", myrepelements_vec[j]), "p_value"], pattern = ".*e"))
 
 }

}


samplesizes[is.na(samplesizes$p_value), "p_value"] <- ""

dat_text <- data.frame(
  label = samplesizes[samplesizes$p_value > -0.001, "p_value"],
  location   = samplesizes[samplesizes$p_value > -0.001, "location"],
  group = samplesizes[samplesizes$p_value > -0.001, "group"]
)

dat_text[, "p_value_parse"] <- paste0(as.character(round(as.numeric(str_remove(dat_text$label, pattern = "e.*")), digits = 2)), " %*% 10^", str_remove(dat_text$label, pattern = ".*e"))
dat_text[!is.na(dat_text$label) & dat_text$label == "0", "p_value_parse"] <- "0"

RNAseq_plot_melt$group <- factor(RNAseq_plot_melt$group, levels = c("H3K9_random", "H3K9_NNMT", "H4K20_random", "H4K20_NNMT"))

pdf("graphics/NCI60_repetitive_elements_boxplot.pdf",
    width = 4.5, height = 6)

ggplot(aes(x = signal, y = group, group = location), data = RNAseq_plot_melt ) + 
  geom_boxplot(aes(group = group), fill = rep(c("grey",
                                            "darkorchid1",
                                            "grey",
                                            "magenta"
                                        ), times = 5), notch = TRUE, outlier.shape = NA) +
  theme_classic() +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey") +
  xlab("Linear model t-value") + 
  ylab(NULL) +
  scale_y_discrete(label = c("H4K20me3 (random genes)",
                             "H4K20me3 (NNMT RNA-seq)",
                                "H3K9me3 (random genes)",
                             "H3K9me3 (NNMT RNA-seq)")
                             ) +
  facet_wrap(~ location, ncol = 1, scales = "free") +
  coord_cartesian(xlim = c(-4, 5.5)) + 
  theme(strip.text = element_text(size = 12, color = "white", face = "bold"),
        strip.background = element_rect(fill = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.title.x = element_text(size = 10, color = "black", face = "bold")) +
  geom_text(aes(x = 4.75, group = location, label = p_value_parse), data = dat_text, parse = TRUE, size = 3) +
  geom_text(aes(x = 4.75, group = location, label = Freq), data = samplesizes, size = 3)

dev.off()

pdf("graphics/NCI60_repetitive_elements_boxplot_nolabels.pdf",
    width = 3.75, height = 6)

ggplot(aes(x = signal, y = group, group = location), data = RNAseq_plot_melt ) + 
  geom_boxplot(aes(group = group), fill = rep(c("grey",
                                                "darkorchid1",
                                                "grey",
                                                "magenta"
  ), times = 5), notch = TRUE, outlier.shape = NA) +
  theme_classic() +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey") +
  xlab("Linear model t-value") + 
  ylab(NULL) +
  facet_wrap(~ location, ncol = 1, scales = "free") +
  coord_cartesian(xlim = c(-4, 5.5)) + 
  theme(strip.text = element_text(size = 12, color = "white", face = "bold"),
        strip.background = element_rect(fill = "black"),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.title.x = element_text(size = 10, color = "black", face = "bold")) +
  geom_text(aes(x = 4.75, group = location, label = p_value_parse), data = dat_text, parse = TRUE, size = 3) +
  geom_text(aes(x = 4.75, group = location, label = Freq), data = samplesizes, size = 3)

dev.off()

#### EXPRESSION FROM HERVs in NCI60 for Fig S9F and Fig S9G ####

H3K9NNMT_RNAseq_repeats_signal_model <- readRDS("output/H3K9NNMT_RNAseq_HERVs_signal_model.rds")
H3K9NNMT_RNAseq_HERVs_signal_vs_t_df <- readRDS("output/H3K9NNMT_RNAseq_HERVs_signal_vs_t_df.rds")

H3K9_HERVs_expressed <- H3K9finaloverlappingcounts[apply(H3K9finaloverlappingcounts, 1, function(x){sum(x == 0) < 30}), ]

H3K9_HERVs_expressed_t <- data.frame(t(H3K9_HERVs_expressed))
H3K9_HERVs_expressed_t[, "NNMT"] <- NCI60_normalised_counts[NNMT_entrezgene_id, colnames(H3K9_HERVs_expressed)]

diseaselookup <- data.frame(disease = str_extract(unique(str_remove(str_remove(computematrix_files, "_bigwigcompare\\.bigWig_.*SAMPLE_computematrix.gz"), "^H[3-4]K[0-9]{1,2}me3_")), "^[A-Z]{2,3}"),
                            name = str_remove(unique(str_remove(str_remove(computematrix_files, "_bigwigcompare\\.bigWig_.*SAMPLE_computematrix.gz"), "^H[3-4]K[0-9]{1,2}me3_")), "^[A-Z]{2,3}_"))

H3K9_HERVs_expressed_t[, "tissue"] <- diseaselookup[match(row.names(H3K9_HERVs_expressed_t), diseaselookup$name), "disease"]

H3K9_HERVs_topbin_expressed_t <- H3K9_HERVs_expressed_t[, colnames(H3K9_HERVs_expressed_t) %in% c(str_remove(row.names(H3K9NNMT_RNAseq_HERVs_signal_vs_t_df[H3K9NNMT_RNAseq_HERVs_signal_vs_t_df$points_bin == 5, ]), "\\.t_value"), "NNMT", "tissue")]
H3K9_HERVs_topbin_expressed_number_of_sites <- ncol(H3K9_HERVs_topbin_expressed_t) - 2

H3K9_HERVs_topbin_expressed_NNMTlm_values <- lapply(colnames(H3K9_HERVs_topbin_expressed_t)[1:H3K9_HERVs_topbin_expressed_number_of_sites], function(x){

  tryCatch(expr = {
    
    templm <- lm(formula = as.formula(paste0("log10(", x, " + 1) ~ log10(NNMT + 1) + tissue + log10(NNMT + 1)*tissue")), data = H3K9_HERVs_topbin_expressed_t)
    templmsum <- summary(templm)
    
    t <- templmsum$coefficients["log10(NNMT + 1)", "t value"]
    
    output_vec <- c(t)
    names(output_vec) <- c("t_value")
    
    return(output_vec)
    
  }, error = function(e){
    
    output_vec <- c(NA)
    names(output_vec) <- c("t_value")
    
    return(output_vec)
    
  }
  )
  
})

names(H3K9_HERVs_topbin_expressed_NNMTlm_values) <- colnames(H3K9_HERVs_topbin_expressed_t)[1:H3K9_HERVs_topbin_expressed_number_of_sites]

# now for random distribution

H3K9_HERVs_topbin_expressed_t_plusrandom <- H3K9_HERVs_topbin_expressed_t

for(i in 1:length(random_expressed_RNAseq_sample)){
  
  H3K9_HERVs_topbin_expressed_t_plusrandom[, paste0(random_expressed_RNAseq_sample[i], "_RNAseq")] <- as.numeric(RNAseq_data_expressed[random_expressed_RNAseq_sample[i], match(row.names(H3K9_HERVs_topbin_expressed_t), colnames(RNAseq_data_expressed))])
      
}

H3K9_HERVs_topbin_expressed_randomlm_values <- sapply(random_expressed_RNAseq_sample[1:100], function(thisgene){

      sapply(colnames(H3K9_HERVs_topbin_expressed_t)[1:H3K9_HERVs_topbin_expressed_number_of_sites], function(x){

        tryCatch(expr = {
          
          glmtemp <- glm.nb(formula = as.formula(paste0(x, " ~ log10(", as.character(thisgene), "_RNAseq + 1) + tissue + log10(", thisgene, "_RNAseq + 1) * tissue")), data = H3K9_HERVs_topbin_expressed_t_plusrandom)
          
          glmsum <- summary(glmtemp)
          
          t <- glmsum$coefficients[paste0("log10(", thisgene, "_RNAseq + 1)"), "z value"]
          
          output_vec <- c(t)
          names(output_vec) <- c("t_value")
          
          return(output_vec)
          
        }, error = function(e){
          
          output_vec <- c(NA)
          names(output_vec) <- c("t_value")
          
          return(output_vec)
          
        }
        )
        
      })
      
    }) 

H3K9_HERVs_topbin_expressed_lm_for_plot <- reshape2::melt(data.frame("random" = rowMeans(H3K9_HERVs_topbin_expressed_randomlm_values, na.rm = TRUE),
                                                      NNMT = unlist(H3K9_HERVs_topbin_expressed_NNMTlm_values)))

colnames(H3K9_HERVs_topbin_expressed_lm_for_plot) <- c("group", "value")
H3K9_HERVs_topbin_expressed_lm_for_plot$group <- factor(H3K9_HERVs_topbin_expressed_lm_for_plot$group, levels = c("random", "NNMT"))

saveRDS(H3K9_HERVs_topbin_expressed_lm_for_plot, "plot_data/H3K9_HERVs_topbin_expressed_lm_for_plot.rds")

pdf("graphics/NCI60_HERVs_expression_vs_NNMT_zvalue_boxplot.pdf",
    width = 2,5,
    height = 2.5)

ggplot(H3K9_HERVs_topbin_expressed_lm_for_plot, aes(x = group, group = group, y = value)) + 
  geom_boxplot(fill = c("grey", "red"),
               outlier.shape = NA) + 
  theme_classic() + 
  theme(axis.text = element_text(size = 10,
                                 colour = "black"),
        axis.title.y = element_text(size = 10,
                                    colour = "black"),
        axis.title.x = element_blank()) + 
  geom_hline(yintercept = 0 , linetype = "dashed", colour = "grey") + 
  coord_cartesian(ylim = c(-4, 3)) + 
  ylab("Model z-value")

dev.off()

# H4K20 - NB results are essentially identical to NNMT (top bin in 70% identical in terms of sites)
H4K20NNMT_RNAseq_repeats_signal_model <- readRDS("output/H4K20NNMT_RNAseq_HERVs_signal_model.rds")
H4K20NNMT_RNAseq_HERVs_signal_vs_t_df <- readRDS("output/H4K20NNMT_RNAseq_HERVs_signal_vs_t_df.rds")

H4K20_HERVs_expressed <- H4K20finaloverlappingcounts[apply(H4K20finaloverlappingcounts, 1, function(x){sum(x == 0) < 30}), ]

H4K20_HERVs_expressed_t <- data.frame(t(H4K20_HERVs_expressed))
H4K20_HERVs_expressed_t[, "NNMT"] <- NCI60_normalised_counts[NNMT_entrezgene_id, colnames(H4K20_HERVs_expressed)]

diseaselookup <- data.frame(disease = str_extract(unique(str_remove(str_remove(computematrix_files, "_bigwigcompare\\.bigWig_.*SAMPLE_computematrix.gz"), "^H[3-4]K[0-9]{1,2}me3_")), "^[A-Z]{2,3}"),
                            name = str_remove(unique(str_remove(str_remove(computematrix_files, "_bigwigcompare\\.bigWig_.*SAMPLE_computematrix.gz"), "^H[3-4]K[0-9]{1,2}me3_")), "^[A-Z]{2,3}_"))

H4K20_HERVs_expressed_t[, "tissue"] <- diseaselookup[match(row.names(H4K20_HERVs_expressed_t), diseaselookup$name), "disease"]

H4K20_HERVs_topbin_expressed_t <- H4K20_HERVs_expressed_t[, colnames(H4K20_HERVs_expressed_t) %in% c(str_remove(row.names(H4K20NNMT_RNAseq_HERVs_signal_vs_t_df[H4K20NNMT_RNAseq_HERVs_signal_vs_t_df$points_bin == 5, ]), "\\.t_value"), "NNMT", "tissue")]
H4K20_HERVs_topbin_expressed_number_of_sites <- ncol(H4K20_HERVs_topbin_expressed_t) - 2

H4K20_HERVs_topbin_expressed_NNMTlm_values <- lapply(colnames(H4K20_HERVs_topbin_expressed_t)[1:H4K20_HERVs_topbin_expressed_number_of_sites], function(x){
  
  tryCatch(expr = {
    
    templm <- glm.nb(formula = as.formula(paste0(x, " ~ log10(NNMT + 1) + tissue + log10(NNMT + 1)*tissue")), data = H4K20_HERVs_topbin_expressed_t)
    templmsum <- summary(templm)
    
    t <- templmsum$coefficients["log10(NNMT + 1)", "z value"]
    
    output_vec <- c(t)
    names(output_vec) <- c("t_value")
    
    return(output_vec)
    
  }, error = function(e){
    
    output_vec <- c(NA)
    names(output_vec) <- c("t_value")
    
    return(output_vec)
    
  }
  )
  
})

names(H4K20_HERVs_topbin_expressed_NNMTlm_values) <- colnames(H4K20_HERVs_topbin_expressed_t)[1:H4K20_HERVs_topbin_expressed_number_of_sites]

# now for random distribution

H4K20_HERVs_topbin_expressed_t_plusrandom <- H4K20_HERVs_topbin_expressed_t

for(i in 1:length(random_expressed_RNAseq_sample)){
  
  H4K20_HERVs_topbin_expressed_t_plusrandom[, paste0(random_expressed_RNAseq_sample[i], "_RNAseq")] <- as.numeric(RNAseq_data_expressed[random_expressed_RNAseq_sample[i], match(row.names(H4K20_HERVs_topbin_expressed_t), colnames(RNAseq_data_expressed))])
  
}

H4K20_HERVs_topbin_expressed_randomlm_values <- sapply(random_expressed_RNAseq_sample[1:10], function(thisgene){
  
  sapply(colnames(H4K20_HERVs_topbin_expressed_t)[1:H4K20_HERVs_topbin_expressed_number_of_sites], function(x){
    
    tryCatch(expr = {
      
      glmtemp <- glm.nb(formula = as.formula(paste0(x, " ~ log10(", as.character(thisgene), "_RNAseq + 1) + tissue + log10(", thisgene, "_RNAseq + 1) * tissue")), data = H4K20_HERVs_topbin_expressed_t_plusrandom)
      
      glmsum <- summary(glmtemp)
      
      t <- glmsum$coefficients[paste0("log10(", thisgene, "_RNAseq + 1)"), "z value"]
      
      output_vec <- c(t)
      names(output_vec) <- c("t_value")
      
      return(output_vec)
      
    }, error = function(e){
      
      output_vec <- c(NA)
      names(output_vec) <- c("t_value")
      
      return(output_vec)
      
    }
    )
    
  })
  
}) 

#### expression by signal ####

H3K9expressing_samples <-apply(H3K9finaloverlappingcounts, 1, function(x){sum(x != 0)})
H3K9countsexpressed_in_half <- H3K9finaloverlappingcounts[names(H3K9expressing_samples[H3K9expressing_samples >30]), ]

H3K9fitmodels <- lapply(row.names(H3K9countsexpressed_in_half), function(x){

  tryCatch({

    model <- glm.nb(unlist(H3K9countsexpressed_in_half[x, match(str_remove(str_remove(row.names(H3K9HERVs_matrix_for_lm), "^H3K9me3_[A-Z]{2,3}_"), "_bigwigcompare\\.bigWig_.*SAMPLE_computematrix.gz"), colnames(H3K9countsexpressed_in_half))]) ~ H3K9HERVs_matrix_for_lm[, x])
    
    summ <- summary(model)
    
    summ$coefficients["H3K9HERVs_matrix_for_lm[, x]", "z value"]
  },
  
  error = function(e){
    
    return(NA)
    
  })
  
})

H4K20expressing_samples <-apply(H4K20finaloverlappingcounts, 1, function(x){sum(x != 0)})
H4K20countsexpressed_in_half <- H4K20finaloverlappingcounts[names(H4K20expressing_samples[H4K20expressing_samples >30]), ]

H4K20fitmodels <- lapply(row.names(H4K20countsexpressed_in_half), function(x){

  tryCatch({
    
    model <- glm.nb(unlist(H4K20countsexpressed_in_half[x, match(str_remove(str_remove(row.names(H4K20HERVs_matrix_for_lm), "^H4K20me3_[A-Z]{2,3}_"), "_bigwigcompare\\.bigWig_.*SAMPLE_computematrix.gz"), colnames(H4K20countsexpressed_in_half))]) ~ H4K20HERVs_matrix_for_lm[, x])
    
    summ <- summary(model)
    
    summ$coefficients["H4K20HERVs_matrix_for_lm[, x]", "z value"]
  },
  
  error = function(e){
    
    return(NA)
    
  })
  
})

HERVs_signal_vs_expression_forplot <- reshape2::melt(list(H3K9 = data.frame(value = unlist(H3K9fitmodels)),
    H4K20 = data.frame(value = unlist(H4K20fitmodels))))

saveRDS(HERVs_signal_vs_expression_forplot, "plot_data/HERVs_signal_vs_expression_forplot.rds")

pdf("graphics/NCI60_HERVexpression_byChIPsignal_boxplot.pdf",
    height = 2.5,
    width = 2.5)

ggplot(data = HERVs_signal_vs_expression_forplot, aes(y = value)) + 
  geom_boxplot(fill = "orange") + 
  facet_wrap(~ L1) + 
  theme_classic() + 
  theme(axis.text.y = element_text(size = 10,
                                 colour = "black"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 10,
                                    colour = "black"),
        axis.title.x = element_blank(),
        strip.text = element_text(size = 12, color = "white", face = "bold"),
        strip.background = element_rect(fill = "black"),
        panel.background = element_rect(fill='transparent')) + 
  geom_hline(yintercept = 0 , linetype = "dashed", colour = "grey") + 
  ylab("Expression ~ signal model z-value")

dev.off()

#### expression of H3K9 marked genes ####

presentmatches <- row.names(NCI60_normalised_counts)[(row.names(NCI60_normalised_counts) %in% str_remove(names(H3K9NNMT_RNAseq_genes_signal_model), "^X"))]

counts_for_H3K9marked <- NCI60_normalised_counts[c(presentmatches, NNMT_entrezgene_id), ]

counts_for_H3K9marked_t <- data.frame(t(counts_for_H3K9marked))

counts_for_H3K9marked_t <- counts_for_H3K9marked_t[apply(counts_for_H3K9marked_t, 2, function(x){sum(as.numeric(x) == 0) < 30})]

counts_for_H3K9marked_t[, "tissue"] <- diseaselookup[match(row.names(counts_for_H3K9marked_t), diseaselookup$name), "disease"]

H3K9marked_genes_NNMT_t_values <- sapply(colnames(counts_for_H3K9marked_t)[1:(ncol(counts_for_H3K9marked_t)-2)], function(x){

  templm <- lm(  as.formula(paste0("log10(", x, "+1) ~ log10(X", NNMT_entrezgene_id, "+1) + tissue + (tissue * log10(X", NNMT_entrezgene_id, "+1))")), data = counts_for_H3K9marked_t)

  tempsumm <- summary(templm)
  
  t <- tempsumm$coefficients[paste0("log10(X", NNMT_entrezgene_id, " + 1)"), "t value"]
  
  return(t)
  
  })

# H4K20

counts_for_H4K20marked <- NCI60_normalised_counts[row.names(NCI60_normalised_counts) %in% str_remove(names(H4K20NNMT_RNAseq_genes_signal_model), "^X"), ]

counts_for_H4K20marked_t <- data.frame(t(counts_for_H4K20marked))

counts_for_H4K20marked_t <- counts_for_H4K20marked_t[apply(counts_for_H4K20marked_t, 2, function(x){sum(as.numeric(x) == 0) < 30})]

counts_for_H4K20marked_t[, "tissue"] <- diseaselookup[match(row.names(counts_for_H4K20marked_t), diseaselookup$name), "disease"]
counts_for_H4K20marked_t[, paste0("X", NNMT_entrezgene_id)] <- NCI60_normalised_counts[NNMT_entrezgene_id, ]

H4K20marked_genes_NNMT_t_values <- sapply(colnames(counts_for_H4K20marked_t)[1:(ncol(counts_for_H4K20marked_t)-2)], function(x){
  
  templm <- lm(  as.formula(paste0("log10(", x, "+1) ~ log10(X", NNMT_entrezgene_id, "+1) + tissue + (tissue * log10(X", NNMT_entrezgene_id, "+1))")), data = counts_for_H4K20marked_t)
  
  tempsumm <- summary(templm)
  
  t <- tempsumm$coefficients[paste0("log10(X", NNMT_entrezgene_id, " + 1)"), "t value"]
  
  return(t)
  
})

#### expression vs signal for genebodies ####

H4K20genes_expression_vs_signal_models <- lapply(colnames(counts_for_H4K20marked_t)[1:(ncol(counts_for_H4K20marked_t) - 2)], function(x){

  tryCatch({

    tempsignal <- H4K20genesmatrix_for_lm[match(row.names(counts_for_H4K20marked_t), str_remove(row.names(H4K20genesmatrix_for_lm), "^H4K20me3_[A-Z]{2,3}_")), x]
    
    model <- lm(unlist(counts_for_H4K20marked_t[, x]) ~ tempsignal)
    
    summ <- summary(model)
    
    summ$coefficients["tempsignal", "t value"]
  },
  
  error = function(e){
    
    return(NA)
    
  })
  
})

H3K9genes_expression_vs_signal_models <- lapply(colnames(counts_for_H3K9marked_t)[1:(ncol(counts_for_H3K9marked_t) - 2)], function(x){

  tryCatch({
    
    tempsignal <- H3K9genesmatrix_for_lm[match(row.names(counts_for_H3K9marked_t), str_remove(row.names(H3K9genesmatrix_for_lm), "^H3K9me3_[A-Z]{2,3}_")), x]
    
    model <- lm(unlist(counts_for_H3K9marked_t[, x]) ~ tempsignal)
    
    summ <- summary(model)
    
    summ$coefficients["tempsignal", "t value"]
  },
  
  error = function(e){
    
    return(NA)
    
  })
  
})

#### NCI60 expression by signal plot ####

NCI60_expression_by_signal_plot <- list(HERVs = list(H3K9 = unlist(H3K9fitmodels),
                                                  H4K20 = unlist(H4K20fitmodels)),
                                     gene_bodies = list(H3K9 = H3K9genes_expression_vs_signal_models,
                                                        H4K20 = H4K20genes_expression_vs_signal_models)
                                     )

NCI60_expr_signal_plot_melt <- melt(NCI60_expression_by_signal_plot)
colnames(NCI60_expr_signal_plot_melt) <- c("signal", "group", "don'tknow", "location")

NCI60_expr_signal_plot_melt$group = paste0(NCI60_expr_signal_plot_melt$group, "me3")

factet_labs <- c("HERVs", "Gene bodies")
names(factet_labs) <- c("HERVs", "gene_bodies")

samplesizes <- table(paste0(NCI60_expr_signal_plot_melt$group, NCI60_expr_signal_plot_melt$location))
samplesizes <- as.data.frame(samplesizes)
samplesizes[, "group"] <- str_remove(str_remove(samplesizes$Var1, pattern = "gene_bodies"), pattern = "HERVs")

for(i in 1:nrow(samplesizes)){
  samplesizes[i, "location"] <- str_remove(samplesizes[i, "Var1"], pattern = samplesizes[i, "group"])
}

samplesizes$Freq = paste0("n = ", samplesizes$Freq)

samplesizes[, "p_value"] <- sapply(levels(samplesizes$Var1), function(x){

  wilcox.test(NCI60_expr_signal_plot_melt[paste0(NCI60_expr_signal_plot_melt$group, NCI60_expr_signal_plot_melt$location) == x, "signal"])$p.value
  
})

NCI60_expr_signal_plot_melt$group <- factor(NCI60_expr_signal_plot_melt$group, levels = c("H3K9me3", "H4K20me3"))

saveRDS(NCI60_expr_signal_plot_melt, "plot_data/NCI60_expr_signal_plot_melt.rds")
# NCI60_expr_signal_plot_melt <- readRDS("plot_data/NCI60_expr_signal_plot_melt.rds")

write.table(NCI60_expr_signal_plot_melt[, c(1:2, 4)],
            file = "plot_data/Fig S9/Fig_S9G_plotdata.txt",
            row.names = FALSE,
            col.names = TRUE,
            sep = "\t")

pdf("graphics/NCI60_expression_vs_signal_boxplot.pdf",
    width = 4, height = 2.8)

ggplot(aes(y = signal, x = group, group = location), data = NCI60_expr_signal_plot_melt ) + 
  geom_boxplot(aes(group = group),
               notch = TRUE,
               outlier.shape = NA,
               fill = c("magenta",
                        "dark orchid",
                        "magenta",
                        "dark orchid")) +
  facet_wrap(~ location, labeller = labeller(location = factet_labs)) +
  theme_classic() +
  theme(strip.text = element_text(size = 12, color = "white", face = "bold"),
                               strip.background = element_rect(fill = "black"),
                               axis.text.y = element_text(size = 10, color = "black"),
                               axis.text.x = element_text(size = 10, color = "black"),
        axis.title.x = element_blank(),
                                panel.border = element_rect(colour = "black",
                                                            fill = "transparent"),
                               axis.title.y = element_text(size = 9, color = "black")) +
  geom_hline(yintercept = 0,
            linetype = "dashed",
            colour = "grey") + 
  geom_text(aes(y = 4.75, group = location, label = Freq),
            data = samplesizes,
            size = 3) + 
  coord_cartesian(ylim = c(-5, 5)) + 
  ylab("Expression ~ ChIP signal linear model t-value")
  
dev.off()

#### expression vs NNMT for gene bodies ####

# since is very similar / partially redundant between marks, will do only for H3K9.

H3K9_40_of_60_modelled_entrez <- str_remove(colnames(H3K9genesmatrix_for_lm)[1:(ncol(H3K9genesmatrix_for_lm)-5)], "^X")

H3K9_40_of_60_counts <- NCI60_normalised_counts[H3K9_40_of_60_modelled_entrez, ]

# need to weed out ones with almost no expression
H3K9_40_of_60_counts <- H3K9_40_of_60_counts[apply(H3K9_40_of_60_counts, 1, function(X){sum(X == 0) < 15}), ]

H3K9_40_of_60_counts_t_df <- data.frame(t(H3K9_40_of_60_counts))

H3K9_40_of_60_counts_t_df[, "NNMT_RNAseq"] <- NCI60_normalised_counts[NNMT_entrezgene_id, row.names(H3K9_40_of_60_counts_t_df)]
H3K9_40_of_60_counts_t_df[, "tissue"] <- diseaselookup[match(row.names(H3K9_40_of_60_counts_t_df), diseaselookup$name), "disease"]

H3K9genes_expression_vs_NNMT_models <- lapply(colnames(H3K9_40_of_60_counts_t_df)[1:(ncol(H3K9_40_of_60_counts_t_df) - 2)], function(x){
  
  tryCatch({

    model <- lm(as.formula(paste0("log10(", x, " + 1) ~ log10(NNMT_RNAseq + 1) + tissue + (tissue * log10(NNMT_RNAseq + 1))")), data = H3K9_40_of_60_counts_t_df)
    
    summ <- summary(model)
    
    summ$coefficients["log10(NNMT_RNAseq + 1)", "t value"]
  },
  
  error = function(e){
    
    return(NA)
    
  })
  
})

names(H3K9genes_expression_vs_NNMT_models) <- colnames(H3K9_40_of_60_counts_t_df)[1:(ncol(H3K9_40_of_60_counts_t_df) - 2)]

# H4K20 too

# since is very similar / partially redundant between marks, will do only for H4K20.

H4K20_40_of_60_modelled_entrez <- str_remove(colnames(H4K20genesmatrix_for_lm)[1:(ncol(H4K20genesmatrix_for_lm)-5)], "^X")

H4K20_40_of_60_counts <- NCI60_normalised_counts[H4K20_40_of_60_modelled_entrez, ]

# need to weed out ones with almost no expression
H4K20_40_of_60_counts <- H4K20_40_of_60_counts[apply(H4K20_40_of_60_counts, 1, function(X){sum(X == 0) < 15}), ]

H4K20_40_of_60_counts_t_df <- data.frame(t(H4K20_40_of_60_counts))

H4K20_40_of_60_counts_t_df[, "NNMT_RNAseq"] <- NCI60_normalised_counts[NNMT_entrezgene_id, row.names(H4K20_40_of_60_counts_t_df)]
H4K20_40_of_60_counts_t_df[, "tissue"] <- diseaselookup[match(row.names(H4K20_40_of_60_counts_t_df), diseaselookup$name), "disease"]

H4K20genes_expression_vs_NNMT_models <- lapply(colnames(H4K20_40_of_60_counts_t_df)[1:(ncol(H4K20_40_of_60_counts_t_df) - 2)], function(x){
  
  tryCatch({
    
    model <- lm(as.formula(paste0("log10(", x, " + 1) ~ log10(NNMT_RNAseq + 1) + tissue + (tissue * log10(NNMT_RNAseq + 1))")), data = H4K20_40_of_60_counts_t_df)
    
    summ <- summary(model)
    
    summ$coefficients["log10(NNMT_RNAseq + 1)", "t value"]
  },
  
  error = function(e){
    
    return(NA)
    
  })
  
})

names(H4K20genes_expression_vs_NNMT_models) <- colnames(H4K20_40_of_60_counts_t_df)[1:(ncol(H4K20_40_of_60_counts_t_df) - 2)]

#### NCI60 expression by NNMT plot ####

NCI60_H4K20expression_by_NNMT_plot <- list(HERVs = list(signal = unlist(H4K20NNMT_RNAseq_HERVs_signal_model[names(H4K20_HERVs_topbin_expressed_NNMTlm_values)]),
                                                     expression = unlist(H4K20_HERVs_topbin_expressed_NNMTlm_values)),
                                        gene_bodies = list(signal = unlist(H4K20NNMT_RNAseq_genes_signal_model[names(H4K20genes_expression_vs_NNMT_models)]),
                                                           expression = unlist(H4K20genes_expression_vs_NNMT_models))
)

NCI60_H4K20expression_by_NNMT_plot_melt <- reshape2::melt(NCI60_H4K20expression_by_NNMT_plot)
colnames(NCI60_H4K20expression_by_NNMT_plot_melt) <- c("tvals", "group", "location")

NCI60_H4K20expression_by_NNMT_plot_melt[NCI60_H4K20expression_by_NNMT_plot_melt$group == "signal", "group"] <- "H4K20me3\nsignal"
NCI60_H4K20expression_by_NNMT_plot_melt[NCI60_H4K20expression_by_NNMT_plot_melt$group == "expression", "group"] <- "Gene\nexpression"

factet_labs <- c("HERVs", "Gene bodies")
names(factet_labs) <- c("HERVs", "gene_bodies")

NCI60_H4K20expression_by_NNMT_plot_melt$group <- factor(NCI60_H4K20expression_by_NNMT_plot_melt$group, levels = c("H4K20me3\nsignal", "Gene\nexpression"))

table(NCI60_H4K20expression_by_NNMT_plot_melt$group)
samplesizes <- table(paste0(NCI60_H4K20expression_by_NNMT_plot_melt$group, NCI60_H4K20expression_by_NNMT_plot_melt$location))

samplesizes_edit <- c("gene_bodies" = max(samplesizes[str_detect(names(samplesizes), "gene_bodies")]),
                      "HERVs" = max(samplesizes[str_detect(names(samplesizes), "HERVs")]))

samplesizes_edit <- data.frame(cbind(samplesizes_edit,
                 names(samplesizes_edit)))

names(samplesizes_edit) <- c("Freq", "location")

samplesizes_edit$Freq = paste0("n = ", samplesizes_edit$Freq)

saveRDS(NCI60_H4K20expression_by_NNMT_plot_melt, "plot_data/NCI60_H4K20expression_by_NNMT_plot_melt.rds")
# NCI60_H4K20expression_by_NNMT_plot_melt <- readRDS("plot_data/NCI60_H4K20expression_by_NNMT_plot_melt.rds")

write.table(NCI60_H4K20expression_by_NNMT_plot_melt,
            file = "plot_data/Fig S9/Fig_S9F_plotdata.txt",
            row.names = FALSE,
            col.names = TRUE,
            sep = "\t")

pdf("graphics/NCI60_H4K20me3_expression_by_NNMT_boxplot.pdf",
    width = 4, height = 3)

ggplot(aes(y = tvals, x = group, group = location), data = NCI60_H4K20expression_by_NNMT_plot_melt) + 
  geom_boxplot(aes(group = group),
               notch = TRUE,
               outlier.shape = NA,
               fill = c("darkorchid1",
                        "red",
                        "darkorchid1",
                        "red")) +
  facet_wrap(~ location, labeller = labeller(location = factet_labs)) +
  theme_classic() +
  theme(strip.text = element_text(size = 12, color = "white", face = "bold"),
        strip.background = element_rect(fill = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.text.x = element_text(size = 10, color = "black"),
        panel.border = element_rect(colour = "black",
                                    fill = "transparent"),
        axis.title.y = element_text(size = 9, color = "black"),
        axis.title.x = element_blank()) +
  geom_hline(yintercept = 0,
             linetype = "dashed",
             colour = "grey") +
  ylab(italic(NNMT)~" Linear model t-value") + 
  geom_text(aes(x = 1.5, y = 3, group = location, label = Freq), data = samplesizes_edit) + 
  coord_cartesian(ylim = c(-4, 3.5)) +
  scale_y_continuous(breaks = c(-3, -1.5, 0, 1.5, 3))

dev.off()

#### CeNDR C. elegans data for Fig 4 and Fig S8 ####

#### process raw data ####

## configure biomaRt to retrieve gene info

# new_config <- httr::config(ssl_verifypeer = FALSE)
# httr::set_config(new_config, override = FALSE)
parasite_mart <- useMart("parasite_mart", dataset = "wbps_gene", host = "https://parasite.wormbase.org", port = 443)

CeNDR_TPMs[, "gene"] <- str_extract(CeNDR_TPMs$transcript, "^[0-9A-Z]+\\.[0-9]{1,2}")

CeNDR_TPMs[is.na(CeNDR_TPMs$gene), "gene"] <- str_extract(CeNDR_TPMs[is.na(CeNDR_TPMs$gene), "transcript"], "^[0-9A-Z]+\\.[a-z]{1,2}")
CeNDR_TPMs[is.na(CeNDR_TPMs$gene), "gene"] <- CeNDR_TPMs[is.na(CeNDR_TPMs$gene), "transcript"]
CeNDR_TPMs <- CeNDR_TPMs[, 2:ncol(CeNDR_TPMs)]

# collapse to gene level by summing transcript counts
CeNDR_TPMs_collapse <- aggregate(. ~ gene, data = CeNDR_TPMs, FUN = sum)
row.names(CeNDR_TPMs_collapse) <- CeNDR_TPMs_collapse$gene
CeNDR_TPMs_collapse <- CeNDR_TPMs_collapse[, 2:ncol(CeNDR_TPMs_collapse)]

saveRDS(CeNDR_TPMs_collapse, "output/Cel_TPM_gene.rds")
# CeNDR_TPMs_collapse <- readRDS("output/Cel_TPM_gene.rds")

CeNDR_TPMgenes <- getBM(mart = parasite_mart,
                  filters = "wormbase_gseqname",
                  value = row.names(CeNDR_TPMs_collapse),
                  attributes = c("wormbase_gene", "wormbase_gseq", "description", "wormbase_locus"))

TPMcollapse_names <- CeNDR_TPMgenes[match(row.names(CeNDR_TPMs_collapse), CeNDR_TPMgenes$wormbase_gseq), "wormbase_gene"]

CeNDR_TPMs_collapse2 <- CeNDR_TPMs_collapse[!is.na(TPMcollapse_names), ]
row.names(CeNDR_TPMs_collapse2) <- TPMcollapse_names[!is.na(TPMcollapse_names)]

# remove genes with zero expression across all samples
CeNDR_TPMs_collapse_expressed <- CeNDR_TPMs_collapse2[apply(CeNDR_TPMs_collapse2, 1, function(x){any(x != 0)}), ]

# transform for later use with RAPToR as per RAPToR vignette
CeNDR_TPMnorm <- limma::normalizeBetweenArrays(CeNDR_TPMs_collapse_expressed, method = "quantile")
CeNDR_TPMnormlog <- log1p(CeNDR_TPMnorm) # log1p(x) = log(x + 1)

# now MRN normalisation of raw counts

# convert raw counts into MRN pseudocounts 

# DESeq2 wants a colData object. Not actually used for the normalisation. Here we can use the sample IDs with tissue.
CeNDR_col_data <- str_extract(colnames(CeNDR_raw_counts)[2:ncol(CeNDR_raw_counts)], "^[0-9A-Z]+")
CeNDR_col_data <- as.matrix(CeNDR_col_data)

rownames(CeNDR_col_data) <- colnames(CeNDR_raw_counts)[2:ncol(CeNDR_raw_counts)]
colnames(CeNDR_col_data) <- c("Strain")

# convert to matrix of integers

CeNDR_raw_counts_mat <- as.matrix(CeNDR_raw_counts[, 2:ncol(CeNDR_raw_counts)])
colnames(CeNDR_raw_counts_mat) <- colnames(CeNDR_raw_counts)[2:ncol(CeNDR_raw_counts)]

CeNDR_raw_counts_mat <- apply(CeNDR_raw_counts_mat, 2, as.integer)
row.names(CeNDR_raw_counts_mat) <- CeNDR_raw_counts[, 1]

# need to create a DESeq2 object. Design set to ~1 allows for use of estimateSizeFactors
tempdds <- DESeqDataSetFromMatrix(countData = CeNDR_raw_counts_mat, colData = CeNDR_col_data, design = ~ 1)

# this function estimates the scaling factors from the samples from the median of ratios wrt to the geometric mean for each gene across samples
tempdds <- estimateSizeFactors(tempdds)

# put the counts normalised by the scaling factors in a new object
CeNDR_normalised_counts <- counts(tempdds, normalized = TRUE)

row.names(CeNDR_normalised_counts) <- row.names(CeNDR_raw_counts_mat)

Cel_genes_extract <- str_extract(row.names(CeNDR_normalised_counts), "^[0-9A-Z]+\\.[0-9]{1,2}")

Cel_CeNDR_normalised_counts <- data.frame(CeNDR_normalised_counts)
Cel_CeNDR_normalised_counts[, "gene"] <- str_extract(row.names(Cel_CeNDR_normalised_counts), "^[0-9A-Z]+\\.[0-9]{1,2}")

# collapse to gene level by summing transcript counts
CeNDR_normalised_counts_collapse <- aggregate(. ~ gene, data = Cel_CeNDR_normalised_counts, FUN = sum)

row.names(CeNDR_normalised_counts_collapse) <- CeNDR_normalised_counts_collapse$gene
CeNDR_normalised_counts_collapse <- CeNDR_normalised_counts_collapse[, 2:ncol(CeNDR_normalised_counts_collapse)]

# remove genes which are not expressed in any sample
CeNDR_normalised_counts_collapse <- CeNDR_normalised_counts_collapse[apply(CeNDR_normalised_counts_collapse, 1, function(x){any(x != 0)}), ]

saveRDS(CeNDR_normalised_counts_collapse, "output/Cel_CeNDR_normalised_counts_gene_level.rds")
# CeNDR_normalised_counts_collapse <- readRDS("output/Cel_CeNDR_normalised_counts_gene_level.rds")

#### define histone methyltransferases ####

worm_hmt_manual <- c("mes-3",     #  H3K27
                     "mes-2",     #  H3K27
                     "mes-4",     #  H3K36
                     "mes-6",     #  H3K27
                     "met-2",     #  H3K9
                     "set-25",    #  H3K9
                     "set-1",     #  predicted
                     "set-2",     #  H3K4
                     "set-17",    #  H3K4
                     "set-16",    #  H3K4
                     "dot-1.5",     #  H3K79
                     "set-27",      #  H3K36
                     "abu-12",      #  H3K79
                     "dot-1.3",     #  H3K79
                     "dot-1.2",     #  H3K79
                     "dot-1.4",     #  H3K79
                     "dot-1.1",     #  H3K79
                     "set-11",      #  ?
                     "set-6",       #  H3K9
                     "set-18",      #  H3K36
                     "set-30",      #  H3K4
                     "set-4",       #  H4K20
                     "set-12",      #  H3K36
                     "set-23",      #  ?
                     "met-1",       # H3K36
                     "set-26",      # H3K9
                     "set-10",      # ?
                     "set-14",
                     "lin-59")

worm_targets <-  c("H3K27",
                   "H3K27",
                   "H3K36",
                   "H3K27",
                   "H3K9",
                   "H3K9",
                   "unknown",     #  predicted
                   "H3K4",     #  H3K4
                   "H3K4",    #  H3K4
                   "H3K4",    #  H3K4
                   "H3K79",     #  H3K79
                   "H3K36",      #  H3K36
                   "H3K79",      #  H3K79
                   "H3K79",     #  H3K79
                   "H3K79",     #  H3K79
                   "H3K79",     #  H3K79
                   "H3K79",     #  H3K79
                   "unknown",      #  ?
                   "H3K9",       #  H3K9
                   "H3K36",      #  H3K36
                   "H3K4",      #  H3K4
                   "H4K20",       #  H4K20
                   "H3K36",      #  H3K36
                   "unknown",      #  ?
                   "H3K36",       # H3K36
                   "H3K9",      # H3K9
                   "unknown",      # ?
                   "unknown",
                   "H3K36")

Cel_HMTgenes <- getBM(mart = parasite_mart, 
               filters = "gene_name",
               value = worm_hmt_manual,
               attributes = c("wormbase_gene",
                              "wormbase_locus",
                              "wormbase_gseq",
                              "description",
                              "start_position",
                              "end_position",
                              "chromosome_name",
                              "strand"))

# eliminate random blank row
Cel_HMTgenes <- Cel_HMTgenes[apply(Cel_HMTgenes, 1, function(x){any(x != "")}), ]
Cel_HMTgenes <- Cel_HMTgenes[!(Cel_HMTgenes$wormbase_gene == ""), ]

Cel_HMTgenes <- Cel_HMTgenes[match(worm_hmt_manual, Cel_HMTgenes$wormbase_locus), ]
Cel_HMTgenes[, "target"] <- worm_targets

write.xlsx(Cel_HMTgenes, "output/Cel_HMTgene.xlsx")

# get positions

View(listAttributes(parasite_mart))

CelHMT_forGR <- Cel_HMTgenes[, c("chromosome_name",
                                 "start_position",
                                 "end_position",
                                 "strand",
                                 "wormbase_locus",
                                 "wormbase_gseq")]

colnames(CelHMT_forGR) <- c("chr", "start", "end", "strand", "wormbase_locus", "wormbase_gseq")
CelHMT_forGR[CelHMT_forGR$strand == 1, "strand"] <- "+"
CelHMT_forGR[CelHMT_forGR$strand == -1, "strand"] <- "-"

CelHMT_forGR$chr <- paste0("chr", CelHMT_forGR$chr) 

CelHMTGR <-  makeGRangesFromDataFrame(CelHMT_forGR, keep.extra.columns = TRUE)

#### anmt and pmt expression versus HMTs no age correction; Fig 11B

# let's try collapsing the expression for the strains

CeNDR_strains <- unique(str_extract(colnames(CeNDR_normalised_counts_collapse), "^[0-9A-Z]+"))

CeNDR_strain_HMTs <- sapply(CeNDR_strains, function(x){

  normalised_counts_this_strain <- CeNDR_normalised_counts_collapse[, str_detect(colnames(CeNDR_normalised_counts_collapse), x)]
  
  if(is.vector(normalised_counts_this_strain)){
    
    names(normalised_counts_this_strain) <- row.names(CeNDR_normalised_counts_collapse)
    
    HMT_reads_t_this_strain <- sapply(Cel_HMTgenes$wormbase_gseq, function(x){
      
      thisgenes <- normalised_counts_this_strain[str_detect(names(normalised_counts_this_strain), x)]
      
    })
    
    return(sapply(HMT_reads_t_this_strain, sum))
    
  } else {
    
    HMT_reads_t_this_strain <- sapply(Cel_HMTgenes$wormbase_gseq, function(x){
      
      thisgenes <- normalised_counts_this_strain[str_detect(row.names(normalised_counts_this_strain), x), ] 
      
      if(is.vector(thisgenes)){
        
        return(thisgenes)
        
      } else {
        
        return(colSums(thisgenes))
        
      }
      
    })
    
  }
  
  HMT_reads_this_strain <- t(HMT_reads_t_this_strain)
  colnames(HMT_reads_this_strain) <- row.names(HMT_reads_t_this_strain)
  row.names(HMT_reads_this_strain) <- colnames(HMT_reads_t_this_strain)
  
  return(apply(HMT_reads_this_strain, 1, function(x){gm_mean(x)}))
  
})

CeNDR_strain_HMTsums <- colSums(CeNDR_strain_HMTs)

CeNDR_strain_anmt1 <- sapply(CeNDR_strains, function(x){

  normalised_counts_this_strain <- CeNDR_normalised_counts_collapse[, str_detect(colnames(CeNDR_normalised_counts_collapse), x)]
  
  if(is.vector(normalised_counts_this_strain)){
    
    names(normalised_counts_this_strain) <- row.names(CeNDR_normalised_counts_collapse)
    
    thisgenes <- normalised_counts_this_strain[str_detect(names(normalised_counts_this_strain), "B0303.2")] 
    
    thisgene_sum <- sum(thisgenes)
    
  } else {
    
    thisgenes <- normalised_counts_this_strain[str_detect(row.names(normalised_counts_this_strain), "B0303.2"), ]
    
    thisgene_sum <- colSums(thisgenes)
    
  }
  
  return(gm_mean(thisgene_sum))
    
  })

CeNDR_strain_anmt3 <- sapply(CeNDR_strains, function(x){
  
  normalised_counts_this_strain <- CeNDR_normalised_counts_collapse[, str_detect(colnames(CeNDR_normalised_counts_collapse), x)]
  
  if(is.vector(normalised_counts_this_strain)){
    
    names(normalised_counts_this_strain) <- row.names(CeNDR_normalised_counts_collapse)
    
    thisgenes <- normalised_counts_this_strain[str_detect(names(normalised_counts_this_strain), "T07C12.9")] 
    
    thisgene_sum <- sum(thisgenes)
    
  } else {
    
    thisgenes <- normalised_counts_this_strain[str_detect(row.names(normalised_counts_this_strain), "T07C12.9"), ]
    
    thisgene_sum <- colSums(thisgenes)
    
  }
  
  return(gm_mean(thisgene_sum))
  
})
  
CeNDR_strain_pmt1 <- sapply(CeNDR_strains, function(x){
  
  normalised_counts_this_strain <- CeNDR_normalised_counts_collapse[, str_detect(colnames(CeNDR_normalised_counts_collapse), x)]
  
  if(is.vector(normalised_counts_this_strain)){
    
    names(normalised_counts_this_strain) <- row.names(CeNDR_normalised_counts_collapse)
    
    thisgenes <- normalised_counts_this_strain[str_detect(names(normalised_counts_this_strain), "ZK622.3")] 
    
    thisgene_sum <- sum(thisgenes)
    
  } else {
    
    thisgenes <- normalised_counts_this_strain[str_detect(row.names(normalised_counts_this_strain), "ZK622.3"), ]
    
    thisgene_sum <- colSums(thisgenes)
    
  }
  
  return(gm_mean(thisgene_sum))
  
})

# stats for anmt plot
format(cor.test(log10(CeNDR_strain_anmt1), log10(CeNDR_strain_HMTsums))$p.value, scientific = T)
format(cor.test(log10(CeNDR_strain_anmt3), log10(CeNDR_strain_HMTsums))$p.value, scientific = T)

strain_uncorrected_for_plots <- data.frame(ANMT1 = CeNDR_strain_anmt1,
                                           ANMT3 = CeNDR_strain_anmt3,
                                           PMT1 = CeNDR_strain_pmt1,
                                           HMT = CeNDR_strain_HMTsums)

saveRDS(strain_uncorrected_for_plots, "plot_data/strain_uncorrected_for_plots.rds")
# strain_uncorrected_for_plots <- readRDS("plot_data/strain_uncorrected_for_plots.rds")

write.table(strain_uncorrected_for_plots,
            file = "plot_data/Fig S11/Fig_S11AB_plotdata.txt",
            col.names = TRUE,
            row.names = TRUE,
            sep = "\t")

pdf("~/NNMT_manuscript/graphics/CeNDR_strains_notcorrected_HMTvsANMT1_scatter.pdf", 
    height = 1.8,
    width = 1.8)

ggplot(data = strain_uncorrected_for_plots, aes(x = HMT, y = ANMT1)) +
  geom_point() +
  scale_y_log10(breaks = c(25, 50, 100, 250),
                expand = c(0, 0),
                limits = c(25, 250)) +
  scale_x_log10(breaks = c(30000, 40000, 50000),
                expand = c(0, 0),
                limits = c(25000, 55000)) +
  theme_classic() + 
  ylab(substitute("Strain"~italic(anmt-1)~"expression")) + 
  xlab(substitute("Strain mean total\n HMT expression (pseudocounts)")) + 
  theme(axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  geom_smooth(method = "lm",
              se = FALSE) +
annotate(geom = "text", 
         x = 29500, 
         y = 36, 
         label = substitute(italic(r)~"= -0.715"),
         size = 2.5) +
  annotate(geom = "text",
           x = 31500,
           y = 30,
           label = substitute(italic(p)~"= 6.94 x"~10^-34),
           size = 2.5)

dev.off()

pdf("~/NNMT_manuscript/graphics/CeNDR_strains_notcorrected_HMTvsANMT3_scatter.pdf", 
    height = 1.8,
    width = 1.8)

ggplot(data = strain_uncorrected_for_plots, aes(x = HMT, y = ANMT3)) +
  geom_point() +
  scale_y_log10(breaks = c(150, 500, 1500),
                expand = c(0, 0),
                limits = c(100, 1800)) +
  scale_x_log10(breaks = c(30000, 40000, 50000),
                expand = c(0, 0),
                limits = c(25000, 55000)) +
  theme_classic() + 
  ylab(substitute("Strain"~italic(anmt-3)~"expression")) + 
  xlab(substitute("Strain mean total\n HMT expression (pseudocounts)")) + 
  theme(axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  geom_smooth(method = "lm",
              se = FALSE) +
  annotate(geom = "text", 
           x = 29500, 
           y = 155, 
           label = substitute(italic(r)~"= -0.732"),
           size = 2.5) +
  annotate(geom = "text",
           x = 31500,
           y = 125,
           label = substitute(italic(p)~"= 3.62 x"~10^-36),
           size = 2.5)

dev.off()

# stats for pmt1 plot
format(cor.test(log10(CeNDR_strain_pmt1), log10(CeNDR_strain_HMTsums))$p.value, scientific = T)

pdf("~/NNMT_manuscript/graphics/CeNDR_strains_notcorrected_HMTvsPMT1_scatter.pdf", 
    height = 1.8,
    width = 1.8)

ggplot(data = strain_uncorrected_for_plots, aes(x = HMT, y = PMT1)) +
  geom_point() +
  scale_y_log10(breaks = c(2500, 5000, 12500),
                expand = c(0, 0),
                limits = c(2500, 12500)) +
  scale_x_log10(breaks = c(30000, 40000, 50000),
                expand = c(0, 0),
                limits = c(25000, 55000)) +
  theme_classic() + 
  ylab(substitute("Strain"~italic(pmt-1)~"expression")) + 
  xlab(substitute("Strain mean total\n HMT expression (pseudocounts)")) + 
  theme(axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  geom_smooth(method = "lm",
              se = FALSE) +
  annotate(geom = "text",
           x = 30000,
           y = 3300,
           label = substitute(italic(r)~"= -0.669"),
           size = 2.5) +
  annotate(geom = "text",
           x = 31500,
           y = 2800,
           label = substitute(italic(p)~"= 2.43 x"~10^-28),
           size = 2.5)

dev.off()

#### age estimates with RAPToR for Fig S11C ####

r_YA2 <- prepare_refdata("Cel_YA_2", 'wormRef', 5000)

andersen_ae2 <- RAPToR::ae(samp = CeNDR_TPMnormlog,                         # input gene expression matrix
                           refdata = r_YA2$interpGE,            # reference gene expression matrix
                           ref.time_series = r_YA2$time.series)

saveRDS(andersen_ae2, "output/andersen_ageestimates_CelYA2.rds")
# andersen_ae2 <- readRDS("output/andersen_ageestimates_CelYA2.rds")

# make df with replicate age and expressions
HMTsums_age2 <- data.frame(HMT = colSums(CeNDR_normalised_counts_collapse[Cel_HMTgenes$wormbase_gseq, ]), 
                           age = andersen_ae2$age.estimates[colnames(CeNDR_normalised_counts_collapse), "age.estimate"],
                           ANMT = unlist(CeNDR_normalised_counts_collapse["B0303.2", ]),
                           ANMT3 = unlist(CeNDR_normalised_counts_collapse["T07C12.9", ]),
                           PMT1 = unlist(CeNDR_normalised_counts_collapse["ZK622.3", ]),
                           batch = str_extract(colnames(CeNDR_normalised_counts_collapse), "[0-9]{2}_[0-9]{4}$"))

# produce plots with splines.
# later want to do analysis with spline residuals (instead of lm)

HMT_spline_df6 <- smooth.spline(x = HMTsums_age2$age,y = log10(HMTsums_age2$HMT), df = 6)
HMTsums_age2[, "HMTsplinefit_df6"] <- predict(HMT_spline_df6, x = HMTsums_age2$age)$y
HMTsums_age2[, "HMTsplineresid_df6"] <- residuals(HMT_spline_df6)

ANMT_spline_df6 <- smooth.spline(x = HMTsums_age2$age,y = log10(HMTsums_age2$ANMT), df = 6)
HMTsums_age2[, "ANMTsplinefit_df6"] <- predict(ANMT_spline_df6, x = HMTsums_age2$age)$y
HMTsums_age2[, "ANMTsplineresid_df6"] <- residuals(ANMT_spline_df6)


ANMT3_spline_df6 <- smooth.spline(x = HMTsums_age2$age,y = log10(HMTsums_age2$ANMT3), df = 6)
HMTsums_age2[, "ANMT3splinefit_df6"] <- predict(ANMT3_spline_df6, x = HMTsums_age2$age)$y
HMTsums_age2[, "ANMT3splineresid_df6"] <- residuals(ANMT3_spline_df6)

PMT1_spline_df6 <- smooth.spline(x = HMTsums_age2$age,y = log10(HMTsums_age2$PMT1), df = 6)
HMTsums_age2[, "PMT1splinefit_df6"] <- predict(PMT1_spline_df6, x = HMTsums_age2$age)$y
HMTsums_age2[, "PMT1splineresid_df6"] <- residuals(PMT1_spline_df6)

saveRDS(HMTsums_age2, "plot_data/HMTsums_age2.rds")

write.table(HMTsums_age2[, 1:5],
            file = "plot_data/Fig S11/Fig_S11C_plotdata.txt",
            col.names = TRUE,
            row.names = TRUE,
            sep = "\t")

# plots with splines

pdf("graphics/HMTvsage_foroverlay.pdf",
    width = 5,
    height = 3.5)

ggplot(data = HMTsums_age2, aes(x = age, y = log10(HMT))) +
  geom_point(colour = "magenta",
             alpha = 0.3) +
  geom_line(aes(x = age, y = HMTsplinefit_df6),
            colour = "magenta", 
            linewidth = 2) + 
  theme_classic() +
  ylab("") + 
  xlab("") + 
  scale_x_continuous(limits = c(58, 73),
                     breaks = c(58, 61, 64, 67, 70, 73)) + 
  theme(axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "magenta"))

dev.off()

pdf("graphics/ANMTvsage_foroverlay.pdf",
    width = 5,
    height = 3.5)

ggplot(data = HMTsums_age2, aes(x = age, y = log10(ANMT))) +
  geom_point(colour = "green",
             alpha = 0.3) +
  geom_line(aes(x = age, y = ANMTsplinefit_df6),
            colour = "green", 
            linewidth = 2) + 
  theme_classic() +
  ylab("") + 
  xlab(" \n") + 
  scale_x_continuous(limits = c(58, 73),
                     breaks = c(58, 61, 64, 67, 70, 73)) + 
  scale_y_continuous(position = "right",
                     breaks = c(1.4, 1.9, 2.4)) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y.right = element_text(colour = "green"),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_rect(fill='transparent'),
        legend.box.background = element_rect(fill='transparent'))

dev.off()

pdf("graphics/PMT1vsage_foroverlay.pdf",
    width = 5,
    height = 3.5)

ggplot(data = HMTsums_age2, aes(x = age, y = log10(PMT1))) +
  geom_point(colour = "orange",
             alpha = 0.3) +
  geom_line(aes(x = age, y = PMT1splinefit_df6),
            colour = "orange",
            linewidth = 2) + 
  theme_classic() +
  ylab("") + 
  xlab(" \n") + 
  scale_x_continuous(limits = c(58, 73),
                     breaks = c(58, 61, 64, 67, 70, 73)) + 
  scale_y_continuous(position = "right",
                     breaks = c(3.5, 3.9)) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y.right = element_text(colour = "orange"),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_rect(fill='transparent'),
        legend.box.background = element_rect(fill='transparent'))

dev.off()

pdf("graphics/ANMT3vsage_foroverlay.pdf",
    width = 5,
    height = 3.5)

ggplot(data = HMTsums_age2, aes(x = age, y = log10(ANMT3))) +
  geom_point(colour = "green",
             alpha = 0.3) +
  geom_line(aes(x = age, y = ANMT3splinefit_df6),
            colour = "green", 
            linewidth = 2) + 
  theme_classic() +
  ylab("") + 
  xlab(" \n") + 
  scale_x_continuous(limits = c(58, 73),
                     breaks = c(58, 61, 64, 67, 70, 73)) + 
  scale_y_continuous(position = "right",
                     breaks = c(2.1, 2.7, 3.3)) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y.right = element_text(colour = "green"),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_rect(fill='transparent'),
        legend.box.background = element_rect(fill='transparent'))

dev.off()

# try correcting by batch

HMTsums_age2[, "HMTsplineresid_batch_df6"] <- residuals(lmer(HMTsplineresid_df6 ~ (1|batch), data = HMTsums_age2))
HMTsums_age2[, "ANMTsplineresid_batch_df6"] <- residuals(lmer(ANMTsplineresid_df6 ~ (1|batch), data = HMTsums_age2))
HMTsums_age2[, "ANMT3splineresid_batch_df6"] <- residuals(lmer(ANMT3splineresid_df6 ~ (1|batch), data = HMTsums_age2))
HMTsums_age2[, "PMT1splineresid_batch_df6"] <- residuals(lmer(PMT1splineresid_df6 ~ (1|batch), data = HMTsums_age2))

HMTsums_age2[, "strain"] <- str_extract(row.names(HMTsums_age2), "^[0-9A-Z]+")

HMTsums_strain_correctedresid <- aggregate(cbind(HMTsplineresid_batch_df6,
                                                 ANMTsplineresid_batch_df6,
                                                 ANMT3splineresid_batch_df6,
                                                 PMT1splineresid_batch_df6) ~ strain,
                                           data = HMTsums_age2,
                                           mean)

cor.test(HMTsums_strain_correctedresid$HMTsplineresid_batch_df6, HMTsums_strain_correctedresid$ANMTsplineresid_batch_df6)
cor.test(HMTsums_strain_correctedresid$HMTsplineresid_batch_df6, HMTsums_strain_correctedresid$ANMTsplineresid_batch_df6, method = "spearman")

cor.test(HMTsums_strain_correctedresid$HMTsplineresid_batch_df6, HMTsums_strain_correctedresid$ANMT3splineresid_batch_df6)
cor.test(HMTsums_strain_correctedresid$HMTsplineresid_batch_df6, HMTsums_strain_correctedresid$ANMT3splineresid_batch_df6, method = "spearman")

cor.test(HMTsums_strain_correctedresid$HMTsplineresid_batch_df6, HMTsums_strain_correctedresid$PMT1splineresid_batch_df6)
cor.test(HMTsums_strain_correctedresid$HMTsplineresid_batch_df6, HMTsums_strain_correctedresid$PMT1splineresid_batch_df6, method = "spearman")

# with splines, whole genome
base_lmdf <- data.frame(batch = HMTsums_age2$batch,
                        strain = HMTsums_age2$strain,
                        age = HMTsums_age2$age,
                        sample = row.names(HMTsums_age2))

allHMTANMTPMTcors <- apply(CeNDR_normalised_counts_collapse[!(row.names(CeNDR_normalised_counts_collapse) %in% Cel_HMTgenes$wormbase_gseq), ], 1, function(x){

  # so for each gene will need to add expression to base table featuring age and batch (or better yet add age and batch to each table with expression)
  # collapse down.
  # then do the linear model for each to regress out age and batch
  # then take the average across the strains
  # then correlate that to HMTs. 
  
  base_lmdf[, "gene_noise"] <- x + rnorm(x, 1, 0.1)
  base_lmdf[, "gene_splineresid"] <- residuals(smooth.spline(base_lmdf$age, log10(base_lmdf$gene_noise), df = 6))
  base_lmdf[, "gene_splinebatchresid"] <- residuals(lmer(gene_splineresid ~ (1|batch), data = base_lmdf))
  
  strain_average <- aggregate(gene_splinebatchresid ~ strain, data = base_lmdf, FUN = mean)
  
  output <- c()
  
  output["HMTcor"] <- cor.test(strain_average[match(CeNDR_strains, strain_average$strain), "gene_splinebatchresid"], HMTsums_strain_correctedresid$HMTsplineresid_batch_df6, method = "spearman")$estimate
  output["HMTpval"] <- cor.test(strain_average[match(CeNDR_strains, strain_average$strain), "gene_splinebatchresid"], HMTsums_strain_correctedresid$HMTsplineresid_batch_df6, method = "spearman")$p.value
 
  output["ANMTcor"] <- cor.test(strain_average[match(CeNDR_strains, strain_average$strain), "gene_splinebatchresid"], HMTsums_strain_correctedresid$ANMTsplineresid_batch_df6, method = "spearman")$estimate
  output["ANMTpval"] <- cor.test(strain_average[match(CeNDR_strains, strain_average$strain), "gene_splinebatchresid"], HMTsums_strain_correctedresid$ANMTsplineresid_batch_df6, method = "spearman")$p.value
  
  output["ANMT3cor"] <- cor.test(strain_average[match(CeNDR_strains, strain_average$strain), "gene_splinebatchresid"], HMTsums_strain_correctedresid$ANMT3splineresid_batch_df6, method = "spearman")$estimate
  output["ANMT3pval"] <- cor.test(strain_average[match(CeNDR_strains, strain_average$strain), "gene_splinebatchresid"], HMTsums_strain_correctedresid$ANMT3splineresid_batch_df6, method = "spearman")$p.value
  
  output["PMT1cor"] <- cor.test(strain_average[match(CeNDR_strains, strain_average$strain), "gene_splinebatchresid"], HMTsums_strain_correctedresid$PMT1splineresid_batch_df6, method = "spearman")$estimate
  output["PMT1pval"] <- cor.test(strain_average[match(CeNDR_strains, strain_average$strain), "gene_splinebatchresid"], HMTsums_strain_correctedresid$PMT1splineresid_batch_df6, method = "spearman")$p.value
  
  return(output)
  
})

which(str_detect(names(allHMTANMTPMTcors["HMTcor", ][order(allHMTANMTPMTcors["HMTcor", ])]),"B0303.2")) / length(allHMTANMTPMTcors["HMTcor", ])
# 7.11% anmt1
which(str_detect(names(allHMTANMTPMTcors["HMTcor", ][order(allHMTANMTPMTcors["HMTcor", ])]),"T07C12.9")) / length(allHMTANMTPMTcors["HMTcor", ])
# 2.40% anmt3

which(str_detect(names(allHMTANMTPMTcors["HMTcor", ][order(allHMTANMTPMTcors["HMTcor", ])]),"ZK622.3")) / length(allHMTANMTPMTcors["HMTcor", ])
# 16.0% pmt1

allANMTcors <- allHMTANMTPMTcors["ANMTcor", ]
allANMTcors["HMTsum"] <- allHMTANMTPMTcors["HMTcor", "B0303.2"]

which(str_detect(names(allANMTcors[order(allANMTcors)]),"HMTsum")) / length(allANMTcors)
# 1.00%

allANMT3cors <- allHMTANMTPMTcors["ANMT3cor", ]
allANMT3cors["HMTsum"] <- allHMTANMTPMTcors["HMTcor", "T07C12.9"]

which(str_detect(names(allANMT3cors[order(allANMT3cors)]),"HMTsum")) / length(allANMT3cors)
# 2.46%

allPMT1cors <- allHMTANMTPMTcors["PMT1cor", ]
allPMT1cors["HMTsum"] <- allHMTANMTPMTcors["HMTcor", "ZK622.3"]

which(str_detect(names(allPMT1cors[order(allPMT1cors)]),"HMTsum")) / length(allPMT1cors)
# 9.46%

# produce scatterplots of residuals

saveRDS(HMTsums_strain_correctedresid, "plot_data/HMTsums_strain_correctedresid.rds")
HMTsums_strain_correctedresid <- readRDS("plot_data/HMTsums_strain_correctedresid.rds")

write.table(HMTsums_strain_correctedresid,
            file = "plot_data/Fig S11/Fig_S11D_plotdata.txt",
            col.names = TRUE,
            row.names = TRUE,
            sep = "\t")

pdf("~/NNMT_manuscript/graphics/CeNDR_strains_corrected_HMTvsANMT1_scatter.pdf", 
    height = 2,
    width = 1.8)

ggplot(data = HMTsums_strain_correctedresid, aes(x = HMTsplineresid_batch_df6, y = ANMTsplineresid_batch_df6)) + 
  geom_point()+
  scale_y_continuous(limits = c(-0.35, 0.25)) +
  scale_x_continuous(limits = c(-0.12, 0.12)) +
  theme_classic() + 
  ylab(substitute("Strain mean"~italic(anmt-1)~"residuals")) + 
  xlab(substitute("Strain mean total HMT residuals")) + 
  theme(axis.text.x = element_text(colour = "black",
                                   angle = 30,
                                   hjust = 1),
        axis.text.y = element_text(colour = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  geom_smooth(method = "lm",
              se = FALSE) +
  annotate(geom = "text",
           x = -0.087,
           y = -0.28,
           label = substitute(italic(r)~"= -0.33"),
           size = 2.5) +
  annotate(geom = "text",
           x = -0.063,
           y = -0.33,
           label = substitute(italic(p)~"= 1.13 x"~10^-7),
           size = 2.5)

dev.off()

format(cor.test(HMTsums_strain_correctedresid$HMTsplineresid_batch_df6, HMTsums_strain_correctedresid$ANMT3splineresid_batch_df6)$p.value, scientific = TRUE)

pdf("~/NNMT_manuscript/graphics/CeNDR_strains_corrected_HMTvsANMT3_scatter.pdf", 
    height = 2,
    width = 1.8)

ggplot(data = HMTsums_strain_correctedresid, aes(x = HMTsplineresid_batch_df6, y = ANMT3splineresid_batch_df6)) + 
  geom_point()+
  # scale_y_continuous(limits = c(-0.35, 0.25)) +
  scale_x_continuous(limits = c(-0.12, 0.12)) +
  theme_classic() + 
  ylab(substitute("Strain mean"~italic(anmt-3)~"residuals")) + 
  xlab(substitute("Strain mean total HMT residuals")) + 
  theme(axis.text.x = element_text(colour = "black",
                                   angle = 30,
                                   hjust = 1),
        axis.text.y = element_text(colour = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  geom_smooth(method = "lm",
              se = FALSE) +
  annotate(geom = "text",
           x = -0.078,
           y = -0.25,
           label = substitute(italic(r)~"= -0.47"),
           size = 2.5) +
  annotate(geom = "text",
           x = -0.053,
           y = -0.3,
           label = substitute(italic(p)~"= 1.48 x"~10^-12),
           size = 2.5)

dev.off()

pdf("~/NNMT_manuscript/graphics/CeNDR_strains_corrected_HMTvsPMT1_scatter.pdf", 
    height = 2,
    width = 1.8)

ggplot(data = HMTsums_strain_correctedresid, aes(x = HMTsplineresid_batch_df6, y = PMT1splineresid_batch_df6)) + 
  geom_point()+
  scale_x_continuous(limits = c(-0.12, 0.12)) +
  theme_classic() + 
  ylab(substitute("Strain mean"~italic(PMT-1)~"residuals")) + 
  xlab(substitute("Strain mean total HMT residuals")) + 
  theme(axis.text.x = element_text(colour = "black",
                                   angle = 30,
                                   hjust = 1),
        axis.text.y = element_text(colour = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  geom_smooth(method = "lm",
              se = FALSE) +
  annotate(geom = "text",
           x = -0.08,
           y = -0.145,
           label = substitute(italic(r)~"= -0.23"),
           size = 2.5) +
  annotate(geom = "text",
           x = -0.06,
           y = -0.18,
           label = substitute(italic(p)~"= 9.08 x"~10^-4),
           size = 2.5)

dev.off()

# volcano plot for ANMT1

# calculate cors for single HMTs
allHMTANMTPMTcors_withsingleHMTs <- apply(CeNDR_normalised_counts_collapse[row.names(CeNDR_normalised_counts_collapse) %in% Cel_HMTgenes$wormbase_gseq, ], 1, function(x){
  
  # so for each gene will need to add expression to base table featuring age and batch (or better yet add age and batch to each table with expression)
  # collapse down.
  # then do the linear model for each to regress out age and batch
  # then take the average across the strains
  # then correlate that to HMTs. 
  
  base_lmdf[, "gene_noise"] <- x + rnorm(x, 1, 0.1)
  base_lmdf[, "gene_splineresid"] <- residuals(smooth.spline(base_lmdf$age, log10(base_lmdf$gene_noise), df = 6))
  base_lmdf[, "gene_splinebatchresid"] <- residuals(lmer(gene_splineresid ~ (1|batch), data = base_lmdf))
  
  strain_average <- aggregate(gene_splinebatchresid ~ strain, data = base_lmdf, FUN = mean)
  
  output <- c()
  
  output["HMTcor"] <- cor.test(strain_average[match(CeNDR_strains, strain_average$strain), "gene_splinebatchresid"], HMTsums_strain_correctedresid$HMTsplineresid_batch_df6, method = "spearman")$estimate
  output["HMTpval"] <- cor.test(strain_average[match(CeNDR_strains, strain_average$strain), "gene_splinebatchresid"], HMTsums_strain_correctedresid$HMTsplineresid_batch_df6, method = "spearman")$p.value
  
  output["ANMTcor"] <- cor.test(strain_average[match(CeNDR_strains, strain_average$strain), "gene_splinebatchresid"], HMTsums_strain_correctedresid$ANMTsplineresid_batch_df6, method = "spearman")$estimate
  output["ANMTpval"] <- cor.test(strain_average[match(CeNDR_strains, strain_average$strain), "gene_splinebatchresid"], HMTsums_strain_correctedresid$ANMTsplineresid_batch_df6, method = "spearman")$p.value
  
  
  output["PMT1cor"] <- cor.test(strain_average[match(CeNDR_strains, strain_average$strain), "gene_splinebatchresid"], HMTsums_strain_correctedresid$PMT1splineresid_batch_df6, method = "spearman")$estimate
  output["PMT1pval"] <- cor.test(strain_average[match(CeNDR_strains, strain_average$strain), "gene_splinebatchresid"], HMTsums_strain_correctedresid$PMT1splineresid_batch_df6, method = "spearman")$p.value
  
  return(output)
  
})

allANMTcor_df <- data.frame(rho = allHMTANMTPMTcors["ANMTcor", ], pval = unlist(allHMTANMTPMTcors["ANMTpval", ]))

# add single HMTs
allANMTcor_df <- rbind(allANMTcor_df, data.frame(rho = allHMTANMTPMTcors_withsingleHMTs["ANMTcor", ], pval = unlist(allHMTANMTPMTcors_withsingleHMTs["ANMTpval", ])))

# remove anmt1 itself
allANMTcor_df <- allANMTcor_df[!(row.names(allANMTcor_df) %in% "B0303.2"), ]

# add total HMTs
allANMTcor_df["total HMTs", "rho"] <- allHMTANMTPMTcors["HMTcor", "B0303.2"]
allANMTcor_df["total HMTs", "pval"] <- allHMTANMTPMTcors["HMTpval", "B0303.2"]

# adjust p values
allANMTcor_df[, "FDR"] <- p.adjust(allANMTcor_df$pval, method = "BH")

allANMTcor_df[, "gene"] <- row.names(allANMTcor_df)

allANMTcor_df$gene <- factor(allANMTcor_df$gene, levels = allANMTcor_df[order(allANMTcor_df$rho, decreasing = FALSE), "gene"])
allANMTcor_df[, "barcolour"] <- "grey"

allANMTcor_df[allANMTcor_df$rho < quantile(allANMTcor_df$rho, 0.025, na.rm = TRUE), "barcolour"] <- "gray55"
allANMTcor_df[allANMTcor_df$rho > quantile(allANMTcor_df$rho, 0.975, na.rm = TRUE), "barcolour"] <- "gray55"

allANMTcor_HMTs_df <- allANMTcor_df[unlist(Cel_HMTgenes$wormbase_gseq), ]

allANMTcor_HMTs_df[, "barcolour"] <- "red"
allANMTcor_df["total HMTs", "barcolour"] <- "black"

allANMTcor_allHMTs_df <- allANMTcor_df["total HMTs", ]

allANMTcor_HMTs_df <- allANMTcor_HMTs_df[!row.names(allANMTcor_HMTs_df) %in% "total HMTs", ]

png("~/NNMT_manuscript/graphics/ANMT1volcano.pdf",
    height = 2.5,
    width = 2.5)

ggplot(aes(y = -log10(FDR), x = rho), data = allANMTcor_df) + 
  geom_jitter(colour = allANMTcor_df$barcolour, width = 0.02, height = 0.02, size = 0.05, alpha = 0.6) +
  theme_classic() +
  theme(axis.text.y = element_text(colour = "black",
                                   size = 10),
        axis.text.x = element_text(colour = "black",
                                   size = 10),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.line.y = element_line(colour = "black"),
        plot.margin = margin(t = 5, r = 10, b = 5, l = 5)) +
  geom_jitter(data = allANMTcor_HMTs_df, width = 0.01, height = 0.06, colour = allANMTcor_HMTs_df$barcolour, size = 0.7) + 
  geom_jitter(data = allANMTcor_allHMTs_df, width = 0.01, height = 0.06, colour = allANMTcor_allHMTs_df$barcolour, size = 2) + 
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(limits = c(-0.55, 0.55)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey") +
  xlab(substitute("Spearman's"~rho)) +
  ylab(substitute(-log[10]~"(FDR)"))

dev.off()

# volcano plot for ANMT3

# calculate cors for single HMTs
allHMTANMTPMTcors_withsingleHMTs <- apply(CeNDR_normalised_counts_collapse[row.names(CeNDR_normalised_counts_collapse) %in% Cel_HMTgenes$wormbase_gseq, ], 1, function(x){
  
  # so for each gene will need to add expression to base table featuring age and batch (or better yet add age and batch to each table with expression)
  # collapse down.
  # then do the linear model for each to regress out age and batch
  # then take the average across the strains
  # then correlate that to HMTs. 
  
  base_lmdf[, "gene_noise"] <- x + rnorm(x, 1, 0.1)
  base_lmdf[, "gene_splineresid"] <- residuals(smooth.spline(base_lmdf$age, log10(base_lmdf$gene_noise), df = 6))
  base_lmdf[, "gene_splinebatchresid"] <- residuals(lmer(gene_splineresid ~ (1|batch), data = base_lmdf))
  
  strain_average <- aggregate(gene_splinebatchresid ~ strain, data = base_lmdf, FUN = mean)
  
  output <- c()
  
  output["HMTcor"] <- cor.test(strain_average[match(CeNDR_strains, strain_average$strain), "gene_splinebatchresid"], HMTsums_strain_correctedresid$HMTsplineresid_batch_df6, method = "spearman")$estimate
  output["HMTpval"] <- cor.test(strain_average[match(CeNDR_strains, strain_average$strain), "gene_splinebatchresid"], HMTsums_strain_correctedresid$HMTsplineresid_batch_df6, method = "spearman")$p.value
  
  output["ANMT3cor"] <- cor.test(strain_average[match(CeNDR_strains, strain_average$strain), "gene_splinebatchresid"], HMTsums_strain_correctedresid$ANMT3splineresid_batch_df6, method = "spearman")$estimate
  output["ANMT3pval"] <- cor.test(strain_average[match(CeNDR_strains, strain_average$strain), "gene_splinebatchresid"], HMTsums_strain_correctedresid$ANMT3splineresid_batch_df6, method = "spearman")$p.value
  
  
  output["PMT1cor"] <- cor.test(strain_average[match(CeNDR_strains, strain_average$strain), "gene_splinebatchresid"], HMTsums_strain_correctedresid$PMT1splineresid_batch_df6, method = "spearman")$estimate
  output["PMT1pval"] <- cor.test(strain_average[match(CeNDR_strains, strain_average$strain), "gene_splinebatchresid"], HMTsums_strain_correctedresid$PMT1splineresid_batch_df6, method = "spearman")$p.value
  
  return(output)
  
})

allANMT3cor_df <- data.frame(rho = allHMTANMTPMTcors["ANMT3cor", ], pval = unlist(allHMTANMTPMTcors["ANMT3pval", ]))

# add single HMTs
allANMT3cor_df <- rbind(allANMT3cor_df, data.frame(rho = allHMTANMTPMTcors_withsingleHMTs["ANMT3cor", ], pval = unlist(allHMTANMTPMTcors_withsingleHMTs["ANMT3pval", ])))

# remove anmt3 itself
allANMT3cor_df <- allANMT3cor_df[!(row.names(allANMT3cor_df) %in% "T07C12.9"), ]

# add total HMTs
allANMT3cor_df["total HMTs", "rho"] <- allHMTANMTPMTcors["HMTcor", "T07C12.9"]
allANMT3cor_df["total HMTs", "pval"] <- allHMTANMTPMTcors["HMTpval", "T07C12.9"]

# adjust p values
allANMT3cor_df[, "FDR"] <- p.adjust(allANMT3cor_df$pval, method = "BH")

allANMT3cor_df[, "gene"] <- row.names(allANMT3cor_df)

allANMT3cor_df$gene <- factor(allANMT3cor_df$gene, levels = allANMT3cor_df[order(allANMT3cor_df$rho, decreasing = FALSE), "gene"])
allANMT3cor_df[, "barcolour"] <- "grey"

allANMT3cor_df[allANMT3cor_df$rho < quantile(allANMT3cor_df$rho, 0.025, na.rm = TRUE), "barcolour"] <- "gray55"
allANMT3cor_df[allANMT3cor_df$rho > quantile(allANMT3cor_df$rho, 0.975, na.rm = TRUE), "barcolour"] <- "gray55"

allANMT3cor_HMTs_df <- allANMT3cor_df[unlist(Cel_HMTgenes$wormbase_gseq), ]

allANMT3cor_HMTs_df[, "barcolour"] <- "red"
allANMT3cor_df["total HMTs", "barcolour"] <- "black"

allANMT3cor_allHMTs_df <- allANMT3cor_df["total HMTs", ]

allANMT3cor_HMTs_df <- allANMT3cor_HMTs_df[!row.names(allANMT3cor_HMTs_df) %in% "total HMTs", ]

write.table(allANMT3cor_df,
            file = "plot_data/Fig S11/Fig_S11E_plotdata.txt",
            row.names = TRUE,
            col.names = TRUE,
            sep = "\t")

pdf("~/NNMT_manuscript/graphics/ANMT3volcano.pdf",
    height = 2.5,
    width = 2.5)

ggplot(aes(y = -log10(FDR), x = rho), data = allANMT3cor_df) + 
  geom_jitter(colour = allANMT3cor_df$barcolour, width = 0.02, height = 0.02, size = 0.05, alpha = 0.6) +
  theme_classic() +
  theme(axis.text.y = element_text(colour = "black",
                                   size = 10),
        axis.text.x = element_text(colour = "black",
                                   size = 10),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.line.y = element_line(colour = "black"),
        plot.margin = margin(t = 5, r = 10, b = 5, l = 5)) +
  geom_jitter(data = allANMT3cor_HMTs_df, width = 0.01, height = 0.06, colour = allANMT3cor_HMTs_df$barcolour, size = 0.7) + 
  geom_jitter(data = allANMT3cor_allHMTs_df, width = 0.01, height = 0.06, colour = allANMT3cor_allHMTs_df$barcolour, size = 2) + 
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(limits = c(-0.8, 0.8)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey") +
  xlab(substitute("Spearman's"~rho)) +
  ylab(substitute(-log[10]~"(FDR)"))

dev.off()

#### heatmap of HMT residuals Fig S11A ####

CeNDR_HMT_reads_tdf <- data.frame(t(CeNDR_normalised_counts_collapse[Cel_HMTgenes$wormbase_gseq, ]))
CeNDR_HMT_reads_tdf[, "age"] <- andersen_ae2$age.estimates[colnames(CeNDR_normalised_counts_collapse), "age.estimate"]
CeNDR_HMT_reads_tdf[, "batch"] <- str_extract(colnames(CeNDR_normalised_counts_collapse), "[0-9]{2}_[0-9]{4}$")

# base_lmdf[, "gene_splineresid"] <- residuals(smooth.spline(base_lmdf$age, log10(base_lmdf$gene_noise)))
# base_lmdf[, "gene_splinebatchresid"] <- residuals(lmer(gene_splineresid ~ (1|batch), data = base_lmdf))

CeNDR_HMT_logreads_corrected <- data.frame(sapply(colnames(CeNDR_HMT_reads_tdf)[1:(ncol(CeNDR_HMT_reads_tdf)-2)], function(thisgene){
  
  thisgene_expression <- CeNDR_HMT_reads_tdf[, paste0(thisgene)]

  tempsplineresid <- residuals(smooth.spline(CeNDR_HMT_reads_tdf$age, log10(thisgene_expression + 1), df = 6))
  residuals(lmer(tempsplineresid ~ (1|CeNDR_HMT_reads_tdf$batch)))
  
}))

row.names(CeNDR_HMT_logreads_corrected) <- row.names(CeNDR_HMT_reads_tdf)

CeNDR_HMT_logreads_corrected[, "strain"] <- str_extract(row.names(CeNDR_HMT_logreads_corrected), "^[0-9A-Z]+")

# try across strains
CeNDR_HMT_logreads_corrected_strain <- aggregate(. ~ strain, data = CeNDR_HMT_logreads_corrected, FUN = mean)
row.names(CeNDR_HMT_logreads_corrected_strain) <- CeNDR_HMT_logreads_corrected_strain$strain
CeNDR_HMT_logreads_corrected_strain <- CeNDR_HMT_logreads_corrected_strain[, 2:ncol(CeNDR_HMT_logreads_corrected_strain)]

Cel_strains_HMTresiduals_cor <- cor(CeNDR_HMT_logreads_corrected_strain)

Cel_strains_HMTresiduals_hclustering <- hclust(dist(Cel_strains_HMTresiduals_cor))

# cut into two clusters
CelHMT_clusters <- cutree(Cel_strains_HMTresiduals_hclustering, k = 2)
names(CelHMT_clusters) <- Cel_HMTgenes[match(names(CelHMT_clusters), Cel_HMTgenes$wormbase_gseq), "wormbase_locus"]

# for Cel, take cluster 2 as HE
Cel_HEcluster <- names(CelHMT_clusters[CelHMT_clusters == 2])
Cel_notHEcluster <- names(CelHMT_clusters[CelHMT_clusters == 1])

pdf("graphics/Cel_strains_HMT_clustering_dendrogram.pdf")

plot(Cel_strains_HMTresiduals_hclustering, 
     hang = -1, # labels same height
     main = NA,
     xlab = NA,
     ylab = NA)

dev.off()

Cel_clust_order <- matrix(nrow = nrow(Cel_strains_HMTresiduals_cor), ncol = ncol(Cel_strains_HMTresiduals_cor))
colnames(Cel_clust_order) <- Cel_strains_HMTresiduals_hclustering$labels[rev(Cel_strains_HMTresiduals_hclustering$order)]
row.names(Cel_clust_order) <- Cel_strains_HMTresiduals_hclustering$labels[rev(Cel_strains_HMTresiduals_hclustering$order)]

for(i in 1:nrow(Cel_strains_HMTresiduals_cor)){
  
  for(j in 1:ncol(Cel_strains_HMTresiduals_cor)) {
    
    Cel_clust_order[Cel_strains_HMTresiduals_hclustering$labels[i], Cel_strains_HMTresiduals_hclustering$labels[j]] <- Cel_strains_HMTresiduals_cor[Cel_strains_HMTresiduals_hclustering$labels[i], Cel_strains_HMTresiduals_hclustering$labels[j]]
    
  }
  
}

row.names(Cel_clust_order) <- Cel_HMTgenes[match(row.names(Cel_clust_order), Cel_HMTgenes$wormbase_gseq), "wormbase_locus"]
colnames(Cel_clust_order) <- Cel_HMTgenes[match(colnames(Cel_clust_order), Cel_HMTgenes$wormbase_gseq), "wormbase_locus"]

# for making expression sidebar to heatmap
Cel_HMT_means_gmmeans <- apply(CeNDR_HMT_reads_tdf[, 1:29], 2, gm_mean)

names(Cel_HMT_means_gmmeans) <- Cel_HMTgenes[match(names(Cel_HMT_means_gmmeans), Cel_HMTgenes$wormbase_gseq), "wormbase_locus"]

Cel_HMT_logmeans <- log10(Cel_HMT_means_gmmeans)

# colours for plot
mainpal <- (colorRampPalette(c("blue", "white", "red"))(100))
side_colors <- viridis::inferno(100)

Cel_HMT_sidebar <- map2color(Cel_HMT_logmeans,
                             side_colors
)

names(Cel_HMT_sidebar) <- names(Cel_HMT_logmeans)

Cel_HMT_sidebar <- Cel_HMT_sidebar[row.names(Cel_clust_order)]

Cel_HMTtarget_sidebar <- Cel_HMTgenes$target
names(Cel_HMTtarget_sidebar) <- Cel_HMTgenes$wormbase_locus

Cel_HMTtarget_sidebar[Cel_HMTtarget_sidebar == "H3K4"] <- "greenyellow"
Cel_HMTtarget_sidebar[Cel_HMTtarget_sidebar == "H3K9"] <- "magenta"
Cel_HMTtarget_sidebar[Cel_HMTtarget_sidebar == "H3K27"] <- "purple"
Cel_HMTtarget_sidebar[Cel_HMTtarget_sidebar == "H3K36"] <- "green4"
Cel_HMTtarget_sidebar[Cel_HMTtarget_sidebar == "H3K79"] <- "deepskyblue"
Cel_HMTtarget_sidebar[Cel_HMTtarget_sidebar == "H4K20"] <- "orchid1"
Cel_HMTtarget_sidebar[Cel_HMTtarget_sidebar == "unknown"] <- "grey"

Cel_HMTtarget_sidebar <- Cel_HMTtarget_sidebar[row.names(Cel_clust_order)]

saveRDS(Cel_clust_order, "plot_data/Cel_clust_order.rds")
# Cel_clust_order <- readRDS("plot_data/Cel_clust_order.rds")

#### human C. elegans orthology ####

# submitted list of worm hmts to Ortholist2 [http://ortholist.shaye-lab.org/]
ortholist_results <- read.csv("input/ortholist_results.csv")

# restrict ortholist results to orthologues in at least 3 / 6 databases
ortholist_results <- ortholist_results[ortholist_results$No..of.Databases > 2, ]

human_orthologues <- lapply(unique(ortholist_results$Common.Name), function(thiscommonname){

  ortholist_results[ortholist_results$Common.Name == thiscommonname, "HGNC.Symbol"]
  
})

names(human_orthologues) <- unique(ortholist_results$Common.Name)

Cel_HEorthologues_logical <- sapply(human_orthologues[Cel_HEcluster], function(thisone){

  any(thisone %in% both_HEclusters)
  
})

names(Cel_HEorthologues_logical) <- Cel_HEcluster

sum(Cel_HEorthologues_logical)

Cel_notHEorthologues_logical <- sapply(human_orthologues[Cel_notHEcluster], function(thisone){
  
  any(thisone %in% both_HEclusters)
  
})

names(Cel_notHEorthologues_logical) <- Cel_notHEcluster

orthology_fisher_df <- data.frame(notinhsHE = c(sum(!Cel_notHEorthologues_logical), sum(!Cel_HEorthologues_logical)),
inhsHE =  c(sum(Cel_notHEorthologues_logical), sum(Cel_HEorthologues_logical)),
row.names = c("notHE", "Cel-HE"))

fisher.test(orthology_fisher_df)

write.table(Cel_clust_order,
            file = "plot_data/Fig S11/Fig_S11A_plotdata.txt",
            col.names = TRUE,
            row.names = TRUE,
            sep = "\t")

pdf("~/NNMT_manuscript/graphics/Celegans_HMT_heatmap_KEY.pdf",
    width = 4,
    height = 4)

heatmap.2(Cel_clust_order,
          col = mainpal,
          dendrogram = "none",
          Rowv = FALSE,
          Colv = FALSE,
          trace = "none",
          RowSideColors = Cel_HMT_sidebar,
          density.info = "none",
          cexRow = 0.8,
          cexCol = 0.8,
          keysize = 1.8,
          key.title = NA)

dev.off()

pdf("~/NNMT_manuscript/graphics/Celegans_HMT_heatmap_NOKEY.pdf",
    width = 4,
    height = 4)

heatmap.2(Cel_clust_order,
          col = mainpal,
          dendrogram = "none",
          Rowv = FALSE,
          Colv = FALSE,
          trace = "none",
          RowSideColors = Cel_HMT_sidebar,
          density.info = "none",
          cexRow = 0.8,
          cexCol = 0.8,
          key = FALSE,
          margins = c(5,5))

dev.off()

pdf("~/NNMT_manuscript/graphics/Celegans_HMTtarget_heatmap_noKEY.pdf",
    width = 4,
    height = 4)

heatmap.2(Cel_clust_order,
          col = mainpal,
          dendrogram = "none",
          Rowv = FALSE,
          Colv = FALSE,
          trace = "none",
          RowSideColors = Cel_HMT_sidebar,
          ColSideColors = Cel_HMTtarget_sidebar,
          density.info = "none",
          cexRow = 0.8,
          cexCol = 0.8,
          key = FALSE)

dev.off()

Cel_HE <- row.names(Cel_clust_order)[20:29]
Cel_notHE <- names(Cel_HMT_means_gmmeans[Cel_HMT_means_gmmeans > 100])[!(names(Cel_HMT_means_gmmeans[Cel_HMT_means_gmmeans > 100]) %in% Cel_HE)]
Cel_HE_gseq <- Cel_HMTgenes[match(Cel_HE, Cel_HMTgenes$wormbase_locus), "wormbase_gseq"]

#### INTRA-HMT CORRELATIONS HUMAN for Fig S10 ####

correlations.within.HMTs <- function(database = c("GTEX", "TCGA"),
                                     corplotmethod = c("color", "circle", "square", "ellipse", "number", "shade", "colour", "pie"),
                                     iterations = 10){
  
  database <- match.arg(database)
  method <- match.arg(corplotmethod)

  if(database == "GTEX"){
    MOR_list_HMTs <- lapply(GTEX_MOR_list, function(x){x[SET_HMTs$ensembl_gene_id, ]})
  }
  
  if(database == "TCGA"){

    MOR_list_HMTs <- lapply(TCGA_MOR_list, function(x){x[SET_HMTs$ensembl_gene_id, ]})
    
    MOR_HMTs_combined <- do.call(cbind, MOR_list_HMTs)
    
  lookuptable <- readRDS("~/NNMT_manuscript/output/TCGA_LUT.rds")
    
  residuals_by_cancer <- lapply(MOR_list_HMTs, function(countdata){

    # restrict to samples present in lookuptable 
    tempcountdata <- countdata[, colnames(countdata) %in% lookuptable$SAMPID]
    
    temp_df <- as.data.frame(matrix(nrow = ncol(tempcountdata), ncol = nrow(tempcountdata) + 5))
    
    colnames(temp_df) <- c(row.names(tempcountdata), "race", "gender", "tumour_stage", "sequencing_centre", "days_to_birth")
    row.names(temp_df) <- colnames(tempcountdata)
    
    # insert tranposed data into new data frame
    temp_df[, 1:nrow(tempcountdata)] <- t(tempcountdata)
    
    # from lookup table, add variables
    temp_df[, "days_to_birth"] <- as.numeric(lookuptable[row.names(temp_df), "days_to_birth"])
    temp_df[, "gender"] <- lookuptable[row.names(temp_df), "gender"]
    temp_df[, "race"] <- lookuptable[row.names(temp_df), "race"]
    temp_df[, "tumour_stage"] <- lookuptable[row.names(temp_df), "tumour_stage"]
    temp_df[, "sequencing_centre"] <- lookuptable[row.names(temp_df), "sequencing_centre"]
    
    # check for any missing values and remove entry if so
    temp_df <- temp_df[!(is.na(temp_df$days_to_birth)), ]
    temp_df <- temp_df[!(is.na(temp_df$gender)), ]
    temp_df <- temp_df[!(is.na(temp_df$race)), ]
    temp_df <- temp_df[!(is.na(temp_df$tumour_stage)), ]
    temp_df <- temp_df[!(is.na(temp_df$sequencing_centre)), ]
    
    # do linear regression across all samples (all tissues combined)
    # perform linear regression with age (bracket midpoint), sex and cause of death
    
    # create dataframe to deposit residuals
    residuals_df <- as.data.frame(matrix(nrow = nrow(temp_df), ncol = ncol(temp_df) - 5))
    row.names(residuals_df) <- row.names(temp_df)
    colnames(residuals_df) <- colnames(temp_df)[1:(ncol(temp_df) - 5)]
    
    # loop over columns (genes); for each gene perform linear regression and deposit residuals in new data frame
    # sex is considered a factor automatically as a character string; sequencing centre is already considered as factor
    # here we note that in some projects, all sequencing centres are equal. Likewise gender for some cancer tpyes (e.g. ovarian)
    # need conditional to ignore seq centre, gender or tumour_stage if only one level exists among considered data in order to avoid an error
    
    if(length(unique(temp_df$race)) == 1){
      
      racevariable <- NULL
      
    } else {
      
      racevariable <- "race"
      
    }
    
    if(length(unique(temp_df$tumour_stage)) == 1){
      
      tumourstagevariable <- NULL
      
    } else {
      
      tumourstagevariable <- "tumour_stage"
      
    }
    
    if(length(unique(temp_df$gender)) == 1) {
      
      gendervariable <- NULL
      
    } else {
      
      gendervariable <- "gender"
    }
    
    if(length(unique(temp_df$sequencing_centre)) == 1) {
      
      seqcentrevariable <- NULL
      
    } else {
      
      seqcentrevariable <- "sequencing_centre"
      
    }
    
    modelvariables <- c("days_to_birth", gendervariable, seqcentrevariable, tumourstagevariable, racevariable)
    
    for(j in 1:(ncol(residuals_df))){
      
      outcome <- paste0("temp_df[, ", j, "]")
      
      f <- as.formula(
        paste(outcome,
              paste(modelvariables, collapse = " + "),
              sep = " ~ "))
      
      residuals_df[, colnames(temp_df)[j]] <- lm(formula = f, data = temp_df)$residuals
      
    } # end of residuals loop
    
    
    return(t(residuals_df))
    
    
  })

  residuals_df <- do.call(cbind, residuals_by_cancer)
      
  }
  
  tissues <- names(MOR_list_HMTs)
  
  if(database == "GTEX"){
    
    lookuptable <- readRDS("~/manuscript/output/GTEX-version8-sampleID-LUT.rds")
   
    if(!exists("GTEX_MOR_across_tissues")){
    GTEX_MOR_across_tissues <- readRDS("output/GTEX_MOR_across_tissues.rds")
    }

    MOR_HMTs_combined <- GTEX_MOR_across_tissues[SET_HMTs$ensembl_gene_id, ]
    
    GTEX_matrix_for_lm_temp <- data.frame(t(MOR_HMTs_combined))

    GTEX_matrix_for_lm_temp[, "age_bracket"] <- sapply(str_split(lookuptable[row.names(GTEX_matrix_for_lm_temp), "age_bracket"], 
                                                 pattern = "-"), function(x){
                                                   
                                                   mean(as.numeric(x))
                                                   
                                                 })    
    
    GTEX_matrix_for_lm_temp[, "sex"] <- lookuptable[row.names(GTEX_matrix_for_lm_temp), "Sex"]
    GTEX_matrix_for_lm_temp[, "death"] <- lookuptable[row.names(GTEX_matrix_for_lm_temp), "DeathScale"]
    GTEX_matrix_for_lm_temp[, "ischemic_time"] <- lookuptable[row.names(GTEX_matrix_for_lm_temp), "IschemicTime"]

    
    # replace death = NA (unknown cause of death) with a 5th factor level
    GTEX_matrix_for_lm_temp[is.na(GTEX_matrix_for_lm_temp$death), "death"] <- 5
    
    # remove samples with NA for ischemic time
    GTEX_matrix_for_lm_temp <- GTEX_matrix_for_lm_temp[!is.na(GTEX_matrix_for_lm_temp$ischemic_time), ]
    
    # ischemic time rescaled so lmer doesn't complain about scale issues
    GTEX_matrix_for_lm_temp[, "ischemic_time"] <- lookuptable[row.names(GTEX_matrix_for_lm_temp), "IschemicTime"]/60 
    
    GTEX_matrix_for_lm_temp[, "batch"] <- lookuptable[row.names(GTEX_matrix_for_lm_temp), "Batch"]
    
    GTEX_matrix_for_lm_temp[, "tissue"] <- lookuptable[row.names(GTEX_matrix_for_lm_temp), "Tissue"]
    
    GTEX_matrix_for_lm_temp[, "donorID"] <- lookuptable[row.names(GTEX_matrix_for_lm_temp), "SUBJID"]
    
    # check for any missing values and remove entry if so

    GTEX_matrix_for_lm_temp <- GTEX_matrix_for_lm_temp[!(is.na(GTEX_matrix_for_lm_temp$ischemic_time)), ]

    # do linear regression across all samples (all tissues combined)
    # perform linear regression with age (bracket midpoint), sex and cause of death
    
    # create dataframe to deposit residuals
    residuals_df <- as.data.frame(matrix(nrow = nrow(GTEX_matrix_for_lm_temp), ncol = ncol(GTEX_matrix_for_lm_temp) - 7))
    row.names(residuals_df) <- row.names(GTEX_matrix_for_lm_temp)
    colnames(residuals_df) <- colnames(GTEX_matrix_for_lm_temp)[1:(ncol(GTEX_matrix_for_lm_temp) - 7)]
    
    # loop over columns (genes); for each gene perform linear regression and deposit residuals in new data frame
    # sex is considered a factor automatically as a character string; as.factor() must be specified for cause of death

    for(j in 1:(ncol(residuals_df))){

      residuals_df[, colnames(GTEX_matrix_for_lm_temp)[j]] <- resid(lmer(GTEX_matrix_for_lm_temp[, j] ~ (sex + as.factor(death) + age_bracket + ischemic_time)*tissue + (1|batch) + (1|donorID), data = GTEX_matrix_for_lm_temp))
      
    }
    
    residuals_df <- data.frame(t(residuals_df))
    
  }

  # get rid of HMTs with very low expression
  retainedHMTs <- row.names(MOR_HMTs_combined[apply(MOR_HMTs_combined, 1, gm_mean) > 100, ])
  MOR_HMTs_combined <- MOR_HMTs_combined[retainedHMTs, ]
  
  residuals_df2 <- data.frame(residuals_df[retainedHMTs, ])
  
  output_list <- list()

  iterations_correlations <- lapply(1:iterations, function(i){
  
  tissuesamples_rankspercent_df <- as.data.frame(matrix(nrow = nrow(residuals_df2), ncol = 0))
  row.names(tissuesamples_rankspercent_df) <- row.names(residuals_df2)
  
  for(j in 1:length(tissues)){
    
    if(database == "GTEX"){
      
      all_tissuesamples <- lookuptable[lookuptable$Tissue == tissues[j], "SAMPID"]
      
      n = 100
      
    }
    
    if(database == "TCGA"){

      all_tissuesamples <- lookuptable[lookuptable$cancer == tissues[j], "SAMPID"]
      all_tissuesamples <- str_replace_all(all_tissuesamples, "-", "\\.")
      
      n = 36
      
    }
    
    present_tissuesamples <- all_tissuesamples[all_tissuesamples %in% colnames(residuals_df2)]
    
    # take 100 random samples from tissue and restrict to OXPHOS genes
    samplecounts <- residuals_df2[, sample(present_tissuesamples, n)]
    
    # rank the samples
    tempdata_ranks <- matrix(nrow = nrow(samplecounts), ncol = ncol(samplecounts))
    row.names(tempdata_ranks) <- row.names(samplecounts)
    colnames(tempdata_ranks) <- colnames(samplecounts)
    
    for (k in 1:nrow(samplecounts)){
      tempdata_ranks[k, ] <- rank(samplecounts[k, ])
    } 
    
    # here we add the ranks for this tissue to those previously calculated
    tissuesamples_rankspercent_df <- merge(tissuesamples_rankspercent_df, tempdata_ranks, all.x = TRUE, all.y = FALSE, by = "row.names")
    row.names(tissuesamples_rankspercent_df) <- tissuesamples_rankspercent_df[, "Row.names"]
    tissuesamples_rankspercent_df <- tissuesamples_rankspercent_df[, 2:ncol(tissuesamples_rankspercent_df)]
    
  }
  
  row.names(tissuesamples_rankspercent_df) <- SET_HMTs[match(row.names(tissuesamples_rankspercent_df), SET_HMTs$ensembl_gene_id), "hgnc_symbol"]
  
  pan_HMT_cor <- Hmisc::rcorr(t(tissuesamples_rankspercent_df), type = "spearman")

  return(pan_HMT_cor)
  
  }) #end of lapply over iterations
  
  # take average correlation
  
  iterations_correlations_alone <- lapply(iterations_correlations, function(x){x[["r"]]})
  iterations_pvals_alone <- lapply(iterations_correlations, function(x){x[["P"]]})
  
  correlations_array <- array(as.numeric(unlist(iterations_correlations_alone)), dim = c(nrow(iterations_correlations_alone[[1]]), ncol(iterations_correlations_alone[[1]]), length(iterations_correlations_alone)))
  average_correlations <- apply(correlations_array , 1:2 , mean )
  
  row.names(average_correlations) <- row.names(iterations_correlations_alone[[1]])
  colnames(average_correlations) <- colnames(iterations_correlations_alone[[1]])
  
  pval_array <- array(as.numeric(unlist(iterations_pvals_alone)), dim = c(nrow(iterations_pvals_alone[[1]]), ncol(iterations_pvals_alone[[1]]), length(iterations_pvals_alone)))
  average_pvals <- apply(pval_array , 1:2 , mean)
  
  pval_row_list <- lapply(2:nrow(average_pvals), function(i){
    
    average_pvals[i, 1:(i-1)]
    
  })
  
  adjust_pval_vec <- p.adjust(unlist(pval_row_list), method = "BH")
  adjust_pval_vec_temp <- adjust_pval_vec
  
  adjusted_matrix <- matrix(nrow = nrow(average_pvals), 
                            ncol = ncol(average_pvals))
  
for(i in nrow(average_pvals):2){

  adjusted_matrix[i, 1:(i-1)] <- tail(adjust_pval_vec_temp, n = length(pval_row_list[[i-1]]))
  
  adjust_pval_vec_temp <- head(adjust_pval_vec_temp, n = (length(adjust_pval_vec_temp) - length(pval_row_list[[i-1]])))
  
}
  
  row.names(adjusted_matrix) <- row.names(iterations_correlations_alone[[1]])
  colnames(adjusted_matrix) <- colnames(iterations_correlations_alone[[1]])
  
  hclustering_temp <- hclust(dist(average_correlations))
  hclustering_temp$labels
  
  pan_HMT_cor_order <- matrix(nrow = nrow(average_correlations), ncol = ncol(average_correlations))
  colnames(pan_HMT_cor_order) <- hclustering_temp$labels[rev(hclustering_temp$order)]
  row.names(pan_HMT_cor_order) <- hclustering_temp$labels[rev(hclustering_temp$order)]
  
  for(i in 1:nrow(average_correlations)){
    
    for(j in 1:ncol(average_correlations)) {
      
      pan_HMT_cor_order[hclustering_temp$labels[i], hclustering_temp$labels[j]] <- average_correlations[hclustering_temp$labels[i], hclustering_temp$labels[j]]
      
    }
    
  }
  
  # for pvals
  pval_HMT_cor_order <- matrix(nrow = nrow(average_correlations), ncol = ncol(average_correlations))
  colnames(pval_HMT_cor_order) <- hclustering_temp$labels[rev(hclustering_temp$order)]
  row.names(pval_HMT_cor_order) <- hclustering_temp$labels[rev(hclustering_temp$order)]
  
  for(i in 1:nrow(average_correlations)){
    
    for(j in 1:ncol(average_correlations)) {

      thisone <- adjusted_matrix[hclustering_temp$labels[i], hclustering_temp$labels[j]]
      thatone <- adjusted_matrix[hclustering_temp$labels[j], hclustering_temp$labels[i]]
      
      together_these <- c(thisone, thatone)
      
      if(sum(is.na(together_these)) == 1){
      together_these <- together_these[!is.na(together_these)]
      pval_HMT_cor_order[hclustering_temp$labels[i], hclustering_temp$labels[j]] <- together_these
      }

      if(sum(is.na(together_these)) == 2){
        pval_HMT_cor_order[hclustering_temp$labels[i], hclustering_temp$labels[j]] <- NA
      }
      
    }
    
  }
  
  pdf(paste0("graphics/", database, "_HMT_corrplots_allcombined.pdf"))
  
  corrplot(pan_HMT_cor_order, 
           method = method,
           tl.col = SET_HMTs[match(row.names(pan_HMT_cor_order),
                                   SET_HMTs$hgnc_symbol), "colour"])
  
  dev.off()
  
  allcombined_list <- list(correlations = pan_HMT_cor_order,
                           FDR = pval_HMT_cor_order)
  
  output_list[["allcombined"]] <- allcombined_list
  
  HMT_cormatrices <- lapply(tissues, function(thistissue){

    if(database == "GTEX"){
    thistissue_residuals <- residuals_df2[, colnames(residuals_df2) %in% lookuptable[lookuptable$Tissue == thistissue, "SAMPID"]]
    }
    
    if(database == "TCGA"){
      thistissue_residuals <- residuals_df2[, colnames(residuals_df2) %in% str_replace_all(lookuptable[lookuptable$cancer == thistissue, "SAMPID"], "-", "\\.")]
    }
    
    row.names(thistissue_residuals) <- SET_HMTs[match(row.names(thistissue_residuals), SET_HMTs$ensembl_gene_id), "hgnc_symbol"]

    cor(t(thistissue_residuals))
    
  })
  
  names(HMT_cormatrices) <- tissues
  
  for(i in 1:length(tissues)){
    
    thistissue_cormatrix <- HMT_cormatrices[[tissues[i]]]
    
    thistissue_order <- matrix(nrow = nrow(pan_HMT_cor_order),
                               ncol = ncol(pan_HMT_cor_order)) 
    
    row.names(thistissue_order) <- row.names(pan_HMT_cor_order)
    colnames(thistissue_order) <- colnames(pan_HMT_cor_order)
    
    for(k in 1:nrow(pan_HMT_cor_order)){
      
      for(j in 1:ncol(pan_HMT_cor_order)) {
        
        thistissue_order[hclustering_temp$labels[k], hclustering_temp$labels[j]] <- thistissue_cormatrix[hclustering_temp$labels[k], hclustering_temp$labels[j]]
        
      }
      
    }
    
    output_list[[paste0(tissues[i])]] <- thistissue_order
    
  }
  
  pdf(paste0("graphics/", database, "_HMT_corrplots_tissue.pdf"))
  
  for(i in 1:length(output_list)){
    
    if(i == 1){
      
      corrplot(output_list[[i]][["correlations"]], tl.col = SET_HMTs[match(row.names(output_list[[i]]), SET_HMTs$hgnc_symbol), "colour"], 
               main = names(output_list)[i],
               mar=c(0,0,1,0),
               method = method)
      
    } else {
    
    corrplot(output_list[[i]], tl.col = SET_HMTs[match(row.names(output_list[[i]]), SET_HMTs$hgnc_symbol), "colour"], 
             main = names(output_list)[i],
             mar=c(0,0,1,0),
             method = method)
      
    }
    
  }
  
  dev.off()
  
  return(output_list)
  
}

GTEX_correlations_within_HMTs <- correlations.within.HMTs(database = "GTEX",
                                     corplotmethod = "color",
                                     iterations = 100)

saveRDS(GTEX_correlations_within_HMTs, "output/GTEX_correlations_within_HMTs.rds")
# GTEX_correlations_within_HMTs <- readRDS("output/GTEX_correlations_within_HMTs.rds")

TCGA_correlations_within_HMTs <- correlations.within.HMTs(database = "TCGA",
                                                          corplotmethod = "color",
                                                          iterations = 100)

saveRDS(TCGA_correlations_within_HMTs, "output/TCGA_correlations_within_HMTs.rds")
# TCGA_correlations_within_HMTs <- readRDS("output/TCGA_correlations_within_HMTs.rds")

#### plot heatmaps ####

GTEX_MOR_across_tissues <- readRDS("output/GTEX_MOR_across_tissues.rds")

allGTEX_HMT_gmmeans <- apply(GTEX_MOR_across_tissues[SET_HMTs$ensembl_gene_id, ], 1, gm_mean)
allGTEX_HMT_logmeans <- log10(allGTEX_HMT_gmmeans)
names(allGTEX_HMT_logmeans) <- SET_HMTs[match(names(allGTEX_HMT_logmeans), SET_HMTs$ensembl_gene_id), "hgnc_symbol"]

# colours for plot
mainpal <- (colorRampPalette(c("blue", "white", "red"))(100))
side_colors <- viridis::inferno(100)

GTEX_act_sidebar <- colnames(GTEX_correlations_within_HMTs$allcombined$correlations)
names(GTEX_act_sidebar) <- colnames(GTEX_correlations_within_HMTs$allcombined$correlations)
selection_colours <- c("green", "magenta", "grey")

for(i in 1:3){

    GTEX_act_sidebar[act_rep_list_hgnc[[i]][act_rep_list_hgnc[[i]] %in% names(GTEX_act_sidebar)]] <- selection_colours[i]
    
  }
  
GTEX_all_sidebar <- map2color(allGTEX_HMT_logmeans[names(allGTEX_HMT_logmeans) %in% colnames(GTEX_correlations_within_HMTs$allcombined$correlations)],
                              side_colors,
                              limits = c(2.3, 3.9))

names(GTEX_all_sidebar) <- names(allGTEX_HMT_logmeans)[names(allGTEX_HMT_logmeans) %in% colnames(GTEX_correlations_within_HMTs$allcombined$correlations)]

GTEX_all_sidebar <- GTEX_all_sidebar[row.names(GTEX_correlations_within_HMTs[["allcombined"]][["correlations"]])]

write.table(GTEX_correlations_within_HMTs[["allcombined"]][["correlations"]],
            file = "plot_data/Fig S10/Fig_S10A_GTEX_plotdata.txt",
            col.names = TRUE,
            row.names = FALSE,
            sep = "\t")

pdf("graphics/GTEX_HMTcorrs_heatmapexpression_combined.pdf",
    height = 4,
    width = 4)

heatmap.2(GTEX_correlations_within_HMTs[["allcombined"]][["correlations"]],
          col = mainpal,
          dendrogram = "none",
          Rowv = FALSE,
          Colv = FALSE,
          trace = "none",
          RowSideColors = GTEX_all_sidebar,
          density.info = "none",
          cexRow = 0.56,
          cexCol = 0.56,
          keysize = 1.8,
          key.title = NA)

dev.off()

pdf("graphics/GTEX_HMTcorrs_heatmapexpression_combined_actsidebar.pdf",
    height = 4,
    width = 4)

heatmap.2(GTEX_correlations_within_HMTs[["allcombined"]][["correlations"]],
          col = mainpal,
          dendrogram = "none",
          Rowv = FALSE,
          Colv = FALSE,
          trace = "none",
          RowSideColors = GTEX_act_sidebar,
          density.info = "none",
          cexRow = 0.56,
          cexCol = 0.56,
          keysize = 1.8,
          key.title = NA)

dev.off()

GTEXhclustering <- hclust(dist(GTEX_correlations_within_HMTs[["allcombined"]][["correlations"]]))

pdf("graphics/GTEX_HMT_clustering_dendrogram.pdf")

plot(GTEXhclustering, 
     hang = -1, # labels same height
     main = NA,
     xlab = NA,
     ylab = NA)

dev.off()

# cut into two or three clusters
GTEXclusters <- cutree(GTEXhclustering, k = 3)

# for GTEX, take cluster 3 as HE
GTEX_HEcluster <- names(GTEXclusters[GTEXclusters == 3])

TCGAhclustering <- hclust(dist(TCGA_correlations_within_HMTs[["allcombined"]][["correlations"]]))

#visualise dendrogram

pdf("graphics/TCGA_HMT_clustering_dendrogram.pdf")

plot(TCGAhclustering, 
     hang = -1, # labels same height
     main = NA,
     xlab = NA,
     ylab = NA)

dev.off()

# cut into 4 clusters. Take middle two
TCGAclusters <- cutree(TCGAhclustering, k = 4)

# for TCGA, take cluster 3 as HE
TCGA_HEcluster <- names(TCGAclusters[TCGAclusters %in% c(2,3)])

# whats in both?
both_HEclusters <- base::union(GTEX_HEcluster, 
                         TCGA_HEcluster)

# what's the same?
consensus_HEclusters <- intersect(GTEX_HEcluster, 
                TCGA_HEcluster)

# what's different?
setdiff(GTEX_HEcluster, 
          TCGA_HEcluster)

TCGA_MOR_across_tissues <- readRDS("~/manuscript/output/TCGA-MOR-normalisation-across-cancers.rds")

allTCGA_HMT_gmmeans <- apply(TCGA_MOR_across_tissues[SET_HMTs$ensembl_gene_id, ], 1, gm_mean)
allTCGA_HMT_logmeans <- log10(allTCGA_HMT_gmmeans)
names(allTCGA_HMT_logmeans) <- SET_HMTs[match(names(allTCGA_HMT_logmeans), SET_HMTs$ensembl_gene_id), "hgnc_symbol"]

TCGA_all_sidebar <- map2color(allTCGA_HMT_logmeans[names(allTCGA_HMT_logmeans) %in% colnames(TCGA_correlations_within_HMTs$allcombined$correlations)],
                              side_colors,
                              limits = c(2.3, 3.9))

names(TCGA_all_sidebar) <- names(allTCGA_HMT_logmeans)[names(allTCGA_HMT_logmeans) %in% colnames(TCGA_correlations_within_HMTs$allcombined$correlations)]

TCGA_all_sidebar <- TCGA_all_sidebar[row.names(TCGA_correlations_within_HMTs[["allcombined"]][["correlations"]])]

TCGA_act_sidebar <- colnames(TCGA_correlations_within_HMTs$allcombined$correlations)
names(TCGA_act_sidebar) <- colnames(TCGA_correlations_within_HMTs$allcombined$correlations)

for(i in 1:3){
  
  TCGA_act_sidebar[act_rep_list_hgnc[[i]][act_rep_list_hgnc[[i]] %in% names(TCGA_act_sidebar)]] <- selection_colours[i]
  
}

write.table(TCGA_correlations_within_HMTs[["allcombined"]][["correlations"]],
            file = "plot_data/Fig S10/Fig_S10A_TCGA_plotdata.txt",
            row.names = TRUE,
            col.names = TRUE,
            sep = "\t")

pdf("graphics/TCGA_HMTcorrs_heatmapexpression_combined.pdf",
    width = 4,
    height = 4)

heatmap.2(TCGA_correlations_within_HMTs[["allcombined"]][["correlations"]],
          col = mainpal,
          dendrogram = "none",
          Rowv = FALSE,
          Colv = FALSE,
          trace = "none",
          RowSideColors = TCGA_all_sidebar,
          density.info = "none",
          cexRow = 0.6,
          cexCol = 0.6,
          keysize = 1.8,
          key.title = NA)

dev.off()

pdf("graphics/TCGA_HMTcorrs_heatmapexpression_combined_actsidebar.pdf",
    width = 4,
    height = 4)

heatmap.2(TCGA_correlations_within_HMTs[["allcombined"]][["correlations"]],
          col = mainpal,
          dendrogram = "none",
          Rowv = FALSE,
          Colv = FALSE,
          trace = "none",
          RowSideColors = TCGA_act_sidebar,
          density.info = "none",
          cexRow = 0.6,
          cexCol = 0.6,
          keysize = 1.8,
          key.title = NA)

dev.off()

pdf("graphics/GTEX_HMTcorrs_heatmapexpression_combined_TCGAorder.pdf",
    width = 4,
    height = 4)

heatmap.2(GTEX_correlations_within_HMTs$allcombined$correlations[row.names(TCGA_correlations_within_HMTs$allcombined$correlations), colnames(TCGA_correlations_within_HMTs$allcombined$correlations)],
col = mainpal,
dendrogram = "none",
Rowv = FALSE,
Colv = FALSE,
trace = "none",
RowSideColors = GTEX_all_sidebar[row.names(TCGA_correlations_within_HMTs$allcombined$correlations)],
density.info = "none",
cexRow = 0.6,
cexCol = 0.6,
keysize = 1.8,
key.title = NA)

dev.off()

# pdf("graphics/TCGA_HMTcorrs_heatmapexpression_combined_GTEXorder.pdf",
#     width = 4,
#     height = 4)
# 
# heatmap.2(TCGA_correlations_within_HMTs$allcombined[row.names(GTEX_correlations_within_HMTs$allcombined), colnames(GTEX_correlations_within_HMTs$allcombined)],
#           col = mainpal,
#           dendrogram = "none",
#           Rowv = FALSE,
#           Colv = FALSE,
#           trace = "none",
#           RowSideColors = TCGA_all_sidebar[row.names(GTEX_correlations_within_HMTs$allcombined)],
#           density.info = "none",
#           cexRow = 0.6,
#           cexCol = 0.6,
#           keysize = 1.8,
#           key.title = NA)
# 
# dev.off()

pdf("graphics/GTEX_HMTcorrelations_tissues_fullplots.pdf")

for (i in 2:length(GTEX_correlations_within_HMTs)){

  print(
        
    heatmap.2(GTEX_correlations_within_HMTs[[i]],
              col = mainpal,
              dendrogram = "none",
              Rowv = FALSE,
              Colv = FALSE,
              trace = "none",
              main = names(GTEX_correlations_within_HMTs)[i],
              RowSideColors = GTEX_all_sidebar,
              ColSideColors = GTEX_act_sidebar,
              density.info = "none",
              cexRow = 0.56,
              cexCol = 0.56,
              keysize = 1.8,
              key.title = NA)
        
        
        )
  
}

dev.off()

pdf("graphics/TCGA_HMTcorrelations_tissues_fullplots.pdf")

for (i in 2:length(TCGA_correlations_within_HMTs)){
  
  print(
    
    heatmap.2(TCGA_correlations_within_HMTs[[i]],
              col = mainpal,
              dendrogram = "none",
              Rowv = FALSE,
              Colv = FALSE,
              trace = "none",
              main = names(TCGA_correlations_within_HMTs)[i],
              RowSideColors = TCGA_all_sidebar,
              ColSideColors = TCGA_act_sidebar,
              density.info = "none",
              cexRow = 0.56,
              cexCol = 0.56,
              keysize = 1.8,
              key.title = NA)
    
    
  )
  
}

dev.off()

#### mean correlations with random genes for Fig S10B ####

# start with TCGA

TCGA_HMT_random_cors <- lapply(colnames(TCGA_correlations_within_HMTs$allcombined$correlations), function(thisHGNCname){
  
  thisensemblname <- SET_HMTs[SET_HMTs$hgnc_symbol == thisHGNCname, "ensembl_gene_id"]
  
  nuclear_notHMT <- nuclear_genes_expressed[!nuclear_genes_expressed %in% SET_HMTs$ensembl_gene_id]
  nuclear_notHMT_inTCGA <- nuclear_notHMT[nuclear_notHMT %in% row.names(TCGA_MOR_list[[1]])]
  
  randomgenes_noHMTs <- sample(nuclear_notHMT_inTCGA, 100)
    
    MOR_list_theseones <- lapply(TCGA_MOR_list, function(x){x[c(thisensemblname, randomgenes_noHMTs), ]})
    
    # MOR_combined <- do.call(cbind, MOR_list_HMTs)
    
    TCGA_LUT <- readRDS("~/NNMT_manuscript/output/TCGA_LUT.rds")
    
    residuals_by_cancer <- lapply(MOR_list_theseones, function(countdata){

      # restrict to samples present in lookuptable 
      tempcountdata <- countdata[, colnames(countdata) %in% TCGA_LUT$SAMPID]
      
      temp_df <- as.data.frame(matrix(nrow = ncol(tempcountdata), ncol = nrow(tempcountdata) + 5))
      
      colnames(temp_df) <- c(row.names(tempcountdata), "race", "gender", "tumour_stage", "sequencing_centre", "days_to_birth")
      row.names(temp_df) <- colnames(tempcountdata)
      
      # insert tranposed data into new data frame
      temp_df[, 1:nrow(tempcountdata)] <- t(tempcountdata)
      
      # from lookup table, add variables
      temp_df[, "days_to_birth"] <- as.numeric(TCGA_LUT[row.names(temp_df), "days_to_birth"])
      temp_df[, "gender"] <- TCGA_LUT[row.names(temp_df), "gender"]
      temp_df[, "race"] <- TCGA_LUT[row.names(temp_df), "race"]
      temp_df[, "tumour_stage"] <- TCGA_LUT[row.names(temp_df), "tumour_stage"]
      temp_df[, "sequencing_centre"] <- TCGA_LUT[row.names(temp_df), "sequencing_centre"]
      
      # check for any missing values and remove entry if so
      temp_df <- temp_df[!(is.na(temp_df$days_to_birth)), ]
      temp_df <- temp_df[!(is.na(temp_df$gender)), ]
      temp_df <- temp_df[!(is.na(temp_df$race)), ]
      temp_df <- temp_df[!(is.na(temp_df$tumour_stage)), ]
      temp_df <- temp_df[!(is.na(temp_df$sequencing_centre)), ]
      
      # do linear regression across all samples (all tissues combined)
      # perform linear regression with age (bracket midpoint), sex and cause of death
      
      # create dataframe to deposit residuals
      residuals_df <- as.data.frame(matrix(nrow = nrow(temp_df), ncol = ncol(temp_df) - 5))
      row.names(residuals_df) <- row.names(temp_df)
      colnames(residuals_df) <- colnames(temp_df)[1:(ncol(temp_df) - 5)]
      
      # loop over columns (genes); for each gene perform linear regression and deposit residuals in new data frame
      # sex is considered a factor automatically as a character string; sequencing centre is already considered as factor
      # here we note that in some projects, all sequencing centres are equal. Likewise gender for some cancer tpyes (e.g. ovarian)
      # need conditional to ignore seq centre, gender or tumour_stage if only one level exists among considered data in order to avoid an error
      
      if(length(unique(temp_df$race)) == 1){
        
        racevariable <- NULL
        
      } else {
        
        racevariable <- "race"
        
      }
      
      if(length(unique(temp_df$tumour_stage)) == 1){
        
        tumourstagevariable <- NULL
        
      } else {
        
        tumourstagevariable <- "tumour_stage"
        
      }
      
      if(length(unique(temp_df$gender)) == 1) {
        
        gendervariable <- NULL
        
      } else {
        
        gendervariable <- "gender"
      }
      
      if(length(unique(temp_df$sequencing_centre)) == 1) {
        
        seqcentrevariable <- NULL
        
      } else {
        
        seqcentrevariable <- "sequencing_centre"
        
      }
      
      modelvariables <- c("days_to_birth", gendervariable, seqcentrevariable, tumourstagevariable, racevariable)
      
      for(j in 1:(ncol(residuals_df))){
        
        outcome <- paste0("temp_df[, ", j, "]")
        
        f <- as.formula(
          paste(outcome,
                paste(modelvariables, collapse = " + "),
                sep = " ~ "))
        
        residuals_df[, colnames(temp_df)[j]] <- lm(formula = f, data = temp_df)$residuals
        
      } # end of residuals loop
      
      return(t(residuals_df))
      
    })
    
    residuals_df <- do.call(cbind, residuals_by_cancer)
  
    HMT_random_corrs <- sapply(randomgenes_noHMTs, function(thisrandom){

      cor.test(residuals_df[thisrandom,], residuals_df[thisensemblname,], method = "spearman")$estimate

    })
    
    return(HMT_random_corrs)

    })

TCGA_HMT_random_cors_forplot_list <- list(random = unlist(TCGA_HMT_random_cors), 
                                          all_HMTs = unlist(TCGA_correlations_within_HMTs$allcombined$correlations)[!unlist(TCGA_correlations_within_HMTs$allcombined$correlations) == 1],
                                          HE_clust = unlist(TCGA_correlations_within_HMTs$allcombined$correlations[TCGA_HEcluster, TCGA_HEcluster])[!unlist(TCGA_correlations_within_HMTs$allcombined$correlations[TCGA_HEcluster, TCGA_HEcluster]) == 1])

TCGA_HMT_random_cors_forplot_df <- reshape2::melt(TCGA_HMT_random_cors_forplot_list)
colnames(TCGA_HMT_random_cors_forplot_df) <- c("corr", "group")

TCGA_HMT_random_cors_forplot_df$group <- factor(TCGA_HMT_random_cors_forplot_df$group, levels = names(TCGA_HMT_random_cors_forplot_list))

TCGA_HMT_random_cors_plot_labs <- c("Random genes", "All HMTs", "High-expressing\ncluster")

saveRDS(TCGA_HMT_random_cors_forplot_df, "plot_data/TCGA_HMT_random_cors_forplot_df.rds")
TCGA_HMT_random_cors_forplot_df <- readRDS("plot_data/TCGA_HMT_random_cors_forplot_df.rds")

write.table(TCGA_HMT_random_cors_forplot_df,
            file = "plot_data/Fig S10/Fig_S10B_TCGA_plotdata.txt",
            row.names = FALSE,
            col.names = TRUE,
            sep = "\t")

pdf("graphics/TCGA_HMT_random_cors_plot.pdf",
    height = 2.5,
    width = 3)

ggplot(data = TCGA_HMT_random_cors_forplot_df, aes(x = group, y = corr)) + 
  geom_boxplot(outlier.size = 0.1,
               fill = c("grey", "dark turquoise", "red")) + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45,
                                   size = 8, 
                                   color = "black",
                                   hjust = 1),
        axis.text.y = element_text(size = 10,
                                   color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10),
        legend.position = "none") + 
  scale_y_continuous(breaks = seq(-0.5, 0.5, 0.5),
                     expand = c(0, 0),
                     limits = c(-0.75, 0.75)) +
  scale_x_discrete(labels = TCGA_HMT_random_cors_plot_labs) +
  geom_hline(yintercept = 0,
             linetype = "dashed",
             color = "grey") +
  ylab(substitute("Spearman's"~rho~"to HMTs"))

dev.off()

if(!exists("GTEX_MOR_across_tissues")){
  GTEX_MOR_across_tissues <- readRDS("output/GTEX_MOR_across_tissues.rds")
}

GTEX_HMT_random_cors <- lapply(colnames(GTEX_correlations_within_HMTs$allcombined$correlations), function(thisHGNCname){
# takes about 50 min per gene. total will take around 30h
  
  thisensemblname <- SET_HMTs[SET_HMTs$hgnc_symbol == thisHGNCname, "ensembl_gene_id"]
  
  randomgenes_noHMTs <- sample(nuclear_genes_expressed[!nuclear_genes_expressed %in% SET_HMTs$ensembl_gene_id], 100)

  # if(!exists("GTEX_MOR_across_tissues")){
  #   GTEX_MOR_across_tissues <- readRDS("output/GTEX_MOR_across_tissues.rds")
  # }
  
  MOR_theseones <- GTEX_MOR_across_tissues[c(thisensemblname, randomgenes_noHMTs), ]
  
  GTEX_LUT <- readRDS("output/GTEX-version8-sampleID-LUT.rds")
  
  GTEX_matrix_for_lm_temp <- data.frame(t(MOR_theseones))
    
  GTEX_matrix_for_lm_temp[, "age_bracket"] <- sapply(str_split(GTEX_LUT[row.names(GTEX_matrix_for_lm_temp), "age_bracket"], 
                                                                 pattern = "-"), function(x){
                                                                   
                                                                   mean(as.numeric(x))
                                                                   
                                                                 })    
    
    GTEX_matrix_for_lm_temp[, "sex"] <- GTEX_LUT[row.names(GTEX_matrix_for_lm_temp), "Sex"]
    GTEX_matrix_for_lm_temp[, "death"] <- GTEX_LUT[row.names(GTEX_matrix_for_lm_temp), "DeathScale"]
    GTEX_matrix_for_lm_temp[, "ischemic_time"] <- GTEX_LUT[row.names(GTEX_matrix_for_lm_temp), "IschemicTime"]
    
    
    # replace death = NA (unknown cause of death) with a 5th factor level
    GTEX_matrix_for_lm_temp[is.na(GTEX_matrix_for_lm_temp$death), "death"] <- 5
    
    # remove samples with NA for ischemic time
    GTEX_matrix_for_lm_temp <- GTEX_matrix_for_lm_temp[!is.na(GTEX_matrix_for_lm_temp$ischemic_time), ]
    
    # ischemic time rescaled so lmer doesn't complain about scale issues
    GTEX_matrix_for_lm_temp[, "ischemic_time"] <- GTEX_LUT[row.names(GTEX_matrix_for_lm_temp), "IschemicTime"]/60 
    
    GTEX_matrix_for_lm_temp[, "batch"] <- GTEX_LUT[row.names(GTEX_matrix_for_lm_temp), "Batch"]
    
    GTEX_matrix_for_lm_temp[, "tissue"] <- GTEX_LUT[row.names(GTEX_matrix_for_lm_temp), "Tissue"]
    
    GTEX_matrix_for_lm_temp[, "donorID"] <- GTEX_LUT[row.names(GTEX_matrix_for_lm_temp), "SUBJID"]
    
    # check for any missing values and remove entry if so
    
    GTEX_matrix_for_lm_temp <- GTEX_matrix_for_lm_temp[!(is.na(GTEX_matrix_for_lm_temp$ischemic_time)), ]
    
    # do linear regression across all samples (all tissues combined)
    # perform linear regression with age (bracket midpoint), sex and cause of death
    
    # create dataframe to deposit residuals
    residuals_df <- as.data.frame(matrix(nrow = nrow(GTEX_matrix_for_lm_temp), ncol = ncol(GTEX_matrix_for_lm_temp) - 7))
    row.names(residuals_df) <- row.names(GTEX_matrix_for_lm_temp)
    colnames(residuals_df) <- colnames(GTEX_matrix_for_lm_temp)[1:(ncol(GTEX_matrix_for_lm_temp) - 7)]
    
    # loop over columns (genes); for each gene perform linear regression and deposit residuals in new data frame
    # sex is considered a factor automatically as a character string; as.factor() must be specified for cause of death
    
    for(j in 1:(ncol(residuals_df))){
      
      residuals_df[, colnames(GTEX_matrix_for_lm_temp)[j]] <- resid(lmer(GTEX_matrix_for_lm_temp[, j] ~ (sex + as.factor(death) + age_bracket + ischemic_time)*tissue + (1|batch) + (1|donorID), data = GTEX_matrix_for_lm_temp))
      
    }
    
    residuals_df <- data.frame(t(residuals_df))

  HMT_random_corrs <- sapply(randomgenes_noHMTs, function(thisrandom){

    cor.test(unlist(residuals_df[thisrandom,]), unlist(residuals_df[thisensemblname,]), method = "spearman")$estimate
    
  })
  
  return(HMT_random_corrs)
  
})

names(GTEX_HMT_random_cors) <- colnames(GTEX_correlations_within_HMTs$allcombined$correlations)
saveRDS(GTEX_HMT_random_cors, "output/GTEX_HMT_random_cors.rds")

# here exploratory boxplot. random corrs vs all HMT-HMT corrs vs high expressed (here defined as 20 with expression > 10^3.25)
GTEX_HMT_random_cors_forplot_list <- list(random = unlist(GTEX_HMT_random_cors), 
        all_HMTs = unlist(GTEX_correlations_within_HMTs$allcombined$correlations)[!unlist(GTEX_correlations_within_HMTs$allcombined$correlations) == 1],
        HE_clust = unlist(GTEX_correlations_within_HMTs$allcombined$correlations[GTEX_HEcluster, GTEX_HEcluster])[!unlist(GTEX_correlations_within_HMTs$allcombined$correlations[GTEX_HEcluster, GTEX_HEcluster]) == 1])

GTEX_HMT_random_cors_forplot_df <- reshape2::melt(GTEX_HMT_random_cors_forplot_list)
colnames(GTEX_HMT_random_cors_forplot_df) <- c("corr", "group")

GTEX_HMT_random_cors_forplot_df$group <- factor(GTEX_HMT_random_cors_forplot_df$group, levels = names(GTEX_HMT_random_cors_forplot_list))

GTEX_HMT_random_cors_plot_labs <- c("Random genes", "All HMTs", "High-expressing\ncluster")

saveRDS(GTEX_HMT_random_cors_forplot_df, "plot_data/GTEX_HMT_random_cors_forplot_df.rds")
GTEX_HMT_random_cors_forplot_df <- readRDS("plot_data/GTEX_HMT_random_cors_forplot_df.rds")

write.table(GTEX_HMT_random_cors_forplot_df,
            file = "plot_data/Fig S10/Fig_S10B_GTEX_plotdata.txt",
            row.names = FALSE,
            col.names = TRUE,
            sep = "\t")

pdf("graphics/GTEX_HMT_random_cors_plot.pdf",
    height = 2.5,
    width = 3)

ggplot(data = GTEX_HMT_random_cors_forplot_df, aes(x = group, y = corr)) + 
  geom_boxplot(outlier.size = 0.1,
               fill = c("grey", "dark turquoise", "red")) + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45,
                                   size = 8, 
                                   color = "black",
                                   hjust = 1),
        axis.text.y = element_text(size = 10,
                                   color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10),
        legend.position = "none") + 
  scale_y_continuous(breaks = seq(-0.5, 0.5, 0.5),
                     expand = c(0, 0),
                     limits = c(-0.75, 0.75)) +
  scale_x_discrete(labels = GTEX_HMT_random_cors_plot_labs) +
  geom_hline(yintercept = 0,
             linetype = "dashed",
             color = "grey") +
  ylab(substitute("Spearman's"~rho~"to HMTs"))
  
dev.off()

#### HMT CORRELATION NETWORK PLOT PREPARATION for S10C ####

## GTEX/TCGA correlation between edges and jaccard index #####

GTEX_HMThm_genes <- colnames(GTEX_correlations_within_HMTs$allcombined$correlations)
TCGA_HMThm_genes <- colnames(TCGA_correlations_within_HMTs$allcombined$correlations)

gene_combos <- data.frame(t(combn(TCGA_HMThm_genes, 2)))

gene_combos_GTEX <- gene_combos
gene_combos_TCGA <- gene_combos

gene_combos_GTEX[, "correlation"] <- apply(gene_combos, 1, function(x){

  GTEX_correlations_within_HMTs$allcombined$correlations[x[1], x[2]]
  
})

gene_combos_TCGA[, "correlation"] <- apply(gene_combos, 1, function(x){
  
  TCGA_correlations_within_HMTs$allcombined$correlations[x[1], x[2]]
  
})

cor.test(gene_combos_GTEX[, "correlation"], 
         gene_combos_TCGA[, "correlation"])

gene_combos_GTEX_cut <- apply(gene_combos_GTEX[abs(gene_combos_GTEX$correlation) > 0.2, c(1,2)], 1, function(x){
  
  gene1 <- x[1]
  gene2 <- x[2]
  
  gene_vec <- c(gene1, gene2)
  gene_vec <- gene_vec[order(gene_vec)]
  
  paste0(gene_vec[1],"_",gene_vec[2] )
  
})

gene_combos_TCGA_cut <- apply(gene_combos_TCGA[abs(gene_combos_TCGA$correlation) > 0.2, c(1,2)], 1, function(x){
  
  gene1 <- x[1]
  gene2 <- x[2]
  
  gene_vec <- c(gene1, gene2)
  gene_vec <- gene_vec[order(gene_vec)]
  
  paste0(gene_vec[1],"_",gene_vec[2] )
  
})

jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}

jaccard(gene_combos_GTEX_cut, gene_combos_TCGA_cut)

#### prepare .gexf file for Gephi import for gtex ####

# Keep only high correlations
GTEX_pan_fornetwork <- GTEX_correlations_within_HMTs$allcombined$correlations

GTEX_pan_fornetwork[GTEX_correlations_within_HMTs$allcombined$FDR > 0.05] <- 0

# Make an Igraph object from this matrix. Use simplify to ensure that there are no duplicated edges or self loops
GTEX_pan_HMT_network <- simplify(graph_from_adjacency_matrix(abs(GTEX_pan_fornetwork),
                                                             weighted = TRUE,
                                                             mode = "undirected",
                                                             diag = FALSE))

dataset.ext.base <- simplify(graph_from_adjacency_matrix(abs(GTEX_pan_fornetwork),
                                                         weighted = TRUE,
                                                         mode = "undirected",
                                                         diag = FALSE))
dataset.ext.base[]
dataset.ext.base <- delete.edges(dataset.ext.base, E(dataset.ext.base)[abs(weight) < 0.2 ])

dataSet.ext <- ddply(as_data_frame(dataset.ext.base), .variables=c("from", "to", "weight"))

dataSet.ext[, "rawweight"] <- apply(dataSet.ext, 1, function(x){GTEX_pan_fornetwork[x[1], x[2]]})
dataSet.ext[, "scaleweight"] <- 0.5 + 0.5*(dataSet.ext$rawweight)

# cull some edges
GTEX_pan_HMT_network <- delete.edges(GTEX_pan_HMT_network, E(GTEX_pan_HMT_network)[ abs(weight) < 0.2 ])

# Calculate degree for all nodes
degAll <- degree(GTEX_pan_HMT_network, v = V(GTEX_pan_HMT_network), mode = "all")

#Calculate Dice similarities between all pairs of nodes
dsAll <- similarity.dice(GTEX_pan_HMT_network, vids = V(GTEX_pan_HMT_network), mode = "all")

GTEX_pan_HMT_network <- set.vertex.attribute(GTEX_pan_HMT_network, "degree", index = V(GTEX_pan_HMT_network), value = degAll)

GTEX_pan_HMT_network <- set.edge.attribute(GTEX_pan_HMT_network, "weight", index = E(GTEX_pan_HMT_network), value = 0)
GTEX_pan_HMT_network <- set.edge.attribute(GTEX_pan_HMT_network, "scaleweight", index = E(GTEX_pan_HMT_network), value = 0)

# The order of interactions in GTEX_pan_HMT_network is not the same as it is in dataSet or as it is in the edge list,
# and for that reason these values cannot be assigned directly

E(GTEX_pan_HMT_network)[dataSet.ext$from %--% dataSet.ext$to]$weight <- as.numeric(dataSet.ext$weight)
E(GTEX_pan_HMT_network)[dataSet.ext$from %--% dataSet.ext$to]$scaleweight <- as.numeric(dataSet.ext$scaleweight)

# Create a dataframe nodes: 1st column - node ID, 2nd column -node name
nodes_df <- data.frame(ID = c(1:vcount(GTEX_pan_HMT_network)),
                       NAME = V(GTEX_pan_HMT_network)$name)

# Create a dataframe edges: 1st column - source node ID, 2nd column -target node ID
edges_df <- as.data.frame(get.edges(GTEX_pan_HMT_network, c(1:ecount(GTEX_pan_HMT_network))))

# Define node and edge attributes - these attributes won't be directly used for network visualization, but they
# may be useful for other network manipulations in Gephi

# Create a dataframe with node attributes: 1st column - attribute 1 (degree), 2nd column - attribute 2 (betweenness)
nodes_att <- data.frame(DEG = V(GTEX_pan_HMT_network)$degree) 

# Create a dataframe with edge attributes: 1st column - attribute 1 (weight), 2nd column - attribute 2 (similarity)
edges_att <- data.frame(WGH = E(GTEX_pan_HMT_network)$weight,
                        SCLWGH = E(GTEX_pan_HMT_network)$scaleweight) 

# Define node/edge visual attributes - these attributes are the ones used for network visualization

# Calculate node coordinate - needs to be 3D

#nodes_coord <- as.data.frame(layout.fruchterman.reingold(GTEX_pan_HMT_network, weights = E(GTEX_pan_HMT_network)$similarity, dim = 3, niter = 10000))
# We'll cheat here, as 2D coordinates result in a better (2D) plot than 3D coordinates
nodes_coord <- as.data.frame(layout.fruchterman.reingold(GTEX_pan_HMT_network, weights = E(GTEX_pan_HMT_network)$scaleweight, dim = 2, niter = 10000))
nodes_coord <- cbind(nodes_coord, rep(0, times = nrow(nodes_coord)))

# Calculate node size
# We'll interpolate node size based on the node degree, using the "approx" function
approxVals <- approx(c(1, 5), n = length(unique(V(GTEX_pan_HMT_network)$degree)))
# And we will assign a node size for each node based on its degree
nodes_size <- sapply(V(GTEX_pan_HMT_network)$degree, function(x) approxVals$y[which(sort(unique(V(GTEX_pan_HMT_network)$degree)) == x)])

# Define node color. Activating / repressing colours

GTEX_HMT_pan_Vactivity <- SET_HMTs_annotated[match(names(V(GTEX_pan_HMT_network)), SET_HMTs_annotated$hgnc_symbol), "Transcription"]
GTEX_HMT_pan_Vactivity[GTEX_HMT_pan_Vactivity == "repressing"] <- "magenta"
GTEX_HMT_pan_Vactivity[GTEX_HMT_pan_Vactivity == "activating"] <- "green"
GTEX_HMT_pan_Vactivity[GTEX_HMT_pan_Vactivity == "unclear"] <- "grey"

nodes_col_df <- as.data.frame(t(col2rgb(GTEX_HMT_pan_Vactivity, alpha = FALSE)))

nodes_att_viz <- list(color = nodes_col_df,
                      position = nodes_coord,
                      size = nodes_size)

# Assign visual attributes to edges using the same approach as we did for nodes
# 
# edges_col <- map2color(E(GTEX_pan_HMT_network)$scaleweight,
#                        mainpal,
#                        limits = c(0,1))

edges_col_vec <- vector(length = length(E(GTEX_pan_HMT_network)$scaleweight))
edges_col_vec[E(GTEX_pan_HMT_network)$scaleweight > 0.5] <- "#FF0000"
edges_col_vec[E(GTEX_pan_HMT_network)$scaleweight < 0.5] <- "#0000FF"

edges_col_df <- as.data.frame(t(col2rgb(edges_col_vec, alpha = FALSE)))
edges_col_df <- cbind(edges_col_df, alpha = rep(1, times = nrow(edges_col_df)))
edges_att_viz <-list(color = edges_col_df)

# Write the network into a gexf (Gephi) file
write.gexf(nodes = nodes_df, edges = edges_df, nodesAtt = nodes_att, edgesWeight = E(GTEX_pan_HMT_network)$weight^2, edgesAtt = edges_att, nodesVizAtt = nodes_att_viz, edgesVizAtt = edges_att_viz, defaultedgetype = "undirected", output = "output/GTEX_HMT_pantissue.gexf")

#### prepare .gexf file for Gephi import for tcga ####

# Keep only high correlations
TCGA_pan_fornetwork <- TCGA_correlations_within_HMTs$allcombined$correlations

TCGA_pan_fornetwork[TCGA_correlations_within_HMTs$allcombined$FDR > 0.05] <- 0

# Make an Igraph object from this matrix. Use simplify to ensure that there are no duplicated edges or self loops
TCGA_pan_HMT_network <- simplify(graph_from_adjacency_matrix(abs(TCGA_pan_fornetwork),
                                                             weighted = TRUE,
                                                             mode = "undirected",
                                                             diag = FALSE))

dataset.ext.base <- simplify(graph_from_adjacency_matrix(abs(TCGA_pan_fornetwork),
                                                         weighted = TRUE,
                                                         mode = "undirected",
                                                         diag = FALSE))

dataset.ext.base <- delete.edges(dataset.ext.base, E(dataset.ext.base)[abs(weight) < 0.2 ])

dataSet.ext <- ddply(as_data_frame(dataset.ext.base), .variables=c("from", "to", "weight"))

dataSet.ext[, "rawweight"] <- apply(dataSet.ext, 1, function(x){TCGA_pan_fornetwork[x[1], x[2]]})
dataSet.ext[, "scaleweight"] <- 0.5 + 0.5*(dataSet.ext$rawweight)

# cull some edges
TCGA_pan_HMT_network <- delete.edges(TCGA_pan_HMT_network, E(TCGA_pan_HMT_network)[ abs(weight) < 0.2 ])

# Calculate degree for all nodes
degAll <- degree(TCGA_pan_HMT_network, v = V(TCGA_pan_HMT_network), mode = "all")

#Calculate Dice similarities between all pairs of nodes
dsAll <- similarity.dice(TCGA_pan_HMT_network, vids = V(TCGA_pan_HMT_network), mode = "all")

TCGA_pan_HMT_network <- set.vertex.attribute(TCGA_pan_HMT_network, "degree", index = V(TCGA_pan_HMT_network), value = degAll)

TCGA_pan_HMT_network <- set.edge.attribute(TCGA_pan_HMT_network, "weight", index = E(TCGA_pan_HMT_network), value = 0)
TCGA_pan_HMT_network <- set.edge.attribute(TCGA_pan_HMT_network, "scaleweight", index = E(TCGA_pan_HMT_network), value = 0)

# The order of interactions in TCGA_pan_HMT_network is not the same as it is in dataSet or as it is in the edge list,
# and for that reason these values cannot be assigned directly

E(TCGA_pan_HMT_network)[dataSet.ext$from %--% dataSet.ext$to]$weight <- as.numeric(dataSet.ext$weight)
E(TCGA_pan_HMT_network)[dataSet.ext$from %--% dataSet.ext$to]$scaleweight <- as.numeric(dataSet.ext$scaleweight)

# Create a dataframe nodes: 1st column - node ID, 2nd column -node name
nodes_df <- data.frame(ID = c(1:vcount(TCGA_pan_HMT_network)),
                       NAME = V(TCGA_pan_HMT_network)$name)

# Create a dataframe edges: 1st column - source node ID, 2nd column -target node ID
edges_df <- as.data.frame(get.edges(TCGA_pan_HMT_network, c(1:ecount(TCGA_pan_HMT_network))))

# Define node and edge attributes - these attributes won't be directly used for network visualization, but they
# may be useful for other network manipulations in Gephi

# Create a dataframe with node attributes: 1st column - attribute 1 (degree), 2nd column - attribute 2 (betweenness)
nodes_att <- data.frame(DEG = V(TCGA_pan_HMT_network)$degree) 

# Create a dataframe with edge attributes: 1st column - attribute 1 (weight), 2nd column - attribute 2 (similarity)
edges_att <- data.frame(WGH = E(TCGA_pan_HMT_network)$weight,
                        SCLWGH = E(TCGA_pan_HMT_network)$scaleweight) 

# Define node/edge visual attributes - these attributes are the ones used for network visualization

# Calculate node coordinate - needs to be 3D

#nodes_coord <- as.data.frame(layout.fruchterman.reingold(TCGA_pan_HMT_network, weights = E(TCGA_pan_HMT_network)$similarity, dim = 3, niter = 10000))
# We'll cheat here, as 2D coordinates result in a better (2D) plot than 3D coordinates
nodes_coord <- as.data.frame(layout.fruchterman.reingold(TCGA_pan_HMT_network, weights = E(TCGA_pan_HMT_network)$scaleweight, dim = 2, niter = 10000))
nodes_coord <- cbind(nodes_coord, rep(0, times = nrow(nodes_coord)))

# Calculate node size
# We'll interpolate node size based on the node degree, using the "approx" function
approxVals <- approx(c(1, 5), n = length(unique(V(TCGA_pan_HMT_network)$degree)))
# And we will assign a node size for each node based on its degree
nodes_size <- sapply(V(TCGA_pan_HMT_network)$degree, function(x) approxVals$y[which(sort(unique(V(TCGA_pan_HMT_network)$degree)) == x)])

# Define node color. Activating / repressing colours

TCGA_HMT_pan_Vactivity <- SET_HMTs_annotated[match(names(V(TCGA_pan_HMT_network)), SET_HMTs_annotated$hgnc_symbol), "Transcription"]
TCGA_HMT_pan_Vactivity[TCGA_HMT_pan_Vactivity == "repressing"] <- "magenta"
TCGA_HMT_pan_Vactivity[TCGA_HMT_pan_Vactivity == "activating"] <- "green"
TCGA_HMT_pan_Vactivity[TCGA_HMT_pan_Vactivity == "unclear"] <- "grey"

nodes_col_df <- as.data.frame(t(col2rgb(TCGA_HMT_pan_Vactivity, alpha = FALSE)))

nodes_att_viz <- list(color = nodes_col_df,
                      position = nodes_coord,
                      size = nodes_size)

# Assign visual attributes to edges using the same approach as we did for nodes

edges_col <- map2color(E(TCGA_pan_HMT_network)$scaleweight,
                       mainpal,
                       limits = c(0,1))

edges_col_df <- as.data.frame(t(col2rgb(edges_col, alpha = FALSE)))
edges_col_df <- cbind(edges_col_df, alpha = rep(1, times = nrow(edges_col_df)))
edges_att_viz <-list(color = edges_col_df)

# Write the network into a gexf (Gephi) file
write.gexf(nodes = nodes_df, edges = edges_df, nodesAtt = nodes_att, edgesWeight = E(TCGA_pan_HMT_network)$weight^2, edgesAtt = edges_att, nodesVizAtt = nodes_att_viz, edgesVizAtt = edges_att_viz, defaultedgetype = "undirected", output = "output/TCGA_HMT_pantissue.gexf")

#### LIN-35 ChIP in C. elegans Fig 4B ####

# get data from Table S1 of Goetsch et al. 2017, PLoS Genetics [https://doi.org/10.1371/journal.pgen.1007088]
download.file("https://doi.org/10.1371/journal.pgen.1007088.s011", "output/Goetsch_TS1.xlsx")

lin35chip <- read_excel("output/Goetsch_TS1.xlsx", sheet = 5)

lin35GR <- makeGRangesFromDataFrame(lin35chip)

## strandedness

# in + strand molecules, the TSS is at the 'start'
Cel_HMTSS_plus <- CelHMTGR[strand(CelHMTGR) == "+"]
end(Cel_HMTSS_plus) <- start(Cel_HMTSS_plus)

# in - strand molecules, the TSS is at the 'end'
Cel_HMTSS_minus <- CelHMTGR[strand(CelHMTGR) == "-"]
start(Cel_HMTSS_minus) <- end(Cel_HMTSS_minus)

# combine the TSSs
CeLHMTTSS <- c(Cel_HMTSS_plus, Cel_HMTSS_minus)
# CeLHMTTSS <- CeLHMTTSS[match(Cel_HMT_locstrand$wormbase_locus, mcols(CeLHMTTSS)$wormbase_locus)]

lin35distances <- distanceToNearest(CeLHMTTSS, lin35GR)
mcols(lin35distances)[, "wormbase_locus"] <- mcols(CeLHMTTSS)$wormbase_locus
# Cel_HMT_locstrand[, "distancetolin35peaks"] <- mcols(lin35distances)$distance
# 
# boxplot(Cel_HMT_locstrand[Cel_HMT_locstrand$wormbase_locus %in% Cel_HE, ]$distancetolin35peaks,
# Cel_HMT_locstrand[!Cel_HMT_locstrand$wormbase_locus %in% Cel_HE, ]$distancetolin35peaks)

# Random genes
if(!exists("CeNDR_normalised_counts_collapse")){
  CeNDR_normalised_counts_collapse <- readRDS("output/Cel_CeNDR_normalised_counts_gene_level.rds")
}

Cel_normalised_expressed <- CeNDR_normalised_counts_collapse[apply(CeNDR_normalised_counts_collapse, 1, function(x){!any(x == 0)}), ]

# gives 11177 genes expressed in all samples
# take random sample of consistently expressed genes
Cel_genes_randomsample <- sample(row.names(Cel_normalised_expressed), 1000)
Cel_genes_randomsampleBM <- getBM(mart = parasite_mart,
                        attributes = c("strand", 
                                       "chromosome_name",
                                       "start_position", 
                                       "end_position", 
                                       "wormbase_locus",
                                       "wormbase_gseq"),
                        values = Cel_genes_randomsample,
                        filters = "wormbase_gseqname")

Cel_genes_randomsampleBM[, "chromosome_name"] <- paste0("chr", Cel_genes_randomsampleBM$chromosome_name)

# need to make colnames 'start' and 'end' for GR - remove '_position'
colnames(Cel_genes_randomsampleBM) <- str_remove(colnames(Cel_genes_randomsampleBM), pattern = "_position")

# convert strand to recognise by GR
Cel_genes_randomsampleBM[Cel_genes_randomsampleBM$strand == 1, "strand"] <- "+"
Cel_genes_randomsampleBM[Cel_genes_randomsampleBM$strand == -1, "strand"] <- "-"

Cel_genes_randomsampleGR <- makeGRangesFromDataFrame(Cel_genes_randomsampleBM, keep.extra.columns = TRUE)

Cel_genes_randomsampleGR_X <- Cel_genes_randomsampleGR[seqnames(Cel_genes_randomsampleGR) == "chrX"]
Cel_genes_randomsampleGR_autosomal <- Cel_genes_randomsampleGR[seqnames(Cel_genes_randomsampleGR) != "chrX"]

Cel_genes_randomsampleXTSS_plus <- Cel_genes_randomsampleGR_X[strand(Cel_genes_randomsampleGR_X) == "+"]
end(Cel_genes_randomsampleXTSS_plus) <- start(Cel_genes_randomsampleXTSS_plus)
Cel_genes_randomsampleXTSS_minus <- Cel_genes_randomsampleGR_X[strand(Cel_genes_randomsampleGR_X) == "-"]
start(Cel_genes_randomsampleXTSS_minus) <- end(Cel_genes_randomsampleXTSS_minus)
Cel_genes_randomsampleXTSS <- c(Cel_genes_randomsampleXTSS_plus, Cel_genes_randomsampleXTSS_minus)

Cel_genes_randomsampleautoTSS_plus <- Cel_genes_randomsampleGR_autosomal[strand(Cel_genes_randomsampleGR_autosomal) == "+"]
end(Cel_genes_randomsampleautoTSS_plus) <- start(Cel_genes_randomsampleautoTSS_plus)
Cel_genes_randomsampleautoTSS_minus <- Cel_genes_randomsampleGR_autosomal[strand(Cel_genes_randomsampleGR_autosomal) == "-"]
start(Cel_genes_randomsampleautoTSS_minus) <- end(Cel_genes_randomsampleautoTSS_minus)
Cel_genes_randomsampleautoTSS <- c(Cel_genes_randomsampleautoTSS_plus, Cel_genes_randomsampleautoTSS_minus)

autosomes_lin35distance <- mcols(lin35distances)$distance[as.vector(seqnames(CeLHMTTSS)) != "chrX"]
autosomes_lin35distance <- mcols(lin35distances)$distance[as.vector(seqnames(CeLHMTTSS)) != "chrX"]


Xchr_lin35distance <- mcols(lin35distances)$distance[as.vector(seqnames(CeLHMTTSS)) == "chrX"]

Cel_genes_random_lin35distance <- distanceToNearest(Cel_genes_randomsampleautoTSS, lin35GR)

Cel_genes_random_autosomes_lin35distance <- mcols(Cel_genes_random_lin35distance)$distance[as.vector(seqnames(CeLHMTTSS)) != "chrX"]
Cel_genes_random_Xchr_lin35distance <- mcols(Cel_genes_random_lin35distance)$distance[as.vector(seqnames(CeLHMTTSS)) == "chrX"]

wilcox.test(Cel_genes_random_autosomes_lin35distance, autosomes_lin35distance)
wilcox.test(Cel_genes_random_Xchr_lin35distance, Xchr_lin35distance)

lin35_distances_for_plot <- reshape2::melt(list("random genes" = c(Cel_genes_random_autosomes_lin35distance, Cel_genes_random_Xchr_lin35distance),
                                 "other HMTs" = mcols(lin35distances)$distance[!mcols(lin35distances)$wormbase_locus %in% Cel_HE],
                                 "HE cluster" = mcols(lin35distances)$distance[mcols(lin35distances)$wormbase_locus %in% Cel_HE]))

lin35_distances_for_plot_copy <- lin35_distances_for_plot
colnames(lin35_distances_for_plot_copy) <- c("lin35_TSSdistance", "group")

write.table(lin35_distances_for_plot_copy,
            file = "plot_data/Fig 4/Fig_4B_plotdata.txt",
            row.names = FALSE,
            col.names = TRUE,
            sep = "\t")

scientific_10 <- function(x) {parse(text=gsub("e\\+*", " %*% 10^", scales::scientific_format()(x)))}
options(scipen = 3)

pdf("graphics/LIN35_binding_distance.pdf",
    height = 2.5,
    width = 3.5)

ggplot(lin35_distances_for_plot, aes(x = value, y = L1, fill = L1)) +
  geom_violin() +
  geom_jitter(width = 0.05,
  height = 0.05, alpha = 0.4) +
  scale_x_continuous(trans = scales::pseudo_log_trans(base = 10),
                     breaks = c(0, 10, 100, 1000, 10000, 100000, 1000000),
                     limits = c(-1, 1000000),
                     labels = trans_format("log10", math_format(10^.x))) +
  theme_classic() + 
  theme(axis.text.x = element_text(colour = "black",
                                   size = 10),
        axis.text.y = element_text(colour = "black",
                                   size = 10),
        legend.position = "none") +
  scale_fill_manual(values = c("red", "blue", "grey")) +
  ylab(NULL) +
  xlab("LIN-35 peak distance to TSS (bp)")
  
dev.off()

wilcox.test(lin35_distances_for_plot[lin35_distances_for_plot$L1 == "other HMTs", "value"],
            lin35_distances_for_plot[lin35_distances_for_plot$L1 == "HE cluster", "value"])

wilcox.test(lin35_distances_for_plot[lin35_distances_for_plot$L1 == "random genes", "value"],
            lin35_distances_for_plot[lin35_distances_for_plot$L1 == "HE cluster", "value"])

wilcox.test(lin35_distances_for_plot[lin35_distances_for_plot$L1 == "random genes", "value"],
            lin35_distances_for_plot[lin35_distances_for_plot$L1 == "other HMTs", "value"])

#### LIN-35 mutations ####

# elsewhere downloaded and realigned raw RNA-seq data from Latorre et al. 2015 [doi:10.1101/gad.255810.114] for lin-35 and N2 L3s and
# Gal et al. 2021 [doi: 10.1016/j.celrep.2021.109835] for N2 and lin-35 L1s
lin35RNAseq_L3 <- readRDS("input/LIN35_normalised_counts.rds")
lin35RNAseq_L1 <- readRDS("input/Petrella_LIN35_normalised_counts.rds")

convertBML3 <- getBM(mart = parasite_mart,
                     values = row.names(lin35RNAseq_L3),
                     filters = "entrezgene_id",
                     attributes = c("entrezgene_id",
                                    "wormbase_locus",
                                    "wormbase_gseq"))

row.names(lin35RNAseq_L3) <- convertBML3[match(row.names(lin35RNAseq_L3), convertBML3$entrezgene_id), "wormbase_gseq"]

convertBML1 <- getBM(mart = parasite_mart,
                     values = row.names(lin35RNAseq_L1),
                     filters = "entrezgene_id",
                     attributes = c("entrezgene_id",
                                    "wormbase_locus",
                                    "wormbase_gseq"))

row.names(lin35RNAseq_L1) <- convertBML1[match(row.names(lin35RNAseq_L1), convertBML1$entrezgene_id), "wormbase_gseq"]

lin35RNAseq_L3 <- lin35RNAseq_L3[!is.na(row.names(lin35RNAseq_L3)), ]
lin35RNAseq_L1 <- lin35RNAseq_L1[!is.na(row.names(lin35RNAseq_L1)), ]

colnames(lin35RNAseq_L3) <- c("wt_r1",
                              "wt_r2",
                              "lin35_r1",
                              "lin35_r2")

Petrella_runinfo <- read.table("input/PetrellaSraRunInfo.txt",
                               sep = ",",
                               header = TRUE)

colnames(lin35RNAseq_L1) <- str_replace_all(paste(Petrella_runinfo[match(str_remove(colnames(lin35RNAseq_L1), pattern = "_sorted.bam"), Petrella_runinfo$Run), "STRAIN"],
                                                  Petrella_runinfo[match(str_remove(colnames(lin35RNAseq_L1), pattern = "_sorted.bam"), Petrella_runinfo$Run), "growth_temperature"], sep = "_"), pattern = " ", replacement = "_")

lin35RNAseq_L1 <- data.frame(lin35RNAseq_L1)

lin35RNAseq_L1_HMT <- colSums(lin35RNAseq_L1[Cel_HMTgenes$wormbase_gseq, ])
lin35RNAseq_L1_HMT <- lin35RNAseq_L1_HMT[c(1:3, 6:9)]

lin35RNAseq_L3_HMT <- colSums(lin35RNAseq_L3[Cel_HMTgenes$wormbase_gseq, ])

foranova_df <- data.frame(HMTs = c(lin35RNAseq_L1_HMT, lin35RNAseq_L3_HMT),
                          stage = c("L1", "L1", "L1", "L1", "L1", "L1", "L1", "L3", "L3", "L3", "L3"))

foranova_df[str_detect(row.names(foranova_df), "N2"), "genotype"] <- "wt"
foranova_df[str_detect(row.names(foranova_df), "wt"), "genotype"] <- "wt"

foranova_df[str_detect(row.names(foranova_df), "lin"), "genotype"] <- "lin35"

my_anova <- aov(log10(HMTs) ~ genotype + stage, data = foranova_df)

# with unbalanced anova design, the order of the factors changes the p value. The appropriate p value is the type II, or in other words, the one where your factor of interest comes second.
Anova(my_anova, type = "II")

# % increase expression
mean(foranova_df[foranova_df$stage == "L1" & foranova_df$genotype == "lin35", ]$HMTs) / mean(foranova_df[foranova_df$stage == "L1" & foranova_df$genotype == "wt", ]$HMTs)
mean(foranova_df[foranova_df$stage == "L3" & foranova_df$genotype == "lin35", ]$HMTs) / mean(foranova_df[foranova_df$stage == "L3" & foranova_df$genotype == "wt", ]$HMTs)

#### Retinoblastoma MUTANT CANCERS TCGA Fig 4D ####

# cancer mafs downloaded with TCGAbiolinks using cluster
cancermaflist <- readRDS("input/cancermaflist.rds")

unique_samples_present <- unique(str_remove(unique(unlist(sapply(cancermaflist, function(x){x$Tumor_Sample_Barcode}))), pattern = "-[:alnum:]+-[:alnum:]+-[:digit:]+$"))

if(!exists("TCGA_LUT")){
  
  TCGA_LUT <- readRDS("~/NNMT_manuscript/output/TCGA_LUT.rds")

}

# subset the TCGA lookuptable according to samples present in the mutation calls
TCGA_LUT_maf <- TCGA_LUT[str_remove(TCGA_LUT[, "SAMPID"], pattern = "-[:alnum:]+-[:alnum:]+-[:digit:]+$") %in% unique_samples_present, ]

# find the mutations in the right place
# looking for mutations that might be deleterious. "Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Ins", "In_Frame_Del", "Frame_Shift_Del"
maf_in_LUT_RB1 <- lapply(cancermaflist, function(x){
  
  x[(str_remove(x$Tumor_Sample_Barcode, pattern = "-[:alnum:]+-[:alnum:]+-[:digit:]+$") %in% str_remove(TCGA_LUT_maf[, "SAMPID"], pattern = "-[:alnum:]+-[:alnum:]+-[:digit:]+$")) & (x$Hugo_Symbol == "RB1") & (x$Variant_Classification %in% c("Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Ins", "In_Frame_Del", "Frame_Shift_Del")), ]
  
})

# find samples with mutation calls but none called in RB1
samples_in_LUT <- unique(unlist(sapply(cancermaflist, function(x){x[str_remove(x$Tumor_Sample_Barcode, pattern = "-[:alnum:]+-[:alnum:]+-[:digit:]+$") %in% str_remove(TCGA_LUT_maf$SAMPID,pattern = "-[:alnum:]+-[:alnum:]+-[:digit:]+$"), "Tumor_Sample_Barcode"]})))
samples_in_LUT_notRBmut <- samples_in_LUT[!(samples_in_LUT %in% unique(unlist(sapply(maf_in_LUT_RB1, function(x){x$Tumor_Sample_Barcode}))))]

# produce table of number of samples in each cancer with no Rb mutations called
wtRB_bycancer <- table(TCGA_LUT_maf[str_remove(TCGA_LUT_maf$SAMPID, pattern = "-[:alnum:]+-[:alnum:]+-[:digit:]+$") %in% str_remove(samples_in_LUT_notRBmut, pattern = "-[:alnum:]+-[:alnum:]+-[:digit:]+$"), "cancer"])

Rbdeleterious_bycancer <- data.frame(Rb_mutant = sapply(maf_in_LUT_RB1, function(x){length(unique(x$Tumor_Sample_Barcode))}),
                                     Rb_wildtype = wtRB_bycancer[names(sapply(maf_in_LUT_RB1, function(x){length(unique(x$Tumor_Sample_Barcode))}))])

Rbdeleterious_SAMPIDs_bycancer <- lapply(maf_in_LUT_RB1, function(x){unique(x$Tumor_Sample_Barcode)})
wtRB_SAMPIDs_bycancer <- lapply(unique(TCGA_LUT$cancer), function(x){TCGA_LUT_maf[(str_remove(TCGA_LUT_maf$SAMPID, pattern = "-[:alnum:]+-[:alnum:]+-[:digit:]+$") %in% str_remove(samples_in_LUT_notRBmut, pattern = "-[:alnum:]+-[:alnum:]+-[:digit:]+$")) & TCGA_LUT_maf$cancer == x, "SAMPID"]})
names(wtRB_SAMPIDs_bycancer) <- unique(TCGA_LUT$cancer)

cancers_with_rbs <- names(Rbdeleterious_SAMPIDs_bycancer)[sapply(Rbdeleterious_SAMPIDs_bycancer, length) > 9]

#### try with linear model correction for confounders ####

if(!(exists("TCGA_MOR_list"))){
  TCGA_MOR_list <- readRDS("output/TCGA_MOR_list.rds")
}

Rb_mut_wt_HMTranks <- lapply(cancers_with_rbs, function(thiscancer){

  thisdata <- TCGA_MOR_list[[thiscancer]]
  HMTsumsforthisone <- colSums(thisdata[SET_HMTs$ensembl_gene_id, ])
  HMTsumsforthisone <- HMTsumsforthisone[names(HMTsumsforthisone) %in% TCGA_LUT
                                         $SAMPID]
  
  HMTsums_for_lm <- data.frame(HMTsums = HMTsumsforthisone)
  HMTsums_for_lm <- cbind(HMTsums_for_lm, TCGA_LUT[match(row.names(HMTsums_for_lm), TCGA_LUT$SAMPID), c("race", "gender", "tumour_stage", "days_to_birth", "sequencing_centre")])
  
  HMTsums_for_lm$days_to_birth <- as.numeric(HMTsums_for_lm$days_to_birth)
  
  HMTsums_for_lm <- HMTsums_for_lm[!apply(HMTsums_for_lm, 1, function(x){any(is.na(x))}), ]
  
  if(length(unique(HMTsums_for_lm$race)) == 1){
    
    racevariable <- NULL
    
  } else {
    
    racevariable <- "race"
    
  }
  
  if(length(unique(HMTsums_for_lm$tumour_stage)) == 1){
    
    tumourstagevariable <- NULL
    
  } else {
    
    tumourstagevariable <- "tumour_stage"
    
  }
  
  if(length(unique(HMTsums_for_lm$gender)) == 1) {
    
    gendervariable <- NULL
    
  } else {
    
    gendervariable <- "gender"
  }
  
  if(length(unique(HMTsums_for_lm$sequencing_centre)) == 1) {
    
    seqcentrevariable <- NULL
    
  } else {
    
    seqcentrevariable <- "sequencing_centre"
    
  }
  
  modelvariables <- c("days_to_birth", gendervariable, seqcentrevariable, tumourstagevariable, racevariable)
  
  f <- as.formula(
    paste("log10(HMTsums + 1)",
          paste(modelvariables, collapse = " + "),
          sep = " ~ "))
  
  HMTsums_for_lm[, "HMTresid"] <- lm(formula = f, data = HMTsums_for_lm)$residuals
  HMTsums_for_lm[, "HMTresid_ranks"] <- rank(HMTsums_for_lm$HMTresid) / (nrow(HMTsums_for_lm) + 1)
  
  output_list <- list()
  
  output_list[["mutant"]] <- HMTsums_for_lm[str_remove(row.names(HMTsums_for_lm), pattern = "-[:alnum:]+-[:alnum:]+-[:digit:]+$") %in%  str_remove(Rbdeleterious_SAMPIDs_bycancer[[thiscancer]], pattern = "-[:alnum:]+-[:alnum:]+-[:digit:]+$"), "HMTresid_ranks"] 
  names(output_list[["mutant"]]) <- row.names(HMTsums_for_lm[str_remove(row.names(HMTsums_for_lm), pattern = "-[:alnum:]+-[:alnum:]+-[:digit:]+$") %in%  str_remove(Rbdeleterious_SAMPIDs_bycancer[[thiscancer]], pattern = "-[:alnum:]+-[:alnum:]+-[:digit:]+$"), ])
  
  output_list[["wildtype"]] <- HMTsums_for_lm[str_remove(row.names(HMTsums_for_lm), pattern = "-[:alnum:]+-[:alnum:]+-[:digit:]+$") %in%  str_remove(wtRB_SAMPIDs_bycancer[[thiscancer]], pattern = "-[:alnum:]+-[:alnum:]+-[:digit:]+$"), "HMTresid_ranks"] 
  names(output_list[["wildtype"]]) <- row.names(HMTsums_for_lm[str_remove(row.names(HMTsums_for_lm), pattern = "-[:alnum:]+-[:alnum:]+-[:digit:]+$") %in%  str_remove(wtRB_SAMPIDs_bycancer[[thiscancer]], pattern = "-[:alnum:]+-[:alnum:]+-[:digit:]+$"), ])
  
  return(output_list)
  
})

names(Rb_mut_wt_HMTranks) <- cancers_with_rbs

Rb_mutant_wt_NNMTranks <- lapply(cancers_with_rbs, function(thiscancer){
  
  thisdata <- TCGA_MOR_list[[thiscancer]]
  
  NNMTforthisone <- thisdata[NNMT_ensembl, ]
  NNMTforthisone <- NNMTforthisone[names(NNMTforthisone) %in% TCGA_LUT$SAMPID]
  
  NNMT_for_lm <- data.frame(NNMT = NNMTforthisone)
  NNMT_for_lm <- cbind(NNMT_for_lm, TCGA_LUT[match(row.names(NNMT_for_lm), TCGA_LUT$SAMPID), c("race", "gender", "tumour_stage", "days_to_birth", "sequencing_centre")])
  
  NNMT_for_lm$days_to_birth <- as.numeric(NNMT_for_lm$days_to_birth)
  
  NNMT_for_lm <- NNMT_for_lm[!apply(NNMT_for_lm, 1, function(x){any(is.na(x))}), ]
  
  if(length(unique(NNMT_for_lm$race)) == 1){
    
    racevariable <- NULL
    
  } else {
    
    racevariable <- "race"
    
  }
  
  if(length(unique(NNMT_for_lm$tumour_stage)) == 1){
    
    tumourstagevariable <- NULL
    
  } else {
    
    tumourstagevariable <- "tumour_stage"
    
  }
  
  if(length(unique(NNMT_for_lm$gender)) == 1) {
    
    gendervariable <- NULL
    
  } else {
    
    gendervariable <- "gender"
  }
  
  if(length(unique(NNMT_for_lm$sequencing_centre)) == 1) {
    
    seqcentrevariable <- NULL
    
  } else {
    
    seqcentrevariable <- "sequencing_centre"
    
  }
  
  modelvariables <- c("days_to_birth", gendervariable, seqcentrevariable, tumourstagevariable, racevariable)
  
  f <- as.formula(
    paste("log10(NNMT + 1)",
          paste(modelvariables, collapse = " + "),
          sep = " ~ "))
  
  NNMT_for_lm[, "NNMTresid"] <- lm(formula = f, data = NNMT_for_lm)$residuals
  NNMT_for_lm[, "NNMTresid_ranks"] <- rank(NNMT_for_lm$NNMTresid) / (nrow(NNMT_for_lm) + 1)
  
  output_list <- list()
  
  output_list[["mutant"]] <- NNMT_for_lm[str_remove(row.names(NNMT_for_lm), pattern = "-[:alnum:]+-[:alnum:]+-[:digit:]+$") %in%  str_remove(Rbdeleterious_SAMPIDs_bycancer[[thiscancer]], pattern = "-[:alnum:]+-[:alnum:]+-[:digit:]+$"), "NNMTresid_ranks"] 
  names(output_list[["mutant"]]) <- row.names(NNMT_for_lm[str_remove(row.names(NNMT_for_lm), pattern = "-[:alnum:]+-[:alnum:]+-[:digit:]+$") %in%  str_remove(Rbdeleterious_SAMPIDs_bycancer[[thiscancer]], pattern = "-[:alnum:]+-[:alnum:]+-[:digit:]+$"), ])
  
  output_list[["wildtype"]] <- NNMT_for_lm[str_remove(row.names(NNMT_for_lm), pattern = "-[:alnum:]+-[:alnum:]+-[:digit:]+$") %in%  str_remove(wtRB_SAMPIDs_bycancer[[thiscancer]], pattern = "-[:alnum:]+-[:alnum:]+-[:digit:]+$"), "NNMTresid_ranks"] 
  names(output_list[["wildtype"]]) <- row.names(NNMT_for_lm[str_remove(row.names(NNMT_for_lm), pattern = "-[:alnum:]+-[:alnum:]+-[:digit:]+$") %in%  str_remove(wtRB_SAMPIDs_bycancer[[thiscancer]], pattern = "-[:alnum:]+-[:alnum:]+-[:digit:]+$"), ])
  
  return(output_list)
  
})

names(Rb_mutant_wt_NNMTranks) <- cancers_with_rbs

# do pan-cancer analysis 
# not all of those with >10 samples in 
RBmutant_select_cancers <- names(sapply(lapply(Rb_mut_wt_HMTranks, function(x){x[["mutant"]]}), length)[sapply(lapply(Rb_mut_wt_HMTranks, function(x){x[["mutant"]]}), length) > 9])

# in no cancer type is the direct comparison significant, although the mean is higher in 13/13 and the median in 11/13
signaltest <- lapply(Rb_mutant_wt_NNMTranks, function(x){
  
  wilcox.test(x[["mutant"]], x[["wildtype"]])
  
})

sapply(Rb_mut_wt_HMTranks, function(x){
  
  outvec <- c(mut_median = mean(x[["mutant"]]), 
              wt_median = mean(x[["wildtype"]]))
  
  return(outvec)
  
})[, RBmutant_select_cancers]

# boxplot(sapply(Rb_mut_wt_HMTranks, function(x){mean(x[["mutant"]])}), sapply(Rb_mut_wt_HMTranks, function(x){mean(x[["wildtype"]])}))

RBmutant_pancancer_HMTmedians <- data.frame(t(sapply(1:1000, function(i){

  Rb_select_vec <- unlist(lapply(RBmutant_select_cancers, function(X){
    
    sample(Rb_mut_wt_HMTranks[[X]][["mutant"]], 10)
    
  }))
  
  WT_select_vec <- unlist(lapply(RBmutant_select_cancers, function(X){
    
    sample(Rb_mut_wt_HMTranks[[X]][["wildtype"]], 10)
    
  }))

  output_vec <- c()
  
  output_vec["mut_rank_median"] <- median(Rb_select_vec)
  output_vec["WT_rank_median"] <- median(WT_select_vec)
  
  output_vec["p_value"] <-   wilcox.test(Rb_select_vec, WT_select_vec)$p.value
  
  return(output_vec)
  
})))

colnames(RBmutant_pancancer_HMTmedians)[1:2] <- c("Deleterious", "Wildtype")

# t.test p value is around 0
t.test(RBmutant_pancancer_HMTmedians$Deleterious, RBmutant_pancancer_HMTmedians$Wildtype)$p.value

RBmutant_pancancer_HMTmedians_melt_forplot <- reshape2::melt(RBmutant_pancancer_HMTmedians[, c("Wildtype", "Deleterious")])

write.table(RBmutant_pancancer_HMTmedians_melt_forplot,
            file = "plot_data/Fig 4/Fig_4D_HMTs_plotdata.txt",
            row.names = FALSE,
            col.names = TRUE,
            sep = "\t")

pdf("graphics/HMTmedianranks_RB1mut.pdf",
    height = 2.5,
    width = 2.5)

ggplot(RBmutant_pancancer_HMTmedians_melt_forplot, aes(x = variable, y = value)) + 
  geom_boxplot(fill = c("gray", "red")) +
  theme_classic() + 
  ylab("Pan-cancer sample\nmedian total HMT percentile") + 
  xlab(italic(RB1)~"mutation status") + 
  theme(axis.text.y = element_text(colour = "black",
                                   size = 10),
        axis.text.x = element_text(colour = "black",
                                   size = 8)) + 
  coord_cartesian(ylim = c(0.3, 0.8)) + 
  geom_hline(yintercept = 0.5,
             linetype = "dashed",
             col = "gray") + 
  geom_text(aes(x = 1.5, 
           y = 0.8, 
           label = "p%~~% 0"),
           parse = TRUE)

dev.off()

## now for NNMT

RBmutant_pancancer_NNMTmedians <- data.frame(t(sapply(1:1000, function(i){
  
  Rb_select_vec <- unlist(lapply(RBmutant_select_cancers, function(X){
    
    sample(Rb_mutant_wt_NNMTranks[[X]][["mutant"]], 10)
    
  }))
  
  WT_select_vec <- unlist(lapply(RBmutant_select_cancers, function(X){
    
    sample(Rb_mutant_wt_NNMTranks[[X]][["wildtype"]], 10)
    
  }))
  
  output_vec <- c()
  
  output_vec["mut_rank_median"] <- median(Rb_select_vec)
  output_vec["WT_rank_median"] <- median(WT_select_vec)
  
  output_vec["p_value"] <-   wilcox.test(Rb_select_vec, WT_select_vec)$p.value
  
  return(output_vec)
  
})))

colnames(RBmutant_pancancer_NNMTmedians)[1:2] <- c("Deleterious", "Wildtype")

# t.test p value is around 0
t.test(RBmutant_pancancer_NNMTmedians$Deleterious, RBmutant_pancancer_NNMTmedians$Wildtype)$p.value

RBmutant_pancancer_NNMTmedians_melt_forplot <- reshape2::melt(RBmutant_pancancer_NNMTmedians[, c("Wildtype", "Deleterious")])

write.table(RBmutant_pancancer_NNMTmedians_melt_forplot,
            file = "plot_data/Fig 4/Fig_4D_NNMT_plotdata.txt",
            row.names = FALSE,
            col.names = TRUE,
            sep = "\t")

pdf("graphics/NNMTmedianranks_RB1mut.pdf",
    height = 2.5,
    width = 2.5)

ggplot(RBmutant_pancancer_NNMTmedians_melt_forplot, aes(x = variable, y = value)) + 
  geom_boxplot(fill = c("gray", "red")) +
  theme_classic() + 
  ylab("Pan-cancer sample\nmedian total NNMT percentile") + 
  xlab(italic(RB1)~"mutation status") + 
  theme(axis.text.y = element_text(colour = "black",
                                   size = 10),
        axis.text.x = element_text(colour = "black",
                                   size = 8)) + 
  coord_cartesian(ylim = c(0.3, 0.8)) + 
  geom_hline(yintercept = 0.5,
             linetype = "dashed",
             col = "gray") + 
  annotate(geom = "text",
           x = 1.5,
           y = 0.8,
           label = substitute("p = 9.76 x"~10^-126))

dev.off()

#### is NNMT downstream of HMTs? Fig 4H ####

all_Rb_ranks <- do.call(c, lapply(Rb_mut_wt_HMTranks, function(X){X[["mutant"]]}))
names(all_Rb_ranks) <- str_remove(names(all_Rb_ranks), pattern = "^TCGA-[A-Z]{2,4}\\.")

all_wt_ranks <- do.call(c, lapply(Rb_mut_wt_HMTranks, function(X){X[["wildtype"]]}))
names(all_wt_ranks) <- str_remove(names(all_wt_ranks), pattern = "^TCGA-[A-Z]{2,4}\\.")

all_Rb_NNMT_ranks <- do.call(c, lapply(Rb_mutant_wt_NNMTranks, function(X){X[["mutant"]]}))
names(all_Rb_NNMT_ranks) <- str_remove(names(all_Rb_NNMT_ranks), pattern = "^TCGA-[A-Z]{2,4}\\.")

all_wt_NNMT_ranks <- do.call(c, lapply(Rb_mutant_wt_NNMTranks, function(X){X[["wildtype"]]}))
names(all_wt_NNMT_ranks) <- str_remove(names(all_wt_NNMT_ranks), pattern = "^TCGA-[A-Z]{2,4}\\.")

HMT_NNMT_rb_tvalues <- lapply(1:1000, function(i){
  
  tryCatch({
    
    Rb_NNMTselect_vec <- unlist(lapply(RBmutant_select_cancers, function(X){
      
      sample(Rb_mutant_wt_NNMTranks[[X]][["mutant"]], 10)
      
    }))
    
    WT_NNMTselect_vec <- unlist(lapply(RBmutant_select_cancers, function(X){
      
      sample(Rb_mutant_wt_NNMTranks[[X]][["wildtype"]], 10)
      
    }))
    
    Rb_HMTselect_vec <- all_Rb_ranks[match(names(Rb_NNMTselect_vec), names(all_Rb_ranks))]
    
    WT_HMTselect_vec <- all_wt_ranks[match(names(WT_NNMTselect_vec), names(all_wt_ranks))]
    
    Rb_lm_df <- data.frame(NNMT = c(WT_NNMTselect_vec, Rb_NNMTselect_vec),
                           HMT = c(WT_HMTselect_vec, Rb_HMTselect_vec))
    
    Rb_lm_df[row.names(Rb_lm_df) %in% names(WT_NNMTselect_vec), "Rbstatus"] <- "wildtype"
    Rb_lm_df[row.names(Rb_lm_df) %in% names(Rb_NNMTselect_vec), "Rbstatus"] <- "deleterious"
    
    NNMTfullmodel <- summary(lm(NNMT ~ HMT + Rbstatus, data = Rb_lm_df))
    NNMTpartmodel <- summary(lm(NNMT ~ Rbstatus, data = Rb_lm_df))
    
    HMTfullmodel <- summary(lm(HMT ~ NNMT + Rbstatus, data = Rb_lm_df))
    HMTpartmodel <- summary(lm(HMT ~ Rbstatus, data = Rb_lm_df))
    
    output_vec <- c(NNMTfullmodel$coefficients["HMT", "t value"],
                    NNMTfullmodel$coefficients["Rbstatuswildtype", "t value"],
                    HMTfullmodel$coefficients["Rbstatuswildtype", "t value"],
                    NNMTpartmodel$coefficients["Rbstatuswildtype", "t value"],
                    HMTpartmodel$coefficients["Rbstatuswildtype", "t value"])
    
    names(output_vec) <- c("HMT_NNMT_t",
                           "NNMT_fullRb_t",
                           "HMT_fullRb_t",
                           "NNMT_partRb_t",
                           "HMT_partRb_t")
    
    return(output_vec)
    
  }, error = function(e){
    
    message("These samples produced a problem")
    
    print(row.names(Rb_lm_df))
    
  })
  
  
})

HMT_NNMT_rb_tvalues <- HMT_NNMT_rb_tvalues[sapply(HMT_NNMT_rb_tvalues, function(x){length(x) == 5})]
HMT_NNMT_rb_tvalues_df <- data.frame(do.call(rbind, HMT_NNMT_rb_tvalues))
HMT_NNMT_rb_tvalues_df <- HMT_NNMT_rb_tvalues_df[, 2:5]

HMT_NNMT_rb_tvalues_meltdf <- reshape2::melt(HMT_NNMT_rb_tvalues_df)

HMT_NNMT_rb_tvalues_meltdf$variable <- factor(HMT_NNMT_rb_tvalues_meltdf$variable, levels = c("NNMT_partRb_t", "NNMT_fullRb_t", "HMT_partRb_t", "HMT_fullRb_t"))

t.test(as.numeric(sapply(HMT_NNMT_rb_tvalues, function(x){x["NNMT_partRb_t"]})),
       as.numeric(sapply(HMT_NNMT_rb_tvalues, function(x){x["NNMT_fullRb_t"]})))$p.value

t.test(as.numeric(sapply(HMT_NNMT_rb_tvalues, function(x){x["HMT_fullRb_t"]})),
       as.numeric(sapply(HMT_NNMT_rb_tvalues, function(x){x["HMT_partRb_t"]})))


ticklabs <- c("RB1 mutation", "RB1 mutation + HMTs", "RB1 mutation", "RB1 mutation + NNMT")

write.table(HMT_NNMT_rb_tvalues_meltdf,
            file = "plot_data/Fig 4/Fig_4H_plotdata.txt",
            row.names = FALSE,
            col.names = TRUE,
            sep = "\t")

pdf("graphics/TCGA_RBmut_HMTNNMTlinearmodel_tvals.pdf",
        width = 5,
        height = 3)

ggplot(HMT_NNMT_rb_tvalues_meltdf, aes(x = variable, y = -value)) + 
  geom_boxplot(fill = c("yellow",
                        "magenta",
                        "yellow",
                        "magenta")) + 
  theme_classic() + 
  geom_hline(yintercept = 0,
             linetype = "dashed",
             color = "gray") + 
  geom_vline(xintercept = 2.5,
             color = "black") + 
  ylab(italic(RB1)~"mutation linear model t-value") + 
  xlab("Linear model terms")  +
  coord_cartesian(ylim = c(-5, 5)) +
  scale_x_discrete(labels= ticklabs) + 
  theme(axis.text.x = element_text(colour = "black",
                                   size = 8),
        axis.text.y = element_text(colour = "black",
                                   size = 10)) + 
  # annotate(x = 1.5,
  #          y = 5,
  #          label = italic(NNMT)~"expression",
  #          geom = "text",
  #          colour = "red",
  #          fontface = 2) + 
  # annotate(x = 3.5,
  #          y = 5,
  #          label = "Total HMTs expression",
  #          geom = "text",
  #          colour = "red",
  #          fontface = 2) + 
  annotate(x = 1.5,
           y = 4.85,
           label = substitute("p = 3.84 x"~10^-67),
           geom = "text",
           colour = "black") +
  annotate(x = 3.5,
           y = 4.85,
           label = "p = 0.526",
           geom = "text",
           colour = "black") 
  # geom_segment(aes(x = 1, xend = 2, y = 4.25, yend = 4.25)) +
  # geom_segment(aes(x = 3, xend = 4, y = 4.25, yend = 4.25))

dev.off()

#### individual HMTs TCGA Supplementary File 7 / File S7 ####

# idenitifes any HMTs which are not expressed in some samples
lowHMTs <- sapply(SET_HMTs$ensembl_gene_id, function(thisone){
  
  any(sapply(TCGA_MOR_list, function(thisdf){
    
    any(thisdf[thisone, ] == 0)
    
  }) == TRUE)
  
})

SET_HMTs[match(names(lowHMTs)[lowHMTs], SET_HMTs$ensembl_gene_id), "hgnc_symbol"]

exprsingleHMTs_RBranks <-  lapply(SET_HMTs[!lowHMTs, "ensembl_gene_id"], function(thisHMT){

  Rb_mutant_ranks_temp <- lapply(RBmutant_select_cancers, function(thiscancer){
    
    thisdata <- TCGA_MOR_list[[thiscancer]]
    HMTexpforthisone <- thisdata[thisHMT, ]
    HMTexpforthisone <- HMTexpforthisone[names(HMTexpforthisone) %in% TCGA_LUT
                                         $SAMPID]
    
    HMTexp_for_lm <- data.frame(HMTexp = HMTexpforthisone)
    HMTexp_for_lm <- cbind(HMTexp_for_lm, TCGA_LUT[match(row.names(HMTexp_for_lm), TCGA_LUT$SAMPID), c("race", "gender", "tumour_stage", "days_to_birth", "sequencing_centre")])
    
    HMTexp_for_lm$days_to_birth <- as.numeric(HMTexp_for_lm$days_to_birth)
    
    HMTexp_for_lm <- HMTexp_for_lm[!apply(HMTexp_for_lm, 1, function(x){any(is.na(x))}), ]
    
    if(length(unique(HMTexp_for_lm$race)) == 1){
      
      racevariable <- NULL
      
    } else {
      
      racevariable <- "race"
      
    }
    
    if(length(unique(HMTexp_for_lm$tumour_stage)) == 1){
      
      tumourstagevariable <- NULL
      
    } else {
      
      tumourstagevariable <- "tumour_stage"
      
    }
    
    if(length(unique(HMTexp_for_lm$gender)) == 1) {
      
      gendervariable <- NULL
      
    } else {
      
      gendervariable <- "gender"
    }
    
    if(length(unique(HMTexp_for_lm$sequencing_centre)) == 1) {
      
      seqcentrevariable <- NULL
      
    } else {
      
      seqcentrevariable <- "sequencing_centre"
      
    }
    
    modelvariables <- c("days_to_birth", gendervariable, seqcentrevariable, tumourstagevariable, racevariable)
    
    f <- as.formula(
      paste("log10(HMTexp)",
            paste(modelvariables, collapse = " + "),
            sep = " ~ "))
    
    HMTexp_for_lm[, "HMTresid"] <- lm(formula = f, data = HMTexp_for_lm)$residuals
    HMTexp_for_lm[, "HMTresid_ranks"] <- rank(HMTexp_for_lm$HMTresid) / (nrow(HMTexp_for_lm) + 1)
    
    HMTexp_for_lm[str_remove(row.names(HMTexp_for_lm), pattern = "-[:alnum:]+-[:alnum:]+-[:digit:]+$") %in%  str_remove(Rbdeleterious_SAMPIDs_bycancer[[thiscancer]], pattern = "-[:alnum:]+-[:alnum:]+-[:digit:]+$"), "HMTresid_ranks"]
    
    output_vec <- HMTexp_for_lm[str_remove(row.names(HMTexp_for_lm), pattern = "-[:alnum:]+-[:alnum:]+-[:digit:]+$") %in%  str_remove(Rbdeleterious_SAMPIDs_bycancer[[thiscancer]], pattern = "-[:alnum:]+-[:alnum:]+-[:digit:]+$"), "HMTresid_ranks"] 
    names(output_vec) <- row.names(HMTexp_for_lm[str_remove(row.names(HMTexp_for_lm), pattern = "-[:alnum:]+-[:alnum:]+-[:digit:]+$") %in%  str_remove(Rbdeleterious_SAMPIDs_bycancer[[thiscancer]], pattern = "-[:alnum:]+-[:alnum:]+-[:digit:]+$"), ])
    
    return(output_vec)
    
  })
  
  names(Rb_mutant_ranks_temp) <- RBmutant_select_cancers
  
  Rb_wt_ranks_temp <- lapply(RBmutant_select_cancers, function(thiscancer){
    
    thisdata <- TCGA_MOR_list[[thiscancer]]
    HMTexpforthisone <- thisdata[thisHMT, ]
    HMTexpforthisone <- HMTexpforthisone[names(HMTexpforthisone) %in% TCGA_LUT
                                         $SAMPID]
    
    HMTexp_for_lm <- data.frame(HMTexp = HMTexpforthisone)
    HMTexp_for_lm <- cbind(HMTexp_for_lm, TCGA_LUT[match(row.names(HMTexp_for_lm), TCGA_LUT$SAMPID), c("race", "gender", "tumour_stage", "days_to_birth", "sequencing_centre")])
    
    HMTexp_for_lm$days_to_birth <- as.numeric(HMTexp_for_lm$days_to_birth)
    
    HMTexp_for_lm <- HMTexp_for_lm[!apply(HMTexp_for_lm, 1, function(x){any(is.na(x))}), ]
    
    if(length(unique(HMTexp_for_lm$race)) == 1){
      
      racevariable <- NULL
      
    } else {
      
      racevariable <- "race"
      
    }
    
    if(length(unique(HMTexp_for_lm$tumour_stage)) == 1){
      
      tumourstagevariable <- NULL
      
    } else {
      
      tumourstagevariable <- "tumour_stage"
      
    }
    
    if(length(unique(HMTexp_for_lm$gender)) == 1) {
      
      gendervariable <- NULL
      
    } else {
      
      gendervariable <- "gender"
    }
    
    if(length(unique(HMTexp_for_lm$sequencing_centre)) == 1) {
      
      seqcentrevariable <- NULL
      
    } else {
      
      seqcentrevariable <- "sequencing_centre"
      
    }
    
    modelvariables <- c("days_to_birth", gendervariable, seqcentrevariable, tumourstagevariable, racevariable)
    
    f <- as.formula(
      paste("HMTexp",
            paste(modelvariables, collapse = " + "),
            sep = " ~ "))
    
    HMTexp_for_lm[, "HMTresid"] <- lm(formula = f, data = HMTexp_for_lm)$residuals
    HMTexp_for_lm[, "HMTresid_ranks"] <- rank(HMTexp_for_lm$HMTresid) / (nrow(HMTexp_for_lm) + 1)
    
    output_vec <- HMTexp_for_lm[str_remove(row.names(HMTexp_for_lm), pattern = "-[:alnum:]+-[:alnum:]+-[:digit:]+$") %in%  str_remove(wtRB_SAMPIDs_bycancer[[thiscancer]], pattern = "-[:alnum:]+-[:alnum:]+-[:digit:]+$"), "HMTresid_ranks"] 
    names(output_vec) <- row.names(HMTexp_for_lm[str_remove(row.names(HMTexp_for_lm), pattern = "-[:alnum:]+-[:alnum:]+-[:digit:]+$") %in%  str_remove(wtRB_SAMPIDs_bycancer[[thiscancer]], pattern = "-[:alnum:]+-[:alnum:]+-[:digit:]+$"), ])
    
    return(output_vec)
    
  })
  
  names(Rb_wt_ranks_temp) <- RBmutant_select_cancers
  
  select_pvals <- sapply(1:1000, function(i){
    
    Rb_select_vec <- unlist(lapply(RBmutant_select_cancers, function(X){

      sample(Rb_mutant_ranks_temp[[X]], 10)
      
    }))
    
    
    WT_select_vec <- unlist(lapply(RBmutant_select_cancers, function(X){
      
      sample(Rb_wt_ranks_temp[[X]], 10)
      
    }))
    
    output_vec <- c()
    
    output_vec["mut_rank_median"] <- median(Rb_select_vec)
    output_vec["WT_rank_median"] <- median(WT_select_vec)
    
    output_vec["p_value"] <-   wilcox.test(Rb_select_vec, WT_select_vec)$p.value
    
    return(output_vec)
    
  })
  
  select_pvals_t <- data.frame(t(select_pvals))
  colnames(select_pvals_t)[1:2] <- c("Deleterious", "Wildtype")
  
  medians_melt_forplot <- reshape2::melt(select_pvals_t[, c("Wildtype", "Deleterious")])
  
  return(medians_melt_forplot)
  
})

names(exprsingleHMTs_RBranks) <- SET_HMTs[!lowHMTs, "hgnc_symbol"]

pdf("graphics/singleHMT_RBmutantplots.pdf")

# lapply(names(exprsingleHMTs_RBranks)[order(names(exprsingleHMTs_RBranks))], function(thisHGNC){
  for (i in 1:length(names(exprsingleHMTs_RBranks))){

  print(
    ggplot(exprsingleHMTs_RBranks[[names(exprsingleHMTs_RBranks)[order(names(exprsingleHMTs_RBranks))][i]]], aes(x = variable, y = value)) +
      geom_boxplot(fill = c("gray", "red")) +
      theme_classic() +
      ylab("Pan-cancer sample median total HMT percentile") +
      xlab(italic(RB1)~"mutation status") +
      ggtitle(paste0(names(exprsingleHMTs_RBranks)[order(names(exprsingleHMTs_RBranks))][i])) +
      theme(axis.text.y = element_text(colour = "black",
                                       size = 10),
            axis.text.x = element_text(colour = "black",
                                       size = 10)
      ) +
      geom_hline(yintercept = 0.5,
                 linetype = "dashed",
                 colour = "grey") +
      coord_cartesian(ylim = c(0.3, 0.8))
  )
  
}

dev.off()

#### Rb mutants in the CCLE for Fig 4E ####

mutcalls <- depmap::depmap_mutationCalls()

RBmutations <- mutcalls[mutcalls$gene_name == "RB1", ]
RBmutations_deleterious <- RBmutations[RBmutations$is_deleterious, ]

# what 
RBmutations_deleterious$depmap_id %in% unlist(sapply(CCLE_MOR_list, colnames))
RBmutations_deleterious_in_expression <- RBmutations_deleterious[RBmutations_deleterious$depmap_id %in% unlist(sapply(CCLE_MOR_list, colnames)), ]

# 90 unique lines with expression data and RB1 mutations
samples_w_muts <- unlist(sapply(CCLE_MOR_list, colnames))[unlist(sapply(CCLE_MOR_list, colnames)) %in% unique(mutcalls$depmap_id)]

samples_w_muts_noRB <- samples_w_muts[!samples_w_muts %in% unique(RBmutations$depmap_id)]
samples_w_muts_noRB_in_expression <- samples_w_muts_noRB[samples_w_muts_noRB %in% unlist(sapply(CCLE_MOR_list, colnames))]

# so 90 with RB deleterious, 1146 without.

if(!exists("metadata")){
  metadata <- depmap_metadata()
}

# View cancer types with mutations
# apart from Brain Cancer with 10, Lung Cancer with 43 is clearly dominant
table(metadata[metadata$depmap_id %in% RBmutations_deleterious_in_expression$depmap_id, "primary_disease"])
table(metadata[metadata$depmap_id %in% RBmutations_deleterious_in_expression$depmap_id, "subtype_disease"])

CCLE_Rb_HMT_wt_ranks <- lapply(names(CCLE_MOR_list), function(thiscancer){
  
  thisdata <- CCLE_MOR_list[[thiscancer]]
  
  HMTforthisone <- colSums(thisdata[SET_HMTs$ensembl_gene_id, ])
  HMTranks <- rank(HMTforthisone) / (length(HMTforthisone) + 1)
  
  output_vec <- HMTranks[names(HMTranks) %in% samples_w_muts_noRB]
  
  return(output_vec)
  
})

names(CCLE_Rb_HMT_wt_ranks) <- names(CCLE_MOR_list)

CCLE_Rb_HMT_mut_ranks <- lapply(names(CCLE_MOR_list), function(thiscancer){
  
  thisdata <- CCLE_MOR_list[[thiscancer]]
  
  HMTforthisone <- colSums(thisdata[SET_HMTs$ensembl_gene_id, ])
  
  HMTranks <- rank(HMTforthisone) / (length(HMTforthisone) + 1)
  
  output_vec <- HMTranks[names(HMTranks) %in% RBmutations_deleterious_in_expression$depmap_id]
  
  return(output_vec)
  
})

names(CCLE_Rb_HMT_mut_ranks) <- names(CCLE_MOR_list)

# note that most mutations are in SCLC and SCLC has higher HMT expression
# so must control for subtype
# build DF with HMTsums and subtypes (wildtype only)

WTraw_HMTsums <- colSums(CCLE_MOR_list[["Lung Cancer"]][, colnames(CCLE_MOR_list[["Lung Cancer"]]) %in% samples_w_muts_noRB][SET_HMTs$ensembl_gene_id, ])
MUTraw_HMTsums <- colSums(CCLE_MOR_list[["Lung Cancer"]][, colnames(CCLE_MOR_list[["Lung Cancer"]]) %in% RBmutations_deleterious_in_expression$depmap_id][SET_HMTs$ensembl_gene_id, ])

Lung_subtypes_df <- data.frame(depmap_id = names(WTraw_HMTsums),
                               HMTsums = WTraw_HMTsums,
                               subtype_disease = metadata[match(names(WTraw_HMTsums), metadata$depmap_id), "subtype_disease"],
                               genotype = "WT")

Lung_subtypes_mut_df <- data.frame(depmap_id = names(MUTraw_HMTsums),
                                   HMTsums = MUTraw_HMTsums,
                                   subtype_disease = metadata[match(names(MUTraw_HMTsums), metadata$depmap_id), "subtype_disease"],
                                   genotype = "MUT")

Lung_subtypes_df <- rbind(Lung_subtypes_df, Lung_subtypes_mut_df)
Lung_subtypes_df <- Lung_subtypes_df[!is.na(Lung_subtypes_df$subtype_disease), ]

Lung_SCLC_df <- Lung_subtypes_df[Lung_subtypes_df$subtype_disease == "Small Cell Lung Cancer (SCLC)", ]

Lung_SCLC_df[Lung_SCLC_df$genotype == "MUT", "genotype"] <- "Deleterious"
Lung_SCLC_df[Lung_SCLC_df$genotype == "WT", "genotype"] <- "Wildtype"

Lung_SCLC_df$genotype <- factor(Lung_SCLC_df$genotype, levels = c("Wildtype", "Deleterious"))

t.test(Lung_SCLC_df[Lung_SCLC_df$genotype == "Wildtype", "HMTsums"],
       Lung_SCLC_df[Lung_SCLC_df$genotype == "Deleterious", "HMTsums"])

scientific_10 <- function(x) {   parse(text=gsub("e\\+*", " %*% 10^", scales::scientific_format()(x))) }
options(scipen = 3)

write.table(Lung_SCLC_df,
            file = "plot_data/Fig 4/Fig_4E_plotdata.txt",
            row.names = TRUE,
            col.names = TRUE,
            sep = "\t")

pdf("graphics/CCLE_SCLC_RB1mut_HMTraw.pdf",
    height = 2.5,
    width = 2.5)

ggplot(Lung_SCLC_df, aes(x = genotype, y = HMTsums)) +
  geom_boxplot(fill = c("gray", "red")) +
  theme_classic() + 
  ylab("Total HMT expression\n(pseudocounts)") + 
  xlab(italic(RB1)~"mutation status") + 
  theme(axis.text.y = element_text(colour = "black",
                                   size = 10),
        axis.text.x = element_text(colour = "black",
                                   size = 8)) +
  coord_cartesian(ylim = c(125000, 250000)) +
  scale_y_continuous(label = scientific_10) +
  annotate(geom = "text",
           x = 1.5,
           y = 247500,
           label = "p = 0.041")

dev.off()

# add single HMTs

Lung_singleHMT_df <- CCLE_MOR_list[["Lung Cancer"]][SET_HMTs$ensembl_gene_id, match(Lung_subtypes_df$depmap_id, colnames(CCLE_MOR_list[["Lung Cancer"]]))]
row.names(Lung_singleHMT_df) <- SET_HMTs$hgnc_symbol

Lung_subtypes_df <- cbind(Lung_subtypes_df, t(Lung_singleHMT_df))

Lung_subtypes_df[, "NNMT"] <- CCLE_MOR_list[["Lung Cancer"]][NNMT_ensembl, match(Lung_subtypes_df$depmap_id, colnames(CCLE_MOR_list[["Lung Cancer"]]))]

# ANOVA

obj <- aov(HMTsums ~ subtype_disease + genotype, data = Lung_subtypes_df)

# significant. p = 0.0344
summary(obj)

# OK how about running the AOV for single HMTs? But it cant really be presented in such a nice way... could still do plots for the SCLC

singleHMTaov <- sapply(SET_HMTs$hgnc_symbol, function(thisHMT){
  
  unlist(summary(aov(as.formula(paste0(thisHMT, " ~ subtype_disease + genotype")), data = Lung_subtypes_df)))["Pr(>F)2"]
  
})

names(singleHMTaov) <- SET_HMTs$hgnc_symbol

# is the SS correct? Since our design is unbalanced we want to use the Type II SS
# since we have put genotype second in the anova model, this is the type II. all good.
adj_p_single <- p.adjust(singleHMTaov, method = "BH")
adj_p_single[adj_p_single < 0.3]

# confirm the direction of effect of signiifcant genes is as expected (upregulation)
singleHMTaov_lm <- sapply(SET_HMTs$hgnc_symbol, function(thisHMT){
  
  kmobj <- lm(as.formula(paste0(thisHMT, " ~ subtype_disease + genotype")), data = Lung_subtypes_df)$coefficients[["genotypeWT"]]
  
})

names(singleHMTaov_lm) <- SET_HMTs$hgnc_symbol

singleHMTaov_lm[names(adj_p_single[adj_p_single < 0.3])]

# all negative (for WT ie expected direction), except SMYD5 which is strongly down, as in TCGA. EZH2 strongest

#### DOROTHEA AND TR ACTIVITY for Fig 4F ####

TCGA_gene_lookup <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                          filters = "ensembl_gene_id",
                          values = row.names(TCGA_MOR_list[[1]]),
                          mart = ensembl)


GTEX_gene_lookup <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                          filters = "ensembl_gene_id",
                          values = row.names(GTEX_MOR_list[[1]]),
                          mart = ensembl)

TCGA_MOR_list_for_decouple <- lapply(TCGA_MOR_list, function(x){

  row.names(x) <- TCGA_gene_lookup[match(row.names(x), TCGA_gene_lookup$ensembl_gene_id), "hgnc_symbol"]

  x <- x[!is.na(row.names(x)), ]
  x <- x[!row.names(x) == "", ]

  return(x)

})

GTEX_MOR_list_for_decouple <- lapply(GTEX_MOR_list, function(x){
  
  row.names(x) <- GTEX_gene_lookup[match(row.names(x), GTEX_gene_lookup$ensembl_gene_id), "hgnc_symbol"]
  
  x <- x[!is.na(row.names(x)), ]
  x <- x[!row.names(x) == "", ]
  
  return(x)
  
})

GTEX_crosstissueMOR_list_for_decouple <- lapply(GTEX_MOR_list, function(x){
  
  thesesamples <- colnames(x)
  tempdata <- GTEX_MOR_across_tissues[, thesesamples]
  
  row.names(tempdata) <- GTEX_gene_lookup[match(row.names(tempdata), GTEX_gene_lookup$ensembl_gene_id), "hgnc_symbol"]
  
  tempdata <- tempdata[!is.na(row.names(tempdata)), ]
  tempdata <- tempdata[!row.names(tempdata) == "", ]
  
  return(tempdata)
  
})

dorothea_tcga <- dorothea_hs_pancancer
dorothea_tcga_highconfidence <- dorothea_tcga[dorothea_tcga$confidence %in% LETTERS[1:4], ]

# Here weight the scores by confidence as in the get_dorothea function of the decoupleR package
# LETTERS is built in vector of capital letters in alphabetical order
for (i in 1:4){

  dorothea_tcga_highconfidence[dorothea_tcga_highconfidence$confidence == LETTERS[i], "mor"] <- unlist(dorothea_tcga_highconfidence[dorothea_tcga_highconfidence$confidence == LETTERS[i], "mor"])/i

}

# exclude TFs with fewer than 10 target genes to judge by

dorothea_tcga_highconfidence <- dorothea_tcga_highconfidence[!(dorothea_tcga_highconfidence$tf %in% names(table(dorothea_tcga_highconfidence$tf)[table(dorothea_tcga_highconfidence$tf) < 10])), ]

saveRDS(dorothea_tcga_highconfidence, "output/dorothea_tcga_highconfidence.rds")
# dorothea_tcga_highconfidence <- readRDS("output/dorothea_tcga_highconfidence.rds")

dorothea_GTEX <- dorothea_hs
dorothea_GTEX_highconfidence <- dorothea_GTEX[dorothea_GTEX$confidence %in% LETTERS[1:4], ]

# Here weight the scores by confidence as in the get_dorothea function of the decoupleR package
# LETTERS is built in vector of capital letters in alphabetical order
for (i in 1:4){
  
  dorothea_GTEX_highconfidence[dorothea_GTEX_highconfidence$confidence == LETTERS[i], "mor"] <- unlist(dorothea_GTEX_highconfidence[dorothea_GTEX_highconfidence$confidence == LETTERS[i], "mor"])/i
  
}

# exclude TFs with fewer than 10 target genes to judge by

dorothea_GTEX_highconfidence <- dorothea_GTEX_highconfidence[!(dorothea_GTEX_highconfidence$tf %in% names(table(dorothea_GTEX_highconfidence$tf)[table(dorothea_GTEX_highconfidence$tf) < 10])), ]

saveRDS(dorothea_GTEX_highconfidence, "output/dorothea_GTEX_highconfidence.rds")
# dorothea_GTEX_highconfidence <- readRDS("output/dorothea_GTEX_highconfidence.rds")

# run through and estimate for TCGA samples

for(i in 1:length(TCGA_MOR_list_for_decouple)){

decouple_temp <- decouple(mat = TCGA_MOR_list_for_decouple[[i]],
                         network = dorothea_tcga_highconfidence,
                         .source = "tf",
                         .target = "target"
)

consensus_temp <- run_consensus(decouple_temp)

tabletemp <- sapply(unique(consensus_temp$condition), function(thissample){

  activities <- unlist(consensus_activities[consensus_activities$condition == thissample, "score"])
  names(activities) <- unlist(consensus_activities[consensus_activities$condition == thissample, "source"])

  return(activities[order(names(activities))])

})

saveRDS(tabletemp, paste0("output/", names(TCGA_MOR_list_for_decouple)[i], "_TFconsensus.rds"))

}

# run through and estimate for GTEX samples

for(i in 1:length(GTEX_MOR_list_for_decouple)){
  
  decouple_temp <- decouple(mat = GTEX_MOR_list_for_decouple[[i]],
                            network = dorothea_GTEX_highconfidence,
                            .source = "tf",
                            .target = "target"
  )
  
  consensus_temp <- run_consensus(decouple_temp)
  
  tabletemp <- sapply(unique(consensus_temp$condition), function(thissample){
    
    activities <- unlist(consensus_activities[consensus_activities$condition == thissample, "score"])
    names(activities) <- unlist(consensus_activities[consensus_activities$condition == thissample, "source"])
    
    return(activities[order(names(activities))])
    
  })
  
  saveRDS(tabletemp, paste0("output/", names(GTEX_MOR_list_for_decouple)[i], "_TFconsensus.rds"))
  
}

# run through and estimate for GTEX samples with cross tissue MOR normalisation

for(i in 1:length(GTEX_crosstissueMOR_list_for_decouple)){
  
  decouple_temp <- decouple(mat = GTEX_crosstissueMOR_list_for_decouple[[i]],
                            network = dorothea_GTEX_highconfidence,
                            .source = "tf",
                            .target = "target"
  )
  
  consensus_temp <- run_consensus(decouple_temp)
  
  tabletemp <- sapply(unique(consensus_temp$condition), function(thissample){
    
    activities <- unlist(consensus_activities[consensus_activities$condition == thissample, "score"])
    names(activities) <- unlist(consensus_activities[consensus_activities$condition == thissample, "source"])
    
    return(activities[order(names(activities))])
    
  })
  
  saveRDS(tabletemp, paste0("output/", names(GTEX_MOR_list_for_decouple)[i], "_TFcrostissueconsensus.rds"))
  
}

# load in estimations

TCGA_activities_files <- list.files("output") 
TCGA_activities_files <- TCGA_activities_files[str_detect(TCGA_activities_files, pattern = "^TCGA-[A-Z]+")]
TCGA_activities_files <- TCGA_activities_files[str_detect(TCGA_activities_files, pattern = "TFconsensus")]

TCGA_activities_list <- lapply(TCGA_activities_files, function(x){
  
  readRDS(paste0("output/", x))
  
})

names(TCGA_activities_list) <- str_remove(TCGA_activities_files, pattern = "_TFconsensus.rds")

GTEX_activities_files <- list.files("output") 
GTEX_activities_files <- GTEX_activities_files[str_detect(GTEX_activities_files, pattern = "TFconsensus")]
GTEX_activities_files <- GTEX_activities_files[!str_detect(GTEX_activities_files, pattern = "TCGA")]

GTEX_activities_list <- lapply(GTEX_activities_files, function(x){
  
  readRDS(paste0("output/", x))
  
})

names(GTEX_activities_list) <- str_remove(GTEX_activities_files, pattern = "_TFconsensus.rds")

GTEX_crosstissueactivities_files <- GTEX_activities_files[str_detect(GTEX_activities_files, pattern = "crostissue")]
GTEX_crosstissueactivities_files <- GTEX_crosstissueactivities_files[!str_detect(GTEX_crosstissueactivities_files, pattern = "TCGA")]

GTEX_crosstissueactivities_list <- lapply(GTEX_crosstissueactivities_files, function(x){
  
  readRDS(paste0("output/", x))
  
})

names(GTEX_crosstissueactivities_list) <- str_remove(GTEX_crosstissueactivities_files, pattern = "_TFcrostissueconsensus.rds")

# restrict TCGA to samples in LUT
TCGA_activities_list <- lapply(TCGA_activities_list, function(x){
  
  x[, colnames(x) %in% TCGA_LUT$SAMPID]
  
})

# 48 samples not in GTEX LUT. Exclude them too
# sum(unlist(!(lapply(GTEX_activities_list, colnames)) %in% GTEX_LUT$SAMPID))

GTEX_activities_list <- lapply(GTEX_activities_list, function(x){
  
  x[, colnames(x) %in% GTEX_LUT$SAMPID]
  
})

GTEX_crosstissueactivities_list <- lapply(GTEX_crosstissueactivities_list, function(x){
  
  x[, colnames(x) %in% GTEX_LUT$SAMPID]
  
})

# limit crosstissue set to tisues in regular set
GTEX_crosstissueactivities_list <- GTEX_crosstissueactivities_list[names(GTEX_activities_list)]

TCGA_MOR_HMT_and_NNMT <- lapply(TCGA_MOR_list, function(x){
  
  totalHMTs <- colSums(x[SET_HMTs$ensembl_gene_id, ])
  
  singlegenes <- x[c(SET_HMTs$ensembl_gene_id, NNMT_ensembl, PEMT_ensembl), ]
  
  output_df <- rbind(singlegenes, totalHMTs)
  
})

TCGA_MOR_HMT_and_NNMT_corrected <- lapply(TCGA_MOR_HMT_and_NNMT, function(x){TCGA.apply.lm.to.counts(log10(x + 1))})

GTEX_MOR_HMT_and_PEMT <- lapply(GTEX_MOR_list, function(x){
  
  totalHMTs <- colSums(x[SET_HMTs$ensembl_gene_id, ])
  
  singlegenes <- x[c(SET_HMTs$ensembl_gene_id, PEMT_ensembl), ]
  
  output_df <- rbind(singlegenes, totalHMTs)
  
})

GTEX_MOR_HMT_and_PEMT_corrected <- lapply(GTEX_MOR_HMT_and_PEMT, function(x){
  
  apply.lm.to.counts(log10(x + 1), lookuptable = GTEX_LUT)
  
})

# correct TR activities too
TCGA_TFactivities_corrected <- lapply(TCGA_activities_list, TCGA.apply.lm.to.counts)

GTEX_TFactivities_corrected <- lapply(GTEX_activities_list, function(x){
  
  apply.lm.to.counts(x, lookuptable = GTEX_LUT)
  
})

# for GTEX
# correlations
GTEX_HMTs_tissue_correlations <- sapply(1:length(GTEX_TFactivities_corrected), function(i){
  
  apply(GTEX_TFactivities_corrected[[i]], 1, function(x){
    
    cor(x, GTEX_MOR_HMT_and_PEMT_corrected[[i]]["totalHMTs", names(x)], method = "spearman")
    
  })
  
})

colnames(GTEX_HMTs_tissue_correlations) <- names(GTEX_TFactivities_corrected)

# FDRs
GTEX_HMTs_tissue_corr_FDRs <- sapply(1:length(GTEX_TFactivities_corrected), function(i){
  
  p.adjust(apply(GTEX_TFactivities_corrected[[i]], 1, function(x){
    
    cor.test(x, GTEX_MOR_HMT_and_PEMT_corrected[[i]]["totalHMTs", names(x)], method = "spearman")$p.value
    
  }), method = "BH")
  
})

colnames(GTEX_HMTs_tissue_corr_FDRs) <- names(GTEX_TFactivities_corrected)

# now for PEMT
GTEX_PEMT_tissue_correlations <- sapply(1:length(GTEX_TFactivities_corrected), function(i){
  
  apply(GTEX_TFactivities_corrected[[i]], 1, function(x){
    
    cor(x, GTEX_MOR_HMT_and_PEMT_corrected[[i]][PEMT_ensembl, names(x)], method = "spearman")
    
  })
  
})

colnames(GTEX_PEMT_tissue_correlations) <- names(GTEX_TFactivities_corrected)

# FDRs
GTEX_PEMT_tissue_corr_FDRs <- sapply(1:length(GTEX_TFactivities_corrected), function(i){
  
  p.adjust(apply(GTEX_TFactivities_corrected[[i]], 1, function(x){
    
    cor.test(x, GTEX_MOR_HMT_and_PEMT_corrected[[i]][PEMT_ensembl, names(x)], method = "spearman")$p.value
    
  }), method = "BH")
  
})

colnames(GTEX_PEMT_tissue_corr_FDRs) <- names(GTEX_TFactivities_corrected)

# now for TCGA

# correlations
TCGA_HMTs_tissue_correlations <- sapply(1:length(TCGA_TFactivities_corrected), function(i){

  apply(TCGA_TFactivities_corrected[[i]], 1, function(x){
    
    cor(x, TCGA_MOR_HMT_and_NNMT_corrected[[i]]["totalHMTs", names(x)], method = "spearman")
    
  })
  
})

colnames(TCGA_HMTs_tissue_correlations) <- names(TCGA_TFactivities_corrected)

# FDRs 
TCGA_HMTs_tissue_corr_FDRs <- sapply(1:length(TCGA_TFactivities_corrected), function(i){
  
  p.adjust(apply(TCGA_TFactivities_corrected[[i]], 1, function(x){
    
    cor.test(x, TCGA_MOR_HMT_and_NNMT_corrected[[i]]["totalHMTs", names(x)], method = "spearman")$p.value
    
  }), method = "BH")
  
})

colnames(TCGA_HMTs_tissue_corr_FDRs) <- names(TCGA_TFactivities_corrected)

# number of TCGA cancers w positive HMT - E2F1 relationship? 32 / 33
sum(TCGA_HMTs_tissue_correlations["E2F1", ]> 0)

# significant? 30/33
sum(TCGA_HMTs_tissue_corr_FDRs["E2F1", ] < 0.1)

# number of GTEX cancers w positive HMT - E2F1 relationship? 37 / 48
sum(GTEX_HMTs_tissue_correlations["E2F1", ]> 0)

# significant? 26/48
sum(GTEX_HMTs_tissue_corr_FDRs["E2F1", ] < 0.1 & GTEX_HMTs_tissue_correlations["E2F1", ] > 0)

# also NNMT 

TCGA_NNMT_tissue_correlations <- sapply(1:length(TCGA_TFactivities_corrected), function(i){
  
  apply(TCGA_TFactivities_corrected[[i]], 1, function(x){
    
    cor(x, TCGA_MOR_HMT_and_NNMT_corrected[[i]][NNMT_ensembl, names(x)], method = "spearman")
    
  })
  
})

colnames(TCGA_NNMT_tissue_correlations) <- names(TCGA_TFactivities_corrected)

# FDRs 
TCGA_NNMT_tissue_corr_FDRs <- sapply(1:length(TCGA_TFactivities_corrected), function(i){
  
  p.adjust(apply(TCGA_TFactivities_corrected[[i]], 1, function(x){
    
    cor.test(x, TCGA_MOR_HMT_and_NNMT_corrected[[i]][NNMT_ensembl, names(x)], method = "spearman")$p.value
    
  }), method = "BH")
  
})

colnames(TCGA_NNMT_tissue_corr_FDRs) <- names(TCGA_TFactivities_corrected)

TCGA_NNMT_tissue_correlations <- TCGA_NNMT_tissue_correlations[!row.names(TCGA_NNMT_tissue_correlations) %in% SET_HMTs$hgnc_symbol, ]
TCGA_NNMT_tissue_corr_FDRs <- TCGA_NNMT_tissue_corr_FDRs[!row.names(TCGA_NNMT_tissue_corr_FDRs) %in% SET_HMTs$hgnc_symbol, ]

TCGA_MOR_HMT_and_NNMT_corrected_rankpercent <- lapply(TCGA_MOR_HMT_and_NNMT_corrected, function(x){
  
  as.data.frame(t(apply(x, 1, rank)/(ncol(x) + 1)))
  
})

TCGA_TFactivities_corrected_rankpercent <- lapply(TCGA_TFactivities_corrected, function(x){
  
  as.data.frame(t(apply(x, 1, rank)/(ncol(x) + 1)))
  
})

TCGA_TFactivities_corrected_rankpercent_df <- do.call(cbind, TCGA_TFactivities_corrected_rankpercent)

iterations_out <- lapply(1:1000, function(i){

  # smallest number of samples of any any cancer type is 36 (CHOL)
  tempsample <- lapply(TCGA_MOR_HMT_and_NNMT_corrected_rankpercent, function(x){
    
    x[, sample(colnames(x), 36)]
    
  })
  
  tempsample_df <- do.call(cbind, tempsample)
  
  HMTout <- apply(TCGA_TFactivities_corrected_rankpercent_df, 1, function(y){
    
    cor.test(unlist(tempsample_df["totalHMTs", ]), y[colnames(tempsample_df)], method = "spearman")$estimate
    
  })
  
  NNMTout <- apply(TCGA_TFactivities_corrected_rankpercent_df, 1, function(y){
    
    cor.test(unlist(tempsample_df[NNMT_ensembl, ]), y[colnames(tempsample_df)], method = "spearman")$estimate
    
  })
  
  list(HMT = HMTout,
       NNMT = NNMTout)
  
})

TCGA_HMT_median_TFcorrs <- matrixStats::rowMedians(as.matrix(data.frame((lapply(iterations_out, function(x){x["HMT"]})))))
TCGA_NNMT_median_TFcorrs <- matrixStats::rowMedians(as.matrix(data.frame((lapply(iterations_out, function(x){x["NNMT"]})))))

names(TCGA_HMT_median_TFcorrs) <- names(iterations_out[[1]][["HMT"]])
names(TCGA_NNMT_median_TFcorrs) <- names(iterations_out[[1]][["HMT"]])

# E2F1 rank?
which(row.names(TCGA_HMT_median_TFcorrs_df) == "E2F1")

NNMTHMTiterations_out <- sapply(1:1000, function(i){
  
  # smallest number of samples of any any cancer type is 36 (CHOL)
  tempsample <- lapply(TCGA_MOR_HMT_and_NNMT_corrected_rankpercent, function(x){
    
    x[, sample(colnames(x), 36)]
    
  })
  
  tempsample_df <- do.call(cbind, tempsample)
  
  cor.test(unlist(tempsample_df["totalHMTs", ]), unlist(tempsample_df[NNMT_ensembl, ]), method = "spearman")$estimate
  
})

TCGA_median_NNMTHMTcorr <- median(NNMTHMTiterations_out)

TCGA_HMT_median_TFcorrs_df <- data.frame(TCGA_HMT_median_TFcorrs[order(TCGA_HMT_median_TFcorrs, decreasing = T)])
TCGA_HMT_median_TFcorrs_df <- data.frame(TCGA_HMT_median_TFcorrs_df[!row.names(TCGA_HMT_median_TFcorrs_df) %in% SET_HMTs$hgnc_symbol, ], 
                                         row.names = row.names(TCGA_HMT_median_TFcorrs_df)[!row.names(TCGA_HMT_median_TFcorrs_df) %in% SET_HMTs$hgnc_symbol])
colnames(TCGA_HMT_median_TFcorrs_df) <- "Spearman_rho"

TCGA_NNMT_median_TFcorrs_df <- data.frame(TCGA_NNMT_median_TFcorrs[order(TCGA_NNMT_median_TFcorrs, decreasing = T)])
TCGA_NNMT_median_TFcorrs_df <- data.frame(TCGA_NNMT_median_TFcorrs_df[!row.names(TCGA_NNMT_median_TFcorrs_df) %in% SET_HMTs$hgnc_symbol, ], 
                                          row.names = row.names(TCGA_NNMT_median_TFcorrs_df)[!row.names(TCGA_NNMT_median_TFcorrs_df) %in% SET_HMTs$hgnc_symbol])
colnames(TCGA_NNMT_median_TFcorrs_df) <- "Spearman_rho"

# Table S6

openxlsx::write.xlsx(list("HMT_correlations" = TCGA_HMTs_tissue_correlations,
                          "HMT_FDRs" = TCGA_HMTs_tissue_corr_FDRs,
                          "HMT_pancancer_correlations" = TCGA_HMT_median_TFcorrs_df,
                          "NNMT_correlations" = TCGA_NNMT_tissue_correlations,
                          "NNMT_FDRs" = TCGA_NNMT_tissue_corr_FDRs,
                          "NNMT_pancancer_correlations" = TCGA_NNMT_median_TFcorrs_df),
                     "output/TCGA_TRactivity_correlations.xlsx",
                     row.names = TRUE)

# GTEX pan

if(!exists("GTEX_MOR_across_tissues")){
GTEX_MOR_across_tissues <- readRDS("output/GTEX_MOR_across_tissues.rds")
}

GTEX_MOR_across_tissues <- data.frame(GTEX_MOR_across_tissues)

GTEX_MOR_across_tissues_HMTs_PEMT <- GTEX_MOR_across_tissues[SET_HMTs$ensembl_gene_id, ]
GTEX_MOR_across_tissues_HMTs_PEMT["totalHMTs", ] <- colSums(GTEX_MOR_across_tissues_HMTs_PEMT)
GTEX_MOR_across_tissues_HMTs_PEMT[PEMT_ensembl, ] <- GTEX_MOR_across_tissues[PEMT_ensembl, ]

rm(GTEX_MOR_across_tissues)

GTEX_MOR_across_tissues_HMTs_PEMT_resid <- GTEX.apply.lm.to.combined.counts(countdata = log10(GTEX_MOR_across_tissues_HMTs_PEMT + 1),
                                                                           lookuptable = GTEX_LUT)
saveRDS(GTEX_MOR_across_tissues_HMTs_PEMT_resid, "output/GTEX_MOR_across_tissues_HMTs_PEMT_resid.rds")
# GTEX_MOR_across_tissues_HMTs_PEMT_resid <- readRDS("output/GTEX_MOR_across_tissues_HMTs_PEMT_resid.rds")

GTEX_crosstissueactivities_df <- do.call(cbind, GTEX_crosstissueactivities_list)
GTEX_activities_across_tissues_resid <- GTEX.apply.lm.to.combined.counts(countdata = GTEX_crosstissueactivities_df,
                                                                        lookuptable = GTEX_LUT)

saveRDS(GTEX_activities_across_tissues_resid, "output/GTEX_activities_across_tissues_resid.rds")
GTEX_activities_across_tissues_resid <- readRDS("output/GTEX_activities_across_tissues_resid.rds")

# rank percentile transformation
tissue_types <- names(GTEX_MOR_HMT_and_PEMT)

GTEX_MOR_across_tissues_HMTs_PEMTresid_rankpercentile_list <- lapply(tissue_types, function(thistissue){
  
  t(apply(GTEX_MOR_across_tissues_HMTs_PEMT_resid[, colnames(GTEX_MOR_across_tissues_HMTs_PEMT_resid) %in% GTEX_LUT[GTEX_LUT$Tissue == thistissue, "SAMPID"]], 1, function(thisgene){
    
    rank(thisgene) / (length(thisgene) + 1)
    
  }))
  
})

names(GTEX_MOR_across_tissues_HMTs_PEMTresid_rankpercentile_list) <- tissue_types

GTEX_activities_across_tissues_rankpercentile_list <- lapply(tissue_types, function(thistissue){
  
  t(apply(GTEX_activities_across_tissues_resid[, colnames(GTEX_activities_across_tissues_resid) %in% GTEX_LUT[GTEX_LUT$Tissue == thistissue, "SAMPID"]], 1, function(thisgene){
    
    rank(thisgene) / (length(thisgene) + 1)
    
  }))
  
})

names(GTEX_activities_across_tissues_rankpercentile_list) <- tissue_types

GTEX_activities_across_tissues_rankpercentile_df <- do.call(cbind, GTEX_activities_across_tissues_rankpercentile_list)

# take samples and run

GTEXiterations_out <- sapply(1:100, function(i){
  
  # smallest number of samples of any any cancer type is 36 (CHOL)
  tempsample <- lapply(GTEX_MOR_across_tissues_HMTs_PEMTresid_rankpercentile_list, function(x){
    
    x[, sample(colnames(x), 100)]
    
  })
  
  tempsample_df <- do.call(cbind, tempsample)
  
  apply(GTEX_activities_across_tissues_rankpercentile_df, 1, function(y){
    
    cor.test(unlist(tempsample_df["totalHMTs", ]), y[colnames(tempsample_df)], method = "spearman")$estimate
    
  })
  
})

GTEX_median_TFcorrs <- rowMedians(GTEXiterations_out)
names(GTEX_median_TFcorrs) <- row.names(GTEXiterations_out)

GTEX_median_TFcorrs_df <- data.frame(GTEX_median_TFcorrs[order(GTEX_median_TFcorrs, decreasing = T)])
colnames(GTEX_median_TFcorrs_df) <- "Spearman_rho"

# take samples and run PEMT

GTEX_PEMTiterations_out <- sapply(1:100, function(i){
  
  tempsample <- lapply(GTEX_MOR_across_tissues_HMTs_PEMTresid_rankpercentile_list, function(x){
    
    x[, sample(colnames(x), 100)]
    
  })
  
  tempsample_df <- do.call(cbind, tempsample)
  
  apply(GTEX_activities_across_tissues_rankpercentile_df, 1, function(y){
    
    cor.test(unlist(tempsample_df[PEMT_ensembl, ]), y[colnames(tempsample_df)], method = "spearman")$estimate
    
  })
  
})

GTEX_median_PEMTTFcorrs <- rowMedians(GTEX_PEMTiterations_out)
names(GTEX_median_PEMTTFcorrs) <- row.names(GTEX_PEMTiterations_out)

GTEX_median_PEMTTFcorrs_df <- data.frame(GTEX_median_PEMTTFcorrs[order(GTEX_median_PEMTTFcorrs, decreasing = T)])
colnames(GTEX_median_PEMTTFcorrs_df) <- "Spearman_rho"

GTEX_median_PEMTTFcorrs_df[] <- GTEX_median_PEMTTFcorrs_df[!GTEX_median_PEMTTFcorrs_df[, 1] %in% SET_HMTs$hgnc_symbol, ]
GTEX_median_PEMTTFcorrs_df_names <- GTEX_median_PEMTTFcorrs_df[, 1]

# Table S7
GTEX_HMTs_tissue_correlations <- read.xlsx("output/GTEX_HMTs_TRactivity_correlations.xlsx",
                                           sheet = 1)
row.names(GTEX_HMTs_tissue_correlations) <- GTEX_HMTs_tissue_correlations[, 1]
GTEX_HMTs_tissue_correlations <- GTEX_HMTs_tissue_correlations[, 2:ncol(GTEX_HMTs_tissue_correlations)]

GTEX_HMTs_tissue_corr_FDRs <- read.xlsx("output/GTEX_HMTs_TRactivity_correlations.xlsx",
                                        sheet = 2)
row.names(GTEX_HMTs_tissue_corr_FDRs) <- GTEX_HMTs_tissue_corr_FDRs[, 1]
GTEX_HMTs_tissue_corr_FDRs <- GTEX_HMTs_tissue_corr_FDRs[, 2:ncol(GTEX_HMTs_tissue_corr_FDRs)]

GTEX_median_TFcorrs_df <- read.xlsx("output/GTEX_HMTs_TRactivity_correlations.xlsx",
                                    sheet = 3)
GTEX_median_TFcorrs_df <- GTEX_median_TFcorrs_df[!GTEX_median_TFcorrs_df[, 1] %in% SET_HMTs$hgnc_symbol, ]
GTEX_median_TFcorrs_df_names <- GTEX_median_TFcorrs_df[, 1]
GTEX_median_TFcorrs_df <- data.frame(GTEX_median_TFcorrs_df[, 2:ncol(GTEX_median_TFcorrs_df)])
row.names(GTEX_median_TFcorrs_df) <- GTEX_median_TFcorrs_df_names
colnames(GTEX_median_TFcorrs_df) <- "Spearman_rho"

openxlsx::write.xlsx(list("HMT_correlations" = GTEX_HMTs_tissue_correlations[!row.names(GTEX_HMTs_tissue_correlations) %in% SET_HMTs$hgnc_symbol, ],
                          "HMT_FDRs" = GTEX_HMTs_tissue_corr_FDRs[!row.names(GTEX_HMTs_tissue_corr_FDRs) %in% SET_HMTs$hgnc_symbol, ],
                          "HMT_pantissue_correlations" = GTEX_median_TFcorrs_df,
                          "PEMT_correlations" = GTEX_PEMT_tissue_correlations[!row.names(GTEX_PEMT_tissue_correlations) %in% SET_HMTs$hgnc_symbol, ],
                          "PEMT_FDRs" = GTEX_PEMT_tissue_corr_FDRs[!row.names(GTEX_PEMT_tissue_corr_FDRs) %in% SET_HMTs$hgnc_symbol, ],
                          "PEMT_pantissue_correlations" = GTEX_median_PEMTTFcorrs_df),
                     "output/GTEX_HMTsPEMT_TRactivity_correlations.xlsx",
                     row.names = TRUE)

#### breast cancer scatterplot Fig 4F ####

TCGA_BRCA_E2F1_plot_df <- data.frame(E2F1 = TCGA_TFactivities_corrected[["TCGA-BRCA"]]["E2F1", ], 
           HMTs = TCGA_MOR_HMT_and_NNMT_corrected[["TCGA-BRCA"]]["totalHMTs", colnames(TCGA_TFactivities_corrected[["TCGA-BRCA"]])], method = "spearman")

cor.test(TCGA_BRCA_E2F1_plot_df[, 1], TCGA_BRCA_E2F1_plot_df[, 2], method = "spearman")$p.value

write.table(TCGA_BRCA_E2F1_plot_df,
            "plot_data/Fig 4/Fig_4F_plotdata.txt",
            row.names = TRUE,
            col.names = TRUE,
            sep = "\t")

pdf("graphics/BRCA-E2F1-HMT-scatter.pdf",
    height = 2.5,
    width = 2.5)

ggplot(TCGA_BRCA_E2F1_plot_df,
       aes(x = E2F1, y = HMTs)) +
  geom_point(size = 0.4, 
             alpha = 0.3) + 
  theme_classic() + 
  ylab("Total HMT expression (corrected)") +
  xlab(substitute(italic(E2F1)~"activity (corrected)")) +
  theme(axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10)) + 
  scale_y_continuous(breaks = seq(-0.25, 0.5, 0.25),
                     limits = c(-0.3, 0.5),
                     expand = c(0,0)) +
  scale_x_continuous(breaks = seq(-1, 1.5, 0.5),
                     limits = c(-1.05, 1.6),
                     expand = c(0,0)) +
  annotate(geom = "text",
           x = -0.16,
           y = 0.475,
           label = substitute(italic(r)~"= 0.607")) +
  annotate(geom = "text",
           x = -0.16,
           y = 0.42,
           label = substitute(italic(p)~"= 1.16 x"~10^-108)) +
  geom_smooth(method = "lm", 
              se = FALSE,
              colour = "red")

dev.off()

#### E2F BINDING IN HUMAN HMTS Fig 4C ####

#code to evaluate potential regulation of HMTs in human
# transcription factor targets downloaded from Harmonizeome [https://maayanlab.cloud/Harmonizome/dataset/ENCODE+Transcription+Factor+Targets]

AttMatrixName <- read.table("input/gene_attribute_matrix.txt",
                            sep = "\t",
                            stringsAsFactors = F,
                            comment.char = "",
                            header = T)

AttMatrixMat <- data.frame(AttMatrixName[3:nrow(AttMatrixName), ])

colnames(AttMatrixMat) <- colnames(AttMatrixName)

AttMatrixMat<-AttMatrixMat[,-2]
colnames(AttMatrixMat)[1]<-"ID"

AttMatrixMat[, 3:ncol(AttMatrixMat)] <- apply(AttMatrixMat[ , 3:ncol(AttMatrixMat)], 2, function(x){
  as.numeric(as.character(x))})

#construct a list with the genes that are regulated by that TF as each entry

RegList <- list()

for(i in 3:ncol(AttMatrixMat)){
  
  RegList[[i - 2]] <- AttMatrixMat[which(AttMatrixMat[, i] == 1), 1]
  
  }


names(RegList) <- colnames(AttMatrixMat)[3:ncol(AttMatrixMat)]

#fisher test for each of the TFs to see whether the SET_HMTs is enriched

#for each TF construct the following table

#SET_HMTs    'SET_HMTs
#TF      ab      a-ab
#TF'     b-ab    tot-a-b+ab

TotGene <- nrow(AttMatrixMat)
#the total number of genes in question is 22819

TotSET_HMTs <- nrow(SET_HMTs)
#the total number of SET_HMTs is 38

TR_ORs <- c()
TR_pvals <- c()

for(i in 1:length(RegList)){

  a <- length(RegList[[i]])
  b <- TotSET_HMTs
  ab <- length(intersect(SET_HMTs$hgnc_symbol, RegList[[i]]))
    
  if(ab > 0){
  
    Fmat <- rbind(c(ab, a-ab), c(b-ab, TotGene-a-b-ab))
    FTest <- fisher.test(Fmat)
    TR_ORs <- c(TR_ORs, FTest$estimate)
    TR_pvals <- c(TR_pvals, FTest$p.value)

  } else {
    
    TR_ORs <- c(TR_ORs, NA)
    TR_pvals <- c(TR_pvals, NA)
    
  }
    
}

names(TR_ORs) <- names(RegList)
names(TR_pvals) <- names(RegList)

TR_ORs_forplot <- data.frame("OR" = TR_ORs,
                             "pval" = TR_pvals)

TR_ORs_forplot[, "colour"] <- "black"
TR_ORs_forplot["E2F1", "colour"] <- "red"
TR_ORs_forplot["E2F4", "colour"] <- "purple"

TR_ORs_forplot[, "size"] <- 1
TR_ORs_forplot["E2F1", "size"] <- 3
TR_ORs_forplot["E2F4", "size"] <- 3

write.table(TR_ORs_forplot,
            file = "plot_data/Fig 4/Fig_4C_plotdata.txt",
            row.names = TRUE,
            col.names = TRUE,
            sep = "\t")

pdf("graphics/TR_binding_to_HMTs_volcano.pdf",
    height = 2.5,
    width = 2.5)

ggplot(data = TR_ORs_forplot, aes(x = log2(OR), y = -log10(pval))) +
  geom_point(col = TR_ORs_forplot$colour,
             size = TR_ORs_forplot$size) +
  theme_classic() +
  theme(axis.text.y = element_text(colour = "black",
                                   size = 10),
        axis.text.x = element_text(colour = "black",
                                   size = 10),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.line.y = element_line(colour = "black"),
        plot.margin = margin(t = 5, r = 10, b = 5, l = 5)) +
  ylab(substitute(-log[10]~"(p value)")) +
  xlab(substitute(-log[2]~"(odds ratio)")) + 
  coord_cartesian(xlim = c(-2, 2)) + 
  annotate(geom = "text",
           x = 0.85,
           y = 1.9,
           label = "E2F1",
           col = "red") +
  annotate(geom = "text",
           x = 1.7,
           y = 1.5,
           label = "E2F4",
           col = "purple")

dev.off()

#### CO-REGULATION CORRELATION MODELLING Fig S2J ####

#let K be the concentration of a TF

#let A1...A10 be activated targets
#let B be a negative target

#let p1...p100 be the elasticities connecting K to each TF; for toy let p1...10 = p1
#let q1 be the elasticity connecting K to B

#each sample has a different K

#question is how many HMTs would have to be regulated?

#single simulation of 1 sample

SimFunc <- function(p1,
                    q1,
                    K,
                    noise){
  
  Pvec <- rnorm(40, mean = p1, sd = noise)
  Avec <- Pvec * K
  
  qSamp <- rnorm(1, mean = q1, sd = noise)
  Bvec <- qSamp * K * (-1)
  
  return(list(Avec,
              Bvec))
  
}

# will repeat 1000 iterations, simulating 500 samples each time
modelingiterations <- lapply(1:1000, function(arbitrary){
print(arbitrary)
  Kvec <- rnormTrunc(500, mean = 1, min = 0)
  
  OUTTotAB <- list()
  
  noisevec <- seq(0, 3, length = 20)

  for(N in 1:20){

    OutputA<-c()
    OutputB<-c()
    
    for(i in 1:500){

      SimFuncT <- SimFunc(p1 = 1, q1 = 1, noise = noisevec[N], K = Kvec[i])
      
      OutputA <- rbind(OutputA, SimFuncT[[1]])
      
      OutputB <-c (OutputB, SimFuncT[[2]])
      
    }
    
    totalCor <- c()
    
    for(i in 2:40){
      
      totalCor <- c(totalCor,
                    cor.test(apply(OutputA[, (1:i)], 1, sum), OutputB, method = "spearman")$estimate)
      
    }
    
    OUTTotAB[[N]] <- totalCor
    
  }
  
  return(OUTTotAB)
  
})

modelingiterations_fordata <- lapply(1:100, function(arbitrary){

  Kvec <- rnormTrunc(500, mean = 1, min = 0)
    
    OutputA<-c()
    OutputB<-c()
    
    # for(i in 1:500){
    for(i in 1:10){
      
      
      SimFuncT <- SimFunc(p1 = 1, q1 = 1, noise = 1.421053, K = Kvec[i])
      
      OutputA <- rbind(OutputA, SimFuncT[[1]])
      
      OutputB <-c (OutputB, SimFuncT[[2]])
      
    }
    
    totalCor <- c()
    
    for(j in 2:40){

      totalCor <- c(totalCor,
                    cor.test(apply(OutputA[, (1:j)], 1, sum), OutputB, method = "spearman")$estimate)
      
    }
  
  return(totalCor)
  
})

model_plotdf <- data.frame(highnoise = do.call(c, lapply(modelingiterations, function(x){x[[15]]})),
                           midnoise = do.call(c, lapply(modelingiterations, function(x){x[[10]]})),
                           lownoise = do.call(c, lapply(modelingiterations, function(x){x[[5]]})),
                           genes = rep(1:length(modelingiterations[[1]][[1]]), length = length(modelingiterations[[1]][[1]]) * length(modelingiterations)))

midnoise_plotdf <- sapply(modelingiterations, function(x){x[[10]]})
midnoise_melt <- melt(midnoise_plotdf)
midnoise_melt <- cbind(midnoise_melt, rep(1:length(modelingiterations[[1]][[1]]), length = length(modelingiterations[[1]][[1]]) * length(modelingiterations)))
colnames(midnoise_melt)[4]  <- "genes"

midnoise_melt <- melt(modelingiterations_fordata)
midnoise_melt <- cbind(midnoise_melt, rep(1:length(modelingiterations_fordata[[1]]), length = length(modelingiterations_fordata[[1]]) * length(modelingiterations_fordata)))
colnames(midnoise_melt)[3]  <- "genes"

saveRDS(midnoise_melt, "plot_data/midnoise_melt.rds")
midnoise_melt <- readRDS("plot_data/midnoise_melt.rds")

write.table(midnoise_melt,
            file = "plot_data/Fig S2/Fig_S2J_plotdata.txt",
            row.names = FALSE,
            col.names = TRUE,
            sep = ",")

pdf("graphics/modelling_coregulation_plot.pdf",
    height = 2.5,
    width = 2.5)

ggplot(aes(x = genes, y = value), data = midnoise_melt) + 
  geom_line(aes(group = Var2), alpha = 0.03) + 
  geom_smooth(col = "red", se = FALSE)  +
  theme_classic() + 
  theme(axis.text = element_text(colour = "black",
                                 size = 10)) + 
  ylab("A-B Spearman's"~rho) +
  xlab("Number of pooled A genes") +
  coord_cartesian(ylim = c(-0.45, 0.05)) +
  geom_hline(yintercept = 0,
             linetype = "dashed",
             colour = "grey")

dev.off()

#### AHCY KO or disruption ####

DZNep_counts_filelist <- list.files("input/Mouse_DZNep_counts/")

Mouse_DZNep_SRA <- read.table("input/Mouse_DZNep_SRA.txt",
                              sep = "\t",
                              header = TRUE)

counts_list <- lapply(DZNep_counts_filelist, function(thisfile){
  
  read.table(paste0("input/Mouse_DZNep_counts/", thisfile),
             sep = "\t",
             header = TRUE)
  
})

names(counts_list) <- str_remove(DZNep_counts_filelist, "_FCcounts.txt")

counts_df <- sapply(counts_list, function(thisone){
  
  thisone$temp_single.sam
  
})

row.names(counts_df) <- counts_list[[1]]$Geneid

AHCY_experiment_samples_WT <- Mouse_DZNep_SRA[Mouse_DZNep_SRA$genotype.variation == "WT", "Run"]
AHCY_experiment_samples_AHCY <- Mouse_DZNep_SRA[Mouse_DZNep_SRA$genotype.variation == "AHCY null", "Run"]

DZNep_experiment_samples_ctrl <- Mouse_DZNep_SRA[Mouse_DZNep_SRA$genotype.variation %in% c("DMSO treated"), "Run"]
DZNep_experiment_samples_dznep <- Mouse_DZNep_SRA[Mouse_DZNep_SRA$genotype.variation %in% c("Dznep treated"), "Run"]

AHCY_counts <- counts_df[, c(AHCY_experiment_samples_WT, AHCY_experiment_samples_AHCY)]

AHCY_coldata <- data.frame(row.names = c(AHCY_experiment_samples_WT, AHCY_experiment_samples_AHCY),
                           group = c(rep("WT", times = length(AHCY_experiment_samples_WT)),
                                     rep("AHCY", times = length(AHCY_experiment_samples_AHCY))))

AHCY_dds <- DESeqDataSetFromMatrix(countData = AHCY_counts,
                                   colData = AHCY_coldata,
                                   design = ~ group)

AHCY_dds$group <- factor(AHCY_dds$group, levels = c("WT","AHCY"))
AHCY_dds <- DESeq(AHCY_dds)                           

AHCY_dds_results <- results(AHCY_dds)

AHCY_statmat <- as.matrix(AHCY_dds_results$stat)
row.names(AHCY_statmat) <- row.names(AHCY_dds_results)

AHCY_normcounts <- counts(AHCY_dds, normalized = TRUE)

# get rid of missing values
AHCY_statmat <- AHCY_statmat[!is.na(AHCY_statmat[, 1]), ]

AHCY_statdf <- data.frame(stat = AHCY_statmat,
                          row.names = names(AHCY_statmat))

mouse_ensembl = useMart("ensembl", dataset = "mmusculus_gene_ensembl")

mouse_BM <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"), 
                  values = row.names(AHCY_dds_results),
                  filters = c("ensembl_gene_id"),
                  mart = mouse_ensembl)

potential_names <- mouse_BM[match(row.names(AHCY_statdf), mouse_BM$ensembl_gene_id), "external_gene_name"]
normcounts_potential_names <- mouse_BM[match(row.names(AHCY_normcounts), mouse_BM$ensembl_gene_id), "external_gene_name"]

# get rid of rows with no unique name
AHCY_statdf <- AHCY_statdf[!(duplicated(potential_names)|duplicated(potential_names, fromLast = TRUE)), ]
AHCY_normcounts_df <- AHCY_normcounts[!(duplicated(normcounts_potential_names)|duplicated(normcounts_potential_names, fromLast = TRUE)), ]

potential_names <- potential_names[!(duplicated(potential_names)|duplicated(potential_names, fromLast = TRUE))]
normcounts_potential_names <- normcounts_potential_names[!(duplicated(normcounts_potential_names)|duplicated(normcounts_potential_names, fromLast = TRUE))]

AHCY_statdf <- data.frame(row.names = potential_names,
                          stat = AHCY_statdf)

AHCY_normcounts_df <-  as.data.frame(AHCY_normcounts_df,
                                     row.names = normcounts_potential_names)

# decoupleR analysis

dorothea_mm_highconfidence <- dorothea_mm[dorothea_mm$confidence %in% c(LETTERS[1:4]), ]

for(i in 1:4){
  
  dorothea_mm_highconfidence[dorothea_mm_highconfidence$confidence == LETTERS[i], "mor"] <- dorothea_mm_highconfidence[dorothea_mm_highconfidence$confidence == LETTERS[i], "mor"] / i
  
}


colinear_TFs <- decoupleR::check_corr(network = dorothea_mm_highconfidence,
                                      .source = "tf")

colinear_to_exclude <- unique(unlist(colinear_TFs[colinear_TFs$correlation > 0.99, c("tf", "tf.2")]))

dorothea_mm_highconfidence <- dorothea_mm_highconfidence[!dorothea_mm_highconfidence$tf %in% colinear_to_exclude, ]

decouple_AHCY <- decouple(mat = AHCY_statdf,
                          network = dorothea_mm_highconfidence,
                          statistic = "mlm",
                          .source = "tf",
                          .target = "target"
)

decouple_AHCY_samples <- decouple(mat = AHCY_normcounts_df,
                                  network = dorothea_mm_highconfidence,
                                  statistic = "mlm",
                                  .source = "tf",
                                  .target = "target"
)

decouple_AHCY_mlm <- decouple_AHCY[decouple_AHCY$statistic == "mlm", ]

# now for DZNep

DZNep_counts <- counts_df[, c(DZNep_experiment_samples_ctrl, DZNep_experiment_samples_dznep)]

DZNep_coldata <- data.frame(row.names = c(DZNep_experiment_samples_ctrl, DZNep_experiment_samples_dznep),
                            group = c(rep("WT", times = length(DZNep_experiment_samples_ctrl)),
                                      rep("DZNep", times = length(DZNep_experiment_samples_dznep))))

DZNep_dds <- DESeqDataSetFromMatrix(countData = DZNep_counts,
                                    colData = DZNep_coldata,
                                    design = ~ group)

DZNep_dds$group <- factor(DZNep_dds$group, levels = c("WT","DZNep"))

DZNep_dds <- DESeq(DZNep_dds)                           

DZNep_dds_results <- results(DZNep_dds)

DZNep_statmat <- as.matrix(DZNep_dds_results$stat)
row.names(DZNep_statmat) <- row.names(DZNep_dds_results)

DZNep_normcounts <- counts(DZNep_dds, normalized = TRUE)

# get rid of missing values
DZNep_statmat <- DZNep_statmat[!is.na(DZNep_statmat[, 1]), ]

DZNep_statdf <- data.frame(stat = DZNep_statmat,
                           row.names = names(DZNep_statmat))

potential_names <- mouse_BM[match(row.names(DZNep_statdf), mouse_BM$ensembl_gene_id), "external_gene_name"]
normcounts_potential_names <- mouse_BM[match(row.names(DZNep_normcounts), mouse_BM$ensembl_gene_id), "external_gene_name"]

DZNep_normcounts_df <- DZNep_normcounts[!(duplicated(normcounts_potential_names)|duplicated(normcounts_potential_names, fromLast = TRUE)), ]

DZNep_statdf <- DZNep_statdf[!(duplicated(potential_names)|duplicated(potential_names, fromLast = TRUE)), ]

# get rid of rows with no unique name
potential_names <- potential_names[!(duplicated(potential_names)|duplicated(potential_names, fromLast = TRUE))]
normcounts_potential_names <- normcounts_potential_names[!(duplicated(normcounts_potential_names)|duplicated(normcounts_potential_names, fromLast = TRUE))]

DZNep_statdf <- data.frame(row.names = potential_names,
                           stat = DZNep_statdf)

DZNep_normcounts_df <-  as.data.frame(DZNep_normcounts_df,
                                      row.names = normcounts_potential_names)

decouple_DZNep <- decouple(mat = DZNep_statdf,
                           network = dorothea_mm_highconfidence,
                           statistic = "mlm",
                           .source = "tf",
                           .target = "target"
)

decouple_DZNep_samples <- decouple(mat = DZNep_normcounts_df,
                                   network = dorothea_mm_highconfidence,
                                   statistic = "mlm",
                                   .source = "tf",
                                   .target = "target"
)

decouple_DZNep_mlm <- decouple_DZNep[decouple_DZNep$statistic == "mlm", ]

decouple_AHCY_mlm_for_table <- decouple_AHCY_mlm[, c(2,3,5,6)]
colnames(decouple_AHCY_mlm_for_table) <- c("method", "transcriptional_regulator_name", "score", "p_value")
decouple_AHCY_mlm_for_table <- decouple_AHCY_mlm_for_table[order(decouple_AHCY_mlm_for_table$score, decreasing = TRUE), ]

decouple_DZNep_mlm_for_table <- decouple_DZNep_mlm[, c(2,3,5,6)]
colnames(decouple_DZNep_mlm_for_table) <- c("method", "transcriptional_regulator_name", "score", "p_value")
decouple_DZNep_mlm_for_table <- decouple_DZNep_mlm_for_table[order(decouple_DZNep_mlm_for_table$score, decreasing = TRUE), ]

# Fig S12 - Ahcy knockout

decouple_AHCY_samples_E2F <- decouple_AHCY_samples[decouple_AHCY_samples$statistic == "mlm" & decouple_AHCY_samples$source == "E2f1", c("source", "score", "condition")]

decouple_AHCY_samples_E2F[decouple_AHCY_samples_E2F$condition %in% AHCY_experiment_samples_WT, "condition"] <- "WT"
decouple_AHCY_samples_E2F[decouple_AHCY_samples_E2F$condition %in% AHCY_experiment_samples_AHCY, "condition"] <- "AHCY"

decouple_AHCY_samples_E2F$condition <- factor(decouple_AHCY_samples_E2F$condition, levels = c("WT", "AHCY"))

pdf("graphics/Mm_AHCYko_E2F1_box.pdf",
    width = 2.5,
    height = 2.5)

ggplot(decouple_AHCY_samples_E2F, aes(x = condition, y = score)) + 
  geom_boxplot(fill = c("grey", "red")) + 
  geom_jitter(alpha = 0.5,
              width = 0.2) +
  theme_classic() + 
  xlab("Genotype") + 
  ylab("E2f1 activity (AU)") 
# annotate(geom = "text",
#          label = substitute("p = 2.2 x"~10^-10),
#          x = 1, 
#          y = 3,
#          size = 3)

dev.off()

wilcox.test(unlist(decouple_AHCY_samples_E2F[decouple_AHCY_samples_E2F$condition == "WT" & decouple_AHCY_samples_E2F$source == "E2f1", "score"]),
            unlist(decouple_AHCY_samples_E2F[decouple_AHCY_samples_E2F$condition == "AHCY" & decouple_AHCY_samples_E2F$source == "E2f1", "score"])) 

wilcox.test(unlist(decouple_AHCY_samples_E2F[decouple_AHCY_samples_E2F$condition == "WT" & decouple_AHCY_samples_E2F$source == "E2f4", "score"]),
            unlist(decouple_AHCY_samples_E2F[decouple_AHCY_samples_E2F$condition == "AHCY" & decouple_AHCY_samples_E2F$source == "E2f4", "score"])) 

decouple_DZNep_samples_E2F <- decouple_DZNep_samples[decouple_DZNep_samples$statistic == "mlm" & decouple_DZNep_samples$source %in% c("E2f1"), c("source", "score", "condition")]

decouple_DZNep_samples_E2F[decouple_DZNep_samples_E2F$condition %in% DZNep_experiment_samples_ctrl, "condition"] <- "vehicle"
decouple_DZNep_samples_E2F[decouple_DZNep_samples_E2F$condition %in% DZNep_experiment_samples_dznep, "condition"] <- "DZNep"

decouple_DZNep_samples_E2F$condition <- factor(decouple_DZNep_samples_E2F$condition, levels = c("vehicle", "DZNep"))

pdf("graphics/Mm_DZNep_E2F1_box.pdf",
    width = 2.5,
    height = 2.5)

ggplot(decouple_DZNep_samples_E2F, aes(x = condition, y = score)) + 
  geom_boxplot(fill = c("grey", "red")) + 
  geom_jitter(alpha = 0.5,
              width = 0.2) +
  theme_classic() + 
  xlab("Treatment") + 
  ylab("E2f1 activity (AU)") +
  annotate(geom = "text",
           label = substitute("p = 2.17 x"~10^-3),
           x = 1, 
           y = 3)

dev.off()

wilcox.test(unlist(decouple_DZNep_samples_E2F[decouple_DZNep_samples_E2F$condition == "vehicle" & decouple_DZNep_samples_E2F$source == "E2f1", "score"]),
            unlist(decouple_DZNep_samples_E2F[decouple_DZNep_samples_E2F$condition == "DZNep" & decouple_DZNep_samples_E2F$source == "E2f1", "score"])) 

wilcox.test(unlist(decouple_DZNep_samples_E2F[decouple_DZNep_samples_E2F$condition == "vehicle" & decouple_DZNep_samples_E2F$source == "E2f4", "score"]),
            unlist(decouple_DZNep_samples_E2F[decouple_DZNep_samples_E2F$condition == "DZNep" & decouple_DZNep_samples_E2F$source == "E2f4", "score"])) 

# now for HMTs totals

mouseHMTs <- read.table("~/NNMT/input/mouse_HMTs.txt",
                        fill = TRUE,
                        header = TRUE,
                        sep = "\t")

AHCY_normcounts_df <- AHCY_normcounts_df[!apply(AHCY_normcounts_df, 1, function(x){any(is.na(x))}), ]

AHCY_HMTsums <- colSums(AHCY_normcounts_df[row.names(AHCY_normcounts_df) %in% mouseHMTs$Mouse.gene.name, ])

AHCY_HMT_plot <- data.frame(AHCY_HMTsums)

AHCY_HMT_plot[row.names(AHCY_HMT_plot) %in% AHCY_experiment_samples_WT, "condition"] <- "WT"
AHCY_HMT_plot[row.names(AHCY_HMT_plot) %in% AHCY_experiment_samples_AHCY, "condition"] <- "AHCY"

AHCY_HMT_plot$condition <- factor(AHCY_HMT_plot$condition, levels = c("WT", "AHCY"))

AHCY_NNMT <- unlist(AHCY_normcounts_df["Nnmt", ])

AHCY_NNMT_plot <- data.frame(AHCY_NNMT)

AHCY_NNMT_plot[row.names(AHCY_NNMT_plot) %in% AHCY_experiment_samples_WT, "condition"] <- "WT"
AHCY_NNMT_plot[row.names(AHCY_NNMT_plot) %in% AHCY_experiment_samples_AHCY, "condition"] <- "AHCY"

AHCY_NNMT_plot$condition <- factor(AHCY_NNMT_plot$condition, levels = c("WT", "AHCY"))

t.test(log10(AHCY_NNMT_plot[AHCY_NNMT_plot$condition == "WT", "AHCY_NNMT"]),
       log10(AHCY_NNMT_plot[AHCY_NNMT_plot$condition == "AHCY", "AHCY_NNMT"]))

t.test(log10(AHCY_HMT_plot[AHCY_NNMT_plot$condition == "WT", "AHCY_HMTsums"]),
       log10(AHCY_HMT_plot[AHCY_NNMT_plot$condition == "AHCY", "AHCY_HMTsums"]))$p.value

write.table(AHCY_NNMT_plot,
            file = "plot_data/Fig S12/Fig_S12A_MmAhcy_plotdata.txt",
            row.names = TRUE,
            col.names = TRUE,
            sep = "\t")

write.table(AHCY_HMT_plot,
            file = "plot_data/Fig S12/Fig_S12B_MmAhcy_plotdata.txt",
            row.names = TRUE,
            col.names = TRUE,
            sep = "\t")

pdf("graphics/Mm_AHCYko_NNMT_box.pdf",
    height = 2.5,
    width = 2.5)

ggplot(AHCY_NNMT_plot, aes(x = condition, y = AHCY_NNMT)) + 
  geom_boxplot(fill = c("grey", "red")) + 
  geom_jitter(alpha = 0.5,
              width = 0.2) +
  theme_classic() + 
  xlab("Genotype") +
  ylab("Nnmt expression\n(pseudocounts)") + 
  annotate(geom = "text",
           label = substitute("p = 7.19 x"~10^-16),
           x = 2,
           y = 75,
           size = 3)

dev.off()

pdf("graphics/Mm_AHCYko_HMT_box.pdf",
    height = 2.5,
    width = 2.5)

ggplot(AHCY_HMT_plot, aes(x = condition, y = AHCY_HMTsums)) + 
  geom_boxplot(fill = c("grey", "red")) + 
  geom_jitter(alpha = 0.5,
              width = 0.2) +
  theme_classic() + 
  xlab("Genotype") +
  ylab("Total HMT expression\n(pseudocounts)") + 
  annotate(geom = "text",
           label = substitute("p = 5.56 x"~10^-18),
           x = 1.1,
           y = 100000,
           size = 3)

dev.off()

DZNep_normcounts_df <- DZNep_normcounts_df[!apply(DZNep_normcounts_df, 1, function(x){any(is.na(x))}), ]

DZNep_HMTsums <- colSums(DZNep_normcounts_df[row.names(DZNep_normcounts_df) %in% mouseHMTs$Mouse.gene.name, ])

DZNep_HMT_plot <- data.frame(DZNep_HMTsums)
DZNep_HMT_plot[row.names(DZNep_HMT_plot) %in% DZNep_experiment_samples_ctrl, "condition"] <- "vehicle"
DZNep_HMT_plot[row.names(DZNep_HMT_plot) %in% DZNep_experiment_samples_dznep, "condition"] <- "DZNep"

DZNep_HMT_plot$condition <- factor(DZNep_HMT_plot$condition, levels = c("vehicle", "DZNep"))

t.test(log10(DZNep_HMT_plot[DZNep_HMT_plot$condition == "vehicle", "DZNep_HMTsums"]),
       log10(DZNep_HMT_plot[DZNep_HMT_plot$condition == "DZNep", "DZNep_HMTsums"]))$p.value

pdf("graphics/Mm_DZNep_HMT_box.pdf",
    height = 2.5,
    width = 2.5)

ggplot(DZNep_HMT_plot, aes(x = condition, y = DZNep_HMTsums)) + 
  geom_boxplot(fill = c("grey", "red")) + 
  geom_jitter(alpha = 0.5,
              width = 0.2) +
  theme_classic() + 
  xlab("Treatment") +
  ylab("Total HMT expression\n(pseudocounts)") +
  annotate(geom = "text",
           label = "p = 0.056",
           x = 2,
           y = 55000,
           size = 3)

dev.off()

DZNep_NNMT <- unlist(DZNep_normcounts_df["Nnmt", ])

DZNep_NNMT_plot <- data.frame(DZNep_NNMT)

DZNep_NNMT_plot[row.names(DZNep_NNMT_plot) %in% DZNep_experiment_samples_ctrl, "condition"] <- "vehicle"
DZNep_NNMT_plot[row.names(DZNep_NNMT_plot) %in% DZNep_experiment_samples_dznep, "condition"] <- "DZNep"

DZNep_NNMT_plot$condition <- factor(DZNep_NNMT_plot$condition, levels = c("vehicle", "DZNep"))

t.test(log10(DZNep_NNMT_plot[DZNep_NNMT_plot$condition == "vehicle", "DZNep_NNMT"]),
       log10(DZNep_NNMT_plot[DZNep_NNMT_plot$condition == "DZNep", "DZNep_NNMT"]))

pdf("graphics/Mm_DZNep_NNMT_box.pdf",
    height = 2.5,
    width = 2.5)

ggplot(DZNep_NNMT_plot, aes(x = condition, y = DZNep_NNMT)) + 
  geom_boxplot(fill = c("grey", "red")) + 
  geom_jitter(alpha = 0.5,
              width = 0.2) +
  theme_classic() + 
  xlab("Treatment") +
  ylab("Nnmt expression\n(pseudocounts)") +
  annotate(geom = "text",
           label = substitute("p = 3.94 x"~10^-7),
           x = 2,
           y = 100,
           size = 3)

dev.off()

write.table(DZNep_NNMT_plot,
            file = "plot_data/Fig S12/Fig_S12A_MmDZNep_plotdata.txt",
            row.names = TRUE,
            col.names = TRUE,
            sep = "\t")

write.table(DZNep_HMT_plot,
            file = "plot_data/Fig S12/Fig_S12B_MmDZNep_plotdata.txt",
            row.names = TRUE,
            col.names = TRUE,
            sep = "\t")

#### MRN normalise RAT DZnep counts ####

RatDZ_counts_filelist <- list.files("input/rat_dznep_counts/")

RatDZ_SRA <- read.table("input/Rat_DZNep_SRA.txt",
                        sep = "\t",
                        header = TRUE)

RatDZ_counts_list <- lapply(RatDZ_counts_filelist, function(thisfile){
  
  read.table(paste0("input/rat_dznep_counts/", thisfile),
             sep = "\t",
             header = TRUE)
  
})

names(RatDZ_counts_list) <- str_remove(RatDZ_counts_filelist, "_FCcounts.txt")

RatDZ_counts_df <- sapply(RatDZ_counts_list, function(thisone){
  
  thisone[, str_detect(colnames(thisone), "\\.sam")]
  
})

row.names(RatDZ_counts_df) <- RatDZ_counts_list[[1]]$Geneid

RatDZ_experiment_samples_ctrl <- RatDZ_SRA[str_detect(RatDZ_SRA$Treatment, "DMSO"), "Run"]
RatDZ_experiment_samples_dznep <- RatDZ_SRA[str_detect(RatDZ_SRA$Treatment, "DZNe"), "Run"]

RatDZ_coldata <- data.frame(row.names = colnames(RatDZ_counts_df),
                            group = RatDZ_SRA[match(colnames(RatDZ_counts_df), RatDZ_SRA$Run), "Treatment"])

RatDZ_coldata[, "group"] <- str_remove(RatDZ_coldata$group, " at 1 uM")

RatDZ_dds <- DESeqDataSetFromMatrix(countData = RatDZ_counts_df,
                                    colData = RatDZ_coldata,
                                    design = ~ group)

RatDZ_dds <- DESeq(RatDZ_dds)                           

RatDZ_dds_results <- results(RatDZ_dds)

RatDZ_statmat <- as.matrix(RatDZ_dds_results$stat)
row.names(RatDZ_statmat) <- row.names(RatDZ_dds_results)

RatDZ_normcounts <- counts(RatDZ_dds, normalized = TRUE)

# get rid of missing values
RatDZ_statmat <- RatDZ_statmat[!is.na(RatDZ_statmat[, 1]), ]

RatDZ_statdf <- data.frame(stat = RatDZ_statmat,
                           row.names = names(RatDZ_statmat))

rat_ensembl = useEnsembl("ensembl", dataset = "rnorvegicus_gene_ensembl", mirror = "useast")

rat_BM <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "mmusculus_homolog_ensembl_gene", "hsapiens_homolog_ensembl_gene", "mmusculus_homolog_associated_gene_name"), 
                values = row.names(RatDZ_normcounts),
                filters = c("ensembl_gene_id"),
                mart = rat_ensembl)

rat_HMTs <- rat_BM[match(mouseHMTs$Mouse.gene.stable.ID, rat_BM$mmusculus_homolog_ensembl_gene), "ensembl_gene_id"]
rat_HMTs <- rat_HMTs[!is.na(rat_HMTs)]

rat_NNMT_ensembl <- rat_BM[rat_BM$external_gene_name == "Nnmt", "ensembl_gene_id"]

RatDZ_HMTsums <- colSums(RatDZ_normcounts[row.names(RatDZ_normcounts) %in% rat_HMTs, ])

RatDZ_HMT_plot <- data.frame(RatDZ_HMTsums)

# the samples are clearly reversed. The figures in the paper show the opposite trends for most notably Fgb and other fibrinogen genes.
# so I will go with the reversal.

# NB here I correct the switched labels
RatDZ_HMT_plot[row.names(RatDZ_HMT_plot) %in% RatDZ_experiment_samples_ctrl, "condition"] <- "DZNep"
RatDZ_HMT_plot[row.names(RatDZ_HMT_plot) %in% RatDZ_experiment_samples_dznep, "condition"] <- "vehicle"

RatDZ_HMT_plot$condition <- factor(RatDZ_HMT_plot$condition, levels = c("vehicle", "DZNep"))

t.test(log10(RatDZ_HMT_plot[RatDZ_HMT_plot$condition == "vehicle", "RatDZ_HMTsums"]),
       log10(RatDZ_HMT_plot[RatDZ_HMT_plot$condition == "DZNep", "RatDZ_HMTsums"]))

RatDZ_NNMT <- unlist(RatDZ_normcounts[rat_NNMT_ensembl, ])

RatDZ_NNMT_plot <- data.frame(RatDZ_NNMT)

# NB here I correct the switched labels
RatDZ_NNMT_plot[row.names(RatDZ_NNMT_plot) %in% RatDZ_experiment_samples_ctrl, "condition"] <- "DZNep"
RatDZ_NNMT_plot[row.names(RatDZ_NNMT_plot) %in% RatDZ_experiment_samples_dznep, "condition"] <- "vehicle"

RatDZ_NNMT_plot$condition <- factor(RatDZ_NNMT_plot$condition, levels = c("vehicle", "DZNep"))

t.test(log10(RatDZ_NNMT_plot[RatDZ_NNMT_plot$condition == "vehicle", "RatDZ_NNMT"]),
       log10(RatDZ_NNMT_plot[RatDZ_NNMT_plot$condition == "DZNep", "RatDZ_NNMT"]))

pdf("graphics/Rn_DZNep_NNMT_box.pdf",
    width = 2.5,
    height = 2.5)

ggplot(RatDZ_NNMT_plot, aes(x = condition, y = RatDZ_NNMT)) + 
  geom_boxplot(fill = c("grey", "red")) + 
  geom_jitter(alpha = 0.5,
              width = 0.2) +
  theme_classic() + 
  xlab("Treatment") +
  ylab("Nnmt expression\n(pseudocounts)") + 
  annotate(geom = "text",
           label = substitute("p = 4.8 x"~10^-04),
           x = 2, 
           y = 75,
           size = 3)

dev.off()

pdf("graphics/Rn_DZNep_HMT_box.pdf",
    width = 2.5,
    height = 2.5)

ggplot(RatDZ_HMT_plot, aes(x = condition, y = RatDZ_HMTsums)) + 
  geom_boxplot(fill = c("grey", "red")) + 
  geom_jitter(alpha = 0.5,
              width = 0.2) +
  theme_classic() + 
  xlab("Treatment") +
  ylab("Total HMT expression\n(pseudocounts)") + 
  annotate(geom = "text",
           label = "p = 0.20",
           x = 1, 
           y = 35000,
           size = 3)

dev.off()

write.table(RatDZ_NNMT_plot,
            file = "plot_data/Fig S12/Fig_S12A_RnDZNep_plotdata.txt",
            row.names = TRUE,
            col.names = TRUE,
            sep = "\t")

write.table(RatDZ_HMT_plot,
            file = "plot_data/Fig S12/Fig_S12B_RnDZNep_plotdata.txt",
            row.names = TRUE,
            col.names = TRUE,
            sep = "\t")

dorothea_rn_highconfidence <- dorothea_mm_highconfidence

dorothea_rn_highconfidence$tf <- rat_BM[match(dorothea_rn_highconfidence$tf, rat_BM$mmusculus_homolog_associated_gene_name), "ensembl_gene_id"]
dorothea_rn_highconfidence$target <- rat_BM[match(dorothea_rn_highconfidence$target, rat_BM$mmusculus_homolog_associated_gene_name), "ensembl_gene_id"]

dorothea_rn_highconfidence <- dorothea_rn_highconfidence[apply(dorothea_rn_highconfidence, 1, function(x){!any(is.na(x))}), ]

dorothea_rn_highconfidence <- dorothea_rn_highconfidence[!(duplicated(paste0(dorothea_rn_highconfidence$tf, "-", dorothea_rn_highconfidence$target))|duplicated(paste0(dorothea_rn_highconfidence$tf, "-", dorothea_rn_highconfidence$target), fromLast = TRUE)), ]

decouple_RatDZ <- decouple(mat = RatDZ_statdf,
                           network = dorothea_rn_highconfidence,
                           statistic = "mlm",
                           .source = "tf",
                           .target = "target"
)

decouple_RatDZ$source <- rat_BM[match(decouple_RatDZ$source, rat_BM$ensembl_gene_id), "external_gene_name"]

decouple_RatDZ_mlm <- decouple_RatDZ[decouple_RatDZ$statistic == "mlm", ]

decouple_RatDZ_samples <- decouple(mat = RatDZ_normcounts,
                                   network = dorothea_rn_highconfidence,
                                   statistic = "mlm",
                                   .source = "tf",
                                   .target = "target"
)

decouple_RatDZ_samples$source <- rat_BM[match(decouple_RatDZ_samples$source, rat_BM$ensembl_gene_id), "external_gene_name"]

decouple_RatDZ_samples_E2F <- decouple_RatDZ_samples[decouple_RatDZ_samples$statistic == "mlm" & decouple_RatDZ_samples$source %in% c("E2f1"), c("source", "score", "condition")]

# NB correcting for swapped samples
decouple_RatDZ_samples_E2F[decouple_RatDZ_samples_E2F$condition %in% RatDZ_experiment_samples_ctrl, "condition"] <- "DZNep"
decouple_RatDZ_samples_E2F[decouple_RatDZ_samples_E2F$condition %in% RatDZ_experiment_samples_dznep, "condition"] <- "vehicle"

decouple_RatDZ_samples_E2F$condition <- factor(decouple_RatDZ_samples_E2F$condition, levels = c("vehicle", "DZNep"))

pdf("graphics/Rn_DZNep_E2f1_box.pdf",
    height = 2.5,
    width = 2.5)

ggplot(decouple_RatDZ_samples_E2F, aes(x = condition, y = score)) + 
  geom_boxplot(fill = c("grey", "red")) + 
  geom_jitter(alpha = 0.5,
              width = 0.2) +
  theme_classic() + 
  xlab("Treatment") + 
  ylab("E2f1 activity (AU)") +
  annotate(geom = "text",
           label = "p = 0.1",
           x = 1,
           y = 0,
           size = 3)

dev.off()

wilcox.test(unlist(decouple_RatDZ_samples_E2F[decouple_RatDZ_samples_E2F$condition == "vehicle" & decouple_RatDZ_samples_E2F$source == "E2f1", "score"]),
            unlist(decouple_RatDZ_samples_E2F[decouple_RatDZ_samples_E2F$condition == "DZNep" & decouple_RatDZ_samples_E2F$source == "E2f1", "score"])) 

decouple_RatDZ_mlm_for_table <- decouple_RatDZ_mlm[, c(2,3,5,6)]
colnames(decouple_RatDZ_mlm_for_table) <- c("method", "transcriptional_regulator_name", "score", "p_value")
decouple_RatDZ_mlm_for_table <- decouple_RatDZ_mlm_for_table[order(decouple_RatDZ_mlm_for_table$score, decreasing = TRUE), ]

# Table S8

openxlsx::write.xlsx(list("MEFs_AHCYko_diffTFactivity" = decouple_AHCY_mlm_for_table,
                          "MEFs_DZNep_diffTFactivity" = decouple_DZNep_mlm_for_table,
                          "RatHSCs_DZNep_diffTFactivity" = decouple_RatDZ_mlm_for_table),
                     "output/AHCY_decouple_tables.xlsx",
                     row.names = FALSE)

#### Glyr1 in mouse study ####

AHCY_GLYR1 <- AHCY_normcounts_df["Glyr1", ]

AHCY_GLYR1_plot <- data.frame(GLYR1 = unlist(AHCY_GLYR1))

AHCY_GLYR1_plot[row.names(AHCY_GLYR1_plot) %in% AHCY_experiment_samples_WT, "condition"] <- "WT"
AHCY_GLYR1_plot[row.names(AHCY_GLYR1_plot) %in% AHCY_experiment_samples_AHCY, "condition"] <- "AHCY"

AHCY_GLYR1_plot$condition <- factor(AHCY_GLYR1_plot$condition, levels = c("WT", "AHCY"))

write.table(AHCY_GLYR1_plot,
            file = "plot_data/Fig S13/Fig_S13A_AhcyKOGlyr1Exp_plotdata.txt",
            col.names = TRUE,
            row.names = TRUE,
            sep = "\t")

# Fig S13A

pdf("graphics/MEF_AHCY_GLYR1exp.pdf",
    height = 2,
    width = 2)

ggplot(AHCY_GLYR1_plot, aes(x = condition, y = GLYR1)) + 
  geom_boxplot(fill = c("grey", "red")) + 
  geom_jitter(alpha = 0.5,
              width = 0.2) +
  theme_classic() + 
  xlab("") +
  ylab("")

dev.off()

DZNep_GLYR1 <- DZNep_normcounts_df["Glyr1", ]

DZNep_GLYR1_plot <- data.frame(GLYR1 = unlist(DZNep_GLYR1))

DZNep_GLYR1_plot[row.names(DZNep_GLYR1_plot) %in% DZNep_experiment_samples_ctrl, "condition"] <- "WT"
DZNep_GLYR1_plot[row.names(DZNep_GLYR1_plot) %in% DZNep_experiment_samples_dznep, "condition"] <- "DZNep"

DZNep_GLYR1_plot$condition <- factor(DZNep_GLYR1_plot$condition, levels = c("WT", "DZNep"))

write.table(DZNep_GLYR1_plot,
            file = "plot_data/Fig S13/Fig_S13A_DZNepGlyr1Exp_plotdata.txt",
            col.names = TRUE,
            row.names = TRUE,
            sep = "\t")

# Fig S13A

pdf("graphics/MEF_DNep_GLYR1exp.pdf",
    height = 2,
    width = 2)

ggplot(DZNep_GLYR1_plot, aes(x = condition, y = GLYR1)) + 
  geom_boxplot(fill = c("grey", "red")) + 
  geom_jitter(alpha = 0.5,
              width = 0.2) +
  theme_classic() + 
  xlab("") +
  ylab("")

dev.off()

decouple_AHCY_samples_GLYR1 <- decouple_AHCY_samples[decouple_AHCY_samples$statistic == "mlm" & decouple_AHCY_samples$source == "Glyr1", c("source", "score", "condition")]

decouple_AHCY_samples_GLYR1[decouple_AHCY_samples_GLYR1$condition %in% AHCY_experiment_samples_WT, "condition"] <- "WT"
decouple_AHCY_samples_GLYR1[decouple_AHCY_samples_GLYR1$condition %in% AHCY_experiment_samples_AHCY, "condition"] <- "AHCY"

decouple_AHCY_samples_GLYR1$condition <- factor(decouple_AHCY_samples_GLYR1$condition, levels = c("WT", "AHCY"))

# Fig S13A

write.table(decouple_AHCY_samples_GLYR1,
            file = "plot_data/Fig S13/Fig_S13A_AhcyKOGlyr1Act_plotdata.txt",
            col.names = TRUE,
            row.names = FALSE,
            sep = "\t")

pdf("graphics/MEF_AHCY_GLYR1act.pdf",
    height = 2,
    width = 2)

ggplot(decouple_AHCY_samples_GLYR1, aes(x = condition, y = score)) + 
  geom_boxplot(fill = c("grey", "red")) + 
  geom_jitter(alpha = 0.5,
              width = 0.2) +
  theme_classic() + 
  xlab("") + 
  ylab("")

dev.off()

wilcox.test(unlist(decouple_AHCY_samples_GLYR1[decouple_AHCY_samples_GLYR1$condition == "WT" & decouple_AHCY_samples_GLYR1$source == "Glyr1", "score"]),
            unlist(decouple_AHCY_samples_GLYR1[decouple_AHCY_samples_GLYR1$condition == "AHCY" & decouple_AHCY_samples_GLYR1$source == "Glyr1", "score"])) 

decouple_DZNep_samples_GLYR1 <- decouple_DZNep_samples[decouple_DZNep_samples$statistic == "mlm" & decouple_DZNep_samples$source == "Glyr1", c("source", "score", "condition")]

decouple_DZNep_samples_GLYR1[decouple_DZNep_samples_GLYR1$condition %in% DZNep_experiment_samples_ctrl, "condition"] <- "vehicle"
decouple_DZNep_samples_GLYR1[decouple_DZNep_samples_GLYR1$condition %in% DZNep_experiment_samples_dznep, "condition"] <- "DZNep"

decouple_DZNep_samples_GLYR1$condition <- factor(decouple_DZNep_samples_GLYR1$condition, levels = c("vehicle", "DZNep"))

write.table(decouple_DZNep_samples_GLYR1,
            file = "plot_data/Fig S13/Fig_S13A_DZNepGlyr1Act_plotdata.txt",
            col.names = TRUE,
            row.names = FALSE,
            sep = "\t")

# Fig S13A

pdf("graphics/MEF_DZNep_GLYR1act.pdf",
    height = 2,
    width = 2)

ggplot(decouple_DZNep_samples_GLYR1, aes(x = condition, y = score)) + 
  geom_boxplot(fill = c("grey", "red")) + 
  geom_jitter(alpha = 0.5,
              width = 0.2) +
  theme_classic() + 
  xlab("") + 
  ylab("")

dev.off()

wilcox.test(unlist(decouple_DZNep_samples_GLYR1[decouple_DZNep_samples_GLYR1$condition == "WT" & decouple_DZNep_samples_GLYR1$source == "Glyr1", "score"]),
            unlist(decouple_DZNep_samples_GLYR1[decouple_DZNep_samples_GLYR1$condition == "DZNep" & decouple_DZNep_samples_GLYR1$source == "Glyr1", "score"])) 

#### nuclear localisations - Table S9 ####

subcellular_localisation <- read.table("~/NNMT_manuscript/REVISION/subcellular_location.tsv",
                                       sep = "\t",
                                       header = TRUE)

transsulfuration_and_GSH_enzymes <- c("AHCY",
                                      "CBS",
                                      "CTH",
                                      "GCLC",
                                      "GCLM",
                                      "GSS")

View(getBM(attributes = c("hgnc_symbol",
                          "ensembl_gene_id",
                          "description"),
           filters = "hgnc_symbol",
           values = transsulfuration_and_GSH_enzymes,
           mart = ensembl))

View(subcellular_localisation[subcellular_localisation$Gene.name %in% transsulfuration_and_GSH_enzymes, c(2:5)])

#### GLYR1 in TCGA #####

## first recreate correlation and partial correlations  ----------------------------

TCGA_activities_files <- list.files("~/TRestimations/TCGA/output/output") 
TCGA_activities_files <- TCGA_activities_files[str_detect(TCGA_activities_files, pattern = "^TCGA-[A-Z]+") & str_detect(TCGA_activities_files, "crosscancer")]

TCGA_activities_list <- lapply(TCGA_activities_files, function(x){
  
  readRDS(paste0("~/TRestimations/TCGA/output/output/", x))
  
})

names(TCGA_activities_list) <- str_remove(TCGA_activities_files, pattern = "_TFcrosscancerconsensus.rds")

TCGA_MOR_HMT_and_NNMT <- lapply(TCGA_MOR_list, function(x){
  
  totalHMTs <- colSums(x[SET_HMTs$ensembl_gene_id, ])
  
  singlegenes <- x[c(SET_HMTs$ensembl_gene_id, NNMT_ensembl), ]
  
  output_df <- rbind(singlegenes, totalHMTs)
  
})

TCGA_MOR_HMT_and_NNMT_corrected <- lapply(TCGA_MOR_HMT_and_NNMT, TCGA.apply.lm.to.counts)

# correct TR activities too for confounders e.g. race, sex, age, cancer stage
TCGA_TFactivities_corrected <- lapply(TCGA_activities_list, TCGA.apply.lm.to.counts)

TCGA_MOR_HMT_and_NNMT_corrected_rankpercent <- lapply(TCGA_MOR_HMT_and_NNMT_corrected, function(x){
  
  as.data.frame(t(apply(x, 1, rank)/(ncol(x) + 1)))
  
})

TCGA_TFactivities_corrected_rankpercent <- lapply(TCGA_TFactivities_corrected, function(x){
  
  as.data.frame(t(apply(x, 1, rank)/(ncol(x) + 1)))
  
})

TCGA_TFactivities_corrected_rankpercent_df <- do.call(cbind, TCGA_TFactivities_corrected_rankpercent)

colnames(TCGA_TFactivities_corrected_rankpercent_df) <- str_remove(colnames(TCGA_TFactivities_corrected_rankpercent_df), 
                                                                   pattern = "^TCGA-[A-Z]{3,4}\\.")

iterations_out2 <- sapply(1:1000, function(i){

  # smallest number of samples of any any cancer type is 36 (CHOL)
  tempsample <- lapply(TCGA_MOR_HMT_and_NNMT_corrected_rankpercent, function(x){
    
    x[, sample(colnames(x), 36)]
    
  })
  
  tempsample_df <- do.call(cbind, tempsample)
  
  apply(TCGA_TFactivities_corrected_rankpercent_df, 1, function(y){
    
    cor.test(unlist(tempsample_df["totalHMTs", ]), y[colnames(tempsample_df)])$estimate
    
  })
  
})

TCGA_median_TFcorrs <- rowMedians(iterations_out2)
names(TCGA_median_TFcorrs) <- row.names(iterations_out2)

TCGANNMTiterations_out <- sapply(1:1000, function(i){
  
  # smallest number of samples of any any cancer type is 36 (CHOL)
  tempsample <- lapply(TCGA_MOR_HMT_and_NNMT_corrected_rankpercent, function(x){
    
    x[, sample(colnames(x), 36)]
    
  })
  
  tempsample_df <- do.call(cbind, tempsample)
  
  apply(TCGA_TFactivities_corrected_rankpercent_df, 1, function(y){
    
    cor.test(unlist(tempsample_df[NNMT_ensembl, ]), y[colnames(tempsample_df)])$estimate
    
  })
  
})

TCGA_NNMT_median_TFcorrs <- rowMedians(TCGANNMTiterations_out)
names(TCGA_NNMT_median_TFcorrs) <- row.names(TCGANNMTiterations_out)

TCGA_NNMT_median_TFcorrs [order(TCGA_NNMT_median_TFcorrs )]

TCGA_ppcor_iterations_out <- lapply(1:1000, function(i){
  
  # smallest number of samples of any any cancer type is 36 (CHOL)
  tempsample <- lapply(TCGA_MOR_HMT_and_NNMT_corrected_rankpercent, function(x){
    
    x[, sample(colnames(x), 36)]
    
  })
  
  tempsample_df <- do.call(cbind, tempsample)
  
  iteration_out <- apply(TCGA_TFactivities_corrected_rankpercent_df, 1, function(y){
    
    temp_mat <- rbind(tempsample_df[c(NNMT_ensembl, "totalHMTs"), ], y[colnames(tempsample_df)])
    temp_pcor <- pcor(t(temp_mat), method = "spearman")
    
    tempoutput <- c(temp_pcor$estimate[NNMT_ensembl, "totalHMTs"], temp_pcor$estimate[NNMT_ensembl, 3], temp_pcor$estimate["totalHMTs", 3])
    names(tempoutput) <- c("HMT-NNMT", "NNMT-TR", "HMT-TR")
    
    return(tempoutput)  
    
  })
  
  return(t(iteration_out))
  
})

saveRDS(TCGA_ppcor_iterations_out, "~/TRestimations/GTEX/output/TCGA_ppcor_iterations_out.rds")
# TCGA_ppcor_iterations_out <- readRDS("~/TRestimations/GTEX/output/TCGA_ppcor_iterations_out.rds")

TCGA_ppcor_array <- array(as.numeric(unlist(TCGA_ppcor_iterations_out)), dim=c(nrow(TCGA_ppcor_iterations_out[[1]]), ncol(TCGA_ppcor_iterations_out[[1]]), length(TCGA_ppcor_iterations_out)))
row.names(TCGA_ppcor_array) <- row.names(TCGA_ppcor_iterations_out[[1]])
colnames(TCGA_ppcor_array) <- colnames(TCGA_ppcor_iterations_out[[1]])

TCGA_ppcor_medians <- data.frame("NNMT_HMT_pcor" = rowMedians(TCGA_ppcor_array[, 1, 1:100]),
                                 "NNMT_TR_pcor" = rowMedians(TCGA_ppcor_array[, 2, 1:100]),
                                 "NNMT_TR_fullcor" = TCGA_NNMT_median_TFcorrs,
                                 "HMT_TR_fullcor" = TCGA_median_TFcorrs,
                                 "HMT_TR_pcor" = rowMedians(TCGA_ppcor_array[, 3, 1:100])
)

TCGA_ppcor_medians[!row.names(TCGA_ppcor_medians) %in% SET_HMTs$hgnc_symbol, ]

write.table(TCGA_ppcor_medians[!row.names(TCGA_ppcor_medians) %in% SET_HMTs$hgnc_symbol, ],
            file = "plot_data/Fig S13/Fig_S13B_plotdata.txt",
            row.names = TRUE,
            col.names = TRUE,
            sep = "\t"
            )

# Fig S13B

pdf("graphics/NNMT_HMT_TFactivity_ppcor_violin.pdf",
    width = 2,
    height = 3)

ggplot(data = TCGA_ppcor_medians[!row.names(TCGA_ppcor_medians) %in% SET_HMTs$hgnc_symbol, ], aes(x = 1, y = NNMT_HMT_pcor)) + 
  geom_violin(fill = "grey") + 
  geom_jitter(width = 0.05, alpha = 0.1) +
  theme_classic() + 
  geom_point(aes(x = 1, y = -0.228), col = "red", size = 3) + 
  coord_cartesian(ylim = c(-0.45, -0.2)) + 
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) + 
  ylab(italic(NNMT)~"-HMT partial correlation") +
  geom_hline(yintercept = -0.365,
             linetype = "dashed")

dev.off()

### GLYR 1 in RB mutants --------------------------------------------------

Rb_mutant_wt_GLYR1actranks <- lapply(cancers_with_rbs, function(thiscancer){
  
  mutant <- TCGA_TFactivities_corrected_rankpercent_df["GLYR1", str_remove(colnames(TCGA_TFactivities_corrected_rankpercent_df), pattern = "-[:alnum:]+-[:alnum:]+-[:digit:]+$") 
                                                       %in%
                                                         str_remove(Rbdeleterious_SAMPIDs_bycancer[[thiscancer]], pattern = "-[:alnum:]+-[:alnum:]+-[:digit:]+$")]
  
  wildtype <- TCGA_TFactivities_corrected_rankpercent_df["GLYR1", str_remove(colnames(TCGA_TFactivities_corrected_rankpercent_df), pattern = "-[:alnum:]+-[:alnum:]+-[:digit:]+$") 
                                                         %in%
                                                           str_remove(wtRB_SAMPIDs_bycancer[[thiscancer]], pattern = "-[:alnum:]+-[:alnum:]+-[:digit:]+$")]
  
  list("mutant" = mutant,
       "wildtype" = wildtype)
  
})

names(Rb_mutant_wt_GLYR1actranks) <- cancers_with_rbs


# do pan-cancer analysis 
# not all of those with >10 samples in 
RBmutant_select_cancers <- names(sapply(lapply(Rb_mutant_wt_GLYR1actranks, function(x){x[["mutant"]]}), length)[sapply(lapply(Rb_mut_wt_HMTranks, function(x){x[["mutant"]]}), length) > 9])

# the direct comparison is significant in UCEC and LUSC. mean and median is higher in 12/13 cancers
# likely that UCEDC wouldnt survive p values correction
signaltest <- lapply(RBmutant_select_cancers, function(x){
  thisname <- x
  x <- Rb_mutant_wt_GLYR1actranks[[x]]
  wt_temp <- wilcox.test(unlist(x[["mutant"]]), unlist(x[["wildtype"]]))
  boxplot(unlist(x[["wildtype"]]), unlist(x[["mutant"]]),
          labels = c("WT", "MUT"),
          main = thisname)
  print(thisname)
  print(paste0("wildtype mean is", mean(unlist(x[["wildtype"]]))))
  print(paste0("mutant mean is", mean(unlist(x[["mutant"]]))))
  print(paste0("wildtype MEDIAN is", median(unlist(x[["wildtype"]]))))
  print(paste0("mutant MEDIAN is", median(unlist(x[["mutant"]]))))
  
  message(paste0("Wilcox p: ", wt_temp$p.value))
  
})

sapply(Rb_mutant_wt_GLYR1actranks, function(x){
  
  outvec <- c(mut_median = mean(unlist(x[["mutant"]])), 
              wt_median = mean(unlist(x[["wildtype"]])))
  
  return(outvec)
  
})[, RBmutant_select_cancers]

RBmutant_pancancer_GLYR1actmedians <- data.frame(t(sapply(1:1000, function(i){
  
  Rb_select_vec <- unlist(lapply(RBmutant_select_cancers, function(X){
    
    sample(Rb_mutant_wt_GLYR1actranks[[X]][["mutant"]], 10)
    
  }))
  
  WT_select_vec <- unlist(lapply(RBmutant_select_cancers, function(X){
    
    sample(Rb_mutant_wt_GLYR1actranks[[X]][["wildtype"]], 10)
    
  }))
  
  output_vec <- c()
  
  output_vec["mut_rank_median"] <- median(Rb_select_vec)
  output_vec["WT_rank_median"] <- median(WT_select_vec)
  
  output_vec["p_value"] <-   wilcox.test(Rb_select_vec, WT_select_vec)$p.value
  
  return(output_vec)
  
})))

colnames(RBmutant_pancancer_GLYR1actmedians)[1:2] <- c("Deleterious", "Wildtype")

# t.test p value is around 0
t.test(RBmutant_pancancer_GLYR1actmedians$Deleterious, RBmutant_pancancer_GLYR1actmedians$Wildtype)$p.value

RBmutant_pancancer_GLYR1actmedians_melt_forplot <- reshape2::melt(RBmutant_pancancer_GLYR1actmedians[, c("Wildtype", "Deleterious")])

saveRDS(RBmutant_pancancer_GLYR1actmedians_melt_forplot,
        "plot_data/RBmutant_pancancer_GLYR1actmedians_melt_forplot.rds")
RBmutant_pancancer_GLYR1actmedians_melt_forplot <- readRDS("plot_data/RBmutant_pancancer_GLYR1actmedians_melt_forplot.rds")

write.table(RBmutant_pancancer_GLYR1actmedians_melt_forplot,
            file = "plot_data/Fig S13/Fig_S13D_plotdata.txt",
            row.names = FALSE,
            col.names = TRUE,
            sep = "\t")

# Fig S13D

pdf("graphics/GLYRactmedianranks_RB1mut.pdf",
    height = 2.5,
    width = 2.5)

ggplot(RBmutant_pancancer_GLYR1actmedians_melt_forplot, aes(x = variable, y = value)) + 
  geom_boxplot(fill = c("gray", "red")) +
  theme_classic() + 
  ylab("Median GLYR1 activity %ile") + 
  xlab(italic(RB1)~"mutation status") + 
  theme(axis.text.y = element_text(colour = "black",
                                   size = 10),
        axis.text.x = element_text(colour = "black",
                                   size = 8)) + 
  coord_cartesian(ylim = c(0.3, 0.8)) + 
  geom_hline(yintercept = 0.5,
             linetype = "dashed",
             col = "gray") + 
  geom_text(aes(x = 1.5, 
                y = 0.8, 
                label = "p%~~% 0"),
            parse = TRUE)

dev.off()

# partial corr ------------------

all_Rb_ranks <- do.call(c, lapply(Rb_mut_wt_HMTranks, function(X){X[["mutant"]]}))
names(all_Rb_ranks) <- str_remove(names(all_Rb_ranks), pattern = "^TCGA-[A-Z]{2,4}\\.")

all_wt_ranks <- do.call(c, lapply(Rb_mut_wt_HMTranks, function(X){X[["wildtype"]]}))
names(all_wt_ranks) <- str_remove(names(all_wt_ranks), pattern = "^TCGA-[A-Z]{2,4}\\.")

all_Rb_NNMT_ranks <- do.call(c, lapply(Rb_mutant_wt_NNMTranks, function(X){X[["mutant"]]}))
names(all_Rb_NNMT_ranks) <- str_remove(names(all_Rb_NNMT_ranks), pattern = "^TCGA-[A-Z]{2,4}\\.")

all_wt_NNMT_ranks <- do.call(c, lapply(Rb_mutant_wt_NNMTranks, function(X){X[["wildtype"]]}))
names(all_wt_NNMT_ranks) <- str_remove(names(all_wt_NNMT_ranks), pattern = "^TCGA-[A-Z]{2,4}\\.")

all_Rb_GLYR1act_ranks <- do.call(c, lapply(Rb_mutant_wt_GLYR1actranks, function(X){unlist(X[["mutant"]])}))
names(all_Rb_GLYR1act_ranks) <- str_remove(names(all_Rb_GLYR1act_ranks), pattern = "^TCGA-[A-Z]{2,4}\\.")

all_wt_GLYR1act_ranks <- do.call(c, lapply(Rb_mutant_wt_GLYR1actranks, function(X){unlist(X[["wildtype"]])}))
names(all_wt_GLYR1act_ranks) <- str_remove(names(all_wt_GLYR1act_ranks), pattern = "^TCGA-[A-Z]{2,4}\\.")

GLYR1_NNMT_rb_tvalues <- lapply(1:1000, function(i){
  
  tryCatch({
    
    Rb_NNMTselect_vec <- unlist(lapply(RBmutant_select_cancers, function(X){
      
      sample(Rb_mutant_wt_NNMTranks[[X]][["mutant"]], 10)
      
    }))
    
    WT_NNMTselect_vec <- unlist(lapply(RBmutant_select_cancers, function(X){
      
      sample(Rb_mutant_wt_NNMTranks[[X]][["wildtype"]], 10)
      
    }))
    
    Rb_GLYR1select_vec <- all_Rb_ranks[match(names(Rb_NNMTselect_vec), names(all_Rb_ranks))]
    
    WT_GLYR1select_vec <- all_wt_ranks[match(names(WT_NNMTselect_vec), names(all_wt_ranks))]
    
    Rb_lm_df <- data.frame(NNMT = c(WT_NNMTselect_vec, Rb_NNMTselect_vec),
                           GLYR1 = c(WT_GLYR1select_vec, Rb_GLYR1select_vec))
    
    Rb_lm_df[row.names(Rb_lm_df) %in% names(WT_NNMTselect_vec), "Rbstatus"] <- "wildtype"
    Rb_lm_df[row.names(Rb_lm_df) %in% names(Rb_NNMTselect_vec), "Rbstatus"] <- "deleterious"
    
    NNMTfullmodel <- summary(lm(NNMT ~ GLYR1 + Rbstatus, data = Rb_lm_df))
    NNMTpartmodel <- summary(lm(NNMT ~ Rbstatus, data = Rb_lm_df))
    
    GLYR1fullmodel <- summary(lm(GLYR1 ~ NNMT + Rbstatus, data = Rb_lm_df))
    GLYR1partmodel <- summary(lm(GLYR1 ~ Rbstatus, data = Rb_lm_df))
    
    output_vec <- c(NNMTfullmodel$coefficients["GLYR1", "t value"],
                    NNMTfullmodel$coefficients["Rbstatuswildtype", "t value"],
                    GLYR1fullmodel$coefficients["Rbstatuswildtype", "t value"],
                    NNMTpartmodel$coefficients["Rbstatuswildtype", "t value"],
                    GLYR1partmodel$coefficients["Rbstatuswildtype", "t value"])
    
    names(output_vec) <- c("GLYR1_NNMT_t",
                           "NNMT_fullRb_t",
                           "GLYR1_fullRb_t",
                           "NNMT_partRb_t",
                           "GLYR1_partRb_t")
    
    return(output_vec)
    
  }, error = function(e){
    
    message("These samples produced a problem")
    
    print(row.names(Rb_lm_df))
    
  })
  
  
})

GLYR1_NNMT_rb_tvalues <- GLYR1_NNMT_rb_tvalues[sapply(GLYR1_NNMT_rb_tvalues, function(x){length(x) == 5})]
GLYR1_NNMT_rb_tvalues_df <- data.frame(do.call(rbind, GLYR1_NNMT_rb_tvalues))
GLYR1_NNMT_rb_tvalues_df <- GLYR1_NNMT_rb_tvalues_df[, 2:5]

GLYR1_NNMT_rb_tvalues_meltdf <- reshape2::melt(GLYR1_NNMT_rb_tvalues_df)

GLYR1_NNMT_rb_tvalues_meltdf$variable <- factor(GLYR1_NNMT_rb_tvalues_meltdf$variable, levels = c("GLYR1_partRb_t", "GLYR1_fullRb_t", "NNMT_partRb_t", "NNMT_fullRb_t"))

t.test(as.numeric(sapply(GLYR1_NNMT_rb_tvalues, function(x){x["NNMT_partRb_t"]})),
       as.numeric(sapply(GLYR1_NNMT_rb_tvalues, function(x){x["NNMT_fullRb_t"]})))$p.value

t.test(as.numeric(sapply(GLYR1_NNMT_rb_tvalues, function(x){x["GLYR1_fullRb_t"]})),
       as.numeric(sapply(GLYR1_NNMT_rb_tvalues, function(x){x["GLYR1_partRb_t"]})))$p.value

ticklabs <- c("RB1 mutation", "RB1 mutation + NNMT", "RB1 mutation", "RB1 mutation + GLYR1")

saveRDS(GLYR1_NNMT_rb_tvalues_meltdf,
        "plot_data/GLYR1_NNMT_rb_tvalues_meltdf.rds")
GLYR1_NNMT_rb_tvalues_meltdf <-readRDS("plot_data/GLYR1_NNMT_rb_tvalues_meltdf.rds")

write.table(GLYR1_NNMT_rb_tvalues_meltdf,
            file = "plot_data/Fig S13/Fig_S13E_plotdata.txt",
            row.names = FALSE,
            col.names = TRUE,
            sep = "\t")

# Fig S13E

pdf("graphics/TCGA_RBmut_GLYR1NNMTlinearmodel_tvals.pdf",
    width = 5,
    height = 3)

print(ggplot(GLYR1_NNMT_rb_tvalues_meltdf, aes(x = variable, y = -value)) + 
        geom_boxplot(fill = c("yellow",
                              "magenta",
                              "yellow",
                              "magenta")) + 
        theme_classic() + 
        geom_hline(yintercept = 0,
                   linetype = "dashed",
                   color = "gray") + 
        geom_vline(xintercept = 2.5,
                   color = "black") + 
        ylab(italic(RB1)~"mutation linear model t-value") + 
        xlab("Linear model terms")  +
        coord_cartesian(ylim = c(-5, 5)) +
        scale_x_discrete(labels= ticklabs) + 
        theme(axis.text.x = element_text(colour = "black",
                                         size = 8),
              axis.text.y = element_text(colour = "black",
                                         size = 10)) + 
        # annotate(x = 3.5,
        #          y = -3,
        #          label = italic(NNMT)~"expression",
        #          geom = "text",
        #          colour = "red",
        #          fontface = 2) +
        # annotate(x = 1.5,
        #          y = -3,
        #          label = "GLYR1 activity",
        #          geom = "text",
        #          colour = "red",
      #          fontface = 2) +
      annotate(x = 3.5,
               y = 4.85,
               label = substitute("p = 1.73 x"~10^-74),
               geom = "text",
               colour = "black") +
        annotate(x = 1.5,
                 y = 4.85,
                 label = "p = 0.315",
                 geom = "text",
                 colour = "black") 
      # geom_segment(aes(x = 1, xend = 2, y = 4.25, yend = 4.25)) +
      # geom_segment(aes(x = 3, xend = 4, y = 4.25, yend = 4.25))
)

dev.off()

# higher expression of GLYR1 in Rb mutant cancers?

GLYR1_ensembl <- "ENSG00000140632"

Rb_mutant_wt_GLYR1ranks <- lapply(cancers_with_rbs, function(thiscancer){
  
  thisdata <- TCGA_MOR_list[[thiscancer]]
  
  GLYR1forthisone <- thisdata[GLYR1_ensembl, ]
  GLYR1forthisone <- GLYR1forthisone[names(GLYR1forthisone) %in% TCGA_LUT$SAMPID]
  
  GLYR1_for_lm <- data.frame(GLYR1 = GLYR1forthisone)
  GLYR1_for_lm <- cbind(GLYR1_for_lm, TCGA_LUT[match(row.names(GLYR1_for_lm), TCGA_LUT$SAMPID), c("race", "gender", "tumour_stage", "days_to_birth", "sequencing_centre")])
  
  GLYR1_for_lm$days_to_birth <- as.numeric(GLYR1_for_lm$days_to_birth)
  
  GLYR1_for_lm <- GLYR1_for_lm[!apply(GLYR1_for_lm, 1, function(x){any(is.na(x))}), ]
  
  if(length(unique(GLYR1_for_lm$race)) == 1){
    
    racevariable <- NULL
    
  } else {
    
    racevariable <- "race"
    
  }
  
  if(length(unique(GLYR1_for_lm$tumour_stage)) == 1){
    
    tumourstagevariable <- NULL
    
  } else {
    
    tumourstagevariable <- "tumour_stage"
    
  }
  
  if(length(unique(GLYR1_for_lm$gender)) == 1) {
    
    gendervariable <- NULL
    
  } else {
    
    gendervariable <- "gender"
  }
  
  if(length(unique(GLYR1_for_lm$sequencing_centre)) == 1) {
    
    seqcentrevariable <- NULL
    
  } else {
    
    seqcentrevariable <- "sequencing_centre"
    
  }
  
  modelvariables <- c("days_to_birth", gendervariable, seqcentrevariable, tumourstagevariable, racevariable)
  
  f <- as.formula(
    paste("log10(GLYR1 + 1)",
          paste(modelvariables, collapse = " + "),
          sep = " ~ "))
  
  GLYR1_for_lm[, "GLYR1resid"] <- lm(formula = f, data = GLYR1_for_lm)$residuals
  GLYR1_for_lm[, "GLYR1resid_ranks"] <- rank(GLYR1_for_lm$GLYR1resid) / (nrow(GLYR1_for_lm) + 1)
  
  output_list <- list()
  
  output_list[["mutant"]] <- GLYR1_for_lm[str_remove(row.names(GLYR1_for_lm), pattern = "-[:alnum:]+-[:alnum:]+-[:digit:]+$") %in%  str_remove(Rbdeleterious_SAMPIDs_bycancer[[thiscancer]], pattern = "-[:alnum:]+-[:alnum:]+-[:digit:]+$"), "GLYR1resid_ranks"] 
  names(output_list[["mutant"]]) <- row.names(GLYR1_for_lm[str_remove(row.names(GLYR1_for_lm), pattern = "-[:alnum:]+-[:alnum:]+-[:digit:]+$") %in%  str_remove(Rbdeleterious_SAMPIDs_bycancer[[thiscancer]], pattern = "-[:alnum:]+-[:alnum:]+-[:digit:]+$"), ])
  
  output_list[["wildtype"]] <- GLYR1_for_lm[str_remove(row.names(GLYR1_for_lm), pattern = "-[:alnum:]+-[:alnum:]+-[:digit:]+$") %in%  str_remove(wtRB_SAMPIDs_bycancer[[thiscancer]], pattern = "-[:alnum:]+-[:alnum:]+-[:digit:]+$"), "GLYR1resid_ranks"] 
  names(output_list[["wildtype"]]) <- row.names(GLYR1_for_lm[str_remove(row.names(GLYR1_for_lm), pattern = "-[:alnum:]+-[:alnum:]+-[:digit:]+$") %in%  str_remove(wtRB_SAMPIDs_bycancer[[thiscancer]], pattern = "-[:alnum:]+-[:alnum:]+-[:digit:]+$"), ])
  
  return(output_list)
  
})

names(Rb_mutant_wt_GLYR1ranks) <- cancers_with_rbs

RBmutant_pancancer_GLYR1medians <- data.frame(t(sapply(1:1000, function(i){
  
  Rb_select_vec <- unlist(lapply(RBmutant_select_cancers, function(X){
    
    sample(Rb_mutant_wt_GLYR1ranks[[X]][["mutant"]], 10)
    
  }))
  
  WT_select_vec <- unlist(lapply(RBmutant_select_cancers, function(X){
    
    sample(Rb_mutant_wt_GLYR1ranks[[X]][["wildtype"]], 10)
    
  }))
  
  output_vec <- c()
  
  output_vec["mut_rank_median"] <- median(Rb_select_vec)
  output_vec["WT_rank_median"] <- median(WT_select_vec)
  
  output_vec["p_value"] <-   wilcox.test(Rb_select_vec, WT_select_vec)$p.value
  
  return(output_vec)
  
})))

colnames(RBmutant_pancancer_GLYR1medians)[1:2] <- c("Deleterious", "Wildtype")

# t.test p value is around 0
t.test(RBmutant_pancancer_GLYR1medians$Deleterious, RBmutant_pancancer_GLYR1medians$Wildtype)$p.value

RBmutant_pancancer_GLYR1medians_melt_forplot <- reshape2::melt(RBmutant_pancancer_GLYR1medians[, c("Wildtype", "Deleterious")])

saveRDS(RBmutant_pancancer_GLYR1medians_melt_forplot,
        "plot_data/RBmutant_pancancer_GLYR1medians_melt_forplot.rds")
RBmutant_pancancer_GLYR1medians_melt_forplot <- readRDS("plot_data/RBmutant_pancancer_GLYR1medians_melt_forplot.rds")

write.table(RBmutant_pancancer_GLYR1medians_melt_forplot,
            file = "plot_data/Fig S13/Fig_S13C_plotdata.txt",
            row.names = FALSE,
            col.names = TRUE,
            sep = "\t")

# Fig S13C

pdf("graphics/GLYR1medianranks_RB1mut.pdf",
    height = 2.5,
    width = 2.5)

ggplot(RBmutant_pancancer_GLYR1medians_melt_forplot, aes(x = variable, y = value)) + 
  geom_boxplot(fill = c("gray", "red")) +
  theme_classic() + 
  ylab("Median"~italic(GLYR1)~"expression %ile") + 
  xlab(italic(RB1)~"mutation status") + 
  theme(axis.text.y = element_text(colour = "black",
                                   size = 10),
        axis.text.x = element_text(colour = "black",
                                   size = 8),
        axis.title.y = element_text(hjust = 0,
                                    size = 10)) + 
  coord_cartesian(ylim = c(0.3, 0.8)) + 
  geom_hline(yintercept = 0.5,
             linetype = "dashed",
             col = "gray") + 
  annotate(geom = "text",
           x = 1.5,
           y = 0.8,
           label = substitute("p = 9.76 x"~10^-126))

dev.off()