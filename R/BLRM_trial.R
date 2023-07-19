library(R2jags)

source("BLRM.bugs")

# Function to calculate the probabilities within different intervals based on MCMC samples
interval_prob <- function(jagsdata, Pint){
  interval_levels <- 1:(length(Pint)-1)
  
  tab <- lapply(2:length(Pint), function(x) jagsdata >= Pint[x])
  total_matrix <- 1 + Reduce("+", tab)
  
  # Crucial step to convert all columns to factors, in case some categories are not present at a certain dose level
  out_fct <- lapply(data.frame(total_matrix), factor, levels = interval_levels)
  tt <- data.frame(do.call(cbind, lapply(out_fct, function(x) {prop.table(table(x, useNA = "no"))})))
  probdt = data.frame(t(tt))
  colnames(probdt) <- c("Punder", "Ptarget", "Pover")
  rownames(probdt) <- NULL
  
  return(probdt)
}

# Function to check stopping criteria for the BLRM design
checkstop_TI <- function(probdt, target.prob = 0.5, min.subj.MTD = 9, max.subj = 40, ewoc = 0.25) {
  probdt$Toxic <- ifelse(probdt$Pover >= ewoc, 1, 0)
  Ntotal <- sum(probdt$Npat)
  currdt <- probdt[probdt$Current == 1, ]
  
  max_target <- max(probdt$Ptarget[probdt$Npat > 0])
  cond1 <- currdt$Ptarget >= target.prob & !currdt$Toxic & currdt$Ptarget == max_target
  cond2 <- currdt$Npat >= min.subj.MTD
  cond3 <- Ntotal >= max.subj
  cond4 <- prod(probdt$Toxic)
  
  stop4mtd <- cond1 & cond2
  stop4tox <- cond4
  stop4ss <- !stop4mtd & cond3
  nostop <- !(stop4mtd | stop4ss | stop4tox)
  stopcode <- 0 * nostop + 1 * stop4ss + 2 * stop4mtd + 3 * stop4tox
  
  return(stopcode)
}

# Function to determine the next action in the BLRM design
action_BLRM <- function(probdt, ewoc = 0.25){
  next_dose <- probdt$Dose[probdt$Pover < ewoc & probdt$Ptarget == max(probdt$Ptarget)]
  if (length(next_dose) > 1) next_dose <- next_dose[length(next_dose)]
  curr_dose <- probdt$Dose[probdt$Current == 1]
  action <- -1 * (next_dose < curr_dose) + 0 * (next_dose == curr_dose) + 1 * (next_dose > curr_dose)
  
  return(action)
}

# Main function that performs the BLRM design
BLRM_design = function(DoseProv = c(10, 25, 50, 100, 200, 400, 800), DoseRef = 100,
                       Tox_prob = c(0.02, 0.07, 0.15, 0.29, 0.48, 0.68, 0.83),
                       init_dose = 1,
                       cohort_size = 3, Nmax = 45, 
                       Pint_BLRM = c(0, 0.16, 0.33, 1), 
                       ewoc = 0.5, target_prob = 0.5,
                       prior_ab = c(-0.693, 0, 2, 1, 0),
                       seed = 2087102){
  
  set.seed(seed)
  
  #toxdt <- data.frame(Dose = DoseProv, Rate = Tox_prob)
  Dose_curr <- DoseProv[init_dose]
  cumdt <- data.frame(Dose = DoseProv, Rate = Tox_prob, Ntox = 0, Npat = 0, Current = 0)
  stopcode <- 0
  cohortdt_s <- NULL 
  cohort_index <- 1
  MTD <- 0
  
  while(stopcode==0){
    # find the current dose
    #toxrate <- subset(toxdt, Dose == Dose_curr)
    
    # Generate toxicity events for the cohort
    Ntox_new <- rbinom(n = 1, size = cohort_size, prob = Tox_prob[DoseProv == Dose_curr])
    
    # Update cohort data frames
    cohortdt_s <- rbind(cohortdt_s, data.frame(Dose = Dose_curr, Ntox = Ntox_new, Npat = cohort_size, cohort = cohort_index))
    cumdt$Ntox <- cumdt$Ntox + ifelse(cumdt$Dose %in% Dose_curr, Ntox_new, 0)
    cumdt$Npat <- cumdt$Npat + ifelse(cumdt$Dose %in% Dose_curr, cohort_size, 0)
    cumdt$Current <- as.integer(cumdt$Dose == Dose_curr)
    
    admdt <- subset(cumdt, Npat > 0)
    
    # Prepare data for JAGS model
    jags_data <- list(
      Ntox = admdt$Ntox,
      Npat = admdt$Npat,
      DoseAdm = admdt$Dose,
      NdoseAdm = nrow(admdt),
      DoseProv = cumdt$Dose,
      NdoseProv = nrow(cumdt),
      DoseRef = DoseRef,
      Prior = prior_ab
    )
    
    # Compile and run the JAGS model
    jags_obj <- jags(model.file = BLRM_orig, data = jags_data, 
                     parameters.to.save = c("Pr_Tox"), 
                     n.chains = 3, n.burnin = 10000, n.iter = 50000, 
                     progress.bar = "none", quiet = T)
    
    # Extract the MCMC samples for the parameter of interest (Pr_Tox in this case)
    PrTox_mcmc <- data.frame(jags_obj$BUGSoutput$sims.matrix[,seq.int(length(DoseProv))])
    
    # Calculate the probability within different intervals
    BLRM_prob <- cbind(cumdt, interval_prob(jagsdata = PrTox_mcmc, Pint = Pint_BLRM))
    BLRM_prob <- data.frame(BLRM_prob)        
    
    stopcode <- checkstop_TI(BLRM_prob, target.prob = target_prob, max.subj = Nmax, ewoc = ewoc)
    
    cohort_index <- cohort_index + 1
    
    #continue to next cohort
    if(stopcode == 0){
      
      action <- action_BLRM(BLRM_prob, ewoc)
      
      Dose_next <- DoseProv[which(DoseProv == Dose_curr) + action]
      cumdt$Current <- as.integer(cumdt$Dose==Dose_next)
      Dose_curr <- as.numeric(cumdt[cumdt$Current == 1, "Dose"])
    } 
    #stop for reaching maximum sample size
    if(stopcode == 1){
      MTD <- 999999
    }
    #stop for declaring MTD
    if(stopcode == 2){
      MTD <- as.numeric(subset(cumdt, Current == 1)$Dose)
    }
    #stop because all doses are toxic
    if(stopcode == 3){
      MTD <- -1
    }
    
  } # End while
  
  return(list(MTD = MTD,
              toxdt = cumdt))
}

# Call the main function to perform the BLRM design
BLRM_design()
