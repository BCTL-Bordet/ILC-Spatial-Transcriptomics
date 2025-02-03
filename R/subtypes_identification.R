#### Assigning ILC TME-Based Subtypes on TPM-normalised RNA-seq HR+/HER2- ILC dataset ####

#### Function to Compute Signature Scores ####
# This function calculates gene signature scores based on predefined gene expression signatures.
# It supports both raw count data and normalized expression data, applying different statistical methods accordingly.

calcSig <- function(d, sig, dropEmpty=TRUE, balanced=FALSE, isCount=FALSE, useNB=FALSE,
                    updateWeights=FALSE, returnW=FALSE, maxN=200)
{ 
  # Check if input is a sparse matrix and set appropriate column operations
  if (is(d, "Matrix")) { colMeans=Matrix::colMeans; colSums=Matrix::colSums; }
  
  # If multiple signatures are provided as a list, recursively apply calcSig to each
  if (is.list(sig) && !is.data.frame(sig))
  { ret = sapply(sig, function(i) calcSig(d, i, balanced=balanced, isCount=isCount, useNB=useNB, updateWeights=updateWeights));
  
    # Drop empty results if required
    if (dropEmpty) { ret = ret[,!colSums(!is.na(ret))==0]; }
    return(ret);
  }
  
  # Extract gene names from the signature file
  s = as.character(sig[,1]);
  s[s==""] = NA;
  x = intersect(rownames(d), s);
  x = x[!is.na(x)];
  
  # Handle cases where no genes in the signature match the dataset
  if (length(x) == 0)
  { if (colnames(sig)[1] == "name") { other = "entrez"; } else { other = "name"; }
    s = as.character(sig[,other]);
    s[s==""] = NA;
    x = intersect(rownames(d), s);
    x = x[!is.na(x)];
    if (length(x) == 0)
    { warning("No gene in common");
      return(rep(NA, ncol(d)));
    }
  }
  
  # Match genes in signature with dataset
  sig = sig[match(x, s),];
  
  # If input is count data, apply negative binomial model
  if (isCount) 
  { 
    if (useNB)
    { 
      if (length(x)==1) { return(log10(1e4*d[x,]/colSums(d)+1)); }
      w = rowSums(d[x,])>1; x = x[w]; sig = sig[w,]; 
      
      # Load required MASS library for negative binomial model
      if(!require(MASS)) { stop("MASS library not installed"); }
      if (any(duplicated(colnames(d)))) { stop("No duplicated colnames"); }
      
      # Prepare dataset for modeling
      y = cbind(as.vector(d[x,]), rep(colSums(d), each=length(x)));
      id = factor(rep(colnames(d), each=length(x)))
      gid = factor(rep(x, ncol(d)));
      a = sig$coefficient;
      Ns = colnames(d);
      
      # Limit computation size if needed
      if (ncol(d)>maxN)
      { Ns = colnames(d)[sample(1:ncol(d), ifelse(length(x)>maxN, maxN/2, maxN))];
        su = id %in% Ns;
        y = y[su,]; id = id[su]; gid = gid[su];
      }
      
      # Fit negative binomial model
      if (all(sig$coefficient==1))
      { fm = glm.nb.noErr(y[,1] ~ id + gid -1 + offset(log(y[,2])));
        co = coef(fm);
        co2 = co[grep("^id", names(co))];
        names(co2) = sub("^id", "", names(co2))
        co2 = co2[Ns]
        
        # Update weights if required
        if (updateWeights)
        { a = tapply(1:nrow(y), gid, function(i) { z=y[i,];
        coef(glm.nb.noErr(z[,1] ~ co2 + offset(log(z[,2]))))["co2"]}) 
        }
      } else { updateWeights=updateWeights+1; }
      
      # Further iterations for weight updating
      for (iter in 1:updateWeights)      
      { g2b = rep(a, length(id)/length(a))
        fm = glm.nb.noErr(y[,1] ~ g2b*id-g2b-id-1 + gid + offset(log(y[,2])));
        co = coef(fm);
        co2 = co[grep("^g2b:id", names(co))]; co2[is.na(co2)] = 0;
        names(co2) = sub("^(g2b)*:id", "", names(co2));
        co2 = co2[Ns]; 
        if (iter<updateWeights)
        {  a = tapply(1:nrow(y), gid, function(i) { z=y[i,];
        coef(glm.nb.noErr(glm.nb(z[,1] ~ co2 + offset(log(z[,2])))))["co2"]} );
        }
      }
      
      # Return weights if required
      if (returnW) { names(a) = x; return(list(sig=co2, weights=a)); }
      return(co2)
    }
    
    # Function to calculate relative expression scores
    fCalc = function(d, x, coef)
    { r = colSums(d[x,,drop=FALSE]*coef, na.rm=TRUE)/colSums(d, na.rm=TRUE);
      sign(r)*log10(1+abs(r)*100);
    }
  } else { fCalc = function(d, x, coef) { colMeans(d[x,,drop=FALSE]*coef); } }
  
  # Handle balanced signature scoring
  if (balanced && any(!(w <- (sig[,"coefficient"] > 0))))
  { 
    val1 = fCalc(d, x[w], sig[w,"coefficient"]);
    val2 = fCalc(d, x[!w], sig[!w,"coefficient"]);
    val = (val1+val2)/2;
  }
  else
  { val = fCalc(d, x, sig[,"coefficient"]);
  }
  return(val);
}

#### Function to Compute Molecular Subtypes ####
# "expression_data" is TPM RNA-seq data (gene names as rownames, sample names al colnames)
# "signature_file" is the list "signatures_subtypes.RDS" containing the subtypes-related gene signatures
library(genefu)
compute_molecular_subtypes <- function(expression_data, signature_file) {
  
  # Load predefined gene signature groups
  sig_groups <- readRDS(signature_file)
  
  # Compute signature scores for each sample
  signature_scores <- calcSig(expression_data, sig_groups)
  
  # Convert results to a dataframe
  signature_scores_df <- as.data.frame(signature_scores)
  
  # Define molecular subtypes to normalize
  subtypes <- c("NSE", "P", "ARE", "MIE")
  
  # Normalize signature scores using rescale function
  normalized_scores <- as.data.frame(
    apply(signature_scores_df[, subtypes], 2, genefu::rescale, q = 0.05)
  )
  
  # Assign each sample to the most dominant subtype
  normalized_scores$group <- colnames(normalized_scores)[max.col(normalized_scores)]
  
  return(normalized_scores)
}

#### Example Usage ####
# expression_data <- readRDS("expression_lobular.RDS")
# signature_file <- readRDS("signatures_subtypes.RDS")
# subtypes_result <- compute_molecular_subtypes(expression_data, signature_file)
# write.csv(subtypes_result, "molecular_subtypes.csv")
