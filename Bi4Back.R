#################################### Genetic expression Binarization ################################################

####################################      Ismail Belgacem 2024       ###############################################   

#################################### The code was produced during this year in 2024  #################################################




#  How to install igraph using "BiocManager"
#  if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("igraph")

require(igraph)
library(igraph)

# Read_Data_File  is a function that Read any CSV file, TSV file or an other type (That contain a header).

rm(Read_Data_File)
Read_Data_File <- function(f = file){
  
  name_data_file <- basename(f)
  
  if(substr(name_data_file, nchar(name_data_file) - 2, nchar(name_data_file)) == "csv"){
    data_dataFrame <- read.csv(f2, sep = ";" )
  }
  else{
    data_dataFrame <- read.delim(f, header=TRUE)
  }
}

# InteracDataFactor is a fuction that:
# fist Read  without a "header"   the the Interaction file (that represent the regulatory graph)  and then
# convert the regulators names in the the interaction file by the corresponding factor in the disc, and 
# keeping only the gene that have their regulators names in the disc, and finally 
# convert the corrected interaction file into a data frame to be treated by using "igraph" package. 

rm(InteracDataFactor)

InteracDataFactor <- function(f = file){
  name_data_file <- basename(f)
  
  if(substr(name_data_file, nchar(name_data_file) - 2, nchar(name_data_file)) == "csv"){
    Dframe <- read.csv(f2, sep = ";" )
  }
  else{
    Dframe <- read.delim(f, header=FALSE)
  }

# Capitalize everything, change the case to simplify the analysis or the detection of factors:
  Dataframe<- data.frame(lapply(Dframe, toupper))
  
# Disc with words that mean an inhibition
  inhibit_disc <- c("INHIBIT", "INHIB", "IN", "DOWN", "DOWNREG",-1)
# Disc with words that mean an activation
  activate_disc <- c("ACTIVATE", "ACTIV", "AC", "UP", "UPREG",1)
  
  Dataframe<-Dataframe[(Dataframe[,2] %in% activate_disc) | (Dataframe[,2] %in% inhibit_disc),]
  
# If none are founded, error
  
  if(nrow(Dataframe)==0){
    tryCatch(
      error = function(cnd) {
        print("Convertion cannot be proceed")
      })
  }
  
  Dataframe[Dataframe[,2] %in% activate_disc,2]<- 1
  Dataframe[Dataframe[,2] %in% inhibit_disc,2]<- -1
  
# Convert the corrected interaction file into a data frame to be treated using "igraph" package.
  
  interaction_data_frame<-data.frame(regulator_gene = Dataframe$V1,
                 regulated_gene = Dataframe$V3,
                 type = Dataframe$V2)
  return(interaction_data_frame)
}

# Regulators is a function that provide all the regulators of a given target using the interaction file.

rm(Regulators)

Regulators <- function(file, gene = gene, regulator_type=""){
  
# Read & Convert the interaction file into a data frame
  
  interaction_data_frame= InteracDataFactor(file)
  
  # If activators are needed
  if(regulator_type == 1){
    interaction_data_frame=interaction_data_frame[interaction_data_frame$type==1,]
    G=graph.data.frame(interaction_data_frame)
  }
  # If inhibitors are needed
  if(regulator_type == -1){
    interaction_data_frame=interaction_data_frame[interaction_data_frame$type==-1,]
    G=graph.data.frame(interaction_data_frame)
  }
  # If the type of regulators is not indicated or and all the regulators are needed
  else{
    G=graph.data.frame(interaction_data_frame)
  }
  
  vertex_var=V(G)$name
  
 # returns the regulators needed
  
  if(gene %in% vertex_var) {
    reg_var=neighbors(G,gene,"in")
    if(length(reg_var) >= 1) {
      return(neighbors(G,gene,"in")$name)
    }
  }
}



# InConsistentRegulators is a function that return the inconsistent regulators of a binarized gene
# using the current binary profile and the interaction data.

rm(InConsistentRegulators)

InConsistentRegulators <- function(f1 = interac_data, data_Frame = exp_data, gene = backward_gene, bg = gene_bin_value){

# Collect the positive and negative regulators
  
data_Frame <- data_Frame[!is.na(data_Frame$Boolean_Profil), ]
  positive_regulator = Regulators(f1, gene, 1)
  negative_regulator = Regulators(f1, gene, -1)
  

# Get the inconsistent regulators
  inconsistent_regulators <- union( data_Frame$genelist[data_Frame$genelist %in% positive_regulator & data_Frame$Boolean_Profil != bg] ,
                                    data_Frame$genelist[data_Frame$genelist %in% negative_regulator & data_Frame$Boolean_Profil == bg])
  
# Returns them
  
  return(inconsistent_regulators)
  
}

# Get_Taux_Values is a function that compute taus values of the consistent regulators of a given gene target (that has a binary value)
# Using interaction data and the expression data

rm(Get_Taux_Values)
Get_Taux_Values <- function(f1 = interac_data, data_Frame = exp_data, gene, Consistent_Reg) {
 
# Initialization of  taus values of the consistent regulators of the target
  
  dataFrame_t_v <- data.frame(CstRegulators= Consistent_Reg, taux_values= NA)
  
# Computing  taus values 
for (reg_var in Consistent_Reg) {
  
  if(reg_var %in% Regulators(f1, gene,1)){
    
    dataFrame_t_v$taux_values[dataFrame_t_v$CstRegulators==reg_var]= ((as.numeric(data_Frame$Boolean_Profil[data_Frame$genelist==gene])*
          (1 - data_Frame$Norm_gene_expression[data_Frame$genelist==reg_var])) 
    + (data_Frame$Norm_gene_expression[data_Frame$genelist==reg_var]*
         (1- as.numeric(data_Frame$Boolean_Profil[data_Frame$genelist==gene]))))
    
  }
  
  if(reg_var %in% Regulators(f1, gene,-1)){
    
    dataFrame_t_v$taux_values[dataFrame_t_v$CstRegulators==reg_var]= ( (as.numeric(data_Frame$Boolean_Profil[data_Frame$genelist==gene]) * 
    data_Frame$Norm_gene_expression[data_Frame$genelist==reg_var]) + ((1 - data_Frame$Norm_gene_expression[data_Frame$genelist==reg_var]) *
      (1- as.numeric(data_Frame$Boolean_Profil[data_Frame$genelist==gene]))) )

  }
  
}
# return them
return(dataFrame_t_v)

}


# Get_Value_Back_Propagation is a function that provide the binary values of regulators in function of the binary value of the target
rm(Get_Value_Back_Propagation) 
Get_Value_Back_Propagation <- function(f1 = interac_data, gene, gene_bin_value, reg_var){
  if(reg_var %in% Regulators(f1, gene,1)){
    return(gene_bin_value)
  }
  else{
    return(!gene_bin_value)
  }
}


# A Normalization function based on min-max methods

Normalize <- function(x) {
  return((x - min(x,na.rm = TRUE)) /(max(x, na.rm = TRUE) - min(x,na.rm = TRUE)))
}


# Bi4Back is the main function of binarisation that  compute of the Boolean profile of gene expression experiments 
# using: the interaction file, the gene expression data file, the gene_name column and the considered experiment column in the gene expression data file

rm(Bi4Back) 

Bi4Back <- function(f1 = Interaction_file, f2 = DATA_file, GeneName_column ,
                         exp_column , delta_real = 0.1, epsilon_real = 0.05){
  
 # Data files treatment where we:
 # Create first  a data frame where 
 # The first column contains all the gene names in the interaction data file (the graph)
 # The second column contains all the corresponding gene expression values in the expression data file (if the measurement of a gene does not exist we assign NA) 
 # The third column contains the normalized gene expression data.
  
  interaction_data_frame <- InteracDataFactor(f1) # Read the interaction data file and convert it to a data frame.
  G3 <-graph.data.frame(interaction_data_frame)
  GeneName_columns <- V(G3)$name
  genelist <- V(G3)$name
  dataFrame_GeneName_columns <- data.frame(genelist=GeneName_columns)
  
  exp_data <- Read_Data_File(f2) # Read the gene Expression Data file
  genelist_in_exp_data  <- genelist[genelist %in% exp_data[,GeneName_column]]
  genelist_not_exp_data <- genelist[!genelist %in% exp_data[,GeneName_column]]
  data_gene_expression  <- exp_data[exp_data[,GeneName_column] %in% genelist_in_exp_data, c(GeneName_column, exp_column)]
  data_Frame <- data.frame(genelist=data_gene_expression[,1], gene_expression=data_gene_expression[,2])
  dataFrame_gene_not<- data.frame(genelist=genelist_not_exp_data, gene_expression=NaN)
  data_Frame <- rbind(data_Frame, dataFrame_gene_not)
  
  # Verify if the data are already normalized, if not do it
  
  if(!all(data_Frame$gene_expression <= 1, na.rm = TRUE) | !all(data_Frame$gene_expression >= 0, na.rm = TRUE)){
    norm_data=Normalize(data_Frame$gene_expression)
  }
  else{
    norm_data <- data_Frame$gene_expression
  }
  
  data_Frame <- cbind(data_Frame, Norm_gene_expression=norm_data)
  
# If some expression data are missing, corrects them by assigning them with neutral values equals to 0.5.
  
  data_Frame$Norm_gene_expression[is.na(data_Frame$Norm_gene_expression)] <- 0.5

  # The initialization of the binary proﬁle to NA and create the forth  column in the previous data frame, that represent Boolean Profil   
  
  data_Frame <- cbind(data_Frame, Boolean_Profil= NA)
  
  # A binary proﬁle is generated with the deﬁnition of the levels of some genes
  # using only extreme cases of gene expression data values (with respect a chosen a Delta value) 
  data_Frame$Boolean_Profil[data_Frame$Norm_gene_expression <= delta_real] <- FALSE
  data_Frame$Boolean_Profil[data_Frame$Norm_gene_expression >= 1 - delta_real] <- TRUE
  
  # Initialization of our condition of loop
  
  existupdate <- TRUE
  
  while (existupdate){
    
    # Stock the previous binary profile
    previous_Boolean_profil <- data_Frame$Boolean_Profil
    
    ############################################ FORWARD PROPAGATION #####################################################
    
    # using the regulatory graph and the underlying Boolean network rules, a completion of NA genes expression levels
    # by forwarding the Boolean values of regulators toward the target.
    
    undefinedvars  <- data_Frame$genelist[is.na(data_Frame$Boolean_Profil)] 
    
    # Loop for each gene without a binarized profil
    
    for(gene in undefinedvars){

      # get its regulators
      positive_regulator = Regulators(f1, gene, 1)
      negative_regulator = Regulators(f1, gene, -1)
    
      # If all the regulators have a binary value, process
      
      if(length(Regulators(f1, gene)) >= 1 && all(!is.na(data_Frame$Boolean_Profil[data_Frame$genelist %in% Regulators(f1, gene)])))
     {
        # If all the activators are on and all the inhibitors are off, affectation a positive binary value
        
        if (all(data_Frame$Boolean_Profil[data_Frame$genelist %in% positive_regulator] == TRUE) &&
            all(data_Frame$Boolean_Profil[data_Frame$genelist %in% negative_regulator] == FALSE)) 
          {
          data_Frame$Boolean_Profil[data_Frame$genelist == gene] <- TRUE
        }
        
        # Vice-versa
        if (all(data_Frame$Boolean_Profil[data_Frame$genelist %in% positive_regulator] == FALSE) &&
            all(data_Frame$Boolean_Profil[data_Frame$genelist %in% negative_regulator] == TRUE)) 
        {
          data_Frame$Boolean_Profil[data_Frame$genelist == gene] <- FALSE
        }
      }
      
    }
    
   ################################################### BACKWARD PROPAGATION  #############################################
    
   #  Using the regulatory graph and the operating rules of Boolean networks, a completion of the not assigned (NA) genes expression levels
   #  by a back propagation of the Boolean values of target toward its regulators
   #  under the condition the target is deﬁned (1 or 0) and all its consistent regulators are NA.
  
   # The already binarized variable
    
    alreadybinarizedvars <- data_Frame$genelist[!is.na(data_Frame$Boolean_Profil)]
    
    # A backward propagation  of the genes with a binary value
    
    for(gene in alreadybinarizedvars) {
      
      if (length(Regulators(f1, gene)) >= 1) {
        ConsistentRegulators = setdiff(Regulators(f1, gene), InConsistentRegulators(f1, data_Frame, gene, data_Frame$Boolean_Profil[data_Frame$genelist==gene]))
       
        
        
        ############################################## INCONSISTENCY TEST #######################################################################################
        
        # If none are found, a re-initialisation of the binary values  of the target (After the gene expression levels are assigned, to avoid confusions and the propagation of inconsistencies an evaluation
        #  using the regulatory graph and the current Boolean values should be performed here)
        
        if (length(ConsistentRegulators) == 0) {
          data_Frame$Boolean_Profil[data_Frame$genelist==gene] <- NA
        }
        
        # If at least one consistent regulators is found and are NA, do Backward propagation
        
        if ((length(ConsistentRegulators) >= 1) && all( is.na (data_Frame$Boolean_Profil[data_Frame$genelist %in% ConsistentRegulators ]) ) ) 
          {
          
          Taux_Consistent_Regulators <- Get_Taux_Values(f1, data_Frame, gene, ConsistentRegulators)
          
          
          chosenregs <- Taux_Consistent_Regulators$CstRegulators[Taux_Consistent_Regulators$taux_values==min(Taux_Consistent_Regulators$taux_values)]
          
         
          
          # Attribution of a binary value of regulators based on taus values
          
          for (var in chosenregs) {
            data_Frame$Boolean_Profil[data_Frame$genelist==var]=Get_Value_Back_Propagation(f1, gene, data_Frame$Boolean_Profil[data_Frame$genelist==gene], var)
          } 
          
          
          #################################################### Harmonisation #############################################################
         
           # Also assign a Boolean value to the regulators that have almost a similar gene expression value compared to regulator that is already assigned using a back propagating.
          
          creg = setdiff(ConsistentRegulators,chosenregs)
          
          reg=Taux_Consistent_Regulators$CstRegulators[Taux_Consistent_Regulators$CstRegulators %in% creg  &  abs(Taux_Consistent_Regulators$taux_values-min(Taux_Consistent_Regulators$taux_values)) <= epsilon_real]
          
          
          for (creg_var in reg) {

            data_Frame$Boolean_Profil[data_Frame$genelist==creg_var]= Get_Value_Back_Propagation(f1, gene, data_Frame$Boolean_Profil[data_Frame$genelist==gene], reg)
            
              }
          
          
          
         }
        
        
    }
    
    }
    
    # Loop condition
    # If a modification was made, loop it
    
    if(identical(previous_Boolean_profil, data_Frame$Boolean_Profil)){
      existupdate <- FALSE
    }
    
  }

  return(data_Frame)
}

#################################################### PROGRAM TEST ###############################################################

Bi4Back("IG.sif","70510d37-dc1a-433b-84ee-d0f7d9ca1c3d.rna_seq.augmented_star_gene_counts.tsv", 2,  5)
