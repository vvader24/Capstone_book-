# Define the function for Bayesian Parallel Analysis --------------------

bayes.pa <- function(data, n.thin = 1, N.draw = 100, n.burnin = 1000, N.CD = 100, quantile = .95, report = "probability"){

  # data refers to the raw data of scores
  # n.thin refers to the thinning parameter for MCMC
  # N.draw refers to the number of draws from the posterior to be used in B-PA
  # n.burnin refers to the number of iterations to burn in MCMC sampler
  # N.CD refers to the number of comparison datasets (per each draw from the posterior)
  # quantile refers to the cutoff of the distribution used to evaluate whether the observed eigenvalue is sufficiently large
  # report refers to whether to report results as a probability or percentage

  # load the required packages
  library(psych)
  library(MCMCpack)



  #function for generating the model syntax to run lavaan
  triag_model <- function(P, num.factor){

    random.rows <- sample(1:P, num.factor, replace=F) #generate random numbers to decide which rows to use

    output_fixed <- matrix(rep(NA,(num.factor)*(num.factor-1)),
                           nrow = (num.factor)*(num.factor-1)/2, ncol = 2)

    output_constraint <- matrix(rep(NA,num.factor*2), ncol = 2)

    count2 <- 1
    for (i1 in 1:(num.factor-1)){
      for (i2 in (i1+1):num.factor){
        output_fixed[count2,1] <- i2
        output_fixed[count2,2] <- random.rows[i1]
        count2 <- count2 + 1
      }
    }

    for (i1 in 1: num.factor){
      output_constraint[i1,1] <- i1
      output_constraint[i1,2] <- random.rows[i1]
    }

    colnames(output_fixed)      <- c("factor","item")
    colnames(output_constraint) <- c("factor","item")


    #specifying the loading matrix
    row_which.factor <- c(rep(NA,num.factor))
    tmp1 <- matrix(c(rep(NA,P)),nrow=num.factor,ncol=P)
    for (which.factor in 1:num.factor){
      for (which.item in 1:P){
        tmp1[which.factor, which.item] <- paste("r",which.factor,which.item, "*V", which.item, sep="")
      }

      row_which.factor[which.factor] <- paste("f", which.factor,"=~", tmp1[which.factor, 1]  ,sep="")
      for (which.item in 2:P){
        row_which.factor[which.factor] <- paste(row_which.factor[which.factor], tmp1[which.factor, which.item], sep="+")
      }
    }

    #facor correlations need to be fixed at 0;
    #given num.factor, I need num.factor-1 rows to make it an identity matrix
    row_factor_cor <- c(rep(NA, (num.factor*(num.factor-1)/2)))
    for (i in 2:num.factor){
      for (j in 1: (i-1))
        row_factor_cor[ (i-1)*(i-2)/2+j ] <- paste("f", i, "~~0*f", j, sep="")
    }

    #next, put the loadings to zero indicated by the randomly selected rows;
    rows_constraint <- rep(NA, (num.factor)*(num.factor-1)/2)

    count<-1

    for (i in 1: (length(random.rows)-1)){
      for (j in (i+1):(num.factor)){
        rows_constraint[count] <- paste("r", j ,random.rows[i],"==0",sep="")
        count<-count+1
      }
    }

    #Last, combine the syntax to form the model specification
    model <- row_which.factor[1]
    for (i in 2:(num.factor)){
      model <- paste(model, row_which.factor[i],sep='\n')
    }


    for (i in 1:length(row_factor_cor)){
      model <- paste(model, row_factor_cor[i],sep='\n')
    }

    for (i in 1:length(rows_constraint)){
      model <- paste(model, rows_constraint[i], sep='\n')
    }
    # print the syntax I wrote using : cat(model)
    #random.rows


    fixed_param <- rep(NA, num.factor)

    for (i in 1:num.factor){

      fixed_param[i] <- paste("r", i, random.rows[i], sep="")

    }

    return(list(model=model,
                fixed_param=fixed_param,
                output_fixed = output_fixed,
                output_constraint = output_constraint))
  }
  #function for generating parallel data set
  gen.pa.data <- function(n_sample, ni, est_lambda_CD, est_unique_CD){

    implied_matrix <- est_lambda_CD %*% t(est_lambda_CD) + diag(est_unique_CD)
    generated.data <- mvrnorm(n = n_sample, mu = rep(0, ni), Sigma = implied_matrix)

  }

  N <- dim(data)[1] #sample size
  P <- dim(data)[2] #n of item

  data <- as.matrix(data)

  c_names <- rep(NA, P)
  for (i in 1:P){c_names[i]<-paste("V",i,sep="")}
  colnames(data) <- c_names


  eigs.data <- fa(data, nfactors=1, rotate="none", SMC=TRUE, max.iter=1, fm="pa", warnings = FALSE)$values # Compute the eigenvalues for observed data
  F.Max <- sum(eigs.data>0) #define the max number of possible factors

  eigs.comp.array <- array(NA, c(P, N.CD, N.draw, F.Max)) # set up array for saving eigenvalues from comparsion datasets

  #create a vector to save the results
  sav_result <- matrix(rep(NA, N.draw * P), nrow = N.draw, P)


  # Generate N.CD comparision datasets based on the kth factor and compare eigen values, test the number of factors
  T.NFactor=1 # Start with testing one factor

  while(T.NFactor <= F.Max){ # Loop from one factor to the max number of factors
    # Testing the first factor
    if(T.NFactor==1){

      print(paste0("Evaluating factor ", T.NFactor, " of ", F.Max))

      for(which.cd in 1:N.CD){ # Loop over N.CD comparsion datasets

        # generate one comparsion dataset using computed loadings
        comp.data <- matrix(rnorm(N*P), nrow=N, ncol=P)
        cor.comp.data <- cor(comp.data) # compute the correlation matrix

        # save the eigenvalues from running PAF on the comparison dataset
        eigs.comp.array[ , which.cd, 1, T.NFactor] <- fa(cor.comp.data, nfactors=1, rotate="none", SMC=TRUE, max.iter=1, fm="pa", warnings = FALSE)$values
      }

      if(eigs.data[T.NFactor]>as.numeric(quantile(eigs.comp.array[T.NFactor, , 1,  T.NFactor], quantile))) T.NFactor <- T.NFactor+1 else break
    }

    # Testing the second factor
    if(T.NFactor==2){

      print(paste0("Evaluating factor ", T.NFactor, " of ", F.Max))

      #step 1, fit the 1 factor CFA model to bayes
      list.of.lambda.constraints <- list(
        V1=list(1,"+") #constrain the first loading to be positive
      )

      posterior.draws <- MCMCfactanal(
        x = data,
        factors = T.NFactor-1,
        lambda.constraints = list.of.lambda.constraints,
        burnin = n.burnin,
        mcmc = N.draw*n.thin,
        thin = n.thin,
        store.scores = FALSE,
        std.var = TRUE
      )

      #step 2, generate CD based on posterior.draws
      #save the estimated loading parameters
      est_lambda <- array(1, c(P, T.NFactor-1, N.draw)) #I don't need N.CD here, because I'm drawing parameters

      for (which.draw in 1: N.draw){
        count <- 1 #a count number to asist entering the parameters

        for (i1 in 1:P){
          for (j1 in 1:(T.NFactor-1)){

            if (est_lambda[i1, j1, which.draw] == 1 ){
              est_lambda[i1, j1, which.draw] <- posterior.draws[which.draw,][count]
              count <- count+1
            }

            if (i1 == P & j1 == (T.NFactor-1) ){ count <- 0 } #reset the count to 0 for the next draw

          }
        }
      }

      #save the estimated unique variances
      est_unique <- matrix(NA, nrow=N.draw, ncol=P)

      for (which.draw in 1: N.draw){
        for (i1 in P:1){
          est_unique[which.draw, i1] <- posterior.draws[which.draw, ncol(posterior.draws)- (P-i1)]
        }
      }

      #get the eigen values for each CD.
      for (which.draw in 1: N.draw){
        for (which.cd2 in 1:N.CD){
          comp.data <- gen.pa.data(N, P, est_lambda[,,which.draw], est_unique[which.draw,])
          cor.comp.data <- cor(comp.data) # compute the correlation matrix
          eigs.comp.array[ , which.cd2, which.draw, T.NFactor] <- fa(cor.comp.data, nfactors=1, rotate="none", SMC=TRUE, max.iter=1, fm="pa", warnings = FALSE)$values
        }

        #decide if we want to move to the next number of factor given this parameter draw
        if(eigs.data[T.NFactor]> as.numeric(quantile(eigs.comp.array[T.NFactor, , which.draw, T.NFactor], quantile))) {
          sav_result[which.draw, T.NFactor] <- 1
        } else {
          sav_result[which.draw, T.NFactor] <- 0
        }
      }

      #decide if we want to move to the next number of factor
      if( sum(sav_result[, T.NFactor], na.rm = TRUE)>0 ) {
        T.NFactor <- T.NFactor+1
      } else break
    }

    # Testing the 3rd, 4th... factor
    if(T.NFactor>2 & T.NFactor<=F.Max){

      print(paste0("Evaluating factor ", T.NFactor, " of ", F.Max))

      Msyntax1 <- triag_model(P,  T.NFactor-1) #get the model syntax to run the CFA in lavaan

      # Create the list of constraints for MCMCpack
      # Initiate file and start sending text
      sink("temp1.R")
      cat("list.of.lambda.constraints <- list(")
      cat("\n") # go to next line

      # loop over constraints that fix loadings to 0
      which.constraint=1
      for(which.constraint in 1:nrow(Msyntax1$output_fixed)){
        cat(
          "V", Msyntax1$output_fixed[which.constraint, 2], "=list(", Msyntax1$output_fixed[which.constraint, 1], ",0),", sep=""
        )
        cat("\n") # go to next line
      }

      cat("\n") # go to next line

      # loop over constraints that fix loadings to be > 0
      which.constraint=1
      for(which.constraint in 1:nrow(Msyntax1$output_constraint)){

        # if it's not the last constraint, include a comma at the end
        if(which.constraint < nrow(Msyntax1$output_constraint)){
          cat(
            "V", Msyntax1$output_constraint[which.constraint, 2], "=list(", Msyntax1$output_constraint[which.constraint, 1], ',"+"),', sep=""
          )
          cat("\n") # go to next line
        }

        # if it's the last constraint, don't include a comma at the end, but do close the whole list
        if(which.constraint == nrow(Msyntax1$output_constraint)){
          cat(
            "V", Msyntax1$output_constraint[which.constraint, 2], "=list(", Msyntax1$output_constraint[which.constraint, 1], ',"+")', sep=""
          )
          cat("\n") # go to next line
          cat(")")
        }

      }

      sink()

      # Run and then remove the file with source code to define the constraints for MCMCpack
      source("temp1.R")
      file.remove("temp1.R")

      #list.of.lambda.constraints

      posterior.draws <- MCMCfactanal(
        x = data,
        factors = T.NFactor-1,
        lambda.constraints = list.of.lambda.constraints,
        burnin = n.burnin,
        mcmc = N.draw*n.thin,
        thin = n.thin,
        store.scores = FALSE,
        std.var = TRUE
      )

      # generate CD based on posterior.draws
      #Msyntax1$output_fixed #indicate which parameters should be fixed = 0

      #save the estimated loading parameters
      est_lambda <- array(1, c(P, T.NFactor-1, N.draw))

      #for fixed parameters, put a value of 0; for freely estimated parameters, keep it 1.
      for (which.draw in 1: N.draw){
        for (i1 in 1: nrow(Msyntax1$output_fixed)){
          est_lambda[Msyntax1$output_fixed[i1,2], Msyntax1$output_fixed[i1,1], which.draw] <- 0
        }
      }

      for (which.draw in 1: N.draw){
        count <- 1 #a count number to asist entering the parameters

        for (i1 in 1:P){
          for (j1 in 1:(T.NFactor-1)){

            if (est_lambda[i1, j1, which.draw] == 1 ){
              est_lambda[i1, j1, which.draw] <- posterior.draws[which.draw,][count]
              count <- count+1
            }

            if (i1 == P & j1 == (T.NFactor-1) ){ count <- 0 } #reset the count to 0 for the next draw

          }
        }
      }


      #save the estimated unique variances
      est_unique <- matrix(NA, nrow=N.draw, ncol=P)

      for (which.draw in 1: N.draw){
        for (i1 in P:1){
          est_unique[which.draw, i1] <- posterior.draws[which.draw, ncol(posterior.draws)- (P-i1)]
        }
      }

      #
      #get the eigen values for each draw.
      for (which.draw in 1: N.draw){
        for (which.cd2 in 1:N.CD){
          comp.data <- gen.pa.data(N, P, est_lambda[,,which.draw], est_unique[which.draw,])
          cor.comp.data <- cor(comp.data) # compute the correlation matrix
          eigs.comp.array[ , which.cd2, which.draw, T.NFactor] <- fa(cor.comp.data, nfactors=1, rotate="none", SMC=TRUE, max.iter=1, fm="pa", warnings = FALSE)$values
        }

        #decide if we want to move to the next number of factor given this parameter draw
        if(eigs.data[T.NFactor]> as.numeric(quantile(eigs.comp.array[T.NFactor, , which.draw, T.NFactor], quantile))) {
          sav_result[which.draw, T.NFactor] <- 1
        } else {
          sav_result[which.draw, T.NFactor] <- 0
        }
      }

      #decide if we want to move to the next number of factor
      if( sum(sav_result[, T.NFactor], na.rm = TRUE)>0 ) {
        T.NFactor <- T.NFactor+1
      } else break
    }
  }

  n.factors = T.NFactor-1 #the max num to be retained


  sav_result <- apply(sav_result, 2, sum)


  if (n.factors == 0){
    sav_result[1] <- 0
  }else{sav_result[1] <- N.draw}

  sav_result[is.na(sav_result)] <- 0

  #vector saving the draws suggesting not moving to the next factor
  sav_result_not <- N.draw - sav_result


  #vector saving the probabilities for each number of factors
  prob_baes <- matrix(rep(1, 1*(1+F.Max)), nrow=1)

  for (i in 0: F.Max){
    if (i == 0){
      prob_baes[(i+1)] <- ifelse(sav_result_not[(i+1)]==0, 0, 100) #prob of keeping 0 factor; either 0 or 100%.
    }

    if (i > 0){
      for (j in 1:i){
        prob_baes[(i+1)] <- prob_baes[(i+1)] * (sav_result[(j)]/N.draw)
      }
      prob_baes[(i+1)] <- prob_baes[(i+1)] *  sav_result_not[(i+1)]/N.draw
    }
  }

  # Specify the column names
  c.names <- rep(NA, ncol(prob_baes))
  for (i in 1: ncol(prob_baes)){
    c.names[i] <- paste("prob_", i-1, "_factors", sep = "")
  }

  # Convert to percentages for reporting, if requested
  if(report=="percentage"){
    prob_baes <- prob_baes*100
    c.names <- rep(NA, ncol(prob_baes))
    for (i in 1: ncol(prob_baes)){
      c.names[i] <- paste("%_", i-1, "_factors", sep = "")
    }
  }


  colnames(prob_baes) <- c.names


  return(prob_baes)
}

# Run the Example --------------------

#read in the example data for illustration
#data <- read.csv("HolzingerSwinefore_Grant-White_19tests.csv")

#run BPA using the bayes.pa function
#bayes.pa(data)






