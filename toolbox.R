#############################################################################
# Functions for the manuscript "The development of ecological systems along #
# paths of least resistance" by Jie Deng, Otto X. Cordero, Tadashi Fukami,  #
# Simon A. Levin, Robert M. Pringle, Ricard Solé, and Serguei Saavedra      #
#############################################################################

# load necessary packages
suppressPackageStartupMessages({
  if(!require(mvtnorm)) {install.packages("mvtnorm"); library(mvtnorm)}
  if(!require(mgcv)) {install.packages("mgcv"); library(mgcv)}
  if(!require(magrittr)) {install.packages("magrittr"); library(magrittr)}
  if(!require(coneproj)) {install.packages("coneproj"); library(coneproj)} # use function `check_irred(mat)` to check 
  if(!require(gtools)) {install.packages("gtools"); library(gtools)} # permutation with repeats
  if(!require(dplyr)) {install.packages("dplyr"); library(dplyr)}
  if(!require(geometry)) {install.packages("geometry"); library(geometry)}
  if(!require(pracma)) {install.packages("pracma"); library(pracma)}
  if(!require(deSolve)) {install.packages("deSolve"); library(deSolve)}
  if(!require(igraph)) {install.packages("igraph"); library(igraph)}
  if(!require(markovchain)) {install.packages("markovchain"); library(markovchain)}
  if(!require(expm)) {install.packages("expm"); library(expm)}
  if(!require(Matrix)) {install.packages("Matrix"); library(Matrix)}
})

# function that computes the feasibility of a system governed by the interaction
# matrix A analytically
# inputs: A = interaction matrix
# output: out = the feasibility of the system
Omega <- function(A) {
  A <- as.matrix(A)
  S <- nrow(A)
  # omega function
  omega <- function(S, Sigma) {
    m <- matrix(0, S, 1)
    a <- matrix(0, S, 1)
    b <- matrix(Inf, S, 1)
    d <- pmvnorm(lower = rep(0, S), upper = rep(Inf, S), mean = rep(0, S), sigma = Sigma)
    out <- d[1] # community level
    return(out)
  }
  # rule out errors
  f <- function(m) class(try(solve(t(m) %*% m), silent = T)) == "matrix"
  if (all(f(A) == FALSE)) {
    return(0)
  }
  else {
    Sigma <- solve(t(A) %*% A)
    return(omega(S, Sigma))
  }
}

# function that computes the feasibility of a system governed by the interaction
# matrix A analytically for many times to remove numerical uncertainty
# inputs: A = interaction matrix
#         num = times of repeat
# output: out = out = the feasibility of the system
OmegaAvg <- function(A, num){
  if(num == FALSE){
    out <- Omega(A)
  }else{
    out <- 0
    for(i in 1:num){
      out <- out + Omega(A = A)
    }
    out <- out/num
  }
  out
}

# function that check the global stability
# inputs: A = the interaction matrix of the entire system under gLV dynamics
# output: TRUE: globally stable
CheckGlobalStability <- function(A){
  all(eigen(A+t(A))$values < 0)
}


# function that generate the interaction matrix by rmax
# inputs: rmax: maximum growth rate of all taxa in the system
#         stable: TRUE = the output interaction matrix is globally stable
# output: interaction_matrix = the generated interaction matrix
ScalingInteractionMatrixRmax <- function(rmax, stable = TRUE, scalar = 0.9){
  num_sp <- length(rmax)
  interactions <- runif(num_sp*num_sp)
  interaction_matrix <- matrix(interactions, ncol = num_sp, nrow = num_sp)
  rmax <- rmax / sum(rmax)
  for (a in 1:num_sp){
    for (b in 1:num_sp){
      interaction_matrix[a,b] <- (rmax[b]/rmax[a])*(1/sqrt(num_sp))
    }
  }
  diag(interaction_matrix) <- 1 # per niche framework
  interaction_matrix <- abs(interaction_matrix)*(-1)  # competition
  if(stable){
    while(!CheckGlobalStability(interaction_matrix)){
      interaction_matrix <- interaction_matrix * scalar
      diag(interaction_matrix) <- -1
    }
  }
  return(interaction_matrix)
}


# function that computes the transition matrix of an one-by-one bottom-up 
# developmental process using feasibility measures
# inputs: A = pairwise interaction matrix of the entire system
# output: transition_matrix["i","j"] = probability of transitioning to system j from system i
TransitionMatrix <- function(A, omega_rep = 30){
  S <- nrow(A) # number of taxa
  index_comm <- c() # index of systems
  for(i in 1:S){
    index_comm <- c(index_comm, combn(S, i, simplify = FALSE))
  }
  index_comm_char <- lapply(index_comm, function(x) paste0(x, collapse = "-"))
  num_comm <- length(index_comm_char)
  # initialize transition matrix
  transition_matrix <- matrix(0, nrow = num_comm, ncol = num_comm)
  colnames(transition_matrix) <- index_comm_char
  rownames(transition_matrix) <- index_comm_char
  # put probabilities into transition matrix by rows
  for(n_comm in 1:(num_comm-1)){
    in_comm <- index_comm[[n_comm]]
    S_in <- length(in_comm) # maximum number of taxa
    out_one_more <- combn(S, S_in+1, simplify = FALSE)
    out_comm <- out_one_more[sapply(out_one_more, function(x) all(in_comm %in% x))]
    prob <- sapply(out_comm, function(x) OmegaAvg(A = A[x,x], num = omega_rep))
    prob <- prob/sum(prob)
    stay <- OmegaAvg(A = A[in_comm,in_comm], num = omega_rep)
    transition_matrix[n_comm, n_comm] <- stay
    out_comm_char <- sapply(out_comm, function(x) paste0(x, collapse = "-"))
    transition_matrix[n_comm, out_comm_char] <- (1-stay)*prob
  }
  err_row <- which(rowSums(transition_matrix)==0)
  transition_matrix[err_row, err_row] <- 1
  transition_matrix/rowSums(transition_matrix)
}


# function that calculate the average probability of "123..(all)" from 
# initial system with different number of taxa by sequential development 
# until reaching the feasibility by simultaneous development
# input: A = pairwise interaction matrix of the entire system
# outputs: avg_prob_per_step = probability of developing "123..(all)" vs # steps
#          critical_step = critical # step of reaching omega by simultaneous development, 
#          omega = feasibility by simultaneous development
DevelopmentalProbability <- function(A, transition_matrix, max_num_step = 20){
  S <- nrow(A)
  # index of initial systems with different number of taxa
  index_in_comm_char <- list()
  initial_prob <- vector("list", S-1)
  num_comm <- c()
  for(i in 1:(S-1)){
    # index of initial systems that we can start with
    index_comm <- combn(S, i, simplify = FALSE)
    num_comm <- c(num_comm, length(index_comm))
    for(comm in index_comm){
      initial_prob[[i]] <- c(initial_prob[[i]], Omega(A[comm, comm]))
    }
    index_in_comm_char[[i]] <- sapply(index_comm, function(x) paste0(x, collapse = "-"))
  }
  # index of the final system that we would like to develop
  index_out_comm_char <- paste0(c(1:S), collapse = "-")
  # initialize the params in loop
  seq_prob <- rep(-1, S-1)
  seq_prob_list <- as.list(seq_prob)
  num_step <- 1
  omega <- Omega(A)
  simu_prob <- c()
  while(num_step <= max_num_step){ # threshold to stop the loop
    mat_prob <- transition_matrix %^% num_step
    for(i in c(1:(S-1))){
      seq_prob_list[[i]] <- c(seq_prob_list[[i]], 
                              initial_prob[[i]] %*% mat_prob[index_in_comm_char[[i]], index_out_comm_char]/num_comm[i])
    }
    simu_prob <- c(simu_prob, 1-(1-omega)^num_step)
    num_step <- num_step + 1
  }
  df_seq_prob <- as.data.frame(seq_prob_list)
  colnames(df_seq_prob) <- paste0(c(1:(S-1)), separate = "-tx")
  df_seq_prob <- df_seq_prob[-1,]
  rownames(df_seq_prob) <- NULL
  seq_prob_list <- as.list(df_seq_prob)
  critical_step <- sapply(seq_prob_list, function(x) which(x > simu_prob)[1])
  avg_initial_prob <- sapply(initial_prob, function(x) mean(x))
  return(list(avg_prob_per_step = df_seq_prob, simu_prob = simu_prob, critical_step = critical_step,
              omega = omega, avg_initial_prob = avg_initial_prob))
}


# function that generate the optimal path of ecological development
# inputs: A: interaction matrix of n taxa, e.g., taxon index 1,2,...,n
#         transition_matrix: the transition matrix of bottom-up assembly using omega
# output: optimal developmental order using taxon index
OptimalPath <- function(A, transition_matrix){
  # index of all systems
  num_sp <- nrow(A)
  index_comm <- c()
  for(i in 1:num_sp){
    index_comm <- c(index_comm, combn(num_sp, i, simplify = FALSE))
  }
  num_comm <- nrow(transition_matrix)
  rownames(transition_matrix) <- NULL
  colnames(transition_matrix) <- NULL
  
  # product of conditional probabilities
  adj_trans_mat <- -log(transition_matrix)
  g_mc <- graph_from_adjacency_matrix(adj_trans_mat, mode="directed", weighted=TRUE)
  
  # start with 1sp systems
  prob <- sapply(1:nrow(A), function(i) {
    exp(-distances(g_mc, v = i, to = num_comm))
  })
  initial_1sp <- which.max(prob) # the initial taxon of the optimal path
  assembly_order_ind <- shortest_paths(g_mc, from = initial_1sp, to = num_comm)$vpath[[1]]
  assembly_order_comm <- index_comm[assembly_order_ind]
  optimal_path <- assembly_order_comm[[1]]
  for(n_comm in 2:length(assembly_order_comm)){
    optimal_path <- c(optimal_path, 
                      setdiff(assembly_order_comm[[n_comm]], 
                              assembly_order_comm[[n_comm-1]]))
  }
  optimal_path
}
######################################################
