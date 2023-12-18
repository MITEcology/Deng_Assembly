##############################################################################
# An example for the manuscript "The development of ecological systems along #
# paths of least resistance" by Jie Deng, Otto X. Cordero, Tadashi Fukami,   #
# Simon A. Levin, Robert M. Pringle, Ricard Solé, and Serguei Saavedra       #
##############################################################################

source("toolbox.R")

# one can repeat this example for 10000 random systems
# with different standard deviation (sdlog_mass) to reproduce Fig 2B

# number of taxa in a random system
num_sp <- 6 

# parameters for the log-normal distribution of body mass
meanlog_mass <- 0 

# parameters for the log-normal distribution of randomly generated body mass
# one can repeat this example using a range of sdlog_mass 
# (e.g., sdlog_mass in [0.01, 2]) to reproduce Fig 2B
sdlog_mass <- 0.5

# scaling relationship between body mass and rmax for bacteria 
# use alpha = -0.25 for metazoans
alpha <- 0.75
  
# generate body mass randomly using log-normal distribution
body_mass <- sort(rlnorm(num_sp, meanlog_mass, sdlog_mass))

# generate rmax using body mass following the scaling relationship
r_max <- body_mass^{alpha}

# calculate the standard deviation of rmax
sd(r_max)

# generate interaction matrix using rmax
A <- ScalingInteractionMatrixRmax(rmax = r_max, stable = TRUE)

# calculate transition matrix
# used in Fig 2 and the pink star of Fig 3
transition_matrix <- TransitionMatrix(A = A, omega_rep = 100)

# calculate developmental probability
# used in Fig 2A
developmental_probability <- DevelopmentalProbability(A, transition_matrix)
sequential_development <- developmental_probability$avg_prob_per_step # sequential developmental probability
simultaneous_development <- developmental_probability$simu_prob # simultaneous developmental probability

# calculate optimal path
# used in the y-axis of Fig 2B and the pink star of Fig 3
optimal_path <- OptimalPath(A = A, transition_matrix = transition_matrix) # 

# calculate the Spearman's rank correlation between optimal path and rmax
# calculated 10000 correlations for the y-axis of Fig 2B and then averaged by rmax groups
cor_rmax_optimal <- cor(r_max, optimal_path, method = "spearman")
