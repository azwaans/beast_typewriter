library(extraDistr)

# Compute the probabilities of observing
nInsertions = 0:5
expectedRate = 0.13
deltaT = 25
lambda =expectedRate * deltaT
maxNInsertions = 5

# Compute the probabilities for obsering 0:5 insertions without truncation
dtpois(nInsertions, lambda = lambda)

# Compute the probabilities for obsering 0:5 insertions with truncation
dtpois(nInsertions, lambda = lambda, b=maxNInsertions)
##NOTE this differs from previous computation on all elements!

# Asserting that the prob of observing 5 insertions is 0, when truncating
# at 4
dtpois(nInsertions, lambda = lambda, b=maxNInsertions-1)


