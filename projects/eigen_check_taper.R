## check the trend of eigenvalues ##
# numerical check for assumption 4&5 #
rm(list = ls())
library(fields)
library(ggplot2)
library(Matrix)

Cov_Matern_eigen_taper <- 
  function(N_side = N_side,  # min within dist is 1/N_side
           dim = d,          # dimension
           nu = nu,          # smoothness
           r = r,            # range in tapering
           phi = 1.0){
  #' return the eigenvalues of the Matern covarianc matrix on a grid
  
  if(d == 1){
    coords = seq(from = 0, to = 1, length.out = N_side + 1)[-1]
  }else if(d == 2){
    coords_side = seq(from = 0, to = 1, length.out = N_side + 1)[-1]
    coords = expand.grid(coords_side, coords_side)
  }
  D <- as.matrix(dist(coords)) 
  # Wendland family of tapered covariance functions
  R <- Matern(D, nu = nu, phi = 1.0, alpha = phi) * 
    Wendland(D, theta = r, dimension = d, k = 2) + diag(N_side^d)
  R <- Matrix(R, sparse = TRUE)    # convert to sparse matrix
  eigen(R)$values - 1.0
}


## check d = 1 ##
d = 1
nu_list = c(0.9, 1.5)#c(0.9, 1.5)
n_list = c(500, 1000, 1500, 2000, 2500, 3000, 3500, 4000)
L = length(n_list)
check_order <- c(0.7, 0.8, 0.9) #c(0.5, 0.75, 0.9) 

lo <- length(check_order)
lnu <- length(nu_list)

check_order_ls <- rep(rep(
  sapply(check_order, f <- function(x){paste0("n^", x)}), L), lnu)
check_N <- rep(rep(n_list, each = length(check_order)), lnu)
nu_ls <- rep(nu_list, each = lo * L)

eigen_check <- rep(NA, lnu * L * lo)
for (k in 1:lnu){
  cat("\n", k, ": \t")
  for(j in 1:L){
    cat(j, "\t")
    N = n_list[j]
    e <- Cov_Matern_eigen_taper(N_side = N, dim = d, nu = nu_list[k], r = 0.5)
    for (i in c(1:(N^d))){
      e[i] <- e[i] / (N^d) * i^(2 * nu_list[k] / d + 1)
    }
    eigen_check[ (k - 1) * (L * lo) + ((j - 1) * lo + 1):(j * lo) ] = 
      c(e[floor((N^d)^check_order)])
  }
}

assump_dat_dim1 <- 
  data.frame(eigen_check = eigen_check,
             order = check_order_ls,
             check_N = check_N, 
             nu_ls = nu_ls,
             group = sapply(1:(lo * L * lnu), 
                            f <- function(x){
                              paste0("nu = ", nu_ls[x], ", i = ", 
                                     check_order_ls[x])}))
assump_dat_dim1$nu_ls <- 
  factor(assump_dat_dim1$nu_ls, levels = check_order,
         labels = sapply(check_order, f <- function(x){
           expression(nu~"="~x)}))

p_dim1 <- ggplot(data = assump_dat_dim1, 
                 aes(x = check_N, y = eigen_check, group = group)) + 
  geom_line(aes(linetype = group)) + 
  geom_point(aes(shape = group)) + 
  scale_shape_discrete(
    name  ="",
    breaks=c("nu = 0.9, i = n^0.7", "nu = 0.9, i = n^0.8", 
             "nu = 0.9, i = n^0.9", "nu = 1.5, i = n^0.7",
             "nu = 1.5, i = n^0.8", "nu = 1.5, i = n^0.9"),
    labels=c(expression(nu ~ "=" ~ 0.9 ~ ","~ alpha ~ "=" ~ 0.7),
             expression(nu ~ "=" ~ 0.9 ~ ","~ alpha ~ "=" ~ 0.8), 
             expression(nu ~ "=" ~ 0.9 ~ ","~ alpha ~ "=" ~ 0.9),
             expression(nu ~ "=" ~ 1.5 ~ ","~ alpha ~ "=" ~ 0.7), 
             expression(nu ~ "=" ~ 1.5 ~ ","~ alpha ~ "=" ~ 0.8),
             expression(nu ~ "=" ~ 1.5 ~ ","~ alpha ~ "=" ~ 0.9))) +
  scale_linetype_discrete(
    name  ="",
    breaks=c("nu = 0.9, i = n^0.7", "nu = 0.9, i = n^0.8", 
             "nu = 0.9, i = n^0.9", "nu = 1.5, i = n^0.7",
             "nu = 1.5, i = n^0.8", "nu = 1.5, i = n^0.9"),
    labels=c(expression(nu ~ "=" ~ 0.9 ~ ","~ alpha ~ "=" ~ 0.7),
             expression(nu ~ "=" ~ 0.9 ~ ","~ alpha ~ "=" ~ 0.8), 
             expression(nu ~ "=" ~ 0.9 ~ ","~ alpha ~ "=" ~ 0.9),
             expression(nu ~ "=" ~ 1.5 ~ ","~ alpha ~ "=" ~ 0.7), 
             expression(nu ~ "=" ~ 1.5 ~ ","~ alpha ~ "=" ~ 0.8),
             expression(nu ~ "=" ~ 1.5 ~ ","~ alpha ~ "=" ~ 0.9))) +
  scale_x_continuous(name = "n") + theme_bw() + ylab(" ") #+ ylim(0.0, NA)

width <- 8.0
height <- 3.0
pointsize <- 16

setEPS()
postscript("./pics/assump3_check.eps", width = width, height = height)
print(p_dim1)
dev.off()
  
####### update plots for jrssb #######
ggsave("./pics/Figure2.pdf", 
       plot = p_dim1, width=width, height=height, dpi = 1000)


## check d = 2 ##
d = 2
nu_list = c(0.5, 1.0)
n_list = c(30, 50, 60, 70, 80, 90, 100, 110)
L = length(n_list)
check_order <- c(0.7, 0.8, 0.9)

lo <- length(check_order)
lnu <- length(nu_list)

check_order_ls <- rep(rep(
  sapply(check_order, f <- function(x){paste0("n^", x)}), L), lnu)
check_N <- rep(rep(n_list^d, each = length(check_order)), lnu)
nu_ls <- rep(nu_list, each = lo * L)

eigen_check <- rep(NA, lnu * L * lo)
for (k in 1:lnu){
  cat("\n", k, ": \t")
  for(j in 1:L){
    cat(j, "\t")
    N = n_list[j]
    e <- Cov_Matern_eigen_taper(N_side = N, dim = d, nu = nu_list[k], r = 0.5, 
                                phi = 18.0)
    e_c <- c()
    for (i in floor((N^d)^check_order)){
      e_c <- c(e_c, e[i] / (N^d) * i^(2 * nu_list[k] / d + 1))
    }
    eigen_check[ (k - 1) * (L * lo) + ((j - 1) * lo + 1):(j * lo) ] = e_c
  }
}

assump_dat_dim2 <- 
  data.frame(eigen_check = eigen_check,
             order = check_order_ls,
             check_N = check_N, 
             nu_ls = nu_ls,
             group = sapply(1:(lo * L * lnu), 
                            f <- function(x){
                              paste0("nu = ", nu_ls[x], ", i = ", 
                                     check_order_ls[x])}))
assump_dat_dim2$nu_ls <- 
  factor(assump_dat_dim2$nu_ls, levels = check_order,
         labels = sapply(check_order, f <- function(x){
           expression(nu~"="~x)}))

p_dim2 <- ggplot(data = assump_dat_dim2, 
                 aes(x = check_N, y = eigen_check, group = group)) + 
  geom_line(aes(linetype = group)) + 
  geom_point(aes(shape = group)) + theme_bw() + 
  ylab(" ") + xlab("n") + theme(legend.position="right", 
                                legend.title = element_blank())

width <- 8.0
height <- 3.0
pointsize <- 16

setEPS()
postscript("./pics/assump3_check2.eps", width = width, height = height)
print(p_dim2)
dev.off()


