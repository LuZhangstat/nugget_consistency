## check the trend of eigenvalues ##
rm(list = ls())
library(fields)
library(ggplot2)

Cov_Matern_eigen <- function(N_side = N_side,  # min within dist is 1/N_side
                             dim = d,          # dimension
                             nu = nu,          # smoothness
                             phi = 1.0){
  #' return the eigenvalues of the Matern covarianc matrix on a grid
  
  if(d == 1){
    coords = seq(from = 0, to = 1, length.out = N_side + 1)[-1]
  }else if(d == 2){
    coords_side = seq(from = 0, to = 1, length.out = N_side + 1)[-1]
    coords = expand.grid(coords_side, coords_side)
  }
  D <- as.matrix(dist(coords))
  R <- Matern(D, nu = nu, phi = 1.0, alpha = phi) + diag(N_side^d)
  eigen(R)$values - 1.0
}

## check d = 1 ##
d = 1
nu_list = c(0.9, 1.5)
n_list = c(100, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200,
           2400, 2600, 2800, 3000)
L = length(n_list)
check_order <- c(0.5, 0.75, 0.9)

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
    e <- Cov_Matern_eigen(N_side = N, dim = d, nu = nu_list[k])
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
    breaks=c("nu = 0.9, i = n^0.5", "nu = 0.9, i = n^0.75", 
             "nu = 0.9, i = n^0.9", "nu = 1.5, i = n^0.5",
             "nu = 1.5, i = n^0.75", "nu = 1.5, i = n^0.9"),
    labels=c(expression(nu ~ "=" ~ 0.9 ~ ","~ alpha ~ "=" ~ 0.5),
             expression(nu ~ "=" ~ 0.9 ~ ","~ alpha ~ "=" ~ 0.75), 
             expression(nu ~ "=" ~ 0.9 ~ ","~ alpha ~ "=" ~ 0.9),
             expression(nu ~ "=" ~ 1.5 ~ ","~ alpha ~ "=" ~ 0.5), 
             expression(nu ~ "=" ~ 1.5 ~ ","~ alpha ~ "=" ~ 0.75),
             expression(nu ~ "=" ~ 1.5 ~ ","~ alpha ~ "=" ~ 0.9))) +
  scale_linetype_discrete(
    name  ="",
    breaks=c("nu = 0.9, i = n^0.5", "nu = 0.9, i = n^0.75", 
             "nu = 0.9, i = n^0.9", "nu = 1.5, i = n^0.5",
             "nu = 1.5, i = n^0.75", "nu = 1.5, i = n^0.9"),
    labels=c(expression(nu ~ "=" ~ 0.9 ~ ","~ alpha ~ "=" ~ 0.5),
             expression(nu ~ "=" ~ 0.9 ~ ","~ alpha ~ "=" ~ 0.75), 
             expression(nu ~ "=" ~ 0.9 ~ ","~ alpha ~ "=" ~ 0.9),
             expression(nu ~ "=" ~ 1.5 ~ ","~ alpha ~ "=" ~ 0.5), 
             expression(nu ~ "=" ~ 1.5 ~ ","~ alpha ~ "=" ~ 0.75),
             expression(nu ~ "=" ~ 1.5 ~ ","~ alpha ~ "=" ~ 0.9))) +
  scale_x_continuous(name = "n") + theme_bw() + ylab(" ") + ylim(0.0, NA)

width <- 8.0
height <- 3.0
pointsize <- 16

setEPS()
postscript("./pics/assump_check.eps", width = width, height = height)
print(p_dim1)
dev.off()



## check d = 2 ##
d = 2
nu_list = c(0.5, 1.5)
n_list = c(10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60)
L = length(n_list)
check_order <- c(0.5, 0.75, 0.95)

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
    e <- Cov_Matern_eigen(N_side = N, dim = d, nu = nu_list[k])
    for (i in c(1:(N^d))){
      e[i] <- e[i] / (N^d) * i^(2 * nu_list[k] / d + 1)
    }
    eigen_check[ (k - 1) * (L * lo) + ((j - 1) * lo + 1):(j * lo) ] = 
      c(e[floor((N^d)^check_order)])
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
postscript("./pics/assump_check2.eps", width = width, height = height)
print(p_dim2)
dev.off()



# 
# 
# threehalf <- function(x) {
#   z <- (1 + abs(x))*exp(-abs(x))
#   return(z)
# }
# 
# covmat <- function(N) {
#   s <- matrix(, nrow = N, ncol = N)
#   for (i in c(1:N)){
#     for (j in c(1:N)){
#       s[i,j] <- threehalf(i/N-j/N)
#     }
#   }
#   return(s)
# }
# 
# 
# a <- rep(0, 10)
# b <- rep(0, 10)
# c <- rep(0, 10)
# for (j in c(1:10)){
#   N = j*100
#   m <- covmat(N)
#   e <- eigen(m)$values
#   for (i in c(1:N)){
#     e[i] <- e[i]/N*i^4
#   }
#   a[j] = e[floor(N^0.5)]
#   b[j] = e[floor(N^0.75)]
#   c[j] = e[floor(N^0.95)]
# }
# 
# dat <-  c(rbind(a,b,c))
