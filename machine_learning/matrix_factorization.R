.libPaths('C:/Users/BrBl1834/R/win-library')
library(tidyverse)

# set.seed(1987)
#if using R 3.6 or later, use `set.seed(1987, sample.kind="Rounding")` instead
set.seed(1987, sample.kind="Rounding")
n <- 100
k <- 8
Sigma <- 64  * matrix(c(1, .75, .5, .75, 1, .5, .5, .5, 1), 3, 3) 
m <- MASS::mvrnorm(n, rep(0, 3), Sigma)
m <- m[order(rowMeans(m), decreasing = TRUE),]
y <- m %x% matrix(rep(1, k), nrow = 1) + matrix(rnorm(matrix(n*k*3)), n, k*3)
colnames(y) <- c(paste(rep("Math",k), 1:k, sep="_"),
                 paste(rep("Science",k), 1:k, sep="_"),
                 paste(rep("Arts",k), 1:k, sep="_"))

# Q1
my_image <- function(x, zlim = range(x), ...){
  colors = rev(RColorBrewer::brewer.pal(9, "RdBu"))
  cols <- 1:ncol(x)
  rows <- 1:nrow(x)
  image(cols, rows, t(x[rev(rows),,drop=FALSE]), xaxt = "n", yaxt = "n",
        xlab="", ylab="",  col = colors, zlim = zlim, ...)
  abline(h=rows + 0.5, v = cols + 0.5)
  axis(side = 1, cols, colnames(x), las = 2)
}

my_image(y)

# Q2
my_image(cor(y), zlim = c(-1,1))
range(cor(y))
axis(side = 2, 1:ncol(y), rev(colnames(y)), las = 2)

# Q3
s <- svd(y)
names(s)

y_svd <- s$u %*% diag(s$d) %*% t(s$v)
max(abs(y - y_svd))

ss_y <- colSums(y^2)
ss_yv <- colSums((y%*%s$v)^2)
sum(ss_y)
sum(ss_yv)

# Q4
plot(1:ncol(y),ss_y)
plot(1:ncol(y),ss_yv)

# Q5
plot(s$d,sqrt(ss_yv))

# Q6

# first 3 columns of YV
ss_yv_first3 <- colSums((y%*%s$v[,1:3])^2)

# total variation explained by first 3 of YV
sum(ss_yv_first3)/sum(ss_y)

# Q7

identical(s$u %*% diag(s$d), sweep(s$u, 2, s$d, FUN = "*"))

# Q8
s$u %*% diag(s$d)

# first component of UD
U1d1.1 <- s$u[,1]*s$d[1]

# average per student
student_avgs <- rowMeans(y)

plot(U1d1.1,student_avgs)

# Q9
cols <- 1:ncol(s$v)
rows <- 1:nrow(s$v)
my_image(s$v)

# Q10

plot(s$u[,1])
range(s$u[,1])

plot(s$v[,1])
plot(t(s$v[,1]))

U1d1.1.tsv1 <- U1d1.1 %*% t(s$v[,1])
my_image(U1d1.1.tsv1)
my_image(y)

# Q11
resid <- y - with(s,(u[, 1, drop=FALSE]*d[1]) %*% t(v[, 1, drop=FALSE]))
my_image(cor(resid), zlim = c(-1,1))
axis(side = 2, 1:ncol(y), rev(colnames(y)), las = 2)

plot(s$u[,2])
range(s$u[,2])

plot(s$v[,2])
plot(t(s$v[,2]))

U2d2.2 <- s$u[,2]*s$d[2]

U2d2.2.tsv2 <- U2d2.2 %*% t(s$v[,2])
my_image(U2d2.2.tsv2)
my_image(y)

# Q12
# variance explained by first two columns
sum(s$d[1:2]^2)/sum(s$d^2) * 100

resid <- y - with(s,sweep(u[, 1:2], 2, d[1:2], FUN="*") %*% t(v[, 1:2]))
my_image(cor(resid), zlim = c(-1,1))
axis(side = 2, 1:ncol(y), rev(colnames(y)), las = 2)

plot(s$u[,3])
range(s$u[,3])

plot(s$v[,3])
plot(t(s$v[,3]))

U3d3.3 <- s$u[,3]*s$d[3]

U3d3.3.tsv3 <- U3d3.3 %*% t(s$v[,3])
my_image(U3d3.3.tsv3)
my_image(y)

# Q13
sum(s$d[1:3]^2)/sum(s$d^2) * 100

resid <- y - with(s,sweep(u[, 1:3], 2, d[1:3], FUN="*") %*% t(v[, 1:3]))
my_image(cor(resid), zlim = c(-1,1))
axis(side = 2, 1:ncol(y), rev(colnames(y)), las = 2)

my_image(U1d1.1.tsv1 + U2d2.2.tsv2 +U3d3.3.tsv3)
my_image(y)
my_image(y-(U1d1.1.tsv1 + U2d2.2.tsv2 +U3d3.3.tsv3))

y_hat <- with(s,sweep(u[, 1:3], 2, d[1:3], FUN="*") %*% t(v[, 1:3]))
my_image(y, zlim = range(y))
my_image(y_hat, zlim = range(y))
my_image(y - y_hat, zlim = range(y))