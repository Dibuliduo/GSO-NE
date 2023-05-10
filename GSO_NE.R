# L_{2,2/3}
GroupThird2 <- function(A, b, c, x1, MAX_ITERS, n, s, p) {
  # 初始化f和error为正无穷
  f <- 0
  error <- 0
  
  # 初始化迭代次数k为1，定义一些变量G、r、lambda、x、NC和v
  k <- 1
  G <- s / p
  r <- n / p
  lambda <- 100
  x <- x1
  NC <- norm(c, "2")
  v <- 0.5
  Va1 <- (.5)^(4 / 3) * 1.5 / v
  
  # 计算Bu1和Bu2
  Bu1 <- 2 * v * t(A) %*% b
  Bu2 <- 2 * v * t(A) %*% A
  
  # 初始化BuG为全零向量
  BuG <- rep(0, r)
  
  while (k < MAX_ITERS) {
    j <- 1
    Bu <- x + Bu1 - Bu2 %*% x
    while (j <= r) {
      BuG[j] <- norm(Bu[((j - 1) * p + 1):(j * p)], "2")
      j <- j + 1
    }
    BuGO <- sort(BuG)
    BuV <- BuGO[r - G]^(4 / 3)
    lambda <- Va1 * BuV
    criterion <- BuGO[r - G]
    
    j <- 1
    while (j <= r) {
      Y1 <- norm(Bu[((j - 1) * p + 1):(j * p)], "2")
      if (Y1 > criterion) {
        q <- 2 * lambda * v
        phi <- acosh(27 * Y1^2 / (16 * q^1.5))
        aa <- 2 * q^0.25 * (cosh(phi / 3))^0.5 / (sqrt(3))
        eta <- 3 * (aa^1.5 + sqrt(2 * Y1 - aa^3))^4
        x[((j - 1) * p + 1):(j * p)] <- (eta / (32 * lambda * v * aa^2 + eta)) * Bu[((j - 1) * p + 1):(j * p)]
      } else {
        x[((j - 1) * p + 1):(j * p)] <- rep(0, p)
      }
      j <- j + 1
    }
    
    f <- c(f, norm(A %*% x - b, "2")^2)
    error <- c(error, norm(x - c, "2") / NC)
    k <- k + 1
  }
  
  hist <- list(f, error)
  return(list(x, hist))
}

# L_{2,1/2}
GroupGP2 <- function(A, b, c, x1, MAX_ITERS, n, s, p) {
  # 初始化f和error为正无穷
  f <- Inf
  error <- Inf
  
  # 初始化迭代次数k为1，定义一些变量G、r、lambda、x、NC和v
  k <- 1
  G <- s / p
  r <- n / p
  lambda <- 100
  x <- x1
  NC <- norm(c, "2")
  v <- 0.5
  
  # 计算Va1
  Va1 <- (2 / 3)^1.5 / v
  
  # 计算Bu1和Bu2
  Bu1 <- 2 * v * t(A) %*% b
  Bu2 <- 2 * v * t(A) %*% A
  
  # 初始化BuG为全零向量
  BuG <- rep(0, r)
  
  # 设置while循环条件，限制迭代次数
  while (k < MAX_ITERS) {
    # 计算Bu、BuG、BuGO和BuV
    j <- 1
    Bu <- x + Bu1 - Bu2 %*% x
    while (j <= r) {
      BuG[j] <- norm(Bu[((j - 1) * p + 1):(j * p)], "2")
      j <- j + 1
    }
    BuGO <- sort(BuG)
    BuV <- BuGO[r - G]^1.5
    lambda <- Va1 * BuV
    criterion <- BuGO[r - G]
    
    j <- 1
    while (j <= r) {
      Y1 <- norm(Bu[((j - 1) * p + 1):(j * p)], "2")
      if (Y1 > criterion) {
        q <- lambda * v / 4
        phi <- acos(q * (3 / Y1)^1.5)
        eta <- 16 * norm(Bu[((j - 1) * p + 1):(j * p)], "2")^(3 / 2) * cos((pi - phi) / 3)^3
        x[((j - 1) * p + 1):(j * p)] <- (eta / (3 * sqrt(3) * lambda * v + eta)) * Bu[((j - 1) * p + 1):(j * p)]
      } else {
        x[((j - 1) * p + 1):(j * p)] <- rep(0, p)
      }
      j <- j + 1
    }
    f <- c(f, norm(A %*% x - b, "2")^2)
    error <- c(error, norm(x - c, "2") / NC)
    k <- k + 1
  }
  
  hist <- list(f, error)
  return(list(x, hist))
}

GroupThird <- function(A, b, c, x1, MAX_ITERS, n, s, p) {
  # 初始化f和error为正无穷
  f <- Inf
  error <- Inf
  
  # 初始化迭代次数k为1，定义一些变量G、r、lambda、x、NC和v
  k <- 1
  G <- s / p
  r <- n / p
  lambda <- 100
  x <- x1
  NC <- norm(c, "2")
  v <- 0.5
  Va1 <- (.5)^(4 / 3) * 1.5 / v
  
  # 计算Bu1和Bu2
  Bu1 <- 2 * v * t(A) %*% b
  Bu2 <- 2 * v * t(A) %*% A
  
  # 初始化BuG为全零向量
  BuG <- rep(0, r)
  
 while (k < MAX_ITERS) {
    j <- 1

    Bu <- x + Bu1 - Bu2 %*% x

    while (j <= r) {
      BuG[j] <- norm(Bu[((j - 1) * p + 1):(j * p)], "2")
      j <- j + 1
    }
    BuGO <- sort(BuG)
    BuV <- BuGO[r - G]^(4 / 3)
    lambda <- Va1 * BuV
    criterion <- BuGO[r - G]
    
    j <- 1
    while (j <= r) {
        
      Y1 <- sum(abs(Bu[((j - 1) * p + 1):(j * p)]))
      Y2 <- norm(Bu[((j - 1) * p + 1):(j * p)], "2")

      if (!is.na(Y2) && !is.na(criterion) && Y2 > criterion) {
        q <- 2 * lambda * v * p
        phi <- acosh(27 * Y1^2 / (16 * q^1.5))
        aa <- 2 * q^0.25 * (cosh(phi / 3))^0.5 / (sqrt(3))
        eta <- 3 * (aa^1.5 + sqrt(2 * Y1 - aa^3))
        x[((j - 1) * p + 1):(j * p)] <- Bu[((j - 1) * p + 1):(j * p)] - (4 * v * lambda * aa^0.5 / eta) * sign(Bu[((j - 1) * p + 1):(j * p)])
      } else {
        x[((j - 1) * p + 1):(j * p)] <- rep(0, p)
      }
      j <- j + 1
    }
    f <- c(f, norm(A %*% x - b, "2")^2)
    error <- c(error, norm(x - c, "2") / NC)
    k <- k + 1
  }
  
  hist <- list(f, error)
  return(list(x, hist))
}

GroupSoft <- function(A, b, c, x1, MAX_ITERS, n, s, p) {
  # 初始化f和error为正无穷
  f <- Inf
  error <- Inf
  
  # 初始化迭代次数k为1，定义一些变量x、v、G、r和NC
  k <- 1
  x <- x1
  v <- 0.5
  G <- s / p
  r <- n / p
  NC <- norm(c, "2")
  
  # 计算Bu1和Bu2
  Bu1 <- 2 * v * t(A) %*% b
  Bu2 <- 2 * v * t(A) %*% A
  
  # 初始化BuG为全零向量
  BuG <- rep(0, r)
  
  while (k < MAX_ITERS) {
    j <- 1
    Bu <- x + Bu1 - Bu2 %*% x
    while (j <= r) {
      BuG[j] <- norm(Bu[((j - 1) * p + 1):(j * p)], "2")
      j <- j + 1
    }
    Bu0 <- sort(BuG)
    criterion <- Bu0[r - G]
    lambda <- criterion / v
    
    j <- 1
    while (j <= r) {
      if (norm(Bu[((j - 1) * p + 1):(j * p)], "2") > criterion) {
        x[((j - 1) * p + 1):(j * p)] <- (1 - v * lambda / norm(Bu[((j - 1) * p + 1):(j * p)], "2")) * Bu[((j - 1) * p + 1):(j * p)]
      } else {
        x[((j - 1) * p + 1):(j * p)] <- rep(0, p)
      }
      j <- j + 1
    }
    
    f <- c(f, norm(A %*% x - b, "2")^2)
    error <- c(error, norm(x - c, "2") / NC)
    k <- k + 1
  }
  
  hist <- list(f, error)
  return(list(x, hist))
}

GroupHard <- function(A, b, c, x1, MAX_ITERS, n, s, p) {
  # 初始化f和error为正无穷
  f <- Inf
  error <- Inf
  
  # 初始化迭代次数k为1，定义一些变量x、v、G、r和NC
  k <- 1
  x <- x1
  v <- 0.5
  G <- s / p
  r <- n / p
  NC <- norm(c, "2")
  
  # 计算Bu1和Bu2
  Bu1 <- 2 * v * t(A) %*% b
  Bu2 <- 2 * v * t(A) %*% A
  
  # 初始化BuG为全零向量
  BuG <- rep(0, r)
  
  while (k < MAX_ITERS) {
    j <- 1
    Bu <- x + Bu1 - Bu2 %*% x
    while (j <= r) {
      BuG[j] <- norm(Bu[((j - 1) * p + 1):(j * p)], "2")
      j <- j + 1
    }
    Bu0 <- sort(BuG)
    criterion <- Bu0[r - G]
    lambda <- criterion^2
    
    j <- 1
    while (j <= r) {
      if (norm(Bu[((j - 1) * p + 1):(j * p)], "2") > criterion) {
        x[((j - 1) * p + 1):(j * p)] <- Bu[((j - 1) * p + 1):(j * p)]
      } else {
        x[((j - 1) * p + 1):(j * p)] <- rep(0, p)
      }
      j <- j + 1
    }
    
    f <- c(f, norm(A %*% x - b, "2")^2)
    error <- c(error, norm(x - c, "2") / NC)
    k <- k + 1
  }
  
  hist <- list(f, error)
  return(list(x, hist))
}


# 定义变量 m, n, sequence, p, gs, s 和 MAX_ITERS
m <- 256
n <- 1024
sequence <- 1:50
p <- 4
gs <- 30
s <- p * gs
MAX_ITERS <- 500

# 创建一个矩阵 A，用正态分布随机数填充，然后进行正交化
A <- matrix(rnorm(m * n), nrow = m, ncol = n)
requireNamespace("ThreeWay")
A <- ThreeWay::orth(t(A))
A <- t(A)

# 定义 ActInd，它是从 1 到 n/p 的随机排列
ActInd <- sample(n / p)

# 初始化向量 c 为全零
c <- rep(0, n)

# 使用 for 循环，根据 ActInd 的值更新向量 c
for (ai in 1:gs) {
  c[((ActInd[ai] - 1) * p + 1):(ActInd[ai] * p)] <- rnorm(p)
}

# 定义噪声系数 sigma，计算向量 b
sigma <- 0.0001
b <- A %*% c + sigma * runif(m)

# 初始化向量 x1 为全零
x1 <- matrix(0,n,1)



results_soft <- GroupSoft(A, b, c, x1, MAX_ITERS, n, s, p)
results_hard <- GroupHard(A, b, c, x1, MAX_ITERS, n, s, p)
results_gp2 <- GroupGP2(A, b, c, x1, MAX_ITERS, n, s, p)
results_third <- GroupThird(A, b, c, x1, MAX_ITERS, n, s, p)
results_third2 <- GroupThird2(A, b, c, x1, MAX_ITERS, n, s, p)


# 加载所需的包
library(ggplot2)

# 从 GroupThird2 函数的结果中提取 error 值
error_values_Third2 <- results_third2[[2]][[2]]
error_values_gp2 <- results_gp2[[2]][[2]]
error_values_Third <- results_third[[2]][[2]]
error_values_soft <- results_soft[[2]][[2]]
error_values_hard <- results_hard[[2]][[2]]



# 计算每个点之间的间隔
interval <- floor(MAX_ITERS / 50)

# 创建一个数据框，包含两组 error 值（从第二个点开始）和对应的迭代次数（仅选择间隔内的点）
error_data_50_third2 <- data.frame(Iteration = seq(1 + interval, MAX_ITERS, by = interval), Error = error_values_Third2[seq(1 + interval, MAX_ITERS, by = interval)], Method = "GroupThird2")
error_data_50_gp2 <- data.frame(Iteration = seq(1 + interval, MAX_ITERS, by = interval), Error = error_values_gp2[seq(1 + interval, MAX_ITERS, by = interval)], Method = "GroupGP2")
error_data_50_soft <- data.frame(Iteration = seq(1 + interval, MAX_ITERS, by = interval), Error = error_values_soft[seq(1 + interval, MAX_ITERS, by = interval)], Method = "GroupSoft")
error_data_50_hard <- data.frame(Iteration = seq(1 + interval, MAX_ITERS, by = interval), Error = error_values_hard[seq(1 + interval, MAX_ITERS, by = interval)], Method = "GroupHard")
error_data_50_third <- data.frame(Iteration = seq(1 + interval, MAX_ITERS, by = interval), Error = error_values_Third[seq(1 + interval, MAX_ITERS, by = interval)], Method = "GroupThird")

# 将两组数据合并为一个数据框
error_data_combined <- rbind(error_data_50_third2, error_data_50_gp2, error_data_50_soft, error_data_50_hard, error_data_50_third)
# error_data_combined <- rbind(error_data_50_third2, error_data_50_gp2, error_data_50_soft, error_data_50_hard)

# error_plot_combined <- ggplot(data = error_data_combined, aes(x = Iteration, y = Error, group = Method, color = Method)) +
#   geom_line() +
#   geom_point(shape = 8, size = 3) +
#   labs(x = "Iteration", y = "Error", title = "Error vs. Iteration (Combined)") +
#   scale_color_manual(values = c("GroupThird2" = "red", "GroupGP2" = "blue"))

# # 打印绘制的图形
# print(error_plot_combined)

# 使用 ggplot2 绘制 error 随迭代次数变化的图形（从第二个点开始，每个方法用不同的颜色和点形状）
error_plot_combined <- ggplot(data = error_data_combined, aes(x = Iteration, y = Error, group = Method, color = Method, shape = Method)) +
  geom_line() +
  geom_point(size = 3) +
  labs(x = "Iteration", y = "Error", title = "Sparsity：10%") +
  scale_color_manual(values = c("GroupThird2" = "red", "GroupGP2" = "blue", "GroupSoft" = "green", "GroupHard" = "black", "GroupThird" = "orange")) +
  scale_shape_manual(values = c("GroupThird2" = 8, "GroupGP2" = 16, "GroupSoft" = 17, "GroupHard" = 18, "GroupThird" = 19))

# 打印绘制的图形
print(error_plot_combined)



















