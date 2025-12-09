library(readxl)
library(forecast)
library(lmtest)
library(LSTS)
library(ggplot2)
library(MASS)

df <- read_excel("stock de cartera consumo normal.xlsx", skip = 3)

df$Fecha <- as.Date(df$Fecha)

df$`Banco de Chile` <- as.numeric(gsub(",", ".", df$`Banco del Estado de Chile`))

serie <- ts(df$`Banco del Estado de Chile`, 
                   start = c(2011, 1),
                   end = c(2025, 5),
                   frequency = 12)

n       <- length(serie)
n_train <- floor(n * 0.9)

serie_train <- subset(serie, end   = n_train)
serie_test  <- subset(serie, start = n_train + 1)


#----------------------------------------------------------
ndiffs(serie_train)

acf(serie_train, lag.max = 30)
pacf(serie_train, lag.max = 30)

modelo <- auto.arima(serie_train, 
                      seasonal = TRUE)
summary(modelo)

run_diagnostics <- function(fit, name, alpha = 0.05) {
  cat("\n==============================\n")
  cat("DIAGNÓSTICOS PARA:", name, "\n")
  cat("==============================\n")
  
  res <- residuals(fit)
  res <- as.numeric(res)
  t   <- 1:length(res)
  
  # Box-Ljung
  cat("\n--- Box-Ljung ---\n")
  bj <- LSTS::Box.Ljung.Test(res, lag = 30)
  print(bj)
  
  # Breusch-Pagan
  cat("\n--- Breusch-Pagan ---\n")
  bp <- lmtest::bptest(res ~ t)
  print(bp)
  
  # Shapiro-Wilk
  cat("\n--- Shapiro-Wilk ---\n")
  # Shapiro no acepta n muy grande; si tu serie es enorme, puedes muestrear:
  if (length(res) <= 5000) {
    sw <- shapiro.test(res)
  } else {
    sw <- shapiro.test(sample(res, 5000))
  }
  print(sw)
  
  # KS vs Normal
  cat("\n--- Kolmogorov-Smirnov ---\n")
  ks <- ks.test(
    res,
    "pnorm",
    mean = mean(res, na.rm = TRUE),
    sd   = sd(res, na.rm = TRUE)
  )
  print(ks)
  
  invisible(list(
    BoxLjung = bj,
    BP       = bp,
    Shapiro  = sw,
    KS       = ks
  ))
}

diag_M1 <- run_diagnostics(modelo, "mati weko")

acf(modelo$residuals^2)

res <- residuals(modelo)
res <- as.numeric(res)

fit_t <- fitdistr(res, densfun = "t")
nu    <- fit_t$estimate["df"]
mu    <- fit_t$estimate["m"]
sigma <- fit_t$estimate["s"]
res_std <- (res - mu) / sigma
ks.test(res_std, "pt", df = nu)
