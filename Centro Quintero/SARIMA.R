###################################
##      Proyecto - EYP3907       ##
## Vicente Garay - Matías Pineda ##
###################################

library(readr)
library(dplyr)
library(lubridate)
library(forecast)
library(LSTS)
library(MASS)
library(purrr)
library(car)
library(moments)
library(imputeTS)

## Data SINCA Estacíon Quintero 

# --- PM25 ---
PM25 <- read_delim(
  "DATA/PM25.csv",
  delim = ";",
  trim_ws = TRUE,
  col_types = cols(.default = "c"),
  show_col_types = FALSE
)

PM25 <- PM25 %>%
  mutate(
    fecha = ymd(sprintf("20%06s", `FECHA (YYMMDD)`)),
    val_valid = as.numeric(gsub(",", ".", `Registros validados`)),
    val_pre   = as.numeric(gsub(",", ".", `Registros preliminares`)),
    val_no    = as.numeric(gsub(",", ".", `Registros no validados`)),
    valor_PM25     = coalesce(val_valid, val_pre, val_no)
  ) %>%
  dplyr::select(fecha, valor_PM25)

summary(PM25)

na_vec <- is.na(PM25$valor_PM25)
rle_na <- rle(na_vec)
rle_na$lengths[rle_na$values == TRUE]
which(na_vec)

## --- Limpiar PM2.5, imputar y JUNTAR TODO ---
# Imputar NA puntuales restantes
PM25$imputado_pm <- is.na(PM25$valor_PM25)

PM25 <- PM25 %>%
  mutate(
    valor_PM25 = na_kalman(valor_PM25, model = "auto.arima")
  )

par(bty="n")
plot(
  PM25$fecha, PM25$valor_PM25, type = "l",
  main = "PM2.5 imputado",
  xlab = "Tiempo",
  ylab = expression("Material particulado 2.5 ( "*mu*"g/m"^3*")")
)
points(x = PM25$fecha[PM25$imputado_pm],
       y = PM25$valor_PM25[PM25$imputado_pm],
       col = "red", pch = 19, cex = 0.5)
legend("topleft",
       legend = c("Valor imputado"),
       col = "red",
       pch = 19,
       bty = "n")

start_year  <- as.numeric(format(min(PM25$fecha, na.rm = TRUE), "%Y"))
start_day   <- as.numeric(format(min(PM25$fecha, na.rm = TRUE), "%j"))
start_ts <- c(start_year, start_day)

PM25_vec <- PM25$valor_PM25 

## Periodograma
LSTS::periodogram(PM25_vec)
per <- LSTS::periodogram(PM25_vec, plot = FALSE)
I_lambda <- per$periodogram
lambda   <- per$lambda

ord <- order(I_lambda, decreasing = TRUE)
lambda_top <- lambda[ord[1:10]]
s_top     <- 2*pi / lambda_top
cbind(lambda_top, s_top)

## Pasar a ts
PM25_ts <- ts(PM25_vec, frequency = 14, start = start_ts)

## separación 90/10
n       <- length(PM25_ts)
n_train <- floor(0.9 * n)

y_train <- window(PM25_ts, end = start(PM25_ts) + c(0, n_train - 1))
y_test  <- window(PM25_ts, start = start(PM25_ts) + c(0, n_train))

## Box Cox
lambda_bc <- BoxCox.lambda(y_train)
lambda_bc

MASS::boxcox(y_train~1)

## ACF y PACF
acf(y_train, lag.max = 600, main = "ACF PM2.5")
pacf(y_train, lag.max = 600, main = "PACF PM2.5")

# Ajuste SARIMA
fit <- auto.arima(
  y_train,
  seasonal   = TRUE,
  lambda = lambda_bc,
  biasadj = TRUE
)

summary(fit)

coefs <- coef(fit)
ses   <- sqrt(diag(fit$var.coef))
tvals <- abs(coefs / ses)
tvals

fit_adj <- Arima(
  y_train,
  order = c(2,1,3),
  seasonal = c(0,0,1),
  lambda = lambda_bc,
  biasadj = TRUE
)

summary(fit_adj)

## Daignóstico
res <- fit$residuals
hist(res)
skewness(res)
kurtosis(res)

# Test
shapiro.test(res)

fit_t <- fitdistr(res, densfun = "t")
nu    <- fit_t$estimate["df"]
mu    <- fit_t$estimate["m"]
sigma <- fit_t$estimate["s"]
res_std <- (res - mu) / sigma
ks.test(res_std, "pt", df = nu)

LSTS::Box.Ljung.Test(res, lag=60)

acf(res,  lag.max=60)
acf(res^2,lag.max=60)

h <- length(y_test)

fc <- forecast(
  fit,
  h    = h
)

accuracy(fc, y_test)

pred <- as.numeric(fc$mean)
lo80 <- as.numeric(fc$lower[,1])
hi80 <- as.numeric(fc$upper[,1])
lo95 <- as.numeric(fc$lower[,2])
hi95 <- as.numeric(fc$upper[,2])

tt <- 1:length(y_test)
yr <- range(c(y_test, lo95, hi95), na.rm = TRUE)

par(bty = "n")
plot(tt, y_test, type = "l", ylim = yr, 
     xlab = "Tiempo (índice en test)", ylab = "PM2.5",
     main = "Validación SARIMA")

polygon(
  c(tt, rev(tt)),
  c(hi95, rev(lo95)),
  col = rgb(0, 0, 1, 0.15), border = NA
)
polygon(
  c(tt, rev(tt)),
  c(hi80, rev(lo80)),
  col = rgb(0, 0, 1, 0.30), border = NA
)

lines(tt, y_test, col = "black")
lines(tt, pred,  col = "red")

legend(
  "topright",
  legend = c("Predicción", "IC 80%", "IC 95%"),
  col    = c("red", rgb(0,0,1,0.30), rgb(0,0,1,0.15)),
  lwd    = c(2,10,10),
  bty    = "n"
)
