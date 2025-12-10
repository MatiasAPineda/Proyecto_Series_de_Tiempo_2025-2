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
library(imputeTS)

## Data SINCA Estacíon Quintero 

# --- PM2,5 ---
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

pm25_vec <- PM25$valor_PM25
pm25_ts <- ts(pm25_vec, frequency = 365, start = start_ts)

## separación 90/10
n       <- length(pm25_ts)
n_train <- floor(0.9 * n)

y_train <- window(pm25_ts, end = start(pm25_ts) + c(0, n_train - 1))
y_test  <- window(pm25_ts, start = start(pm25_ts) + c(0, n_train))


## Box Cox
lambda <- BoxCox.lambda(y_train)
lambda

## ACF y PACF
acf(y_train, lag.max = 600, main = "ACF PM2.5")
pacf(y_train, lag.max = 600, main = "PACF PM2.5")

## Estimador de Dickey
ndiffs(y_train)

## Periodograma
LSTS::periodogram(y_train)

# Ajuste SARIMA
fit <- auto.arima(
  y_train,
  seasonal   = TRUE
)

summary(fit)

coefs <- coef(fit)
ses   <- sqrt(diag(fit$var.coef))
tvals <- abs(coefs / ses)
tvals
