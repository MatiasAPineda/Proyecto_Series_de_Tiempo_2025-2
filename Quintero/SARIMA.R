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
library(car)
library(imputeTS)

## Data SINCA Estacíon Quintero 

# Preparar datos PM2.5
PM25 <- read_delim(
  "PM25.csv",
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
    valor     = coalesce(val_valid, val_pre, val_no)
  ) %>%
  dplyr::select(fecha, valor)

par(bty="n")
plot(PM25$fecha, PM25$valor, type="l", ylab = "PM2.5", xlab = "Time")

summary(PM25)

na_vec <- is.na(PM25$valor)
rle_na <- rle(na_vec)
rle_na$lengths[rle_na$values == TRUE]
which(na_vec)

## --- Limpiar PM2.5, imputar y JUNTAR TODO ---

# Eliminar primer periodo largo de NA y quedarnos desde el primer dato válido
idx <- which(!is.na(PM25$valor))[1]
PM25 <- PM25 %>% slice(idx:n())

# Imputar NA puntuales restantes en PM2.5
PM25$imputado <- is.na(PM25$valor)

PM25 <- PM25 %>%
  mutate(
    valor = na_kalman(valor, model = "auto.arima")
  )

par(bty="n")
plot(PM25$fecha, PM25$valor, type = "l", main = "PM2.5 imputado",
     ylab = "PM2.5", xlab = "Time")
points(x = PM25$fecha[which(PM25$imputado)],     
       y = PM25$valor[PM25$imputado],  
       col = "red", pch = 19, cex = 0.5)
legend("topleft",
       legend = c("Valor imputado"),
       col = "red",
       pch = 19,
       bty = "n")


primer_dia   <- PM25$fecha[1]
anio_inicio  <- year(primer_dia)

pm25_ts7 <- ts(PM25$valor,   start = c(anio_inicio, 1), frequency = 365)

n       <- length(pm25_ts7)
n_train <- floor(0.9 * n)

end_train  <- time(pm25_ts7)[n_train]
start_test <- time(pm25_ts7)[n_train + 1]

y_train <- window(pm25_ts7, end   = end_train)
y_test  <- window(pm25_ts7, start = start_test)

lambda_bc <- BoxCox.lambda(y_train)
lambda_bc

LSTS::periodogram(y_train)
per <- LSTS::periodogram(y_train, plot = FALSE)
I_lambda <- per$periodogram
lambda   <- per$lambda

ord <- order(I_lambda, decreasing = TRUE)
lambda_top <- lambda[ord[1:10]]
s_top     <- 2*pi / lambda_top
cbind(lambda_top, s_top)

par(mfrow = c(1,2))
acf(y_train, lag.max = 600, main = "ACF PM2.5 (freq=7)")
pacf(y_train, lag.max = 600, main = "PACF PM2.5 (freq=7)")
par(mfrow = c(1,1))

fit <- auto.arima(
  y_train,
  seasonal   = TRUE,
  lambda     = lambda_bc,
  biasadj = TRUE
)

summary(fit)

coefs <- coef(fit)
ses   <- sqrt(diag(fit$var.coef))
tvals <- abs(coefs / ses)
tvals

fixed_params <- c(
  NA,   # ar1
  0,   # ar2
  NA, # ar3
  0,   # ma1
  NA,    # ma2
  NA, # ma3
  0 # ma4
)

fit7c <- Arima(
  y_train,
  order = c(3,1,4),
  fixed = fixed_params,
  lambda = lambda_bc,
  biasadj = TRUE
)

summary(fit7c)

res <- fit7c$residuals
hist(res)
shapiro.test(res)

fit_t <- fitdistr(res, densfun = "t")
nu    <- fit_t$estimate["df"]
mu    <- fit_t$estimate["m"]
sigma <- fit_t$estimate["s"]
res_std <- (res - mu) / sigma
ks.test(res_std, "pt", df = nu)

LSTS::Box.Ljung.Test(res, lag=30)

acf(res,  lag.max=30)
acf(res^2,lag.max=30)

Box.test(res^2, lag = 1,  type = "Ljung")
Box.test(res^2, lag = 17, type = "Ljung")

## --- Predicción con SARIMAX (fit_x) ---

h <- length(y_test)

fc <- forecast(
  fit,
  h    = h
)

accuracy(fc, y_test)

# Plot test vs predicción SARIMAX (base R)

pred <- as.numeric(fc$mean)
lo80 <- as.numeric(fc$lower[,1])
hi80 <- as.numeric(fc$upper[,1])
lo95 <- as.numeric(fc$lower[,2])
hi95 <- as.numeric(fc$upper[,2])

tt <- 1:length(y_test)
yr <- range(c(y_test, lo95, hi95), na.rm = TRUE)

par(bty = "n")
plot(tt, y_test, type = "l", lwd = 1,
     ylim = yr, xlab = "Tiempo (índice en test)", ylab = "PM2.5",
     main = "Validación SARIMAX (WS_lag1 + NOX_lag1)")

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

lines(tt, y_test, col = "black", lwd = 1)
lines(tt, pred,  col = "red",   lwd = 1)

legend(
  "topleft",
  legend = c("Predicción", "IC 80%", "IC 95%"),
  col    = c("red", rgb(0,0,1,0.30), rgb(0,0,1,0.15)),
  lwd    = c(2,10,10),
  bty    = "n"
)

res <- residuals(fit7c)
res <- as.numeric(res)
res_train <- res[!is.na(res)]

library(FinTS)
ArchTest(res_train, lags = 12)

library(rugarch)
spec_eg <- ugarchspec(
  variance.model = list(model = "eGARCH", garchOrder = c(1,2)),
  mean.model = list(armaOrder = c(0,0), include.mean = FALSE),
  distribution.model = "std"
)
fit_eg <- ugarchfit(spec_eg, res_train)
show(fit_eg)

# Residuos estandarizados del modelo GARCH
res_std <- residuals(fit_eg, standardize = TRUE)

acf(res_std, main = "ACF residuos GARCH")
acf(res_std^2, main = "ACF residuos^2 GARCH")

Box.test(res_std,   lag = 12, type = "Ljung")
Box.test(res_std^2, lag = 12, type = "Ljung")

fc_test <- forecast(fit7c, h = length(y_test))
mu_test <- as.numeric(fc_test$mean)

h <- length(y_test)

fc_vol <- ugarchforecast(
  fit_eg, 
  n.ahead = h
)

sigma_test <- sigma(fc_vol) 

nu <- coef(fit_eg)["shape"]   # grados de libertad de la t
crit80 <- qt(0.9, df = nu)    # corte para 80%
crit95 <- qt(0.975, df = nu)  # corte para 95%

lo80 <- mu_test - crit80 * sigma_test
hi80 <- mu_test + crit80 * sigma_test

lo95 <- mu_test - crit95 * sigma_test
hi95 <- mu_test + crit95 * sigma_test

pred_test <- mu_test
accuracy(pred_test, y_test)
mean(y_test >= lo95 & y_test <= hi95)

tt <- 1:length(y_test)
yr <- range(c(y_test, lo95, hi95), na.rm = TRUE)

par(bty="n")
plot(tt, y_test, type="l", col="black", lwd=1,
     ylim=yr, xlab="Tiempo", 
     ylab=expression("PM2.5 ("*mu*"g/m"^3*")"),
     main="Validación híbrida SARIMA + EGARCH")

polygon(c(tt, rev(tt)), c(hi95, rev(lo95)),
        col = rgb(0,0,1,0.15), border=NA)

polygon(c(tt, rev(tt)), c(hi80, rev(lo80)),
        col = rgb(0,0,1,0.30), border=NA)

lines(tt, pred_test, col="red", lwd=1)
lines(tt, y_test, col="black", lwd=1)

legend("topleft",
       legend = c("Observado", "Predicción SARIMA", "IC80% híbrido", "IC95% híbrido"),
       col = c("black","red", rgb(0,0,1,0.30), rgb(0,0,1,0.15)),
       lwd = c(1,2,10,10),
       bty="n")

