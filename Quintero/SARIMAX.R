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

summary(PM25)

na_vec <- is.na(PM25$valor)
rle_na <- rle(na_vec)
rle_na$lengths[rle_na$values == TRUE]
which(na_vec)

# --- VIENTO (WS) ---

WS <- read_delim(
  "WS.csv",
  delim = ";",
  trim_ws = TRUE,
  col_types = cols(.default = "c")
)

names(WS) <- c("fecha_raw", "hora_raw", "WS_raw", "col4")

WS <- WS %>%
  mutate(
    fecha = ymd(sprintf("20%06s", fecha_raw)),
    hora  = sprintf("%04s", hora_raw),
    WS    = as.numeric(gsub(",", ".", WS_raw))
  ) %>%
  dplyr::select(fecha, WS)

WS_diario <- WS %>%
  group_by(fecha) %>%
  summarise(WS = mean(WS, na.rm = TRUE))

summary(WS_diario)

## --- NOX ---
NOX <- read_delim(
  "NOX.csv",
  delim = ";",
  col_types = cols(.default = "c"),  # TODO como texto
  trim_ws = TRUE
)

NOX <- NOX %>%
  mutate(
    fecha = ymd(sprintf("20%06s", `FECHA (YYMMDD)`)),
    hora  = sprintf("%04s", `HORA (HHMM)`),
    val_valid = as.numeric(gsub(",", ".", `Registros validados`)),
    val_pre   = as.numeric(gsub(",", ".", `Registros preliminares`)),
    val_no    = as.numeric(gsub(",", ".", `Registros no validados`)),
    NOX_val   = dplyr::coalesce(val_valid, val_pre, val_no)
  ) %>%
  dplyr::select(fecha, NOX_val)

summary(NOX)

na_vec <- is.na(NOX$NOX_val)
rle_na <- rle(na_vec)
rle_na$lengths[rle_na$values == TRUE]
which(na_vec)

## --- Dias de Semana ---
dow <- wday(PM25$fecha, week_start = 1)

dow_fac <- factor(
  dow,
  levels = 1:7,
  labels = c("Lun", "Mar", "Mie", "Jue", "Vie", "Sab", "Dom")
)

dow_mat <- model.matrix(~ dow_fac)[, -1]   # referencia = "Lun"
colnames(dow_mat) <- c("Mar", "Mie", "Jue", "Vie", "Sab", "Dom")


## --- Limpiar PM2.5, imputar y JUNTAR TODO ---

# Eliminar primer periodo largo de NA y quedarnos desde el primer dato válido
idx <- which(!is.na(PM25$valor))[1]
PM25 <- PM25 %>% slice(idx:n())

# Imputar NA puntuales restantes
PM25$imputado_pm <- is.na(PM25$valor)
NOX$imputado_nox <- is.na(NOX$NOX_val)

PM25 <- PM25 %>%
  mutate(
    valor = na_kalman(valor, model = "auto.arima")
  )

NOX <- NOX %>%
  mutate(
    NOX_val = na_kalman(NOX_val, model = "auto.arima")
  )

# Juntar PM2.5 con WS diario y NOX diario
datos <- PM25 %>%
  left_join(WS_diario,  by = "fecha") %>%
  left_join(NOX, by = "fecha")

# Rezagos 1 de WS y NOX
datos <- datos %>%
  mutate(
    WS_lag1   = dplyr::lag(WS, 1),
    NOX_lag1  = dplyr::lag(NOX_val, 1)
  )

# Revisar
summary(datos$WS_lag1)
summary(datos$NOX_lag1)

# Agregar estaciones
datos <- datos %>%
  mutate(
    mes = month(fecha),
    invierno = if_else(mes %in% c(6,7,8), 1, 0),   # JJA
    verano   = if_else(mes %in% c(12,1,2), 1, 0)   # opcional
  )

summary(datos[, c("valor","WS_lag1","NOX_lag1")])

## --- Gráficos ---
par(bty="n")
plot(datos$fecha, datos$WS_lag1, type = "l", main = "Velocidad del Viento (lag 1)",
     ylab = "Velocidad (m/s)", xlab = "Time")

par(bty="n")
plot(datos$fecha, datos$NOX_lag1, type = "l", main = "Óxidos de Nitrógeno (lag 1)",
     ylab = "Óxidos de Nitrógeno (NOX) (ppb)", xlab = "Time")
points(
  x = datos$fecha[datos$imputado_nox],
  y = datos$NOX_lag1[datos$imputado_nox],
  col = "red", pch = 19, cex = 0.5
)
legend("topleft",
       legend = c("Valor imputado"),
       col = "red",
       pch = 19,
       bty = "n")

par(bty="n")
plot(
  datos$fecha, datos$valor, type = "l",
  main = "PM2.5 imputado",
  xlab = "Tiempo",
  ylab = expression("Material particulado 2.5 ( "*mu*"g/m"^3*")")
)
points(x = datos$fecha[datos$imputado_pm],
       y = datos$valor[datos$imputado_pm],
       col = "red", pch = 19, cex = 0.5)
legend("topleft",
       legend = c("Valor imputado"),
       col = "red",
       pch = 19,
       bty = "n")

## --- Serie de tiempo con frecuencia semanal ---

primer_dia   <- datos$fecha[1]
anio_inicio  <- year(primer_dia)
dia_semana   <- wday(primer_dia) 

freq <- 365

pm25_ts <- ts(datos$valor,   start = c(anio_inicio, 1), frequency = freq)
ws_ts   <- ts(datos$WS_lag1, start = c(anio_inicio, 1), frequency = freq)
nox_ts  <- ts(datos$NOX_lag1, start = c(anio_inicio, 1), frequency = freq)
nox2_ts  <- ts(datos$NOX_lag1^2, start = c(anio_inicio, 1), frequency = freq)
dow_ts <- ts(dow_mat, start = start(pm25_ts), frequency = freq)
invierno_ts <- ts(datos$invierno, start=start(pm25_ts), frequency=freq)

# separación 90/10
n       <- length(pm25_ts)
n_train <- floor(0.9 * n)

end_train  <- time(pm25_ts)[n_train]
start_test <- time(pm25_ts)[n_train + 1]

y_train <- window(pm25_ts, end   = end_train)
y_test  <- window(pm25_ts, start = start_test)

# Matriz de regresores (WS_lag1 y NOX_lag1)
X_ts <- cbind(
  ws_ts,
  nox_ts,             # si lo sigues usando
  nox2_ts,            # NOX^2
  dow_ts,
  invierno_ts          # ts(datos$invierno, start=..., freq=365)
)
x_train <- window(X_ts, end   = end_train)
x_test  <- window(X_ts, start = start_test)

## --- Modelación ---

# Box Cox
lambda <- BoxCox.lambda(y_train)
lambda

# ACF y PACF (opcional)
par(mfrow = c(1,2))
acf(y_train, lag.max = 56, main = "ACF PM2.5 (freq=7)")
pacf(y_train, lag.max = 56, main = "PACF PM2.5 (freq=7)")
par(mfrow = c(1,1))

LSTS::periodogram(y_train)

# Ajuste SARIMAX
fit_x <- auto.arima(
  y_train,
  xreg       = x_train,
  seasonal   = TRUE,
  lambda     = lambda,
  biasadj    = TRUE
)

summary(fit_x)

# Significancia
coefs <- coef(fit_x)
ses   <- sqrt(diag(fit_x$var.coef))
tvals <- abs(coefs / ses)
tvals

# Fijamos en 0 solo los coeficientes no significativos
fixed_vec <- rep(0, length(coef(fit_x)))
names(fixed_vec) <- names(coef(fit_x))

fixed_vec["ar1"] <- NA
fixed_vec["ar2"] <- NA
fixed_vec["ma2"] <- NA
fixed_vec["ws_ts"] <- NA
fixed_vec["nox2_ts"] <- NA
fixed_vec["dow_ts.Vie"] <- NA
fixed_vec["dow_ts.Sab"] <- NA

fit7c <-  Arima(
  y_train,
  order    = c(5,1,2), 
  xreg     = x_train,
  lambda   = lambda,
  biasadj  = TRUE,
  fixed    = fixed_vec
)

summary(fit7c)

## --- Diagnóstico (sobre fit7c, univariado) ---

res <- na.omit(fit7c$residuals)
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
Box.test(res^2, lag = 16, type = "Ljung")

## --- Predicción con SARIMAX (fit_x) ---

h <- length(y_test)

fc_x <- forecast(
  fit7c,
  xreg = x_test,
  h    = h
)

accuracy(fc_x, y_test)

# Plot test vs predicción SARIMAX (base R)

pred <- as.numeric(fc_x$mean)
lo80 <- as.numeric(fc_x$lower[,1])
hi80 <- as.numeric(fc_x$upper[,1])
lo95 <- as.numeric(fc_x$lower[,2])
hi95 <- as.numeric(fc_x$upper[,2])

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

# GARCH a los residuos

res <- residuals(fit_x)
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

fc_test <- forecast(fit_x, xreg = x_test, h = length(y_test))
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
