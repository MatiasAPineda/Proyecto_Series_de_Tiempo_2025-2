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

## --- CO ---
CO <- read_delim(
  "DATA/CO.csv",
  delim = ";",
  trim_ws = TRUE,
  col_types = cols(.default = "c"),
  show_col_types = FALSE
)

CO <- CO %>%
  mutate(
    fecha = ymd(sprintf("20%06s", `FECHA (YYMMDD)`)),
    val_valid = as.numeric(gsub(",", ".", `Registros validados`)),
    val_pre   = as.numeric(gsub(",", ".", `Registros preliminares`)),
    val_no    = as.numeric(gsub(",", ".", `Registros no validados`)),
    valor_CO     = coalesce(val_valid, val_pre, val_no)
  ) %>%
  dplyr::select(fecha, valor_CO)

summary(CO)

## --- NO ---
NO <- read_delim(
  "DATA/NO.csv",
  delim = ";",
  trim_ws = TRUE,
  col_types = cols(.default = "c"),
  show_col_types = FALSE
)

NO <- NO %>%
  mutate(
    fecha = ymd(sprintf("20%06s", `FECHA (YYMMDD)`)),
    val_valid = as.numeric(gsub(",", ".", `Registros validados`)),
    val_pre   = as.numeric(gsub(",", ".", `Registros preliminares`)),
    val_no    = as.numeric(gsub(",", ".", `Registros no validados`)),
    valor_NO     = coalesce(val_valid, val_pre, val_no)
  ) %>%
  dplyr::select(fecha, valor_NO)

summary(NO)

## --- NO2 ---
NO2 <- read_delim(
  "DATA/NO2.csv",
  delim = ";",
  trim_ws = TRUE,
  col_types = cols(.default = "c"),
  show_col_types = FALSE
)

NO2 <- NO2 %>%
  mutate(
    fecha = ymd(sprintf("20%06s", `FECHA (YYMMDD)`)),
    val_valid = as.numeric(gsub(",", ".", `Registros validados`)),
    val_pre   = as.numeric(gsub(",", ".", `Registros preliminares`)),
    val_no    = as.numeric(gsub(",", ".", `Registros no validados`)),
    valor_NO2    = coalesce(val_valid, val_pre, val_no)
  ) %>%
  dplyr::select(fecha, valor_NO2)

summary(NO2)

## --- NOX ---
NOX <- read_delim(
  "DATA/NOX.csv",
  delim = ";",
  trim_ws = TRUE,
  col_types = cols(.default = "c"),
  show_col_types = FALSE
)

NOX <- NOX %>%
  mutate(
    fecha = ymd(sprintf("20%06s", `FECHA (YYMMDD)`)),
    val_valid = as.numeric(gsub(",", ".", `Registros validados`)),
    val_pre   = as.numeric(gsub(",", ".", `Registros preliminares`)),
    val_no    = as.numeric(gsub(",", ".", `Registros no validados`)),
    valor_NOX     = coalesce(val_valid, val_pre, val_no)
  ) %>%
  dplyr::select(fecha, valor_NOX)

summary(NOX)

## --- O3 ---
O3 <- read_delim(
  "DATA/O3.csv",
  delim = ";",
  trim_ws = TRUE,
  col_types = cols(.default = "c"),
  show_col_types = FALSE
)

O3 <- O3 %>%
  mutate(
    fecha = ymd(sprintf("20%06s", `FECHA (YYMMDD)`)),
    val_valid = as.numeric(gsub(",", ".", `Registros validados`)),
    val_pre   = as.numeric(gsub(",", ".", `Registros preliminares`)),
    val_no    = as.numeric(gsub(",", ".", `Registros no validados`)),
    valor_O3     = coalesce(val_valid, val_pre, val_no)
  ) %>%
  dplyr::select(fecha, valor_O3)

summary(O3)

## --- PM10 ---
PM10 <- read_delim(
  "DATA/PM10.csv",
  delim = ";",
  trim_ws = TRUE,
  col_types = cols(.default = "c"),
  show_col_types = FALSE
)

PM10 <- PM10 %>%
  mutate(
    fecha = ymd(sprintf("20%06s", `FECHA (YYMMDD)`)),
    val_valid = as.numeric(gsub(",", ".", `Registros validados`)),
    val_pre   = as.numeric(gsub(",", ".", `Registros preliminares`)),
    val_no    = as.numeric(gsub(",", ".", `Registros no validados`)),
    valor_PM10     = coalesce(val_valid, val_pre, val_no)
  ) %>%
  dplyr::select(fecha, valor_PM10)

summary(PM10)

## --- SO2 ---
SO2 <- read_delim(
  "DATA/SO2.csv",
  delim = ";",
  trim_ws = TRUE,
  col_types = cols(.default = "c"),
  show_col_types = FALSE
)

SO2 <- SO2 %>%
  mutate(
    fecha = ymd(sprintf("20%06s", `FECHA (YYMMDD)`)),
    val_valid = as.numeric(gsub(",", ".", `Registros validados`)),
    val_pre   = as.numeric(gsub(",", ".", `Registros preliminares`)),
    val_no    = as.numeric(gsub(",", ".", `Registros no validados`)),
    valor_SO2     = coalesce(val_valid, val_pre, val_no)
  ) %>%
  dplyr::select(fecha, valor_SO2)

summary(SO2)

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

# Juntar datos
datos <- PM25 %>%
  left_join(CO,  by = "fecha") %>%
  left_join(NO, by = "fecha") %>%
  left_join(NO2, by = "fecha") %>%
  left_join(NOX, by = "fecha") %>%
  left_join(O3, by = "fecha") %>%
  left_join(PM10, by = "fecha") %>%
  left_join(SO2, by = "fecha")

head(datos)
summary(datos)

for (columna in colnames(datos)) {
  print(columna)
  na_vec <- is.na(datos[columna])
  rle_na <- rle(as.vector(na_vec))
  print(rle_na$lengths[rle_na$values == TRUE])
}

# Imputar datos faltantes
datos <- datos %>%
  mutate(
    valor_CO = na_kalman(valor_CO, model = "auto.arima"),
    valor_NO = na_kalman(valor_NO, model = "auto.arima"),
    valor_NO2 = na_kalman(valor_NO2, model = "auto.arima"),
    valor_NOX = na_kalman(valor_NOX, model = "auto.arima"),
    valor_O3 = na_kalman(valor_O3, model = "auto.arima"),
    valor_PM10 = na_kalman(valor_PM10, model = "auto.arima"),
    valor_SO2 = na_kalman(valor_SO2, model = "auto.arima")
  )

num_data <- datos %>% dplyr::select(where(is.numeric))

## LAg a todo meno PM25
num_data %>%
  mutate(across(.cols = -all_of("valor_PM25"), ~ dplyr::lag(.x, n = 1)))

## Correlaciones entre datos
library(corrplot)

corrplot(cor(num_data, use = "pairwise.complete.obs"),
         method = "color",
         type = "upper",
         tl.cex = 0.8)

## Pasar a ts
start_year  <- as.numeric(format(min(datos$fecha, na.rm = TRUE), "%Y"))
start_day   <- as.numeric(format(min(datos$fecha, na.rm = TRUE), "%j"))
start_ts <- c(start_year, start_day)

pm25_vec <- num_data[["valor_PM25"]]
pm25_ts <- ts(pm25_vec, frequency = 14, start = start_ts)

X_df <- num_data %>% dplyr::select(-all_of("valor_PM25"))
X_mat <- as.matrix(map_dfc(X_df, as.numeric))

X_ts <- ts(X_mat, frequency = 14, start = start_ts)

## separación 90/10
n       <- length(pm25_ts)
n_train <- floor(0.9 * n)

y_train <- window(pm25_ts, end = start(pm25_ts) + c(0, n_train - 1))
y_test  <- window(pm25_ts, start = start(pm25_ts) + c(0, n_train))

X_train <- window(X_ts, end = start(X_ts) + c(0, n_train - 1))
X_test  <- window(X_ts, start = start(X_ts) + c(0, n_train))

## Modelacion

## Box Cox
lambda <- BoxCox.lambda(y_train)
lambda

MASS::boxcox(y_train~1)

## ACF y PACF
acf(y_train, lag.max = 600, main = "ACF PM2.5")
pacf(y_train, lag.max = 600, main = "PACF PM2.5")

## Estimador de Dickey
ndiffs(y_train)

## Periodograma
LSTS::periodogram(y_train)

# Ajuste SARIMAX 
bruteforce_sarimax <- function(
    y_train, X_train, order, seasonal=NULL, lambda=NULL, verbose=FALSE
) {
  
  # Convertir X a matriz
  X <- as.matrix(X_train)
  p <- ncol(X)
  
  # Todas las combinaciones posibles: 0..255
  combos <- lapply(0:(2^p - 1), function(mask){
    which(as.logical(intToBits(mask)[1:p]))
  })
  
  results <- vector("list", length(combos))
  n_models <- length(combos)
  
  for(i in seq_along(combos)){
    idx <- combos[[i]]
    xreg <- if(length(idx) == 0) NULL else X[, idx, drop = FALSE]
    
    reg_names <- if(length(idx)==0) {
      "ninguno"
    } else {
      paste0(colnames(X)[idx], collapse=", ")
    }
    
    if (verbose) {
      cat("\n==========================================================\n")
      cat(sprintf("📌 Modelo %d / %d\n", i, n_models))
      cat(sprintf("   Regresores: %s\n", reg_names))
      cat("==========================================================\n")
    }
    
    # Ajuste SARIMA/SARIMAX
    fit <- tryCatch({
      Arima(
        y_train,
        order    = order,
        seasonal = seasonal,
        xreg     = xreg,
        lambda   = lambda,
        biasadj = TRUE
      )
    }, error = function(e){
      if(verbose){
        cat("❌ Error al ajustar modelo:", e$message, "\n")
      }
      NULL
    })
    
    if(is.null(fit)){
      # Guardar fallo
      results[[i]] <- list(
        comb_id   = i,
        k_regs    = length(idx),
        regs      = paste(idx, collapse=","),
        AIC       = Inf,
        BIC       = Inf,
        loglik    = NA,
        converged = FALSE
      )
    } else {
      # Guardar modelo OK
      results[[i]] <- list(
        comb_id   = i,
        k_regs    = length(idx),
        regs      = paste(idx, collapse=","),
        AIC       = AIC(fit),
        BIC       = BIC(fit),
        loglik    = as.numeric(logLik(fit)),
        converged = TRUE,
        model     = fit
      )
      
      # Mostrar summary
      if(verbose){
        cat("\n📄 Summary del modelo ajustado:\n")
        print(summary(fit))
      }
    }
  }
  
  # Convertir a tabla ordenada por BIC
  df <- dplyr::bind_rows(results) %>% dplyr::arrange(BIC)
  
  list(
    table  = df,
    models = results
  )
}

res_brute <- bruteforce_sarimax(
  y_train = y_train,
  X_train = X_train,
  order   = c(2,1,3),
  seasonal = c(0,0,1),  
  lambda  = -0.1149271 ,
  verbose = TRUE
)

head(res_brute$table, 5)

fit <- resultado_forward$model

coefs <- coef(fit)
ses   <- sqrt(diag(fit$var.coef))
tvals <- abs(coefs / ses)
tvals

fixed_vec <- rep(NA, length(coef(fit)))
names(fixed_vec) <- names(coef(fit))

fixed_vec["valor_NO2"] <- 0
fixed_vec["valor_NOX"] <- 0

selected <- resultado_forward$selected_regressors
xreg <- as.matrix(X_train[, selected, drop = FALSE])

fit <-  Arima(
  y_train,
  order    = c(3,1,1), 
  xreg     = xreg,
  lambda   = lambda,
  biasadj  = TRUE,
  fixed    = fixed_vec
)
summary(fit)

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
selected <- resultado_forward$selected_regressors
xreg_forecast <- as.matrix(X_test[, selected, drop = FALSE])

fc <- forecast(fit, xreg = xreg_forecast, h = h)

pred_full <- as.numeric(fc$mean)
y_test_vec <- as.numeric(y_test)

accuracy(pred_full, y_test_vec)

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

plot(tail(BoxCox(y_train, -0.1149271), h))

tail(datos, 1)
