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
library(arfima)
library(pracma)
library(fracdiff)

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

## --- Exploración de los Datos ---

# Pasar a TS
primer_dia   <- PM25$fecha[1]
anio_inicio  <- year(primer_dia)

pm25_ts7 <- ts(PM25$valor, start = c(anio_inicio, 1), frequency = 365)

plot(diff(pm25_ts7, lag=365))

# Separar Train/Test

n       <- length(pm25_ts7)
n_train <- floor(0.9 * n)

end_train  <- time(pm25_ts7)[n_train]
start_test <- time(pm25_ts7)[n_train + 1]

y_train <- window(pm25_ts7, end   = end_train)
y_test  <- window(pm25_ts7, start = start_test)

# Box Cox

lambda_bc <- BoxCox.lambda(y_train)
lambda_bc

# Periodograma y frecuencias grandes
LSTS::periodogram(y_train)
per <- LSTS::periodogram(y_train, plot = FALSE)
I_lambda <- per$periodogram
lambda   <- per$lambda

ord <- order(I_lambda, decreasing = TRUE)
lambda_top <- lambda[ord[1:10]]
s_top     <- 2*pi / lambda_top
cbind(lambda_top, s_top)

# ACF y PACF
acf(y_train, lag.max = 600, main = "ACF PM2.5")
pacf(y_train, lag.max = 600, main = "PACF PM2.5")

# A partir del Periodograma y del ACF hay sospecha de Larga Memoria
diff(y_train)

hurstexp(y_train)

# Tenemos que d=0.44 muy significativo
plot(decompose(y_train)) 

# Todavía hay estacionalidad
K <- 3  
X_train <- fourier(y_train, K = K)

fit_fourier_arfima <- arfima(
  z     = y_train,
  order = c(1, 0, 0),   
  xreg  = X_train             
)

summary(fit_fourier_arfima)



