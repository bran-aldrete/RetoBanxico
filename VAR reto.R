# =========================
# Paquetes
# =========================
library(readxl)
library(zoo)
library(lmtest)
library(sandwich)

# =========================
# 1) Leer tu archivo
# =========================
datos <- read_excel("C:/Users/brand/OneDrive/Desktop/ITAM/Banxico Reto/Predicciones/VAR/Datos Var.xlsx")

View(datos)
# Espera columnas: Trim, inf, tasa, Output gap

# =========================
# 2) Limpiar y parsear fecha Trim -> yearqtr
#    (corrige Dic -> Dec y quita espacios)
# =========================
inf<- ts(data=datos$inf,start=c(2006,1), frequency=4)
outgap<-ts(data=datos$`Output gap`,start=c(2006,1), frequency=4)
tiie<-ts(data=datos$tasa,start=c(2006,1), frequency=4)

df <- aggregate(
  cbind(inf, outgap, tiie) ~ time,
  data = df,
  FUN = mean,
  na.rm = TRUE
)

# Verificación rápida (debería ser TRUE)
if (!isTRUE(is.unique(df$time))) {
  stop("Aún hay trimestres duplicados en df$time. Revisa tu columna Trim.")
}

# =========================
# 4) Crear series zoo
# =========================
pi  <- zoo(df$inf,          df$time)  # inflación
r   <- zoo(df$tasa,         df$time)  # tasa
gap <- zoo(df$`Output gap`, df$time)  # output gap

# =========================
# 5) Rezagos
# =========================
L1 <- function(z) lag(z, -1)

pi_1  <- L1(pi)
r_1   <- L1(r)
gap_1 <- L1(gap)

# Dataset para regresiones
dat_z <- na.omit(merge(pi, r, gap, pi_1, r_1, gap_1))
colnames(dat_z) <- c("pi","r","gap","pi_1","r_1","gap_1")

# Pasa a data.frame (lm más estable así)
dat <- data.frame(
  time = index(dat_z),
  coredata(dat_z)
)

# =========================
# 6) Ecuación de inflación (Phillips reducida)
#    pi_t ~ gap_{t-1} + r_{t-1} + pi_{t-1}
# =========================
mod_pi <- lm(pi ~ gap_1 + r_1 + pi_1, data = dat)
print(coeftest(mod_pi, vcov = vcovHC(mod_pi, type = "HC1")))

# =========================
# 7) Ecuación del ciclo (IS reducida)
#    gap_t ~ pi_{t-1} + r_{t-1} + gap_{t-1}
# =========================
mod_gap <- lm(gap ~ pi_1 + r_1 + gap_1, data = dat)
print(coeftest(mod_gap, vcov = vcovHC(mod_gap, type = "HC1")))

# =========================
# 8) Calcular r* (neutral) usando steady state (pi*=3, gap=0)
#    0 = alpha_y + C1*pi* + C2*r* + D1*0
# =========================
b <- coef(mod_gap)
alpha_y <- b["(Intercept)"]
C1 <- b["pi_1"]
C2 <- b["r_1"]

pi_star <- 3
r_star <- - (alpha_y + C1*pi_star) / C2
cat("\nNeutral r* (aprox):", r_star, "\n")

# =========================
# 9) FORZAR phi_y y phi_pi (en vez de estimarlos)
# =========================
phi_y  <- 0.5
phi_pi <- 1.5

# Suavizamiento de tasa
use_smoothing <- TRUE
rho <- 0.8

# =========================
# 10) Simulación hacia adelante (baseline sin shocks)
# =========================
H <- 12  # horizonte (trimestres)

# últimos valores observados (t)
pi_t  <- as.numeric(tail(dat$pi, 1))
r_t   <- as.numeric(tail(dat$r, 1))
gap_t <- as.numeric(tail(dat$gap, 1))

bpi <- coef(mod_pi)
bg  <- coef(mod_gap)

path <- data.frame(h = 1:H, pi = NA_real_, gap = NA_real_, r = NA_real_)

for (h in 1:H) {
  
  # Pronóstico inflación
  pi_next <- bpi["(Intercept)"] +
    bpi["gap_1"] * gap_t +
    bpi["r_1"]   * r_t +
    bpi["pi_1"]  * pi_t
  
  # Pronóstico gap
  gap_next <- bg["(Intercept)"] +
    bg["pi_1"]  * pi_t +
    bg["r_1"]   * r_t +
    bg["gap_1"] * gap_t
  
  # Taylor "deseada" (sin shock)
  r_desired <- r_star + pi_next + phi_y * gap_next + phi_pi * (pi_next - pi_star)
  
  # Suavizamiento opcional
  r_next <- if (use_smoothing) rho * r_t + (1 - rho) * r_desired else r_desired
  
  # Guarda
  path$pi[h]  <- pi_next
  path$gap[h] <- gap_next
  path$r[h]   <- r_next
  
  # Actualiza estado
  pi_t  <- pi_next
  gap_t <- gap_next
  r_t   <- r_next
}

path

