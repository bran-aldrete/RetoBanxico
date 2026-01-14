# =========================
# Paquetes
# =========================
library(readxl)
library(zoo)
library(lmtest)
library(sandwich)
library(ggplot2)
library(tidyr)

# Configurar httpgd para visualización en VS Code
if (requireNamespace("httpgd", quietly = TRUE)) {
  options(vsc.plot = FALSE)
  options(device = function(...) {
    httpgd::hgd(width = 10, height = 7.5, silent = TRUE)
    .Call("hgd_plot_activate", PACKAGE = "httpgd")
  })
}

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

# =========================
# 4) Crear series zoo
# =========================
pi  <- zoo(as.numeric(inf),     time(inf))    # inflación
r   <- zoo(as.numeric(tiie),    time(tiie))   # tasa
gap <- zoo(as.numeric(outgap),  time(outgap)) # output gap

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

# Exporta el pronóstico a CSV
write.csv(
  path,
  "C:/Users/brand/OneDrive/Desktop/ITAM/Banxico Reto/Predicciones/VAR/path_forecast.csv",
  row.names = FALSE
)

# Gráfica estilizada (ggplot2) y exportada a PNG
path_long <- pivot_longer(path, c(pi, gap, r), names_to = "serie", values_to = "valor")

p <- ggplot(path_long, aes(x = h, y = valor, color = serie)) +
  geom_hline(yintercept = 3, linetype = "dashed", color = "gray50", size = 0.7) +
  geom_line(size = 1.1) +
  scale_color_manual(
    values = c(pi = "firebrick", gap = "steelblue", r = "darkgreen"),
    labels = c(pi = "Inflación (pi)", gap = "Output gap", r = "Tasa (r)")
  ) +
  labs(
    title = "Trayectorias simuladas",
    subtitle = "Incluye meta de inflación 3%",
    x = "Horizonte (trimestres)",
    y = "Nivel / %",
    color = "Serie"
  ) +
  annotate("text", x = max(path$h) * 0.98, y = 3, label = "Meta 3%", hjust = 1, vjust = -0.5, color = "gray30", size = 3.5) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.major = element_line(color = "gray85"),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(face = "bold")
  )

ggsave(
  filename = "C:/Users/brand/OneDrive/Desktop/ITAM/Banxico Reto/Predicciones/VAR/path_forecast.png",
  plot = p,
  width = 12,
  height = 9,
  dpi = 200
)

print(p)

