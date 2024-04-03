rm(list=ls())


### Tablas de Mortalidad

# Definición de las funciones para leer las tablas de mortalidad
tabla.PER2020 <- function() {
  read.csv("tablas/PrimerPER2020.csv")
}

tabla.PASEM2020 <- function() {
  read.csv("tablas/PrimerPASEM2020.csv")
}

# Funciones para calcular las tablas de mortalidad por sexo y el ajuste generacional
PER2020.q <- function(anyo = 2012, tantoPOR = 1, prop = 0.5) {
  PER2020 <- tabla.PER2020()
  qbf <- PER2020$qxf * (tantoPOR/1000)
  qbm <- PER2020$qxm * (tantoPOR/1000)
  lambdaf <- PER2020$lambdaf
  lambdam <- PER2020$lambdam
  qf <- qbf * exp(-lambdaf * (anyo + (0:120) - 2012)) 
  qm <- qbm * exp(-lambdam * (anyo + (0:120) - 2012)) 
  q <- prop*qf + (1-prop)*qm
  return(data.frame(qf, qm, q, row.names = PER2020$x))
}

PASEM2020.q <- function(tantoPOR = 1, prop = 0.5) {
  PASEM2020 <- tabla.PASEM2020()
  qf <- PASEM2020$qxf * (tantoPOR/1000)
  qm <- PASEM2020$qxm * (tantoPOR/1000)
  q <- prop*qf + (1-prop)*qm
  return(data.frame(qf, qm, q, row.names = PASEM2020$x))
}

# Funciones para calcular las cohortes
.cohorte <- function(tabla, lx0) {
  lxf <- lxm <- lx <- lx0
  qf <- tabla$qf
  qm <- tabla$qm
  q <- tabla$q
  omega <- nrow(tabla) - 1
  for (ii in 1:omega) {
    lxf[ii+1] <- lxf[ii] * (1-qf[ii])
    lxm[ii+1] <- lxm[ii] * (1-qm[ii])
    lx[ii+1] <- lx[ii] * (1-q[ii])
  }
  return(data.frame(lxf, lxm, lx, row.names = row.names(tabla)))
}

cohorte.PER2020 <- function(anyo = 2012, lx0 = 10^6) {
  tabla <- PER2020.q(anyo)
  return(.cohorte(tabla, lx0))
}

cohorte.PASEM2020 <- function(lx0 = 10^6) {
  tabla <- PASEM2020.q()
  return(.cohorte(tabla, lx0))
}

# Métodos de aproximaciones para edades no enteras
aprox.unif <- function(.lx, .lx.1, t) {
  aprox <- (1-t)*.lx + t*.lx.1
  return(aprox)
}

aprox.hazcons <- function(.lx, .lx.1, t) {
  aprox <- .lx^(1-t) * .lx.1^t
  return(aprox)
}

aprox.balducci <- function(.lx, .lx.1, t) {
  aprox <- ((1-t)/.lx + t/.lx.1)^(-1)
  return(aprox)
}

# Cohortes para edades no enteras
.cohorte.h <- function(tabla, h = 1, tipo = "unif") {
  aprox.fun <- aprox.unif
  if (tipo == "hr") {
    aprox.fun <- aprox.hazcons
  } else if (tipo == "bald") {
    aprox.fun <- aprox.balducci
  }
  omega <- nrow(tabla) - 1
  x <- seq(0, omega, 1/h) 
  x.h <- x * h
  resultados <- lapply(tabla, function(col) {
    sapply(x, function(x) {
      aprox.fun(col[floor(x) + 1], col[ceiling(x) + 1], x - floor(x))
    })
  })
  salida <- data.frame(resultados, row.names = x.h)
  names(salida) <- paste0(names(salida), ".h")
  
  return(salida)
}

cohorte.PER2020.h <- function(anyo = 2012, h = 1, tipo = "unif", lx0 = 10^6) {
  tabla <- cohorte.PER2020(anyo, lx0 = lx0)
  return(.cohorte.h(tabla, h, tipo))
}

cohorte.PASEM2020.h <- function(h = 1, tipo = "unif", lx0 = 10^6) {
  tabla <- cohorte.PASEM2020(lx0 = lx0)
  return(.cohorte.h(tabla, h, tipo))
}
