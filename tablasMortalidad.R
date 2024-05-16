#rm(list=ls())


### Tablas de Mortalidad

print("Funciones basicas:")
print("cohorte.PER2020.h(anyo = 2012, h = 1)")
print("cohorte.PASEM2020.h(h = 1)")

.tabla.PER2020 <- read.csv("https://raw.githubusercontent.com/aozoro/misNotasVida/main/tablas/PrimerPER2020.csv")
.tabla.PASEM2020 <- read.csv("https://raw.githubusercontent.com/aozoro/misNotasVida/main/tablas/PrimerPASEM2020.csv")


# Definición de las funciones para leer las tablas de mortalidad
tabla.PER2020 <- function() {
  .tabla.PER2020 
}

tabla.PASEM2020 <- function() {
  .tabla.PASEM2020
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
  #salida <- data.frame(resultados, row.names = x.h)
  #names(salida) <- paste0(names(salida), ".h")
  salida <- resultados$lx
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

.u <- function(u0, t, m, h,tipo=0, h.prim=h, u1=0){
  #tipo 0: constante, 1: aritmetica , 2: geometrica
  k = h/h.prim
  if (tipo==0) {
    u <- u0 + 0*(t-m*h)
  }  
  
  if (tipo==1) {
    u <- u0 + u1*floor((t-m*h)/k)
  }
  
  if (tipo==2) {
    u <- u0 * u1^floor((t-m*h)/k)
  }
  return(u)
}

.renta <- function(m,n,x,any,h,u,venc,I1){
  lx.h <- cohorte.PER2020.h(anyo = any,h = h)
  t <- (m*h):((m+n)*h-1)
  v <- (1+I1)^(-(t+venc)/h)
  p <- lx.h[x*h+venc+t+1]/lx.h[x*h+1]
  
  return(sum(u*v*p))
}

renta.constante<- function(m,n,x,any.contrato,h,u0,venc,I1,h.prim=h, vitalicio=FALSE){
  omega <- 120
  if (vitalicio) {n <- omega -  x - m}
  
  u <- .u(u0=u0, t=(m*h):((m+n)*h-1), m=m, h=h,tipo=0, h.prim=h.prim)
  
  .renta(m=m,n=n,x=x,any=any.contrato-x,h=h,u=u,venc=venc,I1=I1)
}

renta.aritmetica <- function(m,n,x,any.contrato,h,u0,a,venc,I1,h.prim=h, vitalicio=FALSE){
  omega <- 120
  if (vitalicio) {n <- omega -  x - m}
  
  u <- .u(u0=u0, t=(m*h):((m+n)*h-1), m=m, h=h,tipo=1, h.prim=h.prim,u1 = a)
  
  .renta(m=m,n=n,x=x,any=any.contrato-x,h=h,u=u,venc=venc,I1=I1)
}

renta.geometrica <- function(m,n,x,any.contrato,h,u0,g,venc,I1,h.prim=h, vitalicio=FALSE){
  omega <- 120
  if (vitalicio) {n <- omega -  x - m}
  
  u <- .u(u0=u0, t=(m*h):((m+n)*h-1), m=m, h=h,tipo=2, h.prim=h.prim,u1 = g)
  
  .renta(m=m,n=n,x=x,any=any.contrato-x,h=h,u=u,venc=venc,I1=I1)
}


.seguro <- function(m,n,x,h,u,venc,I1){
  lx.h <- cohorte.PASEM2020.h(h = h)
  t <- (m*h):((m+n)*h-1)
  v <- (1+I1)^(-(t+venc)/h)
  q <- (lx.h[x*h+t+1]-lx.h[x*h+t+2])/ lx.h[x*h+1]
  return(sum(u*v*q))
}


seguro.constante<- function(m,n,x,h,u0,I1,h.prim=h, vitalicio=FALSE){
  venc <- 1 
  omega <- 109
  if (vitalicio) {n <- omega -  x - m}
  
  u <- .u(u0=u0, t=(m*h):((m+n)*h-1), m=m, h=h,tipo=0, h.prim=h.prim)
  
  .seguro(m=m,n=n,x=x,h=h,u=u,venc=venc,I1=I1)
}

seguro.aritmetica <- function(m,n,x,h,u0,a,I1,h.prim=h, vitalicio=FALSE){
  venc <- 1
  omega <- 109
  if (vitalicio) {n <- omega -  x - m}
  
  u <- .u(u0=u0, t=(m*h):((m+n)*h-1), m=m, h=h,tipo=1, h.prim=h.prim,u1 = a)
  
  .seguro(m=m,n=n,x=x,h=h,u=u,venc=venc,I1=I1)
}

seguro.geometrica <- function(m,n,x,h,u0,g,I1,h.prim=h, vitalicio=FALSE){
  venc <- 1
  omega <- 120
  if (vitalicio) {n <- omega -  x - m}
  
  u <- .u(u0=u0, t=(m*h):((m+n)*h-1), m=m, h=h,tipo=2, h.prim=h.prim,u1 = g)
  
  .seguro(m=m,n=n,x=x,h=h,u=u,venc=venc,I1=I1)
}


seguro.constante.cont <- function(m,n,x,h,u0,I1,h.prim=h, vitalicio=FALSE){
  if (vitalicio) {n <- omega -  x - m}
  ((1+I1)^(1/h)-1)*h/log(1+I1)*seguro.constante(m,n,x,h,u0,I1,h.prim=h, vitalicio=FALSE)
}
seguro.aritmetica.cont <- function(m,n,x,h,u0,a,I1,h.prim=h, vitalicio=FALSE){
  if (vitalicio) {n <- omega -  x - m}
  ((1+I1)^(1/h)-1)*h/log(1+I1)*seguro.aritmetica(m,n,x,h,u0,a,I1,h.prim=h, vitalicio=FALSE)
}
seguro.geometrica.cont <- function(m,n,x,h,u0,g,I1,h.prim=h, vitalicio=FALSE){
  if (vitalicio) {n <- omega -  x - m}
  ((1+I1)^(1/h)-1)*h/log(1+I1)*seguro.geometrica(m,n,x,h,u0,g,I1,h.prim=h, vitalicio=FALSE)
}



.seguro.med <- function(m,n,x,h,u,venc,I1){
  lx.h <- cohorte.PASEM2020.h(h = h)
  t <- (m*h):((m+n)*h-1)
  v <- (1+I1)^(-(t+0.5)/h)
  q <- (lx.h[x*h+t+1]-lx.h[x*h+t+2])/ lx.h[x*h+1]
  return(sum(u*v*q))
}


seguro.constante.med<- function(m,n,x,h,u0,I1,h.prim=h, vitalicio=FALSE){
  venc <- 1 
  omega <- 109
  if (vitalicio) {n <- omega -  x - m}
  
  u <- .u(u0=u0, t=(m*h):((m+n)*h-1), m=m, h=h,tipo=0, h.prim=h.prim)
  
  .seguro.med(m=m,n=n,x=x,h=h,u=u,venc=venc,I1=I1)
}

seguro.aritmetica.med <- function(m,n,x,h,u0,a,I1,h.prim=h, vitalicio=FALSE){
  venc <- 1
  omega <- 109
  if (vitalicio) {n <- omega -  x - m}
  
  u <- .u(u0=u0, t=(m*h):((m+n)*h-1), m=m, h=h,tipo=1, h.prim=h.prim,u1 = a)
  
  .seguro.med(m=m,n=n,x=x,h=h,u=u,venc=venc,I1=I1)
}

seguro.geometrica.med <- function(m,n,x,h,u0,g,I1,h.prim=h, vitalicio=FALSE){
  venc <- 1
  omega <- 120
  if (vitalicio) {n <- omega -  x - m}
  
  u <- .u(u0=u0, t=(m*h):((m+n)*h-1), m=m, h=h,tipo=2, h.prim=h.prim,u1 = g)
  
  .seguro.med(m=m,n=n,x=x,h=h,u=u,venc=venc,I1=I1)
}

.renta.varianza <- function(m,n,x,any,h,u,venc,I1){
  p <- cohorte.PER2020.h(anyo = any,h = h)
  t <- (m*h):((m+n)*h-1)
  v <- (1+I1)^(-(t+venc)/h)
  probabilidades <- (p[x*h+t+venc+1]/p[x*h+1])*(1-p[x*h+t+venc+2]/p[x*h+t+venc+1])
  #probabilidades[is.na(probabilidades)] <- 0
  montos <- c(u*v)
  es <- sum(probabilidades*cumsum(montos))
  var <- sum(probabilidades*cumsum(montos)^2)-es^2
  return(list("media"=es,"varianza"=var))
}

renta.constante.varianza<- function(m,n,x,any.contrato,h,u0,venc,I1,h.prim=h, vitalicio=FALSE){
  omega <- 119
  if (vitalicio) {n <- omega -  x - m}
  
  u <- .u(u0=u0, t=(m*h):((m+n)*h-1), m=m, h=h,tipo=0, h.prim=h.prim)
  
  .renta.varianza(m=m,n=n,x=x,any=any.contrato-x,h=h,u=u,venc=venc,I1=I1)
}

renta.aritmetica.varianza <- function(m,n,x,any.contrato,h,u0,a,venc,I1,h.prim=h, vitalicio=FALSE){
  omega <- 119
  if (vitalicio) {n <- omega -  x - m}
  
  u <- .u(u0=u0, t=(m*h):((m+n)*h-1), m=m, h=h,tipo=1, h.prim=h.prim,u1 = a)
  
  .renta.varianza(m=m,n=n,x=x,any=any.contrato-x,h=h,u=u,venc=venc,I1=I1)
}

renta.geometrica.varianza <- function(m,n,x,any.contrato,h,u0,g,venc,I1,h.prim=h, vitalicio=FALSE){
  omega <- 119
  if (vitalicio) {n <- omega -  x - m}
  
  u <- .u(u0=u0, t=(m*h):((m+n)*h-1), m=m, h=h,tipo=2, h.prim=h.prim,u1 = g)
  
  return(.renta.varianza(m=m,n=n,x=x,any=any.contrato-x,h=h,u=u,venc=venc,I1=I1))
}

h.texto <- function(h){
  H <- c(1,2,3,4,6,12)
  per <- c("años","semestres","cuatrimestres",
           "trimestres","semestres","meses")
  return(per[H==h])
}

# Definir la función con los argumentos especificados
esquemaTemporal <- function(x, n, venc, u0, m, h, tipo, h.prim, u1) {
  # Calcular el tiempo y los índices
  t <- (m * h):((m + n) * h - 1)
  
  if (m == 0) {
    cero <- NULL
  } else {
    cero <- 0
  }
  
  df <- data.frame(
    x = c(cero, 7.5, 12.5, 17.5, 27.5, 32.5),
    tiempo = c(cero, m * h, m * h + 1, m * h + 2, (m + n) * h - 1, (m + n) * h)
  )
  
  df$edad <- df$tiempo + x * h
  
  u.idx <- df$tiempo - m * h - venc
  u.idx[u.idx < 0 | u.idx == n * h] <- NA
  
  u <- .u(u0 = u0, t = t, m = m, h = h, tipo = tipo, h.prim = h.prim, u1 = u1)
  u <- round(u[u.idx + 1], 0)
  u[!is.na(u)] <- paste0("u(", u.idx[!is.na(u)], ") = ", u[!is.na(u)])
  
  ggplot(df, aes(x = x)) + 
    annotate("segment", x = min(df$x), xend = max(df$x), y = 0, yend = 0, color = "black") +
    annotate("text", x = df$x, y = 0, label = df$tiempo, vjust = 1.8) +
    annotate("text", x = 35.6, y = 0, label = h.texto(h), vjust = 2.3, size = 2.8) +
    annotate("text", x = df$x, y = 0, label = df$edad, vjust = 3.8) +
    annotate("text", x = 35.6, y = 0, label = "edad en", vjust = 5.8, size = 2) +
    annotate("text", x = 35.6, y = 0, label = h.texto(h), vjust = 7, size = 2) +
    annotate("text", x = df$x[!is.na(u)], y = 0, label = u[!is.na(u)], hjust = -0.2, angle = 90) +
    annotate("text", x = 35.6, y = 0, label = "euros", vjust = -1.8, size = 2.8) +
    annotate("text", x = 22.5, y = 0, label = ". . .", vjust = 1.8) +
    annotate("text", x = 22.5, y = 0, label = ". . .", vjust = 3.8) +
    annotate("text", x = 22.5, y = 0, label = ". . .", vjust = -2.5) +
    geom_point(aes(y = 0), color = "black", size = 2) +
    scale_y_continuous(limits = c(-2, 4)) +
    scale_x_continuous(limits = c(min(df$x), 36.5)) +
    theme_minimal() +
    labs(title = NULL, x = NULL, y = NULL) +
    theme(
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5)
    )
}