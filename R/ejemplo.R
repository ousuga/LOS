#' @name sm
#'
#' @title
#' Separar muestra
#'
#' @description
#' El objetivo es separa las muestra.
#'
#' @param datos conjunto de datos.
#' @param perc percentil.
#'
#' @details
#' La función divide la muestra en dos partes.
#'
#' @return
#' retorna la muestra divididaa
#'
#' @export
#' @examples
#' ## TDivisión de la función uno
#'
sm<-function(datos,perc)
{
  id_muestra<-sample(2,length(datos[,1]),replace = T,prob=c(perc,1-perc))
  fila<-seq(1:length(datos[,1]))
  muestra<-data.frame(id_muestra,fila)

  sub<-subset(muestra[,2],muestra$id_muestra==1)
  test<-subset(muestra$fila,muestra$id_muestra==2)

  datos_train<-datos[sub,]
  datos_test<-datos[test,]

  return(list(train=datos_train,test=datos_test))
}


