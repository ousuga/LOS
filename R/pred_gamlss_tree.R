#' Prediction gamlss tree
#'
#' Prediction gamlss tree.
#'
#' @param newdata	a data.frame.
#' @param objeto gamlss_tree list.
#'
#' @return
#' \code{pred_gamlss_tree} gives a data.frame with data, pred_arbol and pred_gam.
#'
#' @examples
#' require(COUNT)
#' data(azpro)
#' model_amg<-gamlss_tree(los~.,datos=azpro)
#' pred<-pred_gamlss_tree(objeto=model_gam,newdata=azpro)
#' 
#' @export
pred_gamlss_tree<-function(newdata=newdata,objeto)
{
  fo<<-objeto$form
  cyc<<-50
  arbol<<-objeto$arboles$commit1
  new_data_variables<<-data.frame(model.frame(fo,newdata))
  dis_sel<-objeto$dis_sel
  nodos_train<<-objeto$nodos_train
  lista_nodos=NULL
  fo_gamlss=NULL
  datos_nodos_gamlss=NULL
  new_data_pred=NULL
  newdata$pred_arbol<-predict(arbol, new_data_variables)
  nodos<-length(unique(predict(arbol)))
  for(i in 1:nodos )
  {
    assign("n",i,.GlobalEnv)
    lista_nodos[[i]]<-subset(newdata,pred_arbol==objeto$valor_nodo[[1]][i])
    names(lista_nodos)[i]<-paste("nodo",i,sep="")
    fo_gamlss[[i]]<-objeto$modelosxnodo$commit1[[1]]$mu.formula
    datos_nodos_gamlss[[i]]<- data.frame(model.frame(fo_gamlss[[i]],lista_nodos[[i]]))
    model<<-objeto$modelosxnodo$commit1[[i]]
    da<<- datos_nodos_gamlss[[i]]
    if(is.list(objeto$modelosxnodo$commit1[[i]])==T)
    {
      lista_nodos[[i]]$pred_gam<-predict(model, newdata= da)
    }
    else
    {
      lista_nodos[[i]]$pred_gam<-lista_nodos[[i]]$pred_arbol
    }
    new_data_pred<-rbind(new_data_pred, lista_nodos[[i]])
  }
  return(new_data_pred)
}