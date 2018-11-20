#' @name gamlss_tree
#'
#' @title
#' gamlss_tree
#'
#' @description
#' gamlss_tree.
#'
#' @param form percentiles.
#' @param perc percentiles.
#' @param n_dist_mod percentiles.
#' @param var_sel percentiles.
#' @param steps percentiles.
#' @param porc_entre percentiles.
#' @param committess percentiles.
#' @param nom_dist percentiles.
#' @param cyc percentiles.
#' @param prueba_hip percentiles.
#' @param acepta_h percentiles.
#' @param type percentiles.
#' @param arbol_activo percentiles.
#' 
#'
#' @details
#' The gamlss_tree generates two dataset.
#'
#' @return
#' \code{gamlss_tree} gives two columns.

#' @examples
#' ## The gamlss_tree function
#' @export
gamlss_tree<-function(form, datos, n_dist_mod=4,var_sel="aicmodelo",steps=2,
                      porc_entre=0.8,committess=1,
                      nom_dist=c( "exGAUS","GIG","GG","BCCGo","BCPEo","GA", "GB2",
                                  "BCTo","WEI3", "LOGNO", "EXP", "PARETO2","IG",
                                  "IGAMMA","NO"),
                      cyc=50,prueba_hip=TRUE, acepta_h=FALSE,type="counts",arbol_activo=TRUE)
  
{  
 arbol_cart_comm=NULL
  nodos_comm=NULL
  modelonodos_comm=NULL
  variables=NULL
  list_cov=NULL
  nodos_train=NULL
  nodos_test=NULL
  fo=form
  datos_model<-model.frame(fo,datos)
  list_cov<-names(datos_model)
  for (m in 2:length(list_cov))
  {
    if(is.null(variables)==T)
    {
      variables<-list_cov[m]
    }
    else
    {
      variables<-paste(variables,"+",list_cov[m],sep="")
    }  
  }
  datos_sep=split_muestra(datos=datos_model,perc=porc_entre)
  datos_train=datos_sep$train
  datos_test=datos_sep$test
  for(p in 1:committess)
  {      
    if(arbol_activo==T)
    {
      arbol_cart<-rpart(fo,data=datos_train)
      arbol_cart_comm[[p]]<-arbol_cart
      names( arbol_cart_comm)[p]=paste("commit",p,sep="")
      variables_arb<-names(arbol_cart$variable.importance)
      arbol_cart_predict_e<-predict(arbol_cart)
      arbol_cart_predict_t<-predict(arbol_cart,datos_test)
       datos_pred_e<-data.frame(datos_train,arbol_cart_predict_e)
      datos_pred_t<-data.frame(datos_test,arbol_cart_predict_t)
      nodos<-unique(arbol_cart_predict_e)
      nodos_comm[[p]]<-nodos
      n_nodo<-length(nodos)
      for (i in 1:n_nodo)
      {
        nodos_train[[i]]<-subset(datos_pred_e,datos_pred_e$arbol_cart_predict_e==nodos[i])
        names(nodos_train)[i]<-paste("nodo",i,sep="")
        nodos_test[[i]]<-subset(datos_pred_t,datos_pred_t$arbol_cart_predict_t==nodos[i])
        names(nodos_test)[i]<-paste("nodo",i,sep="")
      }
      
    }
    else
    {
      
      n_nodo=1
      nodos_train[[1]]<-datos_train
      names(nodos_train)[1]<-"nodo1"
      nodos_test[[1]]<-datos_test
      names(nodos_test)[1]<-"nodo1"
    }
    n_variables<-length(names(datos_train))-1
    variables_arb<-names(datos_train)[-1]
    nodo_prueba=NULL
    dist_acept=NULL
    n_dist_nodo=NULL
    for (i in 1:n_nodo)
    {
      if(prueba_hip==T)
      {
        nodo_prueba[[i]]<-test_dist(eval(parse(text=paste("nodos_train$nodo",i,sep="")))[,1],type=type,extra=nom_dist,signif=0.05)
        names(nodo_prueba)[i]=paste("nodo",i,sep="")
        nodo_actual<-data.frame(nodo_prueba[[i]])
      }
      else
      {
        nodo_prueba[[i]]<-fitDist(eval(parse(text=paste("nodos_train$nodo",i,sep="")))[,1],type=type,extra=nom_dist)
        nodo_prueba[[i]]<- data.frame("Dist"=names(nodo_prueba[[i]]$fits),"AIC"=nodo_prueba[[i]]$fits)
        row.names(nodo_prueba[[i]]) =NULL  
        names(nodo_prueba)[i]=paste("nodo",i,sep="")
        nodo_actual<-nodo_prueba[[i]]
        } 
      
      if(acepta_h==T)
      {
        dist_acept[[i]]<-nodo_actual%>%
          filter(res_prueba_ks=="Acepta" | res_prueba_ad=="Acepta")%>%
          select(Dist)%>%head(n_dist_mod)
      }
      else
      {
        dist_acept[[i]]<-nodo_actual%>%
          select(Dist)%>%head(n_dist_mod)  
      }  
      n_dist_nodo[[i]]<-length (dist_acept[[i]][,1])
      
    }
    {
       {
        
        nodos_ajuste_cov=NULL
       }
      s=NULL
      for (j in 1:n_nodo)
      {
        formula=""
        for (i in 1:n_variables)
        {
          if (length(unique(eval(parse(text=paste("nodos_train[[",j,"]]$",variables_arb[i],sep="")))))>1)
          {
             if(formula=="") 
            {
              formula<-paste(formula,variables_arb[i],sep="")
            }
            else
            {
              formula<-paste(formula,"+",variables_arb[i],sep="")
            }
          }
        }
        aicmodelo=NULL
        s[j]=1
        mae_gamlss=NULL
        mse_gamlss=NULL
        rmse_gamlss=NULL
        nodo_test_gamlss<<-data.frame(model.frame(paste(form[[2]],"~",formula,sep=""),nodos_test[[j]]))
        rm(da)
        da<<-nodos_train[[j]]
        if (n_dist_nodo[j]>0)
        {
          for (k in 1:n_dist_nodo[j])  
          {
            dist<-as.character(dist_acept[[j]][k,1])
            modelocov<<-
              tryCatch(gamlss(eval(parse(text=paste(fo[[2]],"~",formula,sep=""))),
                              data=da,i.control=glim.control(cyc=cyc), family=dist),error=function(e)-1)
            #modelo<-eval(parse(text=paste("modelonodo",k,"_",j,sep="")))
            aicmodelo=rbind(aicmodelo,tryCatch(modelocov$aic,error=function(e)NA))
            pred<-tryCatch(predict(modelocov,newdata = nodo_test_gamlss),error=function(e)NA)
            mae_gamlss<-rbind(mae_gamlss,mae(nodo_test_gamlss[,1],pred))
            mse_gamlss<-rbind(mse_gamlss,mse(nodo_test_gamlss[,1],pred))
            rmse_gamlss<-rbind(rmse_gamlss,sqrt(mse_gamlss[k]))
          }# cierre for de distribucion
        }  
        {
          mae_gamlss<-c(mae_gamlss)
          mse_gamlss<-c(mse_gamlss)
          rmse_gamlss<-c(rmse_gamlss)
          aicmodelo<-c(aicmodelo)
          nodos_ajuste_cov[[j]]  =  data.frame(dist_acept[[j]],aicmodelo,mae_gamlss,mse_gamlss,rmse_gamlss)
          names(nodos_ajuste_cov)[j]=paste("nodo",j,sep="")
        }
      } ##Cierre for de nodo
    }
    dis_sel=NULL
    #n=1
    for (n in 1:n_nodo) 
    {
      ordenado<-arrange(eval(parse(text=paste("nodos_ajuste_cov$nodo",n,sep=""))),eval(parse(text=var_sel)))
      dis_sel[n]<-as.character(ordenado$Dist[1])
    }
    modelonodos=NULL
    for (n in 1:n_nodo)
    {
      modelonodos[[n]]<- tryCatch(gamlss(eval(parse(text=paste(fo[[2]],"~",formula,sep=""))),
                                         data=nodos_train[[n]],i.control=glim.control(cyc=cyc), family=dis_sel[n]),error=function(e)NA)
      modelonodos[[n]]<-tryCatch(stepGAIC(modelonodos[[n]],steps=steps),error=function(e)modelonodos[[n]])
      
      names(modelonodos)[n]=paste("nodo",n,sep="")
    }
    modelonodos_comm[[p]]<-modelonodos
    names(modelonodos_comm)[[p]]=paste("commit",p,sep="")
    datos_train=NULL
    for (i in 1:n_nodo)
    {  
      if(is.list(modelonodos[[i]])==T)  
      {  
        nodos_train[[i]]$y = nodos_train[[i]][,1] -(predict(modelonodos[[i]])-nodos_train[[i]][,1])
      }
      else
      {
        nodos_train[[i]]$y= nodos_train[[i]][,1] -(mean(nodos_train[[i]][,1])-nodos_train[[i]][,1])
      }
      datos_train<-rbind(datos_train,data.frame(model.frame(eval(parse(text=paste("y~",variables,sep=""))),
                                                            data=nodos_train[[i]])))
    }
    datos_train$y<-case_when(
    datos_train$y<=0~0.1,
    datos_train$y>0~datos_train$y)
    fo<-as.formula(paste('y~',variables,sep=""))
     
  } 
  return(list(arboles=arbol_cart_comm,
              modelosxnodo=modelonodos_comm,
              valor_nodo=nodos_comm,
              dis_sel=dis_sel,nodos_train=
                nodos_train,form=form,datos_test=datos_test,
              nodo_prueba=nodo_prueba,nodos_ajuste_cov=nodos_ajuste_cov))
} 
#' 
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
    datos_nodos_gamlss[[i]]<-
    data.frame(model.frame(fo_gamlss[[i]],lista_nodos[[i]]))
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
#' 
split_sample<-function(datos,perc)
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
#' 
test_dist<-function(y, type="realplus", extra="NO", signif=0.05)
{
  order_fit<-fitDist(y,type=type,extra=extra)
  nom_dist<-names(order_fit$fits)
  n_dist<- length(nom_dist)
  resultado=NULL
  acepta=0
  for (i in 1:n_dist)
  { 
    rm(ks_fit)
    rm(ad_fit)
    fit<-gamlssML(y,family=nom_dist[i])
    n_param=length(fit$parameters)
    pt_dis<-eval(parse(text=paste("p",nom_dist[i],sep="")))
    if (n_param==1)
    {
      ks_fit<- tryCatch(ks.test(x=y, y=pt_dis,
                                mu=fit$mu),error=function(e) -1)
      
      ad_fit<- tryCatch(ad.test(y,null=pt_dis,
                                mu=fit$mu),error=function(e) -1)
    } 
    if (n_param==2)
    {
     ks_fit<- tryCatch(ks.test(x=y, y=pt_dis,
                                mu=fit$mu,
                                sigma=fit$sigma),error=function(e) -1)
      
      ad_fit<- tryCatch(ad.test(y,null=pt_dis,
                                mu=fit$mu,sigma=fit$sigma),error=function(e) -1)
    }
    
    if (n_param==3)
    {
     ks_fit<- tryCatch(ks.test(x=y, y=pt_dis,
                                mu=fit$mu,
                                sigma=fit$sigma,
                                nu=fit$nu),error=function(e) -1)
      
      ad_fit<- tryCatch(ad.test(y,null=pt_dis,
                                mu=fit$mu,
                                sigma=fit$sigma,
                                nu=fit$nu),error=function(e) -1)
    }
    if (n_param==4)
    {
      ks_fit<- tryCatch(ks.test(x=y, y=pt_dis,
                                mu=fit$mu,
                                sigma=fit$sigma,
                                nu=fit$nu,
                                tau=fit$tau),error=function(e) -1)
      
      ad_fit<- tryCatch(ad.test(y,null=pt_dis,
                                mu=fit$mu,
                                sigma=fit$sigma,
                                nu=fit$nu,
                                tau=fit$tau),error=function(e) -1)
    }
    if (n_param>4)
    {
     Print("La distribuci?n ajustada tiene m?s de 4 par?metros")
    }
    if (typeof(ks_fit)!="list"){ks_fit$p.value=-1}
    if (typeof(ad_fit)!="list"){ad_fit$p.value=-1}
    resultado_add<-data.frame(nom_dist[i],ks_fit$p.value,ad_fit$p.value,fit$sbc,fit$aic)
    resultado<-rbind(resultado,resultado_add)
  }
  nom_prueba<-c("Dist","KS", "AD","SBC","AIC")
  names(resultado)<-nom_prueba
  resultado$res_prueba_ks<-case_when(
    resultado$KS==-1 ~"No Aplica",
    resultado$KS<=signif ~"Rechaza",
    resultado$KS>signif~"Acepta"
  )
 resultado$res_prueba_ad<-case_when(
    resultado$AD==-1~"No Aplica",
    resultado$AD<=signif ~"Rechaza",
    resultado$AD>signif~"Acepta")
  return(resultado=resultado)
}