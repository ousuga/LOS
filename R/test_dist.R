#' Test distributions
#'
#' @description
#' Test distributions.
#'
#' @param y	a formula object, with the response on the left of an ~ operator, 
#' and the terms, separated by + operators, on the right. 
#' @param type type of distribution.
#' @param extra additional distributions.
#' @param signif significance level.
#'
#' @return
#' \code{test_dist} gives a data.frame with p-value, AIC, SBC and result of 
#' goodness of fitness test.
#' 
#' 
#' @examples
#' require(COUNT)
#' data(azpro)
#' test<-test_dist(y=azpro$los)
#' test
#' 
#' @export
test_dist<-function(y, type="realplus", extra="NO", signif=0.05)
{
  options(warn=-1)
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