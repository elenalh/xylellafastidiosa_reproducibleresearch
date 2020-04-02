# Function to perform model selection based on WAIC criteria
# Function developed by Joaquín Martínez Minaya (jomartinez@bcamath.org) and Facundo Muñoz Viera (facundo.munoz@cirad.fr)

# Inputs:
# resp -> response variable
# covariates -> model covariates
# database -> database
# n -> number of models to show in the output
# family -> response variable family. Binomial is de default option


Bdiclcpomodel_stack<-function(resp, covariates, database, n, family="binomial",...)
{  
  
  sel.terms <- switch('terms',terms=covariates)
  
  comb.terms <- function(m, v=sel.terms) {
    if(m==0) return('resp ~ 1')
    else {
      combis <- apply(combn(v, m), 2, paste, collapse=' + ')
      return(paste('resp ~ 1', combis, sep=' + '))
    }
  }
  
  # List with all possible models
  f.list <- unlist(sapply(0:length(sel.terms), comb.terms))
  
  dic<-numeric()
  LCPO<-numeric()
  waic<-numeric()
  LCPOfailure<-numeric()
  for(i in 1:length(f.list)){
    res =inla(formula = eval(parse(text=f.list[i])), family=family, data=database, ...)
    dic[i] <- res$dic$dic
    LCPO[i] = -mean(log(res$cpo$cpo))
    waic[i]<-res$waic$waic
    LCPOfailure[i]<-sum(res$cpo$failure)
    print(c(i, dic[i], waic[i], LCPO[i]))
  }
  
  # models according to DIC
  models_dic<-data.frame(f.list[order(dic)[1:n]], dic[order(dic)[1:n]], waic[order(dic)[1:n]], LCPO[order(dic)[1:n]])
  colnames(models_dic)<-c("Models", "Dic", "Waic", "LCPO")
  
  # models according to WAIC
  models_waic<-data.frame(f.list[order(waic)[1:n]], dic[order(waic)[1:n]], waic[order(waic)[1:n]], LCPO[order(waic)[1:n]], LCPOfailure[order(waic)[1:n]])
  colnames(models_waic)<-c("Models", "Dic", "Waic", "LCPO", "LCPO-failure")
  
  
  # models according to CPO
  models_lcpo<-data.frame(f.list[order(LCPO)[1:n]], dic[order(LCPO)[1:n]], waic[order(LCPO)[1:n]], LCPO[order(LCPO)[1:n]], LCPOfailure[order(LCPO)[1:n]])
  colnames(models_lcpo)<-c("Models", "Dic", "Waic", "LCPO", "LCPO-failure")
  
  
  models<-list(models_dic, models_waic, models_lcpo)
  names(models)<-c("Models according to dic", "Models according to waic", "Models according to lcpo")
  models
  
}
