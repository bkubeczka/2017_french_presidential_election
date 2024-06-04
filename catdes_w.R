catdes.w<-function (donnee, num.var, proba = 0.05 , weight=NULL) 
{
  lab.sauv <- lab <- colnames(donnee)
  quali = NULL
  
  for (i in 1:length(lab))
  {
    lab[i] = gsub(" ", ".", lab[i])
    
    if (is.factor(donnee[, i]))
    {
      if (levels(donnee[, i])[1] == "")
      {	
        levels(donnee[, i])[1] = "NA"
      }
      
      if (i != num.var)
      {
        quali = c(quali, i)
      }
    }
  }
  
  quanti = (1:ncol(donnee))[-c(quali, num.var)]
  
  if (length(quanti) == 0)
  {
    quanti = NULL
  }
  
  colnames(donnee) = lab
  res = list()
  nb.modalite <- length(levels(donnee[, num.var]))
  nb.quali = length(quali)
  old.warn = options("warn")
  
  if(length(weight)==0)
  {
    weight<-rep(1,dim(donnee)[1])
  }
  
  if (nb.quali > 0) 
  {
    options(warn = -1)
    
    marge.li = xtabs(weight~donnee[, num.var])
    nom = tri = structure(vector(mode = "list", length = nb.modalite), names = levels(donnee[, num.var]))
    
    for (i in 1:nb.quali) 
    {
      Table <- xtabs(weight~donnee[, num.var] + donnee[, quali[i]])
      marge.col = xtabs(weight~donnee[, quali[i]])
      
      ML<-rowSums(Table)
      
      for (j in 1:nlevels(donnee[, num.var]))
      {
        for (k in 1:nlevels(donnee[, quali[i]]))
        {
          aux2 = Table[j, k]/ML[j]
          
          if(ML[j]==0)
          {
            aux2 = 0
          }
          
          aux3 = marge.col[k]/sum(marge.col)
          
          if (aux2 > aux3) 
          {
            aux4 = phyper(Table[j, k] - 1, ML[j], sum(ML) - ML[j], marge.col[k], lower.tail = FALSE) * 2
          }
          else 
          {
            aux4 = phyper(Table[j, k], ML[j], sum(ML) - ML[j], marge.col[k]) * 2
          }
          
          if (aux4 < proba)
          {
            aux5 = (1 - 2 * as.integer(aux2 > aux3)) * qnorm(aux4/2)
            aux1 = Table[j, k]/marge.col[k]
            
            tri[[j]] = rbind(tri[[j]], c(aux1 * 100, aux2 * 100, aux3 * 100, aux4, aux5))
            nom[[j]] = rbind(nom[[j]], c(levels(donnee[,quali[i]])[k], colnames(donnee)[quali[i]]))
          }
        }
      }
    }
    
    for (j in 1:nb.modalite)
    {
      if (!is.null(tri[[j]]))
      {
        oo = rev(order(tri[[j]][, 5]))
        tri[[j]] = tri[[j]][oo, ]
        nom[[j]] = nom[[j]][oo, ]
        
        if (nrow(matrix(tri[[j]], ncol = 5)) > 1)
        {
          rownames(tri[[j]]) = paste(nom[[j]][, 2], nom[[j]][,1], sep = "=")
        }
        else
        {
          tri[[j]] = matrix(tri[[j]], ncol = 5)
          rownames(tri[[j]]) = paste(nom[[j]][2], nom[[j]][1], sep = "=")
        }
        
        colnames(tri[[j]]) = c("Cla/Mod", "Mod/Cla", "Global", "p.value", "v.test")
      }
      
    }
    
    res$category = tri
  }
  
  if (!is.null(quanti))
  {
    nom = result = structure(vector(mode = "list", length = nb.modalite), names = levels(donnee[, num.var]))
    
    for (i in 1:length(quanti))
    {
      moy.mod = by(donnee[, quanti[i]]*weight, donnee[, num.var], mean, na.rm = TRUE)
      n.mod = summary(donnee[, num.var])
      
      sd.mod = by(donnee[, quanti[i]]*weight, donnee[, num.var], sd, na.rm = TRUE)
      sd.mod = sd.mod * sqrt((n.mod - rep(1, nb.modalite))/n.mod)
      
      moy = mean(donnee[, quanti[i]]*weight, na.rm = TRUE)
      et = sd(donnee[, quanti[i]]*weight, na.rm = TRUE) * sqrt(1 - 1/sum(n.mod))
      
      for (j in 1:nb.modalite)
      {
        v.test = (moy.mod[j] - moy)/et * sqrt(n.mod[j])/sqrt((sum(n.mod) - n.mod[j])/(sum(n.mod) - 1))
        p.value = pnorm(abs(v.test), lower.tail = FALSE) * 2
        
        if (!is.na(v.test))
        {
          if (abs(v.test) > -qnorm(proba/2))
          {
            result[[j]] = rbind(result[[j]], c(v.test, moy.mod[j], moy, sd.mod[j], et, p.value))
            nom[[j]] = c(nom[[j]], colnames(donnee)[quanti[i]])
          }
        }
      }
    }
    
    for (j in 1:nb.modalite)
    {
      if (!is.null(result[[j]]))
      {
        oo = rev(order(result[[j]][, 1]))
        result[[j]] = result[[j]][oo, ]
        nom[[j]] = nom[[j]][oo]
        
        if (nrow(matrix(result[[j]], ncol = 6)) > 1)
        {
          rownames(result[[j]]) = nom[[j]]
          colnames(result[[j]]) = c("v.test", "Mean in category", "Overall mean", "sd in category", "Overall sd", "p.value")
        }
        else
        {
          result[[j]] = matrix(result[[j]], ncol = 6)
          rownames(result[[j]]) = nom[[j]]
          colnames(result[[j]]) = c("v.test", "Mean in category", "Overall mean", "sd in category", "Overall sd", "p.value")
        }
      }
    }
    
    res$quanti = result
  }
  
  options(old.warn)
  class(res) <- c("catdes", "list ")
  return(res)
}