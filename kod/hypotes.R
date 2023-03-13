hypotes2 = function(sigma,matningar,nollH,alfa,riktning,medel){
#  windows()
  layout(matrix(c(1,2),2,1))
  x = seq(nollH - 3.5*sigma/sqrt(matningar),nollH + 3.5*sigma/sqrt(matningar),length.out = 1000)
  y = dnorm(x,nollH,sigma/sqrt(matningar))
  plot(x,y,type = 'l',col = 'red',ylab = '',main = 'Test med fixt alfa-varde', xlab = get.x.lab(sigma,matningar,nollH,alfa,riktning))
  points(c(nollH,nollH),c(0,dnorm(nollH,nollH,sigma/sqrt(matningar))),type = 'l',lty = 2, col = 'red')
  
  if(riktning == '<'){
    k = nollH - qnorm(1 - alfa,0,1)*sigma/sqrt(matningar)
    plot.area(nollH - 3.5*sigma/sqrt(matningar),k,nollH,sigma,matningar,'red')
  }else if(riktning == '>'){
    k = nollH + qnorm(1 - alfa,0,1)*sigma/sqrt(matningar)
    plot.area(k,nollH + 3.5*sigma/sqrt(matningar),nollH,sigma,matningar,'red')
  }else if(riktning == '!='){
    kp = nollH + qnorm(1 - alfa/2,0,1)*sigma/sqrt(matningar)
    km = nollH - qnorm(1 - alfa/2,0,1)*sigma/sqrt(matningar)
    plot.area(nollH - 3.5*sigma/sqrt(matningar), km ,nollH, sigma , matningar ,'red')
    plot.area(kp,nollH + 3.5*sigma/sqrt(matningar),nollH ,sigma,matningar,'red')
  }
  points(c(medel,medel),c(0,dnorm(medel,nollH,sigma/sqrt(matningar))),type = 'l')
  text.fig1.2(sigma,matningar,nollH,riktning,alfa,medel)
  
  plot(x,y,type = 'l',col = 'red',ylab = '',main = 'Test med direktmetoden',xlab = '')
  points(c(nollH,nollH),c(0,dnorm(nollH,nollH,sigma/sqrt(matningar))),type = 'l',lty = 2, col = 'red')
  if(riktning == '!='){
    if(medel < nollH){
      plot.area(nollH - 3.5*sigma/sqrt(matningar),medel,nollH,sigma,matningar,'purple')
    }else{
      plot.area(medel,nollH + 3.5*sigma/sqrt(matningar),nollH,sigma,matningar,'purple')
    }
  }else{
    if(riktning == '<'){
      plot.area(nollH - 3.5*sigma/sqrt(matningar),medel,nollH,sigma,matningar,'purple')
    }else{
      plot.area(medel,nollH + 3.5*sigma/sqrt(matningar),nollH,sigma,matningar,'purple')
    }
  }
  text.fig2.2(sigma,matningar,nollH,riktning,alfa,medel)
  
}


# hypotes1

hypotes1 = function(sigma,matningar,nollH,alfa,riktning,sannaVardet){
#  windows()
  layout(matrix(c(1,2),2,1))
  x = seq(nollH - 3.5*sigma/sqrt(matningar),nollH + 3.5*sigma/sqrt(matningar),length.out = 1000)
  y = dnorm(x,nollH,sigma/sqrt(matningar))
  plot(x,y,type = 'l',col = 'red',ylab = '',main = sprintf('Kritiskt omrade, H0: mu = %4.1f, H1: mu %s %3.1f',nollH,riktning,nollH),xlab = get.x.lab(sigma,matningar,nollH,alfa,riktning))
  points(c(nollH,nollH),c(0,dnorm(nollH,nollH,sigma/sqrt(matningar))),type = 'l',lty = 2, col = 'red')
  
  if(riktning == '<'){
    riktning.mindre(sigma,matningar,nollH,alfa,sannaVardet,x,y)  
  }else if(riktning == '>'){
    riktning.storre(sigma,matningar,nollH,alfa,sannaVardet,x,y)
  }else if(riktning == '!='){
    riktning.skillt(sigma,matningar,nollH,alfa,sannaVardet,x,y)
  }else{
    print('Maste valja en riktningarna: <, > eller !=')
  }
}


get.x.lab = function(sigma,matningar,nollH,alfa,riktning){
  if(riktning == '<'){
    return (sprintf('k = %4.3f',nollH - qnorm(1 - alfa,0,1)*sigma/sqrt(matningar)))
  }else if(riktning == '>'){
    return (sprintf('k = %4.3f',nollH + qnorm(1 - alfa,0,1)*sigma/sqrt(matningar)))
  }else if(riktning == '!='){
    return (sprintf('k1 = %4.3f, k2 = %4.3f',nollH - qnorm(1 - alfa/2,0,1)*sigma/sqrt(matningar),nollH + qnorm(1 - alfa/2,0,1)*sigma/sqrt(matningar)))
  }
}

riktning.mindre = function(sigma,matningar,nollH,alfa,sannaVardet,x,y){
  k = nollH - qnorm(1 - alfa,0,1)*sigma/sqrt(matningar)
  plot.area(nollH - 3.5*sigma/sqrt(matningar),k,nollH,sigma,matningar,'red')
  text.fig1(sigma,matningar,nollH,'<',alfa)
  
  # fig 2
  
  # Text
  beta = 1-pnorm(k,sannaVardet,sigma/sqrt(matningar))
  xlabel = sprintf('rod = alfa = %4.3f  bla = beta = %4.3f  Styrka(%4.2f) = 1 - beta = %4.3f',alfa,beta,sannaVardet,1-beta)
  
  
  xs = seq(sannaVardet - 3.5*sigma/sqrt(matningar),sannaVardet + 3.5*sigma/sqrt(matningar),length.out = 1000)
  ys = dnorm(xs,sannaVardet,sigma/sqrt(matningar))
  yc = cbind(ys,y)
  xc = cbind(xs,x)
  matplot(xc,yc,type = 'l',col = 'red',main = sprintf('Fel av typ 1 och typ 2; Styrka da mu = %4.2f',sannaVardet),xlab = xlabel,ylab = '')
  plot.area(k,sannaVardet + 3.5*sigma/sqrt(matningar),sannaVardet,sigma,matningar,'blue')
  plot.area(nollH - 3.5*sigma/sqrt(matningar),k,nollH,sigma,matningar,'red')
  points(c(nollH,nollH),c(0,dnorm(nollH,nollH,sigma/sqrt(matningar))),type = 'l',lty = 2)
  points(c(sannaVardet,sannaVardet),c(0,dnorm(sannaVardet,sannaVardet,sigma/sqrt(matningar))),type = 'l',lty = 2)
  
  
}



riktning.storre = function(sigma,matningar,nollH,alfa,sannaVardet,x,y){
  k = nollH + qnorm(1 - alfa,0,1)*sigma/sqrt(matningar)
  plot.area(k,nollH + 3.5*sigma/sqrt(matningar),nollH,sigma,matningar,'red')
  text.fig1(sigma,matningar,nollH,'>',alfa)
  
  # fig 2
  
  # Text
  beta = pnorm(k,sannaVardet,sigma/sqrt(matningar))
  xlabel = sprintf('r?d = alfa = %4.3f  bl? = beta = %4.3f  Styrka(%4.2f) = 1 - beta = %4.3f',alfa,beta,sannaVardet,1-beta)
  
  
  
  xs = seq(sannaVardet - 3.5*sigma/sqrt(matningar),sannaVardet + 3.5*sigma/sqrt(matningar),length.out = 1000)
  ys = dnorm(xs,sannaVardet,sigma/sqrt(matningar))
  yc = cbind(ys,y)
  xc = cbind(xs,x)
  matplot(xc,yc,type = 'l',main = sprintf('Fel av typ 1 och typ 2; Styrka d? mu = %4.2f',sannaVardet), xlab = xlabel,ylab = '')
  plot.area(sannaVardet - 3.5*sigma/sqrt(matningar),k,sannaVardet,sigma,matningar,'blue')
  plot.area(k,nollH + 3.5*sigma/sqrt(matningar),nollH,sigma,matningar,'red')
  points(c(nollH,nollH),c(0,dnorm(nollH,nollH,sigma/sqrt(matningar))),type = 'l',lty = 2)
  points(c(sannaVardet,sannaVardet),c(0,dnorm(sannaVardet,sannaVardet,sigma/sqrt(matningar))),type = 'l',lty = 2)
  
}

riktning.skillt = function(sigma,matningar,nollH,alfa,sannaVardet,x,y){
  kp = nollH + qnorm(1 - alfa/2,0,1)*sigma/sqrt(matningar)
  km = nollH - qnorm(1 - alfa/2,0,1)*sigma/sqrt(matningar)
  plot.area(nollH - 3.5*sigma/sqrt(matningar), km ,nollH, sigma , matningar ,'red')
  plot.area(kp,nollH + 3.5*sigma/sqrt(matningar),nollH ,sigma,matningar,'red')
  text.fig1(sigma,matningar,nollH,'!=',alfa)
  
  # fig 2
  
  # Text
  beta = pnorm(kp,sannaVardet,sigma/sqrt(matningar)) - pnorm(km,sannaVardet,sigma/sqrt(matningar)) 
  xlabel = sprintf('r?d = alfa = %4.3f  bl? = beta = %4.3f  Styrka(%4.2f) = 1 - beta = %4.3f',alfa,beta,sannaVardet,1-beta)
  
  
  
  xs = seq(sannaVardet - 3.5*sigma/sqrt(matningar),sannaVardet + 3.5*sigma/sqrt(matningar),length.out = 1000)
  ys = dnorm(xs,sannaVardet,sigma/sqrt(matningar))
  yc = cbind(ys,y)
  xc = cbind(xs,x)
  matplot(xc,yc,type = 'l',main = sprintf('Fel av typ 1 och typ 2; Styrka d? mu = %4.2f',sannaVardet),xlab = xlabel,ylab = '')
  plot.area(nollH - 3.5*sigma/sqrt(matningar), km ,nollH, sigma , matningar ,'red')
  plot.area(kp,nollH + 3.5*sigma/sqrt(matningar),nollH ,sigma,matningar,'red')
  plot.area(km,kp,sannaVardet,sigma,matningar,'blue')
  points(c(nollH,nollH),c(0,dnorm(nollH,nollH,sigma/sqrt(matningar))),type = 'l',lty = 2,col = 'red')
  points(c(sannaVardet,sannaVardet),c(0,dnorm(sannaVardet,sannaVardet,sigma/sqrt(matningar))),type = 'l',lty = 2)
  
  
  
}

plot.area = function(from,to,nollH,sigma,matningar,colour){
  cord.x = c(from,seq(from,to,length.out = 1000),to)
  cord.y = c(0,dnorm(seq(from,to,length.out = 1000),nollH,sigma/sqrt(matningar)),0)
  polygon(cord.x,cord.y,col = colour)
}

text.fig1 = function(sigma,matningar,nollH,riktning,alfa){
  text(nollH+2*sigma/sqrt(matningar),0.9*dnorm(nollH,nollH,sigma/sqrt(matningar)),sprintf('H0: mu = %4.1f\n',nollH))
  text(nollH+2*sigma/sqrt(matningar),0.78*dnorm(nollH,nollH,sigma/sqrt(matningar)),sprintf('H1: mu %s %4.1f\n',riktning,nollH))
  text(nollH+2*sigma/sqrt(matningar),0.73*dnorm(nollH,nollH,sigma/sqrt(matningar)),sprintf('n = %4.0f',matningar))
  text(nollH+2*sigma/sqrt(matningar),0.55*dnorm(nollH,nollH,sigma/sqrt(matningar)),sprintf('sigma = %4.2f\n',sigma))
  text(nollH+2*sigma/sqrt(matningar),0.4*dnorm(nollH,nollH,sigma/sqrt(matningar)),sprintf('alfa = %4.3f\n',alfa))
}

text.fig1.2 = function(sigma,matningar,nollH,riktning,alfa, medel){
  text(nollH-3*sigma/sqrt(matningar),0.9*dnorm(nollH,nollH,sigma/sqrt(matningar)),sprintf('H0: mu = %4.1f\n',nollH))
  text(nollH-3*sigma/sqrt(matningar),0.78*dnorm(nollH,nollH,sigma/sqrt(matningar)),sprintf('H1: mu %s %4.1f\n',riktning,nollH))
  text(nollH-3*sigma/sqrt(matningar),0.73*dnorm(nollH,nollH,sigma/sqrt(matningar)),sprintf('n = %4.0f',matningar))
  text(nollH-3*sigma/sqrt(matningar),0.55*dnorm(nollH,nollH,sigma/sqrt(matningar)),sprintf('sigma = %4.2f\n',sigma))
  text(nollH-3*sigma/sqrt(matningar),0.4*dnorm(nollH,nollH,sigma/sqrt(matningar)),sprintf('alfa = %4.3f\n',alfa))
  
  text(nollH+3*sigma/sqrt(matningar),0.9*dnorm(nollH,nollH,sigma/sqrt(matningar)),sprintf('medel = %4.1f\n',medel))
  if(riktning == '!='){
    if(abs(medel - nollH) < qnorm(1 - alfa/2,0,1)*sigma/sqrt(matningar)){
      text(nollH+2.5*sigma/sqrt(matningar),0.78*dnorm(nollH,nollH,sigma/sqrt(matningar)),sprintf('k1 < medel < k2 <==>'))
      text(nollH+2.5*sigma/sqrt(matningar),0.65*dnorm(nollH,nollH,sigma/sqrt(matningar)),sprintf('H0 forkastas ej'))
    }else{
      if(medel < qnorm(1 - alfa/2,0,1)*sigma/sqrt(matningar)){
        text(nollH+2.5*sigma/sqrt(matningar),0.78*dnorm(nollH,nollH,sigma/sqrt(matningar)),sprintf('medel < k1 <==>'))
        text(nollH+2.5*sigma/sqrt(matningar),0.65*dnorm(nollH,nollH,sigma/sqrt(matningar)),sprintf('H0 forkastas'))
      }else{
        text(nollH+2.5*sigma/sqrt(matningar),0.78*dnorm(nollH,nollH,sigma/sqrt(matningar)),sprintf('medel  > k2 <==>'))
        text(nollH+2.5*sigma/sqrt(matningar),0.65*dnorm(nollH,nollH,sigma/sqrt(matningar)),sprintf('H0 forkastas'))
      }
    }
  }
  if(riktning == '<'){
    if(medel < nollH - qnorm(1 - alfa,0,1)*sigma/sqrt(matningar)){
      text(nollH+2.5*sigma/sqrt(matningar),0.78*dnorm(nollH,nollH,sigma/sqrt(matningar)),sprintf('medel  < k <==>'))
      text(nollH+2.5*sigma/sqrt(matningar),0.65*dnorm(nollH,nollH,sigma/sqrt(matningar)),sprintf('H0 forkastas'))
    }else{
      text(nollH+2.5*sigma/sqrt(matningar),0.78*dnorm(nollH,nollH,sigma/sqrt(matningar)),sprintf('medel  > k <==>'))
      text(nollH+2.5*sigma/sqrt(matningar),0.65*dnorm(nollH,nollH,sigma/sqrt(matningar)),sprintf('H0 forkastas ej'))
    }
  }
  if(riktning == '>'){
    if(medel > nollH + qnorm(1 - alfa,0,1)*sigma/sqrt(matningar)){
      text(nollH+2.5*sigma/sqrt(matningar),0.78*dnorm(nollH,nollH,sigma/sqrt(matningar)),sprintf('medel  > k <==>'))
      text(nollH+2.5*sigma/sqrt(matningar),0.65*dnorm(nollH,nollH,sigma/sqrt(matningar)),sprintf('H0 forkastas'))
    }else{
      text(nollH+2.5*sigma/sqrt(matningar),0.78*dnorm(nollH,nollH,sigma/sqrt(matningar)),sprintf('medel  < k <==>'))
      text(nollH+2.5*sigma/sqrt(matningar),0.65*dnorm(nollH,nollH,sigma/sqrt(matningar)),sprintf('H0 forkastas ej'))
    }
  } 
}

text.fig2.2 = function(sigma,matningar,nollH,riktning,alfa, medel){
  if(riktning == '<'){
    p = pnorm(medel,nollH,sigma/sqrt(matningar))
  }else if(riktning == '>'){
     p = 1 - pnorm(medel,nollH,sigma/sqrt(matningar))
  }else{
    if(medel <= nollH){
      p = pnorm(medel,nollH,sigma/sqrt(matningar))
    }else{
      p = 1 - pnorm(medel,nollH,sigma/sqrt(matningar))
    }
  }
  text(nollH+3*sigma/sqrt(matningar),0.8*dnorm(nollH,nollH,sigma/sqrt(matningar)),sprintf('medel = %4.1f',medel))
  if(riktning == '!='){
    text(nollH+3*sigma/sqrt(matningar),0.6*dnorm(nollH,nollH,sigma/sqrt(matningar)),sprintf('P = 2 * %4.3f\n',p))
    p = 2*p
  }else{
    text(nollH+3*sigma/sqrt(matningar),0.6*dnorm(nollH,nollH,sigma/sqrt(matningar)),sprintf('P = %4.3f\n',p))
  }
  if(p > alfa){
    text(nollH+3*sigma/sqrt(matningar),0.55*dnorm(nollH,nollH,sigma/sqrt(matningar)),sprintf('P > alfa <==>'))
    text(nollH+3*sigma/sqrt(matningar),0.4*dnorm(nollH,nollH,sigma/sqrt(matningar)),sprintf('H0 forkastas ej'))
  }else{
    text(nollH+3*sigma/sqrt(matningar),0.55*dnorm(nollH,nollH,sigma/sqrt(matningar)),sprintf('P < alfa <==>'))
    text(nollH+3*sigma/sqrt(matningar),0.4*dnorm(nollH,nollH,sigma/sqrt(matningar)),sprintf('H0 forkastas'))
  }
}





