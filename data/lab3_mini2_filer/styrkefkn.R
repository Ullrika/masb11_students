styrka1 = function(sigma,matningar,nollH,alfa,riktning,sannaVardet){
#  windows()
  layout(matrix(c(1,2),2,1))
  fig2(sigma,matningar,nollH,alfa,riktning,sannaVardet)
  text(nollH-4*sigma/sqrt(matningar),0.9*dnorm(nollH,nollH,sigma/sqrt(matningar)),sprintf('H0: mu = %4.1f\n',nollH))
  text(nollH-4*sigma/sqrt(matningar),0.78*dnorm(nollH,nollH,sigma/sqrt(matningar)),sprintf('H1: mu %s %4.1f\n',riktning,nollH))
  text(nollH-4*sigma/sqrt(matningar),0.73*dnorm(nollH,nollH,sigma/sqrt(matningar)),sprintf('n = %4.0f',matningar))
  text(nollH-4*sigma/sqrt(matningar),0.55*dnorm(nollH,nollH,sigma/sqrt(matningar)),sprintf('sigma = %4.2f\n',sigma))
  text(nollH-4*sigma/sqrt(matningar),0.4*dnorm(nollH,nollH,sigma/sqrt(matningar)),sprintf('alfa = %4.3f\n',alfa))
  
  st(sigma,matningar,nollH,alfa,riktning,sannaVardet,TRUE)
}

styrka2 = function(sigma,matningar,nollH,alfa,riktning){
 # windows()
  layout(1,1)
  st(sigma,matningar,nollH,alfa,riktning,0,FALSE)
}


fig2 = function(sigma,matningar,nollH,alfa,riktning,sannaVardet){
  x = seq(nollH - 3.5*sigma/sqrt(matningar),nollH + 3.5*sigma/sqrt(matningar),length.out = 1000)
  y = dnorm(x,nollH,sigma/sqrt(matningar))
  if(riktning == '<'){
    k = nollH - qnorm(1 - alfa,0,1)*sigma/sqrt(matningar)
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
    
  }else if(riktning == '>'){
    k = nollH + qnorm(1 - alfa,0,1)*sigma/sqrt(matningar)
    beta = pnorm(k,sannaVardet,sigma/sqrt(matningar))
    xlabel = sprintf('rod = alfa = %4.3f  bla = beta = %4.3f  Styrka(%4.2f) = 1 - beta = %4.3f',alfa,beta,sannaVardet,1-beta)
    
    
    
    xs = seq(sannaVardet - 3.5*sigma/sqrt(matningar),sannaVardet + 3.5*sigma/sqrt(matningar),length.out = 1000)
    ys = dnorm(xs,sannaVardet,sigma/sqrt(matningar))
    yc = cbind(ys,y)
    xc = cbind(xs,x)
    matplot(xc,yc,type = 'l',main = sprintf('Fel av typ 1 och typ 2; Styrka da mu = %4.2f',sannaVardet), xlab = xlabel,ylab = '')
    plot.area(sannaVardet - 3.5*sigma/sqrt(matningar),k,sannaVardet,sigma,matningar,'blue')
    plot.area(k,nollH + 3.5*sigma/sqrt(matningar),nollH,sigma,matningar,'red')
    points(c(nollH,nollH),c(0,dnorm(nollH,nollH,sigma/sqrt(matningar))),type = 'l',lty = 2)
    points(c(sannaVardet,sannaVardet),c(0,dnorm(sannaVardet,sannaVardet,sigma/sqrt(matningar))),type = 'l',lty = 2)
    
  }else{
    kp = nollH + qnorm(1 - alfa/2,0,1)*sigma/sqrt(matningar)
    km = nollH - qnorm(1 - alfa/2,0,1)*sigma/sqrt(matningar)
    beta = pnorm(kp,sannaVardet,sigma/sqrt(matningar)) - pnorm(km,sannaVardet,sigma/sqrt(matningar)) 
    xlabel = sprintf('rod = alfa = %4.3f  bla = beta = %4.3f  Styrka(%4.2f) = 1 - beta = %4.3f',alfa,beta,sannaVardet,1-beta)
    
    
    
    xs = seq(sannaVardet - 3.5*sigma/sqrt(matningar),sannaVardet + 3.5*sigma/sqrt(matningar),length.out = 1000)
    ys = dnorm(xs,sannaVardet,sigma/sqrt(matningar))
    yc = cbind(ys,y)
    xc = cbind(xs,x)
    matplot(xc,yc,type = 'l',main = sprintf('Fel av typ 1 och typ 2; Styrka da mu = %4.2f',sannaVardet),xlab = xlabel,ylab = '')
    plot.area(nollH - 3.5*sigma/sqrt(matningar), km ,nollH, sigma , matningar ,'red')
    plot.area(kp,nollH + 3.5*sigma/sqrt(matningar),nollH ,sigma,matningar,'red')
    plot.area(km,kp,sannaVardet,sigma,matningar,'blue')
    points(c(nollH,nollH),c(0,dnorm(nollH,nollH,sigma/sqrt(matningar))),type = 'l',lty = 2,col = 'red')
    points(c(sannaVardet,sannaVardet),c(0,dnorm(sannaVardet,sannaVardet,sigma/sqrt(matningar))),type = 'l',lty = 2)
    
  }
}


st = function(sigma,n,nollH,alfa,riktning,sannaVardet,line){
  if(riktning == '<'){
    k = nollH - qnorm(1 - alfa,0,1)*sigma/sqrt(n)
    x = seq(k - 3*sigma/sqrt(n),k + 2*sigma/sqrt(n),length.out = 1000)
    plot(x,pnorm(k,x,sigma/sqrt(n)),type = 'l',ylab = '',xlab = 'mu',main = sprintf('S(mu); Sannolikheten att forkasta H0: mu = %3.2f',nollH))
    grid()
  }else if(riktning == '>'){
    k = nollH + qnorm(1 - alfa,0,1)*sigma/sqrt(n)
    x = seq(k - 2*sigma/sqrt(n),k + 3*sigma/sqrt(n),length.out = 1000)
    plot(x,1 - pnorm(k,x,sigma/sqrt(n)), type = 'l',ylab = '',xlab = 'mu',main = sprintf('S(mu); Sannolikheten att forkasta H0: mu = %3.2f',nollH))
    grid()
  }else if(riktning == '!='){
    k1 = nollH - qnorm(1-alfa/2,0,1)*sigma/sqrt(n)
    k2 = nollH + qnorm(1-alfa/2,0,1)*sigma/sqrt(n)
    x = seq(k1 - 3.5*sigma/sqrt(n),k2 + 3.5*sigma/sqrt(n),length.out = 1000)
    plot(x,1-pnorm(k2,x,sigma/sqrt(n)) + pnorm(k1,x,sigma/sqrt(n)),type = 'l',ylab = '',xlab = 'mu',main = sprintf('S(mu); Sannolikheten att forkasta H0: mu = %3.2f',nollH))
    grid()
  }
  if(line){
    if(riktning == '<'){
      f = pnorm(k,sannaVardet,sigma/sqrt(n))
      points(c(sannaVardet,sannaVardet),c(0,f),type = 'l',col = 'red')
    }else if(riktning == '>'){
      f = 1 - pnorm(k,sannaVardet,sigma/sqrt(n))
      points(c(sannaVardet,sannaVardet),c(0,f),type = 'l',col = 'red')
    }else if(riktning == '!='){
      f = 1-pnorm(k2,sannaVardet,sigma/sqrt(n)) + pnorm(k1,sannaVardet,sigma/sqrt(n))
      points(c(sannaVardet,sannaVardet),c(0,f),type = 'l',col = 'red')
    }
    if(riktning == '!='){
      k = k1;
    }
    points(c(k - 4*sigma/sqrt(n),sannaVardet),c(f,f),type = 'l',col = 'red')
  }
}

plot.area = function(from,to,nollH,sigma,matningar,colour){
  cord.x = c(from,seq(from,to,length.out = 1000),to)
  cord.y = c(0,dnorm(seq(from,to,length.out = 1000),nollH,sigma/sqrt(matningar)),0)
  polygon(cord.x,cord.y,col = colour)
}

