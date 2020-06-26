rednoise<-function(phi, len, media = 0, sigma = 1)
{
    ruido<-rnorm(n = len,mean = media, sd = sigma )
    for( i in  1:(len-1) ) {x[i+1]<- phi*x[i]+ruidobranco[i]}
    return(x)
}

rednoise<-function(phi, len, media = 0, sigma = 1)
{
    
    ruido<-rnorm(n = len,mean = media, sd = sigma )
    for( i in  1:999 ) {
        
        x[i+1]<- phi*x[i]+ruidobranco[i]
        
        x[i+1]<- phi*x[i]+ruidobranco[i]
        
    }
    return(x)
    
}



sidemean<-function(matrix, axis = 1) {
    
    sidesize <- dim(matrix)[axis]
    
    finalmean <- numeric(sidesize)
    
    for(i in 1:sidesize)
    {
        if(axis==1)
        {
            finalmean[i]<- mean(matrix[i,])
        }
        else
        {
            finalmean[i]<- mean(matrix[,i])
        }
    }
    
    return(finalmean)
}

lseq <- function(from = 1, to = 100000, length.out = 6) {
    # logarithmic spaced sequence
    # blatantly stolen from library("emdbook"), because needed only this
    exp(seq(log(from), log(to), length.out = length.out))
}

testelele<- function(h){
    print("Hello")
}

Entropia3D<- function(h){
    
    s<- - h$density * log(h$density) * ( h$breaksx[2] - h$breaksx[1] ) * ( h$breaksy[2] - h$breaksy[1] )
    
    s[which(!is.finite(s))]<- 0
    
    shannon<- sum(s)
    
    
    return(shannon)
    
}

incrementos<-function(x) {
    
    x0<-head(x,length(x)-1)
    x1<-tail(x,length(x)-1)
    
    v<-(x1-x0)
    v0<-head(v,length(x)-2)
    v1<-tail(v,length(x)-2)
    
    a<-(v1-v0)
    a0<-head(a,length(x)-3)
    a1<-tail(a,length(x)-3)
    
    da<-(a1-a0)
    da1<-tail(da,length(x)-4)
    da0<-head(da,length(x)-4)
    
    dda<-(da1-da0)
    dda1<-tail(dda,length(x)-5)
    dda0<-head(dda,length(x)-5)
    
    ddda<-(dda1-dda0)
    
    
    
    final<-list("v"= v, "a"= a, "da" = da, "dda" = dda, "ddda" = ddda )
    
    return(final)
}

ShannonS<-function(h)
{
    ss<- - h$density*log(h$density)
    ss[which(!is.finite(ss))]<- 0
    s<-sum(ss)*(1/sum(h$density))
    return(s)
}


FuncaoResposta<-function(qero, Dt, vez){
    
    corrpre <- paste(toString(vez) ,'Bcorrpre', sep="")
    rresposta <- paste(toString(vez) ,'Bresposta', sep="")
    csv<-'.csv'
    
    titulo0<-paste( corrpre , csv , sep="")
    tituloR<-paste( rresposta , csv , sep="")
    
    ns<- qero$ns
    
    kick<- qero$kick
    tkick<- qero$tkick/(qero$dt*qero$nc)
    
    fim<- numeric(Dt)
    
    #mediapre<-   mean( qero$posicao[1:tkick,]  )
    
    corpre<-      (  funcaodecorrelacaoSample(qero$posicao, (tkick - Dt - 1) , Dt , ns)   )
    corpos<-      (  funcaodecorrelacaoSample(qero$posicao, tkick , Dt , ns)       )
    
    
    for (i in 1:Dt)
    {
        mediapos<-  mean(qero$posicao[tkick+i, ]  )
        fim[i]<- mediapos # - mediapre -> mediapre é igual a zero
    }
    
    fim<-fim/kick
    
    SalvaCorrelacao(corpre, titulo0)
    SalvaCorrelacao(fim, tituloR)
    
    final <- list("corpre"= corpre , "resposta" = fim )
    return(final)
    
    
}

#nc = ncoarse

KurchanSilvio<-function(ns, nt, nc = 1, dt = 0.001, massa = 1, gamma, ctelastica, temp = 3 , omega0 = 50 , aa = 0.001966 ){
    
    ntrans <- 10
    
    #lt1 <- 10000000
    
    lt1 <- 1000
    
    dtc <- dt * nc
    
    xx <- matrix(ncol = ns, nrow = nt)
    vv <- matrix(ncol = ns, nrow = nt)
    eo <- matrix(ncol = ns, nrow = nt)
    #ez <- matrix(ncol = ns, nrow = nt)
    
    nsamples<-ns
    ncoarse<- nc
    
    pi <-  3.14159265358979323846264338327950288
    we <- massa
    g <- gamma
    gw <- g / we
    gt <- g * dt
    gtw <- g * dt / we
    gmt <- gt * we
    cl <- 0.1
    c1 <- cl / we
    cm <- 1.0
    c2 <- cm / we
    
    #sigw = 0.0
    tura <- temp
    sigw <- 2.0 * g * tura * dt / we / we
    
    
    ncictrans <- 1000 * ( 1.0 / g / dtc )
    
    
    emts <- 0.0
    emtb <- 0.0
    emtf <- 0.0
    
    
    for(i in 1:nsamples)
    {
        
        print(i)
        print('nsamples')
        
        ww0 <- sqrt(sigw) * rnorm(1)
        x1 <- 0.0
        v1 <- 0.0
        y1 <- 0.1
        u1 <- 0.0
        
        for(jc in 1:ncictrans)
        {
            
            #print(jc)
            #print('ncictrans')
            
            
            for(ic in 1:ncoarse)
            {
                
                
                #print(ic)
                #print('ncoarse')
                
                x0 <- x1
                v0 <- v1
                y0 <- y1
                u0 <- u1
                
                ww <- sqrt(sigw) * rnorm(1)
                
                wwn <- (ww + ww0) / 2.0
                ww0 <- ww
                
                frc <- 0.0
                frc <- frc - c2 * x0
                frc <- frc - gw * v0
                
                trm1 <- 0.0
                trm1 <- dt * frc
                v1 <- v0 + trm1 + wwn
                x1 <- x0 + v1 * dt
                
                emts <- (we*v1*v1 + cm*x1*x1)/2.0
                
                frc <- 0.0
                frc <- - cl * y0 + aa * emts
                trm1 <- dt * frc
                u1 <- u0 + trm1
                y1 <- y0 + (u1+u0) * dt / 2.0
                
                emtb <- (u1*u1 + cl*y1*y1)/2.0
                
            }
            
            
            
        }
        
        #fim do transiente
        
        for(it in 1:lt1)
        {
            for(ic in 1:ncoarse)
            {
                x0 <- x1
                v0 <- v1
                y0 <- y1
                u0 <- u1
                
                sigw <- 0.0
                sigw <- 2.0 * g * tura * dt / we / we
                ww <- sqrt(sigw)  * rnorm(1)   #* gasdev(idum)
                
                wwn <- (ww + ww0) / 2.0
                ww0 <- ww
                
                frc <- 0.0
                frc <- frc - c2 * x0
                frc <- frc - gw * v0
                
                trm1 <- 0.0
                trm1 <- dt * frc
                v1 <- v0 + trm1 + wwn
                
                x1 <- x0 + v1 * dt
                
                emts <- (we*v1*v1 + cm*x1*x1)/2.0
                
                
                frc <- 0.0
                
                ####
                
                frc <- - cl * y0 + aa * emts
                
                trm1 <- dt * frc
                
                u1 <- u0 + trm1
                
                y1 <- y0 + (u1+u0) * dt / 2.0
                emtb <- (u1*u1 + cl*y1*y1)/2.0
                
            }
            
            emtf <- emtf + emtb / (lt1) / (nsamples)
            
            vv[it,i] <- u1
            xx[it,i] <- x1
            eo[it,i] <- emtb
            
        }
        
    }
    
    ez <- emtf
    
    final<-list("v"= vv, "x"= xx, "eo" = eo , "ez" = ez )
    
    return(final)
}

SimulacaoBimodalUOUOUOUOUOU <- function(ns, nt, nc, dt, massa, gamma, ctelastica, alfa, dalfa, xa, dxa){
    
    
    ntc<-nt*nc
    
    xx<-matrix( nrow = nt, ncol = ns)
    vv<-matrix( nrow = nt, ncol = ns)
    ee<-matrix( nrow = nt, ncol = ns)
    
    ########
    
    ji<-matrix( nrow = nt, ncol = ns)
    jd<-matrix( nrow = nt, ncol = ns)
    ruido<-matrix( nrow = nt, ncol = ns)
    
    alfazero<-alfa
    
    ua<- alfa/2
    
    ub<- ua
    
    c1<- ctelastica/massa
    
    gm<- gamma/massa
    
    ie<- 0
    
    it<- 0
    
    
    
    for(i in 1:ns)
    {
        
        print(i)
        #Temp<- ((xa^2)*(massa * alfa))/(gamma * (ctelastica + alfa * (massa * alfa + gamma)  ))
        
        x <- 0
        v <- 0
        
        v12 <- 0
        ji1 <- 0
        jd1 <- 0
        alfa <- alfa + dalfa
        xa <- xa + dxa
        
        xb <- -xa
        ua <-alfa/2
        ub <-ua
        eta <-xa
        
        it <- 0
        ic <- 0
        ir <- 0
        ig <- 0
        ie <- 0
        
        while(it < ntc)
        {
            if(eta==xa)
            {
                
                teta<- - log(1-runif(1))/ub
                neta<- teta%/%dt + 1
                
            } else {
                
                teta<- - log(1-runif(1))/ua
                neta<- teta%/%dt + 1
                
            }
            
            ie<- 0
            
            while(ie<neta & it < ntc)
            {
                it<-it+1
                ie<-ie+1
                
                u<-v
                y<-x
                force<-  -gm*u -c1*y + eta/massa
                v<- u + force*dt
                x<- y+ v*dt
                
                
                #  u<-v
                #y<-x
                #force<-  -gm*u -c1*y + eta/massa
                #v<- u + force*dt
                #x<- y+ v*dt
                
                ji1 <- ji1 + (v+u)*eta*dt/2
                #jd1<-jd1+ gamma*(v^2+u^2)*dt/2
                jd1 <- jd1 + gamma*(((v+u)/2)^2)*dt   #!!!!!
                
                ig<-ig+1
                
                if(ig==nc)
                {
                    
                    xx[it/nc,i]<- x
                    vv[it/nc,i]<- v
                    ee[it/nc,i]<- (ji1-jd1)
                    ji[it/nc,i]<- ji1
                    jd[it/nc,i]<- jd1
                    ruido[it/nc,i]<- eta
                    
                    
                    ig<-0
                }
            }
            if(eta==xa)
            {
                eta<-xb
            } else {
                eta<-xa
            }
        }
        
    }
    
    final<-list( "posicao"=xx, "velocidade"=vv, "ji" = ji , "jd" = jd ,"energia"=ee, "massa"=massa,"gamma"=gamma,"k"=ctelastica, "alfa" = alfazero,"dalfa"= dalfa, "xa"=xa, "dxa"=dxa, "ns" = ns ,"nt"=nt,"nc"=nc,"dt"= dt, "temp"=Temp, "ruido" = ruido)
    
    return(final)
}




KurchanSilvioBIM<-function(ns, nt, nc = 1, dt = 0.000001, massa = 1, gamma, ctelastica, temp = 3 , omega0 = 50 , aa = 0.001966 ){
    
    ntrans <- 10
    #lt1 <- 10000000
    lt1 <- 1000
    
    dtc<- dt* nc
    
    xx <- matrix(ncol = ns, nrow = nt)
    vv <- matrix(ncol = ns, nrow = nt)
    eo <- matrix(ncol = ns, nrow = nt)
    ez <- matrix(ncol = ns, nrow = nt)
    
    nsamples <- ns
    ncoarse <- nc
    
    
    pi <-  3.14159265358979323846264338327950288
    we <- massa
    g <- gamma
    gw <- g / we
    gt <- g * dt
    gtw <- g * dt / we
    gmt <- gt * we
    cl <- 0.1
    c1 <- cl / we
    cm <- 1.0
    c2 <- cm / we
    
    #sigw = 0.0
    tura <- temp
    sigw <- 2.0 * g * tura * dt / we / we
    
    
    ncictrans <- 1000 * ( 1.0 / g / dtc )
    
    
    emts <- 0.0
    emtb <- 0.0
    emtf <- 0.0
    
    
    for(i in 1:nsamples)
    {
        ww0 <- sqrt(sigw) * rnorm(1)
        x1 <- 0.0
        v1 <- 0.0
        y1 <- 0.1
        u1 <- 0.0
        
        alfa<- alfa+ dalfa
        xa<- xa+dxa
        
        xb<- -xa
        ua<- alfa/2
        ub<- ua
        eta<- xa
        
        it<-0
        ic<-0
        ir<-0
        ig<-0
        ie<-0
        
        
        for(jc in 1:ncictrans)
        {
            
            for(ic in 1:ncoarse)
            {
                
                if(eta==xa & ie>neta)
                {
                    eta <- xb
                    
                    teta<- - log(1-runif(1))/ub
                    neta<- teta%/%dt + 1
                    
                    ie <- 0
                    
                } else if(eta==xb & ie>neta) {
                    
                    eta <- xa
                    
                    teta<- - log(1-runif(1))/ua
                    neta<- teta%/%dt + 1
                    
                    ie <- 0
                    
                }
                
                        #it <- it + 1
                        ie <- ie + 1
                        
                        x0 <- x1
                        v0 <- v1
                        y0 <- y1
                        u0 <- u1
                        
                        #force <-  -gm*u -c1*y + eta/massa
                        #frc <- 0.0
                        #frc <- frc - c2 * x0
                        #frc <- frc - gw * v0
                        
                        force <-  -gw*u -c2*y + eta/massa
                        
                        # v <- u + force*dt
                        #x <- y+ v*dt
                        #trm1 <- 0.0
                        
                        trm1 <- dt * force
                        v1 <- v0 + trm1 #+ wwn
                        x1 <- x0 + v1 * dt
                        
                        emts <- (we*v1*v1 + cm*x1*x1)/2.0
                        
                        #frc <- 0.0
                        frc <- - cl * y0 + aa * emts
                        trm1 <- dt * frc
                        u1 <- u0 + trm1
                        y1 <- y0 + (u1+u0) * dt / 2.0
                        
                        emtb <- (u1*u1 + cl*y1*y1)/2.0
                        
                        #ji1 <- ji1 + (v1+v0)*eta*dt/2
                        #jd1 <- jd1 + gamma*(((v1+v0)/2)^2)*dt   #!!!!!
                        
                        #ig<-ig+1
                        
                        emtb <- (u1*u1 + cl*y1*y1)/2.0
                
                
                
                
                #ww <- sqrt(sigw) * rnorm(1)
                #wwn <- (ww + ww0) / 2.0
                #ww0 <- ww
                #trm1 <- 0.0
                #trm1 <- dt * frc
                #v1 <- v0 + trm1 + wwn
                #x1 <- x0 + v1 * dt
                #emts <- (we*v1*v1 + cm*x1*x1)/2.0
                #frc <- 0.0
                #frc <- - cl * y0 + aa * emts
                #trm1 <- dt * frc
                #u1 <- u0 + trm1
                #y1 <- y0 + (u1+u0) * dt / 2.0
                
            }
            
        }
        
        #fim do transiente
        
        for(it in 1:lt1)
        {
            for(ic in 1:ncoarse)
            {
                if(eta==xa & ie>neta)
                {
                    eta <- xb
                    
                    teta<- - log(1-runif(1))/ub
                    neta<- teta%/%dt + 1
                    
                    ie <- 0
                    
                } else if(eta==xb & ie>neta) {
                    
                    eta <- xa
                    
                    teta<- - log(1-runif(1))/ua
                    neta<- teta%/%dt + 1
                    
                    ie <- 0
                    
                }
                
                
                #it <- it + 1
                ie <- ie + 1
                
                #u <- v
                #y <- x
                
                x0 <- x1
                v0 <- v1
                y0 <- y1
                u0 <- u1
                
                
                #force <-  -gm*u -c1*y + eta/massa
                #frc <- 0.0
                #frc <- frc - c2 * x0
                #frc <- frc - gw * v0
                
                force <-  -gw*u -c2*y + eta/massa
                
                
                # v <- u + force*dt
                #x <- y+ v*dt
                
                #trm1 <- 0.0
                
                trm1 <- dt * force
                v1 <- v0 + trm1 #+ wwn
                x1 <- x0 + v1 * dt
                
                emts <- (we*v1*v1 + cm*x1*x1)/2.0
                
                #frc <- 0.0
                frc <- - cl * y0 + aa * emts
                trm1 <- dt * frc
                u1 <- u0 + trm1
                y1 <- y0 + (u1+u0) * dt / 2.0
                
                emtb <- (u1*u1 + cl*y1*y1)/2.0
                
                #ji1 <- ji1 + (v1+v0)*eta*dt/2
                #jd1 <- jd1 + gamma*(((v1+v0)/2)^2)*dt   #!!!!!
                
                #ig<-ig+1
                
                emtb <- (u1*u1 + cl*y1*y1)/2.0
                
                
                
                
                #ww <- sqrt(sigw) * rnorm(1)
                #wwn <- (ww + ww0) / 2.0
                #ww0 <- ww
                #trm1 <- 0.0
                #trm1 <- dt * frc
                #v1 <- v0 + trm1 + wwn
                #x1 <- x0 + v1 * dt
                #emts <- (we*v1*v1 + cm*x1*x1)/2.0
                #frc <- 0.0
                #frc <- - cl * y0 + aa * emts
                #trm1 <- dt * frc
                #u1 <- u0 + trm1
                #y1 <- y0 + (u1+u0) * dt / 2.0
                
            }
            
            emtf <- emtf + emtb / (lt1) / (nsamples)
            
            vv[it,i] <- u1
            xx[it,i] <- x1
            eo[it,i] <- emtb
            
        }
        
    }
    
    ez <- emtf
    
    final<-list("v"= vv, "x"= xx, "eo" = eo )
    
    return(final)
}




Kurchantrans<-function(ns, ntrans, nt, nc=1, dt, massa, gamma, ctelastica, temp , omega0 = 50 , a = 10 ){

# omega0 = 50 , a = 10 ){
    
    #alfa, dalfa, xa, dxa, omega0 = 10 , a = 28 ){
    
    # eh na energia seguir folha
    #<- function(ns, nt, nc, dt, massa, gamma, ctelastica, alfa, dalfa, xa, dxa, deltax, tt){
    
    #qero<-SimulacaoBimodal(ns, nt, nc, dt, massa, gamma, ctelastica, alfa, dalfa, xa, dxa)
    
    tw<-round(nt/10)
    
    qero<-SimulacaoLangevin(ns, nt+ntrans, nc, dt, massa, gamma, ctelastica, temp)

    #wat<-FuncaoResposta(qero,Dt,)
    
    xx <- matrix(ncol = ns, nrow = nt/nc)
    vv <- matrix(ncol = ns, nrow = nt/nc)
    Dx <- matrix(ncol = ns, nrow = nt/nc)
    Dv <- matrix(ncol = ns, nrow = nt/nc)
    #zz <- matrix(ncol = ns, nrow = nt)
    
    
    for(i in 1:ns)
    {
        
        v<- 0
        x<- 0
        
        ic<-0
        
        for(j in 1:(nt+ntrans))
        {
            ic<- ic+1
            
            u<-v
            y<-x
            
            #force<-  -(omega0**2)*x # + a * qero$energia[j,i] # energia de qero é o O no modelo do kurchan
            
            if(j>tw-1)
            {
                force<- -(omega0**2)*x  + a * qero$energia[j,i]
            }
            else
            {
                force<- -(omega0**2)*x
            }
            
            #print(j)
            #print(i)
            
            #v<- u + force*dt
            #x<- y + v*dt
            
            v<- u + force*dt
            x<- y + v*dt
            
            #ji1<-ji1+ (v+u)*eta*dt/2
            #jd1 <- jd1 + gamma*(((v+u)/2)^2)*dt
            
            if(ic==nc && j > ntrans)
            {
                ic <- 0
                
                xx[j,i] <- x
                vv[j,i] <- v
                
                Dx[j,i] <- v*dt
                Dv[j,i] <- force*dt
            
            }
            
        }
    }
    
    Et<-numeric(nt/nc)
    
    Eot<-numeric(nt/nc)
    
    
    for(j in 1:nt/nc)
    {
        
        Et[j]<- (0.5)*(omega0**2) * mean(xx[j,]**2) + (0.5) * mean(vv[j,]**2)
        
        Eot[j]<- - a*mean(qero$energia[j,] * xx[j,])
        
    }
    
    #Energia do oscilador tem que convergir pra temperatura Tzero
    
    E<- (0.5)*(omega0**2)*mean(xx**2) + (0.5) * mean(vv**2)
    
    Eo<- - a*mean(qero$energia * xx)
    
    #    Resposta, começando com Dt arbitrário
    
    Dt<- 100
    
    #resposta <- matrix(nrow = nt, ncol = Dt )
    
    #  resposta <- numeric(Dt)
    
    
    # round(nt/2) é um tempo qualquer onde a principio ja temos estado estacionario
    #
    #  for(j in 1:(nt- Dt ))
    #   {
    #     for(i in 1:Dt)
    #  {
            #resposta[j,i]<-  (  mean(qero$energia[ j +i  , ]) - mean( qero$energia[ j  ,])   )/(  a*Dx[j,2]  )
            
            #       resposta[i]<-  (  mean(qero$energia[ tw +i  , ]) - mean( qero$energia[ tw  ,])   )/(  a* mean(Dx[tw,])  )
            #}
        #  }

    final <- list("Eoscilador"= Et, "Eint" = Eot, "ELangevin" = qero$energia, "ELangemedia" = sidemean(qero$energia)  ,"Eosctot"= E, "Esistot" = Eo,    "posicao"= xx,  "velocidade"= vv, "omega"=omega0, "temp" = temp)
    return(final)
    
}


Kurchan<-function(ns,  nt, nc=1, dt, massa, gamma, ctelastica, temp , omega0 = 50 , a = 10 ){
    

    
    tw<-round(nt/10)
    
    qero<-SimulacaoLangevin(ns, nt*nc, 1, dt, massa, gamma, ctelastica, temp)
    

    print('foi até aqui')
    
    xx <- matrix(ncol = ns, nrow = nt)
    vv <- matrix(ncol = ns, nrow = nt)
    Dx <- matrix(ncol = ns, nrow = nt)
    Dv <- matrix(ncol = ns, nrow = nt)

    print('foi até aqui 2')
    
    for(i in 1:ns)
    {
        
        v<- 0
        x<- 0
        
        ic<-0
        
        for(j in 1:(nt*nc))
        {
            ic<- ic+1
            
            u<-v
            y<-x
            
            #force<-  -(omega0**2)*x # + a * qero$energia[j,i] # energia de qero é o O no modelo do kurchan
            
            if(j>tw-1)
            {
                force<- -(omega0**2)*x  + a * qero$energia[j,i]
            }
            else
            {
                force<- -(omega0**2)*x
            }
            
            #print(j)
            #print(i)
            
            #v<- u + force*dt
            #x<- y + v*dt
            
            v<- u + force*dt
            x<- y + v*dt
            
            #ji1<-ji1+ (v+u)*eta*dt/2
            #jd1 <- jd1 + gamma*(((v+u)/2)^2)*dt
            
            if(ic==nc )
            {
                ic <- 0
                
                xx[j/nc,i] <- x
                vv[j/nc,i] <- v
                
                Dx[j/nc,i] <- v*dt
                Dv[j/nc,i] <- force*dt
                
            }
            
        }
    }
    
    Et<-numeric(nt)
    
    Eot<-numeric(nt)
    print('foi ateh qui')
    
    print(toString(length(xx)))
    print(toString(length(qero$energia)))
    
    for(j in 1:(nt/nc))
    {
        
        Et[j]<- (0.5)*(omega0**2) * mean(xx[j,]**2) + (0.5) * mean(vv[j,]**2)
        
        Eot[j]<- - a*mean(qero$energia[j,] * xx[j,])
        
    }
    
    #Energia do oscilador tem que convergir pra temperatura Tzero
    
    E<- (0.5)*(omega0**2)*mean(xx**2) + (0.5) * mean(vv**2)
    
    Eo<- - a*mean(qero$energia * xx)
    
    #    Resposta, começando com Dt arbitrário
    
    Dt<- 100
    
    #resposta <- matrix(nrow = nt, ncol = Dt )
    
    #  resposta <- numeric(Dt)
    
    
    # round(nt/2) é um tempo qualquer onde a principio ja temos estado estacionario
    #
    #  for(j in 1:(nt- Dt ))
    #   {
    #     for(i in 1:Dt)
    #  {
    #resposta[j,i]<-  (  mean(qero$energia[ j +i  , ]) - mean( qero$energia[ j  ,])   )/(  a*Dx[j,2]  )
    
    #       resposta[i]<-  (  mean(qero$energia[ tw +i  , ]) - mean( qero$energia[ tw  ,])   )/(  a* mean(Dx[tw,])  )
    #}
    #  }
    
    final <- list("Eoscilador"= Et, "Eint" = Eot, "ELangevin" = qero$energia, "ELangemedia" = sidemean(qero$energia)  ,"Eosctot"= E, "Esistot" = Eo,    "posicao"= xx,  "velocidade"= vv, "omega"=omega0, "temp" = temp)
    return(final)
    
}



#SimulacaoLangevin




KurchanLite<-function(qero, omega0 = 0.84 , a = 0.2 ){
    
    #qero<-SimulacaoLangevin(ns, nt*nc, 1, dt, massa, gamma, ctelastica, temp)
    #function(qero , ns,  nt, nc=1, dt, massa, gamma, ctelastica, temp , omega0 = 50 , a = 10 ){
    
    ns<- qero$ns
    nt<- qero$nt
    nc<- qero$nc
    dt<- qero$dt
    temp<- qero$temp
    
    
    tw<-round(nt/10)
    
    # print('foi até aqui')
    
    xx <- matrix(ncol = ns, nrow = nt)
    vv <- matrix(ncol = ns, nrow = nt)
    Dx <- matrix(ncol = ns, nrow = nt)
    Dv <- matrix(ncol = ns, nrow = nt)
    
    
    
    # print('foi até aqui 2')
    
    for(i in 1:ns)
    {
        
        #print(i)
        
        v<- 0
        x<- 0
        ic<- 0
        
        for(j in 1:(nt))
        {
            ic<- ic+1
            
            u<-v
            y<-x
            
            
            if(j>tw-1)
            {
                force<- -(omega0**2)*x  + a * qero$energia[j,i]
            }
            else
            {
                force<- -(omega0**2)*x
            }
            
            
            v<- u + force*dt
            x<- y + v*dt
            
  
            
            #          if(ic==nc )
            #{
            #ic <- 0
                
            xx[j,i] <- x
            vv[j,i] <- v
            
            Dx[j,i] <- v*dt
            Dv[j,i] <- force*dt
            
                #}
            
        }
    }
    
    Et<-numeric(nt)
    
    Eot<-numeric(nt)
    # print('foi ateh qui')
    
    # print(toString(length(xx)))
    #  print(toString(length(qero$energia)))
    
    for(j in 1:(nt/nc))
    {
        
        Et[j]<- (0.5)*(omega0**2) * mean(xx[j,]**2) + (0.5) * mean(vv[j,]**2)
        
        Eot[j]<- - a*mean(qero$energia[j,] * xx[j,])
        
    }
    
    #Energia do oscilador tem que convergir pra temperatura Tzero
    
    E<- (0.5)*(omega0**2)*mean(xx**2) + (0.5) * mean(vv**2)
    
    Eo<- - a*mean(qero$energia * xx)
    
    #    Resposta, começando com Dt arbitrário
    
    Dt<- 100
    
    #resposta <- matrix(nrow = nt, ncol = Dt )
    
    #  resposta <- numeric(Dt)
    
    
    # round(nt/2) é um tempo qualquer onde a principio ja temos estado estacionario
    #
    #  for(j in 1:(nt- Dt ))
    #   {
    #     for(i in 1:Dt)
    #  {
    #resposta[j,i]<-  (  mean(qero$energia[ j +i  , ]) - mean( qero$energia[ j  ,])   )/(  a*Dx[j,2]  )
    
    #       resposta[i]<-  (  mean(qero$energia[ tw +i  , ]) - mean( qero$energia[ tw  ,])   )/(  a* mean(Dx[tw,])  )
    #}
    #  }
    
    final <- list("Eoscilador"= Et, "Eint" = Eot, "ELangevin" = qero$energia, "ELangemedia" = sidemean(qero$energia)  ,"Eosctot"= E, "Esistot" = Eo,    "posicao"= xx,  "velocidade"= vv, "omega"=omega0, "temp" = temp)
    return(final)
    
}

#KurchanLiteNota<-function(qero, omega0 = 0.84 , a = 0.2 ){
KurchanLiteNota<-function(qero, omega0 = 100 , a = 21 ){
    
    wat<-KurchanLite(qero, omega0  , a )
    
    printmeansxv(wat)
    
}

KurchanLiteNotaOUT<-function(qero, omega0 = 100 , a = 21){
    
    wat<-KurchanLite(qero, omega0  , a )
    
    fim<-printmeansxvOUT(wat)
    
    return(fim)
}


Kurchanmk<-function(ns, nt, nc=1, dt, massa, gamma, ctelastica, temp , kk = 50, m = 1 , a = 10 ){
    
    # omega0 = 50 , a = 10 ){
    
    #alfa, dalfa, xa, dxa, omega0 = 10 , a = 28 ){
    
    # eh na energia seguir folha
    #<- function(ns, nt, nc, dt, massa, gamma, ctelastica, alfa, dalfa, xa, dxa, deltax, tt){
    
    #qero<-SimulacaoBimodal(ns, nt, nc, dt, massa, gamma, ctelastica, alfa, dalfa, xa, dxa)
    
    tw<-round(nt/10)
    
    qero<-SimulacaoLangevin(ns, nt, nc, dt, massa, gamma, ctelastica, temp)
    
    #wat<-FuncaoResposta(qero,Dt,)
    
    xx <- matrix(ncol = ns, nrow = nt)
    vv <- matrix(ncol = ns, nrow = nt)
    Dx <- matrix(ncol = ns, nrow = nt)
    Dv <- matrix(ncol = ns, nrow = nt)
    #zz <- matrix(ncol = ns, nrow = nt)
    
    
    for(i in 1:ns)
    {
        
        v<- 0
        x<- 0
        
        
        for(j in 1:nt)
        {
            u<-v
            y<-x
            
            #force<-  -(omega0**2)*x # + a * qero$energia[j,i] # energia de qero é o O no modelo do kurchan
            
            if(j>tw-1)
            {
                force<- - kk * x  + a * qero$energia[j,i]
            }
            else
            {
                force<- -kk*x
            }
            
            #print(j)
            #print(i)
            
            #v<- u + force*dt
            #x<- y + v*dt
            
            v<- u + force*dt
            x<- y + v*dt
            
            #ji1<-ji1+ (v+u)*eta*dt/2
            #jd1 <- jd1 + gamma*(((v+u)/2)^2)*dt
            
            xx[j,i] <- x
            vv[j,i] <- v
            
            Dx[j,i] <- v*dt
            Dv[j,i] <- force*dt
            
        }
    }
    
    Et<-numeric(nt)
    
    Eot<-numeric(nt)
    
    
    for(j in 1:nt)
    {
        
        Et[j]<- (0.5)*(omega0**2) * mean(xx[j,]**2) + (0.5) * mean(vv[j,]**2)
        
        Eot[j]<- - a*mean(qero$energia[j,] * xx[j,])
        
    }
    
    #Energia do oscilador tem que convergir pra temperatura Tzero
    
    E<- (0.5)*(omega0**2)*mean(xx**2) + (0.5) * mean(vv**2)
    
    Eo<- - a*mean(qero$energia * xx)
    
    #    Resposta, começando com Dt arbitrário
    
    Dt<- 100
    
    #resposta <- matrix(nrow = nt, ncol = Dt )
    
    #  resposta <- numeric(Dt)
    
    
    # round(nt/2) é um tempo qualquer onde a principio ja temos estado estacionario
    #
    #  for(j in 1:(nt- Dt ))
    #   {
    #     for(i in 1:Dt)
    #  {
    #resposta[j,i]<-  (  mean(qero$energia[ j +i  , ]) - mean( qero$energia[ j  ,])   )/(  a*Dx[j,2]  )
    
    #       resposta[i]<-  (  mean(qero$energia[ tw +i  , ]) - mean( qero$energia[ tw  ,])   )/(  a* mean(Dx[tw,])  )
    #}
    #  }
    
    final <- list("Eoscilador"= Et, "Eint" = Eot, "ELangevin" = qero$energia, "ELangemedia" = sidemean(qero$energia)  ,"Eosctot"= E, "Esistot" = Eo,    "posicao"= xx,  "velocidade"= vv, "omega"=omega0)
    return(final)
    
}



KurchanEx<-function(ns, nt, nc=1, dt, massa, gamma, ctelastica, temp , omega0 = 50 , a = 10 ){
    
    
    tw<-round(nt/10)
    
    
    xx <- matrix(ncol = ns, nrow = nt)
    vv <- matrix(ncol = ns, nrow = nt)
    Dx <- matrix(ncol = ns, nrow = nt)
    Dv <- matrix(ncol = ns, nrow = nt)
    #zz <- matrix(ncol = ns, nrow = nt)
    
    
    for(i in 1:ns)
    {
        print(i)
        v<- 0
        x<- 0
        
        for(j in 1:nt)
        {
            
            u<-v
            y<-x
            
            if(j>tw-1)
            {
                force<- -(omega0**2)*x  + a * rnorm(mean=0,sd = temp,1)
            }
            else
            {
                force<- -(omega0**2)*x
            }
            
            v<- u + force*dt
            x<- y + v*dt
            
            xx[j,i] <- x
            vv[j,i] <- v
            
            Dx[j,i] <- v*dt
            Dv[j,i] <- force*dt
            
        }
    }
    
    Et<-numeric(nt)
    
    Eot<-numeric(nt)
    
    
    for(j in 1:nt)
    {
        Et[j]<- (0.5)*(omega0**2) * mean(xx[j,]**2) + (0.5) * mean(vv[j,]**2)
    }
    
    E<- (0.5)*(omega0**2)*mean(xx**2) + (0.5) * mean(vv**2)
    
    final <- list("Eoscilador"= Et, "Eint" = Eot, "Eosctot"= E,     "posicao"= xx,  "velocidade"= vv, "omega"=omega0)
    return(final)
    return(final)
    
}


printmeansxv<-function(qero){
    
    
    x2<-( (  (qero$omega**2) * mean(qero$posicao**2)  ) )/2
    v2<-mean(qero$velocidade**2)/2
    
    tot<- x2 + v2
    
    totdif<- (x2-v2)
    totdifesc<- ((x2-v2)/v2)*100
    
    #entot<-
    
    diftemp <- (x2+v2) - qero$temp
    
    diftemp2 <- (x2+v2) - qero$temp**2
    
    diftempesc <- (((x2+v2) - qero$temp)/qero$temp)*100
    
    diftempesc2 <- (((x2+v2) - qero$temp**2)/qero$temp**2)*100
    #qero$omega
    temperatura<- qero$temp
    
    
    print(paste( "temperatura é", toString(temperatura), sep='') )
    cat("\n")
    print(paste( "a medida do termômetro é (total de energia)", toString(tot), sep='') )
    cat("\n")
    print(paste("O valor da energia potencial é =", toString(x2), sep='') )
    cat("\n")
    print(paste("O valor da energia cinética é =", toString(v2), sep='') )
    cat("\n")
    print(paste( "a diferença de energia é ", toString(totdif), sep='') )
    cat("\n")
    print(paste( "a diferença escalada em porcentgem entre kinetic e potential é ", toString(totdifesc), " %  ", sep='') )
    cat("\n")
    print(paste( "a diferença entre a temperatura e a medida do termometro é ", toString(diftemp), sep='') )
    cat("\n")
    print(paste( "a diferença entre a temperatura e a medida do termometro quadrado é ", toString(diftemp2), sep='') )
    cat("\n")
    print(paste( "o erro do termometro linear escalada em porcentgem é ", toString(diftempesc), " %  ", sep='') )
    cat("\n")
    print(paste( "o erro do termometro quadrado escalada em porcentgem é ", toString(diftempesc2), " %  ", sep='') )
    
    
}

printmeansxvOUT<-function(qero){
    
    
    x2<-( (  (qero$omega**2) * mean(qero$posicao**2)  ) )/2
    v2<-mean(qero$velocidade**2)/2
    
    tot<- x2 + v2
    
    totdif<- (x2-v2)
    totdifesc<- ((x2-v2)/v2)*100
    
    #entot<-
    
    diftemp <- (x2+v2) - qero$temp
    
    diftemp2 <- (x2+v2) - qero$temp**2
    
    diftempesc <- (((x2+v2) - qero$temp)/qero$temp)*100
    diftempesc2 <- (((x2+v2) - qero$temp**2)/qero$temp**2)*100
    #qero$omega
    temperatura<- qero$temp
    
    
    # print(paste( "temperatura é", toString(temperatura), sep='') )
    #cat("\n")
    #print(paste( "a medida do termômetro é (total de energia)", toString(tot), sep='') )
    #cat("\n")
    #print(paste("O valor da energia potencial é =", toString(x2), sep='') )
    #cat("\n")
    #print(paste("O valor da energia cinética é =", toString(v2), sep='') )
    #cat("\n")
    #print(paste( "a diferença de energia é ", toString(totdif), sep='') )
    #cat("\n")
    #print(paste( "a diferença escalada em porcentgem entre kinetic e potential é ", toString(totdifesc), " %  ", sep='') )
    #cat("\n")
    #print(paste( "a diferença entre a temperatura e a medida do termometro é ", toString(diftemp), sep='') )
    #cat("\n")
    #print(paste( "a diferença entre a temperatura e a medida do termometro quadrado é ", toString(diftemp2), sep='') )
    #cat("\n")
    #print(paste( "o erro do termometro linear escalada em porcentgem é ", toString(diftempesc), " %  ", sep='') )
    #cat("\n")
    #print(paste( "o erro do termometro quadrado escalada em porcentgem é ", toString(diftempesc2), " %  ", sep='') )
    
    final <- list("temp"= temperatura, "termometro" = tot, "Epot" = x2, "Ekin" = v2  ,"DifE"= totdif, "DifEesc" = totdifesc, "errolinear"= diftemp, "erroquadrado"= diftemp2,  "errolinearesc"= diftempesc, "erroquadradoesc"= diftempesc2)
    return(final)
}


FazGraficoTermometro<-function(){
    
    watb05<-KurchanLiteNotaOUT(qerol05)
    watl05<-KurchanLiteNotaOUT(qero05)
    
    watb10<-KurchanLiteNotaOUT(qerol1)
    watl10<-KurchanLiteNotaOUT(qero)
    
    watb15<-KurchanLiteNotaOUT(qerol15)
    watl15<-KurchanLiteNotaOUT(qero15)
    
    watb20<-KurchanLiteNotaOUT(qerol2)
    watl20<-KurchanLiteNotaOUT(qero2)
    
    watb25<-KurchanLiteNotaOUT(qerol25)
    watl25<-KurchanLiteNotaOUT(qero25)
    
    watb30<-KurchanLiteNotaOUT(qerol3)
    watl30<-KurchanLiteNotaOUT(qero3)
    
    watb35<-KurchanLiteNotaOUT(qerol35)
    watl35<-KurchanLiteNotaOUT(qero35)
    
    watb40<-KurchanLiteNotaOUT(qerol4)
    watl40<-KurchanLiteNotaOUT(qero4)
    
    
    errosb<-numeric(8)
    
    errosb[1]<-watb05$erroquadradoesc
    errosb[2]<-watb10$erroquadradoesc
    errosb[3]<-watb15$erroquadradoesc
    errosb[4]<-watb20$erroquadradoesc
    errosb[5]<-watb25$erroquadradoesc
    errosb[6]<-watb30$erroquadradoesc
    errosb[7]<-watb35$erroquadradoesc
    errosb[8]<-watb40$erroquadradoesc
    
    errosl<-numeric(8)
    
    errosl[1]<-watl05$erroquadradoesc
    errosl[2]<-watl10$erroquadradoesc
    errosl[3]<-watl15$erroquadradoesc
    errosl[4]<-watl20$erroquadradoesc
    errosl[5]<-watl25$erroquadradoesc
    errosl[6]<-watl30$erroquadradoesc
    errosl[7]<-watl35$erroquadradoesc
    errosl[8]<-watl40$erroquadradoesc
    
    
    templ<- seq(0.5,4,0.5)
    
    terml<- numeric(8)
    
    terml[1]<-watl05$termometro
    terml[2]<-watl10$termometro
    terml[3]<-watl15$termometro
    terml[4]<-watl20$termometro
    terml[5]<-watl25$termometro
    terml[6]<-watl30$termometro
    terml[7]<-watl35$termometro
    terml[8]<-watl40$termometro
    
    termb<- numeric(8)
    
    termb[1]<-watb05$termometro
    termb[2]<-watb10$termometro
    termb[3]<-watb15$termometro
    termb[4]<-watb20$termometro
    termb[5]<-watb25$termometro
    termb[6]<-watb30$termometro
    termb[7]<-watb35$termometro
    termb[8]<-watb40$termometro
    
    balancol<- numeric(8)
    
    balancol[1]<-watl05$DifEesc
    balancol[2]<-watl10$DifEesc
    balancol[3]<-watl15$DifEesc
    balancol[4]<-watl20$DifEesc
    balancol[5]<-watl25$DifEesc
    balancol[6]<-watl30$DifEesc
    balancol[7]<-watl35$DifEesc
    balancol[8]<-watl40$DifEesc
    
    balancob<- numeric(8)
    
    balancob[1]<-watb05$DifEesc
    balancob[2]<-watb10$DifEesc
    balancob[3]<-watb15$DifEesc
    balancob[4]<-watb20$DifEesc
    balancob[5]<-watb25$DifEesc
    balancob[6]<-watb30$DifEesc
    balancob[7]<-watb35$DifEesc
    balancob[8]<-watb40$DifEesc
    
    
    
    grafs<-'termometro.jpg'
    jpeg(grafs)
    
    plot(terml~templ, main = "Termometro", xlab='Temperatura', ylab= 'Raiz quadrada da Energia do oscilador', col='blue')
    points(termb~templ, col='red')
    dev.off()
    
    
    grafs<-'errotermometro.jpg'
    jpeg(grafs)
    
    plot(errosb~templ, main = "Erro do termômetro quadrado", xlab='Temperatura', ylab= 'Erro', col='red')
    points(errosl~templ, col='blue')
    dev.off()
    
    grafs<-'balancotermometro.jpg'
    jpeg(grafs)
    
    plot(balancob~templ, main = "Balanço Cinética x Potencial", xlab='Temperatura', ylab= 'Balanço percentual = 100*(Ekin + Epot)/Epot', col='red')
    points(balancol~templ, col='blue')
    dev.off()
    
    
    #grafs<-paste(letra , "TRespostaDcorr", toString(vez) ,".jpg", sep="")
    #grafs<-toString(grafs)
    #jpeg(grafs)
    
    #plot( respostaam ~ tempo , col = 'red' )
    #lines(   - dcorranalitico ~ tempo, col = 'blue' )
    #lines(   - dcorrlange ~ tempo, col = 'green' )
    
    #dev.off()
    
    
    
}

searchparameter<-function(){
    aseq<- seq(0.1,10,0.1)
    ooseq<- seq(0.1,10,0.1)
}

Kurchan2<-function(ns, nt, nc=1, dt, massa, gamma, ctelastica, temp , omega0 = 50 , a = 10 ){
    
    #alfa, dalfa, xa, dxa, omega0 = 10 , a = 28 ){
    
    # eh na energia seguir folha
    #<- function(ns, nt, nc, dt, massa, gamma, ctelastica, alfa, dalfa, xa, dxa, deltax, tt){
    
    #qero<-SimulacaoBimodal(ns, nt, nc, dt, massa, gamma, ctelastica, alfa, dalfa, xa, dxa)
    
    tw<-round(nt/10)
    
    qero<-SimulacaoLangevin(ns, nt, nc, dt, massa, gamma, ctelastica, temp)
    
    #wat<-FuncaoResposta(qero,Dt,)
    
    xx <- matrix(ncol = ns, nrow = nt)
    vv <- matrix(ncol = ns, nrow = nt)
    Dx <- matrix(ncol = ns, nrow = nt)
    Dv <- matrix(ncol = ns, nrow = nt)
    #zz <- matrix(ncol = ns, nrow = nt)
    
    
    for(i in 1:ns)
    {
        
        v<- 0
        x<- 3
        
        
        for(j in 1:nt)
        {
            u<-v
            y<-x
            
            #force<-  -(omega0**2)*x # + a * qero$energia[j,i] # energia de qero é o O no modelo do kurchan
            
            if(j>tw-1)
            {
                force<- -(omega0**2)*x  + a * qero$posicao[j,i]
            }
            else
            {
                force<- -(omega0**2)*x
            }
            
            #print(j)
            #print(i)
            
            v<- u + force*dt
            x<- y + v*dt
            
            #ji1<-ji1+ (v+u)*eta*dt/2
            #jd1 <- jd1 + gamma*(((v+u)/2)^2)*dt
            
            xx[j,i] <- x
            vv[j,i] <- v
            
            Dx[j,i] <- v*dt
            Dv[j,i] <- force*dt
            
        }
    }
    
    Et<-numeric(nt)
    
    Eot<-numeric(nt)
    
    
    for(j in 1:nt)
    {
        
        Et[j]<- (0.5)*(omega0**2) * mean(xx[j,]**2) + (0.5) * mean(vv[j,]**2)
        
        Eot[j]<- - a*mean(qero$posicao[j,] * xx[j,])
        
    }
    
    #Energia do oscilador tem que convergir pra temperatura Tzero
    
    E<- (0.5)*(omega0**2)*mean(xx**2) + (0.5)*mean(vv**2)
    
    Eo<- - a*mean(qero$posicao * xx)
    
    #    Resposta, começando com Dt arbitrário
    
    Dt<- 100
    
    #resposta <- matrix(nrow = nt, ncol = Dt )
    
    resposta <- numeric(Dt)
    
    
    # round(nt/2) é um tempo qualquer onde a principio ja temos estado estacionario
    #
    #  for(j in 1:(nt- Dt ))
    #  {
    for(i in 1:Dt)
    {
        #resposta[j,i]<-  (  mean(qero$energia[ j +i  , ]) - mean( qero$energia[ j  ,])   )/(  a*Dx[j,2]  )
        
        resposta[i]<-  (  mean(qero$posicao[ tw +i  , ]) - mean( qero$posicao[ tw  ,])   )/(  a* mean(Dx[tw,])  )
    }
    #  }
    
    final <- list("Eoscilador"= Et, "Eint" = Eot, "ELangevin" = qero$energia, "ELangemedia" = sidemean(qero$energia)  ,"Eosctot"= E, "Esistot" = Eo,    "resposta"= resposta, "posicao"= xx,  "velocidade"= vv, "omega"= omega0)
    return(final)
    
}



LorenzModel<-function(ns, nt, nc, sigma = 10 , r = 28 , b = (8/3) , deltaz = 8.67, dt=0.001, tkick ){
    
    xx <- matrix(ncol = ns, nrow = nt)
    yy <- matrix(ncol = ns, nrow = nt)
    zz <- matrix(ncol = ns, nrow = nt)
    
    x <- 10
    y <- 10
    z <- 10
    
    for(i in 1:ns){
    
        print(i)
    
    
        for(j in 1:nt){
            
            if (j == tkick)
            {
                 z <- z + deltaz
            }
            
            for(k in 1:nc){
                
                
                dx <-  sigma * dt * (y-x)
                
                dy <- (-x*z +r*x - y) * dt
                
                dz <- (x*y-b*z) * dt
                
                
                x <- x + dx
                
                y <- y + dy
                
                z <- z + dz
            
            }
            
            xx[j,i]<- x
            yy[j,i]<- y
            zz[j,i]<- z
            
            
        }
    
    }
    
    final <- list("x"= xx , "y" = yy , "z" = zz, "tkick" = tkick, "ns" = ns, "nt" = nt, "nc" = nc, "kick" = deltaz )
    return(final)
    
}

TesteLorenzModel<-function(){
    
    # qero<- LorenzModel(ns = 10000, nt = 1000, deltaz = 8.67 , nc = 100, tkick = 5000, dt=)
    
    # 2
    #qero<- LorenzModel(ns = 1000, nt = 1000, deltaz = 8.67 , nc = 100, tkick = 500, dt=0.01)
    
    # 4
    # qero<- LorenzModel(ns = 100, nt = 100000, deltaz = 8.67 , nc = 1, tkick = 5000)
    
    # 5
    
    #qero<- LorenzModel(ns = 100, nt = 1000, deltaz = 8.67 , nc = 1, tkick = 500)
    
    #6
    #qero<- LorenzModel(ns = 100, nt = 1000, deltaz = 8.67 , nc = 100, tkick = 500)
    
    #7
    #qero<- LorenzModel(ns = 100, nt = 1000, deltaz = 8.67 , nc = 10, tkick = 500)
    
    
    #8
    #qero<- LorenzModel(ns = 100, nt = 3000, deltaz = 8.67 , nc = 10, tkick = 2000)
    
    #9
    #qero<- LorenzModel(ns = 200, nt = 3000, deltaz = 8.67 , nc = 10, tkick = 1000)
    
    #9 = 10
    #qero<- LorenzModel(ns = 200, nt = 3000, deltaz = 8.67 , nc = 10, tkick = 1000)
    
    #11 = 10 diferente soh por ns
    #qero<- LorenzModel(ns = 500, nt = 3000, deltaz = 8.67 , nc = 10, tkick = 1000)
    
    #12
    #qero<- LorenzModel(ns = 1000, nt = 2000, deltaz = 8.67 , nc = 10, tkick = 1000)
    
    #13
    
    #qero<- LorenzModel(ns = 10000, nt = 2000, deltaz = 8.67 , nc = 10, tkick = 1000)
    
    #14
    
    qero<- LorenzModel(ns = 100000, nt = 400, deltaz = 0.01*8.67 , nc = 10, tkick = 200)
    
    hz<-hist(qero$z, breaks = 100)
    hhz<-hz$density
    hhb<-hz$mids
    
    csv<-'.csv'
    
    vulpiresposta<- '14lorenzrespostav'
    
    rresposta<- '14lorenzresposta'
    
    corrpre<- '14lorenzcorpre'
    
    corrpos<- '14lorenzcorpos'
    
    
    
    titulo0<-paste( corrpre , csv , sep="")
    
    titulo1<-paste( corrpos ,  csv , sep="")
    
    tituloR<-paste( rresposta , csv , sep="")
    
    ns<- 100000
    
    kick<- 8.67
    tkick<- 200
    
    Dt<- 200
    
    fim<- numeric(Dt)
    
    mediapre<-      mean( qero$z[1:tkick-1,]  )
    
    corpre<-      (  funcaodecorrelacaoSample(qero$z, (tkick-Dt -1) , Dt, ns)   )
    
    corpos<- ( funcaodecorrelacaoSample(qero$z, tkick , Dt ,ns))
    
    
    for (i in 1:Dt)
    {
        
        mediapos<-  mean(qero$z[tkick-1+i, ]  )
        
        fim[i]<- (mediapos-mediapre)
        
    }
    
    fim<-fim/kick
    
    SalvaCorrelacao(corpre, titulo0)
    
    SalvaCorrelacao(corpos, titulo1)
    
    SalvaCorrelacao(fim, tituloR)
    
    
    
    #SalvaCorrelacao(qero$x, 'LorenzX.csv')
    
    #SalvaCorrelacao(qero$y, 'LorenzY.csv')
    
    #SalvaCorrelacao(qero$z, 'LorenzZ.csv')
    
    
}

TesteLorenzModel2<-function(){
    
    # 15 qero<- LorenzModel(ns = 100000, nt = 500, deltaz = 8.67 , nc = 10, tkick = 250)
    
    qero<- LorenzModel(ns = 500000, nt = 500, deltaz = 8.67 , nc = 10, tkick = 250)
    
    hz<-hist(qero$z, breaks = 100)
    hhz<-hz$density
    hhb<-hz$mids
    
    csv<-'.csv'
    
    vulpiresposta<- '16lorenzrespostav'
    
    rresposta<- '16lorenzresposta'
    
    corrpre<- '16lorenzcorpre'
    
    corrpos<- '16lorenzcorpos'
    
    
    
    titulo0<-paste( corrpre , csv , sep="")
    
    titulo1<-paste( corrpos ,  csv , sep="")
    
    tituloR<-paste( rresposta , csv , sep="")
    
    titulovR<-paste( vulpiresposta , csv , sep="")
    
    ns<- 500000
    
    kick<- 8.67
    tkick<- 250
    
    Dt<- 200
    
    breakks<- 60
    
    fim<- numeric(Dt)
    
    mediapre<-      mean( qero$z[1:tkick,]  )
    
    print(1)
    corpre<-      (  funcaodecorrelacaoSample(qero$z, (tkick-Dt -1) , Dt, ns)   )
    
    print(2)
    corpos<- ( funcaodecorrelacaoSample(qero$z, tkick , Dt ,ns)       )
    
    print(3)
    
    respostavulpi<- RespostaLinearSampleVulpiLorenz1(qero, Dt, breakks)
    
    print(4)
    
    for (i in 1:Dt)
    {
        
        mediapos<-  mean(qero$z[tkick+i, ]  )
        
        fim[i]<- mediapos
        
    }
    
    fim<-fim/kick
    
    SalvaCorrelacao(corpre, titulo0)
    
    SalvaCorrelacao(corpos, titulo1)
    
    SalvaCorrelacao(fim, tituloR)
    
    SalvaCorrelacao(respostavulpi, titulovR)
    
    #SalvaCorrelacao(qero$x, 'LorenzX.csv')
    
    #SalvaCorrelacao(qero$y, 'LorenzY.csv')
    
    #SalvaCorrelacao(qero$z, 'LorenzZ.csv')
    
    
}

VulpianiModel<-function(ns, nt, nc, sigma = (1.07238) , D = 1 , deltax = (1.07238)*(7.6)*(0.001), dt = 0.001, tkick ){
    
    xx <- matrix(ncol = ns, nrow = nt)
    
    
    
    for(i in 1:ns){
        
        x<- 0
        
        print(i)
        
        for(j in 1:nt){
            
            if (j == tkick)
            {
                x <- x + deltax
            }
            
            for(k in 1:nc){
                
                
                if(abs(x) <= 1 ){
                    
                    #forceest <- sqrt(2*gamma*temp*dt)*rnorm(1)
                    #forcedet <-  (-gm*u -c1*y)*dt
                
                    dx <- ( -x + sqrt(2*D)*rnorm(1) )*dt
                    x <- x + dx
                    #print(dx)
                 }
                
                 if(abs(x) > 1 ){
                     
                     if(x > 0)
                     {
                     
                         dx <-  ( 1 + sqrt(2*D)*rnorm(1) )*dt
                         x <- x + dx
                         #print(dx)
                         #print('modulo maior que um maior q zero')
                     }
                     
                     if(x <= 0)
                     {
                         dx <- ( -1 + sqrt(2*D)*rnorm(1) )*dt
                         x <- x + dx
                         
                     }
                     
                 }
                
            }
            
            xx[j,i]<- x
            
        }
        
    }
    
    final <- list("z"= xx, "kick" = deltax, "tkick" = tkick, "ns"=ns )
    return(final)
    
}




TesteVulpianiModel<-function(){
    
    
    qero<-VulpianiModel(ns = 100, nt = 1000, nc = 10, tkick = 250)
    
    numteste<- 2
    
    stringera<- toString(numteste)
    
    csv<-'.csv'
    
    rresposta<- 'vulpianiresposta'
    
    rrespostav1<- 'vulpianirespostavum'
    
    rrespostav2<- 'vulpianirespostavdois'
    
    corrpre<- 'vulpianicorpre'
    
    corrpos<- 'vulpianicorpos'
    
    titulo0<-paste(stringera ,corrpre , csv , sep="")
    
    titulo1<-paste(stringera ,corrpos ,  csv , sep="")
    
    tituloR<-paste(stringera ,rresposta , csv , sep="")
    
    tituloRv1<-paste(stringera ,rrespostav1 , csv , sep="")
    
    tituloRv2<-paste(stringera ,rrespostav2 , csv , sep="")
    
    ns<- qero$ns
    
    kick<- qero$kick
    tkick<- qero$tkick
    
    Dt<- 200
    
    breakks<- 80
    
    fim<- numeric(Dt)
    
    mediapre<-      mean( qero$z[1:tkick,]  )
    
    print('aqui foi')
    
    corpre<-      (  funcaodecorrelacaoSample(qero$z, (tkick-Dt -1) , Dt, ns)   )
    
    print('aqui foi 2')
    
    corpos<- ( funcaodecorrelacaoSample(qero$z, tkick , Dt ,ns))
    
    print('aqui foi 3')
    
    resposta1<- RespostaLinearSampleVulpiLorenz1(qero, Dt, breakks)
    
    print('aqui foi 4')
    
    resposta2<- RespostaLinearSampleVulpiLorenz2(qero, Dt, breakks)
    
    print('aqui foi 5')
    
    for (i in 1:Dt)
    {
        
        mediapos<-  mean(qero$z[tkick+i, ]  )
        
        fim[i]<- mediapos - mediapre
        
    }
    
    fim<-fim/kick
    

    
    SalvaCorrelacao(corpre, titulo0)
    
    SalvaCorrelacao(corpos, titulo1)
    
    SalvaCorrelacao(fim, tituloR)
    
    SalvaCorrelacao(resposta1, tituloRv1)
    
    SalvaCorrelacao(resposta2, tituloRv2)
    
    
    
}

rascunhoh<- function(h){

hy<-hist(y,breaks=50)
sy<- - hy$density * log(hy$density) # (hx$breaks[2]-hx$breaks[1])
sy[which(!is.finite(sy))]<- 0
ssy<-sum(sy)*(1/sum(hy$density))

hx<-hist(x,breaks=50)
sx<- - hx$density * log(hx$density)* (hx$breaks[2]-hx$breaks[1])
sx[which(!is.finite(sx))]<- 0
ssx<-sum(sx)

}


Histograma3D<- function(x,y, breaks){
    
    # x e y tem o mesmo comprimento
    
    len<- length(x) # = length(y)
    
    rangeX <- max(x) - min(x)
    
    rangeY <- max(y) - min(y)
    
    intervalox <- rangeX/breaks
    
    intervaloy <- rangeY/breaks
    
    midsxx<- matrix(ncol = breaks, nrow = breaks)
    
    midsyy<- matrix(ncol = breaks, nrow = breaks)
    
    # breaksxx<- matrix(ncol = breaks+1, nrow = breaks)
    
    # breaksyy<- matrix(ncol = breaks, nrow = breaks)
    
    breaksx<- seq( min(x) ,  max(x) , intervalox ) # length BREAKS +1
    
    #print(breaksx)
    
    breaksy<- seq( min(y) ,  max(y) , intervaloy ) # length BREAKS +1
    
    #print(breaksy)
    # for(i in 1:len)
    #{
        
        #    xx[,i]<-rep(breaksx[i], len )
        # yy[i,]<-rep(breaksy[i], len )
        
        #}
    
    pesos<- matrix(ncol = breaks, nrow = breaks)
    pesos<- is.na(pesos) - 0L
    pesos<- pesos - 1

    
    for( i in 1: len)
    {
        #print(i)
        
        for(j in 1:(breaks))
        {
            
            if(breaksx[j] <= x[i] &&  x[i] < breaksx[j+1] )
            {
                indicex<- j
                #print('pegou indicex')
                break
            }
            if(x[i]==max(x))
            {
                indicex<- breaks
                #print('pegou indicex')
                break
            }
        }
        
        for(j in 1:(breaks))
        {
            if(breaksy[j] <= y[i] &&  y[i]  < breaksy[j+1] )
            {
                indicey<- j
                #print('pegou indicey')
                break
            }
            if(y[i]==max(y))
            {
                indicey<- breaks
                #print('pegou indicey')
                break
            }
            
        }
        
        pesos[indicex,indicey] <- pesos[indicex,indicey] + 1
        
    }
    
    pesos <- pesos /len
    
    pesos <- pesos / (( breaksx[2] - breaksx[1] ) * ( breaksy[2] - breaksy[1] ))
    
    #return(pesos)
    
    final<-list("density"=pesos,"breaksx"=breaksx,"breaksy"=breaksy )
    
    return(final)
    
}

cumulantes<- function(a){
    cum<-numeric(6)
    
    cum[1]<-a[2]
    cum[2]<-a[3] - a[2]^2
    cum[3]<-a[4] - 3*a[3]*a[2] + 2*a[2]^3
    cum[4]<-a[5] - 4*a[4]*a[2] - 3*a[3]^2 + 12*a[3]*a[2]^2  - 6*a[2]^4
    cum[5]<-a[6] - 5*a[5]*a[2] - 10*a[4]*a[3] + 20*a[4]*a[2]^2 + 30*a[3]^2 * a[2] - 60 * a[3] * a[2]^3 + 24*a[2]^5
    cum[6]<-a[7] - 6*a[6]*a[2] - 15*a[5]*a[3] + 30*a[5]*a[2]^2 - 10*a[4]^2 + 120 * a[4] * a[3] * a[2] - 120 * a[4] * a[2]^3 + 30*a[3]^3 - 270 * a[3]^2 * a[2]^2 + 360 * a[3] * a[2]^4 - 120* a[2]^6

    return(cum)
}


pegaordemK<-function(a, ordem){
    
    fim<-numeric(length(a))
    
    for(i in 1:length(a))
    {
        fim[i]<-a[[i]][ordem]
    }
    return(fim)
}


pegaordemM<-function(a, ordem){
    
    fim<-numeric(length(a))
    
    for(i in 1:length(a))
    {
        fim[i]<-a[[i]][ordem+1]
    }
    return(fim)
}


derivadazin<-function(a,b){
    
    
    nn<-length(a)-1
    
    fim<-numeric(nn)
    
    for(i in 1: nn )
    {
        fim[i]<-(a[i+1]-a[i])/(b[i+1]-b[i])
        #  derivhard[i]<-(fim$S[i+1]-fim$S[i])/(fim$E[i+1]-fim$E[i])
    }
    
    return(fim)
    
}


SimulacaoLangevin <- function(ns, nt, nc, dt, massa, gamma, ctelastica, temp){
    
    
    ntc<-nt*nc
    
    xx<-matrix( nrow = nt, ncol = ns)
    vv<-matrix( nrow = nt, ncol = ns)
    ee<-matrix( nrow = nt, ncol = ns)
    ji<-matrix( nrow = nt, ncol = ns)
    jd<-matrix( nrow = nt, ncol = ns)
    
    
    c1<- ctelastica/massa
    
    gm<- gamma/massa
    
    ie<- 0
    
    it<- 0
    
    
    
    for(i in 1:ns)
    {
        
        print(i)
        
        x<-0
        v<-0
        u<-0
        y<-0
        
        ji1<-0
        jd1<-0
        
        
        it<-0
        ic<-0
        ir<-0
        ig<-0
        
        forceest0<-0
        forcedet0<-0
        
        while(it < ntc)
        {
            
            it <- it + 1
            
            forceest <- sqrt(2*gamma*temp*dt)*rnorm(1)
            
            #forceeststrato<- (forceest + forceest0)/2
            #forceest0<-forceest
            
            forcedet <-  (-gm*u -c1*y)*dt
            
            #forcedetstrato <- (forcedet + forcedet0)/2
            #forcedet0 <- forcedet
            
            v <- u + forcedet + forceest
            x <- y + v*dt
            
            ji1 <- ji1 + (v+u)*forceest/2
            jd1 <- jd1 + gamma*(((v+u)/2)^2)*dt
            
            # caso bimodal
            # ji1<-ji1+ (v+u)*eta*dt/2
            # jd1<-jd1+ gamma*(v^2+u^2)*dt/2
            
            ig <- ig + 1
            
            u <- v
            y <- x
            
            if(ig==nc)
            {
                
                xx[it/nc,i]<- x
                vv[it/nc,i]<- v
                ee[it/nc,i]<- (ji1-jd1)
                ji[it/nc,i]<- ji1
                jd[it/nc,i]<- jd1
                
                ig<- 0
                
            }
        }
        
    }
    
    
    final<-list("posicao"=xx,"velocidade"=vv, "ji"= ji, "jd"=jd,"energia"=ee,  "massa"=massa,"gamma"=gamma,"k"=ctelastica, "ns" = ns ,"nt"=nt,"nc"=nc,"dt"= dt, "temp"=temp)
    
    return(final)
}

RuidoColorido<- function(alfa,len, media = 0, sigma = 1, dts=0.0001)
{
    
    xi<-0
    xii<-numeric(len)
    ruido<-rnorm(mean = media, sd = sigma,n = len)*dts
    
    for( i in 1:999)
    {
        for(j in 1:999)
        {
            xi<-xi + exp(-alfa*i)*0.5*(ruido[j]+ruido[j+1])
        }
        xii[i]<-xi
    }
    return(xii)
}


SimulacaoLangevinColorida <- function(ns, nt, nc, dt, massa, gamma, ctelastica, temp, alfa){
    
    
    ntc<-nt*nc
    
    xx<-matrix( nrow = nt, ncol = ns)
    vv<-matrix( nrow = nt, ncol = ns)
    ee<-matrix( nrow = nt, ncol = ns)
    ji<-matrix( nrow = nt, ncol = ns)
    jd<-matrix( nrow = nt, ncol = ns)
    
    
    c1<- ctelastica/massa
    
    gm<- gamma/massa
    
    ie<- 0
    
    it<- 0
    
    
    
    for(i in 1:ns)
    {
        #force2<- force2 * exp(- alfa * dt  )  - gm * u * exp(- alfa * dt  ) * dt
        
        
        forceestpre <- numeric(ntc)
        forceestpos <- numeric(ntc)
        #forceest <- sqrt(2*gamma*temp*dt)*rnorm(1)

        print(i)
        
        x<-0
        v<-0
        u<-0
        y<-0
        
        ji1<-0
        jd1<-0
        
        
        it<-0
        ic<-0
        ir<-0
        ig<-0
        
        forceest0<-0
        forcedet0<-0
        
        while(it < ntc)
        {
            
            it <- it + 1
            
            
            
            forceestpre[it] <- sqrt(2*gamma*temp*dt)*rnorm(1)
            
            forceestpos[it] <- forceestpos[it]#xi = xi + DEXP(-REAL(ir)*alfa) *dws
            
            forceestpos<- force2 * exp(- alfa * dt  )  - gm * u * exp(- alfa * dt  ) * dt
            
            forceest <- forcest + exp()
            
            
            #forceeststrato<- (forceest + forceest0)/2
            #forceest0<-forceest
            
            forcedet <-  (-gm*u -c1*y)*dt
            
            #forcedetstrato <- (forcedet + forcedet0)/2
            #forcedet0 <- forcedet
            
            v <- u + forcedet + forceest
            x <- y + v*dt
            
            ji1 <- ji1 + (v+u)*forceest/2
            jd1 <- jd1 + gamma*(((v+u)/2)^2)*dt
            
            # caso bimodal
            # ji1<-ji1+ (v+u)*eta*dt/2
            # jd1<-jd1+ gamma*(v^2+u^2)*dt/2
            
            ig <- ig + 1
            
            u <- v
            y <- x
            
            if(ig==nc)
            {
                
                xx[it/nc,i]<- x
                vv[it/nc,i]<- v
                ee[it/nc,i]<- (ji1-jd1)
                ji[it/nc,i]<- ji1
                jd[it/nc,i]<- jd1
                
                ig<- 0
                
            }
        }
        
    }
    
    
    final<-list("posicao"=xx,"velocidade"=vv, "ji"= ji, "jd"=jd,"energia"=ee,  "massa"=massa,"gamma"=gamma,"k"=ctelastica, "ns" = ns ,"nt"=nt,"nc"=nc,"dt"= dt, "temp"=temp)
    
    return(final)
}

SimulacaoColoridaInterno <- function(ns, nt, nc, dt, massa, gamma, ctelastica, temp){
    
    
    ntc<-nt*nc
    
    xx<-matrix( nrow = nt, ncol = ns)
    vv<-matrix( nrow = nt, ncol = ns)
    ee<-matrix( nrow = nt, ncol = ns)
    ji<-matrix( nrow = nt, ncol = ns)
    jd<-matrix( nrow = nt, ncol = ns)
    
    
    c1<- ctelastica/massa
    
    gm<- gamma/massa
    
    ie<- 0
    
    it<- 0
    
    
    
    for(i in 1:ns)
    {
        
        print(i)
        
        x<-0
        v<-0
        u<-0
        y<-0
        
        ji1<-0
        jd1<-0
        
        
        it<-0
        ic<-0
        ir<-0
        ig<-0
        
        forceest0<-0
        forcedet0<-0
        
        while(it < ntc)
        {
            
            it <- it + 1
            
            forceest <- sqrt(2*gamma*temp*dt)*rnorm(1)
            
            #forceeststrato<- (forceest + forceest0)/2
            #forceest0<-forceest
            
            forcedet <-  (-gm*u -c1*y)*dt
            
            #forcedetstrato <- (forcedet + forcedet0)/2
            #forcedet0 <- forcedet
            
            v <- u + forcedet + forceest
            x <- y + v*dt
            
            ji1 <- ji1 + (v+u)*forceest/2
            jd1 <- jd1 + gamma*(((v+u)/2)^2)*dt
            
            # caso bimodal
            # ji1<-ji1+ (v+u)*eta*dt/2
            # jd1<-jd1+ gamma*(v^2+u^2)*dt/2
            
            ig <- ig + 1
            
            u <- v
            y <- x
            
            if(ig==nc)
            {
                
                xx[it/nc,i]<- x
                vv[it/nc,i]<- v
                ee[it/nc,i]<- (ji1-jd1)
                ji[it/nc,i]<- ji1
                jd[it/nc,i]<- jd1
                
                ig<- 0
                
            }
        }
        
    }
    
    
    final<-list("posicao"=xx,"velocidade"=vv, "ji"= ji, "jd"=jd,"energia"=ee,  "massa"=massa,"gamma"=gamma,"k"=ctelastica, "ns" = ns ,"nt"=nt,"nc"=nc,"dt"= dt, "temp"=temp)
    
    return(final)
}


SalvaLangevin<- function(qero, nome){

    posicao<- paste(nome,'posicao.csv',sep="")
    SalvaCorrelacao(qero$posicao, posicao)
    
    velocidade<- paste(nome,'velocidade.csv',sep="")
    SalvaCorrelacao(qero$velocidade, velocidade)
    
    ji<- paste(nome,'ji.csv',sep="")
    SalvaCorrelacao(qero$ji, ji)
    
    jd<- paste(nome,'jd.csv',sep="")
    SalvaCorrelacao(qero$jd, jd)
    
    energia<- paste(nome,'energia.csv',sep="")
    SalvaCorrelacao(qero$energia, energia)
    
    massa<- paste(nome,'massa.csv',sep="")
    SalvaCorrelacao(qero$massa, massa)
    
    gamma<- paste(nome,'gamma.csv',sep="")
    SalvaCorrelacao(qero$gamma, gamma)
    
    k<- paste(nome,'k.csv',sep="")
    SalvaCorrelacao(qero$k, k)
    
    ns<- paste(nome,'ns.csv',sep="")
    SalvaCorrelacao(qero$ns, ns)
    
    nt<- paste(nome,'nt.csv',sep="")
    SalvaCorrelacao(qero$nt, nt)
    
    nc<- paste(nome,'nc.csv',sep="")
    SalvaCorrelacao(qero$nc, nc)
    
    dt<- paste(nome,'dt.csv',sep="")
    SalvaCorrelacao(qero$dt, dt)
    
    temp<- paste(nome,'temp.csv',sep="")
    SalvaCorrelacao(qero$temp, temp)
    
}

CarregaLange<-function( nome)
{
    
    nt<-read.csv(paste(nome ,'nt.csv',sep=""))
    nt<-nt[1,1]
    
    ns<-read.csv(paste(nome ,'ns.csv',sep=""))
    ns<-ns[1,1]
    
    nc<-read.csv(paste(nome ,'nc.csv',sep=""))
    nc<-nc[1,1]
    
    dt<-read.csv(paste(nome ,'dt.csv',sep=""))
    dt<-dt[1,1]
    
    temp<-read.csv(paste(nome ,'temp.csv',sep=""))
    temp<-temp[1,1]
    
    massa<-read.csv(paste(nome ,'massa.csv',sep=""))
    massa<-massa[1,1]
    
    gamma<-read.csv(paste(nome ,'gamma.csv',sep=""))
    gamma<-gamma[1,1]
    
    ctelastica<-read.csv(paste(nome ,'k.csv',sep=""))
    ctelastica<-ctelastica[1,1]
    
    x<-read.csv(paste(nome ,'posicao.csv',sep=""))
    xx<-matrix( nrow = nt, ncol = ns)
    for(i in 1:ns){for(j in 1:nt){xx[i,j]<-x[i,j] }}
    
    v<-read.csv(paste(nome ,'velocidade.csv',sep=""))
    vv<-matrix( nrow = nt, ncol = ns)
    for(i in 1:ns){for(j in 1:nt){vv[i,j]<-v[i,j] }}
    
    e<-read.csv(paste(nome ,'energia.csv',sep=""))
    ee<-matrix( nrow = nt, ncol = ns)
    for(i in 1:ns){for(j in 1:nt){ee[i,j]<-e[i,j] }}
    
    
    jji<-read.csv(paste(nome ,'ji.csv',sep=""))
    ji<-matrix( nrow = nt, ncol = ns)
    for(i in 1:ns){for(j in 1:nt){ji[i,j]<-jji[i,j] }}
    
    jjd<-read.csv(paste(nome ,'jd.csv',sep=""))
    jd<-matrix( nrow = nt, ncol = ns)
    for(i in 1:ns){for(j in 1:nt){jd[i,j]<-jjd[i,j] }}
    
    
    final<-list("posicao"=xx,"velocidade"=vv, "ji"= ji, "jd"=jd,"energia"=ee,  "massa"=massa,"gamma"=gamma,"k"=ctelastica, "ns" = ns ,"nt"=nt,"nc"=nc,"dt"= dt, "temp"=temp)
    
    return(final)
    
}

SimulacaoLangevinKickF <- function(ns, nt, nc, dt, massa, gamma, ctelastica, temp, deltaf, tt){
    
    
    ntc<-nt*nc
    
    xx<-matrix( nrow = nt, ncol = ns)
    vv<-matrix( nrow = nt, ncol = ns)
    ee<-matrix( nrow = nt, ncol = ns)
    ji<-matrix( nrow = nt, ncol = ns)
    jd<-matrix( nrow = nt, ncol = ns)
    
    c1<- ctelastica/massa
    
    gm<- gamma/massa
    
    ie<- 0
    
    it<- 0
    
    
    
    for(i in 1:ns)
    {
        
        print(i)
        
        x<-0
        v<-0
        u<-0
        y<-0
        
        ji1<-0
        jd1<-0
        
        
        it<-0
        ic<-0
        ir<-0
        ig<-0
        
        forceest0<-0
        forcedet0<-0
        
        indicadorkick<- 0
        
        while(it < ntc)
        {
            
            it <- it + 1
            
            forceest <- sqrt(2*gamma*temp*dt)*rnorm(1)
            
    
            #forceeststrato<- (forceest + forceest0)/2
            #forceest0<-forceest
            
            forcedet <-  (-gm*u -c1*y)*dt
            
            #forcedetstrato <- (forcedet + forcedet0)/2
            #forcedet0 <- forcedet
            
            
            
            
            force <- forcedet + forceest
            
            
            if( it*dt >= tt & indicadorkick==0)
            {
                force <- force + deltaf
                print(i)
                indicadorkick<- 1
            }
            
            v <- u + force
            x <- y + v*dt
            
            ji1 <- ji1 + (v+u)*forceest/2
            jd1 <- jd1 + gamma*(((v+u)/2)^2)*dt
            
            # caso bimodal
            # ji1<-ji1+ (v+u)*eta*dt/2
            # jd1<-jd1+ gamma*(v^2+u^2)*dt/2
            
            ig <- ig + 1
            
            u <- v
            y <- x
            
            if(ig==nc)
            {
                
                xx[it/nc,i]<- x
                vv[it/nc,i]<- v
                ee[it/nc,i]<- (ji1-jd1)
                ji[it/nc,i]<- ji1
                jd[it/nc,i]<- jd1
                
                ig<- 0
                
            }
        }
        
    }
    
    
    final<-list("posicao"=xx,"velocidade"=vv, "ji"= ji, "jd"=jd,"energia"=ee,  "massa"=massa,"gamma"=gamma,"k"=ctelastica, "ns" = ns ,"nt"=nt,"nc"=nc,"dt"= dt, "temp"=temp, "kick" = deltaf, "tkick" = tt)
    
    return(final)
}



analiseentropicaLangevin <- function(listasim ){
    
    #xai<-listasim$xa
    
    mediaE<-numeric(listasim$ns)
    
    
    for(i in 1:listasim$ns){
        mediaE[i]<-mean(listasim$energia[,i])
    }
    
    
    #ENERGIA MÉDIA ANALÍTICA
    
    
    
    
    # histogramas
    stoth<-numeric(listasim$ns)
    for(i in 1:listasim$ns){
        
        hhx<-hist(listasim$posicao[,i], breaks=100, plot= FALSE)
        hhv<-hist(listasim$velocidade[,i], breaks=100, plot= FALSE)
        
        ssx<- - hhx$density*log(hhx$density)*(hhx$breaks[2]-hhx$breaks[1])
        ssx[which(is.na(ssx))]<- 0
        sx<-sum(ssx)*(1/sum(hhx$density))
        ssv<- - hhv$density*log(hhv$density)
        ssv[which(is.na(ssv))]<- 0
        sv<-sum(ssv)*(1/sum(hhv$density))
        
        
        stoth[i]<- sx + sv
    }
    
    
    
    fim<-list("S"=stoth, "E"=mediaE,  "ns"=listasim$ns,"nt"=listasim$nt,"nc"=listasim$nc,"dt"= listasim$dt)
    
    
    return(fim)
    
    
}



evolucaoentropica<-function(listasim){
    
    
    if(listasim$dalfa==0)
    {
        alfai<-rep(listasim$alfa, listasim$ns)
    } else {
        alfai<-seq(listasim$alfa+listasim$dalfa,listasim$alfa+listasim$ns*listasim$dalfa,listasim$dalfa)
    }
    
    if(listasim$dxa==0)
    {
        xai<-rep(listasim$xa, listasim$ns)
    } else {
        xai<-seq(listasim$xa+listasim$dxa,listasim$xa+listasim$ns*listasim$dxa,listasim$dxa)
    }
    
    mediaE<-numeric(listasim$ns)
    sdE<-numeric(listasim$ns)
    
    #ENERGIA MÉDIA ANALÍTICA
    
    anaE<- (((xai)^2 )*(2*listasim$massa*alfai +listasim$gamma))/(2*listasim$gamma*(listasim$k +listasim$massa *( alfai^2) + alfai * listasim$gamma))
    
    
    
    for(i in 1:listasim$ns){
        
        mediaE[i]<-mean(listasim$energia[,i])
        
        sdE[i]<-sd(listasim$energia[,i])
        
    }
    
    erroE<- abs(((anaE-mediaE)/anaE)*100)
    
    # evolucao - NT
    
    nt <- listasim$nt
    
    ns <- listasim$ns
    
    nc <- listasim$nc
    
    dt <- listasim$dt
    
    
    
    tempo<- seq( dt*nc*100    ,   (nt/100)*dt*nc*100  ,  dt*nc*100 )
    
    
    # histogramas
    
    # saida da SimulacalBimodal
    # final<-list("posicao"=xx,"velocidade"=vv, "ji" = ji , "jd" = jd ,"energia"=ee, "massa"=massa,"gamma"=gamma,"k"=ctelastica, "alfa" = alfazero,"dalfa"= dalfa, "xa"=xa, "dxa"=dxa, "ns" = ns ,"nt"=nt,"nc"=nc,"dt"= dt, "temp"=Temp)
    
    
    # stoth<-numeric(listasim$ns)
    
    stoth<- numeric(nt/100)
    
    etoth<- numeric(nt/100)
    
    Ett<- numeric(nt/100)
    
    jji<- numeric(nt/100)
    
    jjd<- numeric(nt/100)
    
    for(j in 1:ns){
        
        for( i in 1:(nt/100) ){
            
                print(i)
                
                hhx<- hist(listasim$posicao[1:(i*100) , j], breaks = 100, plot=FALSE)
                hhv<- hist(listasim$velocidade[1:(i*100) , j], breaks = 100, plot=FALSE)
                
                #a ideia seria mudar a formula da entropia aqui - JAYNES
                
                ssx<- - hhx$density*log(hhx$density)
                ssx[which(!is.finite(ssx))]<- 0
                sx<-sum(ssx)*(1/sum(hhx$density))
                
                ssv<- - hhv$density*log(hhv$density)
                ssv[which(!is.finite(ssv))]<- 0
                sv<-sum(ssv)*(1/sum(hhv$density))
                
                stoth[i]<- stoth[i] + sx + sv
            
                etoth[i]<- etoth[i] + mean(listasim$energia[1:(i*100),j])
                
                Ett[i] <- etoth[i] + listasim$energia[(i*100),j]
                
                jji[i] <- jji[i] + listasim$ji[(i*100),j]
            
                jjd[i] <- jjd[i] + listasim$jd[(i*100),j]
        }
        
    }
    
    jjd <- jjd/ns
    
    jji <- jji/ns
    
    Ett<- Ett/ns
    
    stoth<- stoth/ns
    
    etoth<- etoth/ns
    
    
    
    
    fim<-list("tempo" = tempo  ,"S" = stoth, "Etempo" = Ett , "jd" = jjd   ,"ji" = jji , "Et" = etoth,  "E" = mediaE, "AnaE" = anaE, "erroE" = erroE , "alfa" = alfai, "xa" = xai,  "ns" = listasim$ns , "nt" = listasim$nt, "nc" = listasim$nc,"dt" = listasim$dt)
    
    
    return(fim)
    
}

evolucaoentropica3D<-function(listasim, breakks){
    
    
    if(listasim$dalfa==0)
    {
        alfai<-rep(listasim$alfa, listasim$ns)
    } else {
        alfai<-seq(listasim$alfa+listasim$dalfa,listasim$alfa+listasim$ns*listasim$dalfa,listasim$dalfa)
    }
    
    if(listasim$dxa==0)
    {
        xai<-rep(listasim$xa, listasim$ns)
    } else {
        xai<-seq(listasim$xa+listasim$dxa,listasim$xa+listasim$ns*listasim$dxa,listasim$dxa)
    }
    
    mediaE<-numeric(listasim$ns)
    sdE<-numeric(listasim$ns)
    
    #ENERGIA MÉDIA ANALÍTICA
    
    anaE<- (((xai)^2 )*(2*listasim$massa*alfai +listasim$gamma))/(2*listasim$gamma*(listasim$k +listasim$massa *( alfai^2) + alfai * listasim$gamma))
    
    
    
    for(i in 1:listasim$ns){
        
        mediaE[i]<-mean(listasim$energia[,i])
        
        sdE[i]<-sd(listasim$energia[,i])
        
    }
    
    erroE<- abs(((anaE-mediaE)/anaE)*100)
    
    # evolucao - NT
    
    nt <- listasim$nt
    
    ns <- listasim$ns
    
    nc <- listasim$nc
    
    dt <- listasim$dt
    
    
    
    tempo<- seq( dt*nc*100    ,   (nt/100)*dt*nc*100  ,  dt*nc*100 )
    
    
    
    
    # histogramas
    
    #  hist3d
    
    # saida da SimulacalBimodal
    # final<-list("posicao"=xx,"velocidade"=vv, "ji" = ji , "jd" = jd ,"energia"=ee, "massa"=massa,"gamma"=gamma,"k"=ctelastica, "alfa" = alfazero,"dalfa"= dalfa, "xa"=xa, "dxa"=dxa, "ns" = ns ,"nt"=nt,"nc"=nc,"dt"= dt, "temp"=Temp)
    
    
    # stoth<-numeric(listasim$ns)
    
    stoth<- numeric(nt/100)
    
    stothx<- numeric(nt/100)
    
    stothv<- numeric(nt/100)
    
    etoth<- numeric(nt/100)
    
    Ett<- numeric(nt/100)
    
    jji<- numeric(nt/100)
    
    jjd<- numeric(nt/100)
    
    # for(j in 1:ns){
        
        #print(j)
        
        for( i in 1:(nt/100) ){
            
            print(i)
            
            hhx<- hist(listasim$posicao[1:(i*100) , ], breaks = breakks, plot = FALSE)
            hhv<- hist(listasim$velocidade[1:(i*100) , ], breaks = breakks, plot = FALSE)
            hh<- Histograma3D(listasim$posicao[1:(i*100) , ], listasim$velocidade[1:(i*100) , ], breakks )
            
            ss<- Entropia3D(hh)
            #a ideia seria mudar a formula da entropia aqui - JAYNES
            
            ssx<- - hhx$density*log(hhx$density)
            ssx[which(!is.finite(ssx))]<- 0
            sx<-sum(ssx)*(1/sum(hhx$density))
            #print(sx)
            
            ssv<- - hhv$density*log(hhv$density)
            ssv[which(!is.finite(ssv))]<- 0
            sv<-sum(ssv)*(1/sum(hhv$density))
            #print(sv)
            #print(stothx)
            #print(stothv)
            
            stoth[i]<- stoth[i] + ss
            
            stothx[i]<- stothx[i] + sx
            
            stothv[i]<- stothv[i] + sv
            
            etoth[i]<- etoth[i] + mean(listasim$energia[1:(i*100),])
            
            Ett[i] <- etoth[i] + mean(listasim$energia[(i*100),])
            
            jji[i] <- jji[i] + mean(listasim$ji[(i*100),])
            
            jjd[i] <- jjd[i] + mean(listasim$jd[(i*100),])
            
            # }
        
    }
    
    jjd <- jjd/ns
    
    jji <- jji/ns
    
    Ett <- Ett/ns
    
    stoth <- stoth/ns
    
    stothx <- stothx/ns
    
    stothv <- stothv/ns
    
    etoth <- etoth/ns



    
    fim<-list("tempo" = tempo  ,"Sx" = stothx, "Sv" = stothv ,"S" = stoth, "Etempo" = Ett , "jd" = jjd   ,"ji" = jji , "Et" = etoth,  "E" = mediaE, "AnaE" = anaE, "erroE" = erroE , "alfa" = alfai, "xa" = xai,  "ns" = listasim$ns , "nt" = listasim$nt, "nc" = listasim$nc,"dt" = listasim$dt)
    
    
    return(fim)
}

evolucaoentropica3DLan<-function(listasim, breakks){
    
    

    nt <- listasim$nt
    
    ns <- listasim$ns
    
    nc <- listasim$nc
    
    dt <- listasim$dt
    
    
    
    tempo<- seq( dt*nc*100    ,   (nt/100)*dt*nc*100  ,  dt*nc*100 )

    
    stoth<- numeric(nt/100)
    
    stothx<- numeric(nt/100)
    
    stothv<- numeric(nt/100)
    
    etoth<- numeric(nt/100)
    
    Ett<- numeric(nt/100)
    
    jji<- numeric(nt/100)
    
    jjd<- numeric(nt/100)
    

    
    for( i in 1:(nt/100) ){
        
        print(i)
        
        hhx<- hist(listasim$posicao[1:(i*100) , ], breaks = breakks, plot = FALSE)
        hhv<- hist(listasim$velocidade[1:(i*100) , ], breaks = breakks, plot = FALSE)
        hh<- Histograma3D(listasim$posicao[1:(i*100) , ], listasim$velocidade[1:(i*100) , ], breakks )
        
        ss<- Entropia3D(hh)

        
        ssx<- - hhx$density*log(hhx$density)
        ssx[which(!is.finite(ssx))]<- 0
        sx<-sum(ssx)*(1/sum(hhx$density))

        ssv<- - hhv$density*log(hhv$density)
        ssv[which(!is.finite(ssv))]<- 0
        sv<-sum(ssv)*(1/sum(hhv$density))

        
        stoth[i]<- stoth[i] + ss
        
        stothx[i]<- stothx[i] + sx
        
        stothv[i]<- stothv[i] + sv
        
        etoth[i]<- etoth[i] + mean(listasim$energia[1:(i*100),])
        
        Ett[i] <- etoth[i] + mean(listasim$energia[(i*100),])
        
        jji[i] <- jji[i] + mean(listasim$ji[(i*100),])
        
        jjd[i] <- jjd[i] + mean(listasim$jd[(i*100),])
        

        
    }
    
    jjd <- jjd/ns
    
    jji <- jji/ns
    
    Ett <- Ett/ns
    
    stoth <- stoth/ns
    
    stothx <- stothx/ns
    
    stothv <- stothv/ns
    
    etoth <- etoth/ns
    
    
    fim<-list("tempo" = tempo  ,"Sx" = stothx, "Sv" = stothv ,"S" = stoth, "Etempo" = Ett , "jd" = jjd   ,"ji" = jji , "Et" = etoth,    "ns" = listasim$ns , "nt" = listasim$nt, "nc" = listasim$nc,"dt" = listasim$dt)
    
    
    return(fim)
}




EvolucaoEntropicaQ<-function(listasim, Q)
{
    
    
    
    if(listasim$dalfa==0)
    {
        alfai<-rep(listasim$alfa, listasim$ns)
    } else {
        alfai<-seq(listasim$alfa+listasim$dalfa,listasim$alfa+listasim$ns*listasim$dalfa,listasim$dalfa)
    }
    
    if(listasim$dxa==0)
    {
        xai<-rep(listasim$xa, listasim$ns)
    } else {
        xai<-seq(listasim$xa+listasim$dxa,listasim$xa+listasim$ns*listasim$dxa,listasim$dxa)
    }
    
    mediaE<-numeric(listasim$ns)
    sdE<-numeric(listasim$ns)
    
    #ENERGIA MÉDIA ANALÍTICA
    
    anaE<- (((xai)^2 )*(2*listasim$massa*alfai +listasim$gamma))/(2*listasim$gamma*(listasim$k +listasim$massa *( alfai^2) + alfai * listasim$gamma))
    
    
    
    for(i in 1:listasim$ns){
        
        mediaE[i]<-mean(listasim$energia[,i])
        
        sdE[i]<-sd(listasim$energia[,i])
        
    }
    
    erroE<- abs(((anaE-mediaE)/anaE)*100)
    
    # evolucao - NT
    
    nt <- listasim$nt
    
    ns <- listasim$ns
    
    nc <- listasim$nc
    
    dt <- listasim$dt
    
    
    
    tempo<- seq( dt*nc*100    ,   (nt/100)*dt*nc*100  ,  dt*nc*100 )
    
    
    # histogramas
    
    # saida da SimulacalBimodal
    # final<-list("posicao"=xx,"velocidade"=vv, "ji" = ji , "jd" = jd ,"energia"=ee, "massa"=massa,"gamma"=gamma,"k"=ctelastica, "alfa" = alfazero,"dalfa"= dalfa, "xa"=xa, "dxa"=dxa, "ns" = ns ,"nt"=nt,"nc"=nc,"dt"= dt, "temp"=Temp)
    
    
    # stoth<-numeric(listasim$ns)
    
    stoth<- numeric(nt/100)
    
    etoth<- numeric(nt/100)
    
    Ett<- numeric(nt/100)
    
    jji<- numeric(nt/100)
    
    jjd<- numeric(nt/100)
    
    for(j in 1:ns){
        
        for( i in 1:(nt/100) ){
            
            print(i)
            
            hhx<- hist(listasim$posicao[1:(i*100) , j], breaks = 100, plot=FALSE)
            hhv<- hist(listasim$velocidade[1:(i*100) , j], breaks = 100, plot=FALSE)
            
            #a ideia seria mudar a formula da entropia aqui - JAYNES
            
            #  ssx<- - hhx$density*log(hhx$density)
            #ssx[which(!is.finite(ssx))]<- 0
            #sx<-sum(ssx)*(1/sum(hhx$density))
            
            # ssv<- - hhv$density*log(hhv$density)
            #ssv[which(!is.finite(ssv))]<- 0
            #sv<-sum(ssv)*(1/sum(hhv$density))
            
            sx<- (1- sum( hhx$density^Q) )/(Q-1)
            sv<- (1- sum( hhv$density^Q) )/(Q-1)
            
            stoth[i]<- stoth[i] + sx + sv
            
            etoth[i]<- etoth[i] + mean(listasim$energia[1:(i*100),j])
            
            Ett[i] <- etoth[i] + listasim$energia[(i*100),j]
            
            #
            
            
            jji[i] <- jji[i] + listasim$ji[(i*100),j]
            
            jjd[i] <- jjd[i] + listasim$jd[(i*100),j]
        }
        
    }
    
    jjd <- jjd/ns
    
    jji <- jji/ns
    
    Ett<- Ett/ns
    
    stoth<- stoth/ns
    
    etoth<- etoth/ns
    
    
    
    
    fim<-list("tempo" = tempo  ,"S" = stoth, "Etempo" = Ett , "jd" = jjd   ,"ji" = jji , "Et" = etoth,  "E" = mediaE, "AnaE" = anaE, "erroE" = erroE , "alfa" = alfai, "xa" = xai,  "ns" = listasim$ns , "nt" = listasim$nt, "nc" = listasim$nc,"dt" = listasim$dt)
    
    
    
    
    
  
}




TesteEntropicoLangevin<-function(){
    
    
    qero<-SimulacaoLangevin(1000,20000,100,0.00001,1,1,1,1)
    
    wat<-evolucaoentropicaLangevin(qero)
    
    SalvaCorrelacao(wat$S,'EntropiaLangevin.csv')
    SalvaCorrelacao(wat$Et,'EnergiaLangevin.csv')
    SalvaCorrelacao(wat$injetada,'InjetadaLangevin.csv')
    SalvaCorrelacao(wat$dissipada,'DissipadaLangevin.csv')
    

}

evolucaoentropicaLangevin<-function(listasim){
    
    #   final<-list("posicao"=xx,"velocidade"=vv, "ji"= ji, "jd"=jd,"energia"=ee,  "massa"=massa,"gamma"=gamma,"k"=ctelastica, "ns" = ns ,"nt"=nt,"nc"=nc,"dt"= dt, "temp"=temp)
    
    #SimulacaoLangevin <- function(ns, nt, nc, dt, massa, gamma, ctelastica, temp){
    
    mediaE<-numeric(listasim$ns)
    sdE<-numeric(listasim$ns)
    
    #ENERGIA MÉDIA ANALÍTICA
    
    #anaE<- (((xai)^2 )*(2*listasim$massa*alfai +listasim$gamma))/(2*listasim$gamma*(listasim$k +listasim$massa *( alfai^2) + alfai * listasim$gamma))
    
    #anaE
    
    for(i in 1:listasim$ns){
        
        mediaE[i]<-mean(listasim$energia[,i])
        
        sdE[i]<-sd(listasim$energia[,i])
        
    }
    
    #erroE<- abs(((anaE-mediaE)/anaE)*100)
    
    # evolucao - NT
    
    nt <- listasim$nt
    
    ns<- listasim$ns
    # histogramas
    
    
    
    
    # stoth<-numeric(listasim$ns)
    
    stoth<- numeric(nt/100)
    
    etoth<- numeric(nt/100)
    
    Ett<- numeric(nt/100)
    
    jii<- numeric(nt/100)
    
    jdd<- numeric(nt/100)
    
    for(j in 1:ns){
        
        for( i in 1:(nt/100) ){
            
            print(i)
            
            hhx<- hist(listasim$posicao[1:(i*100) , j], breaks = 100, plot=FALSE)
            hhv<- hist(listasim$velocidade[1:(i*100) , j], breaks = 100, plot=FALSE)
            
            #a ideia seria mudar a formula da entropia aqui - JAYNES
            
            ssx<- - hhx$density*log(hhx$density)
            ssx[which(!is.finite(ssx))]<- 0
            sx<-sum(ssx)*(1/sum(hhx$density))
            
            ssv<- - hhv$density*log(hhv$density)
            ssv[which(!is.finite(ssv))]<- 0
            sv<-sum(ssv)*(1/sum(hhv$density))
            
            stoth[i]<- stoth[i] + sx + sv
            
            etoth[i]<- etoth[i] + mean(listasim$energia[1:(i*100),j])
            
            Ett[i] <- etoth[i] + listasim$energia[(i*100),j]
            
            jii[i]<- jii[i] + listasim$ji[(i*100),j]
            
            jdd[i]<- jdd[i] + listasim$jd[(i*100),j]
            
        }
        
    }
    
    Ett<- Ett/ns
    
    stoth<- stoth/ns
    
    etoth<- etoth/ns
    
    jii<-jii/ns
    
    jdd<-jdd/ns
    
    fim<-list("S" = stoth, "Etempo" = Ett ,"injetada" = jii, "dissipada" = jdd,   "Et" = etoth,  "E" = mediaE,  "ns" = listasim$ns , "nt" = listasim$nt, "nc" = listasim$nc,"dt" = listasim$dt)
    
    
    return(fim)
    
}

#evolucaoentropicaLangevin<-function(listasim){


evolucaoentropicaLangevinQ<-function(listasim, Q){
    
    #   final<-list("posicao"=xx,"velocidade"=vv, "ji"= ji, "jd"=jd,"energia"=ee,  "massa"=massa,"gamma"=gamma,"k"=ctelastica, "ns" = ns ,"nt"=nt,"nc"=nc,"dt"= dt, "temp"=temp)
    
    #SimulacaoLangevin <- function(ns, nt, nc, dt, massa, gamma, ctelastica, temp){
    
    mediaE<-numeric(listasim$ns)
    sdE<-numeric(listasim$ns)
    
    #ENERGIA MÉDIA ANALÍTICA
    
    #anaE<- (((xai)^2 )*(2*listasim$massa*alfai +listasim$gamma))/(2*listasim$gamma*(listasim$k +listasim$massa *( alfai^2) + alfai * listasim$gamma))
    
    #anaE
    
    for(i in 1:listasim$ns){
        
        mediaE[i]<-mean(listasim$energia[,i])
        
        sdE[i]<-sd(listasim$energia[,i])
        
    }
    
    #erroE<- abs(((anaE-mediaE)/anaE)*100)
    
    # evolucao - NT
    
    nt <- listasim$nt
    
    ns<- listasim$ns
    # histogramas
    
    
    
    
    # stoth<-numeric(listasim$ns)
    
    stoth<- numeric(nt/100)
    
    etoth<- numeric(nt/100)
    
    Ett<- numeric(nt/100)
    
    jii<- numeric(nt/100)
    
    jdd<- numeric(nt/100)
    
    for(j in 1:ns){
        
        for( i in 1:(nt/100) ){
            
            print(i)
            
            hhx<- hist(listasim$posicao[1:(i*100) , j], breaks = 100, plot=FALSE)
            hhv<- hist(listasim$velocidade[1:(i*100) , j], breaks = 100, plot=FALSE)

            
            #a ideia seria mudar a formula da entropia aqui - JAYNES
            
            #  ssx<- - hhx$density*log(hhx$density)
            #ssx[which(!is.finite(ssx))]<- 0
            #sx<-sum(ssx)*(1/sum(hhx$density))
            
            # ssv<- - hhv$density*log(hhv$density)
            #ssv[which(!is.finite(ssv))]<- 0
            #sv<-sum(ssv)*(1/sum(hhv$density))
            
            sx<- (1- sum( hhx$density^Q) )/(Q-1)
            sv<- (1- sum( hhv$density^Q) )/(Q-1)
            
            
            stoth[i]<- stoth[i] + sx + sv
            
            etoth[i]<- etoth[i] + mean(listasim$energia[1:(i*100),j])
            
            Ett[i] <- etoth[i] + listasim$energia[(i*100),j]
            
            jii[i]<- jii[i] + listasim$ji[(i*100),j]
            
            jdd[i]<- jdd[i] + listasim$jd[(i*100),j]
            
        }
        
    }
    
    Ett<- Ett/ns
    
    stoth<- stoth/ns
    
    etoth<- etoth/ns
    
    jii<-jii/ns
    
    jdd<-jdd/ns
    
    fim<-list("S" = stoth, "Etempo" = Ett ,"injetada" = jii, "dissipada" = jdd,   "Et" = etoth,  "E" = mediaE,  "ns" = listasim$ns , "nt" = listasim$nt, "nc" = listasim$nc,"dt" = listasim$dt)
    
    
    return(fim)
    
}



TesteFuncaoResposta <- function(){
    
    #source('Entropiaultima.R')
    
    corrpre <- '3corrpre'
    corrpos <- '3corrpos'
    rresposta <- '3rresposta'
    
    csv<-'.csv'
    
    for (i in 1:90)
    {
    
        titulo0<-paste( corrpre , toString(i) , csv , sep="")
    
        titulo1<-paste( corrpos , toString(i) , csv , sep="")
    
        tituloR<-paste( rresposta , toString(i) , csv , sep="")
    
    #SimulacaoBimodalKickF <- function(ns, nt, nc, dt, massa, gamma, ctelastica, alfa, dalfa, xa, dxa, deltax, tt)
    
    # PRMEIRA FOI 10000,  AGOR VOU FAZER 100000
    #1  qero<-SimulacaoBimodalKickF(5000,10000,100,0.001,1,1,1,1,0,1,0,10000, (10*i) + 50 )
    #2  qero<-SimulacaoBimodalKickF(5000,10000,100,0.001,1,1,1,1,0,1,0,100000, (10*i) + 50 )
    #3  qero<-SimulacaoBimodalKickX(5000,10000,100,0.001,1,1,1,1,0,1,0,10, (10*i) + 50 )   -> esse ficou escroto
    #3.2  qero<-SimulacaoBimodalKickF(5000,10000,100,0.001,1,1,1,1,0,1,0,50000, (10*i) + 50 )
    
    
        qero<-SimulacaoBimodalKickF(5000,10000,10000,0.00001,1,1,1,1,0,1,0,50000, (10*i) + 50 )
    
        ns<- qero$ns
    
        kick<- qero$kick
        tkick<- qero$tkick /(qero$dt*qero$nc)
    
        Dt<- 100
    
        fim<- numeric(Dt)
    
        corpre<-      (  funcaodecorrelacaoSample(qero$posicao, (tkick-Dt -1) , Dt, ns)   )
    
        corpos<- ( funcaodecorrelacaoSample(qero$posicao, tkick , Dt ,ns))

    
    
        fim<- fim+ (corpos - corpre)
    
        fim<-fim/kick
    
        fim<- fim/ns
    
    
        SalvaCorrelacao(corpre, titulo0)
    
        SalvaCorrelacao(corpos, titulo1)
    
        SalvaCorrelacao(fim, tituloR)
    
    }
    
}



TesteFormula2<- function(nome='default'){
    
    alpha <- seq(0.5,6,0.5)

    s <- numeric(length(alpha))
    e <- numeric(length(alpha))
    s2 <- numeric(length(alpha))
    e2 <- numeric(length(alpha))

    dsde <- numeric(length(alpha))
    
    for(i in 1:length(alpha))
    {
        print(i)
        qero <- SimulacaoBimodalTcte(10000,1000,1000,0.001,1,1,1, alpha[i],0,1  )
 	qero2 <- SimulacaoBimodalTcte(10000,1000,1000,0.001,1,1,1, alpha[i]+0.1,0,1  )
        wat <- analiseentropicaSIMPLISTT(qero)
        wat2<- analiseentropicaSIMPLISTT(qero2)

        #fim<-list("S"=stoth,"E"=mediaE, "PolyE"=anaE, "alfa"=alfai, "xa"=xai,  "ns"=listasim$ns,"nt"=listasim$nt,"nc"=listasim$nc,"dt"= listasim$dt)
        #SimulacaoBimodalTcte<- function(ns, nt, nc, dt, massa, gamma, ctelastica, alfa, dalfa, Temp){
        
        s[i] <- mean(wat$S)
        e[i] <- mean(wat$E)

	s2[i] <- mean(wat2$S)
        e2[i] <- mean(wat2$E)
    }
    
    nomes <- paste(nome,'entropia.csv', sep='')
    nomes2 <- paste(nome,'entropia2.csv', sep='')

    nomee <- paste(nome,'energia.csv', sep='')
    nomee2 <- paste(nome,'energia2.csv', sep='')

    nomedsde <- paste(nome,'dsde.csv', sep='')
    nomealpha <- paste(nome,'alpha.csv', sep='')
    
    for(i in 1:length(alpha))
    {
        dsde[i] <- (s2[i]-s[i])/(e2[i]-e[i])
    }
    
    SalvaCorrelacao(s, nomes)
    SalvaCorrelacao(e, nomee)

    SalvaCorrelacao(s2, nomes2)
    SalvaCorrelacao(e2, nomee2)

    SalvaCorrelacao(dsde, nomedsde)
    SalvaCorrelacao(alpha, nomealpha)

    
}

TesteFormula<- function(nome='default'){
    
    alpha <- seq(0.5,6,0.5)

    s <- numeric(length(alpha))
    e <- numeric(length(alpha))


    dsde <- numeric(length(alpha))
    
    for(i in 1:length(alpha))
    {
        print(i)
        qero <- SimulacaoBimodalTcte(12000,1000,1000,0.001,1,1,1, alpha[i],0,1  )
        wat <- analiseentropicaSIMPLISTT(qero)


        #fim<-list("S"=stoth,"E"=mediaE, "PolyE"=anaE, "alfa"=alfai, "xa"=xai,  "ns"=listasim$ns,"nt"=listasim$nt,"nc"=listasim$nc,"dt"= listasim$dt)
        #SimulacaoBimodalTcte<- function(ns, nt, nc, dt, massa, gamma, ctelastica, alfa, dalfa, Temp){
        
        s[i] <- mean(wat$S)
        e[i] <- mean(wat$E)

    }
    
    nomes <- paste(nome,'entropia.csv', sep='')
    nomee <- paste(nome,'energia.csv', sep='')
    nomedsde <- paste(nome,'dsde.csv', sep='')
    nomealpha <- paste(nome,'alpha.csv', sep='')
    
    for(i in 2:length(alpha))
    {
        dsde[i] <- (s[i]-s[i-1])/(e[i]-e[i-1])
    }
    
    SalvaCorrelacao(s, nomes)
    SalvaCorrelacao(e, nomee)
    SalvaCorrelacao(dsde, nomedsde)
    SalvaCorrelacao(alpha, nomealpha)
    
}

Resposta2 <- function(qero, Dt){
    
    ns<- qero$ns
    kick<- qero$kick
    
    tkick<- qero$tkick/(qero$dt*qero$nc)
    print(tkick)
    
    respostaa<- numeric(Dt)
    suceptibilidade <- numeric(Dt)
    XTeff<- numeric(Dt)
    
    print('1foi até aqui')
    print('passou o corrsample')
    
    for (i in 1:Dt)
    {
        mediapos <-  mean(qero$posicao[tkick+i, ]  )
        respostaa[i] <- mediapos/kick
        suceptibilidade[i] <- sum(respostaa[1:i]) * qero$dt
    }
    
    #iii <-  iii + 1
    
    assimfdtpre<- AssimetriaFDT(qero, tkick-Dt-1, Dt)
    assimfdtpos<- AssimetriaFDT(qero, tkick, Dt)
    
    assimpre <-Assimetria(qero, tkick-Dt-1, Dt)
    assimpos <-Assimetria(qero, tkick, Dt)
    
    assimtotpre <-Assimetriatot(qero, tkick-Dt-1, Dt)
    assimtotpos <-Assimetriatot(qero, tkick, Dt)
    
    # assim
    
    
    final<-list("assimfdtpre" = assimfdtpre, "assimfdtpos" = assimfdtpos ,"assimpre" = assimpre,"assimpos" = assimpos ,"resposta" = respostaa)
    
    return(final)
    
    #return(assimpre)
    
}


Resposta1 <- function(qero, Dt){
    
    ns<- qero$ns
    kick<- qero$kick
    
    tkick<- qero$tkick/(qero$dt*qero$nc)
    print(tkick)
    
    respostaa<- numeric(Dt)
    suceptibilidade <- numeric(Dt)
    XTeff<- numeric(Dt)
    
    print('1foi até aqui')
    print('passou o corrsample')
    
    for (i in 1:Dt)
    {
        mediapos <-  mean(qero$posicao[tkick+i, ]  )
        respostaa[i] <- mediapos/kick
        #suceptibilidade[i] <- sum(respostaa[1:i]) * qero$dt
    }
    
    #iii <-  iii + 1
    # assim
    
    
    final<-list("resposta" = respostaa)
    
    return(final)
    
    #return(assimpre)
    
}

Resposta3 <- function(qero, Dt){
    
    ns<- qero$ns
    kick<- qero$kick
    
    tkick<- qero$tkick/(qero$dt*qero$nc)
    print(tkick)
    
    respostaa<- numeric(Dt)
    suceptibilidade <- numeric(Dt)
    XTeff<- numeric(Dt)
    
    print('1foi até aqui')
    print('passou o corrsample')
    
    for (i in 1:Dt)
    {
        mediapos <-  mean(qero$posicao[tkick+i, ]  )
        respostaa[i] <- mediapos/kick
        #suceptibilidade[i] <- sum(respostaa[1:i]) * qero$dt
    }
    
    #iii <-  iii + 1
    # assim
    
    
    final<-list("resposta" = respostaa)
    
    return(final)
    
    #return(assimpre)
    
}


TesteFuncaoRespostaF <- function(paralel, nsamples){
    

        
        
        Temp <- 1
        
        alpha <- 1
        
        # xa
        
        # xa<- sqrt( Temp*(1*(1+alpha[kk]*(1+1*alpha[kk])))/(1*alpha[kk])   )
        
        xa<- sqrt( Temp*(1*(1+alpha*(1+1*alpha)))/(1*alpha)   )
        
        # qero<-SimulacaoBimodalKickF(nsamples,2510,10,0.001,1,1,1, alpha[kk] ,0, xa, 0, 0.5,  11 )
        
        qero<-SimulacaoBimodalKickF(nsamples,2510,10,0.001,1,1,1, alpha ,0, xa, 0, 0.5,  11 )
        
        
        Dt<- 1000
        

        
        tempo <- seq(0.01, 10, 0.01)
        
        vez<- 1
        
        letra <- paralel
        
        corrpre <- paste( letra ,  'corrpre', toString(vez)  , sep="")
        rresposta <- paste( letra, 'resposta', toString(vez) , sep="")
        
        sucep <- paste( letra, 'sucep', toString(vez) , sep="")
        Xteff <- paste( letra, 'XTeff', toString(vez) , sep="")
        
        Assimpre <- paste( letra, 'Assimpre', toString(vez) , sep="")
        Assimpos <- paste( letra, 'Assimpos', toString(vez) , sep="")
        
        csv<-'.csv'
        
        titulo0<-paste( corrpre , csv , sep="")
        tituloR<-paste( rresposta , csv , sep="")
        tituloX<-paste( sucep , csv , sep="")
        tituloT<-paste( Xteff , csv , sep="")
        tituloApre<-paste( Assimpre , csv , sep="")
        tituloApos<-paste( Assimpos , csv , sep="")
        
        ns<- qero$ns
        kick<- qero$kick
        
        tkick<- qero$tkick/(qero$dt*qero$nc)
        print(tkick)
        
        respostaa<- numeric(Dt)
        suceptibilidade <- numeric(Dt)
        XTeff<- numeric(Dt)
        
        print('1foi até aqui')
        corpre<-      (  funcaodecorrelacaoSample(qero$posicao, (tkick - Dt - 1) , Dt , ns)   )
        print('passou o corrsample')
        
        for (i in 1:Dt)
        {
            mediapos <-  mean(qero$posicao[tkick+i, ]  )
            respostaa[i] <- mediapos/kick
            #suceptibilidade[i] <- sum(respostaa[1:i])
            suceptibilidade[i] <- sum(respostaa[1:i]) * qero$dt
        }
        
        # iii <-  iii + 1
        
        print('2foi até aqui')
        #assim<- Assimetria(qero, 100, 500)
        
        #NOVIDADE AQUI ASSIMETRIA
        
        print(tkick-Dt-1)
        print(tkick+Dt)
        
        #  assimfdtpre<- AssimetriaFDT(qero, tkick-Dt-1, Dt)
        # assimfdtpos<- AssimetriaFDT(qero, tkick, Dt)
        
        print('chegao no save')
        
        #SalvaCorrelacao(assimfdtpre, tituloApre)
        # SalvaCorrelacao(assimfdtpos, tituloApos)
        
      
        SalvaCorrelacao(corpre, titulo0)
        
        SalvaCorrelacao(suceptibilidade, tituloX)
        
        SalvaCorrelacao(respostaa, tituloR)

        
        
        SalvaCorrelacao(XTeff, tituloT)
        
}

TesteFuncaoRespostaFextensoT <- function(paralel, nsamples){
    
    alpha <- seq(1.0,30.0,1.0)
    
    
    sequencia <- seq(1,30,1)
    
    T <- numeric(  length(sequencia)   )
    Tlan <-     numeric(  length(sequencia)   )
    alppha <- numeric(  length(sequencia)   )
    
    
    rsquared <- numeric(  length(sequencia)   )
    rsquaredlan <- numeric(  length(sequencia)   )
    
    iii <- 0
    
    
    for (kk in sequencia)
    {
        print(alpha[kk])
        print(kk)
    
    
    Temp <- 1
    
    #alpha <- 1
    
    xa<- sqrt( Temp*(1*(1+alpha[kk]*(1+1*alpha[kk])))/(1*alpha[kk])   )
    
    qero<-SimulacaoBimodalKickF(nsamples,2510,10,0.001,1,1,1, alpha[kk] ,0, xa, 0, 0.5,  11 )
    
    Dt<- 1000
    
    tempo <- seq(0.01, 10, 0.01)
    
    vez<- kk
    
    letra <- paralel
    
    corrpre <- paste( letra ,  'corrpre', toString(vez)  , sep="")
    rresposta <- paste( letra, 'resposta', toString(vez) , sep="")
    
    sucep <- paste( letra, 'sucep', toString(vez) , sep="")
    Xteff <- paste( letra, 'XTeff', toString(vez) , sep="")
    
    Assimpre <- paste( letra, 'Assimpre', toString(vez) , sep="")
    Assimpos <- paste( letra, 'Assimpos', toString(vez) , sep="")
    
    csv<-'.csv'
    
    titulo0<-paste( corrpre , csv , sep="")
    tituloR<-paste( rresposta , csv , sep="")
    tituloX<-paste( sucep , csv , sep="")
    tituloT<-paste( Xteff , csv , sep="")
    tituloApre<-paste( Assimpre , csv , sep="")
    tituloApos<-paste( Assimpos , csv , sep="")
    
    ns<- qero$ns
    kick<- qero$kick
    
    tkick<- qero$tkick/(qero$dt*qero$nc)
    print(tkick)
    
    respostaa<- numeric(Dt)
    suceptibilidade <- numeric(Dt)
    XTeff<- numeric(Dt)
    
    print('1foi até aqui')
    corpre<-      (  funcaodecorrelacaoSample(qero$posicao, (tkick - Dt - 1) , Dt , ns)   )
    print('passou o corrsample')
    
    for (i in 1:Dt)
    {
        mediapos <-  mean(qero$posicao[tkick+i, ]  )
        respostaa[i] <- mediapos/kick
        suceptibilidade[i] <- sum(respostaa[1:i]) * qero$dt
    }
    

    
    print('2foi até aqui')

    print(tkick-Dt-1)
    print(tkick+Dt)
    
  
    
    print('chegao no save')
    
    
    SalvaCorrelacao(corpre, titulo0)
    
    SalvaCorrelacao(suceptibilidade, tituloX)
    
    SalvaCorrelacao(respostaa, tituloR)
    
    SalvaCorrelacao(XTeff, tituloT)
    
    }
    
    
    
}


TesteFuncaoRespostaFextensoxa <- function(paralel, nsamples){
    
    #  alpha <- seq(1.0,30.0,1.0)
    
    alpha <- seq(0.05,1.5,0.05)
    
    sequencia <- seq(1,30,1)
    
    T <- numeric(  length(sequencia)   )
    Tlan <-     numeric(  length(sequencia)   )
    alppha <- numeric(  length(sequencia)   )
    
    
    rsquared <- numeric(  length(sequencia)   )
    rsquaredlan <- numeric(  length(sequencia)   )
    
    iii <- 0
    
    xa<- 1
    
    for (kk in sequencia)
    {
        print(alpha[kk])
        print(kk)
        
        
        #Temp <- 1
        
        #alpha <- 1
        
        #xa<- sqrt( Temp*(1*(1+alpha[kk]*(1+1*alpha[kk])))/(1*alpha[kk])   )
        
        
        qero<-SimulacaoBimodalKickF(nsamples,2510,10,0.001,1,1,1, alpha[kk] ,0, xa, 0, 0.5,  11 )
        
        Dt<- 1000
        
        tempo <- seq(0.01, 10, 0.01)
        
        vez<- kk
        
        letra <- paralel
        
        corrpre <- paste( letra ,  'corrpre', toString(vez)  , sep="")
        rresposta <- paste( letra, 'resposta', toString(vez) , sep="")
        
        sucep <- paste( letra, 'sucep', toString(vez) , sep="")
        Xteff <- paste( letra, 'XTeff', toString(vez) , sep="")
        
        Assimpre <- paste( letra, 'Assimpre', toString(vez) , sep="")
        Assimpos <- paste( letra, 'Assimpos', toString(vez) , sep="")
        
        csv<-'.csv'
        
        titulo0<-paste( corrpre , csv , sep="")
        tituloR<-paste( rresposta , csv , sep="")
        tituloX<-paste( sucep , csv , sep="")
        tituloT<-paste( Xteff , csv , sep="")
        tituloApre<-paste( Assimpre , csv , sep="")
        tituloApos<-paste( Assimpos , csv , sep="")
        
        ns<- qero$ns
        kick<- qero$kick
        
        tkick<- qero$tkick/(qero$dt*qero$nc)
        print(tkick)
        
        respostaa<- numeric(Dt)
        suceptibilidade <- numeric(Dt)
        XTeff<- numeric(Dt)
        
        print('1foi até aqui')
        corpre<-      (  funcaodecorrelacaoSample(qero$posicao, (tkick - Dt - 1) , Dt , ns)   )
        print('passou o corrsample')
        
        for (i in 1:Dt)
        {
            mediapos <-  mean(qero$posicao[tkick+i, ]  )
            respostaa[i] <- mediapos/kick
            suceptibilidade[i] <- sum(respostaa[1:i]) * qero$dt
        }
        
        
        
        print('2foi até aqui')
        
        print(tkick-Dt-1)
        print(tkick+Dt)
        
        
        
        print('chegao no save')
        
        
        SalvaCorrelacao(corpre, titulo0)
        
        SalvaCorrelacao(suceptibilidade, tituloX)
        
        SalvaCorrelacao(respostaa, tituloR)
        
        SalvaCorrelacao(XTeff, tituloT)
        
    }
    
    
    
}





TesteFuncaoResposta2Tcte <- function(paralel, nsamples){
    
    #source('Entropiaultima.R')
    
    #rresposta <- 'trueresposta'
    #rresposta <- '2trueresposta'
    #rresposta <- '3trueresposta'
    # rresposta <- '4trueresposta'
    # rresposta <- '5trueresposta'
    # rresposta <- '6trueresposta'
    
    #1   qero<-SimulacaoBimodalKickF(10,10000,100,0.001,1,1,1,1,0,1,0,50000,  50 )
    #2   qero<-SimulacaoBimodalKickF(20,10000,100,0.001,1,1,1,2,0,1,0,50000,  50 )
    #3   qero<-SimulacaoBimodalKickF(20,10000,100,0.001,1,1,1,0.5,0,1,0,50000,  50 )
    
    #4   qero<-SimulacaoBimodalKickF(100,10000,100,0.001,1,1,1,1,0,1,0,500,  50 )
    
    #5
    #qero<-SimulacaoBimodalKickF(1000,10000,100,0.001,1,1,1,1,0,1,0,5000,  100 )
    
    #6
    #qero<-SimulacaoBimodalKickF(1000,10000,100,0.001,1,1,1,1,0,1,0,1000,  50 )
    
    #8   qero<-SimulacaoBimodalKickF(500,10000,100,0.001,1,1,1,1,0,1,0,500,  50 )
    
    #11 qero<-SimulacaoBimodalKickF(1000,10000,100,0.001,1,1,1,0.5,0,1,0,100,  50 )
    
    # 11 qero<-SimulacaoBimodalKickF(1000,10000,100,0.001,1,1,1,0.5,0,1,0,100,  100 )
    
    # 12
    
    # qero<-SimulacaoBimodalKickF(1000,10000,100,0.001,1,1,1,0.5,0,1,0,1000,  100 )
    
    # 14 function(ns, nt, nc, dt, massa, gamma, ctelastica, temp, deltaf, tt){
    
    # 4Lancorrpre'  qero<-SimulacaoLangevinKickF(1000,10000,100,0.001,1,1,1,0.5, 100,  100 )
    
    # 4Lancorrpre'
    
    # 5Lancorrpre'
    
    #qero<-SimulacaoLangevinKickF(500,10000,100,0.001,1,1,1,0.5, 100,  100 )
    
    # '6Lancorrpre'
        #qero<-SimulacaoLangevinKickF(500,10000,100,0.001,1,1,1,2, 100,  100 )
    
    #'Bcorrpre' NOVA SESSION
    
    # 1 qero<-SimulacaoBimodalKickF(500,2000,10,0.001,1,1,1,1,0,1,0, 500,  10 ) DT = 100
    
    #qero<-SimulacaoBimodalKickF(500,2000,10,0.001,1,1,1,1,0,1,0, 500,  10 )
    #Dt<-
    #vez<- 2
    
    #qero<-SimulacaoBimodalKickF(500,2010,10,0.001,1,1,1,1,0,1,0, 500,  10.1 )
    #Dt<-1000
    #vez<- 3
    
    #qero<-SimulacaoBimodalKickF(1000,2010,10,0.001,1,1,1,1,0,1,0, 500,  10.1 )
    #Dt<- 1000
    #vez<- 4
    
    #qero<-SimulacaoBimodalKickF(2000,2010,10,0.001,1,1,1,1,0,1,0, 100,  10.1 )
    #Dt<- 1000
    #vez<- 5
    
    #qero<-SimulacaoBimodalKickF(10000,2010,10,0.001,1,1,1,1,0,1,0, 100,  10.1 )
    #Dt<- 1000
    #vez<- 6
    
    #qero<-SimulacaoBimodalKickF(10000,3010,10,0.001,1,1,1,1,0,1,0, 100,  20.1 )
    #Dt<- 1000
    #vez<- 7
    
    #qero<-SimulacaoBimodalKickF(20000,3010,10,0.001,1,1,1,1,0,1,0, 100,  20.1 )
    #Dt<- 1000
    #vez<- 8
    
    #qero<-SimulacaoBimodalKickF(10000,3010,10,0.001,1,1,1,1,0,1,0, 50,  20.1 )
    #Dt<- 1000
    #vez<- 9
    
    #qero<-SimulacaoBimodalKickF(30000,3010,10,0.001,1,1,1,1,0,1,0, 50,  20.1 )
    #Dt<- 1000
    #vez<- 10
    
    #qero<-SimulacaoBimodalKickF(20000,3010,10,0.001,1,1,1,0.5,0,1,0, 50,  20.1 )
    #Dt<- 1000
    #vez<- 11
    
    #qero<-SimulacaoBimodalKickF(20000,3010,10,0.001,1,1,1,0.3,0,1,0, 50,  20.1 )
    #Dt<- 1000
    #vez<- 12
    
    #qero<-SimulacaoBimodalKickF(20000,3010,10,0.001,1,1,1,1.5,0,1,0, 50,  20.1 )
    ##Dt<- 1000
    #vez<- 13
    
    #qero<-SimulacaoBimodalKickF(100000,3010,10,0.001,1,1,1,1.5,0,1,0, 50,  20.1 )
    #Dt<- 1000
    #vez<- 14
    
    #qero<-SimulacaoBimodalKickF(30000,3010,10,0.001,1,1,1,2,0,1,0, 50,  20.1 )
    #Dt<- 1000
    #vez<- 15
    
    #qero<-SimulacaoBimodalKickF(100000,3010,10,0.001,1,1,1,1,0,1,0, 50,  20.1 )
    #Dt<- 1000
    #vez<- 16
    
    #qero<-SimulacaoBimodalKickF(100000,3010,10,0.001,1,1,1,1,0,1,0, 50,  20.1 )
    #Dt<- 1000
    #vez<- 17
    
    #qero<-SimulacaoBimodalKickF(50000,3010,10,0.001,1,1,1,2,0,1,0, 50,  20.1 )
    #Dt<- 1000
    #vez<- 18
    
    #qero<-SimulacaoBimodalKickF(200000,3010,10,0.001,1,1,1,0.7,0,1,0, 50,  20.1 )
    #Dt<- 1000
    #vez<- 19
    
    #qero<-SimulacaoBimodalKickF(110000,3010,10,0.001,1,1,1,0.9,0,1,0, 50,  20.1 )
    #Dt<- 1000
    #vez<- 20
    
    #qero<-SimulacaoBimodalKickF(40000,3010,10,0.001,1,1,1,0.4,0,1,0, 50,  20.1 )
    #Dt<- 1000
    #vez<- 21
    
    #qero<-SimulacaoBimodalKickF(40000,3010,10,0.001,1,1,1,1,0,1,0, 10,  20.1 ) # aqui mudou o kick!
    #Dt<- 1000
    #vez<- 22
    
    #qero<-SimulacaoBimodalKickF(40000,3010,10,0.001,1,1,1,1,0,1,0, 2,  20.1 ) # aqui mudou o kick!
    #Dt<- 1000
    #vez<- 23
    
    #qero<-SimulacaoBimodalKickF(40000,3010,10,0.001,1,1,1,1,0,1,0, 1,  20.1 ) # aqui mudou o kick!
    #Dt<- 1000
    #vez<- 24
    
    #qero<-SimulacaoBimodalKickF(40000,3010,10,0.001,1,1,1,2,0,1,0, 1,  20.1 ) # aqui mudou o kick!
    #Dt<- 1000
    #vez<- 25
    
    #qero<-SimulacaoBimodalKickF(40000,3010,10,0.001,1,1,1,2,0,1,0, 0.5,  20.1 ) # aqui mudou o kick!
    #Dt<- 1000
    #vez<- 26
    
    #qero<-SimulacaoBimodalKickF(40000,3010,10,0.001,1,1,1,2,0,0.1,0, 5,  20.1 ) # aqui mudou o kick!
    #Dt<- 1000
    #vez<- 27
    
    #qero<-SimulacaoBimodalKickF(100000,3010,10,0.001,1,1,1,1,0,0.1,0, 5,  20.1 ) # aqui mudou o kick!
    #Dt<- 1000
    #vez<- 28
    
    #qero<-SimulacaoBimodalKickF(100000,3010,10,0.001,1,1,1,1,0,0.05,0, 5,  20.1 ) # aqui mudou o kick!
    #Dt<- 1000
    #vez<- 29
    
    #qero<-SimulacaoBimodalKickF(50000,3010,10,0.001,1,1,1,1,0,0.05,0, 5,  20.1 ) # aqui mudou o kick!
    #Dt<- 1000
    #vez<- 30
    
    #qero<-SimulacaoBimodalKickF(50000,3010,10,0.001,1,1,1,0.5,0,0.05,0, 5,  20.1 ) # aqui mudou o kick!
    #Dt<- 1000
    #vez<- 31
    
    #qero<-SimulacaoBimodalKickF(50000,3010,10,0.001,1,1,1,2,0,0.02,0, 1,  20.1 ) # aqui mudou o kick!
    #Dt<- 1000
    #vez<- 32
    
    # qero<-SimulacaoBimodalKickF(50000,3010,10,0.001,1,1,1,5,0,0.1,0, 1,  20.1 ) # aqui mudou o kick!
    #Dt<- 1000
    #vez<- 33
    
    #qero<-SimulacaoBimodalKickX(50000,3010,10,0.001,1,1,1,5,0,1,0, 0.2,  20.1 ) # aqui mudou o kick! e virou kickx
    #Dt<- 1000
    #vez<- 34
    
    #qero<-SimulacaoBimodalKickX(50000,3010,10,0.001,1,1,1,0.5,0,1,0, 0.2,  20.1 ) # aqui mudou o alfa em relacao ao anterior
    #Dt<- 1000
    #vez<- 34
    
    #qero<-SimulacaoBimodalKickX(50000,3010,10,0.001,1,1,1,1,0,1,0, 0.1,  20.1 ) # aqui mudou o alfa em relacao ao anterior
    #Dt<- 1000
    #vez<- 35
    
    #qero<-SimulacaoBimodalKickX(50000,3010,10,0.001,1,1,1,1,0,1,0, 1,  20.1 ) # aqui mudou o alfa em relacao ao anterior
    #Dt<- 1000
    #vez<- 36
    
    #qero<-SimulacaoBimodalKickX(50000,3010,10,0.001,1,1,1,1,0,1,0, 0.5,  20.1 ) # aqui mudou o alfa em relacao ao anterior
    #Dt<- 1000
    #vez<- 37
    
    #qero<-SimulacaoBimodalKickX(50000,3010,10,0.001,1,1,1,1,0,2,0, 0.5,  20.1 ) # aqui mudou o alfa em relacao ao anterior
    #Dt<- 1000
    #vez<- 38
    
    #qero<-SimulacaoBimodalKickX(50000,3010,10,0.001,1,1,1,1.5,0,1.5,0, 0.5,  20.1 ) # aqui mudou o alfa em relacao ao anterior
    #Dt<- 1000
    #vez<- 39
    
    #qero<-SimulacaoBimodalKickX(50000,3010,10,0.001,1,1,1,1.5,0,1.5,0, 0.5,  20.1 ) # aqui mudou o alfa em relacao ao anterior
    #Dt<- 1000
    #vez<- 39
    
    ### salvamos os d corr analiticos previamente nos arquivos Dcorr1 ateh dcorr60 .csv
    #nsamples <- 1000
    
    #alphalto
    #alpha <- seq(5.1,8,0.1)
    
    #alphabemalto
    alpha <- seq(1,30,1)
    
    #sequencia <- seq(4,30,1)
    #alpha <- seq(3,9,0.1)
    
    sequencia <- seq(1,30,1)
    
    T <- numeric(  length(sequencia)   )
    Tlan <-     numeric(  length(sequencia)   )
    
    alppha <- numeric(  length(sequencia)   )
    
    
    rsquared <- numeric(  length(sequencia)   )
    rsquaredlan <- numeric(  length(sequencia)   )
    
    iii <- 0
    
    
    for (kk in sequencia)
    {
    print(alpha[kk])
    print(kk)
    
    # Temp = 0.710526   /. {k -> 1, \[Gamma] -> 1, m -> 1, a -> 1.5, \[Alpha] -> 1.5}   que funcionou
    #vamo tentar manter a temperatura
    
    #Temp <- 0.710526
    
    Temp <- 1
    
    
    #xa<- sqrt( Temp*(gamma*(ctelastica+alpha[kk]*(gamma+massa*alpha[kk])))/(massa*alpha[kk])   )   #-> será que isso garante o encaixe???
    
    xa<- sqrt( Temp*(1*(1+alpha[kk]*(1+1*alpha[kk])))/(1*alpha[kk])   )
    
    #qero<-SimulacaoBimodalKickX(nsamples,2510,2000,0.000005,1,1,1, alpha[kk] ,0, xa, 0, 0.5,  11 ) # aqui mudou o alfa em relacao ao anterior
    
    qero<-SimulacaoBimodalKickXsalvaEta(nsamples,2510,10,0.001,1,1,1, alpha[kk] ,0, xa, 0, 0.5,  11 )
    
    #qero<-SimulacaoBimodalKickV(nsamples,2510,10,0.001,1,1,1, alpha[kk] ,0, xa, 0, 0.5,  11 )
    
    Dt<- 1000
    
    #tempo<-  seq(0.002, 4, 0.002)
    
    tempo <- seq(0.01, 10, 0.01)
    
    vez<- kk
    
    letra <- paralel
    
    #corrpos <- paste(toString(vez) , 'Bcorrpos' , sep = ""  )
    
    corrpre <- paste( letra ,  'corrpre', toString(vez)  , sep="")
    rresposta <- paste( letra, 'resposta', toString(vez) , sep="")
    
    sucep <- paste( letra, 'sucep', toString(vez) , sep="")
    Xteff <- paste( letra, 'XTeff', toString(vez) , sep="")
    
    Assimpre <- paste( letra, 'Assimpre', toString(vez) , sep="")
    Assimpos <- paste( letra, 'Assimpos', toString(vez) , sep="")
    
    csv<-'.csv'
    
    #titulo1<-paste( corrpos ,  csv , sep="")
    
    titulo0<-paste( corrpre , csv , sep="")
    tituloR<-paste( rresposta , csv , sep="")
    tituloX<-paste( sucep , csv , sep="")
    tituloT<-paste( Xteff , csv , sep="")
    tituloApre<-paste( Assimpre , csv , sep="")
    tituloApos<-paste( Assimpos , csv , sep="")
    
    ns<- qero$ns
    kick<- qero$kick
    
    tkick<- qero$tkick/(qero$dt*qero$nc)
    print(tkick)
    
    respostaa<- numeric(Dt)
    suceptibilidade <- numeric(Dt)
    XTeff<- numeric(Dt)
    
    print('1foi até aqui')
    corpre<-      (  funcaodecorrelacaoSample(qero$posicao, (tkick - Dt - 1) , Dt , ns)   )
    print('passou o corrsample')
    
    for (i in 1:Dt)
    {
        mediapos <-  mean(qero$posicao[tkick+i, ]  )
        respostaa[i] <- mediapos/kick
        #suceptibilidade[i] <- sum(respostaa[1:i])
        suceptibilidade[i] <- sum(respostaa[1:i]) * qero$dt
    }
    
    iii <-  iii + 1
    
    print('2foi até aqui')
    #assim<- Assimetria(qero, 100, 500)
    
    #NOVIDADE AQUI ASSIMETRIA
    
    print(tkick-Dt-1)
    print(tkick+Dt)
    
    assimfdtpre<- AssimetriaFDT(qero, tkick-Dt-1, Dt)
    assimfdtpos<- AssimetriaFDT(qero, tkick, Dt)

    print('chegao no save')
    
    SalvaCorrelacao(assimfdtpre, tituloApre)
    SalvaCorrelacao(assimfdtpos, tituloApos)
    
    SalvaCorrelacao(corpre, titulo0)
    SalvaCorrelacao(suceptibilidade, tituloX)
    
    #SalvaCorrelacao(corpos, titulo1)
    
    SalvaCorrelacao(respostaa, tituloR)
    
    #dcorranalitico <-read.table(  paste('TDcorrv', toString(kk) ,'.csv', sep='' )    )
    
    dcorranalitico <-read.table(  paste('TDcorr', toString(kk) ,'.csv', sep='' )    )
    dcorranalitico <- dcorranalitico [[1]]
    
    #dcorrlange <-read.table(  'dcorlangev.csv'   )
    
    dcorrlange <-read.table(  'dcorlange.csv'   )
    dcorrlange <- dcorrlange [[1]]
    lala<-numeric(length(dcorrlange))
    for(i in 1:length(dcorrlange)){lala[i]<- as.numeric(strsplit(toString(dcorrlange[i]), '\\*')[[1]][1])}
    dcorrlange<-lala
    
    corranalitico <-read.table(  paste('corr', toString(kk) ,'.csv', sep='' )    )
    corranalitico <- corranalitico [[1]]
    
    #print(corranalitico)
    #print(dcorranalitico)
    #print(dcorrlange)
    
    T[iii]<-  - min(dcorranalitico)/max(respostaa)
    
    Tlan[iii]<-   0.546293  / max(respostaa)  # parametros mecanicos iguais a um
    
    ####
    
    for(j in 1:Dt)
    {
        XTeff[j]<-(dcorranalitico[j]+min(dcorranalitico))/(respostaa[j]+min(respostaa))
    }
    
    ####
    
    
    SalvaCorrelacao(XTeff, tituloT)
    
    respostaam <- respostaa/max(respostaa)
    dcorranalitico <-  dcorranalitico/ (- min(dcorranalitico) )
    dcorrlange <-  dcorrlange/ (- min(dcorrlange) )
    
    
    
    print('chegao no save 2')
    print(length(respostaam))
    print(length(tempo))
    print(length(dcorranalitico))
    
    #####
    grafs<-paste(letra , "TRespostaDcorr", toString(vez) ,".jpg", sep="")
    grafs<-toString(grafs)
    jpeg(grafs)
    
    plot( respostaam ~ tempo , col = 'red' )
    lines(   - dcorranalitico ~ tempo, col = 'blue' )
    lines(   - dcorrlange ~ tempo, col = 'green' )
    
    dev.off()
    ######
    
    print('chegao no save 3')
    
    
    grafs<-paste(letra ,"XCTkurchan", toString(vez) ,".jpg", sep="")
    grafs<-toString(grafs)
    jpeg(grafs)
    

    plot(  suceptibilidade   ~   corranalitico , col='red' )
    dev.off()
    
    
    #####
    
    grafxteff<-paste(letra , "XTeff", toString(vez) ,".jpg", sep="")
    grafxteff<-toString(grafxteff)
    jpeg(grafxteff)
    
    
    plot( XTeff ~ tempo , col = 'red' )
    
    
    dev.off()
    
    ######
    
    grafassimpre<-paste(letra , "Assimpre", toString(vez) ,".jpg", sep="")
    grafassim<-toString(grafassimpre)
    jpeg(grafassimpre)
    
    plot( assimfdtpre ~ tempo , col = 'red' )
    print('chegao no save 4')
    
    dev.off()
    
    grafassimpos<-paste(letra , "Assimpos", toString(vez) ,".jpg", sep="")
    grafassimpos<-toString(grafassimpos)
    jpeg(grafassimpos)
    
    plot( assimfdtpos ~ tempo , col = 'red' )
    print('chegao no save 5')
    
    dev.off()
    
    ######
    
    
    rsquared[iii] <- cor(dcorranalitico, respostaa  )**2
    rsquaredlan[iii] <- cor(dcorrlange, respostaa  )**2
    alppha[iii] <-  alpha[kk]
    
    }
    
    SalvaCorrelacao(T,       paste(letra ,'temperaturaeff.csv', sep="" )  )
    SalvaCorrelacao(Tlan,       paste(letra ,'temperaturalaneff.csv', sep="" )  )
    SalvaCorrelacao(alppha,       paste(letra ,'alpha.csv' , sep="")  )
    SalvaCorrelacao(rsquaredlan, paste(letra ,'rsquaredlan.csv' , sep="")  )
    SalvaCorrelacao(rsquared, paste(letra ,'rsquared.csv' , sep="")  )
    
    #####
    grafs<-paste(letra , "errobimolan.jpg", sep="")
    grafs<-toString(grafs)
    jpeg(grafs)
    
    plot(  rsquared   ~   alppha , col='red' )
    points(  rsquaredlan   ~  alppha, col='blue' )
    dev.off()
    #####
    
    grafs<-paste(letra , "temperaturaeff.jpg", sep="")
    grafs<-toString(grafs)
    jpeg(grafs)
    
    plot(  T   ~   alppha , col='red' )
    points(  Tlan   ~  alppha, col='blue' )
    dev.off()
    
    
}

TesteFuncaoResposta2xacte <- function(paralel, nsamples){
    
   
    alpha <- seq(1.0,30.0,1.0)
    
    
    sequencia <- seq(1,30,1)
    
    T <- numeric(  length(sequencia)   )
    Tlan <-     numeric(  length(sequencia)   )
    
    alppha <- numeric(  length(sequencia)   )
    
    xa<- 1
    
    rsquared <- numeric(  length(sequencia)   )
    rsquaredlan <- numeric(  length(sequencia)   )
    
    iii <- 0
    
    
    for (kk in sequencia)
    {
        print(alpha[kk])
        print(kk)
        
        # Temp = 0.710526   /. {k -> 1, \[Gamma] -> 1, m -> 1, a -> 1.5, \[Alpha] -> 1.5}   que funcionou
        #vamo tentar manter a temperatura
        
        #Temp <- 0.710526
        
        #Temp <- 1
        
        
        #xa<- sqrt( Temp*(gamma*(ctelastica+alpha[kk]*(gamma+massa*alpha[kk])))/(massa*alpha[kk])   )   #-> será que isso garante o encaixe???
        
        #xa<- sqrt( Temp*(1*(1+alpha[kk]*(1+1*alpha[kk])))/(1*alpha[kk])   )
        
        
        
        #qero<-SimulacaoBimodalKickX(nsamples,2510,2000,0.000005,1,1,1, alpha[kk] ,0, xa, 0, 0.5,  11 ) # aqui mudou o alfa em relacao ao anterior
        
        qero<-SimulacaoBimodalKickXsalvaEta(nsamples,2510,10,0.001,1,1,1, alpha[kk] ,0, xa, 0, 0.5,  11 )
        
        #qero<-SimulacaoBimodalKickV(nsamples,2510,10,0.001,1,1,1, alpha[kk] ,0, xa, 0, 0.5,  11 )
        
        Dt<- 1000
        
        #tempo<-  seq(0.002, 4, 0.002)
        
        tempo <- seq(0.01, 10, 0.01)
        
        vez<- kk
        
        letra <- paralel
        
        #corrpos <- paste(toString(vez) , 'Bcorrpos' , sep = ""  )
        
        corrpre <- paste( letra ,  'corrpre', toString(vez)  , sep="")
        rresposta <- paste( letra, 'resposta', toString(vez) , sep="")
        
        sucep <- paste( letra, 'sucep', toString(vez) , sep="")
        Xteff <- paste( letra, 'XTeff', toString(vez) , sep="")
        
        Assimpre <- paste( letra, 'Assimpre', toString(vez) , sep="")
        Assimpos <- paste( letra, 'Assimpos', toString(vez) , sep="")
        
        csv<-'.csv'
        
        #titulo1<-paste( corrpos ,  csv , sep="")
        
        titulo0<-paste( corrpre , csv , sep="")
        tituloR<-paste( rresposta , csv , sep="")
        tituloX<-paste( sucep , csv , sep="")
        tituloT<-paste( Xteff , csv , sep="")
        tituloApre<-paste( Assimpre , csv , sep="")
        tituloApos<-paste( Assimpos , csv , sep="")
        
        ns<- qero$ns
        kick<- qero$kick
        
        tkick<- qero$tkick/(qero$dt*qero$nc)
        print(tkick)
        
        respostaa<- numeric(Dt)
        suceptibilidade <- numeric(Dt)
        XTeff<- numeric(Dt)
        
        print('1foi até aqui')
        corpre<-      (  funcaodecorrelacaoSample(qero$posicao, (tkick - Dt - 1) , Dt , ns)   )
        print('passou o corrsample')
        
        for (i in 1:Dt)
        {
            mediapos <-  mean(qero$posicao[tkick+i, ]  )
            respostaa[i] <- mediapos/kick
            #suceptibilidade[i] <- sum(respostaa[1:i])
            suceptibilidade[i] <- sum(respostaa[1:i]) * qero$dt
        }
        
        iii <-  iii + 1
        
        print('2foi até aqui')
        #assim<- Assimetria(qero, 100, 500)
        
        #NOVIDADE AQUI ASSIMETRIA
        
        print(tkick-Dt-1)
        print(tkick+Dt)
        
        assimfdtpre<- AssimetriaFDT(qero, tkick-Dt-1, Dt)
        assimfdtpos<- AssimetriaFDT(qero, tkick, Dt)
        
        print('chegao no save')
        
        SalvaCorrelacao(assimfdtpre, tituloApre)
        SalvaCorrelacao(assimfdtpos, tituloApos)
        
        SalvaCorrelacao(corpre, titulo0)
        SalvaCorrelacao(suceptibilidade, tituloX)
        
        #SalvaCorrelacao(corpos, titulo1)
        
        SalvaCorrelacao(respostaa, tituloR)
        
        #dcorranalitico <-read.table(  paste('TDcorrv', toString(kk) ,'.csv', sep='' )    )
        
        dcorranalitico <-read.table(  paste('TDcorr', toString(kk) ,'.csv', sep='' )    )
        dcorranalitico <- dcorranalitico [[1]]
        
        #dcorrlange <-read.table(  'dcorlangev.csv'   )
        
        dcorrlange <-read.table(  'dcorlange.csv'   )
        dcorrlange <- dcorrlange [[1]]
        lala<-numeric(length(dcorrlange))
        for(i in 1:length(dcorrlange)){lala[i]<- as.numeric(strsplit(toString(dcorrlange[i]), '\\*')[[1]][1])}
        dcorrlange<-lala
        
        corranalitico <-read.table(  paste('corr', toString(kk) ,'.csv', sep='' )    )
        corranalitico <- corranalitico [[1]]
        
        #print(corranalitico)
        #print(dcorranalitico)
        #print(dcorrlange)
        
        T[iii]<-  - min(dcorranalitico)/max(respostaa)
        
        Tlan[iii]<-   0.546293  / max(respostaa)  # parametros mecanicos iguais a um
        
        ####
        
        for(j in 1:Dt)
        {
            XTeff[j]<-(dcorranalitico[j]+min(dcorranalitico))/(respostaa[j]+min(respostaa))
        }
        
        ####
        
        
        SalvaCorrelacao(XTeff, tituloT)
        
        respostaam <- respostaa/max(respostaa)
        dcorranalitico <-  dcorranalitico/ (- min(dcorranalitico) )
        dcorrlange <-  dcorrlange/ (- min(dcorrlange) )
        
        
        
        print('chegao no save 2')
        print(length(respostaam))
        print(length(tempo))
        print(length(dcorranalitico))
        
        #####
        grafs<-paste(letra , "TRespostaDcorr", toString(vez) ,".jpg", sep="")
        grafs<-toString(grafs)
        jpeg(grafs)
        
        plot( respostaam ~ tempo , col = 'red' )
        lines(   - dcorranalitico ~ tempo, col = 'blue' )
        lines(   - dcorrlange ~ tempo, col = 'green' )
        
        dev.off()
        ######
        
        print('chegao no save 3')
        
        
        grafs<-paste(letra ,"XCTkurchan", toString(vez) ,".jpg", sep="")
        grafs<-toString(grafs)
        jpeg(grafs)
        
        
        plot(  suceptibilidade   ~   corranalitico , col='red' )
        dev.off()
        
        
        #####
        
        grafxteff<-paste(letra , "XTeff", toString(vez) ,".jpg", sep="")
        grafxteff<-toString(grafxteff)
        jpeg(grafxteff)
        
        
        plot( XTeff ~ tempo , col = 'red' )
        
        
        dev.off()
        
        ######
        
        grafassimpre<-paste(letra , "Assimpre", toString(vez) ,".jpg", sep="")
        grafassim<-toString(grafassimpre)
        jpeg(grafassimpre)
        
        plot( assimfdtpre ~ tempo , col = 'red' )
        print('chegao no save 4')
        
        dev.off()
        
        grafassimpos<-paste(letra , "Assimpos", toString(vez) ,".jpg", sep="")
        grafassimpos<-toString(grafassimpos)
        jpeg(grafassimpos)
        
        plot( assimfdtpos ~ tempo , col = 'red' )
        print('chegao no save 5')
        
        dev.off()
        
        ######
        
        
        rsquared[iii] <- cor(dcorranalitico, respostaa  )**2
        rsquaredlan[iii] <- cor(dcorrlange, respostaa  )**2
        alppha[iii] <-  alpha[kk]
        
    }
    
    SalvaCorrelacao(T,       paste(letra ,'temperaturaeff.csv', sep="" )  )
    SalvaCorrelacao(Tlan,       paste(letra ,'temperaturalaneff.csv', sep="" )  )
    SalvaCorrelacao(alppha,       paste(letra ,'alpha.csv' , sep="")  )
    SalvaCorrelacao(rsquaredlan, paste(letra ,'rsquaredlan.csv' , sep="")  )
    SalvaCorrelacao(rsquared, paste(letra ,'rsquared.csv' , sep="")  )
    
    #####
    grafs<-paste(letra , "errobimolan.jpg", sep="")
    grafs<-toString(grafs)
    jpeg(grafs)
    
    plot(  rsquared   ~   alppha , col='red' )
    points(  rsquaredlan   ~  alppha, col='blue' )
    dev.off()
    #####
    
    grafs<-paste(letra , "temperaturaeff.jpg", sep="")
    grafs<-toString(grafs)
    jpeg(grafs)
    
    plot(  T   ~   alppha , col='red' )
    points(  Tlan   ~  alppha, col='blue' )
    dev.off()
    
    
}




TesteFuncaoResposta3 <- function(paralel, nsamples){
    

    alpha <- seq(0.1,6 ,0.1)

    sequencia <- seq(1,60,1)
    
    T <- numeric(  length(sequencia)   )
    Tlan <-     numeric(  length(sequencia)   )
    alppha <- numeric(  length(sequencia)   )
    
    
    rsquared <- numeric(  length(sequencia)   )
    rsquaredlan <- numeric(  length(sequencia)   )
    
    iii <- 0
    
    
    for (kk in sequencia)
    {
        print(alpha[kk])
        print(kk)
        
        Temp <- 1
        
        xa<- sqrt( Temp*(1*(1+alpha[kk]*(1+1*alpha[kk])))/(1*alpha[kk])   )
        
        #function(ns, nt, nc, dt, massa, gamma, ctelastica, alfa, dalfa, xa, dxa, deltax, tt){
        
        qero<-SimulacaoBimodalKickX(nsamples,2510,10,0.001,1,1,0.3, alpha[kk] ,0, xa, 0, 0.5,  11 )
        
        Dt<- 1000
        
        tempo <- seq(0.01, 10, 0.01)
        
        vez<- kk
        
        letra <- paralel
        
        #corrpos <- paste(toString(vez) , 'Bcorrpos' , sep = ""  )
        
        corrpre <- paste( letra ,  'corrpre', toString(vez)  , sep="")
        rresposta <- paste( letra, 'resposta', toString(vez) , sep="")
        sucep <- paste( letra, 'sucep', toString(vez) , sep="")
        Xteff <- paste( letra, 'XTeff', toString(vez) , sep="")
        Assim <- paste( letra, 'Assim', toString(vez) , sep="")
        
        csv<-'.csv'
        
        #titulo1<-paste( corrpos ,  csv , sep="")
        
        titulo0<-paste( corrpre , csv , sep="")
        tituloR<-paste( rresposta , csv , sep="")
        tituloX<-paste( sucep , csv , sep="")
        tituloT<-paste( Xteff , csv , sep="")
        tituloA<-paste( Assim , csv , sep="")
        
        ns<- qero$ns
        kick<- qero$kick
        
        tkick<- qero$tkick/(qero$dt*qero$nc)
        print(tkick)
        
        respostaa<- numeric(Dt)
        suceptibilidade <- numeric(Dt)
        XTeff<- numeric(Dt)
        
        print('1foi até aqui')
        corpre<-      (  funcaodecorrelacaoSample(qero$posicao, (tkick - Dt - 1) , Dt , ns)   )
        print('passou o corrsample')
        
        for (i in 1:Dt)
        {
            mediapos <-  mean(qero$posicao[tkick+i, ]  )
            respostaa[i] <- mediapos/kick
            #suceptibilidade[i] <- sum(respostaa[1:i])
            suceptibilidade[i] <- sum(respostaa[1:i]) * qero$dt
        }
        
        iii <-  iii + 1
        
        
        assim<- Assimetria(qero, 100, 500)
        
        
        print('chegao no save')
        
        SalvaCorrelacao(assim, tituloA)
        SalvaCorrelacao(corpre, titulo0)
        SalvaCorrelacao(suceptibilidade, tituloX)
        #SalvaCorrelacao(corpos, titulo1)
        SalvaCorrelacao(respostaa, tituloR)
        
        #dcorranalitico <-read.table(  paste('TDcorrv', toString(kk) ,'.csv', sep='' )    )
        dcorranalitico <-read.table(  paste('TDcorr', toString(kk) ,'.csv', sep='' )    )
        dcorranalitico <- dcorranalitico [[1]]
        
        #dcorrlange <-read.table(  'dcorlangev.csv'   )
        dcorrlange <-read.table(  'dcorlange.csv'   )
        dcorrlange <- dcorrlange [[1]]
        lala<-numeric(length(dcorrlange))
        for(i in 1:length(dcorrlange)){lala[i]<- as.numeric(strsplit(toString(dcorrlange[i]), '\\*')[[1]][1])}
        dcorrlange<-lala
        
        corranalitico <-read.table(  paste('corr', toString(kk) ,'.csv', sep='' )    )
        corranalitico <- corranalitico [[1]]
        
        #print(corranalitico)
        #print(dcorranalitico)
        #print(dcorrlange)
        
        T[iii]<-  - min(dcorranalitico)/max(respostaa)
        
        Tlan[iii]<-   0.546293  / max(respostaa)  # parametros mecanicos iguais a um
        
        ####
        
        for(j in 1:Dt)
        {
            XTeff[j]<-(dcorranalitico[j]+min(dcorranalitico))/(respostaa[j]+min(respostaa))
        }
        ####
        
        
        SalvaCorrelacao(XTeff, tituloT)
        
        respostaam <- respostaa/max(respostaa)
        dcorranalitico <-  dcorranalitico/ (- min(dcorranalitico) )
        dcorrlange <-  dcorrlange/ (- min(dcorrlange) )
        
        
        
        print('chegao no save 2')
        print(length(respostaam))
        print(length(tempo))
        print(length(dcorranalitico))
        
        #####
        grafs<-paste(letra , "TRespostaDcorr", toString(vez) ,".jpg", sep="")
        grafs<-toString(grafs)
        jpeg(grafs)
        
        plot( respostaam ~ tempo , col = 'red' )
        lines(   - dcorranalitico ~ tempo, col = 'blue' )
        lines(   - dcorrlange ~ tempo, col = 'green' )
        
        dev.off()
        ######
        
        print('chegao no save 3')
        
        
        grafs<-paste(letra ,"XCTkurchan", toString(vez) ,".jpg", sep="")
        grafs<-toString(grafs)
        jpeg(grafs)
        
        
        plot(  suceptibilidade   ~   corranalitico , col='red' )
        dev.off()
        
        
        #####
        
        grafxteff<-paste(letra , "XTeff", toString(vez) ,".jpg", sep="")
        grafxteff<-toString(grafxteff)
        jpeg(grafxteff)
        
        plot( XTeff ~ tempo , col = 'red' )
        
        
        dev.off()
        
        ######
        
        grafassim<-paste(letra , "Assim", toString(vez) ,".jpg", sep="")
        grafassim<-toString(grafassim)
        jpeg(grafassim)
        
        plot( assim ~ tempo , col = 'red' )
        
        
        dev.off()
        
        ######
        
        
        rsquared[iii] <- cor(dcorranalitico, respostaa  )**2
        rsquaredlan[iii] <- cor(dcorrlange, respostaa  )**2
        alppha[iii] <-  alpha[kk]
        
    }
    
    SalvaCorrelacao(T,       paste(letra ,'temperaturaeff.csv', sep="" )  )
    SalvaCorrelacao(Tlan,       paste(letra ,'temperaturalaneff.csv', sep="" )  )
    SalvaCorrelacao(alppha,       paste(letra ,'alpha.csv' , sep="")  )
    SalvaCorrelacao(rsquaredlan, paste(letra ,'rsquaredlan.csv' , sep="")  )
    SalvaCorrelacao(rsquared, paste(letra ,'rsquared.csv' , sep="")  )
    
    #####
    grafs<-paste(letra , "errobimolan.jpg", sep="")
    grafs<-toString(grafs)
    jpeg(grafs)
    
    plot(  rsquared   ~   alppha , col='red' )
    points(  rsquaredlan   ~  alppha, col='blue' )
    dev.off()
    #####
    
    grafs<-paste(letra , "temperaturaeff.jpg", sep="")
    grafs<-toString(grafs)
    jpeg(grafs)
    
    plot(  T   ~   alppha , col='red' )
    points(  Tlan   ~  alppha, col='blue' )
    dev.off()
    
    
}

AssimetriaFDT<- function(qero, t, lag){
 
 
    ns<- qero$ns
    
    asimean<-numeric(lag)
    
    for(j in 1:lag)
    {
        for(i in 1:ns)
        {
            #asimean[j] <- asimean[j] + (qero$k*qero$posicao[t+j,i]  + qero$gamma*qero$velocidade[t+j,i])*qero$posicao[t,i] - (qero$k*qero$posicao[t,i] + qero$gamma*qero$velocidade[t,i])*qero$posicao[t+j,i]
            
            asimean[j] <- asimean[j] + ( qero$ruido[t+j,i]*qero$posicao[t,i] ) - ( qero$ruido[t,i]*qero$posicao[t+j,i] )
        }
    }
    
    asimean<-asimean/ns
    return(asimean)
    
}

Assimetria <- function(qero, t, lag){
    
    
    ns<- qero$ns
    
    asimean<-numeric(lag)
    
    for(j in 1:lag)
    {
        for(i in 1:ns)
        {
            # asimean[j] <- asimean[j] + (qero$k*qero$posicao[t+j,i]  + qero$gamma*qero$velocidade[t+j,i])*qero$posicao[t,i] - (qero$k*qero$posicao[t,i] + qero$gamma*qero$velocidade[t,i])*qero$posicao[t+j,i]
            
            asimean[j] <- asimean[j] + ( qero$posicao[t+j,i]*qero$posicao[t,i] ) - ( qero$posicao[t,i]*qero$posicao[t+j,i] )
        }

        
    }
    
    asimean<-asimean/ns
    return(asimean)
    
}

Assimetriatot <- function(qero, t, lag){
    
    
    ns<- qero$ns
    
    asimean<-numeric(lag)
    
    for(j in 1:lag)
    {
        for(i in 1:ns)
        {
            asimean[j] <- asimean[j] + (qero$k*qero$posicao[t+j,i]  + qero$gamma*qero$velocidade[t+j,i])*qero$posicao[t,i] - (qero$k*qero$posicao[t,i] + qero$gamma*qero$velocidade[t,i])*qero$posicao[t+j,i]
        }
    }
    
    asimean<-asimean/ns
    return(asimean)
    
}


TesteFuncaoRespostaHard <- function(paralel, nsamples){
    
  
    alpha <- 10
    #sequencia <- seq(1,60,1)
    
  
    
    iii <- 0
    
    

    

    Temp <- 1

    xa<- sqrt( Temp*(1*(1+alpha*(1+1*alpha)))/(1*alpha)   )

    qero<-SimulacaoBimodalKickX(nsamples,25100,10,0.0001,1,1,1, alpha ,0, xa, 0, 0.5,  11 )

    Dt<- 10000

    tempo <- seq(0.001, 10, 0.001)
    
    vez <- 123
    
    letra <- paralel
    
    #corrpos <- paste(toString(vez) , 'Bcorrpos' , sep = ""  )
    
    corrpre <- paste( letra ,  'corrpre', toString(vez)  , sep="")
    rresposta <- paste( letra, 'resposta', toString(vez) , sep="")
    sucep <- paste( letra, 'sucep', toString(vez) , sep="")
    
    csv<-'.csv'
    
    #titulo1<-paste( corrpos ,  csv , sep="")
    
    titulo0<-paste( corrpre , csv , sep="")
    tituloR<-paste( rresposta , csv , sep="")
    tituloX<-paste( sucep , csv , sep="")
    
    ns<- qero$ns
    kick<- qero$kick
    
    tkick<- qero$tkick/(qero$dt*qero$nc)
    print(tkick)
    
    respostaa<- numeric(Dt)
    suceptibilidade <- numeric(Dt)
    
    print('1foi até aqui')
    #corpre<-      (  funcaodecorrelacaoSample(qero$posicao, (tkick - Dt - 1) , Dt , ns)   )
    print('passou o corrsample')
    
    for (i in 1:Dt)
    {
        mediapos <-  mean(qero$posicao[tkick+i, ]  )
        respostaa[i] <- mediapos/kick
        suceptibilidade[i] <- sum(respostaa[1:i])
    }
    
    iii <-  iii + 1
    
    
    print('chegao no save')
    
    #  SalvaCorrelacao(corpre, titulo0)
    SalvaCorrelacao(suceptibilidade, tituloX)
    #SalvaCorrelacao(corpos, titulo1)
    SalvaCorrelacao(respostaa, tituloR)
    
    #dcorranalitico <-read.table(  paste('TDcorrv', toString(kk) ,'.csv', sep='' )    )
    #dcorranalitico <-read.table(  paste('TDcorr', toString(50) ,'.csv', sep='' )    )
    #dcorranalitico <- dcorranalitico [[1]]
    
    #dcorrlange <-read.table(  'dcorlangev.csv'   )
    #dcorrlange <-read.table(  'dcorlange.csv'   )
    #dcorrlange <- dcorrlange [[1]]
    #lala<-numeric(length(dcorrlange))
    #for(i in 1:length(dcorrlange)){lala[i]<- as.numeric(strsplit(toString(dcorrlange[i]), '\\*')[[1]][1])}
    #dcorrlange<-lala
    
    #corranalitico <-read.table(  paste('corr', toString(50) ,'.csv', sep='' )    )
    #corranalitico <- corranalitico [[1]]
    
    #print(corranalitico)
    #print(dcorranalitico)
    #print(dcorrlange)
    
    #T<-  - min(dcorranalitico)/max(respostaa)
    #Tlan<-   0.546293  / max(respostaa)  # parametros mecanicos iguais a um
    
    
    respostaam <- respostaa/max(respostaa)
    # dcorranalitico <-  dcorranalitico/ (- min(dcorranalitico) )
    #dcorrlange <-  dcorrlange/ (- min(dcorrlange) )
    
    
    
    print('chegao no save 2')
    print(length(respostaam))
    print(length(tempo))
    #print(length(dcorranalitico))
    
    #####
    grafs<-paste(letra , "TRespostaDcorr", toString(vez) ,".jpg", sep="")
    grafs<-toString(grafs)
    jpeg(grafs)
    
    plot( respostaam ~ tempo , col = 'red' )
    # lines(   - dcorranalitico ~ tempo, col = 'blue' )
    #lines(   - dcorrlange ~ tempo, col = 'green' )
    
    dev.off()
    ######
    
    print('chegao no save 3')
    
    
    #grafs<-paste(letra ,"XCTkurchan", toString(vez) ,".jpg", sep="")
    #grafs<-toString(grafs)
    #jpeg(grafs)
    
    
    #plot(  suceptibilidade   ~   corranalitico , col='red' )
    #dev.off()
    
    
    
    #rsquared <- cor(dcorranalitico, respostaa  )**2
    #rsquaredlan <- cor(dcorrlange, respostaa  )**2
    

    
    # SalvaCorrelacao(T,       paste(letra ,'temperaturaeff.csv', sep="" )  )
    # SalvaCorrelacao(Tlan,       paste(letra ,'temperaturalaneff.csv', sep="" )  )
    #SalvaCorrelacao(alppha,       paste(letra ,'alpha.csv' , sep="")  )
    #SalvaCorrelacao(rsquaredlan, paste(letra ,'rsquaredlan.csv' , sep="")  )
    #SalvaCorrelacao(rsquared, paste(letra ,'rsquared.csv' , sep="")  )
    
    #####
    #  grafs<-paste(letra , "errobimolan.jpg", sep="")
    # grafs<-toString(grafs)
    #jpeg(grafs)
    
    # plot(  rsquared   ~   alppha , col='red' )
    #points(  rsquaredlan   ~  alppha, col='blue' )
    #dev.off()
    #####
    
    #grafs<-paste(letra , "temperaturaeff.jpg", sep="")
    #grafs<-toString(grafs)
    #jpeg(grafs)
    
    # plot(  T   ~   alppha , col='red' )
    # points(  Tlan   ~  alppha, col='blue' )
    # dev.off()
    
    
}



TesteFuncaoResposta4 <- function(){
    
    #source('Entropiaultima.R')
    
    #rresposta <- 'trueresposta'
    #rresposta <- '2trueresposta'
    #rresposta <- '3trueresposta'
    
    
   
    
    nsamples <- 100
    
    alpha <- seq(0.1,6,0.1)
    

    
    sequencia <- seq(5,60,10)
    T <- numeric(  length(sequencia)   )
    alppha <- numeric(  length(sequencia)   )
    rsquared <- numeric(  length(sequencia)   )
    
    iii <- 0
    
    
    for (kk in sequencia)
    {
        
    print(kk)
    
    qero<-SimulacaoBimodalKickX(nsamples,2510,1000,0.00001,1,1,1, alpha[kk] ,0,1.5,0, 0.5,  11 ) # aqui mudou o alfa em relacao ao anterior
    
    
    Dt<- 1000
    
    #tempo<-  seq(0.002, 4, 0.002)
    
    tempo <- seq(0.01, 10, 0.01)
    
    vez<- kk
    
    letra <- 'D'
    
    #corrpos <- paste(toString(vez) , 'Bcorrpos' , sep = ""  )
    
    corrpre <- paste(toString(vez) , letra  , 'corrpre', sep="")
    rresposta <- paste(toString(vez), letra ,'resposta', sep="")
    sucep <- paste(toString(vez), letra ,'sucep', sep="")
    
    csv<-'.csv'
    
    #titulo1<-paste( corrpos ,  csv , sep="")
    
    titulo0<-paste( corrpre , csv , sep="")

    tituloR<-paste( rresposta , csv , sep="")
    tituloX<-paste( sucep , csv , sep="")
    
    ns<- qero$ns
    kick<- qero$kick
    
    tkick<- qero$tkick/(qero$dt*qero$nc)
    print(tkick)
    
    respostaa<- numeric(Dt)
    
    suceptibilidade <- numeric(Dt)
    
    #mediapre<-   mean( qero$posicao[1:tkick,]  )
    print('1foi até aqui')
    corpre<-      (  funcaodecorrelacaoSample(qero$posicao, (tkick - Dt - 1) , Dt , ns)   )
    #corpos<-      (  funcaodecorrelacaoSample(qero$posicao, tkick , Dt , ns)       )
    print('passou o corrsample')
    
    for (i in 1:Dt)
    {
        mediapos <-  mean(qero$posicao[tkick+i, ]  )
        respostaa[i] <- mediapos/kick # - mediapre -> mediapre é igual a zero
        suceptibilidade[i] <- sum(respostaa[1:i])
    }
    
    iii <-  iii + 1
    
    #fim<-fim/kick
    
    #for (i in 1:Dt)
    #{
    #    mediapos<-  mean(qero$posicao[tkick+i, ]  )
    #    fim[i]<- mediapos # - mediapre -> mediapre é igual a zero
    #}
    print('chegao no save')
    
    SalvaCorrelacao(corpre, titulo0)
    SalvaCorrelacao(suceptibilidade, tituloX)
    #SalvaCorrelacao(corpos, titulo1)
    SalvaCorrelacao(respostaa, tituloR)
    
    dcorranalitico <-read.table(  paste('Dcorr', toString(kk) ,'.csv', sep='' )    )
    dcorranalitico <- dcorranalitico [[1]]
    
    T[iii]<-  - min(dcorranalitico)/max(respostaa)
    
    respostaam <- respostaa/max(respostaa)
    dcorranalitico <-  dcorranalitico/ (- min(dcorranalitico) )
    
    print('chegao no save 2')
    print(length(respostaam))
    print(length(tempo))
    print(length(dcorranalitico))
    
    
    
    #####
    grafs<-paste("RespostaDcorr", toString(vez) ,".jpg", sep="")
    grafs<-toString(grafs)
    jpeg(grafs)
    
    plot(  respostaam   ~   tempo , col='red' )
    lines(   - dcorranalitico   ~  tempo, col='blue' )
    dev.off()
    ######
    
    print('chegao no save 3')
    
    
    grafs<-paste("XCkurchan", toString(vez) ,".jpg", sep="")
    grafs<-toString(grafs)
    jpeg(grafs)
    
    
    plot(  suceptibilidade   ~   corpre , col='red' )
    #lines(   dcorranalitico   ~  tempo, col='blue' )
    dev.off()
    
    #grafs<-paste("RespostaDcorr", ToString(vez) ,".jpg", sep="")
    #grafs<-toString(grafs)
    #jpeg(grafs)
    #plot(  fimm   ~   tempo , col='red' )
    #lines(   dcorranalitico   ~  tempo, col='blue' )
    #dev.off()
    
    
    
    rsquared[iii] <- cor(dcorranalitico, respostaa  )**2
    
    alppha[iii] <-  alpha[kk]
    
    }
    
    SalvaCorrelacao(T, 'temperaturaeff.csv')
    SalvaCorrelacao(alppha, 'alpha.csv')
    SalvaCorrelacao(rsquared, 'rsquared.csv')
    
    
    #grafs<-paste("Rsquared", ToString(vez) ,".jpg", sep="")
    #grafs<-toString(grafs)
    #jpeg(grafs)
    
    #plot(  fimm   ~   tempo , col='red' )
    #lines(   dcorranalitico   ~  tempo, col='blue' )
    #dev.off()
    
    
}



Testinlan<- function(){
    
   
    
    
    sequencia <- seq(5,60,10)
    T <- numeric(  length(sequencia)   )
    alppha <- numeric(  length(sequencia)   )
    rsquared <- numeric(  length(sequencia)   )
    
    iii <- 0
    
    letra <- 'D'
    
    for (kk in sequencia)
    {
        vez<-kk
        
        rresposta <- paste(toString(vez), letra ,'resposta', sep="")
        
        csv<-'.csv'
        
        tituloR<-paste( rresposta , csv , sep="")

        
    
        respostaa <- read.table(  tituloR )
        respostaa <- respostaa [[1]]
        
        print(respostaa)
        
        dcorranalitico <-read.table(  paste('DDcorrLan', toString(kk) ,'.csv', sep='' )    )
        dcorranalitico <- dcorranalitico [[1]]
        
        
        #print(dcorranalitico)
        # T[iii]<-  - min(dcorranalitico)/max(respostaa)
        
        respostaam <- respostaa/max(respostaa)
        dcorranalitico <-  dcorranalitico/ (- min(dcorranalitico) )
        
        print('chegao no save 2')
        #  print(length(respostaam))
        #print(length(tempo))
        #print(length(dcorranalitico))
        
        tempo<-  seq(0.01, 10, 0.01)
        
        #####
        grafs<-paste("LRespostaDcorr", toString(vez) ,".jpg", sep="")
        grafs<-toString(grafs)
        jpeg(grafs)
        
        plot(  respostaam   ~   tempo , col='red' )
        lines(   - dcorranalitico   ~  tempo, col='blue' )
        dev.off()
        ######
        
        print('chegao no save 3')
        
    #  grafs<-paste("LXCkurchan", toString(vez) ,".jpg", sep="")
    #grafs<-toString(grafs)
    #jpeg(grafs)
    
    
    #plot(  suceptibilidade   ~   corpre , col='red' )
    #lines(   dcorranalitico   ~  tempo, col='blue' )
    # dev.off()
    
    #grafs<-paste("RespostaDcorr", ToString(vez) ,".jpg", sep="")
    #grafs<-toString(grafs)
    #jpeg(grafs)
    #plot(  fimm   ~   tempo , col='red' )
    #lines(   dcorranalitico   ~  tempo, col='blue' )
    #dev.off()
    
        rsquared[iii] <- cor(dcorranalitico, respostaa  )**2
        
        alppha[iii] <-  alpha[kk]
    
    }
    
}

TesteFuncaoResposta78 <- function(){
    
   
    nsamples <- 5000
    
    alpha <- seq(0.1,6,0.1)
    
    sequencia <- seq(5,60,10)
    
    T <- numeric(  length(sequencia)   )
    
    alppha <- numeric(  length(sequencia)   )
    
    rsquared <- numeric(  length(sequencia)   )
    
    iii <- 0
    
    
    for (kk in sequencia)
    {
        
        print(kk)
        
        tempo<-  seq(0.01,10,0.01)
        
        qero<-SimulacaoBimodalKickV(nsamples,3010,10,0.001,1,1,1, alpha[kk] ,0,1.5,0, 0.5,  20.1 )
        Dt<- 1000
        vez<- kk
        
        letra <- 'D'
        

        
        corrpre <- paste(toString(vez) , letra  , 'corrpre', sep="")
        rresposta <- paste(toString(vez), letra ,'resposta', sep="")
        sucep <- paste(toString(vez), letra ,'sucep', sep="")
        
        csv<-'.csv'
        
        titulo0<-paste( corrpre , csv , sep="")
        
        tituloR<-paste( rresposta , csv , sep="")
        tituloX<-paste( sucep , csv , sep="")
        
        ns<- qero$ns
        kick<- qero$kick
        tkick<- qero$tkick/(qero$dt*qero$nc)
        
        respostaa<- numeric(Dt)
        
        suceptibilidade <- numeric(Dt)
        
        corpre<-      (  funcaodecorrelacaoSample(qero$velocidade, (tkick - Dt - 1) , Dt , ns)   )
        
        
        for (i in 1:Dt)
        {
            mediapos <-  mean(qero$velocidade[tkick+i, ]  )
            respostaa[i] <- mediapos/kick # - mediapre -> mediapre é igual a zero
            suceptibilidade[i] <- sum(respostaa[1:i])
        }
        
        iii <-  iii + 1
        
        SalvaCorrelacao(corpre, titulo0)
        SalvaCorrelacao(suceptibilidade, tituloX)
        SalvaCorrelacao(respostaa, tituloR)
        
        dcorranalitico <-read.table(  paste('DDcorrv', toString(kk) ,'.csv', sep='' )    )
        dcorranalitico <- dcorranalitico [[1]]
        
        T[iii]<-  - min(dcorranalitico)/max(respostaa)
        
        respostaam <- respostaa/max(respostaa)
        dcorranalitico <-  dcorranalitico/ (- min(dcorranalitico) )
        
        
        #####
        grafs<-paste("VRespostaDcorr", toString(vez) ,".jpg", sep="")
        grafs<-toString(grafs)
        jpeg(grafs)
        
        plot(  respostaam   ~   tempo , col='red' )
        lines(   - dcorranalitico   ~  tempo, col='blue' )
        dev.off()
        ######
        
        grafs<-paste("XVCkurchan", toString(vez) ,".jpg", sep="")
        grafs<-toString(grafs)
        jpeg(grafs)
        
        
        plot(  suceptibilidade   ~   corpre , col='red' )
        #lines(   dcorranalitico   ~  tempo, col='blue' )
        dev.off()
        
        #grafs<-paste("RespostaDcorr", ToString(vez) ,".jpg", sep="")
        #grafs<-toString(grafs)
        #jpeg(grafs)
        #plot(  fimm   ~   tempo , col='red' )
        #lines(   dcorranalitico   ~  tempo, col='blue' )
        #dev.off()
        
        
        
        rsquared[iii] <- cor(dcorranalitico, respostaa  )**2
        
        alppha[iii] <-  alpha[kk]
        
    }
    
    SalvaCorrelacao(T, 'Vtemperaturaeff.csv')
    SalvaCorrelacao(alppha, 'Valpha.csv')
    SalvaCorrelacao(rsquared, 'Vrsquared.csv')
    
    
    #grafs<-paste("Rsquared", ToString(vez) ,".jpg", sep="")
    #grafs<-toString(grafs)
    #jpeg(grafs)
    
    #plot(  fimm   ~   tempo , col='red' )
    #lines(   dcorranalitico   ~  tempo, col='blue' )
    #dev.off()
    
    
}



TesteFuncaoRespostaLangevin <- function(){
    
    #source('Entropiaultima.R')
    
    rresposta <- 'truerespostaLange'
    
    csv<-'.csv'
    
    
    #SimulacaoLangevinKickF <- function(ns, nt, nc, dt, massa, gamma, ctelastica, temp, deltaf){

    
    
    qero<-SimulacaoLangevinKickF(10,10000,100,0.001,1,1,1,1,50000, 50 )
    
    tituloR<-paste( rresposta , csv , sep="")
    
    ns<- qero$ns
    
    kick<- qero$kick
    tkick<- qero$tkick /(qero$dt*qero$nc)
    
    Dt<- 100
    
    fim<- numeric(Dt)
    
    mediapre<-      mean( qero$posicao[1:tkick,]  )
    
    
    for (i in 1:Dt)
    {
        
        mediapos<-  mean(qero$posicao[tkick+i, ]  )
        
        fim[i]<- mediapos - mediapre
        
    }
    
    fim<-fim/kick
    
    
    SalvaCorrelacao(fim, tituloR)
    
}




TesteAgeFuncaoResposta <- function(){
    
    #source('Entropiaultima.R')
    
    corrpre <- '3corrpre'
    corrpos <- '3corrpos'
    rresposta <- '3rresposta'
    
    csv<-'.csv'
    
    for (i in 1:90)
    {
        
        titulo0<-paste( corrpre , toString(i) , csv , sep="")
        
        titulo1<-paste( corrpos , toString(i) , csv , sep="")
        
        tituloR<-paste( rresposta , toString(i) , csv , sep="")
        
        #SimulacaoBimodalKickF <- function(ns, nt, nc, dt, massa, gamma, ctelastica, alfa, dalfa, xa, dxa, deltax, tt)
        
        # PRMEIRA FOI 10000,  AGOR VOU FAZER 100000
        #1  qero<-SimulacaoBimodalKickF(5000,10000,100,0.001,1,1,1,1,0,1,0,10000, (10*i) + 50 )
        #2  qero<-SimulacaoBimodalKickF(5000,10000,100,0.001,1,1,1,1,0,1,0,100000, (10*i) + 50 )
        #3  qero<-SimulacaoBimodalKickX(5000,10000,100,0.001,1,1,1,1,0,1,0,10, (10*i) + 50 )
        #4
        
        
        qero<-SimulacaoBimodalKickF(5000,10000,100,0.001,1,1,1,1,0,1,0,50000, (10*i) + 50 )
        
        ns<- qero$ns
        
        kick<- qero$kick
        tkick<- qero$tkick /(qero$dt*qero$nc)
        
        Dt<- 100
        
        fim<- numeric(Dt)
        
        corpre<-      (  funcaodecorrelacaoSample(qero$posicao, (tkick-Dt -1) , Dt, ns)   )
        
        corpos<- ( funcaodecorrelacaoSample(qero$posicao, tkick , Dt ,ns))
        
        
        
        fim<- fim+ (corpos - corpre)
        
        fim<-fim/kick
        
        fim<- fim/ns
        
        
        SalvaCorrelacao(corpre, titulo0)
        
        SalvaCorrelacao(corpos, titulo1)
        
        SalvaCorrelacao(fim, tituloR)
        
    }
    
}



TesteEvolucaoEntropicaB <- function(){
    
    #SimulacaoBimodal <- function(ns, nt, nc, dt, massa, gamma, ctelastica, alfa, dalfa, xa, dxa){
    
    
    #qero1<-SimulacaoBimodal(100,1000,100,0.0001,1,1,1,1,0,1,0)   indice 2!
    
    
    #qero1<-SimulacaoBimodal(100,10000,100,0.00001,1,1,1,1,0,1,0)  # indice 3!
    
    qero1<-SimulacaoBimodal(50,5000,100,0.00001,1,1,1,1,0,1,0) # indice 4!
    
    wat1<-evolucaoentropica3D(qero1, 60)
    
    vez<- 4
    
    indice<- paste(  toString(vez) , '.csv' , sep="")
    
    
    
    SalvaCorrelacao(wat1$S, paste('evoentropia', indice , sep="")    )
    
    SalvaCorrelacao(wat1$Sx, paste('evoentropiax', indice , sep="")   )
    
    SalvaCorrelacao(wat1$Sv, paste('evoentropiav', indice , sep="")  )
    
    SalvaCorrelacao(wat1$Etempo, paste('evoenergia', indice , sep="")  )
    
    SalvaCorrelacao(wat1$Et, paste('evoenergiaM', indice , sep="")   )
    
    SalvaCorrelacao(wat1$jd, paste('evodissipada', indice , sep="")   )
    
    SalvaCorrelacao(wat1$ji, paste('evoinjetada', indice , sep="")   )
    
    SalvaCorrelacao(wat1$tempo, paste('evotempo', indice , sep="")   )
    

}

TesteEvolucaoEntropicaLan <- function(){
    
    #SimulacaoLangevin <- function(ns, nt, nc, dt, massa, gamma, ctelastica, temp){
    
    #qero1<-SimulacaoBimodal(100,1000,100,0.0001,1,1,1,1,0,1,0)   indice 2!
    
    #qero1<-SimulacaoLangevin(100,10000,100,0.00001,1,1,1,1)  # indice 3!
    
    #qero1<-SimulacaoLangevin(50,5000,100,0.00001,1,1,1,1)  # indice 4!
    
    qero1<-SimulacaoLangevin(50,5000,10,0.001,1,1,1,1) #indice 5 ainda n fiz
    
    wat1<-evolucaoentropica3DLan(qero1, 60)
    
    vez<- 5
    
    indice<- paste(  toString(vez) , '.csv' , sep="")
    
    
    
    SalvaCorrelacao(wat1$S, paste('Levoentropia', indice , sep="")    )
    
    SalvaCorrelacao(wat1$Sx, paste('Levoentropiax', indice , sep="")   )
    
    SalvaCorrelacao(wat1$Sv, paste('Levoentropiav', indice , sep="")  )
    
    SalvaCorrelacao(wat1$Etempo, paste('Levoenergia', indice , sep="")  )
    
    SalvaCorrelacao(wat1$Et, paste('LevoenergiaM', indice , sep="")   )
    
    SalvaCorrelacao(wat1$jd, paste('Levodissipada', indice , sep="")   )
    
    SalvaCorrelacao(wat1$ji, paste('Levoinjetada', indice , sep="")   )
    
    SalvaCorrelacao(wat1$tempo, paste('Levotempo', indice , sep="")   )
    
    
}



TestinEvolucaoS <- function(){
    

    
    #  corpre <- 'corpre'
    #corpos <- 'corpos'
    # resposta <- 'resposta'
    
    # csv<-'.csv'
    
    # for (i in 1:90)
    #{
    
        i<-3
        
        #titulo0<-paste( a , toString(i) , csv , sep="")
        
        #titulo1<-paste( b , toString(i) , csv , sep="")
        
        #tituloR<-paste( c , toString(i) , csv , sep="")
        
        # SimulacaoBimodalKickF <- function(ns, nt, nc, dt, massa, gamma, ctelastica, alfa, dalfa, xa, dxa, deltax, tt)
        
        #temp
        
        qero<-SimulacaoBimodalKickF(10,10000,100,0.001, 1,1,1, 1,0,1,0,1000, (10*i) + 1 )
        
        ns<- qero$ns
        
        kick<- qero$kick
        tkick<- qero$tkick /(qero$dt*qero$nc)
        
        Dt<- 100
        
        fim<- numeric(Dt)
        
        corpre<-      (  funcaodecorrelacaoSample(qero$posicao, (tkick-Dt -1) , Dt, ns)   )
        
        corpos<- ( funcaodecorrelacaoSample(qero$posicao, tkick , Dt ,ns))
        
        #SalvaCorrelacao(cor0, '')
        
        fim<- fim+ (corpos - corpre)
        
        fim<-fim/kick
        
        fim<- fim/ns
        
        # SalvaCorrelacao(corpre, titulo0)
        
        #SalvaCorrelacao(corpos, titulo1)
        
        #SalvaCorrelacao(fim, tituloR)
        
        final<-list("corpre" = corpre , "corpos" = corpos ,  "reposta"=fim)
        
        return(final)
    
}





TPC <- function(){
    
    a <- 'cor1init.csv'
    b <- 'cor2init.csv'
    
    # SimulacaoBimodalTPC1 <- function(ns, nt, nc, dt, massa, gamma, ctelastica, alfa, dalfa, xa, dxa, var)
    
    qero1<-SimulacaoBimodalTPC1(1000, 1000, 100, 0.001, 1, 1, 1, 1, 0, 1, 0, 5)
    
    cor1<-matrix (nrow =100, ncol =10)
    
    #xx<-matrix( nrow = nt, ncol = ns)
    
    
    
    #funcaodecorrelacaoSample<-function(vetor, t1, lagmax, ns)
    
    for (i in 1:9)
    {
        
        print('aaiosnasidnaisudiunasd')
        
        cor1[,i]<-      (  funcaodecorrelacaoSample(qero1$posicao, 1 + ((i-1)*10), 100, 1000)   )
        
        titulo1<- paste(toString(i) , a, sep="")
        
        SalvaCorrelacao(cor1[,i], titulo1)
        
    }
    
    qero2<-SimulacaoBimodalTPC2(1000, 1000, 100, 0.001, 1, 1, 1, 1, 0, 1, 0, 5)
    
    cor2<-matrix (nrow =100, ncol =10)
    
    for (i in 1:8)
    {
        
        print('aaiosnasidnaisudiunasd')
        
        cor2[,i]<-      (  funcaodecorrelacaoSample(qero2$posicao, 1 + ((i-1)*10), 100, 1000)   )
        
        titulo2<- paste(toString(i) , b, sep="")
        
        SalvaCorrelacao(cor2[,i], titulo2)
        
    }
    
    
}
    

#RespostaLinearTradTeste<-function(listasim, Dt){
#
    
    
#    ns<- listasim$ns
    
#   kick<- listasim$kick
#    tkick<- listasim$tkick /(listasim$dt*listasim$nc)
    
    
#    fim<- numeric(Dt)
    
#   cor0<-     - (  funcaodecorrelacaoSample(listasim$posicao, tkick-Dt -1, Dt, ns)   )
    
#   cor1<- ( funcaodecorrelacaoSample(listasim$posicao, tkick , Dt ,ns)
    
    #SalvaCorrelacao(cor0, '')
    
    
#   fim<- fim+ cor1 - cor0
    
    
#   fim<-fim/kick
    
#   fim<- fim/ns
    
    
#   return(fim)
    
#}
deltamove <- function( ns, nt, dt, xa, alfa){
    
    
    
    xx<- matrix( nrow = nt , ncol = ns )
    vv<- matrix( nrow = nt , ncol = ns )
    
    x<- 0
    
    for (j in 1:ns)
    {
        v<- RuidoBimodalSimetrico ( nt,dt,xa,alfa )
        
        for (i in 1:nt){
            x<- x+ v[i]*dt
            xx[i,j]<-x
            vv[i,j]<-v[i]
        }
        
    }
    
    final<-list("posicao"=xx,"velocidade"=vv)
    return(final)
}

RuidoBimodalSimetrico <- function(nt, dt ,xa, alfa){
    
    ruido<-numeric(nt)
    
    
    it <- 0
    eta<-xa
    
    xb<-   -xa
    
    while(it < nt)
    {
        if(eta==xa)
        {
            teta<- - log(1-runif(1))/(alfa/2)
            neta<- teta%/%dt + 1
        } else {
            teta<- - log(1-runif(1))/(alfa/2)
            neta<- teta%/%dt + 1
        }
        ie<-0
        while(ie<neta & it < nt)
        {
            it<- it+1
            ie<- ie+1
    
            ruido[it]<- eta
        }
        
        
        if(eta==xa)
        {
            eta<-xb
        } else {
            eta<-xa
        }
        
    }
    
    return(ruido)
}

RuidoBimodal <- function(nt, dt, xa, xb, mua, mub){
    
    ruido<-numeric(nt)
    
    
    it <- 0
    eta<-xa
    
    while(it < nt)
    {
        if(eta==xa)
        {
            teta<- - log(1-runif(1))/(mub)
            neta<- teta%/%dt + 1
        } else {
            teta<- - log(1-runif(1))/(mua)
            neta<- teta%/%dt + 1
        }
        ie<-0
        while(ie<neta & it < nt)
        {
            it<- it+1
            ie<- ie+1
            
            ruido[it]<- eta
        }
        
        if(eta==xa)
        {
            eta<-xb
        } else {
            eta<-xa
        }
        
        
    }
    
    return(ruido)
}

RuidoBimodalSample <- function(ns,nt, dt, xa, xb, mua, mub){
    
    
    
    
    ruido<- matrix(nrow = nt, ncol = ns)
    
    for(i in 1:ns)
    {
    print(i)
    it <- 0
    eta<-xa
    
    while(it < nt)
    {
        if(eta==xa)
        {
            teta<- - log(1-runif(1))/(mub)
            neta<- teta%/%dt + 1
        } else {
            teta<- - log(1-runif(1))/(mua)
            neta<- teta%/%dt + 1
        }
        ie<-0
        while(ie<neta & it < nt)
        {
            it<- it+1
            ie<- ie+1
            
            ruido[it,i]<- eta
        }
        
        if(eta==xa)
        {
            eta<-xb
        } else {
            eta<-xa
        }
        
        
    }
    
    
    }
    return(ruido)
}


SimulacaoBimodalTcte <- function(ns, nt, nc, dt, massa, gamma, ctelastica, alfa, dalfa, Temp){
    
    tempo<-seq(dt*nc, nt*dt*nc, dt*nc)
    
    ntc<-nt*nc
    
    
    xx<-matrix( nrow = nt, ncol = ns)
    vv<-matrix( nrow = nt, ncol = ns)
    ee<-matrix( nrow = nt, ncol = ns)
    ji<-matrix( nrow = nt, ncol = ns)
    jd<-matrix( nrow = nt, ncol = ns)

    ua<- alfa/2
    
    ub<- ua
    
    c1<- ctelastica/massa
    
    gm<- gamma/massa
    
    alfai<-alfa
    
    ie<-0
    
    it<-0
    
    xai <-numeric(ns)
    

    for(i in 1:ns)
    {
        x<-0
        v<-0
        v12<-0
        ji1<-0
        jd1<-0
        alfai<- alfai+ dalfa
        #AQUI A LINHA QUE GARANTE TEMPERATURA CONSTANTE = 1, para g=1 cm=1 we=1 alfa=1 xa=1
        
        xa<- sqrt( Temp*(gamma*(ctelastica+alfai*(gamma+massa*alfai)))/(massa*alfai)   )
        
        xai[i]<-xa
        
        xb<- -xa
        ua<-alfai/2
        ub<-ua
        eta<-xa
        it<-0
        ic<-0
        ir<-0
        ig<-0
        ie<-0
        while(it < ntc)
        {
            if(eta==xa)
            {
                teta<- - log(1-runif(1))/ub
                neta<- teta%/%dt + 1
            } else {
                teta<- - log(1-runif(1))/ua
                neta<- teta%/%dt + 1
            }
            ie<-0
            while(ie<neta & it < ntc)
            {
                
                it<-it+1
                ie<-ie+1
                u<-v
                y<-x
                force<-  -gm*u -c1*y + eta/massa
                
                v<- u + force*dt
                x<- y+ v*dt
                #print(x)
                ji1<-ji1+ (v+u)*eta*dt/2
                #jd1<-jd1+ gamma*(v^2+u^2)*dt/2
                jd1 <- jd1 + gamma*(((v+u)/2)^2)*dt
                
                ig<-ig+1
                
                if(ig==nc)
                {
                    #print(x)
                    #print(v)
                    xx[it/nc,i]<- x
                    vv[it/nc,i]<- v
                    ee[it/nc,i]<- (ji1-jd1)
                    ji[it/nc,i]<- (ji1)
                    jd[it/nc,i]<- (jd1)
                    
                    ig<-0
                }
            }
            if(eta==xa)
			{
                eta<-xb
			} else {
                eta<-xa
			}
        }
        
    }
    
    
    
    
    
    #function(ns, nt, nc, dt, massa, gamma, ctelastica, alfa, dalfa, Temp)
    final<-list("posicao"=xx,"velocidade"=vv, "energia"=ee, "ji" = ji, "jd" = jd , "massa"=massa,"gamma"=gamma,"k"=ctelastica, "alfa" = alfa, "dalfa"= dalfa, "xa"=xai, "ns"=ns,"nt"=nt,"nc"=nc,"dt"= dt, "temp"=Temp, "t"=tempo)
    return(final)
}


SimulacaoBimodalXActe <- function(ns, nt, nc, dt, massa, gamma, ctelastica, alfa, dalfa, xa){
    
    
    ntc<-nt*nc
    
    
    xx<-matrix( nrow = nt, ncol = ns)
    vv<-matrix( nrow = nt, ncol = ns)
    ee<-matrix( nrow = nt, ncol = ns)
    
    ua<- alfa/2
    
    ub<- ua
    
    c1<- ctelastica/massa
    
    gm<- gamma/massa
    
    alfai<-alfa
    
    ie<-0
    
    it<-0
    
    
    #Temp<- numeric(ns)
    
    for(i in 1:ns)
    {
        
        Temp<- ((xa^2)*(massa * alfa))/(gamma * (ctelastica + alfa * (massa * alfa + gamma)  ))
        
        x<-0
        v<-0
        v12<-0
        ji1<-0
        jd1<-0
        alfai<- alfai+ dalfa

        xb<- -xa
        ua<-alfai/2
        ub<-ua
        eta<-xa
        it<-0
        ic<-0
        ir<-0
        ig<-0
        ie<-0
        while(it < ntc)
        {
            if(eta==xa)
            {
                teta<- - log(1-runif(1))/ub
                neta<- teta%/%dt + 1
            } else {
                teta<- - log(1-runif(1))/ua
                neta<- teta%/%dt + 1
            }
            ie<-0
            while(ie<neta & it < ntc)
            {
                
                it<-it+1
                ie<-ie+1
                u<-v
                y<-x
                force<-  -gm*u -c1*y + eta/massa
                
                v<- u + force*dt
                x<- y+ v*dt
                #print(x)
                ji1<-ji1+ (v+u)*eta*dt/2
                #jd1<-jd1+ gamma*(v^2+u^2)*dt/2
                jd1 <- jd1 + gamma*(((v+u)/2)^2)*dt
                
                
                ig<-ig+1
                
                if(ig==nc)
                {
                    #print(x)
                    #print(v)
                    xx[it/nc,i]<- x
                    vv[it/nc,i]<- v
                    ee[it/nc,i]<- (ji1-jd1)
                    ig<-0
                }
            }
            if(eta==xa)
			{
                eta<-xb
			} else {
                eta<-xa
			}
        }
        
    }
    
    
    
    
    
    #function(ns, nt, nc, dt, massa, gamma, ctelastica, alfa, dalfa, Temp)
    final<-list("posicao"=xx,"velocidade"=vv, "energia"=ee, "massa"=massa,"gamma"=gamma,"k"=ctelastica, "alfa" = alfa,"dalfa"= dalfa, "xa"=xa, "ns"=ns,"nt"=nt,"nc"=nc,"dt"= dt, "temp"=Temp)
    return(final)
}


# hhv<-hist(listasim$velocidade[tkick:tkick+Dt,ns], breaks = Breakks)

# ssx0<- - hhx0$density*log(hhx0$density)
# ssx0[which(is.na(ssx0))]<- 0
# sx0<-sum(ssx0)*(1/sum(hhx0$density))
#ssv0<- - hhv0$density*log(hhv0$density)
#ssv0[which(is.na(ssv0))]<- 0
#sv0<-sum(ssv0)*(1/sum(hhv0$density))

#ssx<- - hhx$density*log(hhx$density)
#ssx[which(is.na(ssx))]<- 0
#sx<-sum(ssx)*(1/sum(hhx$density))
#ssv<- - hhv$density*log(hhv$density)
#ssv[which(is.na(ssv))]<- 0
#sv<-sum(ssv)*(1/sum(hhv$density))


#stoth[i]<- sx + sv

autocorr2<-function(vetor, lag)
{
    media<- 0
    med<- mean(vetor)
    
    for(i in 1:(length(vetor)-lag)  )
    
    {
        media<- media + (vetor[i] - med)*(vetor[i+lag] - med)
    }
    media <- media/(length(vetor) - lag)
    
    
    
    return (media)
}

funcaodecorrelacao2<-function(vetor, lagmax)
{
    media<-numeric(lagmax)
    variancia<- var(vetor)
    
    for(i in 1:lagmax)
    {
        media[i]<-autocorr2(vetor,i)/variancia
    }
    
    return(media)
}



autocorr<-function(vetor, lag)
{
    media<- 0
  
    for(j in 1:  (length(vetor)-lag)  )
    {
        #print(j)
        if (j==0)
        {
            break
        }
        media<- media + vetor[j]*vetor[j+lag]
    }
    media <- media/(length(vetor)-lag)
    
    return (media)
}

funcaodecorrelacao<-function(vetor, lagmax)
{
    fim<-numeric(lagmax)
    
    for(ii in 1:lagmax)
    {
        fim[ii]<-autocorr(vetor,ii)
    }

    return(fim)
}

corrsample<-function(vetor, t1,t2,ns)
{
    
    #s<- abs(t1-t2)
    
    fim<- 0
    
    for(ii in 1:ns)
    {
        fim<- fim + (vetor[t1,ii]*vetor[t2,ii])
    }
    
    fim<-fim/ns
    
    return(fim)
    
}

corrsample2<-function(vetor1, vetor2, t1,t2,ns)
{
    
    #s<- abs(t1-t2)
    
    fim<- 0
    
    for(ii in 1:ns)
    {
        fim<- fim + (vetor1[t1,ii]*vetor2[t2,ii])
    }
    
    fim<-fim/ns
    
    return(fim)
    
}

funcaodecorrelacaoSample<-function(vetor, t1, lagmax, ns)
{
    
    fim<-numeric(lagmax)
    
    for(jj in 1:lagmax)
    {
        fim[jj]<-   corrsample(vetor,t1,(t1+jj), ns)
    }
    
    return(fim)
}


funcaodecorrelacaoSample2<-function(vetor1, vetor2, t1, lagmax, ns)
{
    
    fim<-numeric(lagmax)
    
    for(jj in 1:lagmax)
    {
        fim[jj]<-   corrsample2(vetor1, vetor2,t1,(t1+jj), ns)
    }
    
    return(fim)
}


RespostaLinearTrad<-function(listasim, Dt){
    
    
    
    ns<- listasim$ns
    
    kick<- listasim$kick
    tkick<- listasim$tkick /(listasim$dt*listasim$nc)
    
    
    fim<- numeric(Dt)
    
    cor0<-     - (  funcaodecorrelacaoSample(listasim$posicao, tkick-Dt -1, Dt, ns)   )
    
    cor1<- ( funcaodecorrelacaoSample(listasim$posicao, tkick , Dt ,ns))
    
    #SalvaCorrelacao(cor0, '')
    
    
    fim<- fim+ cor1 - cor0
    
    
    fim<-fim/kick
    
    fim<- fim/ns
    
    
    return(fim)
    
}



Susceptibilidade<-function(resposta, dt){
    
    fim <- sum(resposta)
    fim <- fim*dt
    return(fim)
    
}

TesteKickX<-function(kick, dkick, n)
{
    
    listakick <- seq(kick, kick+ (n-1)*dkick, dkick)

    finis <- numeric(n)
    
    for (i in 1:n)
    {
        
        qero<-SimulacaoBimodalProtocoloX(100,1000,100,0.001,1,1,1,1,0,1,0,listakick[i],50  )
        wat<-RespostaLinearTrad(qero,100)
        finis[i]<- Susceptibilidade(wat,0.001)
        
    }
    
    return(finis)
    
}

TesteKickV<-function(kick, dkick, n)
{
    listakick <- seq(kick, kick+ (n-1)*dkick, dkick)
    # vou deixar os parametros termomecanicos tudo em 1 logo
    finis <- numeric(n)
    
    for (i in 1:n)
    {
        
        qero<-SimulacaoBimodalProtocoloV(100,1000,100,0.001,1,1,1,1,0,1,0,listakick[i],50  )
        wat<-RespostaLinearTrad(qero,100)
        finis[i]<- Susceptibilidade(wat,0.001)
    }
    
    return(finis)
}

RespostaLinearTradv<-function(listasim, Dt){
    
    
    
    ns<- listasim$ns
    
    kick<- listasim$kick
    tkick<- listasim$tkick /(listasim$dt*listasim$nc)
    
    
    fim<- numeric(Dt)
    
    for(i in 1:ns)
    {
        print(i)
        
        fim<- fim+  ( funcaodecorrelacao(listasim$velocidade[1:(tkick+Dt),i], Dt) -   funcaodecorrelacao(listasim$velocidade[1:(tkick-1),i], Dt)   )
        
    }
    
    fim<-fim/kick
    
    fim<- fim/ns
    
    
    return(fim)
    
}

FuncaoRespostaTrad<-function(listasim, Dtmax){
    
}

RespostaLinearMediaSampleVulpiLB <-function(listasim, Dt, Breakks)
{
    
    ns<- listasim$ns
    
    kick<- listasim$kick
    
    tkick<- listasim$tkick   /(listasim$dt*listasim$nc)
    
    #falta o nc aqui
    
    #print(tkick)
    #print(listasim$nt)
    resposta<- 0
    
    for(i in 1:ns)
    {
        hhx0<-hist(listasim$posicao[1:tkick-1 , ns], breaks = Breakks)
        
        #hhx<-hist(listasim$posicao[1:tkick+Dt , ns], breaks = Breakks)
        
        #print(i)
        
           for(j in 1:(length(hhx$breaks)-1) )
        {
           
            if ( hhx$breaks[j] <= listasim$posicao[tkick-1,i] && listasim$posicao[tkick-1,i] < hhx$breaks[j+1] )
            {
               k0<-j
            }
            
            if ( hhx$breaks[j] <= listasim$posicao[tkick,i] && listasim$posicao[tkick,i] < hhx$breaks[j+1] )
            {
                k1<-j
            }
 
        }
        
        
        resposta<- resposta-    listasim$posicao[tkick+Dt,i]*(log(hhx0$density[k1]) - log(hhx0$density[k0]))/kick
        
    }
    
    
    resposta<- resposta/ns
    
    return(resposta)
}

FuncaoRespostaVulpiLB <-function(listasim, Dtmax, Breakks)
{
    resposta<-numeric(Dtmax)
    
    
    for(i in 1:Dtmax)
    {
        print(i)
        resposta[i]<-RespostaLinearMediaSampleVulpiLB(listasim,i,Breakks)
    }
    
    return (resposta)
}


RespostaLinearSampleVulpiLorenz1 <-function(qero, Dt, Breakks)
{
    
    ns<- qero$ns
    
    kick<- qero$kick
    
    tkick<- qero$tkick  # /(listasim$dt*listasim$nc) - pro bimodal
    
    
    resposta<-numeric(Dt)
    
   
    hhx0<-hist(qero$z[1:tkick-1 , ], breaks = Breakks, plot = FALSE)
    
    k1<-numeric(Dt)
    
    for(i in 1:(ns))
    {
        #print(i)
        
    for(j in 1:(length(hhx0$breaks)-1) )
    {
        #print(j)
        if ( hhx0$breaks[j] <= qero$z[tkick-1,i] && qero$z[tkick-1,i] < hhx0$breaks[j+1] )
        {
            k0<-j
            break
        }

    }
    
    if(  qero$z[tkick-1,i] == hhx0$breaks[length(hhx$breaks)] )
    {
        k0<-j
    }
    
    for(k in 1:Dt)
    
    {
        #print(k)
        for(j in 1:(length(hhx0$breaks)-1) )
        {
            if ( hhx0$breaks[j] <= qero$z[tkick-1+k,i] && qero$z[tkick-1+k,i] < hhx0$breaks[j+1] )
            {
                # print('foi')
                k1[k]<-j
                break
            }
        }
        
        #print(k)
        #print(k1)
        #print(qero$z[tkick+k,i])
        #print(log(hhx0$density[k1]))
        #print(log(hhx0$density[k0]))
        resposta[k]<-  resposta[k] - qero$z[tkick+k,i]*(log(hhx0$density[k1[k]]) - log(hhx0$density[k0]))/( qero$z[tkick+k,i] - qero$z[tkick-1,i]  )  #kick
        
    }
    
    }
    
    resposta<- resposta/ns
    
    return(resposta)
    
}

RespostaLinearSampleVulpiLorenz2 <-function(qero, Dt, Breakks) # versao das eq 3.6 e 3.7 no response theory
{
    
    ns<- qero$ns
    kick<- qero$kick
    tkick<- qero$tkick  # /(listasim$dt*listasim$nc) - pro bimodal
    
    resposta<-numeric(Dt)
    
    hhx0<-hist(qero$z[1:(tkick-1) , ], breaks = Breakks, plot = FALSE)
    
    
    for(i in 1:(ns))
    {
        for(j in 1:(length(hhx0$breaks)-1) )
        {
            if ( hhx0$breaks[j] <= qero$z[tkick-1,i] && qero$z[tkick-1,i] < hhx0$breaks[j+1] )
            {
                k0<-j
                break
            }
        }
        
        if(  qero$z[tkick-1,i] == hhx0$breaks[length(hhx$breaks)] )
        {
            k0<-j
        }
        
        for(j in 1:(length(hhx0$breaks)-1) )
        {
            if ( hhx0$breaks[j] <= (qero$z[tkick-1,i] - kick) && (qero$z[tkick-1,i] - kick) < hhx0$breaks[j+1] )
            {
                k1<-j
                break
            }
        }
        
        if(  (qero$z[tkick-1,i] - kick) == hhx0$breaks[length(hhx$breaks)] )
        {
            k1<-j
        }
        
        for(k in 1:Dt)
        {
            resposta[k]<-  resposta[k] - qero$z[tkick+k,i]*(   hhx0$density[k1] - hhx0$density[k0]   )/(hhx0$density[k0] )
        }
    }
    resposta<- resposta/ns
    return(resposta)
    
}

RespostaLinearSampleVulpiLorenz3 <-function(listasim, Dt, Breakks)
{
    
    ns<- listasim$ns
    
    kick<- listasim$kick
    
    tkick<- listasim$tkick  # /(listasim$dt*listasim$nc)
    
    #falta o nc aqui
    
    #print(tkick)
    #print(listasim$nt)
    resposta<-numeric(Dt)
    
    #for(i in 1:ns)
    #{
    
    hhx0<-hist(listasim$z[1:tkick-1 , ], breaks = Breakks, plot = FALSE)
    
    
    for(i in 1:Dt)
    
    {
        hhx<-hist(listasim$z[1:tkick-1 +i , ], breaks = Breakks, plot = FALSE)
        
        
        for(j in 1:(length(hhx$breaks)-1) )
        {
            
            if ( hhx$breaks[j] <= listasim$z[tkick-1,i] && listasim$z[tkick-1,i] < hhx$breaks[j+1] )
            {
                k0<-j
            }
            
            if ( hhx$breaks[j] <= listasim$z[tkick,i] && listasim$z[tkick,i] < hhx$breaks[j+1] )
            {
                k1<-j
            }
            
        }
        
        #}
        
        
        resposta[i]<- resposta-    listasim$z[tkick+i,]*(log(hhx0$density[k1]) - log(hhx0$density[k0]))/kick
        
    }
    
    
    resposta<- resposta/ns
    
    return(resposta)
}



FuncaoRespostaVulpi <-function(listasim, Dtmax, Breakks)
{
    resposta<-numeric(Dtmax)
    
    
    for(i in 1:Dtmax)
    {
        print(i)
        resposta[i]<-RespostaLinearSampleVulpiLorenz(listasim,i,Breakks)
    }
    
    return (resposta)
}








ComparaCorrelacaoFuncaoresposta<-function()
{}


SimulacaoBimodalKickX <- function(ns, nt, nc, dt, massa, gamma, ctelastica, alfa, dalfa, xa, dxa, deltax, tt){
    
    
    ntc<-nt*nc
    
    xx<-matrix( nrow = nt, ncol = ns)
    vv<-matrix( nrow = nt, ncol = ns)
    ee<-matrix( nrow = nt, ncol = ns)
    
    alfazero<-alfa
    
    ua<- alfa/2
    
    ub<- ua
    
    c1<- ctelastica/massa
    
    gm<- gamma/massa
    
    ie<-0
    
    it<-0
    
    tempo<-seq(0,ntc*dt, dt*nc)
    
    for(i in 1:ns)
    {
        
        print(i)
        
        Temp<- ((xa^2)*(massa * alfa))/(gamma * (ctelastica + alfa * (massa * alfa + gamma)  ))
        
        x<-0
        v<-0
        v12<-0
        ji1<-0
        jd1<-0
        alfa<- alfa+ dalfa
        xa<- xa+dxa
        
        xb<- -xa
        ua<-alfa/2
        ub<-ua
        eta<-xa
        it<-0
        ic<-0
        ir<-0
        ig<-0
        ie<-0
        indicadorkick<-0
        
        while(it < ntc)
        {
            
            #print(it*dt)
            if( it*dt >= tt & indicadorkick==0)
            {
                x <- x + deltax
                print(i)
                indicadorkick<- 1
            }
            
            
            if(eta==xa)
            {
                teta<- - log(1-runif(1))/ub
                neta<- teta%/%dt + 1
            } else {
                teta<- - log(1-runif(1))/ua
                neta<- teta%/%dt + 1
            }
            ie<-0
            while(ie<neta & it < ntc)
            {
                it<-it+1
                ie<-ie+1
                
                u<-v
                y<-x
                force<-  -gm*u -c1*y + eta/massa
                v<- u + force*dt
                x<- y+ v*dt
                
                ji1<-ji1+ (v+u)*eta*dt/2
                #jd1<-jd1+ gamma*(v^2+u^2)*dt/2
                jd1 <- jd1 + gamma*(((v+u)/2)^2)*dt
                
                ig<-ig+1
                
                if(ig==nc)
                {
                    xx[it/nc,i]<- x
                    vv[it/nc,i]<- v
                    ee[it/nc,i]<- (ji1-jd1)
                    ig<-0
                }
            }
            if(eta==xa)
            {
                eta<-xb
            } else {
                eta<-xa
            }
        }
        
    }
    final<-list("posicao"=xx,"velocidade"=vv, "energia"=ee, "massa"=massa,"gamma"=gamma,"k"=ctelastica, "alfa" = alfazero,"dalfa"= dalfa, "xa"=xa, "dxa"=dxa, "ns" = ns ,"nt"=nt,"nc"=nc,"dt"= dt, "temp"=Temp, "tkick" = tt, "kick" = deltax, "tempo"= tempo)
    
    return(final)
}



SimulacaoBimodalKickXsalvaEta <- function(ns, nt, nc, dt, massa, gamma, ctelastica, alfa, dalfa, xa, dxa, deltax, tt){
    
    
    ntc<-nt*nc
    
    xx<-matrix( nrow = nt, ncol = ns)
    vv<-matrix( nrow = nt, ncol = ns)
    ee<-matrix( nrow = nt, ncol = ns)
    ruido<-matrix( nrow = nt, ncol = ns)
    
    alfazero<-alfa
    
    ua<- alfa/2
    
    ub<- ua
    
    c1<- ctelastica/massa
    
    gm<- gamma/massa
    
    ie<-0
    
    it<-0
    
    tempo<-seq(0,ntc*dt, dt*nc)
    
    for(i in 1:ns)
    {
        
        print(i)
        
        Temp<- ((xa^2)*(massa * alfa))/(gamma * (ctelastica + alfa * (massa * alfa + gamma)  ))
        
        x<-0
        v<-0
        v12<-0
        ji1<-0
        jd1<-0
        alfa<- alfa+ dalfa
        xa<- xa+dxa
        
        xb<- -xa
        ua<-alfa/2
        ub<-ua
        eta<-xa
        it<-0
        ic<-0
        ir<-0
        ig<-0
        ie<-0
        indicadorkick<-0
        
        while(it < ntc)
        {
            
            #print(it*dt)
            if( it*dt >= tt & indicadorkick==0)
            {
                x <- x + deltax
                print(i)
                indicadorkick<- 1
            }
            
            
            if(eta==xa)
            {
                teta<- - log(1-runif(1))/ub
                neta<- teta%/%dt + 1
            } else {
                teta<- - log(1-runif(1))/ua
                neta<- teta%/%dt + 1
            }
            ie<-0
            while(ie<neta & it < ntc)
            {
                it<-it+1
                ie<-ie+1
                
                u<-v
                y<-x
                force<-  -gm*u -c1*y + eta/massa
                v<- u + force*dt
                x<- y+ v*dt
                
                ji1<-ji1+ (v+u)*eta*dt/2
                #jd1<-jd1+ gamma*(v^2+u^2)*dt/2
                jd1 <- jd1 + gamma*(((v+u)/2)^2)*dt
                
                ig<-ig+1
                
                if(ig==nc)
                {
                    xx[it/nc,i]<- x
                    vv[it/nc,i]<- v
                    ee[it/nc,i]<- (ji1-jd1)
                    ruido[it/nc,i]<-eta
                    ig<-0
                }
            }
            if(eta==xa)
            {
                eta<-xb
            } else {
                eta<-xa
            }
        }
        
    }
    final<-list("posicao" = xx ,"velocidade" = vv, "energia" = ee, "ruido" = ruido, "massa" = massa,"gamma"=gamma,"k"=ctelastica, "alfa" = alfazero,"dalfa"= dalfa, "xa"=xa, "dxa"=dxa, "ns" = ns ,"nt"=nt,"nc"=nc,"dt"= dt, "temp"=Temp, "tkick" = tt, "kick" = deltax, "tempo"= tempo)
    
    return(final)
}

RespostaKickHeta<-function(qero, Dt){
        
        ns<- qero$ns
        kick<- qero$kick
        lambda <- qero$lambda
        tkick<- qero$tkick/(qero$dt*qero$nc)
        print(tkick)
        
        respostaa<- numeric(Dt)
        suceptibilidade <- numeric(Dt)
        XTeff<- numeric(Dt)
        
        print('1foi até aqui')
        print('passou o corrsample')
        
        for (i in 1:Dt)
        {
            mediapos <-  mean(qero$posicao[tkick+i, ]  )
            respostaa[i] <- mediapos/(kick* exp(- lambda * i))
            suceptibilidade[i] <- sum(respostaa[1:i]) * qero$dt
        }
        
        #iii <-  iii + 1
        
        assimfdtpre<- AssimetriaFDT(qero, tkick-Dt-1, Dt)
        assimfdtpos<- AssimetriaFDT(qero, tkick, Dt)
        
        assimpre <-Assimetria(qero, tkick-Dt-1, Dt)
        assimpos <-Assimetria(qero, tkick, Dt)
        
        assimtotpre <-Assimetriatot(qero, tkick-Dt-1, Dt)
        assimtotpos <-Assimetriatot(qero, tkick, Dt)
        
        # assim
        
        
        final<-list("assimfdtpre" = assimfdtpre, "assimfdtpos" = assimfdtpos ,"assimpre" = assimpre,"assimpos" = assimpos ,"resposta" = respostaa)
        
        return(final)
        
    

    
    
    
    
}

PrintRespostaKickH<-function()
{
    
    tituloR<-paste( rresposta , csv , sep="")
    tituloX<-paste( sucep , csv , sep="")
    tituloT<-paste( Xteff , csv , sep="")
    tituloApre<-paste( Assimpre , csv , sep="")
    tituloApos<-paste( Assimpos , csv , sep="")
    
    corrpre <- paste( letra ,  'corrpre', toString(vez)  , sep="")
    rresposta <- paste( letra, 'resposta', toString(vez) , sep="")
    
    sucep <- paste( letra, 'sucep', toString(vez) , sep="")
    Xteff <- paste( letra, 'XTeff', toString(vez) , sep="")
    
    Assimpre <- paste( letra, 'Assimpre', toString(vez) , sep="")
    Assimpos <- paste( letra, 'Assimpos', toString(vez) , sep="")
    
    csv<-'.csv'
}


RespostaKickHeta2<-function(qero, Dt){
    
    ns<- qero$ns
    kick<- qero$kick
    lambda <- qero$lambda
    tkick<- qero$tkick/(qero$dt*qero$nc)
    print(tkick)
    
    respostaa<- numeric(Dt)
    suceptibilidade <- numeric(Dt)
    XTeff<- numeric(Dt)
    
    print('1foi até aqui')
    print('passou o corrsample')
    
    for (i in 1:Dt)
    {
        mediapos <-  mean(qero$posicao[tkick+i, ]  )
        respostaa[i] <- mediapos/(kick)
        suceptibilidade[i] <- sum(respostaa[1:i]) * qero$dt
    }
    
    #iii <-  iii + 1
    
    assimfdtpre<- AssimetriaFDT(qero, tkick-Dt-1, Dt)
    assimfdtpos<- AssimetriaFDT(qero, tkick, Dt)
    
    assimpre <-Assimetria(qero, tkick-Dt-1, Dt)
    assimpos <-Assimetria(qero, tkick, Dt)
    
    assimtotpre <-Assimetriatot(qero, tkick-Dt-1, Dt)
    assimtotpos <-Assimetriatot(qero, tkick, Dt)
    
    # assim
    
    
    final<-list("assimfdtpre" = assimfdtpre, "assimfdtpos" = assimfdtpos ,"assimpre" = assimpre,"assimpos" = assimpos ,"resposta" = respostaa)
    
    return(final)
    
    
    
    
    
    
    
}

SimBimodKickHeta <- function(ns, nt, nc, dt, massa, gamma, ctelastica, alfa, dalfa, xa, dxa, deltaf, tt,lambda, Dt){
    
    
    ntc<-nt*nc
    
    xx<-matrix( nrow = nt, ncol = ns)
    vv<-matrix( nrow = nt, ncol = ns)
    ee<-matrix( nrow = nt, ncol = ns)
    ruido<-matrix( nrow = nt, ncol = ns)
    
    alfazero<-alfa
    
    ua<- alfa/2
    
    ub<- ua
    
    c1<- ctelastica/massa
    
    gm<- gamma/massa
    
    ie<-0
    
    it<-0
    
    tempo<-seq(0,ntc*dt, dt*nc)
    
    for(i in 1:ns)
    {
        
        print(i)
        
        Temp<- ((xa^2)*(massa * alfa))/(gamma * (ctelastica + alfa * (massa * alfa + gamma)  ))
        
        x<-0
        v<-0
        v12<-0
        ji1<-0
        jd1<-0
        alfa<- alfa+ dalfa
        xa<- xa+dxa
        
        xb<- -xa
        ua<-alfa/2
        ub<-ua
        eta<-xa
        it<-0
        ic<-0
        ir<-0
        ig<-0
        ie<-0
        indicadorkick<-0
        
        while(it < ntc)
        {
            
            
            if(eta==xa)
            {
                teta<- - log(1-runif(1))/ub
                neta<- teta%/%dt + 1
            } else {
                teta<- - log(1-runif(1))/ua
                neta<- teta%/%dt + 1
            }
            ie<-0
            while(ie<neta & it < ntc)
            {
                it<-it+1
                ie<-ie+1
                
                u<-v
                y<-x
                
                force<-  -gm*u -c1*y + eta/massa
                
                
                #print(it*dt)
                if( it*dt >= tt & indicadorkick < Dt  )
                {
                    force <- force + deltaf *  exp(- lambda * indicadorkick)
                    print(i)
                    indicadorkick<- indicadorkick + 1
                }
                
                
                v<- u + force*dt
                x<- y+ v*dt
                
                ji1<-ji1+ (v+u)*eta*dt/2
                #jd1<-jd1+ gamma*(v^2+u^2)*dt/2
                jd1 <- jd1 + gamma*(((v+u)/2)^2)*dt
                
                ig<-ig+1
                
                if(ig==nc)
                {
                    xx[it/nc,i]<- x
                    vv[it/nc,i]<- v
                    ee[it/nc,i]<- (ji1-jd1)
                    ruido[it/nc,i]<-eta
                    ig<-0
                }
            }
            if(eta==xa)
            {
                eta<-xb
            } else {
                eta<-xa
            }
        }
        
    }
    final<-list("posicao" = xx ,"velocidade" = vv, "energia" = ee, "ruido" = ruido, "massa" = massa,"gamma"=gamma,"k"=ctelastica, "alfa" = alfazero,"dalfa"= dalfa, "xa"=xa, "dxa"=dxa, "ns" = ns ,"nt"=nt,"nc"=nc,"dt"= dt, "temp"=Temp, "tkick" = tt, "kick" = deltaf, "tempo"= tempo, "lambda" = lambda)
    
    return(final)
}

SimulacaoBimodalKickXsalvaEtaOD <- function(ns, nt, nc, dt, massa, gamma, ctelastica, alfa, dalfa, xa, dxa, deltax, tt){
    
    
    ntc<-nt*nc
    
    xx<-matrix( nrow = nt, ncol = ns)
    vv<-matrix( nrow = nt, ncol = ns)
    ee<-matrix( nrow = nt, ncol = ns)
    ruido<-matrix( nrow = nt, ncol = ns)
    
    alfazero<-alfa
    
    ua<- alfa/2
    
    ub<- ua
    
    c1<- ctelastica/massa
    
    gm<- gamma/massa
    
    ie<-0
    
    it<-0
    
    tempo<-seq(0,ntc*dt, dt*nc)
    
    for(i in 1:ns)
    {
        
        print(i)
        
        Temp<- ((xa^2)*(massa * alfa))/(gamma * (ctelastica + alfa * (massa * alfa + gamma)  ))
        
        x<-0
        v<-0
        v12<-0
        ji1<-0
        jd1<-0
        alfa<- alfa+ dalfa
        xa<- xa+dxa
        
        xb<- -xa
        ua<-alfa/2
        ub<-ua
        eta<-xa
        it<-0
        ic<-0
        ir<-0
        ig<-0
        ie<-0
        indicadorkick<-0
        
        while(it < ntc)
        {
            
            #print(it*dt)
            if( it*dt >= tt & indicadorkick==0)
            {
                x <- x + deltax
                print(i)
                indicadorkick<- 1
            }
            
            
            if(eta==xa)
            {
                teta<- - log(1-runif(1))/ub
                neta<- teta%/%dt + 1
            } else {
                teta<- - log(1-runif(1))/ua
                neta<- teta%/%dt + 1
            }
            ie<-0
            while(ie<neta & it < ntc)
            {
                it<-it+1
                ie<-ie+1
                
                u<-v
                y<-x
                
                
                dx<- (- c1 * y * dt - eta)*(1/gm)
                
                x<-  y + dx
                
                v<-dx/dt
                
                ji1<-ji1+ (v+u)*eta*dt/2
                #jd1<-jd1+ gamma*(v^2+u^2)*dt/2
                jd1 <- jd1 + gamma*(((v+u)/2)^2)*dt
                
                ig<-ig+1
                
                if(ig==nc)
                {
                    xx[it/nc,i]<- x
                    vv[it/nc,i]<- v
                    ee[it/nc,i]<- (ji1-jd1)
                    ruido[it/nc,i]<-eta
                    ig<-0
                }
            }
            if(eta==xa)
            {
                eta<-xb
            } else {
                eta<-xa
            }
        }
        
    }
    final<-list("posicao" = xx ,"velocidade" = vv, "energia" = ee, "ruido" = ruido, "massa" = massa,"gamma"=gamma,"k"=ctelastica, "alfa" = alfazero,"dalfa"= dalfa, "xa"=xa, "dxa"=dxa, "ns" = ns ,"nt"=nt,"nc"=nc,"dt"= dt, "temp"=Temp, "tkick" = tt, "kick" = deltax, "tempo"= tempo)
    
    return(final)
}



SimulacaoBimodalKickF <- function(ns, nt, nc, dt, massa, gamma, ctelastica, alfa, dalfa, xa, dxa, deltaf, tt){
    
    
    ntc<-nt*nc
    
    xx<-matrix( nrow = nt, ncol = ns)
    vv<-matrix( nrow = nt, ncol = ns)
    ee<-matrix( nrow = nt, ncol = ns)
    #ff<-matrix( nrow = nt, ncol = ns)
    
    alfazero<-alfa
    
    ua<- alfa/2
    
    ub<- ua
    
    c1<- ctelastica/massa
    
    gm<- gamma/massa
    
    ie<-0
    
    it<-0
    
    tempo<-seq(0,ntc*dt, dt*nc)
    
    for(i in 1:ns)
    {
        
        print(i)
        
        Temp<- ((xa^2)*(massa * alfa))/(gamma * (ctelastica + alfa * (massa * alfa + gamma)  ))
        
        x<-0
        v<-0
        v12<-0
        ji1<-0
        jd1<-0
        alfa<- alfa+ dalfa
        xa<- xa+dxa
        
        xb<- -xa
        ua<-alfa/2
        ub<-ua
        eta<-xa
        it<-0
        ic<-0
        ir<-0
        ig<-0
        ie<-0
        #### será?
        ii<-0
        indicadorkick<-0
        
        while(it < ntc)
        {
            
            #print(it*dt)
            
            if(eta==xa)
            {
                teta<- - log(1-runif(1))/ub
                neta<- teta%/%dt + 1
            } else {
                teta<- - log(1-runif(1))/ua
                neta<- teta%/%dt + 1
            }
            ie<-0
            while(ie<neta & it < ntc)
            {
                it<-it+1
                ie<-ie+1
                
                u<-v
                y<-x
                force<-  -gm*u -c1*y + eta/massa
                
                
                
                if( it*dt >= tt & indicadorkick==0)
                {
                    force <- force + deltaf
                    #print(i)
                    #indicadorkick<- 1
                }
                
                v<- u + force*dt
                x<- y+ v*dt
                
                ji1<-ji1+ (v+u)*eta*dt/2
                #jd1<-jd1+ gamma*(v^2+u^2)*dt/2
                jd1 <- jd1 + gamma*(((v+u)/2)^2)*dt
                
                ig<-ig+1
                
                if(ig==nc)
                {
                    #ff[it/nc,i]<- force
                    xx[it/nc,i]<- x
                    vv[it/nc,i]<- v
                    ee[it/nc,i]<- (ji1-jd1)
                    ig<-0
                }
            }
            if(eta==xa)
            {
                eta<-xb
            } else {
                eta<-xa
            }
        }
        
    }
    final<-list("posicao"=xx,"velocidade"=vv, "energia"=ee, "massa"=massa,"gamma"=gamma,"k"=ctelastica, "alfa" = alfazero,"dalfa"= dalfa, "xa"=xa, "dxa"=dxa, "ns" = ns ,"nt"=nt,"nc"=nc,"dt"= dt, "temp"=Temp, "tkick" = tt, "kick" = deltaf, "tempo"= tempo) #"forca"= ff
    
    return(final)
}

SimulacaoBimodalKickFdelta <- function(ns, nt, nc, dt, massa, gamma, ctelastica, alfa, dalfa, xa, dxa, deltaf, tt){
    
    
    ntc<-nt*nc
    
    xx<-matrix( nrow = nt, ncol = ns)
    vv<-matrix( nrow = nt, ncol = ns)
    ee<-matrix( nrow = nt, ncol = ns)
    #ff<-matrix( nrow = nt, ncol = ns)
    
    alfazero<-alfa
    
    ua<- alfa/2
    
    ub<- ua
    
    c1<- ctelastica/massa
    
    gm<- gamma/massa
    
    ie<-0
    
    it<-0
    
    tempo<-seq(0,ntc*dt, dt*nc)
    
    for(i in 1:ns)
    {
        
        print(i)
        
        Temp<- ((xa^2)*(massa * alfa))/(gamma * (ctelastica + alfa * (massa * alfa + gamma)  ))
        
        x<-0
        v<-0
        v12<-0
        ji1<-0
        jd1<-0
        alfa<- alfa+ dalfa
        xa<- xa+dxa
        
        xb<- -xa
        ua<-alfa/2
        ub<-ua
        eta<-xa
        it<-0
        ic<-0
        ir<-0
        ig<-0
        ie<-0
        #### será?
        ii<-0
        indicadorkick<-0
        
        while(it < ntc)
        {
            
            #print(it*dt)
            
            if(eta==xa)
            {
                teta<- - log(1-runif(1))/ub
                neta<- teta%/%dt + 1
            } else {
                teta<- - log(1-runif(1))/ua
                neta<- teta%/%dt + 1
            }
            ie<-0
            while(ie<neta & it < ntc)
            {
                it<-it+1
                ie<-ie+1
                
                u<-v
                y<-x
                force<-  -gm*u -c1*y + eta/massa
                
                
                
                if( it*dt >= tt & indicadorkick==0)
                {
                    force <- force + deltaf
                    #print(i)
                    indicadorkick<- 1
                }
                
                v<- u + force*dt
                x<- y+ v*dt
                
                ji1<-ji1+ (v+u)*eta*dt/2
                #jd1<-jd1+ gamma*(v^2+u^2)*dt/2
                jd1 <- jd1 + gamma*(((v+u)/2)^2)*dt
                
                ig<-ig+1
                
                if(ig==nc)
                {
                    #ff[it/nc,i]<- force
                    xx[it/nc,i]<- x
                    vv[it/nc,i]<- v
                    ee[it/nc,i]<- (ji1-jd1)
                    ig<-0
                }
            }
            if(eta==xa)
            {
                eta<-xb
            } else {
                eta<-xa
            }
        }
        
    }
    
    final<-list("posicao"=xx,"velocidade"=vv, "energia"=ee, "massa"=massa,"gamma"=gamma,"k"=ctelastica, "alfa" = alfazero,"dalfa"= dalfa, "xa"=xa, "dxa"=dxa, "ns" = ns ,"nt"=nt,"nc"=nc,"dt"= dt, "temp"=Temp, "tkick" = tt, "kick" = deltaf, "tempo"= tempo) #"forca"= ff
    
    return(final)
}




SimulacaoBimodalKickFrandom <- function(ns, nt, nc, dt, massa, gamma, ctelastica, alfa, dalfa, xa, dxa, deltaf, tt){
    
    
    ntc<-nt*nc
    
    xx<-matrix( nrow = nt, ncol = ns)
    vv<-matrix( nrow = nt, ncol = ns)
    ee<-matrix( nrow = nt, ncol = ns)
    #ff<-matrix( nrow = nt, ncol = ns)
    kick<-numeric(ns)
    
    
    alfazero<-alfa
    
    ua<- alfa/2
    
    ub<- ua
    
    c1<- ctelastica/massa
    
    gm<- gamma/massa
    
    ie<-0
    
    it<-0
    
    tempo<-seq(0,ntc*dt, dt*nc)
    
    for(i in 1:ns)
    {
        
        print(i)
        
        Temp<- ((xa^2)*(massa * alfa))/(gamma * (ctelastica + alfa * (massa * alfa + gamma)  ))
        
        x<-0
        v<-0
        v12<-0
        ji1<-0
        jd1<-0
        alfa<- alfa+ dalfa
        xa<- xa+dxa
        
        xb<- -xa
        ua<-alfa/2
        ub<-ua
        eta<-xa
        it<-0
        ic<-0
        ir<-0
        ig<-0
        ie<-0
        #### será?
        ii<-0
        indicadorkick<-0
        
        while(it < ntc)
        {
            
            #print(it*dt)
            
            if(eta==xa)
            {
                teta<- - log(1-runif(1))/ub
                neta<- teta%/%dt + 1
            } else {
                teta<- - log(1-runif(1))/ua
                neta<- teta%/%dt + 1
            }
            ie<-0
            while(ie<neta & it < ntc)
            {
                it<-it+1
                ie<-ie+1
                
                u<-v
                y<-x
                force<-  -gm*u -c1*y + eta/massa
                
                
                
                if( it*dt >= tt & indicadorkick==0)
                {
                    kick<- deltaf*runif(1,min =-1,max=1)
                    force <- force + kick
                    #print(i)
                    #indicadorkick<- 1
                }
                
                v<- u + force*dt
                x<- y+ v*dt
                
                ji1<-ji1+ (v+u)*eta*dt/2
                #jd1<-jd1+ gamma*(v^2+u^2)*dt/2
                jd1 <- jd1 + gamma*(((v+u)/2)^2)*dt
                
                ig<-ig+1
                
                if(ig==nc)
                {
                    #ff[it/nc,i]<- force
                    xx[it/nc,i]<- x
                    vv[it/nc,i]<- v
                    ee[it/nc,i]<- (ji1-jd1)
                    ig<-0
                }
            }
            if(eta==xa)
            {
                eta<-xb
            } else {
                eta<-xa
            }
        }
        
    }
    final<-list("posicao"=xx,"velocidade"=vv, "energia"=ee, "massa"=massa,"gamma"=gamma,"k"=ctelastica, "alfa" = alfazero,"dalfa"= dalfa, "xa"=xa, "dxa"=dxa, "ns" = ns ,"nt"=nt,"nc"=nc,"dt"= dt, "temp"=Temp, "tkick" = tt, "kick" = deltaf, "tempo"= tempo) #"forca"= ff
    
    return(final)
}


SimulacaoBimodalKickFrandomSameNoise <- function(ns, nt, nc, dt, massa, gamma, ctelastica, alfa, dalfa, xa, dxa, deltaf, tt){
    
    
    ntc<-nt*nc
    
    xx<-matrix( nrow = nt, ncol = ns)
    vv<-matrix( nrow = nt, ncol = ns)
    ee<-matrix( nrow = nt, ncol = ns)
    #ff<-matrix( nrow = nt, ncol = ns)
    kick<-numeric(ns)
    
    
    alfazero<-alfa
    
    ua<- alfa/2
    
    ub<- ua
    
    c1<- ctelastica/massa
    
    gm<- gamma/massa
    
    ie<-0
    
    it<-0
    
    tempo<-seq(0,ntc*dt, dt*nc)
    
    for(i in 1:ns)
    {
        
        print(i)
        
        Temp<- ((xa^2)*(massa * alfa))/(gamma * (ctelastica + alfa * (massa * alfa + gamma)  ))
        
        x<-0
        v<-0
        v12<-0
        ji1<-0
        jd1<-0
        alfa<- alfa+ dalfa
        xa<- xa+dxa
        
        xb<- -xa
        ua<-alfa/2
        ub<-ua
        eta<-xa
        it<-0
        ic<-0
        ir<-0
        ig<-0
        ie<-0
        #### será?
        ii<-0
        indicadorkick<-0
        
        while(it < ntc)
        {
            
            #print(it*dt)
            

                it<-it+1
                ie<-ie+1
                
                u<-v
                y<-x
                force<-  -gm*u -c1*y + ruido[it]/massa
                
                
                
                if( it*dt >= tt & indicadorkick==0)
                {
                    kick<- deltaf*runif(1,min =-1,max=1)
                    force <- force + kick
                    #print(i)
                    #indicadorkick<- 1
                }
                
                v<- u + force*dt
                x<- y+ v*dt
                
                ji1<-ji1+ (v+u)*eta*dt/2
                #jd1<-jd1+ gamma*(v^2+u^2)*dt/2
                jd1 <- jd1 + gamma*(((v+u)/2)^2)*dt
                
                ig<-ig+1
                
                if(ig==nc)
                {
                    #ff[it/nc,i]<- force
                    xx[it/nc,i]<- x
                    vv[it/nc,i]<- v
                    ee[it/nc,i]<- (ji1-jd1)
                    ig<-0
                }
            }

        }
    
    final<-list("posicao"=xx,"velocidade"=vv, "energia"=ee, "massa"=massa,"gamma"=gamma,"k"=ctelastica, "alfa" = alfazero,"dalfa"= dalfa, "xa"=xa, "dxa"=dxa, "ns" = ns ,"nt"=nt,"nc"=nc,"dt"= dt, "temp"=Temp, "tkick" = tt, "kick" = deltaf, "tempo"= tempo) #"forca"= ff
    
    return(final)
}



SimulacaoBimodalKickH <- function(ns, nt, nc, dt, massa, gamma, ctelastica, alfa, dalfa, xa, dxa, deltaf, tt){
    
    
    ntc<-nt*nc
    
    xx<-matrix( nrow = nt, ncol = ns)
    vv<-matrix( nrow = nt, ncol = ns)
    ee<-matrix( nrow = nt, ncol = ns)
    #ff<-matrix( nrow = nt, ncol = ns)
    
    alfazero<-alfa
    
    ua<- alfa/2
    
    ub<- ua
    
    c1<- ctelastica/massa
    
    gm<- gamma/massa
    
    ie<-0
    
    it<-0
    
    tempo<-seq(0,ntc*dt, dt*nc)
    
    for(i in 1:ns)
    {
        
        print(i)
        
        Temp<- ((xa^2)*(massa * alfa))/(gamma * (ctelastica + alfa * (massa * alfa + gamma)  ))
        
        x<-0
        v<-0
        v12<-0
        ji1<-0
        jd1<-0
        alfa<- alfa+ dalfa
        xa<- xa+dxa
        
        xb<- -xa
        ua<-alfa/2
        ub<-ua
        eta<-xa
        it<-0
        ic<-0
        ir<-0
        ig<-0
        ie<-0
        indicadorkick<-0
        
        while(it < ntc)
        {
            
            #print(it*dt)
            
            
            
            if(eta==xa)
            {
                teta<- - log(1-runif(1))/ub
                neta<- teta%/%dt + 1
            } else {
                teta<- - log(1-runif(1))/ua
                neta<- teta%/%dt + 1
            }
            ie<-0
            while(ie<neta & it < ntc)
            {
                it<-it+1
                ie<-ie+1
                
                u<-v
                y<-x
                force<-  -gm*u -c1*y + eta/massa
                
                
                if( it*dt >= tt & indicadorkick==0)
                {
                    force <- force + deltaf*x
                    print(i)
                    indicadorkick<- 1
                }
                
                v<- u + force*dt
                x<- y+ v*dt
                
                ji1<-ji1+ (v+u)*eta*dt/2
                #jd1<-jd1+ gamma*(v^2+u^2)*dt/2
                jd1 <- jd1 + gamma*(((v+u)/2)^2)*dt
                
                ig<-ig+1
                
                if(ig==nc)
                {
                    #ff[it/nc,i]<- force
                    xx[it/nc,i]<- x
                    vv[it/nc,i]<- v
                    ee[it/nc,i]<- (ji1-jd1)
                    ig<-0
                }
            }
            if(eta==xa)
            {
                eta<-xb
            } else {
                eta<-xa
            }
        }
        
    }
    final<-list("posicao"=xx,"velocidade"=vv, "energia"=ee, "massa"=massa,"gamma"=gamma,"k"=ctelastica, "alfa" = alfazero,"dalfa"= dalfa, "xa"=xa, "dxa"=dxa, "ns" = ns ,"nt"=nt,"nc"=nc,"dt"= dt, "temp"=Temp, "tkick" = tt, "kick" = deltaf, "tempo"= tempo) #"forca"= ff
    
    return(final)
}



SimulacaoBimodalKickV <- function(ns, nt, nc, dt, massa, gamma, ctelastica, alfa, dalfa, xa, dxa, deltav, tt){
    
    
    ntc<-nt*nc
    
    xx<-matrix( nrow = nt, ncol = ns)
    vv<-matrix( nrow = nt, ncol = ns)
    ee<-matrix( nrow = nt, ncol = ns)
    
    alfazero<-alfa
    
    ua<- alfa/2
    
    ub<- ua
    
    c1<- ctelastica/massa
    
    gm<- gamma/massa
    
    ie<-0
    
    it<-0
    
    
    
    for(i in 1:ns)
    {
        print(i)
        
        Temp<- ((xa^2)*(massa * alfa))/(gamma * (ctelastica + alfa * (massa * alfa + gamma)  ))
        
        x<-0
        v<-0
        v12<-0
        ji1<-0
        jd1<-0
        alfa<- alfa + dalfa
        xa<- xa + dxa
        
        xb<- -xa
        ua<-alfa/2
        ub<-ua
        eta<-xa
        it<-0
        ic<-0
        ir<-0
        ig<-0
        ie<-0
        indicadorkick<-0
        
        
        while(it < ntc)
        {
            
            if( it*dt >= tt & indicadorkick==0)
            {
                v<- v + deltav
                #x<-  ((x**2 - 2* deltav/ctelastica))
                #print('fooooooi')
                indicadorkick<- 1
            }
            
            if(eta==xa)
            {
                teta<- - log(1-runif(1))/ub
                neta<- teta%/%dt + 1
            } else {
                teta<- - log(1-runif(1))/ua
                neta<- teta%/%dt + 1
            }
            
            ie<-0
            
            while(ie<neta & it < ntc)
            {
                
                
                it<-it+1
                ie<-ie+1
                
                u<-v
                y<-x
                force<-  -gm*u -c1*y + eta/massa
                v<- u + force*dt
                x<- y+ v*dt
                
                ji1<-ji1+ (v+u)*eta*dt/2
                #jd1<-jd1+ gamma*(v^2+u^2)*dt/2
                jd1 <- jd1 + gamma*(((v+u)/2)^2)*dt
                
                ig<-ig+1
                
                if(ig==nc)
                {
                    xx[it/nc,i]<- x
                    vv[it/nc,i]<- v
                    ee[it/nc,i]<- (ji1-jd1)
                    ig<-0
                }
            }
            if(eta==xa)
            {
                eta<-xb
            } else {
                eta<-xa
            }
        }
        
    }
    final<-list("posicao"=xx,"velocidade"=vv, "energia"=ee, "massa"=massa,"gamma"=gamma,"k"=ctelastica, "alfa" = alfazero,"dalfa"= dalfa, "xa"=xa, "dxa"=dxa, "ns" = ns ,"nt"=nt,"nc"=nc,"dt"= dt, "temp"=Temp, "tkick" = tt, "kick" = deltav )
    
    return(final)
}




SimulacaoBimodal <- function(ns, nt, nc, dt, massa, gamma, ctelastica, alfa, dalfa, xa, dxa){

tempo<-seq(dt*nc, nt*dt*nc, dt*nc)

ntc<-nt*nc

xx<-matrix( nrow = nt, ncol = ns)
vv<-matrix( nrow = nt, ncol = ns)
ee<-matrix( nrow = nt, ncol = ns)


ji<-matrix( nrow = nt, ncol = ns)
jd<-matrix( nrow = nt, ncol = ns)
ruido<-matrix( nrow = nt, ncol = ns)

alfazero<-alfa

ua<- alfa/2

ub<- ua

c1<- ctelastica/massa

gm<- gamma/massa

ie<- 0
it<- 0
 
for(i in 1:ns)
{
    
    
    Temp<- ((xa^2)*(massa * alfa))/(gamma * (ctelastica + alfa * (massa * alfa + gamma)  ))
    
	x<-0
	v<-0
	v12<-0
	ji1<-0
	jd1<-0
	alfa<- alfa+ dalfa
    xa<- xa+dxa
	
	xb<- -xa
	ua<-alfa/2
	ub<-ua
	eta<-xa
	it<-0
	ic<-0
	ir<-0
	ig<-0
	ie<-0
	
	while(it < ntc)
	{
		if(eta==xa)
		{
			teta<- - log(1-runif(1))/ub
			neta<- teta%/%dt + 1
		} else {
			teta<- - log(1-runif(1))/ua
			neta<- teta%/%dt + 1
		}
		ie<-0
		while(ie<neta & it < ntc) 
		{
		it<-it+1	
		ie<-ie+1

		u<-v
		y<-x
		force<-  -gm*u -c1*y + eta/massa
		v<- u + force*dt
		x<- y+ v*dt
		
        
        #  u<-v
        #y<-x
        #force<-  -gm*u -c1*y + eta/massa
        #v<- u + force*dt
        #x<- y+ v*dt
        
        
        
        
		ji1 <- ji1 + (v+u)*eta*dt/2
        #jd1<-jd1+ gamma*(v^2+u^2)*dt/2
        jd1 <- jd1 + gamma*(((v+u)/2)^2)*dt   #!!!!!
        
        
        
		ig<-ig+1
		
		if(ig==nc)
			{
                
				xx[it/nc,i]<- x
				vv[it/nc,i]<- v
			 	ee[it/nc,i]<- (ji1-jd1)
                ji[it/nc,i]<- ji1
                jd[it/nc,i]<- jd1
                ruido[it/nc,i]<- eta
                
                
				ig<-0
			}
		}
		if(eta==xa)
			{
			eta<-xb
			} else {
			eta<-xa
			}
	}
	
}

final<-list( "posicao"=xx, "velocidade"=vv, "ji" = ji , "jd" = jd ,"energia"=ee, "massa"=massa,"gamma"=gamma,"k"=ctelastica, "alfa" = alfazero,"dalfa"= dalfa, "xa"=xa, "dxa"=dxa, "ns" = ns ,"nt"=nt,"nc"=nc,"dt"= dt, "temp"=Temp, "ruido" = ruido, "t" = tempo)

return(final)
}



TesteInternal <- function(){
    
    qero1<-SimulacaoBimodalInternal(1000,4000,100,0.001,1,1,1,1,0,1,0)
    
    wat1x<-funcaodecorrelacaoSample(qero1$posicao,3900,100,400)
    wat1v<-funcaodecorrelacaoSample(qero1$velocidade,3900,100,400)
    
    Salva(wat1x,'corr1x.csv')
    Salva(wat1v,'corr1v.csv')
    
    
    qero2k<-SimulacaoBimodalInternal(1000,4000,100,0.001,1,1,2,1,0,1,0)
    
    wat2kx<-funcaodecorrelacaoSample(qero2k$posicao,3900,100,400)
    wat2kv<-funcaodecorrelacaoSample(qero2k$velocidade,3900,100,400)
    
    Salva(wat2kx,'corr2kx.csv')
    Salva(wat2kv,'corr2kv.csv')
    
    
    qero2<-SimulacaoBimodalInternal(1000,4000,100,0.001,1,2,2,2,0,1,0)
    
    wat2x<-funcaodecorrelacaoSample(qero2$posicao,3900,100,400)
    wat2v<-funcaodecorrelacaoSample(qero2$velocidade,3900,100,400)
    
    Salva(wat2x,'corr2x.csv')
    Salva(wat2v,'corr2v.csv')
    
    
    qero3<-SimulacaoBimodalInternal(1000,4000,100,0.001,1,3,3,3,0,1,0)
    
    wat3x<-funcaodecorrelacaoSample(qero3$posicao,3900,100,400)
    wat3v<-funcaodecorrelacaoSample(qero3$velocidade,3900,100,400)
    
    Salva(wat3x,'corr3x.csv')
    Salva(wat3v,'corr3v.csv')
    
    
    qerob<-SimulacaoBimodalInternal(1000,4000,100,0.001,1,2,3,4,0,1,0)
    
    watbx<-funcaodecorrelacaoSample(qerob$posicao,3900,100,400)
    watbv<-funcaodecorrelacaoSample(qerob$velocidade,3900,100,400)
    
    Salva(watbx,'corrbx.csv')
    Salva(watbv,'corrbv.csv')
    
    
    
}

TesteInternal2 <- function(){
    
    qero1<-SimulacaoBimodalInternal(1000,4000,100,0.001,2,1,1,1,0,4,0)
    
    wat1x<-funcaodecorrelacaoSample(qero1$posicao,3900,100,400)
    wat1v<-funcaodecorrelacaoSample(qero1$velocidade,3900,100,400)
    
    Salva(wat1x,'corrax.csv')
    Salva(wat1v,'corrav.csv')
    
    
    
    qero2<-SimulacaoBimodalInternal(1000,4000,100,0.001,3,2,2,2,0,7,0)
    
    wat2x<-funcaodecorrelacaoSample(qero2$posicao,3900,100,400)
    wat2v<-funcaodecorrelacaoSample(qero2$velocidade,3900,100,400)
    
    Salva(wat2x,'corrbx.csv')
    Salva(wat2v,'corrbv.csv')
    
    
    qero3<-SimulacaoBimodalInternal(1000,4000,100,0.001,8,3,3,3,0,9,0)
    
    wat3x<-funcaodecorrelacaoSample(qero3$posicao,3900,100,400)
    wat3v<-funcaodecorrelacaoSample(qero3$velocidade,3900,100,400)
    
    Salva(wat3x,'corrcx.csv')
    Salva(wat3v,'corrcv.csv')
    
    
    
    
    
}

TesteHistogramaInternal<- function()
{
    qero1<-SimulacaoBimodalInternal(1000,4000,100,0.001,2,1,1,1,0,4,0)
    
    #hist(
}


SimulacaoBimodalInternal <- function(ns, nt, nc, dt, massa, gamma, ctelastica, alfa, dalfa, xa, dxa){
    
    tempo<-seq(dt*nc, nt*dt*nc, dt*nc)
    
    ntc<-nt*nc
    
    xx<-matrix( nrow = nt, ncol = ns)
    vv<-matrix( nrow = nt, ncol = ns)
    ee<-matrix( nrow = nt, ncol = ns)
    
    #velocoarse<-numeric(nc*nt)
    
    ji<-matrix( nrow = nt, ncol = ns)
    jd<-matrix( nrow = nt, ncol = ns)
    ruido<-matrix( nrow = nt, ncol = ns)
    
    alfazero<-alfa
    
    ua<- alfa/2
    
    ub<- ua
    
    c1<- ctelastica/massa
    
    gm<- gamma/massa
    
    ie<- 0
    it<- 0
    
    for(i in 1:ns)
    {
        #
        print(i)
        
        
        Temp<- ((xa^2)*(massa * alfa))/(gamma * (ctelastica + alfa * (massa * alfa + gamma)  ))
        
        x<-0
        v<-0
        v12<-0
        ji1<-0
        jd1<-0
        alfa<- alfa+ dalfa
        xa<- xa+dxa
        
        xb<- -xa
        ua<-alfa/2
        ub<-ua
        eta<-xa
        it<-0
        ic<-0
        ir<-0
        ig<-0
        ie<-0
        
        force2<-0
        
        while(it < ntc)
        {
            
            if(eta==xa)
            {
                teta<- - log(1-runif(1))/ub
                neta<- teta%/%dt + 1
            } else {
                teta<- - log(1-runif(1))/ua
                neta<- teta%/%dt + 1
            }
            
            ie<- 0
            
            #force2<-0
            
            while(ie<neta & it < ntc)
            {
                it<-it+1
                ie<-ie+1
                
                #velocoarse[it]<-v
                
                u<-v
                y<-x
                
                # verificar FDT
                force1<-  -c1*y + eta/massa
                
                #force1<-  -c1*y + (eta*sqrt(gamma))/massa
                
                
                # for(ij in 1:it)
                #{
                
                #force2<- force2 * exp(- alfa * dt  )  - gm * u * exp(- alfa * dt  ) * dt
                
                
                force<- force1 + force2
                
                v<- u + force*dt
                
                x<- y+ v*dt
                
                force2<- force2 * exp(- alfa * dt  )  - gm * (v+u) * exp(- alfa * dt  ) * dt/(2* alfa)
                
                # force2<- force2 * exp(- alfa * dt  )  - ( xa**2) *(v+u) * exp(- alfa * dt  ) * dt/(2* alfa)
                
                
                ji1 <- ji1 + (v+u)*eta*dt/2
                
                #jd1<-jd1+ gamma*(v^2+u^2)*dt/2
                
                jd1 <- jd1 + gamma*(((v+u)/2)^2)*dt   #!!!!!
                
                
                
                ig<-ig+1
                
                if(ig==nc)
                {
                    
                    xx[it/nc,i]<- x
                    vv[it/nc,i]<- v
                    ee[it/nc,i]<- (ji1-jd1)
                    ji[it/nc,i]<- ji1
                    jd[it/nc,i]<- jd1
                    ruido[it/nc,i]<- eta
                    
                    
                    ig<-0
                }
            }
            if(eta==xa)
            {
                eta<-xb
            } else {
                eta<-xa
            }
        }
        
    }
    
    final<-list( "posicao"=xx, "velocidade"=vv, "ji" = ji , "jd" = jd ,"energia"=ee, "massa"=massa,"gamma"=gamma,"k"=ctelastica, "alfa" = alfazero,"dalfa"= dalfa, "xa"=xa, "dxa"=dxa, "ns" = ns ,"nt"=nt,"nc"=nc,"dt"= dt, "temp"=Temp, "ruido" = ruido, "t" = tempo)
    
    return(final)
}


HistogramaOD<- function(){
    
    qero<-SimulacaoBimodalODlite(100000,1000,100,0.001,1,1,1,0.02,0,1,0)
    hx<-hist(qero$posicao, breaks=200)
    hv<-hist(qero$velocidade,breaks=200)
    SalvaCorrelacao(hx$density,"002ODhistposicao.csv")
    SalvaCorrelacao(hx$mids,"002ODmidposicao.csv")
    SalvaCorrelacao(hv$density,"002ODhistvelo.csv")
    SalvaCorrelacao(hv$mids,"002ODmidvelo.csv")
    
}

SimulacaoBimodalODlite <- function(ns, nt, nc, dt, massa, gamma, ctelastica, alfa, dalfa, xa, dxa){
    
    tempo<-seq(dt*nc, nt*dt*nc, dt*nc)
    
    ntc<-nt*nc
    
    xx<-matrix( nrow = nt, ncol = ns)
    vv<-matrix( nrow = nt, ncol = ns)
    # ee<-matrix( nrow = nt, ncol = ns)
    
    # ji<-matrix( nrow = nt, ncol = ns)
    #jd<-matrix( nrow = nt, ncol = ns)
    #ruido<-matrix( nrow = nt, ncol = ns)
    
    alfazero<-alfa
    
    ua<- alfa/2
    
    ub<- ua
    
    c1<- ctelastica
    
    gm<- gamma
    
    ie<- 0
    
    it<- 0
    
    for(i in 1:ns)
    {
        print(i)
        
        #  Temp<- ((xa^2)*(massa * alfa))/(gamma * (ctelastica + alfa * (massa * alfa + gamma)  ))
        
        x<- 0
        v<- 0
        v12<- 0
        ji1<- 0
        jd1<- 0
        alfa<- alfa + dalfa
        xa<- xa + dxa
        
        xb<- -xa
        ua<-alfa/2
        ub<-ua
        eta<-xa
        
        it<- 0
        ic<- 0
        ir<- 0
        ig<- 0
        ie<- 0
        
        while(it < ntc)
        {
            if(eta==xa)
            {
                teta<- - log(1-runif(1))/ub
                neta<- teta%/%dt + 1
            } else {
                teta<- - log(1-runif(1))/ua
                neta<- teta%/%dt + 1
            }
            ie<-0
            while(ie<neta & it < ntc)
            {
                it<-it+1
                ie<-ie+1
                
                u<-v
                y<-x
                
                dx<- - ( c1 * y  + eta)*(dt / gm)
                
                x <-  y + dx
                
                v <- dx/dt
                
                
                #ji1 <- ji1 - (v+u)*eta*dt/2 -------> original era com sinal + mas troquei pra acertar (????)
                
                #ji1 <- ji1 - (v+u)*eta*dt/2

                #jd1 <- jd1 + gamma*(((v+u)/2)^2)*dt
                
                ig <- ig+1
                
                if(ig==nc)
                {
                    
                    xx[it/nc,i]<- x
                    vv[it/nc,i]<- v
                    #  ee[it/nc,i]<- (ji1-jd1)
                    #ji[it/nc,i]<- ji1
                    #jd[it/nc,i]<- jd1
                    #ruido[it/nc,i]<- eta
                    
                    ig<-0
                }
            }
            if(eta==xa)
            {
                eta<-xb
            } else {
                eta<-xa
            }
        }
        
    }
    
    #final<-list( "posicao"=xx, "velocidade"=vv, "ji" = ji , "jd" = jd ,"energia"=ee, "massa"=massa,"gamma"=gamma,"k"=ctelastica, "alfa" = alfazero,"dalfa"= dalfa, "xa"=xa, "dxa"=dxa, "ns" = ns ,"nt"=nt,"nc"=nc,"dt"= dt, "temp"=Temp, "ruido" = ruido, "t" = tempo)
    
    final<-list( "posicao"=xx, "velocidade"=vv,  "massa"=massa,"gamma"=gamma,"k"=ctelastica, "alfa" = alfazero,"dalfa"= dalfa, "xa"=xa, "dxa"=dxa, "ns" = ns ,"nt"=nt,"nc"=nc,"dt"= dt,  "t" = tempo)
    
    return(final)
}

SimulacaoBimodalODKickF <- function(ns, nt, nc, dt, massa, gamma, ctelastica, alfa, dalfa, xa, dxa,  deltaf, tt){
    
    tempo<-seq(dt*nc, nt*dt*nc, dt*nc)
    
    ntc<-nt*nc
    
    xx<-matrix( nrow = nt, ncol = ns)
    vv<-matrix( nrow = nt, ncol = ns)
    ee<-matrix( nrow = nt, ncol = ns)
    
    ji<-matrix( nrow = nt, ncol = ns)
    jd<-matrix( nrow = nt, ncol = ns)
    #ruido<-matrix( nrow = nt, ncol = ns)
    
    alfazero<-alfa
    
    ua<- alfa/2
    
    ub<- ua
    
    c1<- ctelastica
    
    gm<- gamma
    
    ie<- 0
    
    it<- 0
    
    for(i in 1:ns)
    {
        print(i)
        
        Temp<- ((xa^2)*(massa * alfa))/(gamma * (ctelastica + alfa * (massa * alfa + gamma)  ))
        
        x<- 0
        v<- 0
        v12<- 0
        ji1<- 0
        jd1<- 0
        alfa<- alfa + dalfa
        xa<- xa + dxa
        
        xb<- -xa
        ua<-alfa/2
        ub<-ua
        eta<-xa
        
        it<- 0
        ic<- 0
        ir<- 0
        ig<- 0
        ie<- 0
        
        indicadorkick<-0
        
        while(it < ntc)
        {
            if(eta==xa)
            {
                teta<- - log(1-runif(1))/ub
                neta<- teta%/%dt + 1
            } else {
                teta<- - log(1-runif(1))/ua
                neta<- teta%/%dt + 1
            }
            ie<-0
            while(ie<neta & it < ntc)
            {
                it<-it+1
                ie<-ie+1
                
                u<-v
                y<-x
                
                force<- - (c1 * y  + eta)
                
                if( it*dt >= tt & indicadorkick==0)
                {
                    force <- force + deltaf
                    #print(i)
                    #indicadorkick<- 1
                }
                
                dx<- force*(dt / gm)
                
                x <-  y + dx
                
                v <- dx/dt
                
                
                #ji1 <- ji1 - (v+u)*eta*dt/2 -------> original era com sinal + mas troquei pra acertar (????)
                
                ji1 <- ji1 - (v+u)*eta*dt/2
                
                jd1 <- jd1 + ctelastica*(((v+u)/2)^2)*dt
                
                ig <- ig+1
                
                if(ig==nc)
                {
                    
                    xx[it/nc,i]<- x
                    vv[it/nc,i]<- v
                    ee[it/nc,i]<- (ji1-jd1)
                    ji[it/nc,i]<- ji1
                    jd[it/nc,i]<- jd1
                    #ruido[it/nc,i]<- eta
                    
                    ig<-0
                }
            }
            if(eta==xa)
            {
                eta<-xb
            } else {
                eta<-xa
            }
        }
        
    }
    
    final<-list( "posicao"=xx, "velocidade"=vv, "ji" = ji , "jd" = jd ,"energia"=ee, "massa"=massa,"gamma"=gamma,"k"=ctelastica, "alfa" = alfazero,"dalfa"= dalfa, "xa"=xa, "dxa"=dxa, "ns" = ns ,"nt"=nt,"nc"=nc,"dt"= dt, "temp"=Temp,  "t" = tempo) # "ruido" = ruido,
    
    return(final)
}

SimulacaoBimodalOD <- function(ns, nt, nc, dt, massa, gamma, ctelastica, alfa, dalfa, xa, dxa){
    
    tempo<-seq(dt*nc, nt*dt*nc, dt*nc)
    
    ntc<-nt*nc
    
    xx<-matrix( nrow = nt, ncol = ns)
    vv<-matrix( nrow = nt, ncol = ns)
    ee<-matrix( nrow = nt, ncol = ns)
    
    ji<-matrix( nrow = nt, ncol = ns)
    jd<-matrix( nrow = nt, ncol = ns)
    #ruido<-matrix( nrow = nt, ncol = ns)
    
    alfazero<-alfa
    
    ua<- alfa/2
    
    ub<- ua
    
    c1<- ctelastica
    
    gm<- gamma
    
    ie<- 0
    
    it<- 0
    
    for(i in 1:ns)
    {
        print(i)
        
        Temp<- ((xa^2)*(massa * alfa))/(gamma * (ctelastica + alfa * (massa * alfa + gamma)  ))
        
        x<- 0
        v<- 0
        v12<- 0
        ji1<- 0
        jd1<- 0
        alfa<- alfa + dalfa
        xa<- xa + dxa
        
        xb<- -xa
        ua<-alfa/2
        ub<-ua
        eta<-xa
        
        it<- 0
        ic<- 0
        ir<- 0
        ig<- 0
        ie<- 0
        
        while(it < ntc)
        {
            if(eta==xa)
            {
                teta<- - log(1-runif(1))/ub
                neta<- teta%/%dt + 1
            } else {
                teta<- - log(1-runif(1))/ua
                neta<- teta%/%dt + 1
            }
            ie<-0
            while(ie<neta & it < ntc)
            {
                it<-it+1
                ie<-ie+1
                
                u<-v
                y<-x
                
                dx<- - ( c1 * y  + eta)*(dt / gm)
                
                x <-  y + dx
                
                v <- dx/dt
                
                
                #ji1 <- ji1 - (v+u)*eta*dt/2 -------> original era com sinal + mas troquei pra acertar (????)
                
                ji1 <- ji1 - (v+u)*eta*dt/2
                
                jd1 <- jd1 + gamma*(((v+u)/2)^2)*dt
                
                #jd1 <- jd1 + ctelastica *(((v+u)/2)^2)*dt
                
                ig <- ig + 1
                
                if(ig==nc)
                {
                    
                    xx[it/nc,i]<- x
                    vv[it/nc,i]<- v
                    ee[it/nc,i]<- (ji1-jd1)
                    ji[it/nc,i]<- ji1
                    jd[it/nc,i]<- jd1
                    #ruido[it/nc,i]<- eta
                    
                    ig<-0
                }
            }
            if(eta==xa)
            {
                eta<-xb
            } else {
                eta<-xa
            }
        }
        
    }
    
    final<-list( "posicao"=xx, "velocidade"=vv, "ji" = ji , "jd" = jd ,"energia"=ee, "massa"=massa,"gamma"=gamma,"k"=ctelastica, "alfa" = alfazero,"dalfa"= dalfa, "xa"=xa, "dxa"=dxa, "ns" = ns ,"nt"=nt,"nc"=nc,"dt"= dt, "temp"=Temp,  "t" = tempo) # "ruido" = ruido,
    
    return(final)
}

SimulacaoLangevinOD <- function(ns, nt, nc, dt, massa, gamma, ctelastica, temp){
    
    tempo<-seq(dt*nc, nt*dt*nc, dt*nc)
    
    ntc<-nt*nc
    
    xx<-matrix( nrow = nt, ncol = ns)
    vv<-matrix( nrow = nt, ncol = ns)
    ee<-matrix( nrow = nt, ncol = ns)
    
    ji<-matrix( nrow = nt, ncol = ns)
    jd<-matrix( nrow = nt, ncol = ns)
    #ruido<-matrix( nrow = nt, ncol = ns)
    
    #alfazero<-alfa
    
    #ua<- alfa/2
    
    #ub<- ua
    
    c1<- ctelastica
    
    gm<- gamma
    
    ie<- 0
    
    it<- 0
    
    for(i in 1:ns)
    {
        print(i)
        
        x<- 0
        v<- 0
        v12<- 0
        ji1<- 0
        jd1<- 0

        
        
        it<- 0
        ic<- 0
        ir<- 0
        ig<- 0
        ie<- 0
        
        while(it < ntc)
        {
            
                it<-it+1
                ie<-ie+1
                
                u<-v
                y<-x
                
                eta<-rnorm(1)
                
                
                forceest <- sqrt(2*gamma*temp*dt)*rnorm(1)
                
   
                
                #forcedet <-  (-gm*u -c1*y)*dt
                
                dx<- - ( c1 * y  )*(dt / gm) - forceest
                
                x <-  y + dx
                
                v <- dx/dt
                
                
                #ji1 <- ji1 - (v+u)*eta*dt/2 -------> original era com sinal + mas troquei pra acertar (????)
                
                #ji1 <- ji1 - (v+u)*eta*dt/2
                
                #jd1 <- jd1 + gamma*(((v+u)/2)^2)*dt
                
                
                ji1 <- ji1 - (v+u)*forceest/2
                jd1 <- jd1 + gamma*(((v+u)/2)^2)*dt
                
                ig <- ig + 1
                
                if(ig==nc)
                {
                    
                    xx[it/nc,i]<- x
                    vv[it/nc,i]<- v
                    ee[it/nc,i]<- (ji1-jd1)
                    ji[it/nc,i]<- ji1
                    jd[it/nc,i]<- jd1
                    #ruido[it/nc,i]<- eta
                    
                    ig<-0
                }
      
        }
        
    }
    
    final<-list( "posicao"=xx, "velocidade"=vv, "ji" = ji , "jd" = jd ,"energia"=ee, "massa"=massa,"gamma"=gamma,"k"=ctelastica,  "ns" = ns ,"nt"=nt,"nc"=nc,"dt"= dt, "temp"=temp,  "t" = tempo) # "ruido" = ruido,
    
    return(final)
}




SimulacaoBimodalDuasParticulas <- function(ns, nt, nc, dt, massa1, gamma1, ctelastica1, alfa1, dalfa1, xa1, dxa1, massa2, gamma2, ctelastica2, alfa2, dalfa2, xa2, dxa2, cacopla ){
    
    
    ntc<-nt*nc
    
    xx1<-matrix( nrow = nt, ncol = ns)
    vv1<-matrix( nrow = nt, ncol = ns)
    ee1<-matrix( nrow = nt, ncol = ns)
    
    xx2<-matrix( nrow = nt, ncol = ns)
    vv2<-matrix( nrow = nt, ncol = ns)
    ee2<-matrix( nrow = nt, ncol = ns)
    
    ########
    
    ji<-matrix( nrow = nt, ncol = ns)
    jd<-matrix( nrow = nt, ncol = ns)
    ruido<-matrix( nrow = nt, ncol = ns)
    
    alfazero<-alfa1
    
    ua<- alfa1/2
    
    ub<- ua
    
    c1<- ctelastica1/massa1
    
    gm1<- gamma1/massa1
    
    alfazero<-alfa2
    
    ua<- alfa2/2
    
    ub<- ua
    
    c2<- ctelastica2/massa2
    
    gm2<- gamma2/massa2
    
    cacop
    
    
    
    
    
    
    ie<-0
    
    it<-0
    
    
    
    for(i in 1:ns)
    {
        
        print(i)
        Temp<- ((xa^2)*(massa * alfa))/(gamma * (ctelastica + alfa * (massa * alfa + gamma)  ))
        
        x<-0
        v<-0
        v12<-0
        ji1<-0
        jd1<-0
        alfa<- alfa+ dalfa
        xa<- xa+dxa
        
        xb<- -xa
        ua<-alfa/2
        ub<-ua
        eta<-xa
        it<-0
        ic<-0
        ir<-0
        ig<-0
        ie<-0
        
        while(it < ntc)
        {
            if(eta==xa)
            {
                teta<- - log(1-runif(1))/ub
                neta<- teta%/%dt + 1
            } else {
                teta<- - log(1-runif(1))/ua
                neta<- teta%/%dt + 1
            }
            ie<-0
            while(ie<neta & it < ntc)
            {
                it<-it+1
                ie<-ie+1
                
                u<-v
                y<-x
                force<-  -gm*u -c1*y + eta/massa  - molaacopla * ( y1-y2 )
                v<- u + force*dt
                x<- y+ v*dt
                
                
                u2<-v2
                y2<-x2
                
                force<-  -gm*u2 -c1*y2 + eta/massa - molaacopla * ( y2 - y1)
                
                v<- u + force*dt
                x<- y+ v*dt
                
                
                
                
                ji1 <- ji1 + (v+u)*eta*dt/2
                #jd1<-jd1+ gamma*(v^2+u^2)*dt/2
                jd1 <- jd1 + gamma*(((v+u)/2)^2)*dt   #!!!!!
                
                
                
                ig<-ig+1
                
                if(ig==nc)
                {
                    
                    xx[it/nc,i]<- x
                    vv[it/nc,i]<- v
                    ee[it/nc,i]<- (ji1-jd1)
                    ji[it/nc,i]<- ji1
                    jd[it/nc,i]<- jd1
                    ruido[it/nc,i]<- eta
                    
                    
                    ig<-0
                }
            }
            if(eta==xa)
            {
                eta<-xb
            } else {
                eta<-xa
            }
        }
        
    }
    
    final<-list( "posicao"=xx, "velocidade"=vv, "ji" = ji , "jd" = jd ,"energia"=ee, "massa"=massa,"gamma"=gamma,"k"=ctelastica, "alfa" = alfazero,"dalfa"= dalfa, "xa"=xa, "dxa"=dxa, "ns" = ns ,"nt"=nt,"nc"=nc,"dt"= dt, "temp"=Temp, "ruido" = ruido)
    
    return(final)
}



SimulacaoBimodalTPC1 <- function(ns, nt, nc, dt, massa, gamma, ctelastica, alfa, dalfa, xa, dxa, var){
    
    
    ntc<-nt*nc
    
    xx<-matrix( nrow = nt, ncol = ns)
    vv<-matrix( nrow = nt, ncol = ns)
    ee<-matrix( nrow = nt, ncol = ns)
    
    alfazero<-alfa
    
    ua<- alfa/2
    
    ub<- ua
    
    c1<- ctelastica/massa
    
    gm<- gamma/massa
    
    ie<-0
    
    it<-0
    
    x00<- numeric(ns)
    
    
    for(i in 1:ns)
    {
        
        print(i)
        Temp<- ((xa^2)*(massa * alfa))/(gamma * (ctelastica + alfa * (massa * alfa + gamma)  ))
        #####
        x<- rnorm(1,0,var)
        #####
        x00[i]<- x
        ####
        v<- 0
        v12<-0
        ji1<-0
        jd1<-0
        alfa<- alfa+ dalfa
        xa<- xa+dxa
        
        xb<- -xa
        ua<-alfa/2
        ub<-ua
        eta<-xa
        it<-0
        ic<-0
        ir<-0
        ig<-0
        ie<-0
        
        while(it < ntc)
        {
            if(eta==xa)
            {
                teta<- - log(1-runif(1))/ub
                neta<- teta%/%dt + 1
            } else {
                teta<- - log(1-runif(1))/ua
                neta<- teta%/%dt + 1
            }
            ie<-0
            while(ie<neta & it < ntc)
            {
                it<-it+1
                ie<-ie+1
                
                u<-v
                y<-x
                force<-  -gm*u -c1*y + eta/massa
                v<- u + force*dt
                x<- y+ v*dt
                
                ji1<-ji1+ (v+u)*eta*dt/2
                #jd1<-jd1+ gamma*(v^2+u^2)*dt/2
                jd1 <- jd1 + gamma*(((v+u)/2)^2)*dt
                
                ig<-ig+1
                
                if(ig==nc)
                {
                    xx[it/nc,i]<- x
                    vv[it/nc,i]<- v
                    ee[it/nc,i]<- (ji1-jd1)
                    ig<-0
                }
            }
            if(eta==xa)
            {
                eta<-xb
            } else {
                eta<-xa
            }
        }
        
    }
    final<-list("posicao"=xx,"velocidade"=vv, "energia"=ee, "massa"=massa,"gamma"=gamma,"k"=ctelastica, "alfa" = alfazero,"dalfa"= dalfa, "xa"=xa, "dxa"=dxa, "ns" = ns ,"nt"=nt,"nc"=nc,"dt"= dt, "temp"=Temp, "x00" = x00 )
    
    return(final)
}

SimulacaoBimodalTPC2 <- function(ns, nt, nc, dt, massa, gamma, ctelastica, alfa, dalfa, xa, dxa, var){
    
    
    ntc<-nt*nc
    
    xx<-matrix( nrow = nt, ncol = ns)
    vv<-matrix( nrow = nt, ncol = ns)
    ee<-matrix( nrow = nt, ncol = ns)
    
    alfazero<-alfa
    
    ua<- alfa/2
    
    ub<- ua
    
    c1<- ctelastica/massa
    
    gm<- gamma/massa
    
    ie<-0
    
    it<-0
    
    v00<- numeric(ns)
    
    
    for(i in 1:ns)
    {
        
        print(i)
        Temp<- ((xa^2)*(massa * alfa))/(gamma * (ctelastica + alfa * (massa * alfa + gamma)  ))
        #####
        v<- rnorm(1,0,var)
        #####
        v00[i]<- v
        ####
        x<- 0
        v12<-0
        ji1<-0
        jd1<-0
        alfa<- alfa+ dalfa
        xa<- xa+dxa
        
        xb<- -xa
        ua<-alfa/2
        ub<-ua
        eta<-xa
        it<-0
        ic<-0
        ir<-0
        ig<-0
        ie<-0
        
        while(it < ntc)
        {
            if(eta==xa)
            {
                teta<- - log(1-runif(1))/ub
                neta<- teta%/%dt + 1
            } else {
                teta<- - log(1-runif(1))/ua
                neta<- teta%/%dt + 1
            }
            ie<-0
            while(ie<neta & it < ntc)
            {
                it<-it+1
                ie<-ie+1
                
                u<-v
                y<-x
                force<-  -gm*u -c1*y + eta/massa
                v<- u + force*dt
                x<- y+ v*dt
                
                ji1<-ji1+ (v+u)*eta*dt/2
                #jd1<-jd1+ gamma*(v^2+u^2)*dt/2
                jd1 <- jd1 + gamma*(((v+u)/2)^2)*dt
                
                ig<-ig+1
                
                if(ig==nc)
                {
                    xx[it/nc,i]<- x
                    vv[it/nc,i]<- v
                    ee[it/nc,i]<- (ji1-jd1)
                    ig<-0
                }
            }
            if(eta==xa)
            {
                eta<-xb
            } else {
                eta<-xa
            }
        }
        
    }
    final<-list("posicao"=xx,"velocidade"=vv, "energia"=ee, "massa"=massa,"gamma"=gamma,"k"=ctelastica, "alfa" = alfazero,"dalfa"= dalfa, "xa"=xa, "dxa"=dxa, "ns" = ns ,"nt"=nt,"nc"=nc,"dt"= dt, "temp"=Temp, "v00" = v00 )
    
    return(final)
}




SimulacaoBimodalProtocolo <- function(ns, nt, nc, dt, massa, gamma, ctelastica, alfa, dalfa, xa, dxa, tau){
    
    
    ntc<-nt*nc
    
    xx<-matrix( nrow = nt, ncol = ns)
    vv<-matrix( nrow = nt, ncol = ns)
    ee<-matrix( nrow = nt, ncol = ns)
    
    alfazero<-alfa
    
    ua<- alfa/2
    
    ub<- ua
    
    c1<- ctelastica/massa
    
    gm<- gamma/massa
    
    ie<-0
    
    it<-0
    
    
    
    for(i in 1:ns)
    {
        
        
        Temp<- ((xa^2)*(massa * alfa))/(gamma * (ctelastica + alfa * (massa * alfa + gamma)  ))
        
        x<-0
        v<-0
        v12<-0
        ji1<-0
        jd1<-0
        alfa<- alfa+ dalfa
        xa<- xa+dxa
        
        xb<- -xa
        ua<-alfa/2
        ub<-ua
        eta<-xa
        it<-0
        ic<-0
        ir<-0
        ig<-0
        ie<-0
        
        while(it < ntc)
        {
            if(eta==xa)
            {
                teta<- - log(1-runif(1))/ub
                neta<- teta%/%dt + 1
            } else {
                teta<- - log(1-runif(1))/ua
                neta<- teta%/%dt + 1
            }
            ie<-0
            while(ie<neta & it < ntc)
            {
                it<-it+1
                ie<-ie+1
                
                u<-v
                y<-x
                force<-  -gm*u -c1*y + eta/massa
                v<- u + force*dt
                x<- y+ v*dt
                
                ji1<-ji1+ (v+u)*eta*dt/2
                #jd1<-jd1+ gamma*(v^2+u^2)*dt/2
                jd1 <- jd1 + gamma*(((v+u)/2)^2)*dt
                
                ig<-ig+1
                
                if(ig==nc)
                {
                    xx[it/nc,i]<- x
                    vv[it/nc,i]<- v
                    ee[it/nc,i]<- (ji1-jd1)
                    ig<-0
                }
            }
            if(eta==xa)
            {
                eta<-xb
            } else {
                eta<-xa
            }
        }
        
    }
    final<-list("posicao"=xx,"velocidade"=vv, "energia"=ee, "massa"=massa,"gamma"=gamma,"k"=ctelastica, "alfa" = alfazero,"dalfa"= dalfa, "xa"=xa, "dxa"=dxa, "ns" = ns ,"nt"=nt,"nc"=nc,"dt"= dt, "temp"=Temp )
    
    return(final)
}


SimulacaoBimodalQuadrada <- function(ns, nt, nc, dt, massa, gamma, ctelastica, alfa, dalfa, xa, dxa){
    
    
    ntc<-nt*nc
    
    xx<-matrix( nrow = nt, ncol = ns)
    vv<-matrix( nrow = nt, ncol = ns)
    ee<-matrix( nrow = nt, ncol = ns)
    
    alfazero<-alfa
    
    ua<- alfa/2
    
    ub<- ua
    
    c1<- ctelastica/massa
    
    gm<- gamma/massa
    
    ie<-0
    
    it<-0
    
    alfa0<-alfa
    xa0<-xa
    
    kkk<-0
    
    for(i in 1:ns)
    {
        
        
        
        
        Temp<- ((xa^2)*(massa * alfa))/(gamma * (ctelastica + alfa * (massa * alfa + gamma)  ))
        
        x<-0
        v<-0
        v12<-0
        ji1<-0
        jd1<-0
        
        
        kkk<-kkk+1
        
        if(kkk<sqrt(ns))
        {
            xa<- xa+dxa
        
        }
        
        else
        {
            
            xa<-xa0
            alfa<-alfa+dalfa
            kkk<-0
        }
        
        xb<- -xa
        ua<-alfa/2
        ub<-ua
        eta<-xa
        it<-0
        ic<-0
        ir<-0
        ig<-0
        ie<-0
        
        while(it < ntc)
        {
            if(eta==xa)
            {
                teta<- - log(1-runif(1))/ub
                neta<- teta%/%dt + 1
            } else {
                teta<- - log(1-runif(1))/ua
                neta<- teta%/%dt + 1
            }
            ie<-0
            while(ie<neta & it < ntc)
            {
                it<-it+1
                ie<-ie+1
                
                u<-v
                y<-x
                force<-  -gm*u -c1*y + eta/massa
                v<- u + force*dt
                x<- y+ v*dt
                
                ji1<-ji1+ (v+u)*eta*dt/2 
                #jd1<-jd1+ gamma*(v^2+u^2)*dt/2
                jd1 <- jd1 + gamma*(((v+u)/2)^2)*dt
                
                ig<-ig+1
                
                if(ig==nc)
                {
                    xx[it/nc,i]<- x
                    vv[it/nc,i]<- v
                    ee[it/nc,i]<- (ji1-jd1)
                    ig<-0
                }
            }
            if(eta==xa)
			{
                eta<-xb
			} else {
                eta<-xa
			}
        }
        
    }
    final<-list("posicao"=xx,"velocidade"=vv, "energia"=ee, "massa"=massa,"gamma"=gamma,"k"=ctelastica, "alfa" = alfazero,"dalfa"= dalfa, "xa"=xa, "dxa"=dxa, "ns" = ns ,"nt"=nt,"nc"=nc,"dt"= dt, "temp"=Temp )
    
    return(final)
}

PrinterQuadrado<- function(listasim)
{
    
    qeromais<-analiseentropicaSIMPLIS(listasim)
    
    for(i in 1:listasim$ns)
    {
    
        titulo<-paste(toString(i), ".csv")

        
        experimento<-data.frame(x,y,ss)
        
        write.table(experimento, "entropia.csv" , sep=",", row.names=FALSE, col.names=FALSE)
        

        
    }
}


analiseentropicaTcte <- function(listasim ){
    
    
if(listasim$dalfa==0)
{
    alfai<-rep(listasim$alfa, listasim$ns)
} else {
    alfai<-seq(listasim$alfa+listasim$dalfa,listasim$alfa+listasim$ns*listasim$dalfa,listasim$dalfa)
}

#MANTTENDO A TEMPERATURA CONSTANTE

xa<- sqrt( listasim$T*(listasim$gamma*(listasim$k+alfai*(listasim$gamma+listasim$massa*alfai)))/(listasim$massa*alfai)   )

mediaE<-numeric(listasim$ns)
for(i in 1:listasim$ns){
	mediaE[i]<-mean(listasim$energia[,i])
}


# histogramas
stoth<-numeric(listasim$ns)
for(i in 1:listasim$ns){
    
    
    
hhx<-hist(listasim$posicao[,i], breaks=100)
hhv<-hist(listasim$velocidade[,i], breaks=100)





ssx<- - hhx$density*log(hhx$density)
ssx[which(is.na(ssx))]<- 0
sx<-sum(ssx)*(1/sum(hhx$density))
ssv<- - hhv$density*log(hhv$density)
ssv[which(is.na(ssv))]<- 0
sv<-sum(ssv)*(1/sum(hhv$density))

stoth[i]<- sx + sv


}


#ENERGIA MÉDIA ANALÍTICA
anaE<- (((xa)^2 )*(2*listasim$massa*alfai +listasim$gamma))/(2*listasim$gamma*(listasim$k +listasim$massa *( alfai^2) + alfai * listasim$gamma))

# APROXIMACAO LINEAR PRA ENERGIA MÉDIA E ENTROPIA CHEGANDO NO MARCONI
#Elm<-lm(mediaE~poly(alfai,4,raw=TRUE))
# A<-Elm$coefficients[1]+alfai*Elm$coefficients[2]+(alfai^2)*Elm$coefficients[3]+(alfai^3)*Elm$coefficients[4]+(alfai^4)*Elm$coefficients[5]

#Slm<-lm(stoth~poly(alfai,5,raw=TRUE))
# B<-Slm$coefficients[1]+alfai*Slm$coefficients[2]+(alfai^2)*Slm$coefficients[3]+(alfai^3)*Slm$coefficients[4]+(alfai^4)*Slm$coefficients[5]+(alfai^5)*Slm$coefficients[6]

###

Slm<-loess(stoth~alfai)
bb<-predict(Slm,alfai, se = TRUE)

B<-bb$fit

Marconi<-numeric(listasim$ns-1)
 for(i in 1:(listasim$ns-1)){Marconi[i]<-(anaE[i+1]-anaE[i])/(B[i+1]-B[i]) }

Marconihard<-numeric(listasim$ns-1)
 for(i in 1:(listasim$ns-1)){Marconihard[i]<-(mediaE[i+1]-mediaE[i])/(stoth[i+1]-stoth[i]) }

alfai<-tail(alfai, listasim$ns-1)
mediaE<-tail(mediaE, listasim$ns-1)
stoth<-tail(stoth, listasim$ns-1)
B<-tail(B, listasim$ns-1)
A<-tail(anaE, listasim$ns-1)

fim<-list("S"=stoth,"E"=mediaE, "polyS"=B,"polyE"=A,"marconi"=Marconi,"marconihard"=Marconihard, "alfa"=alfai,"ns"=listasim$ns,"nt"=listasim$nt,"nc"=listasim$nc,"dt"= listasim$dt )
return(fim)


}


analiseentropicaExai <- function(listasim ){
    
    #alfai<-seq(listasim$alfa+listasim$dalfa,listasim$alfa+listasim$ns*listasim$dalfa,listasim$dalfa)
    
    xai<-listasim$xa
    
    mediaE<-numeric(listasim$ns)
    
    
    
    if(listasim$dalfa==0)
    {
        alfai<-rep(listasim$alfa, listasim$ns)
    } else {
        alfai<-seq(listasim$alfa+listasim$dalfa,listasim$alfa+listasim$ns*listasim$dalfa,listasim$dalfa)
    }
    
    
    

    
    
    
    for(i in 1:listasim$ns){
        mediaE[i]<-mean(listasim$energia[,i])
    }
    
    
    #ENERGIA MÉDIA ANALÍTICA
    anaE<- (((xai)^2 )*(2*listasim$massa*alfai +listasim$gamma))/(2*listasim$gamma*(listasim$k +listasim$massa *( alfai^2) + alfai * listasim$gamma))
    
    
    
    # histogramas
    stoth<-numeric(listasim$ns)
    for(i in 1:listasim$ns){
        
        hhx<-hist(listasim$posicao[,i], breaks=100)
        hhv<-hist(listasim$velocidade[,i], breaks=100)
        
        ssx<- - hhx$density*log(hhx$density)
        ssx[which(is.na(ssx))]<- 0
        sx<-sum(ssx)*(1/sum(hhx$density))
        ssv<- - hhv$density*log(hhv$density)
        ssv[which(is.na(ssv))]<- 0
        sv<-sum(ssv)*(1/sum(hhv$density))
        
        
        stoth[i]<- sx + sv
    }
    
    
    # APROXIMACAO LINEAR PRA ENTROPIA CHEGANDO NO MARCONI

    
    # Slm<-lm(stoth~poly(xai,5,raw=TRUE))
    #B<-Slm$coefficients[1]+xai*Slm$coefficients[2]+(xai^2)*Slm$coefficients[3]+(xai^3)*Slm$coefficients[4]+(xai^4)*Slm$coefficients[5]+(xai^5)*Slm$coefficients[6]
    
    # Marconi<-numeric(listasim$ns-1)
    #for(i in 1:(listasim$ns-1)){Marconi[i]<-(anaE[i+1]-anaE[i])/(B[i+1]-B[i]) }
    
    # Marconihard<-numeric(listasim$ns-1)
    #for(i in 1:(listasim$ns-1)){Marconihard[i]<-(mediaE[i+1]-mediaE[i])/(stoth[i+1]-stoth[i]) }
    
    #alfai<-tail(alfai, listasim$ns-1)
    #xai<-tail(xai, listasim$ns-1)
    #mediaE<-tail(mediaE, listasim$ns-1)
    #anaE<-tail(anaE, listasim$ns-1)
    #stoth<-tail(stoth, listasim$ns-1)
    #B<-tail(B, listasim$ns-1)
    #A<-tail(A, ns-1)
    
    fim<-list("S"=stoth, "E"=mediaE,  "PolyE" = anaE,  "xa"=xai ,"alfa"=alfai,"ns"=listasim$ns,"nt"=listasim$nt,"nc"=listasim$nc,"dt"= listasim$dt)
    

    return(fim)
    
    
}

analiseentropicaEalfai <- function(listasim){
    
    alfai<-seq(listasim$alfa+listasim$dalfa,listasim$alfa+listasim$ns*listasim$dalfa,listasim$dalfa)
    xai<-seq(listasim$xa+listasim$dxa,listasim$xa+listasim$ns*listasim$dxa,listasim$dxa)
    mediaE<-numeric(listasim$ns)
    for(i in 1:listasim$ns){
        mediaE[i]<-mean(listasim$energia[,i])
    }
    
    
    #ENERGIA MÉDIA ANALÍTICA
    anaE<- (((xai)^2 )*(2*listasim$massa*alfai +listasim$gamma))/(2*listasim$gamma*(listasim$k +listasim$massa *( alfai^2) + alfai * listasim$gamma))
    
    # histogramas
    stoth<-numeric(listasim$ns)
    for(i in 1:listasim$ns){
        hhx<-hist(listasim$posicao[,i], breaks=100)
        hhv<-hist(listasim$velocidade[,i], breaks=100)
        
        
        ssx<- - hhx$density*log(hhx$density)
        ssx[which(is.na(ssx))]<-0
        sx<-sum(ssx)*(1/sum(hhx$density))
        ssv<- - hhv$density*log(hhv$density)
        ssv[which(is.na(ssv))]<-0
        sv<-sum(ssv)*(1/sum(hhv$density))


        stoth[i]<- sx + sv
    }
    
    
    # APROXIMACAO LINEAR PRA ENTROPIA CHEGANDO NO MARCONI
    
    
    Slm<-lm(stoth~I(1/alfai))
    B<-Slm$coefficients[1]+alfai*Slm$coefficients[2]+(alfai^2)*Slm$coefficients[3]+(alfai^3)*Slm$coefficients[4]+(alfai^4)*Slm$coefficients[5]+(alfai^5)*Slm$coefficients[6]
    
    Marconi<-numeric(listasim$ns-1)
    for(i in 1:(listasim$ns-1)){Marconi[i]<-(anaE[i+1]-anaE[i])/(B[i+1]-B[i]) }
    
    Marconihard<-numeric(listasim$ns-1)
    for(i in 1:(listasim$ns-1)){Marconihard[i]<-(mediaE[i+1]-mediaE[i])/(stoth[i+1]-stoth[i]) }
    
    alfai<-tail(alfai, listasim$ns-1)
    xai<-tail(xai, listasim$ns-1)
    mediaE<-tail(mediaE, listasim$ns-1)
    anaE<-tail(anaE, listasim$ns-1)
    stoth<-tail(stoth, listasim$ns-1)
    B<-tail(B, listasim$ns-1)
    #A<-tail(A, ns-1)
    
    fim<-list("S"=stoth,"E"=mediaE, "polyS"=B,"polyE"=anaE,"marconi"=Marconi,"marconihard"=Marconihard, "alfa"=alfai ,"xa"=xai,"ns"=listasim$ns,"nt"=listasim$nt,"nc"=listasim$nc,"dt"= listasim$dt)
    return(fim)
    
    
}


analiseentropicaTctepura <- function(listasim ){
    
    
    if(listasim$dalfa==0)
    {
        alfai<-rep(listasim$alfa, listasim$ns)
    } else {
        alfai<-seq(listasim$alfa+listasim$dalfa,listasim$alfa+listasim$ns*listasim$dalfa,listasim$dalfa)
    }
    
    
    
    xa<- sqrt( listasim$T*(listasim$gamma*(listasim$k+alfai*(listasim$gamma+listasim$massa*alfai)))/(listasim$massa*alfai)   )
    
    mediaE<-numeric(listasim$ns)
    for(i in 1:listasim$ns){
        mediaE[i]<-mean(listasim$energia[,i])
    }
    
    
    # histogramas
    stoth<-numeric(listasim$ns)
    for(i in 1:listasim$ns){
        
        
        
        hhx<-hist(listasim$posicao[,i], breaks=100)
        hhv<-hist(listasim$velocidade[,i], breaks=100)
        
        
        
        
        
        ssx<- - hhx$density*log(hhx$density)
        ssx[which(is.na(ssx))]<- 0
        sx<-sum(ssx)*(1/sum(hhx$density))
        ssv<- - hhv$density*log(hhv$density)
        ssv[which(is.na(ssv))]<- 0
        sv<-sum(ssv)*(1/sum(hhv$density))
        
        stoth[i]<- sx + sv
        
        
    }
    
    
    #ENERGIA MÉDIA ANALÍTICA
    anaE<- (((xa)^2 )*(2*listasim$massa*alfai +listasim$gamma))/(2*listasim$gamma*(listasim$k +listasim$massa *( alfai^2) + alfai * listasim$gamma))
    
    # APROXIMACAO LINEAR PRA ENERGIA MÉDIA E ENTROPIA CHEGANDO NO MARCONI
    #Elm<-lm(mediaE~poly(alfai,4,raw=TRUE))
    # A<-Elm$coefficients[1]+alfai*Elm$coefficients[2]+(alfai^2)*Elm$coefficients[3]+(alfai^3)*Elm$coefficients[4]+(alfai^4)*Elm$coefficients[5]
    
    #Slm<-lm(stoth~poly(alfai,5,raw=TRUE))
    # B<-Slm$coefficients[1]+alfai*Slm$coefficients[2]+(alfai^2)*Slm$coefficients[3]+(alfai^3)*Slm$coefficients[4]+(alfai^4)*Slm$coefficients[5]+(alfai^5)*Slm$coefficients[6]
    
    ###
    
    # Slm<-loess(stoth~alfai)
    #bb<-predict(Slm,alfai, se = TRUE)
    
    # B<-bb$fit
    
    #Marconi<-numeric(listasim$ns-1)
    #for(i in 1:(listasim$ns-1)){Marconi[i]<-(anaE[i+1]-anaE[i])/(B[i+1]-B[i]) }
    
    # Marconihard<-numeric(listasim$ns-1)
    #for(i in 1:(listasim$ns-1)){Marconihard[i]<-(mediaE[i+1]-mediaE[i])/(stoth[i+1]-stoth[i]) }
    
    # alfai<-tail(alfai, listasim$ns-1)
    #mediaE<-tail(mediaE, listasim$ns-1)
    #stoth<-tail(stoth, listasim$ns-1)
    
    #B<-tail(B, listasim$ns-1)
    #A<-tail(anaE, listasim$ns-1)
    
    fim<-list("S"=stoth,"E"=mediaE, "PolyE"=anaE, "alfa"=alfai,"ns"=listasim$ns,"nt"=listasim$nt,"nc"=listasim$nc,"dt"= listasim$dt )
    return(fim)
    
    
}

analiseentropicaTcteconte <- function(listasim ){
    
    
    if(listasim$dalfa==0)
    {
        alfai<-rep(listasim$alfa, listasim$ns)
    } else {
        alfai<-seq(listasim$alfa+listasim$dalfa,listasim$alfa+listasim$ns*listasim$dalfa,listasim$dalfa)
    }
    
    
    
    xa<- sqrt( listasim$T*(listasim$gamma*(listasim$k+alfai*(listasim$gamma+listasim$massa*alfai)))/(listasim$massa*alfai)   )
    
    mediaE<-numeric(listasim$ns)
    sdE<-numeric(listasim$ns)
    
    #moments package
    require('moments')
    
    
    
    
    #ENERGIA MÉDIA ANALÍTICA
    
    anaE<- (((xa)^2 )*(2*listasim$massa*alfai +listasim$gamma))/(2*listasim$gamma*(listasim$k +listasim$massa *( alfai^2) + alfai * listasim$gamma))
    
    
    momentationE<- vector("list", listasim$ns)
    
    cumulationE<- vector("list", listasim$ns)
    
    for(i in 1:listasim$ns){
        mediaE[i]<-mean(listasim$energia[,i])
        sdE[i]<-sd(listasim$energia[,i])
        
        momentationE[[i]]<-all.moments(listasim$energia[,i],6)
        
        cumulationE[[i]]<-all.cumulants(momentationE[[i]])
        
        
        
        #edgeworth - soma estandardizada
        
        edgei <- (1/(sqrt(length(listasim$energia[,i])))) * sum(   (listasim$energia[,i] - mediaE[i])/varianE     )
        
        
        edgeii <- (1/(sqrt(length(listasim$energia[,i])))) * sum(   (listasim$energia[,i] - mediaE[i])/varianE     )
        
    }
    
   
    
    # histogramas
    stoth<-numeric(listasim$ns)
    for(i in 1:listasim$ns){
        
        
        
        hhx<-hist(listasim$posicao[,i], breaks=100)
        hhv<-hist(listasim$velocidade[,i], breaks=100)
        
        
        #a ideia seria mudar a formula da entropia aqui - JAYNES
        
        
        
        
        ssx<- log(length(hhx$breaks)) - hhx$density*log(     (hhx$density)/(dunif(1,hhx$breaks[1],hhx$breaks[length(hhx$breaks)]))      )
        ssx[which(is.na(ssx))]<- 0
        
        sx<-sum(ssx)*(1/sum(hhx$density))
        
        ssv<- log(length(hhv$breaks)) - hhv$density*log(     (hhv$density)/(dunif(1,hhv$breaks[1],hhv$breaks[length(hhv$breaks)]))      )
        
        ssv[which(is.na(ssv))]<- 0
        
        sv<-sum(ssv)*(1/sum(hhv$density))
        
        stoth[i]<- sx + sv
        
        
    }
    

    
    # APROXIMACAO LINEAR PRA ENERGIA MÉDIA E ENTROPIA CHEGANDO NO MARCONI
    #Elm<-lm(mediaE~poly(alfai,4,raw=TRUE))
    # A<-Elm$coefficients[1]+alfai*Elm$coefficients[2]+(alfai^2)*Elm$coefficients[3]+(alfai^3)*Elm$coefficients[4]+(alfai^4)*Elm$coefficients[5]
    
    #Slm<-lm(stoth~poly(alfai,5,raw=TRUE))
    # B<-Slm$coefficients[1]+alfai*Slm$coefficients[2]+(alfai^2)*Slm$coefficients[3]+(alfai^3)*Slm$coefficients[4]+(alfai^4)*Slm$coefficients[5]+(alfai^5)*Slm$coefficients[6]
    
    ###
    
    # Slm<-loess(stoth~alfai)
    #bb<-predict(Slm,alfai, se = TRUE)
    
    # B<-bb$fit
    
    #Marconi<-numeric(listasim$ns-1)
    #for(i in 1:(listasim$ns-1)){Marconi[i]<-(anaE[i+1]-anaE[i])/(B[i+1]-B[i]) }
    
    # Marconihard<-numeric(listasim$ns-1)
    #for(i in 1:(listasim$ns-1)){Marconihard[i]<-(mediaE[i+1]-mediaE[i])/(stoth[i+1]-stoth[i]) }
    
    # alfai<-tail(alfai, listasim$ns-1)
    #mediaE<-tail(mediaE, listasim$ns-1)
    #stoth<-tail(stoth, listasim$ns-1)
    
    #B<-tail(B, listasim$ns-1)
    #A<-tail(anaE, listasim$ns-1)
    
    fim<-list("S"=stoth,"E"=mediaE, "PolyE"=anaE, "alfa"=alfai,"ns"=listasim$ns,"nt"=listasim$nt,"nc"=listasim$nc,"dt"= listasim$dt, "momentosE"=momentationE, "cumulantesE"=cumulationE )
    
    return(fim)
    
    
}

analiseentropicaTctecont <- function(listasim ){
    
    
    if(listasim$dalfa==0)
    {
        alfai<-rep(listasim$alfa, listasim$ns)
    } else {
        alfai<-seq(listasim$alfa+listasim$dalfa,listasim$alfa+listasim$ns*listasim$dalfa,listasim$dalfa)
    }
    
    
    
    xa<- sqrt( listasim$T*(listasim$gamma*(listasim$k+alfai*(listasim$gamma+listasim$massa*alfai)))/(listasim$massa*alfai)   )
    
    mediaE<-numeric(listasim$ns)
    sdE<-numeric(listasim$ns)
    
    #moments package
    require('moments')
    
    
    
    
    #ENERGIA MÉDIA ANALÍTICA
    
    anaE<- (((xa)^2 )*(2*listasim$massa*alfai +listasim$gamma))/(2*listasim$gamma*(listasim$k +listasim$massa *( alfai^2) + alfai * listasim$gamma))
    
    #E
    
    momentationE<- vector("list", listasim$ns)
    
    cumulationE<- vector("list", listasim$ns)
    
    #X
    
    # momentationX<- vector("list", listasim$ns)
    
    # cumulationX<- vector("list", listasim$ns)
    
    #V
    
    # momentationV<- vector("list", listasim$ns)
    
    # cumulationV<- vector("list", listasim$ns)
    
    
    for(i in 1:listasim$ns){
        mediaE[i]<-mean(listasim$energia[,i])
        sdE[i]<-sd(listasim$energia[,i])
        
        momentationE[[i]]<-all.moments(listasim$energia[,i],6)
        
        cumulationE[[i]]<-all.cumulants(momentationE[[i]])
        
        #momentationX[[i]]<-all.moments(listasim$posicao[,i],6)
        
        # cumulationX[[i]]<-all.cumulants(momentationX[[i]])
        
        # momentationV[[i]]<-all.moments(listasim$velocidade[,i],6)
        
        # cumulationV[[i]]<-all.cumulants(momentationV[[i]])
        
        #edgeworth - soma estandardizada
        
        #edgei <- (1/(sqrt(length(listasim$energia[,i])))) * sum(   (listasim$energia[,i] - mediaE[i])/varianE     )
        
        
        #  edgeii <- (1/(sqrt(length(listasim$energia[,i])))) * sum(   (listasim$energia[,i] - mediaE[i])/varianE     )
        
    }
    
    
    
    # histogramas
    stoth<-numeric(listasim$ns)
    for(i in 1:listasim$ns){
        
        
        
        hhx<-hist(listasim$posicao[,i], breaks=100)
        hhv<-hist(listasim$velocidade[,i], breaks=100)
        
        
        #a ideia seria mudar a formula da entropia aqui - JAYNES
        
        
        
        
        ssx<- log(length(hhx$breaks)) - hhx$density*log(     (hhx$density)/(dunif(1,hhx$breaks[1],hhx$breaks[length(hhx$breaks)]))      )
        
        #alteracao rolou de is.na pra !is.finite, vamo ver
        
        ssx[which(!is.finite(ssx))]<- 0
        
        sx<-sum(ssx)*(1/sum(hhx$density))
        
        ssv<- log(length(hhv$breaks)) - hhv$density*log(     (hhv$density)/(dunif(1,hhv$breaks[1],hhv$breaks[length(hhv$breaks)]))      )
        
        ssv[which(!is.finite(ssv))]<- 0
        
        sv<-sum(ssv)*(1/sum(hhv$density))
        
        stoth[i]<- sx + sv
        
        
    }
    
    
    
    # APROXIMACAO LINEAR PRA ENERGIA MÉDIA E ENTROPIA CHEGANDO NO MARCONI
    #Elm<-lm(mediaE~poly(alfai,4,raw=TRUE))
    # A<-Elm$coefficients[1]+alfai*Elm$coefficients[2]+(alfai^2)*Elm$coefficients[3]+(alfai^3)*Elm$coefficients[4]+(alfai^4)*Elm$coefficients[5]
    
    #Slm<-lm(stoth~poly(alfai,5,raw=TRUE))
    # B<-Slm$coefficients[1]+alfai*Slm$coefficients[2]+(alfai^2)*Slm$coefficients[3]+(alfai^3)*Slm$coefficients[4]+(alfai^4)*Slm$coefficients[5]+(alfai^5)*Slm$coefficients[6]
    
    ###
    
    # Slm<-loess(stoth~alfai)
    #bb<-predict(Slm,alfai, se = TRUE)
    
    # B<-bb$fit
    
    #Marconi<-numeric(listasim$ns-1)
    #for(i in 1:(listasim$ns-1)){Marconi[i]<-(anaE[i+1]-anaE[i])/(B[i+1]-B[i]) }
    
    # Marconihard<-numeric(listasim$ns-1)
    #for(i in 1:(listasim$ns-1)){Marconihard[i]<-(mediaE[i+1]-mediaE[i])/(stoth[i+1]-stoth[i]) }
    
    # alfai<-tail(alfai, listasim$ns-1)
    #mediaE<-tail(mediaE, listasim$ns-1)
    #stoth<-tail(stoth, listasim$ns-1)
    
    #B<-tail(B, listasim$ns-1)
    #A<-tail(anaE, listasim$ns-1)
    
    fim<-list("S"=stoth,"E"=mediaE, "PolyE"=anaE, "alfa"=alfai,"ns"=listasim$ns,"nt"=listasim$nt,"nc"=listasim$nc,"dt"= listasim$dt, "momentosE"=momentationE, "cumulantesE"=cumulationE )
    
    return(fim)
    
    
}


analiseentropicaTctecontmomen <- function(listasim ){
    
    
    if(listasim$dalfa==0)
    {
        alfai<-rep(listasim$alfa, listasim$ns)
    } else {
        alfai<-seq(listasim$alfa+listasim$dalfa,listasim$alfa+listasim$ns*listasim$dalfa,listasim$dalfa)
    }
    
    
    
    xa<- sqrt( listasim$T*(listasim$gamma*(listasim$k+alfai*(listasim$gamma+listasim$massa*alfai)))/(listasim$massa*alfai)   )
    
    mediaE<-numeric(listasim$ns)
    sdE<-numeric(listasim$ns)
    
    #moments package
    require('moments')
    
    
    
    
    #ENERGIA MÉDIA ANALÍTICA
    
    anaE<- (((xa)^2 )*(2*listasim$massa*alfai +listasim$gamma))/(2*listasim$gamma*(listasim$k +listasim$massa *( alfai^2) + alfai * listasim$gamma))
    
    #E
    
    momentationE<- vector("list", listasim$ns)
    
    cumulationE<- vector("list", listasim$ns)
    
    #X
    
    momentationX<- vector("list", listasim$ns)
    
    cumulationX<- vector("list", listasim$ns)
    
    #V
    
    momentationV<- vector("list", listasim$ns)
    
    cumulationV<- vector("list", listasim$ns)
    
    
    for(i in 1:listasim$ns){
        mediaE[i]<-mean(listasim$energia[,i])
        sdE[i]<-sd(listasim$energia[,i])
        
        momentationE[[i]]<-all.moments(listasim$energia[,i],6)
        
        cumulationE[[i]]<-cumulantes(momentationE[[i]])
        
        momentationX[[i]]<-all.moments(listasim$posicao[,i],6)
        
        cumulationX[[i]]<-cumulantes(momentationX[[i]])
        
        momentationV[[i]]<-all.moments(listasim$velocidade[,i],6)
        
        cumulationV[[i]]<-cumulantes(momentationV[[i]])
        
        #edgeworth - soma estandardizada
        
        #edgei <- (1/(sqrt(length(listasim$energia[,i])))) * sum(   (listasim$energia[,i] - mediaE[i])/varianE     )
        
        
        #  edgeii <- (1/(sqrt(length(listasim$energia[,i])))) * sum(   (listasim$energia[,i] - mediaE[i])/varianE     )
        
    }
    
    
    
    # histogramas
    stoth<-numeric(listasim$ns)
    for(i in 1:listasim$ns){
        
        
        
        hhx<-hist(listasim$posicao[,i], breaks=100)
        hhv<-hist(listasim$velocidade[,i], breaks=100)
        
        
        #a ideia seria mudar a formula da entropia aqui - JAYNES
        
        
        ssx<- - hhx$density*log(hhx$density)
        ssx[which(!is.finite(ssx))]<- 0
        sx<-sum(ssx)*(1/sum(hhx$density))
        ssv<- - hhv$density*log(hhv$density)
        ssv[which(!is.finite(ssv))]<- 0
        sv<-sum(ssv)*(1/sum(hhv$density))

        
        
        
        # ssx<-  - hhx$density*log(     (hhx$density)/(dunif(1,hhx$breaks[1],hhx$breaks[length(hhx$breaks)]))      )# + log(length(hhx$breaks))
        
        #alteracao rolou de is.na pra !is.finite, vamo ver
        
        
        
        # sx<-sum(ssx)*(1/sum(hhx$density))
        
        #ssv<- log(length(hhv$breaks)) - hhv$density*log(     (hhv$density)/(dunif(1,hhv$breaks[1],hhv$breaks[length(hhv$breaks)]))      ) #
        #
        
        
        # sv<-sum(ssv)*(1/sum(hhv$density))
        
        stoth[i]<- sx + sv
        
        
    }
    
    
    
    # APROXIMACAO LINEAR PRA ENERGIA MÉDIA E ENTROPIA CHEGANDO NO MARCONI
    #Elm<-lm(mediaE~poly(alfai,4,raw=TRUE))
    # A<-Elm$coefficients[1]+alfai*Elm$coefficients[2]+(alfai^2)*Elm$coefficients[3]+(alfai^3)*Elm$coefficients[4]+(alfai^4)*Elm$coefficients[5]
    
    #Slm<-lm(stoth~poly(alfai,5,raw=TRUE))
    # B<-Slm$coefficients[1]+alfai*Slm$coefficients[2]+(alfai^2)*Slm$coefficients[3]+(alfai^3)*Slm$coefficients[4]+(alfai^4)*Slm$coefficients[5]+(alfai^5)*Slm$coefficients[6]
    
    ###
    
    # Slm<-loess(stoth~alfai)
    #bb<-predict(Slm,alfai, se = TRUE)
    
    # B<-bb$fit
    
    #Marconi<-numeric(listasim$ns-1)
    #for(i in 1:(listasim$ns-1)){Marconi[i]<-(anaE[i+1]-anaE[i])/(B[i+1]-B[i]) }
    
    # Marconihard<-numeric(listasim$ns-1)
    #for(i in 1:(listasim$ns-1)){Marconihard[i]<-(mediaE[i+1]-mediaE[i])/(stoth[i+1]-stoth[i]) }
    
    # alfai<-tail(alfai, listasim$ns-1)
    #mediaE<-tail(mediaE, listasim$ns-1)
    #stoth<-tail(stoth, listasim$ns-1)
    
    #B<-tail(B, listasim$ns-1)
    #A<-tail(anaE, listasim$ns-1)
    
    fim<-list("S"=stoth,"E"=mediaE, "PolyE"=anaE, "alfa"=alfai,"ns"=listasim$ns,"nt"=listasim$nt,"nc"=listasim$nc,"dt"= listasim$dt, "momentosE"=momentationE, "cumulantesE"=cumulationE, "momentosX"=momentationX, "cumulantesX"=cumulationX, "momentosV"=momentationV, "cumulantesV"=cumulationV )
    
    return(fim)
    
    
}


analiseentropicamomen <- function(listasim ){
    
    
    if(listasim$dalfa==0)
    {
        alfai<-rep(listasim$alfa, listasim$ns)
    } else {
        alfai<-seq(listasim$alfa+listasim$dalfa,listasim$alfa+listasim$ns*listasim$dalfa,listasim$dalfa)
    }
    
    if(listasim$dalfa==0)
    {
        xai<-rep(listasim$xa, listasim$ns)
    } else {
        xai<-seq(listasim$xa+listasim$dxa,listasim$xa+listasim$ns*listasim$dxa,listasim$dxa)
    }
    
    mediaE<-numeric(listasim$ns)
    sdE<-numeric(listasim$ns)
    
    #moments package
    require('moments')
    
    
    
    
    #ENERGIA MÉDIA ANALÍTICA
    
    anaE<- (((xai)^2 )*(2*listasim$massa*alfai +listasim$gamma))/(2*listasim$gamma*(listasim$k +listasim$massa *( alfai^2) + alfai * listasim$gamma))
    
    #E
    
    momentationE<- vector("list", listasim$ns)
    
    cumulationE<- vector("list", listasim$ns)
    
    #X
    
    momentationX<- vector("list", listasim$ns)
    
    cumulationX<- vector("list", listasim$ns)
    
    #V
    
    momentationV<- vector("list", listasim$ns)
    
    cumulationV<- vector("list", listasim$ns)
    
    
    for(i in 1:listasim$ns){
        
        mediaE[i]<-mean(listasim$energia[,i])
        
        sdE[i]<-sd(listasim$energia[,i])
        
        momentationE[[i]]<-all.moments(listasim$energia[,i],6)
        
        cumulationE[[i]]<-cumulantes(momentationE[[i]])
        
        momentationX[[i]]<-all.moments(listasim$posicao[,i],6)
        
        cumulationX[[i]]<-cumulantes(momentationX[[i]])
        
        momentationV[[i]]<-all.moments(listasim$velocidade[,i],6)
        
        cumulationV[[i]]<-cumulantes(momentationV[[i]])
        
        
    }
    
    
    stoth<-numeric(listasim$ns)
    for(i in 1:listasim$ns){
        
        
        
        hhx<-hist(listasim$posicao[,i], breaks=100)
        hhv<-hist(listasim$velocidade[,i], breaks=100)
        
        
        ssx<- - hhx$density*log(hhx$density)
        ssx[which(!is.finite(ssx))]<- 0
        sx<-sum(ssx)*(1/sum(hhx$density))
        ssv<- - hhv$density*log(hhv$density)
        ssv[which(!is.finite(ssv))]<- 0
        sv<-sum(ssv)*(1/sum(hhv$density))
        
        
        stoth[i]<- sx + sv
        
        
    }
    

    
    fim<-list("S"=stoth,"E"=mediaE, "PolyE"=anaE, "alfa"=alfai,"ns"=listasim$ns,"nt"=listasim$nt,"nc"=listasim$nc,"dt"= listasim$dt, "momentosE"=momentationE, "cumulantesE"=cumulationE, "momentosX"=momentationX, "cumulantesX"=cumulationX, "momentosV"=momentationV, "cumulantesV"=cumulationV )
    
    return(fim)
    
    
}

analiseentropicamomenFIT <- function(listasim ){
    
    
    if(listasim$dalfa==0)
    {
        alfai<-rep(listasim$alfa, listasim$ns)
    } else {
        alfai<-seq(listasim$alfa+listasim$dalfa,listasim$alfa+listasim$ns*listasim$dalfa,listasim$dalfa)
    }
    
    if(listasim$dalfa==0)
    {
        xai<-rep(listasim$xa, listasim$ns)
    } else {
        xai<-seq(listasim$xa+listasim$dxa,listasim$xa+listasim$ns*listasim$dxa,listasim$dxa)
    }
    
    mediaE<-numeric(listasim$ns)
    sdE<-numeric(listasim$ns)
    
    #moments package
    require('moments')
    
    
    
    
    #ENERGIA MÉDIA ANALÍTICA
    
    anaE<- (((xai)^2 )*(2*listasim$massa*alfai +listasim$gamma))/(2*listasim$gamma*(listasim$k +listasim$massa *( alfai^2) + alfai * listasim$gamma))
    
    #E
    
    momentationE<- vector("list", listasim$ns)
    
    cumulationE<- vector("list", listasim$ns)
    
    #X
    
    momentationX<- vector("list", listasim$ns)
    
    cumulationX<- vector("list", listasim$ns)
    
    #V
    
    momentationV<- vector("list", listasim$ns)
    
    cumulationV<- vector("list", listasim$ns)
    
    
    for(i in 1:listasim$ns){
        
        mediaE[i]<-mean(listasim$energia[,i])
        
        sdE[i]<-sd(listasim$energia[,i])
        
        momentationE[[i]]<-all.moments(listasim$energia[,i],10)
        
        #cumulationE[[i]]<-cumulantes(momentationE[[i]])
        
        momentationX[[i]]<-all.moments(listasim$posicao[,i],10)
        
        # cumulationX[[i]]<-cumulantes(momentationX[[i]])
        
        momentationV[[i]]<-all.moments(listasim$velocidade[,i],10)
        
        #cumulationV[[i]]<-cumulantes(momentationV[[i]])
        
        
    }
    
    
    stoth<-numeric(listasim$ns)
    for(i in 1:listasim$ns){
        
        
        
        hhx<-hist(listasim$posicao[,i], breaks=100)
        hhv<-hist(listasim$velocidade[,i], breaks=100)
        
        
        ssx<- - hhx$density*log(hhx$density)
        ssx[which(!is.finite(ssx))]<- 0
        sx<-sum(ssx)*(1/sum(hhx$density))
        ssv<- - hhv$density*log(hhv$density)
        ssv[which(!is.finite(ssv))]<- 0
        sv<-sum(ssv)*(1/sum(hhv$density))
        
        
        stoth[i]<- sx + sv
        
        
    }
    
    #experimento<-data.frame(bx,hx)
    #write.table(experimento, "histposicao.csv",sep=',' ,row.names=FALSE, col.names=FALSE)
    
    #experimento<-data.frame(bv,hv)
    #write.table(experimento, "histvelocidade.csv",sep=',' ,row.names=FALSE, col.names=FALSE)

#experimento<-data.frame(seq(0,6),momenV)
#>      write.table(experimento, "momenv.csv",sep=',' ,row.names=FALSE, col.names=FALSE)
#> experimento<-data.frame(seq(0,6),momenX)
#>      write.table(experimento, "momenx.csv",sep=',' ,row.names=FALSE, col.names=FALSE)


    fim<-list("S"=stoth,"E"=mediaE, "PolyE"=anaE, "alfa"=alfai,"ns"=listasim$ns,"nt"=listasim$nt,"nc"=listasim$nc,"dt"= listasim$dt, "momentosE"=momentationE,  "momentosX"=momentationX, "momentosV"=momentationV )
    
    return(fim)
    
    
}


analiseentropicaSIMPLIS <- function(listasim ){
    
    
    if(listasim$dalfa==0)
    {
        alfai<-rep(listasim$alfa, listasim$ns)
    } else {
        alfai<-seq(listasim$alfa+listasim$dalfa,listasim$alfa+listasim$ns*listasim$dalfa,listasim$dalfa)
    }
    
    if(listasim$dxa==0)
    {
        xai<-rep(listasim$xa, listasim$ns)
    } else {
        xai<-seq(listasim$xa+listasim$dxa,listasim$xa+listasim$ns*listasim$dxa,listasim$dxa)
    }
    
    
    
    mediaE<-numeric(listasim$ns)
    sdE<-numeric(listasim$ns)
    
    
    
    
    #ENERGIA MÉDIA ANALÍTICA
    
    anaE<- (((xai)^2 )*(2*listasim$massa*alfai +listasim$gamma))/(2*listasim$gamma*(listasim$k +listasim$massa *( alfai^2) + alfai * listasim$gamma))
    

     for(i in 1:listasim$ns){
         
        mediaE[i]<-mean(listasim$energia[,i])
        
        sdE[i]<-sd(listasim$energia[,i])

    }
    
    
    
    # histogramas
    stoth<-numeric(listasim$ns)
    for(i in 1:listasim$ns){
        
        
        
        hhx<-hist(listasim$posicao[,i], breaks=100, plot = FALSE)
        hhv<-hist(listasim$velocidade[,i], breaks=100, plot = FALSE)
        
        
        #a ideia seria mudar a formula da entropia aqui - JAYNES
        
        
        ssx<- - hhx$density*log(hhx$density)
        ssx[which(!is.finite(ssx))]<- 0
        sx<-sum(ssx)*(1/sum(hhx$density))
        ssv<- - hhv$density*log(hhv$density)
        ssv[which(!is.finite(ssv))]<- 0
        sv<-sum(ssv)*(1/sum(hhv$density))
        
        
        
        
        stoth[i]<- sx + sv
        
        
    }
    

    fim<-list("S"=stoth,"E"=mediaE, "PolyE"=anaE, "alfa"=alfai, "xa"=xai,  "ns"=listasim$ns,"nt"=listasim$nt,"nc"=listasim$nc,"dt"= listasim$dt)
    
    return(fim)
    
    
}

analiseentropicaSIMPLISTT <- function(listasim ){
    
    
    if(listasim$dalfa==0)
    {
        alfai<-rep(listasim$alfa, listasim$ns)
    } else {
        alfai<-seq(listasim$alfa+listasim$dalfa,listasim$alfa+listasim$ns*listasim$dalfa,listasim$dalfa)
    }
    

    xai<- listasim$xa
    
    
    
    mediaE<-numeric(listasim$ns)
    sdE<-numeric(listasim$ns)
    
    
    
    
    #ENERGIA MÉDIA ANALÍTICA
    
    anaE<- (((xai)^2 )*(2*listasim$massa*alfai +listasim$gamma))/(2*listasim$gamma*(listasim$k +listasim$massa *( alfai^2) + alfai * listasim$gamma))
    
    
    for(i in 1:listasim$ns){
        
        mediaE[i]<-mean(listasim$energia[,i])
        
        sdE[i]<-sd(listasim$energia[,i])
        
    }
    
    
    
    # histogramas
    stoth<-numeric(listasim$ns)
    for(i in 1:listasim$ns){
        

        hhx<-hist(listasim$posicao[,i], breaks=100, plot = FALSE)
        hhv<-hist(listasim$velocidade[,i], breaks=100, plot = FALSE)
        
        #a ideia seria mudar a formula da entropia aqui - JAYNES
        
        ssx<- - hhx$density*log(hhx$density)
        ssx[which(!is.finite(ssx))]<- 0
        sx<-sum(ssx)*(1/sum(hhx$density))
        ssv<- - hhv$density*log(hhv$density)
        ssv[which(!is.finite(ssv))]<- 0
        sv<-sum(ssv)*(1/sum(hhv$density))
        
        
        
        stoth[i]<- sx + sv
        
        
    }
    
    
    fim<-list("S"=stoth,"E"=mediaE, "PolyE"=anaE, "alfa"=alfai, "xa"=xai,  "ns"=listasim$ns,"nt"=listasim$nt,"nc"=listasim$nc,"dt"= listasim$dt)
    
    return(fim)
    
    
}


analiseentropicaSIMPLIStemporal <- function(listasim ){
    
    
    if(listasim$dalfa==0)
    {
        alfai<-rep(listasim$alfa, listasim$ns)
    } else {
        alfai<-seq(listasim$alfa+listasim$dalfa,listasim$alfa+listasim$ns*listasim$dalfa,listasim$dalfa)
    }
    
    if(listasim$dxa==0)
    {
        xai<-rep(listasim$xa, listasim$ns)
    } else {
        xai<-seq(listasim$xa+listasim$dxa,listasim$xa+listasim$ns*listasim$dxa,listasim$dxa)
    }
    
    
    
    mediaE<-numeric(listasim$ns)
    sdE<-numeric(listasim$ns)
    
    
    
    
    #ENERGIA MÉDIA ANALÍTICA
    
    anaE<- (((xai)^2 )*(2*listasim$massa*alfai +listasim$gamma))/(2*listasim$gamma*(listasim$k +listasim$massa *( alfai^2) + alfai * listasim$gamma))
    
    
    for(i in 1:listasim$ns){
        
        mediaE[i]<-mean(listasim$energia[,i])
        
        sdE[i]<-sd(listasim$energia[,i])
        
    }
    
    
    
    # histogramas
    stoth<-numeric(listasim$ns)
    for(i in 1:listasim$ns){
        
        
        
        hhx<-hist(listasim$posicao[,i], breaks=100)
        hhv<-hist(listasim$velocidade[,i], breaks=100)
        
        
        #a ideia seria mudar a formula da entropia aqui - JAYNES
        
        
        ssx<- - hhx$density*log(hhx$density)
        ssx[which(!is.finite(ssx))]<- 0
        sx<-sum(ssx)*(1/sum(hhx$density))
        ssv<- - hhv$density*log(hhv$density)
        ssv[which(!is.finite(ssv))]<- 0
        sv<-sum(ssv)*(1/sum(hhv$density))
        
        
        
        
        stoth[i]<- sx + sv
        
        
    }
    
    
    fim<-list("S"=stoth,"E"=mediaE, "PolyE"=anaE, "alfa"=alfai, "xa"=xai,  "ns"=listasim$ns,"nt"=listasim$nt,"nc"=listasim$nc,"dt"= listasim$dt)
    
    return(fim)
    
    
}


analiseentropicaSIMPLISquadrada <- function(listasim ){
    
    
    nns<-sqrt(listasim$ns)

    alfai<-rep(listasim$alfa, nns)

    xai<-seq(listasim$xa+listasim$dxa,listasim$xa+nns*listasim$dxa,listasim$dxa)

    for(i in 1:nns)
    {
        alfai<-c(alfai,rep(listasim$alfa + i*listasim$dalfa , nns))
        xai<-c(xai,  seq(listasim$xa+listasim$dxa,listasim$xa+nns*listasim$dxa,listasim$dxa)   )
    }
    
    
    mediaE<-numeric(listasim$ns)
    sdE<-numeric(listasim$ns)
    
    
    
    
    #ENERGIA MÉDIA ANALÍTICA
    
    anaE<- (((xai)^2 )*(2*listasim$massa*alfai +listasim$gamma))/(2*listasim$gamma*(listasim$k +listasim$massa *( alfai^2) + alfai * listasim$gamma))
    
    
    for(i in 1:listasim$ns){
        
        mediaE[i]<-mean(listasim$energia[,i])
        
        sdE[i]<-sd(listasim$energia[,i])
        
    }
    
    
    
    # histogramas
    stoth<-numeric(listasim$ns)
    for(i in 1:listasim$ns){
        
        
        
        hhx<-hist(listasim$posicao[,i], breaks=100)
        hhv<-hist(listasim$velocidade[,i], breaks=100)
        
        
        
        histv<-data.frame ( hhv$mids  , hhv$density)
        
        histx<-data.frame ( hhx$mids , hhx$density)

        write.table(histv, paste(toString(i), "histv.csv") , sep=",", row.names=FALSE, col.names=FALSE)

        write.table(histx, paste(toString(i), "histx.csv") , sep=",", row.names=FALSE, col.names=FALSE)
        
        
        #a ideia seria mudar a formula da entropia aqui - JAYNES
        
        
        ssx<- - hhx$density*log(hhx$density)
        ssx[which(!is.finite(ssx))]<- 0
        sx<-sum(ssx)*(1/sum(hhx$density))
        ssv<- - hhv$density*log(hhv$density)
        ssv[which(!is.finite(ssv))]<- 0
        sv<-sum(ssv)*(1/sum(hhv$density))
        
        
        
        
        stoth[i]<- sx + sv
        
        
    }
    
    
    fim<-list("S"=stoth,"E"=mediaE, "PolyE"=anaE, "alfa"=alfai, "xa"=xai,  "ns"=listasim$ns,"nt"=listasim$nt,"nc"=listasim$nc,"dt"= listasim$dt)
    
    return(fim)
    
    
}


analiseentropicaQTctecontmomen <- function(listasim, Q){
    
    
    if(listasim$dalfa==0)
    {
        alfai<-rep(listasim$alfa, listasim$ns)
    } else {
        alfai<-seq(listasim$alfa+listasim$dalfa,listasim$alfa+listasim$ns*listasim$dalfa,listasim$dalfa)
    }
    
    
    
    xa<- sqrt( listasim$T*(listasim$gamma*(listasim$k+alfai*(listasim$gamma+listasim$massa*alfai)))/(listasim$massa*alfai)   )
    
    mediaE<-numeric(listasim$ns)
    sdE<-numeric(listasim$ns)
    
    #moments package
    require('moments')
    
    
    
    
    #ENERGIA MÉDIA ANALÍTICA
    
    anaE<- (((xa)^2 )*(2*listasim$massa*alfai +listasim$gamma))/(2*listasim$gamma*(listasim$k +listasim$massa *( alfai^2) + alfai * listasim$gamma))
    
    #E
    
    momentationE<- vector("list", listasim$ns)
    
    cumulationE<- vector("list", listasim$ns)
    
    #X
    
    momentationX<- vector("list", listasim$ns)
    
    cumulationX<- vector("list", listasim$ns)
    
    #V
    
    momentationV<- vector("list", listasim$ns)
    
    cumulationV<- vector("list", listasim$ns)
    
    
    for(i in 1:listasim$ns){
        mediaE[i]<-mean(listasim$energia[,i])
        sdE[i]<-sd(listasim$energia[,i])
        
        momentationE[[i]]<-all.moments(listasim$energia[,i],6)
        
        cumulationE[[i]]<-cumulantes(momentationE[[i]])
        
        momentationX[[i]]<-all.moments(listasim$posicao[,i],6)
        
        cumulationX[[i]]<-cumulantes(momentationX[[i]])
        
        momentationV[[i]]<-all.moments(listasim$velocidade[,i],6)
        
        cumulationV[[i]]<-cumulantes(momentationV[[i]])
        
        #edgeworth - soma estandardizada
        
        #edgei <- (1/(sqrt(length(listasim$energia[,i])))) * sum(   (listasim$energia[,i] - mediaE[i])/varianE     )
        
        
        #  edgeii <- (1/(sqrt(length(listasim$energia[,i])))) * sum(   (listasim$energia[,i] - mediaE[i])/varianE     )
        
    }
    
    
    
    # histogramas
    stoth<-numeric(listasim$ns)
    for(i in 1:listasim$ns){
        
        
        
        hhx<-hist(listasim$posicao[,i], breaks=100)
        hhv<-hist(listasim$velocidade[,i], breaks=100)
        
        
        #a ideia seria mudar a formula da entropia aqui - JAYNES
        
        
        #ssx<- - hhx$density*log(hhx$density)
        
        sx<- (1- sum( hhv$density^Q) )/(Q-1)
        
        
        # ssx[which(!is.finite(ssx))]<- 0
        
        # sx<-sum(ssx)*(1/sum(hhx$density))
        
        
        
        sv<- (1- sum( hhv$density^Q) )/(Q-1)
        
        
        #ssv[which(!is.finite(ssv))]<- 0
        
        #sv<-sum(ssv)*(1/sum(hhv$density))
        
        
        
        
        # ssx<-  - hhx$density*log(     (hhx$density)/(dunif(1,hhx$breaks[1],hhx$breaks[length(hhx$breaks)]))      )# + log(length(hhx$breaks))
        
        #alteracao rolou de is.na pra !is.finite, vamo ver
        
        
        
        # sx<-sum(ssx)*(1/sum(hhx$density))
        
        #ssv<- log(length(hhv$breaks)) - hhv$density*log(     (hhv$density)/(dunif(1,hhv$breaks[1],hhv$breaks[length(hhv$breaks)]))      ) #
        #
        
        
        # sv<-sum(ssv)*(1/sum(hhv$density))
        
        stoth[i]<- sx + sv
        
        
    }
    
    
    
    # APROXIMACAO LINEAR PRA ENERGIA MÉDIA E ENTROPIA CHEGANDO NO MARCONI
    #Elm<-lm(mediaE~poly(alfai,4,raw=TRUE))
    # A<-Elm$coefficients[1]+alfai*Elm$coefficients[2]+(alfai^2)*Elm$coefficients[3]+(alfai^3)*Elm$coefficients[4]+(alfai^4)*Elm$coefficients[5]
    
    #Slm<-lm(stoth~poly(alfai,5,raw=TRUE))
    # B<-Slm$coefficients[1]+alfai*Slm$coefficients[2]+(alfai^2)*Slm$coefficients[3]+(alfai^3)*Slm$coefficients[4]+(alfai^4)*Slm$coefficients[5]+(alfai^5)*Slm$coefficients[6]
    
    ###
    
    # Slm<-loess(stoth~alfai)
    #bb<-predict(Slm,alfai, se = TRUE)
    
    # B<-bb$fit
    
    #Marconi<-numeric(listasim$ns-1)
    #for(i in 1:(listasim$ns-1)){Marconi[i]<-(anaE[i+1]-anaE[i])/(B[i+1]-B[i]) }
    
    # Marconihard<-numeric(listasim$ns-1)
    #for(i in 1:(listasim$ns-1)){Marconihard[i]<-(mediaE[i+1]-mediaE[i])/(stoth[i+1]-stoth[i]) }
    
    # alfai<-tail(alfai, listasim$ns-1)
    #mediaE<-tail(mediaE, listasim$ns-1)
    #stoth<-tail(stoth, listasim$ns-1)
    
    #B<-tail(B, listasim$ns-1)
    #A<-tail(anaE, listasim$ns-1)
    
    fim<-list("S"=stoth,"E"=mediaE, "PolyE"=anaE, "alfa"=alfai,"ns"=listasim$ns,"nt"=listasim$nt,"nc"=listasim$nc,"dt"= listasim$dt, "momentosE"=momentationE, "cumulantesE"=cumulationE, "momentosX"=momentationX, "cumulantesX"=cumulationX, "momentosV"=momentationV, "cumulantesV"=cumulationV )
    
    return(fim)
    
    
}



#output de analiseentropicaTcte (para o output de fim que sao varios quero concatenados)
#1 —— S
#2 —— E
#3 —— POLY S
#4 —— POLY E
#5 — — MARCONI
#6 — — MARCONI HARD
#7 — — ALFA
#8 — — NS
#9 — — NT
#10 — — NC
#11 — — DT


analiseentropicaTctecontXV <- function(listasim ){
    
    
    if(listasim$dalfa==0)
    {
        alfai<-rep(listasim$alfa, listasim$ns)
    } else {
        alfai<-seq(listasim$alfa+listasim$dalfa,listasim$alfa+listasim$ns*listasim$dalfa,listasim$dalfa)
    }
    
    
    
    xa<- sqrt( listasim$T*(listasim$gamma*(listasim$k+alfai*(listasim$gamma+listasim$massa*alfai)))/(listasim$massa*alfai)   )
    
    mediaE<-numeric(listasim$ns)
    sdE<-numeric(listasim$ns)
    
    #moments package
    require('moments')
    
    
    
    
    #ENERGIA MÉDIA ANALÍTICA
    
    anaE<- (((xa)^2 )*(2*listasim$massa*alfai +listasim$gamma))/(2*listasim$gamma*(listasim$k +listasim$massa *( alfai^2) + alfai * listasim$gamma))
    
    
    Epot<- 0.5*var(listasim$posicao)
    Ekin<- 0.5*var(listasim$velocidade)
    
    
    momentationE<- vector("list", listasim$ns)
    
    cumulationE<- vector("list", listasim$ns)
    
    momentationV<- vector("list", listasim$ns)
    
    cumulationV<- vector("list", listasim$ns)
    
    for(i in 1:listasim$ns){
        
        
        mediaE[i]<-mean(listasim$energia[,i])
        
        sdE[i]<-sd(listasim$energia[,i])
        
        momentationE[[i]]<-all.moments(listasim$energia[,i],6)
        
        cumulationE[[i]]<-all.cumulants(momentationE[[i]])
        
        
        
        momentationV[[i]]<-all.moments(listasim$velocidade[,i],6)
        
        cumulationV[[i]]<-all.cumulants(momentationE[[i]])
        
        #edgeworth - soma estandardizada
        
        #edgei <- (1/(sqrt(length(listasim$energia[,i])))) * sum(   (listasim$energia[,i] - mediaE[i])/varianE     )
        
        #  edgeii <- (1/(sqrt(length(listasim$energia[,i])))) * sum(   (listasim$energia[,i] - mediaE[i])/varianE     )
        
        
    }
    
    
    
    # histogramas
    stoth<-numeric(listasim$ns)
    stothx<-numeric(listasim$ns)
    stothv<-numeric(listasim$ns)
    
    for(i in 1:listasim$ns){
        
        
        
        hhx<-hist(listasim$posicao[,i], breaks=100)
        hhv<-hist(listasim$velocidade[,i], breaks=100)
        
        
        #a ideia seria mudar a formula da entropia aqui - JAYNES
        
        
        
        
        ssx<- log(length(hhx$breaks)) - hhx$density*log(     (hhx$density)/(dunif(1,hhx$breaks[1],hhx$breaks[length(hhx$breaks)]))      )
        #alteracao rolou de is.na pra !is.finite, vamo ver
        ssx[which(!is.finite(ssx))]<- 0
        
        stothx[i]<-sum(ssx)*(1/sum(hhx$density))
        
        ssv<- log(length(hhv$breaks)) - hhv$density*log(     (hhv$density)/(dunif(1,hhv$breaks[1],hhv$breaks[length(hhv$breaks)]))      )
        
        ssv[which(!is.finite(ssv))]<- 0
        
        stothv[i]<-sum(ssv)*(1/sum(hhv$density))
        
        stoth[i]<- stothx[i] + stothv[i]
        
        
    }
    
    
    
    # APROXIMACAO LINEAR PRA ENERGIA MÉDIA E ENTROPIA CHEGANDO NO MARCONI
    #Elm<-lm(mediaE~poly(alfai,4,raw=TRUE))
    # A<-Elm$coefficients[1]+alfai*Elm$coefficients[2]+(alfai^2)*Elm$coefficients[3]+(alfai^3)*Elm$coefficients[4]+(alfai^4)*Elm$coefficients[5]
    
    #Slm<-lm(stoth~poly(alfai,5,raw=TRUE))
    # B<-Slm$coefficients[1]+alfai*Slm$coefficients[2]+(alfai^2)*Slm$coefficients[3]+(alfai^3)*Slm$coefficients[4]+(alfai^4)*Slm$coefficients[5]+(alfai^5)*Slm$coefficients[6]
    
    ###
    
    # Slm<-loess(stoth~alfai)
    #bb<-predict(Slm,alfai, se = TRUE)
    
    # B<-bb$fit
    
    #Marconi<-numeric(listasim$ns-1)
    #for(i in 1:(listasim$ns-1)){Marconi[i]<-(anaE[i+1]-anaE[i])/(B[i+1]-B[i]) }
    
    # Marconihard<-numeric(listasim$ns-1)
    #for(i in 1:(listasim$ns-1)){Marconihard[i]<-(mediaE[i+1]-mediaE[i])/(stoth[i+1]-stoth[i]) }
    
    # alfai<-tail(alfai, listasim$ns-1)
    #mediaE<-tail(mediaE, listasim$ns-1)
    #stoth<-tail(stoth, listasim$ns-1)
    
    #B<-tail(B, listasim$ns-1)
    #A<-tail(anaE, listasim$ns-1)
    
    fim<-list("Sx"=stothx, "Sv"=stothv ,"S"=stoth,"E"=mediaE, "PolyE"=anaE,"Epot"=Epot ,"Ekin"=Ekin,  "alfa"=alfai,"ns"=listasim$ns,"nt"=listasim$nt,"nc"=listasim$nc,"dt"= listasim$dt, "momentosE"=momentationE, "cumulantesE"=cumulationE, "momentosV"=momentationV, "cumulantesE"=cumulationV )
    
    return(fim)
    
    
}



salvagraficosalfasimples<-function(quero){
    
    a<-toString(quero$ns)
    b<-toString(quero$nt)
    c<-toString(quero$nc)
    d<-toString(2*quero$alfa[1]-quero$alfa[2])
    e<-toString(quero$alfa[2]-quero$alfa[1])
    
    
    
    grafs<-c("Spoly",a, "ns",b, "nt",c, "nc",d , "alfa",e, "dalfa"," .jpg" )
    grafs<-toString(grafs)
    
    jpeg(grafs)
    
    plot(quero$S ~ quero$alfa )
    lines(quero$polyS~ quero$alfa, col='blue')

    dev.off()
    
    grafe<-c("Epoly",a, "ns",b, "nt",c, "nc",d , "alfa",e, "dalfa"," .jpg" )
    grafe<-toString(grafe)
    
    jpeg(grafe)


    plot(quero$E ~ quero$alfa)
    lines(quero$polyE~ quero$alfa, col='blue')
    
    dev.off()
    
    #    grafm<-c("Marconipoly",a, "ns",b, "nt",c, "nc",d , "alfa",e, "dalfa"," .png" )
    #grafm<-toString(grafm)
    
    # plot(quero$marconi ~ quero$alfa )
    #dev.copy(png, grafm )
    #dev.off()
    
    # grafmm<-c("Marconihard",a, "ns",b, "nt",c, "nc",d , "alfa",e, "dalfa"," .png" )
    #grafmm<-toString(grafmm)
    
    #  plot(quero$marconihard ~ quero$alfa )
    #dev.copy(png, grafmm )
    #dev.off()
    
}


salvagraficosalfa<-function(quero){

a<-toString(quero$ns)
b<-toString(quero$nt)
c<-toString(quero$nc)
d<-toString(2*quero$alfa[1]-quero$alfa[2])
e<-toString(quero$alfa[2]-quero$alfa[1])



grafs<-c("Spoly",a, "ns",b, "nt",c, "nc",d , "alfa",e, "dalfa"," .png" )
grafs<-toString(grafs)

plot(quero$S ~ quero$alfa )
lines(quero$polyS~ quero$alfa, col='blue')
dev.copy(png, grafs )
dev.off()

grafe<-c("Epoly",a, "ns",b, "nt",c, "nc",d , "alfa",e, "dalfa"," .png" )
grafe<-toString(grafe)

plot(quero$E ~ quero$alfa)
lines(quero$polyE~ quero$alfa, col='blue')
dev.copy(png, grafe )
dev.off()

grafm<-c("Marconipoly",a, "ns",b, "nt",c, "nc",d , "alfa",e, "dalfa"," .png" )
grafm<-toString(grafm)

plot(quero$marconi ~ quero$alfa )
dev.copy(png, grafm )
dev.off()

grafmm<-c("Marconihard",a, "ns",b, "nt",c, "nc",d , "alfa",e, "dalfa"," .png" )
grafmm<-toString(grafmm)

plot(quero$marconihard ~ quero$alfa )
dev.copy(png, grafmm )
dev.off()

}

salvagraficosxa<-function(quero, ns, nt, nc , xa, dxa){
    
    a<-toString(ns)
    b<-toString(nt)
    c<-toString(nc)
    d<-toString(xa)
    e<-toString(dxa)
    
    
    
    grafs<-c("Spoly",a, "ns",b, "nt",c, "nc",d , "xa",e, "dxa"," .png" )
    grafs<-toString(grafs)
    
    plot(quero$S ~ quero$xa )
    lines(quero$polyS~ quero$xa, col='blue')
    dev.copy(png, grafs )
    dev.off()
    
    grafe<-c("Epoly",a, "ns",b, "nt",c, "nc",d , "xa",e, "dxa"," .png" )
    grafe<-toString(grafe)
    
    plot(quero$E ~ quero$xa)
    lines(quero$polyE~ quero$xa, col='blue')
    dev.copy(png, grafe )
    dev.off()
    
    grafm<-c("Marconipoly",a, "ns",b, "nt",c, "nc",d , "xa",e, "dxa"," .png" )
    grafm<-toString(grafm)
    
    plot(quero$marconi ~ quero$xa )
    dev.copy(png, grafm )
    dev.off()
    
    grafmm<-c("Marconihard",a, "ns",b, "nt",c, "nc",d , "xa",e, "dxa"," .png" )
    grafmm<-toString(grafmm)
    
    plot(quero$marconihard ~ quero$xa )
    dev.copy(png, grafmm )
    dev.off()
    
}

salvavariosgraficosalfa<-function(fim, temp){
    
    #color<-grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = TRUE)]
    color<-grDevices::rainbow(30)
    
    cores<-tail(color,20)
    color<-rev(color)
    
    a<-toString(fim$ns)
    b<-toString(fim$nt)
    c<-toString(fim$nc)
    d<-toString(2*fim$alfa[1]-fim$alfa[2])
    e<-toString(fim$alfa[2]-fim$alfa[1])
    f<-toString(temp[length(temp)])
    g<-toString(temp[1])


    #grafico de entropia

    grafs<-c("Spoly",a, "ns",b, "nt",c, "nc",d , "alfa",e, "dalfa" ,".jpg" )
    grafs<-toString(grafs)
    jpeg(grafs)
    
    
    #output de analiseentropicaTcte (para o output de fim que sao varios quero concatenados)
    #1 —— S
    #2 —— E
    #3 —— POLY S
    #4 —— POLY E
    #5 — — MARCONI
    #6 — — MARCONI HARD
    #7 — — ALFA
    #8 — — NS
    #9 — — NT
    #10 — — NC
    #11 — — DT
    
    plot(fim$S ~ fim$alfa, col=cores[1],ylim=c(0, 10) )
    lines(fim$polyS~ fim$alfa, col=cores[1])

    
    for (i in 1:9)
    {
        aa<- 11*i +1
        points(fim[[aa]]~fim$alfa,col=cores[i+1] )
        lines(fim[[11*i +3]]~fim$alfa, col=cores[i+1] )
    }
    
    dev.off()
    
    #grafico de energia
    
    
    grafe<-c("Epoly",a, "ns",b, "nt",c, "nc",d , "alfa",e, "dalfa",".jpg" )
    grafe<-toString(grafe)
    jpeg(grafe)
    
    plot(fim$E ~ fim$alfa,col=cores[1], ylim=c(0, 10))
    lines(fim$polyE ~ fim$alfa, col=cores[1])
    
    for (i in 1:9)
    {
        points(fim[[11*i +2]]~fim$alfa,col=cores[i+1] )
        lines(fim[[11*i +4]]~fim$alfa, col=cores[i+1] )
    }
    
    dev.off()
    
    #A DERIVADA QUE QUEREMOS, SUAVIZADA COM LOESS
    
    grafm<-c("Marconipoly",a, "ns",b, "nt",c, "nc",d , "alfa",e, "dalfa",".jpg" )
    grafm<-toString(grafm)
    jpeg(grafm)
    
    plot(fim$marconi ~ fim$alfa )
    
    for (i in 1:9)
    {
        points(fim[[11*i +5]]~fim$alfa,col=cores[i+1] )
    }
    
    dev.off()
    
    grafmm<-c("Marconihard",a, "ns",b, "nt",c, "nc",d , "alfa",e, "dalfa",".jpg" )
    grafmm<-toString(grafmm)
    jpeg(grafmm)
    #A DERIVADA HARD QUE CREIO N VAI ROLAR, MAS BORA
    
    plot(fim$marconihard ~ fim$alfa )
    
    for (i in 1:9)
    {
        points(fim[[11*i +6]]~fim$alfa,col=cores[i+1] )
    }
    
    dev.off()
    
}


variosT<-function(ns, nt, nc, dt, massa, gamma, ctelastica, alfa, dalfa,temp)
{
    #temp<-seq(0.5,5.0,0.5)
    
    fim<-list()
    
    for(i in 1:length(temp))
    {
        listasim<-SimulacaoBimodalTcte(ns, nt, nc, dt, massa, gamma, ctelastica, alfa, dalfa, temp[i])
        quero<-analiseentropicaTcte(listasim)
        fim<-c(fim,quero)
        
    }
    return(fim)
}

variosTamortcont<-function(ns, nt, nc, dt, massa, gamma, ctelastica, alfa, dalfa,temp)
{
    #temp<-seq(0.5,5.0,0.5)
    
    fim<-list()
    
    for(i in 1:length(temp))
    {
        listasim<-SimAmortizadoraTctecont(ns, nt, nc, dt, massa, gamma, ctelastica, alfa, dalfa, temp[i])
        quero<-analiseentropicaTcte(listasim)
        fim<-c(fim,quero)
        
    }
    return(fim)
}





SimAmortizadoraTcte<- function(nns, nx, ns, nt, nc, dt, massa, gamma, ctelastica, alfa, dalfa, temp){

# zoado... nx nao aumentou a precisao do valor da energia, nc aumentou, preciso testar o nns
    
    anaE <- numeric(ns)
    
    energia <- numeric(ns)
    
    entropia <- numeric(ns)
    
    analiticE <- numeric(ns)
    
    erroE <- numeric(ns)
    
    
    
    
    ddalfa <- seq(alfa +dalfa, alfa+ ns * dalfa, dalfa)
    
    
    
    
    for(j in 1:ns)
    {
        
        momentos<-numeric(nx)
        ss<- numeric(nx)
        
        ee<- numeric(nx)
        
        
        for(i in 1:nns)
        {
            listasim<-SimulacaoBimodalTcte(nx, nt, nc, dt, massa, gamma, ctelastica, ddalfa[j], 0, temp)
            keru<-analiseentropicaTctepura(listasim)
            ss<-ss+(keru$S)/nns
            ee<-ee+(keru$E)/nns
        }
        
        sss<-mean(ss)
        eee<-mean(ee)
        
        anaE<-keru$PolyE[1]
        
        erroE[j]<-((eee-anaE)/anaE)*100
        
        entropia[j]<-sss
        
        energia[j]<-eee
        
        analiticE[j]<-anaE
    
    
    }
    
    Slm<-loess(entropia~ddalfa)
    
    bb<-predict(Slm,ddalfa, se = TRUE)
    
    B<-bb$fit

    
    #FIM DO SIMULACAO BIMODAL T CTE
    #   final<-list("posicao"=xx,"velocidade"=vv, "energia"=ee, "massa"=massa,"gamma"=gamma,"k"=ctelastica, "alfa" = alfa,"dalfa"= dalfa, "xa"=xa, "ns"=ns,"nt"=nt,"nc"=nc,"dt"= dt, "T"=Temp)
    #return(final)
    
    final<-list("PolyE"=analiticE, "E"=energia, "S"=entropia, "PolyS"= B , "alfa"=ddalfa, "nns" = nns, "ns"=ns, "nt"=nt, "nc"=nc, "dt"= dt, "temp" = temp, "erroE" = erroE)
    
    return(final)
    
}





SimAmortizadoraTctelog<- function(nns, nx, ns, nt, nc, dt, massa, gamma, ctelastica, alfa, dalfa, temp){
    
    # zoado... nx nao aumentou a precisao do valor da energia, nc aumentou, preciso testar o nns
    
    anaE <- numeric(ns)
    
    energia <- numeric(ns)
    
    entropia <- numeric(ns)
    
    analiticE <- numeric(ns)
    
    erroE <- numeric(ns)
    
    
    
    
    ddalfa <- lseq(dalfa, alfa, ns)
    
    
    
    
    for(j in 1:ns)
    {
        
        
        ss<- numeric(nx)
        
        ee<- numeric(nx)
        
        
        for(i in 1:nns)
        {
            listasim<-SimulacaoBimodalTcte(nx, nt, nc, dt, massa, gamma, ctelastica, ddalfa[j], 0, temp)
            keru<-analiseentropicaTctepura(listasim)
            ss<-ss+(keru$S)/nns
            ee<-ee+(keru$E)/nns
        }
        
        sss<-mean(ss)
        eee<-mean(ee)
        
        anaE<-keru$PolyE[1]
        
        erroE[j]<-((eee-anaE)/anaE)*100
        
        entropia[j]<-sss
        
        energia[j]<-eee
        
        analiticE[j]<-anaE
        
        
    }
    
    Slm<-loess(entropia~ddalfa)
    
    bb<-predict(Slm,ddalfa, se = TRUE)
    
    B<-bb$fit
    
    
    #FIM DO SIMULACAO BIMODAL T CTE
    #   final<-list("posicao"=xx,"velocidade"=vv, "energia"=ee, "massa"=massa,"gamma"=gamma,"k"=ctelastica, "alfa" = alfa,"dalfa"= dalfa, "xa"=xa, "ns"=ns,"nt"=nt,"nc"=nc,"dt"= dt, "T"=Temp)
    #return(final)
    
    final<-list("momentosE" ,"PolyE"=analiticE, "E"=energia, "S"=entropia, "PolyS"= B , "alfa"=ddalfa, "nns" = nns, "ns"=ns, "nt"=nt, "nc"=nc, "dt"= dt, "temp" = temp, "erroE" = erroE)
    
    return(final)
    
}

SimAmortizadoraTctecont<- function(nns, nx, ns, nt, nc, dt, massa, gamma, ctelastica, alfa, dalfa, temp){
    
    
    
    
    
    # nns<-10
    #nx<-10
    #ns<-10
    #nt<-100
    #nc<-100
    #dt<-0.001
    #massa<-1
    #gamma<-1
    #ctelastica<-1
    #alfa<-0
    #dalfa<-0.01
    #temp<-1
    #
    #
    # zoado... nx nao aumentou a precisao do valor da energia, nc aumentou, preciso testar o nns
    
    anaE <- numeric(ns)
    
    energia <- numeric(ns)
    
    entropia <- numeric(ns)
    
    analiticE <- numeric(ns)
    
    erroE <- numeric(ns)
    
    
    
    
    #ddalfa <- lseq(dalfa, alfa, ns)
    ddalfa <- seq(alfa +dalfa, alfa+ ns * dalfa, dalfa)
    
    
    
    for(j in 1:ns)
    {
        
        
        ss<- numeric(nx)
        
        ee<- numeric(nx)
        #momentose<-numeric(nx)
        #cumulantese<-numeric(nx)
        
        for(i in 1:nns)
        {
            listasim<-SimulacaoBimodalTcte(nx, nt, nc, dt, massa, gamma, ctelastica, ddalfa[j], 0, temp)
            keru<-analiseentropicaTctecontmomen(listasim)
            ss<-ss+(keru$S)/nns
            ee<-ee+(keru$E)/nns
            #momentose<-(keru$momentationE)/nns
            #cumulatesE<-(keru$cumulationE)/nns
        }
        
        sss<-mean(ss)
        eee<-mean(ee)
        
        anaE<-keru$PolyE[1]
        
        erroE[j]<-((eee-anaE)/anaE)*100
        
        entropia[j]<-sss
        
        energia[j]<-eee
        
        analiticE[j]<-anaE
        
        
    }
    
    Slm<-loess(entropia~ddalfa)
    
    bb<-predict(Slm,ddalfa, se = TRUE)
    
    B<-bb$fit
    
    
    #FIM DO SIMULACAO BIMODAL T CTE
    #   final<-list("posicao"=xx,"velocidade"=vv, "energia"=ee, "massa"=massa,"gamma"=gamma,"k"=ctelastica, "alfa" = alfa,"dalfa"= dalfa, "xa"=xa, "ns"=ns,"nt"=nt,"nc"=nc,"dt"= dt, "T"=Temp)
    #return(final)
    
    final<-list("PolyE"=analiticE, "E"=energia, "S"=entropia, "PolyS"= B , "alfa"=ddalfa, "nns" = nns, "ns"=ns, "nt"=nt, "nc"=nc, "dt"= dt, "temp" = temp, "erroE" = erroE)
    
    return(final)
    
}





SimAmortizadoraTctecontlogKurt<- function(nns, nx, ns, nt, nc, dt, massa, gamma, ctelastica, alfa, dalfa, temp){
    
    
    anaE <- numeric(ns)
    
    energia <- numeric(ns)
    
    entropia <- numeric(ns)
    
    analiticE <- numeric(ns)
    
    erroE <- numeric(ns)
    
    ddalfa <- lseq(dalfa, alfa, ns)
    
    
    
    #E
    
    momentationE<- vector("list", ns)
    
    cumulationE<- vector("list", ns)
    
    #X
    
    momentationX<- vector("list",ns)
    
    cumulationX<- vector("list", ns)
    
    #V
    
    momentationV<- vector("list", ns)
    
    cumulationV<- vector("list", ns)
    
    
    




    for(j in 1:ns)
    {
        
        momentosE<-numeric(7)
        cumulantesE<-numeric(5)
        momentosX<-numeric(7)
        cumulantesX<-numeric(5)
        momentosV<-numeric(7)
        cumulantesV<-numeric(5)
        
        
        ss<- numeric(nx)
        
        ee<- numeric(nx)
        #momentose<-numeric(nx)
        #cumulantese<-numeric(nx)
        
        for(i in 1:nns)
        {
            listasim<-SimulacaoBimodalTcte(nx, nt, nc, dt, massa, gamma, ctelastica, ddalfa[j], 0, temp)
            keru<-analiseentropicaTctecontmomen(listasim)
            ss<-ss+(keru$S)/nns
            ee<-ee+(keru$E)/nns
            
            #print(momentosE)
            #print(cumulantesE)
            #print(keru$momentosE[[1]])
            
            momE<- numeric(length( keru$momentosE[[1]] ))
            kumE<- numeric(length(keru$cumulantesE[[1]] ))
            momX<- numeric(length(keru$momentosX[[1]] ))
            kumX<- numeric(length(keru$cumulantesX[[1]] ))
            momV<- numeric(length(keru$momentosV[[1]] ))
            kumV<- numeric(length(keru$cumulantesV[[1]] ))
            
            for(l in 1:nx)
            {
                momE<- momE + keru$momentosE[[l]]/nx
                kumE<- kumE + keru$cumulantesE[[l]]/nx
                momX<- momX + keru$momentosX[[l]]/nx
                kumX<- kumX + keru$cumulantesX[[l]]/nx
                momV<- momV + keru$momentosV[[l]]/nx
                kumV<- kumV + keru$cumulantesV[[l]]/nx
            }
            
            
            momentosE <-momentosE+(momE)/nns
            #print(momentosE)
            cumulantesE<-cumulantesE+(kumE)/nns
            momentosX<-momentosX+(momX)/nns
            cumulantesX<-cumulantesX+(kumX)/nns
            momentosV<-momentosV+(momV)/nns
            cumulantesV<-cumulantesV+(kumV)/nns
        }
        
        sss<-mean(ss)
        eee<-mean(ee)
        
        anaE<-keru$PolyE[1]
        
        erroE[j]<-((eee-anaE)/anaE)*100
        
        entropia[j]<-sss
        
        energia[j]<-eee
        
        analiticE[j]<-anaE
        
        momentationE[[j]]<-momentosE
        cumulationE[[j]]<-cumulantesE
        momentationX[[j]]<-momentosX
        cumulationX[[j]]<-cumulantesX
        momentationV[[j]]<-momentosV
        cumulationV[[j]]<-cumulantesV
        
        print(j)
        
    }
    
    Slm<-loess(entropia~ddalfa)
    
    bb<-predict(Slm,ddalfa, se = TRUE)
    
    B<-bb$fit
    
    
    #FIM DO SIMULACAO BIMODAL T CTE
    #   final<-list("posicao"=xx,"velocidade"=vv, "energia"=ee, "massa"=massa,"gamma"=gamma,"k"=ctelastica, "alfa" = alfa,"dalfa"= dalfa, "xa"=xa, "ns"=ns,"nt"=nt,"nc"=nc,"dt"= dt, "T"=Temp)
    #return(final)
    
    final<-list("PolyE"=analiticE, "E"=energia, "S"=entropia, "PolyS"= B , "alfa"=ddalfa, "nns" = nns, "ns"=ns, "nt"=nt, "nc"=nc, "dt"= dt, "temp" = temp, "erroE" = erroE,"momentosE"=momentationE, "cumulantesE"=cumulationE, "momentosV"=momentationV, "cumulantesV"=cumulationV, "momentosX"=momentationX, "cumulantesX"=cumulationX )
    
    return(final)
    
}


#variando T diretamente

SimAmortizadoraAlfactecontKurtT<- function(nns, nx, ns, nt, nc, dt, massa, gamma, ctelastica, alfa,  temp, dtemp){
    
    # nns<-10
    #nx<-10
    #ns<-10
    #nt<-100
    #nc<-100
    #dt<-0.001
    #massa<-1
    #gamma<-1
    #ctelastica<-1
    #alfa<-0
    #dalfa<-0.01
    #temp<-1
    #
    #
    # zoado... nx nao aumentou a precisao do valor da energia, nc aumentou, preciso testar o nns
    
    xai<-numeric(ns)
    
    anaE <- numeric(ns)
    
    energia <- numeric(ns)
    
    entropia <- numeric(ns)
    
    analiticE <- numeric(ns)
    
    erroE <- numeric(ns)
    
    ddtemp <- seq(temp, temp+(ns-1)*dtemp, dtemp)
    
    #ddalfa <- seq(alfa +dalfa, alfa+ ns * dalfa, dalfa)
    
    
    #E
    
    momentationE<- vector("list", ns)
    
    cumulationE<- vector("list", ns)
    
    #X
    
    momentationX<- vector("list",ns)
    
    cumulationX<- vector("list", ns)
    
    #V
    
    momentationV<- vector("list", ns)
    
    cumulationV<- vector("list", ns)
    
    
    
    
    
    
    
    for(j in 1:ns)
    {
        
        alfai<-
        
        momentosE<-numeric(7)
        cumulantesE<-numeric(5)
        momentosX<-numeric(7)
        cumulantesX<-numeric(5)
        momentosV<-numeric(7)
        cumulantesV<-numeric(5)
        
        
        ss<- numeric(nx)
        
        ee<- numeric(nx)
        #momentose<-numeric(nx)
        #cumulantese<-numeric(nx)
        
        for(i in 1:nns)
        {
            listasim<-SimulacaoBimodalTcte(nx, nt, nc, dt, massa, gamma, ctelastica, alfa, 0, ddtemp[j])
            keru<-analiseentropicaTctecontmomen(listasim)
            ss<-ss+(keru$S)/nns
            ee<-ee+(keru$E)/nns
            
            #print(momentosE)
            #print(cumulantesE)
            #print(keru$momentosE[[1]])
            
            momE<- numeric(length( keru$momentosE[[1]] ))
            kumE<- numeric(length(keru$cumulantesE[[1]] ))
            momX<- numeric(length(keru$momentosX[[1]] ))
            kumX<- numeric(length(keru$cumulantesX[[1]] ))
            momV<- numeric(length(keru$momentosV[[1]] ))
            kumV<- numeric(length(keru$cumulantesV[[1]] ))
            
            for(l in 1:nx)
            {
                momE<- momE + keru$momentosE[[l]]/nx
                kumE<- kumE + keru$cumulantesE[[l]]/nx
                momX<- momX + keru$momentosX[[l]]/nx
                kumX<- kumX + keru$cumulantesX[[l]]/nx
                momV<- momV + keru$momentosV[[l]]/nx
                kumV<- kumV + keru$cumulantesV[[l]]/nx
            }
            
            
            momentosE <-momentosE+(momE)/nns
            #print(momentosE)
            cumulantesE<-cumulantesE+(kumE)/nns
            momentosX<-momentosX+(momX)/nns
            cumulantesX<-cumulantesX+(kumX)/nns
            momentosV<-momentosV+(momV)/nns
            cumulantesV<-cumulantesV+(kumV)/nns
        }
        
        xai[j]<-listasim$xa[1]
        
        sss<-mean(ss)
        eee<-mean(ee)
        
        anaE<-keru$PolyE[1]
        
        erroE[j]<-((eee-anaE)/anaE)*100
        
        entropia[j]<-sss
        
        energia[j]<-eee
        
        analiticE[j]<-anaE
        
        momentationE[[j]]<-momentosE
        cumulationE[[j]]<-cumulantesE
        momentationX[[j]]<-momentosX
        cumulationX[[j]]<-cumulantesX
        momentationV[[j]]<-momentosV
        cumulationV[[j]]<-cumulantesV
        
        print(j)
        
    }
    
    Slm<-loess(entropia~ddtemp)
    
    bb<-predict(Slm,ddtemp, se = TRUE)
    
    B<-bb$fit
    
    
    #FIM DO SIMULACAO BIMODAL T CTE
    #   final<-list("posicao"=xx,"velocidade"=vv, "energia"=ee, "massa"=massa,"gamma"=gamma,"k"=ctelastica, "alfa" = alfa,"dalfa"= dalfa, "xa"=xa, "ns"=ns,"nt"=nt,"nc"=nc,"dt"= dt, "T"=Temp)
    #return(final)
    
    final<-list("PolyE"=analiticE, "E"=energia, "S"=entropia, "PolyS"= B , "alfa"=alfa, "nns" = nns, "ns"=ns, "nt"=nt, "nc"=nc, "dt"= dt, "temp" = ddtemp,  "erroE" = erroE,"momentosE"=momentationE, "cumulantesE"=cumulationE, "momentosV"=momentationV, "cumulantesV"=cumulationV, "momentosX"=momentationX, "cumulantesX"=cumulationX, "xa"=xai )
    
    return(final)
    
}






SimAmortizadoraAlfactecontlogKurt<- function(nns, nx, ns, nt, nc, dt, massa, gamma, ctelastica, alfa, xa, dxa){
    
    # nns<-10
    #nx<-10
    #ns<-10
    #nt<-100
    #nc<-100
    #dt<-0.001
    #massa<-1
    #gamma<-1
    #ctelastica<-1
    #alfa<-0
    #dalfa<-0.01
    #temp<-1
    #
    #
    # zoado... nx nao aumentou a precisao do valor da energia, nc aumentou, preciso testar o nns
    
    anaE <- numeric(ns)
    
    energia <- numeric(ns)
    
    entropia <- numeric(ns)
    
    analiticE <- numeric(ns)
    
    erroE <- numeric(ns)
    
    ddxa <- lseq(xa, dxa, ns)
    
    #ddalfa <- seq(alfa +dalfa, alfa+ ns * dalfa, dalfa)
    
    
    #E
    
    momentationE<- vector("list", ns)
    
    cumulationE<- vector("list", ns)
    
    #X
    
    momentationX<- vector("list",ns)
    
    cumulationX<- vector("list", ns)
    
    #V
    
    momentationV<- vector("list", ns)
    
    cumulationV<- vector("list", ns)
    
    
    
    temperatura<-numeric(ns)
    
    
    
    for(j in 1:ns)
    {
        
        momentosE<-numeric(7)
        cumulantesE<-numeric(5)
        momentosX<-numeric(7)
        cumulantesX<-numeric(5)
        momentosV<-numeric(7)
        cumulantesV<-numeric(5)
        
        
        ss<- numeric(nx)
        
        ee<- numeric(nx)
        #momentose<-numeric(nx)
        #cumulantese<-numeric(nx)
        
        for(i in 1:nns)
        {
            listasim<-SimulacaoBimodal(nx, nt, nc, dt, massa, gamma, ctelastica, alfa, 0, ddxa[j], 0)
            keru<-analiseentropicaTctecontmomen(listasim)
            ss<-ss+(keru$S)/nns
            ee<-ee+(keru$E)/nns
            
            #print(momentosE)
            #print(cumulantesE)
            #print(keru$momentosE[[1]])
            
            
            
            momE<- numeric(length( keru$momentosE[[1]] ))
            kumE<- numeric(length(keru$cumulantesE[[1]] ))
            momX<- numeric(length(keru$momentosX[[1]] ))
            kumX<- numeric(length(keru$cumulantesX[[1]] ))
            momV<- numeric(length(keru$momentosV[[1]] ))
            kumV<- numeric(length(keru$cumulantesV[[1]] ))
            
            for(l in 1:nx)
            {
                momE<- momE + keru$momentosE[[l]]/nx
                kumE<- kumE + keru$cumulantesE[[l]]/nx
                momX<- momX + keru$momentosX[[l]]/nx
                kumX<- kumX + keru$cumulantesX[[l]]/nx
                momV<- momV + keru$momentosV[[l]]/nx
                kumV<- kumV + keru$cumulantesV[[l]]/nx
            }
            
            
            
            
            momentosE <-momentosE+(momE)/nns
            #print(momentosE)
            cumulantesE<-cumulantesE+(kumE)/nns
            momentosX<-momentosX+(momX)/nns
            cumulantesX<-cumulantesX+(kumX)/nns
            momentosV<-momentosV+(momV)/nns
            cumulantesV<-cumulantesV+(kumV)/nns
        }
        
        temperatura[j]<-listasim$temp
        
        sss<-mean(ss)
        eee<-mean(ee)
        
        anaE<-keru$PolyE[1]
        
        erroE[j]<-((eee-anaE)/anaE)*100
        
        entropia[j]<-sss
        
        energia[j]<-eee
        
        analiticE[j]<-anaE
        
        momentationE[[j]]<-momentosE
        cumulationE[[j]]<-cumulantesE
        momentationX[[j]]<-momentosX
        cumulationX[[j]]<-cumulantesX
        momentationV[[j]]<-momentosV
        cumulationV[[j]]<-cumulantesV
        
         print(j)
        
    }
    
    Slm<-loess(entropia~ddxa)
    
    bb<-predict(Slm,ddxa, se = TRUE)
    
    B<-bb$fit
    
    
    #FIM DO SIMULACAO BIMODAL T CTE
    #   final<-list("posicao"=xx,"velocidade"=vv, "energia"=ee, "massa"=massa,"gamma"=gamma,"k"=ctelastica, "alfa" = alfa,"dalfa"= dalfa, "xa"=xa, "ns"=ns,"nt"=nt,"nc"=nc,"dt"= dt, "T"=Temp)
    #return(final)
    
    final<-list("PolyE"=analiticE, "E"=energia, "S"=entropia, "PolyS"= B , "alfa"=alfa, "xa"=ddxa , "nns" = nns, "ns"=ns, "nt"=nt, "nc"=nc, "dt"= dt, "temp" = temperatura, "erroE" = erroE,"momentosE"=momentationE, "cumulantesE"=cumulationE, "momentosV"=momentationV, "cumulantesV"=cumulationV, "momentosX"=momentationX, "cumulantesX"=cumulationX )
    
    return(final)
    
}



SimAmortizadoraAlfactecontKurt<- function(nns, nx, ns, nt, nc, dt, massa, gamma, ctelastica, alfa, xa, dxa){
    
    # nns<-10
    #nx<-10
    #ns<-10
    #nt<-100
    #nc<-100
    #dt<-0.001
    #massa<-1
    #gamma<-1
    #ctelastica<-1
    #alfa<-0
    #dalfa<-0.01
    #temp<-1
    #
    #
    # zoado... nx nao aumentou a precisao do valor da energia, nc aumentou, preciso testar o nns
    
    anaE <- numeric(ns)
    
    energia <- numeric(ns)
    
    entropia <- numeric(ns)
    
    analiticE <- numeric(ns)
    
    erroE <- numeric(ns)
    
    ddxa <- seq(xa, xa+ (ns-1)*dxa , dxa)
    
    #ddalfa <- seq(alfa +dalfa, alfa+ ns * dalfa, dalfa)
    
    
    #E
    
    momentationE<- vector("list", ns)
    
    cumulationE<- vector("list", ns)
    
    #X
    
    momentationX<- vector("list",ns)
    
    cumulationX<- vector("list", ns)
    
    #V
    
    momentationV<- vector("list", ns)
    
    cumulationV<- vector("list", ns)
    
    
    
    temperatura<-numeric(ns)
    
    
    
    for(j in 1:ns)
    {
        
        momentosE<-numeric(7)
        cumulantesE<-numeric(5)
        momentosX<-numeric(7)
        cumulantesX<-numeric(5)
        momentosV<-numeric(7)
        cumulantesV<-numeric(5)
        
        
        ss<- numeric(nx)
        
        ee<- numeric(nx)
        #momentose<-numeric(nx)
        #cumulantese<-numeric(nx)
        
        for(i in 1:nns)
        {
            listasim<-SimulacaoBimodal(nx, nt, nc, dt, massa, gamma, ctelastica, alfa, 0, ddxa[j], 0)
            keru<-analiseentropicaTctecontmomen(listasim)
            ss<-ss+(keru$S)/nns
            ee<-ee+(keru$E)/nns
            
            #print(momentosE)
            #print(cumulantesE)
            #print(keru$momentosE[[1]])
            
            
            
            momE<- numeric(length( keru$momentosE[[1]] ))
            kumE<- numeric(length(keru$cumulantesE[[1]] ))
            momX<- numeric(length(keru$momentosX[[1]] ))
            kumX<- numeric(length(keru$cumulantesX[[1]] ))
            momV<- numeric(length(keru$momentosV[[1]] ))
            kumV<- numeric(length(keru$cumulantesV[[1]] ))
            
            for(l in 1:nx)
            {
                momE<- momE + keru$momentosE[[l]]/nx
                kumE<- kumE + keru$cumulantesE[[l]]/nx
                momX<- momX + keru$momentosX[[l]]/nx
                kumX<- kumX + keru$cumulantesX[[l]]/nx
                momV<- momV + keru$momentosV[[l]]/nx
                kumV<- kumV + keru$cumulantesV[[l]]/nx
            }
            
            
            
            
            momentosE <-momentosE+(momE)/nns
            #print(momentosE)
            cumulantesE<-cumulantesE+(kumE)/nns
            momentosX<-momentosX+(momX)/nns
            cumulantesX<-cumulantesX+(kumX)/nns
            momentosV<-momentosV+(momV)/nns
            cumulantesV<-cumulantesV+(kumV)/nns
        }
        
        temperatura[j]<-listasim$temp
        
        sss<-mean(ss)
        eee<-mean(ee)
        
        anaE<-keru$PolyE[1]
        
        erroE[j]<-((eee-anaE)/anaE)*100
        
        entropia[j]<-sss
        
        energia[j]<-eee
        
        analiticE[j]<-anaE
        
        momentationE[[j]]<-momentosE
        cumulationE[[j]]<-cumulantesE
        momentationX[[j]]<-momentosX
        cumulationX[[j]]<-cumulantesX
        momentationV[[j]]<-momentosV
        cumulationV[[j]]<-cumulantesV
        
        print(j)
        
    }
    
    Slm<-loess(entropia~ddxa)
    
    bb<-predict(Slm,ddxa, se = TRUE)
    
    B<-bb$fit
    
    
    #FIM DO SIMULACAO BIMODAL T CTE
    #   final<-list("posicao"=xx,"velocidade"=vv, "energia"=ee, "massa"=massa,"gamma"=gamma,"k"=ctelastica, "alfa" = alfa,"dalfa"= dalfa, "xa"=xa, "ns"=ns,"nt"=nt,"nc"=nc,"dt"= dt, "T"=Temp)
    #return(final)
    
    final<-list("PolyE"=analiticE, "E"=energia, "S"=entropia, "PolyS"= B , "alfa"=alfa, "xa"=ddxa , "nns" = nns, "ns"=ns, "nt"=nt, "nc"=nc, "dt"= dt, "temp" = temperatura, "erroE" = erroE,"momentosE"=momentationE, "cumulantesE"=cumulationE, "momentosV"=momentationV, "cumulantesV"=cumulationV, "momentosX"=momentationX, "cumulantesX"=cumulationX )
    
    return(final)
    
}



SimAmortizadoraAlfactecontKurtpra3D<- function(nns, nx, ns, nt, nc, dt, massa, gamma, ctelastica, alfa, xa, dxa){
    
    # nns<-10
    #nx<-10
    #ns<-10
    #nt<-100
    #nc<-100
    #dt<-0.001
    #massa<-1
    #gamma<-1
    #ctelastica<-1
    #alfa<-0
    #dalfa<-0.01
    #temp<-1
    #
    #
    # zoado... nx nao aumentou a precisao do valor da energia, nc aumentou, preciso testar o nns
    
    anaE <- numeric(ns)
    
    energia <- numeric(ns)
    
    entropia <- numeric(ns)
    
    analiticE <- numeric(ns)
    
    erroE <- numeric(ns)
    
    ddxa <- seq(xa, xa+ (ns-1)*dxa , dxa)
    
    #ddalfa <- seq(alfa +dalfa, alfa+ ns * dalfa, dalfa)
    
    
    #E
    
    momentationE<- vector("list", ns)
    
    cumulationE<- vector("list", ns)
    
    #X
    
    momentationX<- vector("list",ns)
    
    cumulationX<- vector("list", ns)
    
    #V
    
    momentationV<- vector("list", ns)
    
    cumulationV<- vector("list", ns)
    
    
    
    temperatura<-numeric(ns)
    
    
    
    for(j in 1:ns)
    {
        
        momentosE<-numeric(7)
        cumulantesE<-numeric(5)
        momentosX<-numeric(7)
        cumulantesX<-numeric(5)
        momentosV<-numeric(7)
        cumulantesV<-numeric(5)
        
        
        ss<- numeric(nx)
        
        ee<- numeric(nx)
        #momentose<-numeric(nx)
        #cumulantese<-numeric(nx)
        
        for(i in 1:nns)
        {
            listasim<-SimulacaoBimodal(nx, nt, nc, dt, massa, gamma, ctelastica, alfa, 0, ddxa[j], 0)
            keru<-analiseentropicaTctecontmomen(listasim)
            ss<-ss+(keru$S)/nns
            ee<-ee+(keru$E)/nns
            
            #print(momentosE)
            #print(cumulantesE)
            #print(keru$momentosE[[1]])
            
            
            
            momE<- numeric(length( keru$momentosE[[1]] ))
            kumE<- numeric(length(keru$cumulantesE[[1]] ))
            momX<- numeric(length(keru$momentosX[[1]] ))
            kumX<- numeric(length(keru$cumulantesX[[1]] ))
            momV<- numeric(length(keru$momentosV[[1]] ))
            kumV<- numeric(length(keru$cumulantesV[[1]] ))
            
            for(l in 1:nx)
            {
                momE<- momE + keru$momentosE[[l]]/nx
                kumE<- kumE + keru$cumulantesE[[l]]/nx
                momX<- momX + keru$momentosX[[l]]/nx
                kumX<- kumX + keru$cumulantesX[[l]]/nx
                momV<- momV + keru$momentosV[[l]]/nx
                kumV<- kumV + keru$cumulantesV[[l]]/nx
            }
            
            
            
            
            momentosE <-momentosE+(momE)/nns
            #print(momentosE)
            cumulantesE<-cumulantesE+(kumE)/nns
            momentosX<-momentosX+(momX)/nns
            cumulantesX<-cumulantesX+(kumX)/nns
            momentosV<-momentosV+(momV)/nns
            cumulantesV<-cumulantesV+(kumV)/nns
        }
        
        temperatura[j]<-listasim$temp
        
        sss<-mean(ss)
        eee<-mean(ee)
        
        anaE<-keru$PolyE[1]
        
        erroE[j]<-((eee-anaE)/anaE)*100
        
        entropia[j]<-sss
        
        energia[j]<-eee
        
        analiticE[j]<-anaE
        
        momentationE[[j]]<-momentosE
        cumulationE[[j]]<-cumulantesE
        momentationX[[j]]<-momentosX
        cumulationX[[j]]<-cumulantesX
        momentationV[[j]]<-momentosV
        cumulationV[[j]]<-cumulantesV
        
        #print(j)
        
    }
    
    # Slm<-loess(entropia~ddxa)
    
    # bb<-predict(Slm,ddxa, se = TRUE)
    
    # B<-bb$fit
    
    
    #FIM DO SIMULACAO BIMODAL T CTE
    #   final<-list("posicao"=xx,"velocidade"=vv, "energia"=ee, "massa"=massa,"gamma"=gamma,"k"=ctelastica, "alfa" = alfa,"dalfa"= dalfa, "xa"=xa, "ns"=ns,"nt"=nt,"nc"=nc,"dt"= dt, "T"=Temp)
    #return(final)
    
    final<-list("PolyE"=analiticE, "E"=energia, "S"=entropia,  "alfa"=alfa, "xa"=ddxa , "nns" = nns, "ns"=ns, "nt"=nt, "nc"=nc, "dt"= dt, "temp" = temperatura, "erroE" = erroE,"momentosE"=momentationE, "cumulantesE"=cumulationE, "momentosV"=momentationV, "cumulantesV"=cumulationV, "momentosX"=momentationX, "cumulantesX"=cumulationX )
    
    return(final)
    
}

SimSimplisLangevin<- function(nns, nx, ns, nt, nc, dt, massa, gamma, ctelastica, temp){
    
    
    anaE <- numeric(ns)
    
    energia <- numeric(ns)
    
    entropia <- numeric(ns)
    
    analiticE <- numeric(ns)
    
    erroE <- numeric(ns)
    
    #  ddxa <- seq(xa, xa+ (ns-1)*dxa , dxa)
    
    #ddalfa <- seq(alfa +dalfa, alfa+ ns * dalfa, dalfa)
    
    
    
    
    temperatura<-numeric(ns)
    
    
    
    for(j in 1:ns)
    {
        
        
        
        ss<- numeric(nx)
        
        ee<- numeric(nx)
        #momentose<-numeric(nx)
        #cumulantese<-numeric(nx)
        
        for(i in 1:nns)
        {
            
            #SimulacaoLangevin <- function(ns, nt, nc, dt, massa, gamma, ctelastica, temp){

            
            listasim<-SimulacaoLangevin(nx, nt, nc, dt, massa, gamma, ctelastica,temp)
            keru<-analiseentropicaLangevin(listasim)
            
            ss<-ss+(keru$S)/nns
            ee<-ee+(keru$E)/nns
            
        }
        
        temperatura[j]<-listasim$temp
        
        sss<-mean(ss)
        eee<-mean(ee)
        
        anaE<-temp
        
        erroE[j]<-((eee-anaE)/anaE)*100
        
        entropia[j]<-sss
        
        energia[j]<-eee
        
        analiticE[j]<-anaE
        
        
        #print(j)
        
    }
    
    
    
    
    
    
    
    final<-list("PolyE"=analiticE, "E"=energia, "S"=entropia , "nns" = nns, "ns"=ns, "nt"=nt, "nc"=nc, "dt"= dt, "temp" = temperatura, "erroE" = erroE)
    
    return(final)
    
}


SimSimplesAlfactepra3D<- function(nns, nx, ns, nt, nc, dt, massa, gamma, ctelastica, alfa, xa, dxa){
    
    
    anaE <- numeric(ns)
    energia <- numeric(ns)
    entropia <- numeric(ns)
    analiticE <- numeric(ns)
    erroE <- numeric(ns)
    ddxa <- seq(xa, xa+ (ns-1)*dxa , dxa)
    temperatura<-numeric(ns)
    
    for(j in 1:ns)
    {
        
        ss<- numeric(nx)
        ee<- numeric(nx)
        #momentose<-numeric(nx)
        #cumulantese<-numeric(nx)
        
        for(i in 1:nns)
        {
            listasim<-SimulacaoBimodal(nx, nt, nc, dt, massa, gamma, ctelastica, alfa, 0, ddxa[j], 0)
            keru<-analiseentropicaSIMPLIS(listasim)
            
            ss<-ss+(keru$S)/nns
            ee<-ee+(keru$E)/nns
            
        }
        
        temperatura[j]<-listasim$temp
        
        sss<-mean(ss)
        eee<-mean(ee)
        
        anaE<-keru$PolyE[1]
        
        erroE[j]<-((eee-anaE)/anaE)*100
        
        entropia[j]<-sss
        
        energia[j]<-eee
        
        analiticE[j]<-anaE
        
        
        #print(j)
        
    }


    final<-list("PolyE"=analiticE, "E"=energia, "S"=entropia,  "alfa"=alfa, "xa"=ddxa , "nns" = nns, "ns"=ns, "nt"=nt, "nc"=nc, "dt"= dt, "temp" = temperatura, "erroE" = erroE)
    
    return(final)
    
}

SimSimplesAlfactepra3DOD<- function(nns, nx, ns, nt, nc, dt, massa, gamma, ctelastica, alfa, xa, dxa){
    
    anaE <- numeric(ns)
    energia <- numeric(ns)
    entropia <- numeric(ns)
    analiticE <- numeric(ns)
    erroE <- numeric(ns)
    ddxa <- seq(xa, xa+ (ns-1)*dxa , dxa)
    temperatura<-numeric(ns)
    
    for(j in 1:ns)
    {
        ss<- numeric(nx)
        ee<- numeric(nx)
        
        for(i in 1:nns)
        {
            listasim<-SimulacaoBimodalOD(nx, nt, nc, dt, massa, gamma, ctelastica, alfa, 0, ddxa[j], 0)
            keru<-analiseentropicaSIMPLIS(listasim)
            
            ss<-ss+(keru$S)/nns
            ee<-ee+(keru$E)/nns
            
        }
        
        temperatura[j]<-listasim$temp
        
        sss<-mean(ss)
        eee<-mean(ee)
        
        anaE<-keru$PolyE[1]
        
        erroE[j]<-((eee-anaE)/anaE)*100
        
        entropia[j]<-sss
        
        energia[j]<-eee
        
        analiticE[j]<-anaE
        
        
        #print(j)
        
    }
    
    
    final<-list("PolyE"=analiticE, "E"=energia, "S"=entropia,  "alfa"=alfa, "xa"=ddxa , "nns" = nns, "ns"=ns, "nt"=nt, "nc"=nc, "dt"= dt, "temp" = temperatura, "erroE" = erroE)
    
    return(final)
    
}

SimSimplesEctepra3D<- function(nns, nx, ns, nt, nc, dt, massa, gamma, ctelastica, alfa, xa, dxa){
    
    
    anaE <- numeric(ns)
    energia <- numeric(ns)
    entropia <- numeric(ns)
    analiticE <- numeric(ns)
    erroE <- numeric(ns)
    ddxa <- seq(xa, xa+ (ns-1)*dxa , dxa)
    temperatura<-numeric(ns)
    
    
    xa<- (sqrt(2* energia))/sqrt(   1/(ctelastica + massa*alfa**2 + alfa*gamma)   +    (2*massa*alfa)/(ctelastica + massa*alfa**2 + alfa*gamma)   )
    
    for(j in 1:ns)
    {
        ss<- numeric(nx)
        ee<- numeric(nx)
        
        for(i in 1:nns)
        {
            listasim<-SimulacaoBimodal(nx, nt, nc, dt, massa, gamma, ctelastica, alfa, 0, ddxa[j], 0)
            keru<-analiseentropicaSIMPLIS(listasim)
            ss<-ss+(keru$S)/nns
            ee<-ee+(keru$E)/nns
        }
        
        temperatura[j]<-listasim$temp
        
        sss<-mean(ss)
        eee<-mean(ee)
        anaE<-keru$PolyE[1]
        
        erroE[j]<-((eee-anaE)/anaE)*100
        entropia[j]<-sss
        energia[j]<-eee
        analiticE[j]<-anaE
        
    }
    
    final<-list("PolyE"=analiticE, "E"=energia, "S"=entropia,  "alfa"=alfa, "xa"=ddxa , "nns" = nns, "ns"=ns, "nt"=nt, "nc"=nc, "dt"= dt, "temp" = temperatura, "erroE" = erroE)
    
    return(final)
    
}

testeEcte <- function(nome ){
    
    #qero <- SimSimplesEcte(7, 7, 20, 10000, 100, 0.001, 1, 1, 1, 0.2, 0.3, 0.3 )
    qero <- SimSimplesEcte(7, 7, 20, 10000, 100, 0.001, 1, 1, 2, 0.2, 0.3, 0.3 )
    
    #function(nns, nx, ns, nt, nc, dt, massa, gamma, ctelastica, alfa, dalfa, energia ){
    
    #qero2 <- SimSimplesEcte(7, 7, 20, 10000, 100, 0.001, 1, 1, 1, 0.25, 0.3, 0.3 )
    
    
    SalvaCorrelacao(qero$E, paste( nome, "energia.csv", sep = "" ) )
    SalvaCorrelacao(qero$S, paste( nome, "entropia.csv", sep ="" ) )
    
    #SalvaCorrelacao(qero$E, paste( nome, "energia.csv", sep = "" ) )
    #SalvaCorrelacao(qero$S, paste( nome, "entropia.csv", sep ="" ) )
    
    SalvaCorrelacao(qero$alfa, paste( nome, "alfa.csv", sep = "" ) )
    SalvaCorrelacao(qero$erroE, paste( nome, "erroE.csv", sep = "" ) )
    
}

testeTcte <- function(nome ){
    
    #qero <- SimSimplesTcte(7, 7, 5, 10000, 100, 0.001, 1, 1, 1, 0.5, 0.5, 2.2 )
    
    qero <- SimSimplesTcte(7, 7, 15, 10000, 100, 0.001, 1, 1, 2, 0.5, 0.5, 2.2 )
    
    s<-qero$S
    e<-qero$S
    alfa<-qero$alfa
    dsde <- numeric(length(alfa))
    
    for(i in 2:length(alfa))
    {
        dsde[i] <- (s[i]-s[i-1])/(e[i]-e[i-1])
    }
    
    SalvaCorrelacao(qero$E, paste( nome, "energia.csv", sep = "" ) )
    SalvaCorrelacao(qero$S, paste( nome, "entropia.csv", sep ="" ) )
    SalvaCorrelacao(qero$alfa, paste( nome, "alfa.csv", sep = "" ) )
    SalvaCorrelacao(qero$erroE, paste( nome, "erroE.csv", sep = "" ) )
    SalvaCorrelacao(dsde, paste( nome, "dsde.csv", sep = "" ) )
    
}

SimSimplesEcte<- function(nns, nx, ns, nt, nc, dt, massa, gamma, ctelastica, alfa, dalfa, energia ){
    
    
    anaE <- numeric(ns)
    energiaa <- numeric(ns)
    entropia <- numeric(ns)
    analiticE <- numeric(ns)
    erroE <- numeric(ns)
    
    
    alfai <- seq(alfa, alfa+ (ns-1)*dalfa , dalfa)
    temperatura<-numeric(ns)
    
    
    xa<- (sqrt(2* energia))/sqrt(  1/(ctelastica + massa*alfai**2 + alfai*gamma)   +    (2*massa*alfai)/(ctelastica + massa*alfai**2 + alfai*gamma)   )
    
    for(j in 1:ns)
    {
        print(j)
        ss<- numeric(nx)
        ee<- numeric(nx)
        
        for(i in 1:nns)
        {
            listasim<-SimulacaoBimodal(nx, nt, nc, dt, massa, gamma, ctelastica, alfai[j], 0, xa[j], 0)
            keru<-analiseentropicaSIMPLIS(listasim)
            ss<-ss+(keru$S)/nns
            ee<-ee+(keru$E)/nns
        }
        
        temperatura[j]<-listasim$temp
        
        sss<-mean(ss)
        eee<-mean(ee)
        anaE<-keru$PolyE[1]
        
        erroE[j]<-((eee-anaE)/anaE)*100
        entropia[j]<-sss
        energiaa[j]<-eee
        analiticE[j]<-anaE
        
    }
    
    final<-list("PolyE"=analiticE, "E"=energiaa, "S"=entropia,  "alfa"=alfai, "xa"=xa , "nns" = nns, "ns"=ns, "nt"=nt, "nc"=nc, "dt"= dt, "temp" = temperatura, "erroE" = erroE)
    
    return(final)
    
}

SimSimplesTcte<- function(nns, nx, ns, nt, nc, dt, massa, gamma, ctelastica, alfa, dalfa, Temp ){
    
    
    anaE <- numeric(ns)
    energiaa <- numeric(ns)
    entropia <- numeric(ns)
    analiticE <- numeric(ns)
    erroE <- numeric(ns)
    
    
    alfai <- seq(alfa, alfa+ (ns-1)*dalfa , dalfa)
    temperatura<-numeric(ns)
    
    
    #xa<- (sqrt(2* energia))/sqrt(  1/(ctelastica + massa*alfai**2 + alfai*gamma)   +    (2*massa*alfai)/(ctelastica + massa*alfai**2 + alfai*gamma)   )
    
    xa<- sqrt( Temp*(gamma*(ctelastica+alfai*(gamma+massa*alfai)))/(massa*alfai)   )
    
    for(j in 1:ns)
    {
        print(j)
        ss<- numeric(nx)
        ee<- numeric(nx)
        
        for(i in 1:nns)
        {
            listasim<-SimulacaoBimodal(nx, nt, nc, dt, massa, gamma, ctelastica, alfai[j], 0, xa[j], 0)
            keru<-analiseentropicaSIMPLIS(listasim)
            ss<-ss+(keru$S)/nns
            ee<-ee+(keru$E)/nns
        }
        
        temperatura[j]<-listasim$temp
        
        sss<-mean(ss)
        eee<-mean(ee)
        anaE<-keru$PolyE[1]
        
        erroE[j]<-((eee-anaE)/anaE)*100
        entropia[j]<-sss
        energiaa[j]<-eee
        analiticE[j]<-anaE
        
    }
    
    final<-list("PolyE"=analiticE, "E"=energiaa, "S"=entropia,  "alfa"=alfai, "xa"=xa , "nns" = nns, "ns"=ns, "nt"=nt, "nc"=nc, "dt"= dt, "temp" = temperatura, "erroE" = erroE)
    
    return(final)
    
}



SimSimplesXactepra3D<- function(nns, nx, ns, nt, nc, dt, massa, gamma, ctelastica, alfa, dalfa, xa){
    
    
    anaE <- numeric(ns)
    
    energia <- numeric(ns)
    
    entropia <- numeric(ns)
    
    analiticE <- numeric(ns)
    
    erroE <- numeric(ns)
    
    #ddxa <- seq(xa, xa+ (ns-1)*dxa , dxa)
    
    ddalfa <- seq(alfa +dalfa, alfa+ ns * dalfa, dalfa)
    
    
    
    
    temperatura<-numeric(ns)
    
    
    
    for(j in 1:ns)
    {
        
        
        
        ss<- numeric(nx)
        
        ee<- numeric(nx)
        #momentose<-numeric(nx)
        #cumulantese<-numeric(nx)
        
        for(i in 1:nns)
        {
            listasim<-SimulacaoBimodal(nx, nt, nc, dt, massa, gamma, ctelastica, ddalfa[j], 0, xa, 0)
            keru<-analiseentropicaSIMPLIS(listasim)
            
            ss<-ss+(keru$S)/nns
            ee<-ee+(keru$E)/nns
            
        }
        
        temperatura[j]<-listasim$temp
        
        sss<-mean(ss)
        eee<-mean(ee)
        
        anaE<-keru$PolyE[1]
        
        erroE[j]<-((eee-anaE)/anaE)*100
        
        entropia[j]<-sss
        
        energia[j]<-eee
        
        analiticE[j]<-anaE
        
        
        #print(j)
        
    }
    
    
    
    
    
    
    
    final<-list("PolyE"= analiticE, "E"=energia, "S"=entropia,  "alfa"=ddalfa, "xa"=xa , "nns" = nns, "ns"=ns, "nt"=nt, "nc"=nc, "dt"= dt, "temp" = temperatura, "erroE" = erroE)
    
    return(final)
    
}


SimSimplesTctepra3D<- function(nns, nx, ns=3, nt, nc, dt, massa, gamma, ctelastica, alfa, xa, dalfa){
    
    
    anaE <- numeric(ns)
    
    energia <- numeric(ns)
    
    entropia <- numeric(ns)
    
    analiticE <- numeric(ns)
    
    erroE <- numeric(ns)
    
    #ddxa <- seq(xa, xa+ (ns-1)*dxa , dxa)
    
    Temp<- ( alfa*massa* xa**2 )/ (gamma *(ctelastica + alfa * (gamma + massa*alfa) ) )
    
    ddalfa <- seq(alfa +dalfa, alfa+ ns * dalfa, dalfa)
    
    ddxa<- sqrt( Temp*(gamma*(ctelastica+ddalfa*(gamma+massa*ddalfa)))/(massa*ddalfa)   )

    
    temperatura<-numeric(ns)
    
    
    
    for(j in 1:ns)
    {
        #print(j)
        
        
        ss<- numeric(nx)
        
        ee<- numeric(nx)
        #momentose<-numeric(nx)
        #cumulantese<-numeric(nx)
        
        for(i in 1:nns)
        {
            #print(i)
            
            #SimulacaoBimodalTcte <- function(ns, nt, nc, dt, massa, gamma, ctelastica, alfa, dalfa, Temp)
            
            listasim<-SimulacaoBimodalTcte(nx, nt, nc, dt, massa, gamma, ctelastica, ddalfa[j] , 0, Temp)
            keru<-analiseentropicaSIMPLISTT(listasim)
            
            ss<-ss+(keru$S)/nns
            ee<-ee+(keru$E)/nns
            
        }
        
        temperatura[j]<-listasim$T
        
        sss<-mean(ss)
        eee<-mean(ee)
        
        anaE<-keru$PolyE[1]
        
        erroE[j]<-((eee-anaE)/anaE)*100
        
        entropia[j]<-sss
        
        energia[j]<-eee
        
        analiticE[j]<-anaE
        
        
        #print(j)
        
    }
    
    
    final<-list("PolyE"=analiticE, "E"=energia, "S"=entropia,  "alfa"=alfa, "xa"=ddxa , "nns" = nns, "ns"=ns, "nt"=nt, "nc"=nc, "dt"= dt, "temp" = temperatura, "erroE" = erroE)
    
    return(final)
    
}


SimSimplesTctepra3D2<- function(nns, nx, ns=3, nt, nc, dt, massa, gamma, ctelastica, alfa, xa, dxa){
    
    
    anaE <- numeric(ns)
    
    energia <- numeric(ns)
    
    entropia <- numeric(ns)
    
    analiticE <- numeric(ns)
    
    erroE <- numeric(ns)
    
    ddxa <- seq(xa, xa+ (ns-1)*dxa , dxa)
    
    Temp<- ( alfa*massa* xa**2 )/ (gamma *(ctelastica + alfa * (gamma + massa*alfa) ) )
    
    #ddalfa <- seq(alfa +dalfa, alfa+ ns * dalfa, dalfa)
    #ddxa<- sqrt( Temp*(gamma*(ctelastica+ddalfa*(gamma+massa*ddalfa)))/(massa*ddalfa)   )
    
    ddalfa <- ( massa*ddxa**2 - Temp * gamma**2 + sqrt( -4*ctelastica*massa*(Temp **2)*(gamma**2) + ((ddxa**2)*massa - Temp*(gamma**2) )**2 ) )/(2*massa*Temp*gamma)
    
    temperatura<-numeric(ns)
    
    for(j in 1:ns)
    {
        
        ss<- numeric(nx)
        
        ee<- numeric(nx)
        
        for(i in 1:nns)
        {
            
            listasim<-SimulacaoBimodalTcte(nx, nt, nc, dt, massa, gamma, ctelastica, ddalfa[j] , 0, Temp)
            keru<-analiseentropicaSIMPLISTT(listasim)
            
            ss<-ss+(keru$S)/nns
            ee<-ee+(keru$E)/nns
            
        }
        
        temperatura[j]<-listasim$T
        
        sss<-mean(ss)
        eee<-mean(ee)
        
        anaE<-keru$PolyE[1]
        
        erroE[j]<-((eee-anaE)/anaE)*100
        
        entropia[j]<-sss
        
        energia[j]<-eee
        
        analiticE[j]<-anaE
        
        
    }
    
    
    final<-list("PolyE"=analiticE, "E"=energia, "S"=entropia,  "alfa"=alfa, "xa"=ddxa , "nns" = nns, "ns"=ns, "nt"=nt, "nc"=nc, "dt"= dt, "temp" = temperatura, "erroE" = erroE)
    
    return(final)
    
}




D3graf<-function(nns, nx, ns, nsa, nt, nc, dt, massa, gamma, ctelastica, alfa, dalfa, xa, dxa ,deltaxa ){
    
    ddalfa <- seq(alfa, alfa+ (nsa-1)*dalfa , dalfa)
    
    ddxa <- seq(xa, xa+ (ns-1)*deltaxa , deltaxa)
    
    dsde<-matrix( nrow = nsa, ncol = ns)
    
    erroE<-matrix( nrow = nsa, ncol = ns)
    
    tempi<-matrix( nrow = nsa, ncol = ns)
    
    
    tempii<-matrix( nrow = nsa, ncol = ns)
    
    
    si<-matrix( nrow = nsa, ncol = ns)
    
    sii<-matrix( nrow = nsa, ncol = ns)
    
    ei<-matrix( nrow = nsa, ncol = ns)
    
    eii<-matrix( nrow = nsa, ncol = ns)
    
    
    
    # momentationE<- vector("list", ns)
    
    # cumulationE<- vector("list", ns)
    
    #X
    
    #  momentationX<- vector("list",ns)
    
    #  cumulationX<- vector("list", ns)
    
    #V
    
    #  momentationV<- vector("list", ns)
    
    # cumulationV<- vector("list", ns)

    
    
    for(l in 1:nsa){
        
        #print(l)
    
    for(j in 1:ns){
        
        
        qerin<-SimSimplesAlfactepra3D(nns, nx, 2, nt, nc, dt, massa, gamma, ctelastica, ddalfa[l] , ddxa[j] , dxa)
        
        tempi[l,j]<-qerin$temp[1]
        
        tempii[l,j]<-qerin$temp[2]
        
        si[l,j] <- qerin$S[1]
        
        sii[l,j] <- qerin$S[2]
        
        ei[l,j] <- qerin$E[1]
        
        eii[l,j] <- qerin$E[2]
        
        dsde[l,j]<- (qerin$S[2]-qerin$S[1])/(qerin$E[2]-qerin$E[1])
        
        erroE[l,j]<-mean(qerin$erroE)
        
        
    }
    
    
    #momentationE[[j]]<-momentosE
    #cumulationE[[j]]<-cumulantesE
    #momentationX[[j]]<-momentosX
    #cumulationX[[j]]<-cumulantesX
    #momentationV[[j]]<-momentosV
    #cumulationV[[j]]<-cumulantesV



    print(l)
    
    
    }
    
    
    final<-list("ei"=ei, "eii"=eii  ,"si"= si, "sii"= sii , "alfa"=ddalfa, "xa"=ddxa , "nns" = nns, "ns"=ns, "nsa"=nsa, "nt"=nt, "nc"=nc, "dt"= dt, "tempi" = tempi, "tempii" = tempii, "erroE" = erroE, "dsde"=dsde)
    #"momentosE"=momentationE, "cumulantesE"=cumulationE, "momentosV"=momentationV, "cumulantesV"=cumulationV, "momentosX"=momentationX, "cumulantesX"=cumulationX )
    
    return(final)

    
}

#calcula a derivada a esquerda e a direita

D3grafRL<-function(nns, nx, ns, nsa, nt, nc, dt, massa, gamma, ctelastica, alfa, dalfa, xa, dxa ,deltaxa ){
    
    
    
    
    ddalfa <- seq(alfa, alfa+ (nsa-1)*dalfa , dalfa)
    
    ddxa <- seq(xa, xa+ (ns-1)*deltaxa , deltaxa)
    
    dsde<-matrix( nrow = nsa, ncol = ns)
    
    erroE<-matrix( nrow = nsa, ncol = ns)
    
    
    
    tempi<-matrix( nrow = nsa, ncol = ns)
    
    tempr<-matrix( nrow = nsa, ncol = ns)
    
    templ<-matrix( nrow = nsa, ncol = ns)
    
    
    
    si<-matrix( nrow = nsa, ncol = ns)
    
    sr<-matrix( nrow = nsa, ncol = ns)
    
    sl<-matrix( nrow = nsa, ncol = ns)
    
    
    
    ei<-matrix( nrow = nsa, ncol = ns)
    
    er<-matrix( nrow = nsa, ncol = ns)
    
    el<-matrix( nrow = nsa, ncol = ns)
    
    
    
    
    for(l in 1:nsa){
        
        for(j in 1:ns){
            
            qerin <- SimSimplesAlfactepra3D(nns, nx, 3, nt, nc, dt, massa, gamma, ctelastica, ddalfa[l] , ddxa[j] - dxa , dxa)
            
            
            templ[l,j] <- qerin$temp[1]
            
            tempi[l,j] <- qerin$temp[2]
            
            tempr[l,j] <- qerin$temp[2]
            
            
            sl[l,j] <- qerin$S[1]
            
            si[l,j] <- qerin$S[2]
            
            sr[l,j] <- qerin$S[3]
            
            
            el[l,j] <- qerin$E[1]
            
            ei[l,j] <- qerin$E[2]
                        
            er[l,j] <- qerin$E[3]
            
            
            
            dsde[l,j]<- (qerin$S[3]-qerin$S[1])/(qerin$E[3]-qerin$E[1])
            
            erroE[l,j]<-mean(qerin$erroE)
            
            
        }
        
        
        
        
        print(l)
        
        
    }
    
    
    final<-list("ei"=ei, "er"=er  , "el"=el  ,"si"= si, "sii"= sr , "sii"= sl , "alfa"=ddalfa, "xa"=ddxa , "nns" = nns, "ns"=ns, "nsa"=nsa, "nt"=nt, "nc"=nc, "dt"= dt, "tempi" = tempi, "tempr" = tempr, "templ" = templ, "erroE" = erroE, "dsde"=dsde)
    #"momentosE"=momentationE, "cumulantesE"=cumulationE, "momentosV"=momentationV, "cumulantesV"=cumulationV, "momentosX"=momentationX, "cumulantesX"=cumulationX )
    
    return(final)
    
    
}


D3grafLite<-function(nns, nx, ns, nsa, nt, nc, dt, massa, gamma, ctelastica, alfa, dalfa, xa, dxa  ){
    
    ddalfa <- seq(alfa, alfa+ (nsa-1)*dalfa , dalfa)
    
    ddxa <- seq(xa, xa+ (ns-1)*dxa , dxa)
    dsde <- matrix( nrow = nsa, ncol = ns)
    erroE <- matrix( nrow = nsa, ncol = ns)
    
    tempi <- matrix( nrow = nsa, ncol = ns)
    si <- matrix( nrow = nsa, ncol = ns)
    ei <- matrix( nrow = nsa, ncol = ns)

    for(l in 1:nsa){
        for(j in 1:ns){
            
            qerin <- SimSimplesAlfactepra3D(nns, nx, 1, nt, nc, dt, massa, gamma, ctelastica, ddalfa[l] , ddxa[j]  , 0)
            tempi[l,j] <- qerin$temp[1]
            si[l,j] <- qerin$S[1]
            ei[l,j] <- qerin$E[1]
            erroE[l,j]<-mean(qerin$erroE)

        }
        print(l)
    }
    
    final<-list("ei"=ei, "si"= si, "alfa"=ddalfa, "xa"=ddxa , "nns" = nns, "ns"=ns, "nsa"=nsa, "nt"=nt, "nc"=nc, "dt"= dt, "tempi" = tempi, "erroE" = erroE)
    return(final)
    
}

D3grafLiteOD<-function(nns, nx, ns, nsa, nt, nc, dt, massa, gamma, ctelastica, alfa, dalfa, xa, dxa  ){
    
    ddalfa <- seq(alfa, alfa+ (nsa-1)*dalfa , dalfa)
    ddxa <- seq(xa, xa+ (ns-1)*dxa , dxa)
    dsde <- matrix( nrow = nsa, ncol = ns)
    erroE <- matrix( nrow = nsa, ncol = ns)
    
    tempi <- matrix( nrow = nsa, ncol = ns)
    si <- matrix( nrow = nsa, ncol = ns)
    ei <- matrix( nrow = nsa, ncol = ns)
    
    for(l in 1:nsa){
        for(j in 1:ns){
            
            qerin <- SimSimplesAlfactepra3DOD(nns, nx, 1, nt, nc, dt, massa, gamma, ctelastica, ddalfa[l] , ddxa[j]  , 0)
            tempi[l,j] <- qerin$temp[1]
            si[l,j] <- qerin$S[1]
            ei[l,j] <- qerin$E[1]
            erroE[l,j]<-mean(qerin$erroE)
            
        }
        print(l)
    }
    
    final<-list("ei"=ei, "si"= si, "alfa"=ddalfa, "xa"=ddxa , "nns" = nns, "ns"=ns, "nsa"=nsa, "nt"=nt, "nc"=nc, "dt"= dt, "tempi" = tempi, "erroE" = erroE)
    return(final)
    
}


D3grafLiteEcte<-function(nns, nx, ns, nsa, nt, nc, dt, massa, gamma, ctelastica, alfa, dalfa, energia, de  ){
    
    
    ddalfa <- seq(alfa, alfa+ (nsa-1)*dalfa , dalfa)
    #ddxa<- (sqrt(2* energia))/sqrt(   1/(ctelastica + massa*alfa**2 + alfa*gamma)   +    (2*massa*alfa)/(ctelastica + massa*alfa**2 + alfa*gamma)   )
    dde <- seq(energia, energia+ (ns-1)*de , de)
    #ddxa<- matrix( nrow = nsa, ncol = ns)
    ddxa<- numeric( nsa*ns)
    
    
    dsde<-matrix( nrow = nsa, ncol = ns)
    erroE<-matrix( nrow = nsa, ncol = ns)
    
    tempi<-matrix( nrow = nsa, ncol = ns)
    si<-matrix( nrow = nsa, ncol = ns)
    ei<-matrix( nrow = nsa, ncol = ns)
    
    g<-0
    
    
    for(l in 1:nsa){
        for(j in 1:ns){
            
            g<- g +1
            
            xa<- (sqrt(2* dde[j]))/sqrt(   1/(ctelastica + massa* ddalfa[l]**2 +  ddalfa[l]*gamma)   +    (2*massa* ddalfa[l])/(ctelastica + massa* ddalfa[l]**2 +  ddalfa[l]*gamma)   )
            ddxa[g]<-xa
            
            qerin <- SimSimplesAlfactepra3D(nns, nx, 1, nt, nc, dt, massa, gamma, ctelastica, ddalfa[l] , xa  , 0)
            tempi[l,j] <- qerin$temp[1]
            si[l,j] <- qerin$S[1]
            ei[l,j] <- qerin$E[1]
            erroE[l,j]<-mean(qerin$erroE)
            
        }
        print(l)
        
    }
    
    final<-list("dde" = dde,"ei"=ei, "si"= si, "alfa"=ddalfa, "xa" = ddxa , "nns" = nns, "ns"=ns, "nsa"=nsa, "nt"=nt, "nc"=nc, "dt"= dt, "tempi" = tempi, "erroE" = erroE)
    return(final)
    
}


D3grafRLTcte<-function(nns, nx, ns, nsa, nt, nc, dt, massa, gamma, ctelastica, alfa, dalfa, xa, dxa , deltaalfa ){
    
    
    
    ddalfa <- seq(alfa, alfa+ (nsa-1)*dalfa , dalfa)
    
    ddxa <- seq(xa, xa+ (ns-1)*dxa , dxa)
    
    dsde<-matrix( nrow = nsa, ncol = ns)
    
    erroE<-matrix( nrow = nsa, ncol = ns)
    
    
    
    tempi<-matrix( nrow = nsa, ncol = ns)
    
    tempr<-matrix( nrow = nsa, ncol = ns)
    
    templ<-matrix( nrow = nsa, ncol = ns)
    
    
    
    si<-matrix( nrow = nsa, ncol = ns)
    
    sr<-matrix( nrow = nsa, ncol = ns)
    
    sl<-matrix( nrow = nsa, ncol = ns)
    
    
    
    ei<-matrix( nrow = nsa, ncol = ns)
    
    er<-matrix( nrow = nsa, ncol = ns)
    
    el<-matrix( nrow = nsa, ncol = ns)
    
    
    for(l in 1:nsa){
        
        for(j in 1:ns){
            
            # function(nns, nx, ns, nt, nc, dt, massa, gamma, ctelastica, alfa, xa, dalfa){
            
            qerin<-SimSimplesTctepra3D(nns, nx, 3, nt, nc, dt, massa, gamma, ctelastica, ddalfa[l] - deltaalfa , ddxa[j] , deltaalfa)
            
            
            templ[l,j]<-qerin$temp[1]
            
            tempi[l,j]<-qerin$temp[2]
            
            tempr[l,j]<-qerin$temp[2]
            
            
            sl[l,j] <- qerin$S[1]
            
            si[l,j] <- qerin$S[2]
            
            sr[l,j] <- qerin$E[3]
            
            
            el[l,j] <- qerin$E[1]
            
            ei[l,j] <- qerin$E[2]
            
            er[l,j] <- qerin$E[3]
            
            
            
            dsde[l,j]<- (qerin$S[3]-qerin$S[1])/(qerin$E[3]-qerin$E[1])
            
            erroE[l,j]<-mean(qerin$erroE)
            
            
        }
        
        
        print(l)
        
        
    }
    
    
    final<-list("ei"=ei, "er"=er  , "el"=el  ,"si"= si, "sii"= sr , "sii"= sl , "alfa"=ddalfa, "xa"=ddxa , "nns" = nns, "ns"=ns, "nsa"=nsa, "nt"=nt, "nc"=nc, "dt"= dt, "tempi" = tempi, "tempr" = tempr, "templ" = templ, "erroE" = erroE, "dsde"=dsde)
    #"momentosE"=momentationE, "cumulantesE"=cumulationE, "momentosV"=momentationV, "cumulantesV"=cumulationV, "momentosX"=momentationX, "cumulantesX"=cumulationX )
    
    return(final)
    
    
}



D3grafMK<-function(nns, nx, ns, nsa, nt, nc, dt, massa, dmassa, gamma, ctelastica, dctelastica ,alfa, xa ){
    
    
    
    
    ddm <- seq(massa, massa+ (nsa-1)*dmassa , dmassa)
    
    ddk <- seq(ctelastica, ctelastica+ (ns-1)*dctelastica , dctelastica)
    

    
    # dsde<-matrix( nrow = nsa, ncol = ns)
    
    erroE<-matrix( nrow = nsa, ncol = ns)
    
    tempi<-matrix( nrow = nsa, ncol = ns)
    
    si<-matrix( nrow = nsa, ncol = ns)
    
    ei<-matrix( nrow = nsa, ncol = ns)
    
    for(l in 1:nsa){
        
        for(j in 1:ns){
            
            #SimSimplesAlfactepra3D<- function(nns, nx, ns, nt, nc, dt, massa, gamma, ctelastica, alfa, xa, dxa){
            
            qerin<-SimSimplesAlfactepra3D(nns, nx, 1, nt, nc, dt, ddm[j], gamma, ddk[l], alfa , xa , 0)
            
            tempi[l,j]<-qerin$temp[1]

            si[l,j] <- qerin$S[1]
            
            ei[l,j] <- qerin$E[1]
            
            erroE[l,j]<-mean(qerin$erroE)
        
        }
        
        
        
        
        print(l)
        
        
    }
    
    
    final<-list("ei"=ei ,"si"= si , "alfa"=alfa, "xa"=xa , "nns" = nns, "ns"=ns, "nsa"=nsa, "nt"=nt, "nc"=nc, "dt"= dt, "tempi" = tempi, "erroE" = erroE, "m"=ddm, "k"=ddk)
    #"momentosE"=momentationE, "cumulantesE"=cumulationE, "momentosV"=momentationV, "cumulantesV"=cumulationV, "momentosX"=momentationX, "cumulantesX"=cumulationX )
    
    return(final)
    
    
}

D3grafRLalfavaria<-function(nns, nx, ns, nsa, nt, nc, dt, massa, gamma, ctelastica, alfa, dalfa, deltalfa, xa, dxa  ){
    
    
    
    
    ddalfa <- seq(alfa, alfa+ (nsa-1)*dalfa , dalfa)
    
    ddxa <- seq(xa, xa+ (ns-1)*dxa , dxa)
    
    dsde<-matrix( nrow = nsa, ncol = ns)
    
    erroE<-matrix( nrow = nsa, ncol = ns)
    
    
    
    tempi<-matrix( nrow = nsa, ncol = ns)
    
    tempr<-matrix( nrow = nsa, ncol = ns)
    
    templ<-matrix( nrow = nsa, ncol = ns)
    
    
    
    si<-matrix( nrow = nsa, ncol = ns)
    
    sr<-matrix( nrow = nsa, ncol = ns)
    
    sl<-matrix( nrow = nsa, ncol = ns)
    
    
    
    ei<-matrix( nrow = nsa, ncol = ns)
    
    er<-matrix( nrow = nsa, ncol = ns)
    
    el<-matrix( nrow = nsa, ncol = ns)
    
    
    
    
    for(l in 1:nsa){
        
        for(j in 1:ns){
            
            qerin<-SimSimplesXactepra3D(nns, nx, 3, nt, nc, dt, massa, gamma, ctelastica, ddalfa[l] - deltalfa , deltalfa , ddxa[j])
            
            #qerin<-SimSimplesAlfactepra3D(nns, nx, 3, nt, nc, dt, massa, gamma, ctelastica, ddalfa[l] , ddxa[j] - dxa , dxa)
            
            templ[l,j]<-qerin$temp[1]
            
            tempi[l,j]<-qerin$temp[2]
            
            tempr[l,j]<-qerin$temp[2]
            
            
            sl[l,j] <- qerin$S[1]
            
            si[l,j] <- qerin$S[2]
            
            sr[l,j] <- qerin$E[3]
            
            
            el[l,j] <- qerin$E[1]
            
            ei[l,j] <- qerin$E[2]
            
            er[l,j] <- qerin$E[3]
            
            
            
            dsde[l,j]<- (qerin$S[3]-qerin$S[1])/(qerin$E[3]-qerin$E[1])
            
            erroE[l,j]<-mean(qerin$erroE)
            
            
        }
        
        
        
        
        print(l)
        
        
    }
    
    
    final<-list("ei"=ei, "er"=er  , "el"=el  ,"si"= si, "sii"= sr , "sii"= sl , "alfa"=ddalfa, "xa"=ddxa , "nns" = nns, "ns"=ns, "nsa"=nsa, "nt"=nt, "nc"=nc, "dt"= dt, "tempi" = tempi, "tempr" = tempr, "templ" = templ, "erroE" = erroE, "dsde"=dsde)
    #"momentosE"=momentationE, "cumulantesE"=cumulationE, "momentosV"=momentationV, "cumulantesV"=cumulationV, "momentosX"=momentationX, "cumulantesX"=cumulationX )
    
    return(final)
    
    
}

D3grafRLlangevinMK<-function(nns, nx, ns, nsa, nt, nc, dt, massa, dm, gamma, ctelastica, dctelastica, temp){
    
    ddmassa <- seq(massa, massa+ (nsa-1)*dm , dm)
    
    ddk <- seq(ctelastica, ctelastica+ (ns-1)*dctelastica , dctelastica)
    

    erroE<-matrix( nrow = nsa, ncol = ns)
    

    tempi<-matrix( nrow = nsa, ncol = ns)
    
    
    si<-matrix( nrow = nsa, ncol = ns)

    
    ei<-matrix( nrow = nsa, ncol = ns)
    
    for(l in 1:nsa){
        
        for(j in 1:ns){
            
            #
            # SimSimplisLangevin<- function(nns, nx, ns, nt, nc, dt, massa, gamma, ctelastica, temp){
            #
            
            qerin<-SimSimplisLangevin(nns, nx, 1, nt, nc, dt, ddmassa[l], gamma, ddk[j] , temp)
            
            tempi[l,j]<-qerin$temp[1]

            si[l,j] <- qerin$S[1]
            
            ei[l,j] <- qerin$E[1]
    
    
            erroE[l,j]<-mean(qerin$erroE)
            
            
        }
        
        print(l)
        
    }
    
    
    final<-list("ei"=ei  ,"si"= si , "nns" = nns, "ns"=ns, "nsa"=nsa, "nt"=nt, "nc"=nc, "dt"= dt, "tempi" = tempi, "erroE" = erroE, "m"=ddmassa, "k"=ddk)
    #"momentosE"=momentationE, "cumulantesE"=cumulationE, "momentosV"=momentationV, "cumulantesV"=cumulationV, "momentosX"=momentationX, "cumulantesX"=cumulationX )
    
    return(final)
    
    
}

D3grafRLlangevinGT<-function(nns, nx, ns, nsa, nt, nc, dt, massa, gamma, dgamma, ctelastica, temp, dtemp){
    
    ddgamma <- seq(gamma, gama+ (nsa-1)*dgamma , dgamma)
    
    ddtemp <- seq(temp, temp+ (ns-1)*dtemp , dtemp)
    
    
    erroE<-matrix( nrow = nsa, ncol = ns)
    
    
    tempi<-matrix( nrow = nsa, ncol = ns)
    
    
    si<-matrix( nrow = nsa, ncol = ns)
    
    
    ei<-matrix( nrow = nsa, ncol = ns)
    
    for(l in 1:nsa){
        
        for(j in 1:ns){
            
            qerin<-SimSimplisLangevin(nns,nx,1, nt, nc, dt, massa, ddgamma[l], ctelastica, ddtemp[j])
            
            tempi[l,j]<-qerin$temp[1]
            
            si[l,j] <- qerin$S[1]
            
            ei[l,j] <- qerin$E[1]
            
            
            erroE[l,j]<-mean(qerin$erroE)
            
            
        }
        
        print(l)
        
    }
    
    
    final<-list("ei"=ei  ,"si"= si , "nns" = nns, "ns"=ns, "nsa"=nsa, "nt"=nt, "nc"=nc, "dt"= dt, "tempi" = tempi, "erroE" = erroE, "gamma"=ddgamma, "temp"=ddtemp)
    #"momentosE"=momentationE, "cumulantesE"=cumulationE, "momentosV"=momentationV, "cumulantesV"=cumulationV, "momentosX"=momentationX, "cumulantesX"=cumulationX )
    
    return(final)
    
    
}


#salvacsv<-function(final){
    
    
##    experimento<-data.frame(final$ei,final$si,final$dsde,final$tempi )
#    colnames(experimento) <- c('E','S',  'dsde',  'temp')
    
    
##    write.csv(experimento, "experimento.csv" )

#}


PrintD3graf<-function(final){
    
    require('plot3D')
    
    
    for(i in 1:length(final$alfa))
    {
        if(i==1)
        {
            x<-rep(final$alfa[i],length(final$xa))
        }
        else
        {
            x<-c(x,rep(final$alfa[i],length(final$xa)))
        }
    }
    


    y<-rep(final$xa,length(final$alfa))
        
    
    ee<- numeric(length(final$xa)*length(final$alfa))
    
    for(i in 1:length(final$alfa))
    {
        for (j in 1:length(final$xa))
        {
            ee[(i-1)*length(final$xa) + j] <- final$ei[i,j]
        }
        
    }
    
    
    jpeg("energia.jpg")
    
    scatter3D(x,y,ee,main='Energia',xlab='alfa',ylab='amplitude')
    
    dev.off()
    
    
    ss<- numeric(length(final$xa)*length(final$alfa))
    
    for(i in 1:length(final$alfa))
    {
        for (j in 1:length(final$xa))
        {
            ss[(i-1)*length(final$xa) + j] <- final$si[i,j]
        }
        
    }
    
    jpeg("entropia.jpg")
    
    scatter3D(x,y,ss, main='Entropia',xlab='alfa',ylab='amplitude')
    
    dev.off()
    
    
    dd<- numeric(length(final$xa)*length(final$alfa))
    
    for(i in 1:length(final$alfa))
    {
        for (j in 1:length(final$xa))
        {
            dd[(i-1)*length(final$xa) + j] <- final$dsde[i,j]
        }
    }
    
    jpeg("dsde.jpg")
    
    scatter3D(x,y,dd, main='DsDe',xlab='alfa',ylab='amplitude')
    
    dev.off()

    tt<- numeric(length(final$xa)*length(final$alfa))

    for(i in 1:length(final$alfa))
    {
        for (j in 1:length(final$xa))
        {
            tt[(i-1)*length(final$xa) + j] <- final$tempi[i,j]
        }
    }

    jpeg("temp.jpg")

    scatter3D(x,y,tt , main='Temperatura',xlab='alfa',ylab='amplitude')

    dev.off()


    experimento<-data.frame(final$ei,final$si,final$dsde,final$tempi )
    colnames(experimento) <- c('E','S',  'dsde',  'temp')


    write.csv(experimento, "experimento.csv" )


}

PrintD3csv<-function(final, letra){
    
    
    for(i in 1:length(final$alfa))
    {
        if(i==1)
        {
            x<-rep(final$alfa[i],length(final$xa))
        }
        else
        {
            x<-c(x,rep(final$alfa[i],length(final$xa)))
        }
    }
    
    
    
    y<-rep(final$xa,length(final$alfa))
    
    
    ee<- numeric(length(final$xa)*length(final$alfa))
    
    for(i in 1:length(final$alfa))
    {
        for (j in 1:length(final$xa))
        {
            ee[(i-1)*length(final$xa) + j] <- final$ei[i,j]
        }
        
    }
    
    
    
    ss<- numeric(length(final$xa)*length(final$alfa))
    
    for(i in 1:length(final$alfa))
    {
        for (j in 1:length(final$xa))
        {
            ss[(i-1)*length(final$xa) + j] <- final$si[i,j]
        }
        
    }
    
    
    dd<- numeric(length(final$xa)*length(final$alfa))
    
    for(i in 1:length(final$alfa))
    {
        for (j in 1:length(final$xa))
        {
            dd[(i-1)*length(final$xa) + j] <- final$dsde[i,j]
        }
    }
    

    
    tt<- numeric(length(final$xa)*length(final$alfa))
    
    for(i in 1:length(final$alfa))
    {
        for (j in 1:length(final$xa))
        {
            tt[(i-1)*length(final$xa) + j] <- final$tempi[i,j]
        }
    }
    

    
    experimento<-data.frame(x,y,ss)
    write.table(experimento, paste(letra ,"entropia.csv", sep='') , sep=",", row.names=FALSE, col.names=FALSE)
    
    experimento<-data.frame(x,y,ee)
    write.table(experimento, paste(letra,"energia.csv", sep='') , sep=",", row.names=FALSE, col.names=FALSE)
    
    experimento<-data.frame(x,y,tt)
    write.table(experimento, paste(letra,"temp.csv", sep='') , sep=",", row.names=FALSE, col.names=FALSE)
    
    experimento<-data.frame(x,y,dd)
    write.table(experimento, paste(letra,"dsde.csv", sep=''), sep="," ,row.names=FALSE, col.names=FALSE)
    
    
    
    
}


PrintD3csvLite<-function(final, letra){
    
    
    for(i in 1:length(final$alfa))
    {
        if(i==1)
        {
            x<-rep(final$alfa[i],length(final$xa))
        }
        else
        {
            x<-c(x,rep(final$alfa[i],length(final$xa)))
        }
    }
    
    
    
    y<-rep(final$xa,length(final$alfa))
    
    ee<- numeric(length(final$xa)*length(final$alfa))
    
    for(i in 1:length(final$alfa))
    {
        for (j in 1:length(final$xa))
        {
            ee[(i-1)*length(final$xa) + j] <- final$ei[i,j]
        }
        
    }
    
    
    
    ss<- numeric(length(final$xa)*length(final$alfa))
    
    for(i in 1:length(final$alfa))
    {
        for (j in 1:length(final$xa))
        {
            ss[(i-1)*length(final$xa) + j] <- final$si[i,j]
        }
        
    }
    
    
 #   dd<- numeric(length(final$xa)*length(final$alfa))
    
  #  for(i in 1:length(final$alfa))
  #  {
  #      for (j in 1:length(final$xa))
  #      {
  #          dd[(i-1)*length(final$xa) + j] <- final$dsde[i,j]
  #      }
  #  }
    
    
    
    tt<- numeric(length(final$xa)*length(final$alfa))
    
    for(i in 1:length(final$alfa))
    {
        for (j in 1:length(final$xa))
        {
            tt[(i-1)*length(final$xa) + j] <- final$tempi[i,j]
        }
    }
    
    
    
    experimento<-data.frame(x,y,ss)
    write.table(experimento, paste(letra ,"entropia.csv", sep='') , sep=",", row.names=FALSE, col.names=FALSE)
    
    experimento<-data.frame(x,y,ee)
    write.table(experimento, paste(letra,"energia.csv", sep='') , sep=",", row.names=FALSE, col.names=FALSE)
    
    experimento<-data.frame(x,y,tt)
    write.table(experimento, paste(letra,"temp.csv", sep='') , sep=",", row.names=FALSE, col.names=FALSE)
    
 #   experimento<-data.frame(x,y,dd)
 #   write.table(experimento, paste(letra,"dsde.csv", sep=''), sep="," ,row.names=FALSE, col.names=FALSE)
    
}

PrintD3csvLiteEcte<-function(final, letra){
    
    
    nsa<- final$nsa
    ns<- final$ns
    
    for(i in 1:ns)
    {
        if(i==1)
        {
            x<-rep(final$alfa[i],nsa)
        }
        else
        {
            x<-c(x,rep(final$alfa[i],nsa))
        }
    }
    
    y<- final$xa
    
    ee<- numeric(length(final$xa))
    
    for(i in 1:length(final$alfa))
    {
        for (j in 1:nsa)
        {
            ee[(i-1)*nsa + j] <- final$ei[i,j]
        }
        
    }
    
    
    
    ss<- numeric(length(final$xa))
    
    for(i in 1:length(final$alfa))
    {
        for (j in 1:nsa)
        {
            ss[(i-1)*nsa + j] <- final$si[i,j]
        }
        
    }
    
    
    #   dd<- numeric(length(final$xa)*length(final$alfa))
    
    #  for(i in 1:length(final$alfa))
    #  {
    #      for (j in 1:length(final$xa))
    #      {
    #          dd[(i-1)*length(final$xa) + j] <- final$dsde[i,j]
    #      }
    #  }
    
    
    
    tt<- numeric(length(final$xa))
    
    for(i in 1:length(final$alfa))
    {
        for (j in 1:nsa)
        {
            tt[(i-1)*nsa + j] <- final$tempi[i,j]
        }
    }
    
    
    
    experimento<-data.frame(x,y,ss)
    write.table(experimento, paste(letra ,"entropia.csv", sep='') , sep=",", row.names=FALSE, col.names=FALSE)
    
    experimento<-data.frame(x,y,ee)
    write.table(experimento, paste(letra,"energia.csv", sep='') , sep=",", row.names=FALSE, col.names=FALSE)
    
    experimento<-data.frame(x,y,tt)
    write.table(experimento, paste(letra,"temp.csv", sep='') , sep=",", row.names=FALSE, col.names=FALSE)
    
    #   experimento<-data.frame(x,y,dd)
    #   write.table(experimento, paste(letra,"dsde.csv", sep=''), sep="," ,row.names=FALSE, col.names=FALSE)
    
}



PrintD3csvLangeMK<-function(final){
    
    # require('plot3D')
    ###
    #
    #
    #
    #
    #
    #
    
    
    for(i in 1:length(final$m))
    {
        if(i==1)
        {
            x<-rep(final$m[i],length(final$k))
        }
        else
        {
            x<-c(x,rep(final$m[i],length(final$k)))
        }
    }
    
    
    
    y<-rep(final$k,length(final$m))
    
    
    ee<- numeric(length(final$k)*length(final$m))
    
    for(i in 1:length(final$m))
    {
        for (j in 1:length(final$k))
        {
            ee[(i-1)*length(final$k) + j] <- final$ei[i,j]
        }
        
    }
    
    
    
    ss<- numeric(length(final$k)*length(final$m))
    
    for(i in 1:length(final$m))
    {
        for (j in 1:length(final$k))
        {
            ss[(i-1)*length(final$k) + j] <- final$si[i,j]
        }
        
    }
    
    
    #   dd<- numeric(length(final$k)*length(final$m))
    
    #for(i in 1:length(final$m))
    #{
    #    for (j in 1:length(final$k))
    #    {
    #        dd[(i-1)*length(final$k) + j] <- final$dsde[i,j]
    #    }
    #}
    
    
    
    tt<- numeric(length(final$k)*length(final$m))
    
    for(i in 1:length(final$m))
    {
        for (j in 1:length(final$k))
        {
            tt[(i-1)*length(final$k) + j] <- final$tempi[i,j]
        }
    }
    
    
    
    experimento<-data.frame(x,y,ss)
    write.table(experimento, "entropia.csv" , sep=",", row.names=FALSE, col.names=FALSE)
    
    experimento<-data.frame(x,y,ee)
    write.table(experimento, "energia.csv" , sep=",", row.names=FALSE, col.names=FALSE)
    
    experimento<-data.frame(x,y,tt)
    write.table(experimento, "temp.csv" , sep=",", row.names=FALSE, col.names=FALSE)
    
    # experimento<-data.frame(x,y,dd)
    #write.table(experimento, "dsde.csv", sep="," ,row.names=FALSE, col.names=FALSE)
    
    
    
}

PrintD3csvLangeMKERRITO<-function(final){
    
    # require('plot3D')
    ###
    #
    #
    #
    #
    #
    #
    
    
    for(i in 1:length(m))
    {
        if(i==1)
        {
            x<-rep(m[i],length(k))
        }
        else
        {
            x<-c(x,rep(m[i],length(k)))
        }
    }
    
    
    
    y<-rep(k,length(m))
    
    
    ee<- numeric(length(k)*length(m))
    
    for(i in 1:length(m))
    {
        for (j in 1:length(k))
        {
            ee[(i-1)*length(k) + j] <- qero$ei[i,j]
        }
        
    }
    
    
    
    ss<- numeric(length(k)*length(m))
    
    for(i in 1:length(m))
    {
        for (j in 1:length(k))
        {
            ss[(i-1)*length(k) + j] <- qero$si[i,j]
        }
        
    }
    
    
    dd<- numeric(length(k)*length(m))
    
    for(i in 1:length(m))
    {
        for (j in 1:length(k))
        {
            dd[(i-1)*length(K) + j] <- qero$dsde[i,j]
        }
    }
    
    
    
    tt<- numeric(length(k)*length(m))
    
    for(i in 1:length(m))
    {
        for (j in 1:length(k))
        {
            tt[(i-1)*length(k) + j] <- qero$tempi[i,j]
        }
    }
    
    
    
    experimento<-data.frame(x,y,ss)
    write.table(experimento, "entropia.csv" , sep=",", row.names=FALSE, col.names=FALSE)
    
    experimento<-data.frame(x,y,ee)
    write.table(experimento, "energia.csv" , sep=",", row.names=FALSE, col.names=FALSE)
    
    experimento<-data.frame(x,y,tt)
    write.table(experimento, "temp.csv" , sep=",", row.names=FALSE, col.names=FALSE)
    
    experimento<-data.frame(x,y,dd)
    write.table(experimento, "dsde.csv", sep="," ,row.names=FALSE, col.names=FALSE)
    
    
    
}



print3grafagora<-function(final){
    
    ddalfa <- seq(0.1,0.1 + 0.7*(14), 0.7)
    
    ddxa <- seq(0.1,0.1 + 0.7*(14), 0.7)
    
    require("plot3D")
    
    
    
    
    
    for(i in 1:length(ddalfa))
    {
        if(i==1)
        {
            x<-rep(ddalfa[i],length(ddxa))
        }
        else
        {
            x<-c(x,rep(ddalfa[i],length(ddxa)))
        }
    }
    
    
    
    y<-rep(ddxa,length(ddalfa))
    
    
    ee<- numeric(length(ddxa)*length(ddalfa))
    
    for(i in 1:length(ddalfa))
    {
        for (j in 1:length(ddxa))
        {
            z[(i-1)*length(ddxa) + j] <- final$ei[i,j]
        }
        
    }
    
    
    jpeg("energia.jpg")
    
    scatter3D(x,y,ee,main='energia',xlab='alfa',ylab='amplitude')
    
    dev.off()
    
    ss<- numeric(length(ddxa)*length(ddalfa))
    
    for(i in 1:length(ddalfa))
    {
        for (j in 1:length(ddxa))
        {
            ss[(i-1)*length(ddxa) + j] <- final$si[i,j]
        }
        
    }
    
    jpeg("entropia.jpg")
    
    scatter3D(x,y,ss,main='entropia',xlab='alfa',ylab='amplitude')
    
    dev.off()
    
    
    dd<- numeric(length(ddxa)*length(ddalfa))
    
    for(i in 1:length(ddalfa))
    {
        for (j in 1:length(ddxa))
        {
            dd[(i-1)*length(ddxa) + j] <- final$dsde[i,j]
        }
        
    }
    
    jpeg("dsde.jpg")
    
    scatter3D(x,y,dd,main='dsde',xlab='alfa',ylab='amplitude')
    
    dev.off()
    
    tt<- numeric(length(ddxa)*length(ddalfa))
    
    for(i in 1:length(ddalfa))
    {
        for (j in 1:length(ddxa))
        {
            tt[(i-1)*length(ddxa) + j] <- final$tempi[i,j]
        }
        
    }
    jpeg("temp.jpg")
    
    scatter3D(x,y,tt, main='Temperatura',xlab='alfa',ylab='amplitude')
    
    dev.off()



    experimento<-data.frame(x,y,ss)
    write.table(experimento, "entropia.csv" ,row.names=FALSE, col.names=FALSE)

    experimento<-data.frame(x,y,ee)
    write.table(experimento, "energia.csv" ,row.names=FALSE, col.names=FALSE)

    experimento<-data.frame(x,y,tt)
    write.table(experimento, "temp.csv" ,row.names=FALSE, col.names=FALSE)

    experimento<-data.frame(x,y,dd)
    write.table(experimento, "dsde.csv",sep="," ,row.names=FALSE, col.names=FALSE)

    

}



SimAmortizadoraAlfaTctecontlogKurt<- function(nns, nx, ns, nt, nc, dt, massa, gamma, ctelastica, alfa, xa, dxa){
    
    # nns<-10
    #nx<-10
    #ns<-10
    #nt<-100
    #nc<-100
    #dt<-0.001
    #massa<-1
    #gamma<-1
    #ctelastica<-1
    #alfa<-0
    #dalfa<-0.01
    #temp<-1
    #
    #
    # zoado... nx nao aumentou a precisao do valor da energia, nc aumentou, preciso testar o nns
    
    anaE <- numeric(ns)
    
    energia <- numeric(ns)
    
    entropia <- numeric(ns)
    
    analiticE <- numeric(ns)
    
    erroE <- numeric(ns)
    
    ddxa <- seq(xa, xa+ ns*dxa , dxa)
    
    #ddalfa <- seq(alfa +dalfa, alfa+ ns * dalfa, dalfa)
    
    
    #E
    
    momentationE<- vector("list", ns)
    
    cumulationE<- vector("list", ns)
    
    #X
    
    momentationX<- vector("list",ns)
    
    cumulationX<- vector("list", ns)
    
    #V
    
    momentationV<- vector("list", ns)
    
    cumulationV<- vector("list", ns)
    
    
    
    temperatura<-numeric(ns)
    
    
    
    for(j in 1:ns)
    {
        
        momentosE<-numeric(7)
        cumulantesE<-numeric(5)
        momentosX<-numeric(7)
        cumulantesX<-numeric(5)
        momentosV<-numeric(7)
        cumulantesV<-numeric(5)
        
        
        ss<- numeric(nx)
        
        ee<- numeric(nx)
        #momentose<-numeric(nx)
        #cumulantese<-numeric(nx)
        
        for(i in 1:nns)
        {
            listasim<-SimulacaoBimodal(nx, nt, nc, dt, massa, gamma, ctelastica, alfa, 0, ddxa[j], 0)
            keru<-analiseentropicaTctecontmomen(listasim)
            ss<-ss+(keru$S)/nns
            ee<-ee+(keru$E)/nns
            
            #print(momentosE)
            #print(cumulantesE)
            #print(keru$momentosE[[1]])
            
            
            
            momE<- numeric(length( keru$momentosE[[1]] ))
            kumE<- numeric(length(keru$cumulantesE[[1]] ))
            momX<- numeric(length(keru$momentosX[[1]] ))
            kumX<- numeric(length(keru$cumulantesX[[1]] ))
            momV<- numeric(length(keru$momentosV[[1]] ))
            kumV<- numeric(length(keru$cumulantesV[[1]] ))
            
            for(l in 1:nx)
            {
                momE<- momE + keru$momentosE[[l]]/nx
                kumE<- kumE + keru$cumulantesE[[l]]/nx
                momX<- momX + keru$momentosX[[l]]/nx
                kumX<- kumX + keru$cumulantesX[[l]]/nx
                momV<- momV + keru$momentosV[[l]]/nx
                kumV<- kumV + keru$cumulantesV[[l]]/nx
            }
            
            
            
            
            momentosE <-momentosE+(momE)/nns
            #print(momentosE)
            cumulantesE<-cumulantesE+(kumE)/nns
            momentosX<-momentosX+(momX)/nns
            cumulantesX<-cumulantesX+(kumX)/nns
            momentosV<-momentosV+(momV)/nns
            cumulantesV<-cumulantesV+(kumV)/nns
        }
        
        temperatura[j]<-listasim$temp
        
        sss<-mean(ss)
        eee<-mean(ee)
        
        anaE<-keru$PolyE[1]
        
        erroE[j]<-((eee-anaE)/anaE)*100
        
        entropia[j]<-sss
        
        energia[j]<-eee
        
        analiticE[j]<-anaE
        
        momentationE[[j]]<-momentosE
        cumulationE[[j]]<-cumulantesE
        momentationX[[j]]<-momentosX
        cumulationX[[j]]<-cumulantesX
        momentationV[[j]]<-momentosV
        cumulationV[[j]]<-cumulantesV
        
        print(j)
        
    }
    
    Slm<-loess(entropia~ddxa)
    
    bb<-predict(Slm,ddxa, se = TRUE)
    
    B<-bb$fit
    
    
    #FIM DO SIMULACAO BIMODAL T CTE
    #   final<-list("posicao"=xx,"velocidade"=vv, "energia"=ee, "massa"=massa,"gamma"=gamma,"k"=ctelastica, "alfa" = alfa,"dalfa"= dalfa, "xa"=xa, "ns"=ns,"nt"=nt,"nc"=nc,"dt"= dt, "T"=Temp)
    #return(final)
    
    final<-list("PolyE"=analiticE, "E"=energia, "S"=entropia, "PolyS"= B , "alfa"=alfa, "xa"=ddxa , "nns" = nns, "ns"=ns, "nt"=nt, "nc"=nc, "dt"= dt, "temp" = temperatura, "erroE" = erroE,"momentosE"=momentationE, "cumulantesE"=cumulationE, "momentosV"=momentationV, "cumulantesV"=cumulationV, "momentosX"=momentationX, "cumulantesX"=cumulationX )
    
    return(final)
    
}




SimAmortizadoraTctecontlogKurtQ<- function(nns, nx, ns, nt, nc, dt, massa, gamma, ctelastica, alfa, dalfa, temp, Q){
    
    # nns<-10
    #nx<-10
    #ns<-10
    #nt<-100
    #nc<-100
    #dt<-0.001
    #massa<-1
    #gamma<-1
    #ctelastica<-1
    #alfa<-0
    #dalfa<-0.01
    #temp<-1
    #
    #
    # zoado... nx nao aumentou a precisao do valor da energia, nc aumentou, preciso testar o nns
    
    anaE <- numeric(ns)
    
    energia <- numeric(ns)
    
    entropia <- numeric(ns)
    
    analiticE <- numeric(ns)
    
    erroE <- numeric(ns)
    
    ddalfa <- lseq(dalfa, alfa, ns)
    
    #ddalfa <- seq(alfa +dalfa, alfa+ ns * dalfa, dalfa)
    
    
    #E
    
    momentationE<- vector("list", ns)
    
    cumulationE<- vector("list", ns)
    
    #X
    
    momentationX<- vector("list",ns)
    
    cumulationX<- vector("list", ns)
    
    #V
    
    momentationV<- vector("list", ns)
    
    cumulationV<- vector("list", ns)
    
    
    
    
    
    
    
    for(j in 1:ns)
    {
        
        momentosE<-numeric(7)
        cumulantesE<-numeric(5)
        momentosX<-numeric(7)
        cumulantesX<-numeric(5)
        momentosV<-numeric(7)
        cumulantesV<-numeric(5)
        
        
        ss<- numeric(nx)
        
        ee<- numeric(nx)
        #momentose<-numeric(nx)
        #cumulantese<-numeric(nx)
        
        for(i in 1:nns)
        {
            listasim<-SimulacaoBimodalTcte(nx, nt, nc, dt, massa, gamma, ctelastica, ddalfa[j], 0, temp)
            keru<-analiseentropicaQTctecontmomen(listasim, Q)
            ss<-ss+(keru$S)/nns
            ee<-ee+(keru$E)/nns
            
            #print(momentosE)
            #print(cumulantesE)
            #print(keru$momentosE[[1]])
            
            momE<- numeric(length( keru$momentosE[[1]] ))
            kumE<- numeric(length(keru$cumulantesE[[1]] ))
            momX<- numeric(length(keru$momentosX[[1]] ))
            kumX<- numeric(length(keru$cumulantesX[[1]] ))
            momV<- numeric(length(keru$momentosV[[1]] ))
            kumV<- numeric(length(keru$cummulantesV[[1]] ))
            
            for(l in 1:nx)
            {
                momE<- momE + keru$momentosE[[l]]/nx
                kumE<- kumE + keru$cumulantesE[[l]]/nx
                momX<- momX + keru$momentosX[[l]]/nx
                kumX<- kumX + keru$cumulantesX[[l]]/nx
                momV<- momV + keru$momentosV[[l]]/nx
                kumV<- kumV + keru$cummulantesV[[l]]/nx
            }
            
            
            momentosE <-momentosE+(momE)/nns
            #print(momentosE)
            cumulantesE<-cumulantesE+(kumE)/nns
            momentosX<-momentosX+(momX)/nns
            cumulantesX<-cumulantesX+(kumX)/nns
            momentosV<-momentosV+(momV)/nns
            cumulantesV<-cumulantesV+(kumV)/nns
        }
        
        sss<-mean(ss)
        eee<-mean(ee)
        
        anaE<-keru$PolyE[1]
        
        erroE[j]<-((eee-anaE)/anaE)*100
        
        entropia[j]<-sss
        
        energia[j]<-eee
        
        analiticE[j]<-anaE
        
        momentationE[[j]]<-momentosE
        cumulationE[[j]]<-cumulantesE
        momentationX[[j]]<-momentosX
        cumulationX[[j]]<-cumulantesX
        momentationV[[j]]<-momentosV
        cumulationV[[j]]<-cumulantesV
        
        
    }
    
    Slm<-loess(entropia~ddalfa)
    
    bb<-predict(Slm,ddalfa, se = TRUE)
    
    B<-bb$fit
    
    
    #FIM DO SIMULACAO BIMODAL T CTE
    #   final<-list("posicao"=xx,"velocidade"=vv, "energia"=ee, "massa"=massa,"gamma"=gamma,"k"=ctelastica, "alfa" = alfa,"dalfa"= dalfa, "xa"=xa, "ns"=ns,"nt"=nt,"nc"=nc,"dt"= dt, "T"=Temp)
    #return(final)
    
    final<-list("PolyE"=analiticE, "E"=energia, "S"=entropia, "PolyS"= B , "alfa"=ddalfa, "nns" = nns, "ns"=ns, "nt"=nt, "nc"=nc, "dt"= dt, "temp" = temp, "erroE" = erroE,"momentosE"=momentationE, "cumulantesE"=cumulationE, "momentosV"=momentationV, "cumulantesV"=cumulationV, "momentosX"=momentationX, "cumulantesX"=cumulationX )
    
    return(final)
    
}

SimAmortizadoraXActecontlogKurt<- function(nns, nx, ns, nt, nc, dt, massa, gamma, ctelastica, alfa, dalfa, temp){
    
    # nns<-10
    #nx<-10
    #ns<-10
    #nt<-100
    #nc<-100
    #dt<-0.001
    #massa<-1
    #gamma<-1
    #ctelastica<-1
    #alfa<-0
    #dalfa<-0.01
    #temp<-1
    #
    #
    # zoado... nx nao aumentou a precisao do valor da energia, nc aumentou, preciso testar o nns
    
    anaE <- numeric(ns)
    
    energia <- numeric(ns)
    
    entropia <- numeric(ns)
    
    analiticE <- numeric(ns)
    
    erroE <- numeric(ns)
    
    ddalfa <- lseq(dalfa, alfa, ns)
    
    #ddalfa <- seq(alfa +dalfa, alfa+ ns * dalfa, dalfa)
    
    
    #E
    
    momentationE<- vector("list", ns)
    
    cumulationE<- vector("list", ns)
    
    #X
    
    momentationX<- vector("list",ns)
    
    cumulationX<- vector("list", ns)
    
    #V
    
    momentationV<- vector("list", ns)
    
    cumulationV<- vector("list", ns)
    
    
    
    
    
    
    
    for(j in 1:ns)
    {
        
        momentosE<-numeric(7)
        cumulantesE<-numeric(5)
        momentosX<-numeric(7)
        cumulantesX<-numeric(5)
        momentosV<-numeric(7)
        cumulantesV<-numeric(5)
        
        
        ss<- numeric(nx)
        
        ee<- numeric(nx)
        #momentose<-numeric(nx)
        #cumulantese<-numeric(nx)
        
        for(i in 1:nns)
        {
            listasim<-SimulacaoBimodalXActe(nx, nt, nc, dt, massa, gamma, ctelastica, ddalfa[j], 0, temp)
            keru<-analiseentropicaTctecontmomen(listasim)
            ss<-ss+(keru$S)/nns
            ee<-ee+(keru$E)/nns
            
            #print(momentosE)
            #print(cumulantesE)
            #print(keru$momentosE[[1]])
            
            momE<- numeric(length( keru$momentosE[[1]] ))
            kumE<- numeric(length(keru$cumulantesE[[1]] ))
            momX<- numeric(length(keru$momentosX[[1]] ))
            kumX<- numeric(length(keru$cumulantesX[[1]] ))
            momV<- numeric(length(keru$momentosV[[1]] ))
            kumV<- numeric(length(keru$cumulantesV[[1]] ))
            
            for(l in 1:nx)
            {
                momE<- momE + keru$momentosE[[l]]/nx
                kumE<- kumE + keru$cumulantesE[[l]]/nx
                momX<- momX + keru$momentosX[[l]]/nx
                kumX<- kumX + keru$cumulantesX[[l]]/nx
                momV<- momV + keru$momentosV[[l]]/nx
                kumV<- kumV + keru$cumulantesV[[l]]/nx
            }
            
            
            momentosE <-momentosE+(momE)/nns
            #print(momentosE)
            cumulantesE<-cumulantesE+(kumE)/nns
            momentosX<-momentosX+(momX)/nns
            cumulantesX<-cumulantesX+(kumX)/nns
            momentosV<-momentosV+(momV)/nns
            cumulantesV<-cumulantesV+(kumV)/nns
        }
        
        sss<-mean(ss)
        eee<-mean(ee)
        
        anaE<-keru$PolyE[1]
        
        erroE[j]<-((eee-anaE)/anaE)*100
        
        entropia[j]<-sss
        
        energia[j]<-eee
        
        analiticE[j]<-anaE
        
        momentationE[[j]]<-momentosE
        cumulationE[[j]]<-cumulantesE
        momentationX[[j]]<-momentosX
        cumulationX[[j]]<-cumulantesX
        momentationV[[j]]<-momentosV
        cumulationV[[j]]<-cumulantesV
        
        
    }
    
    Slm<-loess(entropia~ddalfa)
    
    bb<-predict(Slm,ddalfa, se = TRUE)
    
    B<-bb$fit
    
    
    #FIM DO SIMULACAO BIMODAL T CTE
    #   final<-list("posicao"=xx,"velocidade"=vv, "energia"=ee, "massa"=massa,"gamma"=gamma,"k"=ctelastica, "alfa" = alfa,"dalfa"= dalfa, "xa"=xa, "ns"=ns,"nt"=nt,"nc"=nc,"dt"= dt, "T"=Temp)
    #return(final)
    
    final<-list("PolyE"=analiticE, "E"=energia, "S"=entropia, "PolyS"= B , "alfa"=ddalfa, "nns" = nns, "ns"=ns, "nt"=nt, "nc"=nc, "dt"= dt, "temp" = temp, "erroE" = erroE,"momentosE"=momentationE, "cumulantesE"=cumulationE, "momentosV"=momentationV, "cumulantesV"=cumulationV, "momentosX"=momentationX, "cumulantesX"=cumulationX )
    
    return(final)
    
}





SimAmortizadoraTctecontlog<- function(nns, nx, ns, nt, nc, dt, massa, gamma, ctelastica, alfa, dalfa, temp){
    
    # nns<-10
    #nx<-10
    #ns<-10
    #nt<-100
    #nc<-100
    #dt<-0.001
    #massa<-1
    #gamma<-1
    #ctelastica<-1
    #alfa<-0
    #dalfa<-0.01
    #temp<-1
    #
    #
    # zoado... nx nao aumentou a precisao do valor da energia, nc aumentou, preciso testar o nns
    
    anaE <- numeric(ns)
    
    energia <- numeric(ns)
    
    entropia <- numeric(ns)
    
    analiticE <- numeric(ns)
    
    erroE <- numeric(ns)
    
    ddalfa <- lseq(dalfa, alfa, ns)
    
    #ddalfa <- seq(alfa +dalfa, alfa+ ns * dalfa, dalfa)
    
    
    for(j in 1:ns)
    {
        
        
        ss<- numeric(nx)
        
        ee<- numeric(nx)
        #momentose<-numeric(nx)
        #cumulantese<-numeric(nx)
        
        for(i in 1:nns)
        {
            listasim<-SimulacaoBimodalTcte(nx, nt, nc, dt, massa, gamma, ctelastica, ddalfa[j], 0, temp)
            keru<-analiseentropicaTctecont(listasim)
            ss<-ss+(keru$S)/nns
            ee<-ee+(keru$E)/nns
            #momentose<-(keru$momentationE)/nns
            #cumulatesE<-(keru$cumulationE)/nns
        }
        
        sss<-mean(ss)
        eee<-mean(ee)
        
        anaE<-keru$PolyE[1]
        
        erroE[j]<-((eee-anaE)/anaE)*100
        
        entropia[j]<-sss
        
        energia[j]<-eee
        
        analiticE[j]<-anaE
        
        
    }
    
    Slm<-loess(entropia~ddalfa)
    
    bb<-predict(Slm,ddalfa, se = TRUE)
    
    B<-bb$fit
    
    
    #FIM DO SIMULACAO BIMODAL T CTE
    #   final<-list("posicao"=xx,"velocidade"=vv, "energia"=ee, "massa"=massa,"gamma"=gamma,"k"=ctelastica, "alfa" = alfa,"dalfa"= dalfa, "xa"=xa, "ns"=ns,"nt"=nt,"nc"=nc,"dt"= dt, "T"=Temp)
    #return(final)
    
    final<-list("PolyE"=analiticE, "E"=energia, "S"=entropia, "PolyS"= B , "alfa"=ddalfa, "nns" = nns, "ns"=ns, "nt"=nt, "nc"=nc, "dt"= dt, "temp" = temp, "erroE" = erroE)
    
    return(final)
    
}




SimAmortizadoraTctecontlogDT<- function(nns, nx, ns, nt, nc, dt, massa, gamma, ctelastica, alfa, dalfa, temp){
    
    # nns<-10
    #nx<-10
    #ns<-10
    #nt<-100
    #nc<-100
    #dt<-0.001
    #massa<-1
    #gamma<-1
    #ctelastica<-1
    #alfa<-0
    #dalfa<-0.01
    #temp<-1
    #
    #
    # zoado... nx nao aumentou a precisao do valor da energia, nc aumentou, preciso testar o nns
    
    anaE <- numeric(ns)
    
    energia <- numeric(ns)
    
    entropia <- numeric(ns)
    
    analiticE <- numeric(ns)
    
    erroE <- numeric(ns)
    
    ddalfa <- lseq(dalfa, alfa, ns)
    
    #ddalfa <- seq(alfa +dalfa, alfa+ ns * dalfa, dalfa)
    
    
    for(j in 1:ns)
    {
        
        
        ss<- numeric(nx)
        
        ee<- numeric(nx)
        #momentose<-numeric(nx)
        #cumulantese<-numeric(nx)
        
        for(i in 1:nns)
        {
            
            if(ddalfa[j]<1)
            {
            listasim<-SimulacaoBimodalTcte(nx, nt, nc, dt, massa, gamma, ctelastica, ddalfa[j], 0, temp)
            keru<-analiseentropicaTctecont(listasim)
            ss<-ss+(keru$S)/nns
            ee<-ee+(keru$E)/nns

            } else {
            
            listasim<-SimulacaoBimodalTcte(nx, nt, nc*10, dt*0.1 , massa, gamma, ctelastica, ddalfa[j], 0, temp)
            keru<-analiseentropicaTctecont(listasim)
            ss<-ss+(keru$S)/nns
            ee<-ee+(keru$E)/nns
            
            }
        }
        
        sss<-mean(ss)
        eee<-mean(ee)
        
        anaE<-keru$PolyE[1]
        
        erroE[j]<-((eee-anaE)/anaE)*100
        
        entropia[j]<-sss
        
        energia[j]<-eee
        
        analiticE[j]<-anaE
        
        
    }
    
    Slm<-loess(entropia~ddalfa)
    
    bb<-predict(Slm,ddalfa, se = TRUE)
    
    B<-bb$fit
    
    
    #FIM DO SIMULACAO BIMODAL T CTE
    #   final<-list("posicao"=xx,"velocidade"=vv, "energia"=ee, "massa"=massa,"gamma"=gamma,"k"=ctelastica, "alfa" = alfa,"dalfa"= dalfa, "xa"=xa, "ns"=ns,"nt"=nt,"nc"=nc,"dt"= dt, "T"=Temp)
    #return(final)
    
    final<-list("PolyE"=analiticE, "E"=energia, "S"=entropia, "PolyS"= B , "alfa"=ddalfa, "nns" = nns, "ns"=ns, "nt"=nt, "nc"=nc, "dt"= dt, "temp" = temp, "erroE" = erroE)
    
    return(final)
    
}





SimAmortizadoraTctecontlogXV<- function(nns, nx, ns, nt, nc, dt, massa, gamma, ctelastica, alfa, dalfa, temp){
    
    # nns<-10
    #nx<-10
    #ns<-10
    #nt<-100
    #nc<-100
    #dt<-0.001
    #massa<-1
    #gamma<-1
    #ctelastica<-1
    #alfa<-0
    #dalfa<-0.01
    #temp<-1
    #
    #
    # zoado... nx nao aumentou a precisao do valor da energia, nc aumentou, preciso testar o nns
    
    anaE <- numeric(ns)
    
    energia <- numeric(ns)
    
    energiapot<- numeric(ns)
    
    energiakin<- numeric(ns)
    
    entropia <- numeric(ns)
    
    entropiax <- numeric(ns)
    
    entropiav <- numeric(ns)
    
    analiticE <- numeric(ns)
    
    erroE <- numeric(ns)
    
    ddalfa <- lseq(dalfa, alfa, ns)
    
    #ddalfa <- seq(alfa +dalfa, alfa+ ns * dalfa, dalfa)
    
    
    for(j in 1:ns)
    {
        
        
        ss<- numeric(nx)
        ssx<-numeric(nx)
        ssv<-numeric(nx)
        
        
        ee<- numeric(nx)
        #momentose<-numeric(nx)
        #cumulantese<-numeric(nx)
        
        for(i in 1:nns)
        {
            listasim<-SimulacaoBimodalTcte(nx, nt, nc, dt, massa, gamma, ctelastica, ddalfa[j], 0, temp)
            keru<-analiseentropicaTctecontXV(listasim)
            ss<-ss+(keru$S)/nns
            ee<-ee+(keru$E)/nns
            
            ssx<-ssx+(keru$Sx)/nns
            ssv<-ssv+(keru$Sv)/nns

            #momentose<-(keru$momentationE)/nns
            #cumulatesE<-(keru$cumulationE)/nns
            
            #momentosV ---- medir o ponto crítico da curtos
        }
        

#sss<-mean(ss)
        eee<-mean(ee)


        anaE<-keru$PolyE[1]
        
        erroE[j]<-((eee-anaE)/anaE)*100
        
        entropia[j]<-mean(ss)
        
        energia[j]<-eee
 
        entropiax[j]<-mean(ssx)
        entropiav[j]<-mean(ssv)
        
        analiticE[j]<-anaE
        
        
    }
    
    Slm<-loess(entropia~ddalfa)
    
    bb<-predict(Slm,ddalfa, se = TRUE)
    
    B<-bb$fit
    
    
    #FIM DO SIMULACAO BIMODAL T CTE
    #   final<-list("posicao"=xx,"velocidade"=vv, "energia"=ee, "massa"=massa,"gamma"=gamma,"k"=ctelastica, "alfa" = alfa,"dalfa"= dalfa, "xa"=xa, "ns"=ns,"nt"=nt,"nc"=nc,"dt"= dt, "T"=Temp)
    #return(final)
    
    
    #function(nns, nx, ns, nt, nc, dt, massa, gamma, ctelastica, alfa, dalfa, temp){
    

    
    
    
    final<-list("PolyE"=analiticE, "E"=energia, "Sx"=entropiax,"Sv"=entropiav  , "S"=entropia, "PolyS"= B , "alfa"=ddalfa, "nns" = nns, "ns"=ns, "nt"=nt, "nc"=nc, "dt"= dt, "temp" = temp, "erroE" = erroE,  "massa"=massa,"gamma"=gamma,"k"=ctelastica, "xa"=xa)
    
    return(final)
    
}




SimAmortizadora<- function(nns, nx, ns, nt, nc, dt, massa, gamma, ctelastica, alfa, dalfa, xa, dxa){
    
    # zoado... nx nao aumentou a precisao do valor da energia, nc aumentou, preciso testar o nns
    
    anaE <- numeric(ns)
    
    energia <- numeric(ns)
    
    entropia <- numeric(ns)
    
    analiticE <- numeric(ns)
    
    erroE <- numeric(ns)
    
    
    
    
    ddalfa <- seq(alfa +dalfa, alfa+ ns * dalfa, dalfa)
    
    ddxa <- seq(xa +dxa, xa+ ns * dxa, dxa)
    
    
    
    for(j in 1:ns)
    {
        print(j)
        
        ss<- numeric(nx)
        
        ee<- numeric(nx)
        
        
        for(i in 1:nns)
        {
            listasim<-SimulacaoBimodal(nx, nt, nc, dt, massa, gamma, ctelastica, ddalfa[j], 0, xa , 0)
            #function(ns, nt, nc, dt, massa, gamma, ctelastica, alfa, dalfa, xa, dxa){
            keru<-analiseentropicaExai(listasim)
            ss<-ss+(keru$S)/nns
            ee<-ee+(keru$E)/nns
        }
        
        sss<-mean(ss)
        eee<-mean(ee)
        
        anaE<-keru$PolyE[1]
        
        erroE[j]<-((eee-anaE)/anaE)*100
        
        entropia[j]<-sss
        
        energia[j]<-eee
        
        analiticE[j]<-anaE
        
        
    }
    
    Slm<-loess(entropia~ddalfa)
    
    bb<-predict(Slm,ddalfa, se = TRUE)
    
    B<-bb$fit
    
    
    temp<- ((ddxa^2)*(massa*ddalfa))/(gamma*(ctelastica + ddalfa*(massa*ddalfa + gamma)))
    
    #FIM DO SIMULACAO BIMODAL T CTE
    #   final<-list("posicao"=xx,"velocidade"=vv, "energia"=ee, "massa"=massa,"gamma"=gamma,"k"=ctelastica, "alfa" = alfa,"dalfa"= dalfa, "xa"=xa, "ns"=ns,"nt"=nt,"nc"=nc,"dt"= dt, "T"=Temp)
    #return(final)
    
    final<-list("PolyE"=analiticE, "E"=energia, "S"=entropia, "PolyS"= B , "alfa"=ddalfa, "nns" = nns, "ns"=ns, "nt"=nt, "nc"=nc, "dt"= dt, "temp" = temp, "erroE" = erroE)
    
    return(final)
    
}






salvavariosgraficoslog<-function(fim){
    
    
    grafe<-"ExAlfa.jpg"
    
    jpeg(grafe)
    
    plot(fim$E ~ fim$alfa, xlab='alfa', ylab='Energia Media')
    lines(fim$PolyE ~ fim$alfa, col='blue')
    dev.off()
    
    grafs<-"SxAlfa.jpg"

    jpeg(grafs)
    
    plot(fim$S  ~ fim$alfa, xlab='alfa', ylab='Entropia')
    lines(fim$PolyS ~ fim$alfa, col='blue')
    dev.off()
    
    grafse<-"SxE.jpg"
    
    
    
    jpeg(grafse)
    
    plot(fim$S  ~ fim$E, xlab='Energia Media', ylab='Entropia')
    #loo<-loess(fim$S~fim$E)
    lom<-lm(fim$S ~poly(fim$E,4,raw=TRUE))
    Blom<-lom$coefficients[1]+fim$E*lom$coefficients[2]+(fim$E^2)*lom$coefficients[3]+(fim$E^3)*lom$coefficients[4]+(fim$E^4)*lom$coefficients[5]
    

    loo<-loess(fim$S  ~ fim$E)
    
    bb<-predict(loo,fim$E, se = TRUE)
    
    B<-bb$fit
    
    
    lines(Blom ~ fim$E, col='blue')
    dev.off()
    
    lines(B ~ fim$E, col='red')
    dev.off()
    
    
    ######  CALCULO DAS DERIVADAS
    
    derivlo<-numeric(length(fim$S)-1)
    
    #  derivhard<-numeric(length(fim$S)-1)

    
    for(i in 1: (length(fim$S)-1))
    {
        derivlo[i]<-(B[i+1]-B[i])/(fim$E[i+1]-fim$E[i])
        #  derivhard[i]<-(fim$S[i+1]-fim$S[i])/(fim$E[i+1]-fim$E[i])
    }
    
    
    ######
    mE<-tail(fim$E,length(fim$E)-1)
    
    grafsee<-"dSdExE.jpg"
    
    jpeg(grafsee)
    
    
    #   plot(derivhard ~ mE, xlab='Energia Media', ylab='dS/dE')
    
    # loh<-loess(derivhard~mE)
    #C<-predict(loh,mE, se = TRUE)
    #cc<-C$fit
    
    low<-loess(derivlo~mE)
    D<-predict(low,mE, se = TRUE)
    dd<-D$fit
    
    #lines(cc ~ mE, col='blue')
    
    plot(derivlo ~ mE, col='red',xlab='Energia Media', ylab='dS/dE')
    lines(dd ~ mE, col='brown')

    
    dev.off()
    
    
    ######
    malfa<-tail(fim$alfa,length(fim$alfa)-1)
    
    
    
    grafses<-"dSdExAlfa.jpg"
    
    jpeg(grafses)
    
    
    
    #plot(derivhard ~ malfa, xlab='alfa', ylab='dS/dE')
    
    # loh<-loess(derivhard~malfa)
    #EE<-predict(loh,malfa, se = TRUE)
    #ee<-EE$fit
    
    low<-loess(derivlo~malfa)
    FF<-predict(low,malfa, se = TRUE)
    ff<-FF$fit
    
    #  lines(ee ~ malfa, col='blue')
    
    
    plot(derivlo ~ malfa, col='red',xlab='alfa', ylab='dS/dE')
    lines(ff ~ malfa, col='brown')
    
    dev.off()
    
    #####
    experimento<-c(fim$E,fim$S,fim$alfa)
    write.csv(experimento, "experimento.csv" )
    
}


salvavariosgraficoslogXV<-function(fim){
    
    
    grafe<-"ExAlfa.jpg"
    
    jpeg(grafe)
    
    plot(fim$E ~ fim$alfa, xlab='alfa', ylab='Energia Media')
    lines(fim$PolyE ~ fim$alfa, col='blue')
    dev.off()
    
    grafs<-"SxAlfa.jpg"
    
    jpeg(grafs)
    
    plot(fim$S  ~ fim$alfa, xlab='alfa', ylab='Entropia')
    lines(fim$PolyS ~ fim$alfa, col='blue')
    dev.off()
    
    grafsv<-"SVxAlfa.jpg"
    
    jpeg(grafsv)
    
    plot(fim$Sx  ~ fim$alfa, xlab='alfa', ylab='Entropia')
    #lines(fim$PolyS ~ fim$alfa, col='blue')
    dev.off()
    
    grafsx<-"SXxAlfa.jpg"
    
    jpeg(grafsx)
    
    plot(fim$Sv  ~ fim$alfa, xlab='alfa', ylab='Entropia')
    #lines(fim$PolyS ~ fim$alfa, col='blue')
    dev.off()
    
    grafse<-"SxE.jpg"
    
    jpeg(grafse)
    
    plot(fim$S  ~ fim$E, xlab='Energia Media', ylab='Entropia')
    #loo<-loess(fim$S~fim$E)
    lom<-lm(fim$S ~poly(fim$E,4,raw=TRUE))
    Blom<-lom$coefficients[1]+fim$E*lom$coefficients[2]+(fim$E^2)*lom$coefficients[3]+(fim$E^3)*lom$coefficients[4]+(fim$E^4)*lom$coefficients[5]
    
    
    loo<-loess(fim$S  ~ fim$E)
    
    bb<-predict(loo,fim$E, se = TRUE)
    
    B<-bb$fit
    
    
    lines(Blom ~ fim$E, col='blue')
    dev.off()
    
    lines(B ~ fim$E, col='red')
    dev.off()
    
    grafse<-"SXxEpot.jpg"
    
    jpeg(grafse)
    
    ##m=1
    #gamma=1
    
    epot<-(fim$temp*fim$gamma)/(2 * fim$massa *fim$alfa )
    
    
    
    plot(fim$Sx  ~ ekin, xlab='Energia Media', ylab='Entropia')
    dev.off()
    
    grafse<-"SVxEkin.jpg"
    
    jpeg(grafse)
    
    ekin<-rep((1/2)*fim$temp,length(fim$S))
    
    
    plot(fim$Sv  ~ ekin, xlab='Energia Media', ylab='Entropia')
    dev.off()

    
    
    ######  CALCULO DAS DERIVADAS
    
    derivlo<-numeric(length(fim$S)-1)
    
    #  derivhard<-numeric(length(fim$S)-1)
    
    
    for(i in 1: (length(fim$S)-1))
    {
        derivlo[i]<-(B[i+1]-B[i])/(fim$E[i+1]-fim$E[i])
        #  derivhard[i]<-(fim$S[i+1]-fim$S[i])/(fim$E[i+1]-fim$E[i])
    }
    
    
    ######
    mE<-tail(fim$E,length(fim$E)-1)
    
    grafsee<-"dSdExE.jpg"
    
    jpeg(grafsee)
    
    
    #   plot(derivhard ~ mE, xlab='Energia Media', ylab='dS/dE')
    
    # loh<-loess(derivhard~mE)
    #C<-predict(loh,mE, se = TRUE)
    #cc<-C$fit
    
    low<-loess(derivlo~mE)
    D<-predict(low,mE, se = TRUE)
    dd<-D$fit
    
    #lines(cc ~ mE, col='blue')
    
    plot(derivlo ~ mE, col='red',xlab='Energia Media', ylab='dS/dE')
    lines(dd ~ mE, col='brown')
    
    
    dev.off()
    
    
    ######
    malfa<-tail(fim$alfa,length(fim$alfa)-1)
    
    
    
    grafses<-"dSdExAlfa.jpg"
    
    jpeg(grafses)
    
    
    
    #plot(derivhard ~ malfa, xlab='alfa', ylab='dS/dE')
    
    # loh<-loess(derivhard~malfa)
    #EE<-predict(loh,malfa, se = TRUE)
    #ee<-EE$fit
    
    low<-loess(derivlo~malfa)
    FF<-predict(low,malfa, se = TRUE)
    ff<-FF$fit
    
    #  lines(ee ~ malfa, col='blue')
    
    
    plot(derivlo ~ malfa, col='red',xlab='alfa', ylab='dS/dE')
    lines(ff ~ malfa, col='brown')
    
    dev.off()
    
    #####
    experimento<-c(fim$E,fim$S,fim$alfa)
    write.csv(experimento, "experimento.csv" )
    
}








salvasimplesXVKurtlog<-function(fim){
    
    
    grafe<-"ExAlfa.jpg"
    
    jpeg(grafe)
    
    plot(fim$E ~ fim$alfa, xlab='alfa', ylab='Energia Media')
    lines(fim$PolyE ~ fim$alfa, col='blue')
    dev.off()
    
    grafs<-"SxAlfa.jpg"
    
    jpeg(grafs)
    
    plot(fim$S  ~ fim$alfa, xlab='alfa', ylab='Entropia')
    lines(fim$PolyS ~ fim$alfa, col='blue')
    dev.off()
    
    
    grafse<-"SxE.jpg"
    
    jpeg(grafse)
    
    plot(fim$S  ~ fim$E, xlab='Energia Media', ylab='Entropia')
    dev.off()
    
    
    
    
    x1cum<-pegaordemK(fim$cumulantesX,1)
    x2cum<-pegaordemK(fim$cumulantesX,2)
    x3cum<-pegaordemK(fim$cumulantesX,3)
    x4cum<-pegaordemK(fim$cumulantesX,4)
    
    v1cum<-pegaordemK(fim$cumulantesV,1)
    v2cum<-pegaordemK(fim$cumulantesV,2)
    v3cum<-pegaordemK(fim$cumulantesV,3)
    v4cum<-pegaordemK(fim$cumulantesV,4)
    
    grafoe<-"K1xXAlfa.jpg"
    
    
    jpeg(grafoe)
    
    plot(x1cum  ~ fim$alfa, ylab='K1 Posicao', xlab = 'Alfa')
    dev.off()
    
    grafoe<-"K2xXAlfa.jpg"
    
    
    jpeg(grafoe)
    
    plot(x2cum  ~ fim$alfa, ylab='K2 Posicao', xlab = 'Alfa')
    dev.off()
    
    grafoe<-"K3xXAlfa.jpg"
    
    
    jpeg(grafoe)
    
    plot(x3cum  ~ fim$alfa, ylab='K3 Posicao', xlab = 'Alfa')
    dev.off()
    
    grafoe<-"K4xXAlfa.jpg"
    
    
    jpeg(grafoe)
    
    plot(x4cum  ~ fim$alfa, ylab='K4 Posicao', xlab = 'Alfa')
    dev.off()
    
    
    grafsu<-"K1vXAlfa.jpg"
    
    
    jpeg(grafsu)
    
    plot(v1cum  ~ fim$alfa, ylab='K1 Velocidade', xlab = 'Alfa')
    dev.off()
    
    grafsu<-"K2vXAlfa.jpg"
    
    
    jpeg(grafsu)
    
    plot(v2cum  ~ fim$alfa, ylab='K2 Velocidade', xlab = 'Alfa')
    dev.off()
    
    grafsu<-"K3vXAlfa.jpg"
    
    
    jpeg(grafsu)
    
    plot(v3cum  ~ fim$alfa, ylab='K3 Velocidade', xlab = 'Alfa')
    dev.off()
    
    grafsu<-"K4vXAlfa.jpg"
    
    
    jpeg(grafsu)
    
    plot(v4cum  ~ fim$alfa, ylab='K4 Velocidade', xlab = 'Alfa')
    dev.off()
    
    
    #grafsu<-"TempxmIDerivada.jpg"
    
    
    #jpeg(grafsu)
    
    olhala<-derivadazin(fim$S ,fim$E )
    
    olhabem<-(-1/olhala)
    
    
    #plot(olhabem ~ head(fim$temp,length(fim$temp)-1) , ylab='Temperatura', xlab = '-dE/dS')
    #dev.off()
    
    
    grafsu<-"V2cumxmIDerivada.jpg"
    
    
    jpeg(grafsu)
    
    
    plot(olhabem  ~ head(v2cum,length(v2cum)-1) , ylab='Temperatura', xlab = '-dE/dS')
    dev.off()
    
    
    
    #experimento<-c(fim$E,fim$S,fim$alfa,x1cum,x2cum,x3cum,x4cum,v1cum,v2cum,v3cum,v4cum, fim$temp)
    
    experimento<-data.frame(fim$E,fim$S,fim$alfa,x1cum,x2cum,x3cum,x4cum,v1cum,v2cum,v3cum,v4cum, fim$temp)
    colnames(experimento) <- c('Energia','Entropia', 'Alfa','x1cum','x2cum','x3cum','x4cum','v1cum','v2cum','v3cum','v4cum', 'temperatura')
    
    
    write.csv(experimento, "experimento.csv" )
    
    
    
}



salvasimplesXVKurtlogAlfacte<-function(fim){
    
    
    grafe<-"ExXA.jpg"
    
    jpeg(grafe)
    
    plot(fim$E ~ fim$xa, xlab='Amplitude', ylab='Energia Media')
    lines(fim$PolyE ~ fim$xa, col='blue')
    dev.off()
    
    grafs<-"SxXA.jpg"
    
    jpeg(grafs)
    
    plot(fim$S  ~ fim$xa, xlab='Amplitude', ylab='Entropia')
    lines(fim$PolyS ~ fim$xa, col='blue')
    dev.off()
    
    
    grafse<-"SxE.jpg"
    
    jpeg(grafse)
    
    plot(fim$S  ~ fim$E, xlab='Energia Media', ylab='Entropia')
    dev.off()
    
    x1cum<-pegaordemK(fim$cumulantesX,1)
    x2cum<-pegaordemK(fim$cumulantesX,2)
    x3cum<-pegaordemK(fim$cumulantesX,3)
    x4cum<-pegaordemK(fim$cumulantesX,4)
    
    v1cum<-pegaordemK(fim$cumulantesV,1)
    v2cum<-pegaordemK(fim$cumulantesV,2)
    v3cum<-pegaordemK(fim$cumulantesV,3)
    v4cum<-pegaordemK(fim$cumulantesV,4)
    
    grafoe<-"K1xXA.jpg"
    
    
    jpeg(grafoe)
    
    plot(x1cum  ~ fim$xa, ylab='K1 Posicao', xlab = 'Amplitude')
    dev.off()
    
    grafoe<-"K2xXA.jpg"
    
    
    jpeg(grafoe)
    
    plot(x2cum  ~ fim$xa, ylab='K2 Posicao', xlab = 'Amplitude')
    dev.off()
    
    grafoe<-"K3xXA.jpg"
    
    
    jpeg(grafoe)
    
    plot(x3cum  ~ fim$xa, ylab='K3 Posicao', xlab = 'Amplitude')
    dev.off()
    
    grafoe<-"K4xXA.jpg"
    
    
    jpeg(grafoe)
    
    plot(x4cum  ~ fim$xa, ylab='K4 Posicao', xlab = 'Amplitude')
    dev.off()
    
    
    grafsu<-"K1vXA.jpg"
    
    
    jpeg(grafsu)
    
    plot(v1cum  ~ fim$xa, ylab='K1 Velocidade', xlab = 'Amplitude')
    dev.off()
    
    grafsu<-"K2vXA.jpg"
    
    
    jpeg(grafsu)
    
    plot(v2cum  ~ fim$xa, ylab='K2 Velocidade', xlab = 'Amplitude')
    dev.off()
    
    grafsu<-"K3vXA.jpg"
    
    
    jpeg(grafsu)
    
    plot(v3cum  ~ fim$xa, ylab='K3 Velocidade', xlab = 'Amplitude')
    dev.off()
    
    grafsu<-"K4vXA.jpg"
    
    
    jpeg(grafsu)
    
    plot(v4cum  ~ fim$xa, ylab='K4 Velocidade', xlab = 'Amplitude')
    dev.off()
    

#grafsu<-"TempxmIDerivada.jpg"
    
    
    #   jpeg(grafsu)
    
    olhala<-derivadazin(fim$S ,fim$E )
    
    olhabem<-(-1/olhala)
    
    
    # plot(olhabem ~ rep(fim$temp,length(fim$ns)-1) , ylab='Temperatura', xlab = '-dE/dS')
    #dev.off()
    
    
    grafsu<-"V2cumxmIDerivada.jpg"
    
    
    jpeg(grafsu)
    
    
    plot(olhabem  ~ head(v2cum,length(v2cum)-1) ,  ylab = '-dE/dS', xlab='Temperatura (<v2>)')
    dev.off()
    
    
    
    #experimento<-c(fim$E,fim$S,fim$alfa,x1cum,x2cum,x3cum,x4cum,v1cum,v2cum,v3cum,v4cum, fim$temp)
    
    experimento<-data.frame(fim$E,fim$S,fim$alfa,x1cum,x2cum,x3cum,x4cum,v1cum,v2cum,v3cum,v4cum, fim$temp)
    colnames(experimento) <- c('E','S', 'alfa','x1cum','x2cum','x3cum','x4cum','v1cum','v2cum','v3cum','v4cum')
    
    
    write.csv(experimento, "experimento.csv" )
    
    
    
}

salvasimplesXVKurtlogAlfacteT<-function(fim){
    
    
    grafe<-"ExTemp.jpg"
    
    jpeg(grafe)
    
    plot(fim$E ~ fim$temp, xlab='Temperatura', ylab='Energia Media')
    lines(fim$PolyE ~ fim$temp, col='blue')
    dev.off()
    
    grafs<-"SxTemp.jpg"
    
    jpeg(grafs)
    
    plot(fim$S  ~ fim$temp, xlab='Temperatura', ylab='Entropia')
    lines(fim$PolyS ~ fim$temp, col='blue')
    dev.off()
    
    
    grafse<-"SxE.jpg"
    
    jpeg(grafse)
    
    plot(fim$S  ~ fim$E, xlab='Energia Media', ylab='Entropia')
    dev.off()
    
    x1cum<-pegaordemK(fim$cumulantesX,1)
    x2cum<-pegaordemK(fim$cumulantesX,2)
    x3cum<-pegaordemK(fim$cumulantesX,3)
    x4cum<-pegaordemK(fim$cumulantesX,4)
    
    v1cum<-pegaordemK(fim$cumulantesV,1)
    v2cum<-pegaordemK(fim$cumulantesV,2)
    v3cum<-pegaordemK(fim$cumulantesV,3)
    v4cum<-pegaordemK(fim$cumulantesV,4)
    
    grafoe<-"K1xTemp.jpg"
    
    
    jpeg(grafoe)
    
    plot(x1cum  ~ fim$temp, ylab='K1 Posicao', xlab = 'Temperatura')
    dev.off()
    
    grafoe<-"K2xTemp.jpg"
    
    
    jpeg(grafoe)
    
    plot(x2cum  ~ fim$temp, ylab='K2 Posicao', xlab = 'Temperatura')
    dev.off()
    
    grafoe<-"K3xTemp.jpg"
    
    
    jpeg(grafoe)
    
    plot(x3cum  ~ fim$temp, ylab='K3 Posicao', xlab = 'Temperatura')
    dev.off()
    
    grafoe<-"K4xTemp.jpg"
    
    
    jpeg(grafoe)
    
    plot(x4cum  ~ fim$temp, ylab='K4 Posicao', xlab = 'Temperatura')
    dev.off()
    
    
    grafsu<-"K1vTemp.jpg"
    
    
    jpeg(grafsu)
    
    plot(v1cum  ~ fim$temp, ylab='K1 Velocidade', xlab = 'Temperatura')
    dev.off()
    
    grafsu<-"K2vTemp.jpg"
    
    
    jpeg(grafsu)
    
    plot(v2cum  ~ fim$temp, ylab='K2 Velocidade', xlab = 'Temperatura')
    dev.off()
    
    grafsu<-"K3vTemp.jpg"
    
    
    jpeg(grafsu)
    
    plot(v3cum  ~ fim$temp, ylab='K3 Velocidade', xlab = 'Temperatura')
    dev.off()
    
    grafsu<-"K4vTemp.jpg"
    
    
    jpeg(grafsu)
    
    plot(v4cum  ~ fim$temp, ylab='K4 Velocidade', xlab = 'Temperatura')
    dev.off()
    
    
    #grafsu<-"TempxmIDerivada.jpg"
    
    
    #   jpeg(grafsu)
    
    olhala<-derivadazin(fim$S ,fim$E )
    
    olhabem<-(-1/olhala)
    
    
    # plot(olhabem ~ rep(fim$temp,length(fim$ns)-1) , ylab='Temperatura', xlab = '-dE/dS')
    #dev.off()
    
    
    grafsu<-"V2cumxmIDerivada.jpg"
    
    
    jpeg(grafsu)
    
    
    plot(olhabem  ~ head(v2cum,length(v2cum)-1) ,  ylab = '-dE/dS', xlab='Temperatura (<v2>)')
    dev.off()
    
    
    
    #experimento<-c(fim$E,fim$S,fim$alfa,x1cum,x2cum,x3cum,x4cum,v1cum,v2cum,v3cum,v4cum, fim$temp)
    
    experimento<-data.frame(fim$E,fim$S,fim$alfa,fim$xa ,fim$temp ,  x1cum,x2cum,x3cum,x4cum,v1cum,v2cum,v3cum,v4cum, fim$temp)
    colnames(experimento) <- c('E','S',  'alfa',  'xa' ,'temp',  'x1cum','x2cum','x3cum','x4cum','v1cum','v2cum','v3cum','v4cum')
    
    
    write.csv(experimento, "experimento.csv" )
    
    
    
}

SalvaCorrelacao<-function(listasim, nome)
{
    ####
    # final<-list("PolyE"=analiticE, "E"=energia, "S"=entropia, "PolyS"= B , "alfa"=ddalfa, "nns" = nns, "ns"=ns, "nt"=nt, "nc"=nc, "dt"= dt, "temp" = temp, "erroE" = erroE)
    experimento<-data.frame(listasim )
    #colnames(experimento) <- c('R' )
    write.csv(experimento, nome , row.names = FALSE)
    #, col.names=FALSE)
    
    #write.table(experimento, file = nome ,row.names=FALSE, na="",col.names=FALSE, sep=",")
}

Salva<-function(listasim, nome)
{
    
    experimento<-data.frame(listasim )

    write.csv(experimento, nome , row.names = FALSE)
   
}


SalvaExperimento<-function(listasim)
{
    # final<-list("PolyE"=analiticE, "E"=energia, "S"=entropia, "PolyS"= B , "alfa"=ddalfa, "nns" = nns, "ns"=ns, "nt"=nt, "nc"=nc, "dt"= dt, "temp" = temp, "erroE" = erroE)
    experimento<-data.frame(listasim$E, listasim$S, listasim$alfa , listasim$temp , listasim$erroE)
    colnames(experimento) <- c('E','S',  'alfa', 'temp',  'erroE' )
    
    
    write.csv(experimento, "experimento.csv" )
}

RetroExperimento<-function()
{
    
    qero<-read.csv(experimento)
    
    final<-list("E"= qero$Energia, "E"= qero$Entropia, "alfa"= qero$Alfa )
    
    return(final)
    
}


TMarconi<-function(a,alfa,m,k,g)
{
    T<-(alfa* m* a^2)/(g *(k + alfa * (g + m * alfa) ) )
    return(T)
}


##x <- 1:10
#y <- c(2,4,6,8,7,8,14,16,18,20)
#lo <- loess(y~x)
#plot(x,y)
#xl <- seq(min(x),max(x), (max(x) - min(x))/1000)
#lines(xl, predict(lo,xl), col='red', lwd=2)

