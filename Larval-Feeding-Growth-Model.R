#Opposed band - Hydroides - Equal Size Distribution
library(deSolve)
parameters<-c(
              #az= proportion of larvae that go from zygote to stage 1
              #a= amount of carbon needed to advance to next stage
              #f=the total amount of organic carbon in the environment
              #d=death rate of larva
              az=0.8,a=11780972.45,f=23561944.9,d=0.2,
              
              #p=the proportion of carbon in each food size
              p.45=0.125,p1=0.125,p3=0.125,p6=0.125,
              p10=0.125,p15=0.125,p20=0.125,p30=0.125,
              
              #x1=the clearance rate of stage 1 larva on different size food
              x1.45=2.53,x11=5.67,x13=22.28,x16=28.44,
              x110=37.51,x115=14.22,x120=0,x130=0,
              
              #x2=the clearance rate of stage 2 larva on different size food
              x2.45=4.10,x21=10.52,x23=37.30,x26=45.60,
              x210=77.43,x215=58.07,x220=11.46,x230=0,
              
              #x3=the clearance rate of stage 3 larva on different size food
              x3.45=3.49,x31=10.91,x33=37.88,x36=60.49,
              x310=66.67,x315=53.34,x320=14.87,x330=0)
state<-c(z=100,s1=0,s2=0,s3=0,m=0)
diff.EQ.model<-function(t,state,parameters){
  with(as.list(c(state,parameters)),{
    
    #the rate of change from nonfeeding zygotes to stage 1 larvae
    dz<--(z*az)-(d*z)
    
    #rate of change from stage 1 to stage 2
    ds1<-(z*az)-(s1*(((f*p.45*x1.45)/a)+((f*p1*x11)/a)+((f*p3*x13)/a)+
      ((f*p6*x16)/a)+((f*p10*x110)/a)+((f*p15*x115)/a)+
      ((f*p20*x120)/a)+((f*p30*x130)/a)))-(d*s1)
    
    #rate of change from stage 2 to stage 3
    ds2<-(s1*(((f*p.45*x1.45)/a)+((f*p1*x11)/a)+((f*p3*x13)/a)+
      ((f*p6*x16)/a)+((f*p10*x110)/a)+((f*p15*x115)/a)+
      ((f*p20*x120)/a)+((f*p30*x130)/a)))-
      
      (s2*(((f*p.45*x2.45)/a)+((f*p1*x21)/a)+((f*p3*x23)/a)+
      ((f*p6*x26)/a)+((f*p10*x210)/a)+((f*p15*x215)/a)+
      ((f*p20*x220)/a)+((f*p30*x230)/a)))-(d*s2)
    
    #rate of change from stage 3 to metamorphosis
    ds3<-(s2*(((f*p.45*x2.45)/a)+((f*p1*x21)/a)+((f*p3*x23)/a)+
      ((f*p6*x26)/a)+((f*p10*x210)/a)+((f*p15*x215)/a)+
      ((f*p20*x220)/a)+((f*p30*x230)/a)))-
      
      (s3*(((f*p.45*x3.45)/a)+((f*p1*x31)/a)+((f*p3*x33)/a)+
      ((f*p6*x36)/a)+((f*p10*x310)/a)+((f*p15*x315)/a)+
      ((f*p20*x320)/a)+((f*p30*x330)/a)))-(d*s3)
    
    #rate of change of the number of metamorphosed adults
    #metamorphosed adults do not die in this model in order to
    #be able to see the final number that do metamorphose
    dm<-(s3*(((f*p.45*x3.45)/a)+((f*p1*x31)/a)+((f*p3*x33)/a)+
      ((f*p6*x36)/a)+((f*p10*x310)/a)+((f*p15*x315)/a)+
      ((f*p20*x320)/a)+((f*p30*x330)/a)))
    
    list(c(dz,ds1,ds2,ds3,dm))
  })
}
times<-seq(0,20,by=1)
output.1<-ode(y=state, times=times, func=diff.EQ.model,parms=parameters)


#plots the changes in different stage population sizes
#orange is zygotes, red is stage one, green is stage two, blue is stage 3,
#purple is metamorphosis
plot((output.1[,1]),(output.1[,2]),xlim=c(1,20),ylim=c(1,100),type="l",col="orange",lwd=3,main="Opposed Band at Equal Size Distribution",xlab="Time",ylab="# of Idividuals")
lines((output.1[,1]),(output.1[,3]),col="red",lwd=3)
lines((output.1[,1]),(output.1[,4]),col="green",lwd=3)
lines((output.1[,1]),(output.1[,5]),col="blue",lwd=3)
lines((output.1[,1]),(output.1[,6]),col="purple",lwd=3)

output.1


##########################################################################

#Opposed band - Hydroides - Gaussian Distribution
library(deSolve)
parameters<-c(
  #az= proportion of larvae that go from zygot to stage 1
  #a= amount of carbon needed to advance to next stage
  #f=the total amount of organic carbon in the environment
  #d=death rate of larva
  az=0.8,a=11780972.45,f=23561944.9,d=0.2,
  
  #p=the proportion of carbon in each food size
  p.45=0.0332,p1=0.0664,p3=0.132,p6=0.265,
  p10=0.265,p15=0.132,p20=0.0664,p30=0.0332,
  
  #x1=the clearance rate of stage 1 larva on different size food
  x1.45=2.53,x11=5.67,x13=22.28,x16=28.44,
  x110=37.51,x115=14.22,x120=0,x130=0,
  
  #x2=the clearance rate of stage 2 larva on different size food
  x2.45=4.10,x21=10.52,x23=37.30,x26=45.60,
  x210=77.43,x215=58.07,x220=11.46,x230=0,
  
  #x3=the clearance rate of stage 3 larva on different size food
  x3.45=3.49,x31=10.91,x33=37.88,x36=60.49,
  x310=66.67,x315=53.34,x320=14.87,x330=0)
state<-c(z=100,s1=0,s2=0,s3=0,m=0)
diff.EQ.model<-function(t,state,parameters){
  with(as.list(c(state,parameters)),{
    
    #the rate of change from nonfeeding zygotes to stage 1 larvae
    dz<--(z*az)-(d*z)
    
    #rate of change from stage 1 to stage 2
    ds1<-(z*az)-(s1*(((f*p.45*x1.45)/a)+((f*p1*x11)/a)+((f*p3*x13)/a)+
                       ((f*p6*x16)/a)+((f*p10*x110)/a)+((f*p15*x115)/a)+
                       ((f*p20*x120)/a)+((f*p30*x130)/a)))-(d*s1)
    
    #rate of change from stage 2 to stage 3
    ds2<-(s1*(((f*p.45*x1.45)/a)+((f*p1*x11)/a)+((f*p3*x13)/a)+
                ((f*p6*x16)/a)+((f*p10*x110)/a)+((f*p15*x115)/a)+
                ((f*p20*x120)/a)+((f*p30*x130)/a)))-
      
      (s2*(((f*p.45*x2.45)/a)+((f*p1*x21)/a)+((f*p3*x23)/a)+
             ((f*p6*x26)/a)+((f*p10*x210)/a)+((f*p15*x215)/a)+
             ((f*p20*x220)/a)+((f*p30*x230)/a)))-(d*s2)
    
    #rate of change from stage 3 to metamorphosis
    ds3<-(s2*(((f*p.45*x2.45)/a)+((f*p1*x21)/a)+((f*p3*x23)/a)+
                ((f*p6*x26)/a)+((f*p10*x210)/a)+((f*p15*x215)/a)+
                ((f*p20*x220)/a)+((f*p30*x230)/a)))-
      
      (s3*(((f*p.45*x3.45)/a)+((f*p1*x31)/a)+((f*p3*x33)/a)+
             ((f*p6*x36)/a)+((f*p10*x310)/a)+((f*p15*x315)/a)+
             ((f*p20*x320)/a)+((f*p30*x330)/a)))-(d*s3)
    
    #rate of change of the number of metamorphosed adults
    #metamorphosed adults do not die in this model in order to
    #be able to see the final number that do metamorphose
    dm<-(s3*(((f*p.45*x3.45)/a)+((f*p1*x31)/a)+((f*p3*x33)/a)+
               ((f*p6*x36)/a)+((f*p10*x310)/a)+((f*p15*x315)/a)+
               ((f*p20*x320)/a)+((f*p30*x330)/a)))
    
    list(c(dz,ds1,ds2,ds3,dm))
  })
}
times<-seq(0,20,by=1)
output.1<-ode(y=state, times=times, func=diff.EQ.model,parms=parameters)


#plots the changes in different stage population sizes
#orange is zygotes, red is stage one, green is stage two, blue is stage 3,
#purple is metamorphosis
plot((output.1[,1]),(output.1[,2]),xlim=c(1,20),ylim=c(1,100),type="l",col="orange",lwd=3,main="Opposed Band at Gaussian Distribution",xlab="Time",ylab="# of Idividuals")
lines((output.1[,1]),(output.1[,3]),col="red",lwd=3)
lines((output.1[,1]),(output.1[,4]),col="green",lwd=3)
lines((output.1[,1]),(output.1[,5]),col="blue",lwd=3)
lines((output.1[,1]),(output.1[,6]),col="purple",lwd=3)

output.1


##########################################################################

#Opposed band - Hydroides - Large Particle Bias
library(deSolve)
parameters<-c(
  #az= proportion of larvae that go from zygot to stage 1
  #a= amount of carbon needed to advance to next stage
  #f=the total amount of organic carbon in the environment
  #d=death rate of larva
  az=0.8,a=11780972.45,f=23561944.9,d=0.2,
  
  #p=the proportion of carbon in each food size
  p.45=0.0332,p1=0.0332,p3=0.0664,p6=0.0664,
  p10=0.132,p15=0.132,p20=0.265,p30=0.265,
  
  #x1=the clearance rate of stage 1 larva on different size food
  x1.45=2.53,x11=5.67,x13=22.28,x16=28.44,
  x110=37.51,x115=14.22,x120=0,x130=0,
  
  #x2=the clearance rate of stage 2 larva on different size food
  x2.45=4.10,x21=10.52,x23=37.30,x26=45.60,
  x210=77.43,x215=58.07,x220=11.46,x230=0,
  
  #x3=the clearance rate of stage 3 larva on different size food
  x3.45=3.49,x31=10.91,x33=37.88,x36=60.49,
  x310=66.67,x315=53.34,x320=14.87,x330=0)
state<-c(z=100,s1=0,s2=0,s3=0,m=0)
diff.EQ.model<-function(t,state,parameters){
  with(as.list(c(state,parameters)),{
    
    #the rate of change from nonfeeding zygotes to stage 1 larvae
    dz<--(z*az)-(d*z)
    
    #rate of change from stage 1 to stage 2
    ds1<-(z*az)-(s1*(((f*p.45*x1.45)/a)+((f*p1*x11)/a)+((f*p3*x13)/a)+
                       ((f*p6*x16)/a)+((f*p10*x110)/a)+((f*p15*x115)/a)+
                       ((f*p20*x120)/a)+((f*p30*x130)/a)))-(d*s1)
    
    #rate of change from stage 2 to stage 3
    ds2<-(s1*(((f*p.45*x1.45)/a)+((f*p1*x11)/a)+((f*p3*x13)/a)+
                ((f*p6*x16)/a)+((f*p10*x110)/a)+((f*p15*x115)/a)+
                ((f*p20*x120)/a)+((f*p30*x130)/a)))-
      
      (s2*(((f*p.45*x2.45)/a)+((f*p1*x21)/a)+((f*p3*x23)/a)+
             ((f*p6*x26)/a)+((f*p10*x210)/a)+((f*p15*x215)/a)+
             ((f*p20*x220)/a)+((f*p30*x230)/a)))-(d*s2)
    
    #rate of change from stage 3 to metamorphosis
    ds3<-(s2*(((f*p.45*x2.45)/a)+((f*p1*x21)/a)+((f*p3*x23)/a)+
                ((f*p6*x26)/a)+((f*p10*x210)/a)+((f*p15*x215)/a)+
                ((f*p20*x220)/a)+((f*p30*x230)/a)))-
      
      (s3*(((f*p.45*x3.45)/a)+((f*p1*x31)/a)+((f*p3*x33)/a)+
             ((f*p6*x36)/a)+((f*p10*x310)/a)+((f*p15*x315)/a)+
             ((f*p20*x320)/a)+((f*p30*x330)/a)))-(d*s3)
    
    #rate of change of the number of metamorphosed adults
    #metamorphosed adults do not die in this model in order to
    #be able to see the final number that do metamorphose
    dm<-(s3*(((f*p.45*x3.45)/a)+((f*p1*x31)/a)+((f*p3*x33)/a)+
               ((f*p6*x36)/a)+((f*p10*x310)/a)+((f*p15*x315)/a)+
               ((f*p20*x320)/a)+((f*p30*x330)/a)))
    
    list(c(dz,ds1,ds2,ds3,dm))
  })
}
times<-seq(0,20,by=1)
output.1<-ode(y=state, times=times, func=diff.EQ.model,parms=parameters)


#plots the changes in different stage population sizes
#orange is zygotes, red is stage one, green is stage two, blue is stage 3,
#purple is metamorphosis
plot((output.1[,1]),(output.1[,2]),xlim=c(1,20),ylim=c(1,100),type="l",col="orange",lwd=3,main="Opposed Band at Large Particle Bias",xlab="Time",ylab="# of Idividuals")
lines((output.1[,1]),(output.1[,3]),col="red",lwd=3)
lines((output.1[,1]),(output.1[,4]),col="green",lwd=3)
lines((output.1[,1]),(output.1[,5]),col="blue",lwd=3)
lines((output.1[,1]),(output.1[,6]),col="purple",lwd=3)

output.1


##########################################################################

#Opposed band - Hydroides - Small Particle Bias
library(deSolve)
parameters<-c(
  #az= proportion of larvae that go from zygot to stage 1
  #a= amount of carbon needed to advance to next stage
  #f=the total amount of organic carbon in the environment
  #d=death rate of larva
  az=0.8,a=11780972.45,f=23561944.9,d=0.2,
  
  #p=the proportion of carbon in each food size
  p.45=0.265,p1=0.265,p3=0.132,p6=0.132,
  p10=0.0664,p15=0.0664,p20=0.0332,p30=0.0332,
  
  #x1=the clearance rate of stage 1 larva on different size food
  x1.45=2.53,x11=5.67,x13=22.28,x16=28.44,
  x110=37.51,x115=14.22,x120=0,x130=0,
  
  #x2=the clearance rate of stage 2 larva on different size food
  x2.45=4.10,x21=10.52,x23=37.30,x26=45.60,
  x210=77.43,x215=58.07,x220=11.46,x230=0,
  
  #x3=the clearance rate of stage 3 larva on different size food
  x3.45=3.49,x31=10.91,x33=37.88,x36=60.49,
  x310=66.67,x315=53.34,x320=14.87,x330=0)
state<-c(z=100,s1=0,s2=0,s3=0,m=0)
diff.EQ.model<-function(t,state,parameters){
  with(as.list(c(state,parameters)),{
    
    #the rate of change from nonfeeding zygotes to stage 1 larvae
    dz<--(z*az)-(d*z)
    
    #rate of change from stage 1 to stage 2
    ds1<-(z*az)-(s1*(((f*p.45*x1.45)/a)+((f*p1*x11)/a)+((f*p3*x13)/a)+
                       ((f*p6*x16)/a)+((f*p10*x110)/a)+((f*p15*x115)/a)+
                       ((f*p20*x120)/a)+((f*p30*x130)/a)))-(d*s1)
    
    #rate of change from stage 2 to stage 3
    ds2<-(s1*(((f*p.45*x1.45)/a)+((f*p1*x11)/a)+((f*p3*x13)/a)+
                ((f*p6*x16)/a)+((f*p10*x110)/a)+((f*p15*x115)/a)+
                ((f*p20*x120)/a)+((f*p30*x130)/a)))-
      
      (s2*(((f*p.45*x2.45)/a)+((f*p1*x21)/a)+((f*p3*x23)/a)+
             ((f*p6*x26)/a)+((f*p10*x210)/a)+((f*p15*x215)/a)+
             ((f*p20*x220)/a)+((f*p30*x230)/a)))-(d*s2)
    
    #rate of change from stage 3 to metamorphosis
    ds3<-(s2*(((f*p.45*x2.45)/a)+((f*p1*x21)/a)+((f*p3*x23)/a)+
                ((f*p6*x26)/a)+((f*p10*x210)/a)+((f*p15*x215)/a)+
                ((f*p20*x220)/a)+((f*p30*x230)/a)))-
      
      (s3*(((f*p.45*x3.45)/a)+((f*p1*x31)/a)+((f*p3*x33)/a)+
             ((f*p6*x36)/a)+((f*p10*x310)/a)+((f*p15*x315)/a)+
             ((f*p20*x320)/a)+((f*p30*x330)/a)))-(d*s3)
    
    #rate of change of the number of metamorphosed adults
    #metamorphosed adults do not die in this model in order to
    #be able to see the final number that do metamorphose
    dm<-(s3*(((f*p.45*x3.45)/a)+((f*p1*x31)/a)+((f*p3*x33)/a)+
               ((f*p6*x36)/a)+((f*p10*x310)/a)+((f*p15*x315)/a)+
               ((f*p20*x320)/a)+((f*p30*x330)/a)))
    
    list(c(dz,ds1,ds2,ds3,dm))
  })
}
times<-seq(0,20,by=1)
output.1<-ode(y=state, times=times, func=diff.EQ.model,parms=parameters)


#plots the changes in different stage population sizes
#orange is zygotes, red is stage one, green is stage two, blue is stage 3,
#purple is metamorphosis
plot((output.1[,1]),(output.1[,2]),xlim=c(1,20),ylim=c(1,100),type="l",col="orange",lwd=3,main="Opposed Band at Small Particle Bias",xlab="Time",ylab="# of Idividuals")
lines((output.1[,1]),(output.1[,3]),col="red",lwd=3)
lines((output.1[,1]),(output.1[,4]),col="green",lwd=3)
lines((output.1[,1]),(output.1[,5]),col="blue",lwd=3)
lines((output.1[,1]),(output.1[,6]),col="purple",lwd=3)

output.1


##########################################################################

#Ciliary Reversal - Dendraster - Equal Size Distribution
library(deSolve)
parameters<-c(
  #az= proportion of larvae that go from zygot to stage 1
  #a= amount of carbon needed to advance to next stage
  #f=the total amount of organic carbon in the environment
  #d=death rate of larva
  az=0.8,a=82446807.16,f=23561944.9,d=0.2,
  
  #p=the proportion of carbon in each food size
  p.45=0.125,p1=0.125,p3=0.125,p6=0.125,
  p10=0.125,p15=0.125,p20=0.125,p30=0.125,
  
  #x1=the clearance rate of stage 1 larva on different size food
  x1.45=0.03,x11=0.03,x13=9.30,x16=47.47,
  x110=87.84,x115=118.55,x120=118.91,x130=94.71,
  
  #x2=the clearance rate of stage 2 larva on different size food
  x2.45=0.04,x21=0.09,x23=17.35,x26=131.79,
  x210=202.70,x215=250.77,x220=317.34,x230=242.28,
  
  #x3=the clearance rate of stage 3 larva on different size food
  x3.45=0.05,x31=0.26,x33=22.06,x36=211.15,
  x310=416.49,x315=563.17,x320=591.63,x330=636.38)

state<-c(z=100,s1=0,s2=0,s3=0,m=0)
diff.EQ.model<-function(t,state,parameters){
  with(as.list(c(state,parameters)),{
    
    #the rate of change from nonfeeding zygotes to stage 1 larvae
    dz<--(z*az)-(d*z)
    
    #rate of change from stage 1 to stage 2
    ds1<-(z*az)-(s1*(((f*p.45*x1.45)/a)+((f*p1*x11)/a)+((f*p3*x13)/a)+
                       ((f*p6*x16)/a)+((f*p10*x110)/a)+((f*p15*x115)/a)+
                       ((f*p20*x120)/a)+((f*p30*x130)/a)))-(d*s1)
    
    #rate of change from stage 2 to stage 3
    ds2<-(s1*(((f*p.45*x1.45)/a)+((f*p1*x11)/a)+((f*p3*x13)/a)+
                ((f*p6*x16)/a)+((f*p10*x110)/a)+((f*p15*x115)/a)+
                ((f*p20*x120)/a)+((f*p30*x130)/a)))-
      
      (s2*(((f*p.45*x2.45)/a)+((f*p1*x21)/a)+((f*p3*x23)/a)+
             ((f*p6*x26)/a)+((f*p10*x210)/a)+((f*p15*x215)/a)+
             ((f*p20*x220)/a)+((f*p30*x230)/a)))-(d*s2)
    
    #rate of change from stage 3 to metamorphosis
    ds3<-(s2*(((f*p.45*x2.45)/a)+((f*p1*x21)/a)+((f*p3*x23)/a)+
                ((f*p6*x26)/a)+((f*p10*x210)/a)+((f*p15*x215)/a)+
                ((f*p20*x220)/a)+((f*p30*x230)/a)))-
      
      (s3*(((f*p.45*x3.45)/a)+((f*p1*x31)/a)+((f*p3*x33)/a)+
             ((f*p6*x36)/a)+((f*p10*x310)/a)+((f*p15*x315)/a)+
             ((f*p20*x320)/a)+((f*p30*x330)/a)))-(d*s3)
    
    #rate of change of the number of metamorphosed adults
    #metamorphosed adults do not die in this model in order to
    #be able to see the final number that do metamorphose
    dm<-(s3*(((f*p.45*x3.45)/a)+((f*p1*x31)/a)+((f*p3*x33)/a)+
               ((f*p6*x36)/a)+((f*p10*x310)/a)+((f*p15*x315)/a)+
               ((f*p20*x320)/a)+((f*p30*x330)/a)))
    
    list(c(dz,ds1,ds2,ds3,dm))
  })
}
times<-seq(0,20,by=1)
output.1<-ode(y=state, times=times, func=diff.EQ.model,parms=parameters)


#plots the changes in different stage population sizes
#orange is zygotes, red is stage one, green is stage two, blue is stage 3,
#purple is metamorphosis
plot((output.1[,1]),(output.1[,2]),xlim=c(1,20),ylim=c(1,100),type="l",col="orange",lwd=3,main="Ciliary Reversal at Equal Size Distribution",xlab="Time",ylab="# of Idividuals")
lines((output.1[,1]),(output.1[,3]),col="red",lwd=3)
lines((output.1[,1]),(output.1[,4]),col="green",lwd=3)
lines((output.1[,1]),(output.1[,5]),col="blue",lwd=3)
lines((output.1[,1]),(output.1[,6]),col="purple",lwd=3)

output.1


##########################################################################

#Ciliary Reversal - Dendraster - Gaussian Distribution
library(deSolve)
parameters<-c(
  #az= proportion of larvae that go from zygot to stage 1
  #a= amount of carbon needed to advance to next stage
  #f=the total amount of organic carbon in the environment
  #d=death rate of larva
  az=0.8,a=82446807.16,f=23561944.9,d=0.2,
  
  #p=the proportion of carbon in each food size
  p.45=0.0332,p1=0.0664,p3=0.132,p6=0.265,
  p10=0.265,p15=0.132,p20=0.0664,p30=0.0332,
  
  #x1=the clearance rate of stage 1 larva on different size food
  x1.45=0.03,x11=0.03,x13=9.30,x16=47.47,
  x110=87.84,x115=118.55,x120=118.91,x130=94.71,
  
  #x2=the clearance rate of stage 2 larva on different size food
  x2.45=0.04,x21=0.09,x23=17.35,x26=131.79,
  x210=202.70,x215=250.77,x220=317.34,x230=242.28,
  
  #x3=the clearance rate of stage 3 larva on different size food
  x3.45=0.05,x31=0.26,x33=22.06,x36=211.15,
  x310=416.49,x315=563.17,x320=591.63,x330=636.38)

state<-c(z=100,s1=0,s2=0,s3=0,m=0)
diff.EQ.model<-function(t,state,parameters){
  with(as.list(c(state,parameters)),{
    
    #the rate of change from nonfeeding zygotes to stage 1 larvae
    dz<--(z*az)-(d*z)
    
    #rate of change from stage 1 to stage 2
    ds1<-(z*az)-(s1*(((f*p.45*x1.45)/a)+((f*p1*x11)/a)+((f*p3*x13)/a)+
                       ((f*p6*x16)/a)+((f*p10*x110)/a)+((f*p15*x115)/a)+
                       ((f*p20*x120)/a)+((f*p30*x130)/a)))-(d*s1)
    
    #rate of change from stage 2 to stage 3
    ds2<-(s1*(((f*p.45*x1.45)/a)+((f*p1*x11)/a)+((f*p3*x13)/a)+
                ((f*p6*x16)/a)+((f*p10*x110)/a)+((f*p15*x115)/a)+
                ((f*p20*x120)/a)+((f*p30*x130)/a)))-
      
      (s2*(((f*p.45*x2.45)/a)+((f*p1*x21)/a)+((f*p3*x23)/a)+
             ((f*p6*x26)/a)+((f*p10*x210)/a)+((f*p15*x215)/a)+
             ((f*p20*x220)/a)+((f*p30*x230)/a)))-(d*s2)
    
    #rate of change from stage 3 to metamorphosis
    ds3<-(s2*(((f*p.45*x2.45)/a)+((f*p1*x21)/a)+((f*p3*x23)/a)+
                ((f*p6*x26)/a)+((f*p10*x210)/a)+((f*p15*x215)/a)+
                ((f*p20*x220)/a)+((f*p30*x230)/a)))-
      
      (s3*(((f*p.45*x3.45)/a)+((f*p1*x31)/a)+((f*p3*x33)/a)+
             ((f*p6*x36)/a)+((f*p10*x310)/a)+((f*p15*x315)/a)+
             ((f*p20*x320)/a)+((f*p30*x330)/a)))-(d*s3)
    
    #rate of change of the number of metamorphosed adults
    #metamorphosed adults do not die in this model in order to
    #be able to see the final number that do metamorphose
    dm<-(s3*(((f*p.45*x3.45)/a)+((f*p1*x31)/a)+((f*p3*x33)/a)+
               ((f*p6*x36)/a)+((f*p10*x310)/a)+((f*p15*x315)/a)+
               ((f*p20*x320)/a)+((f*p30*x330)/a)))
    
    list(c(dz,ds1,ds2,ds3,dm))
  })
}
times<-seq(0,20,by=1)
output.1<-ode(y=state, times=times, func=diff.EQ.model,parms=parameters)


#plots the changes in different stage population sizes
#orange is zygotes, red is stage one, green is stage two, blue is stage 3,
#purple is metamorphosis
plot((output.1[,1]),(output.1[,2]),xlim=c(1,20),ylim=c(1,100),type="l",col="orange",lwd=3,main="Ciliary Reversal at Gaussian Distribution",xlab="Time",ylab="# of Idividuals")
lines((output.1[,1]),(output.1[,3]),col="red",lwd=3)
lines((output.1[,1]),(output.1[,4]),col="green",lwd=3)
lines((output.1[,1]),(output.1[,5]),col="blue",lwd=3)
lines((output.1[,1]),(output.1[,6]),col="purple",lwd=3)

output.1


##########################################################################

#Ciliary Reversal - Dendraster - Large Particle Bias
library(deSolve)
parameters<-c(
  #az= proportion of larvae that go from zygot to stage 1
  #a= amount of carbon needed to advance to next stage
  #f=the total amount of organic carbon in the environment
  #d=death rate of larva
  az=0.8,a=82446807.16,f=23561944.9,d=0.2,
  
  #p=the proportion of carbon in each food size
  p.45=0.0332,p1=0.0332,p3=0.0664,p6=0.0664,
  p10=0.132,p15=0.132,p20=0.265,p30=0.265,
  
  #x1=the clearance rate of stage 1 larva on different size food
  x1.45=0.03,x11=0.03,x13=9.30,x16=47.47,
  x110=87.84,x115=118.55,x120=118.91,x130=94.71,
  
  #x2=the clearance rate of stage 2 larva on different size food
  x2.45=0.04,x21=0.09,x23=17.35,x26=131.79,
  x210=202.70,x215=250.77,x220=317.34,x230=242.28,
  
  #x3=the clearance rate of stage 3 larva on different size food
  x3.45=0.05,x31=0.26,x33=22.06,x36=211.15,
  x310=416.49,x315=563.17,x320=591.63,x330=636.38)

state<-c(z=100,s1=0,s2=0,s3=0,m=0)
diff.EQ.model<-function(t,state,parameters){
  with(as.list(c(state,parameters)),{
    
    #the rate of change from nonfeeding zygotes to stage 1 larvae
    dz<--(z*az)-(d*z)
    
    #rate of change from stage 1 to stage 2
    ds1<-(z*az)-(s1*(((f*p.45*x1.45)/a)+((f*p1*x11)/a)+((f*p3*x13)/a)+
                       ((f*p6*x16)/a)+((f*p10*x110)/a)+((f*p15*x115)/a)+
                       ((f*p20*x120)/a)+((f*p30*x130)/a)))-(d*s1)
    
    #rate of change from stage 2 to stage 3
    ds2<-(s1*(((f*p.45*x1.45)/a)+((f*p1*x11)/a)+((f*p3*x13)/a)+
                ((f*p6*x16)/a)+((f*p10*x110)/a)+((f*p15*x115)/a)+
                ((f*p20*x120)/a)+((f*p30*x130)/a)))-
      
      (s2*(((f*p.45*x2.45)/a)+((f*p1*x21)/a)+((f*p3*x23)/a)+
             ((f*p6*x26)/a)+((f*p10*x210)/a)+((f*p15*x215)/a)+
             ((f*p20*x220)/a)+((f*p30*x230)/a)))-(d*s2)
    
    #rate of change from stage 3 to metamorphosis
    ds3<-(s2*(((f*p.45*x2.45)/a)+((f*p1*x21)/a)+((f*p3*x23)/a)+
                ((f*p6*x26)/a)+((f*p10*x210)/a)+((f*p15*x215)/a)+
                ((f*p20*x220)/a)+((f*p30*x230)/a)))-
      
      (s3*(((f*p.45*x3.45)/a)+((f*p1*x31)/a)+((f*p3*x33)/a)+
             ((f*p6*x36)/a)+((f*p10*x310)/a)+((f*p15*x315)/a)+
             ((f*p20*x320)/a)+((f*p30*x330)/a)))-(d*s3)
    
    #rate of change of the number of metamorphosed adults
    #metamorphosed adults do not die in this model in order to
    #be able to see the final number that do metamorphose
    dm<-(s3*(((f*p.45*x3.45)/a)+((f*p1*x31)/a)+((f*p3*x33)/a)+
               ((f*p6*x36)/a)+((f*p10*x310)/a)+((f*p15*x315)/a)+
               ((f*p20*x320)/a)+((f*p30*x330)/a)))
    
    list(c(dz,ds1,ds2,ds3,dm))
  })
}
times<-seq(0,20,by=1)
output.1<-ode(y=state, times=times, func=diff.EQ.model,parms=parameters)


#plots the changes in different stage population sizes
#orange is zygotes, red is stage one, green is stage two, blue is stage 3,
#purple is metamorphosis
plot((output.1[,1]),(output.1[,2]),xlim=c(1,20),ylim=c(1,100),type="l",col="orange",lwd=3,main="Ciliary Reversal at Large Particle Bias",xlab="Time",ylab="# of Idividuals")
lines((output.1[,1]),(output.1[,3]),col="red",lwd=3)
lines((output.1[,1]),(output.1[,4]),col="green",lwd=3)
lines((output.1[,1]),(output.1[,5]),col="blue",lwd=3)
lines((output.1[,1]),(output.1[,6]),col="purple",lwd=3)

output.1


##########################################################################

#Ciliary Reversal - Dendraster - Small Particle Bias
library(deSolve)
parameters<-c(
  #az= proportion of larvae that go from zygot to stage 1
  #a= amount of carbon needed to advance to next stage
  #f=the total amount of organic carbon in the environment
  #d=death rate of larva
  az=0.8,a=82446807.16,f=23561944.9,d=0.2,
  
  #p=the proportion of carbon in each food size
  p.45=0.265,p1=0.265,p3=0.132,p6=0.132,
  p10=0.0664,p15=0.0664,p20=0.0332,p30=0.0332,
  
  #x1=the clearance rate of stage 1 larva on different size food
  x1.45=0.03,x11=0.03,x13=9.30,x16=47.47,
  x110=87.84,x115=118.55,x120=118.91,x130=94.71,
  
  #x2=the clearance rate of stage 2 larva on different size food
  x2.45=0.04,x21=0.09,x23=17.35,x26=131.79,
  x210=202.70,x215=250.77,x220=317.34,x230=242.28,
  
  #x3=the clearance rate of stage 3 larva on different size food
  x3.45=0.05,x31=0.26,x33=22.06,x36=211.15,
  x310=416.49,x315=563.17,x320=591.63,x330=636.38)

state<-c(z=100,s1=0,s2=0,s3=0,m=0)
diff.EQ.model<-function(t,state,parameters){
  with(as.list(c(state,parameters)),{
    
    #the rate of change from nonfeeding zygotes to stage 1 larvae
    dz<--(z*az)-(d*z)
    
    #rate of change from stage 1 to stage 2
    ds1<-(z*az)-(s1*(((f*p.45*x1.45)/a)+((f*p1*x11)/a)+((f*p3*x13)/a)+
                       ((f*p6*x16)/a)+((f*p10*x110)/a)+((f*p15*x115)/a)+
                       ((f*p20*x120)/a)+((f*p30*x130)/a)))-(d*s1)
    
    #rate of change from stage 2 to stage 3
    ds2<-(s1*(((f*p.45*x1.45)/a)+((f*p1*x11)/a)+((f*p3*x13)/a)+
                ((f*p6*x16)/a)+((f*p10*x110)/a)+((f*p15*x115)/a)+
                ((f*p20*x120)/a)+((f*p30*x130)/a)))-
      
      (s2*(((f*p.45*x2.45)/a)+((f*p1*x21)/a)+((f*p3*x23)/a)+
             ((f*p6*x26)/a)+((f*p10*x210)/a)+((f*p15*x215)/a)+
             ((f*p20*x220)/a)+((f*p30*x230)/a)))-(d*s2)
    
    #rate of change from stage 3 to metamorphosis
    ds3<-(s2*(((f*p.45*x2.45)/a)+((f*p1*x21)/a)+((f*p3*x23)/a)+
                ((f*p6*x26)/a)+((f*p10*x210)/a)+((f*p15*x215)/a)+
                ((f*p20*x220)/a)+((f*p30*x230)/a)))-
      
      (s3*(((f*p.45*x3.45)/a)+((f*p1*x31)/a)+((f*p3*x33)/a)+
             ((f*p6*x36)/a)+((f*p10*x310)/a)+((f*p15*x315)/a)+
             ((f*p20*x320)/a)+((f*p30*x330)/a)))-(d*s3)
    
    #rate of change of the number of metamorphosed adults
    #metamorphosed adults do not die in this model in order to
    #be able to see the final number that do metamorphose
    dm<-(s3*(((f*p.45*x3.45)/a)+((f*p1*x31)/a)+((f*p3*x33)/a)+
               ((f*p6*x36)/a)+((f*p10*x310)/a)+((f*p15*x315)/a)+
               ((f*p20*x320)/a)+((f*p30*x330)/a)))
    
    list(c(dz,ds1,ds2,ds3,dm))
  })
}
times<-seq(0,20,by=1)
output.1<-ode(y=state, times=times, func=diff.EQ.model,parms=parameters)


#plots the changes in different stage population sizes
#orange is zygotes, red is stage one, green is stage two, blue is stage 3,
#purple is metamorphosis
plot((output.1[,1]),(output.1[,2]),xlim=c(1,20),ylim=c(1,100),type="l",col="orange",lwd=3,main="Ciliary Reversal at Small Particle Bias",xlab="Time",ylab="# of Idividuals")
lines((output.1[,1]),(output.1[,3]),col="red",lwd=3)
lines((output.1[,1]),(output.1[,4]),col="green",lwd=3)
lines((output.1[,1]),(output.1[,5]),col="blue",lwd=3)
lines((output.1[,1]),(output.1[,6]),col="purple",lwd=3)

output.1