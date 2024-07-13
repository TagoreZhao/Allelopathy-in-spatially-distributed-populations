###2-species model###


rm(list = ls()) # clears the workspace
dev.off() # clears plots
cat("\014") # clears the console


# Define parameters (I'm taking these from the text, just under equation 2)
a1=3 #intrinsic birth rate of type 1
b1=1 #natural death rate of type 1
a2=4 #intrinsic birth rate of type 2
b2=1 #natural death rate of type 2
c=3 #the rate at which type 1 poisons type 2

dt=0.1 # step size
T=10 # time horizon of the simulation
MaxTime=T/dt # total number of time steps


# You can look at Figure 1 for initial conditions between 0 and 1
# Depending on where you start, you can follow the arrows in the figure to see which equilibrium (the points in the figure) you will end at.
# We haven't really talked about it, but unless you start at the coexistence equilibrium, you should see one species eventually die out.

# initial conditions
u1=0.5
u2=0.5


# Create two matrixes (I do not know how to do that)
# We can create just one matrix to hold the data.

#u1<-matrix(data = NA,nrow = 1, ncol = 1, byrow = FALSE, dimnames = NULL)
#u2<-matrix(data = NA,nrow = 1, ncol = 1, byrow = FALSE, dimnames = NULL)
U<-matrix(0,nrow=MaxTime+1,ncol=3)


#set up the system and storage the value of u1 and u2???
#a1*u1*(1-u1-u2)-b1*u1=0
#a2*u2*(1-u1-u2)-b2*u2-c*u1*u2=0


# Set up the initial conditions for the simulation
U[1,1]<-0 # set initial time step
U[1,2]<-u1 # sets u1
U[1,3]<-u2 # sets u2


# The main iteration
for (t in 1:MaxTime+1){
  U[t,1] <- U[t-1,1] + dt # time value
  
  U[t,2] <- U[t-1,2] + dt*( a1*U[t-1,2]*(1-U[t-1,2]-U[t-1,3]) - b1*U[t-1,2] ) # u1
  U[t,3] <- U[t-1,3] + dt*( a2*U[t-1,3]*(1-U[t-1,2]-U[t-1,3]) - b2*U[t-1,3] - c*U[t-1,2]*U[t-1,3] ) # u2
  
  # extinction thresholds
  if(U[t,2]<0.0001){U[t,2]=0}
  if(U[t,3]<0.0001){U[t,3]=0}
}
plot(U[,1],U[,2],col="black",xlab="time (t)",ylab="population size",ylim=c(0,1))
points(U[,1],U[,3],col="red")


# summary table of our analytical results and the simulation
equil_summary=matrix(0,nrow=4,ncol=2)
equil_summary[1,1]=1 - b1/a1 # u2 dies out
equil_summary[2,2]=1 - b2/a2 # u1 dies out

equil_summary[3,1]=(a2/a1)*(b1/c) - (b2/c) # coexistence, u1
equil_summary[3,2]=1 - (b1/a1) + (b2/a2) - (a2/a1)*(b1/c) # coexistence, u2

equil_summary[4,]=U[MaxTime+1,2:3] # end result of the simulation
print(equil_summary)



###3-species model###

rm(list = ls()) # clears the workspace
dev.off() # clears plots
cat("\014") # clears the console


# Define parameters (I'm taking these from the text, just under equation 2)
a1=3 #intrinsic birth rate of type 1
b1=1 #natural death rate of type 1
a2=4 #intrinsic birth rate of type 2
b2=1 #natural death rate of type 2
a3=5 #intrinsic birth rate of type 3
b3=1 #natural death rate of type 3
c1=3 #the rate at which type 1 poisons type 2
c2=3 #

dt=0.1 # step size
T=20 # time horizon of the simulation
MaxTime=T/dt # total number of time steps

# initial conditions
## make sure that 1-u1-u2-u3 >= 0
u1=0.4
u2=0.4
u3=0.1
u0=1-(u1+u2+u3)


#set the matrix
U<-matrix(0,nrow=MaxTime+1,ncol=5)
## I increased the number of columns. We have u0,u1,u2,u3, and time


# Set up the initial conditions for the simulation
## We needed to create the main data matrix before setting the initial conditions
## I moed that code up
U[1,1]<-0 # set initial time step
U[1,2]<-u1 # sets u1
U[1,3]<-u2 # sets u2
U[1,4]<-u3 # sets u3
U[1,5]<-1-u1-u2-u3 # sets u3



# The main iteration
for (t in 1:MaxTime+1){
  U[t,1] <- U[t-1,1] + dt # time value
  
  U[t,2] <- U[t-1,2] + dt*( a1*U[t-1,2]*u0-b1*U[t-1,2]) # u1
  U[t,3] <- U[t-1,3] + dt*( a2*U[t-1,3]*u0-b2*U[t-1,3] ) # u2
  
  # U[t,4] <- U[t-1,3] + dt*( a3*U[t-1,4]*u0-b3(b3-c1*U[t-1,2]-c2*U[t-1,3]))
  U[t,4] <- U[t-1,3] + dt*( a3*U[t-1,4]*u0-b3*(b3-c1*U[t-1,2]-c2*U[t-1,3]))
  U[t,5] <- 1-U[t-1,2]-U[t-1,3]-U[t-1,4]
  ## The problem was you were missing a *
  ## I added code for u0
  
  # extinction thresholds
  if(U[t,2]<0.0001){U[t,2]=0}
  if(U[t,3]<0.0001){U[t,3]=0}
  if(U[t,4]<0.0001){U[t,4]=0}
}
plot(U[,1],U[,2],col="black",xlab="time (t)",ylab="population size",ylim=c(0,1))
points(U[,1],U[,3],col="red")
points(U[,1],U[,4],col="blue")


equil_summary=matrix(0,nrow=2,ncol=3)
equil_summary[1,1]=b1/a1*(1+a3/c1)-b3/c1
equil_summary[1,2]=0
equil_summary[1,3]=1+b3/c1-b1/a1*(1+a3/c1)

equil_summary[2,]=U[MaxTime+1,2:3] # end result of the simulation
print(equil_summary)

