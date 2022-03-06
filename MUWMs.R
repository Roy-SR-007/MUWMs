rm(list=ls())


# -------------------------------------------------------------------------

# Deletion algorithm for the parallel classes.

delete = function(pi_1,q,deleted_elements)
{
  ct = 1
  flag=1
  x = array(0)
  
  
  for(i in 1:q)
  {
    for(j in 1:q)
    {
      for(k in deleted_elements)
      {
        if(pi_1[i,j] == k)
        {
          flag=0
          break
        } 
      }
      x[ct] = ifelse(flag!=0,pi_1[i,j],NA)
      ct = ct+1
      flag=1
    }
  }
  
  x = na.omit(x)
  
  pi_1_updated = matrix(x,nrow=q,ncol=(q-a),byrow=T)
  return(pi_1_updated)
}

H = matrix(c(1,1,1,-1),nrow=2,byrow=T)
H = (1/sqrt(2))*H

q = 3
a = 1
d = (q-a)*q

X = 1:d
X_bar = 1:(q^2)


# -------------------------------------------------------------------------




# -------------------------------------------------------------------------

# Write a function for this part.

pi_ref = matrix(1:(q^2),ncol=q,byrow=T)
pi_1 = matrix(c(1,5,9,2,6,7,3,4,8),ncol=3,byrow=T)
pi_2 = matrix(c(1,6,8,5,7,3,9,2,4),byrow=T,ncol=3)
pi_3 = matrix(c(1,7,4,6,3,9,8,5,2),ncol=3,byrow=T)

pi_ref_deleted = pi_ref[-q,]

deleted_elements = pi_ref[q,]

pi_1_deleted = delete(pi_1,q,deleted_elements)
pi_2_deleted = delete(pi_2,q,deleted_elements)
pi_3_deleted = delete(pi_3,q,deleted_elements)

pi_list = list(pi_1_deleted,pi_2_deleted,pi_3_deleted)


# -------------------------------------------------------------------------



# -------------------------------------------------------------------------

# Generating the MUWMs.

MUWM = function(d,q,H,pi_list)
{
  I = diag(d)
  MUWM_list = list(0)
  x = array(0)
  x1 = array(0)
  x2 = array(0)
  x3 = array(0)
  
  for(i in 1:q)
  {
    W = matrix(0,ncol=6,nrow=6)
    ct = 1
    M = pi_list[[i]]
    for(j in 1:nrow(M))
    {
      x = M[j,]
      cannonical_bases = I[,x]
      for(j1 in 1:ncol(H))
      {
        x1 = H[,j1]
        for(k in 1:ncol(cannonical_bases))
        {
          x2 = cannonical_bases[,k]
          x3 = x3 + (x1[k]*x2)
        }
        W[,ct] = x3
        ct = ct+1
        x3 = array(0) 
      }
    }
    MUWM_list[[i]] = W
  }
  return(MUWM_list) # Can return the MUWMs directly, by determining,
  # (pi_1)'pi_1 = I, W1 = (pi_1)'pi_2 and W2 = (pi_1)'pi_3
}

M = MUWM(d,q,H,pi_list)

t(M[[1]])%*%M[[1]]

t(M[[1]])%*%M[[2]]

t(M[[1]])%*%M[[3]]


# -------------------------------------------------------------------------


