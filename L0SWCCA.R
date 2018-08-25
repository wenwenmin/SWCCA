SWCCA = function(X, Y, ku, kv, kw, niter=1000, err=0.0001){
  # Input
  # X \in R^{n \times p} (n:samples, p:variables)
  # Y \in R^{n \times q} (n:samples, q:variables)
  
  set.seed(100)
  u0 = matrix(rnorm(ncol(X)),ncol=1);u0 = u0/norm(u0,'E')
  
  set.seed(200)
  v0 = matrix(rnorm(ncol(Y)),ncol=1);v0 = v0/norm(v0,'E')
  
  set.seed(300)
  w0 = matrix(rnorm(nrow(X)),ncol=1);w0 = w0/norm(w0,'E')
  
  u = u0; v = v0; w = w0
  # Iterative algorithm
  for (i in 1:niter){
    z.u = crossprod(X, w*(Y%*%v))
    u = SWCCA.project(z.u, ku)
    
    z.v = crossprod(Y, w*(X%*%u))
    v = SWCCA.project(z.v, kv) 
    
    z.w = (X%*%u)*(Y%*%v)
    w = SWCCA.project(z.w, kw)
    # w = abs(SWCCA.project(z.w, kw))
    # Algorithm termination condition
    if ((norm(u - u0,'E')<= err)&(norm(v - v0,'E')<= err)&(norm(w - w0,'E')<= err)){break}
    else {
      u0=u; v0=v; w0=w}
  }
  obj = sum((X%*%u)*(Y%*%v)*w)
  return (list(u=u, v=v, w=w, obj=obj))
}

# Sparse project function
SWCCA.project = function(z, k){
  #if((length(z)<k)|(sum(z^2)==0)) return(z)
  if(sum(z^2)==0) return(z)
  u = abs(z);
  u[-order(u,decreasing=T)[1:k]] = 0
  u = u/sqrt(sum(u^2))
  u = sign(z)*u
  return(u) 
}
