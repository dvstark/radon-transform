function von_mises_line,x,p,dp

  mu=p[0]
  k=p[1]
  a=p[2]
  b=p[3]
  slope = p[4]

  ;s=p[3]
  ;b=p[3]

  f = a*exp(k*cos(x-mu))/(2*!pi*beseli(k,0))+(b+slope*x);*(1+s*sin(x-mu));+b

  if n_params() gt 2 then begin

     print,'calculating derivatives'
  
     dfdmu = -k*sin(mu-x)*exp(k*cos(x-mu))/(2*!pi*beseli(k,0))
     
     dfdk = beseli(k,1)*exp(k*cos(mu-x))/(2*!pi*beseli(k,0)^2)

     dfda = exp(k*cos(x-mu))/(2*!pi*beseli(k,0))

     dfdb = 0
     
    ; dpds = 0
     
     dp = [[dfdmu],[dfdk],[dfda],[dfdb]];,[dfds]]

  endif

  return,f
  
end


