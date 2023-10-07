function k = comb_mn(m,n)
if n <= m
    k = factorial(m) / (factorial(n) * factorial(m-n));
else
  %  'check the value of m and n. m cannot be less than n.'
    k = 0;
end