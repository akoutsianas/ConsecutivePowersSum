def modular_method_k_5(primes_bound=30):
  """
    INPUT:
      - primes_bounds: an integer which is the bound of primes we use in the elimination step
      
    OUTPUT:
      A list of small primes over the modular method fails
  """

  E = lambda t: EllipticCurve([0, 5*(20*t+3), 0, 5*(5*(20*t + 3)**2 - 1)/4, 0])
  Snew = Newforms(Gamma0(200), 2)
  
  for newf in Snew:
    Pqs = []
    for q in primes(primes_bound+1):
      if q not in [2,5]:
        print("q: {}".format(q))
        prodq = 1
        aqfnew = newf[q]
        for t0 in range(q):
          # print("t0: {}".format(t0))
          y2p = (5*(20*ZZ(t0) + 3)**2 - 1)/4
          if y2p % q != 0:
            # print("Et0: {}".format(E(ZZ(t0))))
            aqE = q + 1 - E(ZZ(t0)).reduction(q).order()
            prodq *= (aqE - aqfnew)
          else:
            prodq *= (aqfnew**2 - (q + 1)**2)
        print("prodq: {}".format(prodq))
        Pqs.append(prodq)
    print("Pqs: {}".format(Pqs))
        