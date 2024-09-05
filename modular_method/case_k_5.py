def modular_method_k_5(primes_bound=30):
  """
    INPUT:
      - primes_bounds: an integer which is the bound of primes we use in the elimination step
      
    OUTPUT:
      A list of small primes over the modular method fails
  """

  y2p = lambda t: QQ(5*(10*t*(t+1) + 3)**2 - 1)/4
  E = lambda t: EllipticCurve([0, 5*(10*t*(t + 1) + 3), 0, 5*y2p(t), 0])
  Snew = Newforms(Gamma0(200), 2)
  
  problematic_fnew = []
  for newf in Snew:
    Pqs = []
    for q in primes(primes_bound+1):
      if q not in [2,5]:
        # print("q: {}".format(q))
        prodq = 1
        aqfnew = newf[q]
        for t0 in range(q):
          # print("t0: {}".format(t0))
          if y2p(ZZ(t0)) % q != 0:
            aqE = q + 1 - E(ZZ(t0)).reduction(q).order()
            prodq *= (aqE - aqfnew)
          else:
            prodq *= (aqfnew**2 - (q + 1)**2)
        # print("prodq: {}\n".format(prodq))
        Pqs.append(prodq)
    if len([c for c in Pqs if c != 0]) == 0:
      problematic_fnew.append(newf)
    print("Pqs: {}".format(Pqs))
  return problematic_fnew