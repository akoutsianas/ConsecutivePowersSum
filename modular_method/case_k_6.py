def modular_method_k_6(primes_bound=30):
    """
      INPUT:
        - primes_bounds: an integer which is the bound of primes we use in the elimination step

      OUTPUT:
        A list of small primes over the modular method fails
    """


    E = lambda x: EllipticCurve([0, 2*χ, 0, (χ**2 + 1), 0])
    N = 2**7 * 3
    Snew = Newforms(Gamma0(N), 2)

    problematic_fnew = []
    for newf in Snew:
        Pqs = []
        for q in primes(primes_bound + 1):
            if q not in [2, 3]:
                # print("q: {}".format(q))
                prodq = 1
                aqfnew = newf[q]
                for x0 in range(q):
                    # print("t0: {}".format(t0))
                    if (x0**2 + 1) % q != 0:
                        aqE = q + 1 - E(x0).reduction(q).order()
                        prodq *= (aqE - aqfnew)
                    else:
                        prodq *= (aqfnew ** 2 - (q + 1) ** 2)
                    if prodq == 0:
                        # print("t0: {}, y2p: {}, q: {}".format(t0, (2**2 * y2p(ZZ(t0)) + 1 - 5*(250*t*(t + 1) + 63)**2) % q, q))
                        break
                # print("prodq: {}\n".format(prodq))
                Pqs.append(prodq)
        if len([c for c in Pqs if c != 0]) == 0:
            problematic_fnew.append(newf)
        print("Pqs: {}".format(gcd(Pqs).factor()))
    return problematic_fnew
