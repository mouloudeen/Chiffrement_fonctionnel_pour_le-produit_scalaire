def GroupGen(K): 
# retourne G, g, p : G groupe cyclique d'ordre p premier de K bits et de générateur g
    
    p = random_prime(2^K, lbound = 2^(K-1))
    G = GF(p)
    g = G.random_element()

    while g == G(0) or g.multiplicative_order() != p-1: 
        g = G.random_element()
    
    return(G, p, g)
