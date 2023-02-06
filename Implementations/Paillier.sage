from sage.stats.distributions.discrete_gaussian_integer \
 import DiscreteGaussianDistributionIntegerSampler

print("Implémentation du schéma basé sur Paillier")

## Setup 

def Setup(K, l, X, Y):
# retourne une clef maître secrète msk et une clef maître publique mpk 
# correspondantes aux paramètres K, l, X et Y

    # pprime : premier nombre premier aléatoire
    pprime = random_prime(2**K,lbound = 2**(K-1), proof=false) 

    
    p = 2*pprime+1
    
    while (not(p.is_prime())): # on veut p premier
        pprime = random_prime(2**K,lbound = 2**(K-1), proof=false) 
        p = 2*pprime+1

    # qprime : premier nombre premier aléatoire	
    qprime = random_prime(2**K,lbound = 2**(K-1), proof=false) 

    
    q = 2*qprime+1

    while (not(q.is_prime())): # on veut q premier
        qprime = random_prime(2**K,lbound = 2**(K-1), proof=false)
        q = 2*qprime+1

    N = p*q # N est le module composite de Paillier
    
    N2 = N**2 # N2 = N**2
    
    gprime = randint(1,N2)
    
    while gcd(gprime,N2)!=1:
        gprime = randint(1,N2) # gprime aléatoire dans (Z/N2Z)*
        
    g = pow(gprime,2*N,N2) # calcul de g
    
    # calcul de sigma tel que sigma > sqrt(landa)*(N)**(5/2)
    sigma = integer_ceil(sqrt(K)*(N)**(5/2)) + 1 

    
    D = DiscreteGaussianDistributionIntegerSampler(sigma=sigma)
    
    s = [D() for _ in range(l)] # s suit une distribution gausienne
    
    h = [g**(s[i]) for i in range(l)]
    
    mpk = (N, g, h, Y)
    msk = (s, X)
    
    return (mpk, msk)


## KeyDer

def KeyDer(msk, y):
# retourne la clef secrète dérivée de y par la clef maître secrète msk 

    sky = vector(msk[0]).inner_product(vector(y))
    return sky

## Encrypt

def Encrypt(mpk, x):
# retourne un chiffré de x par la clé maître publique mpk

    N = mpk[0]
    g = mpk[1]
    h = mpk[2]
    N2 = N**2
    
    Nfloor = integer_floor(N/4)
    r = randint(0, Nfloor) # r aléatoire entre 0 et N/4
    
    C0 = pow(g,r,N2) # C0
    C = [mod((1+x[i]*N)*h[i]**r, N2) for i in range(len(x))] # reste de Ct
    Ct = [C0]+C
    
    return Ct

## Decrypt

def Decrypt(mpk, sky, Ct, y):
# retourne le déchiffré de Ct par la clef secrète dérivée sky
    
    N = mpk[0]
    N2 = N**2
    C0 = Ct[0]
    
    produit = 1; 
    for i in range(l): 
        produit = produit*(Ct[i+1]^y[i]);
    Cy = mod(C0**(-sky)*produit, N2)
    
    return log(Cy,1+N)


def test_paillier(K, l, X, Y, a):
# permet de faire des tests sur les algorithmes implémentés 

    (N, g, h, Y),(s, X)= Setup(K, l, X, Y)
    x = [randint(0, a) for _ in range(l)]
    skx = KeyDer((s, X), x)
    y = [randint(0, a) for _ in range(l)]
    Cy = Encrypt((N, g, h, Y), y)
    print("correct ?", Decrypt((N, g, h, Y), skx, Cy,x) \
== vector(x).inner_product(vector(y)))
    return (Decrypt((N, g, h, Y), skx, Cy,x), vector(x).inner_product(vector(y)))

'''
K = 25
l = 10
X = 100
Y = 100
a = 10
test_paillier(K, l, X, Y, a)
'''

K = 50
l = 20
X = 100
Y = 100
a = 20
test_paillier(K, l, X, Y, a)

'''
K = 100
l = 40
X = 1000
Y = 1000
a = 50
test_paillier(K, l, X, Y, a)
'''
