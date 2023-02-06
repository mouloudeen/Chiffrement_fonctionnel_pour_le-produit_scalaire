load("GroupGen.sage")

print("Implémentation du schéma basé sur DDH du 3e article")

## Setup

def Setup(l, K):
# retourne une clef maître secrète msk et une clef maître publique mpk
# correspondantes aux paramètres K et l

    G, q, g = GroupGen(K)
    Zq = IntegerModRing(q)
    msk = []
    mpk = [] 

    while True:
        g = G.random_element()
        if (g != G(0)) and g.multiplicative_order() == q-1: break
    while True:
        h = G.random_element()
        if (h != G(0)) and h.multiplicative_order() == q-1: break

    for i in range(l):
        (s, t) = (Zq.random_element(), Zq.random_element())
        msk.append((s, t))
        mpk.append(g^s*h^t)

    pub = (G, q, g, h, l) # informations publiques
    
    return(pub, mpk, msk)


## KeyDer

def KeyDer(msk, y):
# retourne la clef secrète dérivée de y par la clef maître secrète msk 
    
    s = vector([ZZ(i[0]) for i in msk])
    t = vector([ZZ(i[1]) for i in msk])
    
    sy = vector(y).inner_product(s)
    ty = vector(y).inner_product(t)
    
    sky = y, sy, ty

    return(sky)

## Encrypt 

def Encrypt(pub, mpk, x):
# retourne un chiffré de x par la clé maître publique mpk

    # definition du contexte
    G, p, g, h, l = pub
    Zq = IntegerModRing(G.order())

    # chiffrement
    r = Zq.random_element()
    C = g^r; D = h^r; E = []
    for i in range(l):
        E.append(g^x[i]*mpk[i]^r)
    Ct = [C]+[D]+E

    return(Ct)


## Decrypt 

def Decrypt(pub, mpk, sky, Ct):
# retourne le déchiffré de Ct par la clef secrète dérivée sky

    # mise du contexte
    G, p, g, h, l = pub
    y,sy,ty = sky
    C,D,E = Ct[0],Ct[1],Ct[2:]
    Ey = 1
    
    # dechiffrement
    for i in range(l):
        Ey *= E[i]^y[i]
    Ey /= C^(sy)*D^(ty)

    return log(Ey,g)


K = 50; l = 20
pub, mpk, msk = Setup(K, l)
G, p, g, h, l = pub

x = [randint(0, p) for _ in range(l)]
Ct = Encrypt(pub, mpk, x)
y = [randint(0, p) for _ in range(l)]
sky = KeyDer(msk, y)
x_y = Decrypt(pub, mpk, sky, Ct)

print(x_y, "correct ?", x_y == (vector(x).inner_product(vector(y)))%(p-1))
