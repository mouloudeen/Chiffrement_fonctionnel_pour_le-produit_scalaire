load("GroupGen.sage")

print("Implémentation du schéma basé sur DDH du 2e article")

## Setup 

def Setup(K, l): 
# retourne une clef maître secrète msk et une clef maître publique mpk
# correspondantes aux paramètres K et l
    
    G, p, g = GroupGen(K)
    msk = []
    mpk = []
    pub = (G, p, g, l) # informations publiques 
    
    for _ in range(l): 
        s = randint(0,p)
        msk.append(s)
        mpk.append(g^s)
        
    return(pub, mpk, msk)


## KeyDer

def KeyDer(msk, y): 
# retourne la clef secrète dérivée de y par la clef maître secrète msk 
    
    return vector(msk).inner_product(vector(y))

## Encrypt 

def Encrypt(pub, mpk, x): 
# retourne un chiffré de x par la clé maître publique mpk
    
    G, p, g, l = pub
    r = randint(0, p)
    
    Ct = [] # chiffré
    Ct.append(g^r) # ct0
    
    for i in range(l): 
        Ct.append((mpk[i]^r)*(g^x[i])) # cti
        
    return(Ct)        


## Decrypt 

def Decrypt(pub, mpk, Ct, sky, y): 
# retourne le déchiffré de Ct par la clef secrète dérivée sky

    G, p, g, l = pub
    
    produit = 1
    for i in range(l): 
        produit = produit*(Ct[i+1]^y[i])
    
    return(log((produit/(Ct[0]^sky)),g))

K = 50
l = 20
pub, mpk, msk = Setup(K, l)
G, p, g, l = pub

B = floor(sqrt(p/l)) # calcul de la borne B telle que l*B^2<p

x = [randint(0, B) for _ in range(l)]
Ct = Encrypt(pub, mpk, x)
y = [randint(0, p) for _ in range(l)]
sky = KeyDer(msk, y)
x_y = Decrypt(pub, mpk, Ct, sky, y)

print(x_y, "correct ?", x_y == (vector(x).inner_product(vector(y)))%(p-1))
