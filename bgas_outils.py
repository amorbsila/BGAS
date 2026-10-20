"""
bgas_outils.py  --  head h_k, k*B(l) et somme des termes, pour phi, tribonacci, plast.

Conventions
  BGASINV2B(a, bits, b) : u <- round(b*u)+bit, en partant de u=a
  B(l) = bnode de longueur l :  0^l  (phi, trib)  ou  1 0^(l-1)  (trib1, plast)
  base 'phi'   : B = 2,3,5,8,13,...            (F_{l+2})
  base 'trib'  : B = 2,4,7,13,24,...           (T_{l+3}, bnodes 0^l)
  base 'trib1' : B = 3,6,11,20,37,...          (bnodes 10^(l-1))
  base 'plast' : B = 2,3,4,5,7,9,12,...        (P_{l+3}, bnodes 10^(l-1))
"""
from decimal import Decimal, getcontext

getcontext().prec = 400


# ---------------- tes fonctions ----------------
def custom_round_decimal(x):
    return int(x.to_integral_value(rounding="ROUND_HALF_UP"))


def BGAS2BB(att, X, B):
    fn = Decimal(X)
    b = Decimal(str(B))
    bits = []
    while fn > att:
        candidate = Decimal(custom_round_decimal(fn / b))
        while Decimal(custom_round_decimal(candidate * b)) > fn:
            candidate -= 1
        bit = 1 if Decimal(custom_round_decimal(candidate * b)) < fn else 0
        bits.append(bit)
        fn = candidate
    bits.reverse()
    return bits


def BGASINV2B(a, bits, b):
    Bse = Decimal(str(b))
    current = Decimal(a)
    for bit in bits:
        current = Decimal(custom_round_decimal(current * Bse)) + Decimal(bit)
    return int(current)


# ---------------- les bases ----------------
def _newton(f, df, x0):
    x = Decimal(x0)
    for _ in range(600):
        x = x - f(x) / df(x)
    return x


PHI = (1 + Decimal(5).sqrt()) / 2
TAU = _newton(lambda x: x**3 - x**2 - x - 1, lambda x: 3 * x**2 - 2 * x - 1, "1.8")   # tribonacci
RHO = _newton(lambda x: x**3 - x - 1, lambda x: 3 * x**2 - 1, "1.3")                    # plast

# nom -> (constante de base, type de bnode : '0' = 0^l ,  '1' = 1 0^(l-1))
BASES = {
    "phi":   (PHI, "0"),
    "trib":  (TAU, "0"),
    "trib1": (TAU, "1"),
    "plast": (RHO, "1"),
}


def _norm(base):
    """
    Accepte :  'phi' | 'trib' | 'trib0' | 'trib1' | 'plast'
               ou directement la constante (Decimal / float / str) : phi, btrib, bplast
               (pour la constante de tribonacci, le regime 0 est pris par defaut ;
                pour le regime 1, ecrire 'trib1').
    """
    if isinstance(base, str) and base in BASES:
        return base
    if base == "trib0":
        return "trib"
    try:
        x = Decimal(str(base))
    except Exception:
        x = None
    if x is not None:
        for nom, cst in (("phi", PHI), ("trib", TAU), ("plast", RHO)):
            if abs(x - cst) < Decimal("1e-12"):
                return nom
    raise ValueError("base inconnue : %r  (valeurs permises : 'phi', 'trib' (ou 'trib0'), "
                     "'trib1', 'plast', ou la constante phi / btrib / bplast)" % (base,))


def bnode_code(l, base):
    """code du bnode de longueur l"""
    base = _norm(base)
    if l < 1:
        raise ValueError("l doit etre >= 1")
    return [0] * l if BASES[base][1] == "0" else [1] + [0] * (l - 1)


_BN = {}          # cache : _BN[base][l] = B(l)


def bnode(l, base):
    """B(l) = BGASINV(1, code du bnode ; base)   (calcule une fois, B(l+1) = round(b*B(l)))"""
    base = _norm(base)
    if l < 1:
        raise ValueError("bnode(l) : l doit etre >= 1 (recu %r)" % (l,))
    t = _BN.setdefault(base, {})
    if not t:
        t[1] = BGASINV2B(1, bnode_code(1, base), BASES[base][0])
    b = BASES[base][0]
    while max(t) < l:
        m = max(t)
        t[m + 1] = custom_round_decimal(b * Decimal(t[m]))    # ajouter un zero
    return t[l]


# ---------------- 1) la tete h_k ----------------
_HEAD = {}


def head(k, base, M=300):
    """
    Retourne (n_k, h_k) :
      n_k = plus petit n tel que, pour tout m >= n, k*|B(m+1) - b*B(m)| < 1/2
            (alors ajouter un zero envoie k*B(m) sur k*B(m+1))
      h_k = code BGAS canonique de k*B(n_k)   (la 'tete')
    """
    base = _norm(base)
    k = int(k)
    if k < 1:
        raise ValueError("k doit etre >= 1")
    if (k, base) in _HEAD:
        return _HEAD[(k, base)]
    b = BASES[base][0]
    half = Decimal("0.5")
    n = 1
    for m in range(M, 0, -1):                      # dernier m ou le critere echoue
        if k * abs(Decimal(bnode(m + 1, base)) - b * Decimal(bnode(m, base))) >= half:
            n = m + 1
            break
    res = (n, BGAS2BB(1, k * bnode(n, base), b))
    _HEAD[(k, base)] = res
    return res


# ---------------- 2) k * B(l) ----------------
def k_bnode(k, l, base):
    """
    k*B(l) par la formule  BGASINV(1, h_k . 0^(l-n_k) ; base)   si l >= n_k,
    sinon (petit l) calcul direct k*B(l).
    """
    base = _norm(base)
    if l < 1:
        raise ValueError("l doit etre >= 1")
    n_k, h = head(k, base)
    if l >= n_k:
        return BGASINV2B(1, h + [0] * (l - n_k), BASES[base][0])
    return k * bnode(l, base)


# ---------------- 3) somme des termes et sequences canoniques ----------------
N_MIN = {"phi": 2, "trib": 3, "trib1": 3, "plast": 10}

def sigma(n, base):
    """sequence sigma_n telle que  S_n = BGASINV(1, sigma_n ; base)"""
    base = _norm(base)
    if n < N_MIN[base]:
        raise ValueError("sigma_n n'est valable que pour n >= %d dans la base %r (recu n=%d)"
                         % (N_MIN[base], base, n))
    if base == "phi":                      # n >= 2
        return [0] + [1] * (n - 2) + [0, 0]
    if base == "trib":                     # n >= 3
        return [0, 0, 0] + [1] * (n - 3) + [0]
    if base == "trib1":                    # n >= 3
        return [1, 0] + [1] * (n - 3) + [0, 0]
    if base == "plast":                    # n >= 10
        return [1, 0, 0, 1] + [0] * 5 + [1] * (n - 9) + [0] * 4
    raise ValueError(base)


def somme_termes(n, base):
    """S_n = B(1)+...+B(n)  via  BGASINV(1, sigma_n ; base)"""
    base = _norm(base)
    return BGASINV2B(1, sigma(n, base), BASES[base][0])


def somme_directe(n, base):
    base = _norm(base)
    return sum(bnode(i, base) for i in range(1, n + 1))


# ---------------- Pell : base 1+sqrt(2) ----------------
BPELL = 1 + Decimal(2).sqrt()


def pell2B(k):
    """P(k) (P0=0,P1=1,P(n)=2P(n-1)+P(n-2)) = BGASINV(0^(k-1)) : BNode 0^l = P(l+1)"""
    if k == 0:
        return 0
    if k == 1:
        return 1
    return BGASINV2B(1, (k - 1) * [0], BPELL)


def pell_lucas2B(k):
    """Q_k (Q0=Q1=2,Q(n)=2Q(n-1)+Q(n-2)) = BGASINV(0 1 0^(k-2)) : motif 01.0^(l-2)"""
    if k <= 1:
        return 2
    return BGASINV2B(1, [0, 1] + (k - 2) * [0], BPELL)


def z2pell2B(k):
    """H_k = Q_k/2 (1,1,3,7,17,41,...) = BGASINV(1 0^(k-2)) : debut de la zone z2 de l'intervalle k-1"""
    if k <= 1:
        return 1
    return BGASINV2B(1, [1] + (k - 2) * [0], BPELL)


def spell2b(n):
    """P(1)+...+P(n) = BGASINV(1^(n-1))   (somme des n premiers nombres de Pell)"""
    if n == 0:
        return 0
    return BGASINV2B(1, [1] * (n - 1), BPELL)


def debuts_zones(k):
    """debuts des quatre zones (deux premiers bits du code) de l'intervalle k = [P(k+1), P(k+2)), k>=3"""
    return {"00": pell2B(k + 1),
            "01": pell_lucas2B(k),
            "10": z2pell2B(k + 1),
            "11": z2pell2B(k + 1) + pell2B(k - 1)}


# ---------------- auto-test ----------------
if __name__ == "__main__":
    print("== somme des termes : BGASINV(sigma_n) == somme directe ==")
    for base in BASES:
        ok = all(somme_termes(n, base) == somme_directe(n, base) for n in range(N_MIN[base], 61))
        print(f"  {base:<6} n={N_MIN[base]}..60 : {ok}")
    print("== k*B(l) : formule tete.0^(l-n_k) == k*B(l) ==")
    for base in BASES:
        ok = all(k_bnode(k, l, base) == k * bnode(l, base) for k in range(1, 31) for l in range(1, 41))
        print(f"  {base:<6} k=1..30, l=1..40 : {ok}")
    print("== exemples ==")
    print("  phi : head(35) =", head(35, "phi"), "| 35*F(10) =", k_bnode(35, 8, "phi"))
    print("  phi : somme_termes(4) =", somme_termes(4, "phi"), "sigma =", sigma(4, "phi"))
    print("== Pell ==")
    Pe = [0, 1]; Qe = [2, 2]; He = [1, 1]
    for _ in range(80):
        Pe.append(2 * Pe[-1] + Pe[-2]); Qe.append(2 * Qe[-1] + Qe[-2]); He.append(2 * He[-1] + He[-2])
    print("  pell2B, pell_lucas2B, z2pell2B (k<=60) :",
          all(pell2B(k) == Pe[k] and pell_lucas2B(k) == Qe[k] and z2pell2B(k) == He[k] for k in range(0, 61)))
    print("  spell2b(n) = P(1)+..+P(n) (n<=60)       :", all(spell2b(n) == sum(Pe[1:n + 1]) for n in range(0, 61)))
    print("  debuts_zones(6) =", debuts_zones(6))
