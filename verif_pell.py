# verif_pell.py : verifie les formules de reconstruction BGASINV en base de Pell (note_pell.tex)
# usage : python3 verif_pell.py  (necessite bgas_outils.py dans le meme dossier)
import sys; sys.path.insert(0,'.')
import bgas_outils as g
from decimal import Decimal as Dc, getcontext
getcontext().prec=300
al=1+Dc(2).sqrt()
inv=lambda bits:g.BGASINV2B(1,bits,al)
code=lambda x:tuple(g.BGAS2BB(1,x,al))
P=[0,1]; Q=[2,2]; H=[1,1]
for _ in range(300): P.append(2*P[-1]+P[-2]); Q.append(2*Q[-1]+Q[-2]); H.append(2*H[-1]+H[-2])
print("P_n = BGASINV(0^(n-1)), n=1..200      :",all(inv([0]*(n-1))==P[n] for n in range(1,201)))
print("Q_n = BGASINV(0 1 0^(n-2)), n=2..200  :",all(inv([0,1]+[0]*(n-2))==Q[n] for n in range(2,201)))
print("H_n = BGASINV(1 0^(n-2)), n=2..200    :",all(inv([1]+[0]*(n-2))==H[n] for n in range(2,201)))
print("Q_n=2H_n et Q_n=(1+r2)^n+(1-r2)^n (n<=50):",all(Q[n]==2*H[n] and abs(Dc(Q[n])-(al**n+(1-Dc(2).sqrt())**n))<Dc(10)**-100 for n in range(0,51)))
print("sum P_j = (P_{n+1}+P_n-1)/2 (n<=200)  :",all(sum(P[1:n+1])*2==P[n+1]+P[n]-1 for n in range(1,201)))
print("BGASINV(1^(n-1)) = sum_{j<=n} P_j (n<=200):",all(inv([1]*(n-1))==sum(P[1:n+1]) for n in range(1,201)))
print("table sqrt2 : n, P_n, P_{n+1}/P_n-1 :",[(n,P[n],round(P[n+1]/P[n]-1,7)) for n in range(3,9)])
print("limite via BGASINV : BGASINV(0^n)/BGASINV(0^(n-1)) - 1 (n=30) =",float(Dc(inv([0]*30))/Dc(inv([0]*29))-1)," sqrt2 =",float(Dc(2).sqrt()))
# zones, k=6
k=6; lo,hi=P[k+1],P[k+2]
first={}
for x in range(lo,hi): first.setdefault(code(x)[:2],x)
print(f"\nk=6, intervalle [{lo},{hi}) : debuts observes 00,01,10,11 =",[first[(0,0)],first[(0,1)],first[(1,0)],first[(1,1)]])
print("formules : P_{k+1}, Q_k, H_{k+1}, H_{k+1}+P_{k-1} =",[P[k+1],Q[k],H[k+1],H[k+1]+P[k-1]],"| avec P_{k-2} on aurait",H[k+1]+P[k-2])
# fibre de 1^9
k=9; lo,hi=P[k+1],P[k+2]
mem=[x for x in range(lo,hi) if code(x)==(1,)*9]
print("\nfibre de 1^9 :",mem)
print("ecarts :",[v-u for u,v in zip(mem,mem[1:])]," = P_1..P_9 :",[v-u for u,v in zip(mem,mem[1:])]==P[1:10])
