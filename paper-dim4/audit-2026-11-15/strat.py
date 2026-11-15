def strat(n,mul_c=1.4):
    if n<=1: return []
    s=[[] for _ in range(n+1)]; c=[0.0]*(n+1)
    for i in range(2,n+1):
        bb=1;bc=float('inf')
        for b in range(1,i):
            cost=c[i-b]+c[b]+b*mul_c+(i-b)*1.0
            if cost<bc: bc=cost;bb=b
        s[i]=[bb]+s[i-bb]+s[bb]; c[i]=bc
    return s[n]
def frac(l):
    st=strat(l); idx=0; level=[0]; n=0
    for k in range(l):
        prev=sum(level); target=l-1-k; first=True
        while prev!=target:
            if first: n+=1; first=False
            level.append(st[idx]); prev+=st[idx]; idx+=1
        level.pop()
    return n
# precomputes: the structure before step k is the gluing codomain (k=0) or the codomain of step k-1
for l in (80,82,83,84,85,86,127,129,130,131,168,170,171,172):
    n=frac(l); print(l, n, '%.1f%% of steps; excluding the start structure: %d of %d codomains = %.1f%%'%(100*n/l, n-1, l, 100*(n-1)/l))
