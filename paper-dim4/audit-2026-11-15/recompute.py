import csv, statistics as st, re, json, sys, math
R='/home/user/work/results/20261114T224715Z-SECENG-2636-clean/session'
def rows(p,name):
    out=[]
    for l in open(f'{R}/pass{p}/raw/{name}.csv'):
        if l.startswith('#') or l.startswith('iteration'): continue
        out.append(l.rstrip('\n').split(','))
    return out
med=st.median
def cell(name, per_sig=False):
    mc=[];ms=[]
    for p in (1,2,3):
        r=rows(p,name)
        if per_sig:
            by={}
            for x in r: by.setdefault(x[3].split(';')[0],[]).append(x)
            mc.append(st.mean(med(int(x[1]) for x in v) for v in by.values())/1e6)
            ms.append(st.mean(med(int(x[2]) for x in v) for v in by.values())/1e6)
        else:
            mc.append(med(int(x[1]) for x in r)/1e6); ms.append(med(int(x[2]) for x in r)/1e6)
    return med(mc), med(ms), mc
def cref(level,op):
    v=[]
    for p in (1,2,3):
        t=open(f'{R}/pass{p}/c-ref-{level}.txt').read()
        m=re.search(op+r'\s*\|.*?median\s+([\d.]+)',t); v.append(float(m.group(1)))
    return med(v), v
C={}
for L in 'I','III','V':
    for op,k in (('keypair','keygen'),('sign','sign'),('verify','verify')):
        C[f'c-{L}-{k}']=cref(L,op)[0]
    for k in 'keygen','sign','verify','verify-compressed':
        C[f'dim2-{L}-{k}']=cell(f'dim2-{L}-{k}')[0]; C[f'ms-dim2-{L}-{k}']=cell(f'dim2-{L}-{k}')[1]
    for k in 'keygen','sign':
        C[f'dim4-{L}-{k}']=cell(f'dim4-{L}-{k}')[0]; C[f'ms-dim4-{L}-{k}']=cell(f'dim4-{L}-{k}')[1]
    for k in 'verify-four','verify-two','verify-two-nocheck':
        ps = L!='I'
        C[f'dim4-{L}-{k}']=cell(f'dim4-{L}-{k}',ps)[0]; C[f'ms-dim4-{L}-{k}']=cell(f'dim4-{L}-{k}',ps)[1]
ref=json.load(open(f'{R}/numbers.json'))['rows']
for k,v in ref.items():
    assert abs(C[k]-v['mc'])<1e-6*max(1,v['mc']) or k.startswith('c-') and abs(C[k]-v['mc'])<1e-3,(k,C[k],v['mc'])
print('numbers.json reproduced from raw for',len(ref),'rows')
f1=lambda x:f'{x:.1f}'
print('\n== Table 3 one-decimal cells')
for k in sorted(C):
    if not k.startswith('ms-'): print(k, f1(C[k]), f'{C[k]:.3f}')
print('\n== ms level I')
for k in sorted(C):
    if k.startswith('ms-') : print(k, f1(C[k]), f'{C[k]:.3f}')
print('\n== ratios (from one-decimal cells | from full precision)')
def r(a,b,nd=1): 
    x=float(f1(C[a]))/float(f1(C[b])); y=C[a]/C[b]; return f'{x:.{nd}f} | {y:.{nd}f}  ({x:.4f} | {y:.4f})'
for L in 'I','III','V':
    print(L,'four/dim2', r(f'dim4-{L}-verify-four',f'dim2-{L}-verify'))
    print(L,'two/dim2 ', r(f'dim4-{L}-verify-two',f'dim2-{L}-verify'))
    print(L,'two/four ', r(f'dim4-{L}-verify-two',f'dim4-{L}-verify-four',2))
    print(L,'sign4/sign2', r(f'dim4-{L}-sign',f'dim2-{L}-sign'))
    print(L,'four/C verify', r(f'dim4-{L}-verify-four',f'c-{L}-verify'))
    print(L,'two/C verify', r(f'dim4-{L}-verify-two',f'c-{L}-verify'))
    print(L,'rust/C verify', r(f'dim2-{L}-verify',f'c-{L}-verify',2),' keygen',r(f'dim2-{L}-keygen',f'c-{L}-keygen',2),' sign',r(f'dim2-{L}-sign',f'c-{L}-sign',2))
    print(L,'compressed/standard verify', r(f'dim2-{L}-verify-compressed',f'dim2-{L}-verify',2))
    print(L,'nocheck/two', r(f'dim4-{L}-verify-two-nocheck',f'dim4-{L}-verify-two',3))
    print(L,'saving pct', 100*(1-C[f'dim4-{L}-verify-two']/C[f'dim4-{L}-verify-four']))
print('per step two points I:', 48.0/174, C['dim4-I-verify-two']/174)
print('\n== profile (per pass medians of each column; shares of total)')
for L in 'I','III','V':
  for mode in 'four','two':
    for p in (1,2,3):
        r_=rows(p,f'dim4-{L}-profile-{mode}')
        cols=list(zip(*[[int(v) for v in x] for x in r_]))
        m=[med(c) for c in cols]
        tot=m[2]; 
        names=['total','decode','s0 basis','s1 corn+kernel','s2 dim2','s3 gluing','s4 pushes','s5 codomain','s6 final','s7 image','hints+pair','goodness']
        print(L,mode,'pass',p,'n',len(r_),' '.join(f'{n}={v/1e6:.3f}({100*v/tot:.1f}%)' for n,v in zip(names,m[2:])))
print('\n== signing spread and loop')
for L in 'I','III','V':
    for p in (1,2,3):
        r_=rows(p,f'dim4-{L}-sign')
        cyc=sorted(int(x[1])/1e6 for x in r_)
        ex=[tuple(int(v) for v in x[3].split(';')) for x in r_]
        s=[e[0] for e in ex]; pt=[e[1] for e in ex]; rs=[e[2] for e in ex]
        q=st.quantiles(cyc,n=4,method='inclusive'); q2=st.quantiles(cyc,n=4,method='exclusive')
        print(L,'pass',p,'n',len(r_),'median %.1f min %.1f max %.1f'%(med(cyc),cyc[0],cyc[-1]),'quartiles incl %.0f-%.0f excl %.0f-%.0f'%(q[0],q[2],q2[0],q2[2]),
          'samples mean %.2f median %s min %d max %d'%(st.mean(s),med(s),min(s),max(s)),'tests mean %.2f'%st.mean(pt),'restarts',sum(rs))
print('\n== heap')
h={}
for l in open(f'{R}/heap.csv'):
    if l.startswith('#') or l.startswith('row'): continue
    a=l.split(','); h.setdefault(a[0],[]).append(int(a[2]))
for k,v in h.items(): print(k,'mean bytes %.0f = %.1f KiB = %.3f MiB'%(st.mean(v),st.mean(v)/1024,st.mean(v)/1048576),'n',len(v), 'min %.0f max %.0f KiB'%(min(v)/1024,max(v)/1024))
print('\n== expectations')
for L,e in (('I',174),('III',264),('V',346)):
    print(L,'e ln2 =',e*math.log(2),'x2 =',2*e*math.log(2))
