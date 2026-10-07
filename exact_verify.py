"""Independent integer-arithmetic orthogonal rank verifier.

No floating angles, no tolerances. All event rays and intervening chambers
are checked. At events, arbitrary assignments of points on cut lines are
allowed for the weak convention, including points on both lines.
"""
import argparse,json,math,itertools,hashlib
from functools import cmp_to_key
from fractions import Fraction
from pathlib import Path
from collections import Counter

def dot(a,b):return a[0]*b[0]+a[1]*b[1]
def primitive(x,y):
    g=math.gcd(x,y)
    return x//g,y//g
def half(v):return 0 if v[1]>0 or (v[1]==0 and v[0]>0) else 1
def compare(a,b):
    if half(a)!=half(b):return -1 if half(a)<half(b) else 1
    z=a[0]*b[1]-a[1]*b[0]
    return -1 if z>0 else 1 if z<0 else 0
def rays(P):
    out=set()
    for i,a in enumerate(P):
        for b in P[:i]:
            x,y=a[0]-b[0],a[1]-b[1]
            if x==y==0:raise ValueError('Points must be distinct')
            for v in [(x,y),(-x,-y),(-y,x),(y,-x)]:out.add(primitive(*v))
    return sorted(out,key=cmp_to_key(compare))

def top_options(P,v,s):
    values=[dot(x,v) for x in P];level=sorted(values,reverse=True)[s-1]
    fixed={i for i,x in enumerate(values) if x>level};ties=[i for i,x in enumerate(values) if x==level]
    return fixed,ties,s-len(fixed)

def event_overlap_range(P,n,k):
    u=(-n[1],n[0]);Af,At,ac=top_options(P,u,len(P)//2);Bf,Bt,bc=top_options(P,n,2*k)
    # Enumerate only one rank-boundary tie group. For fixed A, extrema over B
    # follow by selecting as few/many points of Bt intersecting A as possible.
    if math.comb(len(At),ac)>1000000:raise ValueError('Tie group too large for this verifier')
    low=len(P);high=-1
    for choice in itertools.combinations(At,ac):
        A=Af|set(choice);base=len(A&Bf);inside=sum(i in A for i in Bt)
        low=min(low,base+max(0,bc-(len(Bt)-inside)))
        high=max(high,base+min(bc,inside))
    # The product of the two fixed-cardinality subset graphs is connected;
    # along an exchange the overlap changes by at most one. Thus every integer
    # between the global minimum and maximum is attainable.
    return low,high

def witness(P,n,k):
    u=(-n[1],n[0]);xs=sorted(dot(x,n) for x in P);ys=sorted(dot(x,u) for x in P)
    b=Fraction(xs[-2*k-1]+xs[-2*k],2);a=Fraction(ys[len(P)//2-1]+ys[len(P)//2],2)
    assert xs[-2*k-1]<b<xs[-2*k] and ys[len(P)//2-1]<a<ys[len(P)//2]
    # Direct inequalities, not rank-set intersection, independently recover cells.
    counts=[0]*4
    for x in P:
        X=dot(x,n)-b;Y=dot(x,u)-a
        assert X and Y
        idx=0 if X>0 and Y>0 else 1 if X<0 and Y>0 else 2 if X<0 and Y<0 else 3
        counts[idx]+=1
    assert counts==[k,len(P)//2-k,len(P)//2-k,k]
    return {'cap_normal':n,'median_normal':u,'cap_offset':str(b),'median_offset':str(a),'cyclic_counts':counts}

def verify(P,k,full=False):
    if any(len(x)!=2 or any(type(z) is not int for z in x) for x in P):
        raise ValueError('This verifier requires exact integer coordinates')
    P=[tuple(map(int,x)) for x in P];N=len(P)
    digest=hashlib.sha256(json.dumps(P,separators=(',',':')).encode()).hexdigest()
    assert N%2==0 and 0<k<=N//4
    rr=rays(P);hist=Counter();first=None
    for a,b in zip(rr,rr[1:]+rr[:1]):
        n=(a[0]+b[0],a[1]+b[1]);assert n!=(0,0)
        u=(-n[1],n[0]);ox=sorted(range(N),key=lambda i:dot(P[i],n),reverse=True)
        oy=sorted(range(N),key=lambda i:dot(P[i],u),reverse=True)
        assert all(dot(P[ox[i]],n)>dot(P[ox[i+1]],n) for i in range(N-1))
        assert all(dot(P[oy[i]],u)>dot(P[oy[i+1]],u) for i in range(N-1))
        q=len(set(ox[:2*k])&set(oy[:N//2]));hist[q-k]+=1
        if q==k and first is None:
            first=witness(P,n,k)
            if not full:return {'status':'partition_found','N':N,'k':k,'points_sha256':digest,'witness':first}
    events=[];critical_min=N;critical_max=0
    for n in rr:
        lo,hi=event_overlap_range(P,n,k)
        critical_min=min(critical_min,lo-k);critical_max=max(critical_max,hi-k)
        if lo<=k<=hi:events.append({'normal':n,'overlap_range':[lo,hi]})
    return {'status':'partition_found' if first else 'weak_partition_only' if events else 'counterexample_certified',
            'N':N,'k':k,'points_sha256':digest,'event_ray_count':len(rr),'chamber_count':sum(hist.values()),
            'defect_histogram':dict(sorted(hist.items())),'witness':first,
            'critical_defect_range':[critical_min,critical_max],
            'weak_solution_event_count':len(events),'first_weak_solution_event':events[0] if events else None}

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('file');ap.add_argument('--k',type=int);ap.add_argument('--full',action='store_true');ap.add_argument('--output');args=ap.parse_args()
    data=json.loads(Path(args.file).read_text());P=data['points'];k=args.k or data['k']
    result=verify(P,k,args.full);s=json.dumps(result,indent=2)+'\n'
    if args.output:Path(args.output).write_text(s)
    print(s)
