"""Second verifier: rational-slope charts and direct quadrant counts.

Does not import the ray verifier. Critical directions use an independent
dynamic program for boundary assignments, with fixed target marginal sizes.
"""
import argparse,json,hashlib
from fractions import Fraction as Q
from collections import Counter
from pathlib import Path

def cut_groups(vals,size):
    level=sorted(vals,reverse=True)[size-1]
    return [([1] if v>level else [0] if v<level else [0,1]) for v in vals]

def event_possible(P,n,k):
    N=len(P);xx=[x*n[0]+y*n[1] for x,y in P];yy=[-x*n[1]+y*n[0] for x,y in P]
    bx=cut_groups(xx,2*k);ay=cut_groups(yy,N//2)
    # State = selected A size, selected B size, intersection size.
    fixed=[0,0,0];amb=[]
    for a,b in zip(ay,bx):
        if len(a)==len(b)==1:
            fixed[0]+=a[0];fixed[1]+=b[0];fixed[2]+=a[0]*b[0]
        else:amb.append((a,b))
    states={tuple(fixed)}
    for aa,bb in amb:
        nxt=set()
        for A,B,C in states:
            for a in aa:
                for b in bb:
                    z=(A+a,B+b,C+a*b)
                    if z[0]<=N//2 and z[1]<=2*k and z[2]<=k:nxt.add(z)
        states=nxt
    return (N//2,2*k,k) in states

def verify(P,k):
    if any(len(x)!=2 or any(type(z) is not int for z in x) for x in P):
        raise ValueError('This verifier requires exact integer coordinates')
    N=len(P);slopes=set()
    for i,(x,y) in enumerate(P):
        for xx,yy in P[:i]:
            dx=x-xx;dy=y-yy
            if dy:slopes.add(Q(-dx,dy))
            if dx:slopes.add(Q(dy,dx))
    ss=sorted(slopes)
    samples=[ss[0]-1]+[(a+b)/2 for a,b in zip(ss[:-1],ss[1:])]+[ss[-1]+1]
    hist=Counter();examples=[]
    for t in samples:
        for s in [1,-1]:
            nx=s*t.denominator;ny=s*t.numerator
            X=[x*nx+y*ny for x,y in P];Y=[-x*ny+y*nx for x,y in P]
            sx=sorted(X);sy=sorted(Y)
            bx=sx[-2*k-1]+sx[-2*k];ay=sy[N//2-1]+sy[N//2]
            counts=[0,0,0,0]
            for x,y in zip(X,Y):
                Xp=2*x>bx;Yp=2*y>ay
                counts[(0 if Xp else 1) if Yp else (3 if Xp else 2)]+=1
            assert counts[0]+counts[3]==2*k and counts[0]+counts[1]==N//2
            hist[counts[0]-k]+=1
            if counts[0]==k and not examples:examples.append({'normal':[nx,ny],'double_offsets':[bx,ay],'cyclic_counts':counts})
    chart_count=sum(hist.values())
    vertical_is_event=(len({x for x,y in P})<N or len({y for x,y in P})<N)
    if not vertical_is_event:
        # The two slope charts artificially split one genuine chamber at each
        # vertical ray. Merge these chart pieces; there is no geometric event.
        for s in [1,-1]:
            X=[s*y for x,y in P];Y=[-s*x for x,y in P];sx=sorted(X);sy=sorted(Y)
            bx=sx[-2*k-1]+sx[-2*k];ay=sy[N//2-1]+sy[N//2]
            q=sum(2*x>bx and 2*y>ay for x,y in zip(X,Y));hist[q-k]-=1
    event_hits=0
    for t in ss:
        for s in [1,-1]:event_hits+=event_possible(P,(s*t.denominator,s*t.numerator),k)
    if vertical_is_event:
        for s in [1,-1]:event_hits+=event_possible(P,(0,s),k)
    return {'status':'partition_found' if examples else 'weak_partition_only' if event_hits else 'counterexample_certified',
        'N':N,'k':k,'points_sha256':hashlib.sha256(json.dumps(P,separators=(',',':')).encode()).hexdigest(),
        'chart_cell_count':chart_count,'chamber_count':sum(hist.values()),'defect_histogram':dict(sorted(hist.items())),
        'weak_solution_event_count':event_hits,'witness':examples[0] if examples else None}

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('file');ap.add_argument('--output');a=ap.parse_args();d=json.loads(Path(a.file).read_text());r=verify(d['points'],d['k']);s=json.dumps(r,indent=2)+'\n'
    if a.output:Path(a.output).write_text(s)
    print(s)
