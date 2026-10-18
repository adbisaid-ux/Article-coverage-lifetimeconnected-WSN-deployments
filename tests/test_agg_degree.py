import sys, numpy as np, time
import os; sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))
from model import *
rng=np.random.default_rng(1); bad=0; ds=[]
for it in range(400):
    N=rng.choice([20,30,50,80]); L=rng.choice([60,100,140]); rc=rng.choice([8.,12.,20.])
    sink=np.array([L/2,L/2])
    # cluster around sink to get deep trees
    xy=np.clip(sink+rng.normal(0,L/5,(N,2)),0,L)
    hop,order=bfs_hops(xy,sink,rc*rc)
    d=agg_tree_degree(xy,hop,order,rc*rc); e=exact_agg_degree(xy,sink,rc)
    if (hop>=2).sum()==0: e=0
    ds.append(d)
    if d!=e: bad+=1; print('mismatch',d,e)
print('bad',bad,'dist',np.bincount(np.array(ds,int)))
p=Problem(50,100,100,10,12,20261016,regime="agg"); pop=rng.uniform(0,100,(100,100))
p.evaluate(pop); t=time.perf_counter()
for _ in range(50): F,c,m=p.evaluate(pop)
print('us/eval',1e6*(time.perf_counter()-t)/5000, F[:3], c[:3])
