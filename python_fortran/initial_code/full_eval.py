"""Port of evaluate.cpp adversarial matrix suite for Python bidiag_svd."""
import numpy as np
import sys
import time
from mr3_gk import bidiag_svd

EPS = 2.2204460492503131e-16
SQRT_EPS = 1.4901161193847656e-08
RES_THRESH = 7.0
ORTHO_THRESH = 5.0

class DetRNG:
    def __init__(self, seed=42):
        self.state = np.uint32(seed)
    def next(self):
        self.state = np.uint32(self.state ^ (self.state << np.uint32(13)))
        self.state = np.uint32(self.state ^ (self.state >> np.uint32(17)))
        self.state = np.uint32(self.state ^ (self.state << np.uint32(5)))
        return self.state
    def uniform(self):
        return float(self.next() & np.uint32(0x7FFFFFFF)) / 0x7FFFFFFF

def make(name, n):
    d = np.zeros(n); e = np.zeros(max(n-1,0))
    rng = DetRNG(42)

    # GROUP 1: Original patterns
    if name == 'exponential_graded':
        for i in range(n): d[i] = 10.0**(-(i)*16.0/n)
        for i in range(n-1): e[i] = d[i]*0.5
    elif name == 'glued_repeated':
        block=10
        for b in range(n//block):
            s = b*block
            for i in range(s, min(s+block,n)): d[i] = 1.0+b*1e-14
            for i in range(s, min(s+block-1,n-1)): e[i] = 0.5
            if s>0 and s-1<n-1: e[s-1] = 1e-15
    elif name == 'saw_tooth':
        for i in range(n): d[i] = 1.0 if i%2==0 else 1e-8
        e[:] = 0.5
    elif name == 'stemr_killer':
        for i in range(n): d[i] = 10.0**(-(i)*20.0/n)
        for i in range(n-1): e[i] = np.sqrt(d[i]*d[i+1])
    elif name == 'huge_condition':
        for i in range(n): d[i] = 10.0**(-(i)*15.0/(n-1))
        for i in range(n-1): e[i] = np.sqrt(d[i]*d[i+1])*0.9
    elif name == 'spike':
        d[:] = 0.01; d[n//2] = 100.0; e[:] = 0.01
    elif name == 'wilkinson_like':
        for i in range(n): d[i] = abs(i-n//2)+1.0
        e[:] = 1.0
    elif name == 'two_clusters':
        for i in range(n): d[i] = 1.0 if i<n//2 else 1e-8
        e[:] = 1e-10
        if n//2-1>=0 and n//2-1<n-1: e[n//2-1] = 1e-16
    elif name == 'random_uniform':
        for i in range(n): d[i] = 0.1+rng.uniform()*9.9
        for i in range(n-1): e[i] = 0.01+rng.uniform()*0.99
    elif name == 'diagonal_only':
        for i in range(n): d[i] = 10.0**(-(i)*8.0/(n-1)) if n>1 else 1.0
    elif name == 'constant':
        d[:] = 3.14
    elif name == 'all_equal_nontrivial':
        d[:] = 1.0; e[:] = 1.0
    elif name == 'one_big_cluster':
        for i in range(n): d[i] = 1.0+rng.uniform()*1e-12
        for i in range(n-1): e[i] = rng.uniform()*1e-12
    elif name == 'arithmetic_progression':
        for i in range(n): d[i] = 1.0-(i)/(n-1)*(1.0-1e-8) if n>1 else 1.0
        for i in range(n-1): e[i] = (d[i]+d[i+1])/4.0
    elif name == 'many_near_zero':
        d[:] = 1e-15
        for i in range(min(5,n)): d[i] = 10.0**(-(i))
        e[:] = 1e-16
        for i in range(min(4,n-1)): e[i] = d[i]*0.1
    elif name == 'random_dense_clusters':
        for i in range(n): d[i] = 0.1+rng.uniform()*9.9
        d[:] = np.sort(d)[::-1]
        for i in range(0,n-10,10):
            base = d[i]
            for j in range(i, min(i+10,n)): d[j] = base+rng.uniform()*1e-13
        for i in range(n-1): e[i] = 0.01+rng.uniform()*0.09
    elif name == 'constant_d_graded_e':
        d[:] = 1.0
        for i in range(n-1): e[i] = 10.0**(-(i)*16.0/(n-1)) if n>1 else 1.0
    elif name == 'random_clustered_5':
        centers = [0.5+rng.uniform()*4.5 for _ in range(5)]
        for i in range(n): d[i] = centers[i%5]+rng.uniform()*2e-12-1e-12
        d[:] = np.sort(d)[::-1]
        for i in range(n-1): e[i] = 1e-13+rng.uniform()*(1e-12-1e-13)
    elif name == 'alternating_sign':
        for i in range(n): d[i] = 1.0 if i%2==0 else -1.0
        e[:] = 0.3
    elif name == 'step_function':
        for i in range(n):
            if i<n//3: d[i]=10.0
            elif i<2*n//3: d[i]=1.0
            else: d[i]=0.1
        e[:] = 0.5
    elif name == 'three_clusters':
        for i in range(n):
            if i<n//3: d[i]=1.0
            elif i<2*n//3: d[i]=1e-4
            else: d[i]=1e-8
        e[:] = 1e-10
    elif name == 'random_sparse_e':
        for i in range(n): d[i] = 0.5+rng.uniform()*1.5
        nz = max(1,n//10)
        for k in range(nz):
            idx = int(rng.next() % (n-1))
            e[idx] = 0.1+rng.uniform()*0.9
    # GROUP 2: Marques
    elif name == 'chkbd':
        for i in range(n): d[i] = 10.0**(-(2*i+1))
        for i in range(n-1): e[i] = 10.0**(-(2*(i+1)))
    elif name == 'marques_graded':
        for i in range(n): d[i] = SQRT_EPS**(i/(n-1)) if n>1 else 1.0
        for i in range(n-1): e[i] = np.sqrt(d[i]*d[i+1])*0.3
    # GROUP 3: Grosser-Lang entry-based
    elif name == 'gl_abcon0':
        d[:]=2.0; e[:]=1.0
    elif name == 'gl_abcon1':
        d[:]=1.0; e[:]=2.0
    elif name == 'gl_abcon2':
        for i in range(n): d[i]=float(i+1)
        e[:]=1.0
    elif name == 'gl_abcon3':
        for i in range(n): d[i]=1.0/float(i+1)
        e[:]=1.0
    elif name == 'gl_random':
        for i in range(n): d[i]=rng.uniform()
        for i in range(n-1): e[i]=rng.uniform()
    elif name == 'gl_gradp':
        d[n-1]=1.0
        for i in range(n-2,-1,-1): d[i]=0.5*d[i+1]; e[i]=d[i]
    elif name == 'gl_gradm':
        d[0]=1.0
        for i in range(n-1): d[i+1]=0.5*d[i]; e[i]=d[i+1]
    elif name == 'gl_wilkp':
        L=(n-1)//2; t=float(L)
        for i in range(L): d[i]=t; d[n-1-i]=t; t-=1.0
        d[L]=0.0; e[:]=1.0
    elif name == 'gl_wilkm':
        L=(n-1)//2; t=float(L)
        for i in range(L): d[i]=t; d[n-1-i]=t; t-=1.0
        d[L]=0.0
        for i in range(n-1): e[i]=float(n-1-i)
    elif name == 'gl_wilkw':
        for i in range(n): d[i]=float(i+1)
        for i in range(n-1): e[i]=float(n-1-i)
    elif name == 'gl_wilk2w':
        for i in range(n): d[i]=float(2*(i+1)-1)
        for i in range(n-1): e[i]=float(n-1-i)
    elif name == 'gl_clement':
        for i in range(n-1): e[i]=np.sqrt(float((i+1)*(n-1-i)))
    elif name == 'gl_gro0':
        d[:]=1.0
        for i in range(n-1): e[i]=float(n-1-i)/n
    elif name == 'gl_gro1':
        for i in range(n): d[i]=float(n-i)/n
        for i in range(n-1): e[i]=float(n-1-i)/n
    elif name == 'gl_gro2':
        d[:]=EPS
        for i in range(min(5,n)): d[i]=1.0
        for i in range(n-1): e[i]=np.sqrt(d[i]*d[i+1])
    elif name == 'gl_gro3':
        for i in range(n): d[i]=EPS**(float(i)/(n-1)) if n>1 else 1.0
        for i in range(n-1): e[i]=np.sqrt(d[i]*d[i+1])*0.5
    elif name == 'gl_bwilkp':
        L=(n-1)//2; t=float(L)
        for i in range(L): d[i]=t; d[n-1-i]=t; t-=1.0
        d[L]=0.0
        for i in range(n-1): e[i]=float(min(i+1,n-1-i))
    # GROUP 4: Grosser-Lang spectrum-based (simplified: use diagonal)
    elif name == 'gl_ones':
        d[:]=1.0
    elif name == 'gl_uniform_eps':
        for i in range(n): d[i]=EPS*(1.0+rng.uniform())
    elif name == 'gl_uniform_sqrteps':
        for i in range(n): d[i]=SQRT_EPS*(1.0+rng.uniform())
    elif name == 'gl_uniform_eps_to_1':
        for i in range(n): d[i]=EPS+(1.0-EPS)*rng.uniform()
    elif name == 'gl_geometric_eps_to_1':
        for i in range(n): d[i]=EPS**(1.0-float(i)/(n-1)) if n>1 else 1.0
    elif name == 'gl_random_spectrum':
        for i in range(n): d[i]=rng.uniform()*10.0
        for i in range(n-1): e[i]=rng.uniform()
    elif name == 'gl_clustered_at_1':
        for i in range(n): d[i]=1.0+rng.uniform()*EPS*n
    elif name == 'gl_clustered_at_pm1':
        for i in range(n): d[i]=(1.0 if i%2==0 else -1.0)*(1.0+rng.uniform()*EPS*n)
    elif name == 'gl_clustered_at_eps':
        for i in range(n): d[i]=EPS*(1.0+rng.uniform()*n)
    elif name == 'gl_clustered_at_pmeps':
        for i in range(n): d[i]=(1.0 if i%2==0 else -1.0)*EPS*(1.0+rng.uniform()*n)
    # GROUP 5: Demmel strongly clustered
    elif name == 'demmel_S1pe':
        d[0]=1.0
        for i in range(1,n): d[i]=EPS
        for i in range(n-1): e[i]=np.sqrt(abs(d[i]*d[i+1]))*0.1
    elif name == 'demmel_S1ps':
        d[0]=1.0
        for i in range(1,n): d[i]=SQRT_EPS
        for i in range(n-1): e[i]=np.sqrt(abs(d[i]*d[i+1]))*0.1
    elif name == 'demmel_S2pe':
        for i in range(n): d[i]=1.0 if i<n//2 else EPS
        for i in range(n-1): e[i]=np.sqrt(abs(d[i]*d[i+1]))*0.1
    elif name == 'demmel_S2ps':
        for i in range(n): d[i]=1.0 if i<n//2 else SQRT_EPS
        for i in range(n-1): e[i]=np.sqrt(abs(d[i]*d[i+1]))*0.1
    # GROUP 6: Demmel weakly clustered
    elif name == 'demmel_W1':
        kappa=1.0/EPS
        for i in range(n-1): d[i]=float(i+1)/kappa
        d[n-1]=1.0
        for i in range(n-1): e[i]=np.sqrt(abs(d[i]*d[i+1]))*0.1
    elif name == 'demmel_W2':
        kappa=1.0/SQRT_EPS
        for i in range(n-1): d[i]=float(i+1)/kappa
        d[n-1]=1.0
        for i in range(n-1): e[i]=np.sqrt(abs(d[i]*d[i+1]))*0.1
    elif name == 'demmel_W3':
        for i in range(n): d[i]=1.0/float(n-i)
        for i in range(n-1): e[i]=np.sqrt(abs(d[i]*d[i+1]))*0.1
    # GROUP 7: Demmel geometric/uniform
    elif name in ('demmel_G1','demmel_G1s'):
        kappa=1.0/EPS if name=='demmel_G1' else 1.0/SQRT_EPS
        for i in range(n): d[i]=kappa**(-(i)/(n-1)) if n>1 else 1.0
        for i in range(n-1): e[i]=np.sqrt(d[i]*d[i+1])*0.1
    elif name in ('demmel_U1','demmel_U1s'):
        kappa=1.0/EPS if name=='demmel_U1' else 1.0/SQRT_EPS
        for i in range(n): d[i]=1.0-(1.0-1.0/kappa)*float(i)/(n-1) if n>1 else 1.0
        for i in range(n-1): e[i]=np.sqrt(abs(d[i]*d[i+1]))*0.1
    # GROUP 8: Wilkinson variants
    elif name == 'wilkinson_exact':
        m_w=(n-1)//2
        for i in range(n): d[i]=float(abs(m_w-i))
        e[:]=1.0
    elif name == 'glued_wilkinson':
        m_w=(n-1)//2
        for i in range(n): d[i]=float(abs(m_w-i))+1.0
        for i in range(n-1): e[i]=1.0 if abs(i-n//2)>1 else 1e-10
    elif name == 'glued_wilkinson_tight':
        m_w=(n-1)//2
        for i in range(n): d[i]=float(abs(m_w-i))+1.0
        for i in range(n-1): e[i]=1.0 if abs(i-n//2)>1 else 1e-15
    # GROUP 9: Willems-Lang
    elif name == 'wl_example48':
        for i in range(n): d[i]=0.5+rng.uniform()*0.5
        for i in range(n-1): e[i]=d[i]*1e-8
    # GROUP 10: Parlett-Dhillon
    elif name == 'pd_T0':
        d[:]=1.0
        for i in range(n-1): e[i]=0.5/(i+1)
    elif name == 'pd_T1':
        for i in range(n): d[i]=1.0/(i+1)
        for i in range(n-1): e[i]=np.sqrt(d[i]*d[i+1])
    # GROUP 11: Stress
    elif name == 'zero_diagonal':
        e[:]=1.0
    elif name == 'single_element':
        d[0]=3.14; return d[:1], np.array([])
    elif name == 'near_overflow':
        big=1e150
        for i in range(n): d[i]=big*0.9**i
        for i in range(n-1): e[i]=big*0.9**i*0.5
    elif name == 'near_underflow':
        tiny=1e-150
        for i in range(n): d[i]=tiny*1.1**i
        for i in range(n-1): e[i]=tiny*1.1**i*0.5
    elif name == 'mixed_signs':
        for i in range(n): d[i]=(1.0 if rng.uniform()>0.5 else -1.0)*(0.1+rng.uniform()*9.9)
        for i in range(n-1): e[i]=0.5+rng.uniform()*4.5
    elif name == 'checkerboard':
        for i in range(n): d[i]=1.0 if i%2==0 else EPS
        for i in range(n-1): e[i]=np.sqrt(abs(d[i]*d[i+1]))
    # GROUP 12: Condition variants
    elif name == 'exponential_graded_k4':
        for i in range(n): d[i]=10.0**(-(i)*4.0/n)
        for i in range(n-1): e[i]=d[i]*0.5
    elif name == 'exponential_graded_k8':
        for i in range(n): d[i]=10.0**(-(i)*8.0/n)
        for i in range(n-1): e[i]=d[i]*0.5
    elif name == 'exponential_graded_k12':
        for i in range(n): d[i]=10.0**(-(i)*12.0/n)
        for i in range(n-1): e[i]=d[i]*0.5
    elif name == 'stemr_killer_k5':
        for i in range(n): d[i]=10.0**(-(i)*5.0/n)
        for i in range(n-1): e[i]=np.sqrt(d[i]*d[i+1])
    elif name == 'stemr_killer_k10':
        for i in range(n): d[i]=10.0**(-(i)*10.0/n)
        for i in range(n-1): e[i]=np.sqrt(d[i]*d[i+1])
    elif name == 'huge_condition_k5':
        for i in range(n): d[i]=10.0**(-(i)*5.0/(n-1)) if n>1 else 1.0
        for i in range(n-1): e[i]=np.sqrt(d[i]*d[i+1])*0.9
    elif name == 'huge_condition_k10':
        for i in range(n): d[i]=10.0**(-(i)*10.0/(n-1)) if n>1 else 1.0
        for i in range(n-1): e[i]=np.sqrt(d[i]*d[i+1])*0.9
    elif name == 'demmel_G1_k4':
        kappa=1e4
        for i in range(n): d[i]=kappa**(-(i)/(n-1)) if n>1 else 1.0
        for i in range(n-1): e[i]=np.sqrt(d[i]*d[i+1])*0.1
    elif name == 'demmel_G1_k8':
        kappa=1e8
        for i in range(n): d[i]=kappa**(-(i)/(n-1)) if n>1 else 1.0
        for i in range(n-1): e[i]=np.sqrt(d[i]*d[i+1])*0.1
    elif name == 'demmel_G1_k12':
        kappa=1e12
        for i in range(n): d[i]=kappa**(-(i)/(n-1)) if n>1 else 1.0
        for i in range(n-1): e[i]=np.sqrt(d[i]*d[i+1])*0.1
    elif name == 'demmel_S1pe_k4':
        d[0]=1.0
        for i in range(1,n): d[i]=1e-4
        for i in range(n-1): e[i]=np.sqrt(abs(d[i]*d[i+1]))*0.1
    elif name == 'demmel_S1pe_k8':
        d[0]=1.0
        for i in range(1,n): d[i]=1e-8
        for i in range(n-1): e[i]=np.sqrt(abs(d[i]*d[i+1]))*0.1
    elif name == 'chkbd_4':
        n=4; d=np.zeros(n); e=np.zeros(n-1)
        for i in range(n): d[i]=10.0**(-(2*i+1))
        for i in range(n-1): e[i]=10.0**(-(2*(i+1)))
    elif name == 'chkbd_16':
        n=16; d=np.zeros(n); e=np.zeros(n-1)
        for i in range(n): d[i]=10.0**(-(2*i+1))
        for i in range(n-1): e[i]=10.0**(-(2*(i+1)))
    elif name == 'marques_graded_k4':
        c=1e-4
        for i in range(n): d[i]=c**(float(i)/(n-1)) if n>1 else 1.0
        for i in range(n-1): e[i]=np.sqrt(d[i]*d[i+1])*0.3
    elif name == 'marques_graded_k8':
        c=1e-8
        for i in range(n): d[i]=c**(float(i)/(n-1)) if n>1 else 1.0
        for i in range(n-1): e[i]=np.sqrt(d[i]*d[i+1])*0.3
    else:
        d[:]=1.0
    return d, e

adv_names = [
    'exponential_graded','glued_repeated','saw_tooth','stemr_killer',
    'huge_condition','spike','wilkinson_like','two_clusters',
    'random_uniform','diagonal_only','constant','all_equal_nontrivial',
    'one_big_cluster','arithmetic_progression','many_near_zero',
    'random_dense_clusters','constant_d_graded_e','random_clustered_5',
    'alternating_sign','step_function','three_clusters','random_sparse_e',
    'chkbd','marques_graded',
    'gl_abcon0','gl_abcon1','gl_abcon2','gl_abcon3','gl_random',
    'gl_gradp','gl_gradm','gl_wilkp','gl_wilkm','gl_wilkw','gl_wilk2w',
    'gl_clement','gl_gro0','gl_gro1','gl_gro2','gl_gro3','gl_bwilkp',
    'gl_ones','gl_uniform_eps','gl_uniform_sqrteps','gl_uniform_eps_to_1',
    'gl_geometric_eps_to_1','gl_random_spectrum',
    'gl_clustered_at_1','gl_clustered_at_pm1',
    'gl_clustered_at_eps','gl_clustered_at_pmeps',
    'demmel_S1pe','demmel_S1ps','demmel_S2pe','demmel_S2ps',
    'demmel_W1','demmel_W2','demmel_W3',
    'demmel_G1','demmel_G1s','demmel_U1','demmel_U1s',
    'wilkinson_exact','glued_wilkinson','glued_wilkinson_tight',
    'wl_example48','pd_T0','pd_T1',
    'zero_diagonal','single_element','near_overflow','near_underflow',
    'mixed_signs','checkerboard',
    'exponential_graded_k4','exponential_graded_k8','exponential_graded_k12',
    'stemr_killer_k5','stemr_killer_k10',
    'huge_condition_k5','huge_condition_k10',
    'demmel_G1_k4','demmel_G1_k8','demmel_G1_k12',
    'demmel_S1pe_k4','demmel_S1pe_k8',
    'chkbd_4','chkbd_16','marques_graded_k4','marques_graded_k8',
]

test_sizes = sorted({10, 100, 200})
catastrophic = set()
p=0; t=0

for sz in test_sizes:
    print(f"\n=== Adversarial Matrices (n={sz}) ===")
    for name in adv_names:
        if name in catastrophic:
            t+=1
            print(f"  {name:35s} n={sz:4d}  SKIPPED")
            continue
        try:
            d, e = make(name, sz)
            n = len(d)
            t0 = time.time()
            sigma, U, V, info = bidiag_svd(d, e)
            dt = time.time()-t0
            B = np.zeros((n,n))
            for i in range(n): B[i,i]=d[i]
            for i in range(n-1): B[i,i+1]=e[i]
            bn = max(np.max(np.abs(B)),1e-300)
            res = np.max(np.abs(B-U@np.diag(sigma)@V.T))/(bn*n*EPS)
            ou = np.max(np.abs(U.T@U-np.eye(n)))/(n*EPS)
            ov = np.max(np.abs(V.T@V-np.eye(n)))/(n*EPS)
            ok = res<=RES_THRESH and ou<=ORTHO_THRESH and ov<=ORTHO_THRESH
            if ok: p+=1
            if max(res,ou,ov)>1000: catastrophic.add(name)
            status = 'PASS' if ok else 'FAIL'
            print(f"  {name:35s} n={n:4d}  res={res:8.3f}  ortU={ou:8.3f}  ortV={ov:8.3f}  t={dt:.4f}s  {status}")
        except Exception as ex:
            t+=1
            catastrophic.add(name)
            print(f"  {name:35s} n={sz:4d}  ERROR: {ex}")
        t+=1

print(f"\n{'='*60}")
print(f"TOTAL: {p}/{t} passed ({100*p/max(t,1):.1f}%)")
