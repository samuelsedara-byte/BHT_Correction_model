# Supplementary File S2: BHT correction analysis (Python 3, NumPy, SciPy). Run as a single script.
import datetime as dt
raw = """
1 Idama-04 30-Jan-94 11162 169
2 Idama-04 31-Jan-94 9583 160
3 Idama-04 22-Feb-95 11162 164
4 Idama-04 23-Feb-95 9583 166
5 Inda-08h 15-Apr-12 9560 181
6 Inda-08h 06-Jun-07 9560 185
7 Inda-08h 26-Apr-07 9560 185
8 Inda-08h 18-Feb-07 9560 185
9 Inda-08h 18-Jan-07 9560 185
10 Inda-08h 03-Dec-06 9560 185
11 Inda-08h 22-Jun-06 9560 182
12 Inda-08h 21-Apr-06 9560 180
13 Inda-08h 26-Feb-06 9560 180
14 Rbk-04 29-Sep-85 12950 275
15 Rbk-04 04-Nov-85 12950 274
16 Rbk-04 10-Aug-86 12950 224
17 Rbk-04 17-Aug-86 12950 225
18 Rbk-04 06-Apr-87 12950 230
19 Rbk-04 12-Apr-87 12950 220
20 Rbk-05 30-Oct-85 9850 240
21 Rbk-05 08-Aug-86 9850 175
22 Rbk-05 22-Jul-87 9850 218
23 Rbk-05 20-Aug-88 9850 140
24 Rbk-05 25-Sep-90 9850 190
25 Rbk-05 16-Oct-92 9850 214
26 Rbk-05 23-May-93 9850 240
27 Rbk-08 17-Apr-85 12800 226
28 Rbk-08 14-Jun-86 12800 208
29 Rbk-08 08-Apr-87 12800 226
30 Rbk-08 16-Mar-88 12800 226
31 Rbk-08 05-Aug-90 12800 225
32 Rbk-08 17-May-97 12800 227
33 Rbk-10 12-Aug-86 12800 218
34 Rbk-10 20-Jul-87 12800 215
35 Rbk-10 18-Mar-88 12800 220
36 Rbk-10 08-Aug-90 12800 225
37 Rbk-10 01-May-93 12800 190
38 Rbk-12 08-Aug-95 12800 224
39 Rbk-12 02-Nov-85 10060 190
40 Rbk-12 19-Aug-86 10060 184
41 Rbk-12 21-Jul-87 10060 118
42 Rbk-12 29-Nov-89 10060 218
43 Rbk-12 21-Jan-93 10060 178
44 Rbk-12 16-Aug-95 10060 202
45 Rbk-13 19-Jul-87 9825 200
46 Rbk-13 18-Aug-88 9825 118
47 Rbk-13 23-Jan-93 9825 170
48 Rbk-14 30-Sep-92 14668 250
"""
import csv
rows=[]
for l in raw.strip().splitlines():
    n,w,d,ft,F=l.split()
    rows.append(dict(n=int(n),well=w,date=dt.datetime.strptime(d,'%d-%b-%y').date(),ft=float(ft),F=float(F)))
for r in rows:
    r['z']=r['ft']*0.3048/1000*1.0  # placeholder
    r['z']=r['ft']/3300.0   # file uses 3300 ft/km (3.38242 = 11162/3300)
    r['C']=(r['F']-32)*5/9
    r['C110']=r['C']*1.1
    r['g']=(r['C110']-27)/r['z']
import numpy as np
from scipy import stats, special, optimize
z=np.array([r['z'] for r in rows]); T=np.array([r['C'] for r in rows]); T110=T*1.1
well=np.array([r['well'] for r in rows]); ft=np.array([r['ft'] for r in rows])
Ts=27.0; K=2.4
rw=0.108; alpha=1.0e-6*3600   # m2/h
def x(t): return rw**2/(4*alpha*t)
def E1(t): return special.exp1(x(t))
def R_H(tc,dt): return np.log((tc+dt)/dt)/E1(tc)
def R_LS(tc,dt): return (E1(tc+dt)-E1(dt))/E1(tc)
def R_E(dt,tau): return np.exp(-dt/tau)
def Teq(Tm,R,Tmud=Ts): return Tmud+(Tm-Tmud)/(1-R)
def slope(y): return stats.linregress(z,y).slope
G0=slope(T)
def inv(f,target):
    return optimize.brentq(lambda d: slope(f(d))-target,8,500)
if False:
    print('E1(xc)',E1(4.0),'x_c',x(4.0))
    dtH=inv(lambda d:Teq(T,R_H(4.0,d)),29.0)
    print('Horner implied dt',dtH)
    for tc in (2,4,6): print(tc, slope(Teq(T,R_H(tc,dtH))), slope(Teq(T,R_LS(tc,dtH))))
    print('LS at dtH', slope(Teq(T,R_LS(4.0,dtH))))
    Rneed=1-1/(30.5/G0); tau=-dtH/np.log(Rneed); print('tau for TD=30.5 at dtH',tau, 'Rneed',Rneed)
    print('grid slope Horner/LS: rows dt 6,12,24,48; cols tc 2,4,6')
    for d in (6,12,24,48):
        print(d,[round(slope(Teq(T,R_H(tc,d))),2) for tc in (2,4,6)],[round(slope(Teq(T,R_LS(tc,d))),2) for tc in (2,4,6)])
import json
dtH=inv(lambda d:Teq(T,R_H(4.0,d)),29.0)
Rn=1-1/(30.5/G0); tau=-dtH/np.log(Rneed) if False else -dtH/np.log(1-1/(30.5/G0))
variants={'Measured':T,'110%':T110,'Horner':Teq(T,R_H(4.0,dtH)),'Line source':Teq(T,R_LS(4.0,dtH)),
          'Time-dependent':Teq(T,R_E(dtH,tau)*np.ones_like(T))}
out={'dt_implied':dtH,'tau':tau,'alpha_m2s':1e-6,'rw':rw}
groups={}
for i,(w,f) in enumerate(zip(well,ft)): groups.setdefault((w,f),[]).append(i)
gk=list(groups)
rng=np.random.default_rng(42)
tcrit=stats.t.ppf(0.975,46)
P={}
for n,y in variants.items():
    s=stats.linregress(z,y); res=y-(s.intercept+s.slope*z)
    # cluster bootstrap over well-depth groups
    bs=[]
    for _ in range(5000):
        pick=rng.integers(0,len(gk),len(gk)); idx=np.concatenate([groups[gk[p]] for p in pick])
        if np.ptp(z[idx])<1e-9: continue
        bs.append(stats.linregress(z[idx],y[idx]).slope)
    anch=((y-Ts)*z).sum()/(z**2).sum()
    gper={w:np.mean([(y[i]-Ts)/z[i] for i in range(len(z)) if well[i]==w]) for w in dict.fromkeys(well)}
    P[n]=dict(G=s.slope,se=s.stderr,ci=(s.slope-tcrit*s.stderr,s.slope+tcrit*s.stderr),boot=(np.percentile(bs,2.5),np.percentile(bs,97.5)),
              T0=s.intercept,seT0=s.intercept_stderr,R2=s.rvalue**2,sd=res.std(),maxabs=abs(res).max(),anch=anch,gper=gper,
              hf=K*s.slope,hf_ci=(K*(s.slope-tcrit*s.stderr),K*(s.slope+tcrit*s.stderr)),hf_anch=K*anch)
    print(n,{k:(np.round(v,2) if not isinstance(v,dict) else {a:round(b,2) for a,b in v.items()}) for k,v in P[n].items()})
# residual diagnostics measured
s=stats.linregress(z,T); r=T-(s.intercept+s.slope*z)
print('Shapiro',stats.shapiro(r))
bins=[(2.8,3.1),(3.3,3.5),(3.8,4.0),(4.3,4.5)]
grp=[r[(z>=a)&(z<b)] for a,b in bins]; print([len(g) for g in grp],[round(g.std(ddof=1),1) if len(g)>1 else None for g in grp])
print('Levene (3 bins)',stats.levene(grp[0],grp[1],grp[2]), stats.levene(grp[0],grp[2]))
print('Brown-Forsythe median',stats.levene(grp[0],grp[2],center='median'))
# raw scatter at fixed depth groups
for k,v in groups.items(): 
    t=T[v]; print(k,len(v),round(t.min(),1),round(t.max(),1),round(t.std(ddof=1),1) if len(v)>1 else '')
# per-well tables
for n in variants: print(n,{w:round(v,1) for w,v in P[n]['gper'].items()})
wl=list(dict.fromkeys(well))
for a,b in [('Measured','110%'),('110%','Horner'),('110%','Line source'),('110%','Time-dependent'),('Horner','Time-dependent')]:
    x1=np.array([P[a]['gper'][w] for w in wl]); x2=np.array([P[b]['gper'][w] for w in wl])
    print(a,b,'Wilcoxon',stats.wilcoxon(x1,x2),'mean diff',round((x2-x1).mean(),2),'ttest',stats.ttest_rel(x2,x1))
# t_c and dt sensitivity for pooled (free slope) 
print('tc sens Horner',[round(slope(Teq(T,R_H(tc,dtH))),2) for tc in (2,3,4,5,6)],'LS',[round(slope(Teq(T,R_LS(tc,dtH))),2) for tc in (2,3,4,5,6)],'TD(tau const)',slope(Teq(T,R_E(dtH,tau)*np.ones_like(T))))
print('dt sens Horner(tc4)',[(d,round(slope(Teq(T,R_H(4.0,d))),2)) for d in (12,18,24,30,36,48)])
print('Tm sens: dtH fixed, Tm=27 vs 40',slope(Teq(T,R_H(4.0,dtH),40.0)))
# depth conversion effect
print('3280.84 vs 3300: slope ratio',3300/3280.84)
# max-envelope per group
Tmax=np.array([T[groups[(w,f)]].max() for (w,f) in gk]); zg=np.array([z[groups[(w,f)][0]] for (w,f) in gk])
print('envelope (group max) slope',stats.linregress(zg,Tmax))
json.dump({k:{kk:(vv if not isinstance(vv,np.ndarray) else vv.tolist()) for kk,vv in v.items()} for k,v in P.items()},open('P.json','w'),default=float)
json.dump(out,open('out.json','w'))
import json
dtH=24.25589964; tau=13.637906612
# within-group share of residual variance (measured)
s=stats.linregress(z,T); r=T-(s.intercept+s.slope*z)
groups={}
for i,(w,f) in enumerate(zip(well,ft)): groups.setdefault((w,f),[]).append(i)
ssw=sum(((r[v]-r[v].mean())**2).sum() for v in groups.values()); sst=(r**2).sum()
print('within-group share of residual SS',ssw/sst, 'pure within-group SD pooled', np.sqrt(ssw/(48-11)))
# range of readings within group
for k,v in groups.items(): print(k,round(T[v].max()-T[v].min(),1))
# per-well info
wl=list(dict.fromkeys(well))
for w in wl:
    idx=[i for i in range(48) if well[i]==w]
    print(w,len(idx),sorted(set(np.round(z[idx],2))))
# HF tables
for n in ['Measured','110%','Horner','Line source','Time-dependent']:
    print(n,'pooled HF ranges',round(P[n]['G']*2.2,1),round(P[n]['G']*2.6,1),'1.5-3.5',round(P[n]['G']*1.5,1),round(P[n]['G']*3.5,1),' anch',round(P[n]['anch']*2.2,1),round(P[n]['anch']*2.6,1))
print('perwell HF 110/TD')
for w in wl: 
    g=P['110%']['gper'][w]; h=P['Time-dependent']['gper'][w]
    print(w,round(g*2.4,1),round(g*2.2,1),round(g*2.6,1),'|',round(h*2.4,1),round(h*2.2,1),round(h*2.6,1))
print('mean perwell HF 110',np.mean([P['110%']['gper'][w] for w in wl])*2.4, 'range', min(P['110%']['gper'][w] for w in wl)*2.4,max(P['110%']['gper'][w] for w in wl)*2.4)
# threshold TC for 70 and 60 (free and anchored)
for n in P: print(n,'TC for 70 free',round(70/P[n]['G'],2),'anch TC for 60',round(60/P[n]['anch'],2))
# ratio checks
print('rel to 110 (free)',{n:round((P[n]['hf']/P['110%']['hf']-1)*100,1) for n in P})
print('TD vs Horner per well mean',np.mean([P['Time-dependent']['gper'][w]-P['Horner']['gper'][w] for w in wl]))
# depth 3280.84 effect
print('perwell range 110 g', min(P['110%']['gper'].values()),max(P['110%']['gper'].values()))
