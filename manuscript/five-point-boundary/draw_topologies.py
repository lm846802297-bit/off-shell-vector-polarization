import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm", "font.size": 10, "pdf.fonttype": 42})
fig, axes = plt.subplots(1, 4, figsize=(7.8, 1.88))
fig.subplots_adjust(left=.005, right=.997, bottom=.16, top=.98, wspace=.035)
P = np.array([[0,1.08],[-.98,-.60],[.98,-.60]])
ink = "#101010"
def vector_line(ax, start, stop):
    start, stop = np.asarray(start), np.asarray(stop)
    dr = stop - start
    length = np.linalg.norm(dr)
    normal = np.array([-dr[1], dr[0]])/length
    cycles = max(1, round(length/.16))
    t = np.linspace(0,1,32*cycles+1)
    curve = start + t[:,None]*dr + .025*np.sin(2*np.pi*cycles*t)[:,None]*normal
    ax.plot(*curve.T,c=ink,lw=1.05)

routes = [
 [(0,.30,-.71,.08),(0,.67,.71,.08)],
 [(0,.54,.69,.13),(1,.52,-.48,.65)],
 [(0,.52,-.71,.19),(0,.52,.71,.19)],
 [],
]
for i,ax in enumerate(axes):
 ax.set_aspect("equal")
 ax.set(xlim=(-1.26,1.32),ylim=(-.83,1.32))
 ax.axis("off")
 for arm,(p,label) in enumerate(zip(P, (r"$W^+$",r"$W^-$",r"$W^3$"))):
  ts = sorted(set([0.,1.]+[t for a,t,_,_ in routes[i] if a==arm]+([.38/P[0,1]] if i==3 and arm==0 else [])))
  for t0,t1 in zip(ts,ts[1:]): vector_line(ax,t0*p,t1*p)
  ax.text(*(p*1.13),label,ha="center",va="center")
 ax.plot(0,0,"o",ms=3.4,color="#101010")
 def scalar(a,b,label):
  ax.plot([a[0],b[0]],[a[1],b[1]],ls=(0,(3,2)),c="#101010",lw=1.25)
  ax.text(b[0],b[1]+.06,label,ha="center",va="bottom",fontsize=9)
 for j,(arm,t,dx,dy) in enumerate(routes[i]):
  a=P[arm]*t
  scalar(a,a+[dx,dy],(r"$h_p$",r"$h_r$")[j])
  ax.plot(*a,"o",ms=2.9,color="#101010")
 if i==3:
  a=np.array([0,.38]); b=np.array([.52,.39])
  ax.plot([0,b[0]],[a[1],b[1]],ls=(0,(3,2)),c="#101010",lw=1.25)
  ax.plot(*a,"o",ms=2.9,color="#101010"); ax.plot(*b,"o",ms=2.9,color="#101010")
  ax.text(.25,.48,r"$h^*$",ha="center",fontsize=9)
  scalar(b,[.99,.70],r"$h_p$")
  scalar(b,[.99,.05],r"$h_r$")
 ax.text(.5,-.028,("a","b","c","d")[i],transform=ax.transAxes,ha="center",va="top")
fig.savefig("charged_topologies.pdf",metadata={"CreationDate":None,"ModDate":None})
