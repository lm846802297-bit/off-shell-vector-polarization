import sys,json,numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
blue='#246486'; orange='#a34c2d'; grey='#567080'; black='#1d3043'
plt.rcParams.update({'font.size':10,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
def topo():
    fig,axes=plt.subplots(2,2,figsize=(8.5,6.0))
    fig.subplots_adjust(top=.84,bottom=.11,wspace=.08,hspace=.34)
    fig.text(.12,.936,'Where the Higgs pair attaches',size=19,weight='bold',color=black)
    fig.text(.12,.885,'One cubic gauge vertex in every charged tree. Orange dashed lines carry the shifted Higgs momenta.',size=9.4,color=grey)
    patterns=[('Same arm: ordered',[(0,0.29),(0,.64)],True),('Different arms',[(0,.51),(1,.51)],False),
              ('Same arm: seagull',[(0,.51),(0,.51)],True),('Stationary scalar split',[],True)]
    for idx,(ax,(name,hits,bullet)) in enumerate(zip(axes.flat,patterns)):
        ax.set(xlim=(-1.58,1.58),ylim=(-.97,1.55))
        ax.set_aspect('equal');ax.axis("off")
        el=np.array([[0.,1.16],[-1.04,-.62],[1.04,-.62]])
        for p,l in zip(el,(r'$W^+$',r'$W^-$',r'$W^3$')):
            ax.plot([0,p[0]],[0,p[1]],color=blue,lw=2);q=p*1.16;ax.text(*q,l,ha='center',va='center',color=black)
        ax.plot(0,0,'o',color=black,ms=5)
        for j,(i,f) in enumerate(hits):
            c=el[i]*f
            if idx==0:
                q=c+np.array([-.96 if j==0 else .96,.23])
            elif idx==1:
                q=c+np.array([.88 if j==0 else -.53,.32 if j==0 else .69])
            else:
                q=c+np.array([-.90 if j==0 else .90,.20 if j==0 else .33])
            ax.plot([c[0],q[0]],[c[1],q[1]],ls='--',c=orange,lw=1.7);ax.plot(*c,'o',c=blue,ms=4)
            ax.annotate((r'$h(p+zn)$' if j==0 else r'$h(r-zn)$'),q,xytext=(0,6),textcoords='offset points',ha='center',color=orange,fontsize=8)
        if idx==3:
            c=el[0]*.41;q=c+np.array([.62,.04])
            ax.plot([c[0],q[0]],[c[1],q[1]],ls='--',c=grey,lw=1.5);ax.text(.29,.57,r'$h^*$',color=grey)
            for f,lab in ((.31,r'$p+zn$'),(-.39,r'$r-zn$')):
                t=q+np.array([.43,f]);ax.plot([q[0],t[0]],[q[1],t[1]],ls='--',c=orange,lw=1.7);ax.text(*t,lab,ha='left',va='center',color=orange,fontsize=8)
            ax.plot(*c,'o',c=blue,ms=4)
        ax.text(.5,-.04,name,transform=ax.transAxes,ha='center',weight='bold',color=black)
    fig.text(.12,.048,'Representative labelled topologies. Multiplicities are 6, 6, 3 and 3 (left to right, row by row).',fontsize=9,color=grey)
    return fig
def zc(x):return complex(*x)
def norms(rows,neutral):
    out=[]
    for r in rows:
        A=zc(r['ground_amplitude' if neutral else 'A0']); P=zc(r['physical_pole_sum' if neutral else 'poles_only'])
        B=zc(r['boundary']) if neutral else sum(map(zc,r['terms']))
        assert abs(A-P-B)<1e-8*max(1,abs(A),abs(B))
        Q=max(1,abs(A),abs(B))
        out.append((abs(A-P)/Q,abs(A-P-B)/Q))
    return np.array(out)
def verifyplot(neutral,charged):
    f,axs=plt.subplots(1,2,figsize=(8.5,5.15),sharey=True)
    f.subplots_adjust(top=.78,bottom=.19,left=.15,right=.97,wspace=.16)
    f.text(.15,.926,'The missing term in pole reconstruction',size=19,weight='bold',color=black)
    f.text(.15,.865,'Identical error normalization in the neutral and charged calculations',size=10,color=grey)
    for ax,r,label in zip(axs,(neutral,charged),('Neutral Yukawa six-point','Charged non-Abelian five-point')):
        k=np.arange(len(r))+1
        ax.semilogy(k,r[:,0],'o',ms=3.4,c=orange,label='Moving poles only')
        ax.semilogy(k,np.maximum(r[:,1],1e-16),'v',ms=3.4,c=blue,label='Moving poles + boundary')
        ax.set(xlabel='Deformation index',title=label,ylim=(5e-17,5))
        ax.grid(axis='y',which='major',alpha=.16)
    axs[0].set_ylabel(r'$|A(0)-A_{\rm rec}|/\max(1,|A(0)|,|B_\infty|)$')
    axs[0].legend(loc='center left',frameon=False,fontsize=8.5)
    f.text(.15,.07,'20 phase-space points × 2 null directions per panel. Residuals below the plotting floor are shown at $10^{-16}$.',fontsize=8.8,color=grey)
    return f
if __name__=='__main__':
    neutral=norms(json.load(open(sys.argv[1]))['cases'],True)
    charged=norms(json.load(open('results.json'))['rows'],False)
    with PdfPages('Fixed_Root_Boundary_Figures.pdf',metadata={'CreationDate':None,'ModDate':None}) as pdf:
        pdf.savefig(topo());pdf.savefig(verifyplot(neutral,charged))
