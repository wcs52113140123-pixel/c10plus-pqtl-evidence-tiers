"""Refined flat-vector icons for Fig. 1 (axes in mm: 1 data unit = 1 mm)."""
import numpy as np
from matplotlib.patches import FancyBboxPatch, Circle, Ellipse, Polygon, Rectangle, Wedge
INK='#3A3A3A'; AB='#4A4A4A'; NAVY='#2F5D9E'; ORG='#E07B39'; PROT='#8FA9C9'; PROT2='#6E8FB5'
PUR='#6A5A8C'; TEAL='#2A8C82'; GREY='#8C8C8C'; LG='#D9D9D9'; CREAM='#FFF4D6'
def card(ax,x,y,w,h,fc,ec='none',r=2.2,lw=0.0,z=0,ls='-'):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle=f'round,pad=0,rounding_size={r}',fc=fc,ec=ec,lw=lw,zorder=z,ls=ls))
def arrow(ax,x0,y0,x1,y1,c=GREY,lw=0.8,hl=1.1,hw=0.75,z=3):
    ax.annotate('',xy=(x1,y1),xytext=(x0,y0),arrowprops=dict(arrowstyle=f'-|>,head_length={hl*0.45},head_width={hw*0.45}',lw=lw,color=c,shrinkA=0,shrinkB=0),zorder=z)
def protein(ax,x,y,r=1.7,epi='normal',z=5):
    ang=np.linspace(0,2*np.pi,72,endpoint=False)
    rr=r*(1+0.07*np.sin(3*ang+0.6)+0.05*np.cos(5*ang+1.1)+0.03*np.sin(7*ang))
    pts=np.c_[x+rr*np.cos(ang),y+0.88*rr*np.sin(ang)]
    ax.add_patch(Polygon(pts,closed=True,fc=PROT,ec='none',zorder=z))
    sh=Polygon(pts,closed=True,fc=PROT2,ec='none',zorder=z+0.1); ax.add_patch(sh)
    clip=Rectangle((x-2*r,y-2*r),4*r,2*r-0.15*r,transform=ax.transData); sh.set_clip_path(clip)
    ax.add_patch(Polygon(pts,closed=True,fc='none',ec=INK,lw=0.45,zorder=z+0.2))
    top=y+0.88*r*1.0
    k=r/1.7
    if epi=='normal': ax.add_patch(Wedge((x,top-0.08),0.62*k,0,180,fc=NAVY,ec=INK,lw=0.35,zorder=z+0.3))
    elif epi=='altered': ax.add_patch(Polygon([[x-0.62*k,top-0.08],[x,top+0.95*k],[x+0.62*k,top-0.08]],closed=True,fc=ORG,ec=INK,lw=0.35,zorder=z+0.3))
    return top
def bgprotein(ax,x,y,r=0.9,z=4):
    ang=np.linspace(0,2*np.pi,40,endpoint=False); rr=r*(1+0.1*np.sin(3*ang))
    ax.add_patch(Polygon(np.c_[x+1.25*rr*np.cos(ang),y+0.8*rr*np.sin(ang)],closed=True,fc='#DADDE2',ec='none',zorder=z))
def igg(ax,tip,angle_deg,s=1.0,c=AB,z=7,lw_arm=1.45,lw_stem=2.0):
    """IgG drawn binding with its Fab tips at `tip`; the Fc stem points along angle (deg, 90 = up)."""
    a=np.deg2rad(angle_deg); u=np.array([np.cos(a),np.sin(a)]); n=np.array([-u[1],u[0]])
    tip=np.asarray(tip,float); hinge=tip+u*1.35*s
    l=tip+n*0.95*s; r_=tip-n*0.95*s
    for p in (l,r_): ax.plot([hinge[0],p[0]],[hinge[1],p[1]],color=c,lw=lw_arm,solid_capstyle='round',zorder=z)
    stem=hinge+u*2.1*s; ax.plot([hinge[0],stem[0]],[hinge[1],stem[1]],color=c,lw=lw_stem,solid_capstyle='round',zorder=z)
    return stem,u
def oligo_pair(ax,p1,p2,meet,c=NAVY,z=6):
    for p in (p1,p2):
        t=np.linspace(0,1,30); xs=p[0]+(meet[0]-p[0])*t; ys=p[1]+(meet[1]-p[1])*t+0.18*np.sin(np.pi*t)
        ax.plot(xs,ys,color=c,lw=0.7,zorder=z,solid_capstyle='round')
    # short duplex above the meeting point
    x0,y0=meet
    t=np.linspace(0,1,40); ax.plot(x0-0.75+1.5*t,y0+0.25+0.28*np.sin(2*np.pi*1.5*t),color=c,lw=0.7,zorder=z); ax.plot(x0-0.75+1.5*t,y0+0.25-0.28*np.sin(2*np.pi*1.5*t),color=c,lw=0.7,zorder=z)
def dangling(ax,p,dirx,c=NAVY,z=6):
    t=np.linspace(0,1,30); ax.plot(p[0]+dirx*1.0*t,p[1]+0.7*t-1.2*t**2,color=c,lw=0.7,zorder=z,solid_capstyle='round')
def cross(ax,x,y,s=0.55,c=ORG,z=9):
    ax.plot([x-s,x+s],[y-s,y+s],color=c,lw=0.8,zorder=z,solid_capstyle='round'); ax.plot([x-s,x+s],[y+s,y-s],color=c,lw=0.8,zorder=z,solid_capstyle='round')
def aptamer(ax,x,ybot,c=PUR,z=7,s=1.0):
    # stem rising from (x,ybot), then a hairpin loop on top, with a short 5' overhang
    t=np.linspace(0,1,40)
    ax.plot(x-0.35*s+0.15*np.sin(4*np.pi*t)*s, ybot+2.0*s*t, color=c,lw=1.0,zorder=z,solid_capstyle='round')
    ax.plot(x+0.35*s+0.15*np.sin(4*np.pi*t+np.pi)*s, ybot+2.0*s*t, color=c,lw=1.0,zorder=z,solid_capstyle='round')
    for yy in np.linspace(ybot+0.3*s,ybot+1.8*s,4): ax.plot([x-0.35*s,x+0.35*s],[yy,yy],color=c,lw=0.4,zorder=z)
    th=np.linspace(-0.15*np.pi,1.15*np.pi,60); ax.plot(x+0.95*s*np.cos(th)-0.0,ybot+2.0*s+0.75*s+0.85*s*np.sin(th),color=c,lw=1.0,zorder=z)
    ax.plot([x+0.35*s,x+0.9*s],[ybot+2.0*s,ybot+2.45*s],color=c,lw=1.0,zorder=z); ax.plot([x-0.35*s,x-0.9*s],[ybot+2.0*s,ybot+2.45*s],color=c,lw=1.0,zorder=z)
    ax.plot([x-0.35*s,x-0.9*s],[ybot,ybot-0.25*s],color=c,lw=1.0,zorder=z,solid_capstyle='round')
def scissors(ax,x,y,s=1.0,c=GREY,z=7):
    ax.add_patch(Circle((x-0.55*s,y+0.42*s),0.33*s,fc='white',ec=c,lw=0.55,zorder=z)); ax.add_patch(Circle((x-0.55*s,y-0.42*s),0.33*s,fc='white',ec=c,lw=0.55,zorder=z))
    ax.plot([x-0.3*s,x+0.75*s],[y+0.25*s,y-0.25*s],color=c,lw=0.7,zorder=z); ax.plot([x-0.3*s,x+0.75*s],[y-0.25*s,y+0.25*s],color=c,lw=0.7,zorder=z)
def peptide(ax,x,y,w=1.5,h=0.65,fc=PROT,z=7,dot=False):
    ax.add_patch(FancyBboxPatch((x-w/2,y-h/2),w,h,boxstyle=f'round,pad=0,rounding_size={h/2}',fc=fc,ec=INK,lw=0.35,zorder=z))
    if dot: ax.add_patch(Circle((x+0.15,y),0.32,fc=ORG,ec='none',zorder=z+1))
def spectrum(ax,x,y,w=4.4,h=2.2,c=TEAL,z=7,heights=(0.35,0.9,0.55,1.0,0.3)):
    ax.plot([x,x+w],[y,y],color=GREY,lw=0.45,zorder=z)
    for i,hh in enumerate(heights):
        xx=x+w*(i+0.6)/(len(heights)+0.2); ax.plot([xx,xx],[y,y+h*hh],color=c,lw=0.75,zorder=z,solid_capstyle='butt')
def tube(ax,x,y,s=1.0,z=6):
    ax.add_patch(FancyBboxPatch((x-0.95*s,y),1.9*s,6.0*s,boxstyle=f'round,pad=0,rounding_size={0.9*s}',fc='white',ec=INK,lw=0.45,zorder=z))
    ax.add_patch(FancyBboxPatch((x-0.78*s,y+0.18*s),1.56*s,2.5*s,boxstyle=f'round,pad=0,rounding_size={0.75*s}',fc='#CFCFCF',ec='none',zorder=z+0.1))
    ax.add_patch(Rectangle((x-0.78*s,y+2.68*s),1.56*s,2.0*s,fc=CREAM,ec='none',zorder=z+0.1))
    ax.add_patch(Rectangle((x-1.15*s,y+5.6*s),2.3*s,1.0*s,fc='#7d7d7d',ec=INK,lw=0.4,zorder=z+0.2))
def gene_model(ax,x0,x1,y,var_x,var_c,exons,var_in_exon=None,tss_x=None):
    ax.plot([x0,x1],[y,y],color=INK,lw=0.6,zorder=5)
    for k,ex in enumerate(exons):
        ec=var_c if (var_in_exon==k) else INK; lw=0.8 if var_in_exon==k else 0.4
        ax.add_patch(Rectangle((ex-1.15,y-0.95),2.3,1.9,fc=PROT,ec=ec,lw=lw,zorder=6))
    if tss_x is not None:
        ax.plot([tss_x,tss_x],[y,y+1.9],color=INK,lw=0.6,zorder=6); arrow(ax,tss_x,y+1.9,tss_x+1.6,y+1.9,c=INK,lw=0.6,hl=0.8,hw=0.6,z=6)
    ax.plot([var_x,var_x],[y-1.0,y+2.6],color=var_c,lw=1.1,zorder=8,solid_capstyle='round')
    ax.add_patch(Circle((var_x,y+3.25),0.75,fc=var_c,ec='none',zorder=8))
def helix(ax,x,y,w,h=2.4,c='#9A9A9A',z=4):
    t=np.linspace(0,1,200); ph=2*np.pi*2.5*t
    ax.plot(x+w*t,y+h/2*np.sin(ph),color=c,lw=0.6,zorder=z); ax.plot(x+w*t,y-h/2*np.sin(ph),color=c,lw=0.6,zorder=z)
    for tt in np.linspace(0.03,0.97,16): ax.plot([x+w*tt]*2,[y-h/2*np.sin(2*np.pi*2.5*tt),y+h/2*np.sin(2*np.pi*2.5*tt)],color=LG,lw=0.4,zorder=z-0.1)
