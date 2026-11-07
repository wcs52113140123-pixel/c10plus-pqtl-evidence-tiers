"""Vector icons in mm units for the Fig. 1 schematic (axes with 1 data unit = 1 mm)."""
import numpy as np
from matplotlib.patches import FancyBboxPatch, Circle, Ellipse, Polygon, Rectangle, Wedge, PathPatch
from matplotlib.path import Path
def card(ax,x,y,w,h,fc,ec=None,r=2.0,lw=0.0,z=0):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle=f'round,pad=0,rounding_size={r}',fc=fc,ec=ec or 'none',lw=lw,zorder=z))
def arrow(ax,x0,y0,x1,y1,c='#7a7a7a',lw=0.9,hw=1.6,z=3):
    ax.annotate('',xy=(x1,y1),xytext=(x0,y0),arrowprops=dict(arrowstyle=f'-|>,head_length={hw*0.5},head_width={hw*0.32}',lw=lw,color=c,shrinkA=0,shrinkB=0),zorder=z)
def dna(ax,x,y,w,h=4.0,c='#5a5a5a',mark=None,mark_c='#E07B39',mark_x=0.55,lw=0.8):
    t=np.linspace(0,1,200); ph=2*np.pi*2.0*t
    ax.plot(x+w*t,y+h/2*np.sin(ph),color=c,lw=lw,solid_capstyle='round',zorder=4)
    ax.plot(x+w*t,y-h/2*np.sin(ph),color=c,lw=lw,solid_capstyle='round',zorder=4)
    for tt in np.linspace(0.04,0.96,13):
        ax.plot([x+w*tt]*2,[y-h/2*np.sin(2*np.pi*2*tt),y+h/2*np.sin(2*np.pi*2*tt)],color=c,lw=0.4,alpha=0.6,zorder=3)
    if mark:
        xm=x+w*mark_x; ax.plot([xm,xm],[y-h/2-0.6,y+h/2+0.6],color=mark_c,lw=1.6,zorder=5,solid_capstyle='round')
        ax.add_patch(Circle((xm,y+h/2+1.7),0.9,fc=mark_c,ec='none',zorder=5))
def protein(ax,x,y,r=2.3,c='#6E8FB5',epi=None,epi_c='#E07B39',z=4):
    ang=np.linspace(0,2*np.pi,60,endpoint=False); rr=r*(1+0.10*np.sin(3*ang)+0.05*np.cos(5*ang))
    ax.add_patch(Polygon(np.c_[x+rr*np.cos(ang),y+rr*np.sin(ang)],closed=True,fc=c,ec='none',zorder=z))
    # epitope on top
    if epi=='normal': ax.add_patch(Wedge((x,y+r*1.02),0.95,0,180,fc='#2F5D9E',ec='none',zorder=z+1))
    if epi=='altered': ax.add_patch(Polygon([[x-0.95,y+r*0.98],[x,y+r*0.98+1.25],[x+0.95,y+r*0.98]],closed=True,fc=epi_c,ec='none',zorder=z+1))
def antibody(ax,x,y,s=1.0,c='#3a3a3a',lw=1.1,bound=True):
    # Y shape pointing down at (x,y) = tip of the arms
    ax.plot([x,x],[y+2.2*s,y+4.4*s],color=c,lw=lw,solid_capstyle='round',zorder=6)
    ax.plot([x,x-1.3*s],[y+2.2*s,y+0.4*s],color=c,lw=lw,solid_capstyle='round',zorder=6)
    ax.plot([x,x+1.3*s],[y+2.2*s,y+0.4*s],color=c,lw=lw,solid_capstyle='round',zorder=6)
def aptamer(ax,x,y,s=1.0,c='#6A5A8C',lw=1.1):
    t=np.linspace(0,1,80); ax.plot(x+0.9*s*np.sin(6*np.pi*t),y+0.4*s+4.0*s*t,color=c,lw=lw,solid_capstyle='round',zorder=6)
def spectrum(ax,x,y,w=7,h=4.5,c='#2A8C82',lw=1.0,heights=(0.35,0.8,0.5,1.0,0.3,0.6)):
    ax.plot([x,x+w],[y,y],color='#555555',lw=0.6,zorder=5)
    for i,hh in enumerate(heights):
        xx=x+w*(i+0.5)/len(heights); ax.plot([xx,xx],[y,y+h*hh],color=c,lw=lw,solid_capstyle='butt',zorder=6)
def person(ax,x,y,s=1.0,c='#6f6f6f'):
    ax.add_patch(Circle((x,y+3.0*s),0.95*s,fc=c,ec='none',zorder=5))
    ax.add_patch(FancyBboxPatch((x-1.35*s,y),2.7*s,1.9*s,boxstyle=f'round,pad=0,rounding_size={0.9*s}',fc=c,ec='none',zorder=5))
def tube(ax,x,y,s=1.0,c='#E9C46A'):
    ax.add_patch(FancyBboxPatch((x-0.9*s,y),1.8*s,5.0*s,boxstyle=f'round,pad=0,rounding_size={0.85*s}',fc='white',ec='#6f6f6f',lw=0.6,zorder=5))
    ax.add_patch(FancyBboxPatch((x-0.75*s,y+0.15*s),1.5*s,2.8*s,boxstyle=f'round,pad=0,rounding_size={0.7*s}',fc=c,ec='none',zorder=6))
def magnifier(ax,x,y,s=1.0,c='#3a3a3a'):
    ax.add_patch(Circle((x,y),1.6*s,fc='white',ec=c,lw=1.0,zorder=6))
    ax.plot([x+1.15*s,x+2.6*s],[y-1.15*s,y-2.6*s],color=c,lw=1.6,solid_capstyle='round',zorder=6)
def globe(ax,x,y,r=2.4,c='#9A7B2F'):
    ax.add_patch(Circle((x,y),r,fc='none',ec=c,lw=0.9,zorder=6))
    ax.add_patch(Ellipse((x,y),r*1.0,2*r,fc='none',ec=c,lw=0.7,zorder=6))
    ax.plot([x-r,x+r],[y,y],color=c,lw=0.7,zorder=6)
    for dy in (r*0.55,-r*0.55):
        hw=np.sqrt(r*r-dy*dy); ax.plot([x-hw,x+hw],[y+dy]*2,color=c,lw=0.5,zorder=6)
def doc(ax,x,y,s=1.0,c='#4D4D4D'):
    ax.add_patch(Polygon([[x,y],[x+3.6*s,y],[x+3.6*s,y+3.6*s],[x+2.6*s,y+4.6*s],[x,y+4.6*s]],closed=True,fc='white',ec=c,lw=0.8,zorder=6))
    for k in range(3): ax.plot([x+0.7*s,x+2.9*s],[y+(1.0+0.9*k)*s]*2,color=c,lw=0.6,zorder=7)
def badge(ax,x,y,kind,r=1.55):
    col={'up':'#2F5D9E','down':'#E07B39','eq':'#8C8C8C'}[kind]
    ax.add_patch(Circle((x,y),r,fc=col,ec='none',zorder=7))
    if kind=='up': ax.add_patch(Polygon([[x-0.8,y-0.45],[x+0.8,y-0.45],[x,y+0.75]],closed=True,fc='white',ec='none',zorder=8))
    elif kind=='down': ax.add_patch(Polygon([[x-0.8,y+0.45],[x+0.8,y+0.45],[x,y-0.75]],closed=True,fc='white',ec='none',zorder=8))
    else:
        for dy in (0.35,-0.35): ax.plot([x-0.7,x+0.7],[y+dy]*2,color='white',lw=1.1,zorder=8,solid_capstyle='round')
