import requests,zlib,struct,pandas as pd,sys,json
from concurrent.futures import ThreadPoolExecutor
URL='https://datashare.ed.ac.uk/server/api/core/bitstreams/907d32fe-7513-48e2-a440-643d85384b51/content'
g=pd.read_csv('a4/gs_entry_gene_windows.tsv',sep='\t',dtype={'chrom':str})
def census(name):
    r0=g[g.name==name].iloc[0]; lho,cs,crc_exp=int(r0.lho),int(r0.csize),int(r0.crc)
    hdr=requests.get(URL,headers={'Range':f'bytes={lho}-{lho+29}'},timeout=120).content; fl,el=struct.unpack('<HH',hdr[26:30])
    resp=requests.get(URL,headers={'Range':f'bytes={lho+30+fl+el}-{lho+30+fl+el+cs-1}'},stream=True,timeout=300)
    inf=zlib.decompressobj(-15); gz=zlib.decompressobj(47); ch={}; buf=b''; crc=0
    for c in resp.iter_content(1<<20):
        d=inf.decompress(c); crc=zlib.crc32(d,crc)
        while d:
            t=gz.decompress(d); buf+=t; lines=buf.split(b'\n'); buf=lines.pop()
            for ln in lines:
                k=ln.split(b'\t',1)[0].decode(); ch[k]=ch.get(k,0)+1
            if gz.eof: d=gz.unused_data; gz=zlib.decompressobj(47)
            else: d=b''
    return dict(name=name,gene=r0.symbol,crc_ok=(crc&0xffffffff)==crc_exp,chrom_counts=ch)
names=['P00740.fastGWA.gz','P05543.fastGWA.gz','P27918_E9PAQ1.fastGWA.gz','P27918_C9J7V5_E9PAQ1.fastGWA.gz']
with ThreadPoolExecutor(4) as ex: R=list(ex.map(census,names))
json.dump(R,open('posthoc/gs_x_entries_chrom_census.json','w'),indent=1)
for r in R: print(r['name'],r['gene'],r['crc_ok'],sorted(k for k in r['chrom_counts'] if k!='CHR'),sum(v for k,v in r['chrom_counts'].items() if k!='CHR'))
