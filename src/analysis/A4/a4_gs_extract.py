"""A4: stream GS fastGWA entries from the DataShare ZIP by HTTP Range; keep only GRCh37 cis windows; verify CRC-32."""
import requests, zlib, struct, gzip, hashlib, json, os, sys, time, pandas as pd
from concurrent.futures import ProcessPoolExecutor, as_completed
URL='https://datashare.ed.ac.uk/server/api/core/bitstreams/907d32fe-7513-48e2-a440-643d85384b51/content'
OUT='a4/gs_cis'; os.makedirs(OUT,exist_ok=True)
G=pd.read_csv('hpc/scp-51c9c568cea0821a/a4/gs_entry_gene_windows.tsv',sep='\t',dtype={'chrom':str})
def work(name):
    g=G[G.name==name]; r0=g.iloc[0]; lho,csize,crc_exp=int(r0.lho),int(r0.csize),int(r0.crc)
    wins={}
    for r in g.itertuples(): wins.setdefault(str(r.chrom),[]).append((r.w0_37,r.w1_37))
    outp=f'{OUT}/{name.replace(".fastGWA.gz","")}.cis.tsv.gz'
    for attempt in range(8):
        try:
            rr=requests.get(URL,headers={'Range':f'bytes={lho}-{lho+29}'},timeout=120)
            if rr.status_code!=206: raise IOError(f'HTTP {rr.status_code} on header')
            hdr=rr.content
            sig,ver,flag,meth,mt,md,crc,cs,us,fl,el=struct.unpack('<IHHHHHIIIHH',hdr)
            assert sig==0x04034b50
            start=lho+30+fl+el; end=start+csize-1
            resp=requests.get(URL,headers={'Range':f'bytes={start}-{end}'},stream=True,timeout=300)
            if resp.status_code!=206: raise IOError(f'HTTP {resp.status_code}')
            inf=zlib.decompressobj(-15); crc32=0; sha=hashlib.sha256(); gz=zlib.decompressobj(47)
            buf=b''; nin=nout=0; header=None; nbytes=0
            with gzip.open(outp+'.part','wb',compresslevel=4) as fo:
                def handle(text):
                    nonlocal buf,header,nin,nout
                    buf+=text; lines=buf.split(b'\n'); buf=lines.pop()
                    for ln in lines:
                        if not ln: continue
                        if header is None: header=ln; fo.write(ln+b'\n'); continue
                        nin+=1
                        f=ln.split(None,3)
                        w=wins.get(f[0].decode())
                        if w:
                            p=int(f[2])
                            for a,b in w:
                                if a<=p<=b: fo.write(ln+b'\n'); nout+=1; break
                for chunk in resp.iter_content(1<<20):
                    nbytes+=len(chunk)
                    d=inf.decompress(chunk); crc32=zlib.crc32(d,crc32); sha.update(d)
                    while d:
                        t=gz.decompress(d); handle(t)
                        if gz.eof:
                            d=gz.unused_data; gz=zlib.decompressobj(47)
                        else: d=b''
                d=inf.flush(); crc32=zlib.crc32(d,crc32); sha.update(d)
                if d: handle(gz.decompress(d))
                if buf: handle(b'\n')
            assert nbytes==csize, f'bytes {nbytes}!={csize}'
            ok=(crc32&0xffffffff)==crc_exp
            if not ok: raise ValueError(f'CRC mismatch {crc32} vs {crc_exp}')
            os.replace(outp+'.part',outp)
            osha=hashlib.sha256(open(outp,'rb').read()).hexdigest()
            return dict(name=name,status='ok',attempt=attempt,crc_ok=True,entry_sha256=sha.hexdigest(),rows_in=nin,rows_out=nout,out_sha256=osha,header=header.decode()[:200])
        except Exception as e:
            err=str(e); time.sleep(min(600,60*2**attempt))
    return dict(name=name,status='failed',error=err)
if __name__=='__main__':
    names=sorted(G.name.unique()); done=set()
    if os.path.exists('a4/gs_extract_log.jsonl'):
        for l in open('a4/gs_extract_log.jsonl'):
            j=json.loads(l)
            if j['status']=='ok': done.add(j['name'])
    todo=[n for n in names if n not in done]; print('todo',len(todo),flush=True)
    with ProcessPoolExecutor(int(sys.argv[1]) if len(sys.argv)>1 else 10) as ex:
        fs={ex.submit(work,n):n for n in todo}
        for f in as_completed(fs):
            r=f.result(); open('a4/gs_extract_log.jsonl','a').write(json.dumps(r)+'\n'); print(r['name'],r['status'],r.get('rows_out'),flush=True)
