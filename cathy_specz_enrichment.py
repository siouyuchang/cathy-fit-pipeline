"""Offline server-side DESI FITS / SDSS CSV -> LS CSV enrichment.
No network access. Build a disk spatial index once, then stream photometry.
See Cathy_specz_server_guide.md for scientific rules and commands.
"""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import json
import sqlite3
from collections import Counter
import numpy as np
import pandas as pd
from astropy.io import fits
from scipy.spatial import cKDTree
from astropy.coordinates import SkyCoord, search_around_sky
import astropy.units as u

POLICY = 'desi_sdss_available_columns_v1'
DESI_MAP = dict(desi_targetid='TARGETID', desi_ra='TARGET_RA', desi_dec='TARGET_DEC',
    desi_z='Z', desi_zerr='ZERR', desi_zwarn='ZWARN', desi_type='SPECTYPE',
    desi_deltachi2='DELTACHI2', desi_maskbits='MASKBITS', desi_morphtype='MORPHTYPE')
SDSS_MAP = dict(sdss_specobjid='specobjid', sdss_ra='ra', sdss_dec='dec',
    sdss_z='z', sdss_zerr='zErr', sdss_zwarning='zwarning', sdss_class='class')
STRINGS = {'desi_targetid','sdss_specobjid','desi_type','desi_morphtype','sdss_class'}
CHORD = 2*np.sin(np.deg2rad(1/3600)/2)

def id_dtypes(columns):
    return {c:'string' for c in columns if 'id' in c.lower() or c in ['release','cluster_name','cluster_id']}

def vectors(ra, dec):
    a,d=np.deg2rad(np.asarray(ra,float)),np.deg2rad(np.asarray(dec,float))
    return np.column_stack((np.cos(d)*np.cos(a),np.cos(d)*np.sin(a),np.sin(d)))

def quality_mask(f, survey):
    p=survey.lower()
    good=np.isfinite(f[p+'_z']) & f[p+'_z'].gt(0) & np.isfinite(f[p+'_zerr']) & f[p+'_zerr'].gt(0) & f[p+'_zerr'].lt(.001)
    good &= np.isfinite(f[p+'_ra']) & f[p+'_ra'].between(0,360,inclusive='left') & np.isfinite(f[p+'_dec']) & f[p+'_dec'].between(-90,90)
    if survey=='DESI':
        good &= f.desi_zwarn.eq(0) & f.desi_type.isin(['GALAXY','QSO']) & f.desi_deltachi2.gt(20) & f.desi_maskbits.eq(0) & f.desi_morphtype.notna() & f.desi_morphtype.ne('PSF')
    else:good &= f.sdss_zwarning.eq(0) | ~f.get('sdss_zwarning_applied',pd.Series(True,index=f.index)).astype(bool)
    return good.fillna(False)

def normalize_chunk(frame):
    for c in frame:
        if c in STRINGS:frame[c]=frame[c].astype('string').str.strip()
        else:frame[c]=pd.to_numeric(frame[c],errors='coerce')
    for c in ['desi_type','desi_morphtype','sdss_class']:
        if c in frame:frame[c]=frame[c].str.upper()
    return frame

def desi_chunks(path, chunksize, hdu):
    # Reopen per block so touched memory-map pages are released between blocks.
    with fits.open(path,memmap=True) as hdus:
        names=hdus[hdu].columns.names
        missing=set(DESI_MAP.values())-set(names)
        if missing:raise ValueError(f'DESI missing columns: {missing}')
        total=len(hdus[hdu].data)
    for start in range(0,total,chunksize):
        with fits.open(path,memmap=True) as hdus:
            block=hdus[hdu].data[start:start+chunksize]
            block=block[block["ZCAT_PRIMARY"] == True]
            block=block[block["ZCAT_PRIMARY"] == True]
            data={}
            for canonical,original in DESI_MAP.items():
                arr=np.array(block[original],copy=True)
                if arr.dtype.kind=='S':arr=np.char.decode(arr,'utf-8')
                if canonical in STRINGS:arr=arr.astype(str)
                elif arr.dtype.byteorder not in ('=','|'):arr=arr.astype(arr.dtype.newbyteorder('='))
                data[canonical]=arr
        yield normalize_chunk(pd.DataFrame(data))

def sdss_chunks(path, chunksize, separator):
    names=list(pd.read_csv(path,sep=separator,nrows=0).columns)
    lower={c.lower():c for c in names}
    mapping={k:lower.get(v.lower()) for k,v in SDSS_MAP.items()}
    missing=[k for k,v in mapping.items() if v is None and k!='sdss_zwarning']
    if missing:raise ValueError(f'SDSS missing columns: {missing}; check --sdss-sep')
    cols=[v for v in mapping.values() if v is not None]
    for f in pd.read_csv(path,sep=separator,usecols=cols,dtype='string',chunksize=chunksize):
        f=f.rename(columns={v:k for k,v in mapping.items() if v is not None})
        f['sdss_zwarning_applied']=mapping['sdss_zwarning'] is not None
        if 'sdss_zwarning' not in f:f['sdss_zwarning']=np.nan
        yield normalize_chunk(f[list(SDSS_MAP)+['sdss_zwarning_applied']])

def connect(path):
    db=sqlite3.connect(str(path))
    db.execute('PRAGMA cache_size=-32768')  # ~32 MiB SQLite cache
    db.execute('PRAGMA temp_store=FILE')
    return db

def build_index(desi, sdss, output, chunksize=50000, hdu=1, separator=','):
    output=Path(output)
    if output.exists():raise FileExistsError(f'Index exists: {output}; reuse with match, or choose a new path')
    # Fail fast on SDSS schema before spending hours indexing DESI.
    probe=sdss_chunks(sdss,1,separator)
    next(probe,None);probe.close()
    output.parent.mkdir(parents=True,exist_ok=True)
    partial=output.with_name(output.name+'.partial')
    if partial.exists():raise FileExistsError(f'Previous incomplete index: {partial}; choose a new path')
    db=connect(partial);stats={}
    db.executescript('CREATE TABLE sources (id INTEGER PRIMARY KEY, survey TEXT, catalog_id TEXT, payload TEXT, UNIQUE(survey,catalog_id)); CREATE VIRTUAL TABLE sky USING rtree(id,xmin,xmax,ymin,ymax,zmin,zmax); CREATE TABLE metadata (key TEXT PRIMARY KEY,value TEXT);')
    try:
        for survey,iterator in [('DESI',desi_chunks(desi,chunksize,hdu)),('SDSS',sdss_chunks(sdss,chunksize,separator))]:
            p=survey.lower();idcol='desi_targetid' if survey=='DESI' else 'sdss_specobjid'
            seen=kept=duplicates=0
            for frame in iterator:
                seen+=len(frame);frame=frame.loc[quality_mask(frame,survey)].copy()
                if frame[idcol].isna().any():raise ValueError(f'{survey}: missing ID')
                xyz=vectors(frame[p+'_ra'],frame[p+'_dec'])
                with db:
                    for record,point in zip(frame.to_dict('records'),xyz):
                        payload=json.dumps({k:None if pd.isna(v) else v for k,v in record.items()},sort_keys=True,allow_nan=False)
                        try:
                            cur=db.execute('INSERT INTO sources(survey,catalog_id,payload) VALUES(?,?,?)',(survey,record[idcol],payload))
                        except sqlite3.IntegrityError:
                            old=db.execute('SELECT payload FROM sources WHERE survey=? AND catalog_id=?',(survey,record[idcol])).fetchone()[0]
                            if old!=payload:raise ValueError(f'{survey}: conflicting repeated ID {record[idcol]}; select one primary observation upstream, do not silently deduplicate')
                            duplicates+=1;continue
                        x,y,z=map(float,point)
                        db.execute('INSERT INTO sky VALUES(?,?,?,?,?,?,?)',(cur.lastrowid,x,x,y,y,z,z));kept+=1
                print(f'{survey}: scanned={seen:,}, indexed={kept:,}, identical_duplicates={duplicates:,}',flush=True)
            stats[survey]={'scanned':seen,'indexed':kept,'identical_duplicates':duplicates}
        metadata={'policy':POLICY,'status':'complete','counts':stats,'sdss_zwarning_rule':'apply == 0 if column exists; otherwise skip as user authorized',
            'utc':datetime.now(timezone.utc).isoformat(),'inputs':[{'path':str(Path(p).resolve()),'size':Path(p).stat().st_size,'mtime_ns':Path(p).stat().st_mtime_ns} for p in [desi,sdss]]}
        with db:db.execute('INSERT INTO metadata VALUES(?,?)',('manifest',json.dumps(metadata)))
    finally:db.close()
    partial.replace(output)
    output.with_suffix(output.suffix+'.json').write_text(json.dumps(metadata,indent=2))
    return metadata

def neighbors(db, point, survey=None):
    # 3D unit-vector bounding cube handles RA wrap and poles without exceptions.
    x,y,z=map(float,point);r=CHORD
    sql='SELECT s.id,s.survey,s.payload FROM sky b JOIN sources s ON s.id=b.id WHERE b.xmax>=? AND b.xmin<=? AND b.ymax>=? AND b.ymin<=? AND b.zmax>=? AND b.zmin<=?'
    args=[x-r,x+r,y-r,y+r,z-r,z+r]
    if survey:sql+=' AND s.survey=?';args.append(survey)
    for key,kind,payload in db.execute(sql,args):
        record=json.loads(payload);p=kind.lower()
        v=vectors([record[p+'_ra']],[record[p+'_dec']])[0]
        if np.linalg.norm(v-point)<CHORD:yield key,kind,record,v

def local_catalogs(db, positions, max_sources):
    found={};queue=[]
    def add(items):
        for key,survey,record,v in items:
            if key not in found:
                found[key]=(survey,record);queue.append((survey,v))
                if len(found)>max_sources:raise MemoryError('Local spectrum limit exceeded; reduce --ls-chunksize. Extremely large connected groups need separate review.')
    valid=np.isfinite(positions.ra)&positions.ra.between(0,360,inclusive='left')&np.isfinite(positions.dec)&positions.dec.between(-90,90)
    for v in vectors(positions.loc[valid,'ra'],positions.loc[valid,'dec']):add(neighbors(db,v))
    # Complete cross-survey connected components, even outside the LS footprint.
    # This prevents chunk/footprint edges from concealing conflicts or ambiguity.
    offset=0
    while offset<len(queue):
        survey,v=queue[offset];offset+=1
        add(neighbors(db,v,'SDSS' if survey=='DESI' else 'DESI'))
    return tuple(pd.DataFrame([r for s,r in found.values() if s==survey],columns=list(mapping)+(['sdss_zwarning_applied'] if survey=='SDSS' else [])) for survey,mapping in [('DESI',DESI_MAP),('SDSS',SDSS_MAP)])

def quality_catalog(frame,survey):
    """Revalidate cuts locally, de-duplicate overlap by exact catalog IDs."""
    f=frame.drop(columns=['query_total'],errors='ignore').copy()
    prefix=survey.lower();idcol='desi_targetid' if survey=='DESI' else 'sdss_specobjid'
    strings={idcol,'desi_type','desi_morphtype','desi_survey','desi_program','sdss_class'}
    for c in f:
        if c not in strings:f[c]=pd.to_numeric(f[c],errors='coerce')
    good=np.isfinite(f[prefix+'_z'])&f[prefix+'_z'].gt(0)&np.isfinite(f[prefix+'_zerr'])&f[prefix+'_zerr'].gt(0)&f[prefix+'_zerr'].lt(.001)
    good &= np.isfinite(f[prefix+'_ra'])&f[prefix+'_ra'].between(0,360,inclusive='left')&np.isfinite(f[prefix+'_dec'])&f[prefix+'_dec'].between(-90,90)
    if survey=='DESI':
        good &= f.desi_zwarn.eq(0)&f.desi_type.isin(['GALAXY','QSO'])&f.desi_deltachi2.gt(20)&f.desi_maskbits.eq(0)&f.desi_morphtype.notna()&f.desi_morphtype.ne('PSF')
    else:good &= f.sdss_zwarning.eq(0) | ~f.get('sdss_zwarning_applied',pd.Series(True,index=f.index)).astype(bool)
    # Index construction already applied these cuts; validate again.
    if not good.fillna(False).all():raise ValueError(f'{survey}: indexed rows fail requested quality cuts')
    unique=f.drop_duplicates()
    if unique[idcol].isna().any() or unique.duplicated(idcol).any():raise ValueError(f'{survey}: conflicting duplicate catalog ID')
    return unique.sort_values(idcol).reset_index(drop=True)


def clean_spec_catalog(desi,sdss):
    """Quality-cut spectra -> 1 arcsec consistency check -> minimum-error spectrum.

    Ambiguous many-to-one cross-survey groups are quarantined rather than
    silently choosing which galaxy is duplicated. No fixed DESI preference,
    except an exact error tie, where DESI is a deterministic tie-breaker.
    """
    d=quality_catalog(desi,'DESI');s=quality_catalog(sdss,'SDSS')
    di=np.array([],dtype=int);si=np.array([],dtype=int);seps=np.array([])
    if len(d) and len(s):
        dc=SkyCoord(d.desi_ra.to_numpy()*u.deg,d.desi_dec.to_numpy()*u.deg)
        sc=SkyCoord(s.sdss_ra.to_numpy()*u.deg,s.sdss_dec.to_numpy()*u.deg)
        nd=cKDTree(vectors(d.desi_ra,d.desi_dec));ns=cKDTree(vectors(s.sdss_ra,s.sdss_dec))
        if nd.count_neighbors(ns,CHORD)>1000000:raise MemoryError('Over one million local spec pairs; reduce LS chunk size or review dense duplicates')
        di,si,sep,_=search_around_sky(dc,sc,1*u.arcsec)
        inside=sep.arcsec<1.;di,si,seps=di[inside],si[inside],sep.arcsec[inside]
    degree_d=np.bincount(di,minlength=len(d));degree_s=np.bincount(si,minlength=len(s))
    ambiguous_d=set(np.flatnonzero(degree_d>1));ambiguous_s=set(np.flatnonzero(degree_s>1))
    # Propagate ambiguity through the entire connected cross-survey group.
    changed=True
    while changed:
        before=(len(ambiguous_d),len(ambiguous_s))
        for a,b in zip(di,si):
            if a in ambiguous_d or b in ambiguous_s:ambiguous_d.add(a);ambiguous_s.add(b)
        changed=before!=(len(ambiguous_d),len(ambiguous_s))
    clean=[];rejected=[];pairs=[]
    for a,b,sep in zip(di,si,seps):
        x,y=d.iloc[a],s.iloc[b];delta=abs(x.desi_z-y.sdss_z)
        reason='ambiguous_spec_pair' if a in ambiguous_d or b in ambiguous_s else ('redshift_conflict' if delta>.005+1e-12 else 'accepted')
        pairs.append({**x.to_dict(),**y.to_dict(),'spec_pair_sep_arcsec':float(sep),'spec_delta_z':delta,'pair_status':reason})
        if reason!='accepted':
            for source,r,ra,dec in [('DESI',x,x.desi_ra,x.desi_dec),('SDSS',y,y.sdss_ra,y.sdss_dec)]:
                rejected.append({**x.to_dict(),**y.to_dict(),'rejected_source':source,'ra':ra,'dec':dec,
                                 'spec_pair_sep_arcsec':float(sep),'spec_delta_z':delta,'reject_reason':reason})
            continue
        source='DESI' if x.desi_zerr<=y.sdss_zerr else 'SDSS'
        best=x if source=='DESI' else y;p=source.lower()
        clean.append({**x.to_dict(),**y.to_dict(),'spec_source':source,
            'spec_object_id':('DESI:'+x.desi_targetid) if source=='DESI' else ('SDSS:'+y.sdss_specobjid),
            'spec_ra':best[p+'_ra'],'spec_dec':best[p+'_dec'],'spec_z':best[p+'_z'],'spec_z_err':best[p+'_zerr'],
            'spec_pair_sep_arcsec':float(sep),'spec_delta_z':delta,'spec_selection':'dual_min_error'})
    for survey,frame,degree in [('DESI',d,degree_d),('SDSS',s,degree_s)]:
        p=survey.lower();idcol='desi_targetid' if survey=='DESI' else 'sdss_specobjid'
        for _,r in frame.loc[degree==0].iterrows():
            clean.append({**r.to_dict(),'spec_source':survey,'spec_object_id':survey+':'+r[idcol],
                'spec_ra':r[p+'_ra'],'spec_dec':r[p+'_dec'],'spec_z':r[p+'_z'],'spec_z_err':r[p+'_zerr'],
                'spec_pair_sep_arcsec':np.nan,'spec_delta_z':np.nan,'spec_selection':'single_survey'})
    clean_columns=list(d.columns)+[c for c in s if c not in d]+['spec_source','spec_object_id','spec_ra','spec_dec','spec_z','spec_z_err','spec_pair_sep_arcsec','spec_delta_z','spec_selection']
    cleaned=pd.DataFrame(clean,columns=clean_columns)
    rejects=pd.DataFrame(rejected,columns=list(d.columns)+list(s.columns)+['rejected_source','ra','dec','spec_pair_sep_arcsec','spec_delta_z','reject_reason']).drop_duplicates()
    pairtable=pd.DataFrame(pairs,columns=list(d.columns)+list(s.columns)+['spec_pair_sep_arcsec','spec_delta_z','pair_status'])
    return cleaned,rejects,pairtable,d,s


def match_clean_to_positions(positions,clean,rejected):
    """Photometric object -> nearest clean spectrum at <1 arcsec, as user code.

    Keep rows with no spectra; flagged conflict/ambiguity cannot be rescued by
    a raw DESI/SDSS field. A nearby rejected spectrum blocks spectroscopic use.
    """
    p=positions.reset_index(drop=True).copy()
    for c in ['ra','dec']:p[c]=pd.to_numeric(p[c],errors='coerce')
    columns=list(clean.columns)+['spec_sep_arcsec','spec_clean_policy','spec_clean_status','spec_quality_ok','spec_reused_within_cluster']
    result=pd.DataFrame(index=p.index,columns=columns)
    for c in columns:
        if c.endswith('id') or c in ['spec_source','spec_selection','desi_type','desi_morphtype','desi_survey','desi_program','sdss_class','spec_clean_policy','spec_clean_status']:
            result[c]=pd.Series(pd.NA,index=p.index,dtype='string')
        elif c=='sdss_zwarning_applied':result[c]=pd.Series(pd.NA,index=p.index,dtype='boolean')
        else:result[c]=np.nan
    result['spec_clean_policy']=POLICY;result['spec_clean_status']='no_match'
    result['spec_quality_ok']=False;result['spec_reused_within_cluster']=False
    valid=np.isfinite(p.ra)&p.ra.between(0,360,inclusive='left')&np.isfinite(p.dec)&p.dec.between(-90,90)
    result.loc[~valid,'spec_clean_status']='invalid_position'
    goodidx=p.index[valid]
    if len(goodidx):
        coords=SkyCoord(p.loc[goodidx,'ra'].to_numpy()*u.deg,p.loc[goodidx,'dec'].to_numpy()*u.deg)
        if len(clean):
            cc=SkyCoord(clean.spec_ra.to_numpy(dtype=float)*u.deg,clean.spec_dec.to_numpy(dtype=float)*u.deg)
            nearest,sep,_=coords.match_to_catalog_sky(cc);m=sep.arcsec<1.;rows=goodidx[m]
            for c in clean:
                values=clean.iloc[nearest[m]][c].astype(result[c].dtype).array
                result.loc[rows,c]=values
            result.loc[rows,'spec_sep_arcsec']=sep.arcsec[m]
            result.loc[rows,'spec_clean_status']='matched';result.loc[rows,'spec_quality_ok']=True
        if len(rejected):
            rc=SkyCoord(rejected.ra.to_numpy(dtype=float)*u.deg,rejected.dec.to_numpy(dtype=float)*u.deg)
            nearest,sep,_=coords.match_to_catalog_sky(rc);m=sep.arcsec<1.;rows=goodidx[m]
            # Store raw pairs for audit, but blank the adopted clean pair.
            for c in [c for c in rejected if c.startswith(('desi_','sdss_')) or c in ['spec_delta_z','spec_pair_sep_arcsec']]:
                values=rejected.iloc[nearest[m]][c].astype(result[c].dtype).array
                result.loc[rows,c]=values
            result.loc[rows,['spec_z','spec_z_err']]=np.nan
            result.loc[rows,'spec_object_id']=pd.NA;result.loc[rows,'spec_source']=pd.NA
            result.loc[rows,'spec_clean_status']=rejected.iloc[nearest[m]].reject_reason.to_numpy()
            result.loc[rows,'spec_quality_ok']=False;result.loc[rows,'spec_sep_arcsec']=sep.arcsec[m]
    matched=result.spec_clean_status.eq('matched')
    keys=pd.DataFrame({'cluster_name':p.cluster_name,'id':result.spec_object_id})
    counts=keys.loc[matched].groupby(['cluster_name','id'])['id'].transform('size')
    result.loc[counts.index,'spec_reused_within_cluster']=counts.gt(1).to_numpy()
    assert result.loc[matched,'spec_z_err'].astype(float).between(0,.001,inclusive='neither').all()
    return result


def enrich_csv(index, source, output, chunksize=2000, max_sources=20000):
    source,output=Path(source),Path(output)
    if output.exists() or source.resolve()==output.resolve():raise FileExistsError('Choose a new output path; inputs are never overwritten')
    header=list(pd.read_csv(source,nrows=0).columns)
    if not {'ra','dec'}.issubset(header):raise ValueError('LS CSV must contain ra, dec in degrees')
    output.parent.mkdir(parents=True,exist_ok=True)
    partial=output.with_name(output.name+'.partial')
    if partial.exists():raise FileExistsError(f'Incomplete previous output: {partial}; choose a new output path')
    if not Path(index).is_file():raise FileNotFoundError(index)
    db=connect(index)
    try:
        row=db.execute("SELECT value FROM metadata WHERE key='manifest'").fetchone()
        manifest=json.loads(row[0]) if row else {}
        if manifest.get('status')!='complete' or manifest.get('policy')!=POLICY:raise ValueError('Incomplete or incompatible index')
        total=0;counts=Counter()
        # Audit stores pairs encountered near this LS file; a pair can recur in
        # different chunks. IDs make repeated audit entries identifiable.
        audit=output.with_name(output.stem+'_pairs.csv');audit_partial=audit.with_name(audit.name+'.partial')
        if audit.exists() or audit_partial.exists():raise FileExistsError(audit)
        audit_written=False
        for block_id,chunk in enumerate(pd.read_csv(source,dtype=id_dtypes(header),chunksize=chunksize)):
            chunk=chunk.reset_index(drop=True)
            positions=pd.DataFrame({'ra':pd.to_numeric(chunk.ra,errors='coerce'),'dec':pd.to_numeric(chunk.dec,errors='coerce'),
                'cluster_name':chunk.get('cluster_name',chunk.get('cluster_id',pd.Series('unspecified',index=chunk.index)))})
            desi,sdss=local_catalogs(db,positions,max_sources)
            clean,rejected,pairs,_,_=clean_spec_catalog(desi,sdss)
            matched=match_clean_to_positions(positions,clean,rejected).drop(columns='spec_reused_within_cluster')
            overlap=set(header)&set(matched.columns)
            if any('input_'+c in header for c in overlap):raise ValueError('Use the original LS CSV, not an already enriched file')
            merged=pd.concat([chunk.rename(columns={c:'input_'+c for c in overlap}),matched],axis=1)
            merged.to_csv(partial,index=False,mode='w' if total==0 else 'a',header=total==0)
            if len(pairs):
                pairs.insert(0,'ls_chunk',block_id)
                pairs.to_csv(audit_partial,index=False,mode='a' if audit_written else 'w',header=not audit_written)
                audit_written=True
            total+=len(chunk);counts.update(matched.spec_clean_status)
            print(f'LS: written={total:,}, status={dict(counts)}',flush=True)
        if total==0:
            d=pd.DataFrame(columns=list(DESI_MAP));s=pd.DataFrame(columns=list(SDSS_MAP)+['sdss_zwarning_applied'])
            clean,rejected,*_=clean_spec_catalog(d,s)
            matches=match_clean_to_positions(pd.DataFrame(columns=['ra','dec','cluster_name']),clean,rejected).drop(columns='spec_reused_within_cluster')
            overlap=set(header)&set(matches.columns)
            pd.DataFrame(columns=[('input_'+c if c in overlap else c) for c in header]+list(matches.columns)).to_csv(partial,index=False)
        if audit_written:audit_partial.replace(audit)
        result={'policy':POLICY,'status':'complete','index':str(Path(index).resolve()),'input':str(source.resolve()),'output':str(output.resolve()),'rows':total,'counts':dict(counts),'ls_chunksize':chunksize,'max_local_sources':max_sources,'index_manifest':manifest}
        output.with_suffix(output.suffix+'.json').write_text(json.dumps(result,indent=2))
        partial.replace(output)
        return result
    finally:db.close()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    sub=parser.add_subparsers(dest='command',required=True)
    build=sub.add_parser('build',help='Stream local server catalogs into a reusable disk index')
    build.add_argument('--desi',required=True);build.add_argument('--sdss',required=True)
    build.add_argument('--index',required=True);build.add_argument('--chunksize',type=int,default=50000)
    build.add_argument('--desi-hdu',type=int,default=1)
    build.add_argument('--sdss-sep',choices=['comma','tab'],default='comma')
    match=sub.add_parser('match',help='Stream LS CSV and append cleaned spectroscopic matches')
    match.add_argument('--index',required=True);match.add_argument('--input',required=True);match.add_argument('--output',required=True)
    match.add_argument('--ls-chunksize',type=int,default=2000);match.add_argument('--max-local-sources',type=int,default=20000)
    args=parser.parse_args()
    if args.command=='build':
        if args.chunksize<1:parser.error('chunksize must be positive')
        build_index(args.desi,args.sdss,args.index,args.chunksize,args.desi_hdu,'\t' if args.sdss_sep=='tab' else ',')
    else:
        if args.ls_chunksize<1 or args.max_local_sources<1:parser.error('chunk and source limits must be positive')
        enrich_csv(args.index,args.input,args.output,args.ls_chunksize,args.max_local_sources)


if __name__=='__main__':main()
