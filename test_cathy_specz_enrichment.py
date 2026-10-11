"""Scientific edge cases; run with the same Python environment as the notebook."""
import tempfile
from pathlib import Path
import numpy as np
import pandas as pd
import cathy_specz_enrichment as e


def desi(ra,z=.5,err=.0001,identity=None,**extra):
    identity=str(39627872767312032+int(round(ra-10))) if identity is None else identity
    return dict(desi_targetid=identity,desi_ra=ra,desi_dec=0.,desi_z=z,desi_zerr=err,
                desi_zwarn=0,desi_type='GALAXY',desi_deltachi2=21.,desi_maskbits=0,desi_morphtype='SER',**extra)


def sdss(ra,z=.501,err=.00005,identity=None):
    identity=str(962658753847519232+int(round(ra-10))) if identity is None else identity
    return dict(sdss_specobjid=identity,sdss_ra=ra,sdss_dec=0.,sdss_z=z,sdss_zerr=err,sdss_zwarning=0,sdss_class='GALAXY')


def test_cleaning():
    d=pd.DataFrame([desi(10),desi(20,z=.6),desi(30,z=.7),desi(40,z=.8)])
    s=pd.DataFrame([sdss(10+.2/3600),sdss(20+.2/3600,z=.61),sdss(50,z=.9)])
    clean,rejected,pairs,_,_=e.clean_spec_catalog(d,s)
    both=clean.loc[clean.spec_selection.eq('dual_min_error')].iloc[0]
    assert both.spec_source=='SDSS' and both.spec_z_err==.00005
    assert both.desi_targetid=='39627872767312032' # large ID remains exact string
    assert pairs.pair_status.tolist().count('redshift_conflict')==1
    assert len(clean)==4 and len(rejected)==2
    phot=pd.DataFrame({'cluster_name':['C']*5,'ra':[10,20,30,60,10+.1/3600],'dec':[0.]*5})
    m=e.match_clean_to_positions(phot,clean,rejected)
    assert m.spec_clean_status.tolist()==['matched','redshift_conflict','matched','no_match','matched']
    assert pd.isna(m.spec_z.iloc[1]) and not m.spec_quality_ok.iloc[1]
    assert m.spec_reused_within_cluster.iloc[0] and m.spec_reused_within_cluster.iloc[4]
    assert m.spec_object_id.iloc[0]=='SDSS:962658753847519232'
    # Same-error tie is deterministic; disagreement is still rejected regardless of error.
    tie=e.clean_spec_catalog(pd.DataFrame([desi(10)]),pd.DataFrame([sdss(10,err=.0001)]))[0]
    assert tie.spec_source.iloc[0]=='DESI'
    # Many-to-one counterparts are audited and not silently counted twice.
    amb_s=pd.DataFrame([sdss(10-.1/3600,identity='1'),sdss(10+.1/3600,identity='2')])
    amb=e.clean_spec_catalog(pd.DataFrame([desi(10)]),amb_s)
    assert amb[0].empty and amb[2].pair_status.eq('ambiguous_spec_pair').all()
    # Zero matches and zero catalogs are valid, not exceptions.
    empty=e.clean_spec_catalog(d.iloc[:0],s.iloc[:0])
    assert empty[0].empty
    assert e.match_clean_to_positions(phot,empty[0],empty[1]).spec_clean_status.eq('no_match').all()
    # Exact cut boundaries are strict.
    for column,value in [('desi_zerr',0.),('desi_zerr',.001),('desi_deltachi2',20.),('desi_maskbits',1),('desi_morphtype','PSF'),('desi_type','STAR')]:
        bad=pd.DataFrame([desi(10)]);bad[column]=value
        try:e.quality_catalog(bad,'DESI')
        except ValueError:pass
        else:raise AssertionError((column,value))
    # Photometric match is <1 arcsec; 1.01 arcsec must not match.
    one=e.clean_spec_catalog(pd.DataFrame([desi(10)]),s.iloc[:0])
    far=pd.DataFrame({'cluster_name':['C'],'ra':[10+1.01/3600],'dec':[0.]})
    assert e.match_clean_to_positions(far,one[0],one[1]).spec_clean_status.iloc[0]=='no_match'


def test_missing_survey_string_columns():
    d=pd.DataFrame([desi(10)])
    s=pd.DataFrame([sdss(10)])
    phot=pd.DataFrame({'cluster_name':['C','C'],'ra':[10.,50.],'dec':[0.,0.]})
    for ds,ss,missing,present in [(d,s.iloc[:0],'sdss_specobjid','desi_targetid'),
                                  (d.iloc[:0],s,'desi_targetid','sdss_specobjid')]:
        clean,rejected,*_=e.clean_spec_catalog(ds,ss)
        # Entire opposite-survey ID column is float NaN in the clean table.
        assert clean[missing].isna().all()
        out=e.match_clean_to_positions(phot,clean,rejected)
        assert out[missing].isna().all()
        assert isinstance(out[missing].dtype,pd.StringDtype)
        assert isinstance(out[present].iloc[0],str)
        assert out.spec_clean_status.tolist()==['matched','no_match']
        empty_match=e.match_clean_to_positions(phot.iloc[1:],clean,rejected)
        assert empty_match.spec_clean_status.tolist()==['no_match']


def test_pipeline_adoption():
    import json,os
    os.environ['MPLBACKEND']='Agg';os.environ['MPLCONFIGDIR']='/tmp/cathy_mpl'
    nb=json.loads(Path('Cathy_fit_pipeline.ipynb').read_text());ns={}
    for i in [2,4]:exec(''.join(nb['cells'][i]['source']),ns)
    g=pd.DataFrame({'ra_deg':[1.]*5,'dec_deg':[0.]*5,'photo_z':[.8]*5,'photo_z_err':[.1]*5,
        'desi_z':[.5]*5,'desi_z_err':[.0002]*5,'sdss_z':[.501,.6,.501,.501,.501],
        'sdss_z_err':[.0001]*5,'spec_z':[.501,np.nan,np.nan,np.nan,np.nan],
        'spec_z_err':[.0001,np.nan,np.nan,np.nan,np.nan],
        'spec_clean_policy':[e.POLICY,e.POLICY,e.POLICY,'',''],
        'spec_clean_status':['matched','redshift_conflict','no_match','',''],
        'spec_source':['SDSS',pd.NA,pd.NA,pd.NA,pd.NA],
        'spec_quality_ok':[True,False,False,False,False],
        'desi_quality_ok':[True]*5,'sdss_quality_ok':[True]*5})
    g.loc[4,'sdss_z']=.7
    adopted=ns['adopt_redshifts'](g)
    assert adopted.z_source.tolist()==['SDSS','PHOTO','PHOTO','SDSS','PHOTO']
    assert adopted.adopted_z_err.iloc[0]==.0001 and adopted.spec_blocked.iloc[1]
    assert adopted.spec_blocked.iloc[4]
    # Missing quality columns must not silently bless raw spectroscopy.
    old=g.iloc[3:4].drop(columns=['desi_quality_ok','sdss_quality_ok'])
    assert ns['adopt_redshifts'](old).z_source.iloc[0]=='PHOTO'
    # Scientific boundary: |delta-z| = .005 is accepted, > .005 rejected.
    edge=e.clean_spec_catalog(pd.DataFrame([desi(10,z=.5)]),pd.DataFrame([sdss(10,z=.505)]))
    assert edge[2].pair_status.iloc[0]=='accepted'

def test_offline_streaming():
    from astropy.table import Table
    import json
    with tempfile.TemporaryDirectory() as temp:
        root=Path(temp)
        d=pd.DataFrame([desi(10),desi(20,z=.6),desi(359.99995,identity='99'),desi(40)])
        d.loc[3,'desi_deltachi2']=20 # fail strict quality
        d['desi_targetid']=d.desi_targetid.astype('int64')
        Table.from_pandas(d.rename(columns=e.DESI_MAP).assign(ZCAT_PRIMARY=True)).write(root/'desi.fits')
        # No warning column: explicitly authorized skip, never fabricate zero.
        s=pd.DataFrame([sdss(10+.2/3600),sdss(20+.2/3600,z=.61),sdss(.00005,identity='88')])
        s.rename(columns=e.SDSS_MAP).drop(columns='zwarning').to_csv(root/'sdss.csv',index=False)
        stats=e.build_index(root/'desi.fits',root/'sdss.csv',root/'index.sqlite',chunksize=1)
        assert stats['counts']['DESI']['indexed']==3
        phot=pd.DataFrame({'cluster_name':['C']*6,'ra':[10,20,40,10,0,np.nan],'dec':[0.]*6,
            'object_id':['123456789123456789']*6,'flux':[1,2,3,4,5,6],'desi_z':[.1]*6})
        phot=phot.assign(cluster_ra=10.,cluster_dec=0.,cluster_zphot=.5,ls_photoz=.8,ls_photoz_err=.1,release='9010',brickid='123456',objid=[str(i) for i in range(6)])
        phot.to_csv(root/'ls.csv',index=False);before=(root/'ls.csv').read_bytes()
        for n in [1,4]:e.enrich_csv(root/'index.sqlite',root/'ls.csv',root/f'out{n}.csv',chunksize=n)
        a=pd.read_csv(root/'out1.csv',dtype={'sdss_specobjid':'string','object_id':'string'})
        b=pd.read_csv(root/'out4.csv',dtype={'sdss_specobjid':'string','object_id':'string'})
        pd.testing.assert_frame_equal(a,b)
        assert a.spec_clean_status.tolist()==['matched','redshift_conflict','no_match','matched','matched','invalid_position']
        assert a.sdss_specobjid.iloc[0]=='962658753847519232'
        assert a.object_id.iloc[0]=='123456789123456789'
        assert a.flux.tolist()==[1,2,3,4,5,6] and a.input_desi_z.eq(.1).all()
        assert not a.sdss_zwarning_applied.dropna().astype(bool).any()
        assert (root/'ls.csv').read_bytes()==before
        nb=json.loads(Path('Cathy_fit_pipeline.ipynb').read_text());ns={}
        for i in [2,4]:exec(''.join(nb['cells'][i]['source']),ns)
        clusters,galaxies=ns['normalize_input']({'name':'LS','format':'ls','galaxies':root/'out1.csv','spectra':None})
        assert galaxies.z_source.tolist()==['SDSS','PHOTO','PHOTO','SDSS','SDSS','PHOTO']
        assert galaxies.spec_blocked.iloc[1] and galaxies.adopted_z.iloc[1]==.8
        assert galaxies.sdss_specobjid.iloc[0]=='962658753847519232'
        # Existing zwarning column is honored, including missing/bad values.
        s.rename(columns=e.SDSS_MAP).assign(zwarning=[0,1,np.nan]).to_csv(root/'withwarn.csv',index=False)
        block=next(e.sdss_chunks(root/'withwarn.csv',10,','))
        assert e.quality_mask(block,'SDSS').tolist()==[True,False,False]
        # Follow counterpart chains beyond the initial LS match radius.
        db=e.connect(root/'index.sqlite')
        x,y=e.local_catalogs(db,pd.DataFrame({'ra':[10-.9/3600],'dec':[0.]}),100)
        assert len(x)==1 and len(y)==1
        db.close()
        pd.DataFrame(columns=phot.columns).to_csv(root/'empty.csv',index=False)
        e.enrich_csv(root/'index.sqlite',root/'empty.csv',root/'empty_out.csv')
        assert len(pd.read_csv(root/'empty_out.csv'))==0


if __name__=='__main__':
    test_cleaning();test_missing_survey_string_columns();test_pipeline_adoption();test_offline_streaming()
    print('PASS: offline FITS/CSV streaming, disk index, chunk invariance, cuts, conflict, exact IDs, RA wrap, boundary closure, input preservation, empty input, pipeline adoption')
