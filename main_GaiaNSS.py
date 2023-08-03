# /!\ always do this in the terminal before : export SSL_CERT_FILE=/usr/local/lib/python3.9/site-packages/pip/_vendor/certifi/cacert.pem

def check_mask(t):
    if hasattr(t,'mask'):
        return t.mask
    else:
        return False

def check_len_mask(t):
    if len(t)==0:
        return True
    else:
        return check_mask(t[0])

def check_mask_len(t):
    if check_mask(t):
        return True
    else:
        return len(t)==0

def main():

    import numpy as np
    import warnings
    import astroquery
    from astroquery.gaia import Gaia
    from astroquery.simbad import Simbad
    from astroquery.vizier import Vizier
    from astroquery.esasky import ESASky
    import astropy
    from astropy.io import ascii
    import astropy.units as u
    from astropy.coordinates import SkyCoord
    from astropy.table import Table, vstack, hstack
    from astropy.utils.exceptions import AstropyWarning
    import argparse
        
    print('')
    print('This code will now search for each entry of your target list the Gaia data available '\
          'in the non single star (NSS) two-body orbit catalog. ')
    print('If P and K are available it calculates the corresponding mass. '\
          'The mass of the primary star is retrieved from the gaia catalog as well. ')
    print('It saves everything in a new table **your file**_gaianss.txt.')
    print('')

    # argument parsing
    parser = argparse.ArgumentParser(description="GaiaPMEX-DR3 detection limits map",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("input", help="File with input star parameters",type=str)
    args = parser.parse_args()
    config = vars(args)


    # important parameters to fix
    Gaia.MAIN_GAIA_TABLE = "gaiadr3.gaia_source" # Select Data Release 3
    Gaia_stellar_param = Vizier(catalog="I/355/paramp")
    Gaia_vizier = Vizier(catalog="I/355/gaiadr3")
    TIC_vizier = Vizier(catalog="I/259/tyc2",columns=['HIP'])
    TICsuppl1_vizier = Vizier(catalog="I/259/suppl_1",columns=['HIP'])
    TICsuppl2_vizier = Vizier(catalog="I/259/suppl_2",columns=['HIP'])
    warnings.simplefilter(action="ignore", category=FutureWarning)
    warnings.simplefilter('ignore', category=AstropyWarning)
    warnings.simplefilter('ignore', category=UserWarning)
        
    # file with basic systems data
    filename=config['input']
    data = ascii.read(filename,delimiter='\t',guess=False,data_start=2)
    data.pprint()
    print('')

    # all fields
    Ndata=len(data['name'])
    names_star=data['name'].data
    flaginit=True
    for i in range(Ndata):
        print('{} - {} - '.format(i,names_star[i]),end='',flush=True)
        listname=Table([[data['name'][i]]],names=('name',))
            
        # search simbad
        customSimbad = Simbad()
        customSimbad.add_votable_fields('ra(:;A;ICRS;J2016;2000)', 'dec(:;D;ICRS;J2016;2000)','id(HD)','id(BD)','id(HIP)','id(TYC)','id(Gaia DR3)','sptype')
        customSimbad.remove_votable_fields('coordinates')
        try:
            while True:
                try:
                    simbad_table = customSimbad.query_object(names_star[i])
                    break
                except Exception as e:
                    print(e)
                    print('trying again...')
            if simbad_table is None:
                print('Object not found in SIMBAD'.format(names_star[i]),end='; ',flush=True)
                simbad_table=Table()
            else:
                print('Object found in SIMBAD'.format(names_star[i]),end='; ',flush=True)
                simbad_table=Table([[simbad_table['MAIN_ID'].data[0]],[simbad_table['RA___A_ICRS_J2016_2000'].data[0]],[simbad_table['DEC___D_ICRS_J2016_2000'].data[0]],[simbad_table['ID_HD'].data[0]],[simbad_table['ID_BD'].data[0]],[simbad_table['ID_HIP'].data[0]],[simbad_table['ID_TYC'].data[0]],[simbad_table['ID_Gaia_DR3'].data[0]],[simbad_table['SP_TYPE'].data[0]]],names=('MAIN_ID','RA___A_ICRS_J2016_2000','DEC___D_ICRS_J2016_2000','ID_HD','ID_BD','ID_HIP','ID_TYC','ID_Gaia_DR3','SP_TYPE'))
                if check_len_mask(simbad_table['ID_HIP']):
                    simbad_table['ID_HIP'].mask=False
                    simbad_table['ID_HIP']=['']
                if check_len_mask(simbad_table['ID_HD']):
                    simbad_table['ID_HD'].mask=False
                    simbad_table['ID_HD']=['']
                if check_len_mask(simbad_table['ID_HIP']):
                    simbad_table['ID_BD'].mask=False
                    simbad_table['ID_BD']=['']
                if check_len_mask(simbad_table['ID_HIP']):
                    simbad_table['ID_TYC'].mask=False
                    simbad_table['ID_TYC']=['']
                if check_len_mask(simbad_table['ID_Gaia_DR3']):
                    simbad_table['ID_Gaia_DR3'].mask=False
                    simbad_table['ID_Gaia_DR3']=['']
        except:
            print('error: SIMBAD retrieval failed')
            simbad_table=Table()
            
        # if SIMBAD failed and Gaia DR3 name, search directly the Gaia DR3 to get the coordinates
        if len(simbad_table)==0 and "aia DR3" in names_star[i]:
            print('Using Gaia ID', end='; ',flush=True)
            while True:
                try:
                    gaiatable=Gaia_vizier.query_constraints(Source=names_star[i][9:])[0]
                    break
                except Exception as e:
                    print(e)
                    print('trying again...')
            main_id=str(names_star[i])
            if check_len_mask(gaiatable['TYC2']):
                gaiatable['TYC2'].mask=False
                gaiatable['TYC2']=['']
            elif gaiatable['TYC2'].data[0]!='':
                gaiatable['TYC2']=['TYC {}'.format(gaiatable['TYC2'].data[0])]
                main_id=str('TYC {}'.format(gaiatable['TYC2'].data[0]))
            if check_len_mask(gaiatable['HIP']):
                gaiatable['HIP'].mask=False
                gaiatable['HIP']=['']
            elif gaiatable['HIP'].data[0]!='':
                gaiatable['HIP']=['HIP {}'.format(gaiatable['HIP'].data[0])]
                main_id=str('HIP {}'.format(gaiatable['HIP'].data[0]))
            coord = SkyCoord(ra=gaiatable['RA_ICRS'].data[0], dec=gaiatable['DE_ICRS'].data[0], unit=(u.degree,u.degree), frame='icrs',equinox='J2000',obstime='J2016.0')
            simbad_table=Table([[main_id],[coord.ra.to_string(unit=u.hourangle,sep=':')],[coord.dec.to_string(unit=u.degree,sep=':')],[''],[''],[gaiatable['HIP'].data[0]],[gaiatable['TYC2'].data[0]],[names_star[i]],['']],names=('MAIN_ID','RA___A_ICRS_J2016_2000','DEC___D_ICRS_J2016_2000','ID_HD','ID_BD','ID_HIP','ID_TYC','ID_Gaia_DR3','SP_TYPE'))
        elif len(simbad_table)>0:
            coord = SkyCoord(ra=simbad_table['RA___A_ICRS_J2016_2000'].data[0], dec=simbad_table['DEC___D_ICRS_J2016_2000'].data[0], unit=(u.hourangle, u.deg), frame='icrs',equinox='J2000',obstime='J2016.0')
        else:
            print('Not a Gaia object either.')
            continue
        HIPname=simbad_table['ID_HIP'][0]
        if check_mask_len(HIPname):
            TYCname=simbad_table['ID_TYC'][0]
            if not check_mask_len(TYCname):
                TICliteral=TYCname.split(' ')[-1]
                if len(TICliteral)>0:
                    while True:
                        try:
                            TICcat=TIC_vizier.query_constraints(TYC1=TICliteral.split('-')[0],TYC2=TICliteral.split('-')[1],TYC3=TICliteral.split('-')[2])
                            break
                        except Exception as e:
                            print(e)
                            print('trying again...')
                    if check_len_mask(TICcat):
                        while True:
                            try:
                                TICcat=TICsuppl1_vizier.query_constraints(TYC1=TICliteral.split('-')[0],TYC2=TICliteral.split('-')[1],TYC3=TICliteral.split('-')[2])
                                break
                            except Exception as e:
                                print(e)
                                print('trying again...')
                    if check_len_mask(TICcat):
                        while True:
                            try:
                                TICcat=TICsuppl2_vizier.query_constraints(TYC1=TICliteral.split('-')[0],TYC2=TICliteral.split('-')[1],TYC3=TICliteral.split('-')[2])
                                break
                            except Exception as e:
                                print(e)
                                print('trying again...')
                    if not check_len_mask(TICcat):
                        if TICcat[0]['HIP'].data[0]!='':
                            while True:
                                try:
                                    simbad_table['ID_HIP'][0]='HIP {}'.format(TICcat[0]['HIP'].data[0])
                                    break
                                except Exception as e:
                                    print(e)
                                    print('trying again...')
        
        # now search Gaia catalog
        radius = u.Quantity(5, u.arcsec)
        try:
            queryTable_main = Gaia.query_object(coordinate=coord, radius=radius)
            if len(queryTable_main)==0:
                print('Object not found in Gaia EDR3')
                continue
            else:
                print('Object found in Gaia EDR3',end='; ',flush=True)
        except:
            print('error: Gaia EDR3 retrieval failed')
            continue
        Nfound=len(queryTable_main['source_id'].data)
        if Nfound>1:
            if len(simbad_table['ID_Gaia_DR3'])>0:
                print('{} entries : selecting correct Gaia DR3 ID'.format(Nfound),end='; ',flush=True)
                source_ids=queryTable_main['source_id'].data
                for ii,ids in enumerate(source_ids):
                    if '{}'.format(ids) in simbad_table['ID_Gaia_DR3'].data[0]:
                        queryTable_main=queryTable_main[ii:ii+1]
                        break
            else:
                print('{} entries : selecting brightest'.format(Nfound),end='; ',flush=True)
                ibrightest=np.argsort(queryTable_main['phot_g_mean_mag'].data)[0]
                queryTable_main=queryTable_main[ibrightest:ibrightest+1]
        source_id=queryTable_main['source_id'].data[0]
        
        # now search the NSS two body orbit table
        try:
            job = Gaia.launch_job("SELECT  "
                                  "  source_id, nss_solution_type, ra, dec, parallax, parallax_error, pmra, pmra_error, pmdec, pmdec_error,"
                                  "  period, period_error, t_periastron, t_periastron_error, eccentricity, eccentricity_error, "
                                  "  semi_amplitude_primary, semi_amplitude_primary_error, semi_amplitude_secondary, semi_amplitude_secondary_error,"
                                  "  mass_ratio, mass_ratio_error, inclination, inclination_error, arg_periastron, arg_periastron_error,"
                                  "  astrometric_n_good_obs_al, rv_n_good_obs_primary, rv_n_good_obs_secondary, goodness_of_fit "
                                  "FROM gaiadr3.nss_two_body_orbit "
                                  "WHERE source_id={}".format(source_id))
            gaiansstable = job.get_results()
            gaiansstable['source_id'].name='source_id_2body'
            gaiansstable['nss_solution_type'].name='nss_solution_type_2body'
            gaiansstable['ra'].name='ra_2body'
            gaiansstable['dec'].name='dec_2body'
            gaiansstable['parallax'].name='parallax_2body'
            gaiansstable['parallax_error'].name='parallax_error_2body'
            gaiansstable['pmra'].name='pmra_2body'
            gaiansstable['pmra_error'].name='pmra_error_2body'
            gaiansstable['pmdec'].name='pmdec_2body'
            gaiansstable['pmdec_error'].name='pmdec_error_2body'
            gaiansstable['astrometric_n_good_obs_al'].name='astrometric_n_good_obs_al_2body'
            gaiansstable['goodness_of_fit'].name='goodness_of_fit_2body'
            if len(gaiansstable)>0:
                print('found in NSS 2-body table',end='; ',flush=True)
                flagNSS2body=True
            else:
                print('not found in NSS 2-body table',end='; ',flush=True)
                gaiansstable.add_row(np.ma.array(np.arange(len(gaiansstable.colnames)),mask=True))
                flagNSS2body=False
        except:
            print('error: Gaia NSS 2-body table retrieval failed')
            continue
            
        # now search the NSS acceleration table
        try:
            job = Gaia.launch_job("SELECT  "
                                  "  source_id, nss_solution_type, ra, dec, parallax, parallax_error, pmra, pmra_error, pmdec, pmdec_error,"
                                  "  accel_ra, accel_ra_error, accel_dec, accel_dec_error, deriv_accel_ra, deriv_accel_ra_error, "
                                  "  deriv_accel_dec, deriv_accel_dec_error, "
                                  "  astrometric_n_good_obs_al, goodness_of_fit "
                                  "FROM gaiadr3.nss_acceleration_astro "
                                  "WHERE source_id={}".format(source_id))
            gaianssaccelerationtable = job.get_results()
            gaianssaccelerationtable['source_id'].name='source_id_acc'
            gaianssaccelerationtable['nss_solution_type'].name='nss_solution_type_acc'
            gaianssaccelerationtable['ra'].name='ra_acc'
            gaianssaccelerationtable['dec'].name='dec_acc'
            gaianssaccelerationtable['parallax'].name='parallax_acc'
            gaianssaccelerationtable['parallax_error'].name='parallax_error_acc'
            gaianssaccelerationtable['pmra'].name='pmra_acc'
            gaianssaccelerationtable['pmra_error'].name='pmra_error_acc'
            gaianssaccelerationtable['pmdec'].name='pmdec_acc'
            gaianssaccelerationtable['pmdec_error'].name='pmdec_error_acc'
            gaianssaccelerationtable['astrometric_n_good_obs_al'].name='astrometric_n_good_obs_al_acc'
            gaianssaccelerationtable['goodness_of_fit'].name='goodness_of_fit_acc'
            if len(gaianssaccelerationtable)>0:
                print('found in NSS acceleration table',end='; ',flush=True)
                flagNSSacc=True
            else:
                print('not found in NSS acceleration table',end='; ',flush=True)
                gaianssaccelerationtable.add_row(np.ma.array(np.arange(len(gaianssaccelerationtable.colnames)),mask=True))
                flagNSSacc=False
        except:
            print('error: Gaia NSS acceleration table retrieval failed')
            continue
        
        
        # now retrieve stellar mass from Gaia Collaboration, 2022
        # https://vizier.u-strasbg.fr/viz-bin/VizieR-3?-source=+I%2F355%2Fparamp&-from=nav&-nav=cat%3AI%2F355%26tab%3A%7BI%2F355%2Fparamp%7D%26key%3Asource%3DI%2F355%2Fparamsup%26HTTPPRM%3A%26
        try:
            star_param_table=Gaia_stellar_param.query_constraints(Source=source_id)[0]
            if len(star_param_table)==0:
                print('not found in Gaia star params')
                star_param_table.add_row(np.ma.array(np.arange(len(star_param_table.colnames)),mask=True))
                flagmass=False
            else:
                print('found in Gaia star params')
                flagmass=True
        except:
            print('error: Gaia star params retrieval failed')
            continue
            
        # add a flag when found in the NSS
        NSSflagtable=Table([[queryTable_main['non_single_star'].value[0]]],names=('NSS_flag',))
           
        # finally calculates the mass of the companion
        dataplanet=Table([[],[]], names=('mass_companion_mjup','mass_companion_mjup_error'))
        if flagNSS2body and flagmass:
            for j in range(len(gaiansstable)):
                if not np.ma.isarray(gaiansstable['semi_amplitude_primary'][j]):
                    mstar=star_param_table['Mass-Flame'][0]
                    mstarerr=mstar*0.1 # assuming 10% error
                    K=gaiansstable['semi_amplitude_primary'][j]
                    Kerr=gaiansstable['semi_amplitude_primary_error'][j]
                    P=gaiansstable['period'][j]
                    Perr=gaiansstable['period_error'][j]
                    ecc=gaiansstable['eccentricity'][j]
                    eccerr=gaiansstable['eccentricity_error'][j]
                    if np.ma.isarray(ecc):
                        if ecc.mask:
                            ecc,eccerr=0,0
                    mplanet = K*1e3/28.4329 * np.sqrt(1-ecc**2) * mstar**(2/3) * (P/365.25)**(1/3)
                    mplaneterr = mplanet * np.sqrt( (Kerr/K)**2 + (ecc*eccerr/(1-ecc**2))**2 + (2/3 * mstarerr/mstar)**2 + (1/3 * Perr/P)**2  )
                    dataplanet=vstack([dataplanet,Table([[mplanet],[mplaneterr]], names=('mass_companion_mjup','mass_companion_mjup_error'))])
                else:
                    nullplanet=Table([[0],[0]], names=('mass_companion_mjup','mass_companion_mjup_error'))
                    nullplanet=Table(nullplanet, masked=True, copy=False)
                    nullplanet['mass_companion_mjup'].mask=True
                    nullplanet['mass_companion_mjup_error'].mask=True
                    dataplanet=vstack([dataplanet,nullplanet])
        else:
            nullplanet=Table([[0],[0]], names=('mass_companion_mjup','mass_companion_mjup_error'))
            nullplanet=Table(nullplanet, masked=True, copy=False)
            nullplanet['mass_companion_mjup'].mask=True
            nullplanet['mass_companion_mjup_error'].mask=True
            dataplanet=vstack([dataplanet,nullplanet])
            
            
        
        # group everything in a big Table
        maindata=simbad_table.copy()
        for i in range(1,len(dataplanet)):
            NSSflagtable=vstack([NSSflagtable,NSSflagtable])
            listname=vstack([listname,listname])
            star_param_table=vstack([star_param_table,star_param_table])
            maindata=vstack([maindata,maindata])
        interdata=hstack([listname,maindata,NSSflagtable,dataplanet,gaiansstable,gaianssaccelerationtable,star_param_table])
        if flaginit:
            newdata=interdata.copy()
            flaginit=False
        else:
            newdata=vstack([newdata,interdata])
            
    newdata.write(filename.split('.')[0]+'_gaianss.txt', overwrite=True,format='ascii.fixed_width_two_line')
    newdata.write(filename.split('.')[0]+'_gaianss.csv', overwrite=True,format='ascii.csv',delimiter=',')
    print('')
    newdata.pprint()
    print('')
        

if __name__=='__main__':
    main()


