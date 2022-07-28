# /!\ always do this in the terminal before : export SSL_CERT_FILE=/usr/local/lib/python3.9/site-packages/pip/_vendor/certifi/cacert.pem

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
        interdata=Table([[data['name'][i]]],names=('name',))
            
        # search simbad
        customSimbad = Simbad()
        customSimbad.add_votable_fields('ra(:;A;ICRS;J2016;2000)', 'dec(:;D;ICRS;J2016;2000)', 'id(HIP)')
        customSimbad.remove_votable_fields('coordinates')
        try:
            result_table = customSimbad.query_object(names_star[i])
            if result_table is None:
                print('Object not found in SIMBAD'.format(names_star[i]),end='; ',flush=True)
                continue
            else:
                print('Object found in SIMBAD'.format(names_star[i]),end='; ',flush=True)
        except:
            print('error: SIMBAD retrieval failed')
            continue
        
        # now search Gaia catalog
        coord = SkyCoord(ra=result_table['RA___A_ICRS_J2016_2000'].data[0], dec=result_table['DEC___D_ICRS_J2016_2000'].data[0], unit=(u.hourangle, u.deg), frame='icrs',equinox='J2000')
        radius = u.Quantity(5, u.arcsec)
        try:
            queryTable_main = Gaia.query_object(coordinate=coord, radius=radius)
            if queryTable_main is None:
                print('Object not found in Gaia EDR3')
                continue
            else:
                print('Object found in Gaia EDR3',end='; ',flush=True)
        except:
            print('error: Gaia EDR3 retrieval failed')
            continue
        Nfound=len(queryTable_main['source_id'].data)
        if Nfound>1:
            print('{} entries : selecting brightest'.format(Nfound),end='; ',flush=True)
            ibrightest=np.argsort(queryTable_main['phot_g_mean_mag'].data)[0]
            queryTable_main=queryTable_main[ibrightest:ibrightest+1]
        source_id=queryTable_main['source_id'].data[0]
        designation=Table([[queryTable_main['DESIGNATION'].data[0]]],names=('DESIGNATION',))
        
        # now search the NSS table
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
            if len(gaiansstable)==0:
                print('not found in NSS table',end='; ',flush=True)
                gaiansstable.add_row(np.ma.array(np.arange(len(gaiansstable.colnames)),mask=True))
                flagNSS=False
            else:
                print('found in NSS table',end='; ',flush=True)
                flagNSS=True
        except:
            print('error: Gaia NSS table retrieval failed')
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
           
        # finally calculates the mass of the companion if not SB2
        flagmassplanet=np.zeros((len(gaiansstable),)).astype(bool)
        dataplanet=Table([[],[]], names=('mass_companion_mjup','mass_companion_mjup_error'))
        if flagNSS and flagmass:
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
        for i in range(1,len(dataplanet)):
            interdata=vstack([interdata,interdata])
            star_param_table=vstack([star_param_table,star_param_table])
            designation=vstack([designation,designation])
        interdata=hstack([interdata,designation,dataplanet,gaiansstable,star_param_table])
        if flaginit:
            newdata=interdata.copy()
            flaginit=False
        else:
            newdata=vstack([newdata,interdata])
            
    newdata.write(filename.split('.')[0]+'_gaianss.txt', overwrite=True, format='ascii.tab')
    print('')
    newdata.pprint()
    print('')
        

if __name__=='__main__':
    main()


