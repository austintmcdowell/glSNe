import numpy as np
import matplotlib.pyplot as plt
import os
import glob
import psycopg2
from astropy.io import fits
from mpi4py import MPI


def main():

    names, RA, DEC = read_lens()

    done_list = []
    for i in xrange(len(names)):
        print i
        read_ps1( names[i] )
        
        done_list.append(str(i))

        with open("/scratch1/scratchdirs/amcd/images/checkpoint_file.txt", "w") as f:
            f.write( str(done_list) )
        


def read_ps1( name ):
# read detected star data from ps1

    _split = lambda iterable, n: [iterable[:len(iterable)/n]] + \
        _split(iterable[len(iterable)/n:], n - 1) if n != 0 else []

    os.chdir('/scratch1/scratchdirs/amcd/images/'+name+'_dir/')

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()


    cat_files = glob.glob("*.cat")
    cat_files.sort()

    my_job = _split(cat_files, size)[rank]

    con = psycopg2.connect(host='scidb2.nersc.gov', user='desi_user', password='decals_Fun!', database='desi')

    with con:
        
        cursor = con.cursor()

        
        for job in my_job:
            
            band, fits_file = filter_type( name, job )

            band = band.lower()
            
            table = read_se( name, job )
            query_results = np.zeros_like( table, dtype=[('PS1_MAG', 'f'), ('PS1_MAGERR', 'f'), ('MAG_CHECK', 'int')] )
            #return table

            print 'working on '+job

            for table_row, dbrow in zip(table, query_results):
            # each row is one .cat entry
    
                query = "select ra, dec, %s_median, %s_stdev from dr1.ps1 where q3c_radial_query(ra, dec, %s, %s, 2/3600.) order by q3c_dist(ra, dec, %s, %s) asc;"

                cursor.execute( query % (band, band, table_row['X_WORLD'], 
                                         table_row['Y_WORLD'], table_row['X_WORLD'], 
                                         table_row['Y_WORLD']))
        
                result = cursor.fetchall()
                
                if len(result) != 0:
                    result = result[0]

                    dbrow['PS1_MAG'] = result[2]
                    if (result[2] > 17) and (result[2] < 20):
                        dbrow['MAG_CHECK'] = 1

                    dbrow['PS1_MAGERR'] = result[3]
                    
                    #return dbrow

            ob_zps = query_results['PS1_MAG'][query_results['MAG_CHECK']==1] + 2.5*np.log10(table['FLUX_AUTO'][query_results['MAG_CHECK']==1])
            
            med_zp = np.median(ob_zps)

            absdev = np.abs( ob_zps - med_zp )
            mad = np.median( absdev )
            sigma = 1.418*mad
            
            print 'writing to '+ fits_file

            with fits.open( fits_file, mode='update' ) as hdulist:
            #data, header = fits.getdata( fits_file, header=True )
                hdulist[0].header['AM_ZP'] = str(med_zp)
                hdulist[0].header['AM_ZPERR'] = str(sigma)

            #fits.writeto( fits_file, data, header, clobber=True )

            print 'finishing '+str(job)
            print '============'
            #return 0

            

def filter_type( name, cat_file ):
# gets list of filters types for source images

    os.chdir('/scratch1/scratchdirs/amcd/images/'+name+'_dir/')
    
    start = cat_file.find('_0')+1
    end = cat_file.find('.c')

    fits_file = name+'_image'+cat_file[start:end]+'.fits'
    
    with fits.open( fits_file ) as image:
        filter = image[0].header['filter']

    return filter, fits_file


def read_se( name, cat_file ):
# read in relevant data from SE results for one source, one image

    os.chdir('/scratch1/scratchdirs/amcd/images/'+name+'_dir/')

    dtype = [('X_IMAGE', 'f'),
             ('Y_IMAGE', 'f'),
             ('X_WORLD', 'f'),
             ('Y_WORLD', 'f'),
             ('A_IMAGE', 'f'),
             ('B_IMAGE', 'f'),
             ('FWHM_WORLD', 'f'),
             ('FLAGS', 'int'),
             ('FLUX_AUTO', 'f'),
             ('FLUXERR_AUTO', 'f')]

    array = np.genfromtxt( cat_file, dtype=dtype )
    return array[array['FLAGS']==0]

def read_lens():
# read in names, RA, DEC from masterlens1.txt

    lst = open('/scratch1/scratchdirs/amcd/images/masterlens1.txt')
    l = lst.readlines()

    names = []
    for i in xrange(len(l)):
        names.append( l[i].split()[0] )

    RA = []
    for i in xrange(len(l)):
        RA.append( float(l[i].split()[1]) )

    DEC = []
    for i in xrange(len(l)):
        DEC.append( float(l[i].split()[2]) )

    return names, np.array(RA), np.array(DEC)




main()
