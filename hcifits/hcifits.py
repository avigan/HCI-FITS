import numpy as np
import os

from astroquery.simbad import Simbad
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.io import fits
from astropy.table import Table
from astropy.table import Column, MaskedColumn
from astropy.time import Time
from astropy import wcs
from pathlib import Path


hci_version = 1.0


def resolve_telescope(instrument):    
    telescopes = {'naco':    'VLT',
                  'sphere':  'VLT',
                  'visir':   'VLT',
                  'sinfoni': 'VLT',
                  'niri':    'Gemini-N',
                  'nici':    'Gemini-S',
                  'gpi':     'Gemini-S',
                  'nirc2':   'Keck-II',
                  'nicmos':  'HST',
                  'clio':    'MMT',
                  'irac':    'Spitzer'}

    return telescopes.get(instrument.lower(), None)


class HCIDataset:
    '''
    HCI data set class
    '''

    # header = None
    # data_info = None
    # reduced_data = None
    # snr_map = None
    # sensitivity_map = None
    # detection_limit = None
    # candidates = None

    
    ##################################################
    # Class variables
    ##################################################
    

    ##################################################
    # Constructor
    ##################################################
    
    def __init__(self, target, date, instrument, coordinates=None):
        '''

        '''

        # target name
        self._target = target
        self._simbad = None
        
        # target coordinates
        info = Simbad.query_object(target)
        if info is None:
            if coordinates is None:
                raise ValueError('Target {0} could not be resolved by SIMBAD. '.format(target) +
                                 'Check the target name or pass a (ra, dec) tuple with values in degrees.')
            elif isinstance(coordinates, SkyCoord):
                self._coord = coordinates
            elif isinstance(coordinates, tuple):
                self._coord = SkyCoord(coordinates[0], coordinates[1], frame='icrs', unit=(u.deg, u.deg))
            else:
                raise ValueError('coordinates or not a SkyCoord object or 2-tuple value')
        else:
            self._coord = SkyCoord(info['RA'][0], info['DEC'][0], frame='icrs', unit=(u.hourangle, u.deg))
            self._simbad = info['MAIN_ID'][0].decode('utf-8')

        # instrument
        self._instrument = instrument
        self._telescope = resolve_telescope(instrument)
        
    ##################################################
    # Representation
    ##################################################
    
    def __repr__(self):
        return ''

    def __format__(self):
        return self.__repr__()

    ##################################################
    # Properties
    ##################################################

    @property
    def telescope(self):
        return self._telescope
    
    @property
    def instrument(self):
        return self._instrument

    @property
    def filt(self):
        return self._filter

    @property
    def wave(self):
        return self._wave

    @property
    def bandwidth(self):
        return self._bandwidth

    @property
    def pixscale(self):
        return self._pixscale
    
    @property
    def target(self):
        return self._target

    @property
    def simbad_target(self):
        return self._simbad
    
    @property
    def ra(self):
        return self._coord.ra.deg
    
    @property
    def dec(self):
        return self._coord.dec.deg

    ##################################################
    # Generic class methods
    ##################################################

    def initFromFile(self, filename):
        pass

    

ds = HCIDataset('HIP65426', '2017-03-02', 'SPHERE')

    
##################################################
# DIVA formating
##################################################

# if (resolve_targets is True):
#     print('Resolving targets:')
#     targets = data['target'].unique()
#     for t in targets:
#         print(' * {0}'.format(t))

#         # special case for close binaries
#         if (t == 'HIP21632B') or (t == 'HD77407B') or (t == 'HD125161B') or \
#            (t == 'HIP12787B') or (t == 'HIP16563B') or (t == 'HIP39896B') or \
#            (t == 'HIP63253B') or (t == 'HIP76629B') or (t == 'HIP84586B'):
#             nt = t[:-1]
#             info = Simbad.query_object(nt)
#             resolved = False
#         else:
#             info = Simbad.query_object(t)
#             resolved = True
        
#         coord = SkyCoord(info['RA'][0], info['DEC'][0], frame='icrs', unit=(u.hourangle, u.deg))

# #
# # read, format and save files
# #

# # add high-level science product column
# data['hlsp_file'] = ''

# for index, row in data.iterrows():    
#     #
#     # create main header and fill it
#     #
#     today = Time(Time.now(), format='iso', out_subfmt='date')
    
#     main_hdr = fits.Header()

#     main_hdr.append(end=True)
#     main_hdr.add_blank('FILE INFORMATION')
#     main_hdr.append(end=True)
#     main_hdr.append(('FILETYPE', 'HLSP', 'type of data found in data file'), end=True)
#     main_hdr.append(('ORIGIN', 'DIVA Team', 'FITS file originator'), useblanks=False, bottom=True)
#     main_hdr.append(('DATE', today.value, 'date this file was written (yyyy-mm-dd)'), useblanks=False, bottom=True)
#     main_hdr.append(end=True)
#     main_hdr.add_blank('FITS FILE DATA STRUCTURE')
#     main_hdr.append(end=True)
#     main_hdr.append(('NEXTEND', len(extension_hdus_full), 'number of standard Extensions'), end=True)
#     main_hdr.append(end=True)
#     main_hdr.add_blank('PROGRAM AND INSTRUMENT INFORMATION')
#     main_hdr.append(end=True)
#     main_hdr.append(('TELESCOP', row['telescope'], 'telescope used to acquire data'), end=True)
#     main_hdr.append(('INSTRUME', row['instrument'], 'instrument used to acquire data'), end=True)
#     main_hdr.append(('FILTER', row['filter'], 'filter used during observation'), end=True)
#     main_hdr.append(end=True)
#     main_hdr.add_blank('SURVEY INFORMATION')
#     main_hdr.append(end=True)
#     main_hdr.append(('SURVNAME', row['survey'], 'name of the survey'), end=True)
#     main_hdr.append(('SURVPI', pi, 'last name of the survey principal investigator'), end=True)
#     main_hdr.append(('OBSSTGY', row['obs_strategy'], 'observing strategy'), end=True)
#     main_hdr.append(('BIBREF', row['bibref'], 'bibliographic reference for the survey'), end=True)
#     main_hdr.append(end=True)
#     main_hdr.add_blank('TARGET INFORMATION')
#     main_hdr.append(end=True)
#     main_hdr.append(('TARGNAME', target, 'target name'), end=True)
#     main_hdr.append(('RA_TARG', row['RA2000'], 'RA of target from SIMBAD (deg) (J2000)'), end=True)
#     main_hdr.append(('DEC_TARG', row['DEC2000'], 'DEC of target from SIMBAD (deg) (J2000)'), end=True)
#     main_hdr.append(('EQUINOX', 2000, 'equinox of celestial coord. system'), end=True)
#     main_hdr.append(('RESOLVED', row['resolved'], 'target fully resolved by SIMBAD'), end=True)
#     main_hdr.append(end=True)
#     main_hdr.add_blank('INFORMATION ON OTHER ASTROPHYSICAL SOURCES')
#     main_hdr.append(end=True)
#     main_hdr.append(('CANDNUM', candnum, '# of point sources detections (-1: unknown)'), end=True)
#     main_hdr.append(end=True)
#     main_hdr.add_blank('PHOTOMETRIC INFORMATION')
#     main_hdr.append(end=True)
#     main_hdr.append(('PHOTPLAM', filt_wl, 'pivot wavelength of the photmode (nm)'), end=True)
#     main_hdr.append(('PHOTBW', filt_bw, 'pivot wavelength of the photmode (nm)'), end=True)
#     main_hdr.append(end=True)
#     main_hdr.add_blank('ASTROMETRIC INFORMATION')
#     main_hdr.append(end=True)
#     main_hdr.append(('PIXSCALE', row['pixel_scale'], 'pixel scale [as/pixel]'), end=True)
#     main_hdr.append(end=True)
#     main_hdr.add_blank('ADDITIONAL COMMENTS')
#     main_hdr.append(end=True)
#     main_hdr.append(('COMMENT', row['comment']), end=True)

    
#     #
#     # obs information table
#     #    
    
#     col1  = fits.Column(name='Image_Number', format='K', unit='', array=[1])
#     col2  = fits.Column(name='Orientation', format='D', unit='deg. e of n', array=[pa])
#     col3  = fits.Column(name='Combined_rotation_angle', format='D', unit='deg.', array=[row['field_rotation']])
#     col4  = fits.Column(name='Number_of_Exposures', format='K', unit='', array=[-1])
#     col5  = fits.Column(name='Exposure_Time', format='D', unit='sec', array=[row['exposure_time']])
#     col6  = fits.Column(name='Observation_Start', format='D', unit='Mod. Julian Date', array=[mjd])
#     col7  = fits.Column(name='Observation_End', format='D', unit='Mod. Julian Date', array=[mjd])
#     col8  = fits.Column(name='UT_Midpoint_Date_of_Observation', format='A10', unit='yyyy-mm-dd', array=[row['date']])
#     col9  = fits.Column(name='UT_Midpoint_Time_of_Observation', format='A8', unit='hh:mm:ss', array=['00:00:00'])
#     col10 = fits.Column(name='Wavelength', format='D', unit='nm', array=[filt_wl])
#     col11 = fits.Column(name='Bandwidth', format='D', unit='nm', array=[filt_bw])
#     col12 = fits.Column(name='Polarization', format='A2', unit='', array=['I'])
#     cols  = fits.ColDefs([col1, col2, col3, col4, col5, col6, col7, col8, col9, col10, col11, col12])

#     # extension info
#     extname = 'DATA_INFORMATION'
#     extension_names.append(extname)
    
#     # create HDU
#     nhdu = fits.BinTableHDU.from_columns(cols)
#     nhdu.header.insert('TFIELDS', ('EXTNAME', extname, 'extension name'), after=True)
#     extension_hdus_full.append(nhdu)
    
#     # update main header
#     main_hdr.insert('NEXTEND', ('EXT{0}TYPE'.format(len(extension_hdus_full)), nhdu.header['XTENSION'],
#                                 'Extension {0} type'.format(len(extension_hdus_full))), after=True)
#     main_hdr.insert('NEXTEND', ('EXT{0}NAME'.format(len(extension_hdus_full)), extname,
#                                 'Extension {0} name'.format(len(extension_hdus_full))), after=True)
    
#     #
#     # final reduced file
#     #
#     if isinstance(img_file, str):
#         hdu = fits.open(os.path.join(root_ref, img_file[1:]))
        
#         # extension info
#         extname = 'REDUCED_DATA'
#         extension_names.append(extname)

#         # create HDU
#         nhdu = fits.ImageHDU(data=hdu[0].data)
#         nhdu.header.append(('EXTNAME', extname, 'extension name'))
#         nhdu.header.append(('BUNIT', 'UNKNOWN', 'brightness units'))

#         redstrat = []
#         if (row['adi'] == True):
#             redstrat.append('ADI')
#         if (row['sdi'] == True):
#             redstrat.append('SDI')
#         redstrat = ', '.join(redstrat)
#         nhdu.header.append(('REDSTRAT', redstrat, 'reduction strategy'))
#         nhdu.header.append(('REDALGO', row['algo'], 'reduction strategy'))

#         # Durkan survey has valid WCS
#         if survey == 'Durkan-Spitzer':
#             w = wcs.WCS(hdu[0].header)
#             nhdu.header.extend(w.to_header())
        
#         extension_hdus_full.append(nhdu)
        
#         # update main header
#         main_hdr.insert('EXT{0}TYPE'.format(len(extension_hdus_full)-1),
#                         ('EXT{0}TYPE'.format(len(extension_hdus_full)), nhdu.header['XTENSION'],
#                          'Extension {0} type'.format(len(extension_hdus_full))), after=True)
#         main_hdr.insert('EXT{0}TYPE'.format(len(extension_hdus_full)-1),
#                         ('EXT{0}NAME'.format(len(extension_hdus_full)), extname,
#                          'Extension {0} name'.format(len(extension_hdus_full))), after=True)
    
#     #
#     # psf file
#     #
#     if isinstance(psf_file, str):
#         hdu = fits.open(os.path.join(root_ref, psf_file[1:]))

#         # extension info
#         extname = 'PSF_IMAGE'
#         extension_names.append(extname)

#         # create HDU
#         nhdu = fits.ImageHDU(data=hdu[0].data)
#         nhdu.header.append(('EXTNAME', extname, 'extension name'))
#         nhdu.header.append(('BUNIT', 'UNKNOWN', 'brightness units'))

#         # Durkan survey has valid WCS
#         if survey == 'Durkan-Spitzer':
#             w = wcs.WCS(hdu[0].header)
#             nhdu.header.extend(w.to_header())
                
#         extension_hdus_full.append(nhdu)
        
#         # update main header
#         main_hdr.insert('EXT{0}TYPE'.format(len(extension_hdus_full)-1),
#                         ('EXT{0}TYPE'.format(len(extension_hdus_full)), nhdu.header['XTENSION'],
#                          'Extension {0} type'.format(len(extension_hdus_full))), after=True)
#         main_hdr.insert('EXT{0}TYPE'.format(len(extension_hdus_full)-1),
#                         ('EXT{0}NAME'.format(len(extension_hdus_full)), extname,
#                          'Extension {0} name'.format(len(extension_hdus_full))), after=True)

#     #
#     # sensitivity map file
#     #
#     if isinstance(detmap_file, str):
#         hdu = fits.open(os.path.join(root_ref, detmap_file[1:]))

#         # extension info
#         extname = 'SENSITIVITY_MAP'
#         extension_names.append(extname)

#         # apparent mag for Heinze survey
#         if (survey == 'Heinze-10') or (survey == 'Durkan-Spitzer'):
#             unit = 'mag'
#         else:
#             unit = 'dmag'
        
#         # create HDU
#         nhdu = fits.ImageHDU(data=hdu[0].data)
#         nhdu.header.append(('EXTNAME', extname, 'extension name'))
#         nhdu.header.append(('BUNIT', unit, 'brightness units'))
#         nhdu.header.append(('NSIGMA', 1.0, 'sensitivity limit level'))
#         if (survey == 'Heinze-10') or (survey == 'Durkan-Spitzer'):
#             nhdu.header.append(('COMMENT', row['comment']), end=True)
#         extension_hdus_full.append(nhdu)

#         # update main header
#         main_hdr.insert('EXT{0}TYPE'.format(len(extension_hdus_full)-1),
#                         ('EXT{0}TYPE'.format(len(extension_hdus_full)), nhdu.header['XTENSION'],
#                          'Extension {0} type'.format(len(extension_hdus_full))), after=True)
#         main_hdr.insert('EXT{0}TYPE'.format(len(extension_hdus_full)-1),
#                         ('EXT{0}NAME'.format(len(extension_hdus_full)), extname,
#                          'Extension {0} name'.format(len(extension_hdus_full))), after=True)
            
#     #
#     # detlim file
#     #
#     if isinstance(detlim_file, str):
#         hdu = fits.open(os.path.join(root_ref, detlim_file[1:]))

#         if (survey == 'Durkan-Spitzer'):
#             # sep already in as
#             sep = hdu[0].data[:, 0].squeeze()
            
#             # detlim @ 5-sigma ==> 1-sigma
#             dlim = hdu[0].data[:, 1].squeeze() + np.log10(5)
#         else:
#             # sep in mas
#             sep = hdu[0].data[:, 0].squeeze()/1000.

#             # detlim @ 1-sigma
#             dlim = hdu[0].data[:, 1].squeeze()
        
#         # apparent mag for Heinze survey
#         if (survey == 'Heinze-10') or (survey == 'Durkan-Spitzer'):
#             unit = 'mag'
#         else:
#             unit = 'dmag'
                
#         # create table        
#         col1  = fits.Column(name='Radius', format='D', unit='as', array=sep)
#         col2  = fits.Column(name='Detection_limit', format='D', unit=unit, array=dlim)
#         cols  = fits.ColDefs([col1, col2])

#         # extension info
#         extname = 'DETECTION_LIMIT'
#         extension_names.append(extname)
         
#         # create HDU
#         nhdu = fits.BinTableHDU.from_columns(cols)
#         nhdu.header.insert('TFIELDS', ('NSIGMA', 1.0, 'detection limit level'), after=True)
#         nhdu.header.insert('TFIELDS', ('EXTNAME', extname, 'extension name'), after=True)
#         if (survey == 'Heinze-10') or (survey == 'Durkan-Spitzer'):
#             nhdu.header.append(('COMMENT', row['comment']), end=True)
#         extension_hdus_full.append(nhdu)
    
#         # update main header
#         main_hdr.insert('EXT{0}TYPE'.format(len(extension_hdus_full)-1),
#                         ('EXT{0}TYPE'.format(len(extension_hdus_full)), nhdu.header['XTENSION'],
#                          'Extension {0} type'.format(len(extension_hdus_full))), after=True)
#         main_hdr.insert('EXT{0}TYPE'.format(len(extension_hdus_full)-1),
#                         ('EXT{0}NAME'.format(len(extension_hdus_full)), extname,
#                          'Extension {0} name'.format(len(extension_hdus_full))), after=True)

#     # cc file
#     if isinstance(cc_file, str):
#         cc_table = Table.read(os.path.join(root_ref, cc_file[1:]), format='ascii', delimiter=',', data_start=3,
#                               names=['target', 'cc', 'date', 'filter', 'sep', 'sep_err', 'pa', 'pa_err',
#                                      'dmag', 'dmag_err', 'snr', 'status'])

#         # date test
#         u = np.unique(cc_table['date'])
#         if len(u) > 1:
#             print(' ==> inconsistent number of epochs in cc for target {:s} in survey {:s}'.format(target, survey))

        
#         # target
#         tcol = Column(data=[target for i in range(len(cc_table))], name='target', dtype='str')
#         cc_table.replace_column('target', tcol)

#         # filter
#         cc_table['filter'] = filt

#         # date
#         cc_table['date'] = date
        
#         # unit conversion
#         cc_table['sep'] = cc_table['sep'].astype(float)/1000.
#         cc_table['sep_err'] = cc_table['sep_err'].astype(float)/1000.

#         # masked columns
#         cc_table['dmag'] = MaskedColumn(cc_table['dmag'], mask=(cc_table['dmag'] == 0.0))
#         cc_table['dmag_err'] = MaskedColumn(cc_table['dmag_err'], mask=(cc_table['dmag_err'] == 0.0))
#         cc_table['snr'] = MaskedColumn(cc_table['snr'], mask=(cc_table['snr'] == 0.0))

#         # cc status
#         cc_table['status'][cc_table['status'] == 'PB'] = 'B'

#         # epoch
#         cc_table.add_column(Column(np.repeat(instr, len(cc_table)), name='instrument'), index=3)
#         # cc_table.add_column(Column(np.repeat(mjd, len(cc_table)), name='mjd'), index=3)
#         # cc_table.add_column(Column(np.repeat(epoch, len(cc_table)), name='epoch'), index=3)
        
#         # add dra,ddec columns
#         sep     = cc_table['sep']
#         sep_err = cc_table['sep_err']
#         pa      = np.deg2rad(cc_table['pa'] + 90)
#         pa_err  = np.deg2rad(cc_table['pa_err'])

#         dra      = -sep*np.cos(pa)
#         dra_err  = np.sqrt(np.cos(pa)**2 * sep_err**2 + sep**2*np.sin(pa)**2 * pa_err**2)
#         ddec     = +sep*np.sin(pa)
#         ddec_err = np.sqrt(np.sin(pa)**2 * sep_err**2 + sep**2*np.cos(pa)**2 * pa_err**2)
        
#         cc_table.add_column(Column(ddec_err, name='ddec_err'), index=10)
#         cc_table.add_column(Column(ddec, name='ddec'), index=10)
#         cc_table.add_column(Column(dra_err, name='dra_err'), index=10)
#         cc_table.add_column(Column(dra, name='dra'), index=10)
        
#         # create table
#         col1  = fits.Column(name='Candidate', format='K', unit='', array=cc_table['cc'])
#         col2  = fits.Column(name='SNR', format='D', unit='', array=cc_table['snr'])
#         col3  = fits.Column(name='dRA', format='D', unit='as', array=cc_table['dra'])
#         col4  = fits.Column(name='err_dRA', format='D', unit='as', array=cc_table['dra_err'])
#         col5  = fits.Column(name='dDEC', format='D', unit='as', array=cc_table['ddec'])
#         col6  = fits.Column(name='err_dDEC', format='D', unit='as', array=cc_table['ddec_err'])
#         col7  = fits.Column(name='Sep', format='D', unit='as', array=cc_table['sep'])
#         col8  = fits.Column(name='err_Sep', format='D', unit='as', array=cc_table['sep_err'])
#         col9  = fits.Column(name='PA', format='D', unit='deg', array=cc_table['pa'])
#         col10 = fits.Column(name='err_PA', format='D', unit='deg', array=cc_table['pa_err'])
#         col11 = fits.Column(name='Flux_cs', format='D', unit='c/s', array=np.repeat(np.nan, 1))
#         col12 = fits.Column(name='err_Flux_cs', format='D', unit='c/s', array=np.repeat(np.nan, 1))
#         col13 = fits.Column(name='Flux_mag', format='D', unit='mag', array=np.repeat(np.nan, 1))
#         col14 = fits.Column(name='err_Flux_mag', format='D', unit='mag', array=np.repeat(np.nan, 1))
#         col15 = fits.Column(name='Flux_Jy', format='D', unit='Jy', array=np.repeat(np.nan, 1))
#         col16 = fits.Column(name='err_Flux_Jy', format='D', unit='Jy', array=np.repeat(np.nan, 1))
#         col17 = fits.Column(name='Flux_erg', format='D', unit='erg/cm2/s/A', array=np.repeat(np.nan, 1))
#         col18 = fits.Column(name='err_Flux_erg', format='D', unit='erg/cm2/s/A', array=np.repeat(np.nan, 1))
#         col19 = fits.Column(name='Contrast', format='D', unit='dmag', array=cc_table['dmag'])
#         col20 = fits.Column(name='err_Contrast', format='D', unit='dmag', array=cc_table['dmag_err'])
#         # col21 = fits.Column(name='Status', format='A2', unit='', array=cc_table['status'])
#         cols  = fits.ColDefs([col1, col2, col3, col4, col5, col6, col7, col8, col9, col10, col11,
#                               col12, col13, col14, col15, col16, col17, col18, col19, col20])

#         # extension info
#         extname = 'SOURCE_DETECTION'
#         extension_names.append(extname)
    
#         # create HDU
#         nhdu = fits.BinTableHDU.from_columns(cols)
#         nhdu.header.insert('TFIELDS', ('EXTNAME', extname, 'extension name'), after=True)
#         extension_hdus_full.append(nhdu)

#         # number of candidates
#         candnum = len(cc_table)
        
#         # update main header
#         main_hdr['CANDNUM'] = candnum
#         data.loc[index, 'candidates'] = candnum
#         main_hdr.insert('EXT{0}TYPE'.format(len(extension_hdus_full)-1),
#                         ('EXT{0}TYPE'.format(len(extension_hdus_full)), nhdu.header['XTENSION'],
#                          'Extension {0} type'.format(len(extension_hdus_full))), after=True)
#         main_hdr.insert('EXT{0}TYPE'.format(len(extension_hdus_full)-1),
#                         ('EXT{0}NAME'.format(len(extension_hdus_full)), extname,
#                          'Extension {0} name'.format(len(extension_hdus_full))), after=True)
        
#     # update table
#     data.loc[index, 'candidates'] = candnum
        
#     # empty dates
#     if isinstance(date, float):
#         date = '0000-00-00'

#     # create directory structure
#     path_new = os.path.join(root_new, 'Data', survey)
#     subpath_new = os.path.join('Data', survey)
#     if not os.path.exists(path_new):
#         os.makedirs(path_new)

#     # write HLSP data (full)
#     main_hdr['NEXTEND'] = len(extension_hdus_full)
    
#     main_hdu = fits.PrimaryHDU(header=main_hdr)
#     all_hdus = [main_hdu]
#     all_hdus.extend(extension_hdus_full)
#     list_hdu = fits.HDUList(all_hdus)
    
#     fname = target+'_'+date+'_'+instr+'_'+filt+'_diva_hlsp.fits'
#     list_hdu.writeto(os.path.join(path_new, fname), overwrite=overwrite)
#     data.loc[index, 'hlsp_file'] = os.path.join(subpath_new, fname)

