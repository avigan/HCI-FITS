import numpy as np
import os

import astropy.table as table

from astroquery.simbad import Simbad
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.io import fits
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
    telescope = telescopes.get(instrument.lower(), None)

    if telescope is None:
        print('Warning: telescope cannot be infered from instrument {0}. '.format(instrument) +
              'You will have to set it manually.')
    else:
        print('Warning: automatically set telescope to {0} for instrument {1}'.format(telescope, instrument))
    
    return telescope


class HCIExtension:
    '''
    HCI-FITS extension class
    '''

    ##################################################
    # Constructor
    ##################################################
    
    def __init__(self, exttype, extname):

        self._name   = extname
        self._type   = exttype
        self._header = fits.Header()
        self._data   = None
    
    ##################################################
    # Properties
    ##################################################

    @property
    def name(self):
        return self._name

    @property
    def type(self):
        return self._type

    @property
    def header(self):
        return self._header

    @header.setter
    def header(self, value):
        self._header = value

    @property
    def data(self):
        return self._data

    @data.setter
    def data(self, value):
        if (self._type == 'BINTABLE') and not isinstance(value, table.Table):
            raise ValueError('Extension type is BINTABLE but the provided value is not a table.')
        elif (self._type == 'IMAGE') and not isinstance(value, np.ndarray):
            raise ValueError('Extension type is IMAGE but the provided value is not an numerical array.')
        
        self._data = value
        

class HCIDataset:
    '''
    HCI data set class
    '''
    
    ##################################################
    # Class variables
    ##################################################
    

    ##################################################
    # Constructor
    ##################################################
    
    def __init__(self):
        '''

        '''

        # properties
        self._target          = None
        self._simbad          = None
        self._coord           = None
        self._instrument      = None
        self._telescope       = None
        self._filter          = None
        self._wave            = None
        self._bandwidth       = None
        self._pixscale        = None

        # extensions
        self._data_info       = None
        self._reduced_data    = None
        self._snr_map         = None
        self._sensitivity_map = None
        self._detection_limit = None
        self._detections_info = None

    @classmethod
    def from_file(self, filename):
        ds = HCIDataset()
        # ds.init_from_file(filename)
        return ds

    @classmethod
    def from_observation(self, target, date, instrument, number_of_data, number_of_candidates=0):
        ds = HCIDataset()

        # fill basic information
        ds.target     = target
        ds.date       = date
        ds.instrument = instrument
        
        # setup data information
        c01 = table.Column(name='Image_Number', unit='', dtype='int16', length=number_of_data)
        c02 = table.Column(name='Orientation', unit='deg. e of n', dtype='float64', length=number_of_data)
        c03 = table.Column(name='Combined_rotation_angle', unit='deg.', dtype='float64', length=number_of_data)
        c04 = table.Column(name='Number_of_Exposures', unit='', dtype='int16', length=number_of_data)
        c05 = table.Column(name='Exposure_Time', unit='sec', dtype='float64', length=number_of_data)
        c06 = table.Column(name='Observation_Start', unit='Mod. Julian Date', dtype='float64', length=number_of_data)
        c07 = table.Column(name='Observation_End', unit='Mod. Julian Date', dtype='float64', length=number_of_data)
        c08 = table.Column(name='UT_Midpoint_Date_of_Observation', unit='yyyy-mm-dd', dtype='a20', length=number_of_data)
        c09 = table.Column(name='UT_Midpoint_Time_of_Observation', unit='hh:mm:ss', dtype='a20', length=number_of_data)
        c10 = table.Column(name='Wavelength', unit='nm', dtype='float64', length=number_of_data)
        c11 = table.Column(name='Bandwidth', unit='nm', dtype='float64', length=number_of_data)
        c12 = table.Column(name='Polarization', unit='', dtype='a20', length=number_of_data)

        data_info_table = table.Table()
        data_info_table.add_columns([c01, c02, c03, c04, c05, c06, c07, c08, c09, c10, c11, c12])        

        ds._data_info      = HCIExtension('BINTABLE', 'DATA_INFORMATION')
        ds._data_info.data = data_info_table
        
        # setup source detections
        if number_of_candidates > 0:
            c01 = table.Column(name='Candidate', unit='', dtype='int16', length=number_of_candidates)
            c02 = table.Column(name='SNR', unit='', dtype='float64', shape=(number_of_data, ), length=number_of_candidates)
            c03 = table.Column(name='dRA', unit='as', dtype='float64', shape=(number_of_data, ), length=number_of_candidates)
            c04 = table.Column(name='err_dRA', unit='as', dtype='float64', shape=(number_of_data, ), length=number_of_candidates)
            c05 = table.Column(name='dDEC', unit='as', dtype='float64', shape=(number_of_data, ), length=number_of_candidates)
            c06 = table.Column(name='err_dDEC', unit='as', dtype='float64', shape=(number_of_data, ), length=number_of_candidates)
            c07 = table.Column(name='Sep', unit='as', dtype='float64', shape=(number_of_data, ), length=number_of_candidates)
            c08 = table.Column(name='err_Sep', unit='as', dtype='float64', shape=(number_of_data, ), length=number_of_candidates)
            c09 = table.Column(name='PA', unit='deg', dtype='float64', shape=(number_of_data, ), length=number_of_candidates)
            c10 = table.Column(name='err_PA', unit='deg', dtype='float64', shape=(number_of_data, ), length=number_of_candidates)
            c11 = table.Column(name='Flux_cs', unit='c/s', dtype='float64', shape=(number_of_data, ), length=number_of_candidates)
            c12 = table.Column(name='err_Flux_cs', unit='c/s', dtype='float64', shape=(number_of_data, ), length=number_of_candidates)
            c13 = table.Column(name='Flux_mag', unit='mag', dtype='float64', shape=(number_of_data, ), length=number_of_candidates)
            c14 = table.Column(name='err_Flux_mag', unit='mag', dtype='float64', shape=(number_of_data, ), length=number_of_candidates)
            c15 = table.Column(name='Flux_Jy', unit='Jy', dtype='float64', shape=(number_of_data, ), length=number_of_candidates)
            c16 = table.Column(name='err_Flux_Jy', unit='Jy', dtype='float64', shape=(number_of_data, ), length=number_of_candidates)
            c17 = table.Column(name='Flux_erg', unit='erg/cm2/s/A', dtype='float64', shape=(number_of_data, ), length=number_of_candidates)
            c18 = table.Column(name='err_Flux_erg', unit='erg/cm2/s/A', dtype='float64', shape=(number_of_data, ), length=number_of_candidates)
            c19 = table.Column(name='Contrast', unit='dmag', dtype='float64', shape=(number_of_data, ), length=number_of_candidates)
            c20 = table.Column(name='err_Contrast', unit='dmag', dtype='float64', shape=(number_of_data, ), length=number_of_candidates)

            source_detection_table = table.Table()
            source_detection_table.add_columns([c01, c02, c03, c04, c05, c06, c07, c08, c09, c10,
                                                c11, c12, c13, c14, c15, c16, c17, c18, c19, c20])

            ds._detections_info      = HCIExtension('BINTABLE', 'SOURCE_DETECTION')
            ds._detections_info.data = source_detection_table

        return ds
    
        
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

    @telescope.setter
    def telescope(self, value):
        self._telescope = value
    
    @property
    def instrument(self):
        return self._instrument

    @instrument.setter
    def instrument(self, instrument):
        self._instrument = instrument
        
        # automatically fill telescope if possible
        self._telescope = resolve_telescope(instrument)

    @property
    def date(self):
        return self._date

    @date.setter
    def date(self, date):
        self._date = Time(date)
    
    @property
    def filter(self):
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

    @target.setter
    def target(self, target):
        self._target = target

        # try resolving the target with SIMBAD
        self.resolve_target()

    @property
    def simbad_target(self):
        return self._simbad

    @property
    def coordinates(self):
        return self._coord

    @coordinates.setter
    def coordinates(self, coordinates):
        if isinstance(coordinates, SkyCoord):
            self._coord = coordinates
        elif isinstance(coordinates, (tuple, list, np.ndarray)):
            self._coord = SkyCoord(coordinates[0], coordinates[1], frame='icrs', unit=(u.deg, u.deg))
        else:
            raise ValueError('Coordinates are not a SkyCoord object or 2-elements list/tuple/array')
    
    @property
    def ra(self):
        return self._coord.ra.deg
    
    @property
    def dec(self):
        return self._coord.dec.deg

    ##################################################
    # Generic class methods
    ##################################################

    def resolve_target(self):
        target = self._target
        if target is None:
            print('Warning: target name not set. Cannot resolve coordinates.')
            return
        
        # target coordinates
        info = Simbad.query_object(target)
        if info is None:
            print('Error: target {0} could not be resolved by SIMBAD. '.format(target) +
                  'Target name will be saved but you must enter manual coordinates.')
        else:
            self._coord  = SkyCoord(info['RA'][0], info['DEC'][0], frame='icrs', unit=(u.hourangle, u.deg))
            self._simbad = info['MAIN_ID'][0].decode('utf-8')

    def insert_data_information():
        pass
    
    
ds = HCIDataset.from_observation('HIP65426', '2017-03-01', 'SPHERE', 2, number_of_candidates=3)


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

