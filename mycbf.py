"""
mycbf.py
--------
An augmented version of the pycbf.cbf_handle_struct
"""

__author__ = "Frank Murphy"
__copyright__ = "Copyright 2012 NE-CAT"
__license__ = "BSD"
__version__ = "3clause"
__maintainer__ = "Frank Murphy"
__email__ = "fmurphy@anl.gov"
__status__ = "Development"
__date__ = "2012/12/13"

import pycbf
import re
import pprint

class new_cbf_handle(pycbf.cbf_handle_struct):

    # probes are (name,regex,groups python types, groups user types)
    # following is for NE-CAT Pilatus 6M-F as of 12/13/2012
    __header_probes__ = [('DETECTOR',               'Detector\:\s+(\w[\w\,\s\/\-]+\w)',                         [str],             ['label']),
                         ('TIMESTAMP',              '(\d+\-\d+\-\d+T\d+\:\d+\:\d+\.\d+)',                       [str],             ['timestamp']),
                         ('PIXEL_SIZE',             'Pixel_size\s*([\d\-e]+)\s*m\s*x\s*([\d\e\-]+)\s*(\w+)',    [float,float,str], ['val', 'val', 'units']),
                         ('SENSOR_THICKNESS',       'Silicon sensor\, thickness ([\d\.]+)\s*(\w+)',             [float,str],       ['val','units']),
                         ('EXPOSURE_TIME',          'Exposure_time\s*([\d\.]+)\s*(\w+)',                        [float,str],       ['val', 'units']),
                         ('EXPOSURE_PERIOD',        'Exposure_period\s*([\d\.]+)\s*(\w+)',                      [float,str],       ['val', 'units']),
                         ('TAU',                    'Tau\s*\=\s*([\d\-\.e]+)\s*(\w+)',                          [float,str],       ['val', 'units']),
                         ('COUNT_CUTOFF',           'Count_cutoff\s*([\d]+)\s*(\w+)',                           [int,str],         ['val', 'units']),
                         ('THRESHOLD',              'Threshold_setting\:\s+(\d+)\s*(\w+)',                      [int,str],         ['val', 'units']),
                         ('GAIN',                   'Gain_setting\s*\:\s*([\w\s]*)\s+\(vrf\s\=\s([\-\d\.]+)\)', [str,float],       ['val', 'val']),  
                         ('EXCLUDED_PIXELS_NUMBER', 'N_excluded_pixels\s*\=\s*(\d+)',                           [int],             ['val']), 
                         ('EXCLUDED_PIXELS_FILE',   'Excluded_pixels\:\s+([\w\.]+)',                            [str],             ['val']),
                         ('FLAT_FIELD_FILE',        'Flat_field\:\s+([\(\w\.\)]+)',                             [str],             ['val']),
                         ('TRIM_FILE',              'Trim_file\:\s+([\w\.]+)',                                  [str],             ['val']),
                         ('IMAGE_PATH',             'Image_path\:\s+(\/*[\w\.]+\/*)',                           [str],             ['val']),
                         ('WAVELENGTH',             'Wavelength\s*([\d\.]+)\s*(\w+)',                           [float,str],         ['val', 'units']),
                         ('DETECTOR_DISTANCE',      'Detector_distance\s+([\d\.]+)\s*(\w+)',                    [float,str],       ['val', 'units']),
                         ('BEAM_CENTER',            'Beam_xy\s*\(([\d\.]+)\,\s*([\d\.]+)\)\s*(\w+)',            [float,float,str], ['val', 'val', 'units']),
                         ('TRANSMISSION',           'Filter_transmission\s*([\d\.]+)',                          [float],           ['val']),    
                         ('START_ANGLE',            'Start_angle\s*([\d\.]+)\s*(\w+)\.*',                       [float,str],       ['val', 'units']),
                         ('ANGLE_INCREMENT',        'Angle_increment\s*([\d\.]+)\s*(\w+)\.*',                   [float,str],       ['val', 'units']),
                         ('TWOTHETA',               'Detector_2theta\s*([\d\.]+)\s*(\w+)\.*',                   [float,str],       ['val', 'units'])]


    def __init__(self, filename=False, test=True):
        
        if test:
            print "new_cbf_handle.__init__"
        
        pycbf.cbf_handle_struct.__init__(self)
        
        self.test = test
        
        self.__file_loaded__ = False
        self.__filename__ = None
        self.__handle_structure__ = None
        self.__raw_header__ = None
        self.__parsed_header__ = None
        
        if (filename):
            self.new_load_file(filename)
        
    def new_load_file(self, filename, mode=None):
        """
        Load data from assigned filename
        """
        if self.test:
            print "new_cbf_handle.new_load_file"
        
        self.__filename__ = filename
        
        if (not mode):
            mode = pycbf.CBF
        
        try:
            self.read_file(self.__filename__, mode)
            self.__file_loaded__ = True
        except:
            raise Exception("Error reading %s. Does it exist?" % self.__filename__)
        
    def new_read_header(self):
        """
        Read the raw header
        """
        if self.test:
            print "new_cbf_handle.new_read_header"        
        
        if not self.__file_loaded__:
            raise Exception("File not loaded.")
        
        self.new_read_cbf_structure()
        
    def new_get_raw_header(self):
        """
        Return the raw header
        """
        if self.test:
            print "new_cbf_handle.new_get_raw_header"
        
        if not self.__handle_structure__:
            self.new_read_header()
            
        for cat in self.__handle_structure__.get('contents'):
            for row in cat:
                for col in row:
                    if col.get('column_name') == "header_contents":
                        return(col.get('value'))
      
    def new_get_parsed_header(self):
        """
        Return the dict containing header information
        """
        if self.test:
            print "new_cbf_handle.new_get_parsed_header"
        
        if not self.__parsed_header__:
            self.new_read_cbf_structure()
            
        return self.__parsed_header__
        
    def new_read_cbf_structure(self):
        """
        Get the structure of the cbf in question, with any values to boot,
        but not data
        """
        if self.test:
            print "new_cbf_handle.new_read_cbf_structure"
            
        if not self.__file_loaded__:
            raise Exception("File not loaded.")
        
        handle_structure = {'contents':[]}
        
        self.select_datablock(0)
        handle_structure['datablock_name'] = self.datablock_name()
        self.rewind_category()
        categories = self.count_categories()
        for i in range(categories):
            handle_structure['contents'].append([])
            self.select_category(i)
            rows = self.count_rows()
            cols = self.count_columns()
            self.rewind_column()
            for j in range(rows):
                handle_structure['contents'][i].append([])
                self.select_row(j)
                self.rewind_column()
                for k in range(cols):
                    handle_structure['contents'][i][j].append({})
                    self.select_column(k)
                    handle_structure['contents'][i][j][k]['column_name'] = self.column_name()
                    handle_structure['contents'][i][j][k]['typeofvalue'] = self.get_typeofvalue()
                    if (handle_structure['contents'][i][j][k]['typeofvalue'].find("bnry") == -1):
                        handle_structure['contents'][i][j][k]['value'] = self.get_value()
                        if self.column_name() == "header_contents":
                            self.__raw_header__ = handle_structure['contents'][i][j][k]['value']
                            self.new_parse_header()
                    else:
                        (compression, binaryid, elsize, elsigned, \
                        elunsigned, elements, minelement, maxelement, \
                        byteorder, dimfast, dimmid, dimslow, padding) = \
                        self.get_integerarrayparameters_wdims_fs()
                        if dimfast == 0:
                            dimfast = 1
                        if dimmid == 0:
                            dimmid = 1
                        if dimslow == 0:
                            dimslow = 1
                        handle_structure['contents'][i][j][k]['parameters'] = {'compression': compression,
                                                                               'binaryid': binaryid,
                                                                               'elsize': elsize,
                                                                               'elsigned': elsigned,
                                                                               'elunsigned': elunsigned,
                                                                               'elements': elements,
                                                                               'minelement': minelement,
                                                                               'maxelement': maxelement,
                                                                               'byteorder': byteorder,
                                                                               'dimfast': dimfast,
                                                                               'dimmid': dimmid,
                                                                               'dimslow': dimslow,
                                                                               'padding': padding}
        self.__handle_structure__ = handle_structure
        
        if self.test:
            printer = pprint.PrettyPrinter(indent=2)
            printer.pprint(self.__handle_structure__)
        
    def new_parse_header(self,raw_header=None):
        """
        Parse the raw header into a more usable format
        Depends heavily on __header_probes__
        """    
        if self.test:
            print "new_cbf_handle.new_parse_header", raw_header
        
        if not raw_header:
            raw_header = self.__raw_header__
            
        if not raw_header:
            raise Exception('raw_header is empty - cannot parse.')
            
        header_dict = {}
        for label, pattern, types, descriptions in self.__header_probes__:
            pat = re.compile(pattern)
            result = pat.search(raw_header)
            if result:
                my_res = []
                for i in range(1,len(descriptions)+1):
                    my_res.append(self.__xform_type__(result.group(i), types[i-1]))
                header_dict[label] = (my_res, types, descriptions)
        
        self.__parsed_header__ = header_dict
        
        if self.test:
            printer = pprint.PrettyPrinter(indent=2)
            printer.pprint(self.__parsed_header__)

    def __xform_type__(self,value,type):
        """
        Cast an input value to an input type and return
        """
        #if self.test:
        #    print 'new_cbf_handle.__xform_type__',value,type
        
        try:
            if type == int:
                return int(value)
            elif type == float:
                return float(value)
            elif type == str:
                return str(value)
        except:
            return value
