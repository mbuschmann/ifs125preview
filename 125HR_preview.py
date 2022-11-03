# -*- coding: utf-8 -*-
"""
Standalone QT program to provide live preview of the idle mode measurements of an IFS125HR.

The program includes rudimentary FFT functionality and provides a smoothed (= low pass - filtered) interferogram 
to detect potential detector nonlinearity during alignment.

Author: Matthias Buschmann, IUP Bremen
Date: 2022/11/3
"""

from __future__ import print_function, division
import sys, yaml, requests, struct, io
from PyQt5.QtWidgets import QPushButton
from matplotlib.backends.backend_qt5agg import FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
import numpy as np

from matplotlib.backends.qt_compat import QtWidgets


class ftsreader():
    ''' A striped down version of ftsreader, including support for using a data stream directly from IFS125, 
    without the need to save to disk first.
    
    Full version at: https://github.com/mbuschmann/ftsreader
    '''
    def search_header_par(self, par):
        '''search the header for parameter <par> and return datablock designation '''
        pars = []
        for i in list(self.header.keys()):
            for j in list(self.header[i].keys()):
                if par == j:
                    pars.append(i)
        if len(pars)==1:
            return pars[0]
        elif len(pars)>1:
            if self.verbose: print('Found parameter in multiple datablocks')
            return pars
        else:
            if self.verbose: print('Parameter', par, 'not found in header!')
            return None

    def get_header_par(self, par):
        try:
            return self.header[self.search_header_par(par)][par]
        except:
            print('Parameter not found in header ...')
            return None

    def read_structure(self):
        #t = time.time()
        '''Read the structure of the file and write to ftsreader.fs'''
        # known blocks so far, there is always a block zero, that is still unidentified
        self.__blocknames =    {'160': 'Sample Parameters',
                        '23': 'Data Parameters',
                        '96': 'Optic Parameters',
                        '64': 'FT Parameters',
                        '48': 'Acquisition Parameters',
                        '32': 'Instrument Parameters',
                        '7':  'Data Block',
                        '0':  'something'}
        self.__blocknames2 = {'132': ' ScSm', # another declaration to differentiate blocks between ifg, spc, etc.
                        '4': ' SpSm',
                        '8': ' IgSm',
                        '20': ' TrSm',
                        '12': ' PhSm',
                        b'\x84': ' SpSm/2.Chn.', # some weird stuff going on with python3 decoding here, use binary representation
                        b'\x88': ' IgSm/2.Chn.'}
        self.fs = {}
        fi = self.getfileobject()
        with fi as f: #open(self.path, 'rb') as f:
            f.seek(0)
            self.log.append('Reading structure of file')
            # read beginning of file to assert magic number, total number of blocks and first offset
            # some unidentified numbers in between, do not seem to be necessary for header, spc or ifg blocks
            (magic, something, something, offset1, something, numberofblocks) = struct.unpack('6i', f.read(struct.calcsize('6i')))
            f.seek(offset1) # start at first offset
            for i in range(numberofblocks): # go through all blocks and save all found blocks in self.fs
                s = f.read(struct.calcsize('2BH2i'))
                #read beginning of block, with infos on block types, something yet unidentified/unimportant of size 'H' for now, length and gobal offset of the block
                (blocktype, blocktype2, something, length, offset2) = struct.unpack('2BH2i',s)
                blocktype = str(blocktype)
                blocktype2 = str(blocktype2)
                if blocktype in self.__blocknames.keys():
                    hdrblockname = self.__blocknames[blocktype]
                else:
                    hdrblockname = '[unknown block '+blocktype+']'
                if blocktype2 in self.__blocknames2.keys():
                    hdrblockname += self.__blocknames2[blocktype2]
                else: pass
                self.log.append('Found block '+str(blocktype)+', '+str(blocktype2)+' and identified as '+hdrblockname)
                if blocktype == '0' or blocktype not in self.__blocknames.keys():
                    hdrblockname += ' len %3i' % (length)
                else:
                    pass
                self.fs[hdrblockname] = {'blocktype': blocktype, 'blocktype2': blocktype2, 'length': length, 'offset': offset2}
        fi.close

    def getfileobject(self):
        if self.filemode == 'hdd':
            fi = open(self.path, 'rb')
        elif self.filemode == 'bytesfromfile':
            with open(self.path, 'rb') as f:
                data = f.read(17428)
            fi = io.BytesIO(data)
        elif self.filemode == 'mem':
            #print(streamdata)
            fi = io.BytesIO(self.streamdata)
        else:
            exit('filemode', self.filemode, ' not supported')
        return fi

    def getparamsfromblock(self, offset, length, full=False):
        '''Read all parameters in a block at binary <length> and <offset> and return as dictionary. On request also include binary length and offset of that parameter.'''
        params = {}
        i=0
        test = True
        fullblock = []
        fi = self.getfileobject()
        with fi as f: #with open(self.path, 'rb') as f:
            while test:
                f.seek(offset+i) # goto block offset
                s = f.read(8) # read 8 bytes
                para, thistype, length = struct.unpack('4s2H', s) # unpack to get info on how to unpack block
                if full:
                    fullblocktmp = [para, thistype, length, offset+i]
                i+=8
                if struct.unpack('4c', para)[-1]==b'\x00': #get null terminating string
                    para=para[:-1]
                else: pass
                if para[:3] != b'END' and length>0: # if not empty block
                    f.seek(offset+i)
                    data = f.read(2*length)
                    i+=2*length
                    try:
                        if thistype == 0:
                            val = struct.unpack('%1ii'%(len(data)/4), data)[0]
                        elif thistype == 1:
                            val = struct.unpack('%1id'%(len(data)/8), data)[0]
                        elif thistype >= 2 and thistype <=4:
                            t = struct.unpack('%1is'%(2*length), data)[0].decode('ISO-8859-1')
                            t2 = ''
                            for ji in t: # deal with zeros in byte array
                                if ji!='\x00' and type(ji)==str: # in python2 you might want to add ... or type(ji)=='unicode'):
                                    t2 += ji
                                else:
                                    break
                            val=t2
                        else:
                            val= '[read error]'
                        params[para.decode()] = val
                        if full:
                            fullblocktmp.append(val)
                            fullblock.append(fullblocktmp)
                    except Exception as e:
                        print('Exception in getparamsfromblock')
                        self.log.append(e)
                        print (e)
                else:
                    test = False
        if full:
            return fullblock
        else:
            return params

    def read_header(self):
        '''Read the header and return as a dictionary.'''
        self.log.append('Reading Header ...')
        self.read_structure()
        self.header = {}
        for block in self.fs.keys():
            if block[:10]!='Data Block' and self.fs[block]['length']>0: # if not data block and not empty, try reading header info
                if 'unknown' in block or 'something' in block:
                    pass
                else:
                    try:
                        self.log.append('Reading Header Block: '+block)
                        self.header[block] = self.getparamsfromblock(self.fs[block]['offset'], self.fs[block]['length'], full=False)
                    except Exception as e:
                        print(e)
                        self.log.append(e)
            else: pass
        return 0

    def get_block(self, pointer, length):
        '''Get data block from file object at <pointer> with length <length>.'''
        self.log.append('Getting data block at '+str(pointer)+' with length '+str(length))
        fi = self.getfileobject()
        with fi as f: 
            f.seek(pointer)
            dat = np.array(struct.unpack('%1if'%(length), f.read(length*4)))
        return dat

    def get_datablocks(self, block):
        '''Read a datablock named <block> and retrieve x- and y-axis np.arrays from it.'''
        #t = time.time()
        self.log.append('Getting data blocks')
        yax = np.array(self.get_block(self.search_block(block)['offset'], self.search_block(block)['length']))
        #print(block)
        if block == 'Data Block IgSm':
            self.log.append('Getting ifg data block')
            # crude estimate of opd axis, only for illustratiion purposes, zpd's not included in calculation, and triangular apod. assumption -> 0.9
            xax = np.linspace(0,2*0.9/float(self.header['Acquisition Parameters']['RES']), len(yax))
        if block == 'Data Block SpSm':
            self.log.append('Getting spc data block')
            # calculate wavenumber axis for spectrum from frequencies of first and last point stored in header
            xax = np.linspace(self.header['Data Parameters SpSm']['FXV'], self.header['Data Parameters SpSm']['LXV'], len(yax))
        if block == 'Data Block ScSm':
            self.log.append('Getting spc data block')
            xax = np.linspace(self.header['Data Parameters ScSm']['FXV'], self.header['Data Parameters ScSm']['LXV'], len(yax))
        if block == 'Data Block TrSm':
            self.log.append('Getting trm data block')
            xax = np.linspace(self.header['Data Parameters TrSm']['FXV'], self.header['Data Parameters TrSm']['LXV'], len(yax))
        if block == 'Data Block PhSm':
            self.log.append('Getting pha data block')
            xax = np.linspace(self.header['Data Parameters PhSm']['FXV'], self.header['Data Parameters PhSm']['LXV'], len(yax))
        return xax, yax

    def test_if_ftsfile(self):
        '''Check the initialized filename for FTS magic number.'''
        self.log.append('testing if FTS file')
        # same 4-byte binary representation found on all valid FTS files ... must be magic
        ftsmagicval = b'\n\n\xfe\xfe'
        try:
            fi = self.getfileobject()
            with fi as f: #with open(self.path, 'rb') as f:
                f.seek(0)
                magic = f.read(4)
            if magic==ftsmagicval:
                if self.verbose:
                    self.log.append('Identified '+self.path+' as FTS file ...')
                self.status=True
                self.isftsfile = True
            else:
                self.log.append('Bad Magic found in '+self.path)
                print('Bad Magic in ', self.path)
                self.status=False
                self.isftsfile = False
        except Exception as e:
            self.log.append(e)
            self.status=False
            self.isftsfile = False

    def search_block(self, blockname):
        '''Searches a <blockname> within the identifies FTS file structure. Returns dictionary entry of the block <blockname>.'''
        if blockname in list(self.fs.keys()):
            return self.fs[blockname]
        else:
            self.log.append('Could not find '+str(blockname)+' in self.fs.keys()')

    def has_block(self, blockname):
        '''Check if <blockname> is present in ftsreader.fs'''
        if blockname in self.fs.keys():
            return True
        else:
            return False

    def __init__(self, path, verbose=False, getspc=False, getifg=False, gettrm=False, getpha=False, getslices=False, filemode='hdd', streamdata=None):
        self.log = []
        self.status = True
        self.verbose = verbose
        self.path = path
        self.filemode = filemode
        self.streamdata = streamdata
        if self.verbose:
            print('Initializing ...')
        self.log.append('Initializing')
        try:
            if path.rfind('/')>0:
                self.folder = path[:path.rfind('/')]
                self.filename = path[path.rfind('/')+1:]
            else:
                self.folder = './'
                self.filename = path
            if not getslices:
                self.test_if_ftsfile()
            if self.status:
                if not getslices:
                    self.read_header()
                else: pass
                # get ifg if requested
                if getifg and self.has_block('Data Block IgSm'):
                    self.ifgopd, self.ifg = self.get_datablocks('Data Block IgSm')
                else:
                    self.log.append('No Interferogram requested or not found ... skipping.')
            else: raise(ValueError('Does not seem to be an FTS file ... skipping'))
            if self.verbose and not self.status:
                self.log.append('An error occured.')
                print('An error occured.')
        except Exception as e:
            self.log.append('Problem with '+str(e))
            print('Error while processing '+path+' ... check self.log or do self.print_log()')

class Preview125(QtWidgets.QMainWindow):
    """ A preview of measurements with the IFS125 in idle mode. Similar to the common Check Signal 
    functionality, but with control of the ifg and header data.
    
    ! Before usage: Adjust config.yaml to your situation !"""
    
    def load_yaml(self, yamlfile):
        # load config file
        with open(yamlfile, 'r') as f:
            yamlcontent = yaml.safe_load(f)
        return yamlcontent
     
    def start_measurement(self):
        print('Sending measure command')
        requests.get(self.url_measure)

    def stop_measurement(self):
        print('Sending stop command')
        requests.get(self.url_measurestop)

    def get_status(self):
        # find status info in response to request to ifs
        stat = requests.get(self.stat_htm)
        i1 = stat.text.rfind('ID=MSTCO')
        i2 = stat.text.find('<', i1)
        status = stat.text[i1 + 9:i2]
        return status
    
    def get_preview(self):
        # 
        self.start_measurement()
        status = 'SCN'
        while status != 'IDL':
            # repeat requests until IDL
            status = self.get_status()
        if status=='IDL':
            # find download link
            data = requests.get(self.data_htm)
            i1 = data.text.find('A HREF=')
            i2 = data.text.find('">', i1)
            # download data from ifs
            data = requests.get('/'.join((self.url_ftir,data.text[i1+9:i2])))
            # read in opus format
            self.preview = ftsreader('', verbose=False, getifg=True, filemode='mem', streamdata=data.content)
            self.calc_fft()
        
    def __init__(self):
        # init everyting
        super().__init__()
        # define global variables              
        self.config = self.load_yaml('config.yaml')
        self.site = self.config['selected_site']
        self.siteconfig = self.config[self.site]
        self.url_ftir = 'http://'+self.siteconfig['ip']
        self.url_measure = '/'.join((self.url_ftir, self.siteconfig['preview_commands']))
        self.url_measurestop = '/'.join((self.url_ftir, self.siteconfig['shutdown_commands']))
        self.stat_htm = '/'.join((self.url_ftir,'stat.htm'))
        self.data_htm = '/'.join((self.url_ftir, 'datafile.htm'))
        self.title = '125HR Preview Idle Mode'
        self.setWindowTitle(self.title)
        #
        #
        self._main = QtWidgets.QWidget()
        self.setCentralWidget(self._main)
        layout = QtWidgets.QVBoxLayout(self._main)

        # setup matplotlib canvases
        dynamic_canvas1 = FigureCanvas(Figure(figsize=(5, 3)))
        layout.addWidget(NavigationToolbar(dynamic_canvas1, self))
        layout.addWidget(dynamic_canvas1)
        dynamic_canvas2 = FigureCanvas(Figure(figsize=(5, 3)))
        layout.addWidget(dynamic_canvas2)
        layout.addWidget(NavigationToolbar(dynamic_canvas2, self))
        self._dynamic_ax1 = dynamic_canvas1.figure.subplots()
        self._dynamic_ax1.set_ylim(0,1)
        self._line1, = self._dynamic_ax1.plot(np.arange(4000), [1]*4000, 'b-')
        self._dynamic_ax2 = dynamic_canvas2.figure.subplots()
        self._dynamic_ax1.set_ylim(-1,1)
        self._line2, = self._dynamic_ax2.plot(np.arange(4000), [1]*4000, '-', color='gray')
        self._line3, = self._dynamic_ax2.plot(np.arange(4000), [1]*4000, 'k-')
        # Setup timer to repeat measurement cycle
        self._timer1 = dynamic_canvas1.new_timer(self.config['refreshrate']*1000)
        self._timer1.add_callback(self._update)
        self._timer1.stop()        
        # start button
        self.startButton = QPushButton(self)
        self.startButton.setText('Check Signal')          #text
        self.startButton.setShortcut('Space')  #shortcut key
        self.startButton.clicked.connect(self.startpreview)
        self.startButton.setToolTip('starting timer to perform low-res measurements; Shortcut: [Space]')
        layout.addWidget(self.startButton)
        # stop button
        self.stopButton = QPushButton(self)
        self.stopButton.setText('Stop Measurements')
        self.stopButton.setShortcut('Esc')
        self.stopButton.clicked.connect(self.stoppreview)
        self.stopButton.setToolTip('stop the timer and thus end preview measurements; Shortcut: [Esc]')
        layout.addWidget(self.stopButton)
        
    def startpreview(self):
        print('Check Signal')
        self._timer1.start()
            
    def stoppreview(self):
        print('Stop Measurements')
        self._timer1.stop()
        self.stop_measurement()
        
    def calc_fft(self):
        # get spc via complex fft of ifg
        self.spc = np.fft.fft(self.preview.ifg)
        # calculate wvn axis, LWN is taken from opus header info
        self.wvn = np.fft.fftfreq(int(len(self.spc)),0.5/self.preview.header['Instrument Parameters']['LWN'])[:int(len(self.spc)/2)]
        # determine index of cut-off wavenumber
        l = len(self.wvn[self.wvn<self.config['cutoff']])
        ys = self.spc.copy()
        # set everything in spectrum between larger than cutoff wavenumber to complex 0, same at the end of the array (mirrored spc)
        ys[l:-l] = 0.0+0j
        # calculate inverse fft of low pass filtered spc
        self.ifg_s = np.fft.ifft(ys)
        
    def _update(self):
        # get measurement data
        self.get_preview()
        # update plots
        y = np.abs(self.spc[10:int(len(self.spc)/2)])
        self._line1.set_data((self.wvn[10:], y/np.max(y)))
        self._line1.figure.canvas.draw()
        self._line2.set_data((np.arange(len(self.preview.ifg)), self.preview.ifg/np.max(self.preview.ifg)))
        self._line2.figure.canvas.draw()
        self._line3.set_data((np.arange(len(self.preview.ifg)), self.ifg_s/np.max(self.ifg_s)/2))
        self._line3.figure.canvas.draw()
    
if __name__ == "__main__":
    qapp = QtWidgets.QApplication.instance()
    if not qapp:
        qapp = QtWidgets.QApplication(sys.argv)
    app = Preview125()
    app.show()
    app.activateWindow()
    app.raise_()
    qapp.exec()













