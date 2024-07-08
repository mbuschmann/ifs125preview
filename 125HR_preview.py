# -*- coding: utf-8 -*-
"""
Standalone QT program to provide live preview of the idle mode measurements of an IFS125HR.

The program includes rudimentary FFT functionality and provides a smoothed (= low pass - filtered) interferogram 
to detect potential detector nonlinearity during alignment.

Author: Matthias Buschmann, IUP Bremen
Date: 2023/06/01
"""

from __future__ import print_function, division
import sys, yaml, requests, struct, io
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.qt_compat import QtWidgets
#from PyQt5.QtWidgets import QPushButton, QLabel, QFileDialog
from matplotlib.backends.backend_qtagg import FigureCanvas
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure



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
        if block == 'Data Block IgSm' or block == 'Data Block':
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

def load_yaml(yamlfile):
    # load config file
    with open(yamlfile, 'r') as f:
        yamlcontent = yaml.safe_load(f)
    return yamlcontent

def smooth_ifg(o, lwn=15798.022, cutoff=3700, l0=4000, verbose='no'):
    ''' Taking file from spectrometer, either via commandline input or as data stream from ifs125. Then the following is done:
        - remove offset of ifg by substracting median
        - calculate scale of ifg to apply to smoothed ifg later
        - apodize ifg by hanning window
        - do complex fft
        - set everything in spectrum larger than cuttoff wavenumber to complex 0 by multipying hanning window
        - calculate inverse fft (-> smoothed ifg) and return with scale and offset re-applied.
    different verbosity levels to return different steps of this program for introspection.
    '''
    pkl = o.header['Instrument Parameters']['PKL']
    # zero ifg
    ifg0 = o.ifg[int(pkl-l0/2):int(pkl+l0/2)]
    offset = np.median(o.ifg[int(pkl-l0/20):int(pkl+l0/20)])
    ifgz = ifg0-offset
    # scale the smoothed ifg later with the diff between max and min of original ifg
    scale = 1/(np.max(ifgz[int(l0/2-l0/20):int(l0/2+l0/20)])-np.min(ifgz[int(l0/2-l0/20):int(l0/2+l0/20)]))
    p = 5 # percent apodization region at beginning and end of IFG
    l = int(l0*p/100.0)
    # create hanning apodization function and apply to ifg
    a1 = np.ones(l0)
    a1[:l] = ((np.cos(np.pi*np.arange(l)/l)+1)**2/4)[::-1]
    a1[-l:] = ((np.cos(np.pi*np.arange(l)/l)+1)**2/4)
    ifga = ifgz*a1
    # get spc via complex fft of ifg
    spc = np.fft.fft(ifga)
    # calculate wvn axis, default LWN is taken from opus header info
    #wvn = np.fft.fftfreq(int(len(spc)), 0.5*lwn)[:int(len(spc)/2)]
    #print(wvn.shape, spc.shape)
    wvn = np.fft.fftfreq(int(len(spc)),0.5/lwn)[:int(len(spc)/2)]
    # determine index of cut-off wavenumber
    l = len(wvn[wvn<cutoff])
    ys = spc.copy()
    # set everything in spectrum between larger than cutoff wavenumber to complex 0, same at the end of the array (mirrored spc)
    #ys[l:-l] = 0.0+0j
    # define and apply Hann window function
    sfunc = lambda nu, cutoff: np.cos(np.pi*nu/(2*cutoff))**2
    a2 = np.ones(len(ys))*(0+0j)
    a2[:l] = sfunc(wvn[:l], cutoff)
    # apply to mirrored part of spc as well. careful to use the correct order of wvn here
    a2[-l:] = sfunc(wvn[:l][::-1], cutoff)
    # apply apodization
    ys = a2*ys
    # calculate inverse fft of apodized spc, discarding imaginary part of reverse fft
    if verbose=='ifg':
        return np.fft.ifft(ys).real*scale+offset, a1, ifga, ifgz, ifg0, offset, scale
    elif verbose=='spc':
        return ys, wvn, spc, a2
    else:
        return np.fft.ifft(ys).real*scale+offset, ys, a1, wvn, offset

def calc_dip_from_fit(ifgs, fitwindowsize=200, zpdblock=40, return_fits=False):
    '''Select a subset of points around zpd and fit a linear function to it. 
    The immidiate region around zpd is blocked out as to not fit the dip itself (zpd +- zpdplock/2 points).'''
    from scipy.optimize import curve_fit
    fitfunc = lambda x, a, b: a*x+b
    y = ifg_s
    x = np.arange(len(y))
    center = int(len(y)/2)
    sel1 = (x>center-fitwindowsize/2) & (x<center+fitwindowsize/2)
    x0, y0 = x[sel1], y[sel1]
    sel2 = (x<center-zpdblock/2) & (x>center+zpdblock/2)
    sel = sel1 | sel2
    x1, y1 = x[sel], y[sel]
    popt, pcov = curve_fit(fitfunc, x1, y1, p0=[1.0,1.0])
    x2, y2 = x1, y0-fitfunc(x1, *popt)
    if np.max(y2)>=np.max(np.abs(y2)):
        # 'positive' dip
        dip = np.max(y2)
    else:
        # negative dip
        dip = np.min(y2)
    print(f"DIP amplitude from fit: {dip:.5f} ‰")
    if return_fits:
        return dip, x1, fitfunc(x1, *popt)
    else:
        return dip

def calc_dip_from_minmax(ifgs, offset, fitwindowsize=200):
    ''' Calculate the minimum or maximum of the smoothed ifg +-fitwindowsize points around zpd.'''
    print('\n\nsome error in the dip calculation using minmax... do not trust\n')
    sifgs = ifgs[int(len(ifgs/2)-fitwindowsize/2):int(len(ifgs/2)+fitwindowsize/2)]
    sifgs = sifgs-offset
    if np.abs(np.max(sifgs))>=np.abs(np.min(sifgs)):
        dip = np.max(sifgs)
    else:
        dip = np.min(sifgs)
    print(f"DIP amplitude from minmax: {dip:.5f} ‰")
    return dip

class Preview125(QtWidgets.QMainWindow):
    """ A preview of measurements with the IFS125 in idle mode. Similar to the common Check Signal 
    functionality, but with control of the ifg and header data.
    
    ! Before usage: Adjust config.yaml to your situation !"""
     
    def start_measurement(self):
        if self.running:
            print('Sending measure command')
            self.label.setText("Sending measure command")
            self.label.setStyleSheet("background-color: lightgreen")
            requests.get(self.url_measure)
        else: pass

    def stop_measurement(self):
        print('Stopping measurements')
        self.running = False
        self._timer1.stop()
        self.label.setText("Measurements stopped.")
        self.label.setStyleSheet("background-color: orange")

    def shutdown_measurement(self):
        print('Sending shutdown command')
        requests.get(self.url_measureshutdown)
        self.label.setText("Sent shutdown measure command")
        self.label.setStyleSheet("background-color: orange")

    def save_data_opus(self):
        file , check = QtWidgets.QFileDialog.getSaveFileName(None, "Save last recorded interferogram as OPUS file", "", "All Files (*)")
        if check:
            print('Saving data')
            with open(file, 'wb') as f:
                f.write(self.raw_ifg)
        else:
            print('Not saving anything. Still running = ', self.running)

    def save_data(self):
        file , check = QtWidgets.QFileDialog.getSaveFileName(None, "Save smoothed interferogram as ASCII file", "", "All Files (*)")
        if check and not self.running:
            print('Saving data in ', file)
            with open(file, 'w') as f:
                s = ''
                for l in self.ifg_s:
                    s+=str(l)+'\n'
                f.write(s)
        else: 
            print('Not saving anything. Still running = ', self.running)
        #requests.get(self.url_measureshutdown)
        #self.label.setText("Sent shutdown measure command")
        #self.label.setStyleSheet("background-color: orange")

    def get_status(self):
        # find status info in response to request to ifs
        stat = requests.get(self.stat_htm)
        i1 = stat.text.rfind('ID=MSTCO')
        i2 = stat.text.find('<', i1)
        status = stat.text[i1 + 9:i2]
        return status
    
    def get_preview(self):
        if self.running:
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
                self.raw_ifg = data.content
                # read in opus format
                self.preview = ftsreader('', verbose=False, getifg=True, filemode='mem', streamdata=self.raw_ifg)
                # preselect ifg around zpd
                self.zpd()
                self.ifg = self.preview.ifg[int(self.zpdindex-self.npt/2):int(self.zpdindex+self.npt/2)]
                # get smoothed ifg
                self.ifg_s, self.spc_apodized, self.apo, self.apo_wvn, offset = smooth_ifg(self.preview, lwn=self.config['lwn'],  cutoff=self.config['cutoff'], l0=self.npt)
                #self.calc_spc()
                #print('all 0? ', np.all(self.ifg_s==0))
                # calc spc
                self.spc = np.fft.fft(self.ifg)
                self.wvn = np.fft.fftfreq(int(len(self.spc)),0.5/self.preview.header['Instrument Parameters']['LWN'])[:int(len(self.spc)/2)]
                # calculate dip my min max of smoothed ifg
                #dip1 = calc_dip_from_minmax(self.ifg_s, offset)
                # calculate dip by substracting a fitted linear function around zpd
                dip2 = calc_dip_from_fit(self.ifg_s, fitwindowsize=self.config['fitwindowsize'], zpdblock=self.config['blockzpd'])
                #self.label.setText(f"DIP amplitude (minmax): {dip1:.5f} ‰     ")
                self.label.setText(f"DIP amplitude (fit): {dip2:.5f} ‰")
                self.label.setStyleSheet("background-color: white")
            else: pass
        else:
            pass
        
    def startpreview(self):
        print('Check Signal')
        self.running=True
        self._timer1.start()
            
    def stoppreview(self):
        print('Stop Measurements')
        self.running=False
        self._timer1.stop()
        self.stop_measurement()
    
    def calc_spc(self):
        self.spc = np.fft.fft(self.preview.ifg)
        self.wvn = np.fft.fftfreq(int(len(self.spc)),0.5/self.preview.header['Instrument Parameters']['LWN'])[:int(len(self.spc)/2)]
        self.ifg_s = self.preview.ifg
    
    def zpd(self):
        # use peak location from header
        self.zpdindex = self.preview.header['Instrument Parameters']['PKL']

    def zpd_minmax(self):
        # using mean between indices of max and min values of ifg
        self.zpdindex = int(round(np.mean([np.argmin(self.preview.ifg), np.argmax(self.preview.ifg)]))) 

    def _update(self):
        # get measurement data
        self.get_preview()
        self.run +=1
        # update plots
        y = np.abs(self.spc[10:int(len(self.spc)/2)])
        self._line1.set_data(self.wvn[10:], y)
        #if self.run == 2:
        #    self._dynamic_ax2.set_xlim(self.config['spc_plot_xlim'])
        #    self._dynamic_ax2.set_ylim(np.min(y)*1.2, np.max(y)*1.2)
        self._line1.figure.canvas.draw()
        #if self.scaledifgaxes:
        #    y1 = self.ifg-np.mean(self.ifg)
        #    self._line2.set_data((np.arange(len(self.ifg)), y1/np.max(np.abs(y1))))
        #else:
        self._line2.set_data((np.arange(len(self.ifg)), self.ifg))            
        #if self.run == 2:
        #    #self._dynamic_ax2.set_xlim(0,4000)
        #    self._dynamic_ax2.set_ylim(np.min(self.preview.ifg)*1.2, np.max(self.preview.ifg)*1.2)
        self._line2.figure.canvas.draw()
        #if self.scaledifgaxes:
        #    y2 = self.ifg_s-np.mean(self.ifg_s[self.config['zpd_interval'][0]:self.config['zpd_interval'][1]])
        #    self._line3.set_data((np.arange(len(self.preview.ifg)), y2/np.max(np.abs(y2[self.config['zpd_interval'][0]:self.config['zpd_interval'][1]]))))
        #else:
        self._line3.set_data((np.arange(len(self.ifg_s)), self.ifg_s))
        #if self.run == 2:
        #    #self._dynamic_ax2.set_xlim(0,4000)
        #    self._dynamic_ax3.set_ylim(np.min(self.ifg_s)*1.2, np.max(self.ifg_s)*1.2)
        self._line3.figure.canvas.draw()
        
    #def clickBox(self, b):
    #    if b.isChecked() == True:
    #        self.scaledifgaxes = True
    #        print('IFG y-axis scaling: ON')
    #    else:
    #        self.scaledifgaxes = False
    #        print('IFG y-axis scaling: OFF')
        
    def __init__(self):
        # init everyting
        super().__init__()
        # define global variables              
        self.config = load_yaml('config.yaml')
        self.run = 0
        self.running=False
        self.npt = self.config['npt']
        self.site = self.config['selected_site']
        self.siteconfig = self.config[self.site]
        self.url_ftir = 'http://'+self.siteconfig['ip']
        self.url_measure = '/'.join((self.url_ftir, self.siteconfig['preview_commands']))
        self.url_measureshutdown = '/'.join((self.url_ftir, self.siteconfig['shutdown_commands']))
        self.stat_htm = '/'.join((self.url_ftir,'stat.htm'))
        self.data_htm = '/'.join((self.url_ftir, 'datafile.htm'))
        self.title = '125HR Preview Idle Mode'
        self.setWindowTitle(self.title)
        self.resize(800, 700)
        #
        #
        self._main = QtWidgets.QWidget()
        self.setCentralWidget(self._main)
        layout = QtWidgets.QVBoxLayout(self._main)
        # setup label
        self.label = QtWidgets.QLabel('Press Check Signal to start', self)
        self.label.setStyleSheet("background-color: orange")
        layout.addWidget(self.label)
        # setup matplotlib canvases
        dynamic_canvas1 = FigureCanvas(Figure(figsize=(5, 3)))
        layout.addWidget(NavigationToolbar(dynamic_canvas1, self))
        layout.addWidget(dynamic_canvas1)
        #
        #self.box = QCheckBox('IFG y-axis scaled',self)
        #self.box.setChecked(False)
        self.scaledifgaxes = False
        #self.box.stateChanged.connect(lambda : self.clickBox(self.box))
        #layout.addWidget(self.box)
        #
        dynamic_canvas2 = FigureCanvas(Figure(figsize=(5, 3)))
        layout.addWidget(dynamic_canvas2)
        layout.addWidget(NavigationToolbar(dynamic_canvas2, self))
        self._dynamic_ax1 = dynamic_canvas1.figure.subplots()
        self._dynamic_ax1.set_title('Spectrum preview')
        self._dynamic_ax1.set_xlim(self.config['spc_plot_xlim'])
        self._dynamic_ax1.set_ylim(0,1)
        self._line1, = self._dynamic_ax1.plot(np.linspace(3000, 11000, self.npt), np.zeros(self.npt), 'b-')
        self._dynamic_ax2 = dynamic_canvas2.figure.subplots()
        self._dynamic_ax2.set_title('Raw and smoothed Interferogram preview')
        self._dynamic_ax3 = self._dynamic_ax2.twinx()
        #self._dynamic_ax3.set_title('Smoothed Interferogram preview')
        #self._dynamic_ax2.set_xlim(self.config['zpd_interval'])
        #self._dynamic_ax2.set_ylim(-1,1)
        self._dynamic_ax2.set_ylabel('Original IFG [Bruker units]')
        self._dynamic_ax3.set_ylabel('Smoothed IFG / (max - min)(IFG) * 1e3]')
        self._line2, = self._dynamic_ax2.plot(np.arange(self.npt), np.zeros(self.npt), '-', color='grey')
        self._line3, = self._dynamic_ax3.plot(np.arange(self.npt), np.zeros(self.npt), 'k-')
        # Setup timer to repeat measurement cycle
        self._timer1 = dynamic_canvas1.new_timer(self.config['refreshrate']*1000)
        self._timer1.add_callback(self._update)
        self._timer1.stop()        
        # start button
        self.startButton = QtWidgets.QPushButton(self)
        self.startButton.setText('Check Signal')          #text
        self.startButton.setShortcut('Space')  #shortcut key
        self.startButton.clicked.connect(self.startpreview)
        self.startButton.setToolTip('starting timer to perform low-res measurements; Shortcut: [Space]')
        layout.addWidget(self.startButton)
        # stop button
        self.stopButton = QtWidgets.QPushButton(self)
        self.stopButton.setText('Stop Measurements')
        self.stopButton.setShortcut('Esc')
        self.stopButton.clicked.connect(self.stoppreview)
        self.stopButton.setToolTip('stop the timer and thus end preview measurements; Shortcut: [Esc]')
        layout.addWidget(self.stopButton)
        # save spectra
        self.saveButton = QtWidgets.QPushButton(self)
        self.saveButton.setText('Save smoothed data (stop meas. first!)')
        self.saveButton.setShortcut('Alt+S')
        self.saveButton.clicked.connect(self.save_data)
        self.saveButton.setToolTip('Save current IFG and smoothed spectrum to file; Shortcut: [Alt+S]')
        layout.addWidget(self.saveButton)
        self.saveButtono = QtWidgets.QPushButton(self)
        self.saveButtono.setText('Save interferogram (stop meas. first!)')
        self.saveButtono.clicked.connect(self.save_data_opus)
        self.saveButtono.setToolTip('Save current IFG and smoothed spectrum to file; Shortcut: [Alt+S]')
        layout.addWidget(self.saveButtono)
        # stop button
        self.shutdownButton = QtWidgets.QPushButton(self)
        self.shutdownButton.setText('Send shutdown command (e.g. lamp off)')
        self.shutdownButton.clicked.connect(self.shutdown_measurement)
        self.shutdownButton.setToolTip('Send one last shudown command, defined in config.yaml')
        layout.addWidget(self.shutdownButton)

if __name__ == "__main__":
    if len(sys.argv)==1:
        qapp = QtWidgets.QApplication.instance()
        if not qapp:
            qapp = QtWidgets.QApplication(sys.argv)
        app = Preview125()
        app.show()
        app.activateWindow()
        app.raise_()
        qapp.exec()
    else:
        fname = sys.argv[1]
        config = load_yaml('config.yaml')
        o = ftsreader(fname, getifg=True)
        cutoff = config['cutoff']
        lwn = config['lwn']
        npt = config['npt']
        l0 = npt

        ifg_s, a1, ifga, ifgz, ifg0, offset, scale = smooth_ifg(o, lwn=lwn, cutoff=cutoff, l0=npt, verbose='ifg')

        dip2, xfit, yfit = calc_dip_from_fit(ifg_s, fitwindowsize=200, zpdblock=40, return_fits=True)
        
        # plot ifg
        fig0, ax01 = plt.subplots()
        #ax01.set_title(f"DIP amplitude (minmax): {dip1:.5f} ‰     ")
        ax01.set_title(f"DIP amplitude (fit): {dip2:.5f} ‰")
        ax01.plot(ifg0, label='original ifg')
        ax01.plot(ifg_s, label='smoothed ifg')
        ax01.plot(xfit, yfit, label='background fit')

        # plot spc
        spc = np.abs(np.fft.fft(ifg0))
        wvn = np.fft.fftfreq(int(len(spc)),0.5/lwn)[:int(len(spc)/2)]
        fig, ax21 = plt.subplots(nrows=1)
        ax21.plot(wvn[20:], spc[20:int(len(spc)/2)])
        
        plt.show()










