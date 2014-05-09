#!/usr/bin/env python
"""
GUI tool for inspecting grism spectrum extractions

$Date$
$Rev$

Instructions for working:
- Put gui_3dhst.py in the directory where you plan to work. This is where the output
files will be created. I recommend TEST_SPECTRA_v4.1/

- Open gui_3dhst.py in your favorite text editor.

- Go to line XXX below and edit the paths for RELEASE and RGB_PATH to the correct
location on your disk. Save.

- Open a teminal window. 

- If running Ureka, type "ur_setup"

- Navigate to gui_3dhst.py directory.

Open python (type python or ipython in your terminal).
In Python:

    import gui_3dhst

	### The default behavior will pop out the ds9 window:
	gui_3dhst.go_pointing(pointing = 'goodss-05')
	
	### To avoid opening the ds9 window:
	gui_3dhst.go_pointing(pointing = 'goodss-05',show_ds9=False)
	
	     
"""
import sys
import shutil
import os
import glob

import Tkinter as tk
from PIL import ImageTk, Image

import numpy as np
import pyfits

noNewLine = '\x1b[1A\x1b[1M'

def examples():
    
    import glob
    #import unicorn.inspect
    import pysao
    
    os.chdir(os.getenv('THREEDHST') + '/3DHST_VariableBackgrounds/GS25')
    
    RGB_PATH = os.getenv('RELEASE') + 'v4.0/RGB/All/'
    
    ds9 = pysao.ds9()
    
    #### Objects with fits
    x = unicorn.inspect.ImageClassifier(glob.glob('ZFIT/PNG/*zfit.png'), FITS_PATH='2D/FITS/', logfile='inspect_zfit', RGB_PATH=RGB_PATH, ds9=ds9)
    
    #### Raw 2D files, including objects not fit
    x = unicorn.inspect.ImageClassifier(glob.glob('2D/PNG/*png'), FITS_PATH='2D/FITS/', logfile='inspect_raw2d', RGB_PATH=RGB_PATH, ds9=ds9)
    
    
class TextInput():
    """
    Popup for entering a comment
    """
    def __init__(self, input='---'):
        master = tk.Toplevel()
        master.geometry('250x30') # This work fine
        self.master = master
    
        self.textbox = tk.Entry(self.master)
        self.textbox.delete(0, tk.END)
        self.textbox.insert(0, input)
        
        #self.button_quit = tk.Button(self.master, text = 'Quit', command = self.finish)
        
        self.textbox.grid_rowconfigure(0, minsize=240)
        self.textbox.grid(row=0, column=0, columnspan=5)
        #self.button_quit.grid(row=0, column=5)

        self.master.bind("<Return>", self.finish)

        self.master.mainloop()

    def finish(self, event):
        #self.frame.quit()
        self.text_value = self.textbox.get()
        self.master.destroy()
        self.master.quit()
    
def go_pointing(pointing='goodss-25', RELEASE='/Volumes/Voyager/TEST_SPECTRA_v4.1/GOODS-S/INTERLACE_v4.1/', RGB_PATH='/Volumes/Voyager/TEST_SPECTRA_v4.1/RGB_goodss/',show_ds9=True):
    """
    
    Classify a 3D-HST pointing
        
    """
    import glob
    #from unicorn import inspect

    field = '-'.join(pointing.split('-')[:-1])
    
    ### Example paths for v4.0 and earlier releases
    #PNG_PATH = '%s/%s/%s-WFC3_v4.0_SPECTRA/%s/ZFIT/PNG/' %(RELEASE, field, field, pointing)
    #FITS_PATH = '%s/%s/%s-WFC3_v4.0_SPECTRA/%s/2D/FITS/' %(RELEASE, field, field, pointing)

    #PNG_PATH = '%s/%s/WFC3/%s/ZFIT/PNG/' %(RELEASE, field, pointing)
    #FITS_PATH = '%s/%s/WFC3/%s/2D/FITS/' %(RELEASE, field, pointing)
    
    PNG_PATH = FITS_PATH = RELEASE
    if show_ds9:
    	import pysao
    	ds9=pysao.ds9()
    	x = ImageClassifier(images=glob.glob(PNG_PATH+pointing+'*new_zfit.png'), RGB_PATH=RGB_PATH, FITS_PATH=FITS_PATH, logfile='%s_inspect.info' %(pointing), ds9=ds9)
    else:
    	x = ImageClassifier(images=glob.glob(PNG_PATH+pointing+'*new_zfit.png'), RGB_PATH=RGB_PATH, FITS_PATH=FITS_PATH, logfile='%s_inspect.info' %(pointing))
    return x
    
class myCheckbutton(object):
    """
    Container for CheckButton self + logging variables
    """
    def __init__(self, gui, text="(t)ilt", param=[0], hotkey='xxx'):
        self.gui = gui
        self.param = param
        self.hotkey = hotkey

        self.var = tk.IntVar()
         
        self.thebutton = tk.Checkbutton(self.gui.frame, text=text, variable=self.var, command=self.set_var)
        self.var.set(self.param[0])
        
    def set_var(self):
        self.param[self.gui.i] = self.var.get()
        #print self.gui.i
      
class mySlider(object):
    """
    Container for Scale slider self + logging variables
    """
    def __init__(self, gui, text="(a)bs./break", param=[0], hotkey='xxx', to=3):
        self.gui = gui
        self.param = param
        self.hotkey = hotkey
        self.to = to        
        
        self.var = tk.IntVar()
        
        self.theslider = tk.Scale(self.gui.frame, to=to, resolution=1, orient=tk.HORIZONTAL, variable=self.var, label=text)
        self.var.set(param[0])
        
    def set_var(self):
        self.param[self.gui.i] = self.var.get()
    
    
class ImageClassifier():
    """
    Main classifier tool for 3D-HST fits
    """
    def __init__(self, images = ['UDS_54826.zfit.png', 'UDS_55031.zfit.png'], logfile='inspect_3dhst.info', RGB_PATH='./', FITS_PATH='./', load_log=True, ds9=None):
        """
        GUI tool for inspecting grism redshift fits
        
         x = unicorn.inspect.ImageClassifier(images=glob.glob('Specz/GOODS-S-25*zfit.png'), RGB_PATH=RGB_PATH, FITS_PATH='Specz/')
                  
         """
        if len(images) == 0:
            print 'No images specified'
            return False
            
        if not os.path.exists(images[0]):
            print 'First image not found (%s), is path correct?' %(images[0])
            return False
        
        self.RGB_PATH = RGB_PATH
        self.FITS_PATH = FITS_PATH
        
        #### Check ds9
        self.ds9 = ds9
        if ds9 is True:
            try:
                import pysao
                self.ds9 = pysao.ds9()
            except:
                print 'Couldn\'t import pysao to run DS9'
                self.ds9 = None

        ##### Add .fits to filename and make backup if necessary
        self.logfile = logfile
        if not self.logfile.lower().endswith('.fits'):
            self.logfile += '.fits'
        
        if os.path.exists(self.logfile):
            bk = glob.glob(self.logfile+'.backup*')
            if len(bk) > 0:
                bkup_file = self.logfile + '.backup.%03d' %(len(bk))
            else:
                bkup_file = self.logfile + '.backup'
                
            shutil.copy(self.logfile, bkup_file)
            print 'Made copy of %s -> %s' %(self.logfile, bkup_file)
        
        ####### Initialize parameters
        self.params = {}        
        self.images = images
        
        if os.path.exists(self.logfile) & load_log:
            self.read_fits()
            
        self.N = len(self.images)

        for key in ['line', 'extended', 'absorption', 'unamb', 'misid', 'contam', 'zp', 'tilt', 'investigate', 'star', 'seen', 'bad2d', 'deblend', 'sed']:
            if key not in self.params.keys():
                self.params[key] = np.zeros(self.N, dtype=np.int)
        
        if 'comment' not in self.params.keys():
            self.params['comment'] = ['---' for i in range(self.N)]
                                        
        self.i = 0
        
        ####### Initialize GUI
        master = tk.Toplevel()
        simple_twod = False
        if 'zfit.png' in self.images[0]:
            master.geometry('1050x630')              
        elif '2D.png' in self.images[0]:
            master.geometry('850x750')     
            simple_twod = True
        else:
            master.geometry('1050x630') 
            
        self.master = master
         
        self.frame = tk.Frame(master)
        self.frame.pack()
        
        #### Image Panels
        imageFile = Image.open(self.images[0])
        im = ImageTk.PhotoImage(imageFile)        
        self.panel = tk.Label(self.frame , image=im)
        
        # imageFile2 = Image.open(self.images[0].replace('zfit','zfit.2D')).resize((500,202))
        # im2 = ImageTk.PhotoImage(imageFile2)        
        im2 = ImageTk.PhotoImage(self.get_twod_file(self.images[0]))
        self.panel2 = tk.Label(self.frame , image=im2)
        # self.panel2.configure(image = im2)
        # self.panel2.image = im2
        
        #### RGB Panel        
        im_rgb = ImageTk.PhotoImage(self.get_rgb_file(self.images[0])) 
        self.panel_rgb = tk.Label(self.frame , image=im_rgb)
        
        #### Keypress binding
        self.master.bind("<Key>", self.keypress_event)
        
        ### For logging slider movements
        #self.sliders = {}
        
        self.sliders = {}
        
        ######
        ### Navigation buttons
        ###
        self.button_quit = tk.Button(self.frame, text = '(q)uit', command = self.finish)
        self.button_log = tk.Button(self.frame, text = 'Log to (F)ile', command = self.write_fits)
        
        self.button_prev = tk.Button(self.frame, text = '(p)rev', command = self.img_prev)
        self.button_next = tk.Button(self.frame, text = '(n)ext', command = self.img_next)
        
        self.buttons = {}
        
        #### Emission lne?
        self.sliders['line'] = mySlider(self, text='Emission (l)ine', param=self.params['line'], hotkey='l')
        
        #### Absorption
        self.sliders['absorption'] = mySlider(self, text='(a)bs./break', param=self.params['absorption'], hotkey='a')

        #### Extended line
        self.sliders['extended'] = mySlider(self, text='(e)xtended', param=self.params['extended'], hotkey='e')
        
        #### Unambiguous redshift
        self.buttons['unamb'] = myCheckbutton(self, text='(u)nambiguous', param=self.params['unamb'], hotkey='u')
        
        #### Flag contamination
        self.sliders['contam'] = mySlider(self, text='(b)ad contam', param=self.params['contam'], hotkey='b')

        #### Line misidentification
        self.buttons['misid'] = myCheckbutton(self, text='Line (m)isID', param=self.params['misid'], hotkey='m')
                
        #### Tilt?
        self.buttons['tilt'] = myCheckbutton(self, text='(t)ilt', param=self.params['tilt'], hotkey='t')
                
        #### Bad z_phot / z_spec or both
        self.buttons['zp'] = myCheckbutton(self, text='Bad (z)spec', param=self.params['zp'], hotkey='z')

        #### SED problems
        self.sliders['sed'] = mySlider(self, text='Bad SED (k)', param=self.params['sed'], hotkey='k', to=2)

        #### Flag for manual investigation
        self.sliders['investigate'] = mySlider(self, text='(i)nvestigate', param=self.params['investigate'], hotkey='i', to=4)
        
        ### Object seen already?
        self.buttons['seen'] = myCheckbutton(self, text='SEEN?', param=self.params['seen'], hotkey='xxx')
        
        #### Star?
        self.buttons['star'] = myCheckbutton(self, text='(s)tar', param=self.params['star'], hotkey='s')

        #### Spectrum problem
        self.buttons['bad2d'] = myCheckbutton(self, text='Bad 1D/2D (x)', param=self.params['bad2d'], hotkey='x')

        #### Deblend issue
        self.buttons['deblend'] = myCheckbutton(self, text='(d)eblend', param=self.params['deblend'], hotkey='d')
        
        ### Comment holder
        self.tvar = tk.StringVar()
        self.e_comment = tk.Label(self.frame, textvariable=self.tvar)
        self.e_comment.configure(relief=tk.RIDGE)
        self.tvar.set(self.params['comment'][0])        
        
        #####################        
        ##### Set up the grid for the GUI elements
        self.button_next.grid(row=1, column=4, columnspan=1)
        self.button_prev.grid(row=1, column=5, columnspan=1)
        
        self.button_log.grid(row=2, column=4, columnspan=1)
        self.button_quit.grid(row=2, column=5, columnspan=1)

        self.buttons['seen'].thebutton.grid(row=3, column=5)
        self.buttons['star'].thebutton.grid(row=3, column=4)
        
        self.sliders['line'].theslider.grid(row=4, column=0)
        self.sliders['absorption'].theslider.grid(row=5, column=0)
        
        self.sliders['extended'].theslider.grid(row=4, column=1)
        self.buttons['unamb'].thebutton.grid(row=5, column=1)

        self.sliders['contam'].theslider.grid(row=4, column=2)
        self.buttons['misid'].thebutton.grid(row=5, column=2)

        self.buttons['tilt'].thebutton.grid(row=4, column=3)        
        self.buttons['zp'].thebutton.grid(row=5, column=3)

        self.sliders['sed'].theslider.grid(row=5, column=4)
        
        self.sliders['investigate'].theslider.grid(row=4, column=5)
        self.buttons['bad2d'].thebutton.grid(row=4, column=4)

        self.buttons['deblend'].thebutton.grid(row=5, column=5)
        
        
        self.e_comment.grid(row=6, column=0, columnspan=5)
        
        if simple_twod:
            self.panel.grid(row=0, column=1, rowspan=4, columnspan=3)
        else:
            self.panel.grid(row=0, column=0, columnspan=6)
            
        self.panel2.grid(row=1, column=1, columnspan=3, rowspan=3)
        self.panel_rgb.grid(row=1, column=0, columnspan=1, rowspan=3)
        
        self.master.mainloop()
    
    def get_twod_file(self, image_file):
        """
        Get the zfit.2D file if appropriate
        """
        if ('2D.png' in image_file) | ('_stack.png' in image_file):
            imageFile2 = Image.new('RGB', (1,1), "white")
        else:
            imageFile2 = Image.open(image_file.replace('zfit','zfit.2D')).resize((500,202))
        
        return imageFile2
        
        #im2 = ImageTk.PhotoImage(imageFile2.resize((500,202)))
        #self.panel2.configure(image = im2)
        #self.panel2.image = im2
        
    def get_rgb_file(self, image_file):
        """
        Translate the filename to an RGB thumbnail
        """
        spl = os.path.basename(image_file.replace('new_zfit', 'newzfit')).split('-')
        if len(spl) == 2:
            rgb_file = spl[0]+'_'+spl[1].split('_')[1]
        else:
            rgb_file = ''.join(spl[:2])+'_'+spl[2].split('_')[1]
        
        ### GOODS-S-25_22461 -> goodss_24461_vJH_6    
        rgb_file = os.path.join(self.RGB_PATH, rgb_file.lower().split('.newzfit')[0].split('.zfit')[0].split('_stack')[0].split('.2d.png')[0] + '_vJH_6.png')
        print rgb_file
        
        if not os.path.exists(rgb_file):
            im_rgb = Image.new('RGB', (100,100), "white")
        else:
            im_rgb = Image.open(rgb_file).resize((150,150))
        
        return im_rgb
        
    def keypress_event(self, event):
        key = event.char
        #print 'Keyboard: "%s"' %(key)
        #self.line_flags[self.i] = (not self.params['line'][self.i])*1
        
        #### Buttons
        for id in self.buttons.keys():
            if key == self.buttons[id].hotkey:
                button = self.buttons[id]
                button.param[self.i] = (not button.param[self.i])*1
                button.var.set(button.param[self.i])
                #print button.param[self.i]
                return True
        
        #### Sliders
        for id in self.sliders.keys():
            if key == self.sliders[id].hotkey:
                slider = self.sliders[id]
                slider.param[self.i] = ((slider.param[self.i]+1) % (slider.to+1))
                slider.var.set(slider.param[self.i])
                #print slider.param[self.i]
                return True

        #### Other keys
        if key == 'c':
            #comment = raw_input('Comment: ')
            ctext = TextInput(self.params['comment'][self.i])
            self.tvar.set(ctext.text_value)
            self.params['comment'][self.i] = ctext.text_value
            #print comment

        elif key == '9':
            #### Open DS9
            #import pysao
            #ds9 = pysao.ds9()
            if self.ds9 is None:
                print 'Need to pass a pysao.ds9 object at initialization'
                return False
            
            exts = ['new_zfit.png','_stack.png','2D.png']
            twod_file = '%s' %(self.images[self.i])
            for ext in exts:
                twod_file = twod_file.replace(ext, '2D.fits')
                
            #twod_file = self.images[self.i].replace('zfit.png', '2D.fits')
            twod_file = os.path.join(self.FITS_PATH, os.path.basename(twod_file))
            
            im = pyfits.open(twod_file)
            
            self.ds9.frame(1)
            self.ds9.view(im['DSCI'].data)
            self.ds9.set('scale limits -0.1 1')
            self.ds9.set('cmap value 1 0.173403')

            self.ds9.frame(2)
            self.ds9.view(im['DSEG'].data)
            self.ds9.set('scale limits 0 1')
            
            self.ds9.frame(3)
            self.ds9.view(im['SCI'].data)
            self.ds9.set('scale limits -0.1 1')
            self.ds9.set('cmap value 8.06667 0.12255')

            self.ds9.frame(4)
            self.ds9.view(im['CONTAM'].data)
            self.ds9.set('scale limits -0.1 1')
            self.ds9.set('cmap value 8.06667 0.12255')
            
            self.ds9.set('tile yes')
            
        elif key == 'n':
            ### next
            self.img_next()
        
        elif key == 'p':
            ### previous
            self.img_prev()
        
        elif key == 'f':
            ### Force write the output file
            self.write_fits()
        
        elif key == 'q':
            ### quit
            self.finish()

        elif key == 'N':
            ### Go to next unseen
            self.set_seen()
            self.params['comment'][self.i] = self.tvar.get()
            
            while (self.i < self.N-1):
                #print self.i
                if self.params['seen'][self.i] == 1:
                    self.i += 1
                else:
                    break

            #self.img_next()
            self.load_image()

        elif key == '0':
            ### Go to beginning
            self.set_seen()
            self.params['comment'][self.i] = self.tvar.get()
            
            self.i = 0
            self.load_image()
                    
        elif key == '!':
            ### Go until find '!' in a comment
            self.set_seen()
            self.params['comment'][self.i] = self.tvar.get()

            ### At end
            if self.i > self.N-2:
                return False
            
            self.i += 1
            
            while (self.i < self.N-1):
                print self.i
                if (self.params['seen'][self.i] == 1) & ('!' not in self.params['comment'][self.i]):
                    self.i += 1
                else:
                    break
            
            self.load_image()
            
            # while (self.params['seen'][self.i] == 1) & ('!' not in self.params['comment'][self.i]):
            #     self.img_next()
        
        elif key in ['U','L','E','I','A','S','T','B']:
            ### Go to next non-zero flag for different keys
            key_param = {'U':self.params['unamb'], 'L':self.params['line'], 'E':self.params['extended'], 'I':self.params['investigate'], 'A':self.params['absorption'], 'S':self.params['star'], 'T':self.params['tilt'], 'B':self.params['contam']}
            self.next_nonzero(key_param[key])
            
        elif key == '?':
            print """
    Additional keys:

      'c':  Open comment box, type <tab> to edit and <enter> when done.
      '9':  Open 2D FITS extensions in DS9, if pysao available.
      'N':  skip through through all 'seen'
      '0':  go back to first
      '!':  go to next object with '!' in the comment
      'y':  go to next unambigous
      '?':  this message
      
"""
        else:
            print 'Hotkey (%s) not bound.' %(key)
    
    def next_nonzero(self, param=None):
        """
        Go to next object where self.param is nonzero
        """
        self.set_seen()
        self.params['comment'][self.i] = self.tvar.get()

        ### At end
        if self.i > self.N-2:
            return False
        
        self.i += 1
        
        while (self.i < self.N-1):
            print self.i
            if (self.params['seen'][self.i] == 1) & (param[self.i] == 0):
                self.i += 1
            else:
                break
        
        self.load_image()
            
    def set_seen(self):
        #print 'Seen flag: %d' %(self.svar.get())
        self.params['seen'][self.i] = 1
    
    def img_next(self):
        self.set_seen()
        self.params['comment'][self.i] = self.tvar.get()
        
        if self.i == self.N-1:
            print 'Already at last image'
            return False
                    
        self.i += 1
        self.load_image()
                    
        return True
    #
    def img_prev(self):
        self.set_seen()
        self.params['comment'][self.i] = self.tvar.get() #self.e_comment.get()

        if self.i == 0:
            print 'Already at first image'
            return False
                    
        self.i -= 1
        self.load_image()
        
        return True
        
    def load_image(self):
        
        #print '%d  %s' %(self.i, self.images[self.i])
        #print '%d %d %d' %(self.params['line'][self.i], self.params['unamb'][self.i], self.params['star'][self.i])
        
        print '%s: %d of %d' %(self.images[self.i], self.i+1, self.N)
        
        for key in self.buttons.keys():
            button = self.buttons[key]
            button.var.set(button.param[self.i])
        
        for key in self.sliders.keys():
            slider = self.sliders[key]
            slider.var.set(slider.param[self.i])
        
        #self.e_comment.delete(0, tk.END)
        #self.e_comment.insert(0, self.params['comment'][self.i])
        self.tvar.set(self.params['comment'][self.i])
        
        imageFile = Image.open(self.images[self.i])
        im = ImageTk.PhotoImage(imageFile)
        self.panel.configure(image = im)
        self.panel.image = im

        #imageFile2 = Image.open(self.images[self.i].replace('zfit','zfit.2D'))
        im2 = ImageTk.PhotoImage(self.get_twod_file(self.images[self.i]))
        self.panel2.configure(image = im2)
        self.panel2.image = im2
        
        ### RGB
        im_rgb = ImageTk.PhotoImage(self.get_rgb_file(self.images[self.i])) 
        self.panel_rgb.configure(image = im_rgb)
        self.panel_rgb.image = im_rgb
        
        #self.panel_rgb = tk.Label(self.frame , image=im_rgb)
        
    def write_fits(self):
        """
        Write the FITS log
        """
        
        import time
        import getpass
        
        formats = {}
        formats['bool'] = 'L'
        formats['int16'] = 'I'
        formats['int32'] = 'J'
        formats['int64'] = 'K'
        formats['float32'] = 'E'
        formats['float64'] = 'D'
        
        formats['>i8'] = 'K'
        formats['>f8'] = 'D'
        
        #### Make the table columns, translating numpy data types to "TFORM"
        coldefs = []
        TFORM = 'A'+str(np.array(self.images).dtype).split('S')[1]
        
        coldefs.append(pyfits.Column(name='images', array=np.array(self.images), format=TFORM))
        
        for column in self.params.keys():
            if column == 'comment':
                coldata = np.array(self.params['comment'])
            else:
                coldata = self.params[column]
            #
            dtype = str(coldata.dtype)
            #print column, dtype
            if dtype in formats.keys():
                TFORM=formats[dtype]
            else:
                if 'S' not in dtype:
                    print 'Unrecognized data type in: %s' %(dtype)
                    return False
                #
                TFORM = 'A'+dtype.split('S')[1]
            #
            #data = self.params[column]
            if '>' in dtype:
                cast_types = {'>i8':np.int64, '>f8':np.float64}
                coldata = np.cast[cast_types[dtype]](coldata)
            #
            coldefs.append(pyfits.Column(name=column, array=coldata, format=TFORM))
        
        #### Done, now make the binary table
        tbhdu = pyfits.new_table(coldefs)

        #### Primary HDU
        hdu = pyfits.PrimaryHDU()
        thdulist = pyfits.HDUList([hdu,tbhdu])

        #### Add modification time of "infile" to FITS header
        infile_mod_time = time.strftime("%m/%d/%Y %I:%M:%S %p",
                            time.localtime()) # os.path.getmtime(self.filename)))
        
        thdulist[0].header.update('MODTIME', infile_mod_time)
        thdulist[0].header.update('USER', getpass.getuser())

        thdulist.writeto(self.logfile, clobber=True)
        
        print 'Log to file %s' %(self.logfile)
        
    def read_fits(self):
        """
        Read already saved output
        """
        
        im = pyfits.open(self.logfile)#
        tab = im[1].data
        print "Read log %s from %s (%s)" %(self.logfile, im[0].header['USER'], im[0].header['MODTIME'])
        
        colnames = tab.columns.names
        
        try:
            #### If images in the input list aren't in the file that 
            #### was read, add them to the table
            from astropy.table import Table
            tab = Table.read(self.logfile)
            colnames = tab.columns
            read_images = tab['images']
            for image in self.images:
                if image not in read_images:
                    tab.add_row()
                    tab['images'][-1] = image
                    tab['comment'][-1] = '---'
            
        except:
            #### No astropy?
            print 'No astropy found.  Forcing image list from %s.' %(self.logfile)
            pass
            
        self.images = tab['images']
            
        for c in colnames:
            #print c
            if c == 'images':
                continue

            if c == 'comment':
                self.params[c] = list(tab[c])
                continue
            #    
            self.params[c] = tab[c]
            
    def finish(self):
        #self.frame.quit()
        #self.write_log()
        self.write_fits()
        self.master.destroy()
        self.master.quit()
        
if __name__ == "__main__":
    files = sys.argv[1:]
    widget = ImageClassifier(images=files)
    