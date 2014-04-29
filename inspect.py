#!/usr/bin/env python
"""
GUI tool for inspecting grism spectrum extractions

$Date$

### Copy the test directory from unicorn to your local machine
cd /tmp/
rsync -avz --progress unicorn.astro.yale.edu:~/GUI_Test /tmp/
cd /tmp/

(in Python):

import glob
#### gui_3dhst.py = unicorn.inspect.py
import gui_3dhst

image_list = glob.glob('*zfit.png')

gui_3dhst.ImageClassifier(images=image_list, logfile='test_inspect.info', RGB_PATH='./', FITS_PATH='./', ds9=None)

"""
import sys
import shutil
import os

import Tkinter as tk
from PIL import ImageTk, Image
import numpy as np

import pyfits

noNewLine = '\x1b[1A\x1b[1M'

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
    
        
class ImageClassifier():
    def __init__(self, images = ['UDS_54826_ijh.png', 'UDS_55031_ijh.png'], logfile='inspect_3dhst.info', RGB_PATH='./', FITS_PATH='./', load_log=True, ds9=None):
        """
        GUI tool for inspecting grism redshift fits
        
         x = inspect.ImageClassifier(images=glob.glob('Specz/GOODS-S-25*zfit.png'), RGB_PATH=RGB_PATH, FITS_PATH='Specz/')

         #### Try on unicorn
         pointing='GOODS-S-25'
         field = '-'.join(pointing.split('-')[:-1])
         PNG_PATH = '/3DHST/Spectra/Release/v4.0/%s/%s-WFC3_v4.0_SPECTRA/%s/ZFIT/PNG/' %(field, field, pointing)
         FITS_PATH = '/3DHST/Spectra/Release/v4.0/%s/%s-WFC3_v4.0_SPECTRA/%s/2D/FITS/' %(field, field, pointing)
         RGB_PATH = '/Users/gbrammer/RGB_v4.0_field/'
         import glob
         from unicorn import inspect
         x = inspect.ImageClassifier(images=glob.glob(PNG_PATH+'*_262*zfit.png'), RGB_PATH=RGB_PATH, FITS_PATH=FITS_PATH, logfile='%s_inspect.info' %(pointing))
         
         
         """
        if len(images) == 0:
            print 'No images specified'
            return False
            
        if not os.path.exists(images[0]):
            print 'First image not found (%s), is path correct?' %(images[0])
            return False
            
        self.logfile = logfile
        if os.path.exists(self.logfile):
            shutil.copy(self.logfile, self.logfile+'.backup')
            print 'Made copy of %s -> %s.backup' %(self.logfile, self.logfile)
            
        self.images = images
        self.N = len(images)
        self.line_flags = np.zeros(self.N, dtype=np.int)
        self.unamb_flags = np.zeros(self.N, dtype=np.int)
        self.misid_flags = np.zeros(self.N, dtype=np.int)
        self.star_flags = np.zeros(self.N, dtype=np.int)
        self.contam_flags = np.zeros(self.N, dtype=np.int)
        self.zp_flags = np.zeros(self.N, dtype=np.int)
        self.seen_flags = np.zeros(self.N, dtype=np.int)
        self.comments = ['---' for i in range(self.N)]
        self.RGB_PATH = RGB_PATH
        self.FITS_PATH = FITS_PATH
        
        self.ds9 = ds9
        if ds9 is True:
            try:
                import pysao
                self.ds9 = pysao.ds9()
            except:
                print 'Couldn\'t import pysao to run DS9'
                self.ds9 = None
            
        if os.path.exists(self.logfile) & load_log:
            self.read_log(self.logfile)
            
        self.i = 0
        
        master = tk.Toplevel()
        master.geometry('1050x590')               # This work fine
        self.master = master
         
        #self.frame = self.master     
        self.frame = tk.Frame(master)
        self.frame.pack()
                
        imageFile = Image.open(self.images[0])
        im = ImageTk.PhotoImage(imageFile)        
        self.panel = tk.Label(self.frame , image=im)

        imageFile2 = Image.open(self.images[0].replace('zfit','zfit.2D')).resize((500,202))
        im2 = ImageTk.PhotoImage(imageFile2)        
        self.panel2 = tk.Label(self.frame , image=im2)
        
        ### RGB
        spl = os.path.basename(self.images[0]).split('-')
        if len(spl) == 2:
            rgb_file = spl[0]+'_'+spl[1].split('_')[1]
        else:
            rgb_file = ''.join(spl[:2])+'_'+spl[2].split('_')[1]
        
        ### GOODS-S-25_22461 -> goodss_24461_vJH_6    
        rgb_file = os.path.join(self.RGB_PATH, rgb_file.lower().split('.zfit')[0] + '_vJH_6.png')
        if not os.path.exists(rgb_file):
            im_rgb = Image.new('RGB', (100,100), "white")
        else:
            im_rgb = Image.open(rgb_file).resize((150,150))
        
        im_rgb = ImageTk.PhotoImage(im_rgb)        
        self.panel_rgb = tk.Label(self.frame , image=im_rgb)
            
        self.master.bind("<Key>", self.keypress_event)
        
        
        ### For logging slider movements
        self.sliders = {}
        
        ######
        ### Navigation buttons
        ###
        self.button_quit = tk.Button(self.frame, text = '(q)uit', command = self.finish)
        self.button_log = tk.Button(self.frame, text = 'Log to (F)ile', command = self.write_log)
        
        self.button_prev = tk.Button(self.frame, text = '(p)rev', command = self.img_prev)
        self.button_next = tk.Button(self.frame, text = '(n)ext', command = self.img_next)
        
        ### Object seen already?
        self.fvar = tk.IntVar()
        self.check_seen = tk.Checkbutton(self.frame, text="Seen?", variable=self.fvar, command=None)
        self.fvar.set(self.seen_flags[0])
        
        ### Comment holder
        self.tvar = tk.StringVar()
        self.e_comment = tk.Label(self.frame, textvariable=self.tvar)
        self.e_comment.configure(relief=tk.RIDGE)
        self.tvar.set(self.comments[0])
        
        #### Star?
        self.svar = tk.IntVar()
        self.check_star = tk.Checkbutton(self.frame, text="(S)tar?", variable=self.svar, command=self.set_star)
        self.svar.set(self.star_flags[0])

        #### Bad z_phot / z_spec or both
        self.zp_var = tk.IntVar()
        self.zp_slider = tk.Scale(self.frame, to=3, resolution=1, orient=tk.HORIZONTAL, variable=self.zp_var, label='Bad (z)p/s/p+s')
        self.zp_var.set(self.zp_flags[0])
        self.sliders[self.zp_slider] = (self.zp_var, self.zp_flags)

        #### Line misidentification
        self.misid_var = tk.IntVar()
        self.check_misid = tk.Checkbutton(self.frame, text="Line (m)isID", variable=self.misid_var, command=self.set_misid)
        self.misid_var.set(self.misid_flags[0])

        #### Unambiguous redshift
        self.bvar = tk.IntVar()
        self.check_unamb = tk.Checkbutton(self.frame, text="(U)nambigous", variable=self.bvar, command=self.set_unamb)
        self.bvar.set(self.unamb_flags[0])
        
        #### Flag contamination
        self.contam_var = tk.IntVar()
        self.contam_slider = tk.Scale(self.frame, to=3, resolution=1, orient=tk.HORIZONTAL, variable=self.contam_var, label='(B)ad Contam:')
        self.contam_var.set(self.contam_flags[0])
        self.sliders[self.contam_slider] = (self.contam_var, self.contam_flags)
        
        #### Emission lne?
        self.line_var = tk.IntVar()
        self.line_slider = tk.Scale(self.frame, to=3, resolution=1, orient=tk.HORIZONTAL, variable=self.line_var, label='Emission (L)ine:')
        self.line_var.set(self.line_flags[0])
        self.sliders[self.line_slider] = (self.line_var, self.line_flags)
        
        #####################        
        ##### Set up the grid for the GUI elements
        self.button_next.grid(row=1, column=4, columnspan=1)
        self.button_prev.grid(row=1, column=5, columnspan=1)
        
        self.button_log.grid(row=2, column=4, columnspan=1)
        self.button_quit.grid(row=2, column=5, columnspan=1)
        
        #self.check_disk.grid(row=3, column=0)
        self.line_slider.grid(row=3, column=0) #, columnspan=1)
        self.contam_slider.grid(row=3, column=1)
        
        self.check_unamb.grid(row=3, column=2)
        
        self.check_misid.grid(row=3, column=3)
        self.zp_slider.grid(row=3, column=4)
        
        self.check_star.grid(row=3, column=5)
        self.check_seen.grid(row=4, column=5)
        
        self.e_comment.grid(row=5, column=0, columnspan=5)
        
        self.panel.grid(row=0, column=0, columnspan=6)
        self.panel2.grid(row=1, column=1, columnspan=3, rowspan=2)
        self.panel_rgb.grid(row=1, column=0, columnspan=1, rowspan=2)
        
        self.master.mainloop()
    
    def log_slider(self, log_var, slider_var, loop=4):
        """
        Log results of moving a slider
        """
        log_var[self.i] = (log_var[self.i] + 1) % loop
        slider_var.set(log_var[self.i])
    
    def listen_slider(self, var):
        print var
        
    def keypress_event(self, event):
        key = event.char
        #print 'Keyboard: "%s"' %(key)
        #self.line_flags[self.i] = (not self.line_flags[self.i])*1
        
        if key == 'l':
            ### Has emission line
            self.log_slider(self.line_flags, self.line_var, loop=4)
            #self.dvar.set((self.line_flags[self.i] > 0)*1)
            #self.line_slider['fg'] = 'red'
            
        elif key == 'u':
            ### Unambiguous redshift
            self.unamb_flags[self.i] = (not self.unamb_flags[self.i])*1
            self.bvar.set(self.unamb_flags[self.i])
        
        elif key == 'm':
            ### Line MisIdentification
            self.misid_flags[self.i] = (not self.misid_flags[self.i])*1
            self.misid_var.set(self.misid_flags[self.i])
        
        elif key == 'b':
            ### Bad contamination
            self.log_slider(self.contam_flags, self.contam_var, loop=4)
            #self.contam_flags[self.i] = (not self.contam_flags[self.i])*1
            #self.evar.set(self.contam_flags[self.i])
        
        
        elif key == 'z':
            ### Bad zphot/spec/both
            self.log_slider(self.zp_flags, self.zp_var, loop=4)
        
        elif key == 's':
            ### object is a star
            self.star_flags[self.i] = (not self.star_flags[self.i])*1
            self.svar.set(self.star_flags[self.i])
        
        elif key == 'c':
            #comment = raw_input('Comment: ')
            ctext = TextInput(self.comments[self.i])
            self.tvar.set(ctext.text_value)
            self.comments[self.i] = ctext.text_value
            #print comment

        elif key == 'd':
            ## Open DS9
            #import pysao
            #ds9 = pysao.ds9()
            if self.ds9 is None:
                print 'Need to pass a pysao.ds9 object at initialization'
                return False
            
            twod_file = self.images[self.i].replace('zfit.png', '2D.fits')
            twod_file = os.path.join(FITS_PATH, os.path.basename(twod_file))
            
            im = pyfits.open(twod_file)
            
            self.ds9.frame(1)
            self.ds9.view(im['DSCI'].data)
            self.ds9.set('scale limits -0.1 1')
            self.ds9.set('cmap value 1 0.173403')

            self.ds9.frame(2)
            self.ds9.view(im['DSEG'].data)
            self.ds9.set('scale limits 0 1')
            
            self.ds9.frame(3)
            self.ds9.view_fits(im['SCI'])
            self.ds9.set('scale limits -0.1 1')
            self.ds9.set('cmap value 8.06667 0.12255')

            self.ds9.frame(4)
            self.ds9.view_fits(im['CONTAM'])
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
            self.write_log()
        
        elif key == 'q':
            ### quit
            self.finish()

        elif key == 'x':
            ### Go to next unseen
            self.set_seen()
            self.comments[self.i] = self.tvar.get()
            
            while (self.i < self.N-1):
                #print self.i
                if self.seen_flags[self.i] == 1:
                    self.i += 1
                else:
                    break

            #self.img_next()
            self.load_image()

        elif key == '0':
            ### Go to beginning
            self.set_seen()
            self.comments[self.i] = self.tvar.get()
            
            self.i = 0
            self.load_image()
                    
        elif key == '!':
            ### Go until find '!' in a comment
            self.set_seen()
            self.comments[self.i] = self.tvar.get()

            ### At end
            if self.i > self.N-2:
                return False
            
            self.i += 1
            
            while (self.i < self.N-1):
                print self.i
                if (self.seen_flags[self.i] == 1) & ('!' not in self.comments[self.i]):
                    self.i += 1
                else:
                    break
            
            self.load_image()
            
            # while (self.seen_flags[self.i] == 1) & ('!' not in self.comments[self.i]):
            #     self.img_next()
        
        elif key == 'y':
            ### Go to next unambiguous
            self.img_next()
            while (self.seen_flags[self.i] == 1) & (self.unamb_flags[self.i] == 0):
                self.img_next()
        
        elif key == '?':
            print """
    Additional keys:

      'c':  Open comment box, type <tab> to edit and <enter> when done.
      'x':  skip through through all 'seen'
      '0':  go back to first
      '!':  go to next object with '!' in the comment
      'y':  go to next unambigous
      '?':  this message
"""
    
    # def set_line(self):
    #     #print 'Disk flag: %d' %(self.dvar.get())
    #     self.line_flags[self.i] = self.dvar.get()
    
    def set_unamb(self):
        #print 'Bar flag: %d' %(self.bvar.get())
        self.unamb_flags[self.i] = self.bvar.get()

    def set_misid(self):
        self.misid_flags[self.i] = self.misid_var.get()
    
    def set_star(self):
        #print 'Star flag: %d' %(self.svar.get())
        self.star_flags[self.i] = self.svar.get()
    
    # def set_contam(self):
    #     #print 'Star flag: %d' %(self.svar.get())
    #     self.contam_flags[self.i] = self.evar.get()
    
    def set_seen(self):
        #print 'Seen flag: %d' %(self.svar.get())
        self.seen_flags[self.i] = 1
    
    def img_next(self):
        self.set_seen()
        self.comments[self.i] = self.tvar.get()
        
        if self.i == self.N-1:
            print 'Already at last image'
            return False
        
        for key in self.sliders.keys():
            slider_var, log_var = self.sliders[key]
            log_var[self.i] = slider_var.get()

        self.i += 1
        self.load_image()
        
        for key in self.sliders.keys():
            slider_var, log_var = self.sliders[key]
            slider_var.set(log_var[self.i])
        
        return True
    #
    def img_prev(self):
        self.set_seen()
        self.comments[self.i] = self.tvar.get() #self.e_comment.get()

        if self.i == 0:
            print 'Already at first image'
            return False
        
        for key in self.sliders.keys():
            slider_var, log_var = self.sliders[key]
            log_var[self.i] = slider_var.get()
        
        self.i -= 1
        self.load_image()
        
        for key in self.sliders.keys():
            slider_var, log_var = self.sliders[key]
            slider_var.set(log_var[self.i])
        
        return True
        
    def load_image(self):
        
        #print '%d  %s' %(self.i, self.images[self.i])
        #print '%d %d %d' %(self.line_flags[self.i], self.unamb_flags[self.i], self.star_flags[self.i])
        
        print '%s: %d of %d' %(self.images[self.i], self.i+1, self.N)
        
        #self.dvar.set(self.line_flags[self.i])
        #self.dvar.set((self.line_flags[self.i] > 0)*1)
        self.line_var.set(self.line_flags[self.i])
        
        self.bvar.set(self.unamb_flags[self.i])
        self.misid_var.set(self.misid_flags[self.i])
        self.svar.set(self.star_flags[self.i])
        #self.evar.set(self.contam_flags[self.i])
        self.fvar.set(self.seen_flags[self.i])
        
        #self.e_comment.delete(0, tk.END)
        #self.e_comment.insert(0, self.comments[self.i])
        self.tvar.set(self.comments[self.i])
        
        imageFile = Image.open(self.images[self.i])
        im = ImageTk.PhotoImage(imageFile)
        self.panel.configure(image = im)
        self.panel.image = im

        imageFile2 = Image.open(self.images[self.i].replace('zfit','zfit.2D'))
        
        im2 = ImageTk.PhotoImage(imageFile2.resize((500,202)))
        self.panel2.configure(image = im2)
        self.panel2.image = im2
        
        ### RGB
        spl = os.path.basename(self.images[0]).split('-')
        if len(spl) == 2:
            rgb_file = spl[0]+'_'+spl[1].split('_')[1]
        else:
            rgb_file = ''.join(spl[:2])+'_'+spl[2].split('_')[1]
        
        ### GOODS-S-25_22461 -> goodss_24461_vJH_6    
        rgb_file = os.path.join(self.RGB_PATH, rgb_file.lower().split('.zfit')[0] + '_vJH_6.png')
        if not os.path.exists(rgb_file):
            im_rgb = Image.new('RGB', (100,100), "white")
        else:
            im_rgb = Image.open(rgb_file).resize((150,150))
        
        im_rgb = ImageTk.PhotoImage(im_rgb)        
        self.panel_rgb.configure(image = im_rgb)
        self.panel_rgb.image = im_rgb
        
        #self.panel_rgb = tk.Label(self.frame , image=im_rgb)
    #
    def read_log(self, logfile):
        self.logfile = logfile
        
        print 'Read logfile: %s' %(self.logfile)
        lines = open(self.logfile).readlines()[1:]
        self.N = len(lines)
        
        self.images, dflag, bflag, eflag, sflag, zflag, mflag, fflag, self.comments = [], [], [], [] ,[], [], [], [], []
        for line in lines:
            spl = line.split('"')
            params = spl[0].split()
            #print line, params
            self.images.append(params[0])
            fflag.append(params[1])
            dflag.append(params[2])
            bflag.append(params[3])
            eflag.append(params[4])
            zflag.append(params[5])
            mflag.append(params[6])
            sflag.append(params[7])
            self.comments.append(spl[1])
        #
        self.line_flags = np.cast[int](dflag)
        self.unamb_flags = np.cast[int](bflag)
        self.contam_flags = np.cast[int](eflag)
        self.misid_flags = np.cast[int](mflag)
        self.zp_flags = np.cast[int](zflag)
        self.star_flags = np.cast[int](sflag)
        self.seen_flags = np.cast[int](fflag)
        
    def write_log(self):
        """
        Write the log file.
        """
        
        ### Catch any updates to the sliders on the current object
        for key in self.sliders.keys():
            slider_var, log_var = self.sliders[key]
            log_var[self.i] = slider_var.get()
        
        fp = open(self.logfile, 'w')
        fp.write('#  file            checked  emline  unambiguous  bad_contam bad_zps line_misid star  comment\n')

        for i in range(self.N):
            try:
                fp.write('%s   %d   %d   %d   %d   %d   %d   %d  "%s"\n' %(self.images[i], self.seen_flags[i], self.line_flags[i], self.unamb_flags[i], self.contam_flags[i], self.zp_flags[i], self.misid_flags[i], self.star_flags[i], self.comments[i]))
            except:
                print 'Write failed on line %d (%s)' %(i, self.images[i])
                self.comments[i] = '-'
                fp.write('%s   %d   %d   %d   %d   %d   %d   %d  "%s"\n' %(self.images[i], self.seen_flags[i], self.line_flags[i], self.unamb_flags[i], self.contam_flags[i], self.zp_flags[i], self.misid_flags[i], self.star_flags[i], self.comments[i]))
                #failed = True
        
        fp.close()
        print 'Wrote %s.' %(self.logfile)
        
    def finish(self):
        #self.frame.quit()
        self.write_log()
        self.master.destroy()
        self.master.quit()
        
if __name__ == "__main__":
    files = sys.argv[1:]
    widget = ImageClassifier(images=files)
    