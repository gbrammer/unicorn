#!/usr/bin/env python
"""
GUI tool for inspecting grism spectrum extractions
"""
import sys
import shutil
import os

import Tkinter as tk
from PIL import ImageTk, Image
import numpy as np

noNewLine = '\x1b[1A\x1b[1M'

class TextInput():
    def __init__(self, input='-'):
        master = tk.Toplevel()
        master.geometry('250x30')               # This work fine
        self.master = master
    
        self.textbox = tk.Entry(self.master)
        self.textbox.delete(0, tk.END)
        self.textbox.insert(0, input)
        
        self.button_quit = tk.Button(self.master, text = 'Quit', command = self.finish)
        
        self.textbox.grid_rowconfigure(0, minsize=200)
        self.textbox.grid(row=0, column=0, columnspan=5)
        self.button_quit.grid(row=0, column=5)

        self.master.bind("<Return>", self.finish)

        self.master.mainloop()

    def finish(self, event):
        #self.frame.quit()
        self.text_value = self.textbox.get()
        self.master.destroy()
        self.master.quit()
    
        
class ImageClassifier():
    def __init__(self, images = ['UDS_54826_ijh.png', 'UDS_55031_ijh.png'], logfile='inspect_3dhst.info', load_log=True):
        
        self.logfile = logfile
        if os.path.exists(self.logfile):
            shutil.copy(self.logfile, self.logfile+'.backup')
            print 'Made copy of %s -> %s.backup' %(self.logfile, self.logfile)
            
        self.images = images
        self.N = len(images)
        self.line_flags = np.zeros(self.N, dtype=np.int)
        self.unamb_flags = np.zeros(self.N, dtype=np.int)
        self.star_flags = np.zeros(self.N, dtype=np.int)
        self.contam_flags = np.zeros(self.N, dtype=np.int)
        self.seen_flags = np.zeros(self.N, dtype=np.int)
        self.comments = ['-' for i in range(self.N)]

        if os.path.exists(self.logfile) & load_log:
            self.read_log(self.logfile)
            
        self.i = 0
        
        master = tk.Toplevel()
        master.geometry('1000x550')               # This work fine
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
        spl = self.images[0].split('-')
        if len(spl) == 2:
            rgb_file = spl[0]+'_'+spl[1].split('_')[1]
        else:
            rgb_file = '-'.join(spl[:2])+'_'+spl[2].split('_')[1]
            
        rgb_file = rgb_file.replace('.zfit','_vJH_6')
        if not os.path.exists(rgb_file):
            im_rgb = Image.new('RGB', (100,100), "white")
        else:
            im_rgb = Image.open(rgb_file).resize((150,150))
        
        im_rgb = ImageTk.PhotoImage(im_rgb)        
        self.panel_rgb = tk.Label(self.frame , image=im_rgb)
            
        #self.panel.pack(side=tk.TOP)
        
        self.master.bind("<Key>", self.keypress_event)
        
        #place(relx=0.0, rely=0.0)
        
        self.button_quit = tk.Button(self.frame, text = '(q)uit', command = self.finish)
        #self.button_quit.pack(side=tk.BOTTOM)
        self.button_log = tk.Button(self.frame, text = 'Log to (F)ile', command = self.write_log)
        #self.button_log.pack(side=tk.BOTTOM)
        
        self.button_prev = tk.Button(self.frame, text = '(p)rev', command = self.img_prev)
        #self.button_prev.pack(side=tk.LEFT)
        self.button_next = tk.Button(self.frame, text = '(n)ext', command = self.img_next)
        #self.button_next.pack(side=tk.LEFT)
        
        self.fvar = tk.IntVar()
        self.check_seen = tk.Checkbutton(self.frame, text="Seen?", variable=self.fvar, command=None)
        #self.check_star.var = svar
        #self.check_seen.pack(side=tk.RIGHT)
        self.fvar.set(self.seen_flags[0])

        self.tvar = tk.StringVar()
        self.e_comment = tk.Label(self.frame, textvariable=self.tvar)
        self.e_comment.configure(relief=tk.RIDGE)
        #self.e_comment.grid(row=5, column=0)
        #self.e_comment.delete(0, tk.END)
        #self.e_comment.insert(0, self.comments[0])
        self.tvar.set(self.comments[0])
        
        #self.e_comment.unbind('<Key>')
        
        self.svar = tk.IntVar()
        self.check_star = tk.Checkbutton(self.frame, text="(S)tar?", variable=self.svar, command=self.set_star)
        #self.check_star.var = svar
        #self.check_star.pack(side=tk.RIGHT)
        self.svar.set(self.star_flags[0])

        self.evar = tk.IntVar()
        self.check_edge = tk.Checkbutton(self.frame, text="(B)ad contamination?", variable=self.evar, command=self.set_contam)
        #self.check_star.var = svar
        #self.check_edge.pack(side=tk.RIGHT)
        self.evar.set(self.contam_flags[0])
        
        self.bvar = tk.IntVar()
        self.check_bar = tk.Checkbutton(self.frame, text="(U)nambigous", variable=self.bvar, command=self.set_unamb)
        #self.check_bar.var = dvar
        #self.check_bar.pack(side=tk.RIGHT)
        self.bvar.set(self.unamb_flags[0])

        self.dvar = tk.IntVar()
        self.check_disk = tk.Checkbutton(self.frame, text="Emission (L)ine?", variable=self.dvar, command=self.set_line)
        #self.check_disk.var = self.dvar
        #self.check_disk.pack(side=tk.RIGHT)
        self.dvar.set(self.line_flags[0])
        
        # self.button_next.pack()
        # self.button_prev.pack()
        # self.button_quit.pack()
        # 
        # self.button_next.grid(row=0, column=0)
        # self.button_prev.grid(row=0, column=1)
        # self.button_quit.grid(row=0, column=2)
        
        self.button_next.grid(row=1, column=4, columnspan=1)
        self.button_prev.grid(row=1, column=5, columnspan=1)
        
        self.button_log.grid(row=2, column=4, columnspan=1)
        self.button_quit.grid(row=2, column=5, columnspan=1)
        
        self.check_disk.grid(row=3, column=0)
        self.check_bar.grid(row=3, column=1)
        self.check_edge.grid(row=3, column=2)
        self.check_star.grid(row=3, column=3)
        self.check_seen.grid(row=3, column=4)
        
        self.e_comment.grid(row=4, column=0, columnspan=5)
        
        self.panel.grid(row=0, column=0, columnspan=6)
        self.panel2.grid(row=1, column=1, columnspan=3, rowspan=2)
        self.panel_rgb.grid(row=1, column=0, columnspan=1, rowspan=2)
        
        self.master.mainloop()
    
    def keypress_event(self, event):
        key = event.char
        #print 'Keyboard: "%s"' %(key)
        #self.line_flags[self.i] = (not self.line_flags[self.i])*1
        
        if key == 'l':
            ### Has emission line
            self.line_flags[self.i] = (not self.line_flags[self.i])*1
            self.dvar.set(self.line_flags[self.i])
        
        elif key == 'u':
            ### Unambiguous redshift
            self.unamb_flags[self.i] = (not self.unamb_flags[self.i])*1
            self.bvar.set(self.unamb_flags[self.i])
        
        elif key == 'b':
            ### Bad contamination
            self.contam_flags[self.i] = (not self.contam_flags[self.i])*1
            self.evar.set(self.contam_flags[self.i])
        
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
    
      'x':  skip through through all 'seen'
      '0':  go back to first
      '!':  go to next object with '!' in the comment
      'y':  go to next unambigous
      '?':  this message
"""
    
    def set_line(self):
        #print 'Disk flag: %d' %(self.dvar.get())
        self.line_flags[self.i] = self.dvar.get()
    
    def set_unamb(self):
        #print 'Bar flag: %d' %(self.bvar.get())
        self.unamb_flags[self.i] = self.bvar.get()
    
    def set_star(self):
        #print 'Star flag: %d' %(self.svar.get())
        self.star_flags[self.i] = self.svar.get()
    
    def set_contam(self):
        #print 'Star flag: %d' %(self.svar.get())
        self.contam_flags[self.i] = self.evar.get()
    
    def set_seen(self):
        #print 'Seen flag: %d' %(self.svar.get())
        self.seen_flags[self.i] = 1
    
    def img_next(self):
        self.set_seen()
        self.comments[self.i] = self.tvar.get()

        if self.i == self.N-1:
            print 'Already at last image'
            return False
        
        self.i += 1
        self.load_image()
        
        return True
    #
    def img_prev(self):
        self.set_seen()
        self.comments[self.i] = self.tvar.get() #self.e_comment.get()

        if self.i == 0:
            print 'Already at first image'
            return False
        
        self.i -= 1
        self.load_image()
        return True
        
    def load_image(self):
        
        #print '%d  %s' %(self.i, self.images[self.i])
        #print '%d %d %d' %(self.line_flags[self.i], self.unamb_flags[self.i], self.star_flags[self.i])
        
        print '%s: %d of %d' %(self.images[self.i], self.i+1, self.N)
        
        self.dvar.set(self.line_flags[self.i])
        self.bvar.set(self.unamb_flags[self.i])
        self.svar.set(self.star_flags[self.i])
        self.evar.set(self.contam_flags[self.i])
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
        spl = self.images[self.i].split('-')
        if len(spl) == 2:
            rgb_file = spl[0]+'_'+spl[1].split('_')[1]
        else:
            rgb_file = '-'.join(spl[:2])+'_'+spl[2].split('_')[1]
            
        rgb_file = rgb_file.replace('.zfit','_vJH_6')
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
        
        self.images, dflag, bflag, eflag, sflag, fflag, self.comments = [], [], [], [] ,[], [], []
        for line in lines:
            spl = line.split('"')
            params = spl[0].split()
            #print line, params
            self.images.append(params[0])
            fflag.append(params[1])
            dflag.append(params[2])
            bflag.append(params[3])
            eflag.append(params[4])
            sflag.append(params[5])
            self.comments.append(spl[1])
        #
        self.line_flags = np.cast[int](dflag)
        self.unamb_flags = np.cast[int](bflag)
        self.contam_flags = np.cast[int](eflag)
        self.star_flags = np.cast[int](sflag)
        self.seen_flags = np.cast[int](fflag)
        
    def write_log(self):
        fp = open(self.logfile, 'w')
        # fp.write('#  file            checked  disk  bar  edge star  comment\n')
        fp.write('#  file            checked  emline  unambiguous  bad_contam star  comment\n')

        for i in range(self.N):
            try:
                fp.write('%s   %d   %d   %d   %d   %d  "%s"\n' %(self.images[i], self.seen_flags[i], self.line_flags[i], self.unamb_flags[i], self.contam_flags[i], self.star_flags[i], self.comments[i]))
            except:
                print 'Write failed on line %d (%s)' %(i, self.images[i])
                self.comments[i] = '-'
                fp.write('%s   %d   %d   %d   %d   %d  "%s"\n' %(self.images[i], self.seen_flags[i], self.line_flags[i], self.unamb_flags[i], self.contam_flags[i], self.star_flags[i], self.comments[i]))
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
    