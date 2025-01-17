# -*- coding: utf-8 -*-
import numpy as np #numerical physics package
import tkinter as tk #GUI package
from tkinter import ttk
import cmath #complex math, needed for structure factor calculations
import random #random number generator, draw new basis vector
#import pmw #Python mega widgets, for tooltips. Is not included in Spyder so can be troublesome to install; optional

def drawall(*args): #redraw all interfaces
    realspacedraw(); recipspacedraw(); settingsdraw()
    return

def addbasis(*args): #event for adding an atom in the base
    global basisvec, rad
    newpos = [[tk.DoubleVar(root,round(random.randrange(2,8)*0.1,1)), tk.DoubleVar(root,round(random.randrange(2,8)*0.1,1))]] #choose random new position
    basisvec = np.append(basisvec,newpos,axis=0) #add basis vector at random position
    rad = np.append(rad,tk.DoubleVar(root,1)) #add radius of new atom, default 1
    drawall()
    return

def removebasis(*args): #event for removing an atom in the base
    global basisvec, rad
    if basisvec.shape[0] > 1: #only remove if more than 1 left
        basisvec = np.delete(basisvec,basisvec.shape[0]-1,0)
        rad = np.delete(rad,rad.shape[0]-1)
    drawall()
    return

def calcrecip(real): #calculate reciprocal lattice vectors
    rot = np.array([[0,-1],[1,0]]) #90 deg rotation matrix
    a1 = np.array([real[0,0].get(),real[0,1].get()]) #real space vector 1
    a2 = np.array([real[1,0].get(),real[1,1].get()]) #real space vector 2
    b1 = 2*np.pi*np.matmul(rot,a2)/np.matmul(a1,np.matmul(rot,a2)) #reciprocal vector 1
    b2 = 2*np.pi*np.matmul(rot,a1)/np.matmul(a2,np.matmul(rot,a1)) #reciprocal vector 2
    recip = np.array([[tk.DoubleVar(root,b1[0]), tk.DoubleVar(root,b1[1])], [tk.DoubleVar(root,b2[0]), tk.DoubleVar(root,b2[1])]]) #store calculation to global variable
    return recip

def realspacedraw(*args): #(re)draw real space lattice
    gridsize = 5 #(half)size of grid crosses
    radscale = 10 #pixel amount per unit radius
    num = 20 #min, max no. of h,k
    mid = real_frame.winfo_reqwidth()/2 #middle coordinates of the plot
    
    #start with an empty canvas
    real_frame.delete('all')
    
    #atom colors
    basis_colors = np.array(['red','blue','green','yellow','purple','orange'])
    
    #draw atoms
    for i  in range(basisvec.shape[0]):
        for x_int in range(-num, num+1):
            for y_int in range(-num, num+1):
                x = mid + scale_real.get()*((x_int+basisvec[i,0].get())*realvec[0,0].get() + (y_int+basisvec[i,1].get())*realvec[1,0].get())
                y = mid - scale_real.get()*((x_int+basisvec[i,0].get())*realvec[0,1].get() + (y_int+basisvec[i,1].get())*realvec[1,1].get()) #flip in y, since tk.Canvas draws with topleft positive
                
                if 0 <= x <= 2*mid and 0 <= y <= 2*mid: #only draw if on screen
                    rad_draw = radscale * rad[i].get()
                    real_frame.create_oval(x-rad_draw, y-rad_draw, x+rad_draw, y+rad_draw, fill=basis_colors[i%basis_colors.size]) #draw an atom, looping through color options
    
    #draw real space square grid markers if enabled
    #seperate loop such that it always overlaps atoms. Could probably be made more efficient
    if show_reallattice.get() == True:
        for x_int in range(-num, num+1):
            for y_int in range(-num, num+1):
                x_cr = mid + x_int*scale_real.get()
                y_cr = mid + y_int*scale_real.get()
                real_frame.create_line((x_cr-gridsize, y_cr, x_cr+gridsize, y_cr), fill='black', width='3')
                real_frame.create_line((x_cr, y_cr-gridsize, x_cr, y_cr+gridsize), fill='black', width='3')
        
    #draw real space lattice vectors if enabled
    if show_realvec.get() == True:
        real_frame.create_line((mid, mid, mid+scale_real.get()*realvec[0,0].get(), mid-scale_real.get()*realvec[0,1].get()), fill='blue', arrow='last', width='2')
        real_frame.create_line((mid, mid, mid+scale_real.get()*realvec[1,0].get(), mid-scale_real.get()*realvec[1,1].get()), fill='blue', arrow='last', width='2')
    
    #draw basis vectors if enabled
    if show_basisvec.get() == True:
        for i in range(1,basisvec.shape[0]):
            real_frame.create_line((mid, mid, mid+scale_real.get()*(basisvec[i,0].get()*realvec[0,0].get()+basisvec[i,1].get()*realvec[1,0].get()), mid-scale_real.get()*(basisvec[i,0].get()*realvec[0,1].get()+basisvec[i,1].get()*realvec[1,1].get())), fill='red', arrow='last', width='2')
    
    #draw miller planes if enabled.
    if show_millerplanes.get() == True:
        if millerindices[0].get() == 0: #to avoid dividing by zero, split the three cases where either miller index is zero
            endpt = 10 * mid * np.array([realvec[0,0].get(), realvec[0,1].get()]) #scale by 10*mid, purely so it is large enough to fall outside the frame
            shift = scale_real.get() * np.array([realvec[1,0].get(), realvec[1,1].get()]) / millerindices[1].get() #shift in between planes
        elif millerindices[1].get() == 0:
            endpt = 10 * mid * np.array([realvec[1,0].get(), realvec[1,1].get()])
            shift = scale_real.get() * np.array([realvec[0,0].get(), realvec[0,1].get()]) / millerindices[0].get()
        else:
            vec_dif = np.array([realvec[0,0].get()/millerindices[0].get() - realvec[1,0].get()/millerindices[1].get(), realvec[0,1].get()/millerindices[0].get() - realvec[1,1].get()/millerindices[1].get()]) #miller vector
            endpt = 10 * mid * vec_dif #scale the difference vector (real vector 1 - recip vector 2)
            shift = scale_real.get() * np.array([realvec[0,0].get() / millerindices[0].get(), realvec[0,1].get() / millerindices[0].get()])
        
        for i in range(-50,51): #try to draw 100 lines
            real_frame.create_line((mid-endpt[0]+i*shift[0], mid+endpt[1]-i*shift[1], mid+endpt[0]+i*shift[0], mid-endpt[1]-i*shift[1]), fill='black', width='2')
    return

def recipspacedraw(*args): #(re)draw reciprocal space
    gridsize = 5 #(half)size of grid crosses
    rad_draw = 5 #drawn radius of diffraction spots
    num = 20 #min, max no. of h,k
    mid = recip_frame.winfo_reqwidth()/2 # middle coodinates of the plot
    
    #start with an empty canvas
    recip_frame.delete('all')
    
    #recalculate reciprocal lattice vectors
    recipvec = calcrecip(realvec)
    
    #draw diffraction spots
    for x_int in range(-num, num+1):
        for y_int in range(-num, num+1):
            x_calc = x_int*recipvec[0,0].get()+y_int*recipvec[1,0].get()
            y_calc = x_int*recipvec[0,1].get()+y_int*recipvec[1,1].get() #flip in y, since tk.Canvas draws with topleft positive
            x = mid + scale_recip.get()*x_calc; y = mid - scale_recip.get()*y_calc
            
            struc = 0 #structure factor summation
            for i in range(0,basisvec.shape[0]):
                if include_ff.get() == True: #if form factors are enabled
                    #how quickly the sinc functions dies out depends on too many physical factors so here it's just made up
                    #I found that an extra factor of 0.05/recip_scale works quite well
                    arg = (0.05 / scale_recip.get()) * np.square(rad[i].get())*(np.square(x_calc)+np.square(y_calc))/(2*np.pi)
                    form = np.sinc(arg/np.pi) #np.sinc is normalized to pi, we don't want that
                else:
                    form = 1
                struc += form*cmath.rect(1,2*np.pi*(x_int*basisvec[i,0].get()+y_int*basisvec[i,1].get()))
                
            inten = np.conj(struc)*struc
            
            #calculate greyscale alpha
            maxinten = np.power(basisvec.shape[0],2) #maximum intensity for given base size
            alpha = int(round(inten.real*100/maxinten,0))
            
            color = 'grey%s'%alpha
            
            if show_inten.get() == 'Values':
                recip_frame.create_text(x,y,text=str(alpha),fill='white')
            else:
                recip_frame.create_oval(x-rad_draw, y-rad_draw, x+rad_draw, y+rad_draw,fill=color,width=0) #draw a spot
            
    #draw reciprocal space square grid markers if enabled
    if show_reciplattice.get() == True:     
        for x_int in range(-num, num+1):
            for y_int in range(-num, num+1):
                x_cr = mid + 2*np.pi*x_int*scale_recip.get()
                y_cr = mid + 2*np.pi*y_int*scale_recip.get()
                recip_frame.create_line((x_cr-gridsize, y_cr, x_cr+gridsize, y_cr), fill='red', width='3')
                recip_frame.create_line((x_cr, y_cr-gridsize, x_cr, y_cr+gridsize), fill='red', width='3')
            
    #draw reciprocal lattice vectors if enabled
    if show_recipvec.get() == True:
        recip_frame.create_line((mid, mid, mid+scale_recip.get()*recipvec[0,0].get(), mid-scale_recip.get()*recipvec[0,1].get()), fill='blue', arrow='last', width='2')
        recip_frame.create_line((mid, mid, mid+scale_recip.get()*recipvec[1,0].get(), mid-scale_recip.get()*recipvec[1,1].get()), fill='blue', arrow='last', width='2')

def settingsdraw(*args): 
    
    #start with empty settings
    for widgets in tabcontrol.winfo_children():
        widgets.destroy()
    
    #define tabs
    tab1 = tk.Frame(tabcontrol); tab2 = tk.Frame(tabcontrol)
    tabcontrol.add(tab1, text='Real space'); tabcontrol.add(tab2, text='Reciprocal space')
    
    #--Real space tab--
    
    #show grid, lattice vectors and basis vectors buttons
    B_show_reallattice = tk.Checkbutton(tab1, text='Show square lattice', variable=show_reallattice, command=realspacedraw)
    B_show_reallattice.grid(sticky='w',column=0, row=0)
    #B_show_reallattice_tt = Pmw.Balloon(root); B_show_reallattice_tt.bind(B_show_reallattice,'Visualize a 1 x 1 real space square lattice')
    
    B_show_realvec = tk.Checkbutton(tab1, text='Show lattice vectors', variable=show_realvec, command=realspacedraw)
    B_show_realvec.grid(sticky='w',column=0, row=1)
    #B_show_realvec_tt = Pmw.Balloon(root); B_show_realvec_tt.bind(B_show_realvec,'Show the real space lattice vectors')
    
    B_show_basisvec = tk.Checkbutton(tab1, text='Show basis vectors', variable=show_basisvec, command=realspacedraw)
    B_show_basisvec.grid(sticky='w',column=0, row=2)
    #B_show_basisvec_tt = Pmw.Balloon(root); B_show_basisvec_tt.bind(B_show_basisvec,'Show the basis vectors')
    
    #show Miller planes and edit Miller indices
    B_show_millerplanes = tk.Checkbutton(tab1, text='Show Miller planes [h,k]', variable=show_millerplanes, command=realspacedraw)
    B_show_millerplanes.grid(sticky='w',column=0, row=3)
    #B_show_millerplanes_tt = Pmw.Balloon(root); B_show_millerplanes_tt.bind(B_show_millerplanes,'Show the Miller planes and edit the Miller indices')
    miller1_entry = tk.Entry(tab1,textvariable=millerindices[0],width=5)
    miller1_entry.grid(column=1, row=3); miller1_entry.bind('<Return>',drawall); miller1_entry.bind('<KP_Enter>',drawall)
    miller2_entry = tk.Entry(tab1,textvariable=millerindices[1],width=5)
    miller2_entry.grid(column=2, row=3); miller2_entry.bind('<Return>',drawall); miller2_entry.bind('<KP_Enter>',drawall)
    
    #edit lattice vector values
    x1y1_lab = tk.Label(tab1,text='Lattice vector 1 [x,y]'); x1y1_lab.grid(sticky='w',column=0,row=4)
    
    x1_entry = tk.Entry(tab1,textvariable=realvec[0,0],width=5)
    x1_entry.grid(column=1, row=4); x1_entry.bind('<Return>',drawall); x1_entry.bind('<KP_Enter>',drawall) #pressing enter pushes the changes to the drawing events
    
    y1_entry = tk.Entry(tab1,textvariable=realvec[0,1],width=5)
    y1_entry.grid(column=2, row=4); y1_entry.bind('<Return>',drawall); y1_entry.bind('<KP_Enter>',drawall)
    
    x1y1_lab = tk.Label(tab1,text='Lattice vector 2 [x,y]'); x1y1_lab.grid(sticky='w',column=0,row=5)
    
    x2_entry = tk.Entry(tab1,textvariable=realvec[1,0],width=5)
    x2_entry.grid(column=1, row=5); x2_entry.bind('<Return>',drawall); x2_entry.bind('<KP_Enter>',drawall)
    
    y2_entry = tk.Entry(tab1,textvariable=realvec[1,1],width=5)
    y2_entry.grid(column=2, row=5); y2_entry.bind('<Return>',drawall); y2_entry.bind('<KP_Enter>',drawall)
    
    #add/remove basis vector buttons
    B_add_basis = tk.Button(tab1, text='Add basis vector', command=addbasis)
    B_add_basis.grid(sticky='w',column=0, row=6)
    #B_add_basis_tt = Pmw.Balloon(root)
    #B_add_basis_tt.bind(B_add_basis,'Add a new basis vector')
    
    B_remove_basis = tk.Button(tab1, text='Remove basis vector', command=removebasis)
    B_remove_basis.grid(sticky='w',column=1, row=6)
   #B_remove_basis_tt = Pmw.Balloon(root)
    #B_remove_basis_tt.bind(B_remove_basis,'Remove the last basis vector')
    
    #first basis vector, coordinates cannot be altered
    row_no = 7 #first row to be used
    
    basis_lab = tk.Label(tab1,text='Basis vector 1 [x,y,radius]'); basis_lab.grid(sticky='w',column=0,row=row_no)
    basis_x_entry = tk.Label(tab1,text='0'); basis_x_entry.grid(column=1, row=row_no)
    basis_y_entry = tk.Label(tab1,text='0'); basis_y_entry.grid(column=2, row=row_no)
    
    basis_rad_entry = tk.Entry(tab1,textvariable=rad[0],width=5)
    basis_rad_entry.grid(column=3, row=row_no)
    basis_rad_entry.bind('<Return>',drawall); basis_rad_entry.bind('<KP_Enter>',drawall)
    
    #add entries for every basis vector except the first
    for i in range(1,basisvec[:,0].size):
            basis_no = i+1
            
            basis_lab = tk.Label(tab1,text='Basis vector %s [x,y,radius]'%basis_no)
            basis_lab.grid(sticky='w',column=0,row=row_no+i)
            
            basis_x_entry = tk.Entry(tab1,textvariable=basisvec[i,0],width=5)
            basis_x_entry.grid(column=1, row=row_no+i)
            basis_x_entry.bind('<Return>',drawall); basis_x_entry.bind('<KP_Enter>',drawall) #pressing enter pushes the changes to the drawing events
            
            basis_y_entry = tk.Entry(tab1,textvariable=basisvec[i,1],width=5)
            basis_y_entry.grid(column=2, row=row_no+i)
            basis_y_entry.bind('<Return>',drawall); basis_y_entry.bind('<KP_Enter>',drawall)
            
            basis_rad_entry = tk.Entry(tab1,textvariable=rad[i],width=5)
            basis_rad_entry.grid(column=3, row=row_no+i)
            basis_rad_entry.bind('<Return>',drawall); basis_rad_entry.bind('<KP_Enter>',drawall)
        
    #--Reciprocal space tab--
    
    #show grid and show lattice vectors buttons
    B_show_reciplattice = tk.Checkbutton(tab2, text='Show square lattice', variable=show_reciplattice, command=recipspacedraw)
    B_show_reciplattice.grid(sticky='w', column=0, row=0)
    #B_show_reciplattice_tt = Pmw.Balloon(root)
    #B_show_reciplattice_tt.bind(B_show_reciplattice,'Visualize a 2π x 2π reciprocal space square lattice')
    
    B_show_recipvec = tk.Checkbutton(tab2, text='Show lattice vectors', variable=show_recipvec, command=recipspacedraw)
    B_show_recipvec.grid(sticky='w', column=0, row=1)
    #B_show_recipvec_tt = Pmw.Balloon(root)
    #B_show_recipvec_tt.bind(B_show_recipvec,'Show the reciprocal space lattice vectors')
    
    B_inten_label = tk.Label(tab2,text='Intensity mode')
    B_inten_label.grid(sticky='w',column=0,row=2)
    B_show_inten = tk.OptionMenu(tab2,show_inten,*['Spots','Values'])
    B_show_inten.grid(sticky='w', column=1, row=2)
    #B_show_inten_tt = Pmw.Balloon(root)
    #B_show_inten_tt.bind(B_show_inten,'Represent diffraction intensities as spots or values')
    
    B_include_ff = tk.Checkbutton(tab2, text='Include form factors', variable=include_ff, command=recipspacedraw)
    B_include_ff.grid(sticky='w', column=0, row=3)
    #B_include_ff_tt = Pmw.Balloon(root)
    #B_include_ff_tt.bind(B_include_ff,'Include form factors by taking the atom radius into account')
    
    return

def main(): #main function to loop
    global root
    root=tk.Tk() #define root
    root.title("BRAVAIS") #set title
    root.geometry('1360x506') #set geometry of main frame
    root.configure(bg='grey90')
    #Pmw.initialise(root) #Pmw needed for tooltips
    
    #global toggles
    global show_reallattice, show_realvec, show_basisvec, show_millerplanes, show_reciplattice, show_recipvec, include_ff, show_inten
    show_reallattice = tk.BooleanVar(root,False) #is square lattice shown?
    show_realvec = tk.BooleanVar(root,True) #are real lattice vectors shown?
    show_basisvec = tk.BooleanVar(root,True) #are basis vectors shown?
    show_millerplanes = tk.BooleanVar(root,False) #show the Miller planes
    show_reciplattice = tk.BooleanVar(root,False)
    show_recipvec = tk.BooleanVar(root,True)
    include_ff = tk.BooleanVar(root,True) #calculate form factors from atom radius?
    show_inten = tk.StringVar(root,'Spots'); show_inten.trace('w',recipspacedraw) #show intensities as value or diffraction spot?
    
    #global variables
    global realvec, basisvec, recipvec, millerindices, scale_real, scale_recip, rad
    realvec = np.array([[tk.DoubleVar(root,1),tk.DoubleVar(root,0)] ,[tk.DoubleVar(root,0),tk.DoubleVar(root,1)]]) #real lattice vectors coords [[x1,y1],[x2,y2]]
    basisvec = np.array([[tk.DoubleVar(root,0),tk.DoubleVar(root,0)]]) #real basis vector coords [[0,0],[x2,y2], ... ]
    millerindices = np.array([tk.DoubleVar(root,1),tk.DoubleVar(root,1)]) #indices of the miller planes
    scale_real = tk.DoubleVar(root,75) #scale of plots in real space, pixels per unit distance
    scale_recip = tk.DoubleVar(root,25/(2*np.pi)) #and in reciprocal space
    rad = np.array([tk.DoubleVar(root,1)]) #radii of basis atoms
    
    #initialize frames
    global real_frame, recip_frame, tabcontrol, tab1, tab2
    real_frame = tk.Canvas(root, width=500, height=500, borderwidth=2, relief='solid',bg='grey90'); real_frame.grid(column=0, row=0)
    tabcontrol = ttk.Notebook(root); tabcontrol.grid(column=1, row=0,sticky='nw')
    recip_frame = tk.Canvas(root, width=500, height=500, borderwidth=2, relief='solid', bg='black'); recip_frame.grid(column=2, row=0)
    
    drawall()
    
    root.mainloop()
    
if __name__ == "__main__": #keep looping main() until program is closed
    main()
