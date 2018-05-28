'''
/*
 *             CAPLIB - CAPsid LIBrary (beta-version)
 *
 *  For Calculations on Icosahedrally Symmetric Virus Capsid Structures
 *
 *  The source codes reported in the following article is contained in this directory,
 *  Shigetaka Yoneda, Yukina Hara-Yamada, Aya Kosugi, Maiko Nanao, Takami Saito, Shunsuke Sato, Nozomu Yamada, and Go Watanabe,
 *  "CAPLIB: A New Program Library for the Modeling and Analysis of Icosahedrally Symmetric Viral Capsids",
 *  ACS Omega 2018, 3, 4458−4465.
 *
 *  The purpose of this library is to analyze directions of rotation axes, calculate cell numbers, 
 *  generate the entire structure of capsid from protomer struture, etc.
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, version 2.
 *
 *  Febrary, 2018.
 */
CAPAXIS.py

CAPAXIS.py is composed of Capgen, Caprad, Caprot, Capwat, Capswch, Captwst, get_obj_center and setaxis.

Main functions
Capgen: Generation the entire structures
Caprad: Enlargement of radius
Caprot: Rotation each monomer
Capwat: Addition water molecules
Capswch: Switching around arbitrary axis
Captwst: Twisting around arbitrary axis

Sub functions
get_obj_center: Calculate center of each P-object.
setaxis: Set axis of rotation for Capswch and Captwst.
'''
from pymol import cmd, stored, selector
from ctypes import *
from math import cos, sin, radians, sqrt, acos
import numpy

capsid_center_xyz = [0.,0.,0.]
axis = [0.,0.,0.]
obj_center_init = numpy.zeros([60,3])

def get_obj_center():
    for i in range(60):
        exec('stored.pos%s = []'%i)

    objname = cmd.get_names()
    n = len(objname)
    s = range(n-1)
    for i in s:
        cmd.iterate_state(1,'OBJ%s'%i,'stored.pos%s.append([x,y,z])'%i)

    n2 = len(stored.pos0)
    s2 = range(n2)
    pos_tmp = numpy.zeros([60,n2,3])
    obj_center_xyz = numpy.zeros([60,3])
    for i in range(60):
        for j in s2:
            exec('pos_tmp[i][j][0] = stored.pos%s[j][0]'%i)
            exec('pos_tmp[i][j][1] = stored.pos%s[j][1]'%i)
            exec('pos_tmp[i][j][2] = stored.pos%s[j][2]'%i)
            obj_center_xyz[i][0] = obj_center_xyz[i][0] + pos_tmp[i][j][0]
            obj_center_xyz[i][1] = obj_center_xyz[i][1] + pos_tmp[i][j][1]
            obj_center_xyz[i][2] = obj_center_xyz[i][2] + pos_tmp[i][j][2]
        obj_center_xyz[i][0] = obj_center_xyz[i][0] / n2
        obj_center_xyz[i][1] = obj_center_xyz[i][1] / n2
        obj_center_xyz[i][2] = obj_center_xyz[i][2] / n2
    return obj_center_xyz

'''
Capgen(shift) 

shift = 0: Center of capsid is not moved (default).
shift != 0: Center of capsid is moved to [0,0,0].
'''
def capgen(shift = 0):
    stored.pos = []
    stored.pos2 = []
    stored.j = 0
    rot = c_double * 600; rot_array = rot()
    xyzcent_t = c_double * 200; xyzcent_t_array = xyzcent_t()
    biomt = c_int * 1; nbiomt_t = biomt(); nbiomt_t = 0
    global obj_center_init
    stored.pos_all_obj = []
    stored.k = 0
    sum_x = 0.; sum_y = 0.; sum_z = 0.
    global capsid_center_xyz

    objname = cmd.get_names()
    s = range(60)
    for i in s:
        cmd.copy('OBJ%s'%i,objname[0])

    cmd.iterate_state(1,'OBJ0','stored.pos.append([x,y,z])')
    n = len(stored.pos)
    s2 = range(n)
    stored.pos2 = list(stored.pos)

    cdll.LoadLibrary("/Users/yoneda/Desktop/170131SatoMasterCAPAXIS/C/capaxis_original/src_bak/capgen5.so")
    libc = CDLL("/Users/yoneda/Desktop/170131SatoMasterCAPAXIS/C/capaxis_original/src_bak/capgen5.so")
    file_name = c_char * 20; file_name1 = file_name()
    file_name1 = objname[0] + ".pdb"
    print file_name1

    libc.read_biomt_file(c_char_p(file_name1), byref(rot_array), byref(xyzcent_t_array))
    x = c_double * n; x1 = x(); x2 = x()
    y = c_double * n; y1 = y(); y2 = y()
    z = c_double * n; z1 = z(); z2 = z()
    for i in s2:
        x1[i] = stored.pos2[i][0]
        y1[i] = stored.pos2[i][1]
        z1[i] = stored.pos2[i][2]
        x2[i] = 0.
        y2[i] = 0.
        z2[i] = 0.
    
    for i in s:
        libc.gen_entire(byref(rot_array),byref(xyzcent_t_array),byref(x1),byref(y1),byref(z1),byref(x2),byref(y2),byref(z2),c_int(n),c_int(i))
        for k in s2:
            stored.pos2[k][0] = x2[k]
            stored.pos2[k][1] = y2[k]
            stored.pos2[k][2] = z2[k]
        stored.j = 0
        cmd.alter_state(1,'OBJ%s'%i,'x = stored.pos2[stored.j][0];y = stored.pos2[stored.j][1];z = stored.pos2[stored.j][2];stored.j=stored.j+1')

    objname = cmd.get_names()
    n = len(objname)
    s = range(n-1)
    for i in s:
        cmd.iterate_state(1,'OBJ%s'%i,'stored.pos_all_obj.append([x,y,z])')
    n1 = len(stored.pos_all_obj)
    s1 = range(n1)
    for i in s1:
        sum_x = sum_x + stored.pos_all_obj[i][0]
        sum_y = sum_y + stored.pos_all_obj[i][1]
        sum_z = sum_z + stored.pos_all_obj[i][2]
    center_x = sum_x / n1
    center_y = sum_y / n1
    center_z = sum_z / n1
    capsid_center_xyz = [center_x,center_y,center_z]

    if shift != 0:
        for i in range(60):
            cmd.translate ([-capsid_center_xyz[0],-capsid_center_xyz[1],-capsid_center_xyz[2]],'OBJ%s'%i)

    obj_center_init = get_obj_center()

'''
Caprad(mag,target)

mag: Magnification.
target: Name of operational P-object (default = all).
'''
def caprad(mag,target=all):
    if target == all:
        obj_center_xyz = numpy.zeros([60,3])
        relative_vect_objcenter = numpy.zeros([60,3])
        relative_vect_objcenter_mod = numpy.zeros([60,3])
        d_xyz = numpy.zeros([60,3])
        
        obj_center_xyz = get_obj_center()

        for i in range(60):
            relative_vect_objcenter[i][0] = obj_center_xyz[i][0]
            relative_vect_objcenter[i][1] = obj_center_xyz[i][1]
            relative_vect_objcenter[i][2] = obj_center_xyz[i][2]
            relative_vect_objcenter_mod[i][0] = mag * relative_vect_objcenter[i][0]
            relative_vect_objcenter_mod[i][1] = mag * relative_vect_objcenter[i][1]
            relative_vect_objcenter_mod[i][2] = mag * relative_vect_objcenter[i][2]
            d_xyz[i][0] = relative_vect_objcenter_mod[i][0] - relative_vect_objcenter[i][0]
            d_xyz[i][1] = relative_vect_objcenter_mod[i][1] - relative_vect_objcenter[i][1]
            d_xyz[i][2] = relative_vect_objcenter_mod[i][2] - relative_vect_objcenter[i][2]
        for i in range(60):
            cmd.translate ([d_xyz[i][0],d_xyz[i][1],d_xyz[i][2]],'OBJ%s'%i)

    if target != all:
        stored.pos_target = []
        obj_center_xyz = numpy.zeros([3])
        relative_vect_objcenter = numpy.zeros([3])
        relative_vect_objcenter_mod = numpy.zeros([3])
        d_xyz = numpy.zeros([3])
        cmd.iterate_state(1,target,'stored.pos_target.append([x,y,z])')
        n = len(stored.pos_target)
        s = range(n)
        for i in s:
            obj_center_xyz[0] = obj_center_xyz[0] + stored.pos_target[i][0]
            obj_center_xyz[1] = obj_center_xyz[1] + stored.pos_target[i][1]
            obj_center_xyz[2] = obj_center_xyz[2] + stored.pos_target[i][2]
        obj_center_xyz[0] = obj_center_xyz[0] / n
        obj_center_xyz[1] = obj_center_xyz[1] / n
        obj_center_xyz[2] = obj_center_xyz[2] / n
        relative_vect_objcenter[0] = obj_center_xyz[0]
        relative_vect_objcenter[1] = obj_center_xyz[1]
        relative_vect_objcenter[2] = obj_center_xyz[2]
        relative_vect_objcenter_mod[0] = mag * relative_vect_objcenter[0]
        relative_vect_objcenter_mod[1] = mag * relative_vect_objcenter[1]
        relative_vect_objcenter_mod[2] = mag * relative_vect_objcenter[2]
        d_xyz[0] = relative_vect_objcenter_mod[0] - relative_vect_objcenter[0]
        d_xyz[1] = relative_vect_objcenter_mod[1] - relative_vect_objcenter[1]
        d_xyz[2] = relative_vect_objcenter_mod[2] - relative_vect_objcenter[2]
        cmd.translate ([d_xyz[0],d_xyz[1],d_xyz[2]],target)

'''
Caprot(theta,target,axis)

theta: Angle of rotation.
target: Name of operational P-object (default = all).
axis: Axis of rotation (default = None: axis is set to center of P-object).
'''
def caprot(theta,target=all,axis=None):
    cmd.reset()
    if target == all:
        obj_center_xyz = numpy.zeros([60,3])

        obj_center_xyz = get_obj_center()
        
        cmd.origin (position=[0.,0.,0.])
        for i in range(60):
            cmd.rotate ([obj_center_xyz[i][0],obj_center_xyz[i][1],obj_center_xyz[i][2]], theta, 'OBJ%s'%i)
    
    if target != all:
        obj_center_xyz = numpy.zeros([3])
        stored.pos_target = []
        cmd.iterate_state(1,target,'stored.pos_target.append([x,y,z])')
        n = len(stored.pos_target)
        s = range(n)
        for i in s:
            obj_center_xyz[0] = obj_center_xyz[0] + stored.pos_target[i][0]
            obj_center_xyz[1] = obj_center_xyz[1] + stored.pos_target[i][1]
            obj_center_xyz[2] = obj_center_xyz[2] + stored.pos_target[i][2]
        obj_center_xyz[0] = obj_center_xyz[0] / n
        obj_center_xyz[1] = obj_center_xyz[1] / n
        obj_center_xyz[2] = obj_center_xyz[2] / n
        
        if axis != None:
            obj_center_xyz = axis
        
        cmd.origin (position=[0.,0.,0.])
        cmd.rotate ([obj_center_xyz[0],obj_center_xyz[1],obj_center_xyz[2]], theta, target)

'''
Capwat(infile,wradiu1,wradiu2,fnamex)

infile: Pass to input file.
wradiu1: Outer radius of water.
wradiu2: Inner radius of water.
fnamex: Name of output file.
'''
def capwat(infile,wradiu1,wradiu2,fnamex):
    addmodule = cdll.LoadLibrary("temp26_4.so")
    addmodule.pdbtop_.argtypes = [POINTER(c_char),POINTER(c_float),POINTER(c_float),POINTER(c_char),POINTER(c_int),POINTER(c_int)]
    addmodule.pdbtop_.restype = c_void_p
    
    file_name = c_char * 78; file_name1 = file_name(); file_name2 = file_name()
    file_name1 = infile
    file_name2 = fnamex
    fr1 = len(file_name1)
    fr2 = len(file_name2)
    
    addmodule.pdbtop_(c_char_p(file_name1),byref(c_float(wradiu1)),byref(c_float(wradiu2)),c_char_p(file_name2),byref(c_int(fr1)),byref(c_int(fr2)))
    
    cmd.load(fnamex)

'''
setaxis(sel1,sel2)

sel1: Name of superposing P-object.
sel2: Name of superposed P-object.
'''
def setaxis(sel1,sel2):
    global axis
    cmd.create('tmp', sel1)
    cmd.super('tmp', sel2)
    t = cmd.get_object_matrix('tmp')
    a = (t[6]-t[9])/math.sqrt((t[6]-t[9])**2+(t[8]-t[2])**2+(t[1]-t[4])**2)
    b = (t[8]-t[2])/math.sqrt((t[6]-t[9])**2+(t[8]-t[2])**2+(t[1]-t[4])**2)
    c = (t[1]-t[4])/math.sqrt((t[6]-t[9])**2+(t[8]-t[2])**2+(t[1]-t[4])**2)
    axis = [a, b, c]
    print axis
    cmd.delete('tmp')

'''
Capswch(theta,alpha,center,axis1)

theta: Angle for selecting operational P-object.
alpha: Angle of rotation.
center:  Center of rotation(default = None: center is set to obj_center_init).
axis1: Axis of rotation(default = None: use axis defined by setaxis).
'''
def capswch(theta,alpha,center=None,axis1=None):
    theta_radi = math.radians(theta)
    alpha_radi = math.radians(alpha)
    obj_center_xyz = numpy.zeros([60,3])
    normal_vect = numpy.zeros([3])
    
    cmd.reset()

    if center == None:
        obj_center_xyz = obj_center_init
    else:
        obj_center_xyz = get_obj_center()

    if axis1 == None:
        axis1 = axis

    for i in range(60):
        exec('ring%s = "_1"'%i)
        angle_cos = (axis1[0]*obj_center_xyz[i][0] + axis1[1]*obj_center_xyz[i][1] + axis1[2]*obj_center_xyz[i][2]) / (math.sqrt(axis1[0]**2+axis1[1]**2+axis1[2]**2)*math.sqrt(obj_center_xyz[i][0]**2+obj_center_xyz[i][1]**2+obj_center_xyz[i][2]**2))
        if theta_radi >= math.acos(angle_cos):
            cmd.pseudoatom ('ring%s'%i,pos=[obj_center_xyz[i][0], obj_center_xyz[i][1], obj_center_xyz[i][2]])
            normal_vect[0] = axis1[1]*obj_center_xyz[i][2] - axis1[2]*obj_center_xyz[i][1]
            normal_vect[1] = axis1[2]*obj_center_xyz[i][0] - axis1[0]*obj_center_xyz[i][2]
            normal_vect[2] = axis1[0]*obj_center_xyz[i][1] - axis1[1]*obj_center_xyz[i][0]
            cmd.origin (position=[0.,0.,0.])
            cmd.rotate ([normal_vect[0],normal_vect[1],normal_vect[2]],(theta - math.degrees(math.acos(angle_cos))), 'ring%s'%i)
            cmd.origin ('ring%s'%i)
            cmd.rotate ([normal_vect[0],normal_vect[1],normal_vect[2]],alpha, 'OBJ%s'%i)
            exec('print "Rotate OBJ%s"'%i)
            cmd.delete('ring%s'%i)

'''
Capswch(theta,alpha,center,axis1)

theta: Angle for selecting operational P-object.
alpha: Angle of rotation.
center:  Center of rotation(default = None: center is set to obj_center_init).
axis1: Axis of rotation(default = None: use axis defined by setaxis).
'''
def captwst(theta,alpha,center=None,axis1=None):
    theta_radi = math.radians(theta)
    alpha_radi = math.radians(alpha)
    obj_center_xyz = numpy.zeros([60,3])
    normal_vect = numpy.zeros([3])
    stored.pos_ring = []
    
    cmd.reset()
    
    if center == None:
        obj_center_xyz = obj_center_init
    else:
        obj_center_xyz = get_obj_center()
    
    if axis1 == None:
        axis1 = axis
    
    for i in range(60):
        exec('ring%s = "_1"'%i)
        angle_cos = (axis1[0]*obj_center_xyz[i][0] + axis1[1]*obj_center_xyz[i][1] + axis1[2]*obj_center_xyz[i][2]) / (math.sqrt(axis1[0]**2+axis1[1]**2+axis1[2]**2)*math.sqrt(obj_center_xyz[i][0]**2+obj_center_xyz[i][1]**2+obj_center_xyz[i][2]**2))
        if theta_radi >= math.acos(angle_cos):
            cmd.pseudoatom ('ring%s'%i,pos=[obj_center_xyz[i][0], obj_center_xyz[i][1], obj_center_xyz[i][2]])
            normal_vect[0] = axis1[1]*obj_center_xyz[i][2] - axis1[2]*obj_center_xyz[i][1]
            normal_vect[1] = axis1[2]*obj_center_xyz[i][0] - axis1[0]*obj_center_xyz[i][2]
            normal_vect[2] = axis1[0]*obj_center_xyz[i][1] - axis1[1]*obj_center_xyz[i][0]
            cmd.origin (position=[0.,0.,0.])
            cmd.rotate ([normal_vect[0],normal_vect[1],normal_vect[2]],(theta - math.degrees(math.acos(angle_cos))),'ring%s'%i)
            cmd.origin ('ring%s'%i)
            cmd.iterate_state(1,'ring%s'%i,'stored.pos_ring = [x,y,z]')
            print stored.pos_ring
            cmd.rotate ([stored.pos_ring[0],stored.pos_ring[1],stored.pos_ring[2]],alpha,'OBJ%s'%i)
            exec('print "Rotate OBJ%s"'%i)
            cmd.delete('ring%s'%i)












