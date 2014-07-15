#!/usr/bin/env python

#this program convert movie data to atom data = for each spherocylinder make residue
#consiting of two atoms at begining and end then in vmd use cpk to draw cylinders
#no it also ads a spherocylinder for patch

import os
import sys
import math
import optparse
import commands
import string
import random
import usefulmath

SC = 10
SCN = SC+0
SCA = SC+1
PSC = SC+2
CPSC = SC+3
CHPSC = SC+4
CHCPSC = SC+5
TPSC = SC+6
TCPSC = SC+7
TCHPSC = SC+8
TCHCPSC = SC+9
SP = 30
SPN = SP+0
SPA = SP+1

def convert_geotype(geotype):
    if(geotype == "SC"):
        return SC
    if(geotype == "SCN"):
        return SCN
    if(geotype == "SCA"):
        return SCA
    if(geotype == "PSC"):
        return PSC
    if(geotype == "CPSC"):
        return CPSC
    if(geotype == "CHPSC"):
        return CHPSC
    if(geotype == "CHCPSC"):
        return CHCPSC
    if(geotype == "TPSC"):
        return TPSC
    if(geotype == "TCPSC"):
        return TCPSC
    if(geotype == "TCHPSC"):
        return TCHPSC
    if(geotype == "TCHCPSC"):
        return TCHCPSC
    if(geotype == "SP"):
        return SP
    if(geotype == "SPN"):
        return SPN
    if(geotype == "SPA"):
        return SPA


def write_data(data,box,particles, types, geotypes, params,outfilename,outfilename2, switch, switchtypes, switchfile):
    f = open(outfilename, 'w')
    outstring=""
    if (switch):
	swstring=""
	sf = open(switchfile, 'w')
    #DEBUG
    #print("Length of types")
    #print(len(types))
    #print("Len Data")
    #print(len(data))
    #print(types)
    #print(geotypes)
    for i in range(len(data)):
	newline="CRYST1 %8.3f %8.3f %8.3f  90.00  90.00  90.00 P 1           1\n" % (box[i][0],box[i][1],box[i][2])
	outstring=outstring+newline
	frame=data[i]
	#print("len frame")
	#print(len(frame))
	atm=1
	for j in range(len(frame)):
	    [x,y,z,vx,vy,vz,px,py,pz,sw]=frame[j][:]
	    [x,y,z] = usefulmath.usepbc([x,y,z],box[i])
	    #if (types[j] <SP):
	    #print("j : %d\n" % j)
	    if (convert_geotype(geotypes[particles[j]]) <SP):
		# Handle the switching
		csw = 2 * sw + types[particles[j]] +0
		psw = 2 * sw + types[particles[j]] +1
		tsw = 2 * sw + types[particles[j]] +2
		if (switch):
		    swstring = swstring + "%d %d %d %d " % (csw, csw, psw, psw)
		vec=usefulmath.vec_normalize([vx,vy,vz])
		#print head
		leng = params[particles[j]][6]
		nx=x+leng/2*vec[0]
		ny=y+leng/2*vec[1]
		nz=z+leng/2*vec[2]
		#print("max j: %d; j: %d" % (len(frame), j))
		newline="ATOM  %5d  C%1d  PSC F%4d    % 8.2f% 8.2f% 8.2f  1.00%6.2f\n" % (atm,particles[j],(j+1)%1000,nx,ny,nz,csw)
		outstring=outstring+newline
		atm=atm+1
		#print tail
		nx=x-leng/2*vec[0]
		ny=y-leng/2*vec[1]
		nz=z-leng/2*vec[2]
		newline="ATOM  %5d  C%1d  PSC F%4d    % 8.2f% 8.2f% 8.2f  1.00%6.2f\n" % (atm,particles[j],(j+1)%1000,nx,ny,nz,csw)
		outstring=outstring+newline
		atm=atm+1
		#print patch
		#patch=params[types[j]][3]+2*params[types[j]][4]
		patch=params[particles[j]][4]+2*params[particles[j]][5]
		#print headpatch
		patchmove=math.fabs(0.5*math.cos((patch+10*pow(patch/180,2))/360*math.pi)) # displacement
		#patchmove= params[particles[j]][2] * patch/ 1080 + params[particles[j]][2] / 12 # displacement
		# The displacement varries between sigma / 12 and sigma / 4 depending on the angle
		nx=x+leng/2*vec[0]+patchmove*px
		ny=y+leng/2*vec[1]+patchmove*py
		nz=z+leng/2*vec[2]+patchmove*pz
		newline="ATOM  %5d  P%1d  PSC F%4d    % 8.2f% 8.2f% 8.2f  1.00%6.2f\n" % (atm,particles[j],(j+1)%1000,nx,ny,nz,psw)
		outstring=outstring+newline
		atm=atm+1
		#print tailpatch
		nx=x-leng/2*vec[0]+patchmove*px
		ny=y-leng/2*vec[1]+patchmove*py
		nz=z-leng/2*vec[2]+patchmove*pz
		newline="ATOM  %5d  P%1d  PSC F%4d    % 8.2f% 8.2f% 8.2f  1.00%6.2f\n" % (atm,particles[j],(j+1)%1000,nx,ny,nz,psw)
		outstring=outstring+newline
		atm=atm+1
		if (convert_geotype(geotypes[particles[j]]) >=TPSC):
		    #print second patch
		    if (switch):
			swstring = swstring + "%d %d " % (tsw, tsw)
		    vec2=[px,py,pz]
		    vec2=usefulmath.rotatevec(vec2,vec,params[particles[j]][7])
		    patch=params[particles[j]][8]+2*params[particles[j]][9]
		    #print headpatch
		    patchmove=0.5*math.cos((patch+10*pow(patch/180,2))/360*math.pi) # displacement
		    #patchmove= params[particles[j]][2] * patch/ 1080 + params[particles[j]][2] / 12 # displacement
		    # The displacement varries between sigma / 12 and sigma / 4 depending on the angle
		    nx=x+leng/2*vec[0]+patchmove*vec2[0]
		    ny=y+leng/2*vec[1]+patchmove*vec2[1]
		    nz=z+leng/2*vec[2]+patchmove*vec2[2]
		    newline="ATOM  %5d  T%1d  PSC F%4d    % 8.2f% 8.2f% 8.2f  1.00%6.2f\n" % (atm,particles[j],(j+1)%1000,nx,ny,nz,tsw)
		    outstring=outstring+newline
		    atm=atm+1
		    #print tailpatch
		    nx=x-leng/2*vec[0]+patchmove*vec2[0]
		    ny=y-leng/2*vec[1]+patchmove*vec2[1]
		    nz=z-leng/2*vec[2]+patchmove*vec2[2]
		    newline="ATOM  %5d  T%1d  PSC F%4d    % 8.2f% 8.2f% 8.2f  1.00%6.2f\n" % (atm,particles[j],(j+1)%1000,nx,ny,nz,tsw)
		    outstring=outstring+newline
		    atm=atm+1
		    
	    else:
		#print sphere
		nx=x
		ny=y
		nz=z
		newline="ATOM  %5d  S%1d  PSC F%4d    % 8.2f% 8.2f% 8.2f  1.00%6.2f\n" % (atm,particles[j]%SP,(j+1)%1000,nx,ny,nz,0)
		atm=atm+1
		outstring=outstring+newline
	outstring=outstring+"END\n"
	if (switch):
	    swstring = swstring + "\n"
    f.write(outstring)
    f.close()
    if (switch):
	sf.write(swstring)
	sf.close()

    return 0

def write_psf(outfilename,box,particles, types, geotypes ,params):
    #C-cylinder,P-patch,S-sphere
    f = open(outfilename, 'w')
    outstring="PSF\n"
    npart=0
    for i in range(len(particles)):
	#if (types[i] <SP):
	if (convert_geotype(geotypes[particles[i]]) <SP):
	    if (convert_geotype(geotypes[particles[i]]) <TPSC):
		npart=npart+4
	    else:
		npart=npart+6
	else:
	    npart=npart+1
    newline="%8d !NATOM\n" % (npart)
    outstring=outstring+newline
    bonds=[]
    atm=1
    for i in range(len(particles)):
	#if (particles[i] <SP):
	if (convert_geotype(geotypes[particles[i]]) <SP):
	    newline="%8d C%03d %4d PSC  C%1d   C%1d                   \n" % (atm,particles[i],i+1,particles[i],particles[i])
	    bonds.append(atm)
	    atm=atm+1
	    outstring=outstring+newline
	    newline="%8d C%03d %4d PSC  C%1d   C%1d                   \n" % (atm,particles[i],i+1,particles[i],particles[i])
	    bonds.append(atm)
	    atm=atm+1
	    outstring=outstring+newline
	    newline="%8d P%03d %4d PSC  P%1d   P%1d                   \n" % (atm,particles[i],i+1,particles[i],particles[i])
	    bonds.append(atm)
	    atm=atm+1
	    outstring=outstring+newline
	    newline="%8d P%03d %4d PSC  P%1d   P%1d                   \n" % (atm,particles[i],i+1,particles[i],particles[i])
	    outstring=outstring+newline
	    bonds.append(atm)
	    atm=atm+1
	    if (convert_geotype(geotypes[particles[i]]) >=TPSC):
		newline="%8d T%03d %4d PSC  T%1d   T%1d                   \n" % (atm,particles[i],i+1,particles[i],particles[i])
		bonds.append(atm)
		atm=atm+1
		outstring=outstring+newline
		newline="%8d T%03d %4d PSC  T%1d   T%1d                   \n" % (atm,particles[i],i+1,particles[i],particles[i])
		outstring=outstring+newline
		bonds.append(atm)
		atm=atm+1
	    
	else:
	    newline="%8d S%03d %4d PSC  S%1d   S%1d                   \n" % (atm,particles[i]%SP,i+1,particles[i]%SP,particles[i]%SP)
	    outstring=outstring+newline
	    atm=atm+1
    newline="%8d !NBOND\n" % (len(bonds)/2)
    outstring=outstring+newline
    newline=""
    for i in range(len(bonds)):
	newline+=" %7d" % (bonds[i])
	if (i%8 ==7):
	    outstring=outstring+newline+"\n"
	    newline=""
    outstring=outstring+"%-64s"%(newline)
    f.write(outstring)
    f.close()

    return 0


def read_input(infilename):
    data=[]
    box=[]
    frame=[]
    inp=open(infilename)
    i=0
    atomnum=0
    for line in inp:
	linesplit=line.split()
	if (len(linesplit)==1):
	    atomnum = int(linesplit[0])
	    #print atomnum
	    i=0
	else:
	    if (len(linesplit)==3):
		[bx,by,bz]=linesplit[:]
		box.append([float(bx),float(by),float(bz)])
	    else:
		if (len(linesplit)==6):
		    [sweep,num,boxstr,bx,by,bz]=linesplit[:]
		    box.append([float(bx),float(by),float(bz)])
		else:
		    [x,y,z,vx,vy,vz,px,py,pz,sw]=linesplit[:]
		    frame.append([float(x),float(y),float(z),float(vx),float(vy),float(vz),float(px),float(py),float(pz),float(sw)])
		    i=i+1
	if ( (atomnum!=0) and (i==atomnum) ):
	    #print frame
	    data.append(frame)
	    frame=[]
    if (atomnum==0):
	data.append(frame)
    
    return [box,data]

def write_vmd(outfilename,outfilename2,particles, types, geotypes, params, num_frames, switch, switchtypes):
    f = open("vmd.script", 'w')
    outstring=""
    outstring+="proc setlook {} {\n"
    outstring+="rotate stop\n"
    outstring+="color Display Background white\n"
    outstring+="display projection orthographic\n"
    outstring+="mol delrep 0 0\n"
    outstring+="mol material Edgy\n"
    vmdpckratio=1.50
    sizeratio=1.1

    j=0
    #DEBUG
    #print particles
    #print types
    particles_set = list(set(particles))
    switchtypes_set = list(set(switchtypes))
    t_clean = []
    for t in particles_set:
	if t != []:
	    t_clean.append(int(t))
    for t in switchtypes_set:
	if t != []:
	    if t != -1:
		t_clean.append(int(t))

    type_set = list(set(t_clean))
    num_types = 0

    for i in range(len(type_set)):
        geotype = geotypes[type_set[i]]
        if ((convert_geotype(geotype) == CHCPSC) or (convert_geotype(geotype) == CHPSC) or (convert_geotype(geotype) == PSC) or (convert_geotype(geotype) == CPSC) or (convert_geotype(geotype) == TCHCPSC) or (convert_geotype(geotype) == TCHPSC) or (convert_geotype(geotype) == TPSC) or (convert_geotype(geotype) == TCPSC)):
	    num_types += 2
            radius = params[type_set[i]][1]
            if (switch):
        	outstring+="mol selection \"beta %4.2f\"\n" % float(type_set[i]+0)
            else:
        	outstring+="mol selection \"name C%d\"\n" %(i + 1)
            outstring+="mol representation CPK %f %f 20 20\n" %(radius,radius*vmdpckratio)
            outstring+="mol color ColorID %d\n"%(j)
            #outstring+="mol color Beta\n"
            outstring+="mol addrep 0\n"
            j = j + 1
            patch=params[type_set[i]][4]+2*params[type_set[i]][5]
            angle = patch / 360 * math.pi
            rpatch=math.sin(angle*pow((2-2*angle/math.pi),2))
            #Noah alternative size calculation
            #displ = params[type_set[i]][4] *  params[type_set[i]][2] / 1080 + params[type_set[i]][2] / 12
            #rpatch = radius * math.sqrt(pow(math.sin(angle),2) + pow(math.cos(angle),2) / 4)
            #rpatch = math.sqrt(pow(radius * math.sin(angle),2) + pow(radius * math.cos(angle) - displ, 2))
            
            if (switch):
        	outstring+="mol selection \"beta %4.2f\"\n" % float(type_set[i]+1)
            else:
        	outstring+="mol selection \"name P%d\"\n" %(i + 1)
            outstring+="mol representation CPK %f %f 20 20\n" %(rpatch,rpatch*vmdpckratio)
            outstring+="mol color ColorID %d\n"%(j)
            #outstring+="mol color Beta\n"
            outstring+="mol addrep 0\n"
            j = j + 1
            if ((convert_geotype(geotype) == TCHCPSC) or (convert_geotype(geotype) == TCHPSC) or (convert_geotype(geotype) == TPSC) or (convert_geotype(geotype) == TCPSC)):
        	patch=params[type_set[i]][8]+2*params[type_set[i]][9]
        	#rpatch=sizeratio*math.sin((patch+10)/2/180*math.pi)
        	angle = patch / 360 * math.pi
        	displ = params[type_set[i]][8] *  params[type_set[i]][2] / 1080 + params[type_set[i]][2] / 12
        	rpatch = radius * math.sqrt(pow(math.sin(angle),2) + pow(math.cos(angle),2) / 4)
        	rpatch = math.sqrt(pow(radius * math.sin(angle),2) + pow(radius * math.cos(angle) - displ, 2))
        	if (switch):
        	    outstring+="mol selection \"beta %4.2f\"\n" % float(type_set[i]+2)
        	else:
        	    outstring+="mol selection \"name T%d\"\n" %(i + 1)
        	outstring+="mol representation CPK %f %f 20 20\n" %(rpatch,rpatch*vmdpckratio)
        	outstring+="mol color ColorID %d\n"%(j)
        	#outstring+="mol color Beta\n"
        	outstring+="mol addrep 0\n"
        	j = j + 1
            

        if ((convert_geotype(geotype) == SPN) or (convert_geotype(geotype) == SPA)):
	    num_types += 1
            radius=sizeratio*params[type_set[i]][0]
            outstring+="mol selection \"name S%d\"\n" % type_set[i]
            outstring+="mol representation CPK %f %f 20 20\n" %(radius,0.0)
            outstring+="mol color ColorID %d\n"%(j)
            outstring+="mol addrep 0\n"
            j=j+1


    outstring+="axes location off\n"
    outstring+="}\n"
    outstring+="mol load psf %s \n" %(outfilename2)
    outstring+="mol addfile %s waitfor 2000  0\n"%(outfilename)
    outstring+="setlook\n"

    if (switch):
	outstring+="""set molid 0
    # read in charge data
    set n [molinfo $molid get numframes]; # Somehow that does not really work
    set n """
	outstring+="%d" % (num_frames)
	outstring+="""
puts "reading betas"
set fp [open "switch.dat" r]
for {set i 0} {$i < $n} {incr i} {
    set bt($i) [gets $fp]
}
close $fp

# procedure to change the charge field from the data in $chrg
proc do_charge {args} {
    global bt molid
    set a [molinfo $molid get numatoms]
    set f [molinfo $molid get frame]
    for {set i 0} {$i < $a} {incr i} {
	    set s [atomselect $molid "index $i"]
	    $s set beta [lindex $bt($f) $i]
    }
}

#trace variable vmd_frame($molid) w do_charge
# turn on update of the coloring in each frame.
for {set i 0} {$i < """
	outstring+="%d" % (num_types)
	outstring+="""} {incr i} {
    trace variable vmd_frame($i) w do_charge
    mol selupdate   $i $molid on
    mol colupdate   $i $molid on
}

# set color mapping parameters
color scale method RGB

# go back to the beginning and activate the additional feartures.
animate goto start
do_charge
"""

    f.write(outstring)
    f.close()

    return 0

def read_top(topfilename):
    i=0
    switch=False
    #types=[]
    types=[[]]
    geotypes=[[]]
    params=[]
    for i in range(30):
	params.append([])
	geotypes.append([])
	types.append([])
    #DEBUG
    #print("Geotypes empty")
    #print(geotypes)
    keyword=""
    keystr=""
    molname=""
    sys=[[],[]]
    #0: molname, 1: type, 2: switchtype
    mol=[[],[],[]]
    inp=open(topfilename)
    for line in inp:
	#remove comment
	ind = line.find("#")
	if (ind>-1):
	    line=line[:ind]
	line = line.strip()
	if (len(line)>0):
	    #if not empty continue
	    #print "line:%s"%(line)
	    if (line[0]=="["):
		#add keyword
		line=line[1:]
		[before,sep,after]=line.partition("]")
		if (len(sep)>0):
		    before = before.strip()
		    keyword=before[:]
		    keyword=keyword.upper()
	    else:
		if (keyword == "TYPES"):
		    oneline = line.split()
		    name = oneline[0]
		    type = int(oneline[1])
		    geotype = oneline[2]
		    linesplit = oneline[3:]
		    #[before,sep,after]=line.partition(":")
		    #after = after.strip()
		    #linesplit=after.split()
		    #print "b:%s a:%s"%(before,after)
		    #if ( (long(before) == SCN) or (long(before) == SPN) ):
		    #DEBUG
		    #print("Linesplit")
		    #print(linesplit)
		    #print(type)
		    #print(geotype)
		    if ( (geotype == "SCN") or (geotype == "SPN") ):
			    if (len(linesplit) != 2):
			        print "TOPOLOGY ERROR: wrong number of parameters for given type %s"%(type)
			        return [[],[]]
			    else:
			        params[type]=[float(linesplit[0]),float(linesplit[1])]
			        types[type]=type
			        geotypes[type]=geotype
		    #if ( (long(before) == SCA) or (long(before) == SPA) ):
		    if ( (geotype == "SCA") or (geotype == "SPA") ):
			    if (len(linesplit) != 4):
			        print "TOPOLOGY ERROR: wrong number of parameters for given type %s"%(type)
			        return [[],[]]
			    else: 
			        params[type]=[float(linesplit[0]),float(linesplit[1]),float(linesplit[2]), float(linesplit[3])]
			        types[type]=type
			        geotypes[type]=geotype
		    #if ( (long(before) == PSC) or (long(before) == CPSC) ):
		    if ( (geotype == "PSC") or (geotype == "CPSC") ):
			    if (len(linesplit) != 7):
			        print "TOPOLOGY ERROR: wrong number of parameters for given type %s" %(type)
			        return [[],[]]
			    else: 
			        params[type]=[float(linesplit[0]),float(linesplit[1]),float(linesplit[2]),float(linesplit[3]),float(linesplit[4]), float(linesplit[5]), float(linesplit[6])]
			        types[type]=type
			        geotypes[type]=geotype
		    if ( (geotype == "CHPSC")or (geotype == "CHCPSC") ):
			    if (len(linesplit) != 8):
			        print "TOPOLOGY ERROR: wrong number of parameters for given type %s" %(type)
			        return [[],[]]
			    else: 
			        params[type]=[float(linesplit[0]),float(linesplit[1]),float(linesplit[2]),float(linesplit[3]),float(linesplit[4]), float(linesplit[5]), float(linesplit[6])]
			        types[type]=type
			        geotypes[type]=geotype
		    if ( (geotype == "TPSC") or (geotype == "TCPSC") ):
			    if (len(linesplit) != 10):
			        print "TOPOLOGY ERROR: wrong number of parameters for given type %s" %(type)
			        return [[],[]]
			    else: 
			        params[type]=[float(linesplit[0]),float(linesplit[1]),float(linesplit[2]),float(linesplit[3]),float(linesplit[4]), float(linesplit[5]), float(linesplit[6]),float(linesplit[7]), float(linesplit[8]), float(linesplit[9]) ]
			        types[type]=type
			        geotypes[type]=geotype
		    if ( (geotype == "TCHPSC")or (geotype == "TCHCPSC") ):
			    if (len(linesplit) != 11):
			        print "TOPOLOGY ERROR: wrong number of parameters for given type %s" %(type)
			        return [[],[]]
			    else: 
			        params[type]=[float(linesplit[0]),float(linesplit[1]),float(linesplit[2]),float(linesplit[3]),float(linesplit[4]), float(linesplit[5]), float(linesplit[6]),float(linesplit[7]), float(linesplit[8]), float(linesplit[9]) ]
			        types[type]=type
			        geotypes[type]=geotype
		else:
		    if (keyword == "MOLECULES"):
			#read molecule
			#DEBUG
			#print(line)
			if (molname==""):
			    [before,sep,after]=line.partition("{")
			    before = before.strip()
			    [before,sep,after]=before.partition(":")
			    before = before.strip()
			    if (len(before)>0):
				molname=before[:]
				molname.upper()
			    after = after.strip()
			    if (len(after)>0):
				line=after[:]
				line = line.strip()
			#readmolparams
			ind = line.find("}")
			if (ind > 1):
			    line=line[:ind]
			[before,sep,after]=line.partition(":")
			before=before.strip()
			before=before.upper()
			if ((before=="PARTICLES") and (molname!="")):
			    after=after.strip()
			    mol[0].append(molname)
			    mol[1].append(after.split()[0]) #[0] NB
			    if len(after.split()) > 1:
				mol[2].append(after.split()[1])
				switch=True
			    else:
				mol[2].append(-1)
			    #DEBUG
			    #print("After")
			    #print(after.split())
			    #print(mol[1])
			#closes mol
			if (ind>-1):
			    molname=""
		    else:
			if (keyword == "SYSTEM"):
			    #read system
			    linesplit=line.split()
			    if (len(linesplit) != 2):
				print "TOPOLOGY ERROR: invalid number of input in System"
				return [[],[]]
			    sys[0].append(linesplit[0])
			    sys[1].append(long(linesplit[1]))
			else:
			    if ((keyword == "EXTER") or (keyword == "EXCLUDE")):
				#we do not displaye external potential 
				print ""
			    else:
				print "TOPOLOGY ERROR: invalid keyword %s on line %s"%(keyword,line)
				return [[],[]]
    particles=[]
    switchtype=[]
    for i in range(len(sys[0])):
        for k in range(sys[1][i]):
            for j in range(len(mol[1])):
                if sys[0][i] == mol[0][j]:
                    particles.append(int(mol[1][j]))
                    switchtype.append(int(mol[2][j]))

    #print("Particles")
    #print(particles)

    #DEBUG
    #print("MOL, SYS")
    #print(mol)
    #print(sys)
    #print("PARAMS")
    #print(particles)
    #print(params)
    #print("Types")
    #print(switchtype)
    #print(geotypes)
    #print(types)
    #print(particles)
    #return [types, geotypes, params]
    return [particles, types, geotypes, params, switch, switchtype]


def make(infilename,outfilename,outfilename2,topfilename):
    print "Reading topology..."
    [particles, types,geotypes, params, switch, switchtypes]=read_top(topfilename)
    if (len(types) < 1):
	print "ERROR: types have not been read"
	return 1
    print "Reading coordinates..."
    [box,data]=read_input(infilename)
    if ((len(data) < 1 ) or ( len(box) < 1 )):
	print "ERROR: data has not been read"
	return 1
    if ( (len(data[0]) % len(particles)) != 0):
	print "ERROR: data with length %d does not fit types with length %d" % (len(data[0]),len(particles))
	return 1
    print "Writing data..."
    switchfile = "switch.dat"
    write_data(data,box,particles, types, geotypes, params, outfilename, outfilename2, switch, switchtypes, switchfile)
    print "Writing psf..."
    write_psf(outfilename2,box,particles, types, geotypes, params)
    print "Writing vmd..."
    #print(types)
    write_vmd(outfilename,outfilename2,particles, types, geotypes, params, len(data), switch, switchtypes)
    print "Done"
    return 0

parser=optparse.OptionParser()
help="""Usage:
%prog [options] 
"""
parser.set_usage(help)
parser.add_option(
    "-i",
    "--input",
    help="Set file from which you want to load data",
    dest="infilename",
    default="movie"
    )
parser.add_option(
    "-o",
    "--output",
    help="Set to which file you want to save data",
    dest="outfilename",
    default="movie.pdb"
    )
parser.add_option(
    "--psf",
    help="Set to which file you want to save connectivity - psf",
    dest="outfilename2",
    default="movie.psf"
    )
parser.add_option(
    "-t",
    "--topology",
    help="File where topology is read from",
    dest="topfilename",
    default="top.init"
    )

(options,arguments)=parser.parse_args()
make(options.infilename,options.outfilename,options.outfilename2,options.topfilename)
# vim: set noexpandtab ts=8: 
