import math
import random

PI = 3.141592653589793238462643383279  

def dot_product(A,B):
    dp=0
    for i in range(0,3):
        dp=dp+A[i]*B[i]
    return dp


#  returns a vector perpendicular to A
#  nothing special about the vector except that it's one of the perpendicular options and is normalized

def perp_vec(A):
#    x=A[0]
#    y=A[1]
#    if x == 0:
#       x=1
#    elif y == 0:
#       y=1
#    else:
#       ratio=y/x
#       y=x*ratio*2
#    some_vector= [x,y,A[2]]
#    some_vector= vec_normalize(some_vector)
#    return cross_product(A,some_vector)
    A=vec_normalize(A)
    dp=1
    while ((dp == 1) or (dp==-1)):
        newvec=vec_random()
        dp=dot_product(newvec,A)
    for i in range(len(newvec)):
        newvec[i]=newvec[i]-A[i]*dp

    return newvec

# Perform the multiplication of a matrix A and a vector B where A is the
# first argument and B is the second argument.  The routine will
# return AxB

def matrix_vec_multiply(A,B):
    AB=[]
    side=len(A)
    
    for i in range (0,side):
        sum=0
        # index the row vector from A
        RA=A[i]
        # Now find the dot product of this row with B
        for j in range (0,side):
            sum= sum + RA[j]*B[j]
        AB.append(sum)
    return AB



# Calculate the proportion of each integer in a list of integers
# Arguments:
# ilist: the list of integers
#
def calc_proportions(ilist):
    plist=[]
    # Find all the different integers
    itypelist=ilist[0]
    plist.append([itypelist,0])
    for i in ilist: #tady jsem skoncil
        if itypelist.index(i) == -1:
            itypelist.append(i)
            plist.append([i,0])

    # Now count the number of each type
    for i in ilist:
        vp=plist[itypelist.index(i)]
        v=vp[1] 
        if itypelist.index(i) == -1:
            print("could not calculate proportions")
        plist[itypelist.index(i)][1]=(v + 1)
    return plist

# calculate the average of a list of numerical data
#
def average(data,froom=0,to=0):
        sum= 0
        cnt= 0
        if ((to == 0) or (to >= len(data))):
            to= len(data) - 1
        if ((froom < 0) or (froom >= len(data))):
            return 0

        for i in range(froom,to+1):
            sum=sum + data[i]
            cnt=cnt+1


        return sum/cnt

# calculate the stdev of a list of numerical data
#
def stdev(data,froom=0,to= 0):
        sum= 0
        cnt= 0
        if ((to == 0) or (to >= len(data))):
            to= len(data) - 1 
        if ((froom < 0) or (froom >= len(data))):
            return 0

        av= average(data,froom,to)

        for  i in range(froom,to + 1):
            sum=sum + (data[i] - av)*(data[i] - av)
            cnt=cnt+1


        var= sum/cnt
        return math.sqrt(var)


# calculate an autocorrelation function from a list of numbers
#
def acorr(data):
        tot=len(data)
        to=tot - 1
        acorr=[]

        av= average(data)
        print(av)
        sum= 0
        for x in range(0,tot):
            sum=data[x]*data[x] + sum
        variance= sum/tot - av*av
        print(variance)

        for y in range(0,tot):
            sum1= 0
            sum2= 0
            sum3= 0
            for x in range(0,tot - y):
                sum1= data[x]*data[x+y]+sum1
                sum2= sum2 + data[x]
                sum3= sum3 + data[x+y]
            ct= (sum1-sum2*sum3/(tot-y))/((tot-y)*variance)
            acorr.append(ct)
        return acorr


# Distance between two vectors
#
def vec_distance(vec1,vec2):

    if len(vec1) != len(vec2):
        print("cannot find distance between vectors:", vec1, " and", vec2," because they are different lengths")

    sum= 0
    for i in range(0,len(vec1)):
        diff= vec1[i] - vec2[i]
        sum= sum + diff*diff

    return math.sqrt(sum)


# Normalize a vector
#
def vec_normalize(vec):
    # Find the magnitude
    mag2= 0.0
    for val in vec:
        mag2= mag2 + val*val
    mag= math.sqrt(mag2)
    
    nvec=[]
    # Divide by the magnitude to produce normalized vector
    for val in vec:
        nvec.append(val/mag)
    return nvec

# Scale a vector
#
def vec_scale(vec,scale):
    svec=[]
    for val in vec:
        svec.append(val*scale)
    return svec


# Take a list of integers and construct a list that contains no
# duplication
#
def uniquelist(original):
    short=[]
    short.append(original[0])

    for element in original:
        if  short.index(element) == -1 :
            short.append(element)

    return short


# ROTATE
# Take a {x y z} point and rotate it an angle sigma around either
# the x y or z axis.
#
def rotation_matrix(axis,phi):
    
    x= axis[0]
    y= axis[1]
    z= axis[2]
    x2= x * x
    y2= y * y
    z2= z * z

    rcos= math.cos(phi)
    rsin= math.sin(phi)

    row0= []
    row1= []
    row2= []

    row0.append(x2 + (y2+z2)*rcos)
    row0.append(x*y*(1.0-rcos) - z*rsin)
    row0.append(x*z*(1.0-rcos) + y*rsin)
    row1.append(x*y*(1.0-rcos) + z*rsin)
    row1.append(y2 + (x2+z2)*rcos)
    row1.append(y*z*(1.0-rcos) - x*rsin)
    row2.append(x*z*(1.0-rcos) - y*rsin)
    row2.append(y*z*(1.0-rcos) + x*rsin)
    row2.append(z2 + (x2+y2)*rcos)


    rotation_matrix=[]
    rotation_matrix.append(row0)
    rotation_matrix.append(row1)
    rotation_matrix.append(row2)

    return rotation_matrix

# cross_product
def cross_product(A,B):
    cp=[]
    Ax=A[0]
    Ay=A[1]
    Az=A[2]
    Bx=B[0]
    By=B[1]
    Bz=B[2]

    cp.append( Ay*Bz - Az*By)
    cp.append( -Ax*Bz + Az*Bx)
    cp.append( Ax*By - Ay*Bx)
    return cp

# addition of vectors
def vec_add(A,B):
    C=[]
    C.append(A[0] + B[0])
    C.append(A[1] + B[1])
    C.append(A[2] + B[2])
    return C

# multiplication of matrix
def matrix_multiply(A,B):
    C=[]
    for i in range(0,3):
        row=[]
        for j in range(0,3):
            sum= 0
            for k in range(0,3):
                sum=sum+A[i][k]*B[k][j]
            row.append(sum)
        C.append(row)
    return C

# sum of matrix
def matrix_sum(A,B):
    C=[]
    for i in range(0,3):
        row=[]
        for j in range(0,3):
            sum=A[i][j]+B[i][j]
            row.append(sum)
        C.append(row)
    return C


#generate random unit vector
def vec_random():
    a=10

    while (a > 1.0):
        xi1 = 1.0 - random.random()*2
        xi2 = 1.0 - random.random()*2
        a = xi1*xi1 + xi2*xi2
    
    b = 2.0 * math.sqrt(1.0 - a)
    x = xi1 * b
    y = xi2 * b
    z = 1.0 - 2.0*a
    randomvec=[x,y,z]

    return randomvec

#construct rotation matrix for transformation to given orientation
def rotation_byvec(orient):
    initial_orient=[0,0,1]
    initial_orient=vec_normalize(initial_orient)
    orient= vec_normalize(orient)
    rotation_angle= math.acos(dot_product(initial_orient,orient))
    rotation_axis=cross_product(initial_orient,orient)
    if vec_distance(rotation_axis,[0,0,0]) == 0:
        # if initial_orient and orient are parellel then we  can use any vector perpendicular to this line as the rotation axis 
        rotation_axis=perp_vec(orient)
    rotation_axis= vec_normalize(rotation_axis)
    rotation_matr= rotation_matrix(rotation_axis,rotation_angle)

    return rotation_matr

#construct rotation matrix for transformation to given orientation
def rotation_by2vec(orient,orient2):
    initial_orient=[0,0,1]
    initial_orient2=[1,0,0]
    initial_orient=vec_normalize(initial_orient)
    initial_orient2=vec_normalize(initial_orient2)
    orient= vec_normalize(orient)
    rotation_angle= math.acos(dot_product(initial_orient,orient))
    rotation_axis=cross_product(initial_orient,orient)
    if vec_distance(rotation_axis,[0,0,0]) == 0:
        # if initial_orient and orient are parellel then we  can use any vector perpendicular to this line as the rotation axis 
        rotation_axis=perp_vec(orient)
    rotation_axis= vec_normalize(rotation_axis)
    rotation_matr1= rotation_matrix(rotation_axis,rotation_angle)
    initial_orient2= matrix_vec_multiply(rotation_matr1,initial_orient2)
    orient2= vec_scale(orient2,1 - dot_product(orient2,orient))
    orient2= vec_normalize(orient2)
    rotation_angle= math.acos(dot_product(initial_orient2,orient2))
    rotation_axis= cross_product(initial_orient2,orient2)
    if vec_distance(rotation_axis,[0,0,0]) == 0:
        # if initial_orient and orient are parellel then we  can use any vector perpendicular to this line as the rotation axis
        rotation_axis=perp_vec(orient2)
    rotation_axis= vec_normalize(rotation_axis)
    rotation_matr2= rotation_matrix(rotation_axis,rotation_angle)
    rotation_matr= matrix_multiply(rotation_matr2,rotation_matr1)
    
    return rotation_matr


#construct rotation matrix for transformation of given angles
def rotation_byangles(a,b,c):
    rotation_matr1= rotation_matrix([1,0,0],a)
    rotation_matr2= rotation_matrix([0,1,0],b)
    rotation_matr3= rotation_matrix([0,0,1],c)
    rotation_matr= matrix_multiply(rotation_matr2,rotation_matr1)
    rotation_matr= matrix_multiply(rotation_matr3,rotation_matr)

    return rotation_matr


#transforms vector to new by rotation and translation
def vec_transform(initvec,rotation_matrix,translate):
    rotated_coord= matrix_vec_multiply(rotation_matrix,[initvec[0],initvec[1],initvec[2]])
    x= translate[0]+rotated_coord[0]
    y= translate[1]+rotated_coord[1]
    z= translate[2]+rotated_coord[2]
    return [x,y,z]


# use of periodic boundary conditions
def usepbc(pos,pbc):
    newpos = []
    for coord in pos:
        i=pos.index(coord)
        newcoord=coord
        while newcoord < 0:
            newcoord = newcoord + pbc[i]
        while newcoord > pbc[i]:
            newcoord = newcoord - pbc[i] 
        newpos.append(newcoord)
    return newpos

#random vector on sphere
def ranvecsph():
    while True:
        ranvec[0] = 2.0 * random.random() - 1.0
        ranvec[1] = 2.0 * random.random() - 1.0
        ranvec[2] = 2.0 * random.random() - 1.0
        if (dot_product(ranvec,ranvec) <= 1):
            break
    #print("%f %f %f \n", ranvec[0],ranvec[1],ranvec[2]);
    return ranvec

def quat_random(void):
    newaxis = ranvecsph
    cosv = 1.0 - 2.0*random.random()
    if (random.random() <0.5):
        sinv = math.sqrt(1.0 - cosv*cosv)
    else:
        sinv = -math.sqrt(1.0 - cosv*cosv)
    newquat=[cosv,newaxis[0]*sinv,newaxis[1]*sinv,newaxis[2]*sinv]

    return newquat

def quat_create(vec,vc,vs):
    return [vc,vec[0]*vs,vec[1]*vs,vec[2]*vs]

#rotate vector with quaternion
def rotatevec(vec,axis,angle):
    quat=quat_create(axis,math.cos(angle/360.0*math.pi),math.sin(angle/360.0*math.pi))
    t2 =  quat[0] * quat[1]
    t3 =  quat[0] * quat[2]
    t4 =  quat[0] * quat[3]
    t5 = -quat[1] * quat[1]
    t6 =  quat[1] * quat[2]
    t7 =  quat[1] * quat[3]
    t8 = -quat[2] * quat[2]
    t9 =  quat[2] * quat[3]
    t10 = -quat[3] * quat[3]
    newx = 2.0 * ( (t8+t10)*vec[0] + (t6-t4)*vec[1] + (t3+t7)*vec[2] ) + vec[0]
    newy = 2.0 * ( (t4+t6)*vec[0] + (t5+t10)*vec[1] + (t9-t2)*vec[2] ) + vec[1]
    newz = 2.0 * ( (t7-t3)*vec[0] + (t2+t9)*vec[1] + (t5+t8)*vec[2] ) + vec[2]
    return [newx,newy,newz]
