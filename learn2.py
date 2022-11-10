import os
import numpy as np
import itertools as it
import random
from Onsager_MC_Params import *
from scipy import spatial
from scipy.spatial import distance
from math import cos, sin, tan, atan, sqrt, pi
from matplotlib import collections as mc
from matplotlib.pyplot import cm
import matplotlib.pyplot as plt
from collections import deque

def plot_s_vs_c(o_path,x,y):

    fig, ax = plt.subplots()
    plt.xlabel('Concentration [rods/system surface]')
    plt.ylabel('S', rotation = 0)
    plt.title('Order Parameter with Respect to System Concentration\nHard Wall BC')
    plt.scatter(x, y)
#    plt.xlim(0,max(x))
    plt.ylim(-0.5,1.2)
    plt.yticks([-0.5,0,0.2,0.4,0.6,0.8,1])
    ax.set_xscale('log')
    plt.savefig(o_path+'OP_vs_C.png')
    plt.close()
    

def plot_space(o_path,c,s,rods,r,Lx,Ly):
    
    lines = [[(rods[i,0]+r*cos(rods[i,2]),rods[i,1]+r*sin(rods[i,2])),(rods[i,0]-r*cos(rods[i,2]),rods[i,1]-r*sin(rods[i,2]))] for i in range(len(rods))]
    lc = mc.LineCollection(lines, colors=['black' for i in range(len(rods))], linewidths=0.2)
    
    fig, ax = plt.subplots()
    ax.add_collection(lc)
    ax.autoscale()
    ax.margins(0.1)
    plt.xlim(0,Lx)
    plt.ylim(0,Ly)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.xticks(np.linspace(0,Lx,5))
    plt.yticks(np.linspace(0,Ly,5))
    plt.title('Rods Plot, Hard Wall BC\nConcentration = '+str(c)+' ,Order Parameter = '+str(s))
    plt.axis('scaled')
    
#    plt.plot([-0.2, 0-0.2], [0-0.2, Ly+0.2], color='black', linewidth=0.4)
#    plt.plot([0-0.2, Lx+0.2], [0-0.2, 0-0.2], color='black', linewidth=0.4)
#    plt.plot([Lx+0.2, Lx+0.2], [0-0.2, Ly+0.2], color='black', linewidth=0.4)
#    plt.plot([0-0.2, Lx+0.2], [Ly+0.2, Ly+0.2], color='black', linewidth=0.4)
    plt.plot([0, 0], [0, Ly], color='black', linewidth=0.4)
    plt.plot([0, Lx], [0, 0], color='black', linewidth=0.4)
    plt.plot([Lx, Lx], [0, Ly], color='black', linewidth=0.4)
    plt.plot([0, Lx], [Ly, Ly], color='black', linewidth=0.4)
    
    plt.savefig(o_path+'rods_plot_concentration_'+str(c)+'.png')
#    plt.show()
    plt.close()


def plot_heat_map(o_path,op_mat,c,Lx,Ly,hor_seg=20,ver_seg=20):
    
    fig, ax = plt.subplots()
#    xticks = np.arange(0,Lx*(1 + 1/(hor_seg-1)),hor_seg-1)
#    yticks = np.arange(0,Ly*(1 + 1/(ver_seg-1)),ver_seg-1)
    plt.imshow(op_mat, cmap=cm.plasma, interpolation='nearest', extent=[0,Lx,0,Ly], origin='lower')
#    plt.xticks(xticks)
#    plt.yticks(yticks)
    cbar = plt.colorbar()
    cbar.ax.set_ylabel('S', rotation=0)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.axis('scaled')
    plt.title('Order Parameter, Hard Wall BC\nConcentration = '+str(c))
    plt.savefig(o_path+'OP_heat_map_concentration_'+str(c)+'.png')
#    plt.show()
    plt.close()


def plot_relaxation_time(o_path,mean_angle_arr,c):
    
    
    l = len(mean_angle_arr)
    tau_arr = np.arange(0,l)
    average_director_scal_prod = np.zeros(l)
    
    for tau in tau_arr:
        corelations_sum = 0  
        for j in range(l-tau):
            corelations_sum += sqrt(cos(mean_angle_arr[j]-mean_angle_arr[j+tau]))
        average_director_scal_prod[tau] = corelations_sum/(l-tau)
    
    fig, ax = plt.subplots()
    plt.ylabel('<n(i)*n(i+tau)>')
    plt.xlabel('tau (number of iterations)')
    plt.title('Relexation Time Plot, Hard Wall BC\nConcentration='+str(c))
    plt.plot(tau_arr, average_director_scal_prod)
    plt.xlim(0,max(tau_arr))
    plt.ylim(0,1.1)
    plt.yticks(np.linspace(0,1,6))
    plt.savefig(o_path+'relexation_time_'+str(c)+'.png')        
    plt.close()
    
    
def calc_order_parameter(rods):
    s = np.zeros([len(rods),1])
    for i in range(len(rods)):
        frod_angle =  rods[i,2]
        relative_angles = np.delete(rods,i,axis=0)[:,2] - frod_angle
        s[i] = S(relative_angles).mean()
    return s.mean()
    
#def calc_order_parameter(rods):
#    local_director = rods[:,2].mean()
#    s_for_rods = 1.5*(np.cos(rods[:,2] - local_director))**2 - 0.5
#    
#    return s_for_rods.mean()
    
#def calc_order_parameter(rods):
#    local_director = rods[:,2].mean()
#    relative_angles = rods[:,2] - local_director
#    cos_theta_arr = np.cos(relative_angles)
#    return 1.5*(cos_theta_arr**2).mean() - 0.5


def create_OP_mat(rods,Lx,Ly,hor_seg=20,ver_seg=20):
    
    x_segments = np.linspace(0,Lx,hor_seg)
    y_segments = np.linspace(0,Ly,ver_seg)
    op_mat = np.zeros((ver_seg-1,hor_seg-1))
    op_mat[:] = np.NaN
    
    for i in range(ver_seg-1):
        for j in range(hor_seg-1):
            rods_in_seg = rods[np.logical_and(np.logical_and(np.logical_and(rods[:,1]>=y_segments[i],rods[:,1]<=y_segments[i+1]),rods[:,0]>=x_segments[j]),rods[:,0]<=x_segments[j+1])]
            if len(rods_in_seg)>1:
                op_mat[i,j] = calc_order_parameter(rods_in_seg)
                
    return op_mat
            


def test_all_rods_intersect(rods,r):
    
    d = []
    
    for i in range(len(rods)):
        
        rod = rods[i,:]
        new_rods = np.delete(rods,i,axis=0)
        
        x = (new_rods[:,0]*np.tan(new_rods[:,2]) - rod[0]*tan(rod[2]) - new_rods[:,1] + rod[1])/(np.tan(new_rods[:,2])-tan(rod[2]))
        y = tan(rod[2])*(x - rod[0]) + rod[1]
        
        intersections = np.vstack([x,y])
        intersections = np.transpose(intersections)
        
        distances = np.array([[spatial.distance.euclidean(intersections[j,0:2],rod[0:2]),
                               spatial.distance.euclidean(intersections[j,0:2],new_rods[j,0:2])] for j in range(len(intersections))])
        
        on_rays = [distances[k][0]<r and distances[k][1]<r for k in range(len(distances))]
        
        if True in on_rays:
            print('rod ' + str(i) + ' is in problem')
        else:
            print('rod ' + str(i) + ' is OK')
            
        d.append(new_rods[on_rays])
    
#    return d


def test_rods_intersect(rods,rod,r):
      
    x = (rods[:,0]*np.tan(rods[:,2]) - rod[0]*tan(rod[2]) - rods[:,1] + rod[1])/(np.tan(rods[:,2])-tan(rod[2]))
    y = tan(rod[2])*(x - rod[0]) + rod[1]
        
    intersections = np.vstack([x,y])
#    intersections = np.transpose(intersections)
        
#    distances = np.array([[spatial.distance.euclidean(intersections[j,0:2],rod[0:2]),spatial.distance.euclidean(intersections[j,0:2],rods[j,0:2])] for j in range(len(intersections))])
    distances_rod = distance.cdist([rod[0:2].tolist()],intersections[:,0:2].tolist())
    distances_rods = np.sqrt((intersections[:,0]-rods[:,0])**2 + (intersections[:,1]-rods[:,1])**2)
    
    on_rays = [distances_rod[k][0]<r and distances_rods[k][1]<r for k in range(len(distances_rod))]
    mask_rows = [i for i in range(len(on_rays)) if on_rays[i]]    
            
    return rods[mask_rows,:]
    
        
def S(theta_arr):
    return 2*(np.cos(theta_arr))**2 - 1


def initialize_space(rods_n, rod_l,Lx,Ly):    
    
#    broke = False
    #initializing rods locations array
    rods = np.zeros([1,3])
#    print(str(0))
    #determinig first rod position
    rods[0,0:2] = [np.random.uniform(0,Lx),np.random.uniform(0,Ly)]
    first_angle = draw_allowed_angle(np.array([]),rods[0,0:2],rod_l,Lx,Ly)
#    print(first_angle)
    while first_angle==None:
        rods[0,0:2] = [np.random.uniform(0,Lx),np.random.uniform(0,Ly)]
        first_angle = draw_allowed_angle(np.array([]),rods[0,0:2],rod_l,Lx,Ly)
    rods[0,2] = first_angle
    
    #filling rods_locs array
    for i in range(1,rods_n):
#        print(i)
        new_rod_angle = None
        
        #draw location of new rod's center
        while new_rod_angle == None:
            new_rod_loc = [np.random.uniform(0,Lx),np.random.uniform(0,Ly)]
            while any((rods[:,0:2]==new_rod_loc).all(1)): #rod canot be at the same location of any other rod
                new_rod_loc = [np.random.uniform(0,Lx),np.random.uniform(0,Ly)]
            new_rod_angle = draw_allowed_angle(rods[0:i,:],new_rod_loc,rod_l,Lx,Ly)

        new_rod = np.append(new_rod_loc,new_rod_angle)
#        d = test_rods_intersect(rods,new_rod,rod_l/2)
#        if len(d):
#            broke = True
#            break
        rods = np.vstack([rods,new_rod])
#    if broke:  
#        rods = np.vstack([rods,new_rod])
    return rods
        
    
def draw_allowed_angle(rods,new_rod,rod_l,Lx,Ly):

    r = rod_l/float(2)
    walls = []
    #consider close walls (as rods) in allowed angle range calculation
    if new_rod[0] + r > Lx:
        walls.append([Lx, new_rod[1], np.pi/2])
    if new_rod[0] - r < 0:
        walls.append([0, new_rod[1], np.pi/2])
    if new_rod[1] + r > Ly:
        walls.append([new_rod[0], Ly, 0])
    if new_rod[1] - r < 0:
        walls.append([new_rod[0], 0, 0])
    walls = np.array(walls)
    
    #distances array (between current rod and other rods)
    if len(rods):
        dist = distance.cdist([new_rod],rods[:,0:2].tolist())
#        print(dist)
#        print(dist<rod_l)
        #all close rods including close walls
        mask = dist<rod_l
        neighbor_rods = rods[mask[0]]
        
        if len(neighbor_rods) and len(walls):
            neighbor_rods = np.vstack([neighbor_rods,walls])
        elif len(walls):
            neighbor_rods = np.array(walls)
    else: #for first rod
        neighbor_rods = np.array(walls) 
#        print(neighbor_rods)
    allowed_angle_interval = find_angle_bounderies(neighbor_rods,new_rod,r)
    
    if allowed_angle_interval == None:
        return None
    else:
        return map_num_to_intervals(np.random.uniform(),allowed_angle_interval)
    
"""
find allowed angle range for new rod given one neighbor coordinates 
"""
def find_angle_bounderies(rods, new_rod, r):
    
    forbiden_angles = []
    
    if not rods.size:
        return [[0,pi]]
    else:
        #edges locations of neighbor rod
        rays_starts = np.transpose([rods[:,0]+r*np.cos(rods[:,2]),rods[:,1]+r*np.sin(rods[:,2])])
        rays_ends = np.transpose([rods[:,0]-r*np.cos(rods[:,2]),rods[:,1]-r*np.sin(rods[:,2])])

        points1, points2 = circle_rod_intersections(rays_starts, rays_ends, new_rod, r)
    
        #intersection of rods with two intersection points, then allowed angle ranges
        two_intersecting_points1, two_intersecting_points2 = points1[np.array([p[0][0]!=None and p[1][0]!=None for p in zip(points1,points2)])], points2[np.array([p[0][0]!=None and p[1][0]!=None for p in zip(points1,points2)])]
#        print(str(len(two_intersecting_points1)) + ' neighbors with TWO intersection points')
        for i in range(len(two_intersecting_points1)):
            angle_ranges_TOW_intersections = get_range_TWO_intersections(two_intersecting_points1[i,:], two_intersecting_points2[i,:],new_rod,r)
            [forbiden_angles.append(x) for x in angle_ranges_TOW_intersections]
        
        #intersection of rods with one intersection points, then rays edges for corresponding rods, then allowed angle ranges
        one_intersecting_point = points1[np.array([p[0][0]!=None and p[1][0]==None for p in zip(points1,points2)])]
        rays_starts_one_inter, rays_ends_one_inter = rays_starts[np.array([p[0][0]!=None and p[1][0]==None for p in zip(points1,points2)])], rays_ends[np.array([p[0][0]!=None and p[1][0]==None for p in zip(points1,points2)])]
#        print(str(len(one_intersecting_point)) + ' neighbors with ONE intersection points')
        for i in range(len(one_intersecting_point)):
            angle_ranges_ONE_intersection = get_range_ONE_intersection(one_intersecting_point[i,:],rays_starts_one_inter[i,:], rays_ends_one_inter[i,:], new_rod,r)
            [forbiden_angles.append(x) for x in angle_ranges_ONE_intersection]
        
        forbiden_angles_union = find_forbiden_angles_intervals(forbiden_angles)
        allowd_angle = find_complement_intervals(forbiden_angles_union)
        
#        print(allowed_angles)
        if not len(forbiden_angles):
            return [[0,pi]]
        else:
            return allowd_angle
         
"""
given an array of rays defined by 2 points and a circle defined by center and radius
calculate intersction points relative to rays start point
"""  
def circle_rod_intersections(rays_starts,rays_ends,rod_loc,r):
    
    points1 = np.empty((np.size(rays_ends,0),2),dtype=object)
    points2 = np.empty((np.size(rays_ends,0),2),dtype=object)
    t = np.empty((np.size(rays_ends,0),2),dtype=object)
    #vector from ray_start to ray_end
    d = rays_ends - rays_starts
    #vector from center of circle to ray start
    f = rays_starts - np.tile(rod_loc,(np.size(d,0),1))
    #solving (t^2) * (d DOT d) + 2t*( f DOT d ) + ( f DOT f - r^2 ) = 0 
    #where t is a parameter fo how far the intersection along the ray 
    #this equation can describe 6 different cases:
    #where is a parameter for how long along the ray the intersection(s) happen
    #        -o->             --|-->  |            |  --|->
    #      ->  o               o ->              | -> |
    #only first 3 cases are valid solutions 
    #also point cannot overlap another point so no need to check for that situation as well
    a = np.transpose(np.einsum('ij,ij->i',d,d))
    b = np.transpose(2*np.einsum('ij,ij->i',f,d))
    c = np.transpose(np.einsum('ij,ij->i',f,f) - r**2)
    discriminant = np.transpose(b**2 - 4*a*c)
    
    #filling solutions for t whrer discriminant >=0
    mask = discriminant>=0
    for i, booli in enumerate(mask):
        if booli:
            t[i,0] = (-b[i]-np.sqrt(discriminant[i]))/(2*a[i])
            t[i,1] = (-b[i]+np.sqrt(discriminant[i]))/(2*a[i])
    
    for i, sol in enumerate(t):
#        print(sol[0])
        if not sol[0] == None:
            if (sol[0] <=0 and 0<=sol[1]<=1): #one intersection. ray starts outside circle or outside
                points1[i,:] = np.array([sol[1]*d[i,0] + rays_starts[i,0],sol[1]*d[i,1] + rays_starts[i,1]])
#                print('poke')
            elif (0<=sol[0]<=1 and sol[1]>1): #one intersection. ray starts inside circle or outside
                points1[i,:] = np.array([sol[0]*d[i,0] + rays_starts[i,0],sol[0]*d[i,1] + rays_starts[i,1]])
#                print('exit wound')
            elif (0<=sol[0] and sol[1]<=1) and sol[0]!=sol[1]: #two intersections
                points1[i,:] = np.array([sol[0]*d[i,0] + rays_starts[i,0],sol[0]*d[i,1] + rays_starts[i,1]])
                points2[i,:] = np.array([sol[1]*d[i,0] + rays_starts[i,0],sol[1]*d[i,1] + rays_starts[i,1]])  
#                print('pass through')
            else: #no intersection some other case (or ray is tangent) 
#                print('sol no good')
                pass
            
    return points1, points2
                
    
"""
given a sorted intervals list, generate an interator of mapped intervals (to [0,1] range)
"""
def map_intervals(intervals):
#    print('in map intervals')
#    print(intervals)
    intervals_sizes = [x[1]-x[0] for x in intervals]
    joint_intervals_size = sum(intervals_sizes)
    lower = 0
    for s in intervals_sizes:
        portion = s/joint_intervals_size
        yield [lower,lower+portion]
        lower = lower+portion
        
    
"""
map a number between 0 to 1 to a uniform segmented distribution (many intervals)
"""
def map_num_to_intervals(num,intervals):
#    print('in map_num_to_intervals')
#    print(intervals)
    #map intervals to range [0,1]
    mapped_intervals = map_intervals(intervals)
    #get relevant mapped interval and its index in intervals list (to map back to original intervals)
    for i, inter in enumerate(mapped_intervals):
        if inter[0]<=num<=inter[1]:
            break
    
    mapped_num = (num-inter[0])*(intervals[i][1]-intervals[i][0])/(inter[1]-inter[0]) + intervals[i][0]
    return mapped_num


"""
given a range and a list of sorted intervals in range
find complement intervals
"""
def find_complement_intervals(intervals,lower = 0, upper = pi):
    
    comp_intervals = []
    current_low_val = lower
    
    if len(intervals):
        for i in range(len(intervals)):
            if intervals[i][0] > upper or intervals[i][1] < lower:
                pass
            elif intervals[i][0]<=current_low_val<=intervals[i][1]:
                current_low_val = intervals[i][1]
            elif current_low_val<=intervals[i][0]:
                comp_intervals.append([current_low_val,intervals[i][0]])
                current_low_val = intervals[i][1]
                
        if current_low_val < upper:
            comp_intervals.append([current_low_val,upper])
        if not len(comp_intervals):
            return None
        else:    
            return comp_intervals
        
    else:
        return [[lower,upper]]    


"""
for multiple angle ranges return forbiden angles intervals as union of ranges
"""
def find_forbiden_angles_intervals(intervals):
    
    sorted_by_lower_bound = sorted(intervals, key=lambda tup: tup[0])
    merged = []

    for higher in sorted_by_lower_bound:
        if not merged:
            merged.append(higher)
        else:
            lower = merged[-1]
            # test for intersection between lower and higher:
            # we know via sorting that lower[0] <= higher[0]
            if higher[0] <= lower[1]:
                upper_bound = max(lower[1], higher[1])
                merged[-1] = [lower[0], upper_bound]  # replace by merged interval
            else:
                merged.append(higher)
    
    return merged


"""
get excluded angle range for rodes intersecting new rod's circle twice
"""
def get_range_TWO_intersections(point1,point2,new_rod,r):
    
    a1 = atan((new_rod[1]-point1[1])/(new_rod[0]-point1[0]))
    a2 = atan((new_rod[1]-point2[1])/(new_rod[0]-point2[0]))
    angle1, angle2 = a1, a2
    if angle1 < 0:
        angle1 = pi+angle1
    if angle2 < 0:
        angle2 = pi+angle2   
    angles = [angle1,angle2]
    
    #determinintd in which squer the points are
    p1_vec = point1 - new_rod
    p2_vec = point2 - new_rod
    
    if a1*a2>0: #two points are in the same squere or opposite squeres
        if p1_vec[1]*p2_vec[1]<0: #opposite squeres
            return [[0,min(angles)],[max(angles),pi]]
        elif p1_vec[1]*p2_vec[1]>0: #same squere
            return [[min(angles),max(angles)]]
    elif a1*a2<0: #two points are in adjusent squers
        if p1_vec[1]*p2_vec[1]>0 and p1_vec[0]*p2_vec[0]<0: #first and second or third and fourth 
            return [[min(angles),max(angles)]]
        elif p1_vec[1]*p2_vec[1]<0 and p1_vec[0]*p2_vec[0]>0: #second and third or first and fourth
            return [[0,min(angles)],[max(angles),pi]]
    

"""
get excluded angle range for rodes intersecting new rod's circle once
"""    
def get_range_ONE_intersection(point,ray_start,ray_end,new_rod,r):
    
    #get the closet rod edge to the new rod center
    distances_of_rods_edges = [distance.euclidean(x,new_rod) for x in [ray_start,ray_end]]
    close_edge = [ray_start,ray_end][np.argmin(distances_of_rods_edges)]
    
    a1 = atan((new_rod[1]-close_edge[1])/(new_rod[0]-close_edge[0]))
    a2 = atan((new_rod[1]-point[1])/(new_rod[0]-point[0]))
    angle1, angle2 = a1, a2
    if angle1 < 0:
        angle1 = pi+angle1
    if angle2 < 0:
        angle2 = pi+angle2
    angles = [angle1,angle2]
    
    #determinintd in which squer the points are
    p1_vec = close_edge - new_rod
    p2_vec = point - new_rod
    
    if a1*a2>0: #two points are in the same squere or opposite squeres
        if p1_vec[1]*p2_vec[1]<0: #opposite squeres
            return [[0,min(angles)],[max(angles),pi]]
        elif p1_vec[1]*p2_vec[1]>0: #same squere
            return [[min(angles),max(angles)]]
    elif a1*a2<0: #two points are in adjusent squers
        if p1_vec[1]*p2_vec[1]>0 and p1_vec[0]*p2_vec[0]<0: #first and second or third and fourth 
            return [[min(angles),max(angles)]]
        elif p1_vec[1]*p2_vec[1]<0 and p1_vec[0]*p2_vec[0]>0: #second and third or first and fourth
            return [[0,min(angles)],[max(angles),pi]]
        
        
def change_rod(rods,rod_to_chng,rod_l,Lx,Ly,sd_normal = 1,draw_only_allowed_angle = True):
    
    r = rod_l/float(2)
    walls = []
    gd_stp = 0 #this argument hleps count excepted steps in function iterating on this function
    
    #draw new location (normal distribution, casue its isotropic) and angle
    new_rod_loc = np.array([np.random.normal(loc = rod_to_chng[0], scale = sd_normal),
                            np.random.normal(loc = rod_to_chng[1], scale = sd_normal)])
#    new_rod_angle = np.random.uniform(low = 0, high = pi)
#    new_rod_angle = draw_allowed_angle(rods,new_rod_loc,rod_l,Lx,Ly)
#    chngd_rod = np.append(new_rod_loc,new_rod_angle)

    #if new location is still inside system
    if Lx>new_rod_loc[0]>0 and Ly>new_rod_loc[1]>0 and not(any((rods[:,0:2]==new_rod_loc).all(1))):
    
        if draw_only_allowed_angle:
            new_rod_angle = draw_allowed_angle(rods,new_rod_loc,rod_l,Lx,Ly)
            if new_rod_angle != None:
                return gd_stp+1, np.append(new_rod_loc,new_rod_angle)
            else:
                return gd_stp, rod_to_chng
            
        else:
            new_rod_angle = np.random.uniform(low = 0, high = pi)
            #consider close walls (as rods) in allowed angle range calculation
            if chngd_rod[0] + r > Lx:
                walls.append([Lx, chngd_rod[1], np.pi/2])
            if chngd_rod[0] - r < 0:
                walls.append([0, chngd_rod[1], np.pi/2])
            if chngd_rod[1] + r > Ly:
                walls.append([chngd_rod[0], Ly, 0])
            if chngd_rod[1] - r < 0:
                walls.append([chngd_rod[0], 0, 0])
            walls = np.array(walls)
    
            #check for close neighbors
#            dist = np.array([spatial.distance.euclidean(chngd_rod[0:2],rods[i,0:2]) for i in range(len(rods))])
            dist = distance.cdist([rod_to_chng[0:2].tolist()],rods[:,0:2].tolist())
            #all close rods including close walls
            neighbor_rods = rods[dist<rod_l] 
            if len(neighbor_rods) and len(walls):
                neighbor_rods = np.vstack([neighbor_rods,walls])
            elif len(walls):
                neighbor_rods = np.array(walls)
    
            intersections = test_rods_intersect(neighbor_rods,chngd_rod,r)
    
            if not len(intersections): #rod has been changed
                return gd_stp+1, np.append(new_rod_loc,new_rod_angle)
            else:
                return gd_stp, rod_to_chng
    
    else:
        return gd_stp, rod_to_chng
    
"""
this function iterat over all rods
for each rod it calls change_rod to try and change rod location and orientation
returns new rods array and number of rods for which change was valid and therefor exepted
"""
def move_all_rods(rods,rod_l,Lx,Ly,sd_normal = 1):
    
    good_steps = 0
    rods_ind = np.arange(len(rods))
    indeces = deque(range(len(rods)))

    for i in range(len(rods)): 
        
        rol_i_of_indeces_left = random.choice(range(len(indeces)))
        rod_to_chng_ind = indeces[rol_i_of_indeces_left]
        del(indeces[rol_i_of_indeces_left])
#        print(str(i))
#        print(str(len(rods)))
        rod_to_chng = rods[rod_to_chng_ind] #rod to change
        other_rods = rods[rods_ind != rod_to_chng_ind,:] #all other rods
        gd_stp, chngd_rod = change_rod(other_rods,rod_to_chng,rod_l,Lx,Ly,sd_normal = sd_normal) #try to change rod
        good_steps+=gd_stp #count valid steps
#        print(str(len(other_rods[0:i]))+','+str(len(chngd_rod))+','+str(len(other_rods[i:])))
        rods = np.vstack([other_rods[0:i],chngd_rod,other_rods[i:]]) #place rod back in place
#        print(str(len(rods)))
    return good_steps, rods
# =============================================================================
# def move_all_rods(rods,rod_l,Lx,Ly,sd_normal = 1):
#     
#     good_steps = 0
#     rods_ind = np.arange(len(rods))
#     
#     for i in range(len(rods)): 
# #        print(str(i))
# #        print(str(len(rods)))
#         rod_to_chng = rods[i] #rod to change
#         other_rods = rods[rods_ind != i,:] #all other rods
#         gd_stp, chngd_rod = change_rod(other_rods,rod_to_chng,rod_l,Lx,Ly,sd_normal = sd_normal) #try to change rod
#         good_steps+=gd_stp #count valid steps
# #        print(str(len(other_rods[0:i]))+','+str(len(chngd_rod))+','+str(len(other_rods[i:])))
#         rods = np.vstack([other_rods[0:i],chngd_rod,other_rods[i:]]) #place rod back in place
# #        print(str(len(rods)))
#     return good_steps, rods
# =============================================================================
        
"""
for iter_n iterations, calls move_all_rods to try to change rods
return new rods array after all iterations 
and mean_angle_arr (mean angle for each of the iteration) and op_arr (order parameter for each iteration)
"""
def iter_steps(iter_n,rods,rod_l,Lx,Ly,sd_normal = 1):
    
    #initialize arrays for calculation after each iteration
    mean_angle_arr = np.zeros(iter_n)
    mean_angle_arr[:] = np.NaN
    op_arr = np.zeros(iter_n)
    op_arr[:] = np.NaN
    
    tot_good_steps = 0
    for i in range(iter_n): 
        good_steps, rods = move_all_rods(rods,rod_l,Lx,Ly,sd_normal = sd_normal)
        tot_good_steps += good_steps
        #fill arrays
        mean_angle = rods[:,2].mean()
        mean_angle_arr[i] = mean_angle
        s = calc_order_parameter(rods)
        op_arr[i] = s
        print('iteration number ' + str(i+1) + ', valid steps: ' + str(good_steps))
        print('order parameter: ' + str(round(s,3)) + ', mean_angle: ' + str(round(mean_angle,3)))
        
    print('total valid steps: ' + str(tot_good_steps))
    
    return rods, mean_angle_arr, op_arr
    
    
if __name__ == '__main__':
    
    
    cwd = os.getcwd()
    if not os.path.exists(cwd+'/onsager_mc_hard_walls_graphs/'):
        os.makedirs(cwd+'/onsager_mc_hard_walls_graphs/')
    output_folder = cwd+'/onsager_mc_hard_walls_graphs/'
    
    r= rod_l/float(2)
    concentrations = len(Lx_arr)
    s_arr = np.zeros(concentrations)
    c_arr = np.zeros(concentrations)

    for i in range(concentrations):    
    
        Lx = Lx_arr[i]
        Ly = Ly_arr[i]
        c = float(rods_n)/(Lx*Ly)
        sd_normal = factor_sd_normal*(Lx+Ly)/2
        
        rods = initialize_space(rods_n,rod_l,Lx,Ly)
        initial_average_angle = np.array([rods[:,2].mean()])
        print('\n')
        print(r'Initial space for C = ' + str(round(c,2)) + ':\nOrder Parameter: ' + str(round(calc_order_parameter(rods),2)) + ', Average director angle: ' + str(round(rods[:,2].mean(),2)))
        print('Starting MC process. each iteration tries to move all rods one after another')
        rods, mean_angle_arr, op_arr = iter_steps(iterations,rods,rod_l,Lx,Ly,sd_normal = sd_normal)
        mean_angle_arr = np.append(initial_average_angle,mean_angle_arr)

#        np.save(output_folder+'rods'+str(round(c,3)),rods)
#        np.save(output_folder+'mean_angle_arr'+str(round(c,3)),mean_angle_arr)
#        np.save(output_folder+'op_arr'+str(round(c,3)),op_arr)           
        s_arr[i] = op_arr[-1]
        c_arr[i] = c
        
        plot_space(output_folder,c,round(op_arr[-1],2),rods,r,Lx,Ly)
        op_mat = create_OP_mat(rods,Lx,Ly)
        plot_heat_map(output_folder,op_mat,c,Lx,Ly)
        plot_relaxation_time(output_folder,mean_angle_arr,c)
        
#    np.save(output_folder+'S',s_arr)
#    np.save(output_folder+'C',c_arr)
    
    plot_s_vs_c(output_folder,c_arr,s_arr)
