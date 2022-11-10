"""
This file define the parameters for the Monte Carlo simulation of the Onsager model for the nematic-isotropic phase transition
please read the user guide for more information about the programm and the parameters
you can also contact shoshan.yoav@gmail.com for more information
"""

rods_n = 1000 
rod_l = 5 
iterations = 1000


factor_sd_normal = 0.1 #see description bellow


Lx_arr = [2,4,6,8,10,15,20,25,30,35,40,50,60,70,80,90,100]  
Ly_arr = [2,4,6,8,10,15,20,25,30,35,40,50,60,70,80,90,100] 
"""
IMPORTANT NOTE: Lx_arr and Ly_arr must contain the same number of elements
for the i'th simulation (system with a defined concentration) the system vol will be: 
V = Lx_arr[i]*Ly_arr[i]
the concentration of the system wil be calculated as follows:
C = rods_n/V
the program will iterate over the values of those to arrays together 
and for each iteration will start a the Monte Carlo simulation
"""


