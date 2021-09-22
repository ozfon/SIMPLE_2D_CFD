#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  9 20:29:14 2021

@author: scott
"""

#ME 541
#Project 1

#import libraries---------------------------
import numpy as np
import matplotlib.pyplot as plt
#import xlwings as xw
#Libraries are imported---------------------

#needed for xlwings-------------------------
#wb = xw.Book(r'/Users/scott/Documents/BYU_classes/Winter2021/ME541/Final_Project/SampleOutput_rev_7_Apr_2020.xlsx')
#sht = wb.sheets['Info and Initial Conditions']
#-------------------------------------------

#Define parameters--------------------------
rho  = 1000.#kg/m^3     #Density of water
mu   = 0.001#Pa*s       #viscosity of water
h    = 1/100#m          #height of channel
l    = 5/100#m          #length of channel
w    = 1#m              #arbitrary width of channel
vin  = 0.001#m/s        #uniform inlet velocity
Pout = 0.0#Pa           #uniform outlet pressure
#Relaxation Properties----------------------
alphau = 0.5#u velocity relaxation*****#0.5
alphav = 0.5#v velocity relaxation*****#0.5
alphap = 0.5#p pressure relaxation*****#0.5
#find mass flow rate in
mdotin = rho*vin*h*w
#variable for keeping track of iterations
iterations = 0
#Converging criteria variables
conv_threshold    = 1*10**-8#threshold for iterations to stop******
max_criteria      = 0.99#set initial criteria above conv_threshold
max_allowed_steps = 10000#max steps before iterations is stopped******
#-------------------------------------------

#Define Grid--------------------------------
#setup the nodes
Nux = 40#Number of u velocity nodes in x direction(initially set to 5)****changes grid refinement****
Nuy = Nux-1#number of u veloicty nodes in y direction(initially set to 4)
Nu = np.array([Nuy,Nux])#number of u velocity nodes[y,x]
Nv = np.array([Nuy+1,Nux+1])#number of y velocity nodes [y,x]
Np = np.array([Nuy,Nux-1])#number of p pressure nodes [y,x]
dx = l/(Nu[1]-1)#distance between nodes in x direction
dy = h/(Nv[0]-1)#distance between nodes in y direction
#setup position of u velocity nodes
xu = np.linspace(0,l,Nu[1])
yu = np.linspace(dy/2,h-dy/2,Nu[0])
#setup position of v velocity nodes
xv = np.linspace(-dx/2,l+dx/2,Nv[1])
yv = np.linspace(0,h,Nv[0])
#setup position of p pressure nodes
xp = np.linspace(dx/2,l-dx/2,Np[1])
yp = np.linspace(dy/2,h-dy/2,Np[0])
#-------------------------------------------

#Declare vectors----------------------------
Au = np.zeros((Nu[0]*Nu[1],Nu[0]*Nu[1]))
bu = np.zeros(len(Au))
Av = np.zeros((Nv[0]*Nv[1],Nv[0]*Nv[1]))
bv = np.zeros(len(Av))
Ap = np.zeros(((Np[0]*(Np[1]-1)),(Np[0]*(Np[1]-1))))
bp = np.zeros(len(Ap))
u_new = np.zeros((Nu[0],Nu[1]))
v_new = np.zeros((Nv[0],Nv[1]))
mdotnew = np.zeros(Nu[0])
criteria = np.zeros(3)#[u_crit,v_crit,p_crit]
residualu = np.zeros(len(Au))
residualv = np.zeros(len(Av))
residual = np.zeros(2)
#-------------------------------------------

#STEP A: Setup initial matrices-------------
#u_initial = np.array(sht.range('A17:E20').value)#intial guess for 5,4 nodes
#v_initial = np.array(sht.range('A23:F27').value)#inital guess provided by professor
#p_initial = np.array(sht.range('A30:D33').value)#matrices imported from excel with xlwings
#set the desired initial field for the number of nodes specified above
u_initial = np.full((Nu[0],Nu[1]),vin)#guess vin everywhere for initial guess
v_initial = np.zeros((Nv[0],Nv[1]))
p_initial = np.zeros((Np[0],Np[1]))
#use the desired initial field to populate the velocity and pressure fields initially
u_star = np.array(u_initial[:])#populate first u star matrix
v_star = np.array(v_initial[:])#populate first v star matrix
p_star = np.array(p_initial[:])#populate first p star matrix
#-------------------------------------------

#Begin iterative process--------------------
for j in range(max_allowed_steps):
    
    #Make u velocity A matrix-------------------
    row = -1
    for i in range(len(Au)):#
        if (np.remainder(i,Nu[1]) == 0):#Inlet condition
            row += 1#for statements below to keep track of current row
            Au[i,i] = 1
            
        
        elif (np.remainder(i+1,Nu[1]) == 0):#outlet condition
            Au[i,i] = 1
            Au[i,i-1] = -1
        
        elif (i < Nu[1]):#The first row will be the bottom wall
            Sp = 2*mu*dx/dy#b term(no slip condition) - Part III - eq.B.1a
            De = mu/dx#Part III - eq.B.1g
            Dw = mu/dx#Part III - eq.B.1h
            Dn = mu/dy#Part III - eq.B.1i
            Ds = 2*mu/dy#Part III - eq.B.1j
            Fe = rho/2*(u_star[row,i+1]+u_star[row,i])#Part III - eq.B.1k
            Fw = rho/2*(u_star[row,i]+u_star[row,i-1])#Part III - eq.B.1l
            Fn = rho/2*(v_star[row+1,i+1]+v_star[row+1,i])#Part III - eq.B.1m
            Fs = rho/2*(v_star[row,i+1]+v_star[row,i])#Part III - eq.B.1n
            Au[i,i+1] = -(De*dy+np.max([-Fe,0])*dy)#ae - Part III - eq.B.1b
            Au[i,i-1] = -(Dw*dy+np.max([Fw,0])*dy)#aw - Part III - eq.B.1c
            Au[i,i+Nu[1]] = -(Dn*dx+np.max([-Fn,0])*dx)#an - Part III - eq.B.1d
            Au[i,i] = (-Au[i,i+1] -Au[i,i-1] -Au[i,i+Nu[1]] +(Fe-Fw)*dy +(Fn-Fs)*dx +Sp)/alphau#ap - Part III - eq.B.1f
        
        elif (i < (Nu[1]*Nu[0])-Nu[1]):#Fill in the interior nodes
            column = i-Nu[1]*row
            De = mu/dx#Part III - eq.B.1g
            Dw = mu/dx#Part III - eq.B.1h
            Dn = mu/dy#Part III - eq.B.1i
            Ds = mu/dy#Part III - eq.B.1j
            Fe = rho/2*(u_star[row,column+1]+u_star[row,column])#Part III - eq.B.1k
            Fw = rho/2*(u_star[row,column]+u_star[row,column-1])#Part III - eq.B.1l
            Fn = rho/2*(v_star[row+1,column+1]+v_star[row+1,column])#Part III - eq.B.1m
            Fs = rho/2*(v_star[row,column+1]+v_star[row,column])#Part III - eq.B.1n
            Au[i,i+1] = -(De*dy+np.max([-Fe,0])*dy)#ae - Part III - eq.B.1b
            Au[i,i-1] = -(Dw*dy+np.max([Fw,0])*dy)#aw - Part III - eq.B.1c
            Au[i,i+Nu[1]] = -(Dn*dx+np.max([-Fn,0])*dx)#an - Part III - eq.B.1d
            Au[i,i-Nu[1]] = -(Ds*dx+np.max([Fs,0])*dx)#as - Part III - eq.B.1e
            Au[i,i] = (-Au[i,i+1] -Au[i,i-1] -Au[i,i+Nu[1]] -Au[i,i-Nu[1]] +(Fe-Fw)*dy +(Fn-Fs)*dx)/alphau#ap - Part III - eq.B.1f
        
        else:#The last row or the top wall is left
            column = i-Nu[1]*row
            Sp = 2*mu*dx/dy#b term(no slip condition) - Part III - eq.B.1a
            De = mu/dx#Part III - eq.B.1g
            Dw = mu/dx#Part III - eq.B.1h
            Dn = 2*mu/dy#Part III - eq.B.1i
            Ds = mu/dy#Part III - eq.B.1j
            Fe = rho/2*(u_star[row,column+1]+u_star[row,column])#Part III - eq.B.1k
            Fw = rho/2*(u_star[row,column]+u_star[row,column-1])#Part III - eq.B.1l
            Fn = rho/2*(v_star[row+1,column+1]+v_star[row+1,column])#Part III - eq.B.1m
            Fs = rho/2*(v_star[row,column+1]+v_star[row,column])#Part III - eq.B.1n
            Au[i,i+1] = -(De*dy+np.max([-Fe,0])*dy)#ae - Part III - eq.B.1b
            Au[i,i-1] = -(Dw*dy+np.max([Fw,0])*dy)#aw - Part III - eq.B.1c
            Au[i,i-Nu[1]] = -(Ds*dx+np.max([Fs,0])*dx)#as - Part III - eq.B.1e
            Au[i,i] = (-Au[i,i+1] -Au[i,i-1] -Au[i,i-Nu[1]] +(Fe-Fw)*dy +(Fn-Fs)*dx +Sp)/alphau#ap - Part III - eq.B.1f
    
    #Make u velocity b vector--------------------
    row = -1
    for i in range(len(bu)):
        if (np.remainder(i,Nu[1]) == 0):#Inlet condition
            bu[i] = vin
            row += 1#for statements below to keep track of current row
        
        elif (np.remainder(i+1,Nu[1]) == 0):#outlet condition
            bu[i] = 0
        
        else:#fill in the b vector for everything else
            column = i-Nu[1]*row
            bu[i] = (p_star[row,column-1]-p_star[row,column])*dy+((1-alphau)*Au[i,i])*u_star[row,column]
    #--------------------------------------------
    
    #Calculate new u velocity--------------------
    u_star_new = np.linalg.solve(Au,bu)
    u_star_new.resize((Nu[0],Nu[1]),refcheck=False)#change vector to matrix
    #mass flow correction
    mdotout = 0#reset
    for i in range(Nu[0]):#step through all the rows
        mdotnew[i] = rho*u_star_new[i,-1]*dy*w#integrate mass flow rate over the outlet
        mdotout += mdotnew[i]#sum up the new mass flow rate values at each exit node
    for i in range(Nu[0]):
        u_star_new[i,-1] = u_star_new[i,-1]*(mdotin/mdotout)#calculate new values for last column
    #--------------------------------------------
    
    #Calculate v velocity A matrix---------------
    row = -1
    for i in range(len(Av)):
        if (np.remainder(i,Nv[1]) == 0):#Inlet condition
            Av[i,i] = 1
            row += 1#for statements below to keep track of current row
        
        elif (np.remainder(i+1,Nv[1]) == 0):#outlet condition
            Av[i,i] = 1
            Av[i,i-1] = -1
        
        elif (i < Nv[1]):#First row is the bottom wall
            Av[i,i] = 1#At the bottom wall v velocity is 0
        
        elif (i < Nv[0]*Nv[1]-Nv[1]):#interior nodes
            column = i-Nv[1]*row
            De = mu/dx#Part III - eq.C.1g
            Dw = mu/dx#Part III - eq.C.1h
            Dn = mu/dy#Part III - eq.C.1i
            Ds = mu/dy#Part III - eq.C.1j
            Fe = rho/2*(u_star_new[row,column]+u_star_new[row-1,column])#Part III - eq.C.1k
            Fw = rho/2*(u_star_new[row,column-1]+u_star_new[row-1,column-1])#Part III - eq.C.1l
            Fn = rho/2*(v_star[row,column]+v_star[row+1,column])#Part III - eq.C.1m
            Fs = rho/2*(v_star[row-1,column]+v_star[row,column])#Part III - eq.C.1n
            Av[i,i+1] = -(De*dy+np.max([-Fe,0])*dy)#ae - Part III - eq.C.1b
            Av[i,i-1] = -(Dw*dy+np.max([Fw,0])*dy)#aw - Part III - eq.C.1c
            Av[i,i+Nv[1]] = -(Dn*dx+np.max([-Fn,0])*dx)#an - Part III - eq.C.1d
            Av[i,i-Nv[1]] = -(Ds*dx+np.max([Fs,0])*dx)#as - Part III - eq.C.1e
            Av[i,i] = (-Av[i,i+1] -Av[i,i-1] -Av[i,i+Nv[1]] -Av[i,i-Nv[1]] +(Fe-Fw)*dy +(Fn-Fs)*dx)/alphav#ap - Part III - eq.C.1f
        
        else:#Last row is the top wall
            Av[i,i] = 1#At the top wall v velocity is 0
    #--------------------------------------------
    
    #Calculate v velocity b vector---------------
    row = -1
    for i in range(len(Av)):
        if (np.remainder(i,Nv[1]) == 0):#Inlet condition
            bv[i] = 0
            row += 1#increment the row counter at each new row
        
        elif (np.remainder(i+1,Nv[1]) ==0):#Outlet condition
            bv[i] = 0
        
        elif (i < Nv[1]):#First row is the bottom wall
            bv[i] = 0
        
        elif (i < Nv[0]*Nv[1]-Nv[1]):#interior nodes
            column = i-Nv[1]*row
            bv[i] = (p_star[row-1,column-1]-p_star[row,column-1])*dx+((1-alphav)*Av[i,i])*v_star[row,column]
        
        else:#Last row is the top wall
            bv[i] = 0
    #---------------------------------------------
    
    #Calculate new v velocity---------------------
    v_star_new = np.linalg.solve(Av,bv)
    v_star_new.resize((Nv[0],Nv[1]),refcheck=False)#change vector to matrix
    #---------------------------------------------
    
    #Calculate p pressure A matrix----------------
    row = -1
    for i in range(len(Ap)):#Step through all nodes
        if (i < Np[1]-1):#If statement for bottom row***
            if (np.remainder(i,Np[1]-1) == 0):#Inlet condition
                row += 1
                column = 0#first column
                de = dy/Au[column+1+Nu[1]*(row),column+1+Nu[1]*(row)]
                dw = dy/Au[column+Nu[1]*(row),column+Nu[1]*(row)]*alphau
                dn = dx/Av[column+1+Nv[1]*(row+1),column+1+Nv[1]*(row+1)]
                ds = dx/Av[column+1+Nv[1]*(row),column+1+Nv[1]*(row)]*alphav
                Ap[i,i+1] = -(rho*de*dy)#ae - Part III - eq.D.1h
                Apw = -(rho*dw*dy)
                Ap[i,i+Np[1]-1] = -(rho*dn*dx)
                Aps = -(rho*ds*dx)
                Ap[i,i] = -Ap[i,i+1] -Ap[i,i+Np[1]-1] -Apw -Aps
            elif (np.remainder(i+1,Np[1]-1) == 0):#outlet condition
                column = i#first row so column=i
                de = dy/Au[column+1+Nu[1]*(row),column+1+Nu[1]*(row)]
                dw = dy/Au[column+Nu[1]*(row),column+Nu[1]*(row)]
                dn = dx/Av[column+1+Nv[1]*(row+1),column+1+Nv[1]*(row+1)]
                ds = dx/Av[column+1+Nv[1]*(row),column+1+Nv[1]*(row)]*alphav
                Ape = -(rho*de*dy)#ae - Part III - eq.D.1h
                Ap[i,i-1] = -(rho*dw*dy)
                Ap[i,i+Np[1]-1] = -(rho*dn*dx)
                Aps = -(rho*ds*dx)
                Ap[i,i] = -Ape -Ap[i,i+Np[1]-1] -Ap[i,i-1] -Aps
            
            else:#interior nodes
                column = i#first row so column=i
                de = dy/Au[column+1+Nu[1]*(row),column+1+Nu[1]*(row)]
                dw = dy/Au[column+Nu[1]*(row),column+Nu[1]*(row)]
                dn = dx/Av[column+1+Nv[1]*(row+1),column+1+Nv[1]*(row+1)]
                ds = dx/Av[column+1+Nv[1]*(row),column+1+Nv[1]*(row)]*alphav
                Ap[i,i+1] = -(rho*de*dy)#ae - Part III - eq.D.1h
                Ap[i,i-1] = -(rho*dw*dy)
                Ap[i,i+Np[1]-1] = -(rho*dn*dx)
                Aps = -(rho*ds*dx)
                Ap[i,i] = -Ap[i,i+1] -Ap[i,i+Np[1]-1] -Ap[i,i-1] -Aps
            
        elif (i < len(Ap)-(Np[1]-1)):#Middle rows(excludes top and bottom)
            if (np.remainder(i,Np[1]-1) == 0):#Inlet condition
                row += 1
                column = 0#first column
                de = dy/Au[column+1+Nu[1]*(row),column+1+Nu[1]*(row)]
                dw = dy/Au[column+Nu[1]*(row),column+Nu[1]*(row)]*alphau
                dn = dx/Av[column+1+Nv[1]*(row+1),column+1+Nv[1]*(row+1)]
                ds = dx/Av[column+1+Nv[1]*(row),column+1+Nv[1]*(row)]
                Ap[i,i+1] = -(rho*de*dy)#ae - Part III - eq.D.1h
                Apw = -(rho*dw*dy)
                Ap[i,i+Np[1]-1] = -(rho*dn*dx)
                Ap[i,i-(Np[1]-1)] = -(rho*ds*dx)
                Ap[i,i] = -Ap[i,i+1] -Ap[i,i+Np[1]-1] -Apw -Ap[i,i-(Np[1]-1)]
            
            elif (np.remainder(i+1,Np[1]-1) == 0):#outlet condition
                column = i-(Np[1]-1)*row#first column
                de = dy/Au[column+1+Nu[1]*(row),column+1+Nu[1]*(row)]
                dw = dy/Au[column+Nu[1]*(row),column+Nu[1]*(row)]
                dn = dx/Av[column+1+Nv[1]*(row+1),column+1+Nv[1]*(row+1)]
                ds = dx/Av[column+1+Nv[1]*(row),column+1+Nv[1]*(row)]
                Ape = -(rho*de*dy)#ae - Part III - eq.D.1h
                Ap[i,i-1] = -(rho*dw*dy)
                Ap[i,i+Np[1]-1] = -(rho*dn*dx)
                Ap[i,i-(Np[1]-1)] = -(rho*ds*dx)
                Ap[i,i] = -Ape -Ap[i,i+Np[1]-1] -Ap[i,i-1] -Ap[i,i-(Np[1]-1)]
            
            else:#interior nodes
                column = i-(Np[1]-1)*row#first row so column=i
                de = dy/Au[column+1+Nu[1]*(row),column+1+Nu[1]*(row)]
                dw = dy/Au[column+Nu[1]*(row),column+Nu[1]*(row)]
                dn = dx/Av[column+1+Nv[1]*(row+1),column+1+Nv[1]*(row+1)]
                ds = dx/Av[column+1+Nv[1]*(row),column+1+Nv[1]*(row)]
                Ap[i,i+1] = -(rho*de*dy)#ae - Part III - eq.D.1h
                Ap[i,i-1] = -(rho*dw*dy)
                Ap[i,i+Np[1]-1] = -(rho*dn*dx)
                Ap[i,i-(Np[1]-1)] = -(rho*ds*dx)
                Ap[i,i] = -Ap[i,i+1] -Ap[i,i+Np[1]-1] -Ap[i,i-1] -Ap[i,i-(Np[1]-1)]
            
        else:#Top row
            if (np.remainder(i,Np[1]-1) == 0):#Inlet condition
                row += 1
                column = 0#first column
                de = dy/Au[column+1+Nu[1]*(row),column+1+Nu[1]*(row)]
                dw = dy/Au[column+Nu[1]*(row),column+Nu[1]*(row)]*alphau
                dn = dx/Av[column+1+Nv[1]*(row+1),column+1+Nv[1]*(row+1)]*alphav
                ds = dx/Av[column+1+Nv[1]*(row),column+1+Nv[1]*(row)]
                Ap[i,i+1] = -(rho*de*dy)#ae - Part III - eq.D.1h
                Apw = -(rho*dw*dy)
                Apn = -(rho*dn*dx)
                Ap[i,i-(Np[1]-1)] = -(rho*ds*dx)
                Ap[i,i] = -Ap[i,i+1] -Apn -Apw -Ap[i,i-(Np[1]-1)]
            
            elif (np.remainder(i+1,Np[1]-1) == 0):#outlet condition
                column = i-(Np[1]-1)*row#first column
                de = dy/Au[column+1+Nu[1]*(row),column+1+Nu[1]*(row)]
                dw = dy/Au[column+Nu[1]*(row),column+Nu[1]*(row)]
                dn = dx/Av[column+1+Nv[1]*(row+1),column+1+Nv[1]*(row+1)]*alphav
                ds = dx/Av[column+1+Nv[1]*(row),column+1+Nv[1]*(row)]
                Ape = -(rho*de*dy)#ae - Part III - eq.D.1h
                Ap[i,i-1] = -(rho*dw*dy)
                Apn = -(rho*dn*dx)
                Ap[i,i-(Np[1]-1)] = -(rho*ds*dx)
                Ap[i,i] = -Ape -Apn -Ap[i,i-1] -Ap[i,i-(Np[1]-1)]
            
            else:#interior nodes
                column = i-(Np[1]-1)*row#first row so column=i
                de = dy/Au[column+1+Nu[1]*(row),column+1+Nu[1]*(row)]
                dw = dy/Au[column+Nu[1]*(row),column+Nu[1]*(row)]
                dn = dx/Av[column+1+Nv[1]*(row+1),column+1+Nv[1]*(row+1)]*alphav
                ds = dx/Av[column+1+Nv[1]*(row),column+1+Nv[1]*(row)]
                Ap[i,i+1] = -(rho*de*dy)#ae - Part III - eq.D.1h
                Ap[i,i-1] = -(rho*dw*dy)
                Apn = -(rho*dn*dx)
                Ap[i,i-(Np[1]-1)] = -(rho*ds*dx)
                Ap[i,i] = -Ap[i,i+1] -Apn -Ap[i,i-1] -Ap[i,i-(Np[1]-1)]
    #---------------------------------------------
    
    #Make p pressure b vector---------------------
    row = -1
    for i in range(len(bp)):
        if (np.remainder(i,Np[1]-1) == 0):#Inlet condition
            row += 1#for statements below to keep track of current row
        column = i-(Np[1]-1)*row
        bp[i] = -(rho*u_star_new[row,column+1]*dy)+(rho*u_star_new[row,column]*dy)-(rho*v_star_new[row+1,column+1]*dx)+(rho*v_star_new[row,column+1]*dx)
    #---------------------------------------------
    
    #Calculate new p prime pressure---------------
    p_prime = np.linalg.solve(Ap,bp)
    #resize vector to matrix and append 0 to the end of each row
    for row in range(Np[0]):#stepping through for the number of rows
        p_prime = np.insert(p_prime,[-1+Np[1]+(row*Np[1])],[0])
    p_prime.resize((Np[0],Np[1]),refcheck=False)
    #---------------------------------------------
    
    #Calculate p new------------------------------
    p_new = p_star+alphap*p_prime
    #---------------------------------------------
    
    #Calculate u new------------------------------
    row = -1
    for i in range(len(Au)):
        if (np.remainder(i,Nu[1]) == 0):#Inlet condition
            row += 1#increment row at each inlet
            column = 0#column is 0 at the inlet
            u_new[row,column] = vin
        
        elif (np.remainder(i+1,Nu[1]) == 0):#outlet condition
            column = i-Nu[1]*row
            u_new[row,column] = u_star_new[row,column]+dy*alphau/Au[i,i]*(p_prime[row,column-1]-0)
        
        else:#interior nodes
            column = i-Nu[1]*row
            u_new[row,column] = u_star_new[row,column]+dy/Au[i,i]*(p_prime[row,column-1]-p_prime[row,column])
    #---------------------------------------------
    
    #Calculate v new------------------------------
    row = -1
    for i in range(len(Av)):
        if (np.remainder(i,Nv[1]) == 0):#Inlet condition
            row += 1#increment row at each inlet
            column = 0#column is 0 at the inlet
            v_new[row,column] = 0
        
        elif (np.remainder(i+1,Nv[1]) == 0):#outlet condition
            column = i-Nv[1]*row
            v_new[row,column] = 0
        
        elif (i < Nv[1]):#Bottom wall
            column = i-Nv[1]*row
            v_new[row,column] = 0
        
        elif (i < Nv[0]*Nv[1]-Nv[1]):#interior nodes
            column = i-Nv[1]*row
            v_new[row,column] = v_star_new[row,column]+dx/Av[i,i]*(p_prime[row-1,column-1]-p_prime[row,column-1])
        
        else:#top wall
            column = i-Nv[1]*row
            v_new[row,column] = 0
    #---------------------------------------------
    
    #Calculate convergence criteria---------------
    #criteria[0] = np.max(abs(u_star-u_new))#u_crit
    #criteria[1] = np.max(abs(v_star-v_new))#v_crit
    #criteria[2] = np.max(abs(p_star-p_new))#p_crit
    #max_criteria = np.max(criteria)
    #---------------------------------------------
    
    #Calculate U residuals------------------------
    row = -1
    for i in range(len(Au)):#u velocity residual
        if (np.remainder(i,Nu[1]) == 0):#Inlet condition
            row += 1
            column = i-Nu[1]*row
            residualu[i] = Au[i,i]*u_star[row,column]-bu[i]
    
        elif (np.remainder(i+1,Nu[1]) == 0):#outlet condition
            column = i-Nu[1]*row
            residualu[i] = (Au[i,i]*u_star[row,column]+Au[i,i-1]*u_star[row,column-1])*(mdotout/mdotin)
    
        elif (i < Nu[1]):#First row with no south node
            column = i-Nu[1]*row
            residualu[i] = Au[i,i]*u_star[row,column]+Au[i,i-1]*u_star[row,column-1]+Au[i,i+1]*u_star[row,column+1]+Au[i,i+Nu[1]]*u_star[row+1,column]-bu[i]
    
        elif (i < (Nu[1]*Nu[0])-Nu[1]):#Calculate unbounded nodes
            column = i-Nu[1]*row
            residualu[i] = Au[i,i]*u_star[row,column]+Au[i,i-1]*u_star[row,column-1]+Au[i,i+1]*u_star[row,column+1]+Au[i,i+Nu[1]]*u_star[row+1,column]+Au[i,i-Nu[1]]*u_star[row-1,column]-bu[i]
    
        else:#last row no north node
            column = i-Nu[1]*row
            residualu[i] = Au[i,i]*u_star[row,column]+Au[i,i-1]*u_star[row,column-1]+Au[i,i+1]*u_star[row,column+1]+Au[i,i-Nu[1]]*u_star[row-1,column]-bu[i]
    residual[0] = np.max(abs(residualu))
    #---------------------------------------------
    
    #Calculate V residuals------------------------
    row = -1
    for i in range(len(Av)):
        if (np.remainder(i,Nv[1]) == 0):#Inlet condition
            row += 1
            column = i-Nv[1]*row
            residualv[i] = Av[i,i]*v_star[row,column]-bv[i]
    
        elif (np.remainder(i+1,Nv[1]) == 0):#outlet condition
            column = i-Nv[1]*row
            residualv[i] = Av[i,i]*v_star[row,column]+Av[i,i-1]*v_star[row,column-1]-bv[i]
    
        elif (i < Nv[1]):#First row with no south node
            column = i-Nv[1]*row
            residualv[i] = Av[i,i]*v_star[row,column]+Av[i,i-1]*v_star[row,column-1]+Av[i,i+1]*v_star[row,column+1]+Av[i,i+Nv[1]]*v_star[row+1,column]-bv[i]
    
        elif (i < (Nv[1]*Nv[0])-Nv[1]):#Calculate unbounded nodes
            column = i-Nv[1]*row
            residualv[i] = Av[i,i]*v_star[row,column]+Av[i,i-1]*v_star[row,column-1]+Av[i,i+1]*v_star[row,column+1]+Av[i,i+Nv[1]]*v_star[row+1,column]+Av[i,i-Nv[1]]*v_star[row-1,column]-bv[i]
    
        else:#last row no north node
            column = i-Nv[1]*row
            residualv[i] = Av[i,i]*v_star[row,column]+Av[i,i-1]*v_star[row,column-1]+Av[i,i+1]*v_star[row,column+1]+Av[i,i-Nv[1]]*v_star[row-1,column]-bv[i]
    residual[1] = np.max(abs(residualv))
    max_residual = np.max(residual)
    #---------------------------------------------
    
        #residualu[i] = Au[i,i]*u_new[row,column]-Au[i]*u_star[i-1]-B[i]
    #for i in range(1,len(Av)-1):#v velocity residual
    #    residualv[i] = Av
    #residual_max = np.max(abs(residual))
    #print('momentum residuals = ' + str(residual_max))
    #--------Done calculating residuals-------------
    
    #Set new values for next iteration (Step G)---
    p_star = np.array(p_new)
    u_star = np.array(u_new)
    v_star = np.array(v_new)
    #---------------------------------------------
    #print('max criteria =',str(max_criteria))
    print('max residual =',str(max_residual))
    iterations += 1#keep track of iterations
    #break for loop if the crietia is met---------
    if (max_residual < conv_threshold):
        break
    #---------------------------------------------
#end of iteration process-------------------------
#Printing iterations info
print('Iterations =',str(iterations))
print('Nux =',str(Nux))
#Calculating bottom wall shear stress
tau = np.zeros(Nu[1])
for i in range(Nu[1]):
    tau[i] += mu*(u_new[1,i]-u_new[0,i])/dy#calculate total shear stress T=mu*du/dy
tau_avg = np.sum(tau)/Nu[1]
print('Tau average =',str(tau_avg),'Pa')
#store outlet u velocity
u_outlet = np.array(u_new[:,-1])
u_midline = np.array(u_new[int(Np[0]/2),:])
#Post-Processing----------------------------------
plt.figure(1)#Pressure contour plot
plt.contourf(xp,yp,p_new)
plt.colorbar()
plt.xlabel('Length (m)')
plt.ylabel('Height (m)')
plt.title('Pressure Plot')
plt.show()

plt.figure(2)#u velocity contour plot
plt.contourf(xu,yu,u_new)
plt.colorbar()
plt.xlabel('Length (m)')
plt.ylabel('Height (m)')
plt.title('U Velocity Plot')
plt.show()

plt.figure(3)#v velocity contour plot
plt.contourf(xv,yv,v_new)
plt.colorbar()
plt.xlabel('Length (m)')
plt.ylabel('Height (m)')
plt.title('V Velocity Plot')
plt.show()

plt.figure(4)#Midline pressure plot
plt.plot(xp,p_new[int(Np[0]/2),:])
plt.xlabel('Length (m)')
plt.ylabel('Pressure (Pa)')
plt.title('Midline Pressure Plot')
plt.show()

plt.figure(5)#Bottom wall shear stress
plt.plot(xu,tau)
plt.xlabel('Length (m)')
plt.ylabel('$\\tau$ - Shear Stress (Pa)')
plt.title('Wall Shear Stress')
plt.show()

plt.figure(6)#Outlet velocity profile
plt.plot(yu,u_outlet)
plt.xlabel('Height (m)')
plt.ylabel('U Velocity (m/s)')
plt.title('Outlet velocity profile')
plt.show()

plt.figure(7)#Midline U velocity plot
plt.plot(xu,u_new[int(Np[0]/2),:])
plt.xlabel('Length (m)')
plt.ylabel('U Velocity (m/s)')
plt.title('Midline U Velocity')
plt.show()