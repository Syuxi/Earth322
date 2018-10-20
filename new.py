@@ -0,0 +1,169 @@
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 14 23:24:31 2016

@author: suo
"""
'''
#importing librarys and packages
matplotlib.pyplot for plotting basemap 
package for scientific computing with Python
csv for raw data file
'''
#import matplotlib.pyplot as plt
import numpy as np
import matplotlib.pyplot as plt
import csv 
from mpl_toolkits.basemap import Basemap,cm
'''
Weighted Form of the Obeserve
Create empty lists to store data
'''
observed=[] #store 3 values,lattitudes, longitudes,O3 value
W24hAVG_SurfO3_units= "parts per billion = nmol/mol" #units of the O3 data
#read data 
with open("janfirst.csv") as csvfile:
    #csv.reader() Return a reader object which will iterate over lines in the given csvfile
    reader = csv.reader(csvfile)
    #next() Return the next row of the reader’s iterable object as a list
    next(reader, None)
    #append() extracts elements(lons+lats+O3value) as input and appends them into the list
    for row in reader:
        x,y,z = float(row[5]), float(row[6]),float(row[20])
        #append three elements together
        observed.append([x,y,z])
#Create a class which makes it easiler to store functions that's easier to organize
class coord():
    #init__create empty object with instances customized and passed into the self. parameter.
    #self. so that the object can keep hold of its own data
    def __init__ (self,coordinate=None, grid_x=None, grid_y=None,observed=None):
        '''
        create a theoretical coordinates to store real coordinates 
        #self.ul :upper left coordinate of observed region
        #self.ur :upper right
        #self.ll :lower left
        #self.lr :lower right
        #self.grid_x :number of grids give to latitudes grids =26
        #self.grid_y :number of grids give to longitude grids =60
        #self.observed :OBSERVED data
        '''
        self.ul=coordinate[0][0];self.ur=coordinate[0][1]
        self.ll=coordinate[1][0];self.lr=coordinate[1][1]
        self.grid_x=grid_x; self.grid_y=grid_y
        self.coordinate=coordinate
        self.observed=observed
        
    def result(self, x=None, y=None):
        '''
        x,y(0,1)-->(lat, lon)
        function that divided the theoretical coordinats into numbers we want,
        
        ul-ll: 1, 
        /by grid_x : divided by 26 pieces =1/26
        *x times it with OBSERVED lattitude(1/26* original latitude)
        +ll: add it to the lowest latitude
        
        return the true coordinates put it into the coordinates we created
        '''
        x=(((self.ul-self.ll)/self.grid_x)*x)+self.ll
        y=(((self.ul-self.ur)/self.grid_y)*y)+self.ur
        return (x,y)
    def distance(self,x, y, lon2, lat2):
        '''
        Calculate the great circle distance between two points(in decimal)
        lon1,lats=x,y is the theoretical coordinates.
        lon2,lat2=OBSERVED coordinates.
        
        #np.radians():
            convert decimal degrees to radians 
        #map():
            applies function to every item of iterable yielding the results.
        #radius=6371km
        '''
        #pass the result of x,y from function result
        lon1,lat1=self.result(x,y)
        #np.radians
        #map() 
        lon1, lat1, lon2, lat2 = map(np.radians, [lon1, lat1, lon2, lat2])
        # haversine formula 
        dlon = lon2 - lon1 
        dlat = lat2 - lat1 
        a = np.sin(dlat/2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2)**2
        c = 2. * np.arcsin(np.sqrt(a)) 
        #km=x-xk
        km = 6371. * c
        return km
        
    def z(self,count1=0,count2=0,zdis1=0,zdis2=0,num1=0, num2=0,wk=0):
        beta=-2.5 #parameter concluded
        D=250. #km as a typical size of a urvab areas
        L=500. #km as a typical scale of synoptic meterology that influence surface O3
        olat = [item[0] for item in observed]  
        olon = [item[1] for item in observed]
        oval = [item[2] for item in observed]
        zk=[]
        for xi in range(self.grid_x):#loop over 26 coordinates
            for yi in range(self.grid_y):#loop over 60 coordinates
                xi,yi=self.result(xi,yi)    
                #print (xi,yi)
                count1=0
                zdis1=0
                zdis2=0
                wk=0
                for [olat,olon,oval] in observed:#calculate the part of smaller than D
                    dis=self.distance(olat,olon,xi,yi)
                    #print (dis)
                    if dis<D:
                        count1=count1+1
                        wk_d=np.power(D,beta)
                        z=wk_d*oval
                        zdis1+=z
                        wk+=wk_d                 
                    #calculating the D<X<L part
                    elif dis>=D and dis<=L:
                        wk_dl=np.power(dis, beta)
                        z=wk_dl*oval
                        zdis2+=z
                        wk+=wk_dl

                num1=zdis1
                num2=zdis2
                try:
                    z=(num1+num2)/wk/count1
                except:
                    z=np.NaN
                #store all the weighted sum data
                zk.append(z)
        zk=np.reshape(zk,(self.grid_x,self.grid_y))
        return zk


#creates a new instance of the class and assigns this object to the local variable a.
a=coord(coordinate=[[41.,-124.],[31.,-109.]],grid_x=26,grid_y=60,observed=observed)
zk=a.z()
#经纬度坐标
yy=np.linspace(a.coordinate[0][0],a.coordinate[1][0],a.grid_x)
xx=np.linspace(a.coordinate[0][1],a.coordinate[1][1],a.grid_y)
xx, yy = np.meshgrid(xx,yy)

fig = plt.figure(figsize=(8,8))
ax = fig.add_axes([0.1,0.1,0.8,0.8])
#设定地图，width地图长度，height地图宽度，单位米， lat_0,lon_0地图中心点坐标
m = Basemap(width=2000000,height=2000000,projection='lcc',
                    resolution='c',lat_0=35,lon_0=-121,)
#画海岸线，州边界，国家边界
m.drawcoastlines()
m.drawstates()
m.drawcountries()
#画经纬度曲线
parallels = np.arange(0.,81,10.)
m.drawparallels(parallels,labels=[False,True,True,False])
meridians = np.arange(10.,351.,10.)
m.drawmeridians(meridians,labels=[True,False,False,True])
#将经纬度坐标转换成地图坐标
xx,yy = m(xx,yy)
#投点
m.scatter(xx,yy, c=zk, cmap='seismic_r')
plt.colorbar()
plt.show()
