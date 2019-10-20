# -*- coding: utf-8 -*-
from netCDF4 import Dataset
dataset = Dataset('pentad_20111227_v11l35flk.nc.gz.nc4')

"""Read the `uf` data to numpy array."""

uf = -dataset.variables['uwnd'][0].data[::8] * 0.00152597204
vf = dataset.variables['vwnd'][0].data[::8] * 0.00152597204

uf[uf < -50] = 0
uf[uf > 50] = 0

lon = dataset.variables['lon']
lat = dataset.variables['lat']

vf[vf < -50] = 0
vf[vf > 50] = 0

import matplotlib.pyplot as plt
import numpy as np
import random as rd

x = np.linspace(0.125000000, 359.875000, len(uf))
y = np.linspace(-78.3750000, 78.3750000, len(vf))

factor = 1 # 상수

def findNext(latitude, longitude): # 위도, 경도
    global factor, uf, vf, lon, lat
    lat2 = 139 - (grid(latitude, lat[0]) % len(lat))
    lon2 = grid(longitude, lon[0]) % len(lon)
    u2 = uf[lat2][lon2]
    v2 = vf[lat2][lon2]
    if np.isnan(u2) or np.isnan(v2):
        return False
    finPhi = latitude + factor * v2
    finTheta = longitude + factor * u2
    if finTheta < lon[0]:
        amount = lon[0] - finTheta
        finTheta = lon[-1] - amount
    if finPhi < lat[0]:
        amount = lat[0] - finPhi
        finPhi = lat[-1] - amount
    return (finPhi, finTheta) # 위도, 경도

def grid(A, org):
    return int(round(8*(A % 1)) / 8 + np.floor(A) - org) 
    
def findStart():
    global uf, vf
    longitude = rd.uniform(lon[0], lon[-1])
    latitude = rd.uniform(lat[0], lat[-1])
    Theta = grid(latitude, lat[0]) % len(lat)
    Phi = grid(longitude, lon[0]) % len(lon)
    while np.isnan(uf[Theta][Phi]):
        longitude = rd.uniform(20.5, 379.5)
        latitude = rd.uniform(-69.5, 69.5)
        Theta = 139 - (grid(latitude, lat[0]) % 140)
        Phi = grid(longitude, lon[0]) % 360
    return latitude, longitude # 위도, 경도

def euclidD(A, B):
    return pow((A[0] - B[0])**2 + (A[1] - B[1])**2, 0.5)

X, Y = np.meshgrid(x, y)

plt.figure(figsize=(12.8, 9.6))
plt.quiver(X, Y, uf, vf, scale = 5, scale_units = 'x')
'''
for k in range(1):
    prevposition = findStart() # lat & lon
    routeX = [prevposition[1]] # lon
    routeY = [prevposition[0]] # lat
    print(f"Start Point : {prevposition}")
    for i in range(1000) :
        position = findNext(prevposition[0], prevposition[1])
        if not position or euclidD(position, prevposition) <= 1e-5:
            break
        routeY.append(position[0]) # lat
        routeX.append(position[1]) # lon
        prevposition = position
    plt.scatter(routeX, routeY, 0.1)
'''
plt.show()
print("E")

