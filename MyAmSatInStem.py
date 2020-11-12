#!/usr/bin/env python
# coding: utf-8

# # AmSat to kindle STEM in India

# ## Overview
# 
# We need a LEO Polar Sun Sync orbit AmSat with FM repeater and few sensors data monitor.
# Which inturn can be tracked and listened to by students at heart by using either a suitable FM receiver
# or better still a SDR reciever.
# 
# ### Digital data comm
# 
# In turn use a simple mechanism to transmit the sensor data like
# 
# a) Simple Analog frequency modulation: Different frequencies for different values
# 
# b) Analog frequency modulated tones for different decimal digits from 0 to 9, maybe each 200 Hz appart.
# 
# b) Binary Digital Keying: Presence or absence of Frequency or Toggle btw 2 different sets of frequencies indicate between 0 and 1. 
# 
# ### Orbit
# 
# The Polar orbit should allow the satellite to cover all parts of the earth.
# 
# If a orbit altitude is reached such that the velocity achieved can ensure that N (1 or more) full orbits can be finished in time T, such that T is 24 hours or its multiples, then the satellite will retrace the same path at the same relative time once every T hours i.e T/24 days. NOTE: Not sure if the LEO orbits is sufficiently stable enough over a long period of atleast few weeks, so that it can be used to have a simple tracking wrt the satellite for that period.
# 
# If a sun synchronous orbit is achived or maintained, then depending on the angle wrt earth-sun orbital plane,
# 
# * if it is 90 degrees, then the satellite could be placed always in day light, so that solar power is always available. However one needs to verify if the temperature of the satellite parts can be maintained well within reasonable operating range, in this case. Is the attitude maintaining rotation good enough or not or conductive heat distribution pipes/paths good enough or ...
# 
# * if it is 0 degrees, then the satellite would be always at around roughly mid day or mid night wrt the places on earth, depending on which side of the earth wrt sun, that place is.
# 
# The satellite may pass over a given ground location between 0 to 2 times under normal conditions for a polar sun synchronous orbit. Ideally this would be twice if the satellite (and its RF system) covers all parts of the earth fully in a day.
# 

# ## NOTE
# 
# Below are few simple minded calculations to get a very very rough initial idea of few things
# 

# In[3]:


import math
import numpy


# ## Earth

# In[4]:


# Gravitational Constant (m^3 / (kg * s^2))
iGravitationalConstant = (6.674 / 10**11)

# Earth Radius (m)
iEarthRadius = 6_400_000

# Earth Mass (kg)
iEarthMass = 5.972 * 10**24

# Earth Circumference
iEarthCircumference = 2 * math.pi * iEarthRadius
print("EarthCircumference:", iEarthCircumference)


# ## LEO

# In[5]:


## LowEarthOrbit

iLeoAltitudeAboveGround = 500_000
iLeoAltitudeFromEarthCenter = iEarthRadius + iLeoAltitudeAboveGround
print("LEO orbit altitude wrt Earth Center:", iLeoAltitudeFromEarthCenter)

# LEO Velocity
#F = (GMm)/(r^2) = ma = m(v^2/r) => GM = r(v^2)
iLeoVelocity = math.sqrt((iGravitationalConstant * iEarthMass)/iLeoAltitudeFromEarthCenter)
print("LEO Velocity:", iLeoVelocity)

# Circumference of the orbit
iLeoCircumference = 2 * math.pi * iLeoAltitudeFromEarthCenter
print("LeoCircumference:", iLeoCircumference)

# Time taken for a single orbit
iLeoOrbitTimeTaken = iLeoCircumference/iLeoVelocity
print("LeoOrbitTimeTaken(mins):", iLeoOrbitTimeTaken/60)

# Num of Siteings in a day
iDayInSecs = 24*60*60
iNumOfSiteingsInADay = iDayInSecs/iLeoOrbitTimeTaken
print("NumOfSiteings in a Day:", iNumOfSiteingsInADay)
iDayInSecs/(iDayInSecs%iLeoOrbitTimeTaken)
#0.146*103


# ## Earth coverage

# In[8]:


# Earth Coverage and Altitude

# Sin(Angle) = OppositeSide/Hypotinus
# Cos(Angle) = AdjacentSide/Hypotinus
# Tan(Angle) = Opposite/Adjacent
# Opposite = Tan(Angle) * Adjacent
# AdjacentSide = Altitude from surface
# Angle = FieldOfView
iSatAltitude = numpy.array([200_000, 500_000, 1_000_000, 36_000_000])
iFieldOfView = 106 # Adjusted to sync NumOfSiteings with EarthCircumference given EarthCoverage Radius @500KM
iFieldOfView = 60
iEarthCoverageRadius = numpy.tan(((iFieldOfView/2)/360)*math.pi*2)*iSatAltitude
print("EarthCoverageRadius[@{}]:{}".format(iFieldOfView, iEarthCoverageRadius))
math.atan(iEarthRadius/(iEarthRadius+36_000_000)) * (360/(math.pi*2))
math.atan(iEarthRadius/(iEarthRadius+500_000)) * (360/(math.pi*2))

# Minimum number of orbits to cover the earths equatorial circumference
iMinNumOfOrbits4Equator=iEarthCircumference/(iEarthCoverageRadius*2*2)
print("Minimum number of orbits required to cover earths equatorial circumference fully:", iMinNumOfOrbits4Equator)


# ## Return to Same point in Given time
# 
# ### Find Orbit Radius given a Orbit Time
# 
# iOrbitTime = iOrbitCircumference/iOrbitVelocity
# 
# iOrbitCircumference = iOrbitTime * iOrbitVelocity
# 
# 2 * Pi * iOrbitRadius = iOrbitTime * sqrt((iGravitationalConstant * iEarthMass)/iOrbitRadius)
# 
# 4 * Pi^2 * iOrbitRadius^2 = iOrbitTime^2 * (iGravitationalConstant * iEarthMass) / iOrbitRadius
# 
# iOrbitRadius^3 = (iOrbitTime^2 * iGravitationalConst * iEarthMass) / (4 * Pi^2)

# In[ ]:





# ## Eclipse due to earth

# In[7]:


# EvenMoreCrudeMath: rough amount of time for which satellite will be eclipsed by earth in a orbit
iHalfOfEarthToLeoCircumference = (0.5*iEarthCircumference)/iLeoCircumference
iMaxEclipseTime = iHalfOfEarthToLeoCircumference * iLeoOrbitTimeTaken
print("Satellite could be Eclipsed by Earth for a max of {} mins per Orbit".format(iMaxEclipseTime/60))
# Actual Eclipse time will be smaller than this, as light bends over (among others...)


# ## Util functions

# In[33]:


# dB wrt milliwatt normally
def to_dBm(inWatts):
    return 10 * math.log10(inWatts/0.001)

# dB wrt isotropic antenna
def to_dBi(inPower):
    return 10 * math.log10(inPower/1)

iLightSpeed = 300_000_000 # meters / sec

def wavelen2freq(waveLen):
    return iLightSpeed/waveLen

def freq2wavelen(freq):
    return iLightSpeed/freq


# ## RF Losses

# In[40]:


# Using Friis transmission formula
# Pr = Pt * Dt * Dr * (Lambda/4*Pi*D)^2
# Pr_dBm = Pt_dBm + Dt_dBi + Dr_dBi + 10*log10((Lambda/4*Pi*D)^2)
# Pr_dBm = Pt_dBm + Dt_dBi + Dr_dBi + 20*log10(Lambda/4*Pi*D)
# Pr_dBm = Pt_dBm + Dt_dBi + Dr_dBi + 20*log10(Lambda) - 20*log10(D) - 20*log10(1/4*Pi)
#
powerTransmitter = 0.5 # Watts; this is what is normal for LEO AmSats
powerTransmitterDBm = to_dBm(powerTransmitter)
print("TransmitPower(dBm):", powerTransmitterDBm)
directivityTransmitter = 1
directivityTransmitterDBi = to_dBi(directivityTransmitter)
directivityReciever = 1 # Change based on antenne system used
directivityRecieverDBi = to_dBi(directivityReciever)
frequencies = numpy.array([145_800_000, 435_000_000, 900_000_000, 2_400_000_000 ])
print("Frequencies:", frequencies)
distances = numpy.array([[500_000],[1_000_000]])


# This doesnt consider any losses through_the/due_to medium of propagation, nor multipath, nor ...
# However as the path loss is the major factor given the distances involved, its used to get a rough feel
def freespace_pathloss(Pt_dBm, Dt_dBi, Dr_dBi, frequencies, distances):
    Pr_dBm = Pt_dBm + Dt_dBi + Dr_dBi + 20*numpy.log10(iLightSpeed/(4*numpy.pi*distances*frequencies))
    return Pr_dBm


powerRecieverDBm = freespace_pathloss(powerTransmitterDBm, directivityTransmitterDBi, directivityRecieverDBi,
                                       frequencies, distances)
print("PowerAtReciever_Rough(dBm):\n",powerRecieverDBm)


# ## ToDo
# Doppler effect
# 
