#!/usr/bin/env python
# coding: utf-8

# # AmSat to kindle STEM in India

# ## Overview
# 
# We need a LEO Polar Sun Sync orbit AmSat with FM repeater and few sensors data monitor.
# Which inturn can be tracked and listened to by students at heart by using either a suitable FM receiver
# or better still a SDR reciever.
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
# 

# In[ ]:


import math


# ## LEO

# In[8]:


## Gravitational Constant (m^3 / (kg * s^2))
iGravitationalConstant = (6.674 / 10**11)

## Earth Radius (m)
iEarthRadius = 6_400_000
## Earth Mass (kg)
iEarthMass = 5.972 * 10**24

## LowEarthOrbit
iLeoAltitudeAboveGround = 500_000
iLeoAltitudeFromEarthCenter = iEarthRadius + iLeoAltitudeAboveGround
print("LEO orbit altitude wrt Earth Center:", iLeoAltitudeFromEarthCenter)

## LEO Velocity
#F = (GMm)/(r^2) = ma = m(v^2/r) => GM = r(v^2)
iLeoVelocity = math.sqrt((iGravitationalConstant * iEarthMass)/iLeoAltitudeFromEarthCenter)
print("LEO Velocity:", iLeoVelocity)


# In[18]:


# Circumference of the orbit
iLeoCircumference = 2 * math.pi * iLeoAltitudeFromEarthCenter
print("LeoCircumference:", iLeoCircumference)
# Time taken for a single orbit
iLeoOrbitTimeTaken = iLeoCircumference/iLeoVelocity
print("LeoOrbitTimeTaken(mins):", iLeoOrbitTimeTaken/60)


# In[48]:


iDayInSecs = 24*60*60
iNumOfSiteingsInADay = iDayInSecs/iLeoOrbitTimeTaken
print("NumOfSiteings in a Day:", iNumOfSiteingsInADay)
iDayInSecs/(iDayInSecs%iLeoOrbitTimeTaken)
#0.146*103


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

# In[49]:


# EvenMoreCrudeMath: rough amount of time for which satellite will be eclipsed by earth in a orbit 
iEarthCircumference = 2 * math.pi * iEarthRadius
iHalfOfEarthToLeoCircumference = (0.5*iEarthCircumference)/iLeoCircumference
iMaxEclipseTime = iHalfOfEarthToLeoCircumference * iLeoOrbitTimeTaken
print("Satellite could be Eclipsed by Earth for a max of {} mins per Orbit".format(iMaxEclipseTime/60))
# Actual Eclipse time will be smaller than this, as light bends over (among others...)

