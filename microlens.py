import numpy as np

def magCurve(times,u0,tE,t0,Ftot,fb):
  fluxes=fluxCurve(times,u0,tE,t0,Ftot,fb)
  mags=-2.5*np.log10(fluxes)
  return mags

def fluxCurve(times,u0,tE,t0,Ftot,fb):

  normalizedTime=(times-t0)/tE
#  print 't0 tE',t0,tE
#  print ' normalizedTime',normalizedTime
#  pause(10)
  u=np.sqrt(normalizedTime**2 + u0**2)
  magnification=(u**2+2)/(u*np.sqrt(u**2+4))

  baselineflux=10.**(-Ftot/2.5)
  fluxes=baselineflux*(magnification*fb + 1.-fb)

  return fluxes

def fluxCurveMOA(times,u0,tE,t0,Fs,Fb):

  normalizedTime=(times-t0)/tE
  u=np.sqrt(normalizedTime**2 + u0**2)
  magnification=(u**2+2)/(u*np.sqrt(u**2+4))

  fluxes=Fb + Fs*magnification

  return fluxes
#_________________________________________________

def calcFluxratio(time,flux,fluxerr,u0,t0,tE):

  normalizedTime=(time-t0)/tE
  u=np.sqrt(normalizedTime**2 + u0**2)
  magnification=(u**2+2)/(u*np.sqrt(u**2+4))
  
  fluxerr2=fluxerr**2

  if 0:
    print
    print 'fluxerr values',fluxerr
    for err in fluxerr:
      if err<=0:
        print 'ERROR: bad flux error',err
    print 'u0 tE',u0,tE
    print 'fluxerr',sum(fluxerr)/len(fluxerr),sum(fluxerr)
    print 'fluxerr2',sum(fluxerr2)/len(fluxerr),sum(fluxerr2)
    print 'mag',sum(magnification)/len(fluxerr),sum(magnification)
    print 'flux',sum(flux)/len(flux),sum(flux)
#    pause(1000)
  
  matrixA=np.zeros((2,2))
  matrixA[0,0]=sum(magnification**2/fluxerr2)
  matrixA[0,1]=sum(magnification/fluxerr2)
  matrixA[1,0]=matrixA[0,1]                     
  matrixA[1,1]=sum(1./fluxerr2)
  
  matrixB=np.zeros(2)
  matrixB[0]=sum(magnification*flux/fluxerr2)
  matrixB[1]=sum(flux/fluxerr2)    

  try:
    (Fs,Fb)=np.linalg.solve(matrixA,matrixB)
    if 0:
      print
      print 'u0, t0, tE',u0,t0,tE
      print 'matrixA',matrixA
      print 'matrixB',matrixB
      print 'GOOD ONE: matrix didnt solve',Fs,Fb
  except:
    if 0:
      print
      print 'u0, t0, tE',u0,t0,tE
      print 'matrixA',matrixA
      print 'matrixB',matrixB
      print 'TROUBLE: matrix didnt solve'
    (Fs,Fb)=(1.,0.)
#    return 666.,666.,666.
#    exit('TROUBLE: matrix didnt solve')
    
  offsets=Fs*magnification + Fb -flux
  chi=offsets/fluxerr
  sumdiff2=sum(chi**2)

  return Fs,Fb,sumdiff2
