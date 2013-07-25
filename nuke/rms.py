import math

x=float(raw_input("test perf?"))
low= -91656.1 
high =119600.0
rng = (high-low)/1000.0

mse = 1.0-x
mse *= rng*rng

print "rms:", math.sqrt(mse)
