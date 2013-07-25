

a=open("nuclear.asc").read().split("\n")[:-1]
for k in a:
 columns=k.split()
 a=int(columns[0])
 z_pre=columns[1]
 if(z_pre[-1]!='0'):
  #print "skipping excited state"
  continue

 z=int(columns[1][:-1])
 n=a-z

 mass_excess=columns[3]
 if(mass_excess[-1]=='#'):
  #print "skipping inferred excess mass"
  continue

 mass_excess=float(columns[3]) 
 std_dev=columns[4]
 std_dev=float(std_dev)

 print ",".join([str(k) for k in z,n,mass_excess])
 
