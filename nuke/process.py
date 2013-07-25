import random 

def scale(x,mn,mx):
 return (x-mn)/(mx-mn)

def write_out(fn,data):
 a=open(fn,"w")
 print len(data)
 for k in data:
  a.write(" ".join([str(x) for x in k])+"\n")
 a.close()


lines=open("nuclear.asc").read().split("\n")[:-1]
rows=[]

for k in lines:
 columns=k.split()
 
 #extract atomic number
 a=int(columns[0])

 #get the z with added digit
 z_pre=columns[1]

 if(z_pre[-1]!='0'):
  #print "skipping excited state"
  continue

 z=int(columns[1][:-1])
 
 #calculate n from atomic number and z
 n=a-z

 mass_excess=columns[3]
 if(mass_excess[-1]=='#'):
  #print "skipping inferred excess mass"
  continue

 mass_excess=float(columns[3]) 

 std_dev=columns[4]
 std_dev=float(std_dev)


 #pre-processing data
 mass_excess=scale(mass_excess,-91656.1,119600.0)

 z_parity=z%2
 n_parity=n%2

 z=scale(z,0.0,108.0)
 n=scale(n,0.0,159.0)
 rows.append([z,z_parity,n,n_parity,mass_excess])
 #print " ".join([str(k) for k in z,z_parity,n,n_parity,mass_excess])

random.shuffle(rows)
sz=len(rows)
split2=sz-sz/8
split1=split2-sz/8
training=rows[:split1]
validation=rows[split1:split2]
testing=rows[split2:]

write_out("nuke_train.dat",training)
write_out("nuke_valid.dat",validation)
write_out("nuke_test.dat",testing)
