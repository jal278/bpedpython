#read in rows of data
rows=open("australian.dat").read().split("\n")[:-1]
rows=[[float(x) for x in row.split(" ")] for row in rows]
row_map={1:[1,],4:range(1,4),5:range(1,15),6:range(1,10),8:[1,],9:[1,],11:[1,],12:range(1,4),15:[1,]}

col_count=len(rows[0])

#expand singular categorial inputs into multiple separate inputs by class
nn_rows=[]
for r in rows:
 nn_row=[]
 for x in range(col_count):
  if x+1 in row_map:
   to_append=[]
   for val in row_map[x+1]:
    if r[x]==val:
     to_append.append(1)
    else:
     to_append.append(0)
   nn_row+=to_append
  else:
   nn_row.append(r[x])
 nn_rows.append(nn_row)

col_count=len(nn_rows[0])
print col_count
#normalize inputs between 0 and 1
for x in range(col_count):
 minval = min([row[x] for row in nn_rows])
 maxval = max([row[x] for row in nn_rows])
 if(minval==maxval):
  maxval=minval+1
  print x 
 for k in range(len(nn_rows)):
  v=nn_rows[k][x]
  scaled= (v-minval)/(maxval-minval)
  nn_rows[k][x]=scaled
import random

def write_out(fn,data):
 a=open(fn,"w")
 print len(data)
 for k in data:
  a.write(" ".join([str(x) for x in k])+"\n")
 a.close()

random.shuffle(nn_rows)
splits=10
for k in range(splits):
 split=int(float(len(nn_rows))/splits*k)
 split2=int(float(len(nn_rows))/splits*(k+1))
 train=nn_rows[:split]+nn_rows[split2:]
 test=nn_rows[split:split2]
 write_out("aust%d_train.dat"%k,train)
 write_out("aust%d_test.dat"%k,test)
