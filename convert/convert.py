import sys
from xml.dom import minidom, Node
def extract_heading(d):
   return float(d.getElementsByTagName('robot_heading')[0].childNodes[0].data)

def extract_specific_point(d,tag):
   return extract_point(d.getElementsByTagName(tag)[0])
def extract_lines(d):
   lines=[]
   line_nodes = d.getElementsByTagName('line')
   for k in line_nodes:
      lines.append(extract_line(k))
   return lines

def extract_poi(d):
   pois=[]
   point_nodes = d.getElementsByTagName('Point')
   for k in point_nodes:
      pois.append(extract_Point(k))
   return pois

def extract_line(d):
   line=[]
   p1 = d.getElementsByTagName('p1')[0]
   p2 = d.getElementsByTagName('p2')[0]
   line+=extract_point(p1)
   line+=extract_point(p2)
   return line

def extract_Point(d):
    x=int(d.getElementsByTagName('X')[0].childNodes[0].data)
    y=int(d.getElementsByTagName('Y')[0].childNodes[0].data)
    return [x,y]

def extract_point(d):
    x=int(d.getElementsByTagName('x')[0].childNodes[0].data)
    y=int(d.getElementsByTagName('y')[0].childNodes[0].data)
    return [x,y]


justwalls=True
start_point=[0,0]
goal_point=[0,0]
poi=[[0,0]]

doc = minidom.parse(sys.argv[1])

lines = extract_lines(doc)
heading = extract_heading(doc)
if(not justwalls):
 poi = extract_poi(doc)

 start_point = extract_specific_point(doc,"start_point")
 goal_point = extract_specific_point(doc,"goal_point")

scale=1.0

print len(lines)
print str(int(scale*start_point[0])) + " " + str(int(scale*start_point[1]))
print str(int(heading))
print str(int(scale*goal_point[0])) + " " + str(int(scale*goal_point[1]))
print str(int(scale*poi[0][0])) + " " + str(int(scale*poi[0][1]))
for line in lines:
    print str(int(scale*line[0])) + " " + str(int(scale*line[1])) + " " + str(int(scale*line[2])) + " " + str(int(scale*line[3]))
