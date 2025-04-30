import matplotlib.pyplot as plt

file = open("lagrange_points_and_values.txt","r")
file1 = open("random_points.txt","r")

data2=[]
nodes=[]

for line in file:
    for s3 in [s2 for s2 in line.split('\n')]:
        if 'p' in s3:
            for s in [s4 for s4 in s3.split('\t')]:
                if s!='p':
                    data2.append(float(s))
        if 'n' in s3:
            for s in [s4 for s4 in s3.split('\t')]:
                if s!='n':
                    nodes.append(float(s))


data_random=[]
data_segments=[]

for line in file1:
    for s3 in [s2 for s2 in line.split('\n')]:
        if 'r' in s3:
            for s in [s4 for s4 in s3.split('\t')]:
                if s!='r':
                    data_random.append(float(s))
        if 's' in s3:
            for s in [s4 for s4 in s3.split('\t')]:
                if s!='s':
                    data_segments.append(float(s))
        
lagrange_x=[]
lagrange_y=[]
function=[]

node_x=[]
node_y=[]
for i in range(0,len(data2),3):
    lagrange_x.append(data2[i])
for i in range(1,len(data2),3):
    lagrange_y.append(data2[i])
for i in range(2,len(data2),3):
    function.append(data2[i])

for i in range(0,len(nodes),2):
    node_x.append(nodes[i])
for i in range(1,len(nodes),2):
    node_y.append(nodes[i])

x_random=[]
y_random=[]
x_segments=[]
y_segments=[]

for i in range(0,len(data_random),2):
    x_random.append(data_random[i])
for i in range(1,len(data_random),2):
    y_random.append(data_random[i])
for i in range(0,len(data_segments),2):
    x_segments.append(data_segments[i])
for i in range(1,len(data_segments),2):
    y_segments.append(data_segments[i])



plt.plot(lagrange_x,lagrange_y,label='Approximation')
plt.plot(lagrange_x,function,label='Function')
plt.scatter(node_x,node_y,label='Nodes')
plt.scatter(x_random,y_random,label='Random points')
plt.scatter(x_segments,y_segments,marker='x',color='red',s=120,label='Segment Boundaries')


plt.grid(True)
plt.xlim([lagrange_x[0],lagrange_x[-1]])
plt.ylim([min(function),max(function)])
plt.xlabel('x')
plt.ylabel('y')
plt.legend()
plt.show()
file.close()   

