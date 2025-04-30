import matplotlib.pyplot as plt

file_x=open("points_and_values.txt","r")

data2=[]
nodes=[]
coord=[]

for line in file_x:
    for s3 in [s2 for s2 in line.split('\n')]:
        if 'p' in s3:
            for s in [s4 for s4 in s3.split('\t')]:
                if s!='p':
                    data2.append(float(s))
        if 'n' in s3:
            for s in [s4 for s4 in s3.split('\t')]:
                if s!='n':
                    nodes.append(float(s))
        if 'xy' in s3:
            for s in [s4 for s4 in s3.split('\t')]:
                if s!='xy':
                     coord.append(float(s))


lagrange_x=[]
lagrange_y=[]
func_y=[]
for i in range(0,len(data2),3):
    lagrange_x.append(data2[i])
for i in range(1,len(data2),3):
    lagrange_y.append(data2[i])
for i in range(2,len(data2),3):
    func_y.append(data2[i])


x=[]
y=[]

for i in range(0,len(nodes),2):
    x.append(nodes[i])
for i in range(1,len(nodes),2):
    y.append(nodes[i])


plt.plot(lagrange_x,lagrange_y,label="Approximation")
plt.plot(lagrange_x,func_y,label="Function")
plt.grid(True)
plt.xlim([coord[0],coord[1]])
plt.ylim([min(y)-0.5,max(y)+0.5])
plt.xlabel('x')
plt.ylabel('y')
plt.legend()
plt.scatter(x,y)
plt.show()
file_x.close()

