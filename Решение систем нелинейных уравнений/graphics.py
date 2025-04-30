import matplotlib.pyplot as plt
import pandas as pd

file1 = open("t.txt","r")
x=[]
y=[]
f0=[]
f1=[]
arr=[]
for line in file1:
	arr=((file1.readline().replace('\n','').split('\t')))
	for i in range(0,len(arr)):
		arr[i]=float(arr[i])
	x.append(arr[0])
	y.append(arr[1])
	f0.append(arr[2])
	f1.append(arr[3])


plt.plot(f0,y,label='x=exp[y]')
plt.plot(x,f1,label='y=-sqrt(x)')

plt.grid(True)
plt.xlim([0,2])
plt.ylim([-2,2])
plt.xlabel('x')
plt.ylabel('y')
plt.legend()
plt.show()


fpi=pd.read_csv("fpi.txt",sep=';')

plt.subplot(1,2,1)
plt.plot(fpi['x0'],fpi['norm'])
plt.xlabel('x0')
plt.ylabel('norm')
plt.grid(True)
plt.xlim([min(fpi['x0']),max(fpi['x0'])])
plt.ylim([min(fpi['norm']),max(fpi['norm'])])
plt.title('fixed_point_iteration')
plt.subplot(1,2,2)

plt.plot(fpi['x1'],fpi['norm'])
plt.xlabel('x1')
plt.ylabel('norm')
plt.grid(True)
plt.xlim([min(fpi['x1']),max(fpi['x1'])])
plt.ylim([min(fpi['norm']),max(fpi['norm'])])
plt.show()

newton=pd.read_csv("newton.txt",sep=';')

plt.subplot(1,2,1)
plt.plot(newton['x0'],newton['norm'])
plt.xlabel('x0')
plt.ylabel('norm')
plt.grid(True)
plt.xlim([min(newton['x0']),max(newton['x0'])])
plt.ylim([min(newton['norm']),max(newton['norm'])])
plt.title('newton')
plt.subplot(1,2,2)

plt.plot(newton['x1'],newton['norm'])
plt.xlabel('x1')
plt.ylabel('norm')
plt.grid(True)
plt.xlim([min(newton['x1']),max(newton['x1'])])
plt.ylim([min(newton['norm']),max(newton['norm'])])
plt.show()

newton_diff=pd.read_csv("newton_diff.txt",sep=';')

plt.subplot(1,2,1)
plt.plot(newton_diff['x0'],newton_diff['norm'])
plt.xlabel('x0')
plt.ylabel('norm')
plt.grid(True)
plt.xlim([min(newton_diff['x0']),max(newton_diff['x0'])])
plt.ylim([min(newton_diff['norm']),max(newton_diff['norm'])])
plt.title('newton_diff')
plt.subplot(1,2,2)

plt.plot(newton_diff['x1'],newton_diff['norm'])
plt.xlabel('x1')
plt.ylabel('norm')
plt.grid(True)
plt.xlim([min(newton_diff['x1']),max(newton_diff['x1'])])
plt.ylim([min(newton_diff['norm']),max(newton_diff['norm'])])
plt.show()

newton_mod=pd.read_csv("newton_mod.txt",sep=';')

plt.subplot(1,2,1)
plt.plot(newton_mod['x0'],newton_mod['norm'])
plt.xlabel('x0')
plt.ylabel('norm')
plt.grid(True)
plt.xlim([min(newton_mod['x0']),max(newton_mod['x0'])])
plt.ylim([min(newton_mod['norm']),max(newton_mod['norm'])])
plt.title('newton_mod')
plt.subplot(1,2,2)

plt.plot(newton_mod['x1'],newton_mod['norm'])
plt.xlabel('x1')
plt.ylabel('norm')
plt.grid(True)
plt.xlim([min(newton_mod['x1']),max(newton_mod['x1'])])
plt.ylim([min(newton_mod['norm']),max(newton_mod['norm'])])
plt.show()


file1.close()
