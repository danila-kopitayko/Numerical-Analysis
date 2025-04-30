import matplotlib.pyplot as plt
import numpy as np

file=open('output.txt','r')
file_c=open('output_c.txt','r')

w=[]
theta=[]
m=[]
stress=[]
x=[]

#считываем значения функций
for line in file:
    line_arr=line.split()
    if 'x' in line_arr and '|-->' in line_arr:
        line_arr.remove('x')
        line_arr.remove('|-->')
    if line_arr[0]=='w':
        w.append(float(line_arr[1]))
        x.append(float(line_arr[2]))
    if line_arr[0]=='theta':
        theta.append(float(float(line_arr[1])))
    if line_arr[0]=='m':
        m.append(float(line_arr[1]))
    if line_arr[0]=='stress':
        stress.append(float(line_arr[1])) 

theta_c=[]
q_c=[]
m_c=[]
#stress_c=[]
w_c=[] 
x_c=[]
for line in file_c:
    line_arr=line.split()
    if 'W' in line_arr:
        #print(line_arr)
        line_arr.remove('W')
        w_c.append(float(line_arr[1]))
        x_c.append(float(line_arr[0])) 
    if 'Theta' in line_arr:
        #print(line_arr)
        line_arr.remove('Theta')
        theta_c.append(float(line_arr[1]))
        #x_c.append(float(line_arr[0])) 
    if 'Moment' in line_arr:
        #print(line_arr)
        line_arr.remove('Moment')
        m_c.append(float(line_arr[1]))
        #x_c.append(float(line_arr[0]))
    if 'Q' in line_arr:
        #print(line_arr)
        line_arr.remove('Q')
        q_c.append(float(line_arr[1]))
        #x_c.append(float(line_arr[0])) 
 
           


file.close()
file_c.close()

#считываем входные данные
arr=[]
file=open('values.txt','r')
for line in file:
	if line.replace('\n','').split('#')[0].split('=')!=['']:
		arr.append(line.replace('\n','').split('#')[0].split('='))
a=float(arr[0][1])
b=float(arr[1][1])
c=float(arr[2][1])
d=float(arr[3][1])
e=float(arr[4][1])
f=float(arr[5][1])
g=float(arr[6][1])
h=float(arr[7][1])

file.close()

def plot_graph(title,x_label,y_label,arr1,arr2,arr1_c,arr2_c):
    plt.title(title)
    plt.plot(arr1,arr2,label='IPM')
    plt.plot(arr1_c,arr2_c,label='FEM')
    plt.grid()
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.scatter(0,0,marker='|',color='black')
    plt.scatter([a,c,e],[0,0,0],marker='^',color='blue')
    plt.scatter(d,0,marker='v')
    #plt.scatter([b,f],[0,0],color='red',s=15)
    plt.plot([0,max(x)],[0,0],color='blue')
    plt.plot([b,f],[0,0],color='yellow')
    plt.scatter(g,0)
    plt.scatter(h,0,marker='d')
    plt.xlim(0,max(x))
    plt.legend()
    #plt.figure(figsize=(10,10))


fig=plt.figure()

ax=fig.add_subplot(511)
plt.title('Схема')
ax.plot([0,max(x)],[0,0],color='blue')
ax.grid()
plt.xlabel('x')
plt.ylabel('y')
ax.scatter(0,0,marker='|',color='black')
ax.scatter([a,c,e],[0,0,0],marker='^',color='blue')
ax.scatter(d,0,marker='v')
#plt.scatter([b,f],[0,0],color='red',s=15)
ax.plot([b,f],[0,0],color='yellow')
ax.scatter(g,0)
ax.scatter(h,0,marker='d')
plt.xlim(0,max(x))
#plt.figure(figsize=(10,10))

ax1=fig.add_subplot(512)
plt.title('Перерезывающая сила')
ax1.plot(x,stress,color='blue',label='IPM')
#print(len(x_c),len(q_c))
ax1.plot(x_c,q_c,color='orange',label='FEM')
ax1.plot([0,max(x)],[0,0],color='blue')
ax1.grid()
plt.xlabel('x')
plt.ylabel('Q')
ax1.scatter(0,0,marker='|',color='black')
ax1.scatter([a,c,e],[0,0,0],marker='^',color='blue')
ax1.scatter(d,0,marker='v')
#plt.scatter([b,f],[0,0],color='red',s=15)
ax1.plot([b,f],[0,0],color='yellow')
ax1.scatter(g,0)
ax1.scatter(h,0,marker='d')
plt.xlim(0,max(x))
#fig.legend()
#plt.figure(figsize=(10,10))

ax2=fig.add_subplot(513)
plt.title('Момент')
ax2.plot(x,m,color='blue',label='IPM')
ax2.plot(x_c,m_c,color='orange',label='FEM')
ax2.plot([0,max(x)],[0,0],color='blue')
ax2.grid()
plt.xlabel('x')
plt.ylabel('M')
ax2.scatter(0,0,marker='|',color='black')
ax2.scatter([a,c,e],[0,0,0],marker='^',color='blue')
ax2.scatter(d,0,marker='v')
#plt.scatter([b,f],[0,0],color='red',s=15)
ax2.plot([b,f],[0,0],color='yellow')
ax2.scatter(g,0)
ax2.scatter(h,0,marker='d')
plt.xlim(0,max(x))
#fig.legend()
#plt.figure(figsize=(10,10))

ax3=fig.add_subplot(514)
plt.title('Угол поворота')
ax3.plot(x,theta,color='blue',label='IPM')
ax3.plot(x_c,theta_c,color='orange',label='FEM')
ax3.plot([0,max(x)],[0,0],color='blue')
ax3.grid()
plt.xlabel('x')
plt.ylabel('Theta')
ax3.scatter(0,0,marker='|',color='black')
ax3.scatter([a,c,e],[0,0,0],marker='^',color='blue')
ax3.scatter(d,0,marker='v')
#plt.scatter([b,f],[0,0],color='red',s=15)
ax3.plot([b,f],[0,0],color='yellow')
ax3.scatter(g,0)
ax3.scatter(h,0,marker='d')
plt.xlim(0,max(x))
#fig.legend()
#plt.figure(figsize=(10,10))

ax4=fig.add_subplot(515)
plt.title('Прогиб')
ax4.plot(x,w,color='blue',label='IPM')
ax4.plot(x_c,w_c,color='orange',label='FEM')
ax4.plot([0,max(x)],[0,0],color='blue')
ax4.grid()
plt.xlabel('x')
plt.ylabel('W')
ax4.scatter(0,0,marker='|',color='black')
ax4.scatter([a,c,e],[0,0,0],marker='^',color='blue')
ax4.scatter(d,0,marker='v')
#plt.scatter([b,f],[0,0],color='red',s=15)
ax4.plot([b,f],[0,0],color='yellow')
ax4.scatter(g,0)
ax4.scatter(h,0,marker='d')
plt.xlim(0,max(x))
#fig.legend()
#plt.figure(figsize=(10,10))



plt.figlegend(['a','Заделка','Опора','Сила','Нагрузка','Момент','Пружина','IPM','FEM'])
def onclick_event(event):
    if event.inaxes==ax:
        plt.figure()
        plot_graph('Схема','x','y',[0,max(x)],[0,0])
        plt.show()
    if event.inaxes==ax1:
        plt.figure()
        plot_graph('Перерезывающая сила','x','Q',x,stress,x_c,q_c)
        #plt.plot(x_c, q_c)
        plt.show()
    if event.inaxes==ax2:
        plt.figure()
        plot_graph('Момент','x','M',x,m,x_c,m_c)
        #plt.plot(x_c, m_c)
        plt.show()
    if event.inaxes==ax3:
        plt.figure()
        plot_graph('Угол поворота','x','Theta',x,theta,x_c, theta_c)
        #plt.plot(x_c, theta_c)
        plt.show()
    if event.inaxes==ax4:
        plt.figure()
        plot_graph('Прогиб','x','W',x,w,x_c,w_c)
        #plt.plot(x_c,w_c)
        plt.show()
        
fig.canvas.mpl_connect('button_press_event',onclick_event)
plt.show()
