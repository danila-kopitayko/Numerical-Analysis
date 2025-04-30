import subprocess
import os
var('a,b,c,d,e,f,g,h,R_a,R_c,R_e,q_b,q_f,P,K,M,w_0,theta_0,M_0,Q_0,J,E,k,w_h')

k=(q_f-q_b)/(f-b)
n=500

file=open('values.txt','r')

arr=[]

for line in file:
	if line.replace('\n','').split('#')[0].split('=')!=['']:
		arr.append(line.replace('\n','').split('#')[0].split('='))

for line in arr:
	line[1]=float(line[1].strip())
	print(line)
l=25
file.close()





'''w(x)=1/(E*J)*((M_0*x^2)/2+(Q_0*x^3)/2+(R_a*(x-a)^3)/6*heaviside(x-a)-((q_b+k*(x-b))*(x-b)^4)/24*heaviside(x-b)+((q_b+k*(x-b))*(x-f)^4)/24*heaviside(x-f)+(R_c*(x-c)^3)/6*heaviside(x-c)-(P*(x-d)^3)/6*heaviside(x-d)+(R_e*(x-e)^3)/6*heaviside(x-e)+(M*(x-g)^2)/2*heaviside(x-g)+(K*w_h*(x-h)^3)/6*heaviside(x-h))'''


w(x)=w_0+theta_0*x+1/(E*J)*((M_0*x^2)/2+(Q_0*x^3)/2+(R_a*(x-a)^3)/6*heaviside(x-a)-(q_b*(x-b)^4)/24*heaviside(x-b)-k*(x-b)^5/120+(q_f*(x-f)^4)/24*heaviside(x-f)+k*(x-f)^5/120+(R_c*(x-c)^3)/6*heaviside(x-c)-(P*(x-d)^3)/6*heaviside(x-d)+(R_e*(x-e)^3)/6*heaviside(x-e)+(M*(x-g)^2)/2*heaviside(x-g)-(K*w_h*(x-h)^3)/6*heaviside(x-h))


theta(x)=diff(w(x),x)

f(x)=0
m(x)=diff(E*J*theta.substitute_function(dirac_delta,f),x)


stress=diff(m.substitute_function(dirac_delta,f),x)
g(x)=1
stress=stress.substitute_function(dirac_delta,f)

eq0=w.subs(x=0)==0

eq1=stress.subs(x=l)==0

eq2=m.subs(x=l).substitute_function(dirac_delta,f)==0

eq3=w.substitute(x=a)==0 #ищем R_a

eq4=w.substitute(x=c)==0 #ищем R_c

eq5=w.substitute(x=e)==0 #ищем R_e

eq6=w.substitute(x=h)==w_h

eq7=theta.subs(x=0)==0



ans=solve([eq0.subs(a=arr[0][1],b=arr[1][1],c=arr[2][1],d=arr[3][1],e=arr[4][1],f=arr[5][1],g=arr[6][1],h=arr[7][1],E=arr[8][1],J=arr[9][1],q_b=arr[10][1],q_f=arr[11][1],P=arr[12][1],K=arr[13][1],M=arr[14][1]),eq1.subs(a=arr[0][1],b=arr[1][1],c=arr[2][1],d=arr[3][1],e=arr[4][1],f=arr[5][1],g=arr[6][1],h=arr[7][1],E=arr[8][1],J=arr[9][1],q_b=arr[10][1],q_f=arr[11][1],P=arr[12][1],K=arr[13][1],M=arr[14][1]),eq2.subs(a=arr[0][1],b=arr[1][1],c=arr[2][1],d=arr[3][1],e=arr[4][1],f=arr[5][1],g=arr[6][1],h=arr[7][1],E=arr[8][1],J=arr[9][1],q_b=arr[10][1],q_f=arr[11][1],P=arr[12][1],K=arr[13][1],M=arr[14][1]),eq3.subs(a=arr[0][1],b=arr[1][1],c=arr[2][1],d=arr[3][1],e=arr[4][1],f=arr[5][1],g=arr[6][1],h=arr[7][1],E=arr[8][1],J=arr[9][1],q_b=arr[10][1],q_f=arr[11][1],P=arr[12][1],K=arr[13][1],M=arr[14][1]),eq4.subs(a=arr[0][1],b=arr[1][1],c=arr[2][1],d=arr[3][1],e=arr[4][1],f=arr[5][1],g=arr[6][1],h=arr[7][1],E=arr[8][1],J=arr[9][1],q_b=arr[10][1],q_f=arr[11][1],P=arr[12][1],K=arr[13][1],M=arr[14][1]),eq5.subs(a=arr[0][1],b=arr[1][1],c=arr[2][1],d=arr[3][1],e=arr[4][1],f=arr[5][1],g=arr[6][1],h=arr[7][1],E=arr[8][1],J=arr[9][1],q_b=arr[10][1],q_f=arr[11][1],P=arr[12][1],K=arr[13][1],M=arr[14][1]),eq6.subs(a=arr[0][1],b=arr[1][1],c=arr[2][1],d=arr[3][1],e=arr[4][1],f=arr[5][1],g=arr[6][1],h=arr[7][1],E=arr[8][1],J=arr[9][1],q_b=arr[10][1],q_f=arr[11][1],P=arr[12][1],K=arr[13][1],M=arr[14][1]),eq7.subs(a=arr[0][1],b=arr[1][1],c=arr[2][1],d=arr[3][1],e=arr[4][1],f=arr[5][1],g=arr[6][1],h=arr[7][1],E=arr[8][1],J=arr[9][1],q_b=arr[10][1],q_f=arr[11][1],P=arr[12][1],K=arr[13][1],M=arr[14][1])],w_0,theta_0,Q_0,M_0,R_a,R_c,R_e,w_h,solution_dict=True)
#print('TEST\n',ans[0][w_h])

file=open('output.txt','w')

print(ans[0])

i=0
while i<25:
	file.write('w '+str(w.subs(a=arr[0][1],b=arr[1][1],c=arr[2][1],d=arr[3][1],e=arr[4][1],f=arr[5][1],g=arr[6][1],h=arr[7][1],E=arr[8][1],J=arr[9][1],q_b=arr[10][1],q_f=arr[11][1],P=arr[12][1],K=arr[13][1],M=arr[14][1],w_h=RR(ans[0][w_h]),w_0=RR(ans[0][w_0]),theta_0=RR(ans[0][theta_0]),Q_0=RR(ans[0][Q_0]),M_0=RR(ans[0][M_0]),R_a=RR(ans[0][R_a]),R_c=RR(ans[0][R_c]),R_e=RR(ans[0][R_e]),x=round(i,2)))+' '+str(round(i,2))+'\n')
	i+=l/n

i=0
while i<25:
	file.write('theta '+str(theta.subs(a=arr[0][1],b=arr[1][1],c=arr[2][1],d=arr[3][1],e=arr[4][1],f=arr[5][1],g=arr[6][1],h=arr[7][1],E=arr[8][1],J=arr[9][1],q_b=arr[10][1],q_f=arr[11][1],P=arr[12][1],K=arr[13][1],M=arr[14][1],w_h=RR(ans[0][w_h]),w_0=RR(ans[0][w_0]),theta_0=RR(ans[0][theta_0]),Q_0=RR(ans[0][Q_0]),M_0=RR(ans[0][M_0]),R_a=RR(ans[0][R_a]),R_c=RR(ans[0][R_c]),R_e=RR(ans[0][R_e]),x=round(i,2)))+' '+str(round(i,2))+'\n')
	i+=l/n


i=0
while i<25:
	file.write('m '+str(m.subs(a=arr[0][1],b=arr[1][1],c=arr[2][1],d=arr[3][1],e=arr[4][1],f=arr[5][1],g=arr[6][1],h=arr[7][1],E=arr[8][1],J=arr[9][1],q_b=arr[10][1],q_f=arr[11][1],P=arr[12][1],K=arr[13][1],M=arr[14][1],w_h=RR(ans[0][w_h]),w_0=RR(ans[0][w_0]),theta_0=RR(ans[0][theta_0]),Q_0=RR(ans[0][Q_0]),M_0=RR(ans[0][M_0]),R_a=RR(ans[0][R_a]),R_c=RR(ans[0][R_c]),R_e=RR(ans[0][R_e]),x=round(i,2)).substitute_function(heaviside,g))+' '+str(round(i,2))+'\n')
	i+=l/n

'''хевисайд возникает в слагаемом с моментом, приложенном в точке 21. Поэтому заменяем хевисайд на 1, чтобы учесть влияние момента'''

i=0
while i<25:
	file.write('stress '+str(stress.subs(a=arr[0][1],b=arr[1][1],c=arr[2][1],d=arr[3][1],e=arr[4][1],f=arr[5][1],g=arr[6][1],h=arr[7][1],E=arr[8][1],J=arr[9][1],q_b=arr[10][1],q_f=arr[11][1],P=arr[12][1],K=arr[13][1],M=arr[14][1],w_h=RR(ans[0][w_h]),w_0=RR(ans[0][w_0]),theta_0=RR(ans[0][theta_0]),Q_0=RR(ans[0][Q_0]),M_0=RR(ans[0][M_0]),R_a=RR(ans[0][R_a]),R_c=RR(ans[0][R_c]),R_e=RR(ans[0][R_e]),x=round(i,2)).substitute_function(heaviside,g))+' '+str(round(i,2))+'\n')
	i+=l/n
		
'''хевисайд возникает в слагаемом с точечной силой. Поэтому заменяем хевисайд на 1, чтобы учесть влияние сил'''




file.close()


os.system("gcc -c t_new.c")
os.system("gcc -o t_new.exe t_new.o -lm")
subprocess.call(["./t_new.exe"])

exec(open("graphics.py").read())