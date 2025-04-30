import matplotlib.pyplot as plt

file = open("points_and_values.txt","r")


analytic_der=[]
deriv1=[]
deriv2=[]
runge=[]

for line in file:
    for s3 in [s2 for s2 in line.split('\n')]:
        if 'a' in s3:
            for s in [s4 for s4 in s3.split('\t')]:
                if s!='a':
                    analytic_der.append(float(s))
        if 'b' in s3:
            for s in [s4 for s4 in s3.split('\t')]:
                if s!='b':
                    deriv1.append(float(s))
        if 'c' in s3:
            for s in [s4 for s4 in s3.split('\t')]:
                if s!='c':
                    deriv2.append(float(s))
        if 'd' in s3:
            for s in [s4 for s4 in s3.split('\t')]:
                if s!='d':
                    runge.append(float(s))



analytic_der_x=[]
analytic_der_y=[]

deriv1_x=[]
deriv1_y=[]

deriv2_x=[]
deriv2_y=[]

runge_x=[]
runge_y=[]

for i in range(0,len(analytic_der),2):
    analytic_der_x.append(analytic_der[i])
for i in range(1,len(analytic_der),2):
    analytic_der_y.append(analytic_der[i])

for i in range(0,len(deriv1),2):
    deriv1_x.append(deriv1[i])
for i in range(1,len(deriv1),2):
    deriv1_y.append(deriv1[i])

for i in range(0,len(deriv2),2):
    deriv2_x.append(deriv2[i])
for i in range(1,len(deriv2),2):
    deriv2_y.append(deriv2[i])

for i in range(0,len(runge),2):
    runge_x.append(runge[i])
for i in range(1,len(runge),2):
    runge_y.append(runge[i])


plt.plot(analytic_der_x,analytic_der_y,label='Real derivative',color="black")
plt.scatter(deriv1_x,deriv1_y,label='Derivative 1')
plt.scatter(deriv2_x,deriv2_y,label='Derivative 2')
plt.scatter(runge_x,runge_y,label='Runge')


plt.grid(True)
plt.xlim([analytic_der_x[0],analytic_der_x[-1]])
plt.ylim([min(analytic_der_y),max(analytic_der_y)])
plt.xlabel('x')
plt.ylabel('y')
plt.legend()
plt.show()
file.close()   


