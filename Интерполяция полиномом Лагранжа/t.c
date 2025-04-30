#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

int main(void)
{
FILE *file;
float a1=0, b1=3.14, probe=0.0,psi=0.0;
int N=15, M=1000, i=0, j=0, k=0;
float h=fabs(b1-a1)/(M-1);
float *x, *values1,*values2, *graphic_points, *denom,*lagrange_pol;
x = (float*)malloc(N*sizeof(float));
values1 = (float*)malloc(N*sizeof(float));
values2 = (float*)malloc(M*sizeof(float));
graphic_points = (float*)malloc(M*sizeof(float));
denom = (float*)malloc(N*sizeof(float)); 
lagrange_pol = (float*)malloc(M*sizeof(float));


//генерируем узлы интерпол€ции

srand((unsigned int) time(NULL));

float delta = (b1-a1)/(10*N);
int l=1;

x[0]= ((float) rand() / (float)(RAND_MAX))*(b1-a1)+a1;
while(l < 10*N)
{
x[l]= ((float) rand() / (float)(RAND_MAX))*(b1-a1)+a1;
for(int k=0; k<l; k++) 
 {
 if((fabs(x[k]-x[l])<delta)) {l--;break;}
 }
l++;
if(l==N){break;}
}

//сортировка
float temp;
for(i=0; i<N; i++)
{for(j=0;j<N-1;j++)
 {if(x[j]>x[j+1])
  {temp=x[j];
   x[j]=x[j+1];
   x[j+1]=temp;  
  }
 }
}



//находим значени€ функции в узлах интерпол€ции
for(i=0; i<N; i++)
 {values1[i]=sin(x[i]);
 }

//находим знаменатель
for(i=0; i<N; i++)
{denom[i]=1;
 for(k=0; k<N; k++)
 {
  if(i!=k) {denom[i] *= (x[i]-x[k]);}
 }
}                                  

//определ€ем массив точек, на которых будем строить график
for(j=0;j<M;j++)
{graphic_points[j]=a1 + h * j;
values2[j]=sin(graphic_points[j]);
}

//строим полином лагранжа

float a=0.0;
for(j=0; j<M; j++)
{a=0.0;
 for(i=0; i<N; i++)
 {psi = 1.0;
  for(k=0; k<N; k++)
  {
   if(i!=k){psi*=(graphic_points[j]-x[k]);}
  }
  a+=values1[i]*psi/denom[i];
 }
 lagrange_pol[j]=a;
}

//записываем все данные в файлы

file = fopen("points_and_values.txt","w");

if(file==NULL){printf("Error -1");return -1;}



for(int i=0; i<M; i++)
{fprintf(file,"p\t%f\t%f\t%f\n",graphic_points[i],lagrange_pol[i],values2[i]);
} 

for(int i=0; i<N; i++)
{fprintf(file,"n\t%f\t%f\n",x[i],values1[i]);
}


fprintf(file,"xy\t%f\t%f\n",a1,b1);

fclose(file);

free(x);
free(values1);
free(values2);                  
free(denom);
free(graphic_points);
free(lagrange_pol);
system("Visualization.py");
return 0;
}