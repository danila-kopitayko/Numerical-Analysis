#include <stdio.h>
#include <stdlib.h>
#include <math.h>

float func(float arg)
{return sin(arg);
}

int main(void)
{
float a=0.0;
float b=10.0;
int K=3;//количество интервалов
int N=5;
int M=N*K;//количество узлов на всем отрезке
int P=500;//количество точек, на которых строим график
float h=fabs(b-a)/(M-1);
float h_err=h/100;
float psi=0.0;
float a1=0.0;
float l1_abs=0.0,l2_abs=0.0,l_inf_abs=0.0,l1_rel=0.0,l2_rel=0.0,l_inf_rel=0.0;
int P_new=fabs(b-a)/h_err; //количество точек, по которым будем находить погрешность
float *x, *f, *denom, *lagrange_pol, *graphic_points,*x_new;
x = (float*)malloc(M*sizeof(float));
f = (float*)malloc(M*sizeof(float));
denom = (float*)malloc(N*sizeof(float));
lagrange_pol = (float*)malloc(P*sizeof(float));
graphic_points = (float*)malloc(P*sizeof(float));
x_new = (float*)malloc(P_new*sizeof(float));
FILE *file1,*error;

//генерируем узлы интерполяции
for(int i=0; i<M; i++)
{x[i]=a+i*h; 
}

for(int i=0; i<P_new; i++)
{x_new[i]=a+i*h_err;
}


for(int i=0; i<P; i++)
{graphic_points[i]=a + (b-a)/(P-1) * i; 
}


//получаем значения функции
for(int i=0; i<M; i++)
{f[i]=func(x[i]); 
}


//находим знаменатель
for(int i=0; i<N; i++)
{denom[i]=1.0;
 for(int k=0; k<N; k++)
 {
  if(i!=k) {denom[i] *= (x[i]-x[k]);}
 }                                   
}                                  

/*
int K//количество интервалов
int N//количество узлов на одном интервале
int M=N*K;//количество узлов на всем отрезке
int P//количество точек, на которых строим график
*/
//строим полином Лагранжа
for(int l=0; l<K; l++)
{
 for(int j=l*P/K; j<(l+1)*P/K; j++)
 {a1=0.0;
  for(int i=l*N; i<(l+1)*N; i++)
  {psi = 1.0;
   for(int k=l*N; k<(l+1)*N; k++)
   {
    if(i!=k){psi*=(graphic_points[j]-x[k]);}
   }
   a1+=f[i]*psi/denom[i-l*N];
  }
  lagrange_pol[j]=a1;
 }
}

file1 = fopen("lagrange_points_and_values.txt","w");

if(file1==NULL){printf("Error -1");return -1;}

for(int j=0; j<P; j++){fprintf(file1,"p\t%f\t%f\t%f\n",graphic_points[j],lagrange_pol[j],func(graphic_points[j]));}
for(int j=0; j<M; j++){fprintf(file1,"n\t%f\t%f\n",x[j],f[j]);}

//int P_new количество точек, по которым будем находить погрешность

for(int l=0; l<K; l++)
{
 for(int j=l*P_new/K; j<(l+1)*P_new/K; j++)
 {a1=0.0;
  for(int i=l*N; i<(l+1)*N; i++)
  {psi = 1.0;
   for(int k=l*N; k<(l+1)*N; k++)
   {
    if(i!=k){psi*=(x_new[j]-x[k]);}
   }
   a1+=f[i]*psi/denom[i-l*N];
  }
  l1_abs+=fabs(func(x_new[j])-a1);
  l1_rel+=fabs(func(x_new[j]));
  l2_abs+=(func(x_new[j])-a1)*(func(x_new[j])-a1);
  l2_rel+=(func(x_new[j]))*(func(x_new[j]));
  l_inf_abs=fmax(l_inf_abs,fabs(func(x_new[j])-a1));
  l_inf_rel=fmax(l_inf_rel,fabs(func(x_new[j])));
 }
}


error=fopen("errors.txt","w");

fprintf(error,"     \tABSOLUTE ERROR\t\tRELATIVE ERROR\n");
fprintf(error,"l1   \t%e\t\t%e\n",l1_abs,l1_abs/l1_rel);
fprintf(error,"l2   \t%e\t\t%e\n",sqrt(l2_abs),sqrt(l2_abs)/sqrt(l2_rel));
fprintf(error,"l_inf\t%e\t\t%e\n",l_inf_abs,l_inf_abs/l_inf_rel);



fclose(error);
fclose(file1);
free(x);
free(graphic_points);
free(x_new);
free(lagrange_pol);
free(denom);
free(f);
system("Visualization.py");
return 0;
}