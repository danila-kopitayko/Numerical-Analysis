#include <stdlib.h>                      
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <Eigen/Dense>

float func(float argument)
{return sin(argument);}



using namespace std;
using namespace Eigen;

int main(void)
{
float a=-3.0;
float b=3.0;
int K=4; //количество конечных элементов
int N=5;//количество узлов конечного элемента
int M=10;//количество случайных точек на конечном элементе
int L=N*K-K+1;
int P=1500;
float h=fabs(b-a)/(L-1);
float h1=fabs(b-a)/K;
float h_viz = fabs(b-a)/P;
float psi1=0.0;
float psi2=0.0;
float new_size=100*L-99;
float abs_err1=0.0, rel_err1=0.0, abs_err2=0.0, rel_err2=0.0, abs_err_inf=0.0, rel_err_inf=0.0, norm_l1=0.0, norm_l2=0.0, norm_l_inf=0.0;
float psi;
MatrixXf coeff(L,L); VectorXf b_vect(L); VectorXf C(L);
coeff.setZero();
b_vect.setZero();
C.setZero();
float *x,*x_rand,*graphic_points,*denom,*temp,*interpol_values,*x_new;
interpol_values=(float*)calloc((P+1),sizeof(float));
x_new=(float*)malloc(new_size*sizeof(float));
x=(float*)calloc(L,sizeof(float));
x_rand=(float*)malloc(K*M*sizeof(float));
graphic_points=(float*)malloc((P+1)*sizeof(float));
denom=(float*)malloc(N*sizeof(float));
temp=(float*)malloc(N*sizeof(float));
FILE *file1,*file2, *file_rand;


//заполяем массив x, graphic_points
for(int i=0; i<L; i++)//K*(N-1)+1; i++)
{x[i]=a+i*h;}


for(int i=0; i<=P; i++)
{graphic_points[i]=a+h_viz*i;}

//задаём случайные точки на каждом конечном элементе
srand((unsigned int) time(NULL));
for(int k=0; k<K; k++)
{
 float delta = (b-a)/(K*10*N);
 int l=1;
 x_rand[k*M]= ((float) rand() / (float)(RAND_MAX))*(b-a)/K+a+k*(b-a)/K;


 while(l < 10*N)
 {
 x_rand[k*M+l]= ((float) rand() / (float)(RAND_MAX))*(b - a)/K + a + k*(b-a)/K;
 for(int s=0; s<l; s++) 
  {
  if(fabs(x_rand[k*M+s]-x_rand[k*M+l])<delta) {l--;break;}
  }
 l++;
 if(l==M){break;}
 }
}
//сортировка
float tmp;
for(int i=0; i<K*M; i++)
{for(int j=0;j<K*M-1;j++)
 {if(x_rand[j]>x_rand[(j+1)])
  {tmp=x_rand[j];
   x_rand[j]=x_rand[(j+1)];
   x_rand[(j+1)]=tmp;  
  }
 }
} 


//находим знаменатель многочлена Лагранжа
for(int i=0; i<N; i++)
{denom[i]=1.0;
 for(int k=0; k<N; k++)
 {
  if(i!=k) {denom[i] *= (x[i]-x[k]);}
 }                                   
}                                  

/*
int K//количество конечных элементов
int N//количество узлов конечного элемента
int M//количество случайных точек на конечном элементе
*/

for(int segm=0; segm<K; segm++)
{
 for(int k=0; k<M; k++)
 {for(int i=0; i<N; i++)
  {psi=1.0;
   for(int j=0; j<N; j++)
   {if(i!=j){psi*=(x_rand[k+segm*M]-x[j+segm*(N-1)]);}
   }
   temp[i]=psi/denom[i];
   b_vect(i+segm*(N-1))+=temp[i]*func(x_rand[k+segm*M]);
  } 
  for(int i1=0; i1<N; i1++)
  {
   for(int i2=0; i2<N; i2++)
   {
    coeff(i1+segm*(N-1),i2+segm*(N-1))+=temp[i1]*temp[i2];	
   }
  } 
 }
}


C=coeff.lu().solve(b_vect);

for(int segm=0; segm<K; segm++) 
{for(int j=segm*P/K; j<=(segm+1)*P/K; j++)
 {if((j==segm*P/K)&&(j!=0)){j++;}
  for(int i=0; i<N; i++)
  {psi = 1.0/denom[i];
   for(int k=0; k<N; k++)
   {
    if(i!=k){psi*=(graphic_points[j]-x[k+segm*(N-1)]);}
   }
   interpol_values[j]+=(C(i+segm*(N-1))*psi); 
  }
 }
}



file1=fopen("lagrange_points_and_values.txt","w");

for(int i=0; i<=P; i++)
{fprintf(file1,"p\t%f\t%f\t%f\n",graphic_points[i],interpol_values[i],func(graphic_points[i]));
} 

for(int i=0; i<L; i++)//K*N-K+1; i++)
{fprintf(file1,"n\t%f\t%f\n",x[i],func(x[i]));
}


//считаем погрешности
for(int segm=0; segm<K; segm++) //seg==l
{
 float a1=0.0;
 for(int j=segm*M; j<(segm+1)*M; j++)
 {a1=0.0;
  for(int i=0; i<N; i++)
  {psi = 1.0;
   for(int k=0; k<N; k++)
   {
    if(i!=k){psi*=(x_rand[j]-x[k+segm*(N-1)]);}
   }
   a1+=C(i+segm*(N-1))*psi/denom[i];
  }
  abs_err1+=fabs(a1-func(x_rand[j]));
  norm_l1+=fabs(a1);
  abs_err2+=(a1-func(x_rand[j]))*(a1-func(x_rand[j]));
  norm_l2+=sqrt(a1*a1);
  abs_err_inf=fmax(fabs(a1-func(x_rand[j])),abs_err_inf);
  norm_l_inf=fmax(fabs(norm_l_inf),fabs(a1));
 }
}

file2=fopen("errors.txt","w");
fprintf(file2,"RANDOM POINTS GRID\n");
fprintf(file2,"     \tABSOLUTE ERROR\t\tRELATIVE ERROR\n");
fprintf(file2,"l1   \t%e\t\t%e\n",abs_err1,abs_err1/norm_l1);
fprintf(file2,"l2   \t%e\t\t%e\n",sqrt(abs_err2),sqrt(abs_err2)/sqrt(norm_l2));
fprintf(file2,"l_inf\t%e\t\t%e\n",abs_err_inf,abs_err_inf/norm_l_inf);
fprintf(file2,"\n\n");
h=(b-a)/(100*K*(N-1));

for(int i=0; i<new_size; i++)
{x_new[i]=a+i*h;}

for(int segm=0; segm<K; segm++) 
{
 float a1=0.0;
 for(int j=(int)(segm*new_size/K); j<(int)((segm+1)*new_size/K); j++)
 {a1=0.0;
  for(int i=0; i<N; i++)
  {psi = 1.0;
   for(int k=0; k<N; k++)
   {
    if(i!=k){psi*=(x_new[j]-x[k+segm*(N-1)]);}
   }
   a1+=C(i+segm*(N-1))*psi/denom[i];
  }

  abs_err1+=fabs(a1-func(x_new[j]));
  norm_l1+=fabs(a1);
  abs_err2+=(a1-func(x_new[j]))*(a1-func(x_new[j]));
  norm_l2+=sqrt(a1*a1);
  abs_err_inf=fmax(fabs(a1-func(x_new[j])),abs_err_inf);
  norm_l_inf=fmax(fabs(norm_l_inf),fabs(a1));
 }
}


fprintf(file2,"h/100 STEP GRID\n");
fprintf(file2,"     \tABSOLUTE ERROR\t\tRELATIVE ERROR\n");
fprintf(file2,"l1   \t%e\t\t%e\n",abs_err1,abs_err1/norm_l1);
fprintf(file2,"l2   \t%e\t\t%e\n",sqrt(abs_err2),sqrt(abs_err2)/sqrt(norm_l2));
fprintf(file2,"l_inf\t%e\t\t%e\n",abs_err_inf,abs_err_inf/norm_l_inf);


file_rand=fopen("random_points.txt","w");
for(int i=0; i<K*M; i++){fprintf(file_rand,"r\t%f\t%f\n",x_rand[i],func(x_rand[i]));}

for(int i=0; i<K; i++){fprintf(file_rand,"s\t%f\t%f\n",a+i*fabs(b-a)/(K-1),func(a+i*fabs(b-a)/(K-1)));}


fclose(file1);
fclose(file_rand);

free(temp);
free(graphic_points);
free(x);
free(x_new);
free(x_rand);
free(interpol_values);
free(denom);
system("Visualization.py");
return 0;
}