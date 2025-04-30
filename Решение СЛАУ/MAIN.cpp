#include <stdlib.h>                      
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <Eigen/Dense>


using namespace Eigen;

float* cg(float **mat, float *b, int n, float eps)
{
 float *r, *x, *d, *z;
 float alpha, beta, mem1=0.0, mem2=0.0;
 int cnt=0; 
 r=(float*)calloc(n,sizeof(float));
 x=(float*)calloc(n,sizeof(float));
 d=(float*)malloc(n*sizeof(float));
 z=(float*)malloc(n*sizeof(float));

 for(int i=0; i<n; i++)
 {for(int j=0; j<n; j++)
  {
   r[i]+=mat[i][j]*x[j];
  }
  r[i]=b[i]-r[i];
  d[i]=r[i];
 }


 do{
 cnt++;

 for(int i=0; i<n; i++)
 {z[i]=0.0;
  for(int j=0; j<n; j++)
  {
   z[i]+=mat[i][j]*d[j];
  }
 }
 
 mem1=0.0;
 mem2=0.0;
 beta=0.0; 
 for(int i=0; i<n; i++)
 {
  mem1+=r[i]*r[i];
  mem2+=d[i]*z[i];
 } 
 
 alpha=mem1/mem2; 

 for (int i=0; i<n; i++)
 {
  x[i]=x[i]+alpha*d[i];  //printf("r before:%f\n",r[i]);
  r[i]=r[i]-alpha*z[i];  //printf("r after:%f\n",r[i]);
  beta+=r[i]*r[i];
 }

 beta=beta/mem1; //printf("beta=%f mem=%f\n",beta,mem1); 
 
 for(int i=0; i<n; i++)
 {
  d[i]=r[i]+beta*d[i];
 }

 printf("%d %f\n",cnt,sqrt(beta*mem1));
 
 if(cnt>100)
 {printf("error\n"); break;}

}while(sqrt(beta*mem1)>eps);

free(d);
free(z);
return x;
}

float func(float argument)
{return sin(argument);}

//==============================================================================

float *seidel(float **mat, float *b, int n, float eps)
{
 int cnt=0;
 float res, err, mem, mem1;
 float *y, *x;
 y=(float *)malloc(n*sizeof(float));
 x=(float *)calloc(n,sizeof(float));

 do
 {
 cnt++;
 err=0.0;

 for(int i=0; i<n; i++)
 {
  y[i]=b[i]/mat[i][i];
  mem=0.0;
  
  for(int j=0; j<n; j++)
  {
   if(i!=j) 
   {
    mem+=mat[i][j]*x[j];
   }

  }
  x[i]=y[i]-mem/mat[i][i];
 }

 for(int i=0; i<n; i++)
 {mem1=0.0;
  for(int j=0; j<n; j++)
   {mem1+=mat[i][j]*x[j];}
  err+=(mem1-b[i])*(mem1-b[i]);
 }
 printf("%d %f\n",cnt,sqrt(err));

 if(cnt>100){printf("err\n");break;}
}while(sqrt(err)>eps);

free(y);
return x;
} 

float l_1(float *x,float *y, int n_min, int n_max)
{
 
 float abs_err=0.0,rel_err=0.0,err;
 for(int i=0; i<n_min; i++)
 {abs_err+=fabs(x[i]-y[i*n_max/n_min]);
 
  rel_err+=fabs(x[i]);
 }
err=abs_err/rel_err; 
return err;
}

float l_2(float *x,float *y,int n_min, int n_max)
{ 
 float abs_err=0.0, rel_err=0.0, err;
 for(int i=0; i<n_min; i++)
 {abs_err+=pow((x[i]-y[i*n_max/n_min]),2);
  
  rel_err+=pow(x[i],2);} 
err=sqrt(abs_err)/sqrt(rel_err); 
 
return err;
}

float l_inf(float *x,float *y,int n_min, int n_max)
{ 
 float abs_err=0.0,rel_err=0.0,err;
 for(int i=0; i<n_min; i++)
 {abs_err=fmax(fabs(x[i]-y[i*n_max/n_min]),err);
 
  rel_err=fmax(fabs(x[i]),rel_err);
 } 
 err=fabs(abs_err)/fabs(rel_err);

 return err;
}

float l_1_abs(float *x,float *y, int n_min, int n_max)
{
 
 float abs_err=0.0,rel_err=0.0,err;
 for(int i=0; i<n_min; i++)
 {abs_err+=fabs(x[i]-y[i*n_max/n_min]);
 
 // rel_err+=fabs(x[i]);
 }

return abs_err;
}

float l_2_abs(float *x,float *y,int n_min, int n_max)
{ 
 float abs_err=0.0, rel_err=0.0, err;
 for(int i=0; i<n_min; i++)
 {abs_err+=pow((x[i]-y[i*n_max/n_min]),2);}
  


 
return abs_err;
}

float l_inf_abs(float *x,float *y,int n_min, int n_max)
{ 
 float abs_err=0.0,rel_err=0.0,err;
 for(int i=0; i<n_min; i++)
 {abs_err=fmax(fabs(x[i]-y[i*n_max/n_min]),err);
 
  
 } 


 return abs_err;
}


//using namespace std;
//using namespace Eigen;

int main(void)
{
float a=-3.0;
float b=3.0;
int K=3; //количество конечных элементов
int N=4;//количество узлов конечного элемента
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
float **coeff, *b_vect, *C;
//MatrixXf coeff(L,L); VectorXf b_vect(L); VectorXf C(L);
coeff=(float**)malloc(L*sizeof(float*));
for(int i=0; i<L; i++){coeff[i]=(float*)calloc(L,sizeof(float));}

MatrixXf coeff2(L,L); VectorXf b_vect2(L); VectorXf C5(L);
coeff2.setZero();
b_vect2.setZero();
C5.setZero();


b_vect=(float*)calloc(L,sizeof(float));
C=(float*)calloc(L,sizeof(float));
//coeff.setZero();
//b_vect.setZero();
//C.setZero();
float *x,*x_rand,*graphic_points,*denom,*temp,*interpol_values,*x_new;
interpol_values=(float*)calloc((P+1),sizeof(float));
x_new=(float*)malloc(new_size*sizeof(float));
x=(float*)calloc(L,sizeof(float));
x_rand=(float*)malloc(K*M*sizeof(float));
graphic_points=(float*)malloc((P+1)*sizeof(float));
denom=(float*)malloc(N*sizeof(float));
temp=(float*)malloc(N*sizeof(float));
FILE *file1,*file2, *file_rand, *file;


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
   b_vect[i+segm*(N-1)]+=temp[i]*func(x_rand[k+segm*M]);
   b_vect2(i+segm*(N-1))+=temp[i]*func(x_rand[k+segm*M]);

  } 
  for(int i1=0; i1<N; i1++)
  {
   for(int i2=0; i2<N; i2++)
   {
    coeff[i1+segm*(N-1)][i2+segm*(N-1)]+=temp[i1]*temp[i2];
    coeff2(i1+segm*(N-1),i2+segm*(N-1))+=temp[i1]*temp[i2];	
	
   }
  } 
 }
}

C5=coeff2.lu().solve(b_vect2);


for(int i=0; i<L; i++)printf("%f ",C5(i));
printf("\n");



//===================================================================================================================
//conjugate_gradient

float eps;
eps=pow(10,-6);

C=cg(coeff,b_vect,L,eps);
for(int i=0; i<L; i++){printf("C[%d]=%f\t",i,C[i]);}
printf("\n");

//-------------------------------------------------------------------------------------------------------------------
//cholesky
float **L_mat;
float sum;

L_mat=(float**)malloc(L*sizeof(float*));
for(int i=0; i<L; i++){L_mat[i]=(float*)calloc(L,sizeof(float));}

for(int j=0; j<L; j++)
{
 sum=0.0;
 for(int k=0; k<L; k++)
 {
  sum+=L_mat[j][k]*L_mat[j][k];
 }
 L_mat[j][j]=sqrt(coeff[j][j]-sum);
 for(int i=j+1; i<L; i++)
 {sum=0.0;
  for(int k=0; k<j; k++)
  {sum+=L_mat[i][k]*L_mat[j][k];
  }
 L_mat[i][j]=(1.0/L_mat[j][j]*(coeff[i][j]-sum));
 } 
}


float *sol1,*C1;
sol1=(float*)calloc(L,sizeof(float));
C1=(float*)calloc(L,sizeof(float));



for(int i=0; i<L; i++)
{sum=0.0;
 for(int j=0; j<i; j++)
 {sum+=L_mat[i][j]/L_mat[i][i]*sol1[j]; 
 }
 sol1[i]=b_vect[i]/L_mat[i][i]-sum;     
}


for(int i=L-1; i>=0; i--)
{sum=0.0;
 for(int j=L-1; j>i; j--)
 {sum+=L_mat[j][i]/L_mat[i][i]*C1[j]; 
 }
 C1[i]=sol1[i]/L_mat[i][i]-sum; 
}
for(int i=0; i<L; i++){printf("C1[%d]=%f\t",i,C1[i]);}
printf("\n");
for(int i=0; i<L; i++){free(L_mat[i]);}
free(L_mat);
printf("\n");

//-------------------------------------------------------------------------------------------------------------------     
//seidel


float *C2;
C2=(float*)calloc(L,sizeof(float));
C2=seidel(coeff,b_vect,L,eps);


for(int i=0; i<L; i++){printf("C2[%d]=%f\t",i,C2[i]);}
printf("\n");


//-------------------------------------------------------------------------------------------------------------------     
//gauss

float **coeff1, *C3,*b_vect1;
C3=(float*)calloc(L,sizeof(float));
coeff1=(float**)malloc(L*sizeof(float*));
for(int i=0; i<L; i++){coeff1[i]=(float*)calloc(L,sizeof(float));}
b_vect1=(float*)calloc(L,sizeof(float));


for(int i=0; i<L; i++)
{b_vect1[i]=b_vect[i];
 for(int j=0; j<L; j++)
 {
  coeff1[i][j]=coeff[i][j];
 }
} 


for(int segm=0; segm<K; segm++)
{
 for(int l=0; l<N; l++)
 {
if(segm!=K-1)
{float denom=0.0;
 for(int i=l; i<N-1; i++)
 {denom=coeff1[i+segm*(N-1)][l+segm*(N-1)];//mem*=denom;
  for(int j=l; j<N; j++)
  {//printf("%f/%f\n",coeff[i+segm*(N-1)][j+segm*(N-1)],denom);
   coeff1[i+segm*(N-1)][j+segm*(N-1)]=coeff1[i+segm*(N-1)][j+segm*(N-1)]/denom;
  } 
  b_vect1[i+segm*(N-1)]=b_vect1[i+segm*(N-1)]/denom;
 }

 denom=coeff1[N-1+segm*(N-1)][l+segm*(N-1)];//mem*=denom;
 for(int j=l; j<2*N-1; j++)
 {//printf("%f/%f\n",coeff[N-1+segm*(N-1)][j+segm*(N-1)],denom);
  coeff1[N-1+segm*(N-1)][j+segm*(N-1)]=coeff1[N-1+segm*(N-1)][j+segm*(N-1)]/denom;
 } 
 b_vect1[N-1+segm*(N-1)]=b_vect1[N-1+segm*(N-1)]/denom;
}

else
{float denom=0.0;
 for(int i=l; i<N; i++)
 {denom=coeff1[i+segm*(N-1)][l+segm*(N-1)];//mem*=denom;
  for(int j=l; j<N; j++)
  {//printf("%f/%f\n",coeff[i+segm*(N-1)][j+segm*(N-1)],denom);
   coeff1[i+segm*(N-1)][j+segm*(N-1)]=coeff1[i+segm*(N-1)][j+segm*(N-1)]/denom;
  } 
 b_vect1[i+segm*(N-1)]=b_vect1[i+segm*(N-1)]/denom;
 }

 denom=coeff1[N-1+segm*(N-1)][l+segm*(N-1)];//mem*=denom;
 for(int j=l; j<2*N-1; j++)
 {//printf("%f/%f\n",coeff[N-1+segm*(N-1)][j+segm*(N-1)],denom);
  coeff1[N-1+segm*(N-1)][j+segm*(N-1)]=coeff1[N-1+segm*(N-1)][j+segm*(N-1)]/denom;
 } 
 b_vect1[N-1+segm*(N-1)]=b_vect1[N-1+segm*(N-1)]/denom;
}

  for(int i=l+1; i<N; i++)
  {
   for(int j=l; j<N; j++)
   {
    coeff1[i+segm*(N-1)][j+segm*(N-1)]=coeff1[i+segm*(N-1)][j+segm*(N-1)]-coeff1[l+segm*(N-1)][j+segm*(N-1)];
   }
   b_vect1[i+segm*(N-1)]=b_vect1[i+segm*(N-1)]-b_vect1[l+segm*(N-1)];
  }
 }
} 



for(int i=L-1; i>=0; i--)
{sum=0.0;
 for(int j=L-1; j>i; j--)
 {sum+=coeff1[i][j]*C3[j]; 
 }
 C3[i]=b_vect1[i]-sum; 
}

for(int i=0; i<L; i++){printf("C3[%d]=%f ",i,C3[i]);}
printf("\n");
free(b_vect1);
//-------------------------------------------------------------------------------------------------------------------     
//lu
float **U, **L_mat2, *C4, *sol2;
U=(float**)malloc(L*sizeof(float*));
for(int i=0; i<L; i++){U[i]=(float*)calloc(L,sizeof(float));}
L_mat2=(float**)malloc(L*sizeof(float*));
for(int i=0; i<L; i++){L_mat2[i]=(float*)calloc(L,sizeof(float));}

C4=(float*)calloc(L,sizeof(float));
sol2=(float*)calloc(L,sizeof(float));


for(int i=0; i<L; i++)
{for(int j=0; j<L; j++)
 {
  U[i][j]=coeff[i][j]; 
  if(i==j){L_mat2[j][j]=1.0;}
 }
}

for(int i=0; i<L-1; i++)
{
 for(int k=i+1; k<L; k++)
 {
  L_mat2[k][i]=U[k][i]/U[i][i];
 
  for(int j=i; j<L; j++)
  {
   U[k][j]=U[k][j]-L_mat2[k][i]*U[i][j];  
  }
 
 }
}



for(int i=0; i<L; i++)
{sum=0.0;
 for(int j=0; j<i; j++)
 {sum+=L_mat2[i][j]/L_mat2[i][i]*sol2[j]; 
 }
 sol2[i]=b_vect[i]/L_mat2[i][i]-sum;      
}

for(int i=L-1; i>=0; i--)
{sum=0.0;
 for(int j=L-1; j>i; j--)
 {sum+=U[i][j]/U[i][i]*C4[j]; 
 }
 C4[i]=sol2[i]/U[i][i]-sum; 
}

for(int i=0; i<L; i++){printf("C4[%d]=%f\t",i,C4[i]);}
printf("\n");

//==================================================================================================================
file=fopen("solution.txt","w");
fprintf(file,"     cg\t\t cholesky\t  seidel\t gauss\t\t    lu\t        library\n");
for(int i=0; i<L; i++)
{
 fprintf(file,"%f\t%f\t%f\t%f\t%f\t%f\n",C[i],C1[i],C2[i],C3[i],C4[i],C5[i]);
}

fclose(file);

FILE *res;
res=fopen("res.txt","w");
//float l_1(float *x,float *y, int n_min, int n_max)
float *lib_sol;
lib_sol=(float*)malloc(L*sizeof(float));
for(int i=0; i<L; i++)lib_sol[i]=C5(i);
fprintf(res,"         l1\t\tl2\t\tl_inf\n\n");
fprintf(res,"cg:      %f\t%f\t%f\ncholesky:%f\t%f\t%f\nseidel:  %f\t%f\t%f\ngauss:   %f\t%f\t%f\nlu:      %f\t%f\t%f\n",l_1(C,lib_sol,L,L),l_2(C,lib_sol,L,L),l_inf(C,lib_sol,L,L),l_1(C1,lib_sol,L,L),l_2(C1,lib_sol,L,L),l_inf(C1,lib_sol,L,L),l_1(C2,lib_sol,L,L),l_2(C2,lib_sol,L,L),l_inf(C2,lib_sol,L,L),l_1(C3,lib_sol,L,L),l_2(C3,lib_sol,L,L),l_inf(C3,lib_sol,L,L),l_1(C4,lib_sol,L,L),l_2(C4,lib_sol,L,L),l_inf(C4,lib_sol,L,L));
fprintf(res,"abs_err\n\n");
fprintf(res,"cg:      %f\t%f\t%f\ncholesky:%f\t%f\t%f\nseidel:  %f\t%f\t%f\ngauss:   %f\t%f\t%f\nlu:      %f\t%f\t%f\n",l_1_abs(C,lib_sol,L,L),l_2_abs(C,lib_sol,L,L),l_inf_abs(C,lib_sol,L,L),l_1_abs(C1,lib_sol,L,L),l_2_abs(C1,lib_sol,L,L),l_inf_abs(C1,lib_sol,L,L),l_1_abs(C2,lib_sol,L,L),l_2_abs(C2,lib_sol,L,L),l_inf_abs(C2,lib_sol,L,L),l_1_abs(C3,lib_sol,L,L),l_2_abs(C3,lib_sol,L,L),l_inf_abs(C3,lib_sol,L,L),l_1_abs(C4,lib_sol,L,L),l_2_abs(C4,lib_sol,L,L),l_inf_abs(C4,lib_sol,L,L));
fprintf(res,"res\n\n");
fprintf(res,"cg:      %f\t%f\t%f\ncholesky:%f\t%f\t%f\nseidel:  %f\t%f\t%f\ngauss:   %f\t%f\t%f\nlu:      %f\t%f\t%f\nlib:     %f\t%f\t%f\n",l_1(C,b_vect,L,L),l_2(C,b_vect,L,L),l_inf(C,b_vect,L,L),l_1(C1,b_vect,L,L),l_2(C1,b_vect,L,L),l_inf(C1,b_vect,L,L),l_1(C2,b_vect,L,L),l_2(C2,b_vect,L,L),l_inf(C2,b_vect,L,L),l_1(C3,b_vect,L,L),l_2(C3,b_vect,L,L),l_inf(C3,b_vect,L,L),l_1(C4,b_vect,L,L),l_2(C4,b_vect,L,L),l_inf(C4,b_vect,L,L),l_1(lib_sol,b_vect,L,L),l_2(lib_sol,b_vect,L,L),l_inf(lib_sol,b_vect,L,L));
fprintf(res,"abs_res\n\n");
fprintf(res,"cg:      %f\t%f\t%f\ncholesky:%f\t%f\t%f\nseidel:  %f\t%f\t%f\ngauss:   %f\t%f\t%f\nlu:      %f\t%f\t%f\nlib:     %f\t%f\t%f\n",l_1_abs(C,b_vect,L,L),l_2_abs(C,b_vect,L,L),l_inf_abs(C,b_vect,L,L),l_1_abs(C1,b_vect,L,L),l_2_abs(C1,b_vect,L,L),l_inf_abs(C1,b_vect,L,L),l_1_abs(C2,b_vect,L,L),l_2_abs(C2,b_vect,L,L),l_inf_abs(C2,b_vect,L,L),l_1_abs(C3,b_vect,L,L),l_2_abs(C3,b_vect,L,L),l_inf_abs(C3,b_vect,L,L),l_1_abs(C4,b_vect,L,L),l_2_abs(C4,b_vect,L,L),l_inf_abs(C4,b_vect,L,L),l_1_abs(lib_sol,b_vect,L,L),l_2_abs(lib_sol,b_vect,L,L),l_inf_abs(lib_sol,b_vect,L,L));


fclose(res);
free(lib_sol);
//-------------------------------------------------------------------------------------------------------------------

/*for(int i=0; i<L; i++)
{
 for(int j=0; j<L; j++)
 {
  printf("%f ",coeff[i][j]);
 }
 printf("\t%f\n",b_vect[i]);
} */
//for(int i=0; i<L; i++) printf("b[%d]=%f\n",i, b_vect[i]);
//for(int i=0; i<L; i++) printf("%f\n",C[i]);

for(int segm=0; segm<K; segm++) 
{for(int j=segm*P/K; j<=(segm+1)*P/K; j++)
 {if((j==segm*P/K)&&(j!=0)){j++;}
  for(int i=0; i<N; i++)
  {psi = 1.0/denom[i];
   for(int k=0; k<N; k++)
   {
    if(i!=k){psi*=(graphic_points[j]-x[k+segm*(N-1)]);}
   }
   interpol_values[j]+=(C[i+segm*(N-1)]*psi); 
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
   a1+=C[i+segm*(N-1)]*psi/denom[i];
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
   a1+=C[i+segm*(N-1)]*psi/denom[i];
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
free(C);
free(C1);
free(C2);
free(C3);
free(C4);
free(sol2);
free(b_vect);
for(int i=0; i<L; i++) {free(coeff[i]);}
for(int i=0; i<L; i++){free(coeff1[i]);}
free(coeff1);
for(int i=0; i<L; i++){free(U[i]);}
free(U);
for(int i=0; i<L; i++){free(L_mat2[i]);}
free(L_mat2);


free(coeff);
//system("Visualization.py");
return 0;
}