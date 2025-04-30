#include <stdio.h>
#include <stdlib.h>
#include <math.h>

float f0(float *x)
{
 return exp(x[1]);
}

float f1(float *x)
{      
 return -sqrt(x[0]);
}

//===============================================================================================================

float *fixed_point_iteration(float *x, float eps,FILE *dvg)
{FILE *file;
 float norm=1.0;
 int cnt=0;

file=fopen("fpi.txt","w");
fprintf(file, "x0;x1;norm;cnt\n");
while(norm>eps)
{cnt++;
 x[0]=f0(x);
 x[1]=f1(x);

 norm=pow((f0(x)-x[0]),2)+pow((f1(x)-x[1]),2);//fmax(fabs(f0(x)-x[0]),fabs(f1(x)-x[1]));
 
 
 fprintf(file, "%f;%f;%f;%d\n",x[0],x[1],norm,cnt);
 

 if(cnt==100){printf("cnt>100\n");break;}
}
//dvg=fopen("divergence.txt","w");

fprintf(dvg,"fixed_point_iteration:\t        %d\n",cnt);
fclose(file);
//fclose(dvg);
return x;
}

//===============================================================================================================

float *newton(float *x, float eps,FILE *dvg)
{
 FILE *file;
 float det, mem2, norm, mem;
 int cnt=0, n=2;
 float **J, **J_new, *F;
 F=(float*)malloc(n*sizeof(float));

 J=(float**)malloc(n*sizeof(float));

 for(int i=0; i<n; i++)
 {
  J[i]=(float*)malloc(n*sizeof(float));
 }
 J_new=(float**)malloc(n*sizeof(float));
 
 for(int i=0; i<n; i++)
 {
  J_new[i]=(float*)malloc(n*sizeof(float));
 }

 file=fopen("newton.txt","w");
 fprintf(file,"x0;x1;norm;cnt\n");

 do
 {
  J[0][0]=1;
  J[1][1]=-1;
  J[0][1]=-2*x[1];
  J[1][0]=1/x[0];

  det=J[0][0]*J[1][1]-J[0][1]*J[1][0];
  mem2=J[0][0];
  J[0][0]=J[1][1]/det;
  J[1][1]=mem2/det;

  J[1][0]=-J[1][0]/det;
  J[0][1]=-J[0][1]/det;


  F[0]=x[0]-x[1]*x[1];
  F[1]=log(x[0])-x[1];


  norm=0.0;
  for(int i=0; i<n; i++)
  {mem=0.0;
   for(int j=0; j<n; j++)
   {
    mem+=J[i][j]*F[j]; 
   }
   x[i]=x[i]-mem;
   norm+=mem*mem;
  }
  cnt++;
  fprintf(file, "%f;%f;%f;%d\n",x[0],x[1],norm,cnt);
  if(cnt>100){printf("cnt>100\n");break;}
 }while(eps<sqrt(norm)); 

 for(int i=0; i<n; i++)
 {free(J[i]);
 }
 free(J);

 for(int i=0; i<n; i++)
 {free(J_new[i]);
 }
 free(J_new);
 free(F);
// dvg=fopen("divergence.txt","w");
 //fprintf(dvg,"method                \t    divergence\n");
 fprintf(dvg,"newton:               \t        %d\n",cnt);
// fclose(dvg);

 fclose(file);
 return x;
}

//===============================================================================================================

float *newton_modified(float *x, float eps, int step,FILE *dvg)
{FILE *file;
 int n=2, cnt=0;
 float norm=0.0,det,mem;
 float **J,**J_new,*F;

 F=(float*)malloc(n*sizeof(float));

 J=(float**)malloc(n*sizeof(float));

 for(int i=0; i<n; i++)
 {
  J[i]=(float*)malloc(n*sizeof(float));
 }
 J_new=(float**)malloc(n*sizeof(float));

 for(int i=0; i<n; i++)
 {
  J_new[i]=(float*)malloc(n*sizeof(float));
 }

 float mem2=0.0;
 int period=0;

 file=fopen("newton_mod.txt","w");
 fprintf(file,"x0;x1;norm;cnt\n");

 do
 {period+=1;

   if((period==step)||(cnt==0))
   {period=0;

    J[0][0]=1;
    J[1][1]=-1;
    J[0][1]=-2*x[1];   
    J[1][0]=1/x[0];


    det=J[0][0]*J[1][1]-J[0][1]*J[1][0];
    mem2=J[0][0];
    J[0][0]=J[1][1]/det;
    J[1][1]=mem2/det;
    J[1][0]=-J[1][0]/det;
    J[0][1]=-J[0][1]/det;

   }

  F[0]=x[0]-x[1]*x[1];
  F[1]=log(x[0])-x[1];


  norm=0.0;
  for(int i=0; i<n; i++)
  {mem=0.0;
   for(int j=0; j<n; j++)
   {
    mem+=J[i][j]*F[j];
   }
   x[i]=x[i]-mem;
   norm+=mem*mem;
  }

  cnt++;
  fprintf(file, "%f;%f;%f;%d\n",x[0],x[1],norm,cnt);
  if(cnt>100){printf("cnt>100\n");break;}
 }while(eps<sqrt(norm)); 


 for(int i=0; i<n; i++)
 {free(J[i]);
 }
 free(J);

 for(int i=0; i<n; i++)
 {free(J_new[i]);
 }
 free(J_new);

// dvg=fopen("divergence.txt","w");
 //fprintf(dvg,"method                \t    divergence\n");
 fprintf(dvg,"newton_mod:             \t%d\n",cnt);

// fclose(dvg);


 free(F); 
 fclose(file);
 return x;
}

//===============================================================================================================

float *newton_diff(float *x1,float *x2, float eps,FILE *dvg)
{
 FILE *file;
 int n=2, cnt=0;
 float norm=0.0,det,mem;
 float  **J,*F;

 F=(float*)malloc(n*sizeof(float));
 J=(float**)malloc(n*sizeof(float));


 for(int i=0; i<n; i++)
 {
  J[i]=(float*)malloc(n*sizeof(float));
 }

 float mem2=0.0;
 float h_x=0.0, h_y=0.0;
 
 file=fopen("newton_diff.txt","w");
 fprintf(file,"x0;x1;norm;cnt\n");

 do
 {h_x=fabs(x2[0]-x1[0]);  
  h_y=fabs(x2[1]-x1[1]);

  F[0]=x2[0]-x2[1]*x2[1];
  F[1]=log(x2[0])-x2[1]; 

  J[0][0]=(x2[0]-x1[0])/h_x;
  J[1][1]=-(x2[1]-x1[1])/h_y;
  J[0][1]=-(pow(x2[1],2)-pow(x1[1],2))/h_y;
  J[1][0]=(log(x2[0])-log(x1[0]))/h_x;

  //==========================

  det=J[0][0]*J[1][1]-J[0][1]*J[1][0];
  mem2=J[0][0];
  J[0][0]=J[1][1]/det;
  J[1][1]=mem2/det;
  J[1][0]=-J[1][0]/det;
  J[0][1]=-J[0][1]/det;

  norm=0.0;
  for(int i=0; i<n; i++)
  {mem=0.0;
   for(int j=0; j<n; j++)
   {
    mem+=J[i][j]*F[j];// printf("%f*%f=%f\n",J[i][j],F[j],mem); 
   }
   x2[i]=x2[i]-mem;
   norm+=mem*mem;
//   printf("mem[%d]=%f\n",i,mem);
  }

  cnt++;
  fprintf(file,"%f;%f;%f;%d\n",x2[0],x2[1],norm,cnt);
//  printf("cnt=%d norm=%f\n",cnt,sqrt(norm));
//  printf("x=%f y=%f\n",x2[0],x2[1]);
  if(cnt>100){printf("cnt>100\n");break;}
 }while(eps<sqrt(norm)); 


 for(int i=0; i<n; i++)
 {free(J[i]);
 }
 free(J);

 free(F);
// dvg=fopen("divergence.txt","w");
// fprintf(dvg,"method                \t    divergence\n");
 
 fprintf(dvg,"newton_diff:\t                %d\n",cnt);

// fclose(dvg);

 fclose(file);
 //printf("%f %f\n",x2[0],x2[1]);
 return x2;
}

int main(void)
{
int n=2,cnt=0, P=500;
float norm=1.0, eps=pow(10,-6),tau=0.0001,x1=0.0, x2=2.0, y1=-2.0, y2=1.0,h_x=fabs(x2-x1)/P,h_y=fabs(y2-y1)/P;
float *x,*f;
x=(float *)malloc(n*sizeof(float));
f=(float*)malloc(n*sizeof(float));
FILE *file, *dvg;
dvg=fopen("divergence.txt","w");
fprintf(dvg,"method                \t    divergence\n");

x[0]=0.4;
x[1]=-0.6;

printf("FIXED ITERATION\n");
x=fixed_point_iteration(x,eps,dvg);
for(int i=0; i<2; i++)
 {printf("x[%d]=%f\n",i,x[i]);}

x[0]=0.4;
x[1]=-0.6;

printf("NEWTON\n");
x=newton(x,eps,dvg);
for(int i=0; i<2; i++)
{printf("x[%d]=%f\n",i,x[i]);}

x[0]=0.4;
x[1]=-0.6;

printf("NEWTON_MODIFIED\n");
x=newton_modified(x,eps,4,dvg);
for(int i=0; i<2; i++)
{printf("x[%d]=%f\n",i,x[i]);}

float *x22;
x22=(float*)malloc(n*sizeof(float));
x[0]=0.4;
x[1]=-1.5;
x22[0]=0.45;
x22[1]=-1.4;


printf("NEWTON_DIFF\n");
x=newton_diff(x,x22,eps,dvg);
printf("%f %f\n",x[0],x[1]);
//for(int i=0; i<2; i++)
//{printf("x[%d]=%f\n",i,x[i]);}


//записываем данные в файл
file=fopen("t.txt","w");
float arr[2]={x1,y1};
for(int i=0; i<P; i++)
{
 fprintf(file,"%f\t%f\t%f\t%f\n",arr[0],arr[1],f0(arr),f1(arr));
 arr[0]+=h_x;
 arr[1]+=h_y;
}


fclose(file);
fclose(dvg);
free(x);
free(x22);
free(f);
system("graphics.py");
printf("THIS IS THE END\n");
return 0;
}