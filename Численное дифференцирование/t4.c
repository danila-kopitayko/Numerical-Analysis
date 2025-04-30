#include "stdio.h"
#include "stdlib.h"
#include "math.h"

double func(double arg)
{return (arg)*(arg)*(arg)*(arg);}


double der(double arg)
{return 4*(arg)*(arg)*(arg);}

int main(void)
{
double a=0.0;
double b=3.0;
int M=20; //количество узлов на отрезке 
int P=500; //количество точек, на которых строим графики
double h=fabs(b-a)/(M-1);
double h_viz = fabs(b-a)/P;
int M_new = 2*M-1;
double h_new = h/2;
double *x,*f_der,*derivative,*f_new,*graphic_points,*f,*f_der_new,*error,*second_Runge_formula;
x=(double*)malloc((M)*sizeof(double));
f_der=(double*)malloc((M)*sizeof(double));
f=(double*)malloc((M)*sizeof(double));
derivative=(double*)malloc((M)*sizeof(double));
f_der_new=(double*)malloc((M_new)*sizeof(double));
second_Runge_formula=(double*)malloc((M)*sizeof(double));
graphic_points = (double*)malloc((P+1)*sizeof(double));
error = (double*)malloc((M)*sizeof(double));
FILE *file_der,*file2;

//создаем массив точек, на которых строим производную
for(int i=0; i<M; i++)
{x[i]=a+i*h;
}


//создаем массив значений функций
for(int i=0; i<M; i++)
{f[i]=func(x[i]); 
}



//строим массив точек, на которых строим графики
for(int i=0; i<=P; i++)
{graphic_points[i] = a + i*h_viz;}

//задаем производную

for(int i=1; i<M-1; i++)
{f_der[i] = ((f[i+1]-f[i-1])/(2*h));
}
f_der[0] = (-3*f[0]+4*f[1]-f[2])/(2*h);
f_der[M-1] = (f[M-3]-4*f[M-2]+3*f[M-1])/(2*h);


//измельчаем сетку
f_new=(double*)malloc((M_new)*sizeof(double));
x=(double*)realloc(x,(M_new)*sizeof(double));

for(int i=0; i<M_new; i++){x[i]=a+i*h_new;}
for(int i=0; i<M_new; i++){f_new[i]=func(x[i]);}


for(int i=1; i<M_new-1; i++)
{f_der_new[i] = ((f_new[i+1]-f_new[i-1])/(2*h_new));
}
f_der_new[0] = (-3*f_new[0]+4*f_new[1]-f_new[2])/(2*h_new);
f_der_new[M_new-1] = (f_new[M_new-3]-4*f_new[M_new-2]+3*f_new[M_new-1])/(2*h_new);

//применяем правило Рунге
for(int i=0; i<M; i++)
{error[i]=(f_der[i]-f_der_new[2*i])/(0.5*0.5-1); 
}

//уточняем формулу для производной

for(int i=0; i<M; i++)
{second_Runge_formula[i] = f_der[i]+error[i];
}

double l1_abs=0.0,l2_abs=0.0,l_inf_abs=0.0;
double l1_rel=0.0,l2_rel=0.0,l_inf_rel=0.0;

//погрешности
for(int i=0; i<M; i++)
{l1_abs+=fabs(f_der[i]-der(x[2*i])); 
 l1_rel+=fabs(der(x[2*i]));
 l2_abs+=(f_der[i]-der(x[2*i]))*(f_der[i]-der(x[2*i]));
 l2_rel+=(der(x[2*i]))*(der(x[2*i]));
 l_inf_abs=fmax(l_inf_abs,fabs(der(x[2*i])-f_der[i]));
 l_inf_rel=fmax(l_inf_rel,fabs(der(x[2*i])));
}

file2=fopen("errors.txt","w");

fprintf(file2,"1st GRID\n");
fprintf(file2,"     \tABSOLUTE ERROR\t\tRELATIVE ERROR\n");
fprintf(file2,"l1   \t%e\t\t%e\n",l1_abs,l1_abs/l1_rel);
fprintf(file2,"l2   \t%e\t\t%e\n",sqrt(l2_abs),sqrt(l2_abs)/sqrt(l2_rel));
fprintf(file2,"l_inf\t%e\t\t%e\n",l_inf_abs,l_inf_abs/l_inf_rel);
fprintf(file2,"\n\n");

l1_abs=0.0,l2_abs=0.0,l_inf_abs=0.0;
l1_rel=0.0,l2_rel=0.0,l_inf_rel=0.0;

for(int i=0; i<M_new; i++)
{
 l1_abs+=fabs(f_der_new[i]-der(x[i])); 
 l1_rel+=fabs(der(x[i]));
 l2_abs+=(f_der_new[i]-der(x[i]))*(f_der_new[i]-der(x[i]));
 l2_rel+=(der(x[i]))*(der(x[i]));
 l_inf_abs=fmax(l_inf_abs,fabs(der(x[i])-f_der_new[i]));
 l_inf_rel=fmax(l_inf_rel,fabs(der(x[i])));
}


fprintf(file2,"2nd GRID\n");
fprintf(file2,"     \tABSOLUTE ERROR\t\tRELATIVE ERROR\n");
fprintf(file2,"l1   \t%e\t\t%e\n",l1_abs,l1_abs/l1_rel);
fprintf(file2,"l2   \t%e\t\t%e\n",sqrt(l2_abs),sqrt(l2_abs)/sqrt(l2_rel));
fprintf(file2,"l_inf\t%e\t\t%e\n",l_inf_abs,l_inf_abs/l_inf_rel);
fprintf(file2,"\n\n");

l1_abs=0.0,l2_abs=0.0,l_inf_abs=0.0;
l1_rel=0.0,l2_rel=0.0,l_inf_rel=0.0;


for(int i=0; i<M; i++)
{l1_abs+=fabs(second_Runge_formula[i]-der(x[2*i])); 
 l1_rel+=fabs(der(x[2*i]));
 l2_abs+=(second_Runge_formula[i]-der(x[2*i]))*(second_Runge_formula[i]-der(x[2*i]));
 l2_rel+=(der(x[2*i]))*(der(x[2*i]));
 l_inf_abs=fmax(l_inf_abs,fabs(der(x[2*i])-second_Runge_formula[i]));
 l_inf_rel=fmax(l_inf_rel,fabs(der(x[2*i])));
}


fprintf(file2,"Second Runge Formula\n");
fprintf(file2,"     \tABSOLUTE ERROR\t\tRELATIVE ERROR\n");
fprintf(file2,"l1   \t%e\t\t%e\n",l1_abs,l1_abs/l1_rel);
fprintf(file2,"l2   \t%e\t\t%e\n",sqrt(l2_abs),sqrt(l2_abs)/sqrt(l2_rel));
fprintf(file2,"l_inf\t%e\t\t%e\n",l_inf_abs,l_inf_abs/l_inf_rel);
fprintf(file2,"\n\n");


l1_abs=0.0,l2_abs=0.0,l_inf_abs=0.0,l1_rel=0.0,l2_rel=0.0,l_inf_rel=0.0;

for(int i=0; i<M; i++)
{l1_rel+=fabs(error[i]);
 l2_rel+=(error[i])*(error[i]);
 l_inf_rel=fmax(l_inf_rel,fabs(error[i]));
}

fprintf(file2,"First Runge Formula\n");
fprintf(file2,"     \tABSOLUTE ERROR\n");
fprintf(file2,"l1   \t%e\n",l1_rel);
fprintf(file2,"l2   \t%e\n",sqrt(l2_rel));
fprintf(file2,"l_inf\t%e\n",l_inf_rel);




file_der=fopen("points_and_values.txt","w");
for(int i=0; i<=P; i++){fprintf(file_der,"a\t%f\t%f\n",graphic_points[i],der(graphic_points[i]));}


for(int i=0; i<M; i++){fprintf(file_der,"b\t%f\t%f\n",x[2*i],f_der[i]);}

for(int i=0; i<M_new; i++){fprintf(file_der,"c\t%f\t%f\n",x[i],f_der_new[i]);}
               
for(int i=0; i<M; i++){fprintf(file_der,"d\t%f\t%f\n",x[2*i],second_Runge_formula[i]);}



fclose(file_der);
fclose(file2);
free(graphic_points);
free(f);
free(f_new);
free(f_der_new);
free(f_der);
free(derivative);
free(x);
free(second_Runge_formula);
free(error);
system("Visualization.py");
return 0;
}