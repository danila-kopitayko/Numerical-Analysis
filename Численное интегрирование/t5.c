#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double func(double x)
{return cos(x);}
         
double integral(double x1, double x2)
{return sin(x2)-sin(x1);
}

int main(void)
{
double a=0.0;
double b=3.0;
int K=10; //количество узлов
double h=fabs(b-a)/K; //шаг сетки
double h_new=fabs(b-a)/(2*K);
double mem;
double *x;
double integral_rect, integral_trap, integral_simpson,integral_gauss,integral_newton,error_rect;
double error_rect1,error_rect2,error_trap1,error_trap2,error_simpson1,error_simpson2,error_gauss1,error_gauss2,error_newton1,error_newton2;
x=(double*)malloc((K+1)*sizeof(double));
FILE *f;

for(int i=0; i<=K; i++)
{x[i] = a + i*h;
}

//составная формула прямоугольников
integral_rect=0.0;
for(int i=0; i<K; i++)
{mem=(x[i]+x[i+1])/2.0;
integral_rect += func(mem);
}
integral_rect*=h;


//составная формула трапеций 
integral_trap = func(a)*h/2.0+func(b)*h/2.0;
for(int i=1; i<K; i++) 
{integral_trap += func(x[i]);
}
integral_trap *= h;

//формула Симпсона
integral_simpson=0.0;
for(int i=0; i<K; i++)
{mem=(x[i]+x[i+1])/2.0; 
integral_simpson += (func(x[i]) + 4.0*func(mem) + func(x[i+1]));
}
integral_simpson*=(h/6.0);

//составная формула Гаусса
integral_gauss=0.0;
for(int i=0; i<K; i++)
{double arg1=0.0, arg2=0.0, arg3=0.0;
arg1=-((x[i+1]-x[i])/2.0)*sqrt(0.6)+(x[i]+x[i+1])/2.0;
arg2=(x[i]+x[i+1])/2.0;
arg3=((x[i+1]-x[i])/2.0)*sqrt(0.6)+(x[i+1]+x[i])/2.0;
integral_gauss+=(func(arg1)*5.0/9.0+func(arg2)*8.0/9.0+func(arg3)*5.0/9.0);
}
integral_gauss*=(h/2.0);

double vect[5]={7.0,32.0,12.0,32.0,7.0};
integral_newton=0.0;
double term;

for(int i=0; i<K; i++)
{term=0.0;
 for(int j=0; j<5; j++)
 {
  mem=x[i]+j*(x[i+1]-x[i])/4.0; //делим отрезок на четыре равные части и берем j-ую точку 
  term+=func(mem)*vect[j];
 }
integral_newton+=(h/90.0)*term;
} 


//считаем относительные погрешности
error_rect1 = fabs(integral_rect - integral(a,b))/fabs(integral(a,b));
error_trap1 = fabs(integral_trap - integral(a,b))/fabs(integral(a,b));
error_simpson1 = fabs(integral_simpson - integral(a,b))/fabs(integral(a,b));
error_gauss1 = fabs(integral_gauss - integral(a,b))/fabs(integral(a,b));
error_newton1 = fabs(integral_newton - integral(a,b))/fabs(integral(a,b));



//измельчаем сетку

x=(double*)realloc(x,(2*K+1)*sizeof(double));
for(int i=0; i<=2*K; i++){x[i]=a+h_new*i;}

integral_rect = 0.0;
integral_trap = 0.0;
integral_simpson = 0.0;
integral_gauss = 0.0;
integral_newton = 0.0;


//составная формула прямоугольников
for(int i=0; i<2*K; i++)  //2K-1
{mem=(x[i]+x[i+1])/2.0;
integral_rect += func(mem);
}
integral_rect*=h_new;

//составная формула трапеций
integral_trap = func(a)*h_new/2.0 + func(b)*h_new/2.0;
for(int i=1; i<2*K; i++)
{
integral_trap += func(x[i]);
}
integral_trap *= h_new;

//формула Симпсона
integral_simpson=0.0;
for(int i=0; i<2*K; i++)
{mem=(x[i]+x[i+1])/2.0; 
integral_simpson += (func(x[i]) + 4.0*func(mem) + func(x[i+1]));
}                   
integral_simpson*=h_new/6.0;


for(int i=0; i<2*K; i++)
{double arg1=0.0, arg2=0.0, arg3=0.0;
arg1=((x[i+1]-x[i])/2.0)*(-sqrt(0.6))+(x[i]+x[i+1])/2.0;
arg2=(x[i]+x[i+1])/2.0;
arg3=(x[i+1]-x[i])/2.0*sqrt(0.6)+(x[i+1]+x[i])/2.0;
integral_gauss+=(func(arg1)*(5.0/9.0)+func(arg2)*(8.0/9.0)+func(arg3)*(5.0/9.0));
}
integral_gauss*=h_new/2.0;



for(int i=0; i<2*K; i++)
{term=0.0;
 for(int j=0; j<5; j++)
 {
  mem=x[i]+j*((x[i+1]-x[i])/4.0); //делим отрезок на четыре равные части и берем j-ую точку 
  term+=func(mem)*vect[j];
 }
integral_newton+=term;
} 
integral_newton*=(h_new/90.0);


error_rect2 = fabs(integral_rect - integral(a,b))/fabs(integral(a,b));
error_trap2 = fabs(integral_trap - integral(a,b))/fabs(integral(a,b));
error_simpson2 = fabs(integral_simpson - integral(a,b))/fabs(integral(a,b));
error_gauss2 = fabs(integral_gauss - integral(a,b))/fabs(integral(a,b));
error_newton2 = fabs(integral_newton - integral(a,b))/fabs(integral(a,b));


f=fopen("Errors.txt","w");

fprintf(f,"         \tH STEP GRID\t1/2H STEP GRID\n");
fprintf(f,"Rectangle\t%e\t%e\n",error_rect1,error_rect2);
fprintf(f,"Trapezoid\t%e\t%e\n",error_trap1,error_trap2);
fprintf(f,"Simpson  \t%e\t%e\n",error_simpson1,error_simpson2);
fprintf(f,"Gauss    \t%e\t%e\n",error_gauss1,error_gauss2);
fprintf(f,"Newton   \t%e\t%e\n",error_newton1,error_newton2);
fclose(f);


free(x);
return 0;
}