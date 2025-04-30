#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>

double *read_file(const char *filename, int n)
{
    FILE *fptr;
    char str[50], *s, *s1, *num;
    int i = 0, p = 0, length, flag = 0;

    double *arr;
    arr = (double *)calloc(n, sizeof(double));

    fptr = fopen(filename, "r");

    while (fgets(str, sizeof(str), fptr) && (p < n))
    { // printf("%s\n",str);
        s = str;
        i = 0;
        flag = 0;
        while ((*s) || (i < 50))
        {
            if (*s == '=')
            {
                s1 = s;
                flag++;
            } // printf("%c",*s1);}
            if (!isdigit(*s) && isdigit(*(s - 1)) && (!isdigit(*(s + 1))) && (*s != '-') && (*s != 'e'))
            {
                length = s + 1 - s1;
                num = (char *)malloc(length * sizeof(char));
                strncpy(num, s1 + 1, length - 1);
                arr[p] = atof(s1 + 1);
                // printf("function1 %d %f",p,arr[p]);
                free(num);
            }

            i++;
            s++;
        }
        if (flag++)
        {
            p++;
        }

        // printf("\n");
    }

    // for(int i=0; i<n; i++)printf("function2 %f\n",arr[i]);
    fclose(fptr);
    return arr;
}

double N1(double x, double l)
{
    return 1.0 - 3.0 * pow(x / l, 2) + 2.0 * pow(x / l, 3);
}

double N2(double x, double l)
{
    return 3.0 * pow(x / l, 2) - 2.0 * pow(x / l, 3);
}

double L1(double x, double l)
{
    return x - 2.0 * x * (x / l) + x * pow(x / l, 2);
}

double L2(double x, double l)
{
    return -x * (x / l) + x * pow(x / l, 2);
}

double N1_der1(double x, double l)
{
    return -6.0 * x / pow(l, 2) + 6.0 * pow(x, 2) / pow(l, 3);
}

double N2_der1(double x, double l)
{
    return 6.0 *  x / pow(l, 2) - 6.0 * pow(x, 2) / pow(l, 3);
}

double L1_der1(double x, double l)
{
    return 1 - 4.0 * (x / l) + 3.0 * pow(x / l, 2);
}

double L2_der1(double x, double l)
{
    return -2.0 * x / l + 3.0 * pow(x / l, 2);//-x * (x / l) + x * pow(x / l, 2);
}

double N1_der2(double x, double l)
{
    return -6.0 / pow(l, 2) + 12.0 * x / pow(l, 3);
}

double N2_der2(double x, double l)
{
    return 6.0 / pow(l, 2) - 12.0 * x / pow(l, 3);
}

double L1_der2(double x, double l)
{
    return - 4.0 / l + 6.0 * x / pow(l, 2);
}

double L2_der2(double x, double l)
{
    return -2.0 / l + 6.0 * x / pow(l, 2);
}

double N1_der3(double l)
{
    return 12.0 / pow(l, 3);
}

double N2_der3(double l)
{
    return - 12.0 / pow(l, 3);
}

double L1_der3(double l)
{
    return 6.0 / pow(l, 2);
}

double L2_der3(double l)
{
    return 6.0 / pow(l,2);
}



void print_matrix(double **mat, int m, int n)
{
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            printf("%f ", mat[i][j]);
        }
        printf("\n");
    }
}

double **create_matrix(int m, int n)
{
    double **A;
    A = (double **)calloc(m, sizeof(double *));

    for (int i = 0; i < m; i++)
        A[i] = (double *)calloc(n, sizeof(double));

    return A;
}

double **stiffness_matrix(double **A, double E, double J, double L)
{
    A[0][0] = 12 * E * J / pow(L, 3);
    A[0][1] = 6 * E * J / pow(L, 2);
    A[0][2] = -12 * E * J / pow(L, 3);
    A[0][3] = 6 * E * J / pow(L, 2);
    A[1][0] = 6 * E * J / pow(L, 2);
    A[1][1] = 4 * E * J / L;
    A[1][2] = -6 * E * J / pow(L, 2);
    A[1][3] = 2 * E * J / L;
    A[2][0] = -12 * E * J / pow(L, 3);
    A[2][1] = -6 * E * J / pow(L, 2);
    A[2][2] = 12 * E * J / pow(L, 3);
    A[2][3] = -6 * E * J / pow(L, 2);
    A[3][0] = 6 * E * J / pow(L, 2);
    A[3][1] = 2 * E * J / L;
    A[3][2] = -6 * E * J / pow(L, 2);
    A[3][3] = 4 * E * J / L;

    return A;
}

void clean(double **mat, int m)
{
    for (int i = 0; i < m; i++)
    {
        free(mat[i]);
    }

    free(mat);
}

/*double *force_vect(double P, double l)
{
    double *N;
    N = (double *)malloc(4 * sizeof(double));
    N[0] = 0.5 * l * P;
    N[1] = pow(l, 2) * P / 12;
    N[2] = 0.5 * P * l;
    N[3] = -pow(l, 2) * P / 12;
    return N;
}*/

double *load_vect(double x0, double x1, double q0, double q1, float a, float b, double l)
{
    double *N;
    N = (double *)malloc(4 * sizeof(double));

    if (((a < b) && (b <= x0)) || ((a < b) && (a >= x1)))
    {
        N[0] = 0.0;
        N[1] = 0.0;
        N[2] = 0.0;
        N[3] = 0.0;
    }

    else
    {
        N[0] = l * (7 * q0 + 3 * q1) / 20;
        N[1] = pow(l, 2) * (q0 / 20 + q1 / 30);
        N[2] = l * (3 * q0 + 7 * q1) / 20;
        N[3] = - pow(l, 2) * (q0 / 30 + q1 / 20);
    }
    return N;
}



double *cg_modified(double **A, double *R, int n, float eps, int *index, double *values, int size)
{
    double *G, *x, *P, *Y;
    double alpha, beta, mem1 = 0.0, mem2 = 0.0,res;
    int cnt = 0;
    G = (double *)calloc(n, sizeof(double));
    x = (double *)calloc(n, sizeof(double));
    P = (double *)malloc(n * sizeof(double));
    Y = (double *)malloc(n * sizeof(double));

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            G[i] += A[i][j] * x[j];
        }
        G[i] = G[i] - R[i];
        P[i] = - G[i];
    }

    for (int i = 0; i < size; i++)
    {
        if (index[i] < n)
        {
            P[index[i]] = values[i];
        }
        else
            printf("index[%d]=%d out of range\n", i, index[i]);
    }

    do
    {
        cnt++;

        for (int i = 0; i < n; i++)
        {
            Y[i] = 0.0;
            for (int j = 0; j < n; j++)
            {
                Y[i] += A[i][j] * P[j];
            }
        }

        for (int i = 0; i < size; i++)
        {
            if (index[i] < n)
            {
                Y[index[i]] = values[i];
            }
            else
                printf("index[%d]=%d out of range\n", i, index[i]);
        }

        //for(int i=0; i < n; i++) {printf("Y[%d]=%f\n",i,Y[i]);}

        // printf("cnt=%d  ",cnt);for(int i=0; i<n; i++){printf("cg: z[%d]=%f\t",i,z[i]);}printf("\n");

        mem1 = 0.0;
        mem2 = 0.0;
        beta = 0.0;
        for (int i = 0; i < n; i++)
        {
            mem1 += - P[i] * G[i];
            mem2 += P[i] * Y[i];
        }

        alpha = mem1 / mem2;

        // printf("cnt = %d  mem1=%f   mem2=%f alpha=%f\n",cnt,mem1,mem2,alpha);

        for (int i = 0; i < n; i++)
        {
            x[i] = x[i] + alpha * P[i]; // printf("r before:%f\n",r[i]);
            G[i] = G[i] + alpha * Y[i]; // printf("r after:%f\n",r[i]);
            beta += G[i] * Y[i];
        }

        beta = beta / mem2; // printf("beta=%f mem=%f\n",beta,mem1);

        for (int i = 0; i < n; i++)
        {
            P[i] = - G[i] + beta * P[i];
        }

        for (int i = 0; i < size; i++)
        {
            if (index[i] < n)
            {
                P[index[i]] = values[i];
            }
            else
                printf("index[%d]=%d out of range\n", i, index[i]);
        }

        
        //printf("beta=%f\n", beta);

        if (cnt > 100)
        {
            printf("error\n");
            break;
        }
    res = 0.0;
    for (int i = 0; i< n; i++) {res += G[i] * G[i];}

    printf("%d %f\n", cnt, sqrt(beta*mem2));

    } while (sqrt(beta*mem2) > eps);

    free(P);
    free(Y);
    return x;
}




int find_index(double *arr, int size, double value)
{
    //printf("CALL FROM FIND_INDEX: SIZE = %d\tVALUE=%f\n",size, value);
    for(int i = 0; i < size; i++)
    {   //printf("fabs[i]=%f\n",fabs(arr[i]-value));
        if(arr[i] == value)
        {
             return i;
        }
        
            
	
        
    }
    printf("index not found\n");
    return -1;
  
}

//==================================================================================================================================

int main(void)
{
    int N = 2, P, M, spring, R_a, R_c, R_e, q_b, q_f;
    double **A, *delta, *R, *arr, **A_full, *R_full, **H, *x, *lengths, *W, *Theta, *Moment, *Q;
    double x1, x2, memry;
    arr = (double *)malloc(15 * sizeof(double));
    arr = read_file("values.txt", 15);
    FILE *output;
    int points = 100;
    output = fopen("output_c.txt", "w");
    double x0 = 0.0, a = arr[0], b = arr[1], c = arr[2], d = arr[3], e = arr[4], f = arr[5], g = arr[6], h = arr[7], E = arr[8], J = arr[9], k = arr[13];
    R = (double *)calloc(2 * N, sizeof(double));
    x = (double *)calloc(points, sizeof(double));
    A = create_matrix(2 * N, 2 * N);

    double elements[] = {0.0,1.0,2.0,3.0,a,5.0,6.0,7.0, b,9.0,10.0, c, 12.0, 13.0, 14.0, d, 16.0,17.0,18.0, e, f, g,22.0, h,24.0,25.0};



    int total_nodes, total_elements;// = sizeof(elements) / sizeof(double);
    

    total_nodes = sizeof(elements) / sizeof(double);
    total_elements = total_nodes - 1;
    R_full = (double *)calloc(2 * total_elements + 2, sizeof(double));

    A_full = create_matrix(2 * total_elements + 2, 2 * total_elements + 2);
    delta = (double *)calloc(2 * total_elements + 2, sizeof(double));

    P = find_index(elements, total_elements, d);
    M = find_index(elements, total_elements, g);
    spring = find_index(elements, total_elements, h);
    R_a = find_index(elements, total_elements, a);
    R_c = find_index(elements, total_elements, c);
    R_e = find_index(elements, total_elements, e);
    q_b = find_index(elements, total_elements, b);
    q_f = find_index(elements, total_elements, f);

    //printf("total_nodes = %d  total_elements=%d\n",total_nodes, total_elements);

    for(int i = 0; i < total_nodes; i++)
    {  
        for(int j = 0; j < total_nodes - 1; j++)
        {
            if (elements[j] > elements[j + 1])
            {
                memry = elements[j];
                elements[j] = elements[j + 1];
                elements[j + 1] = memry;               	
            }           
        }         
    } 


   
    lengths = (double *)malloc((total_elements)*sizeof(double));

    for(int i = 0; i < total_elements; i++)
    {
        lengths[i] = elements[i + 1] - elements[i];
        //printf("lengths[%d]=%f\n",i,lengths[i]);
     
    }

    

    double q1 = -arr[10], q11;
    double q2 = -arr[11], q12;

    //-------------------------------------------------------------------------------------------------------------------------------------

    x1 = x0;
    x2 = x0;
    q11 = q1;
    q12 = q1;

    double coeff = (q2 - q1) / (f - b);

    for (int elmnt = 0; elmnt < total_elements; elmnt++)
    {   
        q12 += lengths[elmnt] * coeff;
        A = stiffness_matrix(A, E, J, lengths[elmnt]);
        x2 += lengths[elmnt];
        R = load_vect(b, f, q1, q2, x1, x2, lengths[elmnt]);
        //printf("x1=%f  x2=%f\n",x1,x2);
        x1 += lengths[elmnt];
        q11 += lengths[elmnt] * coeff;
        for (int i = 0; i < 2 * N; i++)
        {
            for (int j = 0; j < 2 * N; j++)
            {
                A_full[i + 2 * elmnt][j + 2 * elmnt] += A[i][j];
            }
            R_full[i + 2 * elmnt] += R[i];
        }
    }


    R_full[2 * M + 1] += -arr[14];
    R_full[2 * P] += -arr[12];
    A_full[2 * spring][2 * spring] += arr[13];

/*    for(int i=0; i<2 * total_elements + 2; i++)
    {printf("THIS IS R_FULL: R[%d]=%f\n",i,R_full[i]);}
*/   

    /*for (int i = 0; i < 2 * total_elements + 2; i++)
        printf("R[%d]=%f\n", i, R_full[i]);
    */


    double values[]={0.0,0.0,0.0,0.0,0.0};
    int index[]={0,1,2 * R_a,2 * R_c,2 * R_e};

    delta = cg_modified(A_full,R_full,2 * total_elements + 2, 1e-3, index, values, 5);                           

/*    for (int i = 0; i < 2 * total_elements + 2; i++)
    {
        printf("delta[%d]=%f\n", i, delta[i]);
    }
*/

    W = (double *)malloc(points * total_elements * sizeof(double));
    Theta = (double *)malloc(points * total_elements * sizeof(double));
    Moment = (double *)malloc(points * total_elements * sizeof(double));
    Q = (double *)malloc(points * total_elements * sizeof(double));



    x1 = x0;
    for (int elmnt = 0; elmnt < total_elements; elmnt++)
    {

        for(int i = 0; i < points; i++)
        {
            x[i] = i * lengths[elmnt] / points;
        }


        for (int i = 0; i < points; i++)
        {
            W[i] = delta[2 * elmnt] * N1(x[i], lengths[elmnt]) + delta[1 + 2 * elmnt] * L1(x[i], lengths[elmnt]) + delta[2 + 2 * elmnt] * N2(x[i], lengths[elmnt]) + delta[3 + 2 * elmnt] * L2(x[i], lengths[elmnt]);
	    fprintf(output, "W\t%1.12e  %1.12e\n", x[i] + elements[elmnt], W[i]);
        }
    }


    x1 = x0;
    for (int elmnt = 0; elmnt < total_elements; elmnt++)
    {

        for(int i = 0; i < points; i++)
        {
            x[i] = i * lengths[elmnt] / points;
        }


        for (int i = 0; i < points; i++)
        {
            Q[i] = E * J * (delta[2 * elmnt] * N1_der3(lengths[elmnt]) + delta[1 + 2 * elmnt] * L1_der3(lengths[elmnt]) + delta[2 + 2 * elmnt] * N2_der3(lengths[elmnt]) + delta[3 + 2 * elmnt] * L2_der3(lengths[elmnt]));	    
	    fprintf(output, "Q\t%1.12e  %1.12e\n", x[i] + elements[elmnt], Q[i]);
        }
    }


    x1 = x0;
    for (int elmnt = 0; elmnt < total_elements; elmnt++)
    {

        for(int i = 0; i < points; i++)
        {
            x[i] = i * lengths[elmnt] / points;
        }


        for (int i = 0; i < points; i++)
        {
            Moment[i] = E * J * (delta[2 * elmnt] * N1_der2(x[i], lengths[elmnt]) + delta[1 + 2 * elmnt] * L1_der2(x[i], lengths[elmnt]) + delta[2 + 2 * elmnt] * N2_der2(x[i], lengths[elmnt]) + delta[3 + 2 * elmnt] * L2_der2(x[i], lengths[elmnt]));
	    fprintf(output, "Moment\t%1.12e  %1.12e\n", x[i] + elmnt * lengths[elmnt], Moment[i]);
        }
    }



    x1 = x0;
    for (int elmnt = 0; elmnt < total_elements; elmnt++)
    {
        for(int i = 0; i < points; i++)
        {
            x[i] = i * lengths[elmnt] / points;
        }

        for (int i = 0; i < points; i++)
        {
            Theta[i] = delta[2 * elmnt] * N1_der1(x[i], lengths[elmnt]) + delta[1 + 2 * elmnt] * L1_der1(x[i], lengths[elmnt]) + delta[2 + 2 * elmnt] * N2_der1(x[i], lengths[elmnt]) + delta[3 + 2 * elmnt] * L2_der1(x[i], lengths[elmnt]);	    
	    fprintf(output, "Theta\t%1.12e  %1.12e\n", x[i] + elmnt * lengths[elmnt], Theta[i]);
        }
    }







    fclose(output);
    free(W);
    free(Moment);
    free(Q);
    free(Theta);
    free(arr);
    free(delta);
    free(R);
    free(R_full);
    free(lengths);
    clean(A, N);
    //clean(H, 2 * N);
    clean(A_full, 2 * total_elements + 2);
    free(x);
    printf("THE END\n");
    return 0;
}