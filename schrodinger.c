/************************************************************************************************
 *												*
 *   COMPILED WITH THE LINE:									*
 *   gcc -Wall -Wextra -o schrodinger.x schrodinger.c -lm -lgsl -lgslcblas && ./schrodinger.x	*
 *												*
 *												*
 ************************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>

#define pi 3.141592653589793
//#define hbar (1.054571800E-34)
#define mass 1 // (6*12+6) // m(benzene)
#define coeff (0.5/(mass)) //  hbar^2/2m in atomic units
#define BasisOrder 3 // for a quadratic polynomial x^0, x^1, x^2 => 3 coefficients 
// !!! only use Lagrange polynomials of even power (odd number of basis function)
#define IntOrder 9 // number of terms in the Gaussian integration

#define MSG(msg) printf( "\n" #msg "\n")
#define DIV printf("\n=======================================================================================\n\n")
#define DEBUG(msg, var, fmt) printf( #msg "\n" #var "=%" #fmt "\n", var)
#define PRINT(var) printf( #var "=%f\n", var)
 #define PRINTF(var, fmt) printf("\n**DEBUG: " #var "=%" #fmt "\n", var)

const double w[9]={0.330239355001260, 0.312347077040003, 0.260610696402935, 0.180648160694857, 0.081274388361574, 0.312347077040003, 0.260610696402935, 0.180648160694857, 0.081274388361574}; // Gaussian weight factors
const double x[9]={0.000000000000, 0.324253423403809, 0.613371432700590, 0.836031107326636, 0.968160239507626, -0.324253423403809, -0.613371432700590, -0.836031107326636, -0.968160239507626}; // x value for each Gaussian weight factor

int xmax; 
int xmin;

// Parameters for the potential
int pottype=0;    // which potential function to use
double cosfactor=1;  // n in cos(n*x)
double boxlength=1;  // L
double boxv1=1;      // height for x < -L/2
double boxv2=1;      // height for x > L/2
double K=1;	// force constant of the harmonic oscillator, K*x^2/2

double XMIN=0;
double XMAX=2*pi;

// Subroutines (defined at the end of the code)
int alloc2d(double **, int, int);
void dealloc(double **,int);
void printArr(double **, int, int);
void printMat(gsl_matrix *, int, int);
void LagrangeCoeff(double **);
double Vbox(double, double, double);
double Vharm(double, double, double);
double Vcos(double, double, double);

/*******************************************************************/

int main()
{

    inline double basis(double **chi, int i, double x)
    {
	double result = chi[i][0] + chi[i][1]*x + chi[i][2]*x*x;
	return result;
    }
    inline double dbasis(double **chi, int i, double x)
    {
	double result = chi[i][1] + 2*chi[i][2]*x;
	return result;
    }

    double (*Vx)(double,double,double); // function pointer for the potential of choice

    FILE *potfile;
    FILE *cfile;

    int Nelem=1; 
    printf("\nEnter number of elements: "); scanf("%d",&Nelem); // number of elements the whole interval (0..2pi) is separated into

    int MatrixOrder=(((BasisOrder-1)*Nelem)+1);
    double **chi = malloc(BasisOrder*sizeof(double *));
//    double **dchi = malloc((BasisOrder)*sizeof(double *));
    double *centre = malloc((Nelem)*sizeof(double));
    double *length = malloc((Nelem)*sizeof(double));
    double *step = malloc((Nelem+1)*sizeof(double));

    double temp = 0; // used for the summation in some loops 
    double debugging = 0; 
    double prevS=0; double prevH=0; // ensure the overlap between two consecutive integration intervals
    int i,j,k,n; // loop indices
    int debug=0;
    int line=0;
    char gnuplot[100] = "gnuplot"; // the string to execute a system call to gnuplot


    gsl_matrix *T, *V, *s; // local matrices for a single element
    gsl_matrix *H, *S, *C; // global matrices for the whole interval
    gsl_matrix *Hb, *Sb; // global matrices after the boundary conditions are applied
    gsl_vector *E;
    gsl_eigen_gensymmv_workspace *wrkEig; // for generalized eigenvalue problem

    T = gsl_matrix_calloc(BasisOrder,BasisOrder);
    V = gsl_matrix_calloc(BasisOrder,BasisOrder);
    s = gsl_matrix_calloc(BasisOrder,BasisOrder);
    H = gsl_matrix_calloc(MatrixOrder,MatrixOrder);
    S = gsl_matrix_calloc(MatrixOrder,MatrixOrder);

    // *wrkEig, *Hb, *Sb, *C, *E are allocated later after the potential type is chosen because the dimension depends on the boundary conditions
     
    /*****************************************************************************************/

    // Gaussian integration parameters
    printf("\nPerforming Gaussian integration with order %d\n", IntOrder);
    printf("\nBasis set polynomial order: %d\n", BasisOrder-1);
    alloc2d(chi,BasisOrder,BasisOrder);
    LagrangeCoeff(chi);
//    printf("\nLocal integration limits: %d, %d\n\n",xmin,xmax);

//    alloc2d(dchi,BasisOrder,BasisOrder-1);

    printf("\nBasis functions:\n A*x^0 + B*x^1 + C*x^2\n");
    printf("%d polynomials of order %d\n\n",BasisOrder,BasisOrder-1); // BasisOrder is the number of coefficients, including the one for x^0
    printArr(chi,BasisOrder,BasisOrder);

    // First derivative of the basis functions
/*    for (i=0;i<BasisOrder;i++)
    {
        for (j=0;j<BasisOrder-1;j++)
        {
	    dchi[i][j]=chi[i][j+1]*(j+1);
	    PRINT(dbasis(chi,i,1));
        }   
    }   
    printf("\nFirst derivative of the basis functions:\n A*x^0 + B*x^1\n");
    printArr(dchi,BasisOrder,BasisOrder-1); // NOT IN USE ANYMORE
  */  

    printf("\nEnter potential type:\n0: cos(x)\n1: Harmonic oscillator\n2: Particle in a box\n");
    scanf("%d",&pottype);
    switch(pottype) {
	case 0:
	    printf("Factor n in cos(n*x) = ");
	    scanf("%lf",&cosfactor);
	    Vx=Vcos;
	    Sb = gsl_matrix_calloc(MatrixOrder,MatrixOrder);
	    Hb = gsl_matrix_calloc(MatrixOrder,MatrixOrder);
	    C = gsl_matrix_calloc(MatrixOrder,MatrixOrder);
	    E = gsl_vector_calloc(MatrixOrder);

	    wrkEig = gsl_eigen_gensymmv_alloc(MatrixOrder); // for NxN matrix wrk is O(4N)
	    break;
	case 1:
	    MSG(Harmonic oscillator);
	    // Harmonic oscillator
	    XMIN=-5;
	    XMAX=5;
	    Vx=Vharm;

	    // Boundary conditions: the two outer rows and outer columns are set to zero
	    Sb = gsl_matrix_calloc(MatrixOrder-2,MatrixOrder-2);
	    Hb = gsl_matrix_calloc(MatrixOrder-2,MatrixOrder-2);
	    C = gsl_matrix_calloc(MatrixOrder-2,MatrixOrder-2);
	    E = gsl_vector_calloc(MatrixOrder-2);

	    wrkEig = gsl_eigen_gensymmv_alloc(MatrixOrder-2); // for NxN matrix wrk is O(4N)
	break;
	case 2:
	    printf("Box length = ");
	    scanf("%lf",&boxlength);
	    printf("Height of the potential for x < -L/2 = ");
	    scanf("%lf",&boxv1);
	    printf("Height of the potential for x > L/2 = ");
	    scanf("%lf",&boxv2);
	    XMIN=-2*boxlength;
	    XMAX=2*boxlength;
	    Vx=Vbox;

	    // Boundary conditions: the two outer rows and outer columns are set to zero
	    Sb = gsl_matrix_calloc(MatrixOrder-2,MatrixOrder-2);
	    Hb = gsl_matrix_calloc(MatrixOrder-2,MatrixOrder-2);
	    C = gsl_matrix_calloc(MatrixOrder-2,MatrixOrder-2);
	    E = gsl_vector_calloc(MatrixOrder-2);

	    wrkEig = gsl_eigen_gensymmv_alloc(MatrixOrder-2); // for NxN matrix wrk is O(4N)
	    break;
    }

    // Define the elements by their end coordinates, centre and length
    printf("\nThe interval [%f,%f] is split into %d subintervals:\n",XMIN,XMAX,Nelem);
    step[0]=XMIN; // always fixed
    step[Nelem]=XMAX; // the number of coordinates are one more than the number of elements
    for (i=1;i<Nelem;i++) // all equal for now
    {
	step[i]=step[i-1]+(XMAX-XMIN)/Nelem;
    }

    for (i=0;i<Nelem;i++)
    {
	centre[i]=(step[i+1]+step[i])*0.5;
	length[i]=(step[i+1]-step[i])*0.5; // half the length of the interval is necessary
    }

    for (i=0;i<Nelem;i++)
    {
	printf("Element %d. Coordinate: %f, length/2 = %f, centre X0 = %f\n",i,step[i],length[i],centre[i]);
    }
    printf("\n");

// Potential

    printf("Plot potential? Y=1\n");
    scanf("%d",&debug);
    if (debug == 1)
    {
	potfile=fopen("potential.dat","w");
	for (i=0;i<Nelem;i++)
	{
	    fprintf(potfile,"Element %d\n", i);
	    for (k=0;k<IntOrder;k++)
	    {
		fprintf(potfile,"%f\t%lf\n", (length[i]*x[k]+centre[i]), (Vx(length[i],centre[i],x[k])));
	    }
	}
	fclose(potfile);
	system("gnuplot -e \"plot \'potential.dat\' u 1:2 notitle; pause -1\"");
    }

/*****************************************************************************************/

//    printf("Preparing the global matrices of order %d\n\n", MatrixOrder);

    for (n=0;n<Nelem;n++) // Inside each element n:
    {
//	PRINT(coeff);
    // Local kinetic energy operator
    for (i=0;i<BasisOrder;i++)
    {
        for (j=0;j<BasisOrder;j++)
        {
	    temp = 0;
//	    PRINTF(i,d); PRINTF(j,d);
	    for (k=0;k<IntOrder;k++)
	    {
//		PRINTF(k,d);
//		PRINT(dbasis(chi,i,x[k])*dbasis(chi,j,x[k]));
		temp+=w[k]*dbasis(chi,i,x[k])*dbasis(chi,j,x[k]);
//		debugging=w[k]*dbasis(chi,i,x[k])*dbasis(chi,j,x[k]);
//		PRINT(debugging);
//		temp+=debugging;
	    }
//	    PRINT(temp);
//	    debugging=temp/length[n];
//	    PRINT(debugging);
//	    PRINTF(coeff,d);
	    gsl_matrix_set(T,i,j,debugging); // coeff=1/2m from the definition of the Hamiltonian
	    gsl_matrix_set(T,i,j,coeff/length[n]*temp); // coeff=1/2m from the definition of the Hamiltonian
        }   
    } // end T[i][j]   
    printf("Local kinetic energy operator in element %d\n",n);
    printMat(T,BasisOrder,BasisOrder);

    // Potential energy operator
    for (i=0;i<BasisOrder;i++)
    {
        for (j=0;j<BasisOrder;j++)
        {
	    temp = 0;
//	    PRINTF(i,d); PRINTF(j,d);
	    for (k=0;k<IntOrder;k++)
	    {
//		PRINTF(k,d);
//		debugging=w[k] * Vx(length[n],centre[n],x[k]) * basis(chi,i,x[k]) * basis(chi,j,x[k]);
//		PRINT(debugging);
		temp+=w[k] * Vx(length[n],centre[n],x[k]) * basis(chi,i,x[k]) * basis(chi,j,x[k]);  // n - the index of the current element
	    }
//	    debugging=temp;
//	    PRINT(debugging);
	    gsl_matrix_set(V,i,j,temp*length[n]);
        }   
    } // end V[i][j]   
    printf("Potential energy operator in element %d\n",n);
    printMat(V,BasisOrder,BasisOrder);

    // Hamiltonian matrix
    gsl_matrix_add(V,T);    // Use the potential energy elements and subtract the kinetic energy ones from them; V holds the local Hamiltonian
    printf("Local Hamiltonian:\n");
//    printMat(V,BasisOrder,BasisOrder);
    temp = gsl_matrix_get(V,0,0);
    temp += prevH; // adds the bottom right matrix element of the previous interval to the top left element of the current interval to ensure the overlap between them
    prevH = gsl_matrix_get(V,2,2); // store the current bottom right matrix element
    gsl_matrix_set(V,0,0,temp);

    printMat(V,BasisOrder,BasisOrder);

    // Build the global Hamiltonian and global overlap matrix
    //
    // Overlap integrals
    for (i=0;i<BasisOrder;i++)
    {   
        for (j=0;j<BasisOrder;j++)
        {
            temp = 0;
            for (k=0;k<IntOrder;k++)
            {
                temp += w[k]*basis(chi,i,x[k])*basis(chi,j,x[k]);
                gsl_matrix_set(s,i,j,temp*length[n]);
            }
        }
    }    // end s[i][j]   
    temp = gsl_matrix_get(s,0,0);
    temp += prevS; // adds the bottom right matrix element of the previous interval to the top left element of the current interval to ensure the overlap between them
    prevS=gsl_matrix_get(s,2,2); // store the current bottom right matrix element
    gsl_matrix_set(s,0,0,temp);
    printf("Overlap matrix in the element %d\n",n);
    printMat(s,BasisOrder,BasisOrder);

    // Building the global matrices
        for (i=0;i<BasisOrder;i++) // looping over the elements of the local matrices
        {
       	    for (j=0;j<BasisOrder;j++)
	    {
		temp = gsl_matrix_get(V,i,j); // Take the local value of the Hamiltonian
		gsl_matrix_set(H,i+2*n,j+2*n,temp); 
		temp = gsl_matrix_get(s,i,j); // Take the local value of the overlap integral
		gsl_matrix_set(S,i+2*n,j+2*n,temp);
	    }
	} // end H[i][j] and S[i][j]
    
    }   // end loop over the elements

/**************************************************************************************/

    // Boundary conditions
   
// Harmonic oscillator
    // Boundary conditions in the overlap and the Hamiltonian matrices
    for (i=1;i<MatrixOrder-1;i++) // looping over the elements of the overlap matrix ignoring the outer columns and rows
    {
	for (j=1;j<MatrixOrder-1;j++)
	{
	    temp=gsl_matrix_get(S,i,j); // take the inner elements
	    gsl_matrix_set(Sb,i-1,j-1,temp); // shift them by one less
	    temp=gsl_matrix_get(H,i,j); // take the inner elements
	    gsl_matrix_set(Hb,i-1,j-1,temp); // shift them by one less
	}
    } // end boundary conditions

    printf("Print Hamiltonian? Y=1\n");
    scanf("%d",&debug);
    if (debug == 1)
    {
	printMat(H,MatrixOrder,MatrixOrder);
	printf("With boundary conditions:\n");
	printMat(Hb,MatrixOrder-2,MatrixOrder-2);
    }

    printf("Print overlap matrix? Y=1\n");
    scanf("%d",&debug);
    if (debug == 1)
    {
	printMat(S,MatrixOrder,MatrixOrder);
	printf("With boundary conditions:\n");
	printMat(Sb,MatrixOrder-2,MatrixOrder-2);
    }

// Solving the generalized eigenvalue problem
    
    gsl_eigen_gensymmv(Hb,Sb,E,C,wrkEig);
    gsl_eigen_symmv_sort(E,C,GSL_EIGEN_SORT_ABS_ASC);

    printf("Print eigenvalues? Y=1\n");
    scanf("%d",&debug);
    if (debug == 1)
    {

    printf("Eigenvalues:\n");
    for (i=0;i<MatrixOrder-2;i++)
    {
	printf("% f\t",gsl_vector_get(E,i));
    }
    printf("\n\n");
    }

    cfile=fopen("coefficients.dat","w");
    for (i=0;i<MatrixOrder-2;i++)
    {   
        fprintf(cfile,"\n%d\t", line); // The new line is printed in order to separate the eigenvectors from each other; Line gives the consecutive number of the line; it will be the x axis on the plot
	line++;
        for (j=0;j<MatrixOrder-2;j++)
        {
            fprintf(cfile,"% f  ",gsl_matrix_get(C,i,j));
        }
    }   
    fclose(cfile);

    printf("Print eigenvectors? Y=1\n");
    scanf("%d",&debug);
    if (debug == 1)
    {
	printf("Coefficients:\n");
	printMat(C,MatrixOrder-2,MatrixOrder-2);
    }

    printf("Plot eigenvectors: (0 - %d). Enter -1 to exit\n", MatrixOrder-3);
    scanf("%d",&line);
    while (line >= 0)
    {
	sprintf(gnuplot,"gnuplot -e \"plot \'coefficients.dat\' u 1:%d w l notitle; pause -1\"", line+2);
	system(gnuplot);
	printf("Plot eigenvectors: (0 - %d). Enter -1 to exit\n", MatrixOrder-3);
	scanf("%d",&line);
    }

/*    printf("First eigenvector:\n");
    for (i=0;i<MatrixOrder-2;i++)
    {
	printf("%d\t%lf\n",i,gsl_matrix_get(C,0,i));
    }
    printf("\n\n");
*/

    dealloc(chi,BasisOrder);
//    dealloc(dchi,BasisOrder-1);

    gsl_matrix_free(T);
    gsl_matrix_free(V);
    gsl_matrix_free(s);
    gsl_matrix_free(H);
    gsl_matrix_free(S);
    gsl_matrix_free(C);

    gsl_matrix_free(Hb);
    gsl_matrix_free(Sb);
    gsl_vector_free(E);

    gsl_eigen_gensymmv_free(wrkEig);

    free(centre);
    free(length);
    free(step);

    return 0;

}
/*****************************************************************/

int alloc2d(double **array, int nrows, int ncols) 
{
    int i,j;
    for (i=0;i<nrows;i++)
    {
        array[i]=malloc(ncols*sizeof(double));
        if ((array[i]==NULL))
        {
            printf("*** Memory allocation error ***");
            exit(1);
        }
    }   
    for (i=0;i<nrows;i++)
    {   
        for (j=0;j<ncols;j++)
        {
           array[i][j]=0;
        }
    }
    return 0;
}

void printMat(gsl_matrix *array, int nrows, int ncols)
{
    int i,j;
    printf("\n");
    for (i=0;i<nrows;i++)
    {   
        printf("\n");
        for (j=0;j<ncols;j++)
        {
            printf("% f  ",gsl_matrix_get(array,i,j));
        }
    }   
    printf("\n\n");
}

void printArr(double **array, int nrows, int ncols)
{
    int i,j;
    printf("\n");
    for (i=0;i<nrows;i++)
    {   
        printf("\n");
        for (j=0;j<ncols;j++)
        {
            printf("%f\t",array[i][j]);
        }
    }   
    printf("\n\n");
}


void dealloc(double **array,int nrows)
{
    int i;
    for (i=0;i<nrows;i++)
    {
        free(array[i]);
    }
    free(array);
}

void LagrangeCoeff(double **chi)
{
    int i,j;
    int x[BasisOrder]; // the values x0, x1, x2, ... of the Lagrange polynomial
    double denominator;
    double coefficient;
    double sum;

    xmax=(BasisOrder-1)*0.5;
    xmin=-xmax;
    x[0]=xmin; // the smallest value of the local integration interval
    for (i=1;i<BasisOrder;i++)
    {
	x[i]=x[i-1]+1;
    }

    printf("\nLocal interval points: \n");
    for (i=0;i<BasisOrder;i++)
	printf("x[%d]=%d\n",i,x[i]);

    for (i=0;i<BasisOrder;i++)
    {
	PRINTF(i,d);
	denominator=1;
	coefficient=1;
	sum=0;
	for (j=0;j<BasisOrder;j++)
	{
	    if (i != j)
	    {
		denominator*=x[i]-x[j];
		coefficient*=x[j];
		sum+=x[j];
	    } // end if
	} // end for j
	PRINT(denominator);
	PRINT(coefficient);
	PRINT(sum);
	chi[i][0]=coefficient/denominator;
	PRINT(chi[i][0]);
	chi[i][1]=-sum/denominator;
	PRINT(chi[i][1]);
	chi[i][2]=1/denominator;
	PRINT(chi[i][2]);
    } // end for i

    return;
}

double Vcos(double length, double centre, double x)
{
    double result = cos(cosfactor*(centre + length*x));
    return result;
}
double Vharm(double length, double centre, double x)
{
    double result = 0.5*K*(centre + length*x)*(centre + length*x); // K -> K*x^2/2
    return result;
}
double Vbox(double length, double centre, double x)
{
    double result;
    if ( (centre + length*x) < (-boxlength*0.5) )
    {
	result = (boxv1);
    }
    else if ( (centre + length*x) > (boxlength*0.5) )
    {
	result = (boxv2);
    }
    else
    {
	result = 0;
    }
    return result;
}

