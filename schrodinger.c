/************************************************************************************************
 *                                                                                              *
 *   COMPILED WITH THE LINE:                                                                    *
 *   gcc -Wall -Wextra -o schrodinger.x schrodinger.c -lm -lgsl -lgslcblas && ./schrodinger.x   *
 *                                                                                              *
 *   MAKEFILE:                                                                                  *
 *   make clean;  make VERSION="-DUSERINPUT" && ./schrodinger                                    *
 *   make clean;  make VERSION="-DUSERINPUT -DDEBUGGING" && ./schrodinger                        *
 *                                                                                              *
 *                                                                                              *
 ************************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include <ctype.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include <time.h>

#define pi              (3.141592653589793)
//#define hartree2kJ    (2625.5) // Hartree to kJ/mol
//#define AU            (1822.888486192)
#define hbar            (1.054571800) // *1E-34 [J.s]
//#define hPlanck       (6.62607004) // *1E-34 [J.s]
#define cLight          (2.99792458) // *1E8 [m/s]
#define NA              (6.02214) // *1E23 [mol^-1] // Avogadro's number
#define kB_cm           (0.69503476) // Boltzmann constant in cm^1/K
#define kB_kJ           (0.0083144621) // Boltzmann constant in kJ/mol/K
#define TempK           3 // The temperature for the Boltzmann population of states

#define wavenum2kJ      0.0119627 // 1 cm^-1 to kJ/mol
#define kJ2wavenum      83.593 // 1 kJ/mol to cm^-1
#define kJ2Kelvin       120.274 // 1 kJ/mol to Kelvin
#define kJ2MHz          (2.50607E6) // 1 kJ/mol to MHz
#define amu2kg          1.660539 // *1E-27 // 1 atomic mass unit to kg

#define I_SI            (amu2kg) // *1E-47  = 1.660539E-47 // The inertia is in Å^2, so multiply by amu and convert to meters
#define B_j             (hbar*hbar/I_SI) // *1E-21 = 0.669735E-21  // The rotation constant in J : B = hbar^2 /amu/1E-20 *0.5/inertia. The inertia is in Å^2
#define B_kJ            (B_j*NA*1E-1) // = 0.403324009  // The rotation constant in kJ/mol : B = hbar^2 /amu/1E-20 *NA *1E-3 *0.5/inertia. The inertia is in Å^2
#define B_cm            (hbar/2/pi/cLight/amu2kg*1E3) // = 33.715259 // The rotation constant in cm^-1 : B = hbar/2/pi/c/amu/1E-20/100 *0.5/inertia. The inertia is in Å^2
#define B_MHz           (hbar/2/pi/amu2kg*1E7) // = 1800704.65 // 1E7 comes from the exponents // The rotation constant in cm^-1 : B = B_MHz / inertia
#define conversion      B_kJ
// E = hbar^2/2/I [J] * NA / 1000  = [kJ/mol]
// E = h*c/lambda
// E[J] /h/c = E [cm^-1]
// nu = E[J] / h = [Hz]
//#define bohr          (0.52977)  // 1 Å = 1/0.52977 bohr
#define mC              (12.0107)
#define mH              (1.00798)
#define mN              (14.0031)
//#define mass (6*mC+6*mH) // m(benzene)
//#define mass (6*mC+6*mH) // m(benzene)

// Parallel benzene
//#define rC (1.217/bohr)
//#define rH (2.170/bohr)
//#define IC (4*mC*rC*rC)
//#define IH (4*mH*rH*rH)

// Normal benzene
#define rC (1.895) 
#define rH (2.985)
#define IH (6*mH*rH*rH)
#define IC (6*mC*rC*rC)

// Ammonia molecule:
// #define rH (0.7135/bohr)
// #define rC 0
// #define IC 0
// #define IH (3*mH*rH*rH)
// #define IH           2.8059/I_SI
// #define I_perpallel  4.4128/I_SI

// Heavy Water molecule:
//#define rH (0.75669/bohr)
//#define IC 0
//#define IH (2*2*mH*rH*rH)

#define Inertia (IC+IH)
//#define Inertia_para (IC_para+IH_para)
//#define Inertia_perp (IC_perp+IH_perp)
#define Inertia_perp  2*Inertia_para
#define IntOrder 9 // number of terms in the Gaussian integration

#define MSG(msg) printf( "\n" #msg "\n")
#define DIV printf("\n=======================================================================================\n\n")
#define DEBUG(msg, var, fmt) printf( #msg "\n" #var " = %" #fmt "\n", var)
#define PRINT(var) printf( #var " = %lf\n", var)
#define PRINTF(var, fmt) printf( #var " = %" #fmt "\n", var)
#define PRINT2(var1, fmt1, var2, fmt2) printf( #var1 " = %"#fmt1 "\t ", var1);  printf( #var2 " = %"#fmt2 "\n", var2);   
#define PRINT3(var1, fmt1, var2, fmt2, var3, fmt3) printf( #var1 " = %"#fmt1 "\t ", var1);  printf( #var2 " = %"#fmt2 "\t ", var2);  printf( #var3 " = %"#fmt3 " \n", var3); 

// Version 
#ifdef USERINPUT
  #define USERVERSION
#else
  #ifdef FILEINPUT
    #define FILEVERSION
  #endif
#endif

//const double w[9]={0.330239355001260, 0.312347077040003, 0.260610696402935, 0.180648160694857, 0.081274388361574, 0.312347077040003, 0.260610696402935, 0.180648160694857, 0.081274388361574};  // Gaussian weight factors
//const double x[9]={0.000000000000, 0.324253423403809, 0.613371432700590, 0.836031107326636, 0.968160239507626, -0.324253423403809, -0.613371432700590, -0.836031107326636, -0.968160239507626};  // x value for each Gaussian weight factor
const double w[9]={ 0.081274388361574, 0.180648160694857,   0.260610696402935, 0.312347077040003, 0.330239355001260, 0.312347077040003, 0.260610696402935, 0.180648160694857, 0.081274388361574};  // Gaussian weight factors
const double x[9]={-0.968160239507626, -0.836031107326636, -0.613371432700590, -0.324253423403809, 0.000000000000, 0.324253423403809, 0.613371432700590, 0.836031107326636, 0.968160239507626};  // x value for each Gaussian weight factor
 
int xmax;  
int xmin; 

unsigned int Nelem=1; 
unsigned int BasisOrder=3; 

unsigned int boundStates=0; 
// Polynomial order: (BasisOrder-1);  for a quadratic polynomial x^0, x^1, x^2 => 3 coefficients 
// !!! only use Lagrange polynomials of even power (odd number of basis function)
double inertia = Inertia;  // by default take the defined value for benzene
double coeff = 0.5;  //  hbar^2/2 in atomic units

// Parameters for the potential
unsigned int pottype=0;     // which potential function to use
int cosfactor=1;   // n in A*cos(n*x)
double amplitude=1;   // A in A*cos(n*x); 
double xmaxCoeff=1; 
double boxlength=1;   // L 
double boxv1=1;       // height for x < -L/2
double boxv2=1;       // height for x > L/2
double K=1;      // force constant of the harmonic oscillator, K*x^2/2

double XMIN = 0; 
double XMAX = 2*pi; 
unsigned int MatrixOrder;  
unsigned int BoundaryMatrixOrder;  // the order of the matrices after the boundary conditions are applied

double **chi;    // Coefficients of the basis set polynomials
double *centre;  // Centre of the elements on the global coordinate axis
double *length;  // Length of the elements on the global coordinate axis
double *coord;    // The coordinate of the beginning of the element on the global coordinate axis

// For calculation time
clock_t t1,t2; 
double tcpu; 

// GSL matrices
gsl_matrix *T, *V, *s;  // local matrices for a single element
gsl_matrix *H, *S, *C;  // global matrices for the whole interval
gsl_matrix *Hb, *Sb;  // global matrices after the boundary conditions are applied
gsl_matrix *d, *D, *Db; // local and global matrices for the dipole moment operator
gsl_vector *E; 
gsl_eigen_gensymmv_workspace *wrkEig;  // for generalized eigenvalue problem

// Subroutines (defined at the end of the code)
int alloc2d(double **, unsigned int, unsigned int); 
void dealloc(double **,unsigned int); 
void printArr(double **, unsigned int, unsigned int); 
void printMat(gsl_matrix *, unsigned int, unsigned int); 
void printVec(gsl_vector *, unsigned int); 
void LagrangeCoeff(double **);  
void readInput(void); 
void printOutput(void); 
void printEV(void); 
void printWF(void); 
void printPot(void); 
void printDipole(void); 
void BoltzmannPop(void); 

int readUserInput(char *); 

double basis(int, double); 
double dbasis(int, double); 
double x2X(int, double); 
double X2x(int, double); 

double Vbox(double, double, double); 
double Vharm(double, double, double); 
double Vcos(double, double, double); 
double muCos(int, double); 
double muX(int, double); 

double dipoleMom(int, int); 

double integration(void);
double intFunction(double);

double (*Vx)(double,double,double);  // function pointer for the potential of choice
double (*Mu)(int, double); 

/*******************************************************************/

// int main(int argc, char *argv[])
int main(void)
{
    // Read input parameters from file
    readInput(); 

#ifdef USERINPUT
    printf("\nEnter number of elements: ");  // number of elements the whole interval (0..2pi) is separated into
    Nelem = readUserInput("%d"); 
#else
    printf("Using %d elements as defined in the input file.\n", Nelem); 
#endif

    MatrixOrder = ((BasisOrder-1)*Nelem) + 1;  // BasisOrder - 1 because the elements overlap;  + 1 comes from the lack of overlap for the first element

    chi = malloc(BasisOrder*sizeof(double *)); 
    centre = malloc((Nelem)*sizeof(double)); 
    length = malloc((Nelem)*sizeof(double)); 
    coord = malloc((Nelem+1)*sizeof(double)); 

    if ((chi == NULL)||(centre == NULL)||(length == NULL)||(coord == NULL))
    {
	puts("Memory allocation error"); 
	exit(EXIT_FAILURE); 
    }

    int userinput = 0; 

    double temp = 0;  // used for the summation in some loops 
    double prevH = 0;  // ensure the overlap between two consecutive integration intervals
    double prevS = 0;  
    double prevD = 0;
    size_t i,j,k,n;  // loop indices
#ifdef USERINPUT
    char gnuplot[500] = "gnuplot";  // the string to execute a system call to gnuplot
#endif

    T = gsl_matrix_calloc(BasisOrder,BasisOrder); 
    V = gsl_matrix_calloc(BasisOrder,BasisOrder); 
    s = gsl_matrix_calloc(BasisOrder,BasisOrder); 

    d = gsl_matrix_calloc(BasisOrder,BasisOrder); 

    H = gsl_matrix_calloc(MatrixOrder,MatrixOrder); 
    S = gsl_matrix_calloc(MatrixOrder,MatrixOrder); 

    D = gsl_matrix_calloc(MatrixOrder,MatrixOrder); 

#ifdef DEBUGGING
    printf("Local %lux%lu matrices allocated\n", s->size1, s->size2); 
    printf("Global %dx%d matrices allocated\n", MatrixOrder, MatrixOrder); 
#endif

    // *wrkEig, *Hb, *Sb, *C, *E are allocated later after the potential type is chosen because the dimension depends on the boundary conditions

    /*****************************************************************************************/

    // Gaussian integration parameters
    alloc2d(chi,BasisOrder,BasisOrder); 
    LagrangeCoeff(chi); 

#ifdef DEBUGGING
    printf("\nPerforming Gaussian integration with order %d\n", IntOrder); 
    printf("\nBasis set polynomial order: %d\n", BasisOrder-1); 
    printf("\nBasis functions:\n A*x^0 + B*x^1 + C*x^2 + ... \n"); 
    printf("%d polynomials of order %d\n\n",BasisOrder,BasisOrder-1);  // BasisOrder is the number of coefficients, including the one for x^0
    printArr(chi,BasisOrder,BasisOrder); 
#endif

#ifdef USERINPUT
    printf("\nEnter potential type:\n0: cos(x)\n1: Harmonic oscillator\n2: Particle in a box\n"); 
    pottype = readUserInput("%d"); 
#else
    printf("Using the potential type defined in the input file.\n"); 
#endif

    switch(pottype) {
	case 0:
	    MSG(Cosine potential); 
#ifdef USERINPUT
	    //	    printf("Factor in cos(n*x)\nn = "); 
	    //	    scanf("%d",&cosfactor); 
	    printf("Amplitude in A*cos(n*x) in kJ/mol\nA = "); 
	    scanf("%lf",&amplitude); 
	    PRINT(amplitude); 
	    PRINT(amplitude*B_kJ); 
	    PRINT(amplitude*B_cm); 
	    //	    printf("Integration range [0:x*2*pi]\nx = "); 
	    //	    scanf("%lf",&xmaxCoeff);  
#endif

	    coeff=0.5/Inertia;  // once the inertia value has been determined, calculate it as 1/2I
	    MSG(); 
	    printf("I = %lf a.u.\n", inertia); 
	    printf("I = %lf E-47 kg m^2\n",inertia*I_SI); 
	    puts(""); 
	    printf("B = %lf E-21 J\n",B_j*0.5/inertia); 
	    printf("B = %lf kJ/mol\n",B_kJ*0.5/inertia); 
	    printf("B = %lf cm^-1\n",B_cm*0.5/inertia); 
	    printf("B = %lf MHz\n",B_MHz*0.5/inertia); 
	    PRINT(coeff); 

	    XMIN=0; 
	    XMAX=2*pi*xmaxCoeff; 
	    //	    xmin=-pi*xmaxCoeff; 
	    //	    xmax=pi*xmaxCoeff; 
	    Vx=Vcos; 
	    Mu = muCos; 
	    // Boundary conditions: the order of the array is one less, then the periodicity is ensured by copying the removed elements to the first row and column respectively, and putting the sum of the [0][0] and [MatrixOrder][MatrixOrder] element in [0][0]
	    BoundaryMatrixOrder=MatrixOrder-1; 
	    break; 
	case 1:
	    MSG(Harmonic oscillator); 
	    XMIN=-10; 
	    XMAX=10; 
	    Vx=Vharm; 
	    Mu = muX; 

	    // Boundary conditions: the two outer rows and outer columns are set to zero
	    BoundaryMatrixOrder=MatrixOrder-2; 
	    break; 
	case 2:
	    MSG(Particle in a box); 
	    printf("Box length = "); 
	    scanf("%lf",&boxlength); 
	    printf("Height of the potential for x < -L/2 = "); 
	    scanf("%lf",&boxv1); 
	    printf("Height of the potential for x > L/2 = "); 
	    scanf("%lf",&boxv2); 
	    XMIN=-2*boxlength; 
	    XMAX=2*boxlength; 
	    Vx=Vbox; 
	    Mu = muX; 

	    // Boundary conditions: the two outer rows and outer columns are set to zero
	    BoundaryMatrixOrder=MatrixOrder-2; 
	    break; 
    } // end switch case

    // Allocate memory for the array after the application of boundary conditions
    Sb = gsl_matrix_calloc(BoundaryMatrixOrder,BoundaryMatrixOrder); 
    Hb = gsl_matrix_calloc(BoundaryMatrixOrder,BoundaryMatrixOrder); 
    Db = gsl_matrix_calloc(BoundaryMatrixOrder,BoundaryMatrixOrder); 
    C = gsl_matrix_calloc(BoundaryMatrixOrder,BoundaryMatrixOrder); 
    E = gsl_vector_calloc(BoundaryMatrixOrder); 

    wrkEig = gsl_eigen_gensymmv_alloc(BoundaryMatrixOrder);  // for NxN matrix wrk is O(4N)

    // Define the elements by their end coordinates, centre and length
    coord[0]=XMIN;  // always fixed
    coord[Nelem]=XMAX;  // the number of coordinates are one more than the number of elements
    for (i=1; i<Nelem; i++) // all equal for now
    {
	coord[i] = coord[i-1] + (XMAX-XMIN)/Nelem; 
    }

    for (i=0; i<Nelem; i++)
    {
	centre[i] = (coord[i+1] + coord[i])*0.5; 
	length[i] = (coord[i+1] - coord[i])*0.5;  // half the length of the interval is necessary;  
    }

#ifdef DEBUGGING
    printf("\nThe interval [%f,%f] is split into %d subintervals\n",XMIN,XMAX,Nelem); 
    for (i=0; i<Nelem; i++)
    {
	printf("Element %lu. Coordinate: %f, length/2 = %f, centre X0 = %f\n",i,coord[i],length[i],centre[i]); 
    }
    printf("\n\n"); 
    printf("Gaussian points on the global coordinate axis: \n"); 
    for (i=0; i<Nelem; i++)
    {
	for (k=0; k<IntOrder; k++)
	{
	    printf("Elem %d   x[%d] = % lf    X = % lf\n", i,k,x[k],x2X(i,x[k])); 
	}  // end for k
    } // end for i
#endif

    // Plot the potential
#ifdef USERINPUT
    printPot(); 
#endif

    /*****************************************************************************************/

#ifdef DEBUGGING
    printf("Preparing the global matrices of order %u\n\n", MatrixOrder); 
#endif

    for (n=0; n<Nelem; n++) // Inside each element n:
    {
	// Local kinetic energy operator
	for (i=0; i<BasisOrder; i++)
	{
	    for (j=0; j<BasisOrder; j++)
	    {
		temp = 0; 
		for (k=0; k<IntOrder; k++)
		{
		    temp += w[k]*dbasis(i,x[k])*dbasis(j,x[k]); 
		}

		gsl_matrix_set(T,i,j,coeff/length[n]*temp);  // coeff = 1/2m from the definition of the Hamiltonian
		//when integration is from -2..2, the length scale difference is two times less than when it is from -1..1 
		// the division by length[n] comes from the derivatives of the basis functions 
	    }      
	} // end T[i][j]   
#ifdef DEBUGGING
	printf("Local kinetic energy operator in element %lu\n",n); 
	printMat(T,BasisOrder,BasisOrder); 
#endif

	// Potential energy operator
	for (i=0; i<BasisOrder; i++)
	{
	    for (j=0; j<BasisOrder; j++)
	    {
		temp = 0; 
		for (k=0;  k<IntOrder;  k++)
		{
		    temp += w[k] * Vx(length[n],centre[n],x[k]) * basis(i,x[k]) * basis(j,x[k]);   // n - the index of the current element
		} // end for k

		gsl_matrix_set(V,i,j,temp*length[n]); 
	    } // end for j  
	} // end V[i][j]   

#ifdef DEBUGGING
	printf("Potential energy operator in element %lu\n",n); 
	printMat(V,BasisOrder,BasisOrder); 
#endif

	// Hamiltonian matrix
	gsl_matrix_add(V,T);     // Use the potential energy elements and subtract the kinetic energy ones from them;  V holds the local Hamiltonian

#ifdef DEBUGGING
	printf("Local Hamiltonian:\n"); 
	printMat(V,BasisOrder,BasisOrder); 
#endif

	temp = gsl_matrix_get(V,0,0); 
	// adds the bottom right matrix element of the previous interval to the top left element of the current interval to ensure the overlap between them
	temp += prevH;  
	prevH = gsl_matrix_get(V,BasisOrder-1,BasisOrder-1);  // store the current bottom right matrix element
	gsl_matrix_set(V,0,0,temp); 

#ifdef DEBUGGING
	printf("Hamiltonian matrix in element %lu\n",n); 
	printMat(V,BasisOrder,BasisOrder); 
#endif

	// Overlap integrals
	for (i=0; i<BasisOrder; i++)
	{   
	    for (j=0; j<BasisOrder; j++)
	    {
		temp = 0; 
		for (k=0; k<IntOrder; k++)
		{
		    temp += w[k]*basis(i,x[k])*basis(j,x[k]); 
		}
		gsl_matrix_set(s,i,j,temp*length[n]); 
	    }
	}    // end s[i][j]   

#ifdef DEBUGGING
	printf("Overlap matrix in the element %lu\n",n); 
	printMat(s,BasisOrder,BasisOrder); 
#endif

	temp = gsl_matrix_get(s,0,0); 
	// add the bottom right matrix element of the previous interval to the top left element of the current interval to ensure the overlap between them
	// prevS = 0 in the beginning, so nothing is added to the [0][0] element
	temp += prevS;  
	prevS = gsl_matrix_get(s,BasisOrder-1,BasisOrder-1);  // store the current bottom right matrix element
	gsl_matrix_set(s,0,0,temp); 

#ifdef DEBUGGING
	printf("Overlap matrix in the element after changing the value of the element\n"); 
	printMat(s,BasisOrder,BasisOrder); 
#endif


	// Dipole moment matrices
	for (i=0; i<BasisOrder; i++)
	{   
	    for (j=0; j<BasisOrder; j++)
	    {
		temp = 0; 
		//PRINT2(i, d, j, d);
		for (k=0; k<IntOrder; k++)
		{
		    temp += w[k]*basis(i,x[k])*basis(j,x[k])*Mu(n, x[k]);
		    //PRINT3(basis(i,x[k]), lf, basis(j,x[k]), lf, Mu(n, x[k]), lf);
		    //PRINT(w[k]*basis(i,x[k])*basis(j,x[k])*Mu(n, x[k]));
		}
		gsl_matrix_set(d,i,j,temp*length[n]); 
	    }
	}    // end d[i][j]   

#ifdef DEBUGGING
	printf("Dipole moment matrix in the element %lu\n",n); 
	printMat(d,BasisOrder,BasisOrder); 
#endif

	temp = gsl_matrix_get(d,0,0); 
	// add the bottom right matrix element of the previous interval to the top left element of the current interval to ensure the overlap between them
	// prevD = 0 in the beginning, so nothing is added to the [0][0] element
	temp += prevD;  
	prevD = gsl_matrix_get(d,BasisOrder-1,BasisOrder-1);  // store the current bottom right matrix element
	gsl_matrix_set(d,0,0,temp); 

#ifdef DEBUGGING
	printf("Dipole moment in the element after changing the value of the element\n"); 
	printMat(d,BasisOrder,BasisOrder); 
#endif

	// Building the global matrices
	for (i=0; i<BasisOrder; i++) // looping over the elements of the local matrices
	{
	    for (j=0; j<BasisOrder; j++)
	    {
		temp = gsl_matrix_get(V,i,j);  // Take the local value of the Hamiltonian
		gsl_matrix_set(H, i+(BasisOrder-1)*n, j+(BasisOrder-1)*n, temp);  
		temp = gsl_matrix_get(s,i,j);  // Take the local value of the overlap integral
		gsl_matrix_set(S, i+(BasisOrder-1)*n, j+(BasisOrder-1)*n, temp); 

		temp = gsl_matrix_get(d,i,j);  // Take the local value of the dipole moment matrix
		gsl_matrix_set(D, i+(BasisOrder-1)*n, j+(BasisOrder-1)*n, temp); 
	    }
	} // end H[i][j] and S[i][j]


    }   // end loop over the elements

#ifdef DEBUGGING
    printf("Global matrices built\n"); 
#endif

    /**************************************************************************************/

    // Boundary conditions

    if ((pottype==1)||(pottype==2))
    {
	// Harmonic oscillator or Particle in a box
	// Boundary conditions in the overlap and the Hamiltonian matrices
	for (i=1; i<MatrixOrder-1; i++) // looping over the elements of the overlap matrix ignoring the outer columns and rows
	{
	    for (j=1; j<MatrixOrder-1; j++)
	    {
		temp = gsl_matrix_get(S,i,j);  // take the inner elements
		gsl_matrix_set(Sb,i-1,j-1,temp);  // shift them by one less
		temp = gsl_matrix_get(H,i,j);  // take the inner elements
		gsl_matrix_set(Hb,i-1,j-1,temp);  // shift them by one less

		temp = gsl_matrix_get(D,i,j);  // take the inner elements
		gsl_matrix_set(Db,i-1,j-1,temp);  // shift them by one less
	    }
	} // end boundary conditions
    } // end if pottype HarmOsc or Box

    else if (pottype==0)
    {
	// cos(x)
	// copy the elements of the matrix to the one with boundary conditions (bottom row and right-most column removed)
	for (i=0; i<MatrixOrder-1; i++) 
	{
	    for (j=0; j<MatrixOrder-1; j++)
	    {
		temp=gsl_matrix_get(S,i,j);  
		gsl_matrix_set(Sb,i,j,temp);  
		temp=gsl_matrix_get(H,i,j);  
		gsl_matrix_set(Hb,i,j,temp);  

		temp=gsl_matrix_get(D,i,j);  
		gsl_matrix_set(Db,i,j,temp);  
	    }  // end for j
	} // end for i

	// sum [0][0] and [MatrixOrder-1][MatrixOrder-1] and put the result in [0][0] of the new matrices
	temp=gsl_matrix_get(S,0,0); 
	temp+=gsl_matrix_get(S,MatrixOrder-1,MatrixOrder-1); 
	gsl_matrix_set(Sb,0,0,temp); 

	temp=gsl_matrix_get(H,0,0); 
	temp+=gsl_matrix_get(H,MatrixOrder-1,MatrixOrder-1); 
	gsl_matrix_set(Hb,0,0,temp); 

	temp=gsl_matrix_get(D,0,0); 
	temp+=gsl_matrix_get(D,MatrixOrder-1,MatrixOrder-1); 
	gsl_matrix_set(Db,0,0,temp); 

	// Add the periodic boundaries
	j = MatrixOrder-BasisOrder;  // The index of the second array element of the last local matrix in the global matrix
	for (i=j; i<MatrixOrder-1; i++)
	{
	    // The elements copied to the first row of the global array:
	    temp=gsl_matrix_get(S,i,MatrixOrder-1);  // The elements located on the bottom row
	    gsl_matrix_set(Sb,0,i,temp);  // Move those elements to the top row, same column
	    temp=gsl_matrix_get(H,i,MatrixOrder-1);  // The elements located on the bottom row
	    gsl_matrix_set(Hb,0,i,temp);  // Move those elements to the top row, same column
	    temp=gsl_matrix_get(D,i,MatrixOrder-1);  // The elements located on the bottom row
	    gsl_matrix_set(Db,0,i,temp);  // Move those elements to the top row, same column

	    temp=gsl_matrix_get(S,MatrixOrder-1,i);  // The elements located on the right-most column
	    gsl_matrix_set(Sb,i,0,temp);  // Move those elements to the first column, same row
	    temp=gsl_matrix_get(H,MatrixOrder-1,i);  // The elements located on the right-most column
	    gsl_matrix_set(Hb,i,0,temp);  // Move those elements to the first column, same row
	    temp=gsl_matrix_get(D,MatrixOrder-1,i);  // The elements located on the right-most column
	    gsl_matrix_set(Db,i,0,temp);  // Move those elements to the first column, same row
	} // end for i
    } // end else if pottype = cos

#ifdef DEBUGGING
    puts("Periodic boundaries done"); 

    printf("Print Hamiltonian? Y=1\n"); 
    userinput = readUserInput("%d"); 
    if (userinput == 1)
    {
	printMat(H,MatrixOrder,MatrixOrder); 
	printf("With boundary conditions:\n"); 
	printMat(Hb,BoundaryMatrixOrder,BoundaryMatrixOrder); 
    }

    printf("Print overlap matrix? Y=1\n"); 
    userinput = readUserInput("%d"); 
    if (userinput == 1)
    {
	printMat(S,MatrixOrder,MatrixOrder); 
	printf("With boundary conditions:\n"); 
	printMat(Sb,BoundaryMatrixOrder,BoundaryMatrixOrder); 
    }

    printf("Print dipole moment matrix? Y=1\n"); 
    userinput = readUserInput("%d"); 
    if (userinput == 1)
    {
	printMat(D,MatrixOrder,MatrixOrder); 
	printf("With boundary conditions:\n"); 
	printMat(Db,BoundaryMatrixOrder,BoundaryMatrixOrder); 
    }

#endif

    printf("Solving the generalized eigenvalue problem\n\n"); 
    // Solving the generalized eigenvalue problem
    // NOTE: The solver destroys array Hb. Sb contains its Cholesky decomposition!!!

    // Start measuring time
    t1=clock(); 
    gsl_eigen_gensymmv(Hb,Sb,E,C,wrkEig); 
    // Stop measuring time
    t2=clock(); 
    tcpu=((double) (t2-t1))/CLOCKS_PER_SEC; 
    printf("\n%s%f%s\n","CPU time: ",tcpu," s"); 

    gsl_eigen_symmv_sort(E,C,GSL_EIGEN_SORT_VAL_ASC);  // Sort the eigenvalues in ascending order

    printf("Eigenvalues, [a.u.]:\t"); 
    for (i=0; (i<8 && i<BoundaryMatrixOrder); i++)
    {
	printf("% f\t",gsl_vector_get(E,i)); //*B_kJ); 
    }
    printf("\n\n"); 

    printWF(); 
    printEV(); 
    
    if (pottype == 0)
    {
	BoltzmannPop(); 
    }

    for (i = 0;  ( (i < 6) && (i < BoundaryMatrixOrder));  i++)
    {
	for (j = 0;  j<=i;  j++)
	{
	    dipoleMom(j,i); 
	}
    } 

    integration();

#ifdef DEBUGGING
    printf("Print eigenvectors? Y=1\n"); 
    userinput = readUserInput("%d"); 
    if (userinput == 1)
    {
	printf("Coefficients:\n"); 
	printMat(C,BoundaryMatrixOrder,BoundaryMatrixOrder); 
    }   
#endif

#ifdef USERINPUT1
    if (pottype == 0)
    {
	printf("\nPlot the energy states?\n"); 
	userinput = readUserInput("%d"); 
	if (userinput == 1)
	{
	    //	    sprintf(gnuplot,"gnuplot -e \"set xrange [%lf:%lf];  unset xtics;  set yrange [%lf:%lf];  set ylabel \'E, [kJ/mol]\';  plot \'eigenvalues.dat\' u 1:2 notitle, \'potential.dat\' u 1:(\\$2*%lf) w l lw 3 notitle;  pause -1\"", XMIN, XMAX, -amplitude*1.1, amplitude*1.1, B_kJ); 
	    sprintf(gnuplot,"gnuplot -e \"set xrange [%lf:%lf];  unset xtics;  set yrange [%lf:%lf];  set ylabel \'E, [kJ/mol]\';  plot (\\\"<awk \'{print( \'0\', \\$3, \'10\')}\' eigenvalues.dat\\\" ) u 1:2:3 with xerrorbars notitle, \'potential.dat\' u 1:(\\$2*%lf) w l lw 3 notitle;  pause -1\"", XMIN, XMAX, -amplitude*1.2, amplitude*1.2, B_kJ); 
	    system(gnuplot); 
	} // end if userinput == 1
    } // end if pottype==1

    userinput = 0; 
    while (userinput >= 0)
    {
	printf("\nPlot eigenvectors: (0 - %d). Enter -1 to exit\n", BoundaryMatrixOrder-1); 
	userinput = readUserInput("%d"); 
	if (userinput < BoundaryMatrixOrder)
	{
	    sprintf(gnuplot,"gnuplot -e \"unset xtics;  set ylabel \'Wavefunction coefficients\';  plot \'coefficients.dat\' u 1:%d w l notitle;  pause -1\"", userinput+2); 
	    system(gnuplot); 
	}
    } // end while

#endif

#ifdef BATCHRUN
    printOutput(); 
#endif

    /*****************************************************************/

    // Clean up memory
    dealloc(chi,BasisOrder); 

    gsl_matrix_free(T); 
    gsl_matrix_free(V); 
    gsl_matrix_free(s); 
    gsl_matrix_free(H); 
    gsl_matrix_free(S); 
    gsl_matrix_free(C); 

    gsl_matrix_free(d); 
    gsl_matrix_free(D); 

    gsl_matrix_free(Hb); 
    gsl_matrix_free(Sb); 
    gsl_matrix_free(Db); 

    gsl_vector_free(E); 

    gsl_eigen_gensymmv_free(wrkEig); 

    free(centre); 
    free(length); 
    free(coord); 

    return 0; 
}

/*****************************************************************/

int alloc2d(double **array, unsigned int nrows, unsigned int ncols) 
{
    size_t i,j; 
    for (i=0; i<nrows; i++)
    {
        array[i]=malloc(ncols*sizeof(double)); 
        if ((array[i]==NULL))
        {
            printf("*** Memory allocation error ***"); 
            exit(EXIT_FAILURE); 
        }
    }   
    for (i=0; i<nrows; i++)
    {   
        for (j=0; j<ncols; j++)
        {
           array[i][j]=0; 
        }
    }
    return 0; 
}

void printMat(gsl_matrix *array, unsigned int nrows, unsigned int ncols)
{
    size_t i,j; 
    printf("\n"); 
    for (i=0; i<nrows; i++)
    {   
        printf("\n"); 
        for (j=0; j<ncols; j++)
        {
            printf("% f  ",gsl_matrix_get(array,i,j)); 
        }
    }   
    printf("\n\n"); 
    return; 
}

void printVec(gsl_vector *vector, unsigned int n)
{
    size_t i; 
    printf("\n"); 
    for (i=0; i<n; i++)
    {   
        printf("\n"); 
            printf("% f  ",gsl_vector_get(vector, i)); 
    }   
    printf("\n\n"); 
    return; 
}

void printArr(double **array, unsigned int nrows, unsigned int ncols)
{
    size_t i,j; 
    printf("\n"); 
    for (i=0; i<nrows; i++)
    {   
        printf("\n"); 
        for (j=0; j<ncols; j++)
        {
            printf("% f\t",array[i][j]); 
        }
    }   
    printf("\n\n"); 
    return; 
}


void dealloc(double **array, unsigned int nrows)
{
    size_t i; 
    for (i=0; i<nrows; i++)
    {
        free(array[i]); 
    }
    free(array); 
    return; 
}

inline double basis(int fn, double x)
{
    int power; 
    double result=0; 
    for (power=0;  power<BasisOrder;  power++)
    {
	result += chi[fn][power]*pow(x,power); 
    }
    return result; 
}
inline double dbasis(int fn, double x)
{
    int power; 
    double result=0; 
    for (power=1;  power < BasisOrder;  power++) // power=0 => the whole term is 0
    {
	result += power*chi[fn][power]*pow(x,power-1); 
    }
    return result; 
}

void LagrangeCoeff(double **chi)
{
    int fn,i,j; 
    int order;  // i has to be signed because in the loop it is compared if equal to zero
    int xIdx=0;  // Keeps track of the current value of x
    int xLagr[BasisOrder];  // the values x0, x1, x2, ... of the Lagrange polynomial
    double denominator; 

    xmax=(BasisOrder-1)*0.5; 
    xmin=-xmax; 
    xLagr[0]=xmin;  // the smallest value of the local integration interval
    for (i=1; i<BasisOrder; i++)
    {
	xLagr[i]=xLagr[i-1]+1; 
    }

#ifdef DEBUGGING
    printf("\nRoots of the Lagrange polynomial: \n"); 
    for (i=0; i<BasisOrder; i++)
	printf("x[%d]=%d\n",i,xLagr[i]); 
#endif

    // Calculating the coefficients
    for (fn=0;  fn<BasisOrder;  fn++) // go along the basis functions
    {
//	puts(""); 
//	PRINTF(fn,d); 

//	set the case for the polynomial order 1: 
	order=1; 
//	PRINTF(order,d); 
	xIdx=fn+order;  // Starting value is 0, but for chi[0] the value xLagr[0] is not present in the numerator, so increase by one 
	// check if it is necessary to overflow the index 
	if (xIdx >= BasisOrder) // there is no element xLagr[BasisOrder]
	{
	    xIdx = xIdx - BasisOrder; 
	} //end if
//	PRINTF(xIdx,d); 

	denominator=xLagr[fn]-xLagr[xIdx]; 

	chi[fn][1]=1; 
	chi[fn][0]=-xLagr[xIdx]; 

//	printf("chi[%d][1] = %lf\n", fn, chi[fn][1]); 
//	printf("chi[%d][0] = %lf\n", fn, chi[fn][0]); 

	for (order=2;  order<BasisOrder;  order++) // go along the powers of the polynomial to scan from power 1 to power BasisOrder-1
	{
//	    PRINTF(order,d); 

	    xIdx=fn+order;  // Starting value is 0, but for chi[0] the value x[0] is not present in the numerator, so increase by one 
	    // check if it is necessary to overflow the index 
	    if (xIdx >= BasisOrder) // there is no element x[BasisOrder]
	    {
		xIdx = xIdx - BasisOrder; 
	    } //end if
//	    PRINTF(xIdx,d); 	    
	    denominator*=xLagr[fn]-xLagr[xIdx]; 
	
	    chi[fn][order]=1;  // The highest order coefficient is always 1
//	    printf("chi[%d][%d] = %lf\n", fn, order, chi[fn][order]); 

	    for (i=order-2;  i>=0;  i-- )
	    {
		chi[fn][i+1] = (-chi[fn][i+1]*xLagr[xIdx] + chi[fn][i]); 
//		printf("chi[%d][%d] = %lf\n", fn, i+1, chi[fn][i+1]); 
	    } // end for i

	    chi[fn][0] = -chi[fn][0]*xLagr[xIdx]; 
//	    printf("chi[%d][0] = %lf\n", fn, chi[fn][0]); 

	} // end for order

//	printf("chi array: \n"); 
//	printArr(chi,BasisOrder,BasisOrder); 

//	PRINT(denominator); 

	// divide the coefficients by the total denominator for the basis function
	for (j=0;  j<BasisOrder;  j++)
	{
	    chi[fn][j]=chi[fn][j]/denominator; 
	} // end for j

    } // end for fn

//    printf("chi array div denominator: \n"); 
	printf("chi array: \n"); 
    printArr(chi,BasisOrder,BasisOrder); 

    return; 
}

double x2X(int element, double x)
{
    double result = length[element]*x + centre[element]; 
    return result; 
}

double X2x(int element, double X)
{
    double result = (X - centre[element]) / length[element]; 
    return result; 
}

double Vcos(double length, double centre, double x)
{
    double result = amplitude/B_kJ*cos(cosfactor*(centre + length*x));  // A*cos(n*X) ;  The unit conversion ensures that the Schrodinger equation is solved in atomic units. Input is read as kJ/mol
    return result; 
}
double Vharm(double length, double centre, double x)
{
    double result = 0.5*K*(centre + length*x)*(centre + length*x);  // K -> K*X^2/2
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

double muCos(int elem, double x)
{
// the dipole moment operator is mu = mu_x + mu_y;  mu_x = r*cos(X) and mu_y = r*sin(X), where r = 1    double result = 
//    double result = cos(centre + length*x) + sin(centre + length*x); 
//    double result = cos(centre + length*x); 
//    double result = cos(x*length[elem] + centre[elem]); 
    double result = cos(x2X(elem, x)); 
    return result;  
}

double muX(int elem, double x)
{
    //double result = (centre[elem] + length[elem]*x)*(centre[elem] + length[elem]*x)/2.0; 
//    double result = (centre[elem] + length[elem]*x);
    //double result = x2X(elem, x)*x2X(elem, x)/2.0; 
    double result = x2X(elem, x); 
    return result; 
}

void readInput(void)
{
    FILE *inpfile; 

    size_t i = 0; 
    size_t nlines = 0;  // counts the number of lines in the input file
    double inputValue = 0; 
    double massC = -1; 
    double massH = -1; 
    char keyword[50] = " "; 
    char junk[10] = " "; 
    char c = 1;  // counts the characters in the file
    printf("Reading input file\n"); 
    
    inpfile=fopen("param.inp","r"); 
    if (inpfile == NULL)
    {   
        fprintf(stderr, "*** Error reading input data\n"); 
    } // end if
    else
    {   
        // count the number of lines in the file
        while (c != EOF)
        {   
            if (c == '\n')
            {   
                nlines++; 
            }   
            c = getc(inpfile); 
        }  // end while
    } // end else

    rewind(inpfile); 

    for (i=0;  i<nlines;  i++)
    {
	fscanf(inpfile,"%s %s %lf", keyword,junk,&inputValue); 
	if (!strcmp(keyword,"polyOrder"))
	{
	    BasisOrder = 1 + (int)inputValue; 
	    PRINTF(BasisOrder,d); 
	}
	else if (!strcmp(keyword,"massC"))
	{
	    massC=inputValue; 
	    PRINT(massC); 
	}
	else if (!strcmp(keyword,"massH"))
	{
	    massH=inputValue; 
	    PRINT(massH); 
	    inertia=6*massC*rC*rC + 6*massH*rH*rH; 
	    PRINT(inertia); 
	}
	else if (!strcmp(keyword,"Nelem")) // if they are not different 
	{
	    Nelem=(int) inputValue; 
	    PRINTF(Nelem,d); 
	}
	else if (!strcmp(keyword,"pottype"))
	{
	    pottype=(int) inputValue; 
	    PRINTF(pottype,d); 
	}
	else if (!strcmp(keyword,"cosfactor"))
	{
	    cosfactor=(int) inputValue; 
	    PRINTF(cosfactor,d); 
	}
	else if (!strcmp(keyword,"amplitude"))
	{
	    amplitude=inputValue; 
	    PRINT(amplitude); 
	}
	else if (!strcmp(keyword,"xmaxCoeff"))
	{
	    xmaxCoeff=inputValue; 
	    PRINT(xmaxCoeff); 
	}
	
    } // end for i
    fclose(inpfile); 
    return; 
} // end readInput


void printOutput(void)
{
    size_t i; 
    FILE *output; 
    output=fopen("output.dat","a"); 
    // fprintf(output,"Nelem\tpottype\tcosfactor\tamplitude\txmaxCoeff\tE0\tE1\tE2\tE3\tE4\tE5\n"); 
    fprintf(output,"%d\t",Nelem); 
    if (pottype==0)
    {
	fprintf(output,"%5.6lf * cos( %d * x)\t",amplitude,cosfactor); 
	fprintf(output,"in [ % lf ;  % lf ]\t", XMIN, XMAX); 
    } // end if
    else if (pottype==1)
    {
	fprintf(output, " %lf * x^2 /2 \t", K); 
	fprintf(output,"in [ % lf ;  % lf ]\t", XMIN, XMAX); 
    } // end else if

    for (i=0;  i<20;  i++)
    {
	fprintf(output,"%lf\t",gsl_vector_get(E,i)*B_cm); 
    }
    fprintf(output,"\n"); 
    fclose(output); 
    return; 
} // end printOutput

void printWF(void)
{
    FILE *cfile; 
    int i,j; 
    int line=0; 
    cfile = fopen("coefficients.dat","w"); 

    // The solver sometimes changes the phase, so the eigenvectors become negative. The absolute values are correct, and so are the eigenvalues
    // when the potential is harmonic oscillator or particle in a box, check if the middle element of the first eigenvector is negative.
    if ( ((pottype==1)||(pottype==2)) && (gsl_matrix_get(C, MatrixOrder>>2, 0) < 0) ) // The bitshift to the right MatrixOrder>>2 is equivalent to division by two
    {
	for (i=0; i<BoundaryMatrixOrder; i++)
	{
	    for (j=0; j<BoundaryMatrixOrder; j++)
	    {
		gsl_matrix_set(C,i,j, -gsl_matrix_get(C,i,j));  // checks the sign of the elements of the eigenvectors;  if they are negative, invert them. The solver may change the phase and the values are of negative sign but correct absolute value
	    } // end for j
	} // end for i
    } // end if

    for (i=0; i<BoundaryMatrixOrder; i++)
    {   
        fprintf(cfile,"\n%d\t", line);  // The new line is printed in order to separate the eigenvectors from each other;  Line gives the consecutive number of the line;  it will be the x axis on the plot
	line++; 
        for (j=0; j<BoundaryMatrixOrder; j++)
        {
            fprintf(cfile,"% f  ", gsl_matrix_get(C,i,j)); 
        }
    } // end for i  
    fclose(cfile); 
} // printWF
    
void printEV(void)
{
    int i=0; 
    int degeneracy = 1; 
    double diff=0; 
    double degeneracyCriterion = 1E-4; 
    FILE *evfile; 
    FILE *degfile; 

    evfile = fopen("eigenvalues.dat","w"); 
    degfile = fopen("degeneracy.tmp","w"); 

    fprintf(evfile,"#\t  E,[H]\t\t\t  E,[kJ/mol]\t\t  E,[cm^-1]\t\tE,[MHz]\t\t\tdE,[H]\t\tdE,[kJ/mol]\tdE,[cm^-1]\tdE,[MHz]\n"); 
    boundStates = 0; 
    for (i=0; i<BoundaryMatrixOrder-1; i++)
    {
	fprintf(evfile,"%d\t% 16.12lf\t",i, gsl_vector_get(E,i));  // In Hartree
	fprintf(evfile,"% 16.12lf\t",gsl_vector_get(E,i)*B_kJ);  // In kJ/mol
	fprintf(evfile,"% 16.12lf\t",gsl_vector_get(E,i)*B_cm);  // In cm^-1
	fprintf(evfile,"% 16.12lf\t",gsl_vector_get(E,i)*B_MHz);  // In MHz
	
	diff = gsl_vector_get(E,i+1)-gsl_vector_get(E,i);  // dE, [H]
	if (diff < degeneracyCriterion)
	{
	    degeneracy++; 
	}
	else
	{
	    fprintf(degfile,"% 22.16lf\t%lf\t%d\t\n",gsl_vector_get(E,i), diff,degeneracy); 
	    degeneracy = 1; 
	    if (gsl_vector_get(E,i) < amplitude/B_kJ)
	    {
		boundStates++; 
	    } // end if energy

	} // end else

	fprintf(evfile,"% lf\t",diff);  // In [H] 
	fprintf(evfile,"% lf\t",diff*B_kJ);  // In kJ/mol
	fprintf(evfile,"% lf\t",diff*B_cm);  // In cm^-1
	fprintf(evfile,"% lf\n",diff*B_MHz);  // In MHz

    } // end for i

    if (pottype == 0)
    {
	printf("%d bound states\n", boundStates); 
    } // end if

    // Print just the energy of the highest energy because energy difference cannot be calculated
    fprintf(evfile,"%d\t% 16.12lf\t",i, gsl_vector_get(E,BoundaryMatrixOrder-1));  // In Hartree
    fprintf(evfile,"% 16.12lf\n",gsl_vector_get(E,BoundaryMatrixOrder-1)*B_kJ);  // In kJ/mol

    fclose(evfile); 
    fclose(degfile); 
    return; 
} // end printEV()

void printPot(void)
{
    FILE *potfile; 
    int i,k; 
    int userinput=0; 

    potfile=fopen("potential.dat","w"); 
    for (i=0; i<Nelem; i++)
    {
	for (k=0; k<IntOrder; k++)
	{
	    fprintf(potfile,"%lf\t%lf\n", x2X(i,x[k]), (Vx(length[i],centre[i],x[k]))); 
	}  // end for k
    } // end for i
    fclose(potfile); 

#ifdef DEBUGGING
    printf("\nPlot potential? Y=1\n"); 
    userinput = readUserInput("%d"); 
    if (userinput == 1)
    {
	system("gnuplot -e \"plot \'potential.dat\' u 1:2 w l lw 2 notitle;  pause -1\""); 
    } // end if
#endif

    return; 
} // end printPot

void printDipole(void)
{
//    FILE *dipfile; 

//    printf("\nCalculating dipole moments of transition\n"); 
//    dipfile=fopen("dipole.dat","w"); 
//    for (l=0;  l<25;  l++)  // transition from state i -> l: 
//    {
//	for (i=l;  i<25;  i++) 
//	{
//	    dipole = 0; 
//	    for (j=0;  j<BoundaryMatrixOrder;  j++)  
//	    {
//		for (k=0;  k<BoundaryMatrixOrder;  k++)
//		{
////		    dipole += gsl_matrix_get(C,l,j) * gsl_matrix_get(C,i,k) * gsl_matrix_get(MxB,j,k);  // overlap between the wavefunctions psi_(i+1) and psi_i
//		    dipole += gsl_matrix_get(MxB,j,k); 
//		} // end for k
//	    } // end for j
//	    fprintf(dipfile,"E%lu->E%lu:\t% lf\t", i, l, dipole); 
//	} // end for i
//	fprintf(dipfile,"\n"); 
//    } // end for l
//    fclose(dipfile); 

//    MSG(\n); 

    return; 
} // end printDipole

void BoltzmannPop()
{
    int i; 
    int degeneracy = 0; 
    int numStates = boundStates + 12; 
    double popSum = 0;  // The sum of the populations
    double energy = 0;  // Read from the file
    double diff = 0;  // Read from file

    double *populations;  // the exponentials exp(-E_i/kB/T)

    FILE *degfile;  // contains the energies and the degeneracies only;  written by printEV()
    FILE *boltzfile;  // contains the above plus N_i/N

    populations=malloc(numStates*sizeof(double));   // the Boltzmann population of states for the bound states and 12 more 
    // Open the file with degeneracies written by the function printEV(). Read the energies and their corresponding degeneracies. Calculate the exp and store the population for each state. After that the file degeneracy.tmp is not needed
    degfile = fopen("degeneracy.tmp","r"); 
    boltzfile = fopen("boltzmann.dat","w"); 

    if ((degfile != NULL)&&(boltzfile != NULL))
    {
	// Read the energies of the states and their degeneracies from degeneracy.tmp. Calculate exp(-E/kB/T) and put it to populations[i]
	popSum = 0; 
	fscanf(degfile,"%lf\t%lf\t%d", &energy, &diff, &degeneracy);  // read the values from the file
	for (i=0;  i< numStates;  i++)
	{
	    populations[i] = degeneracy*exp(-energy*B_cm/kB_cm/TempK); 
	    popSum += populations[i]; 
	    fscanf(degfile,"%lf\t%lf\t%d", &energy, &diff, &degeneracy);  // read the values from the file
	}

	rewind(degfile); 	
	fprintf(boltzfile,"# Rotation barrier: %lf [cm^-1]\n",amplitude/B_kJ*B_cm); 
	fprintf(boltzfile,"# E,[cm^-1]\tdE,[cm^-1]\tg_i\tN_i/N\n#\n"); 

	// Calculate the bound states
	for (i=0;  i<boundStates;  i++)
	{
	    fscanf(degfile,"%lf\t%lf\t%d", &energy, &diff, &degeneracy);  // read the values from the file
	    fprintf(boltzfile,"% 12.6lf\t%10.6lf\t%d\t%7.5lf\n", energy*B_cm, diff*B_cm, degeneracy, populations[i]/popSum); 
	} // end for

	// Print a separating line
	fprintf(boltzfile,"\n"); 

	// Calculate the first 12 states above the rotation barrier
	for (i=boundStates;  i<numStates;  i++)
	{
	    fscanf(degfile,"%lf\t%lf\t%d", &energy, &diff, &degeneracy);  // read the values from the file
	    fprintf(boltzfile,"% 12.6lf\t%10.6lf\t%d\t%7.5lf\n", energy*B_cm, diff*B_cm, degeneracy, populations[i]/popSum); 
	} // end for

	fclose(boltzfile); 
	fclose(degfile); 
	system("rm -rf ./degeneracy.tmp"); 

	free(populations); 
    } // end if file pointers exist
    return; 
} // end BoltzmannPop()

int readUserInput(char * regexp)
{
    int strLength = 64;  // the length of the string that fgets will read from stdin
    char userinput[strLength]; 
    int stat = 0; 
    int num = -1; 
    while (stat != 1)
    {
	fgets(userinput, strLength-1, stdin);  // the length of the string is specified as one less than strLength because of the null character at the end of a string
	stat = sscanf(userinput, regexp, &num);  
    } 
//    PRINTF(num,d); 
    return num; 
}

double dipoleMom(int initial, int final)
{
    int i = 0;
  //  int i,n,k;
    double dipole = 0;
    double temp = 0;
    gsl_vector_view C_init, C_final;
    gsl_vector *product; 
    gsl_vector *d_vec;  

    C_init = gsl_matrix_column(C, initial);
    C_final = gsl_matrix_column(C, final);

#ifdef DEBUGGING
    MSG("C_init");
    printVec(&C_init.vector, BoundaryMatrixOrder);

    MSG("C_final");
    printVec(&C_final.vector, BoundaryMatrixOrder);
#endif

    product = gsl_vector_calloc(BoundaryMatrixOrder); 
    d_vec = gsl_vector_calloc(BoundaryMatrixOrder); 

    for (i = 0; i < BoundaryMatrixOrder; i++)
    {
	gsl_vector_set(d_vec, i, 1.0);
    }
    
    //  int gsl_blas_dgemv (CBLAS_TRANSPOSE_t TransA, double alpha, const gsl_matrix * A, const gsl_vector * x, double beta, gsl_vector * y)
    // These functions compute the matrix-vector product and sum y = alpha op(A) x + beta y, 
    // where op(A) = A, A^T, A^H for TransA = CblasNoTrans, CblasTrans, CblasConjTrans. 
    //
    // Multiply Sb * C_init
    gsl_blas_dgemv(CblasNoTrans, 1.0, Db, &C_init.vector, 0.0, product);
//
//  DEBUGGING: take the overlap matrix only: multiply its elements by a unit vector
//  NOTE: Sb is destroyed by the eigenvalue solver; it contains its Cholesky decomposition
//    gsl_blas_dgemv(CblasNoTrans, 1.0, Sb, d_vec, 0.0, product);
//    gsl_blas_dgemv(CblasNoTrans, 1.0, Db, d_vec, 0.0, product);
//    MSG("Matrix Db");
//    printMat(Db, BoundaryMatrixOrder, BoundaryMatrixOrder);
#ifdef DEBUGGING
    MSG("Product vector");
    printVec(product, BoundaryMatrixOrder);
#endif
//    for (i = 0; i < BoundaryMatrixOrder; i++)
//    {
////	printf("d_vec = %lf\n", gsl_vector_get(d_vec, i));
//	printf("prod = %lf\n", gsl_vector_get(product, i));
//    }

    // Print the product vector between the matrix D and the C_init vector
//    for (i = 0; i < BoundaryMatrixOrder; i++)
//    {
//	PRINT(gsl_vector_get(product, i));
//    }
//
    // Calculate the dot product between C_final and the product:
    // int gsl_blas_ddot (const gsl_vector * x, const gsl_vector * y, double * result)
    // These functions compute the scalar product x^T y for the vectors x and y, returning the result in result. 
    gsl_blas_ddot(&C_final.vector, product, &dipole);

    // DEBUGGING:
//    gsl_blas_ddot(d_vec, product, &dipole);

//    i = 0;
//    prev = 0;
//
//    for (n = 0; n < Nelem; n++)
//    {
//	k = 0;
//	temp = gsl_vector_get(&C_init.vector, i) * Mu(n, x[k]) + prev;
//	gsl_vector_set(d_vec, i, temp);
//	i++;
//
//	for (k = 1; k < IntOrder-1; k++)
//	    //    for (k = 0; k < BoundaryMatrixOrder; k++)
//	{
//	    temp = gsl_vector_get(&C_init.vector, i) * Mu(n, x[k]);
//	    gsl_vector_set(d_vec, i, temp);
//	    i++;
//	}
//
//	k = IntOrder-1;
//	prev = gsl_vector_get(&C_init.vector, i) * Mu(n, x[k]);
//
//    }

//    gsl_blas_dgemv(CblasNoTrans, 1.0, Sb, d_vec, 0.0, product);
//    gsl_blas_ddot(&C_final.vector, product, &dipole);

    PRINT3(initial, d, final, d, dipole, lf);
    
    gsl_vector_free(d_vec);
    gsl_vector_free(product);

    return dipole; 
}


inline double intFunction(double x)
{
    //double result = pow( (-0.998318*x*x + 1.09275) , 2);  
//    double result = -0.998318*x*x + 1.09275;  
    //double result = 0.1675*x*x + 0.825184; 
    double result = exp(-2*x*x)*x*x; 
//    double result = x*x;
    return result; 
}

double integration(void)
{
    int n, k; 
    double temp = 0; 
    double integral = 0; 
    for (n = 0;  n < Nelem;  n++) // Inside each element n:
    {   
        temp = 0; 
        for (k = 0; k < IntOrder; k++)
        {
            temp += w[k]*intFunction(x2X(n,x[k]));   
        } // end for k
        integral += temp*length[n]; 
    } // end for n
    printf("\nInt [%lf, %lf] exp(-2*x*x)*x*x dx = ", XMIN, XMAX); 
    printf("%20.12f\n",integral); 

    return integral; 
}


