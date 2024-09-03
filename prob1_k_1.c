#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct {
    double* v; // Array for velocity
    double* r; // Array for residue
    double* b; // Array for right-hand side of the equation (Au=b)
    int size; // Size of the arrays
} Struct;

double function(int k, double x) {
    int n = 512; // Calculate number of grids at level l
    return sin((k * M_PI * x) / n);
}


void Gauss_Seidel(double V[],double F[])
{
    int i,k,nu;
    double H,alpha;
    H = 1.0/512;
    alpha = 1.0/(2 + H*H);
    k = 0;
    nu = 2;
    do
    {
        for(i=1;i<512;i++)//computing error only at interior points
        {
            V[i] = H*H*alpha*F[i] + alpha*(V[i-1] + V[i+1]);
            //printf("V[%d]=%lf\n",i,V[i]);
        }
     k=k+1;
    }while(k < nu);
}

void relaxation_methods(Struct *l, int i) {
    int j, n;
    double h, alpha;
    n = l[i].size -1;
    h = 1.0 / n;
    alpha = 1.0 / ((h * h) + 2);

    int iteration = 0;
    do {
        for (j = 1; j < n; j++) {
            l[i].v[j] = ((h * h * l[i].b[j]) + l[i].v[j - 1] + l[i].v[j + 1]) * alpha; //Gauss-Seidel
        }
        iteration++;
    } while (iteration < 2); // Simple iteration limit, not a convergence check
}

void restriction_methods(Struct *l, int i) {
    int n = l[i].size -1;
    double h = 1.0 / n;
    
  ///Calculating residue
    for (int j = 1; j < n; j++) {
        l[i].r[j] = l[i].b[j] - ((-l[i].v[j - 1] + (2+(h*h))* l[i].v[j] - l[i].v[j + 1]) / (h * h));
    }
    
  ///Restriction
    for (int j = 1; j <=(n/2) - 1; j++) {
        l[i + 1].b[j] = 0.25 * (l[i].r[(2 * j) - 1] + 2 * l[i].r[2 * j] + l[i].r[(2 * j) + 1]);
    }
    
}

void prolongation_methods(Struct *l, int i) {
    int n = l[i-1].size -1;
    double u[n+1]; 

    for (int j = 0; j <=n; j++) {
        u[j] = 0; // Initialize temporary array
    }

    int m = n/2 -1; 
    for (int j = 0; j <=m; j++) {
        u[2 * j] = l[i].v[j];
        u[2 * j + 1] = 0.5 * (l[i].v[j] + l[i].v[j + 1]);
    }

    // Update the finer grid values
    int g=l[i-1].size -1; // Size of the finer grid
    for (int j = 0; j <=g; j++) {
        l[i-1].v[j] += u[j];
    }

    //free(u); // Free the temporary array
    
}

int main() {  
    int r=512;
	int level = log(r)/log(2) ;
    Struct l[level];
    double h = 1.0 /r ;
    int k = 1;

    for (int i = 0; i < level; i++) {
        
        l[i].size = r / pow(2, i) + 1;

        l[i].v = (double*)malloc(l[i].size * sizeof(double));
        l[i].r = (double*)malloc(l[i].size * sizeof(double));
        l[i].b = (double*)malloc(l[i].size * sizeof(double));

        if (!l[i].v || !l[i].r || !l[i].b) {
            printf("Memory allocation failed at level %d.\n", i);
            return -1; // Early exit
        }
    }
    
    for(int i=0;i<level;i++)
    {	
    	for(int j=0;j<l[i].size;j++)
    	{	
    	l[i].v[j] = 0;
        l[i].r[j] = 0;
        l[i].b[j] = 0;		
    	}
    }


    for (int j = 0; j < l[0].size; j++) {
    	double p=M_PI*k;
        l[0].b[j] = ( pow(p,2)+ 1) * sin((k * M_PI * j) / r );
    }

    int q=l[0].size-1;
    double error=0;
    int iteration=0;
    
FILE *f1= fopen("residue_vs_iteration.dat","w");
fprintf(f1, "Variables = \"iteration\",\"residue\"\n");
do{
    error=0;
    for(int i=1;i<level;i++)
    {	
    	for(int j=0;j<l[i].size;j++)
    	{	
    	l[i].v[j] = 0;
        l[i].r[j] = 0;
      l[i].b[j] = 0;
    			
    	}
    }
    
    for (int i = 0; i < level-1; i++) {
        relaxation_methods(l, i);
        restriction_methods(l, i);
    }
    
    relaxation_methods(l, level-1);
    prolongation_methods(l, level-1);
    
    for(int i=level-2; i>0; i--){
        relaxation_methods(l, i);
        prolongation_methods(l, i);
        
    }
    relaxation_methods(l, 0);
    
    for(int y=1;y<q;y++){
    	error = error + (l[0].r[y])*(l[0].r[y]);
    }
    error = sqrt(error/ ((q) * (q)));
    iteration++;
    
    printf("iteration=%d,residue-norm=%lf\n",iteration,error);
    fprintf(f1,"%d\t%lf\n",iteration,error);
    }while(error>pow(10,-6));
    fclose(f1);

  FILE *f2= fopen("comparision.dat","w");
  fprintf(f2,"\nexact\t computed\n");
  printf("\nexact\t computed\n");
    for (int j = 0; j <=512; j++) {
        double f = function(1,j); 
        fprintf(f2,"exact=%lf, u[%d]=%lf\n", f, j, l[0].v[j]);
        printf("exact=%lf, u[%d]=%lf\n", f, j, l[0].v[j]);
    }
    fclose(f2);
    
    //Solving using Gauss-Seidel
    double V[513], F[513], R[513];
    double H=1.0/512;
    double beta = 1.0 / ((H* H) + 2);
    
    //Initialization
    for(int j=0;j<513;j++){
    	V[j]=0;
    }
    
    for (int j = 0; j < 513; j++) {
    	double p=M_PI*k;
        F[j] = ( pow(p,2)+ 1) * sin((k * M_PI * j) / 512 );
    }
    
    int w=0;
    double error4;
    
    FILE *f3= fopen("Gauss-Seidel.dat","w");
    fprintf(f3,"Variables = \"log(iteration)\",\"log(residue)\"\n");
    do{
    	error4=0;
    	Gauss_Seidel(V,F);
    	//Calculating Residue
        for(int i=1;i<512;i++)
            {
                R[i] = F[i] + ( (V[i+1] + V[i-1]) - (2+H*H)*V[i] )/(H*H);
                //printf("R[%d]=%lf\n",i,R[i]);
            }
            R[0] = 0.0;
            R[512] = 0.0;
        //2-norm residue
    	for(int y=1;y<512;y++){
    		error4 = error4 + (R[y])*(R[y]);
    	}
    	error4 = sqrt(error4/ ((512) * (512)));
        w++;
        //printf("iteration %d\t r_2norm %lf\n",iteration,error4);
        fprintf(f3,"%d\t%lf\n",w,error4);
    }while(error4>pow(10,-6)); //termination condition
 
    // Free allocated memory
    for (int i = 0; i < 8; i++) {
        free(l[i].v);
        free(l[i].r);
        free(l[i].b);
    }

    return 0;
}

