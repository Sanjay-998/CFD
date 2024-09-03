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
}

void v_cycle(Struct *l, int initial, int level){
int k=1;
double p=M_PI*k;
for(int j=0;j<l[initial].size;j++) {
	l[initial].b[j] = ( pow(p,2)+ 1) * sin((k * M_PI * j) / (l[initial].size-1) );
}

int q=l[initial].size-1;
double error;
//do{
    error=0;
    for(int i=initial+1;i<level;i++)
    {	
    	for(int j=0;j<l[i].size;j++)
    	{	
    	l[i].v[j] = 0;
        l[i].r[j] = 0;
      l[i].b[j] = 0;
    			
    	}
    }
    
    for (int j = initial; j < level-1; j++) {
        relaxation_methods(l, j);
        restriction_methods(l, j);
    }
    
    for(int j=level-1; j>initial; j--){
        relaxation_methods(l, j);
        prolongation_methods(l, j);
        
    }
    relaxation_methods(l, initial);
    
    /*for(int y=1;y<q;y++){
    	error = error + (l[initial].r[y])*(l[initial].r[y]);
    }
    error = sqrt(error/ ((q) * (q)));*/
    
    //}while(error>pow(10,-6));
}

int main() {  
    int r=512;
	int level = log(r)/log(2) ;
    Struct l[level];
    //double h = 1.0 /2 ;
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


    for (int j = 0; j < l[level-1].size; j++) {
    	double p=M_PI*k;
        l[level-1].b[j] = ( pow(p,2)+ 1) * sin((k * M_PI * j) / 2 );
    }

  ///FMG
  for(int k=level-1;k>0;k--){
  	relaxation_methods(l,k);
  	prolongation_methods(l,k);
  	v_cycle(l,k-1,level);
  }
  //Here 1 FMG completes
  
  int e=l[0].size;
  double error1=0;
  for(int y=1;y<e;y++){
    		error1 = error1 + (l[0].r[y])*(l[0].r[y]);
    	}
    	error1 = sqrt(error1/ ((e-1) * (e-1)));
  
  FILE *f1= fopen("Residue_vs_iteration.dat","w");
  printf("2-norm residue after 1 FMG=%lf\n",error1);
  fprintf(f1,"2-norm residue after 1 FMG=%lf\n",error1);
  
  //Now apply V-cycle
   double error;
   int iteration=0;
   
   printf("Finding 2-norm residue using iterative V-cycle\n");
   fprintf(f1,"Finding 2-norm residue using iterative V-cycle\n");
   
  do{
  	error=0;
  	v_cycle(l,0,level);
  	
  	//int e=l[0].size;
  	for(int y=1;y<e;y++){
    		error = error + (l[0].r[y])*(l[0].r[y]);
    	}
    	error = sqrt(error/ ((e-1) * (e-1)));
    	iteration++;
    	printf("iteration=%d,residue-norm=%lf\n",iteration,error);
    	fprintf(f1,"iteration=%d,residue-norm=%lf\n",iteration,error);
  }while(error>pow(10,-6));
  
  FILE *f2= fopen("Comparision.dat","w");
  printf("\nexact\t computed\n");
  fprintf(f2,"\nexact\t computed\n");
    for (int j = 0; j <=512; j++) {
        double f = function(1,j); 
        //fprintf(f2,"exact=%lf, u[%d]=%lf\n", f, j, l[0].v[j]);
        printf("exact=%lf, u[%d]=%lf\n", f, j, l[0].v[j]);
        fprintf(f2,"exact=%lf, u[%d]=%lf\n", f, j, l[0].v[j]);
    }
    
    fprintf(f1,"2-norm residue=%lf\n",error);
    
    //fclose(f2);

    /*// Free allocated memory
    for (int i = 0; i < 8; i++) {
        free(l[i].v);
        free(l[i].r);
        free(l[i].b);
    }*/

    return 0;
}

