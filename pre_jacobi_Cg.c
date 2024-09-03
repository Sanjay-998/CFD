#include <stdlib.h>
#include <stdio.h>
#include <math.h>

//2-norm of a vector
double two_norm(int p, double N[p]) {
	double sum = 0.0;
	for (int i = 0; i < p; i++) {
		sum += N[i] * N[i];
	}
	return sqrt(sum);
}

double vectorMulti(int p,double v1[p], double v2[p]){
	double result = 0;
	for (int i = 0; i < p; i++) {
		result += v1[i] * v2[i];
	}
	return result;
}

int main(){
int Ni,Nj,N;
Ni=Nj=N=128;
int p=pow(N,2);
double T[p], Tnew[p], Aw[p], Ae[p],An[p], As[p], Ap[p], b[p], dold[p], dnew[p];

//Initialization
int z=pow(N,2)-N;

for(int i=0;i<pow(N,2);i++){
	T[i]=0;
	Tnew[i]=0;
}

////west co-efficient
for(int i=0;i<pow(N,2);i++){
	Aw[i]=-1;
	for(int j=0;j<=z;j=j+N){
		Aw[j]=0;
	}
}



//east co-efficient
for(int i=0;i<pow(N,2);i++){
	Ae[i]=-1;
	for(int j=N-1;j<p;j=j+N){
		Ae[j]=0;
	}
}

//north co-efficient
for(int i=0;i<pow(N,2);i++){ 
	An[i]=-1;
	for(int j=0;j<N;j++){
		An[j]=0;
	}
}

//south co-efficient
for(int i=0;i<pow(N,2);i++){ 
	As[i]=-1;
	for(int j=z;j<p;j++){
		As[j]=0;
	}
}

//p co-efficient
for(int i=0;i<pow(N,2);i++){
	Ap[i]=4;
	Ap[0]=Ap[N-1]=Ap[z]=Ap[p-1]=6;
	for(int j=1;j<N-1;j++){
		Ap[j]=5;
	}
	for(int j=N;j<=z-N;j=j+N){
		Ap[j]=5;
	}
	for(int j=2*N-1;j<p-N;j=j+N){
		Ap[j]=5;
	}
	for(int j=z+1;j<p-1;j++){
		Ap[j]=5;	
	}
}

//b-vector
for(int i=0;i<pow(N,2);i++){
	b[i]=0;
	for(int j=0;j<N;j++){
		b[j]=2;
	}
}

//printf("Aw\tAe\tAn\tAs\tAp\tRHS\n");

/*for(int i=0;i<pow(N,2);i++){
	printf("%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",Aw[i],Ae[i],An[i],As[i],Ap[i],b[i]);
}*/

double rold[p], rnew[p], alpha, beta;
for(int i=0;i<pow(N,2);i++){
	//rold[i]=b[i]-(Ap[i]*T[i] + Aw[i]*T[i-1] + Ae[i]*T[i+1] + An[i]*T[i-N] + As[i]*T[i+N]);
	
	rold[i]=b[i]-Ap[i]*T[i];
	if(Aw[i]!=0) {rold[i] = rold[i] - Aw[i]*T[i-1]; }
	if(Ae[i]!=0) {rold[i] = rold[i] - Ae[i]*T[i+1]; }
	if(An[i]!=0) {rold[i] = rold[i] - An[i]*T[i-N]; }
	if(As[i]!=0) {rold[i] = rold[i] - As[i]*T[i+N]; }
	
}

for(int i=0;i<pow(N,2);i++){
	dold[i]=rold[i]/Ap[i];
}

double W,Y,error;
double I[p],delta1[p], delta2[p];
error=0;
int iter = 0;
FILE *fi= fopen("iterationvsresidue.dat","w");
fprintf(fi, "Variables = \"iterations\",\"residue\"\n");
do{
	
	
	for(int i=0;i<pow(N,2);i++){
		delta1[i]=rold[i]/Ap[i];
	}

	W=vectorMulti(p,delta1,rold);
	
	for(int i=0;i<pow(N,2);i++){
		I[i]= Ap[i]*dold[i];
		if(Aw[i]!=0) {I[i] = I[i] + Aw[i]*dold[i-1]; }
		if(Ae[i]!=0) {I[i] = I[i] + Ae[i]*dold[i+1]; }
		if(An[i]!=0) {I[i] = I[i] + An[i]*dold[i-N]; }
		if(As[i]!=0) {I[i] = I[i] + As[i]*dold[i+N]; }	
	}
	Y=vectorMulti(p,dold,I);
	alpha=W/Y;

	for(int i=0;i<pow(N,2);i++){
		Tnew[i]=Tnew[i]+alpha*dold[i];
	}
	/*for(int i=0;i<pow(N,2);i++){
	printf("Tnew[%d]=%lf\n",i,Tnew[i]);
	}*/

	for(int i=0;i<pow(N,2);i++){
		rnew[i]=rold[i]-alpha*I[i];
	}
	
	for(int i=0;i<pow(N,2);i++){
		delta2[i]=rnew[i]/Ap[i];
	}

	beta=vectorMulti(p,rnew,delta2)/vectorMulti(p,rold,delta1);

	for(int i=0;i<pow(N,2);i++){
		dnew[i]=delta2[i]+beta*dold[i];
	}

	
	error=sqrt(vectorMulti(p, rnew, rnew));
	
	//update
	for(int i=0;i<pow(N,2);i++){
		dold[i]=dnew[i];
		rold[i]=rnew[i];
	}
	iter += 1;
	//if(iter >=1)	break;
	
	//if(iter%10==0)	printf("\t Iter: %d, Error: %.4f\n", iter, error);
	fprintf(fi,"%d\t%.8f\n", iter, error);
	//printf("error=%lf\n",error);
}while(error > pow(10,-6));

printf("\t Iter: %d\n", iter);

FILE *f1= fopen("T_lexi","w");

for(int i=0; i<=N; i++)	fprintf(f1, "%.6lf\t", 1.0);	// Top Boundary

for(int i=0; i<pow(N,2)-N-1; i++){
	if(i%N==0)	fprintf(f1, "\n%.6lf\t", 0.0);	// Left boundary
	fprintf(f1, "%.6lf\t" , 0.25*(Tnew[i] + Tnew[i+1] + Tnew[i+N] + Tnew[i+N+1]));
	if(i%N==N-2){
		fprintf(f1, "%.6lf", 0.0);	// Right Boundary
		i += 1;
	}
}

fprintf(f1, "\n");
for(int i=0; i<=N; i++)	fprintf(f1, "%.6lf\t", 0.0);	// Bottom Boundary


FILE *f2= fopen("Exact_T_lexi","w");

double dx = (double)1/N;
double dy = dx;

double temp;

for(int i=0; i<=N; i++){
	for(int j=0; j<=N; j++){
		temp = 0.0;
		for(int q=1; q<200; q++)	temp  += (pow(-1, q+1) + 1)/q*sin(q*M_PI*j*dx)*sinh(q*M_PI*(1-i*dy))/sinh(q*M_PI);
		temp = 2.0*temp/M_PI;
		fprintf(f2, "%.6lf\t", temp);
	}
	fprintf(f2, "\n");
}
fclose(f1);
fclose(f2);

double delx=1.0/N;
double dely=1.0/N;


double A[129][129], Aexact[129][129];
FILE *f8= fopen("T_lexi","r");
FILE *f9= fopen("Exact_T_lexi","r");
FILE *f3= fopen("Temperature results.dat","w");
FILE *f4= fopen("Exact Temperature results.dat","w");

for(int j=128;j>=0;j--){
	for(int i=128;i>=0;i--){
		fscanf(f8,"%lf",&A[i][j]);
		fscanf(f9,"%lf",&Aexact[i][j]);	
	}
}

fprintf(f3, "Variables = \"X\",\"Y\", \"Calculated (T)\"\n");
fprintf(f4, "Variables = \"X\",\"Y\", \"Exact (T)\"\n");
for(int i=0;i<129;i++){
	for(int j=0;j<129;j++){
		double o=(double)i*(double)delx;
        	double u=(double)j*dely;
		
		fprintf(f3,"%lf\t%lf\t%lf\n",o,u,A[i][j]);
		fprintf(f4,"%lf\t%lf\t%lf\n",o,u,Aexact[i][j]);
	}
}
fclose(f1);
fclose(f2);
fclose(f3);
fclose(f4);

FILE *f5= fopen("T at centreline(x=0.5).dat","w");
        fprintf(f5, "Variables = \"T at centerline(x=0.5)\",\"Y\"\n");
		for(int j=0;j<=N;j++){
			double q=j*dely;
        		fprintf(f5,"%lf\t%lf\n",A[64][j],q);
        	}
        fclose(f5);

FILE *f6= fopen("T at centreline(y=0.5).dat","w");
fprintf(f6, "Variables = \"X\",\"T at centerline(y=0.5)\"\n");
        for(int i=0;i<=N;i++){ //for v-velocity
			double q=i*delx;
        		fprintf(f6,"%lf\t%lf\n",q,A[i][64]);
        }
fclose(f6);
FILE *fe= fopen("T(exact) at centreline(x=0.5).dat","w");
fprintf(fe, "Variables = \"T at centerline(x=0.5)\",\"Y\"\n");
		for(int j=0;j<=N;j++){
			double q=j*dely;
        		fprintf(fe,"%lf\t%lf\n",Aexact[64][j],q);
        	}
        fclose(fe);
FILE *ff= fopen("T(exact) at centreline(y=0.5).dat","w");
fprintf(ff, "Variables = \"X\",\"T at centerline(y=0.5)\"\n");
        for(int i=0;i<=N;i++){ //for v-velocity
			double q=i*delx;
        		fprintf(ff,"%lf\t%lf\n",q,Aexact[i][64]);
        }
fclose(ff);	
	
}
