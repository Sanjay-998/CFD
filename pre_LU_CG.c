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

void LU(int p,double Aw[p], double Ae[p], double An[p], double Ap[p], double As[p], double H[p], double rold[p]){
	double Lp[p], Lw[p], Ls[p], Un[p], Ue[p], C[p];
	int N=128;
	
	//INItIALIZING
	for(int i=0;i<pow(N,2);i++){
		Ls[i]=0;
		Lp[i]=0;
		Lw[i]=0;
		Un[i]=0;
		Ue[i]=0;
		C[i]=0;
		H[i]=0;
	}

	for(int i=0;i<pow(N,2);i++){
		Lw[i]=0;
		for(int j=N; j<pow(N,2);j++){
			Lw[i]=Aw[i];
		}
	}

	//Ls
	Ls[0]=0;
	for(int i=1;i<pow(N,2);i++){
		Ls[i]=As[i];
	}
	
	//Lp
	Lp[0]=Ap[0];
	for(int i=1;i<N;i++){
		Lp[i]=Ap[i]-Ls[i]*Un[i-1];
	}
	for(int i=N;i<pow(N,2);i++){
		Lp[i]=Ap[i]-Ls[i]*Un[i-1]-Lw[i]*Ue[i-N];
	}
	
	//Un
	int Z=pow(N,2)-1;
	Un[Z]=0;
	for(int i=0; i<pow(N,2)-1;i++){
		Un[i]=An[i]/Lp[i];
	}
	
	//Ue
	for(int i=0; i<pow(N,2)-N;i++){
		Ue[i]=Ae[i]/Lp[i];
	}
	for(int i=pow(N,2)-N; i<pow(N,2); i++){
		Ue[i]=0;
	}
	
	//Forward Substitution
	for(int i=0;i<pow(N,2);i++){
		C[i]=rold[i];
		if(Ls[i]!=0) {C[i] = C[i] - Ls[i]*C[i+N]; }
		if(Lw[i]!=0) {C[i] = C[i] - Lw[i]*C[i-1]; }
		C[i]=C[i]/Lp[i];
		
	}
	
	//Backward Substitution
	for(int i=pow(N,2)-1;i>0;i--){
		H[i]=C[i];
		if(Ue[i]!=0) {H[i] = H[i] - Ue[i]*H[i+1]; }
		if(Un[i]!=0) {H[i] = H[i] - Un[i]*H[i-N]; }
	}
}

int main(){
int Ni,Nj,N;
Ni=Nj=N=128;
int p=pow(N,2);
double T[p], Tnew[p], Aw[p], Ae[p],An[p], As[p], Ap[p], b[p], dold[p], dnew[p], Lp[p], Lw[p], Ls[p], Un[p], Ue[p], H[p], C[p];

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

double rold[p], rnew[p], alpha, beta;
for(int i=0;i<pow(N,2);i++){
	//rold[i]=b[i]-(Ap[i]*T[i] + Aw[i]*T[i-1] + Ae[i]*T[i+1] + An[i]*T[i-N] + As[i]*T[i+N]);
	
	rold[i]=b[i]-Ap[i]*T[i];
	if(Aw[i]!=0) {rold[i] = rold[i] - Aw[i]*T[i-1]; }
	if(Ae[i]!=0) {rold[i] = rold[i] - Ae[i]*T[i+1]; }
	if(An[i]!=0) {rold[i] = rold[i] - An[i]*T[i-N]; }
	if(As[i]!=0) {rold[i] = rold[i] - As[i]*T[i+N]; }
	
}

LU(p,Aw,Ae,An,Ap,As,H,rold);

for(int i=0;i<pow(N,2);i++){
	dold[i]=H[i]*rold[i];
}

double W,Y,error;
double I[p],delta1[p], delta2[p];;
error=0;
int iter = 0;
FILE *fi= fopen("iterationvsresidue.dat","w");
fprintf(fi, "Variables = \"iterations\",\"residue\"\n");
do{
	for(int i=0;i<pow(N,2);i++){
		I[i]= Ap[i]*dold[i];
		if(Aw[i]!=0) {I[i] = I[i] + Aw[i]*dold[i-1]; }
		if(Ae[i]!=0) {I[i] = I[i] + Ae[i]*dold[i+1]; }
		if(An[i]!=0) {I[i] = I[i] + An[i]*dold[i-N]; }
		if(As[i]!=0) {I[i] = I[i] + As[i]*dold[i+N]; }	
	}
	
	for(int i=0;i<pow(N,2);i++){
		delta1[i]=rold[i]*H[i];
	}

	W=vectorMulti(p,rold,delta1);
	Y=vectorMulti(p,dold,I);
	alpha=W/Y;
	//printf("\n\t Alpha: %.4f\n", alpha);

	for(int i=0;i<pow(N,2);i++){
		Tnew[i]=Tnew[i]+alpha*dold[i];
	}

	for(int i=0;i<pow(N,2);i++){
		rnew[i]=rold[i]-alpha*I[i];
	}
	
	for(int i=0;i<pow(N,2);i++){
		delta2[i]=rnew[i]*H[i];
	}

	beta=vectorMulti(p,rnew,rnew)/vectorMulti(p,rold,rold);

	for(int i=0;i<pow(N,2);i++){
		dnew[i]=rnew[i]+beta*dold[i];
	}

	
	error=sqrt(vectorMulti(p, rnew, rnew));
	
	//update
	for(int i=0;i<pow(N,2);i++){
		dold[i]=dnew[i];
		rold[i]=rnew[i];
	}
	iter += 1;
	fprintf(fi,"%d\t%.8f\n", iter, error);
}while(error > pow(10,-6));
fclose(fi);

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
}fclose(f1);
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
	//fscanf(f1,"\n");
	//fscanf(f2,"\n");
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

FILE *fe= fopen("T(exact) at centreline(x=0.5).dat","w");
fprintf(fe, "Variables = \"T at centerline(x=0.5)\",\"Y\"\n");
		for(int j=0;j<=N;j++){
			double q=j*dely;
        		fprintf(fe,"%lf\t%lf\n",Aexact[64][j],q);
        	}
        fclose(fe);

FILE *f6= fopen("T at centreline(y=0.5).dat","w");
fprintf(f6, "Variables = \"X\",\"T at centerline(y=0.5)\"\n");
        for(int i=0;i<=N;i++){ //for v-velocity
			double q=i*delx;
        		fprintf(f6,"%lf\t%lf\n",q,A[i][64]);
        }
fclose(f6);

FILE *ff= fopen("T(exact) at centreline(y=0.5).dat","w");
fprintf(ff, "Variables = \"X\",\"T at centerline(y=0.5)\"\n");
        for(int i=0;i<=N;i++){ //for v-velocity
			double q=i*delx;
        		fprintf(ff,"%lf\t%lf\n",q,Aexact[i][64]);
        }
fclose(ff);	
}
