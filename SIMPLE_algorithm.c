#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Function to find the maximum between two numbers
double maxtwo(double a, double b) {
    return (a > b) ? a : b;
}

// Function to find the maximum among three numbers
double maxthree(double x, double y, double z) {
    return maxtwo(maxtwo(x, y), z);
}

int main(){
	int i,j,g,k,w,Re,m,n,L1,L2;
	k=0;w=0,g=0;
	int iteration =0;
	/*printf("Enter the value of Grid Size\n");
	printf("for m = ");
	scanf("%d",&m);
	printf("\nfor n = ");
	scanf("%d",&n);
	printf("\nEnter The value of Reynold's Number\n");
	scanf("%d",&Re);
	printf("\nEnter the length along x-direction\n");
	scanf("%d",&L1);
	printf("\nEnter the length along y-direction\n");
	scanf("%d",&L2);
	double delx=L1/(m);
	double dely=L2/(n);*/
	n=129;m=129;double delx=1.0/n;double dely=1.0/m; Re=100;
	double Ae,Aw,An,As,De,Dw,Dn,Ds,Fe,Fw,Fn,Fs,F[m+2][n+2],a_e,a_w,a_n,a_s,a_[m+2][n+2];
	double D_e,D_w,D_n,D_s,F_e,F_w,F_n,F_s,F_[m+5][n+5],ae,aw,an,as,a[m+5][n+5];
	double alpha_p=0.5, alpha_u=0.5, alpha_v=0.5;
	
	double oldustar[m+2][n+2], newustar[m+2][n+2], oldvstar[m+2][n+2], newvstar[m+2][n+2], pstar[m+2][n+2], oldpc[m+2][n+2], newpc[m+2][n+2], 	uc[m+2][n+2], vc[m+2][n+2], u[m+2][n+2], uold[m+2][n+2], v[m+2][n+1], vold[m+2][n+2], d_[m+2][n+2], d[m+2][n+2], Ac[m+2][n+2], bc[m+2][n+2], p[m+2][n+2], S_u, Su, S_p, Sp, Ac_e, Ac_w, Ac_n, Ac_s;
	
	double col_u[m+2][n+2], col_v[m+2][n+2], omega[m+2][n+2];
	double error_u=0;
	double error_v=0;
	double error_p=0;
	S_u=0; S_p=0; Su=0; Sp=0;
	
	//Boundary Condition
	for(j=0;j<n;j++){
		newustar[0][j]=0; //left wall
		oldustar[0][j]=0;
		newustar[m][j]=0; //right wall
		oldustar[m][j]=0;
	}
	
	for(i=0;i<m;i++){
		newvstar[i][0]=0; //bottom wall
		oldvstar[i][0]=0;
		newvstar[i][n]=0; //top wall
		oldvstar[i][n]=0;
	}
	
	//Initialization
	for(i=0;i<m;i++){  //for pressure
		for(j=0;j<n;j++){
			pstar[i][j]=0;
			oldpc[i][j]=0;
			newpc[i][j]=0;
			p[i][j]=0;
			omega[i][j]=0;
		}
	}
	
	for(i=1;i<m;i++){ //for u-velocity
		for (j=0;j<n;j++){
			oldustar[i][j]=0;
			newustar[i][j]=0;
			uc[i][j]=0;
			uold[i][j]=0;
			u[i][j]=0;
		}
	}
	
	for(i=0;i<m;i++){ //for v-velocity
		for (j=1;j<n;j++){
			oldvstar[i][j]=0;
			newvstar[i][j]=0;
			vc[i][j]=0;
			vold[i][j]=0;
			v[i][j]=0;
			}
	}
				
	
	
	
	double error1=0;
	double error2=0;
	double error3=0;
	
	do{
		Ae=Aw=delx;
		An=As=dely;
		
		De=Ae/(delx*Re);
		Dn=An/(dely*Re);
		Dw=Aw/(delx*Re);
		Ds=As/(dely*Re);
		
	// u-momentum equation
		for(i=1;i<m;i++){
				for (j=0;j<n;j++){
						Fe=Ae*0.5*(oldustar[i+1][j] + oldustar[i][j]);
						Fw=Aw*0.5*(oldustar[i][j] + oldustar[i-1][j]);
						Fn=An*0.5*(oldvstar[i][j+1] + oldvstar[i-1][j+1]);
						Fs=As*0.5*(oldvstar[i][j] + oldvstar[i-1][j]);
		
						a_w=maxthree(Fw,Dw+(0.5*Fw),0);
						a_e=maxthree(-Fe,De-(0.5*Fe),0);
						a_s=maxthree(Fs,Ds+(0.5*Fs),0);
						a_n=maxthree(-Fn,Dn-(0.5*Fn),0);
						
						S_p = 0;
						S_u=0;
						
						if(j==0){
						Fs=0.0;//As*0.5*(oldvstar[i][j] + oldvstar[i-1][j]);
						a_s=maxthree(Fs,Ds+(0.5*Fs),0);
						S_p=-a_s;}
						
						if(j==n - 1){
						Fn=0;//An*0.5*(oldvstar[i][j+1] + oldvstar[i-1][j+1]);
						a_n=maxthree(-Fn,Dn-(0.5*Fn),0);
						S_p=-a_n;
						S_u=2*a_n;}
						
						a_[i][j] = a_w + a_e + a_n + a_s + (Fe - Fw + Fn - Fs) - S_p;
						d_[i][j]=dely/a_[i][j];
						
						newustar[i][j] = a_w*newustar[i-1][j]+a_e*newustar[i+1][j]+dely*(pstar[i-1][j] - pstar[i][j]) + S_u;
						
						if (j!=0) {newustar[i][j] = newustar[i][j] + a_s*newustar[i][j-1];}
						if (j!=n-1) {newustar[i][j] = newustar[i][j] + a_n*newustar[i][j+1];}
						
						newustar[i][j] = newustar[i][j]/a_[i][j];
				}
			}
			
		
		FILE *file1= fopen("ustar-velocity.dat","w");
        
		//printf("iteration=%d\n",k);
		for(i=1;i<m;i++){
			for (j=0;j<n;j++){
				fprintf(file1,"\nu[%d][%d]=%lf\n",i,j,newustar[i][j]);
			}
		}
		fclose(file1);
		
	//v-momentum equation
		for(i=0;i<m;i++){
				for(j=1;j<n;j++){
						
						F_e=Ae*0.5*(oldustar[i+1][j] + oldustar[i+1][j-1]);
						F_w=Aw*0.5*(oldustar[i][j] + oldustar[i][j-1]);
						F_n=An*0.5*(oldvstar[i][j+1] + oldvstar[i][j]);
						F_s=As*0.5*(oldvstar[i][j] + oldvstar[i][j-1]);
		
						aw=maxthree(F_w,Dw+(0.5*F_w),0);
						ae=maxthree(-F_e,De-(0.5*F_e),0);
						as=maxthree(F_s,Ds+(0.5*F_s),0);
						an=maxthree(-F_n,Dn-(0.5*F_n),0);
						
						Sp=0;
						Su=0;
						
						if(i==0){
						F_w= 0;
						aw=maxthree(F_w,Dw+(0.5*F_w),0);
						Sp=-aw;
						}
						
						if(i==m-1){
						F_e=0;
						ae=maxthree(-F_e,De-(0.5*F_e),0);
						Sp=-ae;
						}
											
						a[i][j] = aw + ae + an + as + (F_e - F_w + F_n - F_s) -Sp;
						d[i][j]=delx/a[i][j];
						
						newvstar[i][j]= an*newvstar[i][j+1] + as*newvstar[i][j-1] + delx*(pstar[i][j-1] - pstar[i][j]);
						if (i!=0) {newvstar[i][j] = newvstar[i][j] + aw*newvstar[i-1][j];}
						if (i!=m-1) {newvstar[i][j] = newvstar[i][j] + ae*newvstar[i+1][j];}
						
						newvstar[i][j]=newvstar[i][j]/a[i][j];
				}
			}
        		
		
		FILE *file2= fopen("vstar-velocity.dat","w");
		for(i=0;i<m;i++){
			for (j=1;j<n;j++){
				fprintf(file2,"\n%lf\n",newvstar[i][j]);
			}
		}
		fclose(file2);
		
	//Solve Pressure Correction Equation
				
		for(i=0;i<m;i++){
			for(j=0;j<n;j++){
				newpc[i][j]=0;
			}
		}
		
		for(i=0;i<m;i++){
			for(j=0;j<n;j++){
						
				Ac_e=Ae*d_[i+1][j];
				Ac_w=Aw*d_[i][j];
				Ac_n=An*d[i][j+1];
				Ac_s=As*d[i][j];
				
				bc[i][j]=(Ae*newustar[i][j])-(Ae*newustar[i+1][j]) + (An*newvstar[i][j])-(An*newvstar[i][j+1]);
				
				if(i==m-1){
					Ac_e=0;
				}
				if(i==0){
					Ac_w=0;
				}
				if(j==n-1){
					Ac_n=0;
				}
				if(j==0){
					Ac_s=0;
				}
				
				Ac[i][j] = Ac_e+Ac_w+Ac_n+Ac_s;
				
				newpc[i][j] = bc[i][j];
				if (Ac_e!=0) {newpc[i][j] = newpc[i][j] + Ac_e*newpc[i+1][j]; }
				if (Ac_w!=0) {newpc[i][j] = newpc[i][j] + Ac_w*newpc[i-1][j]; }
				if (Ac_n!=0) {newpc[i][j] = newpc[i][j] + Ac_n*newpc[i][j+1]; }
				if (Ac_s!=0) {newpc[i][j] = newpc[i][j] + Ac_s*newpc[i][j-1]; }
				newpc[i][j] = newpc[i][j] / Ac[i][j];
			}
		}
		
		
		FILE *file3= fopen("pcorrection-velocity.dat","w");
		for(i=0;i<m;i++){
			for (j=0;j<n;j++){
				fprintf(file3,"\npc[%d][%d]=%lf\n",i,j,newpc[i][j]);
			}
		}
		fclose(file3);
		//printf("iteration=%d\n",w);
		
	///////Improved pressure
		for(i=0;i<m;i++){
			for(j=0;j<n;j++){
				p[i][j]=pstar[i][j] + (alpha_p*newpc[i][j]);
			}
		}
		
		for(i=0;i<m;i++){
			for(j=0;j<n;j++){
				error_p = error_p + (p[i][j] - pstar[i][j])*(p[i][j] - pstar[i][j]);
			}
		}
		error_p = sqrt(error_p/ ((m-2) * (n-2)));
		
	//Solve u-correction equation
		for(i=1;i<m;i++){
			for(j=0;j<n;j++){
				//d_[i][j]=dely/a_[i][j];
				uc[i][j]=newustar[i][j] + d_[i][j]*(newpc[i-1][j]-newpc[i][j]);
			}
		}
		
	//Solve v-correction equation
		for(i=0;i<m;i++){
			for(j=1;j<n;j++){
				//d[i][j]=delx/a[i][j];
				vc[i][j]=newvstar[i][j] + d[i][j]*(newpc[i][j-1]-newpc[i][j]);
			}
		}
		
		
		
	//Improved u-velocity
		for(i=1;i<m;i++){
			for(j=0;j<n;j++){
				u[i][j]=(alpha_u*uc[i][j]) + ((1-alpha_u)*uold[i][j]);
			}
		}
		
        	for(i=1;i<m;i++){
			for(j=0;j<n;j++){
				error_u = error_u + (u[i][j] - uold[i][j])*(u[i][j] - uold[i][j]);
			}
		}
		error_u = sqrt(error_u/ ((m-2) * (n-2)));
		
	//Improved v-velocity
		for(i=0;i<m;i++){
			for(j=1;j<n;j++){
				v[i][j]=(alpha_v*vc[i][j]) + ((1-alpha_v)*vold[i][j]);
			}
		}
		
        	for(i=1;i<m;i++){
			for(int j=0;j<n;j++){
				error_v = error_v + (v[i][j] - vold[i][j])*(v[i][j] - vold[i][j]);
			}
		}
		error_v = sqrt(error_v/ ((m-2) * (n-2)));
        	
        //update
        	for(i=1;i<m;i++){ //for u-velocity
			for(j=0;j<n;j++){
				uold[i][j]=u[i][j];
				oldustar[i][j]=u[i][j];
				newustar[i][j]=u[i][j];
			}
		}
		
		for(i=0;i<m;i++){ //for v-velocity
			for(j=1;j<n;j++){
				vold[i][j]=v[i][j];
				oldvstar[i][j]=v[i][j];
				newvstar[i][j]=v[i][j];
			}
		}
		
		for(i=0;i<m;i++){ //for v-velocity
			for(j=0;j<n;j++){
				pstar[i][j]=p[i][j];
			}
		}
		iteration++;
		printf("error_u=%lf and error_v=%lf\n",error_u,error_v);
        printf("iteration_f=%d\n",iteration);
        }while(error_u>pow(10,-4) || error_v>pow(10,-4)); //&& error_p>pow(10,-5));
        printf("error_u=%lf and error_v=%lf\n",error_u,error_v);
        printf("iteration_f=%d\n",iteration);
        
        //u-centreline vs y
        FILE *file4= fopen("u-velocitycentreline.dat","w");
        fprintf(file4, "Variables = \"U at centerline\",\"Y\"\n");
        //for(i=1;i<m;i++){ //for u-velocity
		for(j=0;j<n;j++){
			double q=j*dely;
        		fprintf(file4,"%lf\t%lf\n",u[65][j],q);
        	}
        //}
        fclose(file4);
        
        //v-centreline vs x
        FILE *file5= fopen("v-velocity.dat","w");
        fprintf(file4, "Variables = \"X\",\"V at centerline\"\n");
        for(i=0;i<m;i++){ //for v-velocity
		//for(j=1;j<n;j++){
			double q=i*delx;
        		fprintf(file5,"%lf\t%lf\n",q,v[i][65]);
        	//}
        }
        fclose(file5);
        
        for(j=0;j<n-1;j++){
        	u[0][j]=0;
        	u[m][j]=0;
        }
        
        for(i=0;i<m-1;i++){
        	v[i][0]=0;
        	v[i][n]=0;
        }
        
        //After Convergence, mapping the Staggered variables to collocated variables
        for(i=0;i<m;i++){
        	for(j=0;j<n;j++){
        		col_u[i][j]=0.5*(u[i][j]+u[i+1][j]);
        		col_v[i][j]=0.5*(v[i][j]+v[i][j+1]);
        	}
        }
        
        FILE *file6= fopen("velocity-magnitude.dat","w");
        fprintf(file6, "Variables = \"X\",\"Y\", \"V\"\n");
        for(i=0;i<m;i++){
        	for(j=0;j<n;j++){
        		double A=i*delx;
        		double B=j*dely;
        		
        		fprintf(file6,"%lf\t%lf\t%lf\n",A,B,pow((col_u[i][j]*col_u[i][j]+col_v[i][j]*col_v[i][j]),0.5));
        	}
        }
        fclose(file6);
        
        //Calculating Omega
        omega[0][0]=(0.25/delx)*(v[1][0]+v[1][1]-u[0][1]-u[1][1]);
        
        omega[m-1][0]=(0.25/delx)*(-v[m-2][0]-v[m-2][1]-u[0][1]-u[1][1]);
        
        omega[0][n-1]=(0.25/delx)*(v[1][n-1]+v[1][n]-2+u[0][n-2]+u[1][n-2]);
        
        omega[m-1][n-1]=(0.25/delx)*(-v[m-2][n-1]-v[m-2][n]-2+u[m-1][n-2]+u[m][n-2]);
        
        for(i=1;i<m-1;i++){
        	omega[i][0]=(0.25/delx)*(v[i+1][0]+v[i+1][1]-v[i-1][0]-v[i-1][1]-u[i][1]-u[i+1][1]);
        }
        
        for(i=1;i<m-1;i++){
        	omega[i][n-1]=(0.25/delx)*(v[i+1][n-1]+v[i+1][n]-v[i-1][n-1]-v[i-1][n]-2+u[i][n-2]+u[i+1][n-2]);
        }
        
        for(j=1;j<n-1;j++){
        	omega[0][j]=(0.25/delx)*(v[1][j]+v[1][j+1]-u[0][j+1]-u[1][j+1]+u[0][j-1]+u[1][j-1]);
        }
        
        for(j=1;j<n-1;j++){
        	omega[m-1][j]=(0.25/delx)*(-v[m-2][j]-v[m-2][j+1]-u[m-1][j+1]-u[m][j+1]+u[m-1][j-1]+u[m][j-1]);
        }
        
        for(i=1;i<m-1;i++){
        	for(j=1;j<n-1;j++){
        		omega[i][j]=(0.25/delx)*(v[i+1][j]+v[i+1][j+1]-v[i-1][j]-v[i-1][j+1]-u[i][j+1]-u[i+1][j+1]+u[i][j-1]+u[i+1][j-1]);
        	}
        }
        
        FILE *file7= fopen("vorticity.dat","w");
        fprintf(file7, "Variables = \"X\",\"Y\", \"VORTICITY\"\n");
        for(i=0;i<m;i++){
        	for(j=0;j<n;j++){
        		double A=i*delx;
        		double B=j*dely;
        		fprintf(file7,"%lf\t%lf\t%lf\n",A,B,omega[i][j]);
        	}
        }
        fclose(file7);
}
