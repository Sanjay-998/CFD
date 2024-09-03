//For  Re=100

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
int main()
{

	int i,j ,m,n,midpts;
	double delx,dely,x,y;
	m=100;
	n=100;
	midpts=(m-2)*(n-2);
	delx=1.0/(m);
	dely=1.0/(n);
	double beta=(delx/dely);
	double Re=100.0;
	double psi[m][n],omega[m][n],u[m][n],v[m][n];
	double psi_prev[m][n],omega_prev[m][n];
	double error_psi=0.0,error_omega=0.0;
	
	int iteration=0;
	for(j=0;j<n;j++)
	{
		for(i=0;i<m;i++)
		{
			if(j==0)//bb
			{
				u[i][j]=0.0;
				v[i][j]=0.0;
				psi[i][j]=0.0;
			}
			else if(j==(n-1))
			{
				u[i][j]=1.0;
				v[i][j]=0.0;
				psi[i][j]=0.0;
					
			}
			else if(i==0)
			{
				u[i][j]=0.0;
				v[i][j]=0.0;
				psi[i][j]=0.0;
					
			}
			else if(i==(m-1))
			{
				u[i][j]=0.0;
				v[i][j]=0.0;
				psi[i][j]=0.0;
					
			}
				else 
			{
				u[i][j]=0.0;
				v[i][j]=0.0;
				psi[i][j]=0.0;
					
			}
		}
	}
	for(j=0;j<n;j++)
	{
		for(i=0;i<m;i++)
		{
			if(j==(n-1))
			{
				omega[i][j]=((2.0/pow(dely,2))*(psi[i][j]-psi[i][j-1]))-(2.0/dely);
			}
			else if(j==0)
			{
				omega[i][j]=(2.0/pow(dely,2))*(psi[i][j]-psi[i][j+1]);
			}
			else if(i==0)
			{
			    omega[i][j]=(2.0/pow(dely,2))*(psi[i][j]-psi[i+1][j]);	
			}
			else if(i==(m-1))
			{
			 omega[i][j]=(2.0/pow(dely,2))*(psi[i][j]-psi[i-1][j]);	
			}
			else
			{
				omega[i][j]=0.0;
			}
		}
	}
	do
	{
	   for(j=0;j<n;j++)
    	{
	    	for(i=0;i<m;i++)
			{
				psi_prev[i][j]=psi[i][j];
				omega_prev[i][j]=omega[i][j];
			}
    	}
       for(j=1;j<(n-1);j++)
    	{
	    	for(i=1;i<(m-1);i++)
			{
			psi[i][j]=(0.5/(1.0+pow(beta,2)))*(psi[i+1][j]+psi[i-1][j]+(pow(beta,2)*(psi[i][j+1]+psi[i][j-1]))+(pow(delx,2)*omega[i][j]));	
	        }
	    }
	   for(j=1;j<(n-1);j++)
    	{
	    	for(i=1;i<(m-1);i++)
			{
				omega[i][j]=(0.5/(1.0+pow(beta,2)))*((1.0-((psi[i][j+1]-psi[i][j-1])*((beta*Re)/4.0)))*omega[i+1][j]
				+(1.0+((psi[i][j+1]-psi[i][j-1])*((beta*Re)/4.0)))*omega[i-1][j]
				+((1.0+((psi[i+1][j]-psi[i-1][j])*(Re/(4.0*beta))))*(pow(beta,2)*omega[i][j+1]))
				+((1.0-((psi[i+1][j]-psi[i-1][j])*(Re/(4.0*beta))))*(pow(beta,2)*omega[i][j-1]))); 
     	    } 
	    }
	    	for(j=0;j<n;j++)
        	{
		        for(i=0;i<m;i++)
	        	{
		    	  if(j==(n-1))
		       	{
			    	omega[i][j]=((2.0/pow(dely,2))*(psi[i][j]-psi[i][j-1]))-(2.0/dely);
		    	}
		    	else if(j==0)
		    	{
			    	omega[i][j]=(2.0/pow(dely,2))*(psi[i][j]-psi[i][j+1]);
		    	}
		    	else if(i==0)
		    	{
			        omega[i][j]=(2.0/pow(dely,2))*(psi[i][j]-psi[i+1][j]);	
		    	}
		    	else if(i==(m-1))
		    	{
			        omega[i][j]=(2.0/pow(dely,2))*(psi[i][j]-psi[i-1][j]);	
		    	}
	        	}
        	}
      error_psi=0.0;
	  error_omega=0.0;
	  
	  	for(j=1;j<(n-1);j++)
        	{
		        for(i=1;i<(m-1);i++)
	        	{
	        		error_psi=error_psi+pow((psi[i][j]-psi_prev[i][j]),2.0);
	        		error_omega=error_omega+pow((omega[i][j]-omega_prev[i][j]),2.0);
	        	}
	        }
	    	error_psi=sqrt(error_psi/midpts);
	    	error_omega=sqrt(error_omega/midpts);
	    	printf("iteration=%d\t",iteration);
	    	printf("error_psi=%.10lf\terror_omega=%.10lf\n",error_psi,error_omega);
	    	iteration++;
	    	for(j=0;j<(n-1);j++)
        	{
		        for(i=0;i<(m-1);i++)
	        	{
	              u[i][j]=(0.5/dely)*(psi[i][j+1]-psi[i][j-1]);
	              v[i][j]=(-0.5/delx)*(psi[i+1][j]-psi[i-1][j]);
	            }
	        }
	    }
	    while(error_psi>1.0e-6||error_omega>1.0e-6);
	    
			
	        FILE*file;
	        file=fopen("(Re100)outputlid.dat","w");
          	
          	FILE*file1;
          	file1=fopen("(Re100)Streamlines.dat","w");
          	
          	FILE*file2;
          	file2=fopen("(Re100)VelocityVectors.dat","w");
          	
          	FILE*file3;
          	file3=fopen("(Re100)Velocity(u)centre.dat","w");
          	
          	FILE*file4;
          	file4=fopen("(Re100)Velocity(v)centre.dat","w");
          	
          	FILE*file5;
          	file5=fopen("(Re100)Vorticitycontour.dat","w");
          	
	        fprintf(file,"Zone I=%d, J=%d\n",m,n);
	        
	    	fprintf(file1, "Variables = \"X\",\"Y\", \"PSI\"\n");
	        
	        fprintf(file2, "Variables = \"X\",\"Y\", \"U\",\"V\"\n");
	        
	   	fprintf(file3, "Variables = \"U at centerline\",\"Y\"\n");
	        
	        fprintf(file4, "Variables = \"X\",\"V at centerline\"\n");
	        
	        fprintf(file5, "Variables = \"X\",\"Y\", \"VORTICITY\"\n");
	        
	        for(j=0;j<n;j++)
           	   {
           	   	y=j*dely;
	    	      for(i=0;i<m;i++)
	    	      {
	    	      	x=i*delx;
	    	      
	    	      	fprintf(file,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",x,y,u[i][j],v[i][j],psi[i][j],omega[i][j]);
	    	      	fprintf(file1,"%lf\t%lf\t%lf\n",x,y,psi[i][j]);
	    	      	fprintf(file2,"%lf\t%lf\t%lf\t%lf\n",x,y,u[i][j],v[i][j]);
	    	      	fprintf(file3,"%lf\t%lf\n",u[i][j],y);
	    	      	fprintf(file4,"%lf\t%lf\n",x,v[i][j]);
	    	      	fprintf(file5,"%lf\t%lf\t%lf\n",x,y,omega[i][j]);
	    	      }
                }
                fclose(file);
}

