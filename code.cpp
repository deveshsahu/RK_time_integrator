#include<iostream>
#include<string>
#include<math.h>
#include<cstdlib>
using namespace std;
#include"readfile.h"
#include"writefile.h"

int main()
{
 readfile obj;
 cout<<"\nCreating Field of size: "<<obj.N<<" time nodes x "<<obj.M<<" grid points"<<endl;
 cout<<"Simulation time "<<obj.delT*obj.N<<"s"<<endl;
 float **u = (float **)malloc(sizeof(float *) * (obj.N));
 for(int i=0;i<obj.N;i++){
	u[i] = (float *)malloc(sizeof(float)*(obj.M+1));
	}

/* float **u;
 u = new float *[obj.N];
 for(int i=0;i<obj.N;i++)
	u[i] = new float [obj.M+1];
*/
 float d = obj.c*obj.delT/obj.delx;			//Courant Number
 cout<<"\nCourant number = "<<d;	
 if (d*d<=1)
 { 
	cout<<"\nThe convergence condition is followed!\nIntializing field...."; 
	 for (int j=0;j<=obj.M;j++) 
			u[0][j]= sin(2*j*3.14/obj.M);		//Initializing the values	

/******************************************************************************************************************************
			Range Kutta Order 2 SOLVER starts here
******************************************************************************************************************************/
	 if(obj.ddt=="rk2")
		{
		cout<<"\nSelected RK-2 scheme for temporal derievative...Solving";
										//For positive x direction
	    if(obj.ddx=="upw"){
		cout<<"\nSelected Upwind Scheme for spatial derievative...Solving";
		if(d>0){	
		for (int n=1;n<obj.N;n++) {
			for (int i=1;i<=obj.M;i++) {
										// Half predictor step
				float uhalf = u[n-1][i] - d/2*(u[n-1][i]-u[n-1][i-1]);
										// Full predictor step
				u[n][i] = uhalf - d/2*(u[n-1][i]-u[n-1][i-1]);
				}
			u[n][0] = u[n][obj.M];
			}
		}
										//For negative x direction
		else {	
			for (int n=1;n<obj.N;n++) {
				for (int i=obj.M-1;i>=0;i--) { 
										// Half predictor step
				float uhalf = u[n-1][i] - d/2*(u[n-1][i+1]-u[n-1][i]);
										// Full predictor step
				u[n][i] = uhalf - d/2*(u[n-1][i+1]-u[n-1][i]);
				}
				u[n][obj.M] = u[n][0];
			}
		 }
		}
	else if(obj.ddx=="ctd"){
		cout<<"\nSelected Central Difference Scheme for Spatial Derievative";
		for (int n=1;n<obj.N;n++) {
			for (int i=1;i<obj.M;i++) {
										// Half predictor step
				float uhalf = u[n-1][i] - d/4*(u[n-1][i+1]-u[n-1][i-1]);
										// Full predictor step
				u[n][i] = uhalf - d/4*(u[n-1][i+1]-u[n-1][i-1]);
				}
										// Half predictor step
			float uhalf = u[n-1][i] - d/4*(u[n-1][1]-u[n-1][obj.M-1]);
										// Full predictor step
			u[n][obj.M] = uhalf - d/4*(u[n-1][1]-u[n-1][obj.M-1]);
			u[n][0] = u[n][obj.M];
		}
	}
/******************************************************************************************************************************
			Range Kutta Order 3 SOLVER starts here
******************************************************************************************************************************/
	 if(obj.ddt=="rk3")
		{
		cout<<"\nSelected RK-3 scheme for temporal derievative...Solving";
										//For positive x direction
	    if(obj.ddx=="upw"){
		cout<<"\nSelected Upwind scheme for spatial derievative...Solving";
		if(d>0){	
		for (int n=1;n<obj.N;n++) {
			for (int i=1;i<=obj.M;i++) {
				float k1 = -d*(u[n-1][i] - u[n-1][i-1]);		// h(f(xi,yi)
				float k2 = -d*(u[n-1][i] + k1/2 - u[n-1][i-1]);
				float k3 = -d*(u[n-1][i] - k1 +2*k2 -u[n-1][i-1]);
										// Full predictor step
				u[n][i] = u[n-1][i] + (k1 + 4 * k2 + k3)/6;
				}
			u[n][0] = u[n][obj.M];
			}
		}
		
		else								//For Negative x direction
			{
			 }


		}
		else if(obj.ddx=="ctd"){
		for (int n=1;n<obj.N;n++) {
			for (int i=1;i<obj.M;i++) {
				float k1 = -d/2*(u[n-1][i+1] - u[n-1][i-1]);		// h(f(xi,yi)
				float k2 = -d/2*(u[n-1][i+1] - u[n-1][i-1]);
				float k3 = -d/2*(u[n-1][i+1] - u[n-1][i-1]);
										// Full predictor step
				u[n][i] = u[n-1][i] + (k1 + 4 * k2 + k3)/6;
				}
			float k1 = -d/2*(u[n-1][1] - u[n-1][obj.M-1]);		// h(f(xi,yi)
			float k2 = -d/2*(u[n-1][1] - u[n-1][obj.M-1]);
			float k3 = -d/2*(u[n-1][1] - u[n-1][obj.M-1]);
			u[n][obj.M] = u[n-1][obj.M] + (k1 + 4 * k2 + k3)/6;
			u[n][0] = u[n][obj.M];
			}
		 
			}
	}

/******************************************************************************************************************************
			Range Kutta Order 4 SOLVER starts here
******************************************************************************************************************************/
	else if(obj.ddt=="rk4")
		{
		cout<<"\nSelected RK-4 scheme for temporal derievative...Solving";
										//For positive x direction
	    if(obj.ddx=="upw"){
		cout<<"\nSelected Upwind scheme for spatial derievative...Solving";
		if(d>0){	
		for (int n=1;n<obj.N;n++) {
			for (int i=1;i<=obj.M;i++) {
				float k1 = -d*(u[n-1][i] - u[n-1][i-1]);		// h(f(xi,yi)
				float k2 = -d*(u[n-1][i] + k1/2 - u[n-1][i-1]);
				float k3 = -d*(u[n-1][i] + k2/2 - u[n-1][i-1]);
				float k4 = -d*(u[n-1][i] + k3   - u[n-1][i-1]);
										// Full predictor step
				u[n][i] = u[n-1][i] + (k1 + 2*k2 + 2*k3 + k4)/6;
				}
			u[n][0] = u[n][obj.M];
			}
		}
		
		else								//For Negative x direction
			{
			 }

		}
		else if(obj.ddx=="ctd"){
		cout<<"\nCentral Difference Scheme...Solving";
		for (int n=1;n<obj.N;n++) {
			for (int i=1;i<obj.M;i++) {
				float k1 = -d/2*(u[n-1][i+1] - u[n-1][i-1]);		// h(f(xi,yi)
				float k2 = -d/2*(u[n-1][i+1] - u[n-1][i-1]);
				float k3 = -d/2*(u[n-1][i+1] - u[n-1][i-1]);
				float k4 = -d/2*(u[n-1][i+1] - u[n-1][i-1]);
										// Full predictor step
				u[n][i] = u[n-1][i] + (k1 + 2*k2 + 2*k3 + k4)/6;
				}
			float k1 = -d/2*(u[n-1][1] - u[n-1][obj.M-1]);		// h(f(xi,yi)
			float k2 = -d/2*(u[n-1][1] - u[n-1][obj.M-1]);
			float k3 = -d/2*(u[n-1][1] - u[n-1][obj.M-1]);
			float k4 = -d/2*(u[n-1][1] - u[n-1][i-1]);
			u[n][obj.M] = u[n-1][obj.M] + (k1 + 2*k2 + 2*k3 + k4)/6;
			u[n][0] = u[n][obj.M];
			}
		 
			}
	}		
	else
		{
		cout<<"\nIncorrect temporal scheme.Terminating....\n";
		return 0;
		}

cout<<endl<<"Writing into file...";

writefile w('A');
if(!w.writing(obj.write_interval,u,obj.N,obj.M))
	{
	cout<<"\nCould not open file... Terminating\n";
	return 0;
	}
cout<<"\nSolution written!!\tFile Closed...Exiting...\n";
return 1;
}

else
{
 cout<<"\nCourant number exceeds 1. Condition Not Satisfied!\nExiting...";
 return 0;
}
}
