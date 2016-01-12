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
			--------------------------------------
For du/dt = f(t,u)
			k1 =  f(t,u)
			k2 =  f(t+h/2,u+h/2*u)
			u_n+1 = u_n + h*k2
******************************************************************************************************************************/
	 if(obj.ddt=="rk2")
		{
		cout<<"\nSelected RK-2 scheme for temporal derievative...Solving";
										//For positive x direction
	    if(obj.ddx=="upw"){
		cout<<"\nSelected Upwind Scheme for spatial derievative...Solving";
		if(d>0){	
		float uhalf[obj.M];
		for (int n=1;n<obj.N;n++) {
			for (int i=1;i<=obj.M;i++) {
				uhalf[i] = u[n-1][i] - d/2*(u[n-1][i] - u[n-1][i-1]);
				}
				uhalf[0] = u[n-1][0] - d/2*(u[n-1][0] - u[n-1][obj.M-1]);
			for (int i=1;i<=obj.M;i++) {
				u[n][i] = u[n-1][i] - d * (uhalf[i] - uhalf[i-1]);
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
		float uhalf[obj.M+1];
		for (int n=1;n<obj.N;n++) {
			for (int i=1;i<obj.M;i++) {
				uhalf[i] = u[n-1][i] - d/4*(u[n-1][i+1] - u[n-1][i-1]);
				}
				uhalf[0] = u[n-1][0] - d/4*(u[n-1][1] - u[n-1][obj.M-1]);
				uhalf[obj.M] = uhalf[0];
			for (int i=1;i<obj.M;i++) {
				u[n][i] = u[n-1][i] - d/2 * (uhalf[i+1] - uhalf[i-1]);
				}
			u[n][obj.M] = u[n-1][obj.M] - d/2 * (uhalf[1] - uhalf[obj.M-1]);
			u[n][0] = u[n][obj.M];
			}
		}
	}
/******************************************************************************************************************************
			Range Kutta Order 3 SOLVER starts here
			--------------------------------------
For du/dt = f(t,u)
				k1 = f(t,u)
				k2 = f(t+dt/2,u + dt*k1/2)
				k3 = f(t+dt,u + dt*(2*k2 - k1))
			u_n+1 = u_n + dt/6*(k1 + 4* k2 + k3)
******************************************************************************************************************************/
	else if(obj.ddt=="rk3")
		{
		cout<<"\nSelected RK-3 scheme for temporal derievative...Solving";
										//For positive x direction
	    if(obj.ddx=="upw"){
		cout<<"\nSelected Upwind scheme for spatial derievative...Solving";
		float uhalf[obj.M+1],ufull[obj.M+1];
		if(d>0){	
		for (int n=1;n<obj.N;n++) {
			for (int i=1;i<=obj.M;i++) {
				uhalf[i] = u[n-1][i] -d*(u[n-1][i] - u[n-1][i-1]);
				}
			uhalf[0] = uhalf[obj.M];
			for (int i=1;i<=obj.M;i++) {
				ufull[i] = u[n-1][i] -d*(2*(uhalf[i] - uhalf[i-1])-(u[n-1][i]-u[n-1][i-1]));
				}
			ufull[0] = ufull[obj.M];
			for (int i=1;i<=obj.M;i++) {
				u[n][i] = u[n-1][i] - 	d*( (u[n-1][i]-u[n-1][i-1])   	// k1
							+ 4 * (uhalf[i] - uhalf[i-1])	// k2
							+ (ufull[i] - ufull[i-1]) )/6;	// k3
				}
			u[n][0] = u[n][obj.M];
			}
		}
		
		else								//For Negative x direction
			{
		for (int n=1;n<obj.N;n++) {
			for (int i=obj.M-1;i>=0;i--) {
				uhalf[i] = u[n-1][i] -d*(u[n-1][i] - u[n-1][i-1]);
				}
			uhalf[0] = uhalf[obj.M];
			for (int i=obj.M-1;i>=0;i--) {
				ufull[i] = u[n-1][i] -d*(2*(uhalf[i] - uhalf[i-1])-(u[n-1][i]-u[n-1][i-1]));
				}
			ufull[0] = ufull[obj.M];
			for (int i=obj.M-1;i>=0;i--) {
				u[n][i] = u[n-1][i] - 	d*( (u[n-1][i]-u[n-1][i-1])   	// k1
							+ 4 * (uhalf[i] - uhalf[i-1])	// k2
							+ (ufull[i] - ufull[i-1]) )/6;	// k3
				}
			u[n][obj.M]=u[n][0];
			}
			 
		}


		}
		else if(obj.ddx=="ctd"){
		cout<<"\nSelected Central Difference Scheme for spatial derievative...Solving";
		float uhalf[obj.M+1],ufull[obj.M+1];
		for (int n=1;n<obj.N;n++) {
			for (int i=1;i<obj.M;i++) {
				uhalf[i] = u[n-1][i] -d/4*(u[n-1][i+1] - u[n-1][i-1]);		//Evaluating U_n+1/2
				}
			uhalf[obj.M] = u[n-1][obj.M] -d/2*(u[n-1][1] - u[n-1][obj.M-1]);
			uhalf[0] = uhalf[obj.M];
			for (int i=1;i<obj.M;i++) {
				ufull[i] = u[n-1][i] -d/2*(2*(uhalf[i+1] - uhalf[i-1])-(u[n-1][i+1]-u[n-1][i-1]));
												//Evaluating U_n+1
				}
			ufull[obj.M] = u[n-1][obj.M] -d/2*(2*(uhalf[1] - uhalf[obj.M-1])-(u[n-1][1]-u[n-1][obj.M-1]));
			ufull[0] = ufull[obj.M];
			for (int i=1;i<obj.M;i++) {
				u[n][i] = u[n-1][i] - 	d/2*( (u[n-1][i+1]-u[n-1][i-1]) 	// k1 
							+ 4 * (uhalf[i+1] - uhalf[i-1])		// k2
							+ (ufull[i+1] - ufull[i-1]) )/6;	// k3
				}
			u[n][obj.M] = u[n-1][obj.M] - 	d/2*( (u[n-1][1]-u[n-1][obj.M-1])  
						+ 4 * (uhalf[1] - uhalf[obj.M-1])
						+ (ufull[1] - ufull[obj.M-1]) )/6;
			u[n][0] = u[n][obj.M];
			}
		}
	}

/******************************************************************************************************************************
			Range Kutta Order 4 SOLVER starts here
			--------------------------------------
For du/dt = f(t,u)
				k1 = f(t,u)
				k2 = f(t+dt/2,u + dt*k1/2)
				k3 = f(t+dt/2,u + dt*k2/2)
				k4 = f(t+dt,u + dt*k3)
			u_n+1 = u_n + dt/6*(k1 + 2*k2 + 2*k3 + k4)
******************************************************************************************************************************/
	else if(obj.ddt=="rk4")
		{
		cout<<"\nSelected RK-4 scheme for temporal derievative...Solving";
										//For positive x direction
	    if(obj.ddx=="upw"){
		cout<<"\nSelected Upwind scheme for spatial derievative...Solving";
		float uhalf1[obj.M+1],uhalf2[obj.M+1],ufull[obj.M+1];
		if(d>0){	
		for (int n=1;n<obj.N;n++) {
			for (int i=1;i<=obj.M;i++) {
				uhalf1[i] = u[n-1][i] -d/2*(u[n-1][i] - u[n-1][i-1]);		//Evaluating U_n+1/2 with k1 as slope
				}
			uhalf1[0] = uhalf1[obj.M];
			for (int i=1;i<=obj.M;i++) {
				uhalf2[i] = u[n-1][i] -d/2*(uhalf1[i] - uhalf1[i-1]);		//Evaluating U_n+1/2 with k2 as slope
				}
			uhalf2[0] = uhalf2[obj.M];
			for (int i=1;i<=obj.M;i++) {
				ufull[i] = u[n-1][i] -d*(uhalf2[i] - uhalf2[i-1]);		//Evaluating U_n+1 with k3 as slope
				}
			ufull[0] = ufull[obj.M];

			for (int i=1;i<=obj.M;i++) {
				u[n][i] = u[n-1][i] - d* ((u[n-1][i] - u[n-1][i-1])		//k1
							 + 2*(uhalf1[i] - uhalf1[i-1])		//k2
							 + 2*(uhalf2[i] - uhalf2[i-1])		//k3
							 + (ufull[i] - ufull[i-1]))/6;		//k4
				}
			u[n][0] = u[n][obj.M];
			}
		}
		
		else								//For Negative x direction
			{
		for (int n=1;n<obj.N;n++) {
			for (int i=obj.M-1;i>=0;i--) {
				uhalf1[i] = u[n-1][i] -d/2*(u[n-1][i] - u[n-1][i-1]);		//Evaluating U_n+1/2 with k1 as slope
				}
			uhalf1[0] = uhalf1[obj.M];
			for (int i=obj.M-1;i>=0;i--) {
				uhalf2[i] = u[n-1][i] -d/2*(uhalf1[i] - uhalf1[i-1]);		//Evaluating U_n+1/2 with k2 as slope
				}
			uhalf2[0] = uhalf2[obj.M];
			for (int i=obj.M-1;i>=0;i--) {
				ufull[i] = u[n-1][i] -d*(uhalf2[i] - uhalf2[i-1]);		//Evaluating U_n+1 with k3 as slope
				}
			ufull[0] = ufull[obj.M];

			for (int i=obj.M-1;i>=0;i--) {
				u[n][i] = u[n-1][i] - d* ((u[n-1][i+1] - u[n-1][i])		//k1
							 + 2*(uhalf1[i+1] - uhalf1[i])		//k2
							 + 2*(uhalf2[i+1] - uhalf2[i])		//k3
							 + (ufull[i+1] - ufull[i]))/6;		//k4
				}
			u[n][obj.M]=u[n][0];
			}
			 }

		}
		else if(obj.ddx=="ctd"){
		cout<<"\nSelected Central Difference Scheme...Solving";
		float uhalf1[obj.M+1],uhalf2[obj.M+1],ufull[obj.M+1];
		for (int n=1;n<obj.N;n++) {
			for (int i=1;i<=obj.M;i++) {
				uhalf1[i] = u[n-1][i] -d/4*(u[n-1][i+1] - u[n-1][i-1]);		//Evaluating U_n+1/2 with k1 as slope
				}
			uhalf1[obj.M] = u[n-1][obj.M] -d/4*(u[n-1][1] - u[n-1][obj.M-1]);
			uhalf1[0] = uhalf1[obj.M];
			for (int i=1;i<=obj.M;i++) {
				uhalf2[i] = u[n-1][i] -d/4*(uhalf1[i+1] - uhalf1[i-1]);		//Evaluating U_n+1/2 with k2 as slope
				}
			uhalf2[obj.M] = u[n-1][obj.M] -d/4*(uhalf1[1] - uhalf1[obj.M-1]);
			uhalf2[0] = uhalf2[obj.M];
			for (int i=1;i<=obj.M;i++) {
				ufull[i] = u[n-1][i] -d/2*(uhalf2[i+1] - uhalf2[i-1]);		//Evaluating U_n+1 with k2 as slope
				}
			ufull[obj.M] = u[n-1][obj.M] -d/2*(uhalf2[1] - uhalf2[obj.M-1]);
			ufull[0] = ufull[obj.M];
			for (int i=1;i<=obj.M;i++) {
				u[n][i] = u[n-1][i] - d/2* ((u[n-1][i+1] - u[n-1][i-1])		// k1
							 + 2*(uhalf1[i+1] - uhalf1[i-1])	// k2	
							 + 2*(uhalf2[i+1] - uhalf2[i-1])	// k3	
							 + (ufull[i+1] - ufull[i-1]))/6;	// k4
				}
			u[n][obj.M] = u[n-1][obj.M] - d/2* ((u[n-1][1] - u[n-1][obj.M-1])
						 + 2*(uhalf1[1] - uhalf1[obj.M-1])
						 + 2*(uhalf2[1] - uhalf2[obj.M-1])
						 + (ufull[1] - ufull[obj.M-1]))/6;
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
