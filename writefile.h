#include<fstream>
using namespace std;
class writefile
{
 int ctr;
 char name,filename[2];
public:
writefile(char a)
{
 name = a;
 filename[0] = name;
}


void increment_filename()
{
 name = name +1;
 filename[0] = name;
}

 int writing(int interval,float *u[],int N,int M)
 {
 ofstream fs;
for(int x=0,i=0;i<N;x=x++,i=i+interval){
 fs.open(filename);
 if(!fs.is_open())
	return 0;
	for(int j=0;j<=M;j++) {
		fs<<j<<" "<<u[i][j]<<"\n";
	}fs<<endl;
 fs.close();
 increment_filename();
 }
 return 1;
}
};
 
