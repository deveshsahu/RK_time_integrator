//Reads Data from controlDict file. Class read stores the settings
#include<fstream>
#include<iostream>
using namespace std;
class readfile
{
public:
 float delT;
 float delx;
 float len; int M;
 int N;
 float c;
 char x[4];
 string ddt;
 string ddx;
 int write_interval;
readfile()
{ 
 ifstream fs;
 fs.open("controlDict");
 string temp;
 //Read delta t
 fs>>temp; fs>>delT;
 //Read delta x 
 fs>>temp>>delx;
 //Read  domain length
 fs>>temp>>len;
 M = len/delx;
  //Read total time steps
 fs>>temp>>N;
  //Read wave speed
 fs>>temp>>c;
  //Read writing interval
 fs>>temp>>write_interval;
  //Read time scheme
 fs>>temp>>ddt;
  //Read space scheme
 fs>>temp>>ddx;
 fs.close();
}
~readfile() {};
};
 

