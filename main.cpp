#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include "arehnus.h"

using namespace std;


int main()
{
    cout<<"bienvenu dans le code de calcul pyrolyse"<<endl;
    int Nx=50;
    double Lx=1;
    double deltax=Lx/((Nx+1)*1.);
    
    int tranche;
    int Ny=50;
    double Ly=1;
    double deltay=Ly/((Ny+1)*1.);

    double tf=10;
    double dt=0.1;
    double t=dt;
    double x,y;
    int coef;
    int nb=int(ceil(tf/dt));


    int kmax=100000;
    double epsilon=0.000000001;

    double erreur=0.;
   
    vector<double> U(Nx*Ny,0);
    vector<double> U1(Nx*Ny,0);
    vector<double> X0(Nx*Ny,1);
    //vector<double> matrice(Nx*Nx*Ny*Ny,0);
    //vector<double> id(Nx*Nx*Ny*Ny,0);
    //vector<double> A(Nx*Nx*Ny*Ny,0);
    vector<double> b(Nx*Ny,0);
    vector<double> exacte(Nx*Ny,0);
    
    

    

    //identite(id,Nx,Ny);
    zero(U,Nx,Ny);
    for(int n=1;n<=nb;n++){
        //zero(matrice,Nx,Ny);
        zero(b,Nx,Ny);
        remplissage_version(b,t,deltax,deltay,Nx,Ny);
        U1=GradienConjugue(X0,plus_vecteur(U,produit_s_vecteur(b,dt,Nx,Ny),Nx,Ny),Nx,Ny,kmax,epsilon,deltax,deltay,dt);
        for(int j=1;j<=Ny;j++){
            for(int i=1;i<=Nx;i++){
                coef=coefficient(i,j,Nx);
                x=i*deltax;
                y=j*deltay;
                exacte[coef]=u(x,y,t);
            }
        }
        //print_1(U1,"numerique",5,Nx,Ny,deltax,t);
        //print_1(exacte,"exacte",5,Nx,Ny,deltax,t);
        print(U1,"numerique",Nx,Ny,deltax,deltay,t);
        U=U1;
        t+=dt;
    }
    return 0;
}