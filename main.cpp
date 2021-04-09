#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include "arehnus.h"
#include <mpi.h>

using namespace std;


int main(int argc, char *argv[])
{
    cout<<"bienvenu dans le code de calcul pyrolyse"<<endl;


    int rank,nproc,iree;
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);

    int i,j,iBeg,iEnd,taille;


    taille=iEnd-iBeg+1;


    int Nx=5;
    double Lx=1;
    double deltax=Lx/((Nx+1)*1.);

    int tranche;
    int Ny=5;
    double Ly=1;
    double deltay=Ly/((Ny+1)*1.);

    double tf=10;
    double dt=0.1;
    double t=dt;
    double x,y;
    int coef;
    int nb=int(ceil(tf/dt));
    charge_3(rank,Ny,nproc,&iBeg,&iEnd);
    taille=iEnd-iBeg+1;
    //printf("Hello,I am rank %d/%d my iBeg is %d et my iEnd is %d\n",rank,nproc,iBeg,iEnd);

    int kmax=100000;
    double epsilon=0.000000001;

    double erreur=0.;

    vector<double> U(taille*Nx,0);
    vector<double> U1(taille*Nx,0);
    vector<double> X0(1,0);

    vector<double> b(taille*Nx,0);


    zero(U,Nx,Ny);

    for(int n=1;n<=nb;n++){
        zero(b,Nx,Ny);
        remplissage_version(b,t,deltax,deltay,Nx,Ny,iBeg,iEnd);
        U1=GradienConjugue(rank,nproc,X0,plus_vecteur(U,produit_s_vecteur(b,dt,Nx,Ny),Nx,Ny),Nx,Ny,iBeg,iEnd,kmax,epsilon,deltax,deltay,dt,n);
        if(n==1&&rank==0){
            //affiche_vecteur(U1,Nx,Ny);
        }

        print(U1,"numerique ",rank,iBeg,iEnd,Nx,Ny,deltax,deltay,t);
        U=U1;
        t+=dt;
    }

    MPI_Finalize();
    return 0;


}
