#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include <vector>
#include "arehnus.h"
#include <cstdlib>

using namespace std;





double u(double x,double y, double t){
    return x*(1-x)*y*(1-y);
    //return sin(x)+cos(y);
    
}

double g(double x,double y, double t){
    //return sin(x)+cos(y);
    //return cos(x)*x*x+cos(y);
    return 0;
}

double h(double x,double y, double t){
    //return sin(x)+cos(y);
    return 1;
}

double g2(double x,double y, double t){
    //return 2*(1-2*x);
    return 0;
}

double h2(double x,double y, double t){
    //return 2*(1-2*y);
    return 0;
}

double f(double x,double y, double t){
    //return 2*(y-y*y+x-x*x);
    return exp(-pow((x-0.5),2))*exp(-pow((y-0.5),2))*cos(M_PI*0.5*t);
    //return ;
    //return sin(x)+cos(y);
}

int coefficient(int i,int j,int Nx){
    return (j-1)*Nx+i-1;
}


void print_1(vector<double> v,string nom,int j,int Nx,int Ny,double deltax,double t){
    ofstream mon_flux;
    string name_file=nom+" la trache "+to_string(j) +" à l'instant "+to_string(t);
    mon_flux.open(name_file, ios::out);
    double x;
    int coef;
    for(int i=1;i<=Nx;i++){
        x=i*deltax;
        coef=coefficient(i,j,Nx);
        mon_flux<<x<<" "<<std::setprecision(16)<< v[coef] << endl;
    }
    mon_flux.close();
}







void identite(vector<double>& id,int Nx,int Ny){
    for (int i = 0; i < Nx*Ny; i++)
    {
        id[i*Nx*Ny+i]=1.;
    }
    
}

void zero(vector<double>& id,int Nx,int Ny){
    for (int i = 0; i < id.size(); i++){
        id[i]=0.;
    }
    
}

void print_matrice(vector<double> matrice,int Nx,int Ny){
    for(int i=0;i<Nx*Ny;i++){
        for(int j=0;j<Nx*Ny;j++){
            printf(" %f",matrice[i*Nx*Ny+j]);
        }
        printf("\n");
    }
}

void remplissage_version(vector<double>& b,double t,double deltax,double deltay,int Nx,int Ny){
    int coef;
    double x;
    double y;
    for (int j=1; j<=Ny;j++){
        for(int i=1; i<=Nx;i++){
            coef=coefficient(i,j,Nx);
            x=i*deltax;
            y=j*deltay;
            if(i==1){
                //matrice.coeffRef(coef,coef)+=2./(deltax*deltax*1.);
                //matrice.coeffRef(coef,coef+1)-=1./(deltax*deltax*1.);
                b[coef]+=f(x,y,t)+g(0,y,0)/(deltax*deltax*1.);

                //neuman
                //matrice.coeffRef(coef,coef)+=1./(deltax*deltax*1.);
                //b(coef)+=f(x,y,t)-g(0,y,t)/(deltax*1.)-g2(0,y,t)/2.;
            }
            if(i==Nx){
                //matrice.coeffRef(coef,coef)+=2./(deltax*deltax*1.);
                //matrice.coeffRef(coef,coef-1)-=1./(deltax*deltax*1.);
               
                b[coef]+=f(x,y,t)+g(1,y,0)/(deltax*deltax*1.);
                //neuman
                //matrice.coeffRef(coef,coef)+=1./(deltax*deltax*1.);
                //b(coef)+=f(x,y,t)+g(1,y,t)/(deltax*1.)-g2(1,y,t)/2.;//1=Lx
            }
            if(i!=Nx&&i!=1){
               
                b[coef]+=f(x,y,t);
            }


            if(j==1){
                //matrice.coeffRef(coef,coef)+=2./(deltay*deltay*1.);
                //matrice.coeffRef(coef,coef+Nx)-=1./(deltay*deltay*1.);
                
                b[coef]+=h(x,0,t)/(deltay*deltay*1.);
                //matrice.coeffRef(coef,coef)+=1./(deltay*deltay*1.);
                //b[coef]-=h(x,0,t)/(deltay*1.)+h2(x,0,t)/2.;
            }
            if(j==Ny){
                //matrice.coeffRef(coef,coef)+=2./(deltay*deltay*1.);
                //matrice.coeffRef(coef,coef-Nx)-=1./(deltay*deltay*1.);
                //b[coef]+=h(x,1,t)/(deltay*deltay*1.);
                
                b[coef]+=h(x,1,t)/(deltay*deltay*1.);
                //matrice.coeffRef(coef,coef)+=1./(deltay*deltay*1.);
                //b[coef]+=h(x,1,t)/(deltay*1.)-h2(x,1,t)/2.;//1=Ly
            }
            if(j!=Ny&&j!=1){
                /*
                matrice.coeffRef(coef,coef)+=2./(deltay*deltay*1.);
                matrice.coeffRef(coef,coef-Nx)-=1./(deltay*deltay*1.);
                matrice.coeffRef(coef,coef+Nx)-=1./(deltay*deltay*1.);
                */
               
            }
        }
    }
}

vector<double> produit(vector<double> A,vector<double> v,int Nx,int Ny){
    vector<double> res(Nx*Ny,0);
    for(int i=0;i<Nx*Ny;i++){
        for(int j=0;j<Nx*Ny;j++){
            res[i]+=A[i*Nx*Ny+j]*v[j];
        }
    }
    return res;
}


vector<double> produit_M(vector<double> v,int Nx,int Ny,double deltax,double deltay,double dt){
    vector<double> res(Nx*Ny,0);
    int coef;
    for(int j=1;j<=Ny;j++){
        for(int i=1;i<=Nx;i++){
            coef=coefficient(i,j,Nx);
            if(i==1){
                res[coef]+=(dt*(2./(deltax*deltax*1.))+1.)*v[coef];
                res[coef]-=dt*(1./(deltax*deltax*1.))*v[coef+1];
            }
            if(i==Nx){
                res[coef]+=(dt*(2./(deltax*deltax*1.))+1.)*v[coef];
                res[coef]-=dt*(1./(deltax*deltax*1.))*v[coef-1];
            }
            if(i!=Nx&&i!=1){
                res[coef]+=(dt*(2./(deltax*deltax*1.))+1)*v[coef];
                res[coef]-=dt*(1./(deltax*deltax*1.))*v[coef+1];
                res[coef]-=dt*(1./(deltax*deltax*1.))*v[coef-1];
            }


            if(j==1){
                res[coef]+=dt*(2./(deltay*deltay*1.))*v[coef];
                res[coef]-=dt*(1./(deltay*deltay*1.))*v[coef+Nx];
            }


            if(j==Ny){
                res[coef]+=dt*(2./(deltay*deltay*1.))*v[coef];
                res[coef]-=dt*(1./(deltay*deltay*1.))*v[coef-Nx];
            }
            if(j!=Ny&&j!=1){
                res[coef]+=dt*(2./(deltay*deltay*1.))*v[coef];
                res[coef]-=dt*(1./(deltay*deltay*1.))*v[coef+Nx];
                res[coef]-=dt*(1./(deltay*deltay*1.))*v[coef-Nx];
            }
        }
    }
    return res;
}


vector<double> moins_vecteur(vector<double> a,vector<double> b,int Nx,int Ny){
    int n=a.size();
    vector<double> res(n,0);
    for(int i=0;i<n;i++){
        res[i]=a[i]-b[i];
    }
    return res;
}

vector<double> plus_vecteur(vector<double> a,vector<double> b,int Nx,int Ny){
    int n=a.size();
    vector<double> res(n,0);
    for(int i=0;i<n;i++){
        res[i]=a[i]+b[i];
    }
    return res;
}

vector<double> produit_s_vecteur(vector<double> a,double scalaire,int Nx,int Ny){
    int n=a.size();
    vector<double> res(n,0);
    for(int i=0;i<n;i++){
        res[i]=scalaire*a[i];
    }
    return res;
}

double produit_scalaire(vector<double> a,vector<double> b,int Nx,int Ny){
    double res=0;
    for(int i=0;i<Nx*Ny;i++){
        res+=a[i]*b[i];
    }
    return res;
}

vector<double> GradienConjugue(vector<double> x0, vector<double> b,int Nx,int Ny, int kmax, double epsilon,double deltax,double deltay,double dt){
  //int n=x0.size();
    int k=0;
    double alpha;
    vector<double> r(Nx*Ny,0),rn(Nx*Ny,0),p(Nx*Ny,0),z(Nx*Ny,0),x(Nx*Ny,0);
    double beta=0,gamma=0;
    x=x0;
    r=moins_vecteur(b,produit_M(x0,Nx,Ny,deltax,deltay,dt),Nx,Ny);
    p=r;
    beta=sqrt(produit_scalaire(r,r,Nx,Ny));
    while((beta>epsilon)&&(k<kmax+1)){
	    //z=produit(A,p,Nx,Ny);
        z=produit_M(p,Nx,Ny,deltax,deltay,dt);
	    alpha=produit_scalaire(r,r,Nx,Ny)/produit_scalaire(z,p,Nx,Ny);
	    x=plus_vecteur(x,produit_s_vecteur(p,alpha,Nx,Ny),Nx,Ny);
        rn=moins_vecteur(r,produit_s_vecteur(z,alpha,Nx,Ny),Nx,Ny);
	    gamma=produit_scalaire(rn,rn,Nx,Ny)/produit_scalaire(r,r,Nx,Ny);
        p=plus_vecteur(rn,produit_s_vecteur(p,gamma,Nx,Ny),Nx,Ny);
	    beta=sqrt(produit_scalaire(r,r,Nx,Ny));
	    k=k+1;
	    r=rn;
        k+=1;
    }
    return x;
}

void affiche_vecteur(vector<double> a,int Nx,int Ny){
    for(int i=0;i<Nx*Ny;i++){
        printf(" %f",a[i]);
    }
    printf("\n");
}

void print(vector<double> v,string nom,int Nx,int Ny,double deltax,double deltay,double t){
    ofstream mon_flux;
    string name_file=nom+" à l'instant "+to_string(t);
    mon_flux.open(name_file, ios::out);
    double x;
    double y;
    int coef;
    for(int j=1;j<=Ny;j++){
        for(int i=1;i<=Nx;i++){
            x=i*deltax;
            y=j*deltay;
            coef=coefficient(i,j,Nx);
            mon_flux<<x<<" "<<y<<" "<<std::setprecision(16)<< v[coef] << endl;
        }
    }
    mon_flux.close();
}

















