#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include <vector>
#include "arehnus.h"
#include <cstdlib>
#include <mpi.h>

using namespace std;



void charge_3(int me,int n,int np,int *iBeg, int *iEnd)
{
    int r=n%np;
    if(me<r){
       *iBeg=(n/np+1)*me;
       *iEnd=(n/np+1)*(me+1)-1;
    }else{
       *iBeg=r+me*(n/np);
       *iEnd=r+(me+1)*(n/np)-1;
    }
}

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

void remplissage_version(vector<double>& b,double t,double deltax,double deltay,int Nx,int Ny,int iBeg,int iEnd){
    int coef;
    int k=0;
    double x;
    double y;
    for (int j=iBeg+1; j<=iEnd+1;j++){
        for(int i=1; i<=Nx;i++){
            coef=coefficient(i,j,Nx);
            x=i*deltax;
            y=j*deltay;
            if(i==1){
                b[k]+=f(x,y,t)+g(0,y,0)/(deltax*deltax*1.);
            }
            if(i==Nx){
                b[k]+=f(x,y,t)+g(1,y,0)/(deltax*deltax*1.);
            }
            if(i!=Nx&&i!=1){

                b[k]+=f(x,y,t);
            }


            if(j==1){
                b[k]+=h(x,0,t)/(deltay*deltay*1.);
            }
            if(j==Ny){
                b[k]+=h(x,1,t)/(deltay*deltay*1.);
            }
            if(j!=Ny&&j!=1){
            }
            k+=1;
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


vector<double> produit_M(vector<double> v,int Nx,int Ny,int iBeg,int iEnd,double deltax,double deltay,double dt){

    int taille=(iEnd-iBeg+1)*Nx;
    vector<double> res(taille,0);
    int coef;
    int k=0;
    for(int j=iBeg+1;j<=iEnd+1;j++){
        for(int i=1;i<=Nx;i++){
            coef=coefficient(i,j,Nx);
            if(i==1){
                res[k]+=(dt*(2./(deltax*deltax*1.))+1.)*v[coef-iBeg*Nx];
                res[k]-=dt*(1./(deltax*deltax*1.))*v[coef-iBeg*Nx+1];
            }
            if(i==Nx){
                res[k]+=(dt*(2./(deltax*deltax*1.))+1.)*v[coef-iBeg*Nx];
                res[k]-=dt*(1./(deltax*deltax*1.))*v[coef-iBeg*Nx-1];
            }
            if(i!=Nx&&i!=1){
                res[k]+=(dt*(2./(deltax*deltax*1.))+1)*v[coef-iBeg*Nx];
                res[k]-=dt*(1./(deltax*deltax*1.))*v[coef-iBeg*Nx+1];
                res[k]-=dt*(1./(deltax*deltax*1.))*v[coef-iBeg*Nx-1];
            }


            if(j==1){
                res[k]+=dt*(2./(deltay*deltay*1.))*v[coef-iBeg*Nx];
                if(0<=coef-iBeg*Nx+Nx<taille){
                  res[k]-=dt*(1./(deltay*deltay*1.))*v[coef-iBeg*Nx+Nx];
                }
            }


            if(j==Ny){
                res[k]+=dt*(2./(deltay*deltay*1.))*v[coef-iBeg*Nx];
                if(0<=coef-iBeg*Nx-Nx<taille){
                  res[k]-=dt*(1./(deltay*deltay*1.))*v[coef-iBeg*Nx-Nx];
                }
            }
            if(j!=Ny&&j!=1){
                res[k]+=dt*(2./(deltay*deltay*1.))*v[coef-iBeg*Nx];
                if(0<=coef-iBeg*Nx+Nx<taille){
                  res[k]-=dt*(1./(deltay*deltay*1.))*v[coef-iBeg*Nx+Nx];
                }
                if(0<=coef-iBeg*Nx-Nx<taille){
                  res[k]-=dt*(1./(deltay*deltay*1.))*v[coef-iBeg*Nx-Nx];
                }
            }
            k+=1;
        }
    }
    return res;
}


vector<double> moins_vecteur(vector<double> a,vector<double> b,int Nx,int Ny){
    int n=b.size();
    vector<double> res(n,0);
    for(int i=0;i<n;i++){
        res[i]=a[i]-b[i];
    }
    return res;
}

void moins_vecteur_haut(vector<double>& a,vector<double> b,int Nx,int Ny,int taille,double valeur){
    for(int i=0;i<Nx;i++){
        a[i]-=valeur*b[i];
    }
    //std::cout<<"moin vecteur haut bon"<<a.size()<<" "<<b.size()<<std::endl;
}


void moins_vecteur_bas(vector<double>& a,vector<double> b,int Nx,int Ny,int taille,double valeur){
    for(int i=1;i<=Nx;i++){
        a[taille-i]-=valeur*b[Nx-i];
    }
    //std::cout<<"moin vecteur bas bon"<<a.size()<<" "<<b.size()<<" taille est"<<taille<<std::endl;
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
    int n=a.size();
    for(int i=0;i<n;i++){
        res+=a[i]*b[i];
    }
    return res;
}

vector<double> GradienConjugue(int rank,int nproc,vector<double> x0, vector<double> b,int Nx,int Ny,int iBeg,int iEnd,int kmax, double epsilon,double deltax,double deltay,double dt,int n){
  //int n=x0.size();
    int k=0;
    double alpha;
    int taille=(iEnd-iBeg+1)*Nx;
    vector<double> r_tranche(taille,0),rn_tranche(taille,0),p_tranche(taille,0),z_tranche(taille,0),x_tranche(taille,0);
    vector<double> p_send_debut(Nx,0),p_send_fin(Nx,0);
    vector<double> p_rec_debut(Nx,0),p_rec_fin(Nx,0);
    vector<double> r_globale(Nx*Ny,0),p_globale(Nx*Ny,0);
    double z_p=0,g_z_p=0;
    double r_r=0,g_r_r=0;
    double rn_rn=0,g_rn_rn=0;
    double beta=1,gamma=0;
    double valeur=dt*(1./(deltay*deltay*1.));
    int iBeg_r;
    int iEnd_r;
    r_tranche=b;
    p_tranche=r_tranche;
    for(int i=1;i<=Nx;i++){
      p_send_debut[i-1]=p_tranche[i-1];
      p_send_fin[Nx-i]=p_tranche[taille-i];
    }



    //std::cout<<n<<std::endl;
    //beta=sqrt(produit_scalaire(r,r,Nx,Ny));
    while((beta>epsilon)&&(k<kmax+1)){
      //std::cout<<rank<<" "<<"gradien"<<std::endl;
      z_tranche=produit_M(p_tranche,Nx,Ny,iBeg,iEnd,deltax,deltay,dt);
      //std::cout<<rank<<" "<<"produit bon"<<std::endl;
      if(rank==0){
        MPI_Send(&p_send_fin[0],Nx,MPI_DOUBLE,1,0,MPI_COMM_WORLD);
        MPI_Recv(&p_rec_fin[0],Nx,MPI_DOUBLE,1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        //std::cout<<rank<<" "<<"send rec bon"<<std::endl;
        moins_vecteur_bas(z_tranche,p_rec_fin,Nx,Ny,taille,valeur);
        //std::cout<<rank<<" "<<"bon z"<<std::endl;
      }

      if(rank==nproc-1){
        MPI_Send(&p_send_debut[0],Nx,MPI_DOUBLE,rank-1,0,MPI_COMM_WORLD);
        MPI_Recv(&p_rec_debut[0],Nx,MPI_DOUBLE,rank-1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        //std::cout<<rank<<" "<<"send rec bon"<<std::endl;
        moins_vecteur_haut(z_tranche,p_rec_debut,Nx,Ny,taille,valeur);
        //std::cout<<rank<<" "<<"bon z"<<std::endl;
      }

      if(rank!=0&&rank!=nproc-1){
        MPI_Send(&p_send_fin[0],Nx,MPI_DOUBLE,rank+1,0,MPI_COMM_WORLD);
        MPI_Send(&p_send_debut[0],Nx,MPI_DOUBLE,rank-1,0,MPI_COMM_WORLD);
        MPI_Recv(&p_rec_fin[0],Nx,MPI_DOUBLE,rank+1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        MPI_Recv(&p_rec_debut[0],Nx,MPI_DOUBLE,rank-1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        std::cout<<rank<<" "<<"send rec bon"<<std::endl;
        moins_vecteur_haut(z_tranche,p_rec_debut,Nx,Ny,taille,valeur);
        moins_vecteur_bas(z_tranche,p_rec_fin,Nx,Ny,taille,valeur);
      }



      z_p=produit_scalaire(z_tranche,p_tranche,Nx,Ny);
      MPI_Allreduce(&z_p, &g_z_p, 1,MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      r_r=produit_scalaire(r_tranche,r_tranche,Nx,Ny);
      MPI_Allreduce(&r_r, &g_r_r, 1,MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	    alpha=g_r_r/(g_z_p*1.);
	    x_tranche=plus_vecteur(x_tranche,produit_s_vecteur(p_tranche,alpha,Nx,Ny),Nx,Ny);
      rn_tranche=moins_vecteur(r_tranche,produit_s_vecteur(z_tranche,alpha,Nx,Ny),Nx,Ny);
      rn_rn=produit_scalaire(rn_tranche,rn_tranche,Nx,Ny);
      MPI_Allreduce(&rn_rn, &g_rn_rn, 1,MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	    gamma=g_rn_rn/(g_r_r*1.);
      p_tranche=plus_vecteur(rn_tranche,produit_s_vecteur(p_tranche,gamma,Nx,Ny),Nx,Ny);
	    beta=sqrt(g_r_r);
      if(n==1&& rank==1&&k==1){
        //std::cout<<z_p<<" "<<z_p<<endl;
        affiche_vecteur(z_tranche,Nx,Ny);
      }
	    r_tranche=rn_tranche;
      for(int l=1;l<=Nx;l++){
        p_send_debut[l-1]=p_tranche[l-1];
        p_send_fin[Nx-l]=p_tranche[taille-l];
      }
      k+=1;
    }
    return x_tranche;
}

void affiche_vecteur(vector<double> a,int Nx,int Ny){
    for(int i=0;i<a.size();i++){
        printf(" %f",a[i]);
    }
    printf("\n");
}

void print(vector<double> v,string nom,int rank,int iBeg,int iEnd,int Nx,int Ny,double deltax,double deltay,double t){
    ofstream mon_flux;
    string name_file=nom+to_string(rank)+" à l'instant "+to_string(t);
    mon_flux.open(name_file, ios::out);
    double x;
    double y;
    int coef;
    int k=0;
    for(int j=iBeg+1;j<=iEnd+1;j++){
        for(int i=1;i<=Nx;i++){
            if(t==0.1){
              //std::cout<<rank<<" en train d'ecrire"<<endl;
            }
            x=i*deltax;
            y=j*deltay;
            coef=coefficient(i,j,Nx);
            mon_flux<<x<<" "<<y<<" "<<std::setprecision(16)<< v[k] << endl;
            k+=1;
        }
    }
    mon_flux.close();
}
