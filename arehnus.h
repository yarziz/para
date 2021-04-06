


double u(double x,double y, double t);
double g(double x,double y, double t);
double h(double x,double y, double t);
double f(double x,double y, double t);
int coefficient(int i,int j,int Nx);

//void remplissage(Eigen::SparseMatrix<double>& matrice,Eigen::VectorXd& b,double t,double deltax,double deltay,int Nx,int Ny);
void print_1(std::vector<double> v,std::string nom,int j,int Nx,int Ny,double deltax,double t);
//void print(Eigen::VectorXd v,std::string nom,int Nx,int Ny,double deltax,double deltay,double t);
void remplissage_version(std::vector<double>& b,double t,double deltax,double deltay,int Nx,int Ny);
void print(std::vector<double> v,std::string nom,int Nx,int Ny,double deltax,double deltay,double t);


void identite(std::vector<double>& id,int Nx,int Ny);
void affiche_vecteur(std::vector<double> a,int Nx,int Ny);
void zero(std::vector<double>& id,int Nx,int Ny);
void print_matrice(std::vector<double> matrice,int Nx,int Ny);
std::vector<double> produit(std::vector<double> A,std::vector<double> v,int Nx,int Ny);
std::vector<double> produit_M(std::vector<double> v,int Nx,int Ny,double deltax,double deltay,double dt);
std::vector<double> moins_vecteur(std::vector<double> a,std::vector<double> b,int Nx,int Ny);
std::vector<double> plus_vecteur(std::vector<double> a,std::vector<double> b,int Nx,int Ny);
double produit_scalaire(std::vector<double> a,std::vector<double> b,int Nx,int Ny);
std::vector<double> produit_s_vecteur(std::vector<double> a,double scalaire,int Nx,int Ny);
std::vector<double> GradienConjugue(std::vector<double> x0,std::vector<double> b,int Nx,int Ny, int kmax, double epsilon,double deltax,double deltay,double dt);