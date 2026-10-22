// Independent numerical reference for the Stan models, not a Stan executable.
// Historical model: centered Gibbs, elliptical slice study effects,
// Metropolis log-SDs and LKJ canonical partial correlations, plus an
// interwoven noncentered SD move. HOPE: Gaussian-prior elliptical slice.
// Compile: g++ -O3 -std=c++17 reference_sampler.cpp -o reference_sampler
// Each chain is independent. No warmup adaptation is retained in sampling.
#include <array>
#include <vector>
#include <random>
#include <fstream>
#include <iostream>
#include <cmath>
#include <string>
#include <algorithm>
#include <iomanip>
using V=std::array<double,6>;
using M=std::array<V,6>;
using C=std::array<double,8>;
struct Row {int j; double y,se,w,active,power;};
struct RNG {
 std::mt19937_64 gen; std::normal_distribution<double> nd{0,1};
 RNG(int s):gen(s){} double n(){return nd(gen);}
 double u(){return std::generate_canonical<double,53>(gen);}
 V normal(){V z;for(auto &a:z)a=n();return z;}
};
M chol(M a) {
 M l{};
 for(int i=0;i<6;i++)for(int j=0;j<=i;j++){
  double v=a[i][j];for(int k=0;k<j;k++)v-=l[i][k]*l[j][k];
  if(i==j){if(v<=0)throw std::runtime_error("Nonpositive covariance");l[i][j]=sqrt(v);}
  else l[i][j]=v/l[j][j];
 } return l;
}
V mul(const M&a,const V&v){V r{};for(int i=0;i<6;i++)for(int j=0;j<6;j++)r[i]+=a[i][j]*v[j];return r;}
V lower_solve(const M&l,V b){for(int i=0;i<6;i++){for(int j=0;j<i;j++)b[i]-=l[i][j]*b[j];b[i]/=l[i][i];}return b;}
V upper_solve(const M&l,V b){for(int i=5;i>=0;i--){for(int j=i+1;j<6;j++)b[i]-=l[j][i]*b[j];b[i]/=l[i][i];}return b;}
M inverse(const M&l){M a{};for(int j=0;j<6;j++){V e{};e[j]=1;V c=upper_solve(l,lower_solve(l,e));for(int i=0;i<6;i++)a[i][j]=c[i];}return a;}
double logistic(double x){return 1/(1+exp(-x));}
C clinical(const V&eta){V u;for(int i=0;i<6;i++)u[i]=logistic(eta[i]);C x;
 x[0]=30*u[0];x[1]=x[0]+(30-x[0])*u[1];x[6]=-x[0]*u[2];
 x[7]=x[6]-(x[1]+x[6])*u[3];x[2]=x[6]*u[4];
 x[3]=x[2]+(x[7]-x[6])*u[5];x[4]=x[6]-x[2];x[5]=x[7]-x[3];return x;}
M partial_L(const std::array<double,15>&v){M l{};l[0][0]=1;int k=0;
 for(int i=1;i<6;i++){double r=1;for(int j=0;j<i;j++){double p=tanh(v[k++]);l[i][j]=p*r;r*=sqrt(1-p*p);}l[i][i]=r;}return l;}
double partial_prior(const std::array<double,15>&v,double shape){int k=0;double out=0;
 for(int i=1;i<6;i++)for(int j=0;j<i;j++){double p=tanh(v[k++]);out+=(shape+(4-j)*.5)*log1p(-p*p);}return out;}
M scale_L(M l,const V&t){for(int i=0;i<6;i++)for(int j=0;j<6;j++)l[i][j]*=t[i];return l;}
double ll_study(const V&e,int study,const std::vector<Row>&ch,const std::vector<Row>&ba){C x=clinical(e);double out=0;
 for(auto r:ch)if(r.j==study){double mu=r.w*(x[2]+r.active*x[4])+(1-r.w)*(x[3]+r.active*x[5]);double d=(r.y-mu)/r.se;out-=.5*r.power*d*d;}
 for(auto r:ba)if(r.j==study){double mu=r.w*x[0]+(1-r.w)*x[1];double d=(r.y-mu)/r.se;out-=.5*d*d;}return out;}
double ll_all(const std::vector<V>&e,const std::vector<Row>&ch,const std::vector<Row>&ba){double out=0;for(int j=0;j<(int)e.size();j++)out+=ll_study(e[j],j,ch,ba);return out;}
void ellipse_study(V&e,const V&mu,const M&l,int j,const std::vector<Row>&ch,const std::vector<Row>&ba,RNG&rng){
 V nu=mul(l,rng.normal()),center;for(int k=0;k<6;k++)center[k]=e[k]-mu[k];
 double cut=ll_study(e,j,ch,ba)+log(rng.u()),theta=rng.u()*2*M_PI,lo=theta-2*M_PI,hi=theta;
 for(int attempt=0;attempt<10000;attempt++){V prop;for(int k=0;k<6;k++)prop[k]=mu[k]+center[k]*cos(theta)+nu[k]*sin(theta);
  if(ll_study(prop,j,ch,ba)>cut){e=prop;return;}if(theta<0)lo=theta;else hi=theta;theta=lo+rng.u()*(hi-lo);}
 throw std::runtime_error("Elliptical slice failed");}
double gaussian_hierarchy(const std::vector<V>&e,const V&mu,const V&t,const M&lo){double logdet=0,quad=0;
 for(int k=0;k<6;k++)logdet+=log(t[k])+log(lo[k][k]);M full=scale_L(lo,t);
 for(auto x:e){for(int k=0;k<6;k++)x[k]-=mu[k];V z=lower_solve(full,x);for(double a:z)quad+=a*a;}
 return -static_cast<double>(e.size())*logdet-.5*quad;}
double tau_prior(const V&t,const V&q){double v=0;for(int k=0;k<6;k++)v+=log(t[k])-.5*pow(t[k]/q[k],2);return v;}
int main(int argc,char**argv){
 if(argc!=8){std::cerr<<"mode input output seed warmup retained thin\n";return 2;}
 std::string mode=argv[1];std::ifstream in(argv[2]);std::ofstream out(argv[3]);out<<std::setprecision(14);
 RNG rng(std::stoi(argv[4]));int warm=std::stoi(argv[5]),keep=std::stoi(argv[6]),thin=std::stoi(argv[7]);
 if(mode=="meta"){
  V m,s,q;double shape;int J,N,NB;for(double &v:m)in>>v;for(double &v:s)in>>v;for(double &v:q)in>>v;in>>shape>>J>>N>>NB;
  std::vector<Row>ch(N),ba(NB);for(auto&r:ch)in>>r.j>>r.y>>r.se>>r.w>>r.active>>r.power;
  for(auto&r:ba)in>>r.j>>r.y>>r.se>>r.w;
  V mu=m,t=q;std::array<double,15>v{};M lo=partial_L(v);std::vector<V>e(J);
  for(int j=0;j<J;j++){e[j]=m;for(int k=0;k<6;k++)e[j][k]+=.2*rng.n();}
  V step;step.fill(.25);std::array<double,15>rstep;rstep.fill(.3);V acc{};std::array<double,15>racc{};
  for(int k=0;k<6;k++)out<<"mu"<<k<<",";for(int k=0;k<6;k++)out<<"tau"<<k<<",";
  for(int k=0;k<15;k++)out<<"partial"<<k<<",";
  for(int j=0;j<J;j++)for(int k=0;k<6;k++)out<<"eta"<<j<<"_"<<k<<",";
  for(int k=0;k<6;k++)out<<"eta_new"<<k<<(k==5?'\n':',');
  for(int it=0;it<warm+keep*thin;it++){
   M full=scale_L(lo,t);
   for(int j=0;j<J;j++)ellipse_study(e[j],mu,full,j,ch,ba,rng);
   // Gaussian conditional for the latent reference-class mean.
   M prec=inverse(full),post_prec=prec;V rhs{},sum{};
   for(auto x:e)for(int k=0;k<6;k++)sum[k]+=x[k];rhs=mul(prec,sum);
   for(int k=0;k<6;k++){post_prec[k][k]*=J;rhs[k]+=m[k]/(s[k]*s[k]);}
   for(int a=0;a<6;a++)for(int b=0;b<6;b++)if(a!=b)post_prec[a][b]*=J;
   for(int k=0;k<6;k++)post_prec[k][k]+=1/(s[k]*s[k]);
   M lp=chol(post_prec);V center=upper_solve(lp,lower_solve(lp,rhs)),noise=upper_solve(lp,rng.normal());
   for(int k=0;k<6;k++)mu[k]=center[k]+noise[k];
   for(int k=0;k<6;k++){
    V prop=t;prop[k]*=exp(step[k]*rng.n());
    double lr=gaussian_hierarchy(e,mu,prop,lo)+tau_prior(prop,q)-gaussian_hierarchy(e,mu,t,lo)-tau_prior(t,q);
    if(log(rng.u())<lr){t=prop;acc[k]++;}
   }
   for(int k=0;k<15;k++){
    auto vp=v;vp[k]+=rstep[k]*rng.n();M lop=partial_L(vp);
    double lr=gaussian_hierarchy(e,mu,t,lop)+partial_prior(vp,shape)-gaussian_hierarchy(e,mu,t,lo)-partial_prior(v,shape);
    if(std::isfinite(lr)&&log(rng.u())<lr){v=vp;lo=lop;racc[k]++;}
   }
   // Interweaving: keep standardized study deviations fixed during SD moves.
   full=scale_L(lo,t);std::vector<V>z(J);
   for(int j=0;j<J;j++){V diff=e[j];for(int k=0;k<6;k++)diff[k]-=mu[k];z[j]=lower_solve(full,diff);}
   for(int k=0;k<6;k++){
    V prop=t;prop[k]*=exp(.35*rng.n());auto ep=e;
    for(int j=0;j<J;j++){ep[j]=mul(scale_L(lo,prop),z[j]);for(int d=0;d<6;d++)ep[j][d]+=mu[d];}
    double lr=ll_all(ep,ch,ba)+tau_prior(prop,q)-ll_all(e,ch,ba)-tau_prior(t,q);
    if(log(rng.u())<lr){t=prop;e=ep;}
   }
   if(it<warm && (it+1)%100==0){double rate=std::min(.05,1/sqrt((it+1)/100.));
    for(int k=0;k<6;k++){step[k]*=exp(rate*(acc[k]/100.-.35));acc[k]=0;}
    for(int k=0;k<15;k++){rstep[k]*=exp(rate*(racc[k]/100.-.35));racc[k]=0;}}
   if(it>=warm && (it-warm)%thin==0){
    for(double a:mu)out<<a<<",";for(double a:t)out<<a<<",";for(double a:v)out<<a<<",";
    for(auto x:e)for(double a:x)out<<a<<",";V en=mul(scale_L(lo,t),rng.normal());
    for(int k=0;k<6;k++)out<<en[k]+mu[k]<<(k==5?'\n':',');
   }
  }
 } else if(mode=="hope") {
  V m;M cov;double em,es;int N,NB;for(double&v:m)in>>v;for(auto&r:cov)for(double&v:r)in>>v;
  in>>em>>es>>N>>NB;std::vector<Row>ch(N),ba(NB);for(auto&r:ch)in>>r.y>>r.se>>r.w>>r.active;
  for(auto&r:ba)in>>r.y>>r.se>>r.w;M l=chol(cov);
  auto ll=[&](const std::array<double,7>&z){V zz;std::copy_n(z.begin(),6,zz.begin());V eta=mul(l,zz);for(int k=0;k<6;k++)eta[k]+=m[k];C x=clinical(eta);double ed=exp(em+es*z[6]),out=0;
   for(auto r:ch){double g=r.active/(r.active+ed),mu=r.w*(x[2]+x[4]*g)+(1-r.w)*(x[3]+x[5]*g);out-=.5*pow((r.y-mu)/r.se,2);}
   for(auto r:ba){double mu=r.w*x[0]+(1-r.w)*x[1];out-=.5*pow((r.y-mu)/r.se,2);}return out;};
  std::array<double,7>z{};for(double&v:z)v=.2*rng.n();out<<"eta0,eta1,eta2,eta3,eta4,eta5,ed50\n";
  for(int it=0;it<warm+keep*thin;it++){
   std::array<double,7>nu,prop;for(double&v:nu)v=rng.n();double cut=ll(z)+log(rng.u()),theta=rng.u()*2*M_PI,lo=theta-2*M_PI,hi=theta;
   for(int attempt=0;attempt<10000;attempt++){
    for(int k=0;k<7;k++)prop[k]=z[k]*cos(theta)+nu[k]*sin(theta);
    if(ll(prop)>cut){z=prop;break;}if(theta<0)lo=theta;else hi=theta;theta=lo+rng.u()*(hi-lo);
    if(attempt==9999)throw std::runtime_error("HOPE elliptical slice failed");
   }
   if(it>=warm && (it-warm)%thin==0){V zz;std::copy_n(z.begin(),6,zz.begin());V eta=mul(l,zz);for(int k=0;k<6;k++)out<<eta[k]+m[k]<<",";out<<exp(em+es*z[6])<<"\n";}
  }
 } else return 3;
 if(!in.good() && !in.eof())return 4;
 return 0;
}
