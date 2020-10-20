#include <random>
#include <array>
#include <algorithm>
#include <functional>
#include <math.h>
#include <iostream>
#include <chrono>
#include <unistd.h>
#include <cmath>
#include <omp.h>


#define NUMTHREADS 12

#define AVERAGES  2000    //numero de sistemas que simularei e posteriormente farei médias dos parametros
#define gamma   5000   // gamma e o quociente entre a massa da parede movel e a massa de uma molecula do gas
#define NL 1000
#define NR 1000
#define TL 20
#define TR 100
#define iters 37000 // numero de medidas que farei, numero de deltas nos quais medirei a Temperatura
#define deltat 0.01 //delta no qual meço a temperatura 
#define XM_ini 0.10

using namespace std;
class Sistema {
public:
  unsigned Nr;
  unsigned Nl;
  double  xM;   // posicao da parede
  double  vM;   // velocidade da parede
  double tempo;
  double dpl; // variações de momento da parede à esquerda, para depois fazer uma média num delta t
  double dpr; 
  double dpM; // variação de momento da parede móvel para posteriormente calcular uma varição total (soma das variações)  num delta t
  double *xr;   // posicoes da direita
  double *vr;   // velocidade da direita
  double *xl;   // posicoes da esquerda
  double *vl;   // velocidade da esquerda
  double eneg_l;
  double eneg_r;
  std::uniform_real_distribution <double> dist;
  std::normal_distribution<double>       gauss;
  std::mt19937 rng;


  Sistema(unsigned nl, unsigned nr) {
    Nr = nr;
    Nl = nl;
    xr = new double[Nr];
    vr = new double[Nr];
    xl = new double[Nl];
    vl = new double[Nl];

    xM = XM_ini;
    vM = 0.;
    tempo = 0;
    dpl=0;
    dpr=0;

    // Iniciar o gerador de numeros aleatorios

    std::random_device r;
    std::array<int, 624> seed_data;
    std::generate(seed_data.begin(), seed_data.end(), std::ref(r));
    std::seed_seq seq(std::begin(seed_data), std::end(seed_data));
    rng.seed(seq);

    distribuir_posicoes_r();
    distribuir_posicoes_l();
  }


  void distribuir_posicoes_r()
  {
    for(unsigned i = 0; i < Nr ; i++)
      xr[i] = xM + dist(rng)*(1. - xM );
  }

  void distribuir_posicoes_l()
  {
    for(unsigned i = 0; i < Nl ; i++)
      xl[i] = xM * dist(rng);
  }


  void distribuir_velocidades_r(double T)
  {
    for(unsigned i = 0; i < Nr ; i++)
      {
	vr[i] = sqrt(T) * gauss(rng);
	eneg_r += 0.5*vr[i]*vr[i];
      }
    
  }
  
  void distribuir_velocidades_l(double T)
  {
    for(unsigned i = 0; i < Nl ; i++)
      {
	vl[i] = sqrt(T) * gauss(rng);
	eneg_l += 0.5*vl[i]*vl[i];
      }
  }

  void reset(double Tl, double Tr)
  {
    xM = XM_ini;
    vM = 0.;
    distribuir_posicoes_r();
    distribuir_posicoes_l();
    distribuir_velocidades_l(Tl);
    distribuir_velocidades_r(Tr);
    tempo = 0;
  }
    
    
  
  void evolucao() {
    
    double tr = 1.E10, tl = 1.E10, tt;
    double vMi=vM;
    int ir, il;
    for(unsigned i = 0; i < Nr ; i++)
      {
        double t0 =  (xr[i]-xM)/(vM-vr[i]); // Colisao com o movel
        double t1 = (1-xr[i])/(vr[i]);        // Colisao com a parede
	
        if(t0 > 0 && t0 < tr)
          {
            tr = t0;
            ir = i;
          }
	
        if(t1 > 0 && t1 < tr )
          {
            tr = t1;
            ir = i;
          }
      }

    for(unsigned i = 0; i < Nl ; i++)
      {
        double t0 = (xl[i]-xM)/(vM-vl[i]); // Colisao com o movel
        double t1 = (-1*xl[i])/(vl[i]); // Colisao com a parede estatica

        if(t0 > 0 && t0 < tl)
          {
            tl = t0;
            il = i;
          }
	
        if(t1 > 0 && t1 < tl )
          {
            tl = t1;
            il = i;
          }
      }

    // Evoluir
    tt = ( tr < tl ? tr : tl )*0.99999999 ;
    tempo += tt;


    xM += vM*tt;
    for(unsigned i = 0; i < Nr ; i++)
      xr[i] += vr[i]*tt;
    for(unsigned i = 0; i < Nl ; i++)
      xl[i] += vl[i]*tt;
    
    // Colidir
    
    if(tr < tl ) //coilidiu no lado direito
      {
        // OU colide na parede em 1 ou na parede movel
        if( fabs(xr[ir] - 1.) <  0.01)
	  {
	    // colidiu com a parede
	    vr[ir] = - vr[ir];
	  
	    dpl=0;
	    dpr=fabs(2*vr[ir]); //variação de momento da parede estatica direita
	    dpM=0;
	  }
	
        else
          {
	    double vk = vr[ir];
	    vr[ir] = (2 * vM * gamma + vr[ir]*(1-gamma))/(gamma + 1);
            eneg_r += 0.5*(pow(vr[ir],2) - vk*vk);
            vM = (2 * vk  + vM*(gamma - 1))/(gamma + 1);
	    
	    dpl=0;
	    dpr=0;
	    dpM=gamma*(vM-vMi);
          }
      }
    else    //colidiu no lado esquerdo
      {
        if( fabs(xl[il] - 0.) <  0.01)	  // colidiu com a parede estatica esquerda
	  {
	    vl[il] = - vl[il];
	    dpl=fabs(2*vl[il]);
	    dpr=0;
	    dpM=0;
	  }
	
        else
          {
	    double vk = vl[il];	
	    vl[il] = (2 * vM*gamma + vl[il] *(1-gamma))/(gamma + 1);
	    eneg_l += 0.5*(pow(vl[il],2) - vk*vk);
	    vM = (2 * vk  + vM*(gamma - 1))/(gamma + 1);
	    dpl=0;
	    dpr=0;
	    dpM=gamma*(vM-vMi);
          }
      }
    
  }


  double energiaR()
  {
    double energia = 0.;
    for(unsigned i = 0; i < Nr ; i++)
      energia += vr[i]*vr[i];
    return energia/Nr;
  }

  double energiaL()
  {
    double energia = 0.;
    for(unsigned i = 0; i < Nl ; i++)
     energia += vl[i]*vl[i];
    return energia/Nl;
  }

  double tempof()
  {
    return tempo;
  }

   double posicaoM()
  {
    return xM;
  }

  double delta_mom_r()
  {
    return dpr;
  }

   double delta_mom_l()
  {
    return dpl;
  }

  double delta_mom_M()
  {
    return dpM;
  }

};  //tenho de por um ; after class definition


int main() {

  unsigned medidas = AVERAGES;         // numero de vezes que farei a mesma medida
  double *tempos   = new double[iters];   // array do tempos Temperaturas
  double *temposP   = new double[iters];   // array do tempos Pressões
  double *Tl       = new double[iters];   // Temperaturas esquerda
  double *Tr       = new double[iters];   // Temperaturas direita
  double *Pl       = new double[iters];   // Pressão esquerda
  double *Pr       = new double[iters];   // Pressão direita
  double *Fm       = new double[iters];   //Força na parede
  double *Xm       = new double[iters];   //Posiçao da parede

  
  for(unsigned i = 1; i < iters ; i++)
    { 
      tempos[i] = i*deltat;
      temposP[i]=i*deltat-0.5*deltat;   //o tempo ao qual atribuo o valor de P vai ser o meio dos dois tempos[i] que calculo
      Tl[i] = 0;
      Tr[i] = 0;
    }
  
  auto start = chrono::steady_clock::now();
  omp_set_num_threads(NUMTHREADS);

  
#pragma omp parallel
  {
    double Temp_old_r;
    double Temp_old_l;
    double X_old;
    
    double *Tr_private = new double[iters]();
    double *Tl_private = new double[iters]();
    double *Pr_private = new double[iters]();
    double *Pl_private = new double[iters]();
    double *Fm_private = new double[iters]();
    double *Xm_private = new double[iters]();
    
    double dpRt;
    dpRt=0;
    double dpLt;
    dpLt=0;
    double DeltapM;
    DeltapM=0;

    Sistema s(NL,NR);
    s.distribuir_velocidades_l(TL);
    s.distribuir_velocidades_r(TR);
    
#pragma omp for
    for(unsigned j = 0; j < medidas; j++)
      {

        s.reset(TL,TR); 
        Tr_private[0] += s.energiaR();
        Tl_private[0] += s.energiaL();
        for(unsigned i = 1;  i < iters; i++)
          {
            while (s.tempof() < tempos[i])
              {
                  
                Temp_old_r = s.energiaR();
                Temp_old_l = s.energiaL();
                X_old = s.posicaoM();
                
                dpRt+=s.delta_mom_r();
                dpLt+=s.delta_mom_l();
                DeltapM+=s.delta_mom_M();
                
                s.evolucao();
              }
            
            Tl_private[i] += Temp_old_l;
            Tr_private[i] += Temp_old_r;
            Pr_private[i] += dpRt;    /// a dividir por deltat;  divido por deltat no fim no critical, fica mais rapido
            Pl_private[i] += dpLt;
            Fm_private[i] += DeltapM;
            Xm_private[i] += X_old;
            
            dpRt=0;
            dpLt=0;
            DeltapM=0;
          }
      }
    
#pragma omp critical
    for (unsigned kk = 0; kk < iters; kk++)
      {
        Tl[kk] += Tl_private[kk];
        Tr[kk] += Tr_private[kk];
        Pl[kk] += Pl_private[kk]/deltat;
        Pr[kk] += Pr_private[kk]/deltat;
        Fm[kk] += Fm_private[kk]/deltat;
        Xm[kk] += Xm_private[kk];
        
      }
    
  }
  for(unsigned jj = 0 ; jj < iters ; jj++)
    {
      Tl[jj] = Tl[jj]/medidas;
      Tr[jj] = Tr[jj]/medidas;
      Pl[jj] = Pl[jj]/medidas;
      Pr[jj] = Pr[jj]/medidas;
      Fm[jj] = Fm[jj]/medidas;
      Xm[jj] = Xm[jj]/medidas;
    }
  
  FILE *output1;
  char name1[100];
  sprintf(name1,"Temperaturas_xmi%.2f_Nr%dNl%dTr%dTl%dgamma%diters%dmedidas%ddeltat%.3f.dat", XM_ini, NR , NL, TR, TL, gamma, iters, medidas, deltat);
  
  FILE *output2;
  char name2[100];
  sprintf(name2,"Pressões_xmi%.2f_Nr%dNl%dTr%dTl%dgamma%diters%dmedidas%ddeltat%.3f.dat", XM_ini , NR , NL, TR, TL, gamma, iters, medidas, deltat);
  
  FILE *output3;
  char name3[100];
  sprintf(name3,"XM_xmi%.2f_Nr%dNl%dTr%dTl%dgamma%diters%dmedidas%ddeltat%.3f.dat", XM_ini, NR , NL, TR, TL, gamma, iters, medidas, deltat);

  output1 = fopen(name1, "a");
  output2 = fopen(name2, "a");
  output3 = fopen(name3, "a");

  for(unsigned k = 0 ; k < iters ; k++)
    {
        
        fprintf(output1, "%f  %f  %f  \n", tempos[k] , Tl[k] , Tr[k] );
        
        fprintf(output2, "%f  %f  %f  %f  \n", temposP[k] , Pl[k] , Pr[k] , Fm[k] );
        
        fprintf(output3, "%f  %f  \n", tempos[k] , Xm[k] );
    }
    
    
  fclose(output1);
  fclose(output2);
  fclose(output3);
  
  auto end = chrono::steady_clock::now();
  
  
  
  
  FILE *output;
  char name[100];
  sprintf(name,"Tempo_xmi%.2f_Nr%dNl%dTr%dTl%dgamma%diters%dmedidas%ddeltat%.3f.dat", XM_ini , NR , NL, TR, TL, gamma, iters, medidas, deltat);
  output = fopen(name,"w");
  fprintf(output, "%e\n" , static_cast<double>(chrono::duration_cast<chrono::seconds>(end - start).count()));
  fclose(output);
  
  
  return 0;
}
