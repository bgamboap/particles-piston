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
#include <cstring>


#define NUMTHREADS 1

#define AVERAGES  1    //numero de sistemas que simularei e posteriormente farei médias dos parametros
#define gamma   5000  // gamma e o quociente entre a massa da parede movel e a massa de uma molecula do gas
#define NL 30
#define NR 0
#define TL 20
#define TR 100
#define iters 10000 // numero de medidas que farei, numero de deltas nos quais medirei a Temperatura
#define deltat 0.01 //delta no qual meço a temperatura 
#define XM_ini 0.10
#define g 1.0
#define len NL*10

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
  double tt;
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
    //rng.seed(seq);
    rng.seed(4);

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
    
    double tr = 1.E10, tl = 1.E10; //,tt
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
        //std::cout << xl[i] << " ";
        
        double a = - 0.5 * g;
        double b = - (vM-vl[i]);
        double c = (xl[i]-xM);
        
        //double t0 =  (xl[i]-xM)/(vM-vl[i]);
        double t1 = (-1*xl[i])/(vl[i]); // Colisao com a parede estatica
        if(t1 > 0 && t1 < tl )
        {
          tl = t1;
          il = i;
        }
          
        if(b*b >= 4*a*c){
          double t0_1 = (-b + std::pow(b * b - 4 * a * c,0.5)) / (2 * a); // Colisao com o movel
          double t0_2 = (-b - std::pow(b * b - 4 * a * c,0.5)) / (2 * a);
          double t0 = t0_1;
          if(t0_2 < t0 && t0_2 > 0){
            t0 = t0_2;
          }

          if(t0 > 0 && t0 < tl)
            {
              tl = t0;
              il = i;
            }
        }
      }


    // Evoluir
    tt = ( tr < tl ? tr : tl )*0.99999999 ;
    tempo += tt;


    xM += vM*tt - 0.5 * g * tt * tt;
    vM -= g * tt;
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
	
      else // colidiu com o movel
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
      if( fabs(xl[il] - 0.) <  0.001)	  // colidiu com a parede estatica esquerda
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
  /*
  double tempoit()
  {
    return tt;
  }
  */
   double posicaoM()
  {
    return xM;
  }

  void posicao_xl(double *posicoes){
    for(unsigned i=0; i < NL; i++){
      posicoes[i] = xl[i];
    }
  }

  void velocidade_vl(double *velocidades)
  {
    for(unsigned i=0; i < NL; i++)
    {
      velocidades[i] = vl[i];
    }
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

  double xl[iters][NL];

  //double **xl, **vl;
  //xl = new 

  double vl[iters][NL];
  //double *timeiter = new double[iters];

  for (unsigned i=0; i<iters; i++)
    {
      for (unsigned j=0; j<NL; j++)
      {
        xl[i][j] = 0.;
        vl[i][j] = 0.;
      }
    }
  
  
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
    //double t_old;
    
    double *Tr_private = new double[iters]();
    double *Tl_private = new double[iters]();
    double *Pr_private = new double[iters]();
    double *Pl_private = new double[iters]();
    double *Fm_private = new double[iters]();
    double *Xm_private = new double[iters]();
    //double *ts_private = new double[iters]();    

    double xl_private[iters][NL];
    double xl_temp[NL];

    double vl_private[iters][NL];
    double vl_temp[NL];
    
    for (unsigned i=0; i<iters; i++)
    {
      for (unsigned j=0; j<NL; j++)
      {
        xl_private[i][j] = 0.;
        vl_private[i][j] = 0.;
      }
    }

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
                

                dpRt+=s.delta_mom_r();
                dpLt+=s.delta_mom_l();
                DeltapM+=s.delta_mom_M();
                
                s.evolucao();
              }

            X_old = s.posicaoM();
            s.posicao_xl(xl_temp);
            s.velocidade_vl(vl_temp);
            //t_old = s.tempoit();

            for(unsigned n = 0; n < NL; n++){
              xl_private[i][n] += xl_temp[n]; 
              vl_private[i][n] += vl_temp[n]; 
            }
            
            Tl_private[i] += Temp_old_l;
            Tr_private[i] += Temp_old_r;
            Pr_private[i] += dpRt;    /// a dividir por deltat;  divido por deltat no fim no critical, fica mais rapido
            Pl_private[i] += dpLt;
            Fm_private[i] += DeltapM;
            Xm_private[i] += X_old;
            //ts_private[i] += t_old;
            
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
        //timeiter[kk] += ts_private[kk];

        for(unsigned n = 0; n<NL; n++){
          xl[kk][n] += xl_private[kk][n]/medidas;
          vl[kk][n] += vl_private[kk][n]/medidas;
        }

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
      //timeiter[jj] = timeiter[jj]/medidas;
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

  FILE *output4;
  char name4[100];
  sprintf(name4,"Posicoes_xl_xmi%.2f_Nr%dNl%dTr%dTl%dgamma%diters%dmedidas%ddeltat%.3f.dat", XM_ini, NR , NL, TR, TL, gamma, iters, medidas, deltat);

  FILE *output5;
  char name5[100];
  sprintf(name5,"Velocidades_vl_xmi%.2f_Nr%dNl%dTr%dTl%dgamma%diters%dmedidas%ddeltat%.3f.dat", XM_ini, NR , NL, TR, TL, gamma, iters, medidas, deltat);
  /*
  FILE *output6;
  char name6[100];
  sprintf(name6,"Tempos_vl_xmi%.2f_Nr%dNl%dTr%dTl%dgamma%diters%dmedidas%ddeltat%.3f.dat", XM_ini, NR , NL, TR, TL, gamma, iters, medidas, deltat);
  */
  output1 = fopen(name1, "w");
  output2 = fopen(name2, "w");
  output3 = fopen(name3, "w");
  output4 = fopen(name4, "w");
  output5 = fopen(name5, "w");
  //output6 = fopen(name6, "w");

  

  char string[len];
  char stringV[len];
  for(unsigned k = 0 ; k < iters ; k++)
    {
      
      strcpy(string, "");
      strcpy(stringV, "");
        
        fprintf(output1, "%f  %f  %f  \n", tempos[k] , Tl[k] , Tr[k] );
        
        fprintf(output2, "%f  %f  %f  %f  \n", temposP[k] , Pl[k] , Pr[k] , Fm[k] );
        
        fprintf(output3, "%f  %f  \n", tempos[k] , Xm[k] );

        //fprintf(output6, "%f  %f  \n", tempos[k] ,  );

        char strin2[10];
        char strin2V[10];

        for(unsigned n = 0; n < NL; n++){
          //posicoes
          sprintf(strin2, "%f ",xl[k][n]);
          strcat(string, strin2);

          //velocidades
          sprintf(strin2V, "%f ",vl[k][n]);
          strcat(stringV, strin2V);
        }
        
        strcat(string, "\n");
        fprintf(output4, string);

        strcat(stringV, "\n");
        fprintf(output5, stringV);
    }
    
    
  fclose(output1);
  fclose(output2);
  fclose(output3);
  fclose(output4);
  fclose(output5);
  
  auto end = chrono::steady_clock::now();
  
  
  
  
  FILE *output;
  char name[100];
  sprintf(name,"Tempo_xmi%.2f_Nr%dNl%dTr%dTl%dgamma%diters%dmedidas%ddeltat%.3f.dat", XM_ini , NR , NL, TR, TL, gamma, iters, medidas, deltat);
  output = fopen(name,"w");
  fprintf(output, "%e\n" , static_cast<double>(chrono::duration_cast<chrono::seconds>(end - start).count()));
  fclose(output);
  
  
  return 0;
}
