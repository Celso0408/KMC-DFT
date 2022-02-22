///-------------------------------------------------------------------------------------///
///----------------------------Potentiostatic KMC for Li-C -----------------------------///
///-------------------------------------------------------------------------------------///

///Developed by Maximiliano Gavilán-Arriazu

///Libraries
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <limits.h>
//#include <conio.h>    ///This may be activated to print in screen, and turned off to run in clusters
#include <time.h>

///----------------------------------------Information----------------------------------------///
/// This code runs a simulation in potentiostatic conditions for Li+/graphite system. The
/// Hamiltonian used in the present version is that taken from the work of Perassi and Leiva,
/// Electrochemistry Communications 65 (2016) 48–52. The equation utilized to calculate the
/// rate of the events is that taken from Mitchell et al. Journal of Electroanalytical
/// Chemistry 493 (2000) 68–74. The graphite substrate is constructed with a lattice-gas
/// composed of fixed sites where Li+ intercalate.
/// In the case the user wants to run simulations without understanding the code in depth,
/// the most important parameters are: dt (simulation time), mui (chemical potential),
/// Nlx (number of sites in the x axis), Nly (number of sites in the y axis), Nlam (number of
/// layers), Nbins (number of time boxes to average the outputs)
/// In the present version, after a run, you'll see 4 output files. More details below.
///-------------------------------------------------------------------------------------------///


///TIME AND POTENTIAL PARAMETERS
const int       PASTEMP=1;
const double    NPTS=PASTEMP;
const double    dt=10.0;                          ///simulation time in seconds
const double    dtequil=0.2e-5;
const double    dtprom=dt-dtequil;
const int       Nbins=60;                         ///number of time bins
const double    pasobin=dt/Nbins;                //dt/Nbins;// ///VELOCIDAD DE BARRIDO
double          mui=-0.095;                       /// chemical potential in eV
double          muf=0.05;                        //eV
double          dmu=(muf-mui)/PASTEMP;

///LATTICE CONSTANTS
const int       Nlx=24;                     ///Number of sites at the x-axis
const int       NLy=24;                     ///Number of sites at the y-axis
const int       Nly=NLy/2;
const double    ladox=Nlx;
const double    ladoy=Nly;
const int       Nlam=4;                     ///Number of layers
const int       NlamCorte=2;
const double    lam=Nlam;
const int       Npt=2*Nlx*Nly*Nlam;         ///total number of sites
const double    Nt=Npt;
const int       NxC=2*Nlx;
const int       NyC=2*Nly;
const int       NClam=NxC*NyC;
const int       NC=NClam*Nlam;     //El numero de C tiene que ser el doble de Li
const double    pi=acos(-1.0);
const double    dcx=sqrt((1.42*1.42)-(0.71*0.71));//0.5; ///cos(60)=0.5
const double    dcy=2.13;                     //0.86602540378443864676372317075294;//1.0;//sqrt((1.42*1.42)-(0.71*0.71));
const double    dvx=1.0;  		            //Distancia vecino en x
const double    dvz=3.35;
const double    Lx=ladox*2*dcx;  		        //tamaño caja en x (esta perfecto)
const double    Ly=ladoy*2*dcy;  		        //tamaño caja en y (esta perfecto)
const double    Lz=lam*dvz;
const double    Area=Lx*Lz*1.0e-16;          //area de la interfase [cm2]
const int       plano=Nlx*Nlam;
const int       Nplanos=2*Nly;
const double    rm=4.26;     		            //distancia minima de LJ
const double    rn=2.47;                     //primeros vecinos
const double    dxC=dcx*2;
const double    dyC=sqrt((1.42*1.42)-((1.0/4.0)*dxC*dxC));
const double    dC=1.42;
const double    dCm=dC*0.5;
const int       frames=0;//408-1; ///frames anteriores al ultimo

///HAMILTONIAN CONSTANTS
const double    BK=0.00008617385; 		        ///Constante de Boltzmann [eV/K]
const double    Eps=BK*296;		                //epsilon para L-J
const double    Kapp=10*BK*296;		            //kappa para Derosa
const double    rb=1.42;		                //rb para Derosa
const double    alfa=4.0;		                //alfa para Derosa
const double    gamm=-1.176*BK*296;
const double    rcortexy=10.0;
const double    rcortexy2=rcortexy+2.5;
const double    rcortez=3.36;
const int       vxy=60;
const int       vxy2=90;
const int       vz=61*(NlamCorte);
const int       vz2=91*(NlamCorte);
const int       Nvec=6;
const int       Nvecij=306;

///RATE EQUATION CONSTANTS
const double    T=296.0;                            ///Temperature [K]
const double    kT=BK*T;
const double    doskT=2*BK*T;
const double    Eadif=0.370/kT;                     ///energy barrier for diffusion [eV]
const double    Eaads=0.655/kT;                     ///energy barrier for intercalation/deintercalation [eV]
const double    kdif=1.0e13;                        ///effective vibrational frequency [s-1]
const double    kads=1.0e13;                        ///effective vibrational frequency [s-1]
const double    kdes=1.0e13;                        ///effective vibrational frequency [s-1]
const int       Nevento=8;                          ///Number of kMC events
double          factor=500.0;                        ///alfa del paper Chaterjee (2010)
const double    sc=2.0;                         ///gama del paper Chaterjee (2010)
const double    delta=0.001;                    ///delta del paper Chaterjee (2010)

///The next sentence is used to read the file "..." from the same folder of the code (if neccesary)
#define lee_info_en     "carga-24x108-SII-2.xyz"

///These define the name of output files
#define grabadif_en     "data.dat"                   ///averaged output data
#define grabaver_en     "checking-data.dat"          ///not-averaged output data
#define grabared_en     "lattice-coordinates.xyz"    ///intercalation lattice coordinates
#define grabavmd_en     "vmd.xyz"                    ///Li+ configuration coordinates
#define grabaC_en 	 "C96x96x4.xyz"		///Create Carbon atoms coordinates (not necessary)	

//#define grabaC_en       "C24x108.xyz"
FILE  *archivo, *archivo1, *archivo2, *archivo3, *archivo7;


///Random number stuffs
#define ULONGMAX	4294917238.0			// unsigned long max + 1
#define _for(a,b,c) 	for(a=b; a<(c); ++a)		// para simplificar
unsigned int seed[256];
unsigned int r;
unsigned char irr;

///DEFINITION OF VARIABLES
int i,j,k,Kk,TT,M,n3,n4,contar[6],Ocup[Npt],Na[Nbins],Nd[Nbins],cont2,numvec[Npt][Nevento],numvec2[Npt][Nvec],Nconst,Npart,ii,jj,J,I,kk,MC,n,n2,KT[PASTEMP],numven[Npt][vxy],numven2[Npt][vz],Nxy[Npt],Nz[Npt];
int NNi[Npt][Nvec][Nvecij],numven3[Npt][vxy2],numven4[Npt][vz2],Nxy2[Npt],Nz2[Npt];
double l,dx,dy,dz,CORD[Npt][3],CORDC[NC][3],Evento[Npt][Nevento],EventoN[Npt][Nevento],ti,sumaV,tempo,dT,HI,R22,lDt,Dt1,DD,RRR,HI3,Wf,W;
double Tiempo,k1,k2,k3,k4,k5,k6,cont,R,Ft,ss1,ss2,Np,H1[Npt][vxy],H2[Npt][vz],Energia[Npt][Nevento];
double TiempoBin[Nbins+1][2],En[Nbins],Timetotal,tita2[Nbins],tita[Nbins],a,TiempoMedio[Nbins],cum;
double Nint, Ndint, Ndif, En2[Nbins],Nf,BinCluster[Nbins],sumaI,It[Nbins],corriente;

///SUB-PROGRAMMS. Here we define the sub-programms that will be used in the code
void Caja();                        ///Simulation box of intercalation sites
void CajaC();                       ///Creation of the graphite substrate (this is not necessary)
void Proceso();                     ///main kMC steps
void Vecinos();                     ///neighbor lists and interaction energies
void Velocidades();                 ///to calculate all rates of a configuration
void VelocidadesAds();              ///to calculate the rates at the interface of a local environment, given a site and a cut-off
void VelocidadesDif();              ///to calculate diffusion rates of a local environment, given a site and a cut-off
void Inicializa_generador();        ///random number inicializer
void vmd();                         ///To print a configuration frame from a kMC step
void Inicial();                     ///Set the initial conditions to run kMC
void GrabaDif();                    ///Output data with instant values of a configuration, given a kMC step, this is used to follow the progress of a run
void GrabaVer();                    ///Output data from the averaged values in the "Nbins" time boxes
void carga();                       ///Read Li+ coordinates from a file (if necessary)


/////////////////////////////////////////////////////////////////////////////////////////////////////////
///	GENERA UN NUMERO REAL PSEUDO-ALEATORIO UNIFORMEMENTE DISTRIBUIDO EN EL INTERVALO [0,1)
/////////////////////////////////////////////////////////////////////////////////////////////////////////

inline double randomm(void)
{
	return (double)(r=seed[irr++]+=seed[r>>24])/ULONGMAX;
}



///-----------------------START OF THE PROGRAMM---------------------------
int main(){
int dd;

Caja();         ///This is how we "call" the sub-programms
Vecinos();
Inicial();

   for(TT=0;TT<PASTEMP;TT++){ ///This is not important by now, change the potential value if neccesary
      Ft=Tiempo=0.0;
      for(i=0;i<Nbins;i++){It[i]=En[i]=En2[i]=tita[i]=tita2[i]=TiempoBin[i][1]=0.0;}

    ///EQUILIBRATION STEP (MUTED BY NOW)
     //while(Tiempo<dtequil){for(M=0;M<Npt;M++){Proceso();}
     // printf("Muestra=%f Tiempo=%f HI=%f\n",(float)(MC),(float)(Tiempo),(float)(HI));
    // }

     ///PROCESSING STEP
     while(Tiempo<dt){          ///The code runs until "Time" reaches "dt"
            Kk++;
            Proceso();          ///This calls to perform a kMC step
            if(Kk%300000==0.0){vmd();GrabaVer();Nint=Ndint=Ndif=0.0;} ///Print coordinates and data every 3e5 (this number can be modified)
      }///END OF PROCESSING STEP

///OUTPUT-AVERAGING STEP
    Timetotal+=Tiempo;
    for(i=0;i<Nbins;i++){
        It[i]/=TiempoBin[i][1];
        En[i]/=TiempoBin[i][1];
        En2[i]/=TiempoBin[i][1];
        tita[i]/=TiempoBin[i][1];
        tita2[i]/=TiempoBin[i][1];
        if(TiempoBin[i][1]==0){En[i]=tita[i]=tita2[i]=0.0;}
        GrabaDif();
    }

  }


}///------------------------END OF THE PROGRAMM---------------------------


///----------------------------SUB-PROGRAMMS AREA---------------------------------------

void Proceso(){
    int ll,dd;
    double Rdos;

        ///First random number
    A:  R=randomm();if(R==0.0){goto A;}

        ///An event is chosen from EventoN[i][k]
        ii=kk=0;cum=0.0;
        for(i=0;i<Npt;i++){for(k=0;k<Nevento;k++){cum+=EventoN[i][k];if(cum>=R){kk=k;break;}}if(cum>=R){ii=i;break;}}if(cum<R){goto A;}

        ///Second random number
   C:   Rdos=randomm();if(Rdos==0.0){goto C;}

        ti=-log(Rdos)/sumaV;Tiempo+=ti;HI+=Energia[ii][kk];

        ///Rates calculation for the new configuration
        switch(kk){case 6:{Ocup[ii]=1;Nconst++;Nint++;}break;
                   case 7:{Ocup[ii]=0;Nconst--;Ndint++;}break;
                   default:{Ocup[ii]=0;Ocup[numvec[ii][kk]]=1;Ndif++;}break;
        }

        sumaV=sumaI=0.0;
        switch(kk){case 6:{VelocidadesAds();for(i=0;i<plano;i++){sumaI-=Evento[i][6];}}break;
                   case 7:{VelocidadesAds();for(i=0;i<plano;i++){sumaI+=Evento[i][7];}}break;
                   default:{VelocidadesDif();}break;
        }

        dd=0;
        while(dd<Nbins){
            if((Tiempo>=TiempoBin[dd][0])&&(Tiempo<TiempoBin[dd+1][0])){
                TiempoBin[dd][1]++;
                It[dd]+=sumaI;
                En[dd]+=HI;
                En2[dd]+=HI*HI;
                tita[dd]+=Nconst;
                tita2[dd]+=Nconst*Nconst;
                if(kk==6){Na[dd]++;}
                if(kk==7){Nd[dd]++;}
                dd=Nbins+1;
            }
            else{dd++;}
        }

}


void Inicial()
{
    double a;

    Tiempo=ti=sumaV=sumaI=corriente=HI=Npart=R22=0.0;

    Inicializa_generador();

    for(i=0;i<=Nbins;i++){for(j=0;j<2;j++){TiempoBin[i][j]=0.0;}}
    a=0.0;//dtequil;
    for(i=0;i<=Nbins;i++){TiempoBin[i][0]=a;a+=pasobin;}
    for(i=0;i<Nbins;i++){tita[i]=TiempoMedio[i]=0.0;Na[i]=Nd[i]=0;}
    for(i=0;i<Nbins;i++){TiempoMedio[i]=TiempoBin[i][0]+((TiempoBin[i+1][0]-TiempoBin[i][0])/2);}

    int e1;
    double HI2,aa,cubr;

    Tiempo=ti=sumaV=sumaI=corriente=Npart=R22=Nconst=HI=0.0;
    for(i=0;i<Nbins;i++){tita[i]=0.0;}
    for(i=0;i<Npt;i++){Ocup[i]=0;}
    for(i=0;i<Npt;i++){for(j=0;j<Nevento;j++){Evento[i][j]=EventoN[i][j]=Energia[i][j]=0.0;}}

    //carga();
    vmd();
    Velocidades();

    ///Creation of the titles of the output files, "//" or "///" mute parts of the codes
    (archivo =fopen (grabadif_en,"a"));fprintf(archivo,"Tiempo(e9) mu[eV] CuentaTiempo <E>[eV] <E>2[eV] <E2>[eV] Entita Na Nd dt[s] Eadif[eV] kdif[s-1] Eads[eV] kads[s-1[ <N> <N>2 <N2> factor Nlx[A] Nly[A] Nlz[A] Muestras\n");fclose(archivo);
    (archivo =fopen (grabaver_en,"a"));fprintf(archivo,"Tiempo(e10) mu[eV] tita NPart Nads Ndes Ndif En Eadif[eV] Eads[eV] factor Nlx Nly Nlz Steps\n");fclose(archivo);


    RRR=0.0;
    Kk=0;
    cont2=0;Timetotal=0.0;
}


void Vecinos(){
    int n,n1,n2,n3,n4,n5,n6,AAB,AA,BB,BBB,CC,CCB,DD,DDB;
    double nn,nn2,nn3,nn4,nn5,cc,r2,dx2,dy2;
    int Nvec=6;

    for(i=0;i<Npt;i++){for(j=0;j<Nvec;j++){numvec[i][j]=0;}}
    for(i=0;i<Npt;i++){for(j=0;j<Nvec;j++){numvec2[i][j]=-1;}}

    ///PRIMER PLANO

    for(i=0;i<plano;i++){n=0;for(j=0;j<Npt;j++){
            dx=fabs(CORD[i][0]-CORD[j][0]);dy=fabs(CORD[i][1]-CORD[j][1]);dz=fabs(CORD[i][2]-CORD[j][2]);
            if(dz>0.5*Lz){dz=Lz-dz;}if(dx>0.5*Lx){dx=Lx-dx;}//if(dy>0.5*Ly){dy=Ly-dy;}
            if(dz==0){r2=sqrt((dx*dx)+(dy*dy));
                if((r2<2.47)&&(r2>2.45)){
                    dx2=CORD[i][0]-CORD[j][0];if(fabs(dx2)>0.5*Lx){if(dx2<0){dx2=Lx-fabs(dx2);}else{dx2=-Lx+dx2;}}
                    dy2=CORD[i][1]-CORD[j][1];//if(fabs(dy2[i][j])>0.5*Ly){if(dy2[i][j]<0){dy2[i][j]=Ly-fabs(dy2[i][j]);}else{dy2[i][j]=-Ly+dy2[i][j];}}
                    //if(i==0){printf("i=%d j=%d dx=%f dy=%f\n",(int)(i),(int)(j),(float)(dx2[i][j]),(float)(dy2[i][j]));}
                                                                            if((dy2==0)&&(dx2<0)){numvec[i][0]=j;}///derecha
                                                                            if((dy2==0)&&(dx2>0)){numvec[i][1]=j;}///izquierda
                                                                            if((dy2<0)&&(dx2<0)){numvec[i][2]=j;}///arriba derecha
                                                                            if((dy2<0)&&(dx2>0)){numvec[i][3]=j;}///arriba izquierda
    //printf("i=%d j=%d 0=%f 1=%f 2=%f 3=%f 4=%f 5=%f\n",(int)(i),(int)(j),(float)(numvec[i][0]),(float)(numvec[i][1]),(float)(numvec[i][2]),(float)(numvec[i][3]),(float)(numvec[i][4]),(float)(numvec[i][5])); getch();
                }
                if((r2<4.27)&&(r2>4.25)){numvec2[i][n]=j;n++;}}}}///printf("i=%f j=%f\n",(float)(i),(float)(j));getch();}}}}


    for(i=Npt-plano;i<Npt;i++){n=0;for(j=0;j<Npt;j++){
            dx=fabs(CORD[i][0]-CORD[j][0]);dy=fabs(CORD[i][1]-CORD[j][1]);dz=fabs(CORD[i][2]-CORD[j][2]);
            if(dz>0.5*Lz){dz=Lz-dz;}if(dx>0.5*Lx){dx=Lx-dx;}//if(dy>0.5*Ly){dy=Ly-dy;}
            if(dz==0){r2=sqrt((dx*dx)+(dy*dy));
                if((r2<2.47)&&(r2>2.45)){
                    dx2=CORD[i][0]-CORD[j][0];if(fabs(dx2)>0.5*Lx){if(dx2<0){dx2=Lx-fabs(dx2);}else{dx2=-Lx+dx2;}}
                    dy2=CORD[i][1]-CORD[j][1];//if(fabs(dy2[i][j])>0.5*Ly){if(dy2[i][j]<0){dy2[i][j]=Ly-fabs(dy2[i][j]);}else{dy2[i][j]=-Ly+dy2[i][j];}}
                    //if(i==0){printf("i=%d j=%d dx=%f dy=%f\n",(int)(i),(int)(j),(float)(dx2[i][j]),(float)(dy2[i][j]));}
                                                                            if((dy2==0)&&(dx2<0)){numvec[i][0]=j;}///derecha
                                                                            if((dy2==0)&&(dx2>0)){numvec[i][1]=j;}///izquierda
                                                                            if((dy2>0)&&(dx2>0)){numvec[i][4]=j;}///abajo izquierda
                                                                            if((dy2>0)&&(dx2<0)){numvec[i][5]=j;}///abajo derecha
    //printf("i=%d j=%d 0=%f 1=%f 2=%f 3=%f 4=%f 5=%f\n",(int)(i),(int)(j),(float)(numvec[i][0]),(float)(numvec[i][1]),(float)(numvec[i][2]),(float)(numvec[i][3]),(float)(numvec[i][4]),(float)(numvec[i][5])); getch();
                }
                if((r2<4.27)&&(r2>4.25)){numvec2[i][n]=j;n++;}}}}


    ///RESTO DE LA RED (6)

    for(i=plano;i<Npt-plano;i++){n=0;for(j=0;j<Npt;j++){
            dx=fabs(CORD[i][0]-CORD[j][0]);dy=fabs(CORD[i][1]-CORD[j][1]);dz=fabs(CORD[i][2]-CORD[j][2]);
            if(dz>0.5*Lz){dz=Lz-dz;}if(dx>0.5*Lx){dx=Lx-dx;}//if(dy>0.5*Ly){dy=Ly-dy;}
            if(dz==0){r2=sqrt((dx*dx)+(dy*dy));
                if((r2<2.47)&&(r2>2.45)){
                    dx2=CORD[i][0]-CORD[j][0];if(fabs(dx2)>0.5*Lx){if(dx2<0){dx2=Lx-fabs(dx2);}else{dx2=-Lx+dx2;}}
                    dy2=CORD[i][1]-CORD[j][1];//if(fabs(dy2[i][j])>0.5*Ly){if(dy2[i][j]<0){dy2[i][j]=Ly-fabs(dy2[i][j]);}else{dy2[i][j]=-Ly+dy2[i][j];}}
                    //if(i==0){printf("i=%d j=%d dx=%f dy=%f\n",(int)(i),(int)(j),(float)(dx2[i][j]),(float)(dy2[i][j]));}
                                                                            if((dy2==0)&&(dx2<0)){numvec[i][0]=j;}///derecha
                                                                            if((dy2==0)&&(dx2>0)){numvec[i][1]=j;}///izquierda
                                                                            if((dy2<0)&&(dx2<0)){numvec[i][2]=j;}///arriba derecha
                                                                            if((dy2<0)&&(dx2>0)){numvec[i][3]=j;}///arriba izquierda
                                                                            if((dy2>0)&&(dx2>0)){numvec[i][4]=j;}///abajo izquierda
                                                                            if((dy2>0)&&(dx2<0)){numvec[i][5]=j;}///abajo derecha
    ///printf("i=%d j=%d 0=%f 1=%f 2=%f 3=%f 4=%f 5=%f\n",(int)(i),(int)(j),(float)(numvec[i][0]),(float)(numvec[i][1]),(float)(numvec[i][2]),(float)(numvec[i][3]),(float)(numvec[i][4]),(float)(numvec[i][5])); getch();
                }
                 if((r2<4.27)&&(r2>4.25)){numvec2[i][n]=j;n++;}}}}




///CALCULO ENERGIA
    double RR,rr2;

    for(i=0;i<Npt;i++){for(j=0;j<vxy;j++){numven[i][j]=0;H1[i][j]=0.0;}}
    for(i=0;i<Npt;i++){for(j=0;j<vz;j++){numven2[i][j]=0;H2[i][j]=0.0;}}
    for(i=0;i<Npt;i++){for(j=0;j<vxy2;j++){numven3[i][j]=0;}}
    for(i=0;i<Npt;i++){for(j=0;j<vz2;j++){numven4[i][j]=0;}}
    for(i=0;i<Npt;i++){Nxy[i]=Nz[i]=Nxy2[i]=Nz2[i]=0;}

    for(i=0;i<Npt;i++){
        n=n2=n3=n4=0;
        for (j=0;j<Npt;j++){
                dx=fabs(CORD[i][0]-CORD[j][0]);if(dx>0.5*Lx){dx=Lx-dx;}
                dy=fabs(CORD[i][1]-CORD[j][1]);//if(dy>0.5*Ly){dy=Ly-dy;}
                dz=fabs(CORD[i][2]-CORD[j][2]);if(dz>0.5*Lz){dz=Lz-dz;}

                RR=sqrt((dx*dx)+(dy*dy)+(dz*dz));
                rr2=sqrt((dx*dx)+(dy*dy));

                ///Lennard-Jones
                if(dz==0.0){if(i!=j){if(rr2<=rcortexy){numven[i][n]=j;
                                                       H1[i][n]=Eps*((pow((rm/rr2), 12.0))-2.0*(pow((rm/rr2), 6.0)));
                                                       n++;
                                                       Nxy[i]++;
                                                       }}}
                ///De-Rosa
                else{if(i!=j){if(rr2<=rcortexy){if(dz<rcortez){numven2[i][n2]=j;
                                                               H2[i][n2]=Kapp*pow(rb/RR, alfa);
                                                               n2++;
                                                               Nz[i]++;
                                                               }}}}

                 ///Lennard-Jones
                if(dz==0.0){if(i!=j){if(rr2<=rcortexy2){numven3[i][n3]=j;
                                                       n3++;
                                                       Nxy2[i]++;
                                                       }}}
                ///De-Rosa
                else{if(i!=j){if(rr2<=rcortexy2){if(dz<rcortez){numven4[i][n4]=j;
                                                               n4++;
                                                               Nz2[i]++;
                                                               }}}}

        } //del for j
        //if(i==0){printf("vi=%d \n",(int)(Nxy[i]+Nz[i]));getch();}
    }//del for i


    ///MATRIZ DE SITIOS COMPARTIDOS
    for(i=0;i<Npt;i++){for(j=0;j<Nvec;j++){for(k=0;k<Nvecij;k++){NNi[i][j][k]=-1;}}}

    for(i=0;i<plano;i++){
            for(j=0;j<Nvec-2;j++){
                n2=0;
                k=numvec[i][j];//printf("vj=%d \n",(int)(Nxy[k]+Nz[k]));getch();
                for(n=0;n<Nxy2[i];n++){
                        n1=numven3[i][n];
                            if((n1!=i)&&(n1!=k)){
                                NNi[i][j][n2]=n1;n2++;}
                }

                for(n=0;n<Nz2[i];n++){
                        n1=numven4[i][n];
                        if((n1!=i)&&(n1!=k)){
                                 NNi[i][j][n2]=n1;n2++;}
                }
                //printf("i=%d j=%d n2=%d\n",(int)(i),(int)(k),(int)(n2));getch();
                for(n=0;n<Nxy2[k];n++){
                    AA: n1=numven3[k][n];if(n==Nxy2[k]){break;}
                    for(n3=0;n3<n2;n3++){
                           if(NNi[i][j][n3]==n1){n++;goto AA;}
                    }
                           if((n1!=i)&&(n1!=k)){
                                    NNi[i][j][n2]=n1;n2++;}
                }//printf("i=%d j=%d n2=%d\n",(int)(i),(int)(k),(int)(n2));getch();
                for(n=0;n<Nz2[k];n++){
                    AAB: n1=numven4[k][n];if(n==Nz2[k]){break;}
                    for(n3=0;n3<n2;n3++){
                           if(NNi[i][j][n3]==n1){n++;goto AAB;}}
                           if((n1!=i)&&(n1!=k)){
                                    NNi[i][j][n2]=n1;n2++;
                }}//if((i==28)&&(k==100)){printf("i=%d j=%d n2=%d\n",(int)(i),(int)(k),(int)(n2));getch();}
                //for(n3=0;n3<n2;n3++){if((i==28)&&(k==100)){printf("i=%d j=%d n2=%d NN=%d\n",(int)(i),(int)(k),(int)(n2),(int)(NNi[i][j][n3]));getch();}}
        }///de j
    }///de i


    for(i=Npt-plano;i<Npt;i++){
            for(j=0;j<Nvec;j++){
                if((j!=2)&&(j!=3)){
                n2=0;

                k=numvec[i][j];//printf("vj=%d \n",(int)(Nxy[k]+Nz[k]));getch();
                for(n=0;n<Nxy2[i];n++){
                        n1=numven3[i][n];
                            if((n1!=i)&&(n1!=k)){
                                NNi[i][j][n2]=n1;n2++;}
                }

                for(n=0;n<Nz2[i];n++){
                        n1=numven4[i][n];
                        if((n1!=i)&&(n1!=k)){
                                 NNi[i][j][n2]=n1;n2++;}
                }
                //printf("i=%d j=%d n2=%d\n",(int)(i),(int)(k),(int)(n2));getch();
                for(n=0;n<Nxy2[k];n++){
                    BB: n1=numven3[k][n];if(n==Nxy2[k]){break;}
                    for(n3=0;n3<n2;n3++){
                           if(NNi[i][j][n3]==n1){n++;goto BB;}
                    }
                           if((n1!=i)&&(n1!=k)){
                                    NNi[i][j][n2]=n1;n2++;}
                }//printf("i=%d j=%d n2=%d\n",(int)(i),(int)(k),(int)(n2));getch();
                for(n=0;n<Nz2[k];n++){
                    BBB: n1=numven4[k][n];if(n==Nz2[k]){break;}
                    for(n3=0;n3<n2;n3++){
                           if(NNi[i][j][n3]==n1){n++;goto BBB;}}
                           if((n1!=i)&&(n1!=k)){
                                    NNi[i][j][n2]=n1;n2++;
                }}//printf("i=%d j=%d n2=%d\n",(int)(i),(int)(k),(int)(n2));getch();
            }///If j
       }///de j
    }///de i


    for(i=plano;i<Npt-plano;i++){
            for(j=0;j<Nvec;j++){
                n2=0;
                k=numvec[i][j];//
                for(n=0;n<Nxy2[i];n++){
                    n1=numven3[i][n];
                            if((n1!=i)&&(n1!=k)){
                                NNi[i][j][n2]=n1;n2++;}
                }

                for(n=0;n<Nz2[i];n++){
                        n1=numven4[i][n];
                        if((n1!=i)&&(n1!=k)){
                                 NNi[i][j][n2]=n1;n2++;}
                }

                for(n=0;n<Nxy2[k];n++){
                    CC: n1=numven3[k][n];if(n==Nxy2[k]){break;}
                    for(n3=0;n3<n2;n3++){
                           if(NNi[i][j][n3]==n1){n++;goto CC;}
                    }
                           if((n1!=i)&&(n1!=k)){
                                    NNi[i][j][n2]=n1;n2++;}
                }
                for(n=0;n<Nz2[k];n++){
                    CCB: n1=numven4[k][n];if(n==Nz2[k]){break;}
                    for(n3=0;n3<n2;n3++){
                           if(NNi[i][j][n3]==n1){n++;goto CCB;}}
                           if((n1!=i)&&(n1!=k)){
                                    NNi[i][j][n2]=n1;n2++;}
                }//if(i==2000){printf("i=%d j=%d n2=%d\n",(int)(i),(int)(k),(int)(n2));getch();}
        }///de j
    }///de i


}///FIN VECINOS


void Caja(){
    /// Coordenadas en x
    for(i=0;i<Nplanos;i++){for(j=0;j<Nlam;j++){for(k=0;k<Nlx;k++){if(i%2==0){CORD[k+j*Nlx+i*plano][0]=k*2*dcx;}else{CORD[k+j*Nlx+i*plano][0]=dcx+k*2*dcx;}}}}
    /// Coordenadas en y
    for(i=0;i<Nplanos;i++){for(j=0;j<Nlam;j++){for(k=0;k<Nlx;k++){if(i%2==0){CORD[k+j*Nlx+i*plano][1]=i*dcy;}else{CORD[k+j*Nlx+i*plano][1]=i*dcy;}}}}
    /// Coordenadas en z
    for(i=0;i<Nplanos;i++){for(j=0;j<Nlam;j++){for(k=0;k<Nlx;k++){if(i%2==0){CORD[k+j*Nlx+i*plano][2]=j*dvz;}else{CORD[k+j*Nlx+i*plano][2]=j*dvz;}}}}

    //for(i=0;i<Npt;i++){CORD[i][2]=0.0;}
    if((archivo2=fopen (grabared_en,"a"))== NULL );
        fprintf(archivo2,"%d  \n \n", (int)(Npt));  //\n graba uno por linea
        for(i=0;i<Npt;i++){fprintf(archivo2,"Li %4.5f %4.5f  %4.5f \n",  (float)(CORD[i][0]), (float)(CORD[i][1]), (float)(CORD[i][2]));}
        fclose(archivo2);
}///FIN CAJA

void CajaC(){
double ll,gg,hh,ff;
    for(i=0;i<Nlam;i++){for(j=0;j<NxC;j++){ll=gg=hh=ff=0.0;for(k=0;k<NyC;k++){CORDC[i*NClam+j*NyC+k][0]=(double)j*dcx;
                                                              if((j+1)%2!=0){if((k+1)%2!=0){CORDC[i*NClam+j*NyC+k][1]=dCm+dyC+ll*(2*dC+2*dyC);ll++;}
                                                                             else{CORDC[i*NClam+j*NyC+k][1]=dCm+dyC+dC+gg*(2*dC+2*dyC);gg++;}
                                                               }
                                                               else{if((k+1)%2!=0){CORDC[i*NClam+j*NyC+k][1]=dCm+hh*(2*dC+2*dyC);hh++;}
                                                                             else{CORDC[i*NClam+j*NyC+k][1]=dCm+(dC+2*dyC)+ff*(2*dC+2*dyC);ff++;}

                                                               }
                                                              CORDC[i*NClam+j*NyC+k][2]=0.5*dvz+(double)i*dvz;
    }}}

    if((archivo2=fopen (grabaC_en,"a"))== NULL );
    fprintf(archivo2,"%d  \n \n", (int)(NC));  //\n graba uno por linea
    for(i=0;i<NC;i++){fprintf(archivo2,"Li %4.5f %4.5f  %4.5f \n",  (float)(CORDC[i][0]), (float)(CORDC[i][1]), (float)(CORDC[i][2]));}
    fclose(archivo2);
}///FIN CAJA
/*
void Graba(){
    (archivo =fopen (graba_en,"a"));
    fprintf(archivo,"%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", (float)(Timetotal*1e9),(float)(mui),
            (float)(3*titaprom/Npt),(float)(titaprom),(float)(titaprom2),(float)(tita2prom),(float)(Ft),(float)(dt),(float)(HI),
            (float)(Eadif),(float)(Eaads),(float)(kdif),(float)(kads),(float)(Nlx),(float)(Nly), (float)(Nlam),(float)(Npt),
            (float)(PASTEMP),(float)(NMUESTRA));
    fclose(archivo);
}*/

void GrabaDif(){
    (archivo =fopen (grabadif_en,"a"));
    fprintf(archivo,"%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f \n", (float)(TiempoMedio[i]*1e9),(float)(mui),(float)(TiempoBin[i][1]),
            (float)(En[i]),(float)(En[i]*En[i]),(float)(En2[i]),(float)(3*tita[i]/(Npt)),(float)(Na[i]),(float)(Nd[i]),(float)(dt),(float)(Eadif*kT),(float)(kdif),(float)(Eaads*kT),
            (float)(kads),(float)(tita[i]),(float)(tita[i]*tita[i]),(float)(tita2[i]),(float)(factor),(float)(Lx),(float)(Ly),(float)(Lz));
    fclose(archivo);
}

void GrabaVer(){
    float NNN;
    NNN=Nconst;
    (archivo =fopen (grabaver_en,"a"));
    fprintf(archivo,"%f %f %f %f %f %f %f %f %f %f %f %d %d %d %d %d\n", (float)(Tiempo*1e10),(float)(mui),
            (float)(3*NNN/Npt),(float)(Nconst),(float)(Nint),(float)(Ndint),(float)(Ndif),(float)(HI),(float)(Eadif*kT),
            (float)(Eaads*kT),(float)(factor),(int)(Nlx),(int)(Nly), (int)(Nlam),(int)(Npt),(int)(Kk));
    fclose(archivo);
}


void Inicializa_generador(void)
{
	int i,j,s;
// inicializo generador pseudorandom: semilla = segundos desde 1970
	srand((unsigned)time(0));
	irr=1;
	_for(i,0,256) seed[i]=rand();
	r=seed[0];
	_for(i,0,70000) r=seed[irr++]+=seed[r>>24];

    /*int j;
	  irr=1;
	  for(j=0;j<256;j++) seed[j]=rand();
	  r=seed[0]; */
}

void vmd(){
    (archivo1 =fopen (grabavmd_en,"a"));fprintf(archivo1,"%d \n", (int)(Npt));fprintf(archivo1,"\n");
    for(i=0;i<Npt;i++){if(Ocup[i]==1){fprintf(archivo1,"Li %4.5f %4.5f %4.5f\n",(float)(CORD[i][0]),(float)(CORD[i][1]),(float)(CORD[i][2]));}
    else{fprintf(archivo1,"Li %4.5f %4.5f -1.0\n",(float)(Lx/2),(float)(Ly/2));}}
    fclose(archivo1);
}

void carga(){
int e1,jj;
int E;
double B,C,D;
char Li;

 if ((archivo1 =fopen (lee_info_en,"r"))== NULL);//{printf("Error al abrir el archivo\n"); exit(0);}

     for(jj=0;jj<frames;jj++){
        fscanf(archivo1,"%d",&E);
        fscanf(archivo1,"\n");
        for(ii=0; ii<Npt; ii++){
        fscanf(archivo1,"%s", &Li);
        fscanf(archivo1,"%le", &B);
        fscanf(archivo1,"%le", &C);
        fscanf(archivo1,"%le", &D);

        }
    } //del for ii
    e1=0;
    fscanf(archivo1,"%d",&E);
    fscanf(archivo1,"\n");
    for(ii=0; ii<Npt; ii++){
        fscanf(archivo1,"%s", &Li);
        fscanf(archivo1,"%le", &B);
        fscanf(archivo1,"%le", &C);
        fscanf(archivo1,"%le", &D);
        if(D<0.0){Ocup[ii]=0;}
        else{Ocup[ii]=1;e1++;}
    }

    fclose(archivo1);

    Nconst=e1;

    for(i=0;i<Npt;i++){
        if(Ocup[i]==1){
	HI+=gamm;
        for (n=0;n<Nxy[i];n++){j=numven[i][n];if(Ocup[j]==1){ HI+=0.5*H1[i][n];}}
        for (n=0;n<Nz[i];n++){j=numven2[i][n];if(Ocup[j]==1){ HI+=0.5*H2[i][n];}}
    }}
}


void Velocidades(){
    int k1,k2,k3,k4,k5,k6,n1,kk1,kk2,kk3,kk4,kk5;
    double Esitio,Ef;

for(i=0;i<Npt;i++){
    Esitio=Ef=0.0;//n1=0;
    for (n=0;n<Nxy[i];n++){j=numven[i][n];if(Ocup[j]==1){ Esitio+=H1[i][n];}}
    for (n=0;n<Nz[i];n++){j=numven2[i][n];if(Ocup[j]==1){ Esitio+=H2[i][n];}}

///Primer plano
    if((i>=0)&&(i<plano)){      ///SIN 4 y 5
        Evento[i][4]=Evento[i][5]=Energia[i][4]=Energia[i][5]=0.0;
        switch(Ocup[i]){
          case 0:{Energia[i][6]=Esitio+gamm;
                  Evento[i][6]=kads*exp(-(Eaads)-((Energia[i][6]-mui)/(doskT)));sumaV+=Evento[i][6];sumaI+=Evento[i][6];
                  Evento[i][0]=Evento[i][1]=Evento[i][2]=Evento[i][3]=Evento[i][7]=Energia[i][0]=Energia[i][1]=Energia[i][2]=Energia[i][3]=Energia[i][7]=0.0;
          }break;///CASE 0
          case 1:{Evento[i][6]=Energia[i][6]=0.0;
            Energia[i][7]=-(Esitio+gamm);
            Evento[i][7]=kads*exp(-(Eaads)-((Energia[i][7]+mui)/(doskT)));sumaV+=Evento[i][7];sumaI+=Evento[i][7];
            j=numvec[i][0];if(Ocup[j]==0){///0 A LA DERECHA, SIN EL VEC 1 PARA LA Ef
                for (n=0;n<Nxy[j];n++){k=numven[j][n];if((k!=i)&&(Ocup[k]==1)){ Ef+=H1[j][n];}}
                for (n=0;n<Nz[j];n++){k=numven2[j][n];if(Ocup[k]==1){ Ef+=H2[j][n];}}
                Energia[i][0]=Ef-Esitio;
                Evento[i][0]=kdif*exp(-(Energia[i][0]/(doskT))-(Eadif))/factor;sumaV+=Evento[i][0];}else{Evento[i][0]=Energia[i][0]=0.0;}
            Ef=0.0;
            j=numvec[i][1];if(Ocup[j]==0){///1 A LA IZDA, SIN EL VEC 0 PARA LA Ef
                for (n=0;n<Nxy[j];n++){k=numven[j][n];if((k!=i)&&(Ocup[k]==1)){ Ef+=H1[j][n];}}
                for (n=0;n<Nz[j];n++){k=numven2[j][n];if(Ocup[k]==1){ Ef+=H2[j][n];}}
                Energia[i][1]=Ef-Esitio;
                Evento[i][1]=kdif*exp((-Energia[i][1]/(doskT))-(Eadif))/factor;sumaV+=Evento[i][1];}else{Evento[i][1]=Energia[i][1]=0.0;}
            Ef=0.0;
            j=numvec[i][2];if(Ocup[j]==0){///2 ARRIBA DCHA, SIN EL 4
                for (n=0;n<Nxy[j];n++){k=numven[j][n];if((k!=i)&&(Ocup[k]==1)){ Ef+=H1[j][n];}}
                for (n=0;n<Nz[j];n++){k=numven2[j][n];if(Ocup[k]==1){ Ef+=H2[j][n];}}
                Energia[i][2]=Ef-Esitio;
                Evento[i][2]=kdif*exp((-Energia[i][2]/(doskT))-(Eadif))/factor;sumaV+=Evento[i][2];}else{Evento[i][2]=Energia[i][2]=0.0;}
            Ef=0.0;
            j=numvec[i][3];if(Ocup[j]==0){///3 ARRIBA IZDA, SIN EL 5
                for (n=0;n<Nxy[j];n++){k=numven[j][n];if((k!=i)&&(Ocup[k]==1)){ Ef+=H1[j][n];}}
                for (n=0;n<Nz[j];n++){k=numven2[j][n];if(Ocup[k]==1){ Ef+=H2[j][n];}}
                Energia[i][3]=Ef-Esitio;
                Evento[i][3]=kdif*exp((-Energia[i][3]/(doskT))-(Eadif))/factor;sumaV+=Evento[i][3];}else{Evento[i][3]=Energia[i][3]=0.0;}
           }break;///CASE 1
        }///SWITCH
    }///PRIMER PLANO

///Ultimo plano
     if((i>=Npt-plano)&&(i<Npt)){  ///SIN 2 y 3
          Evento[i][6]=Evento[i][7]=Evento[i][2]=Evento[i][3]=0.0;
          Energia[i][6]=Energia[i][7]=Energia[i][2]=Energia[i][3]=0.0;
            if(Ocup[i]==1){
            j=numvec[i][0];if(Ocup[j]==0){///0 A LA DERECHA, SIN EL VEC 1 PARA LA Ef
                for (n=0;n<Nxy[j];n++){k=numven[j][n];if((k!=i)&&(Ocup[k]==1)){ Ef+=H1[j][n];}}
                for (n=0;n<Nz[j];n++){k=numven2[j][n];if(Ocup[k]==1){ Ef+=H2[j][n];}}
                Energia[i][0]=Ef-Esitio;
                Evento[i][0]=kdif*exp((-Energia[i][0]/(doskT))-(Eadif))/factor;sumaV+=Evento[i][0];}else{Evento[i][0]=Energia[i][0]=0.0;}
            Ef=0.0;
            j=numvec[i][1];if(Ocup[j]==0){///1 A LA IZDA, SIN EL VEC 0 PARA LA Ef
                for (n=0;n<Nxy[j];n++){k=numven[j][n];if((k!=i)&&(Ocup[k]==1)){ Ef+=H1[j][n];}}
                for (n=0;n<Nz[j];n++){k=numven2[j][n];if(Ocup[k]==1){ Ef+=H2[j][n];}}
                Energia[i][1]=Ef-Esitio;
                Evento[i][1]=kdif*exp((-Energia[i][1]/(doskT))-(Eadif))/factor;sumaV+=Evento[i][1];}else{Evento[i][1]=Energia[i][1]=0.0;}
            Ef=0.0;
            j=numvec[i][4];if(Ocup[j]==0){///4 ABAJO IZDA, SIN EL 2
                for (n=0;n<Nxy[j];n++){k=numven[j][n];if((k!=i)&&(Ocup[k]==1)){ Ef+=H1[j][n];}}
                for (n=0;n<Nz[j];n++){k=numven2[j][n];if(Ocup[k]==1){ Ef+=H2[j][n];}}
                Energia[i][4]=Ef-Esitio;
                Evento[i][4]=kdif*exp((-Energia[i][4]/(doskT))-(Eadif))/factor;sumaV+=Evento[i][4];}else{Evento[i][4]=Energia[i][4]=0.0;}
            Ef=0.0;
            j=numvec[i][5];if(Ocup[j]==0){///5 ABAJO DCHA, SIN EL 3
                for (n=0;n<Nxy[j];n++){k=numven[j][n];if((k!=i)&&(Ocup[k]==1)){ Ef+=H1[j][n];}}
                for (n=0;n<Nz[j];n++){k=numven2[j][n];if(Ocup[k]==1){ Ef+=H2[j][n];}}
                Energia[i][5]=Ef-Esitio;
                Evento[i][5]=kdif*exp((-Energia[i][5]/(doskT))-(Eadif))/factor;sumaV+=Evento[i][5];}else{Evento[i][5]=Energia[i][5]=0.0;}
          }else{Evento[i][0]=Evento[i][1]=Evento[i][4]=Evento[i][5]=Energia[i][0]=Energia[i][1]=Energia[i][4]=Energia[i][5]=0.0;}
    }///DE ULTIMO PLANO


///Resto
        if((i>=plano)&&(i<Npt-plano)){
          Evento[i][6]=Evento[i][7]=Energia[i][6]=Energia[i][7]=0.0;
          if(Ocup[i]==1){
            j=numvec[i][0];if(Ocup[j]==0){
                for (n=0;n<Nxy[j];n++){k=numven[j][n];if((k!=i)&&(Ocup[k]==1)){ Ef+=H1[j][n];}}
                for (n=0;n<Nz[j];n++){k=numven2[j][n];if(Ocup[k]==1){ Ef+=H2[j][n];}}
                Energia[i][0]=Ef-Esitio;
                Evento[i][0]=kdif*exp((-Energia[i][0]/(doskT))-(Eadif))/factor;sumaV+=Evento[i][0];}else{Evento[i][0]=Energia[i][0]=0.0;}
            Ef=0.0;
            j=numvec[i][1];if(Ocup[j]==0){
                for (n=0;n<Nxy[j];n++){k=numven[j][n];if((k!=i)&&(Ocup[k]==1)){ Ef+=H1[j][n];}}
                for (n=0;n<Nz[j];n++){k=numven2[j][n];if(Ocup[k]==1){ Ef+=H2[j][n];}}
                Energia[i][1]=Ef-Esitio;
                Evento[i][1]=kdif*exp((-Energia[i][1]/(doskT))-(Eadif))/factor;sumaV+=Evento[i][1];}else{Evento[i][1]=Energia[i][1]=0.0;}
            Ef=0.0;
            j=numvec[i][2];if(Ocup[j]==0){
                for (n=0;n<Nxy[j];n++){k=numven[j][n];if((k!=i)&&(Ocup[k]==1)){ Ef+=H1[j][n];}}
                for (n=0;n<Nz[j];n++){k=numven2[j][n];if(Ocup[k]==1){ Ef+=H2[j][n];}}
                Energia[i][2]=Ef-Esitio;
                Evento[i][2]=kdif*exp((-Energia[i][2]/(doskT))-(Eadif))/factor;sumaV+=Evento[i][2];}else{Evento[i][2]=Energia[i][2]=0.0;}
            Ef=0.0;
            j=numvec[i][3];if(Ocup[j]==0){
                for (n=0;n<Nxy[j];n++){k=numven[j][n];if((k!=i)&&(Ocup[k]==1)){ Ef+=H1[j][n];}}
                for (n=0;n<Nz[j];n++){k=numven2[j][n];if(Ocup[k]==1){ Ef+=H2[j][n];}}
                Energia[i][3]=Ef-Esitio;
                Evento[i][3]=kdif*exp((-Energia[i][3]/(doskT))-(Eadif))/factor;sumaV+=Evento[i][3];}else{Evento[i][3]=Energia[i][3]=0.0;}
            Ef=0.0;
            j=numvec[i][4];if(Ocup[j]==0){
                for (n=0;n<Nxy[j];n++){k=numven[j][n];if((k!=i)&&(Ocup[k]==1)){ Ef+=H1[j][n];}}
                for (n=0;n<Nz[j];n++){k=numven2[j][n];if(Ocup[k]==1){ Ef+=H2[j][n];}}
                Energia[i][4]=Ef-Esitio;
                Evento[i][4]=kdif*exp((-Energia[i][4]/(doskT))-(Eadif))/factor;sumaV+=Evento[i][4];}else{Evento[i][4]=Energia[i][4]=0.0;}
            Ef=0.0;
            j=numvec[i][5];if(Ocup[j]==0){
                for (n=0;n<Nxy[j];n++){k=numven[j][n];if((k!=i)&&(Ocup[k]==1)){ Ef+=H1[j][n];}}
                for (n=0;n<Nz[j];n++){k=numven2[j][n];if(Ocup[k]==1){ Ef+=H2[j][n];}}
                Energia[i][5]=Ef-Esitio;
                Evento[i][5]=kdif*exp((-Energia[i][5]/(doskT))-(Eadif))/factor;sumaV+=Evento[i][5];}else{Evento[i][5]=Energia[i][5]=0.0;}
          }
          else{Evento[i][0]=Evento[i][1]=Evento[i][2]=Evento[i][3]=Evento[i][4]=Evento[i][5]=Energia[i][0]=Energia[i][1]=Energia[i][2]=Energia[i][3]=Energia[i][4]=Energia[i][5]=0.0;}
      }///DE RESTO
 }///DE FOR NPT

    ///Normalización de velocidades
    //for(i=0;i<Npt;i++){for(j=0;j<Nevento;j++){if(((kk==6)||(kk==7))&&Evento[i][j]!=0.0){printf("i=%d j%d Evento1=%f\n",(int)(i),(int)(j),(float)(Evento[i][j]));getch();}}}
    for(i=0;i<Npt;i++){for(j=0;j<Nevento;j++){EventoN[i][j]=Evento[i][j];EventoN[i][j]/=sumaV;}}
}///-------------------------------FIN VELOCIDADES--------------------------------///





void VelocidadesAds(){
int nnn2,n1;
double Esitio,Ef;
 //printf("VelAds ii=%d kk=%d\n",(int)(ii),(int)(kk));getch();
    i=ii; ///Si se acepta 6 o 7 el evento se modifica en el primer plano
    Esitio=Ef=0.0;
    n1=0;
    for (n=0;n<Nxy[i];n++){j=numven[i][n];if(Ocup[j]==1){ Esitio+=H1[i][n];}}
    for (n=0;n<Nz[i];n++){j=numven2[i][n];if(Ocup[j]==1){ Esitio+=H2[i][n];}}

///Primer plano
        Evento[i][4]=Evento[i][5]=Energia[i][4]=Energia[i][5]=0.0;
        switch(Ocup[i]){
          case 0:{Energia[i][6]=Esitio+gamm;
                  Evento[i][6]=kads*exp(-(Eaads)-((Energia[i][6]-mui)/doskT));
                  Evento[i][0]=Evento[i][1]=Evento[i][2]=Evento[i][3]=Evento[i][7]=Energia[i][0]=Energia[i][1]=Energia[i][2]=Energia[i][3]=Energia[i][7]=0.0;
          }break;///CASE 0
          case 1:{Evento[i][6]=Energia[i][6]=0.0;
            Energia[i][7]=-(Esitio+gamm);
            Evento[i][7]=kads*exp(-(Eaads)-((Energia[i][7]+mui)/(doskT)));
            j=numvec[i][0];if(Ocup[j]==0){///0 A LA DERECHA, SIN EL VEC 1 PARA LA Ef
                for (n=0;n<Nxy[j];n++){k=numven[j][n];if((k!=i)&&(Ocup[k]==1)){ Ef+=H1[j][n];}}
                for (n=0;n<Nz[j];n++){k=numven2[j][n];if(Ocup[k]==1){ Ef+=H2[j][n];}}
                Energia[i][0]=Ef-Esitio;
                Evento[i][0]=kdif*exp(-(Energia[i][0]/(doskT))-(Eadif))/factor;}else{Evento[i][0]=Energia[i][0]=0.0;}
            Ef=0.0;
            j=numvec[i][1];if(Ocup[j]==0){///1 A LA IZDA, SIN EL VEC 0 PARA LA Ef
                for (n=0;n<Nxy[j];n++){k=numven[j][n];if((k!=i)&&(Ocup[k]==1)){ Ef+=H1[j][n];}}
                for (n=0;n<Nz[j];n++){k=numven2[j][n];if(Ocup[k]==1){ Ef+=H2[j][n];}}
                Energia[i][1]=Ef-Esitio;
                Evento[i][1]=kdif*exp((-Energia[i][1]/(doskT))-(Eadif))/factor;}else{Evento[i][1]=Energia[i][1]=0.0;}
            Ef=0.0;
            j=numvec[i][2];if(Ocup[j]==0){///2 ARRIBA DCHA, SIN EL 4
                for (n=0;n<Nxy[j];n++){k=numven[j][n];if((k!=i)&&(Ocup[k]==1)){ Ef+=H1[j][n];}}
                for (n=0;n<Nz[j];n++){k=numven2[j][n];if(Ocup[k]==1){ Ef+=H2[j][n];}}
                Energia[i][2]=Ef-Esitio;
                Evento[i][2]=kdif*exp((-Energia[i][2]/(doskT))-(Eadif))/factor;}else{Evento[i][2]=Energia[i][2]=0.0;}
            Ef=0.0;
            j=numvec[i][3];if(Ocup[j]==0){///3 ARRIBA IZDA, SIN EL 5
                for (n=0;n<Nxy[j];n++){k=numven[j][n];if((k!=i)&&(Ocup[k]==1)){ Ef+=H1[j][n];}}
                for (n=0;n<Nz[j];n++){k=numven2[j][n];if(Ocup[k]==1){ Ef+=H2[j][n];}}
                Energia[i][3]=Ef-Esitio;
                Evento[i][3]=kdif*exp((-Energia[i][3]/(doskT))-(Eadif))/factor;}else{Evento[i][3]=Energia[i][3]=0.0;}
           }break;
        }///SWITCH


///VECINOS IN-PLANE
for(nnn2=0;nnn2<Nxy2[ii];nnn2++){
    i=numven3[ii][nnn2];
    Esitio=Ef=0.0;
    n1=0;
    for (n=0;n<Nxy[i];n++){j=numven[i][n];if(Ocup[j]==1){ Esitio+=H1[i][n];}}
    for (n=0;n<Nz[i];n++){j=numven2[i][n];if(Ocup[j]==1){ Esitio+=H2[i][n];}}

///Primer plano
    if((i>=0)&&(i<plano)){      ///SIN 4 y 5
        Evento[i][4]=Evento[i][5]=Energia[i][4]=Energia[i][5]=0.0;
        switch(Ocup[i]){
          case 0:{Energia[i][6]=Esitio+gamm;
                  Evento[i][6]=kads*exp(-(Eaads)-((Energia[i][6]-mui)/(doskT)));
                  Evento[i][0]=Evento[i][1]=Evento[i][2]=Evento[i][3]=Evento[i][7]=Energia[i][0]=Energia[i][1]=Energia[i][2]=Energia[i][3]=Energia[i][7]=0.0;
          }break;///CASE 0
          case 1:{Evento[i][6]=Energia[i][6]=0.0;
            Energia[i][7]=-(Esitio+gamm);
            Evento[i][7]=kads*exp(-(Eaads)-((Energia[i][7]+mui)/(doskT)));
            j=numvec[i][0];if(Ocup[j]==0){///0 A LA DERECHA, SIN EL VEC 1 PARA LA Ef
                for (n=0;n<Nxy[j];n++){k=numven[j][n];if((k!=i)&&(Ocup[k]==1)){ Ef+=H1[j][n];}}
                for (n=0;n<Nz[j];n++){k=numven2[j][n];if(Ocup[k]==1){ Ef+=H2[j][n];}}
                Energia[i][0]=Ef-Esitio;
                Evento[i][0]=kdif*exp(-(Energia[i][0]/(doskT))-(Eadif))/factor;}else{Evento[i][0]=Energia[i][0]=0.0;}
            Ef=0.0;
            j=numvec[i][1];if(Ocup[j]==0){///1 A LA IZDA, SIN EL VEC 0 PARA LA Ef
                for (n=0;n<Nxy[j];n++){k=numven[j][n];if((k!=i)&&(Ocup[k]==1)){ Ef+=H1[j][n];}}
                for (n=0;n<Nz[j];n++){k=numven2[j][n];if(Ocup[k]==1){ Ef+=H2[j][n];}}
                Energia[i][1]=Ef-Esitio;
                Evento[i][1]=kdif*exp((-Energia[i][1]/(doskT))-(Eadif))/factor;}else{Evento[i][1]=Energia[i][1]=0.0;}
            Ef=0.0;
            j=numvec[i][2];if(Ocup[j]==0){///2 ARRIBA DCHA, SIN EL 4
                for (n=0;n<Nxy[j];n++){k=numven[j][n];if((k!=i)&&(Ocup[k]==1)){ Ef+=H1[j][n];}}
                for (n=0;n<Nz[j];n++){k=numven2[j][n];if(Ocup[k]==1){ Ef+=H2[j][n];}}
                Energia[i][2]=Ef-Esitio;
                Evento[i][2]=kdif*exp((-Energia[i][2]/(doskT))-(Eadif))/factor;}else{Evento[i][2]=Energia[i][2]=0.0;}
            Ef=0.0;
            j=numvec[i][3];if(Ocup[j]==0){///3 ARRIBA IZDA, SIN EL 5
                for (n=0;n<Nxy[j];n++){k=numven[j][n];if((k!=i)&&(Ocup[k]==1)){ Ef+=H1[j][n];}}
                for (n=0;n<Nz[j];n++){k=numven2[j][n];if(Ocup[k]==1){ Ef+=H2[j][n];}}
                Energia[i][3]=Ef-Esitio;
                Evento[i][3]=kdif*exp((-Energia[i][3]/(doskT))-(Eadif))/factor;}else{Evento[i][3]=Energia[i][3]=0.0;}
           }break;///CASE 1
        }///SWITCH
    }///PRIMER PLANO

///Resto
        if((i>=plano)&&(i<Npt-plano)){
          Evento[i][6]=Evento[i][7]=Energia[i][6]=Energia[i][7]=0.0;
          if(Ocup[i]==1){
            n1=0;
            j=numvec[i][0];if(Ocup[j]==0){
                for (n=0;n<Nxy[j];n++){k=numven[j][n];if((k!=i)&&(Ocup[k]==1)){ Ef+=H1[j][n];}}
                for (n=0;n<Nz[j];n++){k=numven2[j][n];if(Ocup[k]==1){ Ef+=H2[j][n];}}
                Energia[i][0]=Ef-Esitio;
                Evento[i][0]=kdif*exp((-Energia[i][0]/(doskT))-(Eadif))/factor;}else{Evento[i][0]=Energia[i][0]=0.0;}
            Ef=0.0;
            j=numvec[i][1];if(Ocup[j]==0){
                for (n=0;n<Nxy[j];n++){k=numven[j][n];if((k!=i)&&(Ocup[k]==1)){ Ef+=H1[j][n];}}
                for (n=0;n<Nz[j];n++){k=numven2[j][n];if(Ocup[k]==1){ Ef+=H2[j][n];}}
                Energia[i][1]=Ef-Esitio;
                Evento[i][1]=kdif*exp((-Energia[i][1]/(doskT))-(Eadif))/factor;}else{Evento[i][1]=Energia[i][1]=0.0;}
            Ef=0.0;
            j=numvec[i][2];if(Ocup[j]==0){
                for (n=0;n<Nxy[j];n++){k=numven[j][n];if((k!=i)&&(Ocup[k]==1)){ Ef+=H1[j][n];}}
                for (n=0;n<Nz[j];n++){k=numven2[j][n];if(Ocup[k]==1){ Ef+=H2[j][n];}}
                Energia[i][2]=Ef-Esitio;
                Evento[i][2]=kdif*exp((-Energia[i][2]/(doskT))-(Eadif))/factor;}else{Evento[i][2]=Energia[i][2]=0.0;}
            Ef=0.0;
            j=numvec[i][3];if(Ocup[j]==0){
                for (n=0;n<Nxy[j];n++){k=numven[j][n];if((k!=i)&&(Ocup[k]==1)){ Ef+=H1[j][n];}}
                for (n=0;n<Nz[j];n++){k=numven2[j][n];if(Ocup[k]==1){ Ef+=H2[j][n];}}
                Energia[i][3]=Ef-Esitio;
                Evento[i][3]=kdif*exp((-Energia[i][3]/(doskT))-(Eadif))/factor;}else{Evento[i][3]=Energia[i][3]=0.0;}
            Ef=0.0;
            j=numvec[i][4];if(Ocup[j]==0){
                for (n=0;n<Nxy[j];n++){k=numven[j][n];if((k!=i)&&(Ocup[k]==1)){ Ef+=H1[j][n];}}
                for (n=0;n<Nz[j];n++){k=numven2[j][n];if(Ocup[k]==1){ Ef+=H2[j][n];}}
                Energia[i][4]=Ef-Esitio;
                Evento[i][4]=kdif*exp((-Energia[i][4]/(doskT))-(Eadif))/factor;}else{Evento[i][4]=Energia[i][4]=0.0;}
            Ef=0.0;
            j=numvec[i][5];if(Ocup[j]==0){
                for (n=0;n<Nxy[j];n++){k=numven[j][n];if((k!=i)&&(Ocup[k]==1)){ Ef+=H1[j][n];}}
                for (n=0;n<Nz[j];n++){k=numven2[j][n];if(Ocup[k]==1){ Ef+=H2[j][n];}}
                Energia[i][5]=Ef-Esitio;
                Evento[i][5]=kdif*exp((-Energia[i][5]/(doskT))-(Eadif))/factor;}else{Evento[i][5]=Energia[i][5]=0.0;}
          }
          else{Evento[i][0]=Evento[i][1]=Evento[i][2]=Evento[i][3]=Evento[i][4]=Evento[i][5]=Energia[i][0]=Energia[i][1]=Energia[i][2]=Energia[i][3]=Energia[i][4]=Energia[i][5]=0.0;}
      }///DE RESTO
}///NXY


///VECINOS OUT-PLANE
for(nnn2=0;nnn2<Nz2[ii];nnn2++){
    i=numven4[ii][nnn2];
    Esitio=Ef=0.0;
    n1=0;
    for (n=0;n<Nxy[i];n++){j=numven[i][n];if(Ocup[j]==1){ Esitio+=H1[i][n];}}
    for (n=0;n<Nz[i];n++){j=numven2[i][n];if(Ocup[j]==1){ Esitio+=H2[i][n];}}

///Primer plano
    if((i>=0)&&(i<plano)){      ///SIN 4 y 5
        Evento[i][4]=Evento[i][5]=Energia[i][4]=Energia[i][5]=0.0;
        switch(Ocup[i]){
          case 0:{Energia[i][6]=Esitio+gamm;
                  Evento[i][6]=kads*exp(-(Eaads)-((Energia[i][6]-mui)/(doskT)));
                  Evento[i][0]=Evento[i][1]=Evento[i][2]=Evento[i][3]=Evento[i][7]=Energia[i][0]=Energia[i][1]=Energia[i][2]=Energia[i][3]=Energia[i][7]=0.0;
          }break;///CASE 0
          case 1:{Evento[i][6]=Energia[i][6]=0.0;
            Energia[i][7]=-(Esitio+gamm);
            Evento[i][7]=kads*exp(-(Eaads)-((Energia[i][7]+mui)/(doskT)));
            j=numvec[i][0];if(Ocup[j]==0){///0 A LA DERECHA, SIN EL VEC 1 PARA LA Ef
                for (n=0;n<Nxy[j];n++){k=numven[j][n];if((k!=i)&&(Ocup[k]==1)){ Ef+=H1[j][n];}}
                for (n=0;n<Nz[j];n++){k=numven2[j][n];if(Ocup[k]==1){ Ef+=H2[j][n];}}
                Energia[i][0]=Ef-Esitio;
                Evento[i][0]=kdif*exp(-(Energia[i][0]/(doskT))-(Eadif))/factor;}else{Evento[i][0]=Energia[i][0]=0.0;}
            Ef=0.0;
            j=numvec[i][1];if(Ocup[j]==0){///1 A LA IZDA, SIN EL VEC 0 PARA LA Ef
                for (n=0;n<Nxy[j];n++){k=numven[j][n];if((k!=i)&&(Ocup[k]==1)){ Ef+=H1[j][n];}}
                for (n=0;n<Nz[j];n++){k=numven2[j][n];if(Ocup[k]==1){ Ef+=H2[j][n];}}
                Energia[i][1]=Ef-Esitio;
                Evento[i][1]=kdif*exp((-Energia[i][1]/(doskT))-(Eadif))/factor;}else{Evento[i][1]=Energia[i][1]=0.0;}
            Ef=0.0;
            j=numvec[i][2];if(Ocup[j]==0){///2 ARRIBA DCHA, SIN EL 4
                for (n=0;n<Nxy[j];n++){k=numven[j][n];if((k!=i)&&(Ocup[k]==1)){ Ef+=H1[j][n];}}
                for (n=0;n<Nz[j];n++){k=numven2[j][n];if(Ocup[k]==1){ Ef+=H2[j][n];}}
                Energia[i][2]=Ef-Esitio;
                Evento[i][2]=kdif*exp((-Energia[i][2]/(doskT))-(Eadif))/factor;}else{Evento[i][2]=Energia[i][2]=0.0;}
            Ef=0.0;
            j=numvec[i][3];if(Ocup[j]==0){///3 ARRIBA IZDA, SIN EL 5
                for (n=0;n<Nxy[j];n++){k=numven[j][n];if((k!=i)&&(Ocup[k]==1)){ Ef+=H1[j][n];}}
                for (n=0;n<Nz[j];n++){k=numven2[j][n];if(Ocup[k]==1){ Ef+=H2[j][n];}}
                Energia[i][3]=Ef-Esitio;
                Evento[i][3]=kdif*exp((-Energia[i][3]/(doskT))-(Eadif))/factor;}else{Evento[i][3]=Energia[i][3]=0.0;}
           }break;///CASE 1
        }///SWITCH
    }///PRIMER PLANO

///Resto
        if((i>=plano)&&(i<Npt-plano)){
          Evento[i][6]=Evento[i][7]=Energia[i][6]=Energia[i][7]=0.0;
          if(Ocup[i]==1){
            n1=0;
            j=numvec[i][0];if(Ocup[j]==0){
                for (n=0;n<Nxy[j];n++){k=numven[j][n];if((k!=i)&&(Ocup[k]==1)){ Ef+=H1[j][n];}}
                for (n=0;n<Nz[j];n++){k=numven2[j][n];if(Ocup[k]==1){ Ef+=H2[j][n];}}
                Energia[i][0]=Ef-Esitio;
                Evento[i][0]=kdif*exp((-Energia[i][0]/(doskT))-(Eadif))/factor;}else{Evento[i][0]=Energia[i][0]=0.0;}
            Ef=0.0;
            j=numvec[i][1];if(Ocup[j]==0){
                for (n=0;n<Nxy[j];n++){k=numven[j][n];if((k!=i)&&(Ocup[k]==1)){ Ef+=H1[j][n];}}
                for (n=0;n<Nz[j];n++){k=numven2[j][n];if(Ocup[k]==1){ Ef+=H2[j][n];}}
                Energia[i][1]=Ef-Esitio;
                Evento[i][1]=kdif*exp((-Energia[i][1]/(doskT))-(Eadif))/factor;}else{Evento[i][1]=Energia[i][1]=0.0;}
            Ef=0.0;
            j=numvec[i][2];if(Ocup[j]==0){
                for (n=0;n<Nxy[j];n++){k=numven[j][n];if((k!=i)&&(Ocup[k]==1)){ Ef+=H1[j][n];}}
                for (n=0;n<Nz[j];n++){k=numven2[j][n];if(Ocup[k]==1){ Ef+=H2[j][n];}}
                Energia[i][2]=Ef-Esitio;
                Evento[i][2]=kdif*exp((-Energia[i][2]/(doskT))-(Eadif))/factor;}else{Evento[i][2]=Energia[i][2]=0.0;}
            Ef=0.0;
            j=numvec[i][3];if(Ocup[j]==0){
                for (n=0;n<Nxy[j];n++){k=numven[j][n];if((k!=i)&&(Ocup[k]==1)){ Ef+=H1[j][n];}}
                for (n=0;n<Nz[j];n++){k=numven2[j][n];if(Ocup[k]==1){ Ef+=H2[j][n];}}
                Energia[i][3]=Ef-Esitio;
                Evento[i][3]=kdif*exp((-Energia[i][3]/(doskT))-(Eadif))/factor;}else{Evento[i][3]=Energia[i][3]=0.0;}
            Ef=0.0;
            j=numvec[i][4];if(Ocup[j]==0){
                for (n=0;n<Nxy[j];n++){k=numven[j][n];if((k!=i)&&(Ocup[k]==1)){ Ef+=H1[j][n];}}
                for (n=0;n<Nz[j];n++){k=numven2[j][n];if(Ocup[k]==1){ Ef+=H2[j][n];}}
                Energia[i][4]=Ef-Esitio;
                Evento[i][4]=kdif*exp((-Energia[i][4]/(doskT))-(Eadif))/factor;}else{Evento[i][4]=Energia[i][4]=0.0;}
            Ef=0.0;
            j=numvec[i][5];if(Ocup[j]==0){
                for (n=0;n<Nxy[j];n++){k=numven[j][n];if((k!=i)&&(Ocup[k]==1)){ Ef+=H1[j][n];}}
                for (n=0;n<Nz[j];n++){k=numven2[j][n];if(Ocup[k]==1){ Ef+=H2[j][n];}}
                Energia[i][5]=Ef-Esitio;
                Evento[i][5]=kdif*exp((-Energia[i][5]/(doskT))-(Eadif))/factor;}else{Evento[i][5]=Energia[i][5]=0.0;}
          }
          else{Evento[i][0]=Evento[i][1]=Evento[i][2]=Evento[i][3]=Evento[i][4]=Evento[i][5]=Energia[i][0]=Energia[i][1]=Energia[i][2]=Energia[i][3]=Energia[i][4]=Energia[i][5]=0.0;}
      }///DE RESTO
}///NZ


    ///Normalización de velocidades
    for(i=0;i<Npt;i++){for(j=0;j<Nevento;j++){sumaV+=Evento[i][j];}}//if(Evento[i][j]!=0.0){printf("i=%d j%d Evento2=%f\n",(int)(i),(int)(j),(float)(Evento[i][j]));getch();}}}
    for(i=0;i<Npt;i++){for(j=0;j<Nevento;j++){EventoN[i][j]=Evento[i][j];EventoN[i][j]/=sumaV;}}

}///-------------------------------FIN VELOCIDADES--------------------------------///




void VelocidadesDif(){
int nnn3,nnn,n1;
double Esitio,Ef;

//printf("VelDif ii=%d kk=%d\n",(int)(ii),(int)(kk));getch();
///ACTUALIZO LAS VELOCIDADES DEL SITIO DE DONDE SE FUE LA PARTICULA
    i=ii;
    Esitio=Ef=0.0;
    for (n=0;n<Nxy[i];n++){j=numven[i][n];if(Ocup[j]==1){ Esitio+=H1[i][n];}}
    for (n=0;n<Nz[i];n++){j=numven2[i][n];if(Ocup[j]==1){ Esitio+=H2[i][n];}}

///Primer plano
    if((i>=0)&&(i<plano)){      ///SIN 4 y 5
        Energia[i][6]=Esitio+gamm;
        Evento[i][6]=kads*exp(-(Eaads)-((Energia[i][6]-mui)/(doskT)));
        Evento[i][4]=Evento[i][5]=Energia[i][4]=Energia[i][5]=Evento[i][0]=Evento[i][1]=Evento[i][2]=Evento[i][3]=Evento[i][7]=Energia[i][0]=Energia[i][1]=Energia[i][2]=Energia[i][3]=Energia[i][7]=0.0;
    }///PRIMER PLANO
    else{
         Evento[i][6]=Evento[i][7]=Energia[i][6]=Energia[i][7]=Evento[i][0]=Energia[i][0]=Evento[i][1]=Energia[i][1]=Evento[i][2]=Energia[i][2]=Evento[i][3]=Energia[i][3]=Evento[i][4]=Energia[i][4]=Evento[i][5]=Energia[i][5]=0.0;
    }


///ACTUALIZO VELOCIDADES DEL SITIO AL QUE FUE A PARAR LA PARTICULA
    i=numvec[ii][kk]; ///Si se acepta Difusion sabemos que no hay Li en el sitio ii, pero si en numvec[ii][kk]
    Esitio=Ef=0.0;
    for (n=0;n<Nxy[i];n++){j=numven[i][n];if(Ocup[j]==1){ Esitio+=H1[i][n];}}
    for (n=0;n<Nz[i];n++){j=numven2[i][n];if(Ocup[j]==1){ Esitio+=H2[i][n];}}

///Primer plano
    if((i>=0)&&(i<plano)){      ///SIN 4 y 5
            Evento[i][6]=Energia[i][6]=Evento[i][4]=Evento[i][5]=Energia[i][4]=Energia[i][5]=0.0;
            Energia[i][7]=-(Esitio+gamm);
            Evento[i][7]=kads*exp(-(Eaads)-((Energia[i][7]+mui)/(doskT)));
            j=numvec[i][0];if(Ocup[j]==0){///0 A LA DERECHA, SIN EL VEC 1 PARA LA Ef
                for (n=0;n<Nxy[j];n++){k=numven[j][n];if((k!=i)&&(Ocup[k]==1)){ Ef+=H1[j][n];}}
                for (n=0;n<Nz[j];n++){k=numven2[j][n];if(Ocup[k]==1){ Ef+=H2[j][n];}}
                Energia[i][0]=Ef-Esitio;
                Evento[i][0]=kdif*exp(-(Energia[i][0]/(doskT))-(Eadif))/factor;}else{Evento[i][0]=Energia[i][0]=0.0;}
            Ef=0.0;
            j=numvec[i][1];if(Ocup[j]==0){///1 A LA IZDA, SIN EL VEC 0 PARA LA Ef
                for (n=0;n<Nxy[j];n++){k=numven[j][n];if((k!=i)&&(Ocup[k]==1)){ Ef+=H1[j][n];}}
                for (n=0;n<Nz[j];n++){k=numven2[j][n];if(Ocup[k]==1){ Ef+=H2[j][n];}}
                Energia[i][1]=Ef-Esitio;
                Evento[i][1]=kdif*exp((-Energia[i][1]/(doskT))-(Eadif))/factor;}else{Evento[i][1]=Energia[i][1]=0.0;}
            Ef=0.0;
            j=numvec[i][2];if(Ocup[j]==0){///2 ARRIBA DCHA, SIN EL 4
                for (n=0;n<Nxy[j];n++){k=numven[j][n];if((k!=i)&&(Ocup[k]==1)){ Ef+=H1[j][n];}}
                for (n=0;n<Nz[j];n++){k=numven2[j][n];if(Ocup[k]==1){ Ef+=H2[j][n];}}
                Energia[i][2]=Ef-Esitio;
                Evento[i][2]=kdif*exp((-Energia[i][2]/(doskT))-(Eadif))/factor;}else{Evento[i][2]=Energia[i][2]=0.0;}
            Ef=0.0;
            j=numvec[i][3];if(Ocup[j]==0){///3 ARRIBA IZDA, SIN EL 5
                for (n=0;n<Nxy[j];n++){k=numven[j][n];if((k!=i)&&(Ocup[k]==1)){ Ef+=H1[j][n];}}
                for (n=0;n<Nz[j];n++){k=numven2[j][n];if(Ocup[k]==1){ Ef+=H2[j][n];}}
                Energia[i][3]=Ef-Esitio;
                Evento[i][3]=kdif*exp((-Energia[i][3]/(doskT))-(Eadif))/factor;}else{Evento[i][3]=Energia[i][3]=0.0;}
    }///PRIMER PLANO

///Ultimo plano
     if((i>=Npt-plano)&&(i<Npt)){  ///SIN 2 y 3
          Evento[i][6]=Evento[i][7]=Evento[i][2]=Evento[i][3]=0.0;
          Energia[i][6]=Energia[i][7]=Energia[i][2]=Energia[i][3]=0.0;
            j=numvec[i][0];if(Ocup[j]==0){///0 A LA DERECHA, SIN EL VEC 1 PARA LA Ef
                for (n=0;n<Nxy[j];n++){k=numven[j][n];if((k!=i)&&(Ocup[k]==1)){ Ef+=H1[j][n];}}
                for (n=0;n<Nz[j];n++){k=numven2[j][n];if(Ocup[k]==1){ Ef+=H2[j][n];}}
                Energia[i][0]=Ef-Esitio;
                Evento[i][0]=kdif*exp((-Energia[i][0]/(doskT))-(Eadif))/factor;}else{Evento[i][0]=Energia[i][0]=0.0;}
            Ef=0.0;
            j=numvec[i][1];if(Ocup[j]==0){///1 A LA IZDA, SIN EL VEC 0 PARA LA Ef
                for (n=0;n<Nxy[j];n++){k=numven[j][n];if((k!=i)&&(Ocup[k]==1)){ Ef+=H1[j][n];}}
                for (n=0;n<Nz[j];n++){k=numven2[j][n];if(Ocup[k]==1){ Ef+=H2[j][n];}}
                Energia[i][1]=Ef-Esitio;
                Evento[i][1]=kdif*exp((-Energia[i][1]/(doskT))-(Eadif))/factor;}else{Evento[i][1]=Energia[i][1]=0.0;}
            Ef=0.0;
            j=numvec[i][4];if(Ocup[j]==0){///4 ABAJO IZDA, SIN EL 2
                for (n=0;n<Nxy[j];n++){k=numven[j][n];if((k!=i)&&(Ocup[k]==1)){ Ef+=H1[j][n];}}
                for (n=0;n<Nz[j];n++){k=numven2[j][n];if(Ocup[k]==1){ Ef+=H2[j][n];}}
                Energia[i][4]=Ef-Esitio;
                Evento[i][4]=kdif*exp((-Energia[i][4]/(doskT))-(Eadif))/factor;}else{Evento[i][4]=Energia[i][4]=0.0;}
            Ef=0.0;
            j=numvec[i][5];if(Ocup[j]==0){///5 ABAJO DCHA, SIN EL 3
                for (n=0;n<Nxy[j];n++){k=numven[j][n];if((k!=i)&&(Ocup[k]==1)){ Ef+=H1[j][n];}}
                for (n=0;n<Nz[j];n++){k=numven2[j][n];if(Ocup[k]==1){ Ef+=H2[j][n];}}
                Energia[i][5]=Ef-Esitio;
                Evento[i][5]=kdif*exp((-Energia[i][5]/(doskT))-(Eadif))/factor;}else{Evento[i][5]=Energia[i][5]=0.0;}
    }///DE ULTIMO PLANO


///Resto
        if((i>=plano)&&(i<Npt-plano)){
          Evento[i][6]=Evento[i][7]=Energia[i][6]=Energia[i][7]=0.0;
            n1=0;
            j=numvec[i][0];if(Ocup[j]==0){
                for (n=0;n<Nxy[j];n++){k=numven[j][n];if((k!=i)&&(Ocup[k]==1)){ Ef+=H1[j][n];}}
                for (n=0;n<Nz[j];n++){k=numven2[j][n];if(Ocup[k]==1){ Ef+=H2[j][n];}}
                Energia[i][0]=Ef-Esitio;
                Evento[i][0]=kdif*exp((-Energia[i][0]/(doskT))-(Eadif))/factor;}else{Evento[i][0]=Energia[i][0]=0.0;}
            Ef=0.0;
            j=numvec[i][1];if(Ocup[j]==0){
                for (n=0;n<Nxy[j];n++){k=numven[j][n];if((k!=i)&&(Ocup[k]==1)){ Ef+=H1[j][n];}}
                for (n=0;n<Nz[j];n++){k=numven2[j][n];if(Ocup[k]==1){ Ef+=H2[j][n];}}
                Energia[i][1]=Ef-Esitio;
                Evento[i][1]=kdif*exp((-Energia[i][1]/(doskT))-(Eadif))/factor;}else{Evento[i][1]=Energia[i][1]=0.0;}
            Ef=0.0;
            j=numvec[i][2];if(Ocup[j]==0){
                for (n=0;n<Nxy[j];n++){k=numven[j][n];if((k!=i)&&(Ocup[k]==1)){ Ef+=H1[j][n];}}
                for (n=0;n<Nz[j];n++){k=numven2[j][n];if(Ocup[k]==1){ Ef+=H2[j][n];}}
                Energia[i][2]=Ef-Esitio;
                Evento[i][2]=kdif*exp((-Energia[i][2]/(doskT))-(Eadif))/factor;}else{Evento[i][2]=Energia[i][2]=0.0;}
            Ef=0.0;
            j=numvec[i][3];if(Ocup[j]==0){
                for (n=0;n<Nxy[j];n++){k=numven[j][n];if((k!=i)&&(Ocup[k]==1)){ Ef+=H1[j][n];}}
                for (n=0;n<Nz[j];n++){k=numven2[j][n];if(Ocup[k]==1){ Ef+=H2[j][n];}}
                Energia[i][3]=Ef-Esitio;
                Evento[i][3]=kdif*exp((-Energia[i][3]/(doskT))-(Eadif))/factor;}else{Evento[i][3]=Energia[i][3]=0.0;}
            Ef=0.0;
            j=numvec[i][4];if(Ocup[j]==0){
                for (n=0;n<Nxy[j];n++){k=numven[j][n];if((k!=i)&&(Ocup[k]==1)){ Ef+=H1[j][n];}}
                for (n=0;n<Nz[j];n++){k=numven2[j][n];if(Ocup[k]==1){ Ef+=H2[j][n];}}
                Energia[i][4]=Ef-Esitio;
                Evento[i][4]=kdif*exp((-Energia[i][4]/(doskT))-(Eadif))/factor;}else{Evento[i][4]=Energia[i][4]=0.0;}
            Ef=0.0;
            j=numvec[i][5];if(Ocup[j]==0){
                for (n=0;n<Nxy[j];n++){k=numven[j][n];if((k!=i)&&(Ocup[k]==1)){ Ef+=H1[j][n];}}
                for (n=0;n<Nz[j];n++){k=numven2[j][n];if(Ocup[k]==1){ Ef+=H2[j][n];}}
                Energia[i][5]=Ef-Esitio;
                Evento[i][5]=kdif*exp((-Energia[i][5]/(doskT))-(Eadif))/factor;}else{Evento[i][5]=Energia[i][5]=0.0;}
      }///DE RESTO



 ///ACTUALIZO VELOCIDADES DE LOS VECINOS DE LOS DOS SITIOS
 nnn=0;
 for(nnn3=0;nnn3<Nvecij;nnn3++){
   if(NNi[ii][kk][nnn3]==-1){break;}
    i=NNi[ii][kk][nnn3]; ///Sitios compartidos por ii y numvec[ii][kk]

   nnn++;
    Esitio=Ef=0.0;

    for (n=0;n<Nxy[i];n++){j=numven[i][n];if(Ocup[j]==1){ Esitio+=H1[i][n];}}
    for (n=0;n<Nz[i];n++){j=numven2[i][n];if(Ocup[j]==1){ Esitio+=H2[i][n];}}

///Primer plano
    if((i>=0)&&(i<plano)){      ///SIN 4 y 5
            // printf("ii=%d jj=%d kk=%d i=%d n1=%d\n",(int)(ii),(int)(numvec[ii][kk]),(int)(kk),(int)(NNi[ii][kk][nnn]),(int)(n1));getch();
   // n1++;
        Evento[i][4]=Evento[i][5]=Energia[i][4]=Energia[i][5]=0.0;
        switch(Ocup[i]){
          case 0:{Energia[i][6]=Esitio+gamm;
                  Evento[i][6]=kads*exp(-(Eaads)-((Energia[i][6]-mui)/(doskT)));
                  Evento[i][0]=Evento[i][1]=Evento[i][2]=Evento[i][3]=Evento[i][7]=Energia[i][0]=Energia[i][1]=Energia[i][2]=Energia[i][3]=Energia[i][7]=0.0;
          }break;///CASE 0
          case 1:{Evento[i][6]=Energia[i][6]=0.0;
            Energia[i][7]=-(Esitio+gamm);
            Evento[i][7]=kads*exp(-(Eaads)-((Energia[i][7]+mui)/(doskT)));
            j=numvec[i][0];if(Ocup[j]==0){///0 A LA DERECHA, SIN EL VEC 1 PARA LA Ef
                for (n=0;n<Nxy[j];n++){k=numven[j][n];if((k!=i)&&(Ocup[k]==1)){ Ef+=H1[j][n];}}
                for (n=0;n<Nz[j];n++){k=numven2[j][n];if(Ocup[k]==1){ Ef+=H2[j][n];}}
                Energia[i][0]=Ef-Esitio;
                Evento[i][0]=kdif*exp(-(Energia[i][0]/(doskT))-(Eadif))/factor;}else{Evento[i][0]=Energia[i][0]=0.0;}
            Ef=0.0;
            j=numvec[i][1];if(Ocup[j]==0){///1 A LA IZDA, SIN EL VEC 0 PARA LA Ef
                for (n=0;n<Nxy[j];n++){k=numven[j][n];if((k!=i)&&(Ocup[k]==1)){ Ef+=H1[j][n];}}
                for (n=0;n<Nz[j];n++){k=numven2[j][n];if(Ocup[k]==1){ Ef+=H2[j][n];}}
                Energia[i][1]=Ef-Esitio;
                Evento[i][1]=kdif*exp((-Energia[i][1]/(doskT))-(Eadif))/factor;}else{Evento[i][1]=Energia[i][1]=0.0;}
            Ef=0.0;
            j=numvec[i][2];if(Ocup[j]==0){///2 ARRIBA DCHA, SIN EL 4
                for (n=0;n<Nxy[j];n++){k=numven[j][n];if((k!=i)&&(Ocup[k]==1)){ Ef+=H1[j][n];}}
                for (n=0;n<Nz[j];n++){k=numven2[j][n];if(Ocup[k]==1){ Ef+=H2[j][n];}}
                Energia[i][2]=Ef-Esitio;
                Evento[i][2]=kdif*exp((-Energia[i][2]/(doskT))-(Eadif))/factor;}else{Evento[i][2]=Energia[i][2]=0.0;}
            Ef=0.0;
            j=numvec[i][3];if(Ocup[j]==0){///3 ARRIBA IZDA, SIN EL 5
                for (n=0;n<Nxy[j];n++){k=numven[j][n];if((k!=i)&&(Ocup[k]==1)){ Ef+=H1[j][n];}}
                for (n=0;n<Nz[j];n++){k=numven2[j][n];if(Ocup[k]==1){ Ef+=H2[j][n];}}
                Energia[i][3]=Ef-Esitio;
                Evento[i][3]=kdif*exp((-Energia[i][3]/(doskT))-(Eadif))/factor;}else{Evento[i][3]=Energia[i][3]=0.0;}
           }break;///CASE 1
        }///SWITCH
    }///PRIMER PLANO

///Ultimo plano
     if((i>=Npt-plano)&&(i<Npt)){  ///SIN 2 y 3
          Evento[i][6]=Evento[i][7]=Evento[i][2]=Evento[i][3]=0.0;
          Energia[i][6]=Energia[i][7]=Energia[i][2]=Energia[i][3]=0.0;
            if(Ocup[i]==1){
            j=numvec[i][0];if(Ocup[j]==0){///0 A LA DERECHA, SIN EL VEC 1 PARA LA Ef
                for (n=0;n<Nxy[j];n++){k=numven[j][n];if((k!=i)&&(Ocup[k]==1)){ Ef+=H1[j][n];}}
                for (n=0;n<Nz[j];n++){k=numven2[j][n];if(Ocup[k]==1){ Ef+=H2[j][n];}}
                Energia[i][0]=Ef-Esitio;
                Evento[i][0]=kdif*exp((-Energia[i][0]/(doskT))-(Eadif))/factor;}else{Evento[i][0]=Energia[i][0]=0.0;}
            Ef=0.0;
            j=numvec[i][1];if(Ocup[j]==0){///1 A LA IZDA, SIN EL VEC 0 PARA LA Ef
                for (n=0;n<Nxy[j];n++){k=numven[j][n];if((k!=i)&&(Ocup[k]==1)){ Ef+=H1[j][n];}}
                for (n=0;n<Nz[j];n++){k=numven2[j][n];if(Ocup[k]==1){ Ef+=H2[j][n];}}
                Energia[i][1]=Ef-Esitio;
                Evento[i][1]=kdif*exp((-Energia[i][1]/(doskT))-(Eadif))/factor;}else{Evento[i][1]=Energia[i][1]=0.0;}
            Ef=0.0;
            j=numvec[i][4];if(Ocup[j]==0){///4 ABAJO IZDA, SIN EL 2
                for (n=0;n<Nxy[j];n++){k=numven[j][n];if((k!=i)&&(Ocup[k]==1)){ Ef+=H1[j][n];}}
                for (n=0;n<Nz[j];n++){k=numven2[j][n];if(Ocup[k]==1){ Ef+=H2[j][n];}}
                Energia[i][4]=Ef-Esitio;
                Evento[i][4]=kdif*exp((-Energia[i][4]/(doskT))-(Eadif))/factor;}else{Evento[i][4]=Energia[i][4]=0.0;}
            Ef=0.0;
            j=numvec[i][5];if(Ocup[j]==0){///5 ABAJO DCHA, SIN EL 3
                for (n=0;n<Nxy[j];n++){k=numven[j][n];if((k!=i)&&(Ocup[k]==1)){ Ef+=H1[j][n];}}
                for (n=0;n<Nz[j];n++){k=numven2[j][n];if(Ocup[k]==1){ Ef+=H2[j][n];}}
                Energia[i][5]=Ef-Esitio;
                Evento[i][5]=kdif*exp((-Energia[i][5]/(doskT))-(Eadif))/factor;}else{Evento[i][5]=Energia[i][5]=0.0;}
          }else{Evento[i][0]=Evento[i][1]=Evento[i][4]=Evento[i][5]=Energia[i][0]=Energia[i][1]=Energia[i][4]=Energia[i][5]=0.0;}
    }///DE ULTIMO PLANO


///Resto
        if((i>=plano)&&(i<Npt-plano)){
          Evento[i][6]=Evento[i][7]=Energia[i][6]=Energia[i][7]=0.0;
          if(Ocup[i]==1){
            n1=0;
            j=numvec[i][0];if(Ocup[j]==0){
                for (n=0;n<Nxy[j];n++){k=numven[j][n];if((k!=i)&&(Ocup[k]==1)){ Ef+=H1[j][n];}}
                for (n=0;n<Nz[j];n++){k=numven2[j][n];if(Ocup[k]==1){ Ef+=H2[j][n];}}
                Energia[i][0]=Ef-Esitio;
                Evento[i][0]=kdif*exp((-Energia[i][0]/(doskT))-(Eadif))/factor;}else{Evento[i][0]=Energia[i][0]=0.0;}
            Ef=0.0;
            j=numvec[i][1];if(Ocup[j]==0){
                for (n=0;n<Nxy[j];n++){k=numven[j][n];if((k!=i)&&(Ocup[k]==1)){ Ef+=H1[j][n];}}
                for (n=0;n<Nz[j];n++){k=numven2[j][n];if(Ocup[k]==1){ Ef+=H2[j][n];}}
                Energia[i][1]=Ef-Esitio;
                Evento[i][1]=kdif*exp((-Energia[i][1]/(doskT))-(Eadif))/factor;}else{Evento[i][1]=Energia[i][1]=0.0;}
            Ef=0.0;
            j=numvec[i][2];if(Ocup[j]==0){
                for (n=0;n<Nxy[j];n++){k=numven[j][n];if((k!=i)&&(Ocup[k]==1)){ Ef+=H1[j][n];}}
                for (n=0;n<Nz[j];n++){k=numven2[j][n];if(Ocup[k]==1){ Ef+=H2[j][n];}}
                Energia[i][2]=Ef-Esitio;
                Evento[i][2]=kdif*exp((-Energia[i][2]/(doskT))-(Eadif))/factor;}else{Evento[i][2]=Energia[i][2]=0.0;}
            Ef=0.0;
            j=numvec[i][3];if(Ocup[j]==0){
                for (n=0;n<Nxy[j];n++){k=numven[j][n];if((k!=i)&&(Ocup[k]==1)){ Ef+=H1[j][n];}}
                for (n=0;n<Nz[j];n++){k=numven2[j][n];if(Ocup[k]==1){ Ef+=H2[j][n];}}
                Energia[i][3]=Ef-Esitio;
                Evento[i][3]=kdif*exp((-Energia[i][3]/(doskT))-(Eadif))/factor;}else{Evento[i][3]=Energia[i][3]=0.0;}
            Ef=0.0;
            j=numvec[i][4];if(Ocup[j]==0){
                for (n=0;n<Nxy[j];n++){k=numven[j][n];if((k!=i)&&(Ocup[k]==1)){ Ef+=H1[j][n];}}
                for (n=0;n<Nz[j];n++){k=numven2[j][n];if(Ocup[k]==1){ Ef+=H2[j][n];}}
                Energia[i][4]=Ef-Esitio;
                Evento[i][4]=kdif*exp((-Energia[i][4]/(doskT))-(Eadif))/factor;}else{Evento[i][4]=Energia[i][4]=0.0;}
            Ef=0.0;
            j=numvec[i][5];if(Ocup[j]==0){
                for (n=0;n<Nxy[j];n++){k=numven[j][n];if((k!=i)&&(Ocup[k]==1)){ Ef+=H1[j][n];}}
                for (n=0;n<Nz[j];n++){k=numven2[j][n];if(Ocup[k]==1){ Ef+=H2[j][n];}}
                Energia[i][5]=Ef-Esitio;
                Evento[i][5]=kdif*exp((-Energia[i][5]/(doskT))-(Eadif))/factor;}else{Evento[i][5]=Energia[i][5]=0.0;}
          }
          else{Evento[i][0]=Evento[i][1]=Evento[i][2]=Evento[i][3]=Evento[i][4]=Evento[i][5]=Energia[i][0]=Energia[i][1]=Energia[i][2]=Energia[i][3]=Energia[i][4]=Energia[i][5]=0.0;}
      }///DE RESTO
}

    ///Normalización de velocidades
    for(i=0;i<Npt;i++){for(j=0;j<Nevento;j++){sumaV+=Evento[i][j];}}
    for(i=0;i<Npt;i++){for(j=0;j<Nevento;j++){EventoN[i][j]=Evento[i][j];EventoN[i][j]/=sumaV;}}

}///-------------------------------FIN VELOCIDADES--------------------------------///
