#include <stdio.h>
#include <math.h>

/*Meta-Constantes*/
#define N 64
#define NL 4
#define Nt 10000
#define L 15.26548385700747

/*Declaracion de variables*/
double T = 100.0;
double dt = 0.001; /*en picosegundos*/
double mass = 39.948; /*en umas*/
double rx[N], ry[N], rz[N];
double vx[N], vy[N], vz[N];

/*Valores del potencial*/
double sigma = 3.4;/*Unidades de longitud en Angstroms*/
double epsilon = 99.3975903614458; /*1.65E-21 / 1.66E-23 en joules, pero no es la unidad que necesitamos, la que queremos seria dividir entre 1.66x10^-23*/

/*Constante de boltzman*/
double KB = 0.8316867469879518; /*1.3806E-23/1.66E-23;*/

/*Declaracion de funciones*/
double Upot(double r);
void calculaFij(double r, int i, int j, double Fij[], int ipx, int ipy, int ipz);
void iniciales(void);

/*Definimos la dinámica del sistema*/
int main(void)
{
    FILE * fp;
    FILE * fpos;

    int i, j, it, ipx, ipy, ipz, tao, icero;
    double Fij[3], Fi[3], Ecin, Ui, r, Epot, Ecinprom;
    double T_t, lambda, doubletao;
    double newrx[N],newry[N],newrz[N];
    double newvx[N],newvy[N],newvz[N];
    tao = 10.0;
    doubletao = (double)tao*dt;
    iniciales();
    fp = fopen ("./MD_64p_T100_tao10.txt","w");
    fpos = fopen ("./pos.txt","w");

    Ecinprom = 0.0;
    icero = 6000;
    for (it=1; it<Nt; it++)
    {
        printf("%i\n",it);
        Ecin = 0.0;
        Epot = 0.0;
        for (i=0; i<N; i++)
        {

            Fi[0] = 0.0;
            Fi[1] = 0.0;
            Fi[2] = 0.0;

            Ui = 0.0;

            for(ipx= -1; ipx <=1; ipx++)
                for(ipy= -1; ipy <=1; ipy++)
                        for(ipz= -1; ipz <=1; ipz++)
                            for (j=0; j<N; j++)
                            {
                                if (!((ipx==0 && ipy==0 && ipz==0)&&(j==i)))
                                {
                                    /*calcular fuerza*/
                                    r = sqrt((rx[j] + ipx*L - rx[i])*(rx[j] + ipx*L - rx[i]) + (ry[j] + ipy*L - ry[i])*(ry[j] + ipy*L - ry[i]) + (rz[j] + ipz*L - rz[i])*(rz[j] + ipz*L - rz[i]));
                                    calculaFij(r, i, j, Fij, ipx, ipy, ipz);
                                    Fi[0] = Fi[0] + Fij[0];
                                    Fi[1] = Fi[1] + Fij[1];
                                    Fi[2] = Fi[2] + Fij[2];
                                    /*r = sqrt((rx[j]-rx[i])*(rx[j]-rx[i]) + (ry[j]-ry[i])*(ry[j]-ry[i]) + (rz[j]-rz[i])*(rz[j]-rz[i]));*/
                                    Ui = Ui + Upot(r);
                                }
                            }
            /*Pasos necesarios para evitar la actualizacion de las variables instantáneas siendo modificadas*/
            /*Se modifica las velocidades y posiciones dependiendo de las velocidades actuales*/
            newvx[i] = vx[i] + Fi[0]/mass*dt;
            newvy[i] = vy[i] + Fi[1]/mass*dt;
            newvz[i] = vz[i] + Fi[2]/mass*dt;


            newrx[i] = rx[i] + vx[i]*dt + 0.5*Fi[0]/mass*dt*dt;
            newry[i] = ry[i] + vy[i]*dt + 0.5*Fi[1]/mass*dt*dt;
            newrz[i] = rz[i] + vz[i]*dt + 0.5*Fi[2]/mass*dt*dt;

            /*
            newrx[i] = rx[i] + vx[i]*dt;
            newry[i] = ry[i] + vy[i]*dt;
            newrz[i] = rz[i] + vz[i]*dt;*/

            /*Se calcula las velocidades actuales*/
            Ecin = Ecin + 0.5 * mass * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
            Epot = Epot + 0.5 * Ui;
        }
        /*Calculo*/
        if (it%tao==0){
            T_t = 2.0*Ecin/(3.0*N*KB);
            lambda = sqrt(1.0+ dt/doubletao*(T/T_t-1.0));
            for (i=0; i<N; i++){
                newvx[i] = newvx[i]*lambda;
                newvy[i] = newvy[i]*lambda;
                newvz[i] = newvz[i]*lambda;
            }
        }
        if (newrx[i]< 0.0) newrx[i] = newrx[i] + L ;
        if (newrx[i]> L) newrx[i] = newrx[i] - L ;
        if (newry[i]< 0.0) newry[i] = newry[i] + L ;
        if (newry[i]> L) newry[i] = newry[i] - L ;
        if (newrz[i]< 0.0) newrz[i] = newrz[i] + L ;
        if (newrz[i]> L) newrz[i] = newrz[i] - L ;

        fprintf (fp,"%i,%e,%e,%e\n",it,Ecin,Epot,Ecin+Epot);
        if (it>icero)
        {
            Ecinprom = Ecinprom + Ecin;
        }
        /*printf ("%i,%e,%e,%e\n",it,Ecin,Epot,Ecin+Epot);*/
        /*Actualizar ahora sí las variables*/
        for (i=0; i<N; i++)
        {
            rx[i] = newrx[i];
            ry[i] = newry[i];
            rz[i] = newrz[i];
            vx[i] = newvx[i];
            vy[i] = newvy[i];
            vz[i] = newvz[i];
        }
        fprintf (fpos,"%i,%e,%e,%e\n",it,rx[1],ry[1],rz[1]);
    }
    Ecinprom = Ecinprom/((double)(Nt-icero));
    printf("%e", Ecinprom);
    fclose (fp);
    fclose (fpos);
    return(1);
}

double Upot(double r)
{
    return(4.0*epsilon*(pow(sigma/r,12)-pow(sigma/r,6)));
}

void calculaFij(double r, int i, int j, double Fij[], int ipx, int ipy, int ipz)
{
    double fij;
    /*r = sqrt((rx[j]-rx[i])*(rx[j]-rx[i]) + (ry[j]-ry[i])*(ry[j]-ry[i]) + (rz[j]-rz[i])*(rz[j]-rz[i]));*/
    fij = -24.0*epsilon/(sigma*sigma)*pow(sigma/r,8)*(2.0*pow(sigma/r,6)-1.0);
    Fij[0] = fij*(rx[j] + L*ipx - rx[i]);
    Fij[1] = fij*(ry[j] + L*ipy - ry[i]);
    Fij[2] = fij*(rz[j] + L*ipz - rz[i]);
}

void iniciales(void)
{
    int i, ix, iy, iz;
    for (i=0; i<N; i++)
    {
        vx[i] = 0.0;
        vy[i] = 0.0;
        vz[i] = 0.0;
    }
    for (iz=0; iz<NL; iz++)
    {
        for (iy=0; iy<NL; iy++)
        {
            for (ix=0; ix<NL; ix++)
            {
            rx[iz*NL*NL+iy*NL+ix] = L /(NL+1.0)*(ix+1);
            ry[iz*NL*NL+iy*NL+ix] = L /(NL+1.0)*(iy+1);
            rz[iz*NL*NL+iy*NL+ix] = L /(NL+1.0)*(iz+1);
            }
        }
    }
}
