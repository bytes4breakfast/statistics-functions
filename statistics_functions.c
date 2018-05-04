
//#define DEBUG

double gen_rand( double a, double b);
                                    /* gen_rand
                                    **Parametros:
                                    ***a: minimo posible
                                    ***b: maximo posible

                                    **gen_rand() utiliza el generador plano de C, rand(), para generar
                                    un número plano en el intervalo [a,b]*/



void histograma(double *dat, double *out, double *m, double *M, double *delta, int N_dat, int N_interval);
                                    /*histograma:

                                    **Parámetros:
                                    *** dat[]: vector con los datos de entrada a procesar.
                                    *** out[]: Salida de la función, con los valores del histograma para cada intervalo.
                                    *** m : (salida)minimo de los datos de entrada
                                    *** M : (salida)Máximo de los datos de entrada
                                    *** delta:(salida) Anchura de los intervalos del histograma.
                                    *** N_dat: Numero de datos de entrada.
                                    *** N_interval: Numero de intervalos que contendrá el histograma

                                    **histograma() calcula la altura de las N_interval columnas de un histograma,
                                    asignando automaticamente la anchura de estas en función de los parametros m, M y N_interval.
                                    */

void histograma_F(double *dat, double *out, double *m, double *M, double *delta, int N_dat, int N_interval, char *nombre_fichero);
                                    /*histograma_F:

                                    **Parámetros:
                                    *** dat[]: vector con los datos de entrada a procesar.
                                    *** out[]: Salida de la función, con los valores del histograma para cada intervalo.
                                    *** m : (salida)minimo de los datos de entrada
                                    *** M : (salida)Máximo de los datos de entrada
                                    *** delta:(salida) Anchura de los intervalos del histograma.
                                    *** N_dat: Numero de datos de entrada.
                                    *** N_interval: Numero de intervalos que contendrá el histograma
                                    *** nombre_fichero: Cadena con el nombre que tendrá el fichero generado

                                    **histograma_F() crea un fichero con las alturas de las N_interval columnas de un histograma (col3),
                                    la posicion de inicio de cada columna(col2) y el numero de cada columna(col1),
                                    asignando automaticamente la anchura de los intervalos en función de los parametros m, M y N_interval.
                                    */

void stat_med_var(double *serie, int N, double *media, double *varianza);

                                    /* stat_med_var:
                                    **Parametros:
                                    *** serie: vector de entrada con los datos.
                                    *** N: numero de datos en serie.
                                    *** media (salida): estimador de la media de la serie.
                                    *** varianza (salida): estimador de la varianza de la serie.

                                    Calcula la media y la varianza de un conjunto de datos.
                                    */

unsigned int PR_rand(void);  /*PR_rand:
			 	 
				 Llamar a srand() antes de usar la función.
					

                                 Utiliza el algoritmo de Parisi-Rapuano (Shift-register) para
                                 generar números aleatorios planos en el intervalo [0,2^32-1].
				 Utiliza rand() como fuente de entropía.
                                 */




//*******HISTOGRAMA*****************************************************************


 void histograma(double *dat, double *out, double *m, double *M, double *delta, int N_dat, int N_interval)
{
    int i, indicehist;
    double norm;

    *(m)=10000000;  //Inicializo el máximo y el mínimo a la primera componente del array que contiene los datos.
    *(M)=-10000000;

    for(i=0; i<N_interval; i++) //Inicializo el vector de salida a 0
        out[i]=0;

    for(i=0; i<N_dat; i++)
    {
        if(*(m) > dat[i]){*(m)=dat[i];}   //Encontramos el maximo y el minimo de la función
        if(*(M) < dat[i]){*(M)=dat[i];}
    }

    #ifdef DEBUG
    printf("maximo: %lf    minimo: %lf    \n", *M, *m);
    #endif

    *(delta)= (*(M)-*(m))/N_interval;   //Calculamos delta y nos aseguramos de que no es nulo para evitar dividir por 0
    if(*(delta)==0)
        {
        printf("Error: No se puede calcular el intervalo (Delta=0) .");

        exit(1);
        }

    for(i=0; i<N_dat; i++)
        {
        indicehist=(int)(dat[i]-(*m))/(*delta);
        out[indicehist]++; //Calculo del histograma (cada elemento de out, será una colunmna)
        }

    norm=N_dat* (*delta);

    #ifdef DEBUG
    printf("tamaño del intervalo: %lf   factor normalizacion: %lf", *delta, norm);
    #endif
    for(i=0; i<N_interval; i++)
        out[i]*= (1.0/norm);
}


//*******HISTOGRAMA_F*****************************************************************


void histograma_F(double *dat, double *out, double *m, double *M, double *delta, int N_dat, int N_interval, char *nombre_fichero)
{
    int i, indicehist;
    double norm;
    FILE *f;

    *(m)=10000000;  //Inicializo el máximo y el mínimo a la primera componente del array que contiene los datos.
    *(M)=-10000000;

    for(i=0; i<N_interval; i++) //Inicializo el vector de salida a 0
        out[i]=0;

    for(i=0; i<N_dat; i++)
    {
        if(*(m) > dat[i]){*(m)=dat[i];}   //Encontramos el maximo y el minimo de la función
        if(*(M) < dat[i]){*(M)=dat[i];}
    }

    #ifdef DEBUG
    printf("maximo: %lf    minimo: %lf    \n", *M, *m);
    #endif

    *(delta)= (*(M)-*(m))/(double)N_interval;   //Calculamos delta y nos aseguramos de que no es nulo para evitar dividir por 0

    printf("DELTA HIST %lf \n", *delta);

    if(*(delta)==0)
        {
        printf("Error: No se puede calcular el intervalo (Delta=0) .");

        exit(1);
        }

    for(i=0; i<N_dat; i++)
        {
        indicehist=(dat[i]-(*m))/(*delta);
        out[indicehist]++; //Calculo del histograma (cada elemento de out, será la altura de una colunmna)
        }

    norm=(double)N_dat * (*delta);

    #ifdef DEBUG
    printf("tamaño del intervalo: %lf   factor normalizacion: %lf", *delta, norm);
    #endif
    for(i=0; i<N_interval; i++)
        out[i]*= (1.0/norm);

    f=fopen(nombre_fichero,"w");

    for(i=0; i<N_interval; i++)
        fprintf(f,"%d %lf %lf \n", i, i*(*delta)+ (*m), out[i]);

    fclose(f);
}



//********GEN_RAND***********************************************************************


double gen_rand( double a, double b)
    {

        return (a+(b-a)*(((double)rand())/((double)RAND_MAX+1)));
    }


//***************MEDIA Y VARIANZA********************************************************

void stat_med_var(double *serie, int N, double *media, double *varianza)
{
    int i;
    (*media)=(*varianza)=0;

    for(i=0; i<N; i++)
        {
            (*media)+=serie[i];

        }
    (*media)*=1.0/((double)N);

    for(i=0; i<N; i++)
        {
            (*varianza)+= (serie[i]-(*media))*(serie[i]-(*media));
        }
    (*varianza)*=(1.0/(N-1));

}

//************GENERADOR DE PARISI-RAPUANO**********************************

int PR_rand(void)
    {
        static unsigned int  v[256];
        static unsigned char wa,wb,wc, j=0;  //Utiliza unsigned char porque sólo puede adquirir valores en [0,255]
        int r;

        for (i=0; i<256; i++)
                *(v+i)=rand()+(rand()<<16);

        wa=30; wb=70; wc=97;

        }

      // printf("%d   %d  %d \n", wa,wb, wc);
        v[j]=v[wa]+v[wb];

        r=v[j]^v[wc];
        j++;
        wa++;
        wb++;
        wc++;


        return r;
    }
