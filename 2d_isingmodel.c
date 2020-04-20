#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<time.h>

//Global Variables:
double En = 0; //energy
double Mag = 0;//magnetization 

double spin_arr1[1000][1000] = {1.}; //array for 2d spin array
double sp_arr1[1000][1000] = {1.};

//take a random number in range of N
int rand_int(int N){
    int rn = rand() % N;
    return rn; 
}
double rand_number(){
    double rn = ((double)rand())/((double) RAND_MAX);//going from 0 to 1
    return rn; 
}

//update the 2d array
void update2d(int N, double kT){
   
    int i = rand_int(N);
    int j = rand_int(N);
    int flip = 0;
    
    //change in energy from nearest neighbors
    double dE = 2*sp_arr1[i][j]*( sp_arr1[i-1][j] + sp_arr1[(i+1)%N][j]
                      + sp_arr1[i][j-1] + sp_arr1[i][(j+1)%N]);
    
    if (dE < 0.0){//predicting the likelihood of a dipole flipping
        flip=1;
    }else{
        double p = exp(-1*dE/kT);
        if(rand_number() < p){
            flip = 1;
        }
    }
    if(flip == 1){
        En = En + dE;
        Mag = Mag - 2*sp_arr1[i][j];
        sp_arr1[i][j] = -1*sp_arr1[i][j];
    }
}

void initialize(int N){ //initialize the array
    double p = .45; //this value dictates how if there will be an up or down flip- 0.5 is neutral
    //double E = 0;
    //double M = 0;
    for(int i = 1;i<N;i++){
        for(int j = 1; j<N;j++){
            if (rand_number() < p){
                sp_arr1[i][j]=-1;
            }
            En = En - sp_arr1[i-1][j-1]*sp_arr1[i][j];
            Mag = Mag + sp_arr1[i][j];
        }
    }

    En = En - sp_arr1[N-1][N-1]*sp_arr1[0][0];
    Mag = Mag + sp_arr1[0][0];
}

//reset the array for multiple runs
void reset(int N){
    En = 0;
    Mag = 0;
    for(long int x = 0; x<N;x++){
        for(long int y = 0; y<N;y++){
            sp_arr1[x][y] = 1.;
        }
    }
}

int main(){
    FILE *fp; //opens file
    fp = fopen("/home/jpl/testdata_2d.dat","w" ); //writes files
    srand(time(NULL));
    //printf("%ld\n", sizeof(sp_arr1));
    int N = 1000;
    int passes = 10;
    long int Nmc = passes*pow(N,2);
    printf("#Temp  Energy  Mag\n");
    for(int i=1; i<=60;i++){ //run this code multiple times and take the average

        double kT = 0.1*i;
        reset(N);
        initialize(N);
        for(long int k = 0; k < Nmc; k++){
            update2d(N, kT);
        }
        
        double E1 = 0.;
        double M1 = 0.;
        for(long int k = 0; k < Nmc; k++){
            update2d(N, kT);
            E1 = E1 + En;
            M1 = M1 + Mag;
        }
        E1 = E1/Nmc;
        M1 = M1/Nmc;
        E1 = E1/pow(N,2);
        //M1 = fabs(M1/pow(N,2));
        fprintf(fp,"%2.2f  %2.5e %2.5e\n", kT, E1, M1);
        printf("%2.2f  %2.5e %2.5e\n", kT, E1, M1);
    }
    fclose(fp); 
    return 0;
}
