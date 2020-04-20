#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<time.h>

//Global Variables:
double En = 0; //energy
double Mag = 0;//magnetization
double spin_arr[100][100][100] = {1.};//array for 3d spin array

//take a random number in range of N
int rand_int(int N){
    int rn = rand() % N;
    return rn; 
}
double rand_number(){
    double rn = ((double)rand())/((double) RAND_MAX);//going from 0 to 1
    return rn; 
}

//update the 3d array
void update3d(int N, double kT){
    //the position of the dipole (index of the array)
    int i = rand_int(N);
    int j = rand_int(N);
    int k = rand_int(N);
    int flip = 0;
    
    //change in energy from nearest neighbors
    double dE = 2*spin_arr[i][j][k]*( spin_arr[i-1][j][k] + spin_arr[(i+1)%N][j][k]
                      + spin_arr[i][j-1][k] + spin_arr[i][(j+1)%N][k] 
		      +spin_arr[i][j][k-1] + spin_arr[i][j][(k+1)%N] );
    
    if (dE < 0.0){//predicting the likelihood of a dipole flipping
        flip=1;
    }else{
        double p = exp(-1*dE/kT);
        if(rand_number(N) < p){
            flip = 1;
        }
    }
    if(flip == 1){
        En = En + dE;
        Mag = Mag - 2*spin_arr[i][j][k];
        spin_arr[i][j][k] = -1*spin_arr[i][j][k];
    }
}

void initialize(int N){//initialize the array
    double p = .5;//this value dictates how if there will be an up or down flip- 0.5 is neutral
    //double E = 0;
    //double M = 0;
    for(int i = 1;i<N;i++){
        for(int j = 1; j<N;j++){
            for(int k = 1; k<N;k++){
                if (rand_number() < p){
                    spin_arr[i][j][k]=-1;
                }
                En = En - spin_arr[i-1][j-1][k-1]*spin_arr[i][j][k];
                Mag = Mag + spin_arr[i][j][k];
                
            }
        }
    }

    En = En - spin_arr[N-1][N-1][N-1]*spin_arr[0][0][0];
    Mag = Mag + spin_arr[0][0][0];
}

int main(){
    FILE *fp; //opens file
    fp = fopen("/home/jpl/testdata_3d.dat","w" );
    srand(time(NULL));
    int N = 100;
    int passes = 100;
    long int Nmc = passes*pow(N,3);
    printf("#Temp Energy       Mag\n");
    for(int i=1; i<=60;i++){//run the code multiple times
        En = 0;
        Mag = 0;
        for(long int x = 0; x<N;x++){
            for(long int y = 0; y<N;y++){
                for(long int z = 0; z<N;z++){
                    spin_arr[x][y][z] = 1;
                }
            }
        }
        double kT = 0.1*i;
        initialize(N); 
        for(long int k = 0; k < Nmc; k++){
            //printf("%ld\n",k);
            update3d(N, kT);
        }
        
        //reset the array
        double E1 = 0.;
        double M1 = 0.;
        for(long int k = 0; k < Nmc; k++){
            update3d(N, kT);
            E1 = E1 + En;
            M1 = M1 + Mag;
        }
        
        //getting final values
        E1 = E1/Nmc;
        M1 = M1/Nmc;
        E1 = E1/pow(N,3);
        M1 = M1/pow(N,3);
        fprintf(fp,"%2.2f  %2.5e %2.5e\n", kT, E1, M1);//printing values
        printf("%2.2f  %2.5e %2.5e\n", kT, E1, M1);
    }
    
    fclose(fp);
    return 0;
}
