
#include <string.h>

void SHMDistance(int *Nstrings, int *Nnucs, char **S,  double *MATRIX)
{
     double SymmetricDistanceArray[36] = { 0.000000, 2.864539, 1.000000, 2.138698, 0.000000, 0.000000, 2.864539, 0.000000, 2.138698, 1.000000, 0.000000, 0.000000, 1.000000, 2.138698, 0.000000, 2.864539, 0.000000, 0.000000, 2.138698, 1.000000, 2.864539, 0.000000, 0.000000, 0.000000 ,0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000 ,0.000000, 0.000000, 0.000000, 0.000000 };
    int n = Nstrings[0];
    int N[Nnucs[0]][Nstrings[0]];
    char S1[Nnucs[0]],S2[Nnucs[0]];
    char C[1];
    int i,j,k;
    double d12;    
     for(j=0;j<Nstrings[0];j++){
        strcpy(S1,S[j]);
        for(i=0;i<Nnucs[0];i++){
            memcpy(C,S1+i,1);
            N[i][j]=atoi(C);
        }          
     }

     for(j=0;j<Nstrings[0];j++){
        for(i=j+1;i<Nstrings[0];i++){
            d12=0.0;
            for(k=0;k<Nnucs[0];k++){
            d12+=SymmetricDistanceArray[(N[k][i]-1)*6+N[k][j]-1];
            }
         MATRIX[(j)*Nstrings[0]+i]=MATRIX[(i)*Nstrings[0]+j]=d12;   
        }
     } 
}




