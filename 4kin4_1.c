#include "WolframLibrary.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "mkl_lapacke.h"
#include "mkl_solvers_ee.h"
#include <time.h>
#include "omp.h"
void Err(const char* str) {
    FILE *log=fopen("/home/urko/CdSe/1.log","a");
    fprintf(log,"%s\n",str);//,out_dims[0],out_dims[1]);
//    fflush(log);
    fclose(log);
    return;
}
static struct hop {int x,y,z,o1,o2; double re,im; long double ekin;} H[3000000];
static int nHop,nWan;
void load(char* hamFile) {
    char str[100];
    FILE *in=fopen(hamFile,"r");
    if (in==NULL) {printf("file not found: %s\n",str); exit(0);}
    fgets(str,90,in);
    int skip,i,tmp;
    fscanf(in,"%d %d",&nWan,&skip);
    for (i=0;i<skip;i++) fscanf(in,"%d",&tmp);
    struct hop *p;
    for (p=H;!feof(in);p++) {
        fscanf(in,"%d %d %d %d %d %lf %lf",&(p->x),&(p->y),&(p->z),&(p->o1),&(p->o2),&(p->re),&(p->im));
        p->ekin=0.;
    }
    nHop=p-H-1;
//    sprintf(str,"%d %d %d %d %d %lf %lf %d\n",x,y,z,o1,o2,re,im,rows); Err(str);
    fclose(in);
    return;
}
complex double* Hk(complex double *hk, double kx, double ky, double kz) {
    int i,j;
    for(i=0;i<nWan*nWan;i++)
        hk[i]=0.;
    for(i=0;i<nHop;i++)
        hk[(H[i].o1-1)*nWan+(H[i].o2-1)]+=(H[i].re + H[i].im*I)*cexp(-(kx*H[i].x+ky*H[i].y+kz*H[i].z)*I);
    return hk;
}
// Auxiliary routine: printing a matrix
extern void print_matrix( char* desc, MKL_INT m, MKL_INT n, double* a, MKL_INT lda );
void print_matrix( char* desc, MKL_INT m, MKL_INT n, double* a, MKL_INT lda ) {
        MKL_INT i, j;
        printf( "\n %s\n", desc );
        for( i = 0; i < m; i++ ) {
                for( j = 0; j < n; j++ ) printf( " %6.2f", a[i*lda+j] );
                printf( "\n" );
        }
}
void writeRho(complex double *rho) {
    int i,j;
    FILE *out=fopen("/home/urko/CdSe/1.log","w");
    for (i=0;i<nWan;i++) {
        for (j=0;j<nWan;j++)
            fprintf(out,"%lf %+lfI ",creal(rho[i*nWan+j]),cimag(rho[i*nWan+j]));
//            fprintf(out,"% 4.2lf%+4.2lfI ",creal(rho[i*nWan+j]),cimag(rho[i*nWan+j]));
        fprintf(out,"\n");
    }
    fclose(out);
}
int main(int argc, char **argv) {
    int image;
    if (argc<2) {printf("error>specify the Hamiltonian file name as a first argument\n");exit(0);}
    printf("Image file=%s\n", argv[1]);
    load(argv[1]);
    printf("sizeof(struct hop)=%d, nHop=%d\n",sizeof(struct hop),nHop); int nh,i,j,k;
//    for (nh=nHop-1;nh>nHop-10;nh--)
//	printf("%d %d %d %d %d %lf %lf\n",H[nh].x,H[nh].y,H[nh].z,H[nh].o1,H[nh].o2,H[nh].re,H[nh].im);
    MKL_INT info;
    double *w=malloc(sizeof(double)*nWan);
    double mu=7.7482, eT=0.01;
//    printf("argc=%d\n",argc);
    char str[100];
    if (argc>2) sscanf(argv[2],"%lf",&mu);
    printf("E_Fermi= %lf eV\n",mu);
    complex double *hk = (complex double *)malloc(sizeof(struct hop)*nWan*nWan);
//    complex double *rho=malloc(sizeof(complex double)*nWan*nWan);
//    int nK[]={10,10,10};
    int nK[]={5,5,5};
    if (argc>3) {
        sscanf(argv[3],"%d",&nK[0]);
        sscanf(argv[4],"%d",&nK[1]);
        sscanf(argv[5],"%d",&nK[2]);
    }
    printf("k_grid= %d %d %d \n",nK[0],nK[1],nK[2]);
    int nKs=1;
    for (i=0;i<3;i++) nKs*=nK[i];
    double kx,ky,kz; time_t runtime;
    runtime=-clock();

/*    printf("k=(0,0,0.25pi):");
    Hk(hk,0.,0.,0.25*M_PI);
    info = LAPACKE_zheev( LAPACK_COL_MAJOR, 'V', 'U', nWan, (MKL_Complex16 *)hk, nWan, w );
    printf("\neigenenergies:\n");
    for (i=0;i<nWan;i++) {
        printf("% 8.4lf ",w[i]);
        if((i+1)%8==0) printf("\n");
    }
    printf("\noccupations:\n");
    for (i=0;i<nWan;i++) {
        printf("% 8.4lf ",1./(exp((w[i]-mu)/eT)+1.));
        if((i+1)%8==0) printf("\n");
    }*/
/*    
// band structure Z(0,0,1)-GM-M(1,0,0)
    FILE *out=fopen("4kin.bands","w");
    kx=0;ky=0;
    for (kz=M_PI;kz>0;kz-=M_PI/20) {
        fprintf(out,"k=% 8.4lf %8.4lf %8.4lf w[eV]= ",kx,ky,kz);
        Hk(hk,kx,ky,kz);
        LAPACKE_zheev( LAPACK_COL_MAJOR, 'V', 'U', nWan, (MKL_Complex16 *)hk, nWan, w );
        for (i=0;i<nWan;i++)
            fprintf(out,"% 8.4lf ",w[i]);
        fprintf(out,"\n");
    }
    kx=0;ky=0;kz=0;
    for (kx=0;kx<=M_PI;kx+=M_PI/20) {
        fprintf(out,"k=% 8.4lf %8.4lf %8.4lf w[eV]= ",kx,ky,kz);
        Hk(hk,kx,ky,kz);
        LAPACKE_zheev( LAPACK_COL_MAJOR, 'V', 'U', nWan, (MKL_Complex16 *)hk, nWan, w );
        for (i=0;i<nWan;i++)
            fprintf(out,"% 8.4lf ",w[i]);
        fprintf(out,"\n");
    }
    fclose(out);
*/
    double totocc=0,ekin1=0,occ;
    for (kx=M_PI/nK[0];kx<2.*M_PI;kx+=2.*M_PI/nK[0])
        for (ky=M_PI/nK[1];ky<2.*M_PI;ky+=2.*M_PI/nK[1])
            for (kz=M_PI/nK[2];kz<2.*M_PI;kz+=2.*M_PI/nK[2]) {
//kx=M_PI/2;ky=0.2345*M_PI;kz=0*M_PI;
                Hk(hk,kx,ky,kz);
//                printf("k= %lf %lf %lf, H[%d,%d]=%lf%+lfI\n",kx,ky,kz,nWan,nWan,creal(hk[nWan*nWan-1]),cimag(hk[nWan*nWan-1]));
//              vzexp(&nWan,hk,hko);
                info = LAPACKE_zheev( LAPACK_COL_MAJOR, 'V', 'U', nWan, (MKL_Complex16 *)hk, nWan, w );
                //eigenvectors are contigious in a
                if( info > 0 ) {
                    printf( "LAPACKE_zheev failed to compute eigenvalues.\n" );
                    exit( 1 );
                }
// #pragma omp parallel for
                for (i=0;i<nWan;i++) {
                    //printf("eig %d=%lf\t",i,w[i]);
                    occ=1./(exp((w[i]-mu)/eT)+1.);
                    ekin1+=w[i]*occ;
                    totocc+=occ;

                    int j;
#pragma omp parallel for
                    for(j=0;j<nHop;j++)
                        H[j].ekin+=2*creal(
                            conj(hk[i*nWan+H[j].o1-1])*  hk[i*nWan+H[j].o2-1]*  (H[j].re + H[j].im*I)*  cexp((kx*H[j].x+ky*H[j].y+kz*H[j].z)*I)
                            )*occ/nKs;
//printf("Re t=%lf, Im t=%lf, Re v[i]=%lf, Im v[i]=%lf occ=%lf\n",H[i].re,H[i].im,creal(hk[i*nWan+H[j].o1-1]),cimag(hk[i*nWan+H[j].o1-1]),occ);
                }
//                printf("k=%lf %lf %lf\n",kx,ky,kz);
//                for (i=0;i<nWan;i++) printf("%lf %lf   ",i,creal(hk[i]),cimag(hk[i]));
//                exit(0);
            }
    ekin1/=0.5*nKs;
    totocc/=0.5*nKs;
    runtime+=clock();
    long double ekin=0;
    for (i=0;i<nHop;i++) ekin+=H[i].ekin;
    printf("elapsed time %lf sec.\n Total band energy= %lf eV= %lf Ry,\n total occ= %lf,\n sum of Kin contributions= %lf, \nfirst hopping contributes %lf eV\n",
        (double)runtime/CLOCKS_PER_SEC, ekin1,       ekin1*0.07349864,          totocc,                 (double)ekin,                   (double)H[0].ekin);
        
    FILE* out;
    out=fopen("contrib.out","w");
    for(j=0;j<nHop;j++)
        fprintf(out,"%d %d %d %d %d %lf %lf %Lg\n",(H[j].x),(H[j].y),(H[j].z),(H[j].o1),(H[j].o2),(H[j].re),(H[j].im),H[j].ekin);
    fclose(out);
    
//    writeRho(rho);
//    complex double count;
//    for (i=0,count=0;i<nWan;i++) count+=rho[i*(nWan+1)];
//    printf("total occ=%lf%+lfI\n",creal(count),cimag(count));
    /* Print eigenvalues */
//    print_matrix( "Eigenvalues", 1, nWan, w, 1 );
    /* Print eigenvectors */
//    print_matrix( "Eigenvectors (stored columnwise)", nWan, nWan, rho, nWan );
    free(w); free(hk); //free(rho);
    exit( 0 );
}
void printEig(complex double *hk,double *w) {
    int i,j;
    for (i=0;i<nWan;i++)
            printf("eig %d=%lf+%lfI\n",i,creal(w[i]),cimag(w[i]));
    for (i=0;i<nWan;i++)
        for (j=0;j<nWan;j++) printf("v[%d]=%lf%+lfI\n",i*nWan+j,creal(hk[i*nWan+j]),cimag(hk[i*nWan+j]));

    complex double *count=(complex double *)malloc(sizeof(complex double)*nWan);
    for (j=0;j<nWan;j++) {
        count[j]=0;
        for (i=0;i<nWan;i++) count[j]+=conj(hk[i+nWan*j])*hk[i+nWan*j];
        printf("norm 1=%lf+%lfI\n",creal(count[j]),cimag(count[j]));
    }
    for (i=0;i<nWan;i++) {
        count[i]=0;
        for (j=0;j<nWan;j++) count[i]+=conj(hk[i+nWan*j])*hk[i+nWan*j];
        printf("norm column 1=%lf+%lfI\n",creal(count[i]),cimag(count[i]));
    }
}
