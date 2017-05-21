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
static struct hop {int x,y,z,o1,o2;double re,im;long double ekin;} H[3000000];
static int nHop,nWan;
void load(char* hamFile) {
    char str[100];
    FILE *in=fopen(hamFile,"r");
    if (in==NULL) { printf("file not found: %s\n",str); exit(0); }
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
    double eT=0.01;
    const int nMu=100;
    double mu[nMu];
    for (i=0; i<nMu; i++) mu[i]=7.25+(7.85-7.25)*i/nMu;
//    printf("argc=%d\n",argc);
    char str[100];
//    if (argc>2) sscanf(argv[2],"%lf",&mu);
//    printf("E_Fermi= %lf eV\n",mu);
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
    double occ[nMu];
    for (kx=M_PI/nK[0];kx<2.*M_PI;kx+=2.*M_PI/nK[0])
        for (ky=M_PI/nK[1];ky<2.*M_PI;ky+=2.*M_PI/nK[1])
            for (kz=M_PI/nK[2];kz<2.*M_PI;kz+=2.*M_PI/nK[2]) {
                Hk(hk,kx,ky,kz);
                info = LAPACKE_zheev( LAPACK_COL_MAJOR, 'V', 'U', nWan, (MKL_Complex16 *)hk, nWan, w );
                //eigenvectors are contigious in a
                for (i=0;i<nWan;i++) {
                    //printf("eig %d=%lf\t",i,w[i]);
                    int j;
#pragma omp parallel for
                    for (j=0;j<nMu;j++)
                        occ[j]+=1./(exp((w[i]-mu[j])/eT)+1.);
                }
            }
    runtime+=clock();
    for (j=0;j<nMu-1;j++)
        if (    (2.*occ[j]/nKs<155.)&&  (2.*occ[j+1]/nKs>=155.)&&  
                    (fabs(2.*occ[j+1]/nKs-2.*occ[j]/nKs)<1e-1)) {
            printf("j=%d Fermi level= %lf, total occupation= %g\n",j,mu[j],2.*occ[j]/nKs); 
            break; }
    free(w); free(hk);
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
