// Implémenté par NGUYEN Van-Khoa(van-khoa.nguyen@imt-atlantique.net)

void intercovariance(double input_1[],double input_2[], int nb_samples, int step, int time_delay, double *R_av){
// well-done
    int I = 0;
    int J = 0;

    J = I + time_delay;

    // Matrice intercovariance
    double R[time_delay][time_delay];
    double X[time_delay][time_delay];


    matrix_initialize(&R[0][0],time_delay);

    // Number of loops.
    double n = 0;

    while (J <= nb_samples){

          // Etract the windows of the input_1, the input_2.
          double window1[time_delay], window2[time_delay];
          int i;

          for (i = I; i < J; i++){

              window1[i - I] = input_1[i];
              window2[i-I] = input_2[i];
              }

          // Calculate the vector multiplication.
          vector_dot(window1,window2,&X[0][0],time_delay);


          // Addition operator between two matrix
          matrix_add(&R[0][0],&X[0][0],time_delay);

          n = n+1;

          I = I + step;

          J = I + time_delay;

          }
    // Division a matrix for a scalar
    matrix_scalar_division(&R[0][0],&R_av,n,time_delay);

    }


void off(double* R, int time_delay, double *F){
// well-done
     double sum = 0;

     int i,j;
     for(i = 0; i < time_delay; i++){
           for(j = 0; j< time_delay; j++){
                 if(i != j){
                      sum = sum + *(R  + i*time_delay +j);
                      }
                 }
           }

     *F  = sum/(time_delay*(time_delay - 1));
     }


void tr(double *R, int time_delay, double *F){
// well done

    double sum = 0;

    int i,j;
    for(i=0; i < time_delay; i++){

        for(j=0; j< time_delay; j++){

            if(i==j){

                sum = sum + *(R + i*time_delay + j);
            }
        }
    }

    *F = sum/time_delay;
}


void AFTabc(double x1[], double x2[],int nb_samples,int step,int time_delay,double gamma, double *A){
// well done

    double R11[time_delay][time_delay];
    double R22[time_delay][time_delay];
    double R12[time_delay][time_delay];

    double F1, F2, F12;
    double T1, T2, T12;

    double alpha, beta;
    double chi2;
    double d1, d2;
    double Ain[2][2];

    intercovariance(x1,x1,nb_samples,step,time_delay,&R11[0][0]);
    intercovariance(x1,x2,nb_samples,step,time_delay,&R12[0][0]);
    intercovariance(x2,x2,nb_samples,step,time_delay,&R22[0][0]);

    off(&R11[0][0], time_delay,&F1);
    off(&R12[0][0], time_delay,&F12);
    off(&R22[0][0], time_delay,&F2);

    tr(&R11[0][0], time_delay,&T1);
    tr(&R12[0][0], time_delay,&T12);
    tr(&R22[0][0], time_delay,&T2);

    alpha = 2*F12*T12 -(F1*(T2 - pow(gamma,2)) + F2*(T1 - pow(gamma,2)));
    beta = 2*(pow(T12,2) - (T1 - pow(gamma,2))*(T2 - pow(gamma,2)));

    chi2 = pow(F1*(T2 - pow(gamma,2)) - F2*(T1 - pow(gamma,2)),2) + 4*(F12*(T2 - pow(gamma,2)) - T12*F2)*(F12*(T1 -pow(gamma,2))-T12*F1);

    d1 = alpha - pow(chi2,0.5);
    d2 = alpha + pow(chi2,0.5);

    Ain[0][0] = beta*F1 - (T1 - pow(gamma,2))*d1;
    Ain[0][1] = beta*F12 - T12*d2;
    Ain[1][0] = beta*F12 -T12*d1;
    Ain[1][1] = beta*F2 - (T2 -pow(gamma,2))*d2;

    int i,j;
    for(i = 0; i< 2; i++){
        for(j=0; j<2; j++){
            *(A + 2*i + j) = Ain[i][j];
        }
    }
}


void separation_v0(double input_1[],double input_2[],double output_1[],double output_2[],int nb_samples,int time_delay){
    // well done
    double A[2][2];

    AFTabc(input_1,input_2,nb_samples,time_delay,time_delay,0,&A);
    normalized_matrix(&A);
    matrix_invers_2(&A);

    int i;
    for(i = 0; i< nb_samples; i++){
        output_1[i] = A[0][0]*input_1[i] + A[0][1]*input_2[i];
        output_2[i] = A[1][0]*input_1[i] + A[1][1]*input_2[i];
    }
    }
