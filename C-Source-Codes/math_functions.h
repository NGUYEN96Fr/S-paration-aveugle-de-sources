
void vector_dot(double v1[], double v2[], double *Xin, int time_delay){
// well-done
    int i,j;

    for(i = 0; i < time_delay; i++){

        for(j =0; j < time_delay; j++){
            *(Xin + time_delay*i + j) = v1[i]*v2[j];
        }
    }

}

void matrix_add(double *R,double *X, int time_delay){
// well-done
    int i,j;

    for (i = 0; i < time_delay; i++){

        for(j = 0; j < time_delay; j++){
              *(R + i*time_delay + j) =  *(R + i*time_delay + j) + *(X + i*time_delay + j);
              }
        }
}

void matrix_scalar_division(double *R,double **R_av, double n, int time_delay){
// well-done
     int i,j;

     for (i = 0; i < time_delay; i++){

         for(j=0; j< time_delay; j++){

                  *(*R_av + time_delay*i + j) = *(R + time_delay*i + j)/n;
                  }
         }
     }

void matrix_initialize(double *Ri,int n){
// well-done
    int i,j;

    for(i = 0; i < n; i++){
        for(j= 0; j<n; j++){
            *(Ri + n*i + j) = 0;
        }
    }
}


double dot_product(double v1[], double v2[], int nb_samples){
// well done
    int i;
    double v = 0;
    for(i = 0; i < nb_samples; i++){
           v = v + v1[i]*v2[i];
    }
    return v;
}


double norm(double v1[], int nb_samples){
// well done
    int i;
    double v = 0;
    for(i = 0; i < nb_samples; i++){
           v = v + v1[i]*v1[i];
    }
    v =  sqrt(v);

    return v;
}

void normalized_matrix(double* A){
// well done
    int i,j;
    double v = 0;

    // norm 1
    for(i = 0; i < 2; i++){
        for(j = 0; j< 2; j++)
           v = v + fabs(*(A + 2*i + j));
    }
    // normalize
    for(i = 0; i < 2; i++){
        for(j = 0; j< 2; j++){
             *(A + 2*i + j) = *(A + 2*i + j)/v;
        }
    }
}

void matrix_invers_2(double *A){
// well done
    double det;
    double a,b,c,d;
    a = *A;
    b = *(A + 1);
    c = *(A + 2);
    d = *(A + 3);

    det = a*d - b*c;

    *A = d/det;
    *(A+1) = -b/det;
    *(A + 2) = -c/det;
    *(A + 3) = a/det;
}
