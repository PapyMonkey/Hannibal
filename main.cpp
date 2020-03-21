#include "functions.h"

double pi = acos(-1.);
double G = 6.67408e-11;
int N_body = 9;


int main(){

    vector <object> body(N_body); //bodies
    vector <object> C(N_body); //copies (for the integrator)
    object COM(0.); //center of mass

    body[0] = object();//object(1.989e30); //sun
    body[1] = object(0.33011e24,0.38709893,0.20563069,7.00487,77.45645,48.33167,252.25084); //mercury
    body[2] = object(4.8675e24,0.72333199,0.00677323,3.39471,131.53298,76.68069,181.97973); //venus
    body[3] = object(5.972e24,1.00000011,0.01671022,0.00005,102.94719,-11.26064,100.46435); //earth
    body[4] = object(0.64171e24,1.52366231,0.09341233,1.85061,336.04084,49.57854,355.45332); //mars
    body[5] = object(1.89813e27,5.20336301,0.04839266,1.30530,14.75385,100.55615,34.40438); //jupiter
    body[6] = object(568.34e24,9.53707032,0.05415060,2.48446,92.43194,113.71504,49.94432); //saturn
    body[7] = object(86.813e24,19.19126393,0.04716771,0.76986,170.96424,74.22988,313.23218); //uranus
    body[8] = object(102.413e24,30.06896348,0.00858587,1.76917,44.97135,131.72169,304.88003); //neptune

    double GM = compute_GM(body);
    double T = 1.*orbital_period(GM,body[8].a); //simulation time
    int it = 0;
    double dt = orbital_period(GM,body[1].a)/1000.; //time step
    int snapshot = 100; //snapshot in time iterations
    double t = 0.; //secs
    double t0 = 0.; //J2000 epoch

    ofstream fileH("../CPP_OUTPUTS/H.txt");

    for(int i=0;i<N_body;i++){
        check_body_mass(body[i]);
        get_cartesian(GM,t,t0,body[i]);
        ofstream file;
        file.open("../CPP_OUTPUTS/Trajectories/m"+to_string(i)+".txt");
    }

    while(t < T){

        if((it%snapshot) == 0){
            cout << "Simulation Progress: " << floor(t*100./T)+1 << "%\r";
            update_COM(COM,body);
            for(int i=0;i<N_body;i++){
                ofstream file;
                file.open("../CPP_OUTPUTS/Trajectories/m"+to_string(i)+".txt",fstream::app);
                write_in_OBS(body[i],COM,file);
            }
            write_hamiltonian(body,t,fileH);
        }
        RK4(body,C,dt);
        t+=dt;
        it+=1;
    }

    return 0;
}

