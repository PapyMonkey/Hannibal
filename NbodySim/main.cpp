#include "functions.h"

double pi = acos(-1.);
double G = 6.67408e-11;
int N_body = 6;


int main(){

    vector <object> body(N_body); //bodies
    vector <object> C(N_body); //copies (for the integrator)
    object COM(0.); //center of mass

    body[0] = object(1.89813e27); //jupiter
    body[1] = object(893.2e20,421.8e3,0.004,0.04); //io
    body[2] = object(480.0e20,671.1e3,0.009,0.47); //europa
    body[3] = object(1481.9e20,1070.4e3,0.001,0.18); //ganymede
    body[4] = object(1075.9e20,1882.7e3,0.007,0.19); //callisto
    body[5] = object(104.,500.,.7,20.); //juice

    double GM = compute_GM(body);
    double T = 1.*orbital_period(GM,body[3].a); //simulation time
    int it = 0;
    double dt = orbital_period(G*body[3].m,body[5].a)/1000.; //time step
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

    get_relative_cartesian(G*body[3].m,body[3],t,t0,body[5]);

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

