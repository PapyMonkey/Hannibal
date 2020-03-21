#include "functions.h"

using namespace std;
extern double G;
extern double pi;
extern int N_body;

double rij(object &mi, object &mj){
    double x2 = (mj.x-mi.x)*(mj.x-mi.x);
    double y2 = (mj.y-mi.y)*(mj.y-mi.y);
    double z2 = (mj.z-mi.z)*(mj.z-mi.z);
    return sqrt(x2+y2+z2);
}

double gravity_x(object &mi, object &mj,double &rij_3){
    //Returns gravitational force between mi and mj along x axis. rij_3 is rij^3
    return G*mi.m*mj.m*(mj.x-mi.x)/rij_3;
}

double gravity_y(object &mi, object &mj,double &rij_3){
    //Returns gravitational force between mi and mj along y axis. rij_3 is rij^3
    return G*mi.m*mj.m*(mj.y-mi.y)/rij_3;
}

double gravity_z(object &mi, object &mj,double rij_3){
    //Returns gravitational force between mi and mj along z axis. rij_3 is rij^3
    return G*mi.m*mj.m*(mj.z-mi.z)/rij_3;
}

void compute_forces(vector<vector<double>> &Fx,vector<vector<double>> &Fy,vector<vector<double>> &Fz,vector<object> &B){
    //Computes all forces along x,y,z for all N_bodies in matrix of forces Fx,y,z (whose dimensions are N_body*N_body).
    double r,rij_3;
    for(int i=0;i<N_body-1;i++){
        for(int j=i+1;j<N_body;j++){
            r = rij(B[j],B[i]);rij_3 = r*r*r;
            Fx[i][j] = gravity_x(B[i],B[j],rij_3);
            Fx[j][i] = -Fx[i][j];
            Fy[i][j] = gravity_y(B[i],B[j],rij_3);
            Fy[j][i] = -Fy[i][j];
            Fz[i][j] = gravity_z(B[i],B[j],rij_3);
            Fz[j][i] = -Fz[i][j];
        }
    }
}

void RK4(vector<object> &bodies,vector<object> &C,double dt){
    //Order 4 Runge-Kutta N body integrator with timestep dt. Uses copies (C) of bodies for my convenience.

    double k1x[N_body],k2x[N_body],k3x[N_body],k4x[N_body];
    double k1y[N_body],k2y[N_body],k3y[N_body],k4y[N_body];
    double k1z[N_body],k2z[N_body],k3z[N_body],k4z[N_body];
    double l1x[N_body],l2x[N_body],l3x[N_body],l4x[N_body];
    double l1y[N_body],l2y[N_body],l3y[N_body],l4y[N_body];
    double l1z[N_body],l2z[N_body],l3z[N_body],l4z[N_body];

    vector<vector<double>> Fx;
    Fx.resize(N_body, vector<double>(N_body, .0)); //create N_body*N_body matrix with 0 as values
    vector<vector<double>> Fy;
    Fy.resize(N_body, vector<double>(N_body, .0));
    vector<vector<double>> Fz;
    Fz.resize(N_body, vector<double>(N_body, .0));

    for(int i=0;i<N_body;i++){
        C[i] = bodies[i];
    }

    compute_forces(Fx,Fy,Fz,C);

    for(int i=0;i<N_body;i++){
        k1x[i]=0.;
        k1y[i]=0.;
        k1z[i]=0.;
        //computing RK coefficients:
        for(int j=0;j<N_body;j++){
            k1x[i] += Fx[i][j]/C[i].m;
            k1y[i] += Fy[i][j]/C[i].m;
            k1z[i] += Fz[i][j]/C[i].m;
        }
        k1x[i] *= dt;
        k1y[i] *= dt;
        k1z[i] *= dt;
        l1x[i] = dt*C[i].vx;
        l1y[i] = dt*C[i].vy;
        l1z[i] = dt*C[i].vz;
        //Updating copies for next coefficients:
        C[i].x = bodies[i].x+.5*l1x[i];
        C[i].y = bodies[i].y+.5*l1y[i];
        C[i].z = bodies[i].z+.5*l1z[i];
        C[i].vx = bodies[i].vx+.5*k1x[i];
        C[i].vy = bodies[i].vy+.5*k1y[i];
        C[i].vz = bodies[i].vz+.5*k1z[i];
    }

    compute_forces(Fx,Fy,Fz,C);

    for(int i=0;i<N_body;i++){
        k2x[i]=0.;
        k2y[i]=0.;
        k2z[i]=0.;
        //computing RK coefficients:
        for(int j=0;j<N_body;j++){
            k2x[i] += Fx[i][j]/C[i].m;
            k2y[i] += Fy[i][j]/C[i].m;
            k2z[i] += Fz[i][j]/C[i].m;
        }
        k2x[i] *= dt;
        k2y[i] *= dt;
        k2z[i] *= dt;
        l2x[i] = dt*C[i].vx;
        l2y[i] = dt*C[i].vy;
        l2z[i] = dt*C[i].vz;
        //Updating copies for next coefficients:
        C[i].x = bodies[i].x+.5*l2x[i];
        C[i].y = bodies[i].y+.5*l2y[i];
        C[i].z = bodies[i].z+.5*l2z[i];
        C[i].vx = bodies[i].vx+.5*k2x[i];
        C[i].vy = bodies[i].vy+.5*k2y[i];
        C[i].vz = bodies[i].vz+.5*k2z[i];
    }

    compute_forces(Fx,Fy,Fz,C);

    for(int i=0;i<N_body;i++){
        k3x[i]=0.;
        k3y[i]=0.;
        k3z[i]=0.;
        //computing RK coefficients:
        for(int j=0;j<N_body;j++){
            k3x[i] += Fx[i][j]/C[i].m;
            k3y[i] += Fy[i][j]/C[i].m;
            k3z[i] += Fz[i][j]/C[i].m;
        }
        k3x[i] *= dt;
        k3y[i] *= dt;
        k3z[i] *= dt;
        l3x[i] = dt*C[i].vx;
        l3y[i] = dt*C[i].vy;
        l3z[i] = dt*C[i].vz;
        //Updating copies for next coefficients:
        C[i].x = bodies[i].x+l3x[i];
        C[i].y = bodies[i].y+l3y[i];
        C[i].z = bodies[i].z+l3z[i];
        C[i].vx = bodies[i].vx+k3x[i];
        C[i].vy = bodies[i].vy+k3y[i];
        C[i].vz = bodies[i].vz+k3z[i];
    }

    compute_forces(Fx,Fy,Fz,C);

    for(int i=0;i<N_body;i++){
        k4x[i]=0.;
        k4y[i]=0.;
        k4z[i]=0.;
        //computing RK coefficients:
        for(int j=0;j<N_body;j++){
            k4x[i] += Fx[i][j]/C[i].m;
            k4y[i] += Fy[i][j]/C[i].m;
            k4z[i] += Fz[i][j]/C[i].m;
        }
        k4x[i] *= dt;
        k4y[i] *= dt;
        k4z[i] *= dt;
        l4x[i] = dt*C[i].vx;
        l4y[i] = dt*C[i].vy;
        l4z[i] = dt*C[i].vz;
        //UPDATING POSITIONS AND VELOCITIES:
        bodies[i].vx += 1./6. * (k1x[i]+2.*k2x[i]+2.*k3x[i]+k4x[i]);
        bodies[i].vy += 1./6. * (k1y[i]+2.*k2y[i]+2.*k3y[i]+k4y[i]);
        bodies[i].vz += 1./6. * (k1z[i]+2.*k2z[i]+2.*k3z[i]+k4z[i]);
        bodies[i].x += 1./6. * (l1x[i]+2.*l2x[i]+2.*l3x[i]+l4x[i]);
        bodies[i].y += 1./6. * (l1y[i]+2.*l2y[i]+2.*l3y[i]+l4y[i]);
        bodies[i].z += 1./6. * (l1z[i]+2.*l2z[i]+2.*l3z[i]+l4z[i]);
    }
}

void get_cartesian(double GM,double t,double t0,object &m){
//get cartesian elements from keplerian elements
//MAKE SURE THAT t AND t0 ARE IN SECONDS
    double a = m.a*1.496e11;
    double e = m.e;
    //1: setting M(t)
    double M;
    if(t==t0){
        M = m.M*pi/180.;
    }
    else {
        double dt = t-t0;
        if (a == 0.){
            M = 0.;
        } else {
            M = m.M*pi/180.+dt*sqrt(GM/(a*a*a));
        }
    }
    //2: Computing eccentric anomaly with Newton-Raphson
    double E = M; double E0 = E;
    for(int i=0;i<500;i++){
        E0 = E;
        E = E - (E-e*sin(E)-M)/(1-e*cos(E));
        if(fabs(E-E0) < 1e-8){
            break;
        }
    }
    //3: Obtaining true anomaly
    double mu = 2.*atan2(sqrt(1.+e)*sin(0.5*E),sqrt(1.-e)*cos(0.5*E));
    //4: Obtaining distance to central body
    double rc = a*(1.-e*cos(E));
    //5: Obtaining x,y,z,vx,vy,vz in the orbital frame
    m.x = rc*(cos(mu)); m.vx = sqrt(GM*a)/rc *-sin(E);
    m.y = rc*(sin(mu)); m.vy = sqrt(GM*a)/rc *sqrt(1-e*e)*cos(E);
    m.z = 0.; m.vz = 0.;
    if(rc == 0.){
        m.vx = 0.;
        m.vy = 0.;
        m.vz = 0.;
    }
    //6: Obtaining x,y,z,vx,vy,vz in the inertial frame
    double w = m.w*pi/180.;
    double Om = m.Om*pi/180.;
    double i = m.i*pi/180.;
    double x = m.x; double y = m.y;
    double vx = m.vx; double vy = m.vy;

    m.x = x*(cos(w)*cos(Om)-sin(w)*cos(i)*sin(Om)) - y*(sin(w)*cos(Om)+cos(w)*cos(i)*sin(Om));
    m.y = x*(cos(w)*sin(Om)+sin(w)*cos(i)*cos(Om)) + y*(cos(w)*cos(i)*cos(Om)-sin(w)*sin(Om));
    m.z = x*(sin(w)*sin(i)) + y*(cos(w)*sin(i));
    m.vx = vx*(cos(w)*cos(Om)-sin(w)*cos(i)*sin(Om)) - vy*(sin(w)*cos(Om)+cos(w)*cos(i)*sin(Om));
    m.vy = vx*(cos(w)*sin(Om)+sin(w)*cos(i)*cos(Om)) + vy*(cos(w)*cos(i)*cos(Om)-sin(w)*sin(Om));
    m.vz = vx*(sin(w)*sin(i)) + vy*(cos(w)*sin(i));

}

void get_keplerian(double GM,object &m,object COM){
    //get keplerian elements from cartesian elements
    //1: Preparations
    double x = m.x-COM.x;
    double y = m.y-COM.y;
    double z = m.z-COM.z;
    double vx = m.vx-COM.vx;
    double vy = m.vy-COM.vy;
    double vz = m.vz-COM.vz;
    double r = sqrt(x*x+y*y+z*z);
    double v = sqrt(vx*vx+vy*vy+vz*vz);
    //angular momentum:
    double hx = y*vz-z*vy;
	double hy = z*vx-x*vz;
	double hz = x*vy-y*vx;
	double h = sqrt(hx*hx+hy*hy+hz*hz);
	//obtaining eccentricity:
	double ex = (vy*hz-vz*hy)/GM-x/r;
	double ey = (vz*hx-vx*hz)/GM-y/r;
	double ez = (vx*hy-vy*hx)/GM-z/r;
	m.e = sqrt(ex*ex+ey*ey+ez*ez);
	//obtaining node vector:
	double nx = hy;
	double ny = -hx;
	double n=sqrt(nx*nx+ny*ny);
	//obtaining semimajor axis:
	m.a = -GM/2./(v*v/2. - GM/r);
	//2: Obtaining inclination
	m.i = acos(hz/h);
	//3: Obtaining longitude of the ascending node and argument of periapsis
	if(fabs(sin(m.i)) < 1e-8){
        m.Om = 0.;
	} else {
        m.Om = acos(-nx/n);
        if(m.Om && ny>0.){
            m.Om = 2.*pi-m.Om;
        }
	}
    if(m.e == 0.){
        m.w = 0.;
    }
	else {
		m.w = acos((cos(m.Om)*ex+sin(m.Om)*ey)/m.e);
		if(m.w && ez<0.){
            m.w = 2.*pi-m.w;
		}
	}
	//4: Obtaining true anomaly and eccentric anomaly
	double cosi = cos(hz/h);
	double comega = cos(m.Om);
	double somega = sin(m.Om);
	double cosR = x*comega+y*somega;
	double sinR = (y*comega-x*somega)*cosi+z*sin(m.i);
	double norm = sqrt(cosR*cosR+sinR*sinR);
    cosR /= norm;
    sinR /= norm;
	double T = acos(cosR);
	if(sinR < 0){
        T = 2*pi-T;
	}
	T -= m.w;
	double E = acos((cos(T)+m.e)/(1.+m.e*cos(T)));
	if(sin(T) < 0.){
        E = 2.*pi-E;
	}
	//5: Obtaining Mean anomaly
	m.M = E-m.e*sin(E);

	//6: Converting
	m.a /= 1.496e11;
	m.i *= 180./pi;
	m.w *= 180./pi;
	m.Om *= 180./pi;
	m.M *= 180./pi;
}

void write_in_OBS(object m,object OBS,ofstream& file){
    //writes coordinates of m in the frame of the observer
    m.x -= OBS.x;
    m.y -= OBS.y;
    m.z -= OBS.z;
    m.vx -= OBS.vx;
    m.vy -= OBS.vy;
    m.vz -= OBS.vz;
    m.write_out(file);
}

double orbital_period(double GM, double a){
    //computes orbital period using Kepler law
    a *= 1.496e11;
    double T = sqrt(4*pi*pi*a*a*a/GM);
    return T;
}

void write_hamiltonian(vector<object> &bodies,double t,ofstream& file){
    double T = 0.; //kinetic energy
    double U = 0.; //potential energy
    double vi_2 = 0.; //speed squared of body i

    for(int i=0;i<N_body;i++){
        vi_2 = bodies[i].vx*bodies[i].vx+bodies[i].vy*bodies[i].vy+bodies[i].vz*bodies[i].vz;
        T += .5*bodies[i].m*vi_2;
        for(int j=i+1;j<N_body;j++){
            U -= G*bodies[i].m*bodies[j].m/rij(bodies[i],bodies[j]);
        }
    }
    double H = T+U;
    file << t << "\t" << H << endl;
}

void update_COM(object &COM,vector<object> &bodies){
    double mi_xi = 0.; //products of mi and xi
    double mi_yi = 0.; //products of mi and yi
    double mi_zi = 0.; //products of mi and zi
    double mi_vxi = 0.; //products of mi and vxi
    double mi_vyi = 0.; //products of mi and vyi
    double mi_vzi = 0.; //products of mi and vzi
    double mi = 0.; //sum of all masses
    for(int i=0;i<N_body;i++){
        mi_xi += bodies[i].m*bodies[i].x;
        mi_yi += bodies[i].m*bodies[i].y;
        mi_zi += bodies[i].m*bodies[i].z;
        mi_vxi += bodies[i].m*bodies[i].vx;
        mi_vyi += bodies[i].m*bodies[i].vy;
        mi_vzi += bodies[i].m*bodies[i].vz;
        mi += bodies[i].m;
    }
    COM.x = mi_xi/mi;
    COM.y = mi_yi/mi;
    COM.z = mi_zi/mi;
    COM.vx = mi_vxi/mi;
    COM.vy = mi_vyi/mi;
    COM.vz = mi_vzi/mi;
}

double fRand(double fMin, double fMax){
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

double compute_GM(vector<object> &bodies){
    //returns G*(sum of masses)
    double M = 0.; //total mass of the N bodies
    for(int i=0;i<N_body;i++){
        M += bodies[i].m;
    }
    return G*M;
}

void check_body_mass(object &b){
//Checks if current body has null mass (i.e m=0) raises error if not.
    try{
        if(b.m == 0.){
            throw "ERROR: NULL MASS BODY ENCOUNTERED";
        }
    } catch (const char* err) {
     cerr << err << endl;
     exit(1);
   }
}
