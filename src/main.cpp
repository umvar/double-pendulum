#include <iostream>
#include <stdio.h>
#include <cmath>
#include "cpgplot.h"

// Gravitational acceleration
float g = 9.81;

// Pendulum 1 mass, length
float m1 = 1.0;
float l1 = 1.0;

// Pendulum 2 mass, length
float m2 = 1.0;
float l2 = 1.0;

// Derivative of double pendulum equations of motion
void derive(float t, float state_in[], float dydx[]) {
    float theta1 = state_in[0];
    float theta2 = state_in[1];
    float omega1 = state_in[2];
    float omega2 = state_in[3];

    float theta1_ = omega1;
    float theta2_ = omega2;

    float top;
    float bottom;
    float delta = theta2 - theta1;

    top = m2*l1*omega1*omega1*sin(delta)*cos(delta) + m2*(g*sin(theta2)*cos(delta) + l2*omega2*omega2*sin(delta)) - (m1 + m2)*g*sin(theta1);
    bottom = l1*(m1 + m2*(1 - cos(delta)*cos(delta)));

    float omega1_ = top/bottom;

    top = -m2*l2*omega2*omega2*sin(delta)*cos(delta) + (m1 + m2)*(g*sin(theta1)*cos(delta) - l1*omega1*omega1*sin(delta)) - (m1 + m2)*g*sin(theta2);
    bottom = l2*(m1 + m2*(1 - cos(delta)*cos(delta)));

    float omega2_ = top/bottom;

    dydx[0] = theta1_;
    dydx[1] = theta2_;
    dydx[2] = omega1_;
    dydx[3] = omega2_;
}

// Euler integration
void euler(float t, float state_in[], float state_out[], float dt) {
  float dydx[4];
  derive(t, state_in, dydx);
  for (int i = 0; i < 4; i++) {
    state_out[i] = state_in[i] + dydx[i]*dt;
  }
}



// Runge-Kutta 4 integration
void rk4(float t, float state_in[], float state_out[], float dt) {
    float k1[4], k2[4], k3[4], k4[4], dydx[4], dydx_temp[4], y_temp[4];

    derive(t, state_in, dydx);
    for (int i = 0; i < 4; i++) {
        k1[i] = dt*dydx[i];
        y_temp[i] = state_in[i] + 0.5*k1[i];
    }

    derive(t + 0.5*dt, y_temp, dydx_temp);
    for (int i = 0; i < 4; i++) {
        k2[i] = dt*dydx_temp[i];
        y_temp[i] = state_in[i] + 0.5*k2[i];
    }

    derive(t + 0.5*dt, y_temp, dydx_temp);
    for (int i = 0; i < 4; i++) {
        k3[i] = dt*dydx_temp[i];
        y_temp[i] = state_in[i] + k3[i];
    }

    derive(t + dt, y_temp, dydx_temp);
    for (int i = 0; i < 4; i++) {
        k4[i] = dt*dydx_temp[i];
        state_out[i] = state_in[i] + (k1[i] + 2.0*(k2[i] + k3[i]) + k4[i])/6.0;
    }
}

int main() {
    // Initial angles
    float theta1 = 77*3.14/180;
    float theta2 = 135*3.14/180;

    // Initial angular velocities
    float omega1 = 0.0;
    float omega2 = 0.0;

    // Time range
    float t_end = 25.0;
    float t_start = 0.0;

    // Number of steps
    int steps = 2500;

    // Time step
    float dt = (t_end - t_start)/steps;

    // Initialize arrays for storing time, theta1, theta2, omega1, omega2
    float t[steps], th1[steps], th2[steps], om1[steps], om2[steps];
    th1[0] = theta1;
    th2[0] = theta2;
    om1[0] = omega1;
    om2[0] = omega2;

    // Initialize arrays for storing x and y coordinates of pendulum 1 and 2
    float x1[steps], y1[steps], x2[steps], y2[steps];

    // Set time array given time step and number of steps
    for (int i = 0; i < steps; i++) {
        t[i] = i * dt;
    }

    // Initialize array vectors for storing state variables
    float state_in[4]; // state_in = [theta1, theta2, omega1, omega2]
    float state_out[4]; // state_out = [theta1, theta2, omega1, omega2] after integrating


    // Integrate and find theta1, theta2, omega1, omega2
    for (int i = 0; i < steps - 1; i++) {
        // set state_in
        state_in[0] = th1[i];
        state_in[1] = th2[i];
        state_in[2] = om1[i];
        state_in[3] = om2[i];

        // integrate
        //euler(t[i], state_in, state_out, dt)
        rk4(t[i], state_in, state_out, dt); 

        // set state_out after integrating
        th1[i+1] = state_out[0];
        th2[i+1] = state_out[1];
        om1[i+1] = state_out[2];
        om2[i+1] = state_out[3]; 
    }

    // Find trajectories of pendulums and store in x1, y1, x2, y2
    for (int i = 0; i < steps; i++) {
        x1[i] = l1*sin(th1[i]);
        y1[i] = -l1*cos(th1[i]);
        x2[i] = x1[i] + l2*sin(th2[i]);
        y2[i] = y1[i] - l2*cos(th2[i]);
    }

    // PLOTTING

    // Found that when using /XSERVE we can have multiple windows open.
    if (!cpgopen("/XSERVE")) return 1;

    // Plot angles of pendulums vs time
    cpgopen("/XSERVE");
    cpgenv(0., 25., -10., 40., 0, 1);
    cpglab("t", "theta1, theta2 (radians)", "Angles of Pendulums vs Time");
    cpgsci(4);
    cpgline(steps, t, th1);
    cpgsci(5);
    cpgline(steps, t, th2);

    // Plot angular velocities of pendulums vs time
    cpgopen("/XSERVE");
    cpgenv(0., 25., -10., 40., 0, 1);
    cpglab("t", "omega1, omega2 (radians)", "Angular Velocities of Pendulums vs Time");
    cpgsci(4);
    cpgline(steps, t, om1);
    cpgsci(5);
    cpgline(steps, t, om2);

    // Plot trajectories of pendulums
    cpgopen("/XSERVE");
    cpgenv(-2., 2., -2., 2., 0, 0);
    cpglab("x", "y", "Trajectories of Pendulums");
    cpgsci(4);
    cpgline(steps, x1, y1);
    cpgsci(5);
    cpgline(steps, x2, y2);
}