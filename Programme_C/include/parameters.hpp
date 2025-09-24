#pragma once

//Paramètres Physiques
constexpr double D_COEFF=0.01;
// 1D
constexpr double C_COEFF=0.03;
// 2D
constexpr double C_COEFF_X=0.03;
constexpr double C_COEFF_Y=0.02;

//Paramètres de la Simulation
constexpr double SIMULATION_TIME=2.0;

//Paramètres par défaut du Maillage 1D
constexpr double DEFAULT_LENGTH=1.0;
constexpr int DEFAULT_NX=100;

//Paramètres par défaut pour le Maillage 2D
constexpr double DEFAULT_LENGTH_X=1.0;
constexpr double DEFAULT_LENGTH_Y=1.0;
constexpr int DEFAULT_NX_2D=50;
constexpr int DEFAULT_NY_2D=50;