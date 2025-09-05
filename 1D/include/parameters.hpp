#pragma once

// Paramètres physiques
constexpr double L = 1.0;
constexpr double T_total = 2.0;
constexpr int nx = 100;
constexpr double D = 0.01;
constexpr double C = 0.03;

// Paramètres numériques
constexpr double dx = L / (nx - 1);
constexpr double dt = dx * dx / (2.0 * D);
constexpr int nt = static_cast<int>(T_total / dt);
