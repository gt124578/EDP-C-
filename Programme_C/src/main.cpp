//sudo apt-get install libboost-all-dev gnuplot

//g++ -I../include main.cpp boundary.cpp solver.cpp source_term.cpp solver2D.cpp Mesh2D.cpp Mesh1D.cpp -o programme -lboost_iostreams -lboost_system -std=c++17

#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <numeric>

#include "parameters.hpp"
#include "boundary.hpp"
#include "source_term.hpp"
#include "Mesh1D.hpp"
#include "solver.hpp"
#include "Mesh2D.hpp"
#include "solver2D.hpp"
#include "gnuplot-iostream.h"

using namespace std;

vector<vector<double>> transpose(const vector<vector<double>>& data) {
    if (data.empty()) return {};
    size_t rows=data.size();
    size_t cols=data[0].size();
    vector<vector<double>> transposed_data(cols, vector<double>(rows));
    for (size_t i=0; i<cols; ++i) {
        for (size_t j=0; j<rows; ++j) {
            transposed_data[i][j]=data[j][i];
        }
    }
    return transposed_data;
}




void run_simulation_1d() {
    cout<<"\n=================================================="<<endl;
    cout<<"          LANCEMENT DE LA SIMULATION 1D"<<endl;
    cout<<"=======================================================\n"<<endl;

    //1-Création du maillage et des paramètres
    Mesh1D mesh(DEFAULT_LENGTH,DEFAULT_NX);
    const double dt=mesh.dx()*mesh.dx()/(2.0*D_COEFF);
    const int nt=static_cast<int>(SIMULATION_TIME/dt);
    vector<vector<double>> u(nt,vector<double>(mesh.size(),0.0));

    //2-Condition initiale
    for (int i=0;i<mesh.size();++i) {
        u[0][i]=sin(M_PI*mesh.x()[i]);
    }

    //3-Résolution
    cout<<"Calcul en cours... (nt="<<nt<<", nx="<<mesh.size()<<")"<<endl;
    solve_heat_equation(u,mesh,nt,dt);
    cout<<"Simulation 1D terminee."<<endl;

    //4-Visualisation
    static Gnuplot gp1,gp2;

    //Visualisation 1D (Courbe finale)
    cout<<"Affichage du graphique 1D..."<<endl;
    gp1<<"set title 'Equation de la chaleur 1D (t="<<SIMULATION_TIME<<")'\n";
    gp1<<"set xlabel 'Position x'\n";
    gp1<<"set ylabel 'u(x, t)'\n";
    gp1<<"set grid\n";
    gp1<<"plot '-' with lines title 'Solution 1D'\n";
    gp1.send1d(make_tuple(mesh.x(),u.back()));
    gp1.flush();

    //Diagramme Spatio-Temporel 1D
    cout<<"Affichage de la heatmap 1D (x,t)..."<<endl;
    gp2<<"set title 'Evolution Spatio-Temporelle 1D'\n";
    gp2<<"set xlabel 'Position x'\n";
    gp2<<"set ylabel 'Temps t'\n";
    gp2<<"set xrange [0:1]\n";
    gp2<<"set yrange [0:2.1]\n";
    gp2<<"set cbrange [0:3.0]\n";
    gp2<<"set palette defined ( 0 'black', 0.25 'dark-red', 0.5 'red', 0.75 'yellow', 1 'white' )\n";
    gp2<<"set colorbox\n";
    double pixel_w=DEFAULT_LENGTH/mesh.size();
    double pixel_h=SIMULATION_TIME/nt;
    gp2<<"plot '-' binary array=("<<mesh.size()<<","<<nt<<") format='%double' dx="<<pixel_w<<" dy="<<pixel_h<<" with image notitle\n";
    vector<vector<double>> transposed_u=transpose(u);
    gp2.sendBinary(transposed_u);
    gp2.flush();
}


void run_simulation_2d() {
    cout<<"\n================================================"<<endl;
    cout<<"          LANCEMENT DE LA SIMULATION 2D"<<endl;
    cout<<"================================================\n"<<endl;

    //1-Création du maillage et des paramètres
    Mesh2D mesh(DEFAULT_LENGTH_X,DEFAULT_LENGTH_Y,DEFAULT_NX_2D,DEFAULT_NY_2D);
    double dt_2d=min(mesh.dx(),mesh.dy())*min(mesh.dx(),mesh.dy())/(4.0*D_COEFF);
    int nt_2d=static_cast<int>(SIMULATION_TIME/dt_2d);
    Solver2D::SolutionHistory u(nt_2d,vector<vector<double>>(mesh.size_x(),vector<double>(mesh.size_y(),0.0)));

    //2-Condition initiale
    for (int i=0;i<mesh.size_x();++i) {
        for (int j=0;j<mesh.size_y();++j) {
            u[0][i][j]=sin(M_PI*mesh.x()[i])*sin(M_PI*mesh.y()[j]);
        }
    }

    //3-Résolution
    cout<<"Calcul en cours... (nt="<<nt_2d<<", nx="<<mesh.size_x()<<", ny="<<mesh.size_y()<<")"<<endl;
    Solver2D solver;
    solver.solve(u,mesh,nt_2d,dt_2d);
    cout<<"Simulation 2D terminee."<<endl;

    //4-Visualisation
    static Gnuplot gp_2d;
    cout<<"Affichage de la heatmap 2D (x,y)..."<<endl;
    gp_2d<<"set title 'Equation de la chaleur 2D (t="<<SIMULATION_TIME<<")'\n";
    gp_2d<<"set xlabel 'Position x'\n";
    gp_2d<<"set ylabel 'Position y'\n";
    gp_2d<<"set view map\n";
    gp_2d<<"set size ratio -1\n";
    gp_2d<<"set xrange [-1:50]\n";
    gp_2d<<"set yrange [-1:50]\n";
    gp_2d<<"set cbrange [0:3.0]\n";
    gp_2d<<"set palette defined ( 0 'black', 0.25 'dark-red', 0.5 'red', 0.75 'yellow', 1 'white' )\n";
    gp_2d<<"set colorbox\n";
    gp_2d<<"splot '-' matrix with image notitle\n";
    const auto& final_solution_2d=u.back();
    for (int j=0;j<mesh.size_y();++j) {
        string row_data;
        row_data.reserve(mesh.size_x()*10);
        for (int i=0;i<mesh.size_x();++i) {
            row_data+=to_string(final_solution_2d[i][j]);
            if (i<mesh.size_x()-1) row_data+=" ";
        }
        gp_2d<<row_data<<"\n";
    }
    gp_2d<<"e\n";
    gp_2d.flush();
}



int main() {
    locale::global(locale("C"));

    //Lancer la première simulation
    run_simulation_1d();

    //Lancer la seconde simulation
    run_simulation_2d();

    //Pause générale pour garder toutes les fenêtres ouvertes
    cout<<"\n\nToutes les simulations sont terminees."<<endl;
    cout<<"Toutes les fenetres de graphiques sont affichees."<<endl;
    cout<<"Appuyez sur Entree pour quitter le programme."<<endl;
    cin.get();

    return 0;
}