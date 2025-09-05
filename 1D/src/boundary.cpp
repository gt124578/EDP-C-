#include "boundary.hpp"

void apply_boundary_conditions(std::vector<double>& u_n) {
    u_n.front() = 0.0;
    u_n.back() = 0.0;
}
