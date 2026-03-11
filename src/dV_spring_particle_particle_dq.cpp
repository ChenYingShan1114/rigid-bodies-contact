#include <dV_spring_particle_particle_dq.h>

void dV_spring_particle_particle_dq(Eigen::Ref<Eigen::Vector6d> f, Eigen::Ref<const Eigen::Vector3d> q0,  Eigen::Ref<const Eigen::Vector3d>     q1, double l0, double stiffness) {
    f.setZero();

    double dist = (q1 - q0).norm();
    double force = -stiffness * (dist - l0) / dist;
    f << force * (q1 - q0), force * (q0 - q1);
    
}