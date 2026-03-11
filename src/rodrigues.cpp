#include <rodrigues.h>
#include <cmath>

void rodrigues(Eigen::Matrix3d &R, Eigen::Ref<const Eigen::Vector3d> omega, double dt) {
    R.setZero();

    double theta = omega.norm() * dt;
    
    Eigen::Vector3d omega_nor = omega.normalized();
    Eigen::Matrix3d omega_cross = Eigen::Matrix3d::Zero();
    omega_cross << 0, -omega_nor(2), omega_nor(1),
                   omega_nor(2), 0, -omega_nor(0),
                   -omega_nor(1), omega_nor(0), 0;

    Eigen::Matrix3d I = Eigen::Matrix3d::Zero();
    I.setIdentity();

    R = I + sin(theta) * omega_cross + (1 - cos(theta)) * omega_cross * omega_cross;
}