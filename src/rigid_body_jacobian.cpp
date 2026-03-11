#include <rigid_body_jacobian.h>

void rigid_body_jacobian(Eigen::Matrix36d &J, 
                         Eigen::Ref<const Eigen::Matrix3d> R, Eigen::Ref<const Eigen::Vector3d> p, 
                         Eigen::Ref<const Eigen::Vector3d> x) {

    J.setZero();
           
    Eigen::Vector3d X_bar = x - p;
    Eigen::Matrix3d X_cross = Eigen::Matrix3d::Zero();
    X_cross << 0, -X_bar(2), X_bar(1),
               X_bar(2), 0, -X_bar(0),
               -X_bar(1), X_bar(0), 0;

    Eigen::Matrix3d I = Eigen::Matrix3d::Zero();
    I.setIdentity();

    Eigen::Matrix36d J1 = Eigen::Matrix36d::Zero();
    J1 << X_cross.transpose(), I;

    Eigen::Matrix66d J2 = Eigen::Matrix66d::Zero();
    J2.setIdentity();
    J2.block(0, 0, 3, 3) = R.transpose();
    J2.block(3, 3, 3, 3) = R.transpose();

    J = R * J1 * J2;

}