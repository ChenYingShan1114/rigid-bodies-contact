#include <Eigen/Dense>
#include <EigenTypes.h>
#include <rodrigues.h>
#include <iostream>
#include <rigid_body_jacobian.h>
#include <inverse_rigid_body.h>

//Input:
//  q - 12n vector where n is the number of rigid bodies. Each rigid body is stored as 12 doubles. 
//      The first 9 doubles are the columns of a 3x3 rotation matrix and the final 3 doubles are the world space position of the object's center of mass.
//  qdot - 6n vector of generalied velocities. The first 3 doubles of each body are the world space angular velocity and 
//         the second 3 are the world space linear velocity.
//  dt - the integration time step
//  masses - a vector to mass matrices for each rigid body
//  forces - a 6n vector of generalized forces for n rigid bodies. The first 3 doubles of each rigid body are the torques acting on the object
//           while the second 3 doubles are the linear forces.
//  n - list of collision normals
//  x - list of world space collision points
//  obj - list of collision object ids 
//Output:
//  q - updated generalized coordinates 
//  qdot - updated generalized velocities 
inline void exponential_euler_lcp_contact(Eigen::VectorXd &q, Eigen::VectorXd &qdot, double dt, 
                            std::vector<Eigen::Matrix66d> &masses, Eigen::Ref<const Eigen::VectorXd> forces,
                            std::vector<Eigen::Vector3d> &n, std::vector<Eigen::Vector3d> &x, std::vector<std::pair<int,int> > &obj) {
    for (int i = 0; i < masses.size(); i++) {

        // old time step (t) data
        // inertia and mass
        Eigen::Matrix3d I = Eigen::Matrix3d::Zero();
        I = masses[i].block(0, 0, 3, 3);
        double mass = masses[i](3, 3);

        // torque and force
        Eigen::Vector3d torque = Eigen::Map<const Eigen::Vector3d>(forces.segment<3>(6 * i).data());
        Eigen::Vector3d force = Eigen::Map<const Eigen::Vector3d>(forces.segment<3>(6 * i + 3).data());
        
        // rotation matrix and center
        Eigen::Matrix3d R = Eigen::Map<const Eigen::Matrix3d>(q.segment<9>(12 * i).data());
        Eigen::Matrix3d R2 = R * I * R.transpose();
        Eigen::Vector3d center = Eigen::Map<const Eigen::Vector3d>(q.segment<3>(12 * i + 9).data());

        // anguler velocity and velocity
        Eigen::Vector3d omega = Eigen::Map<const Eigen::Vector3d>(qdot.segment<3>(6 * i).data());
        Eigen::Matrix3d omega_cross = Eigen::Matrix3d::Zero();
        omega_cross << 0, -omega(2), omega(1),
                    omega(2), 0, -omega(0),
                    -omega(1), omega(0), 0;
        Eigen::Vector3d velocity = Eigen::Map<const Eigen::Vector3d>(qdot.segment<3>(6 * i + 3).data());

        // rodrigues rotation matrix
        Eigen::Matrix3d R_tmp = Eigen::Matrix3d::Zero();
        rodrigues(R_tmp, omega, dt);

        // calculate new time step (t + 1) value
        Eigen::Vector3d omega_new = R2.inverse() * (R2 * omega + dt * omega_cross * (R2 * omega) + dt * torque);
        Eigen::Vector3d velocity_new = velocity + dt * force / mass;
        Eigen::Matrix3d R_new = R_tmp * R;
        Eigen::Vector3d center_new = center + dt * velocity;

        // put new value back to generalized coordinate
        qdot.segment(6 * i    , 3) = omega_new;
        qdot.segment(6 * i + 3, 3) = velocity_new;
        q.segment(12 * i    , 9) = Eigen::Map<Eigen::Matrix<double, 1, 9> >(R_new.data());
        q.segment(12 * i + 9, 3) = center_new;

    }
}