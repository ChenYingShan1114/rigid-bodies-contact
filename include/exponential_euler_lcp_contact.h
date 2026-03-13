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
    for (int irb = 0; irb < masses.size(); irb++) {

        // old time step (t) data
        // inertia and mass
        Eigen::Matrix3d I = masses[irb].block(0, 0, 3, 3);
        Eigen::Matrix3d M = masses[irb].block(3, 3, 3, 3);

        // rotation matrix and center
        Eigen::Matrix3d R = Eigen::Map<const Eigen::Matrix3d>(q.segment<9>(12 * irb).data());
        Eigen::Matrix3d R2 = R * I * R.transpose();
        Eigen::Vector3d p = Eigen::Map<const Eigen::Vector3d>(q.segment<3>(12 * irb + 9).data());

        Eigen::Matrix66d masses_inv = Eigen::Matrix66d::Zero();
        masses_inv.block(0, 0, 3, 3) = R2.inverse();
        masses_inv.block(3, 3, 3, 3) = M.inverse();

        // anguler velocity
        Eigen::Vector3d omega = Eigen::Map<const Eigen::Vector3d>(qdot.segment<3>(6 * irb).data());

        // torque and force
        Eigen::Vector6d force_all = forces.segment(6 * irb, 6);
        force_all.segment(6 * irb, 3) += omega.cross(R2 * omega);
        
        // angular velocity and velocity with no contact
        qdot.segment(6 * irb, 6) += dt * masses_inv * force_all;

        // find contact points of object irb (this part not test yet)
        std::vector<Eigen::Vector3d> x_irb;
        std::vector<Eigen::Vector3d> n_irb; 
        for (int i = 0; i < obj.size(); i++) {
            if (obj[i].first == irb) {
                x_irb.push_back(x[i]);
                n_irb.push_back(n[i]);
            } else if (obj[i].second == irb) {
                x_irb.push_back(x[i]);
                n_irb.push_back(-n[i]);
            }
        }

        // projected Gauss-Seidel
        Eigen::VectorXd alpha = Eigen::VectorXd::Zero(x_irb.size());
        Eigen::Matrix36d J;
        Eigen::Vector3d x_body;
        Eigen::Vector6d fA;
        std::vector<Eigen::Vector6d> gAs(x_irb.size());
        std::vector<double> deltas(x_irb.size());

        for (int i = 0; i < x_irb.size(); i++) {
            inverse_rigid_body(x_body, x_irb[i], R, p);
            rigid_body_jacobian(J, R, p, x_body);
            gAs[i] = J.transpose() * n_irb[i];
            deltas[i] = dt * gAs[i].transpose() * masses_inv * gAs[i];
        }

        for (int iter = 0; iter < 1; iter++) {  
            for (int i = 0; i < x_irb.size(); i++) {
                fA.setZero();
                for (int j = 0; j < x_irb.size(); j++) {
                    if (i != j) {
                        fA = fA + alpha(j) * gAs[j];
                    }
                }
                double gamma = gAs[i].transpose() * (qdot + dt * masses_inv * fA);
                alpha(i) = std::max(0.0, -gamma / deltas[i]);
                qdot = qdot + dt * masses_inv * alpha(i) * gAs[i];

                // More efficient way?
                // double gamma = gAs[i].transpose() * qdot;
                // double alpha_old = alpha(i);
                // alpha(i) = std::max(0.0, alpha_old - gamma / deltas[i]);
                // qdot = qdot + dt * masses_inv * (alpha(i) - alpha_old) * gAs[i];
            }
        }

        // calculate new time step (t + 1) value
        Eigen::Vector3d omega_new = Eigen::Map<const Eigen::Vector3d>(qdot.segment<3>(6 * irb).data());
        Eigen::Vector3d velocity_new = Eigen::Map<const Eigen::Vector3d>(qdot.segment<3>(6 * irb + 3).data());

        // rodrigues rotation matrix
        Eigen::Matrix3d R_tmp = Eigen::Matrix3d::Zero();
        rodrigues(R_tmp, omega_new, dt);

        Eigen::Matrix3d R_new = R_tmp * R;
        Eigen::Vector3d p_new = p + dt * velocity_new;

        // put new value back to generalized coordinate
        q.segment(12 * irb    , 9) = Eigen::Map<Eigen::Matrix<double, 1, 9> >(R_new.data());
        q.segment(12 * irb + 9, 3) = p_new;
        
    }
}