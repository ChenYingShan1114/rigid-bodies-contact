#include <inertia_matrix.h>
#include <cassert>

//compute inertia matrix and volume by integrating on surfaces
void inertia_matrix(Eigen::Matrix3d &I, Eigen::Vector3d & center, double &mass, Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> F, double density) {
    mass = 0;
    center.setZero();
    I.setZero();

    double center_x = 0, center_y = 0, center_z = 0;

    double P1 = 0, P2 = 0, P3 = 0, P4 = 0, P5 = 0, P6 = 0, P7 = 0, P8 = 0, P9 = 0;
    double I_xy = 0, I_yz = 0, I_zx = 0, I_xx = 0, I_yy = 0, I_zz = 0;

    for (int i = 0; i < F.rows(); i++) {

        double x0 = V(F(i, 0), 0), y0 = V(F(i, 0), 1), z0 = V(F(i, 0), 2);
        double x1 = V(F(i, 1), 0), y1 = V(F(i, 1), 1), z1 = V(F(i, 1), 2);
        double x2 = V(F(i, 2), 0), y2 = V(F(i, 2), 1), z2 = V(F(i, 2), 2);
        double dx1 = x1 - x0, dy1 = y1 - y0, dz1 = z1 - z0;
        double dx2 = x2 - x0, dy2 = y2 - y0, dz2 = z2 - z0;
        
        // normal vector and area
        Eigen::Vector3d dX1 = Eigen::Vector3d::Zero();
        dX1 << dx1, dy1, dz1;
        Eigen::Vector3d dX2 = Eigen::Vector3d::Zero();
        dX2 << dx2, dy2, dz2;
        Eigen::Vector3d n = Eigen::Vector3d::Zero();
        n = dX1.cross(dX2);
        double area = 0.5 * n.norm();
        n.normalize();

        // mass
        mass += 2 * area * density * n(0) / 6 * (x0 + x1 + x2);
    
        // center of mass
        center_x += density / 12.0 * n(0) * area * (x0 * x0 + x1 * x1 + x2 * x2 + x0 * x1 + x1 * x2 + x2 * x0);
        center_y += density / 12.0 * n(1) * area * (y0 * y0 + y1 * y1 + y2 * y2 + y0 * y1 + y1 * y2 + y2 * y0);
        center_z += density / 12.0 * n(2) * area * (z0 * z0 + z1 * z1 + z2 * z2 + z0 * z1 + z1 * z2 + z2 * z0);

        // inertia matrix        
        // I_xy integral by z; I_yz integral by x; I_zx integral by y
        P1 = (x0 * dy1 + y0 * dx1) * z0 + x0 * y0 * dz1, P2 = (x0 * dy2 + y0 * dx2) * z0 + x0 * y0 * dz2;  P3 = dx1 * dy1 * z0 + dx1 * y0 * dz1 + x0 * dy1 * dz1; P4 = dx2 * dy2 * z0 + dx2 * y0 * dz2 + x0 * dy2 * dz2; P5 = (x0 * dy1 + y0 * dx1) * dz2 + (x0 * dy2 + y0 * dx2) * dz1 + (dx1 * dy2 + dy1 * dx2) * z0; P6 = dx1 * dy1 * dz1; P7 = dx2 * dy2 * dz2; P8 = dx1 * dy1 * dz2 + dx1 * dy2 * dz1 + dx2 * dy1 * dz1; P9 = dx2 * dy2 * dz1 + dx2 * dy1 * dz2 + dx1 * dy2 * dz2;
        I_xy += -density * n(2) * area * (x0 * y0 * z0 + (P1 + P2) / 3.0 + (P3 + P4) / 6.0 + P5 / 12.0 + (P6 + P7) / 10.0 + (P8 + P9) / 30.0);
        I_yz += -density * n(0) * area * (x0 * y0 * z0 + (P1 + P2) / 3.0 + (P3 + P4) / 6.0 + P5 / 12.0 + (P6 + P7) / 10.0 + (P8 + P9) / 30.0);
        I_zx += -density * n(1) * area * (x0 * y0 * z0 + (P1 + P2) / 3.0 + (P3 + P4) / 6.0 + P5 / 12.0 + (P6 + P7) / 10.0 + (P8 + P9) / 30.0);

        // I_xx integral by x
        P1 = 3 * dx1 * x0 * x0; P2 = 3 * dx2 * x0 * x0; P3 = 3 * x0 * dx1 * dx1; P4 = 3 * x0 * dx2 * dx2; P5 = 6 * x0 * dx1 * dx2; P6 = dx1 * dx1 * dx1; P7 = dx2 * dx2 * dx2; P8 = 3 * dx1 * dx1 * dx2; P9 = 3 * dx1 * dx2 * dx2;
        I_xx += density / 3.0 * n(0) * area * (x0 * x0 * x0 + (P1 + P2) / 3.0 + (P3 + P4) / 6.0 + P5 / 12.0 + (P6 + P7) / 10.0 + (P8 + P9) / 30.0);

        // I_yy integral by y
        P1 = 3 * dy1 * y0 * y0; P2 = 3 * dy2 * y0 * y0; P3 = 3 * y0 * dy1 * dy1; P4 = 3 * y0 * dy2 * dy2; P5 = 6 * y0 * dy1 * dy2; P6 = dy1 * dy1 * dy1; P7 = dy2 * dy2 * dy2; P8 = 3 * dy1 * dy1 * dy2; P9 = 3 * dy1 * dy2 * dy2;
        I_yy += density / 3.0 * n(1) * area * (y0 * y0 * y0 + (P1 + P2) / 3.0 + (P3 + P4) / 6.0 + P5 / 12.0 + (P6 + P7) / 10.0 + (P8 + P9) / 30.0);
        
        // I_zz integral by z
        P1 = 3 * dz1 * z0 * z0; P2 = 3 * dz2 * z0 * z0; P3 = 3 * z0 * dz1 * dz1; P4 = 3 * z0 * dz2 * dz2; P5 = 6 * z0 * dz1 * dz2; P6 = dz1 * dz1 * dz1; P7 = dz2 * dz2 * dz2; P8 = 3 * dz1 * dz1 * dz2; P9 = 3 * dz1 * dz2 * dz2;
        I_zz += density / 3.0 * n(2) * area * (z0 * z0 * z0 + (P1 + P2) / 3.0 + (P3 + P4) / 6.0 + P5 / 12.0 + (P6 + P7) / 10.0 + (P8 + P9) / 30.0);
        
    }

    center << center_x, center_y, center_z;
    center = center / mass;
    I << I_yy + I_zz,        I_xy,        I_zx,
                I_xy, I_xx + I_zz,        I_yz,
                I_zx,        I_yz, I_xx + I_yy;
    
}