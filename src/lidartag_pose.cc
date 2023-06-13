#include "lidartag.h"
#include "utils.h"

#include <nlopt.hpp>
#include <iostream>

using namespace std;

namespace BipedLab
{
double checkCost(double point, double cons1, double cons2)
{
  if (point >= cons1 && point <= cons2)
    return 0;
  else
  {
    return std::min(std::abs(point - cons2), std::abs(point - cons1));
  }
}
Eigen::VectorXd d_px_euler(double x11, double y11, double z11, double rpy11, double rpy12, double rpy13)
{
  double t2 = (rpy12 * M_PI) / 1.8e+2;
  double t3 = (rpy13 * M_PI) / 1.8e+2;
  double t4 = std::cos(t2);
  double t5 = std::cos(t3);
  double t6 = std::sin(t2);
  double t7 = std::sin(t3);

  Eigen::VectorXd d_px(6);
  d_px << 1.0, 0.0, 0.0, 0.0,
      (z11 * t4 * M_PI) / 1.8e+2 - (x11 * t5 * t6 * M_PI) / 1.8e+2 + (y11 * t6 * t7 * M_PI) / 1.8e+2,
      t4 * M_PI * (x11 * t7 + y11 * t5) * (-1.0 / 1.8e+2);
  return d_px;
}

Eigen::VectorXd d_py_euler(double x11, double y11, double z11, double rpy11, double rpy12, double rpy13)
{
  double t2 = (rpy11 * M_PI) / 1.8e+2;
  double t3 = (rpy12 * M_PI) / 1.8e+2;
  double t4 = (rpy13 * M_PI) / 1.8e+2;
  double t5 = std::cos(t2);
  double t6 = std::cos(t3);
  double t7 = std::cos(t4);
  double t8 = std::sin(t2);
  double t9 = std::sin(t3);
  double t10 = std::sin(t4);

  Eigen::VectorXd d_py(6);
  d_py << 0.0, 1.0, 0.0,
      -x11 * ((t8 * t10 * M_PI) / 1.8e+2 - (t5 * t7 * t9 * M_PI) / 1.8e+2) -
          y11 * ((t7 * t8 * M_PI) / 1.8e+2 + (t5 * t9 * t10 * M_PI) / 1.8e+2) - (z11 * t5 * t6 * M_PI) / 1.8e+2,
      (t8 * M_PI * (z11 * t9 + x11 * t6 * t7 - y11 * t6 * t10)) / 1.8e+2,
      x11 * ((t5 * t7 * M_PI) / 1.8e+2 - (t8 * t9 * t10 * M_PI) / 1.8e+2) -
          y11 * ((t5 * t10 * M_PI) / 1.8e+2 + (t7 * t8 * t9 * M_PI) / 1.8e+2);
  return d_py;
}

Eigen::VectorXd d_pz_euler(double x11, double y11, double z11, double rpy11, double rpy12, double rpy13)
{
  double t2 = (rpy11 * M_PI) / 1.8e+2;
  double t3 = (rpy12 * M_PI) / 1.8e+2;
  double t4 = (rpy13 * M_PI) / 1.8e+2;
  double t5 = std::cos(t2);
  double t6 = std::cos(t3);
  double t7 = std::cos(t4);
  double t8 = std::sin(t2);
  double t9 = std::sin(t3);
  double t10 = std::sin(t4);

  Eigen::VectorXd d_pz(6);
  d_pz << 0.0, 0.0, 1.0,
      x11 * ((t5 * t10 * M_PI) / 1.8e+2 + (t7 * t8 * t9 * M_PI) / 1.8e+2) +
          y11 * ((t5 * t7 * M_PI) / 1.8e+2 - (t8 * t9 * t10 * M_PI) / 1.8e+2) - (z11 * t6 * t8 * M_PI) / 1.8e+2,
      t5 * M_PI * (z11 * t9 + x11 * t6 * t7 - y11 * t6 * t10) * (-1.0 / 1.8e+2),
      x11 * ((t7 * t8 * M_PI) / 1.8e+2 + (t5 * t9 * t10 * M_PI) / 1.8e+2) -
          y11 * ((t8 * t10 * M_PI) / 1.8e+2 - (t5 * t7 * t9 * M_PI) / 1.8e+2);
  return d_pz;
}

void d_px_lie_mat(const Eigen::MatrixXf& X1, const Eigen::MatrixXf& Y1, const Eigen::MatrixXf& Z1, const double& r11,
                  const double& r12, const double& r13, Eigen::MatrixXf& dx_mat)
{
  double t2 = std::abs(r11);
  double t3 = std::abs(r12);
  double t4 = std::abs(r13);
  double t5 = utils::get_sign(r11);
  double t6 = utils::get_sign(r12);
  double t7 = utils::get_sign(r13);
  double t8 = std::pow(r12, 2);
  double t9 = std::pow(r13, 2);
  double t10 = std::pow(t2, 2);
  double t11 = std::pow(t3, 2);
  double t12 = std::pow(t4, 2);
  double t13 = t8 + t9;
  double t14 = t10 + t11 + t12;
  double t15 = 1 / t14;
  double t17 = std::sqrt(t14);
  double t16 = std::pow(t15, 2);
  double t18 = 1 / t17;
  double t20 = std::cos(t17);
  double t21 = std::sin(t17);
  double t19 = std::pow(t18, 3);
  double t22 = t20 - 1;
  double t23 = t18 * t21;
  double t24 = r11 * t15 * t22;
  double t25 = -t24;

  dx_mat.setZero();
  dx_mat.row(0).setOnes();
  dx_mat.row(3) = Y1 * (-r12 * t15 * t22 - r13 * t2 * t5 * t15 * t20 + r13 * t2 * t5 * t19 * t21 +
                        r11 * r12 * t2 * t5 * t16 * t22 * 2 + r11 * r12 * t2 * t5 * t19 * t21) +
                  Z1 * (-r13 * t15 * t22 + r12 * t2 * t5 * t15 * t20 - r12 * t2 * t5 * t19 * t21 +
                        r11 * r13 * t2 * t5 * t16 * t22 * 2 + r11 * r13 * t2 * t5 * t19 * t21) -
                  X1 * (t2 * t5 * t13 * t16 * t22 * 2 + t2 * t5 * t13 * t19 * t21);

  dx_mat.row(4) = -X1 * (r12 * t15 * t22 * (-2) + t3 * t6 * t13 * t16 * t22 * 2 + t3 * t6 * t13 * t19 * t21) +
                  Y1 * (t25 - r13 * t3 * t6 * t15 * t20 + r13 * t3 * t6 * t19 * t21 +
                        r11 * r12 * t3 * t6 * t16 * t22 * 2 + r11 * r12 * t3 * t6 * t19 * t21) +
                  Z1 * (t23 + r12 * t3 * t6 * t15 * t20 - r12 * t3 * t6 * t19 * t21 +
                        r11 * r13 * t3 * t6 * t16 * t22 * 2 + r11 * r13 * t3 * t6 * t19 * t21);

  dx_mat.row(5) = -X1 * (r13 * t15 * t22 * (-2) + t4 * t7 * t13 * t16 * t22 * 2 + t4 * t7 * t13 * t19 * t21) +
                  Y1 * (-t23 - r13 * t4 * t7 * t15 * t20 + r13 * t4 * t7 * t19 * t21 +
                        r11 * r12 * t4 * t7 * t16 * t22 * 2 + r11 * r12 * t4 * t7 * t19 * t21) +
                  Z1 * (t25 + r12 * t4 * t7 * t15 * t20 - r12 * t4 * t7 * t19 * t21 +
                        r11 * r13 * t4 * t7 * t16 * t22 * 2 + r11 * r13 * t4 * t7 * t19 * t21);
}

void d_py_lie_mat(const Eigen::MatrixXf& X1, const Eigen::MatrixXf& Y1, const Eigen::MatrixXf& Z1, const double& r11,
                  const double& r12, const double& r13, Eigen::MatrixXf& dy_mat)
{
  double t2 = std::abs(r11);
  double t3 = std::abs(r12);
  double t4 = std::abs(r13);
  double t5 = utils::get_sign(r11);
  double t6 = utils::get_sign(r12);
  double t7 = utils::get_sign(r13);
  double t8 = std::pow(r11, 2);
  double t9 = std::pow(r13, 2);
  double t10 = std::pow(t2, 2);
  double t11 = std::pow(t3, 2);
  double t12 = std::pow(t4, 2);
  double t13 = t8 + t9;
  double t14 = t10 + t11 + t12;
  double t15 = 1 / t14;
  double t17 = std::sqrt(t14);
  double t16 = std::pow(t15, 2);
  double t18 = 1 / t17;
  double t20 = std::cos(t17);
  double t21 = std::sin(t17);
  double t19 = std::pow(t18, 3);
  double t22 = t20 - 1;
  double t23 = t18 * t21;
  double t24 = r12 * t15 * t22;
  double t25 = -t24;

  dy_mat.setZero();
  dy_mat.row(1).setOnes();
  dy_mat.row(3) = -Y1 * (r11 * t15 * t22 * (-2) + t2 * t5 * t13 * t16 * t22 * 2 + t2 * t5 * t13 * t19 * t21) +
                  Z1 * (-t23 - r11 * t2 * t5 * t15 * t20 + r11 * t2 * t5 * t19 * t21 +
                        r12 * r13 * t2 * t5 * t16 * t22 * 2 + r12 * r13 * t2 * t5 * t19 * t21) +
                  X1 * (t25 + r13 * t2 * t5 * t15 * t20 - r13 * t2 * t5 * t19 * t21 +
                        r11 * r12 * t2 * t5 * t16 * t22 * 2 + r11 * r12 * t2 * t5 * t19 * t21);

  dy_mat.row(4) = X1 * (-r11 * t15 * t22 + r13 * t3 * t6 * t15 * t20 - r13 * t3 * t6 * t19 * t21 +
                        r11 * r12 * t3 * t6 * t16 * t22 * 2 + r11 * r12 * t3 * t6 * t19 * t21) +
                  Z1 * (-r13 * t15 * t22 - r11 * t3 * t6 * t15 * t20 + r11 * t3 * t6 * t19 * t21 +
                        r12 * r13 * t3 * t6 * t16 * t22 * 2 + r12 * r13 * t3 * t6 * t19 * t21) -
                  Y1 * (t3 * t6 * t13 * t16 * t22 * 2 + t3 * t6 * t13 * t19 * t21);

  dy_mat.row(5) = -Y1 * (r13 * t15 * t22 * (-2) + t4 * t7 * t13 * t16 * t22 * 2 + t4 * t7 * t13 * t19 * t21) +
                  X1 * (t23 + r13 * t4 * t7 * t15 * t20 - r13 * t4 * t7 * t19 * t21 +
                        r11 * r12 * t4 * t7 * t16 * t22 * 2 + r11 * r12 * t4 * t7 * t19 * t21) +
                  Z1 * (t25 - r11 * t4 * t7 * t15 * t20 + r11 * t4 * t7 * t19 * t21 +
                        r12 * r13 * t4 * t7 * t16 * t22 * 2 + r12 * r13 * t4 * t7 * t19 * t21);
}

void d_pz_lie_mat(const Eigen::MatrixXf& X1, const Eigen::MatrixXf& Y1, const Eigen::MatrixXf& Z1, const double& r11,
                  const double& r12, const double& r13, Eigen::MatrixXf& dz_mat)
{
  double t2 = std::abs(r11);
  double t3 = std::abs(r12);
  double t4 = std::abs(r13);
  double t5 = utils::get_sign(r11);
  double t6 = utils::get_sign(r12);
  double t7 = utils::get_sign(r13);
  double t8 = std::pow(r11, 2);
  double t9 = std::pow(r12, 2);
  double t10 = std::pow(t2, 2);
  double t11 = std::pow(t3, 2);
  double t12 = std::pow(t4, 2);
  double t13 = t8 + t9;
  double t14 = t10 + t11 + t12;
  double t15 = 1 / t14;
  double t17 = std::sqrt(t14);
  double t16 = std::pow(t15, 2);
  double t18 = 1 / t17;
  double t20 = std::cos(t17);
  double t21 = std::sin(t17);
  double t19 = std::pow(t18, 3);
  double t22 = t20 - 1;
  double t23 = t18 * t21;
  double t24 = r13 * t15 * t22;
  double t25 = -t24;

  dz_mat.setZero();
  dz_mat.row(2).setOnes();
  dz_mat.row(3) = -Z1 * (r11 * t15 * t22 * (-2) + t2 * t5 * t13 * t16 * t22 * 2 + t2 * t5 * t13 * t19 * t21) +
                  X1 * (t25 - r12 * t2 * t5 * t15 * t20 + r12 * t2 * t5 * t19 * t21 +
                        r11 * r13 * t2 * t5 * t16 * t22 * 2 + r11 * r13 * t2 * t5 * t19 * t21) +
                  Y1 * (t23 + r11 * t2 * t5 * t15 * t20 - r11 * t2 * t5 * t19 * t21 +
                        r12 * r13 * t2 * t5 * t16 * t22 * 2 + r12 * r13 * t2 * t5 * t19 * t21);

  dz_mat.row(4) = -Z1 * (r12 * t15 * t22 * (-2) + t3 * t6 * t13 * t16 * t22 * 2 + t3 * t6 * t13 * t19 * t21) +
                  X1 * (-t23 - r12 * t3 * t6 * t15 * t20 + r12 * t3 * t6 * t19 * t21 +
                        r11 * r13 * t3 * t6 * t16 * t22 * 2 + r11 * r13 * t3 * t6 * t19 * t21) +
                  Y1 * (t25 + r11 * t3 * t6 * t15 * t20 - r11 * t3 * t6 * t19 * t21 +
                        r12 * r13 * t3 * t6 * t16 * t22 * 2 + r12 * r13 * t3 * t6 * t19 * t21);
  dz_mat.row(5) = X1 * (-r11 * t15 * t22 - r12 * t4 * t7 * t15 * t20 + r12 * t4 * t7 * t19 * t21 +
                        r11 * r13 * t4 * t7 * t16 * t22 * 2 + r11 * r13 * t4 * t7 * t19 * t21) +
                  Y1 * (-r12 * t15 * t22 + r11 * t4 * t7 * t15 * t20 - r11 * t4 * t7 * t19 * t21 +
                        r12 * r13 * t4 * t7 * t16 * t22 * 2 + r12 * r13 * t4 * t7 * t19 * t21) -
                  Z1 * (t4 * t7 * t13 * t16 * t22 * 2 + t4 * t7 * t13 * t19 * t21);
}

void d_px_euler_mat(const Eigen::MatrixXf& x11, const Eigen::MatrixXf& y11, const Eigen::MatrixXf& z11,
                    const double& rpy11, const double& rpy12, const double& rpy13, Eigen::MatrixXf& dx_mat)
{
  double t4 = std::cos(rpy12);
  double t5 = std::cos(rpy13);
  double t6 = std::sin(rpy12);
  double t7 = std::sin(rpy13);

  dx_mat.setZero();
  dx_mat.row(0).setOnes();
  dx_mat.row(4) = (z11 * t4 - x11 * t5 * t6 + y11 * t6 * t7);
  dx_mat.row(5) = -t4 * (x11 * t7 + y11 * t5);
}

void d_py_euler_mat(const Eigen::MatrixXf& x11, const Eigen::MatrixXf& y11, const Eigen::MatrixXf& z11,
                    const double& rpy11, const double& rpy12, const double& rpy13, Eigen::MatrixXf& dy_mat)
{
  double t5 = std::cos(rpy11);
  double t6 = std::cos(rpy12);
  double t7 = std::cos(rpy13);
  double t8 = std::sin(rpy11);
  double t9 = std::sin(rpy12);
  double t10 = std::sin(rpy13);

  dy_mat.setZero();
  dy_mat.row(1).setOnes();
  dy_mat.row(3) = (M_PI / 180) * (-x11 * (t8 * t10 - t5 * t7 * t9) - y11 * (t7 * t8 + t5 * t9 * t10) - z11 * t5 * t6);
  dy_mat.row(4) = t8 * (z11 * t9 + x11 * t6 * t7 - y11 * t6 * t10);
  dy_mat.row(5) = (x11 * (t5 * t7 - t8 * t9 * t10) - y11 * (t5 * t10 + t7 * t8 * t9));
}

void d_pz_euler_mat(const Eigen::MatrixXf& x11, const Eigen::MatrixXf& y11, const Eigen::MatrixXf& z11,
                    const double& rpy11, const double& rpy12, const double& rpy13, Eigen::MatrixXf& dz_mat)
{
  double t5 = std::cos(rpy11);
  double t6 = std::cos(rpy12);
  double t7 = std::cos(rpy13);
  double t8 = std::sin(rpy11);
  double t9 = std::sin(rpy12);
  double t10 = std::sin(rpy13);
  dz_mat.setZero();
  dz_mat.row(2).setOnes();
  dz_mat.row(3) = (x11 * (t5 * t10 + t7 * t8 * t9) + y11 * (t5 * t7 - t8 * t9 * t10) - z11 * t6 * t8);
  dz_mat.row(4) = -t5 * (z11 * t9 + x11 * t6 * t7 - y11 * t6 * t10);
  dz_mat.row(5) = (x11 * (t7 * t8 + t5 * t9 * t10) - y11 * (t8 * t10 - t5 * t7 * t9));
}

double getCost(const Eigen::Vector4f& template_bound, Eigen::MatrixXf& transformed_points)
{
  Eigen::MatrixXf transfomed_cost_ary = (transformed_points.cwiseAbs()).colwise() - template_bound;
  Eigen::Vector3f cost_vec = transfomed_cost_ary.cwiseMax(0).rowwise().sum();
  double cost = cost_vec[0] + cost_vec[1] + cost_vec[2];

  return cost;
}

double evaluateCost(const Eigen::Matrix4f& H, const Eigen::MatrixXf& points, const Eigen::Vector4f& template_bound,
                    Eigen::MatrixXf& transformed_points)
{
  transformed_points = H * points;  // 4 x n
  return getCost(template_bound, transformed_points);
}

double computeCost_lie(const std::vector<double>& x, std::vector<double>& grad, void* func_data)
{
  Eigen::MatrixXf* d = reinterpret_cast<Eigen::MatrixXf*>(func_data);
  int num_points = d->cols() - 1;  // last column is passed as boundary

  const double box_width = (*d)(0, num_points);
  const double tag_size = (*d)(0, num_points) * 2;
  Eigen::Vector3d q(x[3], x[4], x[5]);
  Eigen::Matrix4f H = Eigen::Matrix4f::Identity(4, 4);
  H.topLeftCorner(3, 3) = utils::Exp_SO3(q);
  H.topRightCorner(3, 1) << x[0], x[1], x[2];

  Eigen::MatrixXf transfomed_points_mat;

  double cost = evaluateCost(H, (*d).topLeftCorner(4, num_points), (*d).col(num_points), transfomed_points_mat);

  Eigen::Matrix<float, 1, Eigen::Dynamic> transformed_x_ind =
      (transfomed_points_mat.row(0).array() > box_width).cast<float>() -
      (transfomed_points_mat.row(0).array() < -box_width).cast<float>();  // 1 x n

  Eigen::Matrix<float, 1, Eigen::Dynamic> transformed_y_ind =
      (transfomed_points_mat.row(1).array() > tag_size / 2).cast<float>() -
      (transfomed_points_mat.row(1).array() < -tag_size / 2).cast<float>();  // 1 x n

  Eigen::Matrix<float, 1, Eigen::Dynamic> transformed_z_ind =
      (transfomed_points_mat.row(2).array() > tag_size / 2).cast<float>() -
      (transfomed_points_mat.row(2).array() < -tag_size / 2).cast<float>();  // 1 x n

  Eigen::Matrix<float, 1, Eigen::Dynamic> transformed_points_ind(1, 3 * (num_points));  // 1 x 3n
  transformed_points_ind << transformed_x_ind, transformed_y_ind, transformed_z_ind;

  Eigen::VectorXf grad_eig = Eigen::VectorXf::Zero(6);
  if (!grad.empty())
  {
    Eigen::MatrixXf dx_mat = Eigen::MatrixXf::Zero(6, num_points);
    d_px_lie_mat(transfomed_points_mat.row(0), transfomed_points_mat.row(1), transfomed_points_mat.row(2), x[3], x[4],
                 x[5],
                 dx_mat);  // 6 x n

    Eigen::MatrixXf dy_mat = Eigen::MatrixXf::Zero(6, num_points);
    d_py_lie_mat(transfomed_points_mat.row(0), transfomed_points_mat.row(1), transfomed_points_mat.row(2), x[3], x[4],
                 x[5], dy_mat);

    Eigen::MatrixXf dz_mat = Eigen::MatrixXf::Zero(6, num_points);
    d_pz_lie_mat(transfomed_points_mat.row(0), transfomed_points_mat.row(1), transfomed_points_mat.row(2), x[3], x[4],
                 x[5], dz_mat);

    Eigen::MatrixXf d_mat(6, 3 * (num_points));
    d_mat << dx_mat, dy_mat, dz_mat;

    grad_eig = d_mat * transformed_points_ind.transpose();
    grad[0] = grad_eig[0];
    grad[1] = grad_eig[1];
    grad[2] = grad_eig[2];
    grad[3] = grad_eig[3];
    grad[4] = grad_eig[4];
    grad[5] = grad_eig[5];
  }
  return cost;
}

double computeCost_euler(const std::vector<double>& x, std::vector<double>& grad, void* func_data)
{
  Eigen::MatrixXf* d = reinterpret_cast<Eigen::MatrixXf*>(func_data);
  int num_points = d->cols() - 1;  // last column is passed as boundary

  const double box_width = (*d)(0, num_points);
  const double tag_size = (*d)(0, num_points) * 2;

  Eigen::Quaternion<float> q = Eigen::AngleAxisf(x[3], Eigen::Vector3f::UnitX()) *
                               Eigen::AngleAxisf(x[4], Eigen::Vector3f::UnitY()) *
                               Eigen::AngleAxisf(x[5], Eigen::Vector3f::UnitZ());
  Eigen::Matrix4f H = Eigen::Matrix4f::Identity(4, 4);
  H.topLeftCorner(3, 3) = q.matrix();
  H.topRightCorner(3, 1) << x[0], x[1], x[2];

  Eigen::MatrixXf transfomed_points_mat;

  double cost = evaluateCost(H, (*d).topLeftCorner(4, num_points), (*d).col(num_points), transfomed_points_mat);

  Eigen::Matrix<float, 1, Eigen::Dynamic> transformed_x_ind =
      (transfomed_points_mat.row(0).array() > box_width).cast<float>() -
      (transfomed_points_mat.row(0).array() < -box_width).cast<float>();  // 1 x n

  Eigen::Matrix<float, 1, Eigen::Dynamic> transformed_y_ind =
      (transfomed_points_mat.row(1).array() > tag_size / 2).cast<float>() -
      (transfomed_points_mat.row(1).array() < -tag_size / 2).cast<float>();  // 1 x n

  Eigen::Matrix<float, 1, Eigen::Dynamic> transformed_z_ind =
      (transfomed_points_mat.row(2).array() > tag_size / 2).cast<float>() -
      (transfomed_points_mat.row(2).array() < -tag_size / 2).cast<float>();  // 1 x n

  Eigen::Matrix<float, 1, Eigen::Dynamic> transformed_points_ind(1, 3 * (num_points));  // 1 x 3n
  transformed_points_ind << transformed_x_ind, transformed_y_ind, transformed_z_ind;

  Eigen::VectorXf grad_eig = Eigen::VectorXf::Zero(6);
  if (!grad.empty())
  {
    Eigen::MatrixXf dx_mat = Eigen::MatrixXf::Zero(6, num_points);
    d_px_euler_mat(transfomed_points_mat.row(0), transfomed_points_mat.row(1), transfomed_points_mat.row(2), x[3], x[4],
                   x[5],
                   dx_mat);  // 6 x n

    Eigen::MatrixXf dy_mat = Eigen::MatrixXf::Zero(6, num_points);
    d_py_euler_mat(transfomed_points_mat.row(0), transfomed_points_mat.row(1), transfomed_points_mat.row(2), x[3], x[4],
                   x[5], dy_mat);

    Eigen::MatrixXf dz_mat = Eigen::MatrixXf::Zero(6, num_points);
    d_pz_euler_mat(transfomed_points_mat.row(0), transfomed_points_mat.row(1), transfomed_points_mat.row(2), x[3], x[4],
                   x[5], dz_mat);

    Eigen::MatrixXf d_mat(6, 3 * (num_points));
    d_mat << dx_mat, dy_mat, dz_mat;
    grad_eig = d_mat * transformed_points_ind.transpose();
    grad[0] = grad_eig[0];
    grad[1] = grad_eig[1];
    grad[2] = grad_eig[2];
    grad[3] = grad_eig[3];
    grad[4] = grad_eig[4];
    grad[5] = grad_eig[5];
  }
  return cost;
}
//  3 if optimization failed, use initial pose
//  2 if optimization diverged, use initial pose
//  1 if successful
// -1 if rejected by initial guess
// -2 if rejected by coverage of area
// -3 initial guess is bad and optimization diverged
// -4 initial guess is bad and optimization failed
int LiDARTag::_optimizePose(ClusterFamily_t& cluster)
{
  int num_points = cluster.merged_data.cols();
  // box_width, tag size
  Eigen::Vector4f template_bound(0.02, cluster.tag_size / 2, cluster.tag_size / 2, 0);  // 4 x 1

  Eigen::MatrixXf transformed_points;  // 4 x n
  double initial_cost =
      evaluateCost(cluster.initial_pose.homogeneous, cluster.merged_data_h, template_bound, transformed_points);
  int status = 1;

  if (_debug_info)
  {
    ROS_DEBUG_STREAM("==== _optimizePose ====");
    float distance = std::sqrt(pow(cluster.average.x, 2) + pow(cluster.average.y, 2) + pow(cluster.average.z, 2));
    ROS_DEBUG_STREAM("Distance : " << distance);
    ROS_DEBUG_STREAM("Actual Points: " << cluster.data.size() + cluster.edge_points.size());
    ROS_DEBUG_STREAM("Inital Cost: " << initial_cost);
    ROS_DEBUG_STREAM("Cost Threshold: " << _optimization_percent * cluster.inliers / 100);
  }

  if (initial_cost > _optimization_percent * cluster.inliers / 100)
  {
    status = -1;
    ROS_DEBUG_STREAM("Status: " << status);

    return status;
  }

  // compute convex hull and its area
  Eigen::MatrixXf convex_hull;
  utils::constructConvexHull(transformed_points, convex_hull);
  float initial_area = utils::computePolygonArea(convex_hull);
  float coverage_area = initial_area / (cluster.tag_size * cluster.tag_size);

  if (coverage_area < _coa_tunable)
  {
    status = -2;
    ROS_DEBUG_STREAM("Status: " << false);

    return status;
  }

  // initial guess
  std::vector<double> x(6);
  x[0] = cluster.initial_pose.translation[0];
  x[1] = cluster.initial_pose.translation[1];
  x[2] = cluster.initial_pose.translation[2];
  if (_derivative_method)
  {
    x[3] = cluster.initial_pose.roll;
    x[4] = cluster.initial_pose.pitch;
    x[5] = cluster.initial_pose.yaw;
  }
  else
  {
    Eigen::Vector3d lie_algebra = utils::Log_SO3(cluster.initial_pose.rotation);
    x[3] = lie_algebra[0];
    x[4] = lie_algebra[1];
    x[5] = lie_algebra[2];
  }

  // x, y, z, r, p, y
  std::vector<double> lb(6), ub(6);
  lb[0] = x[0] - 0.2;  // in meters
  lb[1] = x[1] - 0.2;
  lb[2] = x[2] - 0.2;
  lb[3] = x[3] - _opt_lb;  // 1.57 in rad (180 deg)
  lb[4] = x[4] - _opt_lb;
  lb[5] = x[5] - _opt_lb;

  ub[0] = x[0] + 0.2;
  ub[1] = x[1] + 0.2;
  ub[2] = x[2] + 0.2;
  ub[3] = x[3] + _opt_ub;
  ub[4] = x[4] + _opt_ub;
  ub[5] = x[5] + _opt_ub;

  // Numerical gradient methods
  // good: LN_PRAXIS   LN_NEWUOA_BOUND  LN_SBPLX
  // medium: LN_BOBYQA LN_NELDERMEAD
  // bad: LN_COBYLA
  // nlopt::opt opt(nlopt::LN_SBPLX, 6);

  // Analytical gradient methods
  // LD_SLSQP  LD_MMA
  // LD_TNEWTON_PRECOND_RESTART
  // LD_TNEWTON_PRECOND
  // LD_TNEWTON_RESTART
  // LD_TNEWTON

  nlopt::algorithm opt_method;
  switch (_optimization_solver)
  {
    // Below is numerical gradient-based methods
    case 1:
      opt_method = nlopt::LN_PRAXIS;
      break;
    case 2:
      opt_method = nlopt::LN_NEWUOA_BOUND;
      break;
    case 3:
      opt_method = nlopt::LN_SBPLX;  // recommended
      break;
    case 4:
      opt_method = nlopt::LN_BOBYQA;
      break;
    case 5:
      opt_method = nlopt::LN_NELDERMEAD;
      break;
    case 6:
      opt_method = nlopt::LN_COBYLA;
      break;
    // Below is analytical gradient-based methods
    case 7:
      opt_method = nlopt::LD_SLSQP;  // recommended   200Hz
      break;
    case 8:
      opt_method = nlopt::LD_MMA;  // recommended   120Hz
      break;
    case 9:
      opt_method = nlopt::LD_TNEWTON_PRECOND_RESTART;  // fail 90%
      break;
    case 10:
      opt_method = nlopt::LD_TNEWTON_PRECOND;  // fail 90%
      break;
    case 11:
      opt_method = nlopt::LD_TNEWTON_RESTART;  // fail 80%
      break;
    case 12:
      opt_method = nlopt::LD_TNEWTON;  // fail 90%
      break;
    case 13:
      opt_method = nlopt::LD_LBFGS;  // fail 90%
      break;
    case 14:
      opt_method = nlopt::LD_VAR1;  // fail 90%
      break;
    case 15:
      opt_method = nlopt::LD_VAR2;  // fail 90%
      break;
    default:
      opt_method = nlopt::LD_SLSQP;  // recommended
  }

  // last column saved as other information such as tagsizes
  // Eigen::Vector4f info_for_cost_function;
  // info_for_cost_function << cluster.tag_size, 0, 0, 0;
  Eigen::MatrixXf data(4, num_points + 1);
  data << cluster.merged_data_h, template_bound;

  // XXX: tolerance and timeout
  float x_tol = 1e-5;
  double max_time = 1.0 / 15;  // 15 Hz
  double minf;

  nlopt::opt opt(opt_method, 6);
  opt.set_lower_bounds(lb);
  opt.set_upper_bounds(ub);
  opt.set_xtol_rel(x_tol);
  std::vector<double> steps = { 0.01, 0.01, 0.01, 0.01, 0.01, 0.01 };  // tx, ty, tz, r, p, y
  opt.set_default_initial_step(steps);
  if (_derivative_method)
  {
    opt.set_min_objective(computeCost_euler, &data);
  }
  else
  {
    opt.set_min_objective(computeCost_lie, &data);
  }

  // opt.set_maxtime(max_time);

  // [Error Code]
  // https://github.com/stevengj/nlopt/blob/
  // 2f1fa1cc5f8a3f79834de9d155d5b6ef98aa7e50/src/api/nlopt.h#L162
  // NLOPT_FAILURE = -1,         /* generic failure code */
  // NLOPT_INVALID_ARGS = -2,
  // NLOPT_OUT_OF_MEMORY = -3,
  // NLOPT_ROUNDOFF_LIMITED = -4,
  // NLOPT_FORCED_STOP = -5,
  // NLOPT_SUCCESS = 1,          /* generic success code */
  // NLOPT_STOPVAL_REACHED = 2,
  // NLOPT_FTOL_REACHED = 3,
  // NLOPT_XTOL_REACHED = 4,
  // NLOPT_MAXEVAL_REACHED = 5,
  // NLOPT_MAXTIME_REACHED = 6,

  try
  {
    nlopt::result result = opt.optimize(x, minf);

    if (minf > _optimization_percent * cluster.inliers / 1000)
    {
      status = -3;
      if (_debug_info)
      {
        ROS_WARN_STREAM("Optimized Cost too large: " << std::setprecision(3) << minf);
        ROS_WARN_STREAM("Inital Cost: " << initial_cost);
      }
      if (initial_cost < 0.1 * _optimization_percent * cluster.inliers / 1000)
      {
        status = 2;
        if (_debug_info)
        {
          ROS_WARN_STREAM("Use initial pose.");
        }

        cluster.pose_tag_to_lidar.homogeneous = cluster.initial_pose.homogeneous;
        cluster.pose_tag_to_lidar.translation = cluster.initial_pose.homogeneous.topRightCorner(3, 1);
        cluster.pose_tag_to_lidar.rotation = cluster.initial_pose.homogeneous.topLeftCorner(3, 3);
      }

      return status;
    }

    Eigen::Matrix3f rotation;
    if (_derivative_method)
    {
      Eigen::Quaternion<float> q = Eigen::AngleAxisf(x[3], Eigen::Vector3f::UnitX()) *
                                   Eigen::AngleAxisf(x[4], Eigen::Vector3f::UnitY()) *
                                   Eigen::AngleAxisf(x[5], Eigen::Vector3f::UnitZ());
      rotation = q.matrix();
    }
    else
    {
      Eigen::Vector3d q(x[3], x[4], x[5]);
      rotation = utils::Exp_SO3(q);
    }

    if (_debug_info)
    {
      if (_derivative_method)
      {
        ROS_DEBUG_STREAM("Optimzed euler angle : " << x[3] * 180 / M_PI << ", " << x[4] * 180 / M_PI << ", "
                                                   << x[5] * 180 / M_PI);
      }
      else
      {
        ROS_DEBUG_STREAM("Optimzed lie algebra : " << x[3] << ", " << x[4] << ", " << x[5]);
      }
    }
    Eigen::Vector3f translation(x[0], x[1], x[2]);
    Eigen::Matrix4f homogeneous;
    homogeneous.topLeftCorner(3, 3) = rotation;
    homogeneous.topRightCorner(3, 1) = translation;
    homogeneous.row(3) << 0, 0, 0, 1;

    cluster.pose_tag_to_lidar.homogeneous = homogeneous;
    cluster.pose_tag_to_lidar.translation = homogeneous.topRightCorner(3, 1);
    cluster.pose_tag_to_lidar.rotation = homogeneous.topLeftCorner(3, 3);

    if (_debug_info)
    {
      nlopt_result results_nlopt = static_cast<nlopt_result>(result);
      ROS_DEBUG_STREAM("Optimization result: " << std::string(nlopt_result_to_string(results_nlopt)));
      ROS_DEBUG_STREAM("Optimized cost is: " << std::setprecision(3) << minf);
      ROS_DEBUG_STREAM("Found minimum at \n" << homogeneous);
    }
  }
  catch (std::exception& e)
  {
    status = -4;
    if (_debug_info)
      ROS_WARN_STREAM("Pose optimization failed: " << e.what());
    if (initial_cost < 0.1 * _optimization_percent * cluster.inliers / 1000)
    {
      status = 3;
      cluster.pose_tag_to_lidar.homogeneous = cluster.initial_pose.homogeneous;
      cluster.pose_tag_to_lidar.translation = cluster.initial_pose.homogeneous.topRightCorner(3, 1);
      cluster.pose_tag_to_lidar.rotation = cluster.initial_pose.homogeneous.topLeftCorner(3, 3);
      if (_debug_info)
        ROS_WARN_STREAM("Use initial pose.");
    }
  }

  if (_debug_info)
  {
    ROS_DEBUG_STREAM("Status: " << status);
  }

  return status;
}
}  // namespace BipedLab
