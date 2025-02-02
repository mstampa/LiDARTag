/* Copyright (C) 2013-2020, The Regents of The University of Michigan.
 * All rights reserved.
 * This software was developed in the Biped Lab (https://www.biped.solutions/)
 * under the direction of Jessy Grizzle, grizzle@umich.edu. This software may
 * be available under alternative licensing terms; contact the address above.
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 1. Redistributions of source code must retain the above copyright notice, this
 *    list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * The views and conclusions contained in the software and documentation are those
 * of the authors and should not be interpreted as representing official policies,
 * either expressed or implied, of the Regents of The University of Michigan.
 *
 * AUTHOR: Bruce JK Huang (bjhuang@umich.edu)
 * WEBSITE: https://www.brucerobot.com/
 */

#include "utils.h"

#include <algorithm>  // sort, stable_sort
#include <chrono>
#include <iostream>
#include <math.h>
#include <numeric>  // iota

using namespace std;

namespace BipedLab
{
namespace utils
{
// second
double spendCPUTime(const clock_t& t_end, const clock_t& t_start)
{
  return (((double)(t_end - t_start)) / CLOCKS_PER_SEC);
}

double spendCPUHz(const clock_t& t_end, const clock_t& t_start)
{
  return 1.0 / spendCPUTime(t_end, t_start);
}

double printSpendCPUHz(const clock_t& t_end, const clock_t& t_start, string txt)
{
  auto ans = spendCPUHz(t_end, t_start);
  cout << fixed << setprecision(2) << txt << ans << " [Hz]" << endl;
  return ans;
}

double printSpendCPUHz(const clock_t& t_end, const clock_t& t_start)
{
  string text = "CPU time used: ";
  return printSpendCPUHz(t_end, t_start, text);
}

double spendElapsedTime(const chrono::steady_clock::time_point& t_end, const chrono::steady_clock::time_point& t_start)
{
  chrono::duration<double> duration = t_end - t_start;
  return duration.count();
}

double spendElapsedTimeMilli(const chrono::steady_clock::time_point& t_end,
                             const chrono::steady_clock::time_point& t_start)
{
  return 1e3 * spendElapsedTime(t_end, t_start);
}

double spendElapsedHz(const chrono::steady_clock::time_point& t_end, const chrono::steady_clock::time_point& t_start)
{
  return 1.0 / spendElapsedTime(t_end, t_start);
}

double printSpendElapsedHz(const chrono::steady_clock::time_point& t_end,
                           const chrono::steady_clock::time_point& t_start, string txt)
{
  auto ans = spendElapsedHz(t_end, t_start);
  cout << fixed << setprecision(2) << txt << ans << " [Hz]" << endl;
  return ans;
}

double printSpendElapsedHz(const chrono::steady_clock::time_point& t_end,
                           const chrono::steady_clock::time_point& t_start)
{
  return printSpendElapsedHz(t_end, t_start, "Elapsed time: ");
}

string tranferToLowercase(string& t_data)
{
  transform(t_data.begin(), t_data.end(), t_data.begin(), ::tolower);
  return t_data;
}

void pressEnterToContinue()
{
  int c;
  printf("Press [Enter] key to continue.\n");
  while (getchar() != '\n')
    ;         // option TWO to clean stdin
  getchar();  // wait for ENTER
}

// Checks if a matrix is a valid rotation matrix.
bool isRotationMatrix(Eigen::Matrix3f& t_R)
{
  Eigen::Matrix3f should_be_identity = t_R * t_R.transpose();
  return (should_be_identity - Eigen::Matrix3f::Identity()).norm() < 1e-6;
}

Eigen::Vector3f rotationMatrixToEulerAngles(Eigen::Matrix3f& t_R)
{
  float sy = sqrt(t_R(0, 0) * (0, 0) + t_R(1, 0) * (1, 0));
  bool singular = sy < 1e-6;
  float x, y, z;

  if (!singular)
  {
    x = RAD2DEG(atan(t_R(2, 1) / t_R(2, 2)));
    y = RAD2DEG(atan(-t_R(2, 0) / sy));
    z = RAD2DEG(atan(t_R(1, 0) / t_R(0, 0)));
  }
  else
  {
    x = RAD2DEG(atan(-t_R(1, 2) / t_R(1, 1)));
    y = RAD2DEG(atan(-t_R(2, 0) / sy));
    z = 0;
  }

  return Eigen::Vector3f(x, y, z);
}

/*
 * A function to check if get all parameters
 */
bool checkParameters(int t_n, ...)
{
  va_list vl_num;
  va_start(vl_num, t_n);
  bool Pass = true;

  for (int i = 0; i < t_n; ++i)
  {
    bool Got = va_arg(vl_num, int);
    if (Got)
      continue;
    else
    {
      cout << "didn't get i: " << i << " in the launch file" << endl;
      Pass = false;
      break;
    }
  }
  va_end(vl_num);
  return Pass;
}

// Overload operator <<
// DO NOT want to change their stucture
void COUT(const PointXYZIRT& t_p)
{
  cout << "x: " << t_p.x << ", y: " << t_p.y << ", z: " << t_p.z << ", ring: " << t_p.ring
       << ", intensity: " << t_p.intensity << endl;
}

bool compareIndex(LiDARPoints_t* A, LiDARPoints_t* B)
{
  return A->index < B->index;
}

uint64_t bitShift(string const& t_value)
{
  uint64_t result = 0;

  char const* p = t_value.c_str();
  char const* q = p + t_value.size();
  while (p < q)
  {
    result = (result << 1) + (result << 3) + *(p++) - '0';
  }

  return result;
}

void normalize(vector<float>& x, vector<float>& y, vector<float>& z, vector<float>& I,
               const pcl::PointCloud<LiDARPoints_t*> t_payload)
{
  // normlize the y,z so the top left is (0,0) and bottom right is (1,1)
  // as well as x axis
  //                                           o
  // top left                                 /|
  //        o----_o         LiDAR ---> front o |  back
  //        |     |                          | o
  //        |     |                          |/
  //        o-----o                          o
  //               bottom right
  float front_x = 1e8;
  float back_x = -1e8;
  float bottom_right_y = 1e8;
  float top_left_y = -1e8;
  float bottom_right_z = 1e8;
  float top_left_z = -1e8;

  float max_intensity = -1e8;

  for (int i = 0; i < t_payload.size(); ++i)
  {
    if (t_payload[i]->point.x > back_x)
      back_x = t_payload[i]->point.x;
    if (t_payload[i]->point.x < front_x)
      front_x = t_payload[i]->point.x;

    if (t_payload[i]->point.y > top_left_y)
      top_left_y = t_payload[i]->point.y;
    if (t_payload[i]->point.y < bottom_right_y)
      bottom_right_y = t_payload[i]->point.y;

    if (t_payload[i]->point.z > top_left_z)
      top_left_z = t_payload[i]->point.z;
    if (t_payload[i]->point.z < bottom_right_z)
      bottom_right_z = t_payload[i]->point.z;
    if (t_payload[i]->point.intensity > max_intensity)
      max_intensity = t_payload[i]->point.intensity;
  }

  float dx = abs(front_x - back_x);
  float dy = abs(top_left_y - bottom_right_y);
  float dz = abs(top_left_z - bottom_right_z);
  for (int i = 0; i < t_payload.size(); ++i)
  {
    x[i] = (back_x - t_payload[i]->point.x) / 8;
    y[i] = (top_left_y - t_payload[i]->point.y) / 8;
    z[i] = (top_left_z - t_payload[i]->point.z) / 8;
    I[i] = (t_payload[i]->point.intensity) / 1.5;
  }
}

void normalizeByAve(vector<float>& x, vector<float>& y, vector<float>& z, vector<float>& I,
                    const pcl::PointCloud<LiDARPoints_t*> t_payload)
{
  float ave_x = 0;
  float ave_y = 0;
  float ave_z = 0;

  for (int i = 0; i < t_payload.size(); ++i)
  {
    ave_x += t_payload[i]->point.x;
    ave_y += t_payload[i]->point.y;
    ave_z += t_payload[i]->point.z;
    x[i] = t_payload[i]->point.x;
    y[i] = t_payload[i]->point.y;
    z[i] = t_payload[i]->point.z;
    I[i] = t_payload[i]->point.intensity;
  }
  ave_x /= t_payload.size();
  ave_y /= t_payload.size();
  ave_z /= t_payload.size();

  for (int i = 0; i < t_payload.size(); ++i)
  {
    x[i] = (x[i] - ave_x) / 5;
    y[i] = (y[i] - ave_y) / 5;
    z[i] = (z[i] - ave_z) / 5;
    I[i] /= 1.5;
  }
}

PointXYZIRT pointsAddDivide(const PointXYZIRT& t_p1, const PointXYZIRT& t_p2, float t_d)
{
  assert(t_d != 0);
  PointXYZIRT tmp;
  tmp.x = (t_p1.x + t_p2.x) / t_d;
  tmp.y = (t_p1.y + t_p2.y) / t_d;
  tmp.z = (t_p1.z + t_p2.z) / t_d;
  tmp.intensity = (t_p1.intensity + t_p2.intensity) / t_d;
  return tmp;
}

// form vector from p1 to p2. ie p2-p1
PointXYZIRT vectorize(const PointXYZIRT& t_p1, const PointXYZIRT& t_p2)
{
  PointXYZIRT tmp;
  tmp.x = (t_p2.x - t_p1.x);
  tmp.y = (t_p2.y - t_p1.y);
  tmp.z = (t_p2.z - t_p1.z);
  tmp.intensity = (t_p2.intensity - t_p1.intensity);
  return tmp;
}

float dot(const PointXYZIRT& t_p1, const PointXYZIRT& t_p2)
{
  return t_p1.y * t_p2.y + t_p1.z * t_p2.z;
}

float Norm(const PointXYZIRT& t_p)
{
  return hypot(t_p.y, t_p.z);
}

double MVN(const float& t_tag_size, const int& t_d, const Eigen::Vector2f& t_X, const Eigen::Vector2f t_mean)
{
  Eigen::Matrix2f Sigma;
  Sigma << t_tag_size / t_d / 2, 0, 0, t_tag_size / t_d / 2;
  double sqrt2pi = sqrt(2 * M_PI);
  double QuadForm = (t_X - t_mean).transpose() * Sigma.inverse() * (t_X - t_mean);
  double norm = pow(sqrt2pi, -2) * pow(Sigma.determinant(), -0.5);
  return norm * exp(-0.5 * QuadForm);
}

// step between p1 and p2
float getStep(const PointXYZIRT& t_p1, const PointXYZIRT& t_p2, const int t_d)
{
  return hypot(t_p2.y - t_p1.y, t_p2.z - t_p1.z) / t_d;
}

// To get the t where p1 + t * v12 is the point that p projects onto line p12
void getProjection(const PointXYZIRT& t_p1, const PointXYZIRT& t_p2, const PointXYZIRT& t_p, float& k,
                   Eigen::Vector2f& t_v)
{
  // form vector from p1 to p2 and p1 to p
  PointXYZIRT v12 = vectorize(t_p1, t_p2);
  PointXYZIRT v1p = vectorize(t_p1, t_p);

  k = abs(dot(v12, v1p) / Norm(v12));
}

void assignCellIndex(const float& t_tag_size, const Eigen::Matrix3f& t_R, PointXYZIRT& t_p_reference,
                     const PointXYZIRT& t_average, const int t_d, PayloadVoting_t& t_vote)
{
  // R: Payload p -> reference x
  // prepare for Gaussian
  float xOffset = t_vote.p->x - t_average.x;
  float yOffset = t_vote.p->y - t_average.y;
  float zOffset = t_vote.p->z - t_average.z;

  float x = xOffset * t_R(0, 0) + yOffset * t_R(0, 1) + zOffset * t_R(0, 2);
  float y = xOffset * t_R(1, 0) + yOffset * t_R(1, 1) + zOffset * t_R(1, 2);
  float z = xOffset * t_R(2, 0) + yOffset * t_R(2, 1) + zOffset * t_R(2, 2);

  // y,z should range int_ between -3s and 3s
  t_p_reference.x = x;
  t_p_reference.y = y;
  t_p_reference.z = z;
  float ss = t_tag_size / t_d;                                    // scale back to the unit square
  y = max(min(y, t_d / 2 * ss), (-t_d / 2 * ss + (float)0.001));  // don't match to 6
  z = max(min(z, t_d / 2 * ss), (-t_d / 2 * ss + (float)0.001));  // don't match to 6
  int cellIndexT = t_d / 2 + floor(-y / ss);
  int cellIndexK = t_d / 2 + floor(-z / ss);

  float cy = (ceil(y / ss) - 0.5) * ss;  // offset to center of each ceil
  float cz = (ceil(z / ss) - 0.5) * ss;

  // which grid it belongs to (in 1-16 vector form)?
  Eigen::Vector2f X(y, z);
  Eigen::Vector2f Mean(cy, cz);
  t_vote.centroid.x = 0;
  t_vote.centroid.y = cy;
  t_vote.centroid.z = cz;
  t_vote.cell = t_d * cellIndexK + cellIndexT;
  t_vote.weight = MVN(t_tag_size, t_d, X, Mean);
}

// normalize weight and classify them into grid
void sortPointsToGrid(vector<vector<PayloadVoting_t*>>& t_grid, vector<PayloadVoting_t>& t_votes)
{
  for (int i = 0; i < t_votes.size(); ++i)
    t_grid[t_votes[i].cell].push_back(&t_votes[i]);
}

void formGrid(Eigen::MatrixXf& t_vertices, float x, float y, float z, float t_tag_size)
{
  // define 5 points in reference coord frame: x0,...x4 (x0==(0,0,0))
  Eigen::Vector3f tmp;
  tmp << x, y, z;  // 0,0,0

  // center
  t_vertices.col(0) = tmp;  // center of ref model

  // p1
  tmp[1] = y + t_tag_size / 2;
  tmp[2] = z + t_tag_size / 2;
  t_vertices.col(1) = tmp;

  // p2
  tmp[1] = y + t_tag_size / 2;
  tmp[2] = z - t_tag_size / 2;
  t_vertices.col(2) = tmp;

  // p3
  tmp[1] = y - t_tag_size / 2;
  tmp[2] = z - t_tag_size / 2;
  t_vertices.col(3) = tmp;

  // p4
  tmp[1] = y - t_tag_size / 2;
  tmp[2] = z + t_tag_size / 2;
  t_vertices.col(4) = tmp;
}

void fitGrid(Eigen::MatrixXf& GridVertices, Eigen::Matrix3f& H, const PointXYZIRT& t_p1, const PointXYZIRT& t_p2,
             const PointXYZIRT& t_p3, const PointXYZIRT& t_p4)
{
  Eigen::MatrixXf payload_vertices(3, 4);
  payload_vertices(0, 0) = t_p1.x;
  payload_vertices(1, 0) = t_p1.y;
  payload_vertices(2, 0) = t_p1.z;

  payload_vertices(0, 1) = t_p2.x;
  payload_vertices(1, 1) = t_p2.y;
  payload_vertices(2, 1) = t_p2.z;

  payload_vertices(0, 2) = t_p3.x;
  payload_vertices(1, 2) = t_p3.y;
  payload_vertices(2, 2) = t_p3.z;

  payload_vertices(0, 3) = t_p4.x;
  payload_vertices(1, 3) = t_p4.y;
  payload_vertices(2, 3) = t_p4.z;

  Eigen::Matrix3f M = GridVertices.rightCols(4) * payload_vertices.transpose();
  Eigen::JacobiSVD<Eigen::MatrixXf> svd(M, Eigen::ComputeFullU | Eigen::ComputeFullV);

  Eigen::Matrix<float, 3, 3, Eigen::DontAlign> R = svd.matrixU() * svd.matrixV().transpose();
  H = R;  // H: payload -> ref
}

void fitGrid_new(Eigen::MatrixXf& GridVertices, Eigen::Matrix3f& H, Eigen::MatrixXf& payload_vertices)
{
  Eigen::Matrix3f M = GridVertices.rightCols(4) * payload_vertices.transpose();
  Eigen::JacobiSVD<Eigen::MatrixXf> svd(M, Eigen::ComputeFullU | Eigen::ComputeFullV);
  Eigen::Matrix<float, 3, 3, Eigen::DontAlign> R = svd.matrixU() * svd.matrixV().transpose();
  H = R;  // H: payload -> ref
}
PointXYZIRT toVelodyne(const Eigen::Vector3f& t_p)
{
  PointXYZIRT point;
  point.x = t_p[0];
  point.y = t_p[1];
  point.z = t_p[2];
  return point;
}

Eigen::Vector3f toEigen(const PointXYZIRT& t_point)
{
  Eigen::Vector3f tmp;
  tmp[0] = t_point.x;
  tmp[1] = t_point.y;
  tmp[2] = t_point.z;
  return tmp;
}

void minus(PointXYZIRT& t_p1, const PointXYZIRT& t_p2)
{
  t_p1.x = t_p1.x - t_p2.x;
  t_p1.y = t_p1.y - t_p2.y;
  t_p1.z = t_p1.z - t_p2.z;
}

float distance(const PointXYZIRT& t_p1, const PointXYZIRT& t_p2)
{
  return hypot(t_p1.x - t_p2.x, t_p1.y - t_p2.y, t_p1.z - t_p2.z);
}

/*
 * A function to calculate angle between va and vb
 * return: angle in degree
 */
template <class T, class U>
float getAngle(T a, U b)
{
  return RAD2DEG(acos(dot(a, b) / (Norm(a) * Norm(b))));
}

/*
 * Check if 4 four corners are valid
 * return  0: valid corners
 * return -1: incorrect distance
 * return -2: incorrect angle
 */
int checkCorners(const float Tagsize, const PointXYZIRT& t_p1, const PointXYZIRT& t_p2, const PointXYZIRT& t_p3,
                 const PointXYZIRT& t_p4)
{
  // XXX tunable
  float AngleLowerBound = 75;
  float AngleUpperBound = 105;

  auto dists = vector({ distance(t_p1, t_p2), distance(t_p1, t_p3), distance(t_p1, t_p4), distance(t_p2, t_p3),
                        distance(t_p2, t_p4), distance(t_p3, t_p4) });
  auto cmp = Tagsize / 3.0;
  if (any_of(dists.begin(), dists.end(), [cmp](float d) { return d < cmp; }))
    return -1;

  // angle between p12 and p14
  float Angle1 = getAngle<PointXYZIRT, PointXYZIRT>(vectorize(t_p1, t_p2), vectorize(t_p1, t_p4));
  if ((Angle1 < AngleLowerBound) || (AngleUpperBound < Angle1))
    return -2;

  // angle between p21 and p23
  float Angle2 = getAngle<PointXYZIRT, PointXYZIRT>(vectorize(t_p2, t_p1), vectorize(t_p2, t_p3));
  if ((Angle2 < AngleLowerBound) || (AngleUpperBound < Angle2))
    return -2;

  // angle between p32 and p34
  float Angle3 = getAngle<PointXYZIRT, PointXYZIRT>(vectorize(t_p3, t_p2), vectorize(t_p3, t_p4));
  if ((Angle3 < AngleLowerBound) || (AngleUpperBound < Angle3))
    return -2;

  // angle between p43 and p41
  float Angle4 = getAngle<PointXYZIRT, PointXYZIRT>(vectorize(t_p4, t_p3), vectorize(t_p4, t_p1));
  if ((Angle4 < AngleLowerBound) || (AngleUpperBound < Angle4))
    return -2;

  return 0;
}

/* Function for creating blockdiagonal given arbitrary number of arguments.  */
template <class T>
T blockMatrix(int t_n, ...)
{
  va_list vl_num;
  va_start(vl_num, t_n);
  int cols_now = 0;
  int rows_now = 0;

  for (int i = 0; i < t_n; ++i)
  {
    T matrix = va_arg(vl_num, T);
    cols_now = cols_now + matrix.cols();
    rows_now = rows_now + matrix.rows();
  }
  va_end(vl_num);
  T Mblock = T::Zero(rows_now, cols_now);
  va_list vl;
  va_start(vl, t_n);
  int rows = 0;
  int cols = 0;
  for (int i = 0; i < t_n; ++i)
  {
    T matrix = va_arg(vl, T);
    Mblock.block(rows, cols, matrix.rows(), matrix.cols()) = matrix;
    rows += matrix.rows();
    cols += matrix.cols();
  }
  return Mblock;
}

// pose is geometry_msgs pose
// template <class T>
// Eigen::Matrix4d poseToEigenMatrix(const T &pose){
Eigen::Matrix4d poseToEigenMatrix(const geometry_msgs::Pose& t_pose)
{
  Eigen::Matrix4d matrix_pose = Eigen::Matrix4d::Identity();
  matrix_pose(0, 3) = t_pose.position.x;
  matrix_pose(1, 3) = t_pose.position.y;
  matrix_pose(2, 3) = t_pose.position.z;
  matrix_pose.topLeftCorner(3, 3) << qToR(t_pose);

  return matrix_pose;
}

// pose is geometry_msgs pose
template <class T>
Eigen::Matrix3d qToR(const T& t_pose)
{
  Eigen::Matrix3d R = Eigen::Matrix3d::Identity();
  double a = t_pose.orientation.w;
  double b = t_pose.orientation.x;
  double c = t_pose.orientation.y;
  double d = t_pose.orientation.z;
  R(0, 0) = a * a + b * b - c * c - d * d;
  R(0, 1) = 2 * (b * c - a * d);
  R(0, 2) = 2 * (b * d + a * c);

  R(1, 0) = 2 * (b * c + a * d);
  R(1, 1) = a * a - b * b + c * c - d * d;
  R(1, 2) = 2 * (c * d - a * b);

  R(2, 0) = 2 * (b * d - a * c);
  R(2, 1) = 2 * (c * d + a * b);
  R(2, 2) = a * a - b * b - c * c + d * d;

  return R;
}

Eigen::Matrix3d qToR(const Eigen::Vector3f& t_pose)
{
  Eigen::Matrix3d R = Eigen::Matrix3d::Identity();
  double a = 0;          // w
  double b = t_pose(0);  // x
  double c = t_pose(1);  // y
  double d = t_pose(2);  // z
  R(0, 0) = a * a + b * b - c * c - d * d;
  R(0, 1) = 2 * (b * c - a * d);
  R(0, 2) = 2 * (b * d + a * c);

  R(1, 0) = 2 * (b * c + a * d);
  R(1, 1) = a * a - b * b + c * c - d * d;
  R(1, 2) = 2 * (c * d - a * b);

  R(2, 0) = 2 * (b * d - a * c);
  R(2, 1) = 2 * (c * d + a * b);
  R(2, 2) = a * a - b * b - c * c + d * d;

  return R;
}

/*
 * Returns all numbers not in set, where the total set is [0,n)
 */
vector<int> complementOfSet(const vector<int>& set, size_t n)
{
  size_t curr_idx = 0;
  vector<int> complement;

  for (auto i = 0; i < n; i++)
  {
    if (curr_idx >= set.size())
    {
      complement.push_back(i);
    }
    else if (i < set[curr_idx])
    {
      complement.push_back(i);
    }
    else if (i == set[curr_idx])
    {
      curr_idx++;
    }  // Inlier
  }

  return complement;
}

Eigen::Vector3f cross_product(Eigen::Vector3f v1, Eigen::Vector3f v2)
{
  Eigen::Vector3f res;
  res[0] = v1[1] * v2[2] - v1[2] * v2[1];
  res[1] = v1[2] * v2[0] - v1[0] * v2[2];
  res[2] = v1[0] * v2[1] - v1[1] * v2[0];
  return res;
}

void readDirectory(const string& name, vector<string>& v)
{
  boost::filesystem::path p(name);
  boost::filesystem::directory_iterator start(p);
  boost::filesystem::directory_iterator end;
  transform(start, end, back_inserter(v), PathLeafString());
  sort(v.begin(), v.end());
}

float computeMedian(const Eigen::VectorXf& eigen_vec)
{
  assert(eigen_vec.size() != 0);
  vector<float> vec(eigen_vec.data(), eigen_vec.data() + eigen_vec.size());
  assert(vec.size() != 0);
  if (vec.size() % 2 == 0)
  {
    const auto median_it1 = vec.begin() + vec.size() / 2 - 1;
    const auto median_it2 = vec.begin() + vec.size() / 2;

    nth_element(vec.begin(), median_it1, vec.end());
    const auto e1 = *median_it1;
    nth_element(vec.begin(), median_it2, vec.end());
    const auto e2 = *median_it2;

    return (e1 + e2) / 2;
  }
  else
  {
    const auto median_it = vec.begin() + vec.size() / 2;
    nth_element(vec.begin(), median_it, vec.end());

    return *median_it;
  }
}

Eigen::MatrixXf convertXYZIToHomogeneous(const Eigen::MatrixXf& mat_xyzi)
{
  assert(mat_xyzi.rows() == 4 && "The input dimension is wrong, it should be four");
  Eigen::MatrixXf mat_h = mat_xyzi;
  mat_h.row(3).setOnes();
  return mat_h;
}

/* A function takes in rot_v as r, p, y in deg and trans_v as x, y, z in meter.
 * It returns a rigid-body transformation.
 * [Note] The rotation matrix from rpy follows the "XYZ" convention
 */
Eigen::Matrix4f computeTransformation(Eigen::Vector3f rot_v, Eigen::Vector3f trans_v)
{
  Eigen::Matrix4f H = Eigen::Matrix4f::Identity(4, 4);
  H.topLeftCorner(3, 3) = computeRotX(rot_v(0)) * computeRotY(rot_v(1)) * computeRotZ(rot_v(2));
  H.topRightCorner(3, 1) = trans_v;
  return H;
}

double get_sign(double x)
{
  return (x >= 0) ? 1 : -1;
}

Eigen::Matrix3f skew(const Eigen::Vector3d v)
{
  Eigen::Matrix3f m;
  m << 0, -v[2], v[1], v[2], 0, -v[0], -v[1], v[0], 0;
  return m;
}

Eigen::Vector3d unskew(const Eigen::Matrix3f Ax)
{
  Eigen::Vector3d v(Ax(2, 1), Ax(0, 2), Ax(1, 0));
  return v;
}

Eigen::Matrix3f Exp_SO3(const Eigen::Vector3d w)
{
  double theta = w.norm();
  Eigen::Matrix3f A = skew(w);
  Eigen::Matrix3f output = Eigen::Matrix3f::Identity();

  if (theta == 0)
    output += +(sin(theta) / theta) * A + ((1 - cos(theta)) / (theta * theta)) * A * A;

  return output;
}

Eigen::Vector3d Log_SO3(const Eigen::Matrix3f A)
{
  double theta = acos((A(0, 0) + A(1, 1) + A(2, 2) - 1) / 2);
  Eigen::Matrix3f A_transpose = A.transpose();
  Eigen::Vector3d output;
  if (theta == 0)
  {
    output = Eigen::Vector3d::Zero(3);
  }
  else
  {
    output = unskew(theta * (A - A_transpose) / (2 * sin(theta)));
  }
  return output;
}

// 3D cross product of OA and OB vectors,
// (i.e x-component of their "2D" cross product,
// but remember that it is not defined in "2D").
// Returns a positive value, if OAB makes a counter-clockwise turn,
// negative for clockwise turn, and zero if the points are collinear.
// compute on y-z plane and direction on x or -x axis
double cross(const Eigen::Vector4f& O, const Eigen::Vector4f& A, const Eigen::Vector4f& B)
{
  return (A[1] - O[1]) * (B[2] - O[2]) - (A[2] - O[2]) * (B[1] - O[1]);
}

// comparator of transformed points, used for convex hull
void sortEigenVectorIndices(const Eigen::MatrixXf& mat, Eigen::VectorXi& indices)
{
  // initialize original index locations
  int num_elements = mat.cols();
  vector<int> idx(num_elements);
  iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing values in v
  // using stable_sort instead of sort
  // to avoid unnecessary index re-orderings
  // when v contains elements of equal values
  stable_sort(idx.begin(), idx.end(), [&mat](int i1, int i2) {
    return (mat(1, i1) < mat(1, i2)) || (mat(1, i1) == mat(1, i2) && mat(2, i1) < mat(2, i2));
  });

  int* ptr = &idx[0];
  Eigen::Map<Eigen::VectorXi> tmp(ptr, num_elements);
  indices = tmp;
}

// Vertices are in a 4xn matrix as [x,y,z,1]
// This function computes the area on y-z plane
float computePolygonArea(const Eigen::MatrixXf& vertices)
{
  // Initialze area
  float area = 0.0;

  int num_pts = vertices.cols();

  // Calculate value of shoelace formula
  int j = num_pts - 1;
  for (int i = 0; i < num_pts; i++)
  {
    area += (vertices(1, j) + vertices(1, i)) * (vertices(2, j) - vertices(2, i));
    j = i;  // j is previous vertex to i
  }

  // Return absolute value
  return abs(area / 2.0);
}

// Returns a list of points on the convex hull in counter-clockwise order.
void constructConvexHull(const Eigen::MatrixXf& P, Eigen::MatrixXf& convex_hull)
{
  size_t n = P.cols();
  if (n <= 3)
    convex_hull = P;

  // Sort points lexicographically
  Eigen::VectorXi indices;
  sortEigenVectorIndices(P, indices);
  Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> perm(indices);
  Eigen::MatrixXf sorted_points = P * perm;
  Eigen::Vector4f left_most = sorted_points.col(0);
  Eigen::Vector4f right_most = sorted_points.col(n - 1);
  Eigen::MatrixXf up = Eigen::MatrixXf::Zero(4, n);
  Eigen::MatrixXf down = Eigen::MatrixXf::Zero(4, n);
  up.col(0) = left_most;
  int k = 0;
  int j = 0;

  for (int i = 0; i < sorted_points.cols(); i++)
  {
    Eigen::Vector4f cur_pt = sorted_points.col(i);
    if (i == sorted_points.cols() - 1 || cross(left_most, cur_pt, right_most) > 0)
    {
      while (k >= 2 && cross(up.col(k - 2), up.col(k - 1), cur_pt) <= 0)
        k--;
      up.col(k++) = cur_pt;
    }

    if (i == sorted_points.cols() - 1 || cross(left_most, cur_pt, right_most) <= 0)
    {
      while (j >= 2 && cross(down.col(j - 2), down.col(j - 1), cur_pt) >= 0)
        j--;
      down.col(j++) = cur_pt;
    }
  }
  convex_hull.resize(up.rows(), k + j - 1);
  convex_hull << up.leftCols(k), down.leftCols(j - 1).rowwise().reverse();
}

}  // namespace utils
}  // namespace BipedLab
