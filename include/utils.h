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

#pragma once
#ifndef UTILS_H
#define UTILS_H

#include "lidartag.h"

#include <boost/filesystem.hpp>
#include <velodyne_pcl/point_types.h>
#include <velodyne_pointcloud/pointcloudXYZIRT.h>

#include <Eigen/Core>

#include <algorithm>  // to use to_lower function
#include <math.h>
#include <stdarg.h>  // for variadic functions
#include <string>
#include <vector>

namespace BipedLab
{
namespace utils
{
using PointXYZIRT = BipedLab::PointXYZIRT;

double spendCPUTime(const std::clock_t& t_end, const std::clock_t& t_start);
double spendCPUHz(const std::clock_t& t_end, const std::clock_t& t_start);
double printSpendCPUHz(const std::clock_t& t_end, const std::clock_t& t_start);
double printSpendCPUHz(const std::clock_t& t_end, const std::clock_t& t_start, std::string txt);

double spendElapsedTimeMilli(const std::chrono::steady_clock::time_point& t_end,
                             const std::chrono::steady_clock::time_point& t_start);

double spendElapsedTime(const std::chrono::steady_clock::time_point& t_end,
                        const std::chrono::steady_clock::time_point& t_start);

double spendElapsedHz(const std::chrono::steady_clock::time_point& t_end,
                      const std::chrono::steady_clock::time_point& t_start);

double printSpendElapsedHz(const std::chrono::steady_clock::time_point& t_end,
                           const std::chrono::steady_clock::time_point& t_start, std::string txt);

double printSpendElapsedHz(const std::chrono::steady_clock::time_point& t_end,
                           const std::chrono::steady_clock::time_point& t_start);

/*
 * Generic function to find an element in vector and also its position.
 * It returns a pair of bool & int i.e.
 * bool : Represents if element is present in vector or not.
 * int : Represents the index of element in vector if its found else -1
 */
template <typename T>
std::pair<bool, int> findInVector(const std::vector<T>& vecOfElements, const T& element)
{
  std::pair<bool, int> result;
  // Find given element in vector
  auto it = std::find(vecOfElements.begin(), vecOfElements.end(), element);
  if (it != vecOfElements.end())
  {
    result.second = distance(vecOfElements.begin(), it);
    result.first = true;
  }
  else
  {
    result.first = false;
    result.second = -1;
  }
  return result;
}

std::string tranferToLowercase(std::string& t_data);

void pressEnterToContinue();

bool isRotationMatrix(Eigen::Matrix3f& t_R);
Eigen::Vector3f rotationMatrixToEulerAngles(Eigen::Matrix3f& t_R);

bool checkParameters(int t_n, ...);
void COUT(const PointXYZIRT& t_p);
bool compareIndex(LiDARPoints_t* A, LiDARPoints_t* B);
uint64_t bitShift(std::string const& t_value);

void normalizeByAve(std::vector<float>& t_x, std::vector<float>& t_y, std::vector<float>& t_z, std::vector<float>& t_I,
                    const pcl::PointCloud<LiDARPoints_t*> payload);
void normalize(std::vector<float>& t_x, std::vector<float>& t_y, std::vector<float>& t_z, std::vector<float>& t_I,
               const pcl::PointCloud<LiDARPoints_t*> t_payload);

PointXYZIRT pointsAddDivide(const PointXYZIRT& t_p1, const PointXYZIRT& t_p2, float t_d = 1);

PointXYZIRT vectorize(const PointXYZIRT& t_p1, const PointXYZIRT& t_p2);

float dot(const PointXYZIRT& t_p1, const PointXYZIRT& t_p2);

float Norm(const PointXYZIRT& t_p);

// a function to determine the step of given two points
float getStep(const PointXYZIRT& t_p1, const PointXYZIRT& t_p2, const int t_d);

void getProjection(const PointXYZIRT& t_p1, const PointXYZIRT& t_p2, const PointXYZIRT& t_p, float& t_k,
                   Eigen::Vector2f& t_v);

double MVN(const float& t_tag_size, const int& t_d, const Eigen::Vector2f& t_X, const Eigen::Vector2f t_mean);

void assignCellIndex(const float& t_tag_size, const Eigen::Matrix3f& t_R, PointXYZIRT& t_p_reference,
                     const PointXYZIRT& t_average, const int t_d, PayloadVoting_t& t_vote);

void sortPointsToGrid(std::vector<std::vector<PayloadVoting_t*>>& t_grid, std::vector<PayloadVoting_t>& t_votes);

Eigen::Vector2f pointToLine(const PointXYZIRT& t_p1, const PointXYZIRT& t_p2, const PointXYZIRT& t_p);

void formGrid(Eigen::MatrixXf& t_vertices, float t_x, float t_y, float t_z, float t_tag_size);

void fitGrid(Eigen::MatrixXf& t_vertices, Eigen::Matrix3f& t_R, const PointXYZIRT& t_p1, const PointXYZIRT& t_p2,
             const PointXYZIRT& t_p3, const PointXYZIRT& t_p4);

void fitGrid_new(Eigen::MatrixXf& t_vertices, Eigen::Matrix3f& H, Eigen::MatrixXf& t_payload_vertices);

float distance(const PointXYZIRT& t_p1, const PointXYZIRT& t_p2);
template <class T, class U>
float getAngle(T a, U b);
double get_sign(double x);
Eigen::Matrix3f skew(const Eigen::Vector3d t_v);

Eigen::Matrix3f Exp_SO3(const Eigen::Vector3d t_w);

Eigen::Vector3d unskew(const Eigen::Matrix3f t_Ax);

Eigen::Vector3d Log_SO3(const Eigen::Matrix3f t_A);

int checkCorners(const float t_tag_size, const PointXYZIRT& t_p1, const PointXYZIRT& t_p2, const PointXYZIRT& t_p3,
                 const PointXYZIRT& t_p4);

PointXYZIRT toVelodyne(const Eigen::Vector3f& t_p);
Eigen::Vector3f toEigen(const PointXYZIRT& t_point);
void minus(PointXYZIRT& t_p1, const PointXYZIRT& t_p2);

template <class T>
T blockMatrix(int t_n, ...);

// template <class T>
Eigen::Matrix4d poseToEigenMatrix(const geometry_msgs::Pose& t_pose);

template <class T>
Eigen::Matrix3d qToR(const T& t_pose);

Eigen::Matrix3d qToR(const Eigen::Vector3f& t_pose);

template <class T>
bool goodNumber(T t_number)
{
  return !std::isinf(t_number) && !std::isnan(t_number);
}

/*
 * Removes all undesired elements in vec listed by index in rm
 */
template <typename T, template <class> class C>
void removeIndicesFromVector(C<T>& c, std::vector<int>& rm)
{
  using std::begin;
  using std::end;

  // If indices in rm are not sorted, sort
  if (!std::is_sorted(rm.begin(), rm.end()))
  {
    std::sort(rm.begin(), rm.end());
  }

  auto rm_iter = rm.begin();
  std::size_t curr_idx = 0;

  const auto pred = [&](const T& /*val*/) {
    // Any more to remove?
    if (rm_iter == rm.end())
    {
      return false;
    }

    // Does current index match remove index
    if (*rm_iter == curr_idx++)
    {
      ++rm_iter;
      return true;
    }
    return false;
  };

  c.erase(std::remove_if(begin(c), end(c), pred), end(c));
}

std::vector<int> complementOfSet(const std::vector<int>& set, std::size_t n);

Eigen::Vector3f cross_product(Eigen::Vector3f v1, Eigen::Vector3f v2);

/* A function to read data in a csv file
 * and load the data to an eigen matrix
 */
template <typename M>
M loadCSV(const std::string& path)
{
  std::ifstream indata;
  indata.open(path);
  std::string line;
  std::vector<typename M::Scalar> values;
  uint rows = 0;
  while (std::getline(indata, line))
  {
    std::stringstream lineStream(line);
    std::string cell;
    while (std::getline(lineStream, cell, ','))
    {
      values.push_back(std::stod(cell));
    }
    ++rows;
  }

  return Eigen::Map<
      const Eigen::Matrix<typename M::Scalar, M::RowsAtCompileTime, M::ColsAtCompileTime, Eigen::RowMajor>>(
      values.data(), rows, values.size() / rows);
}

/* A function to read files from a folder
 *
 */
void readDirectory(const std::string& name, std::vector<std::string>& v);

/* A function to compute the median of a eigen vector
 */
float computeMedian(const Eigen::VectorXf& eigen_vec);

template <typename T>
void printVector(const std::vector<T>& vec)
{
  std::copy(vec.begin(), vec.end(), std::ostream_iterator<T>(std::cout, "\n"));
}

template <typename T>
std::vector<typename T::Scalar> convertEigenToSTDVector(const T& mat)
{
  return std::vector<typename T::Scalar>(mat.data(), mat.data() + mat.rows() * mat.cols());
}

Eigen::MatrixXf convertXYZIToHomogeneous(const Eigen::MatrixXf& mat_xyzi);

template <typename T>
Eigen::Matrix3f computeRotX(T deg)
{
  T rad = DEG2RAD(deg);
  Eigen::Matrix3f rotx;
  rotx << 1, 0, 0, 0, std::cos(rad), -std::sin(rad), 0, std::sin(rad), std::cos(rad);
  return rotx;
}

template <typename T>
Eigen::Matrix3f computeRotY(T deg)
{
  T rad = DEG2RAD(deg);
  Eigen::Matrix3f roty;
  roty << std::cos(rad), 0, std::sin(rad), 0, 1, 0, -std::sin(rad), 0, std::cos(rad);
  return roty;
}

template <typename T>
Eigen::Matrix3f computeRotZ(T deg)
{
  T rad = DEG2RAD(deg);
  Eigen::Matrix3f rotz;
  rotz << std::cos(rad), -std::sin(rad), 0, std::sin(rad), std::cos(rad), 0, 0, 0, 1;
  return rotz;
}

double cross(const Eigen::Vector4f& O, const Eigen::Vector4f& A, const Eigen::Vector4f& B);

Eigen::Matrix4f computeTransformation(Eigen::Vector3f rot_v, Eigen::Vector3f trans_v);

void sortEigenVectorIndices(const Eigen::VectorXf& v, Eigen::VectorXi& indices);
void constructConvexHull(const Eigen::MatrixXf& P, Eigen::MatrixXf& convex_hull);

float computePolygonArea(const Eigen::MatrixXf& vertices);

}  // namespace utils

}  // namespace BipedLab

#endif
