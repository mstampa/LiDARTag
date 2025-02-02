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

#include "nanoflann.hpp"

#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>
#include <pcl_conversions/pcl_conversions.h>
#include <pcl/point_types.h>
#include <pcl/PCLPointCloud2.h>
#include <pcl/features/normal_3d.h>
#include <pcl/ModelCoefficients.h>
#include <tf/LinearMath/Transform.h>
#include <velodyne_pcl/point_types.h>

#include <Eigen/Core>
#include <Eigen/Sparse>

#include <chrono>  // high_resolution_clock

namespace BipedLab
{
using PointXYZIRT = velodyne_pcl::PointXYZIRT;

typedef struct QuickDecodeEntry
{
  uint64_t rcode;     // the queried code
  uint16_t id;        // the tag id (a small integer)
  uint16_t hamming;   // how many errors corrected?
  uint16_t rotation;  // number of rotations [0, 3]
} QuickDecodeEntry_t;

typedef struct QuickDecode
{
  int nentries;
  QuickDecodeEntry_t* entries;
} QuickDecode_t;

typedef struct PayloadVoting
{
  PointXYZIRT* p;
  float weight;
  int cell;
  PointXYZIRT centroid;
} PayloadVoting_t;

typedef struct MaxMin
{
  int min;
  int average;
  int max;
} MaxMin_t;

struct angleComparision
{
  bool operator()(const float& i, const float& j) const
  {
    float threshold = 0.3;
    return (std::abs(i - j) < threshold ? false : i < j);
  }
};

// Structure for LiDAR system
typedef struct LiDARSystem
{
  std::vector<std::vector<int>> point_count_table;  // point per ring  PointCountTable[Scan][ring]
  std::vector<MaxMin_t> max_min_table;              // max min points in a scan
  std::vector<MaxMin_t> ring_average_table;      // max, min, average points in a ring, examed through out a few seconds
  std::set<float, angleComparision> angle_list;  // store the angle of each point

  double points_per_square_meter_at_one_meter;  // TODO: only assume place the tag at dense-point area
  double beam_per_vertical_radian;
  double point_per_horizontal_radian;
} LiDARSystem_t;

// Struture for LiDAR PointCloud with index
typedef struct LiDARPoints
{
  PointXYZIRT point;
  int index;
  int valid;
  double tag_size;   // only take abs value due to uncertain direction
  double box_width;  // Also account for direction by knowing tag is white to black
  double threshold_intensity;
} LiDARPoints_t;

typedef struct TagLines
{
  int upper_ring;
  int lower_ring;
  std::vector<LiDARPoints_t*> upper_line;  // basically just a specific ring, just point to it should be fine
  std::vector<LiDARPoints_t*> lower_line;  // same above
  std::vector<LiDARPoints_t*> left_line;   // same
  std::vector<LiDARPoints_t*> right_line;  // same above

  std::vector<LiDARPoints_t*> bottom_left;   // basically just a specific ring, just point to it should be fine
  std::vector<LiDARPoints_t*> bottom_right;  // same above
  std::vector<LiDARPoints_t*> top_left;      // same
  std::vector<LiDARPoints_t*> top_right;     // same above
} TagLines_t;

typedef struct TagBoundaries
{
  int status;                              // 0 is up right, 1 is tilted
  std::vector<LiDARPoints_t*> line_one;    // basically just a specific ring, just point to it should be fine
  std::vector<LiDARPoints_t*> line_two;    // same above
  std::vector<LiDARPoints_t*> line_three;  // same
  std::vector<LiDARPoints_t*> line_four;   // same above
} TagBoundaries_t;

typedef struct Homogeneous
{
  // EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  float roll;
  float pitch;
  float yaw;
  Eigen::Matrix<float, 3, 1, Eigen::DontAlign> translation;
  Eigen::Matrix<float, 3, 3, Eigen::DontAlign> rotation;
  Eigen::Matrix<float, 4, 4, Eigen::DontAlign> homogeneous;
} Homogeneous_t;

typedef struct Grid
{
  float cx;
  float cz;
  float cy;
} Grid_t;

typedef struct RKHSDecoding
{
  Eigen::MatrixXf initial_template_points;
  Eigen::MatrixXf template_points;
  Eigen::MatrixXf template_points_xyz;
  Eigen::VectorXf template_points_feat;
  Eigen::MatrixXf template_points_3d;
  Eigen::MatrixXf* associated_pattern_3d;
  std::vector<float> score;
  int num_points;
  int size_num;
  int rotation_angle;
  double ell;
  double ave_intensity;
  int id;
  float id_score;
} RKHSDecoding_t;

typedef struct ClusterFamily
{
  // EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  int cluster_id;
  int valid;
  int top_ring;
  int bottom_ring;
  PointXYZIRT top_most_point;
  PointXYZIRT bottom_most_point;
  PointXYZIRT front_most_point;
  PointXYZIRT back_most_point;
  PointXYZIRT right_most_point;
  PointXYZIRT left_most_point;
  PointXYZIRT average;                  // Average point
  PointXYZIRT max_intensity;            // Maximux intensity point
  PointXYZIRT min_intensity;            // Minimum intensity point
  pcl::PointCloud<LiDARPoints_t> data;  // data doesn't have edge points
  pcl::PointCloud<LiDARPoints_t> edge_points;
  pcl::PointCloud<LiDARPoints_t> transformed_edge_points;

  // If the first point of the ring is the cluster.
  // If so, the the indices fo the two sides will be far away
  int special_case;
  Eigen::MatrixXf merged_data;    // this includes edge and filled-in points
  Eigen::MatrixXf merged_data_h;  // this includes edge and filled-in points

  std::vector<MaxMin_t> max_min_index_of_each_ring;             // to fill in points between end points in this cluster
  std::vector<std::vector<LiDARPoints_t*>> ordered_points_ptr;  // of the cluster (to find black margin of the tag)
  std::vector<double> accumulate_intensity_of_each_ring;        // to find the upper/lower lines of the tag
  TagLines_t tag_edges;                                         // store line segment from points
  TagBoundaries_t tag_boundaries;

  std::vector<LiDARPoints_t*> payload_right_boundary_ptr;  // of the cluster (to find black margin of the tag)
  std::vector<LiDARPoints_t*> payload_left_boundary_ptr;   // of the cluster (to find black margin of the tag)
  std::vector<LiDARPoints_t*> payload_boundary_ptr;        // of the cluster (to find black margin of the tag)
  int data_inliers;
  int edge_inliers;
  int inliers;
  double percentages_inliers;
  int boundary_pts;
  int boundary_rings;
  pcl::PointCloud<LiDARPoints_t*> payload;        // payload points with boundary
  pcl::PointCloud<LiDARPoints_t*> RLHS_decoding;  // payload points transformed
  int payload_without_boundary;                   // size of payload points without boundary
  double tag_size;
  double box_width;

  pcl::PointCloud<LiDARPoints_t> edge_group1;
  pcl::PointCloud<LiDARPoints_t> edge_group2;
  pcl::PointCloud<LiDARPoints_t> edge_group3;
  pcl::PointCloud<LiDARPoints_t> edge_group4;

  Eigen::Matrix<float, 3, 1, Eigen::DontAlign> normal_vector;
  Eigen::Matrix<float, 3, 3, Eigen::DontAlign> principal_axes;
  QuickDecodeEntry_t entry;
  Homogeneous_t pose_tag_to_lidar;
  Homogeneous_t pose;
  Homogeneous_t initial_pose;
  tf::Transform transform;
  RKHSDecoding_t rkhs_decoding;

  /* VectorXf:
   *          point_on_line.x : the X coordinate of a point on the line
   *          point_on_line.y : the Y coordinate of a point on the line
   *          point_on_line.z : the Z coordinate of a point on the line
   *          line_direction.x : the X coordinate of a line's direction
   *          line_direction.y : the Y coordinate of a line's direction
   *          line_direction.z : the Z coordinate of a line's direction
   */
  std::vector<Eigen::VectorXf> line_coeff;  // Upper, left, bottom, right line (count-clockwise)
} ClusterFamily_t;

typedef struct GrizTagFamily
{
  // How many codes are there in this tag family?
  uint32_t ncodes;

  // The codes in the family.
  uint64_t* codes;

  // how wide (in bit-sizes) is the black border? (usually 1)
  uint32_t black_border;

  // how many bits tall and wide is it? (e.g. 36bit tag ==> 6)
  uint32_t d;

  // minimum hamming distance between any two codes. (e.g. 36h11 => 11)
  uint32_t h;

  // a human-readable name, e.g., "tag36h11"
  char* name;

  // some detector implementations may preprocess codes in order to
  // accelerate decoding.  They put their data here. (Do not use the
  // same apriltag_family instance in more than one implementation)
  void* impl;
} GrizTagFamily_t;

typedef struct ClusterRemoval
{
  int minimum_return;
  int maximum_return;
  int plane_fitting;
  int plane_outliers;
  int boundary_point_check;
  int minimum_ring_points;
  int no_edge_check;
  int line_fitting;
  int pose_optimization;
  int decoding_failure;

  // for weighted gaussian
  int decoder_not_return;
  int decoder_fail_corner;
} ClusterRemoval_t;

typedef struct Statistics
{
  ClusterRemoval_t cluster_removal;
  int original_cluster_size;
  int remaining_cluster_size;
  int point_cloud_size;
  int edge_cloud_size;
} Statistics_t;

typedef struct Timing
{
  // in ms
  std::chrono::steady_clock::time_point start_total_time;
  std::chrono::steady_clock::time_point start_computation_time;
  std::chrono::steady_clock::time_point timing;

  double duration;
  double total_time;
  double total_duration;
  double edging_and_clustering_time;
  double to_pcl_vector_time;
  double fill_in_time;
  double point_check_time;
  double plane_fitting_removal_time;
  double line_fitting_time;
  double organize_points_time;
  double pca_time;
  double split_edge_time;
  double pose_optimization_time;
  double store_template_time;
  double payload_decoding_time;
  double normal_vector_time;
  double tag_to_robot_time;
} Timing_t;

typedef struct TimeDecoding
{
  std::chrono::steady_clock::time_point timing;  // ms
  double original;
  double matrix;
  double vectorization;
  double tbb_original;
  double tbb_vectorization;
  double manual_scheduling_tbb_vectorization;
  double tbb_scheduling_tbb_vectorization;
  double tbb_kd_tree;
} TimeDecoding_t;

typedef struct TestCluster
{
  int flag;
  ClusterFamily_t new_cluster;
} TestCluster_t;

typedef struct Debug
{
  std::vector<ClusterFamily_t*> point_check;
  std::vector<ClusterFamily_t*> boundary_point;
  std::vector<ClusterFamily_t*> no_edge;
  std::vector<ClusterFamily_t*> extract_payload;
} Debug_t;

typedef struct PathLeafString
{
  std::string operator()(const boost::filesystem::directory_entry& entry) const
  {
    return entry.path().leaf().string();
  }
} PathLeafString_t;

typedef nanoflann::KDTreeEigenMatrixAdaptor<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>, -1,
                                            nanoflann::metric_L2, false>
    kd_tree_t;

typedef Eigen::Triplet<float> Trip_t;

}  // namespace BipedLab
