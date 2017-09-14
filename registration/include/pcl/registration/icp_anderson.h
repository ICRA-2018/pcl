/*
 * Software License Agreement (BSD License)
 *
 *  Point Cloud Library (PCL) - www.pointclouds.org
 *  Copyright (c) 2010, Willow Garage, Inc.
 *  Copyright (c) 2012-, Open Perception, Inc.
 *
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions
 *  are met:
 *
 *   * Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *   * Redistributions in binary form must reproduce the above
 *     copyright notice, this list of conditions and the following
 *     disclaimer in the documentation and/or other materials provided
 *     with the distribution.
 *   * Neither the name of the copyright holder(s) nor the names of its
 *     contributors may be used to endorse or promote products derived
 *     from this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 *  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 *  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 *  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 *  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 *  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 *  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 *  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 *
 * $Id$
 *
 */

#ifndef PCL_ICP_ANDERSON_H_
#define PCL_ICP_ANDERSON_H_

#include <pcl/registration/icp.h>

namespace pcl
{
  /** \brief TODO.
    * \ingroup registration
    */
  template <typename PointSource, typename PointTarget, typename Scalar = float>
  class AndersonIterativeClosestPoint : public IterativeClosestPoint<PointSource, PointTarget>
  {
    public:
      using IterativeClosestPoint<PointSource, PointTarget>::reg_name_;
      using IterativeClosestPoint<PointSource, PointTarget>::getClassName;
      using IterativeClosestPoint<PointSource, PointTarget>::indices_;
      using IterativeClosestPoint<PointSource, PointTarget>::target_;
      using IterativeClosestPoint<PointSource, PointTarget>::input_;
      using IterativeClosestPoint<PointSource, PointTarget>::tree_;
      using IterativeClosestPoint<PointSource, PointTarget>::tree_reciprocal_;
      using IterativeClosestPoint<PointSource, PointTarget>::nr_iterations_;
      using IterativeClosestPoint<PointSource, PointTarget>::max_iterations_;
      using IterativeClosestPoint<PointSource, PointTarget>::previous_transformation_;
      using IterativeClosestPoint<PointSource, PointTarget>::final_transformation_;
      using IterativeClosestPoint<PointSource, PointTarget>::transformation_;
      using IterativeClosestPoint<PointSource, PointTarget>::transformation_epsilon_;
      using IterativeClosestPoint<PointSource, PointTarget>::converged_;
      using IterativeClosestPoint<PointSource, PointTarget>::corr_dist_threshold_;
      using IterativeClosestPoint<PointSource, PointTarget>::inlier_threshold_;
      using IterativeClosestPoint<PointSource, PointTarget>::min_number_correspondences_;
      using IterativeClosestPoint<PointSource, PointTarget>::update_visualizer_;

      using IterativeClosestPoint<PointSource, PointTarget>::need_target_blob_;
      using IterativeClosestPoint<PointSource, PointTarget>::correspondence_estimation_;
      using IterativeClosestPoint<PointSource, PointTarget>::target_has_normals_;
      using IterativeClosestPoint<PointSource, PointTarget>::convergence_criteria_;
      using IterativeClosestPoint<PointSource, PointTarget>::euclidean_fitness_epsilon_;
      using IterativeClosestPoint<PointSource, PointTarget>::transformation_rotation_epsilon_;
      using IterativeClosestPoint<PointSource, PointTarget>::need_source_blob_;
      using IterativeClosestPoint<PointSource, PointTarget>::use_reciprocal_correspondence_;
      using IterativeClosestPoint<PointSource, PointTarget>::correspondences_;
      using IterativeClosestPoint<PointSource, PointTarget>::correspondence_rejectors_;
      using IterativeClosestPoint<PointSource, PointTarget>::source_has_normals_;
      using IterativeClosestPoint<PointSource, PointTarget>::transformation_estimation_;
      using IterativeClosestPoint<PointSource, PointTarget>::determineRequiredBlobData;

      using IterativeClosestPoint<PointSource, PointTarget>::transformCloud;

      typedef Eigen::Matrix<Scalar, 3, 1> Vector3;
      typedef Eigen::Matrix<Scalar, 6, 1> Vector6;
      typedef Eigen::Matrix<Scalar, 3, 3> Matrix3;
      typedef Eigen::Matrix<Scalar, 4, 4> Matrix4;
      typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> VectorX;
      typedef Eigen::Matrix<Scalar, 6, Eigen::Dynamic> Matrix6X;

      typedef pcl::PointCloud<PointSource> PointCloudSource;
      typedef typename PointCloudSource::Ptr PointCloudSourcePtr;
      typedef typename PointCloudSource::ConstPtr PointCloudSourceConstPtr;

      /** \brief Empty constructor. */

      AndersonIterativeClosestPoint () 
        : alpha_limit_min_(-10)
        , alpha_limit_max_(10)
        , beta_(1.0)
        , small_step_threshold_(3)
        , error_overflow_threshold_(0.05)
      {
      }
      
      /** \brief Rigid transformation computation method  with initial guess.
        * \param output the transformed input point cloud dataset using the rigid transformation found
        * \param guess the initial guess of the transformation to compute
        */
      void 
      computeTransformation (PointCloudSource &output, const Eigen::Matrix4f &guess);

      void
      setAlphaThreshold (double min, double max) { alpha_limit_min_ = min; alpha_limit_max_ = max; }
      void
      getAlphaThreshold (double *min, double *max) { *min = alpha_limit_min_; *max = alpha_limit_max_; }

      void
      setBeta (double beta) { beta_ = beta; }
      double
      getBeta () { return beta_; }

      void
      setSmallStepThreshold (int n) { small_step_threshold_ = n; }
      int
      getSmallStepThreshold () { return small_step_threshold_; }

      void
      setErrorOverflowThreshold (double e) { error_overflow_threshold_ = e; }
      double
      getErrorOverflowThreshold () { return error_overflow_threshold_; }

      double alpha_limit_min_;
      double alpha_limit_max_;
      double beta_;
      int small_step_threshold_;
      double error_overflow_threshold_;

    private:
      inline double 
      getMSE (Matrix4 transformation_);

      typename pcl::AndersonIterativeClosestPoint<PointSource, PointTarget>::Vector6
      Matrix42Vector6 (const Matrix4 m);

      typename pcl::AndersonIterativeClosestPoint<PointSource, PointTarget>::Matrix4
      Vector62Matrix4 (const Vector6 v);

      typename pcl::AndersonIterativeClosestPoint<PointSource, PointTarget>::VectorX
      get_alphas_lstsq (const Matrix6X f);

      int
      alphas_cond (const VectorX alphas);

      typename pcl::AndersonIterativeClosestPoint<PointSource, PointTarget>::Vector6
      get_next_u (const Matrix6X u, const Matrix6X g, const Matrix6X f);

  };
}

#include <pcl/registration/impl/icp_anderson.hpp>

#endif  //#ifndef PCL_ICP_ANDERSON_H_
